#!/usr/bin/env python3

"""
START EXPERIMENT 182

TWO-SIDED RESIDUE / QUOTIENT INTERVAL SEARCH

Experiment 181 showed that pure p-residue CRT does not help enough:

    p ≡ a mod r

still leaves approximately sqrt(n)/r candidates for each residue,
and enumerating all residues essentially recreates the same search.

This experiment uses BOTH factor sides.

For selected radices r,s choose:

    p ≡ a (mod r)
    q ≡ b (mod s)

Since pq=n:

    p ≡ n*b^{-1} (mod s)
    q ≡ n*a^{-1} (mod r)

Thus a single (a,b) pair gives:

    p ≡ P (mod r*s)
    q ≡ Q (mod r*s)

and both must satisfy:

    p*q = n
    p <= sqrt(n)
    q >= sqrt(n).

Instead of scanning every p in the progression, parameterize:

    p = P + M*k
    q = Q + M*l

with

    M = r*s.

Substitution gives:

    (P+Mk)(Q+Ml)=n.

Because n is fixed, k and l are tightly coupled.

The experiment derives the exact relation:

    M*P*l + M*Q*k + M^2*k*l
        =
    n - P*Q.

For a fixed k:

    l =
      (n - P*Q - M*Q*k)
      / (M*P + M^2*k)

whenever integral.

Symmetrically for fixed l.

We test whether these two-sided constraints significantly reduce the
candidate count compared with Experiment 181.

No C-values are used.

The hidden p,q are used only for benchmark verification.

"""


import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

MAX_PAIR_ASSIGNMENTS = 2_000_000

MAX_EXACT_TESTS = 2_000_000


# ============================================================
# NUMBER THEORY
# ============================================================

def is_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    ]

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = 41
    step = 2

    while d * d <= n:

        if n % d == 0:
            return False

        d += step
        step = 6 - step

    return True


def egcd(a, b):

    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(
        b,
        a % b,
    )

    return (
        g,
        y1,
        x1 - (a // b) * y1,
    )


def inv_mod(a, m):

    g, x, _ = egcd(
        a,
        m,
    )

    if g != 1:

        raise ValueError(
            f"{a} not invertible mod {m}"
        )

    return x % m


# ============================================================
# SEMIPRIME GENERATOR
# ============================================================

def make_semiprime(bits):

    low = 1 << (
        bits // 2 - 1
    )

    high = 1 << (
        bits // 2 + 1
    )

    for p in range(
        low | 1,
        high,
        2,
    ):

        if not is_prime(p):
            continue

        target = (
            (1 << bits)
            // p
        )

        for delta in range(
            -1000,
            1001,
            2,
        ):

            q = target + delta

            if q <= p:
                continue

            if not is_prime(q):
                continue

            n = p * q

            if n.bit_length() == bits:

                return p, q, n

    raise RuntimeError(
        "semiprime generation failed"
    )


# ============================================================
# CRT
# ============================================================

def crt_pair(
    a1,
    m1,
    a2,
    m2,
):

    g = math.gcd(
        m1,
        m2,
    )

    if (a2 - a1) % g != 0:
        return None

    m1r = m1 // g
    m2r = m2 // g

    rhs = (
        (a2 - a1)
        // g
    ) % m2r

    if m2r == 1:

        t = 0

    else:

        t = (
            rhs
            * inv_mod(
                m1r % m2r,
                m2r,
            )
        ) % m2r

    x = (
        a1
        + m1 * t
    )

    return (
        x % (m1 * m2r),
        m1 * m2r,
    )


# ============================================================
# TWO-SIDED RESIDUE CONSTRAINT
# ============================================================

def build_constraint(
    n,
    r,
    s,
    a,
    b,
):
    """
    Given:

        p ≡ a mod r
        q ≡ b mod s

    derive:

        p ≡ n*b^-1 mod s
        q ≡ n*a^-1 mod r

    Then build combined residue classes for p and q.
    """

    p_mod_s = (
        n
        * inv_mod(
            b,
            s,
        )
    ) % s

    q_mod_r = (
        n
        * inv_mod(
            a,
            r,
        )
    ) % r

    p_class = crt_pair(
        a,
        r,
        p_mod_s,
        s,
    )

    q_class = crt_pair(
        b,
        s,
        q_mod_r,
        r,
    )

    if p_class is None:
        return None

    if q_class is None:
        return None

    P, M1 = p_class
    Q, M2 = q_class

    if M1 != M2:

        raise RuntimeError(
            "expected equal moduli"
        )

    return P, Q, M1


# ============================================================
# DIRECT TWO-SIDED SEARCH
# ============================================================

def search_constraint(
    n,
    r,
    s,
    a,
    b,
):

    constraint = build_constraint(
        n,
        r,
        s,
        a,
        b,
    )

    if constraint is None:
        return None

    P, Q, M = constraint

    sqrt_n = math.isqrt(n)

    exact_tests = 0

    generated = 0

    start = time.perf_counter()

    # --------------------------------------------------------
    # Normalize P and Q to positive representatives.
    # --------------------------------------------------------

    if P == 0:
        P = M

    if Q == 0:
        Q = M

    # --------------------------------------------------------
    # p = P + M*k
    #
    # p <= sqrt(n)
    # --------------------------------------------------------

    k_max = (
        sqrt_n - P
    ) // M

    if k_max < 0:

        return {
            "found": None,
            "generated": 0,
            "exact_tests": 0,
            "time":
                time.perf_counter()
                - start,
        }

    # --------------------------------------------------------
    # q = Q + M*l
    #
    # Since p <= sqrt(n) <= q for the smaller factor:
    #
    # q >= ceil(n / sqrt(n)).
    # --------------------------------------------------------

    min_q = (
        n + sqrt_n - 1
    ) // sqrt_n

    if Q < min_q:

        l_min = (
            min_q - Q + M - 1
        ) // M

    else:

        l_min = 0

    # --------------------------------------------------------
    # Exact bilinear equation:
    #
    # (P+Mk)(Q+Ml)=n
    #
    # For fixed k:
    #
    # numerator =
    #     n - P*Q - M*Q*k
    #
    # denominator =
    #     M*P + M^2*k
    #
    # l = numerator / denominator
    # --------------------------------------------------------

    for k in range(
        k_max + 1
    ):

        p_candidate = (
            P + M * k
        )

        if p_candidate < 2:
            continue

        if p_candidate > sqrt_n:
            break

        numerator = (
            n
            - P * Q
            - M * Q * k
        )

        denominator = (
            M * P
            + M * M * k
        )

        if numerator < 0:
            continue

        if denominator <= 0:
            continue

        if numerator % denominator != 0:
            continue

        l = (
            numerator
            // denominator
        )

        if l < l_min:
            continue

        q_candidate = (
            Q + M * l
        )

        if q_candidate < min_q:
            continue

        generated += 1

        exact_tests += 1

        if exact_tests > MAX_EXACT_TESTS:

            return {
                "found": None,
                "generated": generated,
                "exact_tests":
                    exact_tests,
                "time":
                    time.perf_counter()
                    - start,
                "aborted":
                    "MAX_EXACT_TESTS",
            }

        if (
            p_candidate
            * q_candidate
            == n
        ):

            return {
                "found": (
                    p_candidate,
                    q_candidate,
                ),
                "generated": generated,
                "exact_tests":
                    exact_tests,
                "time":
                    time.perf_counter()
                    - start,
                "aborted": None,
            }

    return {
        "found": None,
        "generated": generated,
        "exact_tests":
            exact_tests,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# ENUMERATE ALL (a,b) FOR A RADIX PAIR
# ============================================================

def search_radix_pair(
    n,
    r,
    s,
):

    start = time.perf_counter()

    tested_pairs = 0

    exact_tests = 0

    for a in range(
        1,
        r,
    ):

        for b in range(
            1,
            s,
        ):

            tested_pairs += 1

            if (
                tested_pairs
                > MAX_PAIR_ASSIGNMENTS
            ):

                return {
                    "found": None,
                    "tested_pairs":
                        tested_pairs,
                    "exact_tests":
                        exact_tests,
                    "time":
                        time.perf_counter()
                        - start,
                    "aborted":
                        "MAX_PAIR_ASSIGNMENTS",
                }

            result = search_constraint(
                n,
                r,
                s,
                a,
                b,
            )

            exact_tests += (
                result["exact_tests"]
            )

            if result.get(
                "found"
            ):

                return {
                    "found":
                        result["found"],
                    "tested_pairs":
                        tested_pairs,
                    "exact_tests":
                        exact_tests,
                    "time":
                        time.perf_counter()
                        - start,
                    "aborted":
                        None,
                }

            if result.get(
                "aborted"
            ):

                return {
                    "found": None,
                    "tested_pairs":
                        tested_pairs,
                    "exact_tests":
                        exact_tests,
                    "time":
                        time.perf_counter()
                        - start,
                    "aborted":
                        result["aborted"],
                }

    return {
        "found": None,
        "tested_pairs":
            tested_pairs,
        "exact_tests":
            exact_tests,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# SINGLE-SIDE BASELINE
# ============================================================

def single_side_candidate_count(
    n,
    r,
):

    sqrt_n = math.isqrt(n)

    total = 0

    for a in range(
        1,
        r,
    ):

        P = a

        if P > sqrt_n:
            continue

        total += (
            (sqrt_n - P)
            // r
            + 1
        )

    return total


# ============================================================
# TRUE RESIDUES
# ============================================================

def true_residue_pair(
    p,
    q,
    r,
    s,
):

    return (
        p % r,
        q % s,
    )


# ============================================================
# RUN INSTANCE
# ============================================================

def run_instance(bits):

    print()
    print("=" * 72)
    print(
        f"START INSTANCE {bits}-BIT"
    )
    print("=" * 72)

    p, q, n = make_semiprime(
        bits
    )

    sqrt_n = math.isqrt(n)

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"sqrt(n) = {sqrt_n}"
    )

    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    # --------------------------------------------------------
    # True residue pairs.
    # --------------------------------------------------------

    print()
    print(
        "TRUE RESIDUES"
    )

    for r, s in [
        (17, 19),
        (17, 37),
        (43, 19),
        (43, 37),
        (59, 61),
    ]:

        a, b = true_residue_pair(
            p,
            q,
            r,
            s,
        )

        print(
            f"  ({r},{s}) "
            f"= ({a},{b})"
        )

    # --------------------------------------------------------
    # Baseline.
    # --------------------------------------------------------

    print()
    print(
        "ONE-SIDED CANDIDATE COUNTS"
    )

    for r in [
        17,
        43,
        59,
        71,
        83,
    ]:

        count = (
            single_side_candidate_count(
                n,
                r,
            )
        )

        print(
            f"  r={r} "
            f"candidates={count}"
        )

    # --------------------------------------------------------
    # Two-sided pair searches.
    # --------------------------------------------------------

    print()
    print(
        "TWO-SIDED RESIDUE SEARCH"
    )

    radix_pairs = [
        (17, 19),
        (17, 37),
        (17, 61),
        (17, 73),
        (17, 89),
        (43, 19),
        (43, 37),
        (43, 61),
        (59, 19),
        (59, 37),
        (59, 61),
    ]

    for r, s in radix_pairs:

        a, b = true_residue_pair(
            p,
            q,
            r,
            s,
        )

        result = search_radix_pair(
            n,
            r,
            s,
        )

        M = r * s

        print(
            f"  ({r},{s}) "
            f"M={M} "
            f"true=({a},{b}) "
            f"pairs={result['tested_pairs']} "
            f"exact={result['exact_tests']} "
            f"time={result['time']:.6f}s"
        )

        if result.get(
            "found"
        ):

            fp, fq = result[
                "found"
            ]

            print()
            print(
                "*** FACTOR FOUND ***"
            )

            print(
                f"    p = {fp}"
            )

            print(
                f"    q = {fq}"
            )

            print(
                f"    correct = "
                f"{fp * fq == n}"
            )

            print()
            print(
                f"FINISHED INSTANCE "
                f"{bits}-BIT"
            )

            return

        if result.get(
            "aborted"
        ):

            print(
                f"    status="
                f"{result['aborted']}"
            )

    print()
    print(
        "No factor recovered."
    )

    print()
    print(
        f"FINISHED INSTANCE {bits}-BIT"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print(
        "START EXPERIMENT 182"
    )

    print()

    print(
        f"R1 = {R1}"
    )

    print(
        f"R2 = {R2}"
    )

    print(
        f"BITS = {BITS}"
    )

    print(
        f"MAX_PAIR_ASSIGNMENTS = "
        f"{MAX_PAIR_ASSIGNMENTS}"
    )

    print(
        f"MAX_EXACT_TESTS = "
        f"{MAX_EXACT_TESTS}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 182"
    )


if __name__ == "__main__":
    main()
