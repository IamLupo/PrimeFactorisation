#!/usr/bin/env python3

"""
START EXPERIMENT 181

MULTI-RADIX RESIDUE -> RATIONAL APPROXIMATION

Experiment 180 established that the C-derived relations are tautological
when C is calculated from a guessed residue vector.

Therefore we now discard C completely.

The only useful information is:

    p ≡ a_i (mod r_i)

for small prime radices r_i.

Since

    q = n / p

and p,q are balanced, q is approximately sqrt(n).

If

    p = a + r*k

then

    n = p*q
      = (a+r*k)q.

For a fixed residue a:

    q ≡ n*a^{-1} (mod r).

Thus every candidate residue a gives simultaneous congruence information
for p and q.

Instead of enumerating all A/B vectors, this experiment reconstructs
candidate p values from a SINGLE residue using a quotient parameter:

    p = a + r*k

with

    k ≈ sqrt(n)/r.

We then use the opposite congruence:

    q = q0 + r*t

and the exact equation

    (a+r*k)(q0+r*t)=n.

This gives a bilinear Diophantine equation.

The experiment measures how many k values survive when the two-sided
congruence is enforced.

A second stage combines two independent radices using generalized CRT.

No C-values are used.

No hidden residues are used during the search.

The hidden p,q are used only for verification.

"""


import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

MAX_EXACT_TESTS = 2_000_000

# Maximum number of residue assignments examined per instance.
MAX_RESIDUE_COMBINATIONS = 2_000_000

# Trial division baseline cap.
MAX_TRIALS = 2_000_000


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
# GENERALIZED CRT
# ============================================================

def crt_pair_general(
    a1,
    m1,
    a2,
    m2,
):
    """
    Solve:

        x = a1 mod m1
        x = a2 mod m2

    Returns:

        (x, lcm(m1,m2))

    or None.
    """

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

    x = a1 + m1 * t

    modulus = m1 * m2r

    return (
        x % modulus,
        modulus,
    )


# ============================================================
# SEMIPRIME GENERATION
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
# SINGLE-RADIX RESIDUE SEARCH
# ============================================================

def search_single_radix(
    n,
    r,
):
    """
    Try every nonzero a mod r.

    For each:

        p = a mod r

    and therefore

        q = n*a^-1 mod r.

    Parameterize:

        p = a + r*k
        q = q0 + r*t

    and search the balanced region.

    This is NOT intended to be asymptotically better than trial
    division. The goal is to measure the constraint density.
    """

    sqrt_n = math.isqrt(n)

    exact_tests = 0

    candidate_p_count = 0

    surviving_residue_count = 0

    t0 = time.perf_counter()

    for a in range(
        1,
        r,
    ):

        q0 = (
            n
            * inv_mod(
                a,
                r,
            )
        ) % r

        # Smallest positive q in this class.
        if q0 == 0:
            q0 = r

        # p <= sqrt(n)
        #
        # p = a + r*k

        k_max = (
            sqrt_n - a
        ) // r

        if k_max < 0:
            continue

        # Balanced q condition.
        #
        # q is at least p for the first factor.
        #
        # We only need to inspect q around n/p.
        #
        # Instead of enumerating all q classes, use the exact
        # equation to calculate whether q=n/p is integral.

        local_candidates = 0

        for k in range(
            k_max + 1
        ):

            candidate_p = (
                a
                + r * k
            )

            if candidate_p < 2:
                continue

            candidate_p_count += 1

            if (
                candidate_p
                > sqrt_n
            ):
                break

            # q must satisfy its residue class.
            if (
                n % candidate_p
                != 0
            ):
                continue

            exact_tests += 1

            q = (
                n
                // candidate_p
            )

            if q % r != q0:
                raise RuntimeError(
                    "internal congruence error"
                )

            if (
                candidate_p * q
                == n
            ):

                return {
                    "found": (
                        candidate_p,
                        q,
                    ),
                    "candidate_p_count":
                        candidate_p_count,
                    "exact_tests":
                        exact_tests,
                    "surviving_residue_count":
                        surviving_residue_count + 1,
                    "time":
                        time.perf_counter()
                        - t0,
                }

        surviving_residue_count += (
            local_candidates > 0
        )

        if (
            exact_tests
            > MAX_EXACT_TESTS
        ):
            break

    return {
        "found": None,
        "candidate_p_count":
            candidate_p_count,
        "exact_tests":
            exact_tests,
        "surviving_residue_count":
            surviving_residue_count,
        "time":
            time.perf_counter()
            - t0,
    }


# ============================================================
# TWO-SIDED CRT SEARCH
# ============================================================

def search_two_radices(
    n,
    r1,
    r2,
):
    """
    Search all nonzero residue pairs

        a1 = p mod r1
        a2 = p mod r2

    and combine them directly.

    Since r1,r2 are coprime:

        p = P mod (r1*r2)

    This creates a sparse arithmetic progression of p candidates.

    """

    sqrt_n = math.isqrt(n)

    modulus = r1 * r2

    residue_classes = []

    for a1 in range(
        1,
        r1,
    ):

        for a2 in range(
            1,
            r2,
        ):

            merged = crt_pair_general(
                a1,
                r1,
                a2,
                r2,
            )

            if merged is None:
                continue

            P, M = merged

            residue_classes.append(
                P
            )

    candidate_classes = len(
        residue_classes
    )

    exact_tests = 0

    t0 = time.perf_counter()

    for P in residue_classes:

        if P == 0:
            P = modulus

        if P > sqrt_n:
            continue

        count = (
            sqrt_n - P
        ) // modulus

        for t in range(
            count + 1
        ):

            candidate_p = (
                P
                + t * modulus
            )

            if candidate_p < 2:
                continue

            exact_tests += 1

            if (
                exact_tests
                > MAX_EXACT_TESTS
            ):

                return {
                    "found": None,
                    "candidate_classes":
                        candidate_classes,
                    "exact_tests":
                        exact_tests,
                    "time":
                        time.perf_counter()
                        - t0,
                    "aborted":
                        "MAX_EXACT_TESTS",
                }

            if n % candidate_p == 0:

                q = (
                    n
                    // candidate_p
                )

                return {
                    "found": (
                        candidate_p,
                        q,
                    ),
                    "candidate_classes":
                        candidate_classes,
                    "exact_tests":
                        exact_tests,
                    "time":
                        time.perf_counter()
                        - t0,
                    "aborted": None,
                }

    return {
        "found": None,
        "candidate_classes":
            candidate_classes,
        "exact_tests":
            exact_tests,
        "time":
            time.perf_counter()
            - t0,
        "aborted": None,
    }


# ============================================================
# MULTI-RADIX CRT SEARCH
# ============================================================

def search_multi_radix(
    n,
    radices,
):
    """
    Combine residue classes of p over several radices.

    This is intentionally pure residue search:

        p ≡ a_i mod r_i

    with every nonzero residue considered.

    The goal is to measure the number of resulting p candidates
    below sqrt(n).
    """

    sqrt_n = math.isqrt(n)

    states = {
        (0, 1)
    }

    generated = 0

    for r in radices:

        new_states = set()

        for P, M in states:

            for a in range(
                1,
                r,
            ):

                merged = crt_pair_general(
                    P,
                    M,
                    a,
                    r,
                )

                if merged is None:
                    continue

                P2, M2 = merged

                if P2 == 0:
                    smallest = M2
                else:
                    smallest = P2

                if smallest > sqrt_n:
                    continue

                new_states.add(
                    (
                        P2,
                        M2,
                    )
                )

                generated += 1

                if (
                    len(new_states)
                    > MAX_RESIDUE_COMBINATIONS
                ):

                    return {
                        "found": None,
                        "states": len(
                            new_states
                        ),
                        "generated":
                            generated,
                        "exact_tests": 0,
                        "time": 0.0,
                        "aborted":
                            "MAX_RESIDUE_COMBINATIONS",
                    }

        states = new_states

        if not states:

            return {
                "found": None,
                "states": 0,
                "generated": generated,
                "exact_tests": 0,
                "time": 0.0,
                "aborted": None,
            }

    exact_tests = 0

    t0 = time.perf_counter()

    for P, M in states:

        if P == 0:
            P = M

        if P > sqrt_n:
            continue

        count = (
            sqrt_n - P
        ) // M

        for t in range(
            count + 1
        ):

            candidate_p = (
                P
                + t * M
            )

            if candidate_p < 2:
                continue

            exact_tests += 1

            if (
                exact_tests
                > MAX_EXACT_TESTS
            ):

                return {
                    "found": None,
                    "states": len(states),
                    "generated": generated,
                    "exact_tests":
                        exact_tests,
                    "time":
                        time.perf_counter()
                        - t0,
                    "aborted":
                        "MAX_EXACT_TESTS",
                }

            if n % candidate_p == 0:

                q = (
                    n
                    // candidate_p
                )

                return {
                    "found": (
                        candidate_p,
                        q,
                    ),
                    "states": len(states),
                    "generated": generated,
                    "exact_tests":
                        exact_tests,
                    "time":
                        time.perf_counter()
                        - t0,
                    "aborted": None,
                }

    return {
        "found": None,
        "states": len(states),
        "generated": generated,
        "exact_tests":
            exact_tests,
        "time":
            time.perf_counter()
            - t0,
        "aborted": None,
    }


# ============================================================
# TRUE RESIDUES
# ============================================================

def true_residue(
    p,
    q,
):
    return (
        [p % r for r in R1],
        [q % s for s in R2],
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

    # Hidden values only for verification.
    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    true_A, true_B = true_residue(
        p,
        q,
    )

    print(
        f"true A = {true_A}"
    )

    print(
        f"true B = {true_B}"
    )

    # --------------------------------------------------------
    # Single-radix baseline.
    # --------------------------------------------------------

    print()
    print(
        "SINGLE-RADIX TESTS"
    )

    for r in (
        R1[0],
        R1[1],
        R2[0],
        R2[1],
    ):

        if r in R1:

            true_a = p % r

        else:

            true_a = q % r

        result = search_single_radix(
            n,
            r,
        )

        print(
            f"  r={r} "
            f"true_residue={true_a} "
            f"candidate_p={result['candidate_p_count']} "
            f"exact={result['exact_tests']} "
            f"time={result['time']:.6f}s"
        )

        if result.get("found"):

            fp, fq = result["found"]

            print(
                f"    FOUND "
                f"({fp}, {fq})"
            )

    # --------------------------------------------------------
    # Two-radix combinations.
    # --------------------------------------------------------

    print()
    print(
        "TWO-RADIX P-CRT TESTS"
    )

    pairs = [
        (17, 43),
        (17, 59),
        (17, 71),
        (17, 83),
        (43, 59),
        (43, 71),
    ]

    for r1, r2 in pairs:

        result = search_two_radices(
            n,
            r1,
            r2,
        )

        print(
            f"  ({r1},{r2}) "
            f"M={r1*r2} "
            f"classes={result['candidate_classes']} "
            f"exact={result['exact_tests']} "
            f"time={result['time']:.6f}s"
        )

        if result.get("found"):

            fp, fq = result["found"]

            print(
                f"    FOUND "
                f"({fp}, {fq})"
            )

            print()
            print(
                f"FINISHED INSTANCE "
                f"{bits}-BIT"
            )

            return

        if result.get("aborted"):

            print(
                f"    status="
                f"{result['aborted']}"
            )

    # --------------------------------------------------------
    # Multi-radix tests.
    # --------------------------------------------------------

    print()
    print(
        "MULTI-RADIX TESTS"
    )

    radix_sets = [
        [17, 43],
        [17, 43, 59],
        [17, 43, 59, 71],
        [17, 43, 59, 71, 83],
    ]

    for radices in radix_sets:

        result = search_multi_radix(
            n,
            radices,
        )

        product = math.prod(
            radices
        )

        print(
            f"  {radices} "
            f"M={product} "
            f"states={result['states']} "
            f"generated={result['generated']} "
            f"exact={result['exact_tests']} "
            f"time={result['time']:.6f}s"
        )

        if result.get("found"):

            fp, fq = result["found"]

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

        if result.get("aborted"):

            print(
                f"    status="
                f"{result['aborted']}"
            )

    print()
    print(
        "No factor recovered by pure "
        "multi-radix residue search."
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
        "START EXPERIMENT 181"
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
        f"MAX_EXACT_TESTS = "
        f"{MAX_EXACT_TESTS}"
    )

    print(
        f"MAX_RESIDUE_COMBINATIONS = "
        f"{MAX_RESIDUE_COMBINATIONS}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 181"
    )


if __name__ == "__main__":
    main()
