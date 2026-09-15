#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 120
# Quotient-lattice congruence solve
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20
R10 = r11 * r20
R01 = r10 * r21
R11 = r11 * r21

RANDOM_SEED = 120


# ============================================================
# PRIME TEST
# ============================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    )

    for p in small:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in (
        2, 3, 5, 7,
        11, 13, 17
    ):

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(bits):

    while True:

        p = random.getrandbits(bits)

        p |= 1 << (bits - 1)
        p |= 1

        if is_probable_prime(p):
            return p


def make_semiprime(bits):

    fb = bits // 2

    while True:

        p = random_prime(fb)
        q = random_prime(fb)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n.bit_length() >= bits - 1:
            return n, p, q


# ============================================================
# CARRY
# ============================================================

def carry_E(r1, r2, k, l, a, b):

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    return c1 + c2 + c3


# ============================================================
# CROSS-SCALE QUOTIENT RESIDUES
# ============================================================

def quotient_residue_states():

    """
    Enumerate all possible residue states:

        a0 = p mod 3
        a1 = p mod 29
        b0 = q mod 5
        b1 = q mod 31

    and derive:

        k0 = p//3  mod 29
        l0 = q//5  mod 31.

    No n is used here.
    """

    states = []

    inv3_mod29 = pow(
        3,
        -1,
        29
    )

    inv5_mod31 = pow(
        5,
        -1,
        31
    )

    for a0 in range(r10):

        for a1 in range(r11):

            kr = (
                (a1 - a0)
                * inv3_mod29
            ) % r11

            for b0 in range(r20):

                for b1 in range(r21):

                    lr = (
                        (b1 - b0)
                        * inv5_mod31
                    ) % r21

                    states.append(
                        {
                            "a0": a0,
                            "a1": a1,
                            "b0": b0,
                            "b1": b1,
                            "kr": kr,
                            "lr": lr,
                        }
                    )

    return states


# ============================================================
# E BOUND
# ============================================================

def E_upper_bound(k, l):

    """
    Exact carry structure gives

        c1 <= k-1
        c2 <= l-1

    and c3 is small.

    We use a safe bound of:

        E <= k + l + 10.
    """

    return k + l + 10


# ============================================================
# SOLVE QUOTIENT LATTICE
# ============================================================

def solve_state(n, state):

    """
    Let

        k = kr + 29u
        l = lr + 31v.

    We need

        Q00 - E = k*l.

    Instead of scanning all k/l pairs, use the fact that

        E <= k+l+10.

    Therefore:

        Q00 - (k+l+10) <= k*l <= Q00.

    Rearranging:

        (k+1)(l+1) <= Q00 + 11.

    This gives a finite hyperbolic region.

    The experiment enumerates the smaller u dimension and
    solves the resulting quadratic inequality for v.

    It deliberately does NOT assume k ~= l.
    """

    Q = n // R00

    kr = state["kr"]
    lr = state["lr"]

    # --------------------------------------------------------
    # Maximum k from p <= sqrt(n).
    # --------------------------------------------------------

    p_limit = math.isqrt(n)

    k_max = p_limit // r10

    if kr > k_max:
        return []

    # First valid k.

    k_first = kr

    if k_first == 0:
        k_first = r11

    if k_first > k_max:
        return []

    # --------------------------------------------------------
    # Rather than enumerate every l, use:
    #
    # k*l <= Q
    #
    # and
    #
    # k+l <= Q-k*l+10
    #
    # to get a narrow interval.
    #
    # We still scan u, but importantly this is one-dimensional.
    # --------------------------------------------------------

    results = []

    k = k_first

    while k <= k_max:

        # ----------------------------------------------------
        # Upper bound from k*l <= Q.
        # ----------------------------------------------------

        l_max = Q // k

        if l_max <= 0:
            break

        # ----------------------------------------------------
        # Lower bound from
        #
        # Q - k*l <= k+l+10
        #
        # => Q-k-10 <= l*(k+1)
        # ----------------------------------------------------

        numerator = (
            Q
            - k
            - 10
        )

        if numerator <= 0:
            l_min = 1
        else:
            l_min = (
                numerator
                + k
            ) // (
                k + 1
            )

        # ----------------------------------------------------
        # Enforce l == lr mod 31.
        # ----------------------------------------------------

        if lr == 0:

            first_l = 31

            if l_min <= first_l:
                l = first_l
            else:
                l = (
                    l_min
                    + (
                        first_l - l_min
                    ) % 31
                )

        else:

            l = (
                l_min
                + (
                    lr - l_min
                ) % 31
            )

        # ----------------------------------------------------
        # There should normally be very few l values.
        # ----------------------------------------------------

        while l <= l_max:

            E = Q - k * l

            # Exact carry bound.

            if 0 <= E <= E_upper_bound(k, l):

                a0 = state["a0"]
                b0 = state["b0"]

                p = (
                    r10 * k
                    + a0
                )

                q = (
                    r20 * l
                    + b0
                )

                # Verify cross-scale residues.

                if (
                    p % r11
                    == state["a1"]
                    and
                    q % r21
                    == state["b1"]
                ):

                    # Exact carry equation.

                    got = carry_E(
                        r10,
                        r20,
                        k,
                        l,
                        a0,
                        b0
                    )

                    if got == E:

                        results.append(
                            {
                                "p": p,
                                "q": q,
                                "k0": k,
                                "l0": l,
                                "a0": a0,
                                "a1": state["a1"],
                                "b0": b0,
                                "b1": state["b1"],
                                "E00": E,
                                "kr": kr,
                                "lr": lr,
                            }
                        )

            l += r21

        k += r11

    return results


# ============================================================
# FULL SEARCH
# ============================================================

def search(n):

    states = quotient_residue_states()

    state_tests = 0
    quotient_candidates = 0

    survivors = []
    exact = []

    start = time.perf_counter()

    for state in states:

        state_tests += 1

        results = solve_state(
            n,
            state
        )

        quotient_candidates += len(
            results
        )

        for result in results:

            survivors.append(
                result
            )

            if (
                result["p"]
                * result["q"]
                == n
            ):

                exact.append(
                    result
                )

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "states": len(states),
        "state_tests": state_tests,
        "quotient_candidates":
            quotient_candidates,
        "survivors":
            survivors,
        "exact":
            exact,
        "time":
            elapsed,
    }


# ============================================================
# DIRECT p-SCAN CONTROL
# ============================================================

def direct_scan(n):

    limit = math.isqrt(n)

    tested = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        tested += 1

        if n % p == 0:

            q = n // p

            return {
                "tested": tested,
                "p": p,
                "q": q,
                "time":
                    time.perf_counter()
                    - start,
            }

    return {
        "tested": tested,
        "p": None,
        "q": None,
        "time":
            time.perf_counter()
            - start,
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    result,
    direct
):

    print()
    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("GRID")

    print(
        "R00 =",
        R00
    )

    print(
        "R10 =",
        R10
    )

    print(
        "R01 =",
        R01
    )

    print(
        "R11 =",
        R11
    )

    print()
    print("DIRECT CONTROL")

    print(
        "p tested =",
        direct["tested"]
    )

    print(
        "runtime  = %.6f s"
        % direct["time"]
    )

    print()
    print("QUOTIENT-LATTICE SEARCH")

    print(
        "residue states       =",
        result["states"]
    )

    print(
        "state tests          =",
        result["state_tests"]
    )

    print(
        "quotient candidates  =",
        result["quotient_candidates"]
    )

    print(
        "survivors            =",
        len(result["survivors"])
    )

    print(
        "exact                =",
        len(result["exact"])
    )

    print(
        "runtime              = %.6f s"
        % result["time"]
    )

    print()
    print("SURVIVORS")

    for s in result["survivors"][:25]:

        print(
            "  "
            f"p={s['p']} "
            f"q={s['q']} | "
            f"k0={s['k0']} "
            f"l0={s['l0']} | "
            f"a=({s['a0']},{s['a1']}) "
            f"b=({s['b0']},{s['b1']}) | "
            f"E00={s['E00']} | "
            f"k0mod29={s['kr']} "
            f"l0mod31={s['lr']} | "
            f"exact={s['p']*s['q']==n}"
        )

    if len(result["survivors"]) > 25:

        print(
            "  ..."
            f"{len(result['survivors'])-25}"
            " more"
        )

    print()
    print(
        "TRUE FOUND =",
        any(
            s["p"] == true_p
            and
            s["q"] == true_q
            for s in result["survivors"]
        )
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(
        RANDOM_SEED
    )

    print("=" * 72)
    print("START EXPERIMENT 120")
    print("Quotient-lattice congruence solve")
    print("=" * 72)

    print()

    print(
        "r1 =",
        R1
    )

    print(
        "r2 =",
        R2
    )

    print(
        "R00 =",
        R00
    )

    print(
        "R10 =",
        R10
    )

    print(
        "R01 =",
        R01
    )

    print(
        "R11 =",
        R11
    )

    print()

    for bits in (
        30,
        36,
        42,
    ):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(
            bits
        )

        direct = direct_scan(
            n
        )

        result = search(
            n
        )

        report(
            n,
            p,
            q,
            result,
            direct
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 120")
    print("=" * 72)


if __name__ == "__main__":
    main()
