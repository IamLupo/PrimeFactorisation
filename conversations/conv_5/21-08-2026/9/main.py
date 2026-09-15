#!/usr/bin/env python3

import sympy as sp
from collections import Counter, defaultdict


# ============================================================
# GENERAL 2^z AFFINE FACTOR EXPERIMENT
#
# k = 2^(z-1)
#
# p = k(y-x) + a
# q = k(y+x) - a
#
# pq = k^2(y^2-x^2) + 2akx - a^2
#
# Since k = 2^(z-1):
#
#   k^2 ≡ 0 (mod 2^z)
#   2ak = a * 2^z ≡ 0 (mod 2^z)
#
# therefore:
#
#   pq ≡ -a^2 (mod 2^z)
#
# ============================================================


# ------------------------------------------------------------
# SEARCH SETTINGS
# ------------------------------------------------------------

Z_VALUES = range(2, 9)       # z = 2 ... 8
A_VALUES = range(1, 16)      # test a = 1 ... 15

X_MIN = -50
X_MAX = 50

Y_MIN = -50
Y_MAX = 50

SHOW_PRIME_CASES = 20


# ------------------------------------------------------------
# THEORETICAL FORMULA
# ------------------------------------------------------------

def predicted_residue(z, a):
    modulus = 1 << z
    return (-a * a) % modulus


# ------------------------------------------------------------
# MAIN EXPERIMENT
# ------------------------------------------------------------

def run_experiment(z, a):

    modulus = 1 << z
    k = 1 << (z - 1)

    expected = predicted_residue(z, a)

    total = 0
    identity_failures = 0
    congruence_failures = 0

    residue_counts = Counter()
    factor_residue_counts = Counter()

    prime_cases = []

    parity_classes = defaultdict(Counter)

    for x in range(X_MIN, X_MAX + 1):
        for y in range(Y_MIN, Y_MAX + 1):

            # ------------------------------------------------
            # Parameterization
            # ------------------------------------------------

            p = k * (y - x) + a
            q = k * (y + x) - a

            # ------------------------------------------------
            # Exact polynomial
            # ------------------------------------------------

            pq = p * q

            polynomial = (
                k * k * (y * y - x * x)
                + 2 * a * k * x
                - a * a
            )

            total += 1

            # ------------------------------------------------
            # Identity check
            # ------------------------------------------------

            if pq != polynomial:
                identity_failures += 1

            # ------------------------------------------------
            # Mod 2^z check
            # ------------------------------------------------

            actual_residue = pq % modulus

            residue_counts[actual_residue] += 1

            if actual_residue != expected:
                congruence_failures += 1

            # ------------------------------------------------
            # Factor residues
            # ------------------------------------------------

            rp = p % modulus
            rq = q % modulus

            factor_residue_counts[(rp, rq)] += 1

            # ------------------------------------------------
            # Parity structure
            # ------------------------------------------------

            parity_classes[(x % 2, y % 2)][(rp, rq)] += 1

            # ------------------------------------------------
            # Prime × prime
            # ------------------------------------------------

            if (
                p > 1
                and q > 1
                and sp.isprime(p)
                and sp.isprime(q)
            ):
                prime_cases.append(
                    {
                        "x": x,
                        "y": y,
                        "p": p,
                        "q": q,
                        "n": pq,
                        "p_mod": rp,
                        "q_mod": rq,
                        "n_mod": actual_residue,
                    }
                )

    # ========================================================
    # REPORT
    # ========================================================

    print()
    print("=" * 90)
    print(f"z={z}, a={a}")
    print("=" * 90)

    print()
    print("PARAMETERS")
    print(f"  z        = {z}")
    print(f"  modulus  = 2^{z} = {modulus}")
    print(f"  k        = 2^{z-1} = {k}")
    print(f"  a        = {a}")

    print()
    print("FORMULAS")
    print(f"  p = {k}(y-x) + {a}")
    print(f"  q = {k}(y+x) - {a}")

    print()
    print("THEORETICAL CONGRUENCE")
    print(f"  pq ≡ -a² (mod 2^{z})")
    print(f"  -a² ≡ {expected} (mod {modulus})")

    print()
    print("CHECKS")
    print(f"  total cases          = {total}")
    print(f"  identity failures    = {identity_failures}")
    print(f"  congruence failures  = {congruence_failures}")
    print(
        f"  identity PASS        = "
        f"{identity_failures == 0}"
    )
    print(
        f"  congruence PASS      = "
        f"{congruence_failures == 0}"
    )

    print()
    print("OBSERVED n RESIDUES")
    for residue, count in sorted(residue_counts.items()):
        print(
            f"  n ≡ {residue:>3} (mod {modulus}) : "
            f"{count}"
        )

    print()
    print("FACTOR RESIDUE PAIRS")
    for (rp, rq), count in sorted(factor_residue_counts.items()):
        print(
            f"  p ≡ {rp:>3}, "
            f"q ≡ {rq:>3} : "
            f"{count}"
        )

    print()
    print("PARITY STRUCTURE")
    print("  x%2 y%2 | possible (p,q) residues")
    print("  -----------------------------------")

    for (xp, yp), pairs in sorted(parity_classes.items()):

        pair_text = ", ".join(
            f"({rp},{rq})={count}"
            for (rp, rq), count
            in sorted(pairs.items())
        )

        print(
            f"   {xp}   {yp}  | {pair_text}"
        )

    print()
    print(f"PRIME × PRIME CASES: {len(prime_cases)}")

    print("-" * 90)

    for case in prime_cases[:SHOW_PRIME_CASES]:

        print(
            f"x={case['x']:4d} "
            f"y={case['y']:4d} "
            f"p={case['p']:6d} "
            f"q={case['q']:6d} "
            f"n={case['n']:10d} "
            f"| "
            f"p%{modulus}={case['p_mod']:3d} "
            f"q%{modulus}={case['q_mod']:3d} "
            f"n%{modulus}={case['n_mod']:3d}"
        )


# ------------------------------------------------------------
# SUMMARY ACROSS ALL a
# ------------------------------------------------------------

def summary(z):

    modulus = 1 << z
    k = 1 << (z - 1)

    print()
    print("#" * 90)
    print(f"SUMMARY FOR z={z}")
    print(f"modulus = {modulus}")
    print(f"k       = {k}")
    print("#" * 90)

    print()
    print(
        f"{'a':>3} "
        f"{'a² mod 2^z':>12} "
        f"{'-a² mod 2^z':>14}"
    )

    print("-" * 40)

    for a in A_VALUES:

        a2 = (a * a) % modulus
        target = (-a * a) % modulus

        print(
            f"{a:3d} "
            f"{a2:12d} "
            f"{target:14d}"
        )


# ------------------------------------------------------------
# RUN
# ------------------------------------------------------------

def main():

    print("=" * 90)
    print("GENERALIZED 2^z PARAMETERIZATION EXPLORER")
    print("=" * 90)

    print()
    print("General family:")
    print()
    print("  k = 2^(z-1)")
    print("  p = k(y-x) + a")
    print("  q = k(y+x) - a")
    print()
    print("Therefore:")
    print()
    print("  pq = k²(y²-x²) + 2akx - a²")
    print()
    print("and modulo 2^z:")
    print()
    print("  pq ≡ -a² (mod 2^z)")
    print()

    # --------------------------------------------------------
    # First show residue summary
    # --------------------------------------------------------

    for z in Z_VALUES:
        summary(z)

    # --------------------------------------------------------
    # Detailed experiments
    # --------------------------------------------------------

    # Keep this small enough to read.
    # Change A_VALUES / Z_VALUES above for larger searches.

    for z in Z_VALUES:

        for a in A_VALUES:

            run_experiment(z, a)


if __name__ == "__main__":
    main()
