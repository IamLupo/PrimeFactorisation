#!/usr/bin/env python3

import sympy as sp


# ============================================================
# EXPERIMENT: SCALED FACTOR PARAMETERIZATION
#
# p = 4y - 4x + 3
# q = 4y + 4x - 3
#
# Therefore:
#
# pq = 16y^2 - 16x^2 + 24x - 9
#
# ============================================================

LIMIT_X = 50
LIMIT_Y = 50


def experiment():
    print("=" * 80)
    print("EXPERIMENT: MOD-8 SCALED PARAMETERIZATION")
    print("=" * 80)

    print()
    print("PARAMETERIZATION")
    print("  p = 4y - 4x + 3")
    print("  q = 4y + 4x - 3")
    print()
    print("POLYNOMIAL")
    print("  n = 16y^2 - 16x^2 + 24x - 9")
    print()

    counts = {
        1: 0,
        3: 0,
        5: 0,
        7: 0,
    }

    prime_cases = []
    failures = []

    for x in range(-LIMIT_X, LIMIT_X + 1):
        for y in range(-LIMIT_Y, LIMIT_Y + 1):

            p = 4 * y - 4 * x + 3
            q = 4 * y + 4 * x - 3

            n_factor = p * q

            n_formula = (
                16 * y * y
                - 16 * x * x
                + 24 * x
                - 9
            )

            # ------------------------------------------------
            # Identity check
            # ------------------------------------------------

            if n_factor != n_formula:
                failures.append(
                    ("IDENTITY", x, y, p, q, n_factor, n_formula)
                )
                continue

            # ------------------------------------------------
            # Mod 8
            # ------------------------------------------------

            r = n_factor % 8

            if r in counts:
                counts[r] += 1

            # ------------------------------------------------
            # Prime factor cases
            # ------------------------------------------------

            if p > 1 and q > 1 and sp.isprime(p) and sp.isprime(q):
                prime_cases.append(
                    (x, y, n_factor, p, q, p % 8, q % 8)
                )

    # ========================================================
    # RESULTS
    # ========================================================

    print("IDENTITY CHECK")
    print(f"  failures = {len(failures)}")
    print(f"  PASS     = {len(failures) == 0}")
    print()

    print("MOD-8 DISTRIBUTION")
    for r in [1, 3, 5, 7]:
        print(f"  n ≡ {r} (mod 8): {counts[r]}")
    print()

    print("THEORETICAL FACTOR RESIDUES")
    print("  p = 4(y-x) + 3")
    print("  q = 4(y+x) - 3")
    print()
    print("  Therefore:")
    print("    p ≡ 3 or 7 (mod 8)")
    print("    q ≡ 1 or 5 (mod 8)")
    print()

    print("POSSIBLE PRIME FACTOR RESIDUE PAIRS")
    residue_pairs = {}

    for x, y, n, p, q, rp, rq in prime_cases:
        key = (rp, rq)
        residue_pairs[key] = residue_pairs.get(key, 0) + 1

    for key, count in sorted(residue_pairs.items()):
        rp, rq = key
        print(
            f"  p ≡ {rp} (mod 8), "
            f"q ≡ {rq} (mod 8): {count}"
        )

    print()

    print("FIRST PRIME × PRIME CASES")
    print("-" * 80)

    for case in prime_cases[:30]:
        x, y, n, p, q, rp, rq = case

        print(
            f"x={x:3d} "
            f"y={y:3d} "
            f"-> p={p:5d} "
            f"q={q:5d} "
            f"n={n:8d} "
            f"| p%8={rp} "
            f"q%8={rq} "
            f"n%8={n % 8}"
        )

    print()
    print("=" * 80)

    # ========================================================
    # PARITY ANALYSIS
    # ========================================================

    print()
    print("PARITY ANALYSIS")
    print("-" * 80)

    parity_counts = {}

    for x in range(-LIMIT_X, LIMIT_X + 1):
        for y in range(-LIMIT_Y, LIMIT_Y + 1):

            p = 4 * y - 4 * x + 3
            q = 4 * y + 4 * x - 3
            n = p * q

            key = (
                x % 2,
                y % 2,
                p % 8,
                q % 8,
                n % 8,
            )

            parity_counts[key] = parity_counts.get(key, 0) + 1

    print(
        "  x%2  y%2 | p%8  q%8 | n%8"
    )

    seen = set()

    for key in sorted(parity_counts):
        xpar, ypar, rp, rq, rn = key

        summary = (xpar, ypar, rp, rq, rn)

        if summary not in seen:
            seen.add(summary)

            print(
                f"   {xpar}    {ypar}  |  "
                f" {rp}    {rq}  |  {rn}"
            )

    print()
    print("=" * 80)


if __name__ == "__main__":
    experiment()
