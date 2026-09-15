#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 147 — EXACT 2-ADIC UNIT-GAUGE / RANK-1 INVARIANCE AUDIT
==============================================================================

Experiment 146 showed:

    raw R:
        rank-one through 2^3;

    row-normalized R:
        rank-one only through 2^2;

    row+column-normalized R:
        rank-one already fails modulo 2.

The key mathematical distinction is that dividing by powers of 2 is NOT
a unit change modulo 2^e.

Experiment 147 therefore removes that ambiguity.

We test only diagonal transformations whose entries are odd:

    R' = D_L R D_R

with D_L and D_R invertible modulo 2^e.

Such transformations preserve determinantal rank over Z/2^e Z up to
equivalent row/column operations.

We investigate whether the observed 2-adic rank-one structure is stable
under several exact unit gauges:

    1. odd parts of row gcds;
    2. odd parts of column gcds;
    3. random deterministic odd diagonal units;
    4. combinations of the above.

For every gauge we test:

    rank mod 2,
    rank-one determinantal condition mod 4,
    mod 8,
    mod 16.

We also compare the first 2-adic defect rank.

The purpose is to distinguish:

    true congruence structure
        from
    artifacts introduced by non-invertible scaling.

No floating point.
No SymPy.
No extrapolation.
No connection to the original (p,q)-kernel is assumed.
"""

from __future__ import annotations

from fractions import Fraction
import sys


# ============================================================================
# DATA
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],
    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],
    5: [
        -62398,
        4771718,
        16027881,
    ],
    7: [
        1,
    ],
}

D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}

TRANSITIONS = [
    (1, 3),
    (3, 5),
    (5, 7),
]

DELTA_11 = 9512681472


# ============================================================================
# HELPERS
# ============================================================================

def q_value(p: int, r: int) -> int:
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def monomial_basis(p: int, d: int):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def valuation(n: int, p: int) -> int:
    n = abs(int(n))

    if n == 0:
        return 10**9

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def odd_part(n: int) -> int:
    n = abs(int(n))

    if n == 0:
        return 1

    while n % 2 == 0:
        n //= 2

    return n


# ============================================================================
# EXACT BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    M = [[int(x) for x in row] for row in A]
    n = len(M)

    if n == 1:
        return M[0][0]

    previous = 1
    sign = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):

            if M[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return 0

        if pivot_row != k:
            M[k], M[pivot_row] = M[pivot_row], M[k]
            sign *= -1

        pivot = M[k][k]

        for i in range(k + 1, n):

            for j in range(k + 1, n):

                value = (
                    M[i][j] * pivot
                    - M[i][k] * M[k][j]
                )

                if k > 0:

                    if value % previous != 0:
                        raise ArithmeticError(
                            "Bareiss division failed."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# BUILD SYSTEM
# ============================================================================

def build_full_system():

    M = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            basis = monomial_basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = basis[j] * q0

            for j in range(6):
                row[6 + j] = basis[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


# ============================================================================
# BUILD F CORE
# ============================================================================

def build_core(M):

    b0 = M[5]
    b1 = M[10]

    H = [
        [b0[12], b0[13]],
        [b1[12], b1[13]],
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_h == 0:
        raise ArithmeticError(
            "Boundary block is singular."
        )

    H_inv = [
        [
            Fraction(H[1][1], det_h),
            Fraction(-H[0][1], det_h),
        ],
        [
            Fraction(-H[1][0], det_h),
            Fraction(H[0][0], det_h),
        ],
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * b0[j]
                + H_inv[i][1] * b1[j]
            )

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        reduced = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral core entry."
                )

            reduced.append(value.numerator)

        core.append(reduced)

    return core


# ============================================================================
# NORMALIZED COFACTORS
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [i for i in range(len(A)) if i != omit_row]
    cols = [j for j in range(len(A[0])) if j != omit_col]

    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


def normalized_cofactors(A):

    n = len(A)

    R = []

    for i in range(n):

        row = []

        for j in range(n):

            mdet = bareiss_det(
                maximal_minor(A, i, j)
            )

            if mdet % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failed."
                )

            row.append(
                mdet // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# DIAGONAL UNIT GAUGES
# ============================================================================

def row_gcds(A):

    result = []

    for row in A:

        g = 0

        for x in row:

            g = __import__("math").gcd(
                g,
                abs(int(x)),
            )

        result.append(g)

    return result


def col_gcds(A):

    n = len(A[0])
    result = []

    for j in range(n):

        g = 0

        for i in range(len(A)):

            g = __import__("math").gcd(
                g,
                abs(int(A[i][j])),
            )

        result.append(g)

    return result


def make_unit_gauge_row_from_gcd(A):

    return [
        odd_part(g)
        for g in row_gcds(A)
    ]


def make_unit_gauge_col_from_gcd(A):

    return [
        odd_part(g)
        for g in col_gcds(A)
    ]


def deterministic_units(seed, n):

    # All values are odd.
    # This is deterministic and contains no randomness dependency.

    x = seed
    result = []

    for _ in range(n):

        x = (1103515245 * x + 12345) & 0x7FFFFFFF

        result.append(
            1 + 2 * (x % 15)
        )

    return result


def apply_unit_gauge(A, row_units, col_units):

    return [
        [
            row_units[i] * A[i][j] * col_units[j]
            for j in range(len(A[0]))
        ]
        for i in range(len(A))
    ]


# ============================================================================
# MODULAR RANK OVER F_2
# ============================================================================

def rank_mod_prime(A, p):

    M = [
        [x % p for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])
    rank = 0

    for col in range(cols):

        pivot = None

        for r in range(rank, rows):

            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = (
            M[pivot],
            M[rank],
        )

        inv = pow(
            M[rank][col],
            -1,
            p,
        )

        for j in range(col, cols):

            M[rank][j] = (
                M[rank][j] * inv
            ) % p

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, cols):

                M[r][j] = (
                    M[r][j]
                    - factor * M[rank][j]
                ) % p

        rank += 1

    return rank


# ============================================================================
# 2x2 MINOR VALUATION
# ============================================================================

def minimum_2x2_v2(A):

    n = len(A)

    minimum = 10**9
    positions = []

    from itertools import combinations

    for i, j in combinations(range(n), 2):

        for a, b in combinations(range(n), 2):

            value = (
                A[i][a] * A[j][b]
                - A[i][b] * A[j][a]
            )

            if value == 0:
                continue

            v = valuation(value, 2)

            if v < minimum:

                minimum = v
                positions = [
                    (i, j, a, b)
                ]

            elif v == minimum:

                positions.append(
                    (i, j, a, b)
                )

    return minimum, positions


# ============================================================================
# FIRST 2-ADIC DEFECT
# ============================================================================

def rank1_factorization_mod8(A):

    modulus = 8
    n = len(A)

    pivot_position = None
    pivot_inverse = None

    for i in range(n):

        for j in range(n):

            value = A[i][j] % modulus

            if value == 0:
                continue

            try:
                inverse = pow(
                    value,
                    -1,
                    modulus,
                )
            except ValueError:
                continue

            pivot_position = (i, j)
            pivot_inverse = inverse
            break

        if pivot_position is not None:
            break

    if pivot_position is None:
        return None

    i0, j0 = pivot_position

    u = [
        A[i][j0] % modulus
        for i in range(n)
    ]

    v = [
        A[i0][j] * pivot_inverse % modulus
        for j in range(n)
    ]

    exact = True

    for i in range(n):

        for j in range(n):

            if (
                u[i] * v[j] % modulus
                != A[i][j] % modulus
            ):

                exact = False
                break

        if not exact:
            break

    if not exact:
        return None

    E = []

    for i in range(n):

        row = []

        for j in range(n):

            difference = (
                A[i][j]
                - u[i] * v[j]
            )

            if difference % modulus != 0:
                raise ArithmeticError(
                    "Invalid mod-8 factorization."
                )

            row.append(
                difference // modulus
            )

        E.append(row)

    return {
        "position": pivot_position,
        "u": u,
        "v": v,
        "E": E,
        "defect_rank_mod2": rank_mod_prime(E, 2),
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 147 — EXACT 2-ADIC UNIT-GAUGE / "
        "RANK-1 INVARIANCE AUDIT"
    )
    print("=" * 78)

    M = build_full_system()
    core = build_core(M)
    R = normalized_cofactors(core)

    # ------------------------------------------------------------------
    # 1. RAW BASELINE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. RAW BASELINE")
    print("=" * 78)

    raw_rank = rank_mod_prime(R, 2)
    raw_min_v, raw_positions = minimum_2x2_v2(R)
    raw_defect = rank1_factorization_mod8(R)

    print(
        f"  rank_mod_2={raw_rank}"
    )

    print(
        f"  minimum_2x2_v2={raw_min_v}"
    )

    print(
        f"  rank1_mod_4={raw_min_v >= 2}"
    )

    print(
        f"  rank1_mod_8={raw_min_v >= 3}"
    )

    print(
        f"  rank1_mod_16={raw_min_v >= 4}"
    )

    print(
        f"  defect_rank_mod2="
        f"{raw_defect['defect_rank_mod2'] if raw_defect else None}"
    )

    # ------------------------------------------------------------------
    # 2. UNIT GAUGES
    # ------------------------------------------------------------------

    row_unit_gauge = make_unit_gauge_row_from_gcd(R)
    col_unit_gauge = make_unit_gauge_col_from_gcd(R)

    gauges = {
        "raw": (
            [1] * 12,
            [1] * 12,
        ),

        "odd_row_gcd_units": (
            row_unit_gauge,
            [1] * 12,
        ),

        "odd_col_gcd_units": (
            [1] * 12,
            col_unit_gauge,
        ),

        "odd_row_col_gcd_units": (
            row_unit_gauge,
            col_unit_gauge,
        ),

        "deterministic_random_units": (
            deterministic_units(17, 12),
            deterministic_units(91, 12),
        ),
    }

    print()
    print("=" * 78)
    print("2. UNIT-GAUGE DEFINITIONS")
    print("=" * 78)

    print(
        f"  odd_row_gcd_units={row_unit_gauge}"
    )

    print(
        f"  odd_col_gcd_units={col_unit_gauge}"
    )

    for name in (
        "deterministic_random_units",
    ):

        print(
            f"  {name}_rows="
            f"{gauges[name][0]}"
        )

        print(
            f"  {name}_cols="
            f"{gauges[name][1]}"
        )

    # ------------------------------------------------------------------
    # 3. UNIT GAUGE INVARIANCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT UNIT-GAUGE RANK AUDIT")
    print("=" * 78)

    results = {}

    for name, (rows, cols) in gauges.items():

        A = apply_unit_gauge(
            R,
            rows,
            cols,
        )

        rank2 = rank_mod_prime(A, 2)

        min_v, positions = minimum_2x2_v2(A)

        defect = rank1_factorization_mod8(A)

        results[name] = {
            "rank2": rank2,
            "min_v": min_v,
            "defect_rank": (
                defect["defect_rank_mod2"]
                if defect is not None
                else None
            ),
            "mod4": min_v >= 2,
            "mod8": min_v >= 3,
            "mod16": min_v >= 4,
        }

        print()
        print(
            f"  {name}"
        )

        print(
            f"    rank_mod_2={rank2}"
        )

        print(
            f"    minimum_2x2_v2={min_v}"
        )

        print(
            f"    rank1_mod_4={min_v >= 2}"
        )

        print(
            f"    rank1_mod_8={min_v >= 3}"
        )

        print(
            f"    rank1_mod_16={min_v >= 4}"
        )

        print(
            f"    defect_rank_mod2="
            f"{results[name]['defect_rank']}"
        )

    # ------------------------------------------------------------------
    # 4. INVARIANCE CHECK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. UNIT-GAUGE INVARIANCE CHECK")
    print("=" * 78)

    reference = results["raw"]

    invariance_ok = True

    for name, result in results.items():

        same = (
            result["rank2"]
            == reference["rank2"]
            and result["mod4"]
            == reference["mod4"]
            and result["mod8"]
            == reference["mod8"]
            and result["mod16"]
            == reference["mod16"]
        )

        print(
            f"  {name}: invariant={same}"
        )

        invariance_ok = (
            invariance_ok
            and same
        )

    # ------------------------------------------------------------------
    # 5. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 146 used non-unit 2-adic scalings. Those transformations
are not invertible modulo powers of 2, so changes in rank persistence
cannot be interpreted as genuine coordinate invariance.

Experiment 147 uses only odd diagonal factors.

Every odd diagonal factor is a unit in Z/2^e Z for every e >= 1.

Therefore these transformations preserve the relevant rank conditions
under invertible row and column operations.

The key diagnostic is whether:

    rank mod 2,
    rank-one modulo 4,
    rank-one modulo 8,
    rank-one modulo 16,

and the first-defect rank remain unchanged.

If they do, the 2-adic collapse is intrinsic under this class of
diagonal gauges.

If they differ, then the implementation or the chosen gauge has to be
examined carefully because a unit gauge should not change the
determinantal divisibility thresholds.

No claim is made about the origin of the matrix R.
"""
    )

    # ------------------------------------------------------------------
    # 6. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        raw_rank == 1
        and raw_min_v == 3
        and invariance_ok
    )

    print(
        f"  raw_baseline_exact={raw_rank == 1 and raw_min_v == 3}"
    )

    print(
        f"  unit_gauge_invariance={invariance_ok}"
    )

    print(
        f"  all_unit_gauges_completed={len(results) == 5}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 147 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

