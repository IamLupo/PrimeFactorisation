#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 144 — EXACT p-ADIC LOW-RANK FACTOR / LIFT-PERSISTENCE AUDIT
==============================================================================

Experiment 143 established:

    R mod 2 : rank = 1
    R mod 3 : rank = 1
    R mod 7 : rank = 2

and, from minor valuations,

    2-adic rank-1 condition survives through 2^3,
    3-adic rank-1 condition survives only modulo 3,
    7-adic rank-2 condition survives only modulo 7.

Experiment 144 recovers explicit low-rank factorizations modulo prime
powers.

For rank 1:

    R_ij = u_i v_j  (mod p^e)

is tested exactly by choosing a nonzero pivot R_ab and constructing

    u_i = R_ib
    v_j = R_aj / R_ab

when the pivot is invertible modulo p^e.

A rank-1 factorization is verified pointwise.

For p=7, rank 2 is treated via an exact 2-dimensional row-space basis
over the finite ring when possible.

The experiment also searches for the maximum e for which a rank-1
factorization exists for p=2,3 and the maximum e for which a rank-2
factorization exists for p=7.

The purpose is descriptive:

    finite-field rank collapse
        versus
    genuine low-rank p-adic lifting.

No claim is made about the source of the phenomenon.

No floating point.
No SymPy.
No extrapolation.
No connection to the original (p,q)-kernel is assumed.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import sys


# ============================================================================
# EXACT DATA
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
# EXACT HELPERS
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


def inv_mod(a: int, m: int):
    a %= m

    if a == 0:
        return None

    # Python's exact modular inverse.
    try:
        return pow(a, -1, m)
    except ValueError:
        return None


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
                            "Bareiss division failure."
                        )
                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# BUILD FULL SYSTEM
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
# BUILD CORE
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
        raise ArithmeticError("Singular boundary block.")

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

        out = []

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

            out.append(value.numerator)

        core.append(out)

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

            det_minor = bareiss_det(
                maximal_minor(A, i, j)
            )

            if det_minor % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failure."
                )

            row.append(
                det_minor // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# RANK-1 FACTOR TEST MOD p^e
# ============================================================================

def rank1_factorization(A, modulus):

    n = len(A)

    # We require a unit pivot.
    pivot = None
    pivot_pos = None

    for i in range(n):

        for j in range(n):

            if A[i][j] % modulus != 0:

                inv = inv_mod(
                    A[i][j],
                    modulus
                )

                if inv is not None:

                    pivot = A[i][j] % modulus
                    pivot_pos = (i, j)
                    break

        if pivot_pos is not None:
            break

    if pivot_pos is None:
        return None

    i0, j0 = pivot_pos
    pivot_inv = inv_mod(pivot, modulus)

    u = [
        A[i][j0] % modulus
        for i in range(n)
    ]

    v = [
        (A[i0][j] * pivot_inv) % modulus
        for j in range(n)
    ]

    exact = True

    for i in range(n):

        for j in range(n):

            lhs = A[i][j] % modulus
            rhs = (u[i] * v[j]) % modulus

            if lhs != rhs:
                exact = False
                break

        if not exact:
            break

    return {
        "pivot": pivot,
        "pivot_position": pivot_pos,
        "u": u,
        "v": v,
        "exact": exact,
    }


# ============================================================================
# RANK-2 FACTORIZATION MOD p^e
# ============================================================================

def rank2_factorization_mod_prime(A, prime):
    """
    Exact rank-2 factorization over F_p.

    This is deliberately only used modulo the prime itself.
    We do not pretend it is a factorization over Z/p^e when lifting
    is not established.
    """

    n = len(A)

    M = [
        [x % prime for x in row]
        for row in A
    ]

    # Find two independent rows.
    basis_rows = []

    def row_rank(rows):

        if not rows:
            return 0

        B = [r[:] for r in rows]
        rnk = 0
        cols = len(B[0])

        for c in range(cols):

            pivot = None

            for rr in range(rnk, len(B)):

                if B[rr][c] % prime != 0:
                    pivot = rr
                    break

            if pivot is None:
                continue

            B[rnk], B[pivot] = B[pivot], B[rnk]

            inv = pow(
                B[rnk][c] % prime,
                -1,
                prime,
            )

            for jj in range(c, cols):
                B[rnk][jj] = (
                    B[rnk][jj] * inv
                ) % prime

            for rr in range(len(B)):

                if rr == rnk:
                    continue

                f = B[rr][c] % prime

                if f == 0:
                    continue

                for jj in range(c, cols):
                    B[rr][jj] = (
                        B[rr][jj]
                        - f * B[rnk][jj]
                    ) % prime

            rnk += 1

        return rnk

    for i in range(n):

        candidate = basis_rows + [M[i]]

        if row_rank(candidate) > len(basis_rows):

            basis_rows.append(M[i])

            if len(basis_rows) == 2:
                break

    if len(basis_rows) != 2:
        return None

    # Express every row as alpha*u + beta*v.
    U = [
        [
            basis_rows[0][j],
            basis_rows[1][j],
        ]
        for j in range(n)
    ]

    coeffs = []

    # Find two pivot columns for the 2x2 basis.
    pivot_cols = None

    for a, b in combinations(range(n), 2):

        det = (
            basis_rows[0][a] * basis_rows[1][b]
            - basis_rows[0][b] * basis_rows[1][a]
        ) % prime

        if det != 0:

            pivot_cols = (a, b)
            inv_det = pow(det, -1, prime)
            break

    if pivot_cols is None:
        return None

    a, b = pivot_cols

    for i in range(n):

        x = M[i][a]
        y = M[i][b]

        # Solve
        # alpha*u_a + beta*v_a = x
        # alpha*u_b + beta*v_b = y
        #
        # Cramer's rule mod p.

        u_a = basis_rows[0][a]
        u_b = basis_rows[0][b]
        v_a = basis_rows[1][a]
        v_b = basis_rows[1][b]

        alpha = (
            x * v_b
            - y * v_a
        ) * inv_det % prime

        beta = (
            u_a * y
            - u_b * x
        ) * inv_det % prime

        coeffs.append(
            (alpha, beta)
        )

    exact = True

    for i in range(n):

        for j in range(n):

            predicted = (
                coeffs[i][0] * basis_rows[0][j]
                + coeffs[i][1] * basis_rows[1][j]
            ) % prime

            if predicted != M[i][j]:
                exact = False
                break

        if not exact:
            break

    return {
        "basis_rows": basis_rows,
        "pivot_cols": pivot_cols,
        "coeffs": coeffs,
        "exact": exact,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 144 — EXACT p-ADIC LOW-RANK FACTOR / "
        "LIFT-PERSISTENCE AUDIT"
    )
    print("=" * 78)

    M = build_full_system()
    core = build_core(M)
    R = normalized_cofactors(core)

    print()
    print("=" * 78)
    print("1. EXACT DATA")
    print("=" * 78)

    print("  R_shape=(12,12)")
    print(f"  Delta_11={DELTA_11}")

    # ------------------------------------------------------------------
    # 2. RANK-1 LIFT FOR 2
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. 2-ADIC RANK-1 FACTORIZATION")
    print("=" * 78)

    rank2_power_results = {}

    for e in range(1, 6):

        modulus = 2 ** e

        result = rank1_factorization(
            R,
            modulus,
        )

        exact = (
            result is not None
            and result["exact"]
        )

        rank2_power_results[e] = exact

        print()
        print(
            f"  modulus=2^{e}={modulus}"
        )

        print(
            f"    rank1_factorization={exact}"
        )

        if exact:

            print(
                f"    pivot_position="
                f"{result['pivot_position']}"
            )

            print(
                f"    u={result['u']}"
            )

            print(
                f"    v={result['v']}"
            )

    # ------------------------------------------------------------------
    # 3. 3-ADIC RANK-1 LIFT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. 3-ADIC RANK-1 FACTORIZATION")
    print("=" * 78)

    rank3_power_results = {}

    for e in range(1, 4):

        modulus = 3 ** e

        result = rank1_factorization(
            R,
            modulus,
        )

        exact = (
            result is not None
            and result["exact"]
        )

        rank3_power_results[e] = exact

        print()
        print(
            f"  modulus=3^{e}={modulus}"
        )

        print(
            f"    rank1_factorization={exact}"
        )

        if exact:

            print(
                f"    pivot_position="
                f"{result['pivot_position']}"
            )

            print(
                f"    u={result['u']}"
            )

            print(
                f"    v={result['v']}"
            )

    # ------------------------------------------------------------------
    # 4. 7-ADIC RANK-2 BASELINE + ATTEMPTED LIFT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. 7-ADIC RANK-2 FACTORIZATION")
    print("=" * 78)

    mod7 = rank2_factorization_mod_prime(
        R,
        7,
    )

    rank7_exact = (
        mod7 is not None
        and mod7["exact"]
    )

    print(
        f"  modulo=7"
    )

    print(
        f"    rank2_factorization={rank7_exact}"
    )

    if rank7_exact:

        print(
            f"    pivot_columns="
            f"{mod7['pivot_cols']}"
        )

        print(
            f"    basis_row_0="
            f"{mod7['basis_rows'][0]}"
        )

        print(
            f"    basis_row_1="
            f"{mod7['basis_rows'][1]}"
        )

        print(
            f"    coefficient_rows="
            f"{mod7['coeffs']}"
        )

    print()
    print(
        "  Full rank-2 lifting modulo 49 is deliberately "
        "not assumed."
    )

    # ------------------------------------------------------------------
    # 5. LIFT SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. LIFT-PERSISTENCE SUMMARY")
    print("=" * 78)

    max2 = max(
        [
            e
            for e, ok in rank2_power_results.items()
            if ok
        ],
        default=0,
    )

    max3 = max(
        [
            e
            for e, ok in rank3_power_results.items()
            if ok
        ],
        default=0,
    )

    print(
        f"  2-adic maximum tested rank1 precision="
        f"{max2}"
    )

    print(
        f"  3-adic maximum tested rank1 precision="
        f"{max3}"
    )

    print(
        f"  7-adic rank2 modulo-7 exact="
        f"{rank7_exact}"
    )

    # ------------------------------------------------------------------
    # 6. IMPORTANT INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
This experiment distinguishes two kinds of low-rank behavior.

A rank statement modulo p is first-order information.

A factorization modulo p^e is stronger: it asks whether the same
low-rank structure survives at higher p-adic precision.

For p=2 and p=3, rank one is tested by an explicit outer-product
factorization.

For p=7, the observed rank two is reconstructed over F_7 and the
resulting basis is checked pointwise.

The important output is therefore not the particular vectors u and v,
but the maximum prime-power precision at which an exact low-rank
factorization survives.

No causal interpretation is attached to this arithmetic pattern.
In particular, this experiment does not identify a hidden recurrence,
a determinant mechanism, or the original (p,q)-kernel.
"""
    )

    # ------------------------------------------------------------------
    # 7. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    baseline_2 = rank2_power_results.get(1, False)
    baseline_3 = rank3_power_results.get(1, False)

    final_ok = (
        baseline_2
        and baseline_3
        and rank7_exact
    )

    print(
        f"  rank1_mod_2_exact={baseline_2}"
    )

    print(
        f"  rank1_mod_3_exact={baseline_3}"
    )

    print(
        f"  rank2_mod_7_exact={rank7_exact}"
    )

    print(
        f"  lift_audit_completed=True"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 144 COMPLETE")


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

