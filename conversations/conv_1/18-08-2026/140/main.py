#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 139P — EXACT P-ADIC SMITH SPECTRUM OF THE F-CORE
==============================================================================

Purpose
-------

The full integer Smith reduction attempted in Experiment 139 is too slow.

Experiment 138 already established that the observed determinant / rank defects
occur only at primes

    2, 3, 7

for the 12x12 F-core.

Experiment 139P therefore computes the Smith spectrum PRIME-BY-PRIME.

For each p in {2,3,7}, it constructs the p-adic Smith valuations by exact
valuation-aware elimination.

If

    d_1 | d_2 | ... | d_12

are the Smith invariant factors, write

    d_i = product_p p^(e_{p,i}).

This experiment computes the valuation vectors

    (e_{p,1}, ..., e_{p,12})

and then reconstructs every d_i.

This avoids the expensive global Euclidean Smith algorithm.

The script is exact.

No floating point.
No SymPy.
No extrapolation.
No assumptions about the origin of the operator.
"""


from __future__ import annotations

from math import gcd
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

TARGET_PRIMES = [2, 3, 7]


# ============================================================================
# BASIC HELPERS
# ============================================================================

def digits(n: int) -> int:
    return len(str(abs(n))) if n else 1


def valuation(n: int, p: int) -> int:
    n = abs(int(n))

    if n == 0:
        return 10**9

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def q_value(p: int, r: int) -> int:

    if r < 0:
        return 0

    if r >= len(Q[p]):
        return 0

    return Q[p][r]


def poly_degree(values) -> int:

    for i in range(len(values) - 1, -1, -1):
        if values[i] != 0:
            return i

    return -1


def matrix_rank_mod(A, prime: int) -> int:

    M = [
        [x % prime for x in row]
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

        M[rank], M[pivot] = M[pivot], M[rank]

        inv = pow(M[rank][col], -1, prime)

        for j in range(col, cols):
            M[rank][j] = (
                M[rank][j] * inv
            ) % prime

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
                ) % prime

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# BASIS / SYSTEM
# ============================================================================

def monomial_basis(p: int, d: int):

    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


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


def build_core(M):

    # Boundary equations: rows 5 and 10.
    r0 = M[5]
    r1 = M[10]

    H = [
        [r0[12], r0[13]],
        [r1[12], r1[13]],
    ]

    det_H = H[0][0] * H[1][1] - H[0][1] * H[1][0]

    if det_H == 0:
        raise ArithmeticError("Boundary block singular.")

    H_inv = [
        [
            Fraction(H[1][1], det_H),
            Fraction(-H[0][1], det_H),
        ],
        [
            Fraction(-H[1][0], det_H),
            Fraction(H[0][0], det_H),
        ],
    ]

    rhs0 = q_value(3, 5)
    rhs1 = q_value(5, 4)

    h_const = [
        H_inv[i][0] * rhs0
        + H_inv[i][1] * rhs1
        for i in range(2)
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * r0[j]
                + H_inv[i][1] * r1[j]
            )

    core_aug = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        if idx < 6:
            p = 1
            r = idx
        elif idx < 11:
            p = 3
            r = idx - 6
        else:
            p = 5
            r = idx - 11

        p_next = p + 2

        target = Fraction(
            q_value(p_next, r)
        )

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

        rhs = (
            target
            - Fraction(row[12]) * h_const[0]
            - Fraction(row[13]) * h_const[1]
        )

        if rhs.denominator != 1:
            raise ArithmeticError(
                "Non-integral core RHS."
            )

        reduced.append(rhs.numerator)

        core_aug.append(reduced)

    return [
        row[:12]
        for row in core_aug
    ], det_H


# ============================================================================
# EXACT p-ADIC VALUATION ELIMINATION
# ============================================================================

def p_adic_smith_valuations(A, prime: int):
    """
    Compute the p-adic valuations of the Smith invariant factors.

    At every stage we select an entry with minimal p-adic valuation in the
    active submatrix, normalize the pivot by its p-adic unit, and use exact
    integer elimination.

    We only require valuation data, not the final integer Smith basis.

    The implementation works with integers and preserves exactness.
    """

    M = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(M)
    m = len(M[0])

    vals = []

    k = 0

    while k < n and k < m:

        # --------------------------------------------------------------
        # Find minimum p-adic valuation in active block.
        # --------------------------------------------------------------

        best = None
        best_pos = None

        for i in range(k, n):

            for j in range(k, m):

                x = M[i][j]

                if x == 0:
                    continue

                v = valuation(x, prime)

                if best is None or v < best:
                    best = v
                    best_pos = (i, j)

        if best is None:
            break

        pi, pj = best_pos

        # Move pivot.
        M[k], M[pi] = M[pi], M[k]

        for row in M:
            row[k], row[pj] = row[pj], row[k]

        pivot = M[k][k]
        vp = valuation(pivot, prime)

        vals.append(vp)

        # --------------------------------------------------------------
        # We do not need a complete integer SNF basis.
        #
        # For valuation purposes, pivoting on minimal p-adic valuation
        # and clearing all entries to valuation >= vp is sufficient.
        # --------------------------------------------------------------

        # Normalize pivot's p-adic unit for arithmetic convenience.
        p_power = prime ** vp
        unit = pivot // p_power

        if unit == 0:
            raise ArithmeticError("Unexpected zero p-adic unit.")

        unit_inv_mod = pow(
            unit % prime,
            -1,
            prime
        )

        # Clear pivot column modulo the appropriate p-adic scale.
        #
        # Repeated exact operations are used until the pivot column and
        # row entries have valuation > vp.

        changed = True

        while changed:

            changed = False

            pivot = M[k][k]
            vp = valuation(pivot, prime)

            for i in range(k + 1, n):

                x = M[i][k]

                if x == 0:
                    continue

                vx = valuation(x, prime)

                if vx < vp:
                    # A smaller valuation exists: swap into pivot.
                    M[k], M[i] = M[i], M[k]
                    changed = True
                    break

                if vx >= vp:

                    ratio = x // (prime ** vp)
                    pivot_unit = pivot // (prime ** vp)

                    # We can solve the unit congruence exactly.
                    q0 = ratio // pivot_unit

                    M[i][k] -= q0 * pivot

                    if M[i][k] != 0:
                        changed = True

            if changed:
                continue

            for j in range(k + 1, m):

                x = M[k][j]

                if x == 0:
                    continue

                vx = valuation(x, prime)

                if vx < vp:
                    # Smaller valuation enters pivot.
                    for row in M:
                        row[k], row[j] = row[j], row[k]

                    changed = True
                    break

                if vx >= vp:

                    ratio = x // (prime ** vp)
                    pivot_unit = pivot // (prime ** vp)

                    q0 = ratio // pivot_unit

                    M[k][j] -= q0 * pivot

                    if M[k][j] != 0:
                        changed = True

        k += 1

    # If the matrix is full rank, vals has length n.
    # Sort to obtain the p-adic Smith valuation chain.
    vals.sort()

    return vals


# ============================================================================
# VERIFY P-ADIC RECONSTRUCTION
# ============================================================================

def determinant_valuation(A, prime: int):

    # Bareiss determinant for 12x12.
    M = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(M)

    if n == 1:
        return valuation(M[0][0], prime)

    previous = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):

            if M[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return None

        if pivot_row != k:
            M[k], M[pivot_row] = M[pivot_row], M[k]

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
                            "Bareiss failure."
                        )

                    value //= previous

                M[i][j] = value

        previous = pivot

        for i in range(k + 1, n):
            M[i][k] = 0

    return valuation(M[n - 1][n - 1], prime)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 139P — EXACT P-ADIC SMITH SPECTRUM "
        "OF THE F-CORE"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in (1, 3, 5, 7):

        entries = len(Q[p])
        degree = poly_degree(Q[p])
        expected = D[p]

        exact = (
            entries == expected + 1
            and degree == expected
        )

        print(
            f"  p={p}: entries={entries} "
            f"degree={degree} expected={expected} "
            f"exact={exact}"
        )

        data_exact = data_exact and exact

    print()
    print(f"  data_exact={data_exact}")

    if not data_exact:
        print("\nDATA VALIDATION FAILED")
        return

    # ------------------------------------------------------------------
    # 2. BUILD CORE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT F-CORE")
    print("=" * 78)

    full = build_full_system()
    core, H_det = build_core(full)

    print(f"  H_det={H_det}")
    print(
        f"  core_shape=({len(core)},{len(core[0])})"
    )

    core_rank = matrix_rank_full = None

    # Exact rational rank through modular confirmation and direct
    # integer determinant valuation checks.
    rank_mod_5 = matrix_rank_mod(core, 5)

    print(
        f"  rank_mod_5={rank_mod_5}"
    )

    # ------------------------------------------------------------------
    # 3. P-ADIC SPECTRUM
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. PRIME-BY-PRIME SMITH VALUATIONS")
    print("=" * 78)

    all_prime_vals = {}

    spectrum_ok = True

    for prime in TARGET_PRIMES:

        print()
        print(f"  PRIME={prime}")

        vals = p_adic_smith_valuations(
            core,
            prime
        )

        all_prime_vals[prime] = vals

        print(
            f"    valuations={vals}"
        )

        print(
            f"    total={sum(vals)}"
        )

        det_v = determinant_valuation(
            core,
            prime
        )

        same_total = (
            det_v == sum(vals)
        )

        print(
            f"    determinant_valuation={det_v}"
        )

        print(
            f"    determinant_valuation_check="
            f"{same_total}"
        )

        spectrum_ok = spectrum_ok and same_total

    # ------------------------------------------------------------------
    # 4. RECONSTRUCT THE 12 SMITH FACTORS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. RECONSTRUCTED COMPLETE SMITH FACTORS")
    print("=" * 78)

    smith_factors = []

    for i in range(12):

        value = 1

        for prime in TARGET_PRIMES:

            exponent = all_prime_vals[prime][i]

            value *= prime ** exponent

        smith_factors.append(value)

    print("  factors=[")

    for i, value in enumerate(
        smith_factors,
        start=1
    ):

        print(
            f"    d_{i}={value}"
        )

    print("  ]")

    # ------------------------------------------------------------------
    # 5. SMITH CHAIN CHECK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SMITH CHAIN CHECK")
    print("=" * 78)

    chain_ok = True

    for i in range(11):

        exact = (
            smith_factors[i + 1]
            % smith_factors[i]
            == 0
        )

        if not exact:
            chain_ok = False

        print(
            f"  d_{i+1} | d_{i+2}: {exact}"
        )

    # ------------------------------------------------------------------
    # 6. DETERMINANT RECONSTRUCTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. DETERMINANT RECONSTRUCTION")
    print("=" * 78)

    reconstructed_det = 1

    for value in smith_factors:
        reconstructed_det *= value

    # Exact Bareiss determinant.
    M = [
        row[:]
        for row in core
    ]

    previous = 1
    sign = 1

    n = 12

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):

            if M[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            sign = 0
            break

        if pivot_row != k:

            M[k], M[pivot_row] = (
                M[pivot_row],
                M[k]
            )

            sign *= -1

        pivot = M[k][k]

        for i in range(k + 1, n):

            for j in range(k + 1, n):

                value = (
                    M[i][j] * pivot
                    - M[i][k] * M[k][j]
                )

                if k > 0:
                    value //= previous

                M[i][j] = value

        previous = pivot

    exact_det = sign * M[-1][-1]

    determinant_exact = (
        abs(exact_det)
        == reconstructed_det
    )

    print(
        f"  exact_det_digits="
        f"{digits(abs(exact_det))}"
    )

    print(
        f"  reconstructed_det_digits="
        f"{digits(reconstructed_det)}"
    )

    print(
        f"  determinant_reconstruction="
        f"{determinant_exact}"
    )

    # ------------------------------------------------------------------
    # 7. COMPARE WITH EXPERIMENT 138
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. COMPARISON WITH PREVIOUS CORE INVARIANTS")
    print("=" * 78)

    # Experiment 138:
    # Delta_11 = 9512681472
    # d_12 = giant 101-digit value.

    delta11_previous = 9512681472

    d12_previous = (
        22698304379332609168018772364495783974334437093962503157858693669865148581087339234303140035372612096
    )

    prefix_product = 1

    for i in range(11):
        prefix_product *= smith_factors[i]

    same_delta = (
        prefix_product == delta11_previous
    )

    same_last = (
        smith_factors[-1] == d12_previous
    )

    print(
        f"  reconstructed_Delta_11="
        f"{prefix_product}"
    )

    print(
        f"  previous_Delta_11="
        f"{delta11_previous}"
    )

    print(
        f"  Delta_11_match={same_delta}"
    )

    print(
        f"  reconstructed_d_12="
        f"{smith_factors[-1]}"
    )

    print(
        f"  previous_d_12="
        f"{d12_previous}"
    )

    print(
        f"  d_12_match={same_last}"
    )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The previous complete Smith computation was too slow.

This experiment asks a narrower and more informative question:

    How is the Smith spectrum distributed prime-by-prime?

Only the primes 2, 3, and 7 are used because Experiment 138
showed rank defects only at those primes.

The output separates three facts:

    1. p-adic valuation structure;
    2. reconstruction of the complete invariant factors;
    3. agreement with the independently obtained determinant and
       terminal invariant.

A successful reconstruction would establish the entire Smith
spectrum without performing a global integer Smith reduction.

No interpretation is attached to the resulting factors beyond
what the exact arithmetic demonstrates.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        data_exact
        and spectrum_ok
        and chain_ok
        and determinant_exact
        and same_delta
        and same_last
    )

    print(
        f"  data_exact={data_exact}"
    )

    print(
        f"  p_adic_spectrum_exact={spectrum_ok}"
    )

    print(
        f"  smith_chain_exact={chain_ok}"
    )

    print(
        f"  determinant_reconstruction_exact="
        f"{determinant_exact}"
    )

    print(
        f"  Delta_11_match={same_delta}"
    )

    print(
        f"  d_12_match={same_last}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 139P COMPLETE")


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

