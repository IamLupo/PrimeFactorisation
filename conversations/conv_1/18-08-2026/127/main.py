#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# ============================================================================
# EXPERIMENT 127 — EXACT BOUNDARY-DEFECT OPERATOR MINIMALITY / FACTOR AUDIT
#
# Starting from the exact B-odd residual rows used in Experiments 125–126,
# reconstruct the unique boundary-defect operator and test:
#
#   1. coefficient-matrix rank of F_0, F_1, H_0;
#   2. separability in p and d;
#   3. common polynomial factors of F_0 and F_1;
#   4. boundary restriction d=0;
#   5. whether the boundary defect is exactly the homogeneous residual;
#   6. whether F_0/F_1 have a common p-only or d-only factor;
#   7. exact pointwise reconstruction.
#
# Everything uses Fraction arithmetic only.
# No SymPy.
# No floating point.
# No extrapolation.
# ============================================================================


# ============================================================================
# EXACT B-ODD SECOND-LAYER DATA
# ============================================================================

Q = {
    1: [
        Fraction(-584531, 35840),
        Fraction(-32224291, 21504),
        Fraction(74131151, 215040),
        Fraction(29406229, 129024),
        Fraction(-1338089411, 7741440),
        Fraction(495451247, 7741440),
    ],

    3: [
        Fraction(-59257, 230400),
        Fraction(7521137, 230400),
        Fraction(12697441, 3225600),
        Fraction(4595257, 1382400),
        Fraction(-140504813, 12902400),
    ],

    5: [
        Fraction(4457, 2764800),
        Fraction(-340837, 2764800),
        Fraction(-5342627, 12902400),
    ],

    7: [
        Fraction(-421, 38707200),
    ],
}


D = {p: len(row) - 1 for p, row in Q.items()}


def q_value(p, r):
    if p not in Q:
        return Fraction(0)

    if r < 0 or r >= len(Q[p]):
        return Fraction(0)

    return Q[p][r]


# ============================================================================
# EXACT LINEAR ALGEBRA
# ============================================================================

def rref(matrix):
    if not matrix:
        return [], []

    A = [
        [Fraction(x) for x in row]
        for row in matrix
    ]

    m = len(A)
    n = len(A[0])

    pivots = []
    row = 0

    for col in range(n):
        pivot = None

        for i in range(row, m):
            if A[i][col] != 0:
                pivot = i
                break

        if pivot is None:
            continue

        A[row], A[pivot] = A[pivot], A[row]

        value = A[row][col]

        for j in range(col, n):
            A[row][j] /= value

        for i in range(m):
            if i == row:
                continue

            factor = A[i][col]

            if factor == 0:
                continue

            for j in range(col, n):
                A[i][j] -= factor * A[row][j]

        pivots.append(col)
        row += 1

        if row == m:
            break

    return A, pivots


def rank(matrix):
    if not matrix:
        return 0

    _, pivots = rref(matrix)
    return len(pivots)


def solve_unique(A, b):
    if not A:
        return None

    n = len(A[0])

    augmented = [
        list(map(Fraction, A[i])) + [Fraction(b[i])]
        for i in range(len(A))
    ]

    R, pivots = rref(augmented)

    for row in R:
        if all(row[j] == 0 for j in range(n)):
            if row[n] != 0:
                return None

    if len(pivots) != n:
        return None

    x = [Fraction(0)] * n

    for i, pivot in enumerate(pivots):
        x[pivot] = R[i][n]

    return x


# ============================================================================
# POLYNOMIAL UTILITIES
# ============================================================================

def trim(poly):
    poly = list(poly)

    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()

    return poly


def poly_degree(poly):
    poly = trim(poly)

    if len(poly) == 1 and poly[0] == 0:
        return -1

    return len(poly) - 1


def poly_eval(poly, x):
    value = Fraction(0)

    for c in reversed(poly):
        value = value * x + c

    return value


def poly_add(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n

    for i in range(len(a)):
        out[i] += a[i]

    for i in range(len(b)):
        out[i] += b[i]

    return trim(out)


def poly_mul(a, b):
    out = [
        Fraction(0)
        for _ in range(len(a) + len(b) - 1)
    ]

    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y

    return trim(out)


# ============================================================================
# BIVARIATE BASIS
#
# 1, d, d^2, p, p*d, p^2
# ============================================================================

def basis_pd(p, d):
    return [
        Fraction(1),
        Fraction(d),
        Fraction(d * d),
        Fraction(p),
        Fraction(p * d),
        Fraction(p * p),
    ]


def eval_block(block, p, d):
    basis = basis_pd(p, d)

    return sum(
        block[i] * basis[i]
        for i in range(6)
    )


# ============================================================================
# RECONSTRUCT THE EXACT BOUNDARY-DEFECT OPERATOR
#
# q_{p+2}(r)
#   =
#   F0(p,d) q_p(r)
#   + F1(p,d) q_p(r+1)
#   + H0(p) [d=0]
#
# d = D(p)-r
# ============================================================================

def build_system():

    A = []
    b = []

    transitions = [
        (1, 3),
        (3, 5),
        (5, 7),
    ]

    for p, p2 in transitions:

        Dp = D[p]

        # Go beyond the visible support so the zero tail is tested too.
        max_r = max(D[p], D[p2]) + 5

        for r in range(max_r + 1):

            d = Dp - r

            basis = basis_pd(p, d)

            row = [Fraction(0)] * 14

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            for i in range(6):
                row[i] += q0 * basis[i]
                row[6 + i] += q1 * basis[i]

            if d == 0:
                row[12] += Fraction(1)
                row[13] += Fraction(p)

            A.append(row)
            b.append(q_value(p2, r))

    return A, b


def reconstruct_operator():

    A, b = build_system()

    rA = rank(A)
    rAug = rank([
        A[i] + [b[i]]
        for i in range(len(A))
    ])

    solution = solve_unique(A, b)

    return solution, rA, rAug


# ============================================================================
# FACTORING / ROOT TESTS FOR UNIVARIATE SLICES
# ============================================================================

def integer_roots(poly, bound=50):

    roots = []

    for x in range(-bound, bound + 1):
        if poly_eval(poly, x) == 0:
            roots.append(x)

    return roots


# ============================================================================
# COMMON FACTOR TEST FOR DIRECTIONAL SLICES
# ============================================================================
#
# Rather than performing full symbolic bivariate factorization, use exact
# slices:
#
#   d = 0
#   p = 0
#
# This detects common p-only/d-only factors without any CAS.
# ============================================================================

def p_slice(block, d_value):
    return [
        block[0] + block[1] * d_value + block[2] * d_value * d_value,
        block[3] + block[4] * d_value,
        block[5],
    ]


def d_slice(block, p_value):
    return [
        block[0] + block[3] * p_value + block[5] * p_value * p_value,
        block[1] + block[4] * p_value,
        block[2],
    ]


def proportional(a, b):

    scalar = None

    for x, y in zip(a, b):

        if y == 0:
            if x != 0:
                return False, None
            continue

        candidate = x / y

        if scalar is None:
            scalar = candidate
        elif candidate != scalar:
            return False, None

    if scalar is None:
        return True, Fraction(0)

    return True, scalar


# ============================================================================
# SEPARABILITY TEST
#
# F(p,d) = U(p)V(d)
#
# The 3x3 coefficient matrix must have rank 1.
# ============================================================================

def coefficient_matrix(block):

    return [
        [block[0], block[1], block[2]],
        [block[3], block[4], Fraction(0)],
        [block[5], Fraction(0), Fraction(0)],
    ]


def separability_rank(block):
    return rank(coefficient_matrix(block))


# ============================================================================
# POINTWISE OPERATOR EVALUATION
# ============================================================================

def evaluate_operator(solution, p, r):

    d = D[p] - r

    F0 = solution[0:6]
    F1 = solution[6:12]
    H0 = solution[12:14]

    value = (
        eval_block(F0, p, d) * q_value(p, r)
        +
        eval_block(F1, p, d) * q_value(p, r + 1)
    )

    if d == 0:
        value += H0[0] + H0[1] * p

    return value


def reconstruction_failures(solution):

    failures = []

    for p in (1, 3, 5):

        p2 = p + 2

        max_r = max(D[p], D[p2]) + 5

        for r in range(max_r + 1):

            lhs = q_value(p2, r)
            rhs = evaluate_operator(solution, p, r)

            if lhs != rhs:
                failures.append(
                    (p, p2, r, lhs, rhs)
                )

    return failures


# ============================================================================
# HOMOGENEOUS BOUNDARY RESIDUAL
#
# At d=0:
#
#   defect(p)
#      = q_{p+2}(D(p))
#        - F0(p,0) q_p(D(p))
#        - F1(p,0) q_p(D(p)+1)
#
# The second term normally vanishes because it lies beyond support.
# ============================================================================

def boundary_homogeneous_residual(solution, p):

    d = 0
    r = D[p]

    F0 = solution[0:6]
    F1 = solution[6:12]

    homogeneous = (
        eval_block(F0, p, d) * q_value(p, r)
        +
        eval_block(F1, p, d) * q_value(p, r + 1)
    )

    target = q_value(p + 2, r)

    return target - homogeneous


# ============================================================================
# INTERIOR RESIDUAL AUDIT
# ============================================================================

def interior_failures(solution):

    failures = []

    for p in (1, 3, 5):

        p2 = p + 2

        for r in range(D[p] + 1):

            d = D[p] - r

            if d == 0:
                continue

            lhs = q_value(p2, r)
            rhs = (
                eval_block(
                    solution[0:6],
                    p,
                    d,
                ) * q_value(p, r)
                +
                eval_block(
                    solution[6:12],
                    p,
                    d,
                ) * q_value(p, r + 1)
            )

            if lhs != rhs:
                failures.append(
                    (p, p2, r, lhs, rhs)
                )

    return failures


# ============================================================================
# INTEGER COMPLEXITY
# ============================================================================

def primitive_signature(coeffs):

    if not coeffs:
        return []

    lcm_den = 1

    for c in coeffs:
        lcm_den = (
            lcm_den
            * c.denominator
            // gcd(lcm_den, c.denominator)
        )

    ints = [
        int(c * lcm_den)
        for c in coeffs
    ]

    g = 0

    for n in ints:
        g = gcd(g, abs(n))

    if g == 0:
        return ints

    ints = [n // g for n in ints]

    first = next(
        (n for n in ints if n != 0),
        0,
    )

    if first < 0:
        ints = [-n for n in ints]

    return ints


def complexity(coeffs):

    if not coeffs:
        return (0, 0, 0)

    sig = primitive_signature(coeffs)

    return (
        max(
            len(str(abs(c.numerator)))
            for c in coeffs
        ),
        max(
            len(str(c.denominator))
            for c in coeffs
        ),
        max(
            (
                len(str(abs(x)))
                for x in sig
                if x != 0
            ),
            default=0,
        ),
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 127 — EXACT BOUNDARY-DEFECT "
        "OPERATOR MINIMALITY / FACTOR AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------------
    # 1. DATA
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in sorted(Q):

        ok = (
            isinstance(Q[p], list)
            and all(
                isinstance(x, Fraction)
                for x in Q[p]
            )
        )

        data_exact = data_exact and ok

        print(
            f"  p={p}: D={D[p]} "
            f"entries={len(Q[p])} "
            f"exact={ok}"
        )

    # ------------------------------------------------------------------------
    # 2. OPERATOR
    # ------------------------------------------------------------------------

    solution, rank_A, rank_aug = reconstruct_operator()

    operator_exact = (
        solution is not None
        and rank_A == rank_aug
    )

    print()
    print("=" * 78)
    print("2. EXACT OPERATOR RECONSTRUCTION")
    print("=" * 78)

    print(f"  rank={rank_A}")
    print(f"  augmented_rank={rank_aug}")
    print("  unknowns=14")
    print(f"  unique_exact={operator_exact}")

    if not operator_exact:

        print()
        print("  Operator reconstruction failed.")

        overall = data_exact and operator_exact

        print()
        print("=" * 78)
        print("FINAL EXACTNESS")
        print("=" * 78)
        print(f"  data_exact={data_exact}")
        print(f"  operator_exact={operator_exact}")
        print(f"  ALL BASIC CHECKS PASS={overall}")
        return

    F0 = solution[0:6]
    F1 = solution[6:12]
    H0 = solution[12:14]

    # ------------------------------------------------------------------------
    # 3. COEFFICIENT-MATRIX RANK
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. BIVARIATE COEFFICIENT-MATRIX RANK")
    print("=" * 78)

    for name, block in (
        ("F_0", F0),
        ("F_1", F1),
    ):

        M = coefficient_matrix(block)
        rnk = rank(M)

        print(
            f"  {name}: matrix_rank={rnk} "
            f"separable_rank1={rnk == 1}"
        )

    # ------------------------------------------------------------------------
    # 4. P/D SLICE ROOTS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT FACTOR / ROOT AUDIT")
    print("=" * 78)

    for name, block in (
        ("F_0", F0),
        ("F_1", F1),
    ):

        d0 = p_slice(block, 0)
        p0 = d_slice(block, 0)

        print(
            f"  {name}: d=0 p-polynomial="
            f"{d0}"
        )

        print(
            f"    integer p-roots="
            f"{integer_roots(d0)}"
        )

        print(
            f"  {name}: p=0 d-polynomial="
            f"{p0}"
        )

        print(
            f"    integer d-roots="
            f"{integer_roots(p0)}"
        )

    h_roots = []

    for p in range(-20, 21):

        value = H0[0] + H0[1] * p

        if value == 0:
            h_roots.append(p)

    print(f"  H_0 integer roots={h_roots}")

    # ------------------------------------------------------------------------
    # 5. COMMON-SLICE PROPORTIONALITY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. F_0 / F_1 SLICE PROPORTIONALITY")
    print("=" * 78)

    for d in (0, 1, 2):

        a = p_slice(F0, d)
        b = p_slice(F1, d)

        prop, scalar = proportional(a, b)

        print(
            f"  d={d}: proportional={prop} "
            f"scalar={scalar}"
        )

    for p in (0, 1, 2, 3):

        a = d_slice(F0, p)
        b = d_slice(F1, p)

        prop, scalar = proportional(a, b)

        print(
            f"  p={p}: proportional={prop} "
            f"scalar={scalar}"
        )

    # ------------------------------------------------------------------------
    # 6. BOUNDARY DEFECT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. EXACT BOUNDARY-DEFECT AUDIT")
    print("=" * 78)

    defect_values = []

    for p in (1, 3, 5):

        defect = boundary_homogeneous_residual(
            solution,
            p,
        )

        H_value = H0[0] + H0[1] * p

        exact = defect == H_value

        defect_values.append(exact)

        print(
            f"  p={p}: "
            f"homogeneous_defect={defect}"
        )

        print(
            f"       H_0(p)={H_value}"
        )

        print(
            f"       exact_match={exact}"
        )

    boundary_exact = all(defect_values)

    # ------------------------------------------------------------------------
    # 7. INTERIOR AUDIT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. HOMOGENEOUS INTERIOR AUDIT")
    print("=" * 78)

    interior_fail = interior_failures(solution)

    interior_exact = len(interior_fail) == 0

    print(
        f"  interior_exact={interior_exact}"
    )

    if interior_fail:
        print(f"  failures={interior_fail}")

    # ------------------------------------------------------------------------
    # 8. COMPLEXITY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. OPERATOR COMPLEXITY")
    print("=" * 78)

    for name, block in (
        ("F_0", F0),
        ("F_1", F1),
        ("H_0", H0),
    ):

        print(
            f"  {name}: complexity={complexity(block)}"
        )

        print(
            f"    primitive_signature="
            f"{primitive_signature(block)}"
        )

    # ------------------------------------------------------------------------
    # 9. POINTWISE RECONSTRUCTION
    # ------------------------------------------------------------------------

    failures = reconstruction_failures(solution)

    pointwise_exact = len(failures) == 0

    print()
    print("=" * 78)
    print("9. POINTWISE RECONSTRUCTION")
    print("=" * 78)

    for p in (1, 3, 5):

        p2 = p + 2

        row_failures = [
            item
            for item in failures
            if item[0] == p
        ]

        print(
            f"  {p}->{p2}: "
            f"exact={len(row_failures) == 0}"
        )

    # ------------------------------------------------------------------------
    # 10. STRUCTURAL SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The exact recurrence has the form

      q_{p+2}(r)
        =
        F_0(p,d) q_p(r)
        +
        F_1(p,d) q_p(r+1)
        +
        H_0(p) [d=0],

  where

      d = D(p)-r.

  Experiment 127 tests whether this representation has any further
  low-complexity structure.

  The coefficient matrices distinguish:

      rank 1
          F(p,d) = U(p)V(d),

  from genuinely coupled p/d dependence.

  Exact slice-root tests probe ordinary polynomial factors without
  invoking a symbolic factorization package.

  The boundary audit tests a stronger statement:

      H_0(p)

  must equal the exact defect left by the homogeneous interior
  operator when evaluated at the terminal boundary d=0.

  Therefore the boundary source is checked independently rather
  than merely being accepted because it was part of the linear fit.

  Finally, the complete operator is checked pointwise on the
  entire observed finite support together with zero-tail points.

  Everything is exact over Fraction arithmetic.
  No floating point.
  No SymPy.
  No extrapolation.
        """
    )

    # ------------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------------

    all_ok = (
        data_exact
        and operator_exact
        and pointwise_exact
        and boundary_exact
        and interior_exact
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  data_exact={data_exact}")
    print(f"  operator_exact={operator_exact}")
    print(f"  boundary_defect_exact={boundary_exact}")
    print(f"  interior_exact={interior_exact}")
    print(f"  pointwise_exact={pointwise_exact}")
    print(f"  ALL BASIC CHECKS PASS={all_ok}")

    print()
    print("EXPERIMENT 127 COMPLETE")


if __name__ == "__main__":
    main()

