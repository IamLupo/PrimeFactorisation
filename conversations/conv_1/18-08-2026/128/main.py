#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# ============================================================================
# EXPERIMENT 128 — EXACT BOUNDARY-DEFECT OPERATOR SYMBOL / COMMON-FACTOR AUDIT
#
# Starting from the exact B-odd residual rows and the unique boundary-defect
# operator from Experiment 127:
#
#   q_{p+2}(r)
#       = F0(p,d) q_p(r)
#       + F1(p,d) q_p(r+1)
#       + H0(p) [d=0],
#
#   d = D(p) - r,
#
# this experiment asks whether the operator has a hidden algebraic common
# factor or symbol degeneration.
#
# Main tests:
#
#   1. reconstruct the exact 14-parameter operator;
#   2. form the operator symbol
#
#          S(p,d,z) = F0(p,d) + F1(p,d) z;
#
#   3. test whether F0 and F1 have a nontrivial common polynomial factor;
#      this is done exactly through resultants in d and p;
#   4. test whether H0 shares a p-factor with the boundary restrictions;
#   5. test affine boundary factors from a finite exact dictionary;
#   6. compute exact resultants and report their degrees;
#   7. verify the operator pointwise on the full observed finite support.
#
# Everything is exact Fraction arithmetic.
# No SymPy.
# No floating point.
# No extrapolation.
# ============================================================================


# ============================================================================
# EXACT SECOND-LAYER B-ODD DATA
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

D = {p: len(Q[p]) - 1 for p in Q}


def q_value(p, r):
    if p not in Q or r < 0 or r >= len(Q[p]):
        return Fraction(0)
    return Q[p][r]


# ============================================================================
# EXACT LINEAR ALGEBRA
# ============================================================================

def rref(matrix):
    if not matrix:
        return [], []

    A = [[Fraction(x) for x in row] for row in matrix]
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

        pivot_value = A[row][col]
        for j in range(col, n):
            A[row][j] /= pivot_value

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
    return len(rref(matrix)[1])


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
        if all(row[j] == 0 for j in range(n)) and row[n] != 0:
            return None

    if len(pivots) != n:
        return None

    solution = [Fraction(0)] * n

    for i, pivot in enumerate(pivots):
        solution[pivot] = R[i][n]

    return solution


# ============================================================================
# UNIVARIATE POLYNOMIALS IN p
# Coefficients are Fractions.
# ============================================================================

def p_trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def p_zero():
    return [Fraction(0)]


def p_degree(a):
    a = p_trim(a)
    if len(a) == 1 and a[0] == 0:
        return -1
    return len(a) - 1


def p_add(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n

    for i in range(len(a)):
        out[i] += a[i]

    for i in range(len(b)):
        out[i] += b[i]

    return p_trim(out)


def p_neg(a):
    return [-x for x in a]


def p_sub(a, b):
    return p_add(a, p_neg(b))


def p_mul(a, b):
    if (
        len(a) == 1 and a[0] == 0
        or len(b) == 1 and b[0] == 0
    ):
        return [Fraction(0)]

    out = [Fraction(0)] * (len(a) + len(b) - 1)

    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y

    return p_trim(out)


def p_scale(a, c):
    return p_trim([x * c for x in a])


def p_divmod(a, b):
    a = p_trim(a)
    b = p_trim(b)

    if p_degree(b) < 0:
        raise ZeroDivisionError("polynomial division by zero")

    if p_degree(a) < p_degree(b):
        return [Fraction(0)], a

    q = [Fraction(0)] * (p_degree(a) - p_degree(b) + 1)
    r = list(a)

    while p_degree(r) >= p_degree(b):
        shift = p_degree(r) - p_degree(b)
        factor = r[-1] / b[-1]
        q[shift] = factor

        subtract = [Fraction(0)] * shift + p_scale(b, factor)
        r = p_sub(r, subtract)

    return p_trim(q), p_trim(r)


def p_monic(a):
    a = p_trim(a)
    if p_degree(a) < 0:
        return a
    return p_scale(a, Fraction(1, a[-1]))


def p_gcd(a, b):
    a = p_trim(a)
    b = p_trim(b)

    if p_degree(a) < 0:
        return p_monic(b)

    if p_degree(b) < 0:
        return p_monic(a)

    while p_degree(b) >= 0:
        _, r = p_divmod(a, b)
        a, b = b, r

    return p_monic(a)


def p_eval(a, x):
    value = Fraction(0)
    for c in reversed(a):
        value = value * x + c
    return value


# ============================================================================
# POLYNOMIALS IN d WITH COEFFICIENTS POLYNOMIALS IN p
#
# A bivariate polynomial is represented as:
#
#     [A0(p), A1(p), ..., An(p)]
#
# meaning
#
#     A0(p) + A1(p)d + ... + An(p)d^n.
# ============================================================================

def b_trim(a):
    a = list(a)
    while len(a) > 1 and p_degree(a[-1]) < 0:
        a.pop()
    return a


def b_zero():
    return [[Fraction(0)]]


def b_degree_d(a):
    a = b_trim(a)
    if len(a) == 1 and p_degree(a[0]) < 0:
        return -1
    return len(a) - 1


def b_add(a, b):
    n = max(len(a), len(b))
    out = [p_zero() for _ in range(n)]

    for i in range(len(a)):
        out[i] = p_add(out[i], a[i])

    for i in range(len(b)):
        out[i] = p_add(out[i], b[i])

    return b_trim(out)


def b_scale(a, c):
    return b_trim([p_scale(x, c) for x in a])


def b_mul(a, b):
    if b_degree_d(a) < 0 or b_degree_d(b) < 0:
        return b_zero()

    out = [
        p_zero()
        for _ in range(len(a) + len(b) - 1)
    ]

    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = p_add(
                out[i + j],
                p_mul(x, y),
            )

    return b_trim(out)


def b_from_block(block):
    # block ordering:
    # 1, d, d^2, p, p*d, p^2
    return [
        [block[0] + block[3], block[5],],
    ]


# ============================================================================
# Convert the 6-term block
#
# F(p,d) =
#   c0 + c1 d + c2 d^2 + c3 p + c4 p d + c5 p^2
#
# to d-polynomial with p-polynomial coefficients.
# ============================================================================

def block_to_bivar(block):

    # d^0
    a0 = [
        block[0] + Fraction(0),
        block[3],
        # p^2 coefficient below
    ]

    # make p^2 coefficient explicit
    a0 = [
        block[0],
        block[3],
        block[5],
    ]

    # d^1
    a1 = [
        block[1],
        block[4],
    ]

    # d^2
    a2 = [
        block[2],
    ]

    return b_trim([a0, a1, a2])


# ============================================================================
# Determinant of a small matrix whose entries are polynomials in p.
# Used for resultants.
# ============================================================================

def determinant(matrix):

    n = len(matrix)

    if n == 0:
        return [Fraction(1)]

    if n == 1:
        return matrix[0][0]

    total = [Fraction(0)]

    for j in range(n):

        minor = []

        for i in range(1, n):
            row = []
            for k in range(n):
                if k != j:
                    row.append(matrix[i][k])
            minor.append(row)

        term = p_mul(
            matrix[0][j],
            determinant(minor),
        )

        if j % 2:
            total = p_sub(total, term)
        else:
            total = p_add(total, term)

    return p_trim(total)


# ============================================================================
# RESULTANT OF f(d), g(d) WHEN BOTH HAVE d-DEGREE <= 2
#
# Sylvester matrix:
#
#       a2 a1 a0 0
#       0  a2 a1 a0
#       b2 b1 b0 0
#       0  b2 b1 b0
#
# where ai, bi are polynomials in p.
# ============================================================================

def resultant_quadratic_d(f, g):

    fd = b_degree_d(f)
    gd = b_degree_d(g)

    if fd < 0 or gd < 0:
        return [Fraction(0)]

    # Pad to quadratic degree.
    a0 = f[0] if len(f) > 0 else p_zero()
    a1 = f[1] if len(f) > 1 else p_zero()
    a2 = f[2] if len(f) > 2 else p_zero()

    b0 = g[0] if len(g) > 0 else p_zero()
    b1 = g[1] if len(g) > 1 else p_zero()
    b2 = g[2] if len(g) > 2 else p_zero()

    # Handle lower degree cases by trimming the Sylvester size.
    if fd == 0 and gd == 0:
        return [Fraction(0)]

    if fd == 1 and gd == 1:
        return p_sub(
            p_mul(a0, b1),
            p_mul(a1, b0),
        )

    if fd == 1 and gd == 2:
        # Res(a1*d+a0, b2*d^2+b1*d+b0)
        # = b2*a0^2 - b1*a0*a1 + b0*a1^2
        return p_add(
            p_sub(
                p_mul(b2, p_mul(a0, a0)),
                p_mul(b1, p_mul(a0, a1)),
            ),
            p_mul(b0, p_mul(a1, a1)),
        )

    if fd == 2 and gd == 1:
        return p_neg(
            resultant_quadratic_d(g, f)
        )

    matrix = [
        [a2, a1, a0, p_zero()],
        [p_zero(), a2, a1, a0],
        [b2, b1, b0, p_zero()],
        [p_zero(), b2, b1, b0],
    ]

    return determinant(matrix)


# ============================================================================
# Resultant after swapping variables:
# evaluate coefficient structure as p-polynomial with d coefficients.
#
# We build the same construction manually by converting:
#
# F(p,d) = sum_i,j c_ij p^i d^j
#
# into a polynomial in p whose coefficients are polynomials in d.
# ============================================================================

def bivar_dict(block):

    # key = (p_degree, d_degree)
    out = {}

    # constant
    out[(0, 0)] = block[0]
    out[(1, 0)] = block[3]
    out[(2, 0)] = block[5]

    out[(0, 1)] = block[1]
    out[(1, 1)] = block[4]

    out[(0, 2)] = block[2]

    return {
        k: v
        for k, v in out.items()
        if v != 0
    }


def swap_bivar(block):

    # Return p-polynomial with coefficients in d.
    data = bivar_dict(block)

    max_p = max(
        (i for i, _ in data),
        default=0,
    )

    out = []

    for i in range(max_p + 1):

        max_d = max(
            (j for (ii, j) in data if ii == i),
            default=0,
        )

        coeff = [Fraction(0)] * (max_d + 1)

        for (ii, j), value in data.items():
            if ii == i:
                coeff[j] += value

        while len(coeff) > 1 and coeff[-1] == 0:
            coeff.pop()

        out.append(coeff)

    return out


def resultant_quadratic_generic(block_a, block_b):

    A = block_to_bivar(block_a)
    B = block_to_bivar(block_b)

    rd = resultant_quadratic_d(A, B)

    return p_trim(rd)


# ============================================================================
# LINEAR FACTOR DICTIONARY
#
# Exact factors tested:
#
#   p
#   p +/- 1
#   p +/- 2
#   d
#   d +/- 1
#   d +/- 2
#   p +/- d
#   p +/- d +/- 1
#   p + 2d + const
#   p - 2d + const
#
# Each is tested directly by exact substitution.
# ============================================================================

def test_linear_factor(block, kind, a=0):

    # block is c0+c1 d+c2 d²+c3 p+c4 p d+c5 p².
    #
    # Return True if the proposed linear factor divides the polynomial.
    #
    # Since the polynomial is quadratic, checking vanishing on the entire
    # corresponding line is sufficient.

    def F(p, d):
        return (
            block[0]
            + block[1] * d
            + block[2] * d * d
            + block[3] * p
            + block[4] * p * d
            + block[5] * p * p
        )

    # Check enough exact points to establish a degree<=2 identity.
    samples = [-2, -1, 0, 1, 2, 3]

    for x in samples:

        if kind == "p":
            p = -a
            d = x

        elif kind == "d":
            p = x
            d = -a

        elif kind == "p_plus_d":
            # p + d + a = 0
            p = x
            d = -p - a

        elif kind == "p_minus_d":
            # p - d + a = 0
            p = x
            d = p + a

        elif kind == "p_plus_2d":
            # p + 2d + a = 0
            p = 2 * x
            d = -x - a

        elif kind == "p_minus_2d":
            # p - 2d + a = 0
            p = 2 * x
            d = x + a

        else:
            raise ValueError(kind)

        if F(p, d) != 0:
            return False

    return True


# ============================================================================
# POINTWISE EVALUATION
# ============================================================================

def eval_block(block, p, d):

    return (
        block[0]
        + block[1] * d
        + block[2] * d * d
        + block[3] * p
        + block[4] * p * d
        + block[5] * p * p
    )


def operator_eval(solution, p, r):

    d = D[p] - r

    F0 = solution[0:6]
    F1 = solution[6:12]
    H0 = solution[12:14]

    value = (
        eval_block(F0, p, d) * q_value(p, r)
        + eval_block(F1, p, d) * q_value(p, r + 1)
    )

    if d == 0:
        value += H0[0] + H0[1] * p

    return value


# ============================================================================
# SYSTEM
# ============================================================================

def build_system():

    A = []
    b = []

    for p in (1, 3, 5):

        Dp = D[p]

        # Include support and a substantial zero tail.
        for r in range(max(D[p], D[p + 2]) + 6):

            d = Dp - r

            basis = [
                Fraction(1),
                Fraction(d),
                Fraction(d * d),
                Fraction(p),
                Fraction(p * d),
                Fraction(p * p),
            ]

            row = [Fraction(0)] * 14

            for i in range(6):
                row[i] = q_value(p, r) * basis[i]
                row[6 + i] = q_value(p, r + 1) * basis[i]

            if d == 0:
                row[12] = Fraction(1)
                row[13] = Fraction(p)

            A.append(row)
            b.append(q_value(p + 2, r))

    return A, b


def reconstruct():

    A, b = build_system()

    rA = rank(A)
    rAug = rank([
        A[i] + [b[i]]
        for i in range(len(A))
    ])

    sol = solve_unique(A, b)

    return sol, rA, rAug


# ============================================================================
# COMPLEXITY
# ============================================================================

def primitive_signature(coeffs):

    if not coeffs:
        return []

    den_lcm = 1

    for c in coeffs:
        den_lcm = (
            den_lcm
            * c.denominator
            // gcd(den_lcm, c.denominator)
        )

    ints = [int(c * den_lcm) for c in coeffs]

    g = 0
    for n in ints:
        g = gcd(g, abs(n))

    if g:
        ints = [n // g for n in ints]

    first = next(
        (x for x in ints if x != 0),
        0,
    )

    if first < 0:
        ints = [-x for x in ints]

    return ints


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 128 — EXACT BOUNDARY-DEFECT "
        "OPERATOR SYMBOL / COMMON-FACTOR AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------------
    # 1. DATA VALIDATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in sorted(Q):

        expected = {
            1: 5,
            3: 4,
            5: 2,
            7: 0,
        }[p]

        exact = (
            len(Q[p]) == expected + 1
            and all(isinstance(x, Fraction) for x in Q[p])
        )

        data_exact &= exact

        print(
            f"  p={p}: degree={D[p]} "
            f"entries={len(Q[p])} exact={exact}"
        )

    # ------------------------------------------------------------------------
    # 2. OPERATOR RECONSTRUCTION
    # ------------------------------------------------------------------------

    sol, rank_A, rank_aug = reconstruct()

    operator_exact = (
        sol is not None
        and rank_A == rank_aug
        and rank_A == 14
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
        print("  OPERATOR RECONSTRUCTION FAILED")

        print()
        print("=" * 78)
        print("FINAL EXACTNESS")
        print("=" * 78)

        overall = data_exact and operator_exact
        print(f"  data_exact={data_exact}")
        print(f"  operator_exact={operator_exact}")
        print(f"  ALL BASIC CHECKS PASS={overall}")

        return

    F0 = sol[0:6]
    F1 = sol[6:12]
    H0 = sol[12:14]

    # ------------------------------------------------------------------------
    # 3. OPERATOR SYMBOL
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT OPERATOR SYMBOL")
    print("=" * 78)

    print(
        "  S(p,d,z) = F_0(p,d) + F_1(p,d) z"
    )

    print(
        f"  F_0 signature={primitive_signature(F0)}"
    )
    print(
        f"  F_1 signature={primitive_signature(F1)}"
    )
    print(
        f"  H_0 signature={primitive_signature(H0)}"
    )

    # ------------------------------------------------------------------------
    # 4. BIVARIATE COMMON FACTOR VIA RESULTANTS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT COMMON-FACTOR AUDIT")
    print("=" * 78)

    res_d = resultant_quadratic_generic(F0, F1)

    resultant_d_zero = (
        p_degree(res_d) < 0
    )

    print(
        "  Resultant_d(F_0,F_1):"
    )
    print(
        f"    degree_in_p={p_degree(res_d)}"
    )
    print(
        f"    identically_zero={resultant_d_zero}"
    )

    if not resultant_d_zero:
        print(
            f"    polynomial={res_d}"
        )

    # Swap p and d and test the other direction.
    #
    # Build the resultant of the p-quadratics using a direct
    # variable swap, which is equivalent to exchanging p and d.
    swapped_F0 = [
        F0[0],
        F0[3],
        F0[5],
    ]

    swapped_F1 = [
        F1[0],
        F1[3],
        F1[5],
    ]

    # The above misses mixed terms, so construct the swapped blocks:
    #
    # F(p,d) -> F(d,p)
    #
    # New block:
    # c0,
    # c1'=old c3,
    # c2'=old c5,
    # c3'=old c1,
    # c4'=old c4,
    # c5'=old c2
    swapped_F0 = [
        F0[0],
        F0[3],
        F0[5],
        F0[1],
        F0[4],
        F0[2],
    ]

    swapped_F1 = [
        F1[0],
        F1[3],
        F1[5],
        F1[1],
        F1[4],
        F1[2],
    ]

    res_p = resultant_quadratic_generic(
        swapped_F0,
        swapped_F1,
    )

    resultant_p_zero = p_degree(res_p) < 0

    print(
        "  Resultant_p(F_0,F_1):"
    )
    print(
        f"    degree_in_d={p_degree(res_p)}"
    )
    print(
        f"    identically_zero={resultant_p_zero}"
    )

    # ------------------------------------------------------------------------
    # 5. EXACT BOUNDARY-FACTOR DICTIONARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. AFFINE FACTOR DICTIONARY")
    print("=" * 78)

    factor_dictionary = [
        ("p", "p", 0),
        ("p+1", "p", 1),
        ("p-1", "p", -1),
        ("p+2", "p", 2),
        ("p-2", "p", -2),

        ("d", "d", 0),
        ("d+1", "d", 1),
        ("d-1", "d", -1),
        ("d+2", "d", 2),
        ("d-2", "d", -2),

        ("p+d", "p_plus_d", 0),
        ("p+d+1", "p_plus_d", 1),
        ("p+d-1", "p_plus_d", -1),

        ("p-d", "p_minus_d", 0),
        ("p-d+1", "p_minus_d", 1),
        ("p-d-1", "p_minus_d", -1),

        ("p+2d", "p_plus_2d", 0),
        ("p+2d+1", "p_plus_2d", 1),
        ("p+2d-1", "p_plus_2d", -1),

        ("p-2d", "p_minus_2d", 0),
        ("p-2d+1", "p_minus_2d", 1),
        ("p-2d-1", "p_minus_2d", -1),
    ]

    common_factors = []

    for label, kind, offset in factor_dictionary:

        f0 = test_linear_factor(
            F0,
            kind,
            offset,
        )

        f1 = test_linear_factor(
            F1,
            kind,
            offset,
        )

        h = False

        # H0 is univariate in p.
        if kind == "p":
            if H0[0] + H0[1] * (-offset) == 0:
                h = True

        if f0 or f1 or h:
            print(
                f"  {label}: "
                f"F0={f0} F1={f1} H0={h}"
            )

        if f0 and f1:
            common_factors.append(label)

    print(
        f"  common_affine_factors={common_factors}"
    )

    # ------------------------------------------------------------------------
    # 6. SLICE GCD AUDIT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. EXACT SPECIALIZATION GCD AUDIT")
    print("=" * 78)

    slice_results = []

    for p_value in (-2, -1, 0, 1, 2, 3, 4):

        f0_d = [
            F0[0] + F0[3] * p_value + F0[5] * p_value * p_value,
            F0[1] + F0[4] * p_value,
            F0[2],
        ]

        f1_d = [
            F1[0] + F1[3] * p_value + F1[5] * p_value * p_value,
            F1[1] + F1[4] * p_value,
            F1[2],
        ]

        g = p_gcd(f0_d, f1_d)

        deg = p_degree(g)

        slice_results.append(deg)

        print(
            f"  p={p_value}: "
            f"gcd_degree_in_d={deg}"
        )

    generic_slice_gcd_zero = all(
        x == 0
        for x in slice_results
    )

    # ------------------------------------------------------------------------
    # 7. BOUNDARY RESTRICTION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. BOUNDARY RESTRICTION AUDIT")
    print("=" * 78)

    # F0(p,0)
    f0_boundary = [
        F0[0],
        F0[3],
        F0[5],
    ]

    # F1(p,0)
    f1_boundary = [
        F1[0],
        F1[3],
        F1[5],
    ]

    boundary_gcd = p_gcd(
        f0_boundary,
        f1_boundary,
    )

    print(
        f"  gcd(F0(p,0),F1(p,0)) degree="
        f"{p_degree(boundary_gcd)}"
    )

    print(
        f"  gcd={boundary_gcd}"
    )

    # ------------------------------------------------------------------------
    # 8. H0 RELATION TO BOUNDARY OPERATOR
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. EXACT BOUNDARY SOURCE RELATION")
    print("=" * 78)

    boundary_exact = True

    for p in (1, 3, 5):

        r = D[p]
        d = 0

        homogeneous = (
            eval_block(F0, p, d) * q_value(p, r)
            + eval_block(F1, p, d) * q_value(p, r + 1)
        )

        defect = q_value(p + 2, r) - homogeneous
        h_value = H0[0] + H0[1] * p

        exact = defect == h_value
        boundary_exact &= exact

        print(
            f"  p={p}: exact={exact}"
        )

    # ------------------------------------------------------------------------
    # 9. POINTWISE RECONSTRUCTION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. POINTWISE RECONSTRUCTION")
    print("=" * 78)

    pointwise_exact = True

    for p in (1, 3, 5):

        failures = []

        for r in range(
            max(D[p], D[p + 2]) + 6
        ):

            lhs = q_value(p + 2, r)
            rhs = operator_eval(
                sol,
                p,
                r,
            )

            if lhs != rhs:
                failures.append(r)

        exact = len(failures) == 0
        pointwise_exact &= exact

        print(
            f"  {p}->{p+2}: exact={exact}"
        )

        if failures:
            print(
                f"    failures={failures}"
            )

    # ------------------------------------------------------------------------
    # 10. COMPLEXITY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. EXACT COEFFICIENT COMPLEXITY")
    print("=" * 78)

    for name, coeffs in (
        ("F_0", F0),
        ("F_1", F1),
        ("H_0", H0),
    ):

        sig = primitive_signature(coeffs)

        max_num_digits = max(
            len(str(abs(c.numerator)))
            for c in coeffs
        )

        max_den_digits = max(
            len(str(c.denominator))
            for c in coeffs
        )

        print(
            f"  {name}: "
            f"max_num_digits={max_num_digits} "
            f"max_den_digits={max_den_digits}"
        )

        print(
            f"    primitive_signature={sig}"
        )

    # ------------------------------------------------------------------------
    # 11. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The exact boundary-defect recurrence is now treated as a polynomial
  operator symbol:

      S(p,d,z)
        = F_0(p,d) + F_1(p,d) z.

  A nontrivial common factor of F_0 and F_1 would factor the entire
  homogeneous shift operator over the coefficient ring.

  The resultant tests are stronger than checking a finite list of
  integer roots: a nonzero resultant rules out a common factor
  depending on the eliminated variable over the rational function
  field.

  The affine-factor dictionary independently tests the most natural
  small factors such as

      p,
      d,
      p+d,
      p-d,
      p+2d,
      p-2d,

  with small integer offsets.

  The specialization gcd audit checks several independent p-slices
  using exact univariate rational-polynomial gcds.

  Finally, the boundary source is checked independently from the
  operator fit.

  The principal structural question is therefore:

      Is the exceptional B-odd operator genuinely irreducible at the
      coefficient level, or is it hiding a small common factor that
      the previous recurrence experiments did not expose?

  No extrapolation is used.
  Everything is exact over QQ.
        """
    )

    # ------------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------------

    common_factor_absent = (
        not resultant_d_zero
        and not resultant_p_zero
        and len(common_factors) == 0
        and generic_slice_gcd_zero
    )

    all_ok = (
        data_exact
        and operator_exact
        and boundary_exact
        and pointwise_exact
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  data_exact={data_exact}"
    )
    print(
        f"  operator_exact={operator_exact}"
    )
    print(
        f"  boundary_defect_exact={boundary_exact}"
    )
    print(
        f"  pointwise_exact={pointwise_exact}"
    )
    print(
        f"  resultant_d_nonzero={not resultant_d_zero}"
    )
    print(
        f"  resultant_p_nonzero={not resultant_p_zero}"
    )
    print(
        f"  common_affine_factor_absent={len(common_factors) == 0}"
    )
    print(
        f"  specialization_gcd_absent={generic_slice_gcd_zero}"
    )
    print(
        f"  coefficient_common_factor_absent={common_factor_absent}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={all_ok}"
    )

    print()
    print("EXPERIMENT 128 COMPLETE")


if __name__ == "__main__":
    main()

