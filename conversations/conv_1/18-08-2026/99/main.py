#!/usr/bin/env python3

import math
import sympy as sp

# ==============================================================================
# EXPERIMENT 101 — EXACT RESIDUAL-KERNEL STRUCTURE AUDIT
#
# Established:
#
#   P_k(j) = j_(k) (m-j)_(s(k)) R_k(j)
#
# with
#
#   s(k) = max(0, k-4)
#
# and
#
#   m_A = 8
#   m_B = 7
#
# Exact QQ arithmetic only.
# No floating point.
# No recurrence search.
# ==============================================================================

j = sp.symbols("j")


# ==============================================================================
# EXACT TRIANGULAR DATA
# ==============================================================================

A = [
    [
        -2, -154, sp.Rational(-818), sp.Rational(-1360, 3),
        sp.Rational(8435, 24), sp.Rational(-5851, 120),
        sp.Rational(-13373, 720), sp.Rational(51773, 5040),
        sp.Rational(-4913, 1920),
    ],
    [
        -550, -5015, -4734, sp.Rational(7879, 3),
        sp.Rational(-5147, 120), sp.Rational(-210877, 720),
        sp.Rational(83651, 720), sp.Rational(-1028053, 40320),
    ],
    [
        -7125, -14567, 4635, sp.Rational(25508, 15),
        sp.Rational(-198919, 144), sp.Rational(427555, 1008),
        sp.Rational(-3174439, 40320),
    ],
    [
        -11900, -1711, sp.Rational(28949, 6),
        sp.Rational(-31711, 15), sp.Rational(404513, 840),
        sp.Rational(-234707, 4032),
    ],
    [
        sp.Rational(-17875, 6), sp.Rational(51337, 30),
        sp.Rational(-52447, 180), sp.Rational(-14333, 210),
        sp.Rational(710501, 13440),
    ],
    [
        sp.Rational(-1001, 12), sp.Rational(26687, 360),
        sp.Rational(-40921, 1260), sp.Rational(9389, 1008),
    ],
    [
        sp.Rational(-5, 72), sp.Rational(5, 72),
        sp.Rational(-5, 144),
    ],
]

B = [
    [
        25, 619, sp.Rational(3231, 2), sp.Rational(-33, 2),
        sp.Rational(-1675, 4), sp.Rational(3363, 20),
        sp.Rational(-9991, 360), sp.Rational(-421, 2520),
    ],
    [
        1750, 8624, sp.Rational(6829, 3), sp.Rational(-27341, 8),
        sp.Rational(10551, 10), sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        9690, 10234, sp.Rational(-57829, 8),
        sp.Rational(148151, 120), sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3), sp.Rational(-19045, 12),
        sp.Rational(-5577, 4), sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24), sp.Rational(-22061, 40),
        sp.Rational(132343, 720), sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12), sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# ==============================================================================
# EXACT POLYNOMIAL HELPERS
# ==============================================================================

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def poly(expr):
    return sp.Poly(clean(expr), j, domain=sp.QQ)


def degree(expr):
    expr = clean(expr)

    if expr == 0:
        return -sp.oo

    return poly(expr).degree()


def falling(n):
    result = sp.Integer(1)

    for q in range(n):
        result *= j - q

    return sp.expand(result)


def upper_falling(m, s):
    result = sp.Integer(1)

    for q in range(s):
        result *= m - j - q

    return sp.expand(result)


def reconstruct(row, k):
    """
    Offset-row convention:

        row[offset] = C[k, k+offset]

    Therefore

        P_k(j) = sum_offset C[k,k+offset] * j_(k+offset).
    """

    total = sp.Integer(0)

    for offset, coefficient in enumerate(row):
        r = k + offset
        total += sp.Rational(coefficient) * falling(r)

    return clean(total)


def exact_quotient(numerator, denominator):
    pn = poly(numerator)
    pd = poly(denominator)

    if pd.is_zero:
        raise ZeroDivisionError("Polynomial denominator is zero.")

    quotient, remainder = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if remainder.is_zero:
        return clean(quotient.as_expr())

    return None


def scalar_multiple(f, g):
    f = clean(f)
    g = clean(g)

    if f == 0 and g == 0:
        return True, sp.Integer(1)

    if f == 0 or g == 0:
        return False, None

    pf = poly(f)
    pg = poly(g)

    if pf.degree() != pg.degree():
        return False, None

    scalar = sp.cancel(pf.LC() / pg.LC())

    if clean(f - scalar * g) == 0:
        return True, scalar

    return False, None


def primitive_integer_coefficients(expr):
    """
    Convert a QQ polynomial to a primitive integer coefficient vector.

    Correctly handles:
      * zero polynomial
      * constant polynomial
      * singleton denominator list
      * ordinary multi-term polynomials
    """

    p = poly(expr)

    if p.is_zero:
        return [0]

    coeffs = [
        sp.Rational(c)
        for c in p.all_coeffs()
    ]

    denominators = [
        int(c.q)
        for c in coeffs
    ]

    # SymPy ilcm() requires at least two arguments.
    if len(denominators) == 1:
        lcm_den = denominators[0]
    else:
        lcm_den = 1
        for d in denominators:
            lcm_den = math.lcm(lcm_den, d)

    integer_coeffs = [
        int(c * lcm_den)
        for c in coeffs
    ]

    common_gcd = 0

    for value in integer_coeffs:
        common_gcd = math.gcd(
            common_gcd,
            abs(value),
        )

    if common_gcd > 1:
        integer_coeffs = [
            value // common_gcd
            for value in integer_coeffs
        ]

    return integer_coeffs


def coefficient_vector(expr):
    p = poly(expr)

    if p.is_zero:
        return [sp.Integer(0)]

    return [
        sp.Rational(c)
        for c in p.all_coeffs()
    ]


def coefficient_matrix(expressions):
    if not expressions:
        return sp.Matrix([])

    max_degree = max(
        int(degree(expr))
        for expr in expressions
    )

    rows = []

    for expr in expressions:
        p = poly(expr)

        rows.append([
            p.nth(d)
            for d in range(max_degree, -1, -1)
        ])

    return sp.Matrix(rows)


# ==============================================================================
# BUILD ORIGINAL POLYNOMIALS
# ==============================================================================

A_poly = [
    reconstruct(row, k)
    for k, row in enumerate(A)
]

B_poly = [
    reconstruct(row, k)
    for k, row in enumerate(B)
]


# ==============================================================================
# EXTRACT RESIDUALS
# ==============================================================================

def residual_family(channel, m):
    residuals = []

    for k, P in enumerate(channel):
        s = max(0, k - 4)

        divisor = clean(
            falling(k) * upper_falling(m, s)
        )

        R = exact_quotient(
            P,
            divisor,
        )

        if R is None:
            raise RuntimeError(
                f"Residual extraction failed: k={k}, m={m}"
            )

        residuals.append(R)

    return residuals


A_residuals = residual_family(
    A_poly,
    8,
)

B_residuals = residual_family(
    B_poly,
    7,
)


# ==============================================================================
# 1. HEADER
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 101 — EXACT RESIDUAL-KERNEL STRUCTURE AUDIT")
print("=" * 78)


# ==============================================================================
# 2. RESIDUAL DEGREE PROFILE
# ==============================================================================

print()
print("=" * 78)
print("1. RESIDUAL DEGREE PROFILE")
print("=" * 78)

for name, residuals, m in [
    ("A", A_residuals, 8),
    ("B", B_residuals, 7),
]:
    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals):
        s = max(0, k - 4)

        observed = degree(R)
        predicted = m - k - s

        print(
            f"    k={k}:"
            f" s={s}"
            f" degree={observed}"
            f" predicted={predicted}"
            f" exact={observed == predicted}"
        )


# ==============================================================================
# 3. EXACT RESIDUAL POLYNOMIALS
# ==============================================================================

print()
print("=" * 78)
print("2. EXACT RESIDUAL POLYNOMIALS")
print("=" * 78)

for name, residuals in [
    ("A", A_residuals),
    ("B", B_residuals),
]:
    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals):
        print(
            f"    R_{k}(j) = {sp.factor(R)}"
        )


# ==============================================================================
# 4. PRIMITIVE INTEGER SIGNATURES
# ==============================================================================

print()
print("=" * 78)
print("3. PRIMITIVE INTEGER SIGNATURES")
print("=" * 78)

for name, residuals in [
    ("A", A_residuals),
    ("B", B_residuals),
]:
    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals):
        print(
            f"    R_{k}: "
            f"{primitive_integer_coefficients(R)}"
        )


# ==============================================================================
# 5. CONSTANT-TERM NORMALIZATION
# ==============================================================================

print()
print("=" * 78)
print("4. CONSTANT-TERM NORMALIZATION")
print("=" * 78)

for name, residuals in [
    ("A", A_residuals),
    ("B", B_residuals),
]:
    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals):
        p = poly(R)
        constant = p.nth(0)

        if constant == 0:
            print(
                f"    k={k}: constant=0"
            )
            continue

        normalized = clean(
            R / constant
        )

        print(
            f"    k={k}: {sp.factor(normalized)}"
        )


# ==============================================================================
# 6. LEADING-COEFFICIENT NORMALIZATION
# ==============================================================================

print()
print("=" * 78)
print("5. LEADING-COEFFICIENT NORMALIZATION")
print("=" * 78)

for name, residuals in [
    ("A", A_residuals),
    ("B", B_residuals),
]:
    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals):
        p = poly(R)
        lead = p.LC()

        normalized = clean(
            R / lead
        )

        print(
            f"    k={k}: {sp.factor(normalized)}"
        )


# ==============================================================================
# 7. REFLECTION AUDIT
# ==============================================================================

print()
print("=" * 78)
print("6. RESIDUAL REFLECTION AUDIT")
print("=" * 78)

for name, residuals, m in [
    ("A", A_residuals, 8),
    ("B", B_residuals, 7),
]:
    print()
    print(
        f"  {name}: j -> {m}-j"
    )

    for k, R in enumerate(residuals):

        reflected = clean(
            R.subs(
                j,
                m - j,
            )
        )

        g = sp.gcd(
            poly(R),
            poly(reflected),
        )

        proportional, scalar = scalar_multiple(
            reflected,
            R,
        )

        print(
            f"    k={k}:"
            f" gcd_degree={g.degree()}"
            f" proportional={proportional}"
            f" scalar={scalar}"
        )


# ==============================================================================
# 8. CONSECUTIVE GCD AUDIT
# ==============================================================================

print()
print("=" * 78)
print("7. CONSECUTIVE RESIDUAL-GCD AUDIT")
print("=" * 78)

for name, residuals in [
    ("A", A_residuals),
    ("B", B_residuals),
]:
    print()
    print(f"  {name} channel")

    for k in range(len(residuals) - 1):

        g = sp.gcd(
            poly(residuals[k]),
            poly(residuals[k + 1]),
        )

        print(
            f"    {k}->{k+1}:"
            f" gcd_degree={g.degree()}"
            f" gcd={sp.factor(g.as_expr())}"
        )


# ==============================================================================
# 9. CONSECUTIVE RATIO PROFILE
# ==============================================================================

print()
print("=" * 78)
print("8. CONSECUTIVE RESIDUAL RATIO AUDIT")
print("=" * 78)

for name, residuals in [
    ("A", A_residuals),
    ("B", B_residuals),
]:
    print()
    print(f"  {name} channel")

    for k in range(len(residuals) - 1):

        R0 = residuals[k]
        R1 = residuals[k + 1]

        ratio = clean(
            R1 / R0
        )

        numerator, denominator = sp.fraction(
            sp.cancel(ratio)
        )

        print(
            f"    {k}->{k+1}:"
            f" num_degree={degree(numerator)}"
            f" den_degree={degree(denominator)}"
        )

        print(
            f"      ratio={ratio}"
        )


# ==============================================================================
# 10. CROSS-CHANNEL RESIDUAL COMPARISON
# ==============================================================================

print()
print("=" * 78)
print("9. CROSS-CHANNEL RESIDUAL COMPARISON")
print("=" * 78)

overlap = min(
    len(A_residuals),
    len(B_residuals),
)

for k in range(overlap):

    RA = A_residuals[k]
    RB = B_residuals[k]

    g = sp.gcd(
        poly(RA),
        poly(RB),
    )

    proportional, scalar = scalar_multiple(
        RA,
        RB,
    )

    ratio = clean(
        RA / RB
    )

    print()
    print(f"  k={k}")
    print(
        f"    gcd_degree={g.degree()}"
    )
    print(
        f"    proportional={proportional}"
        f" scalar={scalar}"
    )
    print(
        f"    A/B={ratio}"
    )


# ==============================================================================
# 11. SHIFTED CROSS-CHANNEL COMPARISON
# ==============================================================================

print()
print("=" * 78)
print("10. SHIFTED CROSS-CHANNEL RESIDUAL AUDIT")
print("=" * 78)

for k in range(overlap):

    RA = A_residuals[k]
    RB = B_residuals[k]

    print()
    print(f"  k={k}")

    for shift in [-1, 0, 1]:

        shifted = clean(
            RB.subs(
                j,
                j + shift,
            )
        )

        g = sp.gcd(
            poly(RA),
            poly(shifted),
        )

        proportional, scalar = scalar_multiple(
            RA,
            shifted,
        )

        print(
            f"    B(j{shift:+d}) vs A(j):"
            f" gcd_degree={g.degree()}"
            f" proportional={proportional}"
            f" scalar={scalar}"
        )


# ==============================================================================
# 12. RESIDUAL COEFFICIENT-MATRIX RANK
# ==============================================================================

print()
print("=" * 78)
print("11. RESIDUAL COEFFICIENT-MATRIX RANK")
print("=" * 78)

A_matrix = coefficient_matrix(
    A_residuals
)

B_matrix = coefficient_matrix(
    B_residuals
)

print(
    f"  A:"
    f" shape={A_matrix.shape}"
    f" rank={A_matrix.rank()}"
)

print(
    f"  B:"
    f" shape={B_matrix.shape}"
    f" rank={B_matrix.rank()}"
)


# ==============================================================================
# 13. SMALL-MINOR / RANK-ONE TEST
# ==============================================================================

print()
print("=" * 78)
print("12. RESIDUAL MATRIX SEPARABILITY TEST")
print("=" * 78)


def rank_one_test(M):
    if M.rows == 0 or M.cols == 0:
        return True

    # A matrix is rank one iff every 2x2 minor vanishes.
    # For this small finite matrix, an exact rank test is simplest.
    return M.rank() <= 1


print(
    "  A rank-1 =",
    rank_one_test(A_matrix)
)

print(
    "  B rank-1 =",
    rank_one_test(B_matrix)
)


# ==============================================================================
# 14. FRESH RECONSTRUCTION
# ==============================================================================

print()
print("=" * 78)
print("13. FULL TWO-SIDED RECONSTRUCTION")
print("=" * 78)

reconstruction_ok = True

for name, channel, residuals, m in [
    ("A", A_poly, A_residuals, 8),
    ("B", B_poly, B_residuals, 7),
]:

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s = max(0, k - 4)

        rebuilt = clean(
            falling(k)
            * upper_falling(m, s)
            * residuals[k]
        )

        ok = clean(
            rebuilt - P
        ) == 0

        reconstruction_ok = (
            reconstruction_ok and ok
        )

        print(
            f"    k={k}: exact={ok}"
        )


# ==============================================================================
# 15. FRESH NUMERIC/EXACT SANITY POINT
# ==============================================================================

# This does not use floating point. It simply evaluates every original
# polynomial and its reconstructed two-sided form at several exact integer
# j-values.

print()
print("=" * 78)
print("14. EXACT POINTWISE SANITY CHECK")
print("=" * 78)

pointwise_ok = True

test_points = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]

for name, channel, residuals, m in [
    ("A", A_poly, A_residuals, 8),
    ("B", B_poly, B_residuals, 7),
]:

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s = max(0, k - 4)

        rebuilt = clean(
            falling(k)
            * upper_falling(m, s)
            * residuals[k]
        )

        local_ok = True

        for value in test_points:

            lhs = clean(
                P.subs(j, value)
            )

            rhs = clean(
                rebuilt.subs(j, value)
            )

            if lhs != rhs:
                local_ok = False
                break

        pointwise_ok = (
            pointwise_ok and local_ok
        )

        print(
            f"    k={k}: pointwise_exact={local_ok}"
        )


# ==============================================================================
# 16. STRUCTURAL SUMMARY
# ==============================================================================

print()
print("=" * 78)
print("15. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  Experiment 100 established the exact unified endpoint law

      s(k) = max(0,k-4)

  for both channels:

      A:
        P_k(j)
          = j_(k) (8-j)_(s(k)) R_A,k(j)

      B:
        P_k(j)
          = j_(k) (7-j)_(s(k)) R_B,k(j).

  This experiment removes both known endpoint factors.

  The residual kernels are therefore the part of the coefficient
  structure that has not yet been explained by endpoint combinatorics.

  The main tests are:

      * residual degree;
      * residual factorization;
      * consecutive residual gcds;
      * consecutive residual ratios;
      * endpoint reflection;
      * A/B residual alignment;
      * shifted A/B alignment;
      * coefficient-matrix rank.

  A low-rank or common-factor result would indicate that the
  endpoint factors account for only part of the structure and
  that the residuals themselves have a reusable kernel.

  A persistently full-rank / irregular residual family would indicate
  that the two endpoint factors are the main simple structure currently
  visible in the finite data.

  All computations are exact over QQ.
  No floating point.
  No extrapolation.
  No large recurrence search.
"""
)


# ==============================================================================
# 17. FINAL EXACTNESS
# ==============================================================================

failures = 0

if not reconstruction_ok:
    failures += 1

if not pointwise_ok:
    failures += 1

print()
print("=" * 78)
print("16. FINAL EXACTNESS")
print("=" * 78)

print(
    "  two_sided_reconstruction =",
    reconstruction_ok,
)

print(
    "  pointwise_exactness =",
    pointwise_ok,
)

print(
    "  failures =",
    failures,
)

print(
    "  ALL BASIC CHECKS PASS =",
    failures == 0,
)

print()
print("EXPERIMENT 101 COMPLETE")