#!/usr/bin/env python3

import math
import sympy as sp

j, y, ksym = sp.symbols("j y k")


# =============================================================================
# EXPERIMENT 104
# EXACT CENTERED PARITY-COEFFICIENT LAW AUDIT
# =============================================================================
#
# Exact QQ arithmetic only.
# No floating point.
# No recurrence search.
#
# Established structure:
#
#   P_k(j) = j_(k) (m-j)_(max(0,k-4)) R_k(j)
#
# with
#
#   m = 8 for A
#   m = 7 for B.
#
# Center:
#
#   y = 2j-m.
#
# Then
#
#   R_k(j) = E_k(y) + O_k(y)
#
# where E is even and O is odd.
# =============================================================================


# =============================================================================
# 1. EXACT TRIANGULAR DATA
# =============================================================================

A = [
    [
        -2,
        -154,
        sp.Rational(-818),
        sp.Rational(-1360, 3),
        sp.Rational(8435, 24),
        sp.Rational(-5851, 120),
        sp.Rational(-13373, 720),
        sp.Rational(51773, 5040),
        sp.Rational(-4913, 1920),
    ],
    [
        -550,
        -5015,
        -4734,
        sp.Rational(7879, 3),
        sp.Rational(-5147, 120),
        sp.Rational(-210877, 720),
        sp.Rational(83651, 720),
        sp.Rational(-1028053, 40320),
    ],
    [
        -7125,
        -14567,
        4635,
        sp.Rational(25508, 15),
        sp.Rational(-198919, 144),
        sp.Rational(427555, 1008),
        sp.Rational(-3174439, 40320),
    ],
    [
        -11900,
        -1711,
        sp.Rational(28949, 6),
        sp.Rational(-31711, 15),
        sp.Rational(404513, 840),
        sp.Rational(-234707, 4032),
    ],
    [
        sp.Rational(-17875, 6),
        sp.Rational(51337, 30),
        sp.Rational(-52447, 180),
        sp.Rational(-14333, 210),
        sp.Rational(710501, 13440),
    ],
    [
        sp.Rational(-1001, 12),
        sp.Rational(26687, 360),
        sp.Rational(-40921, 1260),
        sp.Rational(9389, 1008),
    ],
    [
        sp.Rational(-5, 72),
        sp.Rational(5, 72),
        sp.Rational(-5, 144),
    ],
]

B = [
    [
        25,
        619,
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        1750,
        8624,
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        9690,
        10234,
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# =============================================================================
# 2. SAFE EXACT HELPERS
# =============================================================================

def clean(expr):
    """Exact simplification while keeping everything symbolic."""
    return sp.factor(sp.cancel(sp.expand(sp.sympify(expr))))


def is_zero(expr):
    """Exact zero test."""
    return clean(expr) == 0


def poly(expr, var):
    """
    Convert to an exact QQ polynomial.

    The caller should not use this for a zero expression when it needs
    an ordinary finite degree; use safe_degree() for that.
    """
    return sp.Poly(
        sp.expand(sp.sympify(expr)),
        var,
        domain=sp.QQ,
    )


def safe_degree(expr, var):
    """
    Return -1 for the zero polynomial instead of SymPy's -oo.
    """
    expr = clean(expr)

    if expr == 0:
        return -1

    return int(poly(expr, var).degree())


def safe_terms(expr, var):
    """Return polynomial terms, empty for zero."""
    expr = clean(expr)

    if expr == 0:
        return []

    return poly(expr, var).terms()


def coefficient(expr, var, power):
    """
    Exact coefficient. Returns 0 for zero expressions and missing powers.
    """
    expr = clean(expr)

    if expr == 0:
        return sp.Integer(0)

    p = poly(expr, var)
    return sp.Rational(p.nth(power))


def falling(x, n):
    out = sp.Integer(1)

    for q in range(n):
        out *= x - q

    return sp.expand(out)


def upper_falling(m, x, n):
    out = sp.Integer(1)

    for q in range(n):
        out *= m - x - q

    return sp.expand(out)


def exact_quotient(num, den, var):
    """
    Exact polynomial quotient over QQ.
    Returns None when division is not exact.
    """
    num = clean(num)
    den = clean(den)

    if den == 0:
        return None

    pn = poly(num, var)
    pd = poly(den, var)

    q, r = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(q.as_expr())


# =============================================================================
# 3. RECONSTRUCT P_k(j)
# =============================================================================

def reconstruct(row, k):
    out = sp.Integer(0)

    for offset, coeff in enumerate(row):
        r = k + offset
        out += sp.sympify(coeff) * falling(j, r)

    return clean(out)


# =============================================================================
# 4. RESIDUAL EXTRACTION
# =============================================================================

def residuals(channel, m):
    result = []

    for k, row in enumerate(channel):

        P = reconstruct(row, k)

        s = max(0, k - 4)

        divisor = clean(
            falling(j, k) *
            upper_falling(m, j, s)
        )

        R = exact_quotient(
            P,
            divisor,
            j,
        )

        if R is None:
            raise RuntimeError(
                f"Residual extraction failed for k={k}"
            )

        result.append(R)

    return result


# =============================================================================
# 5. CENTERING
# =============================================================================

def centered(R, m):
    """
    y = 2j-m
    hence j = (y+m)/2.
    """
    substitution = (y + sp.Integer(m)) / 2

    return clean(
        R.subs(j, substitution)
    )


# =============================================================================
# 6. EVEN / ODD PARTS
# =============================================================================

def parity_parts(F):
    """
    F(y) = E(y) + O(y)
    E(-y)=E(y)
    O(-y)=-O(y)
    """
    E = clean(
        (F + F.subs(y, -y)) / 2
    )

    O = clean(
        (F - F.subs(y, -y)) / 2
    )

    return E, O


# =============================================================================
# 7. SAFE PARITY COEFFICIENT ROW
# =============================================================================

def parity_coefficient_row(expr, parity):
    """
    Return coefficients in descending powers of the requested parity.

    parity = "even" or "odd".

    Zero polynomials return [].
    """
    expr = clean(expr)

    if expr == 0:
        return []

    d = safe_degree(expr, y)

    if d < 0:
        return []

    if parity == "even":
        powers = [
            q
            for q in range(d, -1, -1)
            if q % 2 == 0
        ]
    elif parity == "odd":
        powers = [
            q
            for q in range(d, -1, -1)
            if q % 2 == 1
        ]
    else:
        raise ValueError("parity must be 'even' or 'odd'")

    return [
        sp.Rational(
            coefficient(expr, y, q)
        )
        for q in powers
    ]


# =============================================================================
# 8. INTEGER SIGNATURE
# =============================================================================

def primitive_integer_signature(values):
    """
    Convert a rational row to primitive integer coefficients.

    Unlike sp.ilcm(*args), this explicitly handles a one-element row.
    """
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    denominator_lcm = 1

    for v in values:
        denominator_lcm = math.lcm(
            denominator_lcm,
            int(v.q),
        )

    ints = [
        int(v * denominator_lcm)
        for v in values
    ]

    g = 0

    for z in ints:
        g = math.gcd(
            g,
            abs(z),
        )

    if g:
        ints = [
            z // g
            for z in ints
        ]

    # Normalize sign.
    for z in ints:
        if z != 0:
            if z < 0:
                ints = [-q for q in ints]
            break

    return ints


# =============================================================================
# 9. EXACT INDEX POLYNOMIAL
# =============================================================================

def interpolate_index(values):
    points = [
        (
            sp.Integer(idx),
            sp.Rational(value),
        )
        for idx, value in enumerate(values)
    ]

    if not points:
        return sp.Integer(0)

    return clean(
        sp.interpolate(
            points,
            ksym,
        )
    )


# =============================================================================
# 10. EXACT DEGREE <= 2 RATIO SEARCH
# =============================================================================

def ratio_law_degree_two(values):
    """
    Search exactly for

        r(k) =
        (a2*k^2+a1*k+a0) /
        (b2*k^2+b1*k+b0).

    Returns None unless all ratios are represented exactly.
    """
    ratios = []

    for idx in range(len(values) - 1):

        a = sp.Rational(values[idx])
        b = sp.Rational(values[idx + 1])

        if a == 0:
            return None

        ratios.append(
            sp.cancel(b / a)
        )

    if len(ratios) < 3:
        return None

    equations = []

    for idx, ratio in enumerate(ratios):
        x = sp.Integer(idx)

        equations.append([
            x**2,
            x,
            1,
            -ratio * x**2,
            -ratio * x,
            -ratio,
        ])

    M = sp.Matrix(equations)

    nullspace = M.nullspace()

    if not nullspace:
        return None

    v = nullspace[0]

    numerator = clean(
        v[0] * ksym**2 +
        v[1] * ksym +
        v[2]
    )

    denominator = clean(
        v[3] * ksym**2 +
        v[4] * ksym +
        v[5]
    )

    if denominator == 0:
        return None

    # Verify exactly.
    for idx, ratio in enumerate(ratios):

        lhs = clean(
            numerator.subs(
                ksym,
                idx,
            )
        )

        rhs = clean(
            ratio *
            denominator.subs(
                ksym,
                idx,
            )
        )

        if lhs != rhs:
            return None

    return clean(
        numerator / denominator
    )


# =============================================================================
# 11. SAFE PARITY MATRIX
# =============================================================================

def parity_matrix(rows, parity):
    """
    Build a rectangular exact QQ coefficient matrix.

    Zero parity polynomials contribute all-zero rows.
    """
    if parity not in ("even", "odd"):
        raise ValueError("invalid parity")

    powers = set()

    for R in rows:

        d = safe_degree(R, y)

        if d < 0:
            continue

        if parity == "even":
            candidate = range(
                0,
                d + 1,
                2,
            )
        else:
            candidate = range(
                1,
                d + 1,
                2,
            )

        powers.update(candidate)

    powers = sorted(powers)

    if not powers:
        # Preserve the number of rows.
        return sp.zeros(
            len(rows),
            0,
        ), []

    data = []

    for R in rows:
        data.append([
            coefficient(
                R,
                y,
                q,
            )
            for q in powers
        ])

    return sp.Matrix(data), powers


# =============================================================================
# 12. BUILD DATA
# =============================================================================

A_R = residuals(A, 8)
B_R = residuals(B, 7)

A_centered = [
    centered(R, 8)
    for R in A_R
]

B_centered = [
    centered(R, 7)
    for R in B_R
]

A_parity = [
    parity_parts(F)
    for F in A_centered
]

B_parity = [
    parity_parts(F)
    for F in B_centered
]

A_E = [
    item[0]
    for item in A_parity
]

A_O = [
    item[1]
    for item in A_parity
]

B_E = [
    item[0]
    for item in B_parity
]

B_O = [
    item[1]
    for item in B_parity
]


# =============================================================================
# 13. OUTPUT
# =============================================================================

print("=" * 78)
print("EXPERIMENT 104 — EXACT CENTERED PARITY-COEFFICIENT LAW AUDIT")
print("=" * 78)


# =============================================================================
# 1. PARITY DEGREE PROFILE
# =============================================================================

print()
print("=" * 78)
print("1. PARITY DEGREE PROFILE")
print("=" * 78)

for name, rows in [
    ("A-even", A_E),
    ("A-odd", A_O),
    ("B-even", B_E),
    ("B-odd", B_O),
]:
    print()
    print(name)

    for kk, R in enumerate(rows):
        print(
            f"  k={kk}: degree={safe_degree(R, y)}"
        )


# =============================================================================
# 2. PARITY COEFFICIENT TABLES
# =============================================================================

print()
print("=" * 78)
print("2. PARITY COEFFICIENT TABLES")
print("=" * 78)

for name, rows, parity in [
    ("A-even", A_E, "even"),
    ("A-odd", A_O, "odd"),
    ("B-even", B_E, "even"),
    ("B-odd", B_O, "odd"),
]:
    print()
    print(name)

    for kk, R in enumerate(rows):
        print(
            f"  k={kk}: "
            f"{parity_coefficient_row(R, parity)}"
        )


# =============================================================================
# 3. PRIMITIVE INTEGER SIGNATURES
# =============================================================================

print()
print("=" * 78)
print("3. PRIMITIVE INTEGER PARITY SIGNATURES")
print("=" * 78)

for name, rows, parity in [
    ("A-even", A_E, "even"),
    ("A-odd", A_O, "odd"),
    ("B-even", B_E, "even"),
    ("B-odd", B_O, "odd"),
]:
    print()
    print(name)

    for kk, R in enumerate(rows):

        row = parity_coefficient_row(
            R,
            parity,
        )

        print(
            f"  k={kk}: "
            f"{primitive_integer_signature(row)}"
        )


# =============================================================================
# 4. EXACT INDEX POLYNOMIALS
# =============================================================================

print()
print("=" * 78)
print("4. EXACT INDEX POLYNOMIALS FOR INDIVIDUAL PARITY COEFFICIENTS")
print("=" * 78)

for name, rows, parity in [
    ("A-even", A_E, "even"),
    ("A-odd", A_O, "odd"),
    ("B-even", B_E, "even"),
    ("B-odd", B_O, "odd"),
]:
    print()
    print(name)

    max_degree = max(
        [
            safe_degree(R, y)
            for R in rows
        ] + [-1]
    )

    if max_degree < 0:
        print("  all rows are zero")
        continue

    if parity == "even":
        powers = range(
            0,
            max_degree + 1,
            2,
        )
    else:
        powers = range(
            1,
            max_degree + 1,
            2,
        )

    for power in powers:

        values = [
            coefficient(
                R,
                y,
                power,
            )
            for R in rows
        ]

        if all(
            value == 0
            for value in values
        ):
            continue

        P = interpolate_index(
            values
        )

        print()
        print(
            f"  y^{power}: "
            f"values={values}"
        )

        print(
            f"    degree_k={safe_degree(P, ksym)}"
        )

        print(
            f"    P(k)={P}"
        )


# =============================================================================
# 5. NORMALIZED PARITY COEFFICIENTS
# =============================================================================

print()
print("=" * 78)
print("5. NORMALIZED PARITY COEFFICIENTS")
print("=" * 78)

for name, rows, parity in [
    ("A-even", A_E, "even"),
    ("A-odd", A_O, "odd"),
    ("B-even", B_E, "even"),
    ("B-odd", B_O, "odd"),
]:
    print()
    print(name)

    max_degree = max(
        [
            safe_degree(R, y)
            for R in rows
        ] + [-1]
    )

    if max_degree < 0:
        print("  all rows are zero")
        continue

    powers = (
        range(0, max_degree + 1, 2)
        if parity == "even"
        else range(1, max_degree + 1, 2)
    )

    for power in powers:

        values = [
            coefficient(
                R,
                y,
                power,
            )
            for R in rows
        ]

        nonzero = [
            value
            for value in values
            if value != 0
        ]

        if not nonzero:
            continue

        base = nonzero[0]

        normalized = [
            sp.cancel(
                value / base
            )
            if value != 0
            else sp.Integer(0)
            for value in values
        ]

        print(
            f"  y^{power}: "
            f"{normalized}"
        )


# =============================================================================
# 6. SUCCESSIVE-k RATIO SEARCH
# =============================================================================

print()
print("=" * 78)
print("6. SUCCESSIVE-k RATIO SEARCH")
print("=" * 78)

for name, rows, parity in [
    ("A-even", A_E, "even"),
    ("A-odd", A_O, "odd"),
    ("B-even", B_E, "even"),
    ("B-odd", B_O, "odd"),
]:
    print()
    print(name)

    max_degree = max(
        [
            safe_degree(R, y)
            for R in rows
        ] + [-1]
    )

    if max_degree < 0:
        print("  all rows are zero")
        continue

    powers = (
        range(0, max_degree + 1, 2)
        if parity == "even"
        else range(1, max_degree + 1, 2)
    )

    for power in powers:

        values = [
            coefficient(
                R,
                y,
                power,
            )
            for R in rows
        ]

        if all(
            value == 0
            for value in values
        ):
            continue

        law = ratio_law_degree_two(
            values
        )

        print(
            f"  y^{power}: "
            f"degree<=2_ratio={law}"
        )


# =============================================================================
# 7. PARITY COEFFICIENT-MATRIX RANK
# =============================================================================

print()
print("=" * 78)
print("7. PARITY COEFFICIENT-MATRIX RANK")
print("=" * 78)

matrix_data = {}

for name, rows, parity in [
    ("A-even", A_E, "even"),
    ("A-odd", A_O, "odd"),
    ("B-even", B_E, "even"),
    ("B-odd", B_O, "odd"),
]:

    M, powers = parity_matrix(
        rows,
        parity,
    )

    matrix_data[name] = (
        M,
        powers,
    )

    print()
    print(
        f"  {name}: "
        f"shape={M.shape} "
        f"rank={M.rank()}"
    )


# =============================================================================
# 8. CROSS-CHANNEL PARITY COMPARISON
# =============================================================================

print()
print("=" * 78)
print("8. CROSS-CHANNEL PARITY COMPARISON")
print("=" * 78)

overlap = min(
    len(A_E),
    len(B_E),
)

for kk in range(overlap):

    print()
    print(
        f"  k={kk}"
    )

    print("    EVEN")

    powers = sorted(
        set(
            q
            for R in [
                A_E[kk],
                B_E[kk],
            ]
            for q in range(
                0,
                max(
                    0,
                    safe_degree(R, y),
                ) + 1,
                2,
            )
        )
    )

    for q in powers:

        av = coefficient(
            A_E[kk],
            y,
            q,
        )

        bv = coefficient(
            B_E[kk],
            y,
            q,
        )

        ratio = (
            "undefined"
            if bv == 0
            else sp.cancel(av / bv)
        )

        print(
            f"      y^{q}: "
            f"A/B={ratio}"
        )

    print("    ODD")

    powers = sorted(
        set(
            q
            for R in [
                A_O[kk],
                B_O[kk],
            ]
            for q in range(
                1,
                max(
                    0,
                    safe_degree(R, y),
                ) + 1,
                2,
            )
        )
    )

    for q in powers:

        av = coefficient(
            A_O[kk],
            y,
            q,
        )

        bv = coefficient(
            B_O[kk],
            y,
            q,
        )

        ratio = (
            "undefined"
            if bv == 0
            else sp.cancel(av / bv)
        )

        print(
            f"      y^{q}: "
            f"A/B={ratio}"
        )


# =============================================================================
# 9. ZERO-POLYNOMIAL SANITY
# =============================================================================

print()
print("=" * 78)
print("9. ZERO-POLYNOMIAL SANITY AUDIT")
print("=" * 78)

for name, rows in [
    ("A-even", A_E),
    ("A-odd", A_O),
    ("B-even", B_E),
    ("B-odd", B_O),
]:
    print()
    print(name)

    for kk, R in enumerate(rows):
        print(
            f"  k={kk}: "
            f"is_zero={is_zero(R)} "
            f"degree={safe_degree(R, y)}"
        )


# =============================================================================
# 10. EXACT CENTERED RECONSTRUCTION
# =============================================================================

print()
print("=" * 78)
print("10. EXACT CENTERED RECONSTRUCTION")
print("=" * 78)

reconstruction_ok = True

for name, originals, centered_rows, m in [
    ("A", A_R, A_centered, 8),
    ("B", B_R, B_centered, 7),
]:

    print()
    print(name)

    for kk, (R, F) in enumerate(
        zip(
            originals,
            centered_rows,
        )
    ):

        rebuilt = clean(
            F.subs(
                y,
                2 * j - m,
            )
        )

        ok = (
            clean(
                rebuilt - R
            ) == 0
        )

        reconstruction_ok &= ok

        print(
            f"  k={kk}: exact={ok}"
        )


# =============================================================================
# 11. INTERPRETATION
# =============================================================================

print()
print("=" * 78)
print("11. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
r"""
  The centered residuals satisfy

      R_k(j) = E_k(y) + O_k(y),

      y = 2j-m,

  with

      E_k(-y) = E_k(y),
      O_k(-y) = -O_k(y).

  This experiment studies the parity sectors coefficient-by-coefficient.

  The main objects are

      E_k(y) = sum_r e[k,r] y^(2r)

      O_k(y) = y sum_r o[k,r] y^(2r).

  For every nonzero coefficient family we test:

      * exact polynomial dependence on k;
      * exact degree <= 2 rational transition laws;
      * normalization across k;
      * coefficient-matrix rank;
      * A/B cross-channel ratios.

  A zero parity polynomial is treated as an exact zero object rather
  than being assigned SymPy's -oo polynomial degree.

  The purpose is to determine whether centering really separates the
  residual kernel into simpler indexed components.

  Everything is exact over QQ.
  No floating point.
  No recurrence search.
  No extrapolation.
"""
)


# =============================================================================
# 12. FINAL EXACTNESS
# =============================================================================

print()
print("=" * 78)
print("12. FINAL EXACTNESS")
print("=" * 78)

print(
    f"  centered_reconstruction = {reconstruction_ok}"
)

print(
    f"  failures = {0 if reconstruction_ok else 1}"
)

print(
    f"  ALL BASIC CHECKS PASS = {reconstruction_ok}"
)

print()
print("EXPERIMENT 104 COMPLETE")