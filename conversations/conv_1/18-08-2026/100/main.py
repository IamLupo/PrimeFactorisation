#!/usr/bin/env python3

import math
import sympy as sp

# ==============================================================================
# EXPERIMENT 103 — EXACT CENTERED/PARITY RESIDUAL-KERNEL AUDIT
# ==============================================================================

j, y = sp.symbols("j y")


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
# EXACT SYMBOLIC HELPERS
# ==============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def poly(expr, var):
    expr = sp.sympify(expr)
    return sp.Poly(
        sp.expand(expr),
        var,
        domain=sp.QQ,
    )


def degree(expr, var):
    expr = sp.expand(sp.sympify(expr))

    if expr == 0:
        return sp.S.NegativeInfinity

    return poly(expr, var).degree()


def degree_int(expr, var):
    d = degree(expr, var)

    if d == sp.S.NegativeInfinity:
        return -1

    return int(d)


def falling(var, n):
    out = sp.Integer(1)

    for q in range(n):
        out *= var - q

    return sp.expand(out)


def upper_falling(m, var, s):
    out = sp.Integer(1)

    for q in range(s):
        out *= m - var - q

    return sp.expand(out)


def exact_quotient(numerator, denominator, var):
    numerator = sp.expand(sp.sympify(numerator))
    denominator = sp.expand(sp.sympify(denominator))

    if denominator == 0:
        raise ZeroDivisionError("zero polynomial denominator")

    pn = poly(numerator, var)
    pd = poly(denominator, var)

    q, r = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(q.as_expr())


# ==============================================================================
# ORIGINAL TRIANGLE RECONSTRUCTION
# ==============================================================================

def reconstruct(row, k):
    out = sp.Integer(0)

    for offset, coefficient in enumerate(row):
        r = k + offset

        out += (
            sp.sympify(coefficient)
            * falling(j, r)
        )

    return clean(out)


def build_residuals(channel, m):
    residuals = []

    for k, row in enumerate(channel):
        P = reconstruct(row, k)

        s = max(0, k - 4)

        lower = falling(j, k)
        upper = upper_falling(m, j, s)

        divisor = clean(lower * upper)

        R = exact_quotient(
            P,
            divisor,
            j,
        )

        if R is None:
            raise RuntimeError(
                f"Residual extraction failed "
                f"for k={k}, m={m}"
            )

        residuals.append(R)

    return residuals


# ==============================================================================
# CORRECT CENTERING
# ==============================================================================

def centered_form(expr, m):
    """
    y = 2j - m
    therefore
        j = (y+m)/2.

    IMPORTANT:
    Do NOT use sp.Rational(y+m, 2).
    Rational() is for numeric rationals only.
    """
    substitution = (
        y + sp.Integer(m)
    ) / sp.Integer(2)

    return clean(
        sp.expand(
            sp.sympify(expr).subs(
                j,
                substitution,
            )
        )
    )


def parity_parts(expr):
    expr = clean(expr)

    even = clean(
        (
            expr
            + expr.subs(y, -y)
        ) / sp.Integer(2)
    )

    odd = clean(
        (
            expr
            - expr.subs(y, -y)
        ) / sp.Integer(2)
    )

    return even, odd


# ==============================================================================
# SAFE GCD
# ==============================================================================

def gcd_poly(a, b, var):
    a = clean(a)
    b = clean(b)

    if a == 0 and b == 0:
        return sp.Poly(
            0,
            var,
            domain=sp.QQ,
        )

    if a == 0:
        return poly(b, var)

    if b == 0:
        return poly(a, var)

    return sp.gcd(
        poly(a, var),
        poly(b, var),
    )


def gcd_degree(a, b, var):
    g = gcd_poly(a, b, var)

    if g.is_zero:
        return -1

    return int(g.degree())


def proportional(a, b, var):
    a = clean(a)
    b = clean(b)

    if a == 0 and b == 0:
        return True, sp.Integer(1)

    if a == 0 or b == 0:
        return False, None

    pa = poly(a, var)
    pb = poly(b, var)

    if pa.degree() != pb.degree():
        return False, None

    scale = sp.cancel(
        pa.LC() / pb.LC()
    )

    ok = clean(
        a - scale * b
    ) == 0

    return ok, (
        scale if ok else None
    )


# ==============================================================================
# PRIMITIVE INTEGER SIGNATURE
# ==============================================================================

def primitive_integer_coefficients(expr, var):
    expr = clean(expr)

    if expr == 0:
        return [0]

    p = poly(expr, var)

    coeffs = [
        sp.Rational(c)
        for c in p.all_coeffs()
    ]

    denominator_lcm = 1

    for c in coeffs:
        denominator_lcm = math.lcm(
            denominator_lcm,
            int(c.q),
        )

    integers = [
        int(c * denominator_lcm)
        for c in coeffs
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    if g > 1:
        integers = [
            value // g
            for value in integers
        ]

    # Stable sign normalization.
    for value in integers:
        if value != 0:
            if value < 0:
                integers = [
                    -v for v in integers
                ]
            break

    return integers


# ==============================================================================
# SAFE COEFFICIENT MATRIX
# ==============================================================================

def coefficient_matrix(expressions, var):
    """
    Handles zero expressions explicitly.

    For a zero polynomial the degree is -oo,
    so it must never be passed through int().
    """
    if not expressions:
        return sp.Matrix([])

    nonzero = [
        clean(expr)
        for expr in expressions
        if clean(expr) != 0
    ]

    if not nonzero:
        return sp.zeros(
            len(expressions),
            1,
        )

    max_degree = max(
        degree_int(expr, var)
        for expr in nonzero
    )

    matrix_rows = []

    for expr in expressions:
        expr = clean(expr)

        if expr == 0:
            matrix_rows.append(
                [
                    sp.Integer(0)
                    for _ in range(
                        max_degree + 1
                    )
                ]
            )
            continue

        p = poly(expr, var)

        matrix_rows.append(
            [
                sp.Rational(
                    p.nth(power)
                )
                for power in range(
                    max_degree,
                    -1,
                    -1,
                )
            ]
        )

    return sp.Matrix(matrix_rows)


# ==============================================================================
# BUILD OBJECTS
# ==============================================================================

A_residuals = build_residuals(
    A,
    8,
)

B_residuals = build_residuals(
    B,
    7,
)

A_centered = [
    centered_form(R, 8)
    for R in A_residuals
]

B_centered = [
    centered_form(R, 7)
    for R in B_residuals
]

A_even = []
A_odd = []

for R in A_centered:
    E, O = parity_parts(R)
    A_even.append(E)
    A_odd.append(O)

B_even = []
B_odd = []

for R in B_centered:
    E, O = parity_parts(R)
    B_even.append(E)
    B_odd.append(O)


# ==============================================================================
# OUTPUT
# ==============================================================================

print("=" * 78)
print(
    "EXPERIMENT 103 — "
    "EXACT CENTERED/PARITY RESIDUAL-KERNEL AUDIT"
)
print("=" * 78)


# ------------------------------------------------------------------------------
# 1
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("1. CENTERED RESIDUAL POLYNOMIALS")
print("=" * 78)

print()
print("A channel: y = 2*j - 8")

for k, R in enumerate(A_centered):
    print(
        f"  k={k}: "
        f"degree={degree_int(R, y)}"
    )
    print(
        f"    F_{k}(y) = "
        f"{sp.factor(R)}"
    )

print()
print("B channel: y = 2*j - 7")

for k, R in enumerate(B_centered):
    print(
        f"  k={k}: "
        f"degree={degree_int(R, y)}"
    )
    print(
        f"    F_{k}(y) = "
        f"{sp.factor(R)}"
    )


# ------------------------------------------------------------------------------
# 2
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("2. EVEN / ODD DECOMPOSITION")
print("=" * 78)

for name, even_list, odd_list in [
    ("A", A_even, A_odd),
    ("B", B_even, B_odd),
]:

    print()
    print(f"{name} channel")

    for k, (E, O) in enumerate(
        zip(even_list, odd_list)
    ):

        print()
        print(f"  k={k}")

        print(
            f"    even_degree="
            f"{degree_int(E, y)}"
        )

        print(
            f"    odd_degree="
            f"{degree_int(O, y)}"
        )

        print(
            f"    E_{k}(y) = "
            f"{sp.factor(E)}"
        )

        print(
            f"    O_{k}(y) = "
            f"{sp.factor(O)}"
        )


# ------------------------------------------------------------------------------
# 3
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("3. PARITY-DOMINANCE AUDIT")
print("=" * 78)

for name, even_list, odd_list in [
    ("A", A_even, A_odd),
    ("B", B_even, B_odd),
]:

    print()
    print(f"{name} channel")

    for k, (E, O) in enumerate(
        zip(even_list, odd_list)
    ):

        E_zero = clean(E) == 0
        O_zero = clean(O) == 0

        if E_zero and O_zero:
            status = "ZERO"
        elif O_zero:
            status = "PURE_EVEN"
        elif E_zero:
            status = "PURE_ODD"
        else:
            de = degree_int(E, y)
            do = degree_int(O, y)

            if de > do:
                status = "EVEN_DOMINANT"
            elif do > de:
                status = "ODD_DOMINANT"
            else:
                status = "BALANCED"

        print(
            f"  k={k}: "
            f"even_degree={degree_int(E, y)} "
            f"odd_degree={degree_int(O, y)} "
            f"status={status}"
        )


# ------------------------------------------------------------------------------
# 4
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("4. CENTERED INTEGER ROOT AUDIT")
print("=" * 78)

for name, centered in [
    ("A", A_centered),
    ("B", B_centered),
]:

    print()
    print(f"{name} channel")

    for k, F in enumerate(centered):

        roots = []

        for value in range(
            -12,
            13,
        ):
            if clean(
                F.subs(
                    y,
                    value,
                )
            ) == 0:
                roots.append(value)

        print(
            f"  k={k}: "
            f"integer_y_roots={roots}"
        )


# ------------------------------------------------------------------------------
# 5
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("5. CENTERED PRIMITIVE INTEGER SIGNATURES")
print("=" * 78)

for name, centered in [
    ("A", A_centered),
    ("B", B_centered),
]:

    print()
    print(f"{name} channel")

    for k, F in enumerate(centered):

        print(
            f"  k={k}: "
            f"{primitive_integer_coefficients(F, y)}"
        )


# ------------------------------------------------------------------------------
# 6
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("6. EVEN / ODD PRIMITIVE SIGNATURES")
print("=" * 78)

for name, even_list, odd_list in [
    ("A", A_even, A_odd),
    ("B", B_even, B_odd),
]:

    print()
    print(f"{name} channel")

    for k, (E, O) in enumerate(
        zip(even_list, odd_list)
    ):

        print()
        print(f"  k={k}")

        print(
            "    even =",
            primitive_integer_coefficients(
                E,
                y,
            ),
        )

        print(
            "    odd  =",
            primitive_integer_coefficients(
                O,
                y,
            ),
        )


# ------------------------------------------------------------------------------
# 7
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("7. EVEN / ODD COEFFICIENT-MATRIX RANKS")
print("=" * 78)

for name, even_list, odd_list in [
    ("A", A_even, A_odd),
    ("B", B_even, B_odd),
]:

    ME = coefficient_matrix(
        even_list,
        y,
    )

    MO = coefficient_matrix(
        odd_list,
        y,
    )

    print()
    print(f"{name} channel")

    print(
        f"  even_shape={ME.shape} "
        f"even_rank={ME.rank()}"
    )

    print(
        f"  odd_shape={MO.shape} "
        f"odd_rank={MO.rank()}"
    )


# ------------------------------------------------------------------------------
# 8
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("8. A/B CENTERED PARITY ALIGNMENT")
print("=" * 78)

overlap = min(
    len(A_even),
    len(B_even),
)

for k in range(overlap):

    print()
    print(f"  k={k}")

    eprop, escale = proportional(
        A_even[k],
        B_even[k],
        y,
    )

    oprop, oscale = proportional(
        A_odd[k],
        B_odd[k],
        y,
    )

    print(
        f"    even-even:"
        f" proportional={eprop}"
        f" scalar={escale}"
        f" gcd_degree="
        f"{gcd_degree(A_even[k], B_even[k], y)}"
    )

    print(
        f"    odd-odd:"
        f" proportional={oprop}"
        f" scalar={oscale}"
        f" gcd_degree="
        f"{gcd_degree(A_odd[k], B_odd[k], y)}"
    )


# ------------------------------------------------------------------------------
# 9
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("9. SHIFTED A/B CENTERED COMPARISON")
print("=" * 78)

for k in range(overlap):

    print()
    print(f"  k={k}")

    for shift in [-2, -1, 0, 1, 2]:

        shifted_even = clean(
            B_even[k].subs(
                y,
                y + shift,
            )
        )

        shifted_odd = clean(
            B_odd[k].subs(
                y,
                y + shift,
            )
        )

        eprop, escale = proportional(
            A_even[k],
            shifted_even,
            y,
        )

        oprop, oscale = proportional(
            A_odd[k],
            shifted_odd,
            y,
        )

        print(
            f"  shift={shift:+d}: "
            f"even_prop={eprop} "
            f"even_scale={escale} "
            f"odd_prop={oprop} "
            f"odd_scale={oscale}"
        )


# ------------------------------------------------------------------------------
# 10
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("10. CENTERED DERIVATIVE GCD AUDIT")
print("=" * 78)

for name, centered in [
    ("A", A_centered),
    ("B", B_centered),
]:

    print()
    print(f"{name} channel")

    for k in range(
        len(centered) - 1
    ):

        F = centered[k]
        G = centered[k + 1]

        dF = clean(
            sp.diff(F, y)
        )

        print(
            f"  {k}->{k+1}: "
            f"gcd_degree(G,dF)="
            f"{gcd_degree(G, dF, y)}"
        )


# ------------------------------------------------------------------------------
# 11
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("11. ODD-PART / y REDUCTION")
print("=" * 78)

for name, odd_list in [
    ("A", A_odd),
    ("B", B_odd),
]:

    print()
    print(f"{name} channel")

    for k, O in enumerate(odd_list):

        if clean(O) == 0:
            reduced = sp.Integer(0)
        else:
            reduced = exact_quotient(
                O,
                y,
                y,
            )

            if reduced is None:
                raise RuntimeError(
                    f"{name} odd part at k={k} "
                    "is not divisible by y."
                )

        print()
        print(f"  k={k}")

        print(
            f"    O/y degree="
            f"{degree_int(reduced, y)}"
        )

        print(
            f"    O/y = "
            f"{sp.factor(reduced)}"
        )


# ------------------------------------------------------------------------------
# 12
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("12. EXACT CENTERED RECONSTRUCTION")
print("=" * 78)

reconstruction_ok = True

for name, originals, centered, m in [
    ("A", A_residuals, A_centered, 8),
    ("B", B_residuals, B_centered, 7),
]:

    print()
    print(f"{name} channel")

    for k, (R, F) in enumerate(
        zip(originals, centered)
    ):

        rebuilt = clean(
            F.subs(
                y,
                2 * j - m,
            )
        )

        ok = clean(
            rebuilt - R
        ) == 0

        reconstruction_ok = (
            reconstruction_ok and ok
        )

        print(
            f"  k={k}: exact={ok}"
        )


# ------------------------------------------------------------------------------
# 13
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("13. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
r"""
  The established exact residual structure is

      P_k(j)
        = j_(k)
          (m-j)_(max(0,k-4))
          R_k(j).

  This experiment centers each residual at the midpoint:

      y = 2j-m.

  Thus

      j = (y+m)/2.

  Each centered residual is decomposed exactly as

      R_k = E_k + O_k,

  with

      E_k(-y) = E_k(y),
      O_k(-y) = -O_k(y).

  The main questions are whether the midpoint exposes:

      * pure parity;
      * lower-dimensional even/odd sectors;
      * low-rank coefficient matrices;
      * A/B alignment;
      * shifted channel equivalences;
      * derivative/common-factor structure.

  All tests are finite exact identities over QQ.

  Zero polynomials are handled explicitly.
  Symbolic substitutions are never passed through Rational().
  No floating point.
  No extrapolation.
  No recurrence search.
"""
)


# ------------------------------------------------------------------------------
# 14
# ------------------------------------------------------------------------------

print()
print("=" * 78)
print("14. FINAL EXACTNESS")
print("=" * 78)

print(
    f"  centered_reconstruction = "
    f"{reconstruction_ok}"
)

print(
    f"  ALL BASIC CHECKS PASS = "
    f"{reconstruction_ok}"
)

print()
print("EXPERIMENT 103 COMPLETE")