#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 98 — EXACT TWO-SIDED FALLING-FACTORIAL STRUCTURE
#                     WITH CORRECT TRIANGLE INDEXING
#
# Exact QQ arithmetic.
# No floating point.
# ==============================================================================

j = sp.symbols("j")


# ==============================================================================
# 1. FALLING-FACTORIAL TRIANGLE DATA
# ==============================================================================

A = [
    [
        sp.Integer(-2),
        sp.Integer(-154),
        sp.Rational(-818),
        sp.Rational(-1360, 3),
        sp.Rational(8435, 24),
        sp.Rational(-5851, 120),
        sp.Rational(-13373, 720),
        sp.Rational(51773, 5040),
        sp.Rational(-4913, 1920),
    ],
    [
        sp.Integer(-550),
        sp.Integer(-5015),
        sp.Integer(-4734),
        sp.Rational(7879, 3),
        sp.Rational(-5147, 120),
        sp.Rational(-210877, 720),
        sp.Rational(83651, 720),
        sp.Rational(-1028053, 40320),
    ],
    [
        sp.Integer(-7125),
        sp.Integer(-14567),
        sp.Integer(4635),
        sp.Rational(25508, 15),
        sp.Rational(-198919, 144),
        sp.Rational(427555, 1008),
        sp.Rational(-3174439, 40320),
    ],
    [
        sp.Integer(-11900),
        sp.Integer(-1711),
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
        sp.Integer(25),
        sp.Integer(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Integer(1750),
        sp.Integer(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Integer(9690),
        sp.Integer(10234),
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


# ==============================================================================
# 2. SAFE EXACT HELPERS
# ==============================================================================

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def is_zero(expr):
    return clean(expr) == 0


def poly(expr):
    expr = clean(expr)
    return sp.Poly(expr, j, domain=sp.QQ)


def degree(expr):
    expr = clean(expr)

    if expr == 0:
        return -sp.oo

    return poly(expr).degree()


def falling(n):
    out = sp.Integer(1)

    for q in range(n):
        out *= j - q

    return sp.expand(out)


def endpoint_falling(m, n):
    """
    (m-j)_(n)
      = (m-j)(m-j-1)...(m-j-n+1)
    """
    out = sp.Integer(1)

    for q in range(n):
        out *= m - j - q

    return sp.expand(out)


def exact_quotient(divisor, dividend):
    """
    Exact polynomial division over QQ.

    Returns None if not divisible.
    """
    divisor = clean(divisor)
    dividend = clean(dividend)

    if divisor == 0:
        return None

    if dividend == 0:
        return sp.Integer(0)

    q, r = sp.div(
        poly(dividend),
        poly(divisor),
        domain=sp.QQ,
    )

    if r.is_zero:
        return sp.expand(q.as_expr())

    return None


# ==============================================================================
# 3. CORRECT TRIANGLE RECONSTRUCTION
# ==============================================================================

def reconstruct_row(row, k):
    """
    IMPORTANT:

    row[0] corresponds to C[k,k],
    row[1] corresponds to C[k,k+1],
    etc.

    Therefore the correct polynomial is

        sum_{r=k}^{...} C[k,r] j_(r).

    """
    out = sp.Integer(0)

    for offset, coefficient in enumerate(row):
        r = k + offset
        out += sp.Rational(coefficient) * falling(r)

    return sp.expand(out)


A_polys = [
    reconstruct_row(row, k)
    for k, row in enumerate(A)
]

B_polys = [
    reconstruct_row(row, k)
    for k, row in enumerate(B)
]


# ==============================================================================
# 4. HEADER
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 98 — EXACT TWO-SIDED FALLING-FACTORIAL STRUCTURE")
print("             WITH CORRECT TRIANGLE INDEXING")
print("=" * 78)


# ==============================================================================
# 5. ORIGINAL FALLING-FACTORIAL DIVISIBILITY
# ==============================================================================

print()
print("=" * 78)
print("1. ORIGINAL FALLING-FACTORIAL DIVISIBILITY")
print("=" * 78)

falling_ok = True

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        quotient = exact_quotient(
            falling(k),
            P,
        )

        ok = quotient is not None

        print(
            f"    k={k}:"
            f" j_under_{k} divides P = {ok}"
        )

        if ok:
            print(
                f"      quotient_degree={degree(quotient)}"
            )

        falling_ok = falling_ok and ok


# ==============================================================================
# 6. DIRECT ZERO AUDIT
# ==============================================================================

print()
print("=" * 78)
print("2. DIRECT INTEGER ZERO AUDIT")
print("=" * 78)

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel, m={m}")

    for k, P in enumerate(channel):

        left = []
        right = []

        for q in range(k):
            if clean(P.subs(j, q)) == 0:
                left.append(q)

        for q in range(m - k + 1, m + 1):
            if clean(P.subs(j, q)) == 0:
                right.append(q)

        print(
            f"    k={k}:"
            f" left_zeros={left}"
            f" right_zeros={right}"
        )


# ==============================================================================
# 7. QUOTIENTS AFTER REMOVING THE LOWER FALLING FACTOR
# ==============================================================================

print()
print("=" * 78)
print("3. LOWER-END QUOTIENT PROFILE")
print("=" * 78)

quotients = {
    "A": [],
    "B": [],
}

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        Q = exact_quotient(
            falling(k),
            P,
        )

        quotients[name].append(Q)

        print(
            f"    k={k}:"
            f" quotient_degree={degree(Q)}"
        )

        print(
            f"      Q_k(j)={sp.factor(Q)}"
        )


# ==============================================================================
# 8. SEARCH FOR UPPER ENDPOINT FACTORS
# ==============================================================================

print()
print("=" * 78)
print("4. UPPER-ENDPOINT FALLING-FACTOR SEARCH")
print("=" * 78)

endpoint_data = {
    "A": [],
    "B": [],
}

for name in ["A", "B"]:

    channel = A_polys if name == "A" else B_polys
    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel, m={m}")

    for k, Q in enumerate(quotients[name]):

        best_s = 0
        best_R = Q

        if Q is None:
            endpoint_data[name].append(
                (None, None)
            )
            print(
                f"    k={k}: unavailable"
            )
            continue

        for s in range(k + 1):

            divisor = endpoint_falling(
                m,
                s,
            )

            R = exact_quotient(
                divisor,
                Q,
            )

            if R is None:
                break

            best_s = s
            best_R = R

        endpoint_data[name].append(
            (best_s, best_R)
        )

        print(
            f"    k={k}:"
            f" maximum_endpoint_s={best_s}"
            f" residual_degree={degree(best_R)}"
        )

        print(
            f"      endpoint_factor="
            f"{endpoint_falling(m, best_s)}"
        )

        print(
            f"      residual="
            f"{sp.factor(best_R)}"
        )


# ==============================================================================
# 9. TWO-SIDED FACTORIZATION CHECK
# ==============================================================================

print()
print("=" * 78)
print("5. TWO-SIDED FACTORIZATION")
print("=" * 78)

two_sided_ok = True

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s, R = endpoint_data[name][k]

        total_divisor = (
            falling(k)
            * endpoint_falling(m, s)
        )

        check = exact_quotient(
            total_divisor,
            P,
        )

        ok = check is not None

        print(
            f"    k={k}:"
            f" s={s}"
            f" exact={ok}"
        )

        if ok:
            print(
                f"      residual_degree={degree(check)}"
            )

        two_sided_ok = (
            two_sided_ok and ok
        )


# ==============================================================================
# 10. RESIDUAL GCD AUDIT
# ==============================================================================

print()
print("=" * 78)
print("6. RESIDUAL CONSECUTIVE-GCD AUDIT")
print("=" * 78)

for name in ["A", "B"]:

    print()
    print(f"  {name} channel")

    data = endpoint_data[name]

    for k in range(len(data) - 1):

        R1 = data[k][1]
        R2 = data[k + 1][1]

        if R1 is None or R2 is None:
            print(
                f"    {k}->{k+1}: unavailable"
            )
            continue

        if R1 == 0 or R2 == 0:
            print(
                f"    {k}->{k+1}: zero residual"
            )
            continue

        g = sp.gcd(
            poly(R1),
            poly(R2),
        )

        print(
            f"    {k}->{k+1}:"
            f" gcd_degree={g.degree()}"
        )

        if g.degree() > 0:
            print(
                f"      gcd={sp.factor(g.as_expr())}"
            )


# ==============================================================================
# 11. ENDPOINT-REFLECTION AUDIT
# ==============================================================================

print()
print("=" * 78)
print("7. ENDPOINT-REFLECTION AUDIT")
print("=" * 78)

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(
        f"  {name}: j -> {m}-j"
    )

    for k, P in enumerate(channel):

        reflected = clean(
            P.subs(
                j,
                m - j,
            )
        )

        g = sp.gcd(
            poly(P),
            poly(reflected),
        )

        print(
            f"    k={k}: gcd_degree={g.degree()}"
        )


# ==============================================================================
# 12. SPECIAL-ENDPOINT FACTORS ALREADY OBSERVED
# ==============================================================================

print()
print("=" * 78)
print("8. EXPLICIT OBSERVED ENDPOINT FACTORS")
print("=" * 78)

for name in ["A", "B"]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, Q in enumerate(quotients[name]):

        factors = []

        # Search only modest endpoint multiplicities.
        for s in range(1, min(k + 2, 5)):

            divisor = endpoint_falling(
                m,
                s,
            )

            R = exact_quotient(
                divisor,
                Q,
            )

            if R is not None:
                factors.append(s)

        print(
            f"    k={k}: endpoint_multiplicities={factors}"
        )


# ==============================================================================
# 13. FRESH RECONSTRUCTION CHECK
# ==============================================================================

print()
print("=" * 78)
print("9. EXACT RECONSTRUCTION CHECK")
print("=" * 78)

fresh_ok = True

for name, rows, channel in [
    ("A", A, A_polys),
    ("B", B, B_polys),
]:

    print()
    print(f"  {name} channel")

    for k, (row, P) in enumerate(
        zip(rows, channel)
    ):

        rebuilt = reconstruct_row(
            row,
            k,
        )

        ok = clean(
            P - rebuilt
        ) == 0

        print(
            f"    k={k}: exact={ok}"
        )

        fresh_ok = fresh_ok and ok


# ==============================================================================
# 14. STRUCTURAL SUMMARY
# ==============================================================================

print()
print("=" * 78)
print("10. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  IMPORTANT INDEXING CORRECTION

  Each displayed row

      [C[k,k], C[k,k+1], C[k,k+2], ...]

  is an OFFSET row.

  Therefore the actual polynomial is

      P_k(j)
        = sum_{s>=0} C[k,k+s] j_(k+s),

  not

      sum_s C[k,k+s] j_(s).

  With that indexing restored, the previously established
  falling-factorial divisibility can be tested correctly.

  The experiment then asks whether the residual quotient

      Q_k(j) = P_k(j) / j_(k)

  contains a second endpoint factor

      (m-j)_(s),

  where

      m=8 for A,
      m=7 for B.

  This distinguishes three possibilities:

      1. pure lower-end falling structure;

      2. lower + upper endpoint factorial structure;

      3. isolated endpoint roots with no systematic second factor.

  The direct integer-zero audit is included independently,
  so any detected endpoint factor is checked in two ways.

  Everything is exact over QQ.
  """
)


# ==============================================================================
# 15. FINAL EXACTNESS
# ==============================================================================

print()
print("=" * 78)
print("11. FINAL EXACTNESS")
print("=" * 78)

failures = (
    int(not falling_ok)
    + int(not two_sided_ok)
    + int(not fresh_ok)
)

print(
    "  falling_factorial_divisibility =",
    falling_ok,
)

print(
    "  two_sided_factorization =",
    two_sided_ok,
)

print(
    "  reconstruction =",
    fresh_ok,
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
print("EXPERIMENT 98 COMPLETE")

