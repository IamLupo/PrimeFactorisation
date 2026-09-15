#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 100 — EXACT UNIFIED TWO-SIDED ENDPOINT LAW
#
# Candidate law:
#
#     s(k) = max(0, k-4)
#
# for BOTH channels.
#
# A:
#     P_k(j) = j_(k) (8-j)_(s(k)) R_k(j)
#
# B:
#     P_k(j) = j_(k) (7-j)_(s(k)) R_k(j)
#
# Exact QQ arithmetic only.
# No floating point.
# No recurrence search.
# ==============================================================================

j = sp.symbols("j")


# ==============================================================================
# 1. EXACT FALLING-BASIS DATA
# ==============================================================================

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


# ==============================================================================
# 2. HELPERS
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


def endpoint_falling(m, s):
    result = sp.Integer(1)

    for q in range(s):
        result *= m - j - q

    return sp.expand(result)


def exact_quotient(dividend, divisor):
    dividend = poly(dividend)
    divisor = poly(divisor)

    q, r = sp.div(
        dividend,
        divisor,
        domain=sp.QQ,
    )

    if r.is_zero:
        return sp.expand(q.as_expr())

    return None


def reconstruct(row, k):
    result = sp.Integer(0)

    for offset, coefficient in enumerate(row):
        r = k + offset
        result += sp.Rational(coefficient) * falling(r)

    return sp.expand(result)


def evaluate_zero(expr, value):
    return clean(expr.subs(j, value)) == 0


def exact_scalar_multiple(f, g):
    """
    Return (True, c) if f = c*g exactly.
    """
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

    c = sp.cancel(pf.LC() / pg.LC())

    if clean(f - c * g) == 0:
        return True, c

    return False, None


# ==============================================================================
# 3. RECONSTRUCT P_k(j)
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
# 4. HEADER
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 100 — EXACT UNIFIED TWO-SIDED ENDPOINT LAW")
print("=" * 78)


# ==============================================================================
# 5. CANDIDATE LAW
# ==============================================================================

print()
print("=" * 78)
print("1. CANDIDATE UNIFIED MULTIPLICITY LAW")
print("=" * 78)

def s_law(k):
    return max(0, k - 4)


print()
print("  s(k) = max(0, k-4)")
print("  A endpoint m = 8")
print("  B endpoint m = 7")

for name, channel in [
    ("A", A_poly),
    ("B", B_poly),
]:
    print()
    print(f"  {name}:")
    for k in range(len(channel)):
        print(
            f"    k={k}: s(k)={s_law(k)}"
        )


# ==============================================================================
# 6. EXACT TWO-SIDED DIVISIBILITY
# ==============================================================================

print()
print("=" * 78)
print("2. EXACT TWO-SIDED DIVISIBILITY")
print("=" * 78)

divisibility_ok = True

residuals = {
    "A": [],
    "B": [],
}

for name, channel in [
    ("A", A_poly),
    ("B", B_poly),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel, m={m}")

    for k, P in enumerate(channel):

        s = s_law(k)

        divisor = clean(
            falling(k)
            * endpoint_falling(m, s)
        )

        R = exact_quotient(P, divisor)

        ok = R is not None

        divisibility_ok = divisibility_ok and ok

        print(
            f"    k={k}:"
            f" lower={k}"
            f" upper={s}"
            f" exact={ok}"
        )

        residuals[name].append(R)


# ==============================================================================
# 7. MAXIMALITY OF s(k)
# ==============================================================================

print()
print("=" * 78)
print("3. EXACT MAXIMALITY AUDIT")
print("=" * 78)

maximality_ok = True

for name, channel in [
    ("A", A_poly),
    ("B", B_poly),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s = s_law(k)

        current_divisor = clean(
            falling(k)
            * endpoint_falling(m, s)
        )

        current = exact_quotient(
            P,
            current_divisor,
        )

        stronger_divisor = clean(
            falling(k)
            * endpoint_falling(m, s + 1)
        )

        stronger = exact_quotient(
            P,
            stronger_divisor,
        )

        maximal = (
            current is not None
            and stronger is None
        )

        maximality_ok = (
            maximality_ok and maximal
        )

        print(
            f"    k={k}:"
            f" s={s}"
            f" maximal={maximal}"
        )


# ==============================================================================
# 8. RESIDUAL DEGREE LAW
# ==============================================================================

print()
print("=" * 78)
print("4. RESIDUAL DEGREE LAW")
print("=" * 78)

residual_degree_ok = True

for name in ["A", "B"]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals[name]):

        if R is None:
            residual_degree_ok = False
            print(
                f"    k={k}: unavailable"
            )
            continue

        observed = degree(R)

        predicted = (
            m
            - k
            - s_law(k)
        )

        ok = observed == predicted

        residual_degree_ok = (
            residual_degree_ok and ok
        )

        print(
            f"    k={k}:"
            f" observed={observed}"
            f" predicted={predicted}"
            f" exact={ok}"
        )


# ==============================================================================
# 9. RESIDUAL FACTORIZATION
# ==============================================================================

print()
print("=" * 78)
print("5. RESIDUAL FACTORIZATION")
print("=" * 78)

for name in ["A", "B"]:

    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals[name]):

        if R is None:
            continue

        print(
            f"    R_{k}(j) = {sp.factor(R)}"
        )


# ==============================================================================
# 10. RESIDUAL INTEGER ZERO AUDIT
# ==============================================================================

print()
print("=" * 78)
print("6. RESIDUAL INTEGER-ZERO AUDIT")
print("=" * 78)

for name in ["A", "B"]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals[name]):

        if R is None:
            continue

        zeros = [
            q
            for q in range(m + 1)
            if evaluate_zero(R, q)
        ]

        print(
            f"    k={k}: zeros={zeros}"
        )


# ==============================================================================
# 11. RESIDUAL REFLECTION STRUCTURE
# ==============================================================================

print()
print("=" * 78)
print("7. RESIDUAL ENDPOINT-REFLECTION STRUCTURE")
print("=" * 78)

for name in ["A", "B"]:

    m = 8 if name == "A" else 7

    print()
    print(
        f"  {name}: reflection j -> {m}-j"
    )

    for k, R in enumerate(residuals[name]):

        if R is None:
            continue

        reflected = clean(
            R.subs(
                j,
                m - j,
            )
        )

        proportional, scalar = exact_scalar_multiple(
            reflected,
            R,
        )

        gcd = sp.gcd(
            poly(R),
            poly(reflected),
        )

        print(
            f"    k={k}:"
            f" gcd_degree={gcd.degree()}"
            f" proportional={proportional}"
            f" scalar={scalar}"
        )


# ==============================================================================
# 12. RESIDUAL CONSECUTIVE GCD
# ==============================================================================

print()
print("=" * 78)
print("8. RESIDUAL CONSECUTIVE-GCD AUDIT")
print("=" * 78)

for name in ["A", "B"]:

    channel = residuals[name]

    print()
    print(f"  {name} channel")

    for k in range(len(channel) - 1):

        R1 = channel[k]
        R2 = channel[k + 1]

        if R1 is None or R2 is None:
            print(
                f"    {k}->{k+1}: unavailable"
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


# ==============================================================================
# 13. RESIDUAL NORMALIZED POLYNOMIALS
# ==============================================================================

print()
print("=" * 78)
print("9. RESIDUAL NORMALIZATION")
print("=" * 78)

for name in ["A", "B"]:

    print()
    print(f"  {name} channel")

    for k, R in enumerate(residuals[name]):

        if R is None:
            continue

        p = poly(R)

        print(
            f"    k={k}:"
            f" degree={p.degree()}"
            f" leading={p.LC()}"
            f" constant={p.TC()}"
        )


# ==============================================================================
# 14. RECONSTRUCTION CHECK
# ==============================================================================

print()
print("=" * 78)
print("10. FULL RECONSTRUCTION")
print("=" * 78)

reconstruction_ok = True

for name, channel in [
    ("A", A_poly),
    ("B", B_poly),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s = s_law(k)
        R = residuals[name][k]

        if R is None:
            reconstruction_ok = False
            print(
                f"    k={k}: missing residual"
            )
            continue

        rebuilt = clean(
            falling(k)
            * endpoint_falling(m, s)
            * R
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
# 15. COMPARE AGAINST OLD WRONG B LAW
# ==============================================================================

print()
print("=" * 78)
print("11. OLD-vs-CORRECTED B INDEX LAW")
print("=" * 78)

print()
print("  Old tested law:")
print("      s_B(k) = max(0,k-3)")
print()
print("  Corrected law:")
print("      s_B(k) = max(0,k-4)")
print()

old_B_ok = True

for k, P in enumerate(B_poly):

    old_s = max(0, k - 3)

    divisor = clean(
        falling(k)
        * endpoint_falling(7, old_s)
    )

    ok = exact_quotient(
        P,
        divisor,
    ) is not None

    old_B_ok = old_B_ok and ok

    print(
        f"    k={k}:"
        f" old_s={old_s}"
        f" exact={ok}"
    )


# ==============================================================================
# 16. STRUCTURAL INTERPRETATION
# ==============================================================================

print()
print("=" * 78)
print("12. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  Experiment 98 established

      P_k(j) = j_(k) Q_k(j).

  Experiment 99 tested the candidate

      A: s_A(k)=max(0,k-4)
      B: s_B(k)=max(0,k-3)

  and the B law failed at k=4 and k=5.

  The direct zero data instead supports the unified law

      s(k)=max(0,k-4)

  for BOTH channels.

  Therefore the corrected structural candidate is

      A:
        P_k(j)
          = j_(k)
            (8-j)_(max(0,k-4))
            R_k(j)

      B:
        P_k(j)
          = j_(k)
            (7-j)_(max(0,k-4))
            R_k(j).

  This experiment tests:

      * exact divisibility;
      * maximality;
      * residual degree;
      * residual factorization;
      * endpoint reflection;
      * consecutive residual gcds;
      * exact reconstruction.

  The most useful outcome is not another recurrence.
  It is establishing whether BOTH channels share the same
  endpoint-transition law.

  Everything is exact over QQ.
  No floating point.
  No extrapolation.
  """
)


# ==============================================================================
# 17. FINAL EXACTNESS
# ==============================================================================

failures = (
    int(not divisibility_ok)
    + int(not maximality_ok)
    + int(not residual_degree_ok)
    + int(not reconstruction_ok)
)

print()
print("=" * 78)
print("13. FINAL EXACTNESS")
print("=" * 78)

print(
    "  corrected_endpoint_law =",
    divisibility_ok,
)

print(
    "  endpoint_maximality =",
    maximality_ok,
)

print(
    "  residual_degree_law =",
    residual_degree_ok,
)

print(
    "  reconstruction =",
    reconstruction_ok,
)

print(
    "  old_B_law_valid =",
    old_B_ok,
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
print("EXPERIMENT 100 COMPLETE")

