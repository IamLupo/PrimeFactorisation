#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 99 — EXACT TWO-SIDED ENDPOINT LAW / RESIDUAL KERNEL AUDIT
#
# Target:
#
#   P_k(j) = j_(k) (m-j)_(s_k) R_k(j)
#
# with the candidate laws
#
#   A: m=8, s_k=max(0,k-4)
#   B: m=7, s_k=max(0,k-3)
#
# Exact QQ arithmetic only.
# Floating point forbidden.
# No large recurrence search.
# ==============================================================================

j = sp.symbols("j")


# ==============================================================================
# 1. EXACT FALLING-BASIS DATA
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
# 2. EXACT HELPERS
# ==============================================================================

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def poly(expr):
    return sp.Poly(
        clean(expr),
        j,
        domain=sp.QQ,
    )


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


def endpoint_falling(m, s):
    out = sp.Integer(1)

    for q in range(s):
        out *= m - j - q

    return sp.expand(out)


def exact_divide(divisor, dividend):
    divisor = clean(divisor)
    dividend = clean(dividend)

    q, r = sp.div(
        poly(dividend),
        poly(divisor),
        domain=sp.QQ,
    )

    if r.is_zero:
        return sp.expand(q.as_expr())

    return None


def reconstruct_row(row, k):
    """
    row[0] = C[k,k]
    row[1] = C[k,k+1]
    etc.
    """
    out = sp.Integer(0)

    for offset, coefficient in enumerate(row):
        r = k + offset
        out += sp.Rational(coefficient) * falling(r)

    return sp.expand(out)


# ==============================================================================
# 3. RECONSTRUCT ORIGINAL POLYNOMIALS
# ==============================================================================

A_poly = [
    reconstruct_row(row, k)
    for k, row in enumerate(A)
]

B_poly = [
    reconstruct_row(row, k)
    for k, row in enumerate(B)
]


# ==============================================================================
# 4. HEADER
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 99 — EXACT TWO-SIDED ENDPOINT LAW / RESIDUAL KERNEL AUDIT")
print("=" * 78)


# ==============================================================================
# 5. CANDIDATE MULTIPLICITY LAW
# ==============================================================================

print()
print("=" * 78)
print("1. CANDIDATE ENDPOINT MULTIPLICITY LAW")
print("=" * 78)

candidate_s = {
    "A": lambda k: max(0, k - 4),
    "B": lambda k: max(0, k - 3),
}

for name, channel in [
    ("A", A_poly),
    ("B", B_poly),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel, m={m}")

    for k in range(len(channel)):
        s = candidate_s[name](k)

        print(
            f"    k={k}: predicted s(k)={s}"
        )


# ==============================================================================
# 6. EXACT VALIDATION OF THE MULTIPLICITY LAW
# ==============================================================================

print()
print("=" * 78)
print("2. EXACT VALIDATION OF s(k)")
print("=" * 78)

law_ok = True
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
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s = candidate_s[name](k)

        lower = falling(k)
        upper = endpoint_falling(m, s)

        divisor = sp.expand(
            lower * upper
        )

        R = exact_divide(
            divisor,
            P,
        )

        ok = R is not None

        law_ok = law_ok and ok

        print(
            f"    k={k}:"
            f" lower_s={k}"
            f" upper_s={s}"
            f" exact={ok}"
        )

        if ok:
            residuals[name].append(R)

            print(
                f"      residual_degree={degree(R)}"
            )
        else:
            residuals[name].append(None)


# ==============================================================================
# 7. MAXIMALITY CHECK
# ==============================================================================

print()
print("=" * 78)
print("3. MAXIMALITY OF THE ENDPOINT FACTOR")
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

        s = candidate_s[name](k)

        lower = falling(k)
        upper = endpoint_falling(m, s)

        current = exact_divide(
            lower * upper,
            P,
        )

        if current is None:
            print(
                f"    k={k}: candidate decomposition FAILED"
            )
            maximality_ok = False
            continue

        next_upper = endpoint_falling(
            m,
            s + 1,
        )

        stronger = exact_divide(
            lower * next_upper,
            P,
        )

        maximal = stronger is None

        print(
            f"    k={k}:"
            f" s={s}"
            f" maximal={maximal}"
        )

        maximality_ok = (
            maximality_ok and maximal
        )


# ==============================================================================
# 8. RESIDUAL DEGREE LAW
# ==============================================================================

print()
print("=" * 78)
print("4. RESIDUAL DEGREE LAW")
print("=" * 78)

degree_ok = True

for name in ["A", "B"]:

    channel = residuals[name]

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    observed = []

    for k, R in enumerate(channel):

        if R is None:
            degree_ok = False
            continue

        d = degree(R)

        predicted = (
            m
            - k
            - candidate_s[name](k)
        )

        observed.append(d)

        ok = d == predicted

        print(
            f"    k={k}:"
            f" observed={d}"
            f" predicted={predicted}"
            f" exact={ok}"
        )

        degree_ok = degree_ok and ok


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
            f"    R_{k}(j) ="
            f" {sp.factor(R)}"
        )


# ==============================================================================
# 10. RESIDUAL INTEGER-ZERO AUDIT
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

        zeros = []

        for q in range(m + 1):
            if clean(R.subs(j, q)) == 0:
                zeros.append(q)

        print(
            f"    k={k}: residual_zeros={zeros}"
        )


# ==============================================================================
# 11. RESIDUAL REFLECTION AUDIT
# ==============================================================================

print()
print("=" * 78)
print("7. RESIDUAL REFLECTION AUDIT")
print("=" * 78)

reflection_ok = True

for name in ["A", "B"]:

    m = 8 if name == "A" else 7

    print()
    print(
        f"  {name}: reflection j -> {m}-j"
    )

    for k, R in enumerate(residuals[name]):

        if R is None:
            reflection_ok = False
            continue

        reflected = clean(
            R.subs(
                j,
                m - j,
            )
        )

        p = poly(R)
        q = poly(reflected)

        g = sp.gcd(p, q)

        proportional = False
        scalar = None

        if degree(R) == degree(reflected):
            lc1 = p.LC()
            lc2 = q.LC()

            if lc1 != 0:
                scalar = sp.cancel(
                    lc2 / lc1
                )

                proportional = (
                    clean(
                        reflected
                        - scalar * R
                    ) == 0
                )

        print(
            f"    k={k}:"
            f" gcd_degree={g.degree()}"
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

    print()
    print(f"  {name} channel")

    channel = residuals[name]

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
# 13. RESIDUAL VALUE NORMALIZATION
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

        lc = p.LC()
        const = p.TC()

        print(
            f"    k={k}:"
            f" leading={lc}"
            f" constant={const}"
        )


# ==============================================================================
# 14. FULL RECONSTRUCTION
# ==============================================================================

print()
print("=" * 78)
print("10. FULL TWO-SIDED RECONSTRUCTION")
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

        s = candidate_s[name](k)
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

        ok = (
            clean(
                rebuilt - P
            ) == 0
        )

        print(
            f"    k={k}: exact={ok}"
        )

        reconstruction_ok = (
            reconstruction_ok and ok
        )


# ==============================================================================
# 15. STRUCTURAL SUMMARY
# ==============================================================================

print()
print("=" * 78)
print("11. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  Experiment 98 established the corrected decomposition

      P_k(j)
        = j_(k) Q_k(j).

  The new observation is that the high-index rows admit
  additional endpoint factors.

  The candidate law is

      A: s_A(k) = max(0, k-4),
      B: s_B(k) = max(0, k-3).

  Thus the proposed two-sided form is

      A:
        P_k(j)
          = j_(k)
            (8-j)_(max(0,k-4))
            R_k(j)

      B:
        P_k(j)
          = j_(k)
            (7-j)_(max(0,k-3))
            R_k(j).

  The experiment tests this law exactly, including maximality.

  A particularly strong outcome would be:

      1. the multiplicity law is exact;
      2. the residual degree obeys
             deg R_k = m-k-s_k;
      3. residuals acquire a simpler common structure.

  That would mean the coefficient triangle is controlled by
  BOTH endpoints, with a transition near the middle of the
  finite layer range.

  We deliberately do not search for a recurrence here.

  No extrapolation is performed.
  No factorization algorithm is inferred.
  """
)


# ==============================================================================
# 16. FINAL EXACTNESS
# ==============================================================================

failures = (
    int(not law_ok)
    + int(not maximality_ok)
    + int(not degree_ok)
    + int(not reconstruction_ok)
)

print()
print("=" * 78)
print("12. FINAL EXACTNESS")
print("=" * 78)

print(
    "  endpoint_index_law =",
    law_ok,
)

print(
    "  endpoint_maximality =",
    maximality_ok,
)

print(
    "  residual_degree_law =",
    degree_ok,
)

print(
    "  full_reconstruction =",
    reconstruction_ok,
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
print("EXPERIMENT 99 COMPLETE")

