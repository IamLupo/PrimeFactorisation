#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 97 — EXACT TWO-SIDED FALLING-FACTORIAL / ENDPOINT AUDIT
#
# Exact QQ arithmetic only.
# No floating point.
# ==============================================================================

j = sp.symbols("j")


# ==============================================================================
# 1. EXACT FALLING-BASIS TRIANGLES
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
# 2. SAFE EXACT HELPERS
# ==============================================================================

def clean(expr):
    if expr is None:
        return None
    return sp.cancel(sp.expand(sp.sympify(expr)))


def is_zero(expr):
    if expr is None:
        return False
    return clean(expr) == 0


def poly_qq(expr):
    if expr is None:
        raise ValueError("Cannot create polynomial from None.")

    return sp.Poly(
        sp.expand(sp.sympify(expr)),
        j,
        domain=sp.QQ,
    )


def degree(expr):
    if expr is None:
        return -sp.oo

    expr = clean(expr)

    if expr == 0:
        return -sp.oo

    return poly_qq(expr).degree()


def falling(n):
    if n == 0:
        return sp.Integer(1)

    out = sp.Integer(1)

    for r in range(n):
        out *= j - r

    return sp.expand(out)


def endpoint_falling(m, s):
    """
    (m-j)_under_s
      = (m-j)(m-j-1)...(m-j-s+1)
    """
    if s == 0:
        return sp.Integer(1)

    out = sp.Integer(1)

    for r in range(s):
        out *= m - j - r

    return sp.expand(out)


def exact_division(divisor, dividend):
    """
    Returns:
        quotient polynomial if exact,
        0 if dividend is identically zero,
        None otherwise.
    """
    divisor = clean(divisor)
    dividend = clean(dividend)

    if divisor is None or dividend is None:
        return None

    if divisor == 0:
        return None

    if dividend == 0:
        return sp.Integer(0)

    p_divisor = poly_qq(divisor)
    p_dividend = poly_qq(dividend)

    q, r = sp.div(
        p_dividend,
        p_divisor,
        domain=sp.QQ,
    )

    if r.is_zero:
        return sp.expand(q.as_expr())

    return None


def triangle_row_to_poly(row):
    out = sp.Integer(0)

    for r, coeff in enumerate(row):
        out += sp.Rational(coeff) * falling(r)

    return sp.expand(out)


def original_falling_quotient(P, k):
    return exact_division(
        falling(k),
        P,
    )


# ==============================================================================
# 3. BUILD POLYNOMIALS
# ==============================================================================

A_polys = [
    triangle_row_to_poly(row)
    for row in A
]

B_polys = [
    triangle_row_to_poly(row)
    for row in B
]


# ==============================================================================
# 4. MAXIMUM ENDPOINT FACTOR
# ==============================================================================

def max_endpoint_factor(P, k, m):
    """
    Find the largest s for which

        j_under_k * (m-j)_under_s

    divides P.

    This routine is robust to zero residuals.
    """
    Q = original_falling_quotient(P, k)

    if Q is None:
        return None, None

    Q = clean(Q)

    # Zero polynomial is divisible by every polynomial.
    # There is no meaningful finite "maximum" multiplicity.
    # We report the full remaining degree as the usable bound.
    if Q == 0:
        return 0, sp.Integer(0)

    max_possible = degree(Q)

    best_s = 0
    best_R = Q

    for s in range(max_possible + 1):

        divisor = endpoint_falling(m, s)

        candidate = exact_division(
            divisor,
            Q,
        )

        if candidate is None:
            break

        best_s = s
        best_R = clean(candidate)

    return best_s, best_R


# ==============================================================================
# 5. HEADER
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 97 — EXACT TWO-SIDED FALLING-FACTORIAL / ENDPOINT AUDIT")
print("=" * 78)


# ==============================================================================
# 6. ORIGINAL FALLING VALIDATION
# ==============================================================================

print()
print("=" * 78)
print("1. ORIGINAL FALLING-FACTORIAL VALIDATION")
print("=" * 78)

all_falling = True

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel, endpoint m={m}")

    for k, P in enumerate(channel):

        Q = original_falling_quotient(P, k)
        ok = Q is not None

        print(
            f"    k={k}: divisible={ok}"
        )

        all_falling = all_falling and ok


# ==============================================================================
# 7. ENDPOINT MULTIPLICITY
# ==============================================================================

print()
print("=" * 78)
print("2. MAXIMUM ENDPOINT FALLING-FACTOR MULTIPLICITIES")
print("=" * 78)

endpoint_data = {
    "A": [],
    "B": [],
}

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s, R = max_endpoint_factor(
            P,
            k,
            m,
        )

        endpoint_data[name].append(
            (s, R)
        )

        if R is None:
            print(
                f"    k={k}: endpoint analysis unavailable"
            )
            continue

        if is_zero(R):
            print(
                f"    k={k}: residual = 0"
                f"  endpoint_s={s}"
            )
            continue

        print(
            f"    k={k}:"
            f" endpoint_s={s}"
            f" residual_degree={degree(R)}"
        )

        print(
            f"      endpoint_factor="
            f"{endpoint_falling(m, s)}"
        )

        print(
            f"      residual={sp.factor(R)}"
        )


# ==============================================================================
# 8. TWO-SIDED FACTORIZATION
# ==============================================================================

print()
print("=" * 78)
print("3. TWO-SIDED FACTORIZATION PROFILE")
print("=" * 78)

all_two_sided = True

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        s, R = endpoint_data[name][k]

        if R is None:
            print(
                f"    k={k}: unavailable"
            )
            all_two_sided = False
            continue

        total_factor = (
            falling(k)
            * endpoint_falling(m, s)
        )

        Q = exact_division(
            total_factor,
            P,
        )

        ok = Q is not None

        print(
            f"    k={k}:"
            f" j_under_{k} * (m-j)_under_{s}"
            f" exact={ok}"
        )

        all_two_sided = (
            all_two_sided and ok
        )

        if ok:

            print(
                f"      quotient_degree={degree(Q)}"
            )

            print(
                f"      quotient={sp.factor(Q)}"
            )


# ==============================================================================
# 9. DIRECT INTEGER ZERO AUDIT
# ==============================================================================

print()
print("=" * 78)
print("4. DIRECT INTEGER-POINT ZERO AUDIT")
print("=" * 78)

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        left_zeros = []
        right_zeros = []

        for q in range(0, k + 1):

            value = clean(
                P.subs(
                    j,
                    sp.Integer(q),
                )
            )

            if value == 0:
                left_zeros.append(q)

        for q in range(m, max(m - k, -1) - 1, -1):

            value = clean(
                P.subs(
                    j,
                    sp.Integer(q),
                )
            )

            if value == 0:
                right_zeros.append(q)

        print(
            f"    k={k}:"
            f" left_zeros={left_zeros}"
            f" right_zeros={right_zeros}"
        )


# ==============================================================================
# 10. RESIDUAL DEGREE PROFILE
# ==============================================================================

print()
print("=" * 78)
print("5. RESIDUAL DEGREE PROFILE")
print("=" * 78)

for name in ["A", "B"]:

    print()
    print(f"  {name} channel")

    for k, (s, R) in enumerate(
        endpoint_data[name]
    ):

        if R is None:
            print(
                f"    k={k}: unavailable"
            )
            continue

        print(
            f"    k={k}:"
            f" endpoint_s={s}"
            f" residual_degree={degree(R)}"
        )


# ==============================================================================
# 11. CONSECUTIVE RESIDUAL GCD AUDIT
# ==============================================================================

print()
print("=" * 78)
print("6. RESIDUAL FACTORIZATION / COMMON-GCD AUDIT")
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
                f"    residuals {k},{k+1}: unavailable"
            )
            continue

        if is_zero(R1) or is_zero(R2):

            print(
                f"    residuals {k},{k+1}:"
                f" one residual is zero"
            )
            continue

        g = sp.gcd(
            poly_qq(R1),
            poly_qq(R2),
        )

        print(
            f"    residuals {k},{k+1}:"
            f" gcd_degree={g.degree()}"
        )

        if g.degree() > 0:
            print(
                f"      gcd={sp.factor(g.as_expr())}"
            )


# ==============================================================================
# 12. ENDPOINT-COMPLEMENT REFLECTION
# ==============================================================================

print()
print("=" * 78)
print("7. ENDPOINT-COMPLEMENT REFLECTION TEST")
print("=" * 78)

for name, channel in [
    ("A", A_polys),
    ("B", B_polys),
]:

    m = 8 if name == "A" else 7

    print()
    print(
        f"  {name} channel:"
        f" reflection j -> {m}-j"
    )

    for k, P in enumerate(channel):

        reflected = clean(
            P.subs(
                j,
                sp.Integer(m) - j,
            )
        )

        g = sp.gcd(
            poly_qq(P),
            poly_qq(reflected),
        )

        print(
            f"    k={k}:"
            f" gcd_degree={g.degree()}"
        )

        if g.degree() > 0:
            print(
                f"      gcd={sp.factor(g.as_expr())}"
            )


# ==============================================================================
# 13. NORMALIZED RESIDUAL REFLECTION
# ==============================================================================

print()
print("=" * 78)
print("8. NORMALIZED ENDPOINT REFLECTION")
print("=" * 78)


def proportional_reflection(R, m):
    if R is None:
        return False, None

    R = clean(R)

    if R == 0:
        return True, sp.Integer(0)

    reflected = clean(
        R.subs(
            j,
            sp.Integer(m) - j,
        )
    )

    p = poly_qq(R)
    q = poly_qq(reflected)

    ratio = None

    for power in range(
        max(p.degree(), q.degree()) + 1
    ):

        a = sp.Rational(p.nth(power))
        b = sp.Rational(q.nth(power))

        if a == 0 and b == 0:
            continue

        if a == 0 or b == 0:
            return False, None

        current = sp.cancel(
            b / a
        )

        if ratio is None:
            ratio = current
        elif current != ratio:
            return False, None

    return True, ratio


for name in ["A", "B"]:

    m = 8 if name == "A" else 7

    print()
    print(f"  {name} channel")

    for k, (s, R) in enumerate(
        endpoint_data[name]
    ):

        ok, scalar = proportional_reflection(
            R,
            m,
        )

        print(
            f"    k={k}:"
            f" proportional={ok}"
            f" scalar={scalar}"
        )


# ==============================================================================
# 14. TEST WHETHER ENDPOINT MULTIPLICITY HAS A SIMPLE LAW
# ==============================================================================

print()
print("=" * 78)
print("9. ENDPOINT-MULTIPLICITY INDEX PROFILE")
print("=" * 78)

for name in ["A", "B"]:

    data = endpoint_data[name]

    s_values = [
        s
        for s, _ in data
    ]

    print(
        f"  {name}: s(k) = {s_values}"
    )

    # Test simple candidate laws.
    candidates = {
        "0": lambda k, m: 0,
        "k": lambda k, m: k,
        "m-k": lambda k, m: m - k,
        "max(0,k-m)": lambda k, m: max(0, k - m),
        "max(0,2*k-m)": lambda k, m: max(0, 2*k - m),
    }

    m = 8 if name == "A" else 7

    for label, fn in candidates.items():

        matches = all(
            s == fn(k, m)
            for k, (s, _) in enumerate(data)
        )

        print(
            f"    candidate {label}:"
            f" exact={matches}"
        )


# ==============================================================================
# 15. FRESH EXACT CONSISTENCY
# ==============================================================================

print()
print("=" * 78)
print("10. FRESH EXACT SYMBOLIC CONSISTENCY")
print("=" * 78)

all_fresh = True

for name, channel, source_rows in [
    ("A", A_polys, A),
    ("B", B_polys, B),
]:

    print()
    print(f"  {name} channel")

    for k, P in enumerate(channel):

        rebuilt = triangle_row_to_poly(
            source_rows[k]
        )

        ok = clean(
            P - rebuilt
        ) == 0

        print(
            f"    k={k}: exact={ok}"
        )

        all_fresh = (
            all_fresh and ok
        )


# ==============================================================================
# 16. STRUCTURAL INTERPRETATION
# ==============================================================================

print()
print("=" * 78)
print("11. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  The established exact structure is

      P_k(j) = j_under_k Q_k(j).

  Earlier experiments found isolated upper-endpoint factors such as

      A:
          (j-8),
          (j-8)(j-7),

      B:
          (j-7).

  This experiment tests whether those factors belong to a systematic
  second falling-factorial structure

      (m-j)_under_s,

  with

      m = 8  for A,
      m = 7  for B.

  A positive result would mean the original coefficient triangle has
  two-sided endpoint combinatorics:

      P_k(j)
        = j_under_k
          (m-j)_under_s
          R_k(j).

  A simple law for s(k) would be particularly significant.

  A reflection relation for R_k under

      j -> m-j

  would provide an additional symmetry.

  A negative result is also useful: it would indicate that the visible
  endpoint factors are boundary artifacts rather than a complete
  second factorial kernel.

  Everything above is exact over QQ.

  No floating point.
  No large recurrence search.
  No extrapolation.
  No factorization algorithm is inferred.
  """
)


# ==============================================================================
# 17. FINAL EXACTNESS
# ==============================================================================

print()
print("=" * 78)
print("12. FINAL EXACTNESS")
print("=" * 78)

print(
    "  original_falling_divisibility =",
    all_falling,
)

print(
    "  two_sided_factorization_checks =",
    all_two_sided,
)

print(
    "  fresh_consistency =",
    all_fresh,
)

failures = (
    int(not all_falling)
    + int(not all_two_sided)
    + int(not all_fresh)
)

print(
    "  failures =",
    failures,
)

print(
    "  ALL EXACT CHECKS PASS =",
    failures == 0,
)

print()
print("EXPERIMENT 97 COMPLETE")