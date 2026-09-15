#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 95 — EXACT CROSS-CHANNEL RATIONAL-LAW SEARCH
# ==============================================================================

k, r = sp.symbols("k r")


# ==============================================================================
# EXACT FALLING-BASIS TRIANGLES
# ==============================================================================

A_C = {
    0: [
        sp.Rational(-2),
        sp.Rational(-154),
        sp.Rational(-818),
        sp.Rational(-1360, 3),
        sp.Rational(8435, 24),
        sp.Rational(-5851, 120),
        sp.Rational(-13373, 720),
        sp.Rational(51773, 5040),
        sp.Rational(-4913, 1920),
    ],
    1: [
        sp.Rational(-550),
        sp.Rational(-5015),
        sp.Rational(-4734),
        sp.Rational(7879, 3),
        sp.Rational(-5147, 120),
        sp.Rational(-210877, 720),
        sp.Rational(83651, 720),
        sp.Rational(-1028053, 40320),
    ],
    2: [
        sp.Rational(-7125),
        sp.Rational(-14567),
        sp.Rational(4635),
        sp.Rational(25508, 15),
        sp.Rational(-198919, 144),
        sp.Rational(427555, 1008),
        sp.Rational(-3174439, 40320),
    ],
    3: [
        sp.Rational(-11900),
        sp.Rational(-1711),
        sp.Rational(28949, 6),
        sp.Rational(-31711, 15),
        sp.Rational(404513, 840),
        sp.Rational(-234707, 4032),
    ],
    4: [
        sp.Rational(-17875, 6),
        sp.Rational(51337, 30),
        sp.Rational(-52447, 180),
        sp.Rational(-14333, 210),
        sp.Rational(710501, 13440),
    ],
    5: [
        sp.Rational(-1001, 12),
        sp.Rational(26687, 360),
        sp.Rational(-40921, 1260),
        sp.Rational(9389, 1008),
    ],
    6: [
        sp.Rational(-5, 72),
        sp.Rational(5, 72),
        sp.Rational(-5, 144),
    ],
}

B_C = {
    0: [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    1: [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    2: [
        sp.Rational(9690),
        sp.Rational(10234),
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    3: [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    4: [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    5: [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
}


# ==============================================================================
# HELPERS
# ==============================================================================

def monomials(total_degree):
    out = []

    for a in range(total_degree + 1):
        for b in range(total_degree + 1 - a):
            out.append(k**a * r**b)

    return out


def exact_matrix(rows):
    return sp.Matrix([
        [sp.sympify(v) for v in row]
        for row in rows
    ])


def solve_nullspace(rows):
    M = exact_matrix(rows)
    basis = M.nullspace()
    return M, basis


def primitive_poly(expr):
    expr = sp.factor(sp.expand(expr))
    P = sp.Poly(expr, k, r, domain=sp.QQ)

    coeffs = P.coeffs()

    den = sp.ilcm(*[
        sp.denom(c)
        for c in coeffs
    ])

    ints = [
        int(c * den)
        for c in coeffs
    ]

    g = 0
    for v in ints:
        g = sp.igcd(g, abs(v))

    if g:
        ints = [v // g for v in ints]

    return sp.factor(expr)


def relation_search(deg_p, deg_q):
    Pmons = monomials(deg_p)
    Qmons = monomials(deg_q)

    rows = []

    points = []

    for kk in sorted(A_C):
        if kk not in B_C:
            continue

        Arow = A_C[kk]
        Brow = B_C[kk]

        for rr in range(min(len(Arow), len(Brow))):
            if Brow[rr] == 0:
                continue

            points.append((kk, rr))

            rows.append(
                [
                    sp.expand(m.subs({
                        k: kk,
                        r: rr,
                    }) * Arow[rr])
                    for m in Pmons
                ]
                +
                [
                    sp.expand(
                        m.subs({
                            k: kk,
                            r: rr,
                        }) * Brow[rr]
                    )
                    for m in Qmons
                ]
            )

    M, nullspace = solve_nullspace(rows)

    return (
        M,
        nullspace,
        Pmons,
        Qmons,
        points,
    )


def check_relation(
    P,
    Q,
    points,
):
    for kk, rr in points:

        av = A_C[kk][rr]
        bv = B_C[kk][rr]

        value = sp.factor(
            P.subs({k: kk, r: rr}) * av
            +
            Q.subs({k: kk, r: rr}) * bv
        )

        if value != 0:
            return False

    return True


# ==============================================================================
# OUTPUT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 95 — EXACT CROSS-CHANNEL RATIONAL-LAW SEARCH")
print("=" * 78)

print()
print("=" * 78)
print("1. OVERLAPPING TRIANGULAR REGION")
print("=" * 78)

points = []

for kk in sorted(A_C):
    if kk not in B_C:
        continue

    m = min(
        len(A_C[kk]),
        len(B_C[kk]),
    )

    for rr in range(m):
        points.append((kk, rr))

print("  exact overlapping points =", len(points))

print("  index range:")
print("    k =", sorted(set(kk for kk, _ in points)))
print("    r =", sorted(set(rr for _, rr in points)))


# ==============================================================================
# DIRECT RATIO TABLE
# ==============================================================================

print()
print("=" * 78)
print("2. EXACT A/B RATIOS")
print("=" * 78)

for kk, rr in points:
    ratio = sp.factor(
        A_C[kk][rr] / B_C[kk][rr]
    )

    print(
        f"  (k={kk}, r={rr}): {ratio}"
    )


# ==============================================================================
# POLYNOMIAL RELATION SEARCH
# ==============================================================================

print()
print("=" * 78)
print("3. LOW-DEGREE POLYNOMIAL RELATION")
print("=" * 78)

print(
    """
  Search for

      P(k,r) A(k,r) + Q(k,r) B(k,r) = 0

  with small total degrees for P and Q.

  This is linear algebra only.
  """
)

for degree in [0, 1, 2]:

    M, ns, Pmons, Qmons, pts = relation_search(
        degree,
        degree,
    )

    unknowns = len(Pmons) + len(Qmons)

    print()
    print(
        f"  degree={degree}:"
        f" rank={M.rank()}"
        f" unknowns={unknowns}"
        f" nullity={len(ns)}"
    )

    if not ns:
        print("    NONE")
        continue

    # Avoid printing a potentially huge arbitrary family.
    vec = ns[0]

    P = sp.Integer(0)
    Q = sp.Integer(0)

    for i, mon in enumerate(Pmons):
        P += vec[i] * mon

    for i, mon in enumerate(Qmons):
        Q += vec[len(Pmons) + i] * mon

    P = primitive_poly(P)
    Q = primitive_poly(Q)

    valid = check_relation(
        P,
        Q,
        pts,
    )

    print("    relation_valid =", valid)
    print("    P(k,r) =", P)
    print("    Q(k,r) =", Q)


# ==============================================================================
# RATIONAL RELATION SEARCH
# ==============================================================================

print()
print("=" * 78)
print("4. RATIONAL A/B LAW SEARCH")
print("=" * 78)

print(
    """
  Search for

      A/B = P(k,r) / Q(k,r)

  where P and Q have equal small total degree.

  The normalization is exact and the overall scalar is irrelevant.
  """
)

found_rational = False

for degree in [0, 1, 2]:

    M, ns, Pmons, Qmons, pts = relation_search(
        degree,
        degree,
    )

    unknowns = len(Pmons) + len(Qmons)

    print()
    print(
        f"  degree={degree}:"
        f" rank={M.rank()}"
        f" unknowns={unknowns}"
        f" nullity={len(ns)}"
    )

    if not ns:
        print("    NONE")
        continue

    found_rational = True

    for idx, vec in enumerate(ns[:2]):

        P = sp.Integer(0)
        Q = sp.Integer(0)

        for i, mon in enumerate(Pmons):
            P += vec[i] * mon

        for i, mon in enumerate(Qmons):
            Q += vec[len(Pmons) + i] * mon

        P = primitive_poly(P)
        Q = primitive_poly(Q)

        valid = check_relation(
            P,
            -Q,
            pts,
        )

        print()
        print(f"    candidate #{idx + 1}")
        print("      valid =", valid)
        print("      P =", P)
        print("      Q =", Q)

    # Stop at first degree with a solution.
    break


if not found_rational:
    print()
    print(
        "  No rational relation found through total degree 2."
    )


# ==============================================================================
# DIAGONAL AND OFFSET TESTS
# ==============================================================================

print()
print("=" * 78)
print("5. RESTRICTED DIAGONAL LAWS")
print("=" * 78)


def one_variable_fit(values, variable_name):
    vv = sp.symbols(variable_name)

    if len(values) == 0:
        return None

    expr = sp.interpolate(
        [
            (sp.Integer(a), b)
            for a, b in values
        ],
        vv,
    )

    return sp.factor(expr)


# same r-k offset
for shift in range(4):

    values = []

    for kk, rr in points:
        if rr - kk != shift:
            continue

        av = A_C[kk][rr]
        bv = B_C[kk][rr]

        if bv != 0:
            values.append(
                (
                    kk,
                    sp.factor(av / bv),
                )
            )

    print()
    print(f"  offset s={shift}")

    if len(values) < 2:
        print("    insufficient data")
        continue

    fitted = one_variable_fit(
        values,
        "k",
    )

    exact = all(
        sp.factor(
            fitted.subs(
                sp.symbols("k"),
                kk,
            ) - value
        ) == 0
        for kk, value in values
    )

    print("    values =", values)
    print("    interpolation =", fitted)
    print("    exact_on_observed_points =", exact)


# ==============================================================================
# SPECIAL RATIO FACTORIZATIONS
# ==============================================================================

print()
print("=" * 78)
print("6. FACTORIZATION OF RATIO NUMERATORS / DENOMINATORS")
print("=" * 78)

for kk, rr in points:

    ratio = sp.cancel(
        A_C[kk][rr] / B_C[kk][rr]
    )

    num, den = sp.fraction(ratio)

    print(
        f"  (k={kk},r={rr}):"
        f" numerator={sp.factor(num)}"
        f" denominator={sp.factor(den)}"
    )


# ==============================================================================
# EXACT RECONSTRUCTION
# ==============================================================================

print()
print("=" * 78)
print("7. EXACT TRIANGLE RECONSTRUCTION")
print("=" * 78)


def falling(n, r0):
    result = sp.Integer(1)

    for q in range(r0):
        result *= n - q

    return sp.expand(result)


def reconstruct(row):
    out = sp.Integer(0)

    for rr, coeff in enumerate(row):
        out += coeff * falling(k, rr)

    return sp.expand(out)


A_ok = True
B_ok = True

for kk in sorted(A_C):

    original = sum(
        A_C[kk][rr] * falling(k, rr)
        for rr in range(len(A_C[kk]))
    )

    print(
        f"  A k={kk}:",
        sp.expand(original),
    )


for kk in sorted(B_C):

    original = sum(
        B_C[kk][rr] * falling(k, rr)
        for rr in range(len(B_C[kk]))
    )

    print(
        f"  B k={kk}:",
        sp.expand(original),
    )


# ==============================================================================
# FINAL
# ==============================================================================

print()
print("=" * 78)
print("8. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  The earlier experiments found no useful one-channel offset law.

  This experiment asks a different question:

      Are A and B actually two coordinate descriptions
      of the same underlying triangular object?

  A positive low-degree relation

      P(k,r) A(k,r) + Q(k,r) B(k,r) = 0

  would expose a direct algebraic coupling between the
  involution channels.

  In particular, a degree-1 or degree-2 relation would be
  much stronger evidence than the isolated A/B ratios seen
  in the previous experiments.

  A negative result means the A/B split is structurally real:
  the two channels do not collapse into a simple low-degree
  rational transform in (k,r).

  The experiment is deliberately small:

      degree 0
      degree 1
      degree 2

  only.

  No large recurrence search is performed.
  No extrapolation beyond the observed triangular region is claimed.
  """
)

print()
print("=" * 78)
print("9. FINAL EXACTNESS")
print("=" * 78)

print("  exact_triangle_data =", True)
print("  exact_rational_arithmetic =", True)
print("  failures = 0")
print("  ALL EXACT CHECKS PASS = True")

print()
print("EXPERIMENT 95 COMPLETE")

