#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 94 — EXACT FULL OFFSET-GENERATOR FACTORIZATION AUDIT
# ==============================================================================

u, j, s, x = sp.symbols("u j s x")


# ==============================================================================
# HELPERS
# ==============================================================================

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def poly(expr, var):
    return sp.Poly(clean(expr), var, domain=sp.QQ)


def degree(expr, var):
    p = poly(expr, var)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def falling(n, r):
    out = sp.Integer(1)
    for q in range(r):
        out *= n - q
    return clean(out)


def coefficient_values(sequence, k):
    points = []

    for jj, expr in enumerate(sequence):
        p = poly(expr, u)
        value = p.coeff_monomial(u ** k)
        points.append((sp.Integer(jj), value))

    return clean(sp.interpolate(points, j))


def falling_coefficients(P):
    deg = int(degree(P, j))

    values = [
        clean(P.subs(j, q))
        for q in range(deg + 1)
    ]

    coeffs = {}
    current = values

    for r in range(deg + 1):
        coeffs[r] = clean(current[0] / sp.factorial(r))

        if len(current) == 1:
            break

        current = [
            clean(current[q + 1] - current[q])
            for q in range(len(current) - 1)
        ]

    return coeffs


def build_triangle(sequence, max_k):
    out = {}

    for k in range(max_k + 1):
        P = coefficient_values(sequence, k)
        out[k] = falling_coefficients(P)

    return out


def primitive_integer_coeffs(expr, var):
    p = poly(expr, var)

    coeffs = p.all_coeffs()

    denoms = [sp.denom(c) for c in coeffs]
    lcm_den = sp.ilcm(*denoms)

    ints = [
        int(c * lcm_den)
        for c in coeffs
    ]

    g = 0
    for q in ints:
        g = sp.igcd(g, abs(q))

    if g == 0:
        return ints

    ints = [q // g for q in ints]

    # normalize sign
    for q in ints:
        if q != 0:
            if q < 0:
                ints = [-v for v in ints]
            break

    return ints


def integer_roots_in_range(expr, var, lo, hi):
    out = []

    for n in range(lo, hi + 1):
        if clean(expr.subs(var, n)) == 0:
            out.append(n)

    return out


# ==============================================================================
# EXACT CHANNEL DATA
# ==============================================================================

A = {
    0: -2,

    1: -2 * (275*u + 78),

    2: -2 * (7125*u**2 + 5565*u + 973),

    3: -4 * (
        17850*u**3
        + 32538*u**2
        + 15036*u
        + 2023
    ),

    4: -(
        71500*u**4
        + 326664*u**3
        + 323868*u**2
        + 112964*u
        + 12879
    ),

    5: -1001 * (
        10*u**5
        + 152*u**4
        + 340*u**3
        + 256*u**2
        + 77*u
        + 8
    ),

    6: -(
        50*u**6
        + 6686*u**5
        + 50200*u**4
        + 92208*u**3
        + 63401*u**2
        + 18109*u
        + 1820
    ),

    7: -(
        276*u**5
        + 2700*u**4
        + 5478*u**3
        + 3966*u**2
        + 1169*u
        + 120
    ),

    8: -(
        9*u**4
        + 30*u**3
        + 27*u**2
        + 9*u
        + 1
    ),
}


B = {
    0: 25,

    1: 14 * (125*u + 46),

    2: 6 * (
        3230*u**2
        + 3458*u
        + 749
    ),

    3: 4 * (
        11050*u**3
        + 29886*u**2
        + 17663*u
        + 2869
    ),

    4: (
        17875*u**4
        + 138710*u**3
        + 188409*u**2
        + 83097*u
        + 11441
    ),

    5: 13 * (
        50*u**5
        + 1784*u**4
        + 6480*u**3
        + 6812*u**2
        + 2639*u
        + 336
    ),

    6: (
        144*u**5
        + 3370*u**4
        + 11332*u**3
        + 11589*u**2
        + 4431*u
        + 560
    ),

    7: (
        44*u**4
        + 230*u**3
        + 282*u**2
        + 119*u
        + 16
    ),
}


# ==============================================================================
# BUILD EXACT FALLING-BASIS TRIANGLES
# ==============================================================================

Aseq = [A[q] for q in range(9)]
Bseq = [B[q] for q in range(8)]

A_C = build_triangle(Aseq, 6)
B_C = build_triangle(Bseq, 5)


# ==============================================================================
# FULL OFFSET GENERATORS
# ==============================================================================

def diagonal_normalized_generator(C, k, max_r):
    """
    G_k(x) =
        sum_s C[k,k+s]/C[k,k] * x^s
    """
    d = clean(C[k][k])

    terms = []

    for r in range(k, max_r + 1):
        c = clean(C[k].get(r, 0))

        if d == 0:
            value = sp.Integer(0)
        else:
            value = clean(c / d)

        terms.append(value * x ** (r - k))

    return clean(sum(terms))


def binomial_removed_generator(C, k, max_r):
    """
    F_k(x) =
      sum_s
        C[k,k+s]
        -------------------------------- x^s
        C[k,k] * binom(k+s,k)
    """
    d = clean(C[k][k])

    terms = []

    for r in range(k, max_r + 1):
        c = clean(C[k].get(r, 0))

        if d == 0:
            value = sp.Integer(0)
        else:
            value = clean(
                c / (
                    d
                    * sp.binomial(r, k)
                )
            )

        terms.append(value * x ** (r - k))

    return clean(sum(terms))


# ==============================================================================
# 1. FULL GENERATORS
# ==============================================================================

print("=" * 78)
print(
    "EXPERIMENT 94 — EXACT FULL OFFSET-GENERATOR "
    "FACTORIZATION AUDIT"
)
print("=" * 78)

print()
print("=" * 78)
print("1. FULL DIAGONAL-NORMALIZED OFFSET GENERATORS")
print("=" * 78)

A_diag_generators = {}
B_diag_generators = {}

for k in range(7):
    G = diagonal_normalized_generator(
        A_C,
        k,
        8,
    )

    A_diag_generators[k] = G

    print()
    print(f"  A k={k}")
    print("    degree =", degree(G, x))
    print("    generator =", sp.factor(G))
    print(
        "    primitive_integer_coefficients =",
        primitive_integer_coeffs(G, x)
    )
    print(
        "    integer_roots_in_0..8 =",
        integer_roots_in_range(G, x, 0, 8)
    )


for k in range(6):
    G = diagonal_normalized_generator(
        B_C,
        k,
        7,
    )

    B_diag_generators[k] = G

    print()
    print(f"  B k={k}")
    print("    degree =", degree(G, x))
    print("    generator =", sp.factor(G))
    print(
        "    primitive_integer_coefficients =",
        primitive_integer_coeffs(G, x)
    )
    print(
        "    integer_roots_in_0..7 =",
        integer_roots_in_range(G, x, 0, 7)
    )


# ==============================================================================
# 2. BINOMIAL-REMOVED GENERATORS
# ==============================================================================

print()
print("=" * 78)
print("2. FULL BINOMIAL-REMOVED OFFSET GENERATORS")
print("=" * 78)

A_bin_generators = {}
B_bin_generators = {}

for k in range(7):
    G = binomial_removed_generator(
        A_C,
        k,
        8,
    )

    A_bin_generators[k] = G

    print()
    print(f"  A k={k}")
    print("    degree =", degree(G, x))
    print("    generator =", sp.factor(G))
    print(
        "    primitive_integer_coefficients =",
        primitive_integer_coeffs(G, x)
    )
    print(
        "    integer_roots_in_0..8 =",
        integer_roots_in_range(G, x, 0, 8)
    )


for k in range(6):
    G = binomial_removed_generator(
        B_C,
        k,
        7,
    )

    B_bin_generators[k] = G

    print()
    print(f"  B k={k}")
    print("    degree =", degree(G, x))
    print("    generator =", sp.factor(G))
    print(
        "    primitive_integer_coefficients =",
        primitive_integer_coeffs(G, x)
    )
    print(
        "    integer_roots_in_0..7 =",
        integer_roots_in_range(G, x, 0, 7)
    )


# ==============================================================================
# 3. COMMON FACTOR SEARCH BETWEEN ROW GENERATORS
# ==============================================================================

print()
print("=" * 78)
print("3. COMMON FACTORS BETWEEN CONSECUTIVE ROW GENERATORS")
print("=" * 78)


def print_gcd_chain(generators, name):
    print()
    print(f"  {name}")

    for k in range(len(generators) - 1):

        g = sp.gcd(
            poly(generators[k], x),
            poly(generators[k + 1], x),
        ).as_expr()

        g = clean(g)

        print(
            f"    rows {k},{k+1}: gcd = {sp.factor(g)}"
        )


print_gcd_chain(
    A_diag_generators,
    "A diagonal-normalized"
)

print_gcd_chain(
    B_diag_generators,
    "B diagonal-normalized"
)

print_gcd_chain(
    A_bin_generators,
    "A binomial-removed"
)

print_gcd_chain(
    B_bin_generators,
    "B binomial-removed"
)


# ==============================================================================
# 4. SAME-SHAPE / AFFINE-COMPOSITION SEARCH
# ==============================================================================

print()
print("=" * 78)
print("4. AFFINE-COMPOSITION SEARCH")
print("=" * 78)

print(
    """
  Test whether neighboring generators have the form

      G_{k+1}(x)
        = c_k * G_k(a_k*x + b_k)

  for simple rational a_k,b_k,c_k.

  This is intentionally restricted to a tiny finite candidate set.
  """
)


candidate_values = [
    sp.Rational(-2),
    sp.Rational(-1),
    sp.Rational(-1, 2),
    sp.Rational(1, 2),
    sp.Rational(1),
    sp.Rational(2),
]


def affine_match(G1, G2):
    matches = []

    for a in candidate_values:
        for b in candidate_values:

            transformed = clean(
                G1.subs(x, a*x + b)
            )

            p2 = poly(G2, x)

            if p2.is_zero:
                continue

            lc1 = poly(transformed, x).LC()
            lc2 = p2.LC()

            if lc1 == 0:
                continue

            c = clean(lc2 / lc1)

            test = clean(
                G2 - c * transformed
            )

            if test == 0:
                matches.append(
                    (a, b, c)
                )

    return matches


for name, generators in (
    ("A diagonal", A_diag_generators),
    ("B diagonal", B_diag_generators),
    ("A binomial", A_bin_generators),
    ("B binomial", B_bin_generators),
):

    print()
    print(f"  {name}")

    for k in range(len(generators) - 1):

        matches = affine_match(
            generators[k],
            generators[k + 1],
        )

        print(
            f"    {k}->{k+1}: {matches}"
        )


# ==============================================================================
# 5. REVERSAL / PALINDROME AUDIT
# ==============================================================================

print()
print("=" * 78)
print("5. RECIPROCAL / PALINDROMIC AUDIT")
print("=" * 78)


def reciprocal_status(G):
    p = poly(G, x)

    d = p.degree()

    rev = clean(
        x ** d * G.subs(x, 1/x)
    )

    if clean(rev - G) == 0:
        return "reciprocal"

    if clean(rev + G) == 0:
        return "anti-reciprocal"

    return "none"


for name, generators in (
    ("A diagonal", A_diag_generators),
    ("B diagonal", B_diag_generators),
    ("A binomial", A_bin_generators),
    ("B binomial", B_bin_generators),
):

    print()
    print(f"  {name}")

    for k, G in generators.items():
        print(
            f"    k={k}: {reciprocal_status(G)}"
        )


# ==============================================================================
# 6. EVALUATION AT SPECIAL OFFSET VALUES
# ==============================================================================

print()
print("=" * 78)
print("6. SPECIAL-OFFSET EVALUATION")
print("=" * 78)


special_points = [
    sp.Integer(-2),
    sp.Integer(-1),
    sp.Integer(0),
    sp.Integer(1),
    sp.Integer(2),
]


for name, generators in (
    ("A diagonal", A_diag_generators),
    ("B diagonal", B_diag_generators),
    ("A binomial", A_bin_generators),
    ("B binomial", B_bin_generators),
):

    print()
    print(f"  {name}")

    for k, G in generators.items():

        values = [
            (q, clean(G.subs(x, q)))
            for q in special_points
        ]

        print(
            f"    k={k}: {values}"
        )


# ==============================================================================
# 7. CROSS-CHANNEL GENERATOR GCD
# ==============================================================================

print()
print("=" * 78)
print("7. CROSS-CHANNEL GENERATOR GCD AUDIT")
print("=" * 78)


for k in range(6):

    GA = A_diag_generators[k]
    GB = B_diag_generators[k]

    g = sp.gcd(
        poly(GA, x),
        poly(GB, x),
    ).as_expr()

    print(
        f"  diagonal k={k}: gcd={sp.factor(clean(g))}"
    )


for k in range(6):

    GA = A_bin_generators[k]
    GB = B_bin_generators[k]

    g = sp.gcd(
        poly(GA, x),
        poly(GB, x),
    ).as_expr()

    print(
        f"  binomial-removed k={k}: gcd={sp.factor(clean(g))}"
    )


# ==============================================================================
# 8. EXACT RECONSTRUCTION
# ==============================================================================

print()
print("=" * 78)
print("8. EXACT ORIGINAL TRIANGLE RECONSTRUCTION")
print("=" * 78)


def reconstruct_from_triangle(C, k, max_r):
    return clean(
        sum(
            C[k].get(r, sp.Integer(0))
            * falling(j, r)
            for r in range(max_r + 1)
        )
    )


all_ok = True


for k in range(7):

    P = coefficient_values(
        Aseq,
        k,
    )

    R = reconstruct_from_triangle(
        A_C,
        k,
        8,
    )

    ok = clean(P - R) == 0
    all_ok &= ok

    print(f"  A k={k}: {ok}")


for k in range(6):

    P = coefficient_values(
        Bseq,
        k,
    )

    R = reconstruct_from_triangle(
        B_C,
        k,
        7,
    )

    ok = clean(P - R) == 0
    all_ok &= ok

    print(f"  B k={k}: {ok}")


# ==============================================================================
# 9. INTERPRETATION
# ==============================================================================

print()
print("=" * 78)
print("9. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  The preceding experiments found

      P_k(j)
        = sum_{r>=k} C[k,r] j_(r).

  Experiment 93 showed that even after dividing by

      C[k,k] * binom(k+s,k),

  the remaining offset matrix still has rank 4 in both channels.

  The present experiment therefore changes the question.

  Instead of asking for separability, we form the complete polynomial

      G_k(x)
        = sum_s C[k,k+s]/C[k,k] * x^s

  and its binomial-removed version

      F_k(x)
        = sum_s
            C[k,k+s]
            ----------------------------- x^s
            C[k,k] * binom(k+s,k).

  We then test for:

      * factorization in x;
      * common factors between neighboring rows;
      * affine composition between rows;
      * reciprocal/palindromic structure;
      * special-offset values;
      * cross-channel common factors.

  A positive result here would be substantially more informative
  than another generic interpolation identity.

  In particular, a relation such as

      G_{k+1}(x)
        = c_k G_k(a_k x + b_k)

  or a persistent polynomial factor would indicate that the
  triangular coefficient family is generated by a genuine
  transformation in the offset variable.

  If these tests are all negative, the evidence increasingly points
  toward the falling-factorial triangle itself being the natural
  finite construction object rather than a disguised one-variable
  recurrence.

  Everything remains exact over QQ.
  No extrapolation beyond the known range is performed.
  """
)

print()
print("=" * 78)
print("10. FINAL EXACTNESS")
print("=" * 78)

print(
    "  triangle reconstruction =",
    all_ok
)

print(
    "  all checks pass =",
    all_ok
)

print()
print("EXPERIMENT 94 COMPLETE")
