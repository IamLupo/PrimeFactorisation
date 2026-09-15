#!/usr/bin/env python3

"""
EXPERIMENT 82 — EXACT BIVARIATE INDEX RECONSTRUCTION

Goal
----
Replace the expensive unrestricted indexed-recurrence search by an exact
finite-dimensional reconstruction of the channel coefficients as
polynomials in the layer index j.

We study

    A_j(u) = A_{2j}(u),     j = 0,...,8
    B_j(u) = B_{2j+1}(u),   j = 0,...,7

and reconstruct each coefficient sequence

    [u^k] A_j(u)
    [u^k] B_j(u)

as an exact polynomial in j.

This is NOT interpreted as proof of a universal recurrence beyond the
observed range.  It is an exact interpolation/reconstruction of the
available symbolic data.

No floating point is used.
"""

import sympy as sp
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

u = sp.Symbol("u")
j = sp.Symbol("j")

t = sp.Symbol("t")
N = sp.Symbol("N")
X = sp.Symbol("X")


# ---------------------------------------------------------------------------
# Exact channel data
# ---------------------------------------------------------------------------

# Actual channel polynomials, not primitive-normalized versions.

A = {
    0: sp.Integer(-2),

    1: -2 * (275*u + 78),

    2: -2 * (7125*u**2 + 5565*u + 973),

    3: -4 * (17850*u**3 + 32538*u**2 + 15036*u + 2023),

    4: -(71500*u**4 + 326664*u**3 + 323868*u**2
         + 112964*u + 12879),

    5: -1001 * (10*u**5 + 152*u**4 + 340*u**3
                 + 256*u**2 + 77*u + 8),

    6: -(50*u**6 + 6686*u**5 + 50200*u**4
         + 92208*u**3 + 63401*u**2
         + 18109*u + 1820),

    7: -(276*u**5 + 2700*u**4 + 5478*u**3
         + 3966*u**2 + 1169*u + 120),

    8: -(9*u**4 + 30*u**3 + 27*u**2 + 9*u + 1),
}


B = {
    0: sp.Integer(25),

    1: 14 * (125*u + 46),

    2: 6 * (3230*u**2 + 3458*u + 749),

    3: 4 * (11050*u**3 + 29886*u**2
            + 17663*u + 2869),

    4: (17875*u**4 + 138710*u**3 + 188409*u**2
        + 83097*u + 11441),

    5: 13 * (50*u**5 + 1784*u**4 + 6480*u**3
             + 6812*u**2 + 2639*u + 336),

    6: (144*u**5 + 3370*u**4 + 11332*u**3
        + 11589*u**2 + 4431*u + 560),

    7: (44*u**4 + 230*u**3 + 282*u**2 + 119*u + 16),
}


# ---------------------------------------------------------------------------
# Full normalized layer polynomials h_d(t)
# ---------------------------------------------------------------------------

H = {
    0: sp.Integer(-2),

    1: 25*(2*t + 1),

    2: -2*(275*t**2 + 275*t + 78),

    3: 14*(2*t + 1)*(125*t**2 + 125*t + 46),

    4: -2*(7125*t**4 + 14250*t**3
           + 12690*t**2 + 5565*t + 973),

    5: 6*(2*t + 1)*(3230*t**4 + 6460*t**3
                     + 6688*t**2 + 3458*t + 749),

    6: -4*(17850*t**6 + 53550*t**5
           + 86088*t**4 + 82926*t**3
           + 47574*t**2 + 15036*t + 2023),

    7: 4*(2*t + 1)*(11050*t**6 + 33150*t**5
                     + 63036*t**4 + 70822*t**3
                     + 47549*t**2 + 17663*t + 2869),

    8: -(71500*t**8 + 286000*t**7
         + 755664*t**6 + 1265992*t**5
         + 1375360*t**4 + 974400*t**3
         + 436832*t**2 + 112964*t + 12879),

    9: (2*t + 1)*(17875*t**8 + 71500*t**7
                   + 245960*t**6 + 487630*t**5
                   + 622414*t**4 + 515528*t**3
                   + 271506*t**2 + 83097*t + 11441),

    10: -1001*(10*t**10 + 50*t**9 + 252*t**8
               + 708*t**7 + 1302*t**6 + 1638*t**5
               + 1428*t**4 + 852*t**3 + 333*t**2
               + 77*t + 8),

    11: 13*(2*t + 1)*(50*t**10 + 250*t**9
                      + 2284*t**8 + 7636*t**7
                      + 17434*t**6 + 26626*t**5
                      + 28036*t**4 + 20104*t**3
                      + 9451*t**2 + 2639*t + 336),

    12: -(50*t**12 + 300*t**11 + 7436*t**10
          + 34430*t**9 + 117810*t**8
          + 267960*t**7 + 426888*t**6
          + 484110*t**5 + 390225*t**4
          + 219010*t**3 + 81510*t**2
          + 18109*t + 1820),

    13: (2*t + 1)*(144*t**10 + 720*t**9
                    + 4810*t**8 + 14920*t**7
                    + 32272*t**6 + 47620*t**5
                    + 48955*t**4 + 34510*t**3
                    + 16020*t**2 + 4431*t + 560),

    14: -(276*t**10 + 1380*t**9 + 5460*t**8
          + 13560*t**7 + 23058*t**6
          + 27510*t**5 + 23100*t**4
          + 13410*t**3 + 5135*t**2
          + 1169*t + 120),

    15: (2*t + 1)*(44*t**8 + 176*t**7
                    + 494*t**6 + 866*t**5
                    + 1016*t**4 + 794*t**3
                    + 401*t**2 + 119*t + 16),

    16: -(3*t**2 + 3*t + 1) * (
        3*t**6 + 9*t**5 + 18*t**4
        + 21*t**3 + 15*t**2 + 6*t + 1
    ),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def clean(expr):
    return sp.cancel(sp.expand(expr))


def poly_u(expr):
    return sp.Poly(clean(expr), u, domain=sp.QQ)


def poly_t(expr):
    return sp.Poly(clean(expr), t, domain=sp.QQ)


def degree_u(expr):
    p = poly_u(expr)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def degree_t(expr):
    p = poly_t(expr)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def coefficient_u(expr, k):
    return sp.expand(poly_u(expr).nth(k))


def interpolate_in_j(values):
    """
    Exact interpolation of values indexed by consecutive integers.

    values = [v_0, ..., v_m]
    returns exact QQ polynomial P(j) satisfying P(r)=v_r.
    """
    pts = [(sp.Integer(idx), sp.Rational(val)) for idx, val in enumerate(values)]
    return sp.interpolate(pts, j)


def verify_index_polynomial(poly_j, values):
    for idx, val in enumerate(values):
        lhs = clean(poly_j.subs(j, idx))
        rhs = sp.Rational(val)
        if lhs != rhs:
            return False
    return True


def primitive_integer_coefficients(expr):
    p = poly_u(expr)

    if p.is_zero:
        return []

    coeffs = [sp.Rational(c) for c in p.all_coeffs()]

    den = 1
    for c in coeffs:
        den = sp.ilcm(den, int(c.q))

    ints = [int(c * den) for c in coeffs]

    g = 0
    for x in ints:
        g = sp.igcd(g, abs(x))

    if g == 0:
        return [0] * len(ints)

    ints = [x // g for x in ints]

    if ints[0] < 0:
        ints = [-x for x in ints]

    return ints


def content_qq(expr):
    p = poly_u(expr)
    if p.is_zero:
        return sp.Integer(0)

    coeffs = [sp.Rational(c) for c in p.all_coeffs()]

    den = 1
    for c in coeffs:
        den = sp.ilcm(den, int(c.q))

    ints = [int(c * den) for c in coeffs]

    g = 0
    for x in ints:
        g = sp.igcd(g, abs(x))

    return sp.Rational(g, den)


def finite_difference_degree(values):
    vals = [sp.Rational(v) for v in values]

    if len(vals) <= 1:
        return 0

    current = vals

    for degree in range(len(vals)):
        if len(current) <= 1:
            return degree

        if all(x == 0 for x in current):
            return degree

        current = [
            sp.simplify(current[i + 1] - current[i])
            for i in range(len(current) - 1)
        ]

    return len(vals) - 1


def factor_j(expr):
    return sp.factor(clean(expr))


def factor_ju(expr):
    return sp.factor(
        sp.Poly(clean(expr), j, u, domain=sp.QQ).as_expr()
    )


# ---------------------------------------------------------------------------
# Exact involution decomposition
# ---------------------------------------------------------------------------

def involution_check(h):
    partner = clean(h.subs(t, -1 - t))
    even_part = clean((h + partner) / 2)
    odd_part = clean((h - partner) / 2)

    invariant = clean(partner - h) == 0
    anti = clean(partner + h) == 0

    if invariant:
        return True, False, even_part, sp.Integer(0)

    if anti:
        divisor = 2*t + 1
        q, r = sp.div(
            poly_t(odd_part),
            poly_t(divisor),
            domain=sp.QQ,
        )

        divisible = r == 0
        Bpoly = q.as_expr() if divisible else None

        return False, True, sp.Integer(0), clean(Bpoly)

    return False, False, sp.Integer(0), sp.Integer(0)


def recover_A_B():
    recovered = {}

    for d, h in H.items():
        invariant, anti, A_d, B_d = involution_check(h)

        # The expected channel follows parity.
        if d % 2 == 0:
            A_d = clean(A_d)
            B_d = sp.Integer(0)
        else:
            A_d = sp.Integer(0)
            B_d = clean(B_d)

        reconstructed = clean(
            A_d + (2*t + 1)*B_d - h
        ) == 0

        recovered[d] = {
            "A": A_d,
            "B": B_d,
            "invariant": invariant,
            "anti": anti,
            "reconstruction": reconstructed,
        }

    return recovered


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 82 — EXACT BIVARIATE INDEX RECONSTRUCTION")
    print("=" * 78)

    print()
    print("=" * 78)
    print("0. EXACT SYMBOLIC SETUP")
    print("=" * 78)
    print("  t = N/X")
    print("  u = t(t+1)")
    print("  j = channel index")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    # -----------------------------------------------------------------------
    # 1. Involution
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT INVOLUTION VALIDATION")
    print("=" * 78)

    recovered = recover_A_B()

    all_involution = True

    for d in range(16, -1, -1):
        info = recovered[d]
        print(
            f"  d={d:2d}: "
            f"invariant={info['invariant']} "
            f"anti_invariant={info['anti']} "
            f"reconstruction={info['reconstruction']}"
        )

        all_involution = (
            all_involution
            and info["reconstruction"]
            and (
                info["invariant"] if d % 2 == 0
                else info["anti"]
            )
        )

    print()
    print(f"  ALL INVOLUTION CHECKS = {all_involution}")

    # -----------------------------------------------------------------------
    # 2. Build actual A_j/B_j sequences
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SCALAR CHANNEL DATA")
    print("=" * 78)

    A_seq = [clean(A[k]) for k in range(9)]
    B_seq = [clean(B[k]) for k in range(8)]

    print("  A_j = A_{2j}, j=0,...,8")
    for idx, expr in enumerate(A_seq):
        print(
            f"    j={idx}: "
            f"degree_u={degree_u(expr)} "
            f"coefficients={poly_u(expr).all_coeffs()}"
        )

    print()
    print("  B_j = B_{2j+1}, j=0,...,7")
    for idx, expr in enumerate(B_seq):
        print(
            f"    j={idx}: "
            f"degree_u={degree_u(expr)} "
            f"coefficients={poly_u(expr).all_coeffs()}"
        )

    # -----------------------------------------------------------------------
    # 3. Coefficient extraction and exact j interpolation
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT COEFFICIENT POLYNOMIALS IN j")
    print("=" * 78)

    A_j_polys = {}
    B_j_polys = {}

    max_A_degree = max(degree_u(x) for x in A_seq)
    max_B_degree = max(degree_u(x) for x in B_seq)

    print()
    print("  A-channel coefficient laws")

    for k in range(max_A_degree + 1):
        values = [
            coefficient_u(expr, k)
            for expr in A_seq
        ]

        Pj = clean(interpolate_in_j(values))

        ok = verify_index_polynomial(Pj, values)
        fd = finite_difference_degree(values)

        A_j_polys[k] = Pj

        print()
        print(f"    coefficient [u^{k}]")
        print(f"      finite_difference_degree = {fd}")
        print(f"      interpolated_j_degree    = {sp.degree(Pj, j)}")
        print(f"      verified_all_points      = {ok}")
        print(f"      P_A_{k}(j) = {factor_j(Pj)}")

    print()
    print("  B-channel coefficient laws")

    for k in range(max_B_degree + 1):
        values = [
            coefficient_u(expr, k)
            for expr in B_seq
        ]

        Pj = clean(interpolate_in_j(values))

        ok = verify_index_polynomial(Pj, values)
        fd = finite_difference_degree(values)

        B_j_polys[k] = Pj

        print()
        print(f"    coefficient [u^{k}]")
        print(f"      finite_difference_degree = {fd}")
        print(f"      interpolated_j_degree    = {sp.degree(Pj, j)}")
        print(f"      verified_all_points      = {ok}")
        print(f"      P_B_{k}(j) = {factor_j(Pj)}")

    # -----------------------------------------------------------------------
    # 4. Build the bivariate reconstructed channel families
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. RECONSTRUCTED BIVARIATE CHANNEL LAWS")
    print("=" * 78)

    A_bivariate = clean(
        sum(Pj * u**k for k, Pj in A_j_polys.items())
    )

    B_bivariate = clean(
        sum(Pj * u**k for k, Pj in B_j_polys.items())
    )

    print()
    print("  A(j,u) =")
    print(f"    {factor_ju(A_bivariate)}")

    print()
    print("  B(j,u) =")
    print(f"    {factor_ju(B_bivariate)}")

    # -----------------------------------------------------------------------
    # 5. Exact reconstruction check
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT BIVARIATE RECONSTRUCTION")
    print("=" * 78)

    A_all_ok = True
    B_all_ok = True

    for idx, expr in enumerate(A_seq):
        reconstructed = clean(
            A_bivariate.subs(j, idx) - expr
        ) == 0

        print(
            f"  A_{idx}: reconstructed={reconstructed}"
        )

        A_all_ok = A_all_ok and reconstructed

    print()

    for idx, expr in enumerate(B_seq):
        reconstructed = clean(
            B_bivariate.subs(j, idx) - expr
        ) == 0

        print(
            f"  B_{idx}: reconstructed={reconstructed}"
        )

        B_all_ok = B_all_ok and reconstructed

    print()
    print(f"  A reconstruction = {A_all_ok}")
    print(f"  B reconstruction = {B_all_ok}")

    # -----------------------------------------------------------------------
    # 6. Degree table
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. BIVARIATE DEGREE PROFILE")
    print("=" * 78)

    print()
    print("  A(j,u):")
    print(f"    degree_j = {sp.degree(A_bivariate, j)}")
    print(f"    degree_u = {sp.degree(A_bivariate, u)}")
    print(
        f"    total_degree = "
        f"{sp.Poly(A_bivariate, j, u, domain=sp.QQ).total_degree()}"
    )

    print()
    print("  B(j,u):")
    print(f"    degree_j = {sp.degree(B_bivariate, j)}")
    print(f"    degree_u = {sp.degree(B_bivariate, u)}")
    print(
        f"    total_degree = "
        f"{sp.Poly(B_bivariate, j, u, domain=sp.QQ).total_degree()}"
    )

    # -----------------------------------------------------------------------
    # 7. Factorization audit
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. BIVARIATE FACTORIZATION AUDIT")
    print("=" * 78)

    A_factor = sp.factor(A_bivariate)
    B_factor = sp.factor(B_bivariate)

    print()
    print("  factor(A(j,u)) =")
    print(f"    {A_factor}")

    print()
    print("  factor(B(j,u)) =")
    print(f"    {B_factor}")

    print()
    print("  factor over QQ[j,u] complete =")
    print(
        "    A:",
        sp.factor_list(A_bivariate),
    )
    print(
        "    B:",
        sp.factor_list(B_bivariate),
    )

    # -----------------------------------------------------------------------
    # 8. Common-factor audit
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. COMMON FACTOR AUDIT")
    print("=" * 78)

    gcd_AB = sp.gcd(
        sp.Poly(A_bivariate, j, u, domain=sp.QQ),
        sp.Poly(B_bivariate, j, u, domain=sp.QQ),
    ).as_expr()

    print(f"  gcd(A(j,u), B(j,u)) = {sp.factor(gcd_AB)}")

    # -----------------------------------------------------------------------
    # 9. Check whether coefficient polynomial degrees have a simple law
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. INDEX-DEGREE PROFILE OF INDIVIDUAL COEFFICIENTS")
    print("=" * 78)

    print()
    print("  A coefficients:")

    for k, Pj in A_j_polys.items():
        print(
            f"    [u^{k}]: degree_j={sp.degree(Pj, j)} "
            f"leading={sp.LC(sp.Poly(Pj, j, domain=sp.QQ))}"
        )

    print()
    print("  B coefficients:")

    for k, Pj in B_j_polys.items():
        print(
            f"    [u^{k}]: degree_j={sp.degree(Pj, j)} "
            f"leading={sp.LC(sp.Poly(Pj, j, domain=sp.QQ))}"
        )

    # -----------------------------------------------------------------------
    # 10. Newton forward-difference representation
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. NEWTON FORWARD-DIFFERENCE REPRESENTATION")
    print("=" * 78)

    def newton_forward(values):
        vals = [sp.Rational(v) for v in values]

        coeffs = []
        current = vals[:]

        while current:
            coeffs.append(current[0])
            current = [
                sp.simplify(current[i + 1] - current[i])
                for i in range(len(current) - 1)
            ]

        result = sp.Integer(0)

        for r, c in enumerate(coeffs):
            if c == 0:
                continue

            falling = sp.Integer(1)

            for m in range(r):
                falling *= (j - m)

            result += c * falling / sp.factorial(r)

        return sp.expand(result)

    print()
    print("  A-channel Newton forms")

    for k in range(max_A_degree + 1):
        values = [coefficient_u(expr, k) for expr in A_seq]
        P_newton = clean(newton_forward(values))
        P_direct = clean(A_j_polys[k])

        print(
            f"    [u^{k}] matches direct interpolation = "
            f"{P_newton == P_direct}"
        )

    print()
    print("  B-channel Newton forms")

    for k in range(max_B_degree + 1):
        values = [coefficient_u(expr, k) for expr in B_seq]
        P_newton = clean(newton_forward(values))
        P_direct = clean(B_j_polys[k])

        print(
            f"    [u^{k}] matches direct interpolation = "
            f"{P_newton == P_direct}"
        )

    # -----------------------------------------------------------------------
    # 11. Fresh exact parameter substitution
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. FRESH EXACT SEMIPRIME CONSISTENCY TEST")
    print("=" * 78)

    p = sp.Integer(106621)
    q = sp.Integer(246473)

    N_val = p*q
    S_val = p+q
    X_val = S_val + 1

    t_val = sp.Rational(N_val, X_val)
    u_val = clean(t_val * (t_val + 1))

    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {N_val}")
    print(f"  S = {S_val}")
    print(f"  X = {X_val}")
    print(f"  t = {t_val}")
    print(f"  u = {u_val}")

    fresh_A_ok = True
    fresh_B_ok = True

    # Test every known j.
    for idx, expr in enumerate(A_seq):
        actual = clean(expr.subs(u, u_val))
        predicted = clean(
            A_bivariate.subs({
                j: idx,
                u: u_val,
            })
        )

        ok = actual == predicted
        fresh_A_ok = fresh_A_ok and ok

        print(
            f"  fresh A_{idx}: {ok}"
        )

    for idx, expr in enumerate(B_seq):
        actual = clean(expr.subs(u, u_val))
        predicted = clean(
            B_bivariate.subs({
                j: idx,
                u: u_val,
            })
        )

        ok = actual == predicted
        fresh_B_ok = fresh_B_ok and ok

        print(
            f"  fresh B_{idx}: {ok}"
        )

    # -----------------------------------------------------------------------
    # 12. Reconstruct full h_d(t) from indexed laws
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. FULL h_d(t) RECONSTRUCTION FROM INDEXED LAWS")
    print("=" * 78)

    full_h_ok = True

    for d in range(17):

        if d % 2 == 0:
            idx = d // 2

            generated_A = clean(
                A_bivariate.subs({
                    j: idx,
                    u: clean(t*(t+1)),
                })
            )

            generated = clean(generated_A)

        else:
            idx = (d - 1) // 2

            generated_B = clean(
                B_bivariate.subs({
                    j: idx,
                    u: clean(t*(t+1)),
                })
            )

            generated = clean(
                (2*t + 1) * generated_B
            )

        ok = clean(generated - H[d]) == 0
        full_h_ok = full_h_ok and ok

        print(
            f"  d={d:2d}: exact={ok}"
        )

    # -----------------------------------------------------------------------
    # 13. Final interpretation
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The expensive recurrence search has been removed.

  Instead we construct the exact finite bivariate representation

      A_j(u) = A(j,u)
      B_j(u) = B(j,u)

  by reconstructing every coefficient [u^k] as a polynomial in j.

  This answers a different and much more controlled question:

      Does the observed layer family admit a compact
      polynomial law in the layer index?

  Important:

      This is an exact identity on the available
      index range.

  It is NOT automatically a proof that the same
  polynomial law holds for arbitrary future j.

  The important structural possibilities are now:

    1. A(j,u) and B(j,u) factor substantially.
       This would expose the hidden construction blocks.

    2. Their coefficient polynomials factor in j.
       This could reveal binomial/falling-factorial structure.

    3. The Newton forward-difference coefficients
       display simple arithmetic patterns.
       That would suggest a combinatorial construction.

    4. No simplification occurs.
       Then the complicated coefficient structure is
       probably intrinsic to this kernel.

  This experiment therefore targets the construction law
  directly while avoiding the expensive unrestricted
  recurrence search.

  The inverse implication remains:

      observed layers
          -> t
          -> X
          -> N
          -> p,q

  The unresolved implication remains:

      N
          -> layer information
          -> t
          -> X
          -> p,q.
        """
    )

    # -----------------------------------------------------------------------
    # 14. Exactness
    # -----------------------------------------------------------------------

    checks = {
        "involution": all_involution,
        "A_index_reconstruction": A_all_ok,
        "B_index_reconstruction": B_all_ok,
        "fresh_A": fresh_A_ok,
        "fresh_B": fresh_B_ok,
        "full_h": full_h_ok,
    }

    failures = sum(not bool(v) for v in checks.values())

    print()
    print("=" * 78)
    print("FINAL EXACTNESS")
    print("=" * 78)

    for name, value in checks.items():
        print(f"  {name} = {value}")

    print(f"  failures = {failures}")
    print(f"  ALL EXACT CHECKS PASS = {failures == 0}")

    print()
    print("=" * 78)
    print("EXPERIMENT 82 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()
