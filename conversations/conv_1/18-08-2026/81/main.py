#!/usr/bin/env python3

import sympy as sp

# ==============================================================
# EXPERIMENT 83
# EXACT FALLING-FACTORIAL / BINOMIAL CONSTRUCTION AUDIT
# ==============================================================

u = sp.Symbol("u")
j = sp.Symbol("j")
t = sp.Symbol("t")


# --------------------------------------------------------------
# Exact channel data
# --------------------------------------------------------------

A = {
    0: -2,
    1: -2 * (275*u + 78),
    2: -2 * (7125*u**2 + 5565*u + 973),
    3: -4 * (17850*u**3 + 32538*u**2 + 15036*u + 2023),
    4: -(71500*u**4 + 326664*u**3 + 323868*u**2 + 112964*u + 12879),
    5: -1001 * (10*u**5 + 152*u**4 + 340*u**3 + 256*u**2 + 77*u + 8),
    6: -(50*u**6 + 6686*u**5 + 50200*u**4 + 92208*u**3
         + 63401*u**2 + 18109*u + 1820),
    7: -(276*u**5 + 2700*u**4 + 5478*u**3 + 3966*u**2 + 1169*u + 120),
    8: -(9*u**4 + 30*u**3 + 27*u**2 + 9*u + 1),
}

B = {
    0: 25,
    1: 14 * (125*u + 46),
    2: 6 * (3230*u**2 + 3458*u + 749),
    3: 4 * (11050*u**3 + 29886*u**2 + 17663*u + 2869),
    4: 17875*u**4 + 138710*u**3 + 188409*u**2 + 83097*u + 11441,
    5: 13 * (50*u**5 + 1784*u**4 + 6480*u**3
             + 6812*u**2 + 2639*u + 336),
    6: 144*u**5 + 3370*u**4 + 11332*u**3
       + 11589*u**2 + 4431*u + 560,
    7: 44*u**4 + 230*u**3 + 282*u**2 + 119*u + 16,
}


# --------------------------------------------------------------
# Helpers
# --------------------------------------------------------------

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def poly_in(expr, variable):
    expr = clean(expr)
    return sp.Poly(expr, variable, domain=sp.QQ)


def degree_in(expr, variable):
    p = poly_in(expr, variable)
    return -sp.oo if p.is_zero else p.degree()


def coeff_u(expr, k):
    return sp.expand(poly_in(expr, u).nth(k))


def falling_factorial(k):
    result = sp.Integer(1)
    for r in range(k):
        result *= (j - r)
    return sp.expand(result)


def binomial_polynomial(k):
    return sp.expand(falling_factorial(k) / sp.factorial(k))


def divide_in_j(expr, divisor):
    p = poly_in(expr, j)
    q = poly_in(divisor, j)
    quotient, remainder = sp.div(p, q, domain=sp.QQ)
    return clean(quotient.as_expr()), clean(remainder.as_expr())


def interpolate_values(values):
    pts = [(sp.Integer(idx), sp.Rational(value))
           for idx, value in enumerate(values)]
    return clean(sp.interpolate(pts, j))


def factor(expr):
    return sp.factor(clean(expr))


def integer_content(expr):
    p = poly_in(expr, j)

    if p.is_zero:
        return sp.Integer(0)

    coeffs = [sp.Rational(c) for c in p.all_coeffs()]

    denominator_lcm = 1
    for c in coeffs:
        denominator_lcm = sp.ilcm(denominator_lcm, int(c.q))

    integers = [int(c * denominator_lcm) for c in coeffs]

    g = 0
    for value in integers:
        g = sp.igcd(g, abs(value))

    if g == 0:
        return sp.Integer(0)

    return sp.Rational(g, denominator_lcm)


def newton_coefficients(values):
    current = [sp.Rational(v) for v in values]
    result = []

    while current:
        result.append(clean(current[0]))

        if len(current) == 1:
            break

        current = [
            clean(current[r + 1] - current[r])
            for r in range(len(current) - 1)
        ]

    return result


# --------------------------------------------------------------
# Main
# --------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 83 — EXACT FALLING-FACTORIAL / BINOMIAL CONSTRUCTION AUDIT")
    print("=" * 78)

    print()
    print("=" * 78)
    print("0. EXACT SYMBOLIC SETUP")
    print("=" * 78)
    print("  A_j = A_{2j}(u), j=0,...,8")
    print("  B_j = B_{2j+1}(u), j=0,...,7")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    # ----------------------------------------------------------
    # 1. Channel validation
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("1. CHANNEL DATA VALIDATION")
    print("=" * 78)

    a_loaded = True
    b_loaded = True

    for idx, expr in A.items():
        ok = clean(expr - A[idx]) == 0
        a_loaded = a_loaded and ok

    for idx, expr in B.items():
        ok = clean(expr - B[idx]) == 0
        b_loaded = b_loaded and ok

    print(f"  A channel loaded exactly = {a_loaded}")
    print(f"  B channel loaded exactly = {b_loaded}")

    # ----------------------------------------------------------
    # 2. Build coefficient polynomials in j
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT COEFFICIENT POLYNOMIALS IN j")
    print("=" * 78)

    max_a_k = max(
        int(degree_in(expr, u))
        for expr in A.values()
    )

    max_b_k = max(
        int(degree_in(expr, u))
        for expr in B.values()
    )

    A_poly = {}
    B_poly = {}

    print()
    print("  A-channel")

    for k in range(max_a_k + 1):

        values = [
            coeff_u(A[idx], k)
            for idx in range(9)
        ]

        P = interpolate_values(values)
        A_poly[k] = P

        print(
            f"    [u^{k}] "
            f"degree_j={degree_in(P, j)} "
            f"verified={all(clean(P.subs(j, idx) - values[idx]) == 0 for idx in range(9))}"
        )

    print()
    print("  B-channel")

    for k in range(max_b_k + 1):

        values = [
            coeff_u(B[idx], k)
            for idx in range(8)
        ]

        P = interpolate_values(values)
        B_poly[k] = P

        print(
            f"    [u^{k}] "
            f"degree_j={degree_in(P, j)} "
            f"verified={all(clean(P.subs(j, idx) - values[idx]) == 0 for idx in range(8))}"
        )

    # ----------------------------------------------------------
    # 3. Falling factorial divisibility
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FALLING-FACTORIAL DIVISIBILITY")
    print("=" * 78)

    A_ff = {}
    B_ff = {}

    print()
    print("  A-channel")

    for k, P in A_poly.items():

        divisor = falling_factorial(k)
        quotient, remainder = divide_in_j(P, divisor)
        exact = remainder == 0

        A_ff[k] = (quotient, remainder, exact)

        print()
        print(f"    [u^{k}]")
        print(f"      divisor = {divisor}")
        print(f"      exact = {exact}")

        if exact:
            print(f"      quotient = {factor(quotient)}")
            print(f"      quotient_degree = {degree_in(quotient, j)}")
        else:
            print(f"      remainder = {factor(remainder)}")

    print()
    print("  B-channel")

    for k, P in B_poly.items():

        divisor = falling_factorial(k)
        quotient, remainder = divide_in_j(P, divisor)
        exact = remainder == 0

        B_ff[k] = (quotient, remainder, exact)

        print()
        print(f"    [u^{k}]")
        print(f"      divisor = {divisor}")
        print(f"      exact = {exact}")

        if exact:
            print(f"      quotient = {factor(quotient)}")
            print(f"      quotient_degree = {degree_in(quotient, j)}")
        else:
            print(f"      remainder = {factor(remainder)}")

    # ----------------------------------------------------------
    # 4. Binomial normalization
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("4. BINOMIAL-BASIS NORMALIZATION")
    print("=" * 78)

    print()
    print("  A-channel")

    A_binom = {}

    for k, P in A_poly.items():

        divisor = binomial_polynomial(k)
        quotient, remainder = divide_in_j(P, divisor)
        exact = remainder == 0

        if exact:
            A_binom[k] = quotient
            print(
                f"    [u^{k}] "
                f"exact=True "
                f"degree={degree_in(quotient, j)}"
            )
            print(f"      quotient = {factor(quotient)}")
        else:
            print(
                f"    [u^{k}] "
                f"exact=False"
            )
            print(f"      remainder = {factor(remainder)}")

    print()
    print("  B-channel")

    B_binom = {}

    for k, P in B_poly.items():

        divisor = binomial_polynomial(k)
        quotient, remainder = divide_in_j(P, divisor)
        exact = remainder == 0

        if exact:
            B_binom[k] = quotient
            print(
                f"    [u^{k}] "
                f"exact=True "
                f"degree={degree_in(quotient, j)}"
            )
            print(f"      quotient = {factor(quotient)}")
        else:
            print(
                f"    [u^{k}] "
                f"exact=False"
            )
            print(f"      remainder = {factor(remainder)}")

    # ----------------------------------------------------------
    # 5. Quotient complexity
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("5. QUOTIENT COMPLEXITY")
    print("=" * 78)

    print()
    print("  A-channel")

    for k, qexpr in A_binom.items():
        p = poly_in(qexpr, j)
        print(
            f"    [u^{k}] "
            f"degree={p.degree()} "
            f"terms={len(p.terms())} "
            f"content={integer_content(qexpr)}"
        )

    print()
    print("  B-channel")

    for k, qexpr in B_binom.items():
        p = poly_in(qexpr, j)
        print(
            f"    [u^{k}] "
            f"degree={p.degree()} "
            f"terms={len(p.terms())} "
            f"content={integer_content(qexpr)}"
        )

    # ----------------------------------------------------------
    # 6. Cross-channel comparison
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CROSS-CHANNEL BINOMIAL QUOTIENT COMPARISON")
    print("=" * 78)

    common = sorted(
        set(A_binom.keys()) &
        set(B_binom.keys())
    )

    if not common:
        print("  no common coefficient indices")
    else:
        for k in common:

            qa = clean(A_binom[k])
            qb = clean(B_binom[k])
            ratio = clean(qa / qb)

            print()
            print(f"    k={k}")
            print(f"      A quotient = {factor(qa)}")
            print(f"      B quotient = {factor(qb)}")
            print(f"      ratio      = {factor(ratio)}")

    # ----------------------------------------------------------
    # 7. Explicit Newton coefficients
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("7. NEWTON FORWARD-DIFFERENCE COEFFICIENTS")
    print("=" * 78)

    print()
    print("  A-channel")

    for k in range(max_a_k + 1):

        values = [
            coeff_u(A[idx], k)
            for idx in range(9)
        ]

        coeffs = newton_coefficients(values)

        print(f"    [u^{k}] = {coeffs}")

    print()
    print("  B-channel")

    for k in range(max_b_k + 1):

        values = [
            coeff_u(B[idx], k)
            for idx in range(8)
        ]

        coeffs = newton_coefficients(values)

        print(f"    [u^{k}] = {coeffs}")

    # ----------------------------------------------------------
    # 8. Root / zero audit
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("8. INTEGER-ZERO AUDIT OF INDEX POLYNOMIALS")
    print("=" * 78)

    print()
    print("  A-channel")

    for k, P in A_poly.items():

        factored = factor(P)
        roots = []

        for r in range(0, 9):
            if clean(P.subs(j, r)) == 0:
                roots.append(r)

        print(
            f"    [u^{k}] roots within 0..8 = {roots}"
        )
        print(
            f"      factorization = {factored}"
        )

    print()
    print("  B-channel")

    for k, P in B_poly.items():

        factored = factor(P)
        roots = []

        for r in range(0, 8):
            if clean(P.subs(j, r)) == 0:
                roots.append(r)

        print(
            f"    [u^{k}] roots within 0..7 = {roots}"
        )
        print(
            f"      factorization = {factored}"
        )

    # ----------------------------------------------------------
    # 9. Strong construction test
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRONG CONSTRUCTION TEST")
    print("=" * 78)

    a_ff_count = sum(
        1 for _, _, exact in A_ff.values()
        if exact
    )

    b_ff_count = sum(
        1 for _, _, exact in B_ff.values()
        if exact
    )

    print(
        "  falling-factorial exact divisions:"
    )
    print(f"    A = {a_ff_count}/{len(A_ff)}")
    print(f"    B = {b_ff_count}/{len(B_ff)}")

    print()
    print(
        "  Interpretation:"
    )
    print(
        "    A positive result means the coefficient naturally contains"
    )
    print(
        "    the combinatorial factor j(j-1)...(j-k+1)."
    )
    print(
        "    Any endpoint factors beyond the natural falling factor"
    )
    print(
        "    must be treated cautiously because the kernel has only"
    )
    print(
        "    finitely many observed channel indices."
    )

    # ----------------------------------------------------------
    # 10. Fresh exact semiprime
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FRESH EXACT SEMIPRIME CONSISTENCY")
    print("=" * 78)

    p = sp.Integer(106621)
    q = sp.Integer(246473)

    N_val = p * q
    S_val = p + q
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

    fresh_ok = True

    for idx, expr in A.items():

        generated = sp.Integer(0)

        for k, P in A_poly.items():
            generated += P.subs(j, idx) * u_val**k

        actual = clean(expr.subs(u, u_val))
        generated = clean(generated)

        ok = actual == generated
        fresh_ok = fresh_ok and ok

        print(f"  A_{idx}: {ok}")

    for idx, expr in B.items():

        generated = sp.Integer(0)

        for k, P in B_poly.items():
            generated += P.subs(j, idx) * u_val**k

        actual = clean(expr.subs(u, u_val))
        generated = clean(generated)

        ok = actual == generated
        fresh_ok = fresh_ok and ok

        print(f"  B_{idx}: {ok}")

    # ----------------------------------------------------------
    # 11. Final interpretation
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The previous experiment established exact polynomial dependence
  of every channel coefficient on the channel index j.

  This experiment asks whether that dependence is naturally expressed
  in the falling-factorial or binomial basis.

  The important quantities are

      j(j-1)...(j-k+1)

  and

      binomial(j,k).

  A strong positive result would indicate a combinatorial construction
  mechanism behind the layer coefficients.

  A weaker result, where the quotients remain high-degree, would mean
  that polynomial interpolation is an exact description of the finite
  data but does not reveal a simple construction law.

  We do not extrapolate beyond the available j-range.

  We do not claim a factorization algorithm.

  The current inverse direction remains

      observed layers -> t -> X -> N -> p,q.

  The unresolved direction remains

      N -> layer information -> t -> X -> p,q.
        """
    )

    # ----------------------------------------------------------
    # 12. Exactness
    # ----------------------------------------------------------

    checks = {
        "channel_data": a_loaded and b_loaded,
        "falling_factorial_A": a_ff_count > 0,
        "falling_factorial_B": b_ff_count > 0,
        "fresh_validation": fresh_ok,
    }

    failures = sum(
        1 for value in checks.values()
        if not value
    )

    print()
    print("=" * 78)
    print("FINAL EXACTNESS")
    print("=" * 78)

    for name, value in checks.items():
        print(f"  {name} = {value}")

    print(f"  failures = {failures}")
    print(f"  ALL BASIC CHECKS PASS = {failures == 0}")

    print()
    print("=" * 78)
    print("EXPERIMENT 83 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()