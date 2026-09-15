#!/usr/bin/env python3

import sympy as sp

# ==============================================================
# EXPERIMENT 84
# EXACT FALLING-FACTORIAL QUOTIENT / COMBINATORIAL KERNEL AUDIT
# ==============================================================

u = sp.Symbol("u")
j = sp.Symbol("j")


# --------------------------------------------------------------
# Exact scalar channel data
# --------------------------------------------------------------

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
         + 92208*u**3 + 63401*u**2 + 18109*u + 1820),
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
    4: 17875*u**4 + 138710*u**3 + 188409*u**2
       + 83097*u + 11441,
    5: 13 * (50*u**5 + 1784*u**4 + 6480*u**3
             + 6812*u**2 + 2639*u + 336),
    6: 144*u**5 + 3370*u**4 + 11332*u**3
       + 11589*u**2 + 4431*u + 560,
    7: 44*u**4 + 230*u**3 + 282*u**2 + 119*u + 16,
}


# --------------------------------------------------------------
# Safe exact helpers
# --------------------------------------------------------------

def S(x):
    return sp.sympify(x)


def clean(x):
    return sp.cancel(sp.expand(S(x)))


def poly(expr, var):
    return sp.Poly(clean(expr), var, domain=sp.QQ)


def degree(expr, var):
    p = poly(expr, var)
    return -sp.oo if p.is_zero else p.degree()


def coeff_u(expr, k):
    return poly(expr, u).nth(k)


def falling(k):
    result = sp.Integer(1)

    for r in range(k):
        result *= (j - r)

    return sp.expand(result)


def binom_poly(k):
    return sp.expand(falling(k) / sp.factorial(k))


def divide_j(expr, divisor):
    p = poly(expr, j)
    q = poly(divisor, j)
    quotient, remainder = sp.div(p, q, domain=sp.QQ)

    return (
        clean(quotient.as_expr()),
        clean(remainder.as_expr()),
    )


def interpolate(values):
    points = [
        (sp.Integer(idx), S(value))
        for idx, value in enumerate(values)
    ]

    return clean(sp.interpolate(points, j))


def evaluate_integer(expr, n):
    return clean(S(expr).subs(j, sp.Integer(n)))


def finite_differences(values):
    current = [S(v) for v in values]
    rows = [current]

    while len(current) > 1:
        current = [
            clean(current[r + 1] - current[r])
            for r in range(len(current) - 1)
        ]
        rows.append(current)

    return rows


def gcd_list(expressions):
    expressions = [clean(x) for x in expressions if clean(x) != 0]

    if not expressions:
        return sp.Integer(0)

    g = expressions[0]

    for expr in expressions[1:]:
        g = sp.gcd(g, expr)

    return clean(g)


def content(expr):
    p = poly(expr, j)

    if p.is_zero:
        return sp.Integer(0)

    coeffs = [
        sp.Rational(c)
        for c in p.all_coeffs()
    ]

    denominator_lcm = sp.Integer(1)

    for c in coeffs:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(c.q),
        )

    integer_coeffs = [
        int(c * denominator_lcm)
        for c in coeffs
    ]

    g = 0

    for value in integer_coeffs:
        g = sp.igcd(g, abs(value))

    return sp.Rational(g, denominator_lcm)


def primitive_part(expr):
    expr = clean(expr)
    c = content(expr)

    if c == 0:
        return sp.Integer(0)

    return clean(expr / c)


# --------------------------------------------------------------
# Main
# --------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 84 — EXACT FALLING-FACTORIAL QUOTIENT / COMBINATORIAL KERNEL AUDIT")
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
    # 1. Construct coefficient polynomials in j
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT COEFFICIENT LAWS")
    print("=" * 78)

    A_coeff = {}
    B_coeff = {}

    max_A_k = max(
        int(degree(expr, u))
        for expr in A.values()
    )

    max_B_k = max(
        int(degree(expr, u))
        for expr in B.values()
    )

    print()
    print("  A-channel")

    for k in range(max_A_k + 1):

        values = [
            coeff_u(A[idx], k)
            for idx in range(9)
        ]

        P = interpolate(values)
        A_coeff[k] = P

        verified = all(
            evaluate_integer(P, idx) == values[idx]
            for idx in range(9)
        )

        print(
            f"    [u^{k}] degree_j={degree(P, j)} "
            f"verified={verified}"
        )

    print()
    print("  B-channel")

    for k in range(max_B_k + 1):

        values = [
            coeff_u(B[idx], k)
            for idx in range(8)
        ]

        P = interpolate(values)
        B_coeff[k] = P

        verified = all(
            evaluate_integer(P, idx) == values[idx]
            for idx in range(8)
        )

        print(
            f"    [u^{k}] degree_j={degree(P, j)} "
            f"verified={verified}"
        )

    # ----------------------------------------------------------
    # 2. Remove natural falling factorial
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FALLING-FACTORIAL QUOTIENTS")
    print("=" * 78)

    A_quot = {}
    B_quot = {}

    print()
    print("  A-channel")

    for k, P in A_coeff.items():

        divisor = falling(k)
        Q, R = divide_j(P, divisor)

        exact = R == 0

        if exact:
            A_quot[k] = Q

        print()
        print(f"    [u^{k}]")
        print(f"      P(j)      = {sp.factor(P)}")
        print(f"      divisor   = {divisor}")
        print(f"      exact     = {exact}")

        if exact:
            print(f"      Q(j)      = {sp.factor(Q)}")
            print(f"      degree(Q) = {degree(Q, j)}")
            print(f"      content   = {content(Q)}")
        else:
            print(f"      remainder = {sp.factor(R)}")

    print()
    print("  B-channel")

    for k, P in B_coeff.items():

        divisor = falling(k)
        Q, R = divide_j(P, divisor)

        exact = R == 0

        if exact:
            B_quot[k] = Q

        print()
        print(f"    [u^{k}]")
        print(f"      P(j)      = {sp.factor(P)}")
        print(f"      divisor   = {divisor}")
        print(f"      exact     = {exact}")

        if exact:
            print(f"      Q(j)      = {sp.factor(Q)}")
            print(f"      degree(Q) = {degree(Q, j)}")
            print(f"      content   = {content(Q)}")
        else:
            print(f"      remainder = {sp.factor(R)}")

    # ----------------------------------------------------------
    # 3. Search for additional endpoint factors
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("3. ADDITIONAL ENDPOINT-FACTOR AUDIT")
    print("=" * 78)

    print()
    print("  Candidate factors:")
    print("    j-k")
    print("    j-(m-k)")
    print("    j")
    print("    j-1")
    print("    j-2")
    print("    j-3")

    def endpoint_audit(quotients, max_index):

        for k, Q in quotients.items():

            candidates = []

            for r in range(0, max_index + 1):

                factor_expr = j - r
                QQ, RR = divide_j(Q, factor_expr)

                if RR == 0:
                    candidates.append(r)

            print(
                f"    [u^{k}] additional roots/factors = {candidates}"
            )

    print()
    print("  A-channel")
    endpoint_audit(A_quot, 8)

    print()
    print("  B-channel")
    endpoint_audit(B_quot, 7)

    # ----------------------------------------------------------
    # 4. Quotient finite differences
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FINITE-DIFFERENCE DEGREE AFTER FALLING-FACTORIAL REMOVAL")
    print("=" * 78)

    A_fd_degree = {}
    B_fd_degree = {}

    print()
    print("  A-channel")

    for k, Q in A_quot.items():

        values = [
            evaluate_integer(Q, idx)
            for idx in range(9)
        ]

        rows = finite_differences(values)

        deg = None

        for r, row in enumerate(rows):

            if all(v == 0 for v in row):
                deg = r - 1
                break

        if deg is None:
            deg = len(rows) - 1

        A_fd_degree[k] = deg

        print(
            f"    [u^{k}] quotient finite-difference degree = {deg}"
        )

    print()
    print("  B-channel")

    for k, Q in B_quot.items():

        values = [
            evaluate_integer(Q, idx)
            for idx in range(8)
        ]

        rows = finite_differences(values)

        deg = None

        for r, row in enumerate(rows):

            if all(v == 0 for v in row):
                deg = r - 1
                break

        if deg is None:
            deg = len(rows) - 1

        B_fd_degree[k] = deg

        print(
            f"    [u^{k}] quotient finite-difference degree = {deg}"
        )

    # ----------------------------------------------------------
    # 5. Binomial-basis conversion of quotient
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("5. BINOMIAL-BASIS STRUCTURE OF QUOTIENTS")
    print("=" * 78)

    print()
    print("  A-channel")

    A_binom_quot = {}

    for k, Q in A_quot.items():

        vals = [
            evaluate_integer(Q, idx)
            for idx in range(9)
        ]

        diffs = finite_differences(vals)

        coeffs = []

        for r in range(len(diffs)):
            if diffs[r]:
                coeffs.append((r, diffs[r][0]))

        A_binom_quot[k] = coeffs

        print(
            f"    [u^{k}] Newton coefficients = {coeffs}"
        )

    print()
    print("  B-channel")

    B_binom_quot = {}

    for k, Q in B_quot.items():

        vals = [
            evaluate_integer(Q, idx)
            for idx in range(8)
        ]

        diffs = finite_differences(vals)

        coeffs = []

        for r in range(len(diffs)):
            if diffs[r]:
                coeffs.append((r, diffs[r][0]))

        B_binom_quot[k] = coeffs

        print(
            f"    [u^{k}] Newton coefficients = {coeffs}"
        )

    # ----------------------------------------------------------
    # 6. Cross-channel quotient gcd
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CROSS-CHANNEL QUOTIENT GCD AUDIT")
    print("=" * 78)

    quotient_gcds = []

    common = sorted(
        set(A_quot.keys()) &
        set(B_quot.keys())
    )

    for k in common:

        Qa = primitive_part(A_quot[k])
        Qb = primitive_part(B_quot[k])

        g = clean(sp.gcd(Qa, Qb))

        quotient_gcds.append(g)

        print()
        print(f"    k={k}")
        print(f"      primitive A quotient = {sp.factor(Qa)}")
        print(f"      primitive B quotient = {sp.factor(Qb)}")
        print(f"      gcd = {sp.factor(g)}")

    global_gcd = gcd_list(quotient_gcds)

    print()
    print(f"  global cross-channel gcd = {sp.factor(global_gcd)}")

    # ----------------------------------------------------------
    # 7. Search for simple common kernels
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("7. SIMPLE COMMON-KERNEL SEARCH")
    print("=" * 78)

    candidates = [
        j,
        j - 1,
        j - 2,
        j - 3,
        j - 4,
        8 - j,
        7 - j,
        8 - 2*j,
        7 - 2*j,
        2*j + 1,
        2*j + 2,
    ]

    for candidate in candidates:

        a_hits = []
        b_hits = []

        for k, Q in A_quot.items():

            _, R = divide_j(Q, candidate)

            if R == 0:
                a_hits.append(k)

        for k, Q in B_quot.items():

            _, R = divide_j(Q, candidate)

            if R == 0:
                b_hits.append(k)

        if a_hits or b_hits:

            print()
            print(f"  candidate = {candidate}")
            print(f"    A hits = {a_hits}")
            print(f"    B hits = {b_hits}")

    # ----------------------------------------------------------
    # 8. Fresh exact consistency
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("8. FRESH EXACT SEMIPRIME CONSISTENCY")
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

        expr = S(expr)

        generated = sp.Integer(0)

        for k, P in A_coeff.items():
            generated += P.subs(j, idx) * u_val**k

        actual = clean(expr.subs(u, u_val))
        generated = clean(generated)

        ok = actual == generated
        fresh_ok = fresh_ok and ok

        print(f"  A_{idx}: {ok}")

    for idx, expr in B.items():

        expr = S(expr)

        generated = sp.Integer(0)

        for k, P in B_coeff.items():
            generated += P.subs(j, idx) * u_val**k

        actual = clean(expr.subs(u, u_val))
        generated = clean(generated)

        ok = actual == generated
        fresh_ok = fresh_ok and ok

        print(f"  B_{idx}: {ok}")

    # ----------------------------------------------------------
    # 9. Structural summary
    # ----------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The previous experiment established that every coefficient [u^k]
  contains the exact falling-factorial factor

      j(j-1)...(j-k+1).

  This experiment removes that factor and studies what remains.

  The key object is therefore

      P_k(j) = [u^k] / j_under_k.

  Three possibilities are being separated.

    1. P_k(j) is simple.
       Then the coefficient may come from a compact combinatorial
       construction.

    2. P_k(j) has additional endpoint factors.
       Then the finite layer range may encode boundary conditions.

    3. P_k(j) remains complicated.
       Then interpolation may describe the data exactly without
       exposing the underlying construction.

  In particular, if the quotient degrees fall sharply after removal
  of the falling factorial, that is strong evidence that the
  coefficient family is organized combinatorially.

  If A and B quotients share factors, those factors may represent
  common construction blocks.

  This experiment does not extrapolate beyond the observed index
  ranges and does not claim a factorization algorithm.

  The inverse direction remains

      observed layers -> t -> X -> N -> p,q.

  The unresolved direction remains

      N -> layer information -> t -> X -> p,q.
        """
    )

    # ----------------------------------------------------------
    # 10. Exactness audit
    # ----------------------------------------------------------

    checks = {
        "A_coefficients_exact": all(
            degree(P, j) >= 0
            for P in A_coeff.values()
        ),
        "B_coefficients_exact": all(
            degree(P, j) >= 0
            for P in B_coeff.values()
        ),
        "A_falling_factorials": len(A_quot) == len(A_coeff),
        "B_falling_factorials": len(B_quot) == len(B_coeff),
        "fresh_validation": fresh_ok,
    }

    failures = sum(
        1
        for value in checks.values()
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
    print("EXPERIMENT 84 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

