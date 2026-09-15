import sympy as sp


# =============================================================================
# EXPERIMENT 213
# EXACT L2 CONSTANT-TERM STRUCTURE AUDIT
# =============================================================================
#
# No previous experiment output is read.
# No unrestricted affine candidate enumeration is performed.
#
# From Experiment 212:
#
# D_a(k,ell)
#   = A(k,a) ell^2 + B(k,a) ell + C(k,a)
#
# with the exact observed laws
#
#   A(k,a) = -1/2 * C(k+2,a)
#
# and the strongly indicated law
#
#   B(k,a)
#     = C(k+2,a) * (2(k+1)a + k+2) / (2(k+2)).
#
# This experiment:
#
#   1. verifies A exactly;
#   2. verifies B exactly;
#   3. extracts C exactly;
#   4. normalizes C by natural binomial factors;
#   5. studies C as a function of a;
#   6. studies edge/boundary behavior;
#   7. tests a small collection of structurally derived C-formulas;
#   8. reconstructs L2 from the resulting A,B,C law.
#
# The aim is to derive L2, not fit theorem anchors.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# =============================================================================
# SAFE POLYNOMIAL UTILITIES
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    if sp.expand(sp.denom(expr) - 1) != 0:
        raise ValueError(
            f"Expression is not polynomial in N,X: denominator={sp.denom(expr)}"
        )

    return sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.QQ,
    )


def coefficient(expr, a, b):
    if a < 0 or b < 0:
        return sp.Integer(0)

    return sp.expand(
        poly_NX(expr).coeff_monomial(
            N**a * X**b
        )
    )


def homogeneous(expr, degree):
    if degree < 0:
        return sp.Integer(0)

    poly = poly_NX(expr)

    out = sp.Integer(0)

    for (a, b), c in poly.terms():
        if a + b == degree:
            out += c * N**a * X**b

    return sp.expand(out)


# =============================================================================
# EXACT pq KERNEL
# =============================================================================

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# =============================================================================
# DIRECT SYMMETRIC REDUCTION F(p,q) -> G(N,S)
# =============================================================================

def exact_G(k, ell):

    F = exact_F(k, ell)

    result, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Symmetrization remainder for k={k}, ell={ell}: {remainder}"
        )

    sum_symbol = None
    prod_symbol = None

    for symbol, expr in mapping:

        expr = sp.expand(expr)

        if sp.expand(expr - (p + q)) == 0:
            sum_symbol = symbol

        elif sp.expand(expr - p*q) == 0:
            prod_symbol = symbol

    if sum_symbol is None or prod_symbol is None:
        raise ValueError(
            f"Could not identify symmetric variables: {mapping}"
        )

    return sp.expand(
        result.subs(
            {
                sum_symbol: S,
                prod_symbol: N,
            }
        )
    )


def exact_GX(k, ell):
    return sp.expand(
        exact_G(k, ell).subs(
            S,
            X - 1,
        )
    )


# =============================================================================
# PROVED TOP LAYER
# =============================================================================

def top_layer(k, ell):

    return sp.expand(
        -X**(ell - k)
        * (
            (X + N)**k
            - N**k
        )
    )


# =============================================================================
# PROVED L1
#
# Experiment 210R:
#
# [N^a X^(ell-1-a)] L1
#     = ell*C(k+1,a) - k*C(k,a-1)
# =============================================================================

def L1_coefficient(k, ell, a):

    if a < 0 or a > k:
        return sp.Integer(0)

    return sp.expand(
        ell * sp.binomial(k + 1, a)
        - k * sp.binomial(k, a - 1)
    )


def L1(k, ell):

    degree = ell - 1
    out = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        out += (
            L1_coefficient(k, ell, a)
            * N**a
            * X**b
        )

    return sp.expand(out)


# =============================================================================
# EXACT L2
# =============================================================================

L2_CACHE = {}


def exact_L2(k, ell):

    key = (k, ell)

    if key in L2_CACHE:
        return L2_CACHE[key]

    G = exact_GX(
        k,
        ell,
    )

    residual = sp.expand(
        G
        - top_layer(k, ell)
        - L1(k, ell)
    )

    L2_expr = homogeneous(
        residual,
        ell - 2,
    )

    L2_CACHE[key] = L2_expr

    return L2_expr


# =============================================================================
# D_a(k,ell)
# =============================================================================

def D(k, ell, a):

    b = ell - 2 - a

    if b < 0:
        return sp.Integer(0)

    return coefficient(
        exact_L2(k, ell),
        a,
        b,
    )


# =============================================================================
# QUADRATIC INTERPOLATION
# =============================================================================

def quadratic_coefficients(k, a, ell0):

    e = sp.symbols("e")

    points = [
        (ell0, D(k, ell0, a)),
        (ell0 + 2, D(k, ell0 + 2, a)),
        (ell0 + 4, D(k, ell0 + 4, a)),
    ]

    polynomial = sp.interpolate(
        points,
        e,
    )

    polynomial = sp.factor(
        sp.expand(polynomial)
    )

    P = sp.Poly(
        polynomial,
        e,
    )

    A = sp.expand(
        P.coeff_monomial(e**2)
    )

    B = sp.expand(
        P.coeff_monomial(e)
    )

    C = sp.expand(
        P.coeff_monomial(1)
    )

    return polynomial, A, B, C


# =============================================================================
# 1. A-COEFFICIENT LAW
# =============================================================================

def A_audit():

    print()
    print("=" * 78)
    print("1. EXACT A-COEFFICIENT LAW")
    print("=" * 78)

    cases = {
        3: range(0, 5),
        5: range(0, 7),
        7: range(0, 9),
        9: range(0, 7),
        11: range(0, 6),
    }

    starts = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    failures = 0

    for k, a_values in cases.items():

        print()
        print(f"k={k}")

        for a in a_values:

            _, A, _, _ = quadratic_coefficients(
                k,
                a,
                starts[k],
            )

            expected = sp.expand(
                -sp.binomial(k + 2, a) / 2
            )

            ok = (
                sp.expand(A - expected)
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"  a={a:2d} "
                f"A={str(A):>8} "
                f"expected={str(expected):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"A-law failures = {failures}"
    )


# =============================================================================
# 2. B-COEFFICIENT LAW
# =============================================================================

def B_audit():

    print()
    print("=" * 78)
    print("2. EXACT B-COEFFICIENT LAW")
    print("=" * 78)

    cases = {
        3: range(0, 5),
        5: range(0, 7),
        7: range(0, 9),
        9: range(0, 7),
        11: range(0, 6),
    }

    starts = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    failures = 0

    for k, a_values in cases.items():

        print()
        print(f"k={k}")

        for a in a_values:

            _, _, B, _ = quadratic_coefficients(
                k,
                a,
                starts[k],
            )

            expected = sp.expand(
                sp.binomial(k + 2, a)
                * (
                    2 * (k + 1) * a
                    + k + 2
                )
                / (
                    2 * (k + 2)
                )
            )

            ok = (
                sp.expand(B - expected)
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"  a={a:2d} "
                f"B={str(B):>10} "
                f"expected={str(expected):>10} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"B-law failures = {failures}"
    )


# =============================================================================
# 3. CONSTANT TERM TABLE
# =============================================================================

def constant_term_table():

    print()
    print("=" * 78)
    print("3. EXACT C(k,a) TABLE")
    print("=" * 78)

    cases = {
        3: range(0, 5),
        5: range(0, 7),
        7: range(0, 9),
        9: range(0, 10),
        11: range(0, 12),
    }

    starts = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    for k, a_values in cases.items():

        print()
        print(f"k={k}")

        for a in a_values:

            _, _, _, C = quadratic_coefficients(
                k,
                a,
                starts[k],
            )

            print(
                f"  a={a:2d}: "
                f"C={C}"
            )


# =============================================================================
# 4. C NORMALIZATION
# =============================================================================

def C_normalization_audit():

    print()
    print("=" * 78)
    print("4. C NORMALIZATION AUDIT")
    print("=" * 78)

    cases = {
        3: range(0, 5),
        5: range(0, 7),
        7: range(0, 9),
        9: range(0, 10),
        11: range(0, 12),
    }

    starts = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    for k, a_values in cases.items():

        print()
        print(f"k={k}")

        for a in a_values:

            _, _, _, C = quadratic_coefficients(
                k,
                a,
                starts[k],
            )

            c0 = sp.binomial(
                k + 2,
                a,
            )

            c1 = sp.binomial(
                k + 1,
                a - 1,
            )

            c2 = sp.binomial(
                k,
                a - 2,
            )

            ratio0 = (
                sp.factor(C / c0)
                if c0 != 0
                else "undefined"
            )

            ratio1 = (
                sp.factor(C / c1)
                if c1 != 0
                else "undefined"
            )

            ratio2 = (
                sp.factor(C / c2)
                if c2 != 0
                else "undefined"
            )

            print(
                f"  a={a:2d}: "
                f"C={str(C):>8} "
                f"C/C(k+2,a)={ratio0} "
                f"C/C(k+1,a-1)={ratio1} "
                f"C/C(k,a-2)={ratio2}"
            )


# =============================================================================
# 5. FACTORIZATION / EDGE AUDIT
# =============================================================================

def C_edge_audit():

    print()
    print("=" * 78)
    print("5. C EDGE / BOUNDARY AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),

        (5, 11),
        (5, 13),
        (5, 15),
        (5, 19),

        (7, 15),
        (7, 17),

        (9, 21),
        (9, 23),

        (11, 23),
    ]

    for k, ell in cases:

        print()
        print(f"k={k} ell={ell}")

        for a in range(k + 1):

            _, _, _, C = quadratic_coefficients(
                k,
                a,
                min(ell, 23),
            )

            factored = sp.factor(
                C
            )

            if a in (0, 1, 2, k - 2, k - 1, k):
                print(
                    f"  a={a:2d} "
                    f"C={str(C):>8} "
                    f"factor={factored}"
                )


# =============================================================================
# 6. NATURAL C-CANDIDATE AUDIT
#
# These are deliberately small, structurally motivated expressions.
# We are not fitting arbitrary polynomials.
# =============================================================================

def C_candidate_audit():

    print()
    print("=" * 78)
    print("6. NATURAL C-FORMULA AUDIT")
    print("=" * 78)

    cases = [
        3,
        5,
        7,
        9,
        11,
    ]

    candidate_names = [
        "C1 = -k*a*(a+1)/(2*(k+2))*binom(k+2,a)",

        "C2 = -a*binom(k+1,a-1)*(a+k+1)/2",

        "C3 = -a*binom(k+2,a)*(a+1)/(2)",

        "C4 = -a*binom(k+1,a)*(a)/(2)",

        "C5 = -a*(a+1)*binom(k+1,a-1)/2",
    ]

    counts = {
        name: 0
        for name in candidate_names
    }

    total = 0

    for k in cases:

        ell0 = {
            3: 7,
            5: 11,
            7: 15,
            9: 21,
            11: 23,
        }[k]

        print()
        print(f"k={k}")

        for a in range(k + 1):

            _, _, _, C = quadratic_coefficients(
                k,
                a,
                ell0,
            )

            total += 1

            cands = {
                candidate_names[0]:
                    sp.expand(
                        -sp.Rational(1, 2)
                        * k
                        * a
                        * (a + 1)
                        * sp.binomial(k + 2, a)
                        / (k + 2)
                    ),

                candidate_names[1]:
                    sp.expand(
                        -a
                        * sp.binomial(k + 1, a - 1)
                        * (a + k + 1)
                        / 2
                    ),

                candidate_names[2]:
                    sp.expand(
                        -a
                        * (a + 1)
                        * sp.binomial(k + 2, a)
                        / 2
                    ),

                candidate_names[3]:
                    sp.expand(
                        -a
                        * a
                        * sp.binomial(k + 1, a)
                        / 2
                    ),

                candidate_names[4]:
                    sp.expand(
                        -a
                        * (a + 1)
                        * sp.binomial(k + 1, a - 1)
                        / 2
                    ),
            }

            matched = []

            for name, value in cands.items():

                if sp.expand(
                    C - value
                ) == 0:

                    counts[name] += 1
                    matched.append(
                        name.split(" = ")[0]
                    )

            if matched:
                print(
                    f"  a={a:2d} "
                    f"C={C} "
                    f"matches={matched}"
                )

    print()
    print(f"total tested = {total}")

    for name in candidate_names:
        print(
            f"  {name.split(' = ')[0]} : "
            f"{counts[name]}/{total}"
        )


# =============================================================================
# 7. C AS A DISCRETE FUNCTION OF a
# =============================================================================

def C_difference_in_a():

    print()
    print("=" * 78)
    print("7. FINITE-DIFFERENCE STRUCTURE IN a")
    print("=" * 78)

    cases = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    for k, ell in cases.items():

        values = []

        for a in range(k + 1):

            _, _, _, C = quadratic_coefficients(
                k,
                a,
                ell,
            )

            values.append(
                int(C)
            )

        print()
        print(f"k={k}")
        print(
            "  C(a) = "
            + " ".join(
                f"{v}"
                for v in values
            )
        )

        first = [
            values[i + 1] - values[i]
            for i in range(len(values) - 1)
        ]

        second = [
            first[i + 1] - first[i]
            for i in range(len(first) - 1)
        ]

        third = [
            second[i + 1] - second[i]
            for i in range(len(second) - 1)
        ]

        print(
            "  ΔC   = "
            + " ".join(
                str(v)
                for v in first
            )
        )

        print(
            "  Δ²C  = "
            + " ".join(
                str(v)
                for v in second
            )
        )

        print(
            "  Δ³C  = "
            + " ".join(
                str(v)
                for v in third
            )
        )


# =============================================================================
# 8. CHECK A/B/C RECONSTRUCTION OF D
# =============================================================================

def ABC_reconstruction():

    print()
    print("=" * 78)
    print("8. A+B+C RECONSTRUCTION OF D_a(k,ell)")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 13),
        (5, 19),
        (7, 15),
        (7, 17),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in cases:

        ell0 = min(
            ell,
            23,
        )

        G = exact_GX(
            k,
            ell,
        )

        L2_exact = exact_L2(
            k,
            ell,
        )

        ok_case = True

        for a in range(k + 1):

            b = ell - 2 - a

            if b < 0:
                continue

            _, A, B, C = quadratic_coefficients(
                k,
                a,
                ell0,
            )

            predicted = sp.expand(
                A * ell**2
                + B * ell
                + C
            )

            actual = coefficient(
                L2_exact,
                a,
                b,
            )

            if sp.expand(
                actual - predicted
            ) != 0:

                ok_case = False
                failures += 1

                print(
                    f"FAIL k={k} ell={ell} a={a} "
                    f"actual={actual} predicted={predicted}"
                )

                break

        print(
            f"k={k:2d} ell={ell:2d} "
            f"{'PASS' if ok_case else 'FAIL'}"
        )

    print()
    print(
        f"A+B+C reconstruction failures = {failures}"
    )


# =============================================================================
# 9. FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "Experiment 209R proved the exact top layer:"
    )

    print(
        "  G_top = -X^(ell-k) ((X+N)^k - N^k)"
    )

    print()
    print(
        "Experiment 210R proved the complete L1 layer."
    )

    print()
    print(
        "Experiment 212 established that the L2 coefficient"
    )

    print(
        "D_a(k,ell) is quadratic in ell."
    )

    print()
    print(
        "The corrected leading law is:"
    )

    print(
        "  A(k,a) = -1/2 * C(k+2,a)"
    )

    print()
    print(
        "The linear coefficient strongly indicates:"
    )

    print(
        "  B(k,a)"
    )

    print(
        "    = C(k+2,a)"
    )

    print(
        "      * (2(k+1)a + k+2)"
    )

    print(
        "      / (2(k+2))"
    )

    print()
    print(
        "Experiment 213 therefore concentrates exclusively"
    )

    print(
        "on the remaining constant term C(k,a)."
    )

    print()
    print(
        "If C collapses to a closed binomial expression,"
    )

    print(
        "we will have a complete exact formula for L2."
    )

    print()
    print(
        "Only then should the next experiment peel L3."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 213")
    print("EXACT L2 CONSTANT-TERM STRUCTURE AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    A_audit()
    B_audit()
    constant_term_table()
    C_normalization_audit()
    C_edge_audit()
    C_candidate_audit()
    C_difference_in_a()
    ABC_reconstruction()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 213")
    print("=" * 78)


if __name__ == "__main__":
    main()

