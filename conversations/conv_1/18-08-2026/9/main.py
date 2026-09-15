import sympy as sp
from functools import lru_cache

# ==============================================================================
# EXPERIMENT 221
# EXACT L3 INTERIOR VALIDATION + BOUNDARY EXTRACTION
# ==============================================================================
#
# No previous experiment output is read.
# No unrestricted candidate enumeration is performed.
# No floating-point arithmetic is used.
#
# Goal:
#
#   exact pq kernel
#       -> G(N,S)
#       -> G(N,X-1)
#       -> top layer
#       -> L1
#       -> L2
#       -> exact L3
#       -> exact interior formula
#       -> boundary coefficients
#
# The previous run exposed two implementation problems:
#
#   (1) Python expressions such as 1/3 created floats.
#   (2) SymPy Integers were formatted using Python numeric format specs.
#
# Both are eliminated here.
# ==============================================================================

p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# ==============================================================================
# EXACT pq KERNEL
# ==============================================================================

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ==============================================================================
# EXACT SYMMETRIC REDUCTION
# ==============================================================================

@lru_cache(maxsize=None)
def exact_G(k, ell):

    F = exact_F(k, ell)

    sym_expr, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Symmetric reduction failed: "
            f"k={k}, ell={ell}, remainder={remainder}"
        )

    sum_symbol = None
    prod_symbol = None

    for symbol, replacement in mapping:

        replacement = sp.expand(replacement)

        if sp.expand(replacement - (p + q)) == 0:
            sum_symbol = symbol

        elif sp.expand(replacement - p*q) == 0:
            prod_symbol = symbol

    if sum_symbol is None or prod_symbol is None:
        raise ValueError(
            f"Could not identify symmetric variables: mapping={mapping}"
        )

    return sp.expand(
        sym_expr.subs(
            {
                sum_symbol: S,
                prod_symbol: N,
            }
        )
    )


@lru_cache(maxsize=None)
def exact_GX(k, ell):
    return sp.expand(
        exact_G(k, ell).subs(
            S,
            X - 1,
        )
    )


# ==============================================================================
# EXACT COEFFICIENT
# ==============================================================================

def coefficient(expr, a, b):

    if a < 0 or b < 0:
        return sp.Integer(0)

    expr = sp.expand(expr)

    poly = sp.Poly(
        expr,
        N,
        X,
        domain="QQ",
    )

    return sp.expand(
        poly.coeff_monomial(
            N**a * X**b
        )
    )


def homogeneous(expr, degree):

    poly = sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain="QQ",
    )

    result = sp.Integer(0)

    for monom, value in poly.terms():

        a, b = monom

        if a + b == degree:
            result += (
                value
                * N**a
                * X**b
            )

    return sp.expand(result)


# ==============================================================================
# TOP LAYER
# ==============================================================================

def top_layer(k, ell):
    return sp.expand(
        -X**(ell - k)
        * (
            (X + N)**k
            - N**k
        )
    )


# ==============================================================================
# EXACT L1
# ==============================================================================

def L1(k, ell):

    degree = ell - 1
    result = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        value = sp.expand(
            ell * sp.binomial(k + 1, a)
            - k * sp.binomial(k, a - 1)
        )

        result += (
            value
            * N**a
            * X**b
        )

    return sp.expand(result)


# ==============================================================================
# EXACT L2
# ==============================================================================

def L2_interior(k, ell, a):

    A = (
        -sp.Rational(1, 2)
        * sp.binomial(k + 2, a)
    )

    B = (
        sp.binomial(k + 2, a)
        * (
            2 * (k + 1) * a
            + k + 2
        )
        / (
            2 * (k + 2)
        )
    )

    C = (
        -k
        * a
        * (a + 1)
        * sp.binomial(k + 2, a)
        / (
            2 * (k + 2)
        )
    )

    return sp.expand(
        A * ell**2
        + B * ell
        + C
    )


def L2_endpoint_k(k, ell):

    return sp.expand(
        L2_interior(k, ell, k)
        + sp.Rational(
            k * (k + 1),
            2,
        )
    )


def L2_endpoint_kplus1(k, ell):

    return sp.expand(
        sp.Rational(
            2,
            k + 1,
        )
        * L2_endpoint_k(
            k,
            ell,
        )
    )


def L2(k, ell):

    degree = ell - 2
    result = sp.Integer(0)

    # Interior: a = 0,...,k-1
    for a in range(k):

        b = degree - a

        if b >= 0:
            result += (
                L2_interior(k, ell, a)
                * N**a
                * X**b
            )

    # Boundary a=k
    a = k
    b = degree - a

    if b >= 0:
        result += (
            L2_endpoint_k(k, ell)
            * N**a
            * X**b
        )

    # Boundary a=k+1
    a = k + 1
    b = degree - a

    if b >= 0:
        result += (
            L2_endpoint_kplus1(k, ell)
            * N**a
            * X**b
        )

    return sp.expand(result)


# ==============================================================================
# EXACT L3
# ==============================================================================

@lru_cache(maxsize=None)
def exact_L3(k, ell):

    G = exact_GX(k, ell)

    residual = sp.expand(
        G
        - top_layer(k, ell)
        - L1(k, ell)
        - L2(k, ell)
    )

    return homogeneous(
        residual,
        ell - 3,
    )


# ==============================================================================
# EXACT INTERIOR L3 FORMULA
# ==============================================================================

def predicted_components(k, a):

    choose = sp.binomial(
        k + 3,
        a,
    )

    A = (
        sp.Rational(1, 6)
        * choose
    )

    B = (
        -choose
        * (
            a * (k + 2)
            + k + 3
        )
        / (
            2 * (k + 3)
        )
    )

    C = (
        choose
        * (
            sp.Rational(1, 3)
            + (
                a
                * (
                    (k + 1) * a
                    + 2 * k
                    + 3
                )
                / (
                    2 * (k + 3)
                )
            )
        )
    )

    D = (
        -k
        * a
        * (a + 1)
        * (a + 2)
        * choose
        / (
            6 * (k + 3)
        )
    )

    return (
        sp.expand(A),
        sp.expand(B),
        sp.expand(C),
        sp.expand(D),
    )


def predicted_L3(k, ell, a):

    A, B, C, D = predicted_components(
        k,
        a,
    )

    return sp.expand(
        A * ell**3
        + B * ell**2
        + C * ell
        + D
    )


# ==============================================================================
# SAFE STRING FORMATTING
# ==============================================================================

def fmt(value):
    return str(sp.expand(value))


# ==============================================================================
# 1. DIRECT EXACT VALIDATION
# ==============================================================================

def direct_validation():

    print()
    print("=" * 78)
    print("1. EXACT INTERIOR L3 VALIDATION")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (3, 15),

        (5, 11),
        (5, 13),
        (5, 15),
        (5, 17),
        (5, 19),

        (7, 15),
        (7, 17),
        (7, 19),

        (9, 21),
        (9, 23),

        (11, 23),
        (11, 25),
    ]

    tested = 0
    failures = 0

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3
        local_failures = 0

        for a in range(k):

            b = degree - a

            if b < 0:
                continue

            actual = coefficient(
                exact,
                a,
                b,
            )

            predicted = predicted_L3(
                k,
                ell,
                a,
            )

            tested += 1

            if sp.expand(
                actual - predicted
            ) != 0:

                failures += 1
                local_failures += 1

                print()
                print(
                    f"FAIL k={k} ell={ell} a={a}"
                )
                print(
                    f"  actual    = {fmt(actual)}"
                )
                print(
                    f"  predicted = {fmt(predicted)}"
                )
                print(
                    f"  difference= "
                    f"{fmt(actual - predicted)}"
                )

        status = (
            "PASS"
            if local_failures == 0
            else "FAIL"
        )

        print(
            f"k={k:2d} "
            f"ell={ell:2d} "
            f"tested={k:2d} "
            f"failures={local_failures:2d} "
            f"{status}"
        )

    print()
    print(
        f"tested coefficients = {tested}"
    )
    print(
        f"exact formula failures = {failures}"
    )


# ==============================================================================
# 2. EXACT COMPONENT AUDIT
# ==============================================================================

def component_validation():

    print()
    print("=" * 78)
    print("2. EXACT A/B/C/D COMPONENT AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 13),
        (5, 11),
        (5, 17),
        (7, 15),
        (7, 19),
        (9, 21),
        (9, 23),
        (11, 23),
    ]

    failures = 0

    e = sp.symbols(
        "e",
        integer=True,
    )

    for k, ell0 in cases:

        print()
        print(
            f"k={k} ell0={ell0}"
        )

        for a in range(k):

            values = []

            # Need four exact ell samples.
            for j in range(4):

                ell = ell0 + 2 * j

                exact = exact_L3(
                    k,
                    ell,
                )

                b = ell - 3 - a

                if b < 0:
                    values.append(
                        sp.Integer(0)
                    )
                else:
                    values.append(
                        coefficient(
                            exact,
                            a,
                            b,
                        )
                    )

            interpolated = sp.interpolate(
                [
                    (
                        ell0 + 2 * j,
                        values[j],
                    )
                    for j in range(4)
                ],
                e,
            )

            interpolated = sp.expand(
                interpolated
            )

            poly = sp.Poly(
                interpolated,
                e,
                domain="QQ",
            )

            actual_A = poly.coeff_monomial(
                e**3
            )
            actual_B = poly.coeff_monomial(
                e**2
            )
            actual_C = poly.coeff_monomial(
                e
            )
            actual_D = poly.coeff_monomial(
                1
            )

            expected_A, expected_B, expected_C, expected_D = (
                predicted_components(
                    k,
                    a,
                )
            )

            checks = [
                ("A", actual_A, expected_A),
                ("B", actual_B, expected_B),
                ("C", actual_C, expected_C),
                ("D", actual_D, expected_D),
            ]

            row_failed = False

            for name, got, want in checks:

                ok = (
                    sp.expand(
                        got - want
                    )
                    == 0
                )

                if not ok:
                    failures += 1
                    row_failed = True

                print(
                    f"  a={a:2d} "
                    f"{name}: "
                    f"{fmt(got)} "
                    f"expected={fmt(want)} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

            if row_failed:
                print(
                    "    COMPONENT FAILURE"
                )

    print()
    print(
        f"A/B/C/D component failures = {failures}"
    )


# ==============================================================================
# 3. THIRD-DIFFERENCE AUDIT
# ==============================================================================

def finite_difference_audit():

    print()
    print("=" * 78)
    print("3. THIRD-DIFFERENCE AUDIT IN ELL")
    print("=" * 78)

    failures = 0

    cases = [
        (3, 7),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    for k, ell0 in cases:

        print()
        print(
            f"k={k}"
        )

        for a in range(k):

            values = []

            for j in range(4):

                ell = ell0 + 2 * j

                exact = exact_L3(
                    k,
                    ell,
                )

                b = ell - 3 - a

                values.append(
                    coefficient(
                        exact,
                        a,
                        b,
                    )
                )

            d1 = [
                sp.expand(
                    values[i + 1]
                    - values[i]
                )
                for i in range(3)
            ]

            d2 = [
                sp.expand(
                    d1[i + 1]
                    - d1[i]
                )
                for i in range(2)
            ]

            d3 = sp.expand(
                d2[1]
                - d2[0]
            )

            A = sp.expand(
                d3 / sp.Integer(48)
            )

            expected_A = sp.expand(
                sp.Rational(1, 6)
                * sp.binomial(
                    k + 3,
                    a,
                )
            )

            ok = (
                sp.expand(
                    A - expected_A
                )
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"  a={a:2d} "
                f"Delta3={fmt(d3):>8} "
                f"A={fmt(A)} "
                f"expected={fmt(expected_A)} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"third-difference failures = {failures}"
    )


# ==============================================================================
# 4. BOUNDARY EXTRACTION
# ==============================================================================

def boundary_extraction():

    print()
    print("=" * 78)
    print("4. L3 BOUNDARY EXTRACTION")
    print("=" * 78)

    cases = [
        (3, 9),
        (5, 13),
        (7, 17),
        (9, 23),
        (11, 25),
    ]

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3

        print()
        print(
            f"k={k} ell={ell}"
        )

        for a in range(
            max(0, k - 2),
            k + 6,
        ):

            b = degree - a

            if b < 0:
                continue

            value = coefficient(
                exact,
                a,
                b,
            )

            print(
                f"  a={a:2d} "
                f"[N^{a}X^{b}] = {fmt(value)}"
            )


# ==============================================================================
# 5. ELL-SEQUENCE OF THE TRUE BOUNDARY
# ==============================================================================

def boundary_sequences():

    print()
    print("=" * 78)
    print("5. L3 BOUNDARY ELL-SEQUENCES")
    print("=" * 78)

    configurations = [
        (3, [7, 9, 11, 13, 15]),
        (5, [11, 13, 15, 17, 19]),
        (7, [15, 17, 19]),
        (9, [21, 23, 25]),
        (11, [23, 25]),
    ]

    for k, ells in configurations:

        print()
        print(
            f"k={k}"
        )

        for a in range(
            max(0, k - 1),
            k + 5,
        ):

            values = []

            for ell in ells:

                exact = exact_L3(
                    k,
                    ell,
                )

                b = ell - 3 - a

                if b < 0:
                    values.append(
                        sp.Integer(0)
                    )
                else:
                    values.append(
                        coefficient(
                            exact,
                            a,
                            b,
                        )
                    )

            print(
                f"  a={a:2d}: "
                f"{[fmt(v) for v in values]}"
            )

            if len(values) >= 4:

                d1 = [
                    sp.expand(
                        values[i + 1]
                        - values[i]
                    )
                    for i in range(
                        len(values) - 1
                    )
                ]

                d2 = [
                    sp.expand(
                        d1[i + 1]
                        - d1[i]
                    )
                    for i in range(
                        len(d1) - 1
                    )
                ]

                print(
                    "       Delta2="
                    + str(
                        [fmt(v) for v in d2]
                    )
                )

                if len(d2) >= 2:

                    d3 = [
                        sp.expand(
                            d2[i + 1]
                            - d2[i]
                        )
                        for i in range(
                            len(d2) - 1
                        )
                    ]

                    print(
                        "       Delta3="
                        + str(
                            [fmt(v) for v in d3]
                        )
                    )


# ==============================================================================
# 6. A-LAW ONLY: FAST LARGE-K CHECK
# ==============================================================================

def a_law_fast_check():

    print()
    print("=" * 78)
    print("6. A-COEFFICIENT EXACT BINOMIAL LAW")
    print("=" * 78)

    failures = 0
    tested = 0

    # This section only evaluates the exact A coefficient from
    # four direct L3 samples. It does not print full polynomials.

    for k in range(
        3,
        14,
        2,
    ):

        ell0 = k + 4

        for a in range(k):

            e = sp.symbols(
                "e",
                integer=True,
            )

            points = []

            for j in range(4):

                ell = ell0 + 2 * j

                exact = exact_L3(
                    k,
                    ell,
                )

                b = ell - 3 - a

                value = coefficient(
                    exact,
                    a,
                    b,
                )

                points.append(
                    (
                        ell,
                        value,
                    )
                )

            polynomial = sp.expand(
                sp.interpolate(
                    points,
                    e,
                )
            )

            poly = sp.Poly(
                polynomial,
                e,
                domain="QQ",
            )

            A = poly.coeff_monomial(
                e**3
            )

            expected = sp.expand(
                sp.Rational(1, 6)
                * sp.binomial(
                    k + 3,
                    a,
                )
            )

            tested += 1

            if sp.expand(
                A - expected
            ) != 0:

                failures += 1

                print(
                    f"FAIL k={k} a={a} "
                    f"A={fmt(A)} "
                    f"expected={fmt(expected)}"
                )

    print()
    print(
        f"A-law tests = {tested}"
    )
    print(
        f"A-law failures = {failures}"
    )


# ==============================================================================
# 7. FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "All fractional constants are constructed with "
        "sp.Rational; no Python float division is used."
    )

    print(
        "All SymPy values are converted with str() before "
        "width formatting."
    )

    print()
    print(
        "The interior L3 candidate being tested is:"
    )

    print(
        "  A = C(k+3,a)/6"
    )

    print(
        "  B = -C(k+3,a)"
        " * (a(k+2)+k+3)"
        " / (2(k+3))"
    )

    print(
        "  C = C(k+3,a)"
        " * [1/3 + a((k+1)a+2k+3)/(2(k+3))]"
    )

    print(
        "  D = -k*a(a+1)(a+2)"
        " * C(k+3,a)/(6(k+3))"
    )

    print()
    print(
        "The only genuinely unresolved part is the boundary"
    )

    print(
        "  a >= k."
    )

    print()
    print(
        "The next derivation should therefore not alter the"
    )

    print(
        "interior formula merely because of boundary failures."
    )

    print(
        "Instead, derive the boundary rows directly from the"
    )

    print(
        "exact pq kernel."
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 221")
    print("EXACT L3 INTERIOR VALIDATION + BOUNDARY EXTRACTION")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted candidate enumeration is performed."
    )

    direct_validation()
    component_validation()
    finite_difference_audit()
    boundary_extraction()
    boundary_sequences()
    a_law_fast_check()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 221")
    print("=" * 78)


if __name__ == "__main__":
    main()