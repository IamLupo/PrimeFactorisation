import sympy as sp
from functools import lru_cache

# ==============================================================================
# EXPERIMENT 223
# COMPLETE EXACT L3 CLOSED FORM
# ==============================================================================

# No previous experiment output is read.
# No unrestricted candidate enumeration is performed.
#
# Established:
#
#   G_top = -X^(ell-k) ((X+N)^k - N^k)
#
#   L1 is exact.
#
#   L2 is exact.
#
# The remaining task is to prove the COMPLETE L3 structure:
#
#   0 <= a < k:
#       interior formula
#
#   a = k:
#       interior formula + boundary correction
#
#   a = k+1:
#       same A,B as interior extrapolation,
#       corrected C,D
#
#   a >= k+2:
#       zero
#
# Everything is exact over QQ.
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

    remainder = sp.expand(remainder)

    if remainder != 0:
        raise ValueError(
            f"Symmetric reduction failed "
            f"for k={k}, ell={ell}: {remainder}"
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
            f"Could not identify symmetric variables: "
            f"mapping={mapping}"
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
# EXACT COEFFICIENT EXTRACTION
# ==============================================================================

def coefficient(expr, a, b):

    if a < 0 or b < 0:
        return sp.Integer(0)

    poly = sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain="QQ",
    )

    return sp.expand(
        poly.coeff_monomial(
            N**a * X**b
        )
    )


# ==============================================================================
# HOMOGENEOUS PART
# ==============================================================================

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
# L1
# ==============================================================================

def L1(k, ell):

    degree = ell - 1
    result = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        value = (
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
# L2 INTERIOR
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


def L2_k(k, ell):

    return sp.expand(
        L2_interior(k, ell, k)
        + sp.binomial(k + 1, 2)
    )


def L2_kplus1(k, ell):

    return sp.expand(
        sp.Rational(2, k + 1)
        * L2_k(k, ell)
    )


def L2(k, ell):

    degree = ell - 2
    result = sp.Integer(0)

    # 0 <= a < k
    for a in range(k):

        b = degree - a

        if b >= 0:
            result += (
                L2_interior(k, ell, a)
                * N**a
                * X**b
            )

    # a = k
    a = k
    b = degree - a

    if b >= 0:
        result += (
            L2_k(k, ell)
            * N**a
            * X**b
        )

    # a = k+1
    a = k + 1
    b = degree - a

    if b >= 0:
        result += (
            L2_kplus1(k, ell)
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
# L3 INTERIOR COEFFICIENT
# ==============================================================================

def L3_interior_coeff(k, ell, a):

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

    return sp.expand(
        A * ell**3
        + B * ell**2
        + C * ell
        + D
    )


# ==============================================================================
# L3 a=k EXACT BOUNDARY
# ==============================================================================

def L3_k_coeff(k, ell):

    return sp.expand(
        L3_interior_coeff(
            k,
            ell,
            k,
        )
        + sp.binomial(
            k + 2,
            3,
        )
    )


# ==============================================================================
# L3 a=k+1 EXACT BOUNDARY
# ==============================================================================

def L3_kplus1_coeff(k, ell):

    a = k + 1

    interior = L3_interior_coeff(
        k,
        ell,
        a,
    )

    correction_C = -sp.binomial(
        k + 2,
        2,
    )

    correction_D = (
        sp.Rational(1, 2)
        * k
        * (k + 2)
        * (k + 3)
    )

    # A and B stay unchanged.
    # C gets -(k+1)(k+2)/2.
    # D gets +k(k+2)(k+3)/2.

    return sp.expand(
        interior
        + correction_C * ell
        + correction_D
    )


# ==============================================================================
# 1. A=k CORRECTION AUDIT
# ==============================================================================

def audit_k_boundary():

    print()
    print("=" * 78)
    print("1. a=k EXACT BOUNDARY CORRECTION")
    print("=" * 78)

    cases = {
        3: [7, 9, 11, 13, 15],
        5: [11, 13, 15, 17, 19],
        7: [15, 17, 19],
        9: [21, 23, 25],
        11: [23, 25],
        13: [25, 27],
    }

    failures = 0

    for k, ells in cases.items():

        expected_correction = sp.binomial(
            k + 2,
            3,
        )

        print()
        print(f"k={k}")

        for ell in ells:

            exact = exact_L3(k, ell)

            b = ell - 3 - k

            actual = coefficient(
                exact,
                k,
                b,
            )

            naive = L3_interior_coeff(
                k,
                ell,
                k,
            )

            correction = sp.expand(
                actual - naive
            )

            ok = (
                sp.expand(
                    correction
                    - expected_correction
                )
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"  ell={ell:2d} "
                f"actual={str(actual):>8} "
                f"correction={str(correction):>8} "
                f"expected={str(expected_correction):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"a=k correction failures = {failures}"
    )


# ==============================================================================
# 2. a=k+1 COMPONENT AUDIT
# ==============================================================================

def audit_kplus1_components():

    print()
    print("=" * 78)
    print("2. a=k+1 COMPONENT STRUCTURE")
    print("=" * 78)

    cases = {
        3: [7, 9, 11, 13, 15],
        5: [11, 13, 15, 17, 19],
        7: [15, 17, 19],
        9: [21, 23, 25],
        11: [23, 25],
        13: [25, 27],
    }

    for k, ells in cases.items():

        a = k + 1

        print()
        print(
            f"k={k} a={a}"
        )

        poly_values = []

        for ell in ells:

            exact = exact_L3(k, ell)

            b = ell - 3 - a

            actual = coefficient(
                exact,
                a,
                b,
            )

            poly_values.append(
                (sp.Integer(ell), actual)
            )

        polynomial = sp.factor(
            sp.interpolate(
                poly_values,
                e,
            )
        )

        print(
            f"  exact cubic = {polynomial}"
        )

        poly = sp.Poly(
            sp.expand(polynomial),
            e,
            domain="QQ",
        )

        A = poly.coeff_monomial(e**3)
        B = poly.coeff_monomial(e**2)
        C = poly.coeff_monomial(e)
        D = poly.coeff_monomial(1)

        interior_poly = sp.Poly(
            sp.expand(
                L3_interior_coeff(
                    k,
                    e,
                    a,
                )
            ),
            e,
            domain="QQ",
        )

        Ai = interior_poly.coeff_monomial(e**3)
        Bi = interior_poly.coeff_monomial(e**2)
        Ci = interior_poly.coeff_monomial(e)
        Di = interior_poly.coeff_monomial(1)

        print(
            f"  A actual={str(A):>10} "
            f"interior={str(Ai):>10} "
            f"{'PASS' if A == Ai else 'FAIL'}"
        )

        print(
            f"  B actual={str(B):>10} "
            f"interior={str(Bi):>10} "
            f"{'PASS' if B == Bi else 'FAIL'}"
        )

        expected_C_correction = -sp.binomial(
            k + 2,
            2,
        )

        expected_D_correction = (
            sp.Rational(1, 2)
            * k
            * (k + 2)
            * (k + 3)
        )

        actual_C_correction = sp.expand(
            C - Ci
        )

        actual_D_correction = sp.expand(
            D - Di
        )

        print(
            "  C correction="
            f"{str(actual_C_correction):>8} "
            "expected="
            f"{str(expected_C_correction):>8} "
            f"{'PASS' if actual_C_correction == expected_C_correction else 'FAIL'}"
        )

        print(
            "  D correction="
            f"{str(actual_D_correction):>8} "
            "expected="
            f"{str(expected_D_correction):>8} "
            f"{'PASS' if actual_D_correction == expected_D_correction else 'FAIL'}"
        )


# ==============================================================================
# 3. ZERO TAIL
# ==============================================================================

def audit_zero_tail():

    print()
    print("=" * 78)
    print("3. ZERO TAIL: a >= k+2")
    print("=" * 78)

    cases = [
        (3, 9),
        (5, 13),
        (7, 17),
        (9, 23),
        (11, 25),
        (13, 27),
    ]

    failures = 0

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3

        row = []

        for a in range(
            k + 2,
            degree + 1,
        ):

            b = degree - a

            value = coefficient(
                exact,
                a,
                b,
            )

            row.append(
                (a, value)
            )

            if value != 0:
                failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"tail="
            f"{[(a, str(v)) for a, v in row]}"
        )

    print()
    print(
        f"zero-tail failures = {failures}"
    )


# ==============================================================================
# 4. COMPLETE L3 FORMULA
# ==============================================================================

def predicted_L3(k, ell, a):

    if a < 0:
        return sp.Integer(0)

    if a < k:
        return L3_interior_coeff(
            k,
            ell,
            a,
        )

    if a == k:
        return L3_k_coeff(
            k,
            ell,
        )

    if a == k + 1:
        return L3_kplus1_coeff(
            k,
            ell,
        )

    return sp.Integer(0)


# ==============================================================================
# 5. COMPLETE COEFFICIENT AUDIT
# ==============================================================================

def complete_coefficient_audit():

    print()
    print("=" * 78)
    print("4. COMPLETE L3 COEFFICIENT AUDIT")
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
        (9, 25),

        (11, 23),
        (11, 25),

        (13, 25),
        (13, 27),
    ]

    failures = 0
    tested = 0

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3

        local_failures = 0

        for a in range(
            degree + 1
        ):

            b = degree - a

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

                local_failures += 1
                failures += 1

                print(
                    f"FAIL k={k} ell={ell} "
                    f"a={a}"
                )
                print(
                    f"  actual    = {actual}"
                )
                print(
                    f"  predicted = {predicted}"
                )
                print(
                    f"  difference= "
                    f"{sp.expand(actual - predicted)}"
                )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"tested={degree + 1:2d} "
            f"failures={local_failures:2d} "
            f"{'PASS' if local_failures == 0 else 'FAIL'}"
        )

    print()
    print(
        f"tested coefficients = {tested}"
    )
    print(
        f"complete formula failures = {failures}"
    )


# ==============================================================================
# 6. EXACT HOMOGENEOUS RECONSTRUCTION
# ==============================================================================

def reconstruction_audit():

    print()
    print("=" * 78)
    print("5. EXACT L3 RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 13),
        (5, 15),
        (7, 15),
        (7, 17),
        (9, 21),
        (9, 23),
        (11, 23),
        (11, 25),
        (13, 25),
    ]

    failures = 0

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3

        predicted = sp.Integer(0)

        for a in range(
            degree + 1
        ):

            b = degree - a

            value = predicted_L3(
                k,
                ell,
                a,
            )

            if value != 0:
                predicted += (
                    value
                    * N**a
                    * X**b
                )

        difference = sp.expand(
            exact - predicted
        )

        ok = (
            difference == 0
        )

        if not ok:
            failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction="
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            print(
                f"  difference = {difference}"
            )

    print()
    print(
        f"L3 reconstruction failures = {failures}"
    )


# ==============================================================================
# 7. SYMBOLIC BOUNDARY FORMULAS
# ==============================================================================

def symbolic_boundary_formulas():

    print()
    print("=" * 78)
    print("6. SYMBOLIC L3 BOUNDARY FORMULAS")
    print("=" * 78)

    K, L = sp.symbols(
        "K L",
        integer=True,
        nonnegative=True,
    )

    # a = K
    interior_k = sp.expand(
        L3_interior_coeff(
            K,
            L,
            K,
        )
    )

    exact_k = sp.expand(
        interior_k
        + sp.binomial(
            K + 2,
            3,
        )
    )

    # a = K+1
    interior_k1 = sp.expand(
        L3_interior_coeff(
            K,
            L,
            K + 1,
        )
    )

    exact_k1 = sp.expand(
        interior_k1
        - sp.binomial(
            K + 2,
            2,
        )
        * L
        + sp.Rational(
            1,
            2,
        )
        * K
        * (K + 2)
        * (K + 3)
    )

    print()
    print("a = k:")
    print(
        sp.factor(exact_k)
    )

    print()
    print("a = k+1:")
    print(
        sp.factor(exact_k1)
    )

    print()
    print(
        "a >= k+2: 0"
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 223")
    print("COMPLETE EXACT L3 CLOSED FORM")
    print("=" * 78)

    print()
    print(
        "All arithmetic is exact over QQ."
    )
    print(
        "No floating-point equality tests are used."
    )
    print(
        "No unrestricted candidate enumeration is performed."
    )

    audit_k_boundary()
    audit_kplus1_components()
    audit_zero_tail()
    complete_coefficient_audit()
    reconstruction_audit()
    symbolic_boundary_formulas()

    print()
    print("=" * 78)
    print("END EXPERIMENT 223")
    print("=" * 78)


if __name__ == "__main__":
    main()
