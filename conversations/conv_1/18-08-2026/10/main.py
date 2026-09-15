import sympy as sp
from functools import lru_cache

# ==============================================================================
# EXPERIMENT 222
# EXACT L3 BOUNDARY DERIVATION: a=k AND a=k+1
# ==============================================================================
#
# No previous experiment output is read.
# No unrestricted candidate enumeration is performed.
#
# Established facts being used:
#
#   G_top = -X^(ell-k) ((X+N)^k - N^k)
#
#   L1 is exact.
#
#   L2 is exact.
#
#   L3 is the homogeneous degree-(ell-3) part of
#
#       G - G_top - L1 - L2.
#
# The current goal is ONLY:
#
#   1. prove the interior formula at a=k differs by a boundary correction;
#   2. derive the exact a=k correction;
#   3. derive the exact a=k+1 row;
#   4. prove a>=k+2 vanishes;
#   5. reconstruct L3 completely for the tested cases.
#
# Everything is exact over QQ.
# ==============================================================================

p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")
e = sp.symbols("e")


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
        + sp.Rational(
            k * (k + 1),
            2,
        )
    )


def L2_kplus1(k, ell):

    return sp.expand(
        sp.Rational(
            2,
            k + 1,
        )
        * L2_k(k, ell)
    )


def L2(k, ell):

    degree = ell - 2
    result = sp.Integer(0)

    # a = 0,...,k-1
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
# INTERIOR L3 FORMULA
# ==============================================================================

def L3_interior(k, ell, a):

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
# a = k BOUNDARY:
# COMPARE EXACT ROW TO INTERIOR EXTRAPOLATION
# ==============================================================================

def audit_a_k():

    print()
    print("=" * 78)
    print("1. a=k BOUNDARY CORRECTION")
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

        corrections = []

        print()
        print(f"k={k}")

        for ell in ells:

            exact = exact_L3(
                k,
                ell,
            )

            b = ell - 3 - k

            actual = coefficient(
                exact,
                k,
                b,
            )

            naive = L3_interior(
                k,
                ell,
                k,
            )

            correction = sp.expand(
                actual - naive
            )

            corrections.append(
                correction
            )

            print(
                f"  ell={ell:2d}"
                f"  actual={str(actual):>8}"
                f"  naive={str(naive):>8}"
                f"  correction={str(correction):>8}"
            )

        unique = sorted(
            set(corrections),
            key=str,
        )

        predicted = sp.binomial(
            k + 2,
            3,
        )

        print(
            f"  unique corrections = "
            f"{[str(x) for x in unique]}"
        )

        print(
            f"  C(k+2,3) = {str(predicted)}"
        )

        print(
            "  correction law = "
            + (
                "PASS"
                if all(
                    sp.expand(
                        c - predicted
                    ) == 0
                    for c in corrections
                )
                else "FAIL"
            )
        )


# ==============================================================================
# a = k+1 BOUNDARY
# ==============================================================================
#
# We do NOT assume a formula.
#
# First extract the exact cubic in ell from four or five samples.
# Then inspect its coefficients.
# ==============================================================================

def interpolate_boundary_row(k, a, ells):

    points = []

    for ell in ells:

        exact = exact_L3(
            k,
            ell,
        )

        b = ell - 3 - a

        if b < 0:
            value = sp.Integer(0)
        else:
            value = coefficient(
                exact,
                a,
                b,
            )

        points.append(
            (ell, value)
        )

    polynomial = sp.interpolate(
        points,
        e,
    )

    return sp.factor(
        sp.expand(polynomial)
    )


def audit_a_kplus1():

    print()
    print("=" * 78)
    print("2. a=k+1 EXACT BOUNDARY")
    print("=" * 78)

    cases = {
        3: [7, 9, 11, 13, 15],
        5: [11, 13, 15, 17, 19],
        7: [15, 17, 19, 21],
        9: [21, 23, 25, 27],
        11: [23, 25, 27, 29],
    }

    for k, ells in cases.items():

        a = k + 1

        polynomial = interpolate_boundary_row(
            k,
            a,
            ells,
        )

        degree = sp.degree(
            polynomial,
            e,
        )

        print()
        print(
            f"k={k} a={a}"
        )
        print(
            f"  polynomial = {polynomial}"
        )
        print(
            f"  degree     = {degree}"
        )

        if degree <= 3:

            poly = sp.Poly(
                polynomial,
                e,
                domain="QQ",
            )

            A = poly.coeff_monomial(e**3)
            B = poly.coeff_monomial(e**2)
            C = poly.coeff_monomial(e)
            D = poly.coeff_monomial(1)

            print(
                f"  A = {str(A)}"
            )
            print(
                f"  B = {str(B)}"
            )
            print(
                f"  C = {str(C)}"
            )
            print(
                f"  D = {str(D)}"
            )


# ==============================================================================
# SUPPORT: TEST a >= k+2
# ==============================================================================

def audit_zero_tail():

    print()
    print("=" * 78)
    print("3. ZERO-TAIL AUDIT: a >= k+2")
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
            min(
                degree + 1,
                k + 7,
            ),
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
        f"nonzero-tail failures = {failures}"
    )


# ==============================================================================
# a=k RELATION TO INTERIOR FORMULA
# ==============================================================================

def prove_a_k_formula():

    print()
    print("=" * 78)
    print("4. a=k CLOSED FORM")
    print("=" * 78)

    ksym, L = sp.symbols(
        "k L",
        integer=True,
        nonnegative=True,
    )

    # Interior formula extrapolated to a=k.
    naive = sp.expand(
        L3_interior(
            ksym,
            L,
            ksym,
        )
    )

    correction = sp.binomial(
        ksym + 2,
        3,
    )

    corrected = sp.expand(
        naive + correction
    )

    print()
    print(
        "interior extrapolation at a=k:"
    )
    print(
        sp.factor(naive)
    )

    print()
    print(
        "boundary correction:"
    )
    print(
        sp.factor(correction)
    )

    print()
    print(
        "candidate exact a=k:"
    )
    print(
        sp.factor(corrected)
    )


# ==============================================================================
# CROSS-CHECK COMPLETE L3
# ==============================================================================

def complete_reconstruction():

    print()
    print("=" * 78)
    print("5. COMPLETE L3 RECONSTRUCTION")
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
    ]

    failures = 0

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3
        reconstructed = sp.Integer(0)

        # Interior 0 <= a < k
        for a in range(k):

            b = degree - a

            if b >= 0:
                reconstructed += (
                    L3_interior(
                        k,
                        ell,
                        a,
                    )
                    * N**a
                    * X**b
                )

        # a=k
        a = k
        b = degree - a

        if b >= 0:
            reconstructed += (
                (
                    L3_interior(
                        k,
                        ell,
                        k,
                    )
                    + sp.binomial(
                        k + 2,
                        3,
                    )
                )
                * N**a
                * X**b
            )

        # a=k+1 is deliberately NOT inserted yet.
        # We test how much remains.
        difference = sp.expand(
            exact - reconstructed
        )

        residual_terms = []

        poly = sp.Poly(
            difference,
            N,
            X,
            domain="QQ",
        )

        for monom, value in poly.terms():

            a, b = monom

            if a + b == degree:
                residual_terms.append(
                    (
                        a,
                        b,
                        sp.expand(value),
                    )
                )

        print()
        print(
            f"k={k:2d} ell={ell:2d}"
        )
        print(
            "  residual boundary terms = "
            + str(
                [
                    (
                        a,
                        b,
                        str(v),
                    )
                    for a, b, v in residual_terms
                ]
            )
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 222")
    print("EXACT L3 BOUNDARY DERIVATION")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted candidate enumeration is performed."
    )

    audit_a_k()
    audit_a_kplus1()
    audit_zero_tail()
    prove_a_k_formula()
    complete_reconstruction()

    print()
    print("=" * 78)
    print("END EXPERIMENT 222")
    print("=" * 78)


if __name__ == "__main__":
    main()

