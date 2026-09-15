import sympy as sp
from functools import lru_cache

# ==============================================================================
# EXPERIMENT 224
# EXACT L3 BOUNDARY AUDIT — PURE QQ
# ==============================================================================

p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")
e = sp.symbols("e")


# ==============================================================================
# EXACTNESS GUARD
# ==============================================================================

def assert_QQ(expr, label="expression"):
    expr = sp.expand(expr)

    bad = expr.atoms(sp.Float)

    if bad:
        raise TypeError(
            f"FLOAT CONTAMINATION in {label}: {bad}"
        )

    return expr


def QQ(num, den=1):
    """
    Exact rational quotient.
    Never use Python '/' for mathematical coefficients.
    """
    return sp.cancel(
        sp.sympify(num) * sp.Rational(1, den)
    )


# ==============================================================================
# EXACT ORIGINAL KERNEL
# ==============================================================================

def exact_F(k, ell):
    k = sp.Integer(k)
    ell = sp.Integer(ell)

    return assert_QQ(
        sp.expand(
            p**k * (1 + q)**ell
            + q**k * (1 + p)**ell
            - p**ell * (1 + q)**k
            - q**ell * (1 + p)**k
        ),
        "F",
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
            f"non-symmetric remainder for "
            f"k={k}, ell={ell}: {remainder}"
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
            f"could not identify S,N mapping: {mapping}"
        )

    result = sp.expand(
        sym_expr.subs(
            {
                sum_symbol: S,
                prod_symbol: N,
            }
        )
    )

    return assert_QQ(
        result,
        "G(N,S)",
    )


@lru_cache(maxsize=None)
def exact_GX(k, ell):

    result = sp.expand(
        exact_G(k, ell).subs(
            S,
            X - 1,
        )
    )

    return assert_QQ(
        result,
        "G(N,X-1)",
    )


# ==============================================================================
# EXACT COEFFICIENT
# ==============================================================================

def coefficient(expr, a, b):

    if a < 0 or b < 0:
        return sp.Integer(0)

    expr = assert_QQ(
        expr,
        "coefficient input",
    )

    poly = sp.Poly(
        expr,
        N,
        X,
        domain="QQ",
    )

    value = poly.coeff_monomial(
        N**a * X**b
    )

    return assert_QQ(
        value,
        f"[N^{a}X^{b}]",
    )


# ==============================================================================
# HOMOGENEOUS PART
# ==============================================================================

def homogeneous(expr, degree):

    expr = assert_QQ(
        expr,
        "homogeneous input",
    )

    poly = sp.Poly(
        expr,
        N,
        X,
        domain="QQ",
    )

    result = sp.Integer(0)

    for (a, b), value in poly.terms():

        if a + b == degree:
            result += (
                value
                * N**a
                * X**b
            )

    return assert_QQ(
        sp.expand(result),
        f"homogeneous degree {degree}",
    )


# ==============================================================================
# TOP LAYER
# ==============================================================================

def top_layer(k, ell):

    k = sp.Integer(k)
    ell = sp.Integer(ell)

    result = sp.expand(
        -X**(ell - k)
        * (
            (X + N)**k
            - N**k
        )
    )

    return assert_QQ(
        result,
        "L0",
    )


# ==============================================================================
# L1
# ==============================================================================

def L1(k, ell):

    k = sp.Integer(k)
    ell = sp.Integer(ell)

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

    return assert_QQ(
        sp.expand(result),
        "L1",
    )


# ==============================================================================
# L2 INTERIOR
# ==============================================================================

def L2_interior(k, ell, a):

    k = sp.Integer(k)
    ell = sp.Integer(ell)
    a = sp.Integer(a)

    choose = sp.binomial(k + 2, a)

    A = QQ(
        -choose,
        2,
    )

    B = sp.cancel(
        choose
        * (
            2 * (k + 1) * a
            + k + 2
        )
        / (
            2 * (k + 2)
        )
    )

    C = sp.cancel(
        -k
        * a
        * (a + 1)
        * choose
        / (
            2 * (k + 2)
        )
    )

    result = (
        A * ell**2
        + B * ell
        + C
    )

    return assert_QQ(
        sp.expand(result),
        "L2 interior coefficient",
    )


def L2_k(k, ell):

    result = sp.expand(
        L2_interior(k, ell, k)
        + sp.binomial(k + 1, 2)
    )

    return assert_QQ(
        result,
        "L2 a=k",
    )


def L2_kplus1(k, ell):

    result = sp.cancel(
        QQ(2, k + 1)
        * L2_k(k, ell)
    )

    return assert_QQ(
        result,
        "L2 a=k+1",
    )


def L2(k, ell):

    ell = sp.Integer(ell)

    degree = ell - 2
    result = sp.Integer(0)

    for a in range(
        min(int(degree), k + 1) + 1
    ):

        b = degree - a

        if b < 0:
            continue

        if a < k:
            value = L2_interior(
                k,
                ell,
                a,
            )

        elif a == k:
            value = L2_k(
                k,
                ell,
            )

        else:
            value = L2_kplus1(
                k,
                ell,
            )

        result += (
            value
            * N**a
            * X**b
        )

    return assert_QQ(
        sp.expand(result),
        "L2",
    )


# ==============================================================================
# EXACT L3
# ==============================================================================

@lru_cache(maxsize=None)
def exact_L3(k, ell):

    result = sp.expand(
        exact_GX(k, ell)
        - top_layer(k, ell)
        - L1(k, ell)
        - L2(k, ell)
    )

    result = homogeneous(
        result,
        ell - 3,
    )

    return assert_QQ(
        result,
        "L3",
    )


# ==============================================================================
# INTERIOR L3
# ==============================================================================

def L3_interior(k, ell, a):

    k = sp.Integer(k)
    ell = sp.Integer(ell)
    a = sp.Integer(a)

    choose = sp.binomial(
        k + 3,
        a,
    )

    A = QQ(
        choose,
        6,
    )

    B = sp.cancel(
        -choose
        * (
            a * (k + 2)
            + k + 3
        )
        / (
            2 * (k + 3)
        )
    )

    C = sp.cancel(
        choose
        * (
            QQ(1, 3)
            + sp.cancel(
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

    D = sp.cancel(
        -k
        * a
        * (a + 1)
        * (a + 2)
        * choose
        / (
            6 * (k + 3)
        )
    )

    result = (
        A * ell**3
        + B * ell**2
        + C * ell
        + D
    )

    return assert_QQ(
        sp.expand(result),
        "L3 interior coefficient",
    )


# ==============================================================================
# BOUNDARY CANDIDATES
# ==============================================================================

def L3_k(k, ell):

    result = sp.expand(
        L3_interior(k, ell, k)
        + sp.binomial(k + 2, 3)
    )

    return assert_QQ(
        result,
        "L3 a=k",
    )


def L3_kplus1(k, ell):

    a = k + 1

    result = sp.expand(
        L3_interior(k, ell, a)
        - sp.binomial(k + 2, 2) * ell
        + QQ(
            k * (k + 2) * (k + 3),
            2,
        )
    )

    return assert_QQ(
        result,
        "L3 a=k+1",
    )


# ==============================================================================
# COMPLETE PREDICTION
# ==============================================================================

def predicted_L3(k, ell, a):

    k = int(k)
    ell = int(ell)
    a = int(a)

    if a < k:
        return L3_interior(
            k,
            ell,
            a,
        )

    if a == k:
        return L3_k(
            k,
            ell,
        )

    if a == k + 1:
        return L3_kplus1(
            k,
            ell,
        )

    return sp.Integer(0)


# ==============================================================================
# 1. FLOAT CONTAMINATION AUDIT
# ==============================================================================

def float_audit():

    print()
    print("=" * 78)
    print("1. FLOAT CONTAMINATION AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 13),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in cases:

        objects = {
            "G": exact_GX(k, ell),
            "L0": top_layer(k, ell),
            "L1": L1(k, ell),
            "L2": L2(k, ell),
            "L3": exact_L3(k, ell),
        }

        bad = {}

        for name, expr in objects.items():
            floats = expr.atoms(sp.Float)

            if floats:
                bad[name] = floats

        if bad:
            failures += 1
            print(
                f"k={k} ell={ell} "
                f"FLOATS={bad}"
            )
        else:
            print(
                f"k={k} ell={ell} "
                f"FLOATS=NONE PASS"
            )

    print()
    print(
        f"float-contamination failures = {failures}"
    )


# ==============================================================================
# 2. a=k CORRECTION
# ==============================================================================

def audit_k_boundary():

    print()
    print("=" * 78)
    print("2. a=k EXACT BOUNDARY CORRECTION")
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

        expected = sp.binomial(
            k + 2,
            3,
        )

        print()
        print(
            f"k={k} expected correction={expected}"
        )

        for ell in ells:

            exact = exact_L3(k, ell)

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

            ok = (
                correction
                == expected
            )

            if not ok:
                failures += 1

            print(
                f"  ell={ell:2d} "
                f"actual={str(actual):>8} "
                f"correction={str(correction):>8} "
                f"expected={str(expected):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"a=k correction failures = {failures}"
    )


# ==============================================================================
# 3. a=k+1 EXACT COMPONENT AUDIT
# ==============================================================================

def audit_kplus1():

    print()
    print("=" * 78)
    print("3. a=k+1 EXACT COMPONENT AUDIT")
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

        a = k + 1

        values = []

        for ell in ells:

            exact = exact_L3(
                k,
                ell,
            )

            b = ell - 3 - a

            values.append(
                (
                    sp.Integer(ell),
                    coefficient(
                        exact,
                        a,
                        b,
                    ),
                )
            )

        poly = sp.interpolate(
            values,
            e,
        )

        poly = sp.Poly(
            sp.expand(poly),
            e,
            domain="QQ",
        )

        actual_A = poly.coeff_monomial(e**3)
        actual_B = poly.coeff_monomial(e**2)
        actual_C = poly.coeff_monomial(e)
        actual_D = poly.coeff_monomial(1)

        interior = sp.Poly(
            L3_interior(
                k,
                e,
                a,
            ),
            e,
            domain="QQ",
        )

        int_A = interior.coeff_monomial(e**3)
        int_B = interior.coeff_monomial(e**2)
        int_C = interior.coeff_monomial(e)
        int_D = interior.coeff_monomial(1)

        expected_C = -sp.binomial(
            k + 2,
            2,
        )

        expected_D = QQ(
            k * (k + 2) * (k + 3),
            2,
        )

        tests = [
            ("A", actual_A, int_A),
            ("B", actual_B, int_B),
            (
                "C correction",
                sp.expand(actual_C - int_C),
                expected_C,
            ),
            (
                "D correction",
                sp.expand(actual_D - int_D),
                expected_D,
            ),
        ]

        print()
        print(
            f"k={k} a={a}"
        )

        for name, actual, expected in tests:

            ok = (
                sp.expand(
                    actual - expected
                )
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"  {name:14s} "
                f"actual={str(actual):>12} "
                f"expected={str(expected):>12} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"a=k+1 component failures = {failures}"
    )


# ==============================================================================
# 4. ZERO TAIL
# ==============================================================================

def audit_zero_tail():

    print()
    print("=" * 78)
    print("4. ZERO TAIL AUDIT")
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

        tail = []

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

            tail.append(
                (a, value)
            )

            if value != 0:
                failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"tail="
            f"{[(a, str(v)) for a, v in tail]}"
        )

    print()
    print(
        f"zero-tail failures = {failures}"
    )


# ==============================================================================
# 5. COMPLETE L3 AUDIT
# ==============================================================================

def complete_audit():

    print()
    print("=" * 78)
    print("5. COMPLETE L3 COEFFICIENT AUDIT")
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

    total = 0
    failures = 0

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

            total += 1

            if sp.expand(
                actual - predicted
            ) != 0:

                local_failures += 1
                failures += 1

                print(
                    f"FAIL k={k} ell={ell} a={a}"
                )

                print(
                    f"  actual   = {actual}"
                )

                print(
                    f"  predicted= {predicted}"
                )

                print(
                    f"  diff     = "
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
        f"tested coefficients = {total}"
    )
    print(
        f"complete formula failures = {failures}"
    )


# ==============================================================================
# 6. SYMBOLIC BOUNDARY FORMULAS
# ==============================================================================

def symbolic_boundary():

    print()
    print("=" * 78)
    print("6. SYMBOLIC BOUNDARY FORMULAS")
    print("=" * 78)

    K, L = sp.symbols(
        "K L",
        integer=True,
        nonnegative=True,
    )

    int_k = L3_interior(
        K,
        L,
        K,
    )

    exact_k = sp.expand(
        int_k
        + sp.binomial(
            K + 2,
            3,
        )
    )

    int_k1 = L3_interior(
        K,
        L,
        K + 1,
    )

    exact_k1 = sp.expand(
        int_k1
        - sp.binomial(
            K + 2,
            2,
        ) * L
        + QQ(
            K * (K + 2) * (K + 3),
            2,
        )
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
        "a >= k+2 = 0"
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 224")
    print("EXACT L3 BOUNDARY AUDIT — PURE QQ")
    print("=" * 78)

    print()
    print(
        "Every rational coefficient is constructed "
        "without Python floating-point division."
    )

    float_audit()
    audit_k_boundary()
    audit_kplus1()
    audit_zero_tail()
    complete_audit()
    symbolic_boundary()

    print()
    print("=" * 78)
    print("END EXPERIMENT 224")
    print("=" * 78)


if __name__ == "__main__":
    main()
