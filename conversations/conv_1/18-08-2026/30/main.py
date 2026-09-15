# ==============================================================================
# EXPERIMENT 242
# EXACT L2/L3 BOUNDARY LADDER AND GENERALIZATION AUDIT
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No exact_F supplied externally
# No L4 analysis
#
# Purpose:
#
#   1. Reconstruct the exact pq kernel.
#   2. Extract exact L2 and L3 homogeneous layers.
#   3. Verify the now-observed boundary corrections exactly.
#   4. Verify that support terminates at:
#
#        L2: a <= k+1
#        L3: a <= k+2
#
#   5. Test the boundary degree ladder:
#
#        Delta_r(k+j) has degree j in ell.
#
#   6. Test whether the terminal boundary row is simply
#      minus the interior extrapolation:
#
#        Delta_2(k+2) = -P_2(k,k+2,ell)
#        Delta_3(k+3) = -P_3(k,k+3,ell)
#
#   7. Compare the L2 and L3 ladders structurally.
#
# No L4 work is performed.
# ==============================================================================

import sympy as sp


# ==============================================================================
# SYMBOLS
# ==============================================================================

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")

K, A, L = sp.symbols(
    "K A L",
    integer=True,
    nonnegative=True,
)


# ==============================================================================
# EXACT HELPERS
# ==============================================================================

def simp(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(expr)
        )
    )


def exact_equal(a, b):
    return sp.expand(a - b) == 0


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
# EXACT SYMMETRIC TRANSFORMATION
#
# X = p + q + 1
# N = p q
# ==============================================================================

def exact_G(k, ell):
    F = exact_F(k, ell)

    sym_expr, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Non-symmetric remainder for k={k}, ell={ell}: "
            f"{remainder}"
        )

    s1 = mapping[0][0]
    s2 = mapping[1][0]

    G = sym_expr.subs(
        {
            s1: X - 1,
            s2: N,
        }
    )

    return sp.expand(G)


# ==============================================================================
# HOMOGENEOUS PART
# ==============================================================================

def homogeneous_part(expr, degree):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        X,
    )

    result = sp.Integer(0)

    for monom, coeff in poly.terms():
        a, b = monom

        if a + b == degree:
            result += coeff * N**a * X**b

    return sp.expand(result)


# ==============================================================================
# EXACT LAYERS
# ==============================================================================

def exact_L2(k, ell):
    G = exact_G(k, ell)
    return homogeneous_part(
        G,
        ell - 2,
    )


def exact_L3(k, ell):
    G = exact_G(k, ell)
    return homogeneous_part(
        G,
        ell - 3,
    )


# ==============================================================================
# COEFFICIENT
# ==============================================================================

def coeff_layer(layer, a, degree):
    b = degree - a

    if a < 0 or b < 0:
        return sp.Integer(0)

    return sp.expand(layer).coeff(
        N,
        a,
    ).coeff(
        X,
        b,
    )


# ==============================================================================
# INTERIOR PRODUCT LAW
# ==============================================================================

def P_interior(r, k, a, ell):
    product_part = sp.Integer(1)

    for j in range(1, r):
        product_part *= (
            ell - a - j
        )

    return simp(
        sp.Rational(
            (-1) ** (r + 1),
            sp.factorial(r),
        )
        * sp.binomial(
            k + r,
            a,
        )
        * product_part
        * (
            ((k + r) * ell - k * a)
            / (k + r)
        )
    )


# ==============================================================================
# OBSERVED EXACT L2 BOUNDARY FORMULAS
# ==============================================================================

def delta2_j0(k, ell):
    # a = k
    return sp.binomial(
        k + 1,
        2,
    )


def delta2_j1(k, ell):
    # a = k+1
    return (
        k * (k + 2)
        - (k + 1) * ell
    )


def delta2_j2(k, ell):
    # a = k+2
    return (
        (ell - k - 3)
        * (ell - k)
        / 2
    )


def delta2_expected(j, k, ell):
    if j == 0:
        return simp(
            delta2_j0(k, ell)
        )

    if j == 1:
        return simp(
            delta2_j1(k, ell)
        )

    if j == 2:
        return simp(
            delta2_j2(k, ell)
        )

    raise ValueError(
        f"No L2 formula for j={j}"
    )


# ==============================================================================
# OBSERVED EXACT L3 BOUNDARY FORMULAS
# ==============================================================================

def delta3_j0(k, ell):
    # a = k
    return sp.binomial(
        k + 2,
        3,
    )


def delta3_j1(k, ell):
    # a = k+1
    return simp(
        sp.Rational(
            k + 2,
            2,
        )
        * (
            k * (k + 3)
            - (k + 1) * ell
        )
    )


def delta3_j2(k, ell):
    # a = k+2
    return simp(
        (
            (ell - k - 4)
            * (
                (2 * k - 1) * ell
                - k * (k + 3)
            )
            / 2
        )
    )


def delta3_j3(k, ell):
    # a = k+3
    return simp(
        -(
            (ell - k)
            * (ell - k - 4)
            * (ell - k - 5)
            / 6
        )
    )


def delta3_expected(j, k, ell):
    if j == 0:
        return simp(
            delta3_j0(k, ell)
        )

    if j == 1:
        return simp(
            delta3_j1(k, ell)
        )

    if j == 2:
        return simp(
            delta3_j2(k, ell)
        )

    if j == 3:
        return simp(
            delta3_j3(k, ell)
        )

    raise ValueError(
        f"No L3 formula for j={j}"
    )


# ==============================================================================
# CASES
# ==============================================================================

CASES = {
    3:  [7, 9, 11, 13, 15, 17],
    5:  [11, 13, 15, 17, 19, 21],
    7:  [15, 17, 19, 21],
    9:  [21, 23, 25, 27],
    11: [23, 25, 27, 29],
    13: [25, 27, 29, 31],
    15: [27, 29, 31],
    17: [29, 31, 33],
}


# ==============================================================================
# CACHE
# ==============================================================================

L2_CACHE = {}
L3_CACHE = {}


def get_L2(k, ell):
    key = (k, ell)

    if key not in L2_CACHE:
        L2_CACHE[key] = exact_L2(
            k,
            ell,
        )

    return L2_CACHE[key]


def get_L3(k, ell):
    key = (k, ell)

    if key not in L3_CACHE:
        L3_CACHE[key] = exact_L3(
            k,
            ell,
        )

    return L3_CACHE[key]


# ==============================================================================
# EXTRACT ACTUAL DELTA
# ==============================================================================

def actual_delta2(k, ell, j):
    a = k + j

    L2 = get_L2(
        k,
        ell,
    )

    actual = coeff_layer(
        L2,
        a,
        ell - 2,
    )

    interior = P_interior(
        2,
        k,
        a,
        ell,
    )

    return simp(
        actual - interior
    )


def actual_delta3(k, ell, j):
    a = k + j

    L3 = get_L3(
        k,
        ell,
    )

    actual = coeff_layer(
        L3,
        a,
        ell - 3,
    )

    interior = P_interior(
        3,
        k,
        a,
        ell,
    )

    return simp(
        actual - interior
    )


# ==============================================================================
# 1. L2 BOUNDARY LAW AUDIT
# ==============================================================================

def audit_L2():
    print("=" * 78)
    print("1. EXACT L2 BOUNDARY LAW AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    for k, ells in CASES.items():
        for ell in ells:
            print(
                f"k={k:2d} ell={ell:2d}"
            )

            for j in range(3):
                actual = actual_delta2(
                    k,
                    ell,
                    j,
                )

                expected = delta2_expected(
                    j,
                    k,
                    ell,
                )

                residual = simp(
                    actual - expected
                )

                ok = residual == 0

                print(
                    f"  j={j} "
                    f"a={k+j:2d} "
                    f"actual={str(actual):>8} "
                    f"expected={str(expected):>8} "
                    f"residual={str(residual):>8} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

                tested += 1

                if not ok:
                    failures += 1

            print()

    print(
        f"L2 tested = {tested}"
    )
    print(
        f"L2 failures = {failures}"
    )
    print()


# ==============================================================================
# 2. L3 BOUNDARY LAW AUDIT
# ==============================================================================

def audit_L3():
    print("=" * 78)
    print("2. EXACT L3 BOUNDARY LAW AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    for k, ells in CASES.items():
        for ell in ells:
            print(
                f"k={k:2d} ell={ell:2d}"
            )

            for j in range(4):
                actual = actual_delta3(
                    k,
                    ell,
                    j,
                )

                expected = delta3_expected(
                    j,
                    k,
                    ell,
                )

                residual = simp(
                    actual - expected
                )

                ok = residual == 0

                print(
                    f"  j={j} "
                    f"a={k+j:2d} "
                    f"actual={str(actual):>8} "
                    f"expected={str(expected):>8} "
                    f"residual={str(residual):>8} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

                tested += 1

                if not ok:
                    failures += 1

            print()

    print(
        f"L3 tested = {tested}"
    )
    print(
        f"L3 failures = {failures}"
    )
    print()


# ==============================================================================
# 3. SUPPORT AUDIT
# ==============================================================================

def support_audit():
    print("=" * 78)
    print("3. EXACT BOUNDARY SUPPORT AUDIT")
    print("=" * 78)

    failures = 0

    for k, ells in CASES.items():
        for ell in ells:
            L2 = get_L2(
                k,
                ell,
            )

            L3 = get_L3(
                k,
                ell,
            )

            support2 = []

            for a in range(0, k + 5):
                if coeff_layer(
                    L2,
                    a,
                    ell - 2,
                ) != 0:
                    support2.append(a)

            support3 = []

            for a in range(0, k + 6):
                if coeff_layer(
                    L3,
                    a,
                    ell - 3,
                ) != 0:
                    support3.append(a)

            expected2 = list(
                range(
                    0,
                    min(k + 1, ell - 2) + 1,
                )
            )

            expected3 = list(
                range(
                    0,
                    min(k + 2, ell - 3) + 1,
                )
            )

            ok2 = (
                support2 == expected2
            )

            ok3 = (
                support3 == expected3
            )

            print(
                f"k={k:2d} ell={ell:2d} "
                f"L2={support2} "
                f"{'PASS' if ok2 else 'FAIL'}"
            )

            print(
                f"                 "
                f"L3={support3} "
                f"{'PASS' if ok3 else 'FAIL'}"
            )

            if not ok2:
                failures += 1

            if not ok3:
                failures += 1

    print()
    print(
        f"support failures = {failures}"
    )
    print()


# ==============================================================================
# 4. DEGREE LADDER AUDIT
# ==============================================================================

def finite_differences(values):
    rows = [
        list(values)
    ]

    current = list(values)

    while len(current) >= 2:
        current = [
            sp.expand(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        rows.append(current)

    return rows


def difference_order(values):
    rows = finite_differences(
        values
    )

    for order in range(
        1,
        len(rows),
    ):
        row = rows[order]

        if row and all(
            x == row[0]
            for x in row
        ):
            return order

    return None


def degree_ladder_audit():
    print("=" * 78)
    print("4. BOUNDARY DEGREE LADDER")
    print("=" * 78)

    for r in [2, 3]:
        print(
            f"L{r}"
        )

        for k, ells in CASES.items():
            print(
                f"  k={k}"
            )

            for j in range(r + 1):
                values = []

                for ell in ells:
                    if r == 2:
                        values.append(
                            actual_delta2(
                                k,
                                ell,
                                j,
                            )
                        )
                    else:
                        values.append(
                            actual_delta3(
                                k,
                                ell,
                                j,
                            )
                        )

                order = difference_order(
                    values
                )

                print(
                    f"    j={j} "
                    f"values={values} "
                    f"degree={order} "
                    f"expected={j} "
                    f"{'PASS' if order == j else 'FAIL'}"
                )

            print()


# ==============================================================================
# 5. TERMINAL-ROW CANCELLATION AUDIT
# ==============================================================================

def terminal_cancellation_audit():
    print("=" * 78)
    print("5. TERMINAL BOUNDARY CANCELLATION")
    print("=" * 78)

    failures = 0

    # L2 terminal j=2
    for k, ells in CASES.items():
        for ell in ells:
            actual = actual_delta2(
                k,
                ell,
                2,
            )

            expected = simp(
                -P_interior(
                    2,
                    k,
                    k + 2,
                    ell,
                )
            )

            residual = simp(
                actual - expected
            )

            ok = residual == 0

            print(
                f"L2 k={k:2d} ell={ell:2d} "
                f"residual={str(residual):>5} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()

    # L3 terminal j=3
    for k, ells in CASES.items():
        for ell in ells:
            actual = actual_delta3(
                k,
                ell,
                3,
            )

            expected = simp(
                -P_interior(
                    3,
                    k,
                    k + 3,
                    ell,
                )
            )

            residual = simp(
                actual - expected
            )

            ok = residual == 0

            print(
                f"L3 k={k:2d} ell={ell:2d} "
                f"residual={str(residual):>5} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(
        f"terminal cancellation failures = {failures}"
    )
    print()


# ==============================================================================
# 6. SYMBOLIC BOUNDARY LADDER IDENTITIES
# ==============================================================================

def symbolic_boundary_audit():
    print("=" * 78)
    print("6. SYMBOLIC BOUNDARY LADDER IDENTITIES")
    print("=" * 78)

    formulas2 = {
        0: delta2_j0(K, L),
        1: delta2_j1(K, L),
        2: delta2_j2(K, L),
    }

    formulas3 = {
        0: delta3_j0(K, L),
        1: delta3_j1(K, L),
        2: delta3_j2(K, L),
        3: delta3_j3(K, L),
    }

    for j, expr in formulas2.items():
        print(
            f"L2 Delta_(k+{j}) ="
        )
        print(
            f"  {sp.factor(expr)}"
        )
        print()

    for j, expr in formulas3.items():
        print(
            f"L3 Delta_(k+{j}) ="
        )
        print(
            f"  {sp.factor(expr)}"
        )
        print()

    print()


# ==============================================================================
# 7. L2 VS L3 COMPARISON
# ==============================================================================

def compare_ladders():
    print("=" * 78)
    print("7. L2 / L3 BOUNDARY LADDER COMPARISON")
    print("=" * 78)

    print(
        "L2:"
    )

    for j in range(3):
        print(
            f"  j={j}: "
            f"{sp.factor(delta2_expected(j, K, L))}"
        )

    print()

    print(
        "L3:"
    )

    for j in range(4):
        print(
            f"  j={j}: "
            f"{sp.factor(delta3_expected(j, K, L))}"
        )

    print()


# ==============================================================================
# 8. CROSS-k SYMBOLIC CONSISTENCY
# ==============================================================================

def cross_k_audit():
    print("=" * 78)
    print("8. CROSS-k CONSISTENCY")
    print("=" * 78)

    # Test symbolic formulas over a larger odd-k range
    failures = 0

    test_k = [
        3, 5, 7, 9,
        11, 13, 15, 17
    ]

    for k in test_k:
        for ell in CASES[k]:
            for j in range(3):
                actual = actual_delta2(
                    k,
                    ell,
                    j,
                )

                expected = delta2_expected(
                    j,
                    k,
                    ell,
                )

                if actual != expected:
                    failures += 1

            for j in range(4):
                actual = actual_delta3(
                    k,
                    ell,
                    j,
                )

                expected = delta3_expected(
                    j,
                    k,
                    ell,
                )

                if actual != expected:
                    failures += 1

    print(
        f"cross-k exact failures = {failures}"
    )
    print()


# ==============================================================================
# 9. FINAL STRUCTURAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():
    print("=" * 78)
    print("9. FINAL STRUCTURAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The exact kernel is now treated as the sole source
of the coefficients.

The established interior law is:

  P_r(k,a,L)
    = (-1)^(r+1)/r!
      * C(k+r,a)
      * product_{j=1}^{r-1}(L-a-j)
      * ((k+r)L-ka)/(k+r)

For the boundary layers currently known:

  L2:
    Delta(k)   = C(k+1,2)

    Delta(k+1)
      = k(k+2) - (k+1)L

    Delta(k+2)
      = (L-k-3)(L-k)/2

  L3:
    Delta(k)   = C(k+2,3)

    Delta(k+1)
      = (k+2)/2 * [k(k+3)-(k+1)L]

    Delta(k+2)
      = (L-k-4)
        * [(2k-1)L-k(k+3)]/2

    Delta(k+3)
      = -(L-k)(L-k-4)(L-k-5)/6

The important structural observation is:

  boundary offset j
      -> polynomial degree j in L

for every currently established case.

The terminal row also satisfies:

  Delta_2(k+2) = -P_2(k,k+2,L)

  Delta_3(k+3) = -P_3(k,k+3,L)

The next mathematical question is therefore no longer
whether these individual formulas are correct.

The next question is:

  WHAT IS THE GENERAL FORMULA FOR
  Delta_r(k+j,L)?

In particular, determine whether there is a universal
closed boundary product analogous to the interior product

  product_{j=1}^{r-1}(L-a-j).

Do not begin L4 until a genuine general boundary
mechanism has been identified or the derivation shows
exactly why the pattern terminates.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 242")
    print("EXACT L2/L3 BOUNDARY LADDER AND GENERALIZATION AUDIT")
    print("=" * 78)
    print()
    print("Standalone main.py")
    print("Exact arithmetic over QQ")
    print("No previous experiment imported")
    print("No filesystem access")
    print("No L4 analysis")
    print()

    audit_L2()
    audit_L3()
    support_audit()
    degree_ladder_audit()
    terminal_cancellation_audit()
    symbolic_boundary_audit()
    compare_ladders()
    cross_k_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()

