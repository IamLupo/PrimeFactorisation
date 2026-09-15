# ==============================================================================
# EXPERIMENT 241
# EXACT KERNEL DERIVATION OF L2/L3 BOUNDARY CORRECTIONS
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No exact_F supplied from another file
#
# Exact kernel:
#
#   F =
#       p^k (1+q)^ell
#     + q^k (1+p)^ell
#     - p^ell (1+q)^k
#     - q^ell (1+p)^k
#
# Variables:
#
#   N = p q
#   X = p + q + 1
#
# Objective:
#
#   1. Reconstruct exact G(N,X).
#   2. Extract L2 and L3 directly.
#   3. For the boundary rows, compute
#
#        Delta_2(a) = actual L2 coefficient
#                     - interior P2 extrapolation
#
#        Delta_3(a) = actual L3 coefficient
#                     - interior P3 extrapolation
#
#   4. Determine the degree in ell from exact finite differences.
#   5. Factor the resulting boundary polynomials.
#   6. Test cross-k structure.
#
# No boundary formula is assumed to be correct.
# No L4 analysis is performed.
# ==============================================================================

import sympy as sp


# ==============================================================================
# SYMBOLS
# ==============================================================================

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")

K, A, L = sp.symbols("K A L", integer=True, nonnegative=True)


# ==============================================================================
# BASIC EXACT UTILITIES
# ==============================================================================

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


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
# pq -> (N,X) TRANSFORMATION
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
            f"Non-symmetric remainder for k={k}, ell={ell}: {remainder}"
        )

    s1 = mapping[0][0]   # p+q
    s2 = mapping[1][0]   # pq

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
    poly = sp.Poly(sp.expand(expr), N, X)

    result = sp.Integer(0)

    for (a, b), coeff in poly.terms():
        if a + b == degree:
            result += coeff * N**a * X**b

    return sp.expand(result)


# ==============================================================================
# EXACT LAYERS
# ==============================================================================

def exact_layers(k, ell):
    G = exact_G(k, ell)

    L2 = homogeneous_part(G, ell - 2)
    L3 = homogeneous_part(G, ell - 3)

    return G, L2, L3


# ==============================================================================
# COEFFICIENT EXTRACTION
# ==============================================================================

def coefficient(layer, a, degree):
    b = degree - a

    if a < 0 or b < 0:
        return sp.Integer(0)

    return sp.expand(layer).coeff(N, a).coeff(X, b)


# ==============================================================================
# ESTABLISHED INTERIOR PRODUCT LAW
# ==============================================================================

def P_interior(r, k, a, ell):
    k = sp.sympify(k)
    a = sp.sympify(a)
    ell = sp.sympify(ell)

    product_part = sp.Integer(1)

    for j in range(1, r):
        product_part *= ell - a - j

    return simp(
        sp.Rational((-1) ** (r + 1), sp.factorial(r))
        * sp.binomial(k + r, a)
        * product_part
        * (((k + r) * ell - k * a) / (k + r))
    )


# ==============================================================================
# SAMPLE POINTS
# ==============================================================================

CASES = {
    3:  [7, 9, 11, 13, 15, 17],
    5:  [11, 13, 15, 17, 19, 21],
    7:  [15, 17, 19, 21],
    9:  [21, 23, 25, 27],
    11: [23, 25, 27, 29],
    13: [25, 27, 29, 31],
}


# ==============================================================================
# CACHE
# ==============================================================================

LAYER_CACHE = {}


def get_layers(k, ell):
    key = (k, ell)

    if key not in LAYER_CACHE:
        LAYER_CACHE[key] = exact_layers(k, ell)

    return LAYER_CACHE[key]


# ==============================================================================
# 1. KERNEL / TOP SANITY
# ==============================================================================

def sanity_audit():
    print("=" * 78)
    print("1. EXACT KERNEL / TOP SANITY")
    print("=" * 78)

    failures = 0

    for k, ells in CASES.items():
        ell = ells[0]

        F = exact_F(k, ell)

        symmetric = exact_equal(
            F,
            F.xreplace({p: q, q: p}),
        )

        _, _, _ = get_layers(k, ell)

        G = get_layers(k, ell)[0]

        L0 = homogeneous_part(G, ell)

        expected_top = (
            -X**(ell - k)
            * ((X + N)**k - N**k)
        )

        top_ok = exact_equal(L0, expected_top)

        ok = symmetric and top_ok

        print(
            f"k={k:2d} ell={ell:2d} "
            f"symmetric={symmetric} "
            f"top={top_ok} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"sanity failures = {failures}")
    print()


# ==============================================================================
# 2. INTERIOR L2/L3 CONTROL
# ==============================================================================

def interior_control():
    print("=" * 78)
    print("2. INTERIOR L2/L3 CONTROL")
    print("=" * 78)

    failures = 0
    tested = 0

    for k, ells in CASES.items():
        for ell in ells[:4]:
            _, L2, L3 = get_layers(k, ell)

            local = 0

            for r, layer in [(2, L2), (3, L3)]:
                degree = ell - r

                for a in range(k):
                    actual = coefficient(layer, a, degree)
                    expected = P_interior(r, k, a, ell)

                    tested += 1

                    if not exact_equal(actual, expected):
                        local += 1
                        failures += 1

            print(
                f"k={k:2d} ell={ell:2d} "
                f"interior failures={local} "
                f"{'PASS' if local == 0 else 'FAIL'}"
            )

    print()
    print(f"tested = {tested}")
    print(f"interior failures = {failures}")
    print()


# ==============================================================================
# 3. EXTRACT L2 BOUNDARY CORRECTIONS
# ==============================================================================

def extract_L2_boundary():
    print("=" * 78)
    print("3. EXACT L2 BOUNDARY CORRECTIONS")
    print("=" * 78)

    for k, ells in CASES.items():
        print(f"k={k}")

        for ell in ells:
            _, L2, _ = get_layers(k, ell)
            degree = ell - 2

            print(f"  ell={ell}")

            for offset in [0, 1, 2]:
                a = k + offset

                actual = coefficient(L2, a, degree)
                interior = P_interior(2, k, a, ell)

                delta = simp(actual - interior)

                print(
                    f"    a={a:2d} "
                    f"actual={str(actual):>8} "
                    f"interior={str(interior):>8} "
                    f"Delta={str(delta):>8}"
                )

        print()


# ==============================================================================
# 4. EXTRACT L3 BOUNDARY CORRECTIONS
# ==============================================================================

def extract_L3_boundary():
    print("=" * 78)
    print("4. EXACT L3 BOUNDARY CORRECTIONS")
    print("=" * 78)

    for k, ells in CASES.items():
        print(f"k={k}")

        for ell in ells:
            _, _, L3 = get_layers(k, ell)
            degree = ell - 3

            print(f"  ell={ell}")

            for offset in [0, 1, 2, 3]:
                a = k + offset

                actual = coefficient(L3, a, degree)
                interior = P_interior(3, k, a, ell)

                delta = simp(actual - interior)

                print(
                    f"    a={a:2d} "
                    f"actual={str(actual):>8} "
                    f"interior={str(interior):>8} "
                    f"Delta={str(delta):>8}"
                )

        print()


# ==============================================================================
# 5. FINITE-DIFFERENCE HELPERS
# ==============================================================================

def finite_differences(values):
    rows = [list(values)]

    current = list(values)

    while len(current) >= 2:
        current = [
            sp.expand(current[i + 1] - current[i])
            for i in range(len(current) - 1)
        ]
        rows.append(current)

    return rows


def first_constant_difference_order(values):
    rows = finite_differences(values)

    for order, row in enumerate(rows[1:], start=1):
        if row and all(x == row[0] for x in row):
            return order

    if len(values) <= 1:
        return 0

    return None


# ==============================================================================
# 6. L2 DEGREE AUDIT
# ==============================================================================

def L2_degree_audit():
    print("=" * 78)
    print("6. L2 BOUNDARY DEGREE IN ell")
    print("=" * 78)

    for k, ells in CASES.items():
        print(f"k={k}")

        for offset in [0, 1, 2]:
            a = k + offset

            values = []

            for ell in ells:
                _, L2, _ = get_layers(k, ell)

                actual = coefficient(L2, a, ell - 2)
                interior = P_interior(2, k, a, ell)

                values.append(simp(actual - interior))

            order = first_constant_difference_order(values)

            print(
                f"  a=k+{offset} "
                f"Delta={values} "
                f"constant-difference-order={order}"
            )

        print()


# ==============================================================================
# 7. L3 DEGREE AUDIT
# ==============================================================================

def L3_degree_audit():
    print("=" * 78)
    print("7. L3 BOUNDARY DEGREE IN ell")
    print("=" * 78)

    for k, ells in CASES.items():
        print(f"k={k}")

        for offset in [0, 1, 2, 3]:
            a = k + offset

            values = []

            for ell in ells:
                _, _, L3 = get_layers(k, ell)

                actual = coefficient(L3, a, ell - 3)
                interior = P_interior(3, k, a, ell)

                values.append(simp(actual - interior))

            order = first_constant_difference_order(values)

            print(
                f"  a=k+{offset} "
                f"Delta={values} "
                f"constant-difference-order={order}"
            )

        print()


# ==============================================================================
# 8. INTERPOLATE ONLY AFTER DIFFERENCE TEST
# ==============================================================================

def interpolate_boundary(k, offset, layer_number):
    ell_symbol = L
    a = k + offset

    data = []

    for ell in CASES[k]:
        _, L2, L3 = get_layers(k, ell)

        if layer_number == 2:
            layer = L2
            degree = ell - 2
        else:
            layer = L3
            degree = ell - 3

        actual = coefficient(layer, a, degree)
        interior = P_interior(layer_number, k, a, ell)

        delta = simp(actual - interior)

        data.append((sp.Integer(ell), delta))

    if len(data) < 3:
        return None

    poly = sp.interpolate(data, ell_symbol)

    return simp(poly)


# ==============================================================================
# 9. L2 FACTORIZATION
# ==============================================================================

def L2_factorization_audit():
    print("=" * 78)
    print("9. L2 BOUNDARY FACTORIZATION")
    print("=" * 78)

    for k in CASES:
        print(f"k={k}")

        for offset in [0, 1]:
            poly = interpolate_boundary(k, offset, 2)

            if poly is None:
                print(
                    f"  a=k+{offset}: "
                    "insufficient data"
                )
                continue

            print(
                f"  a=k+{offset}: "
                f"{sp.factor(poly)}"
            )

        print()


# ==============================================================================
# 10. L3 FACTORIZATION
# ==============================================================================

def L3_factorization_audit():
    print("=" * 78)
    print("10. L3 BOUNDARY FACTORIZATION")
    print("=" * 78)

    for k in CASES:
        print(f"k={k}")

        for offset in [0, 1, 2, 3]:
            poly = interpolate_boundary(k, offset, 3)

            if poly is None:
                print(
                    f"  a=k+{offset}: "
                    "insufficient data"
                )
                continue

            print(
                f"  a=k+{offset}: "
                f"{sp.factor(poly)}"
            )

        print()


# ==============================================================================
# 11. CROSS-k L3 a=k+2 TEST
# ==============================================================================

def L3_cross_k_a_k2():
    print("=" * 78)
    print("11. CROSS-k L3 a=k+2 STRUCTURE")
    print("=" * 78)

    print(
        "For each k, Delta is computed directly from the exact kernel."
    )
    print()

    for k in CASES:
        values = []

        for ell in CASES[k]:
            _, _, L3 = get_layers(k, ell)

            a = k + 2

            actual = coefficient(L3, a, ell - 3)
            interior = P_interior(3, k, a, ell)

            delta = simp(actual - interior)

            values.append((ell, delta))

        print(f"k={k}:")
        print(f"  {values}")

        if len(values) >= 3:
            poly = interpolate_boundary(k, 2, 3)
            print(f"  interpolated Delta(L) = {sp.factor(poly)}")

        print()


# ==============================================================================
# 12. DIRECT TEST OF THE LOW-k a=k+2 CANDIDATE
# ==============================================================================

def low_k_candidate_a_k2():
    print("=" * 78)
    print("12. LOW-k a=k+2 CANDIDATE TEST")
    print("=" * 78)

    failures = 0

    candidate = (
        (L - K - 4)
        * (
            (2 * K - 1) * L
            - K * (K + 3)
        )
        / 2
    )

    print("candidate:")
    print(f"  {sp.factor(candidate)}")
    print()

    for k, ells in CASES.items():
        for ell in ells:
            _, _, L3 = get_layers(k, ell)

            actual = coefficient(
                L3,
                k + 2,
                ell - 3,
            )

            interior = P_interior(
                3,
                k,
                k + 2,
                ell,
            )

            delta = simp(actual - interior)

            expected = candidate.subs(
                {
                    K: k,
                    L: ell,
                }
            )

            residual = simp(delta - expected)

            ok = residual == 0

            print(
                f"k={k:2d} ell={ell:2d} "
                f"actual={str(delta):>8} "
                f"candidate={str(expected):>8} "
                f"residual={str(residual):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(f"candidate failures = {failures}")
    print()


# ==============================================================================
# 13. SYMBOLIC INTERIOR LAW
# ==============================================================================

def symbolic_interior_audit():
    print("=" * 78)
    print("13. SYMBOLIC INTERIOR PRODUCT LAW")
    print("=" * 78)

    for r in [1, 2, 3]:
        expected = (
            (-1) ** (r + 1)
            / sp.factorial(r)
            * sp.binomial(K + r, A)
        )

        product_part = sp.Integer(1)

        for j in range(1, r):
            product_part *= L - A - j

        expected *= product_part
        expected *= ((K + r) * L - K * A) / (K + r)

        expected = simp(expected)

        print(f"r={r}")
        print(f"  P_r = {expected}")
        print()


# ==============================================================================
# 14. FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():
    print("=" * 78)
    print("14. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
This experiment deliberately does NOT assume the earlier
L2 or L3 boundary formulas.

Every boundary correction is computed as:

  Delta_r(k,a,L)
    = exact kernel coefficient
      - interior product-law extrapolation.

The research questions are now:

  1. What is the exact L2 boundary ladder?
       a=k
       a=k+1
       a=k+2
       ...

  2. What is the exact L3 boundary ladder?
       a=k
       a=k+1
       a=k+2
       a=k+3
       ...

  3. What degree in ell does Delta_r(k,k+j,ell)
     actually have?

  4. Does the degree equal j?

  5. Do the factorizations reveal a general boundary
     product analogous to the interior product law?

The experiment never upgrades an interpolation to a theorem
without displaying the exact finite-difference evidence.

The most important output sections are:

  3. EXACT L2 BOUNDARY CORRECTIONS
  4. EXACT L3 BOUNDARY CORRECTIONS
  6. L2 DEGREE AUDIT
  7. L3 DEGREE AUDIT
  9. L2 BOUNDARY FACTORIZATION
  10. L3 BOUNDARY FACTORIZATION
  11. CROSS-k L3 a=k+2 STRUCTURE
  12. LOW-k a=k+2 CANDIDATE TEST

No L4 analysis is performed.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 241")
    print("EXACT KERNEL DERIVATION OF L2/L3 BOUNDARY CORRECTIONS")
    print("=" * 78)
    print()
    print("Standalone main.py")
    print("Exact arithmetic over QQ")
    print("No previous experiment imported")
    print("No filesystem access")
    print("No L4 analysis")
    print()

    sanity_audit()
    interior_control()
    extract_L2_boundary()
    extract_L3_boundary()
    L2_degree_audit()
    L3_degree_audit()
    L2_factorization_audit()
    L3_factorization_audit()
    L3_cross_k_a_k2()
    low_k_candidate_a_k2()
    symbolic_interior_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()

