import sympy as sp


# =============================================================================
# EXPERIMENT 243
# EXACT L4 INTERIOR + BOUNDARY LADDER DERIVATION
# =============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No floating-point arithmetic
#
# Purpose:
#   1. Reconstruct the exact G(N,X) polynomial directly from the pq kernel.
#   2. Extract homogeneous layers L0,L1,L2,L3,L4.
#   3. Verify the established general interior product law through r=4.
#   4. Extract exact boundary corrections for r=2,3,4.
#   5. Test whether boundary correction degree in ell equals j.
#   6. Factor the exact boundary corrections.
#   7. Look for a genuine general boundary-ladder mechanism.
#
# No proposed L4 boundary formula is assumed to be true.
# =============================================================================


p, q = sp.symbols("p q")
N, X = sp.symbols("N X")
K, A, L = sp.symbols("K A L")


# =============================================================================
# EXACT pq KERNEL
# =============================================================================

def exact_F(k, ell):
    k = sp.Integer(k)
    ell = sp.Integer(ell)

    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# =============================================================================
# EXACT CONVERSION TO G(N,X)
#
# We first symmetrize in p,q:
#
#   s1 = p+q
#   s2 = pq
#
# and then use
#
#   X = p+q+1
#   N = pq
#
# so s1 = X-1.
# =============================================================================

def exact_G(k, ell):
    F = exact_F(k, ell)

    sym_expr, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True
    )

    if sp.expand(remainder) != 0:
        raise AssertionError(
            f"Symmetrization remainder is nonzero: {remainder}"
        )

    # mapping is [(s1, p+q), (s2, pq)]
    s1 = mapping[0][0]
    s2 = mapping[1][0]

    G = sym_expr.subs({
        s1: X - 1,
        s2: N,
    })

    return sp.expand(G)


# =============================================================================
# HOMOGENEOUS LAYER EXTRACTION
#
# G has total degree ell.
#
# Layer r is the homogeneous piece of total degree ell-r.
#
# r=0 -> L0
# r=1 -> L1
# r=2 -> L2
# r=3 -> L3
# r=4 -> L4
# =============================================================================

def homogeneous_layer(G, degree):
    poly = sp.Poly(sp.expand(G), N, X)

    result = sp.Integer(0)

    for (n_exp, x_exp), coeff in poly.terms():
        if n_exp + x_exp == degree:
            result += coeff * N**n_exp * X**x_exp

    return sp.expand(result)


def exact_layers(k, ell, max_r=4):
    G = exact_G(k, ell)

    layers = []

    for r in range(max_r + 1):
        degree = ell - r
        layers.append(homogeneous_layer(G, degree))

    return G, layers


# =============================================================================
# COEFFICIENT EXTRACTION
#
# [N^a X^(ell-r-a)] L_r
# =============================================================================

def layer_coefficient(layer, a, ell, r):
    x_exp = ell - r - a

    if x_exp < 0:
        return sp.Integer(0)

    return sp.expand(
        sp.Poly(layer, N, X).coeff_monomial(
            N**a * X**x_exp
        )
    )


# =============================================================================
# ESTABLISHED INTERIOR PRODUCT LAW
#
# P_r(k,a,L)
# =
# (-1)^(r+1)/r!
# * C(k+r,a)
# * product_{j=1}^{r-1}(L-a-j)
# * ((k+r)L-ka)/(k+r)
#
# Validated experimentally for the interior range a < k.
# =============================================================================

def interior_P(r, k, a, ell):
    k = sp.sympify(k)
    a = sp.sympify(a)
    ell = sp.sympify(ell)

    product_part = sp.Integer(1)

    for j in range(1, r):
        product_part *= ell - a - j

    return sp.cancel(
        sp.Rational((-1) ** (r + 1), 1)
        * sp.binomial(k + r, a)
        * product_part
        * (
            ((k + r) * ell - k * a)
            / (k + r)
        )
        / sp.factorial(r)
    )


# =============================================================================
# SAFE EXACT EQUALITY
# =============================================================================

def exact_equal(a, b):
    return sp.expand(a - b) == 0


# =============================================================================
# NUMERICAL INTERIOR AUDIT FOR r=1..4
# =============================================================================

def interior_audit():
    print("=" * 78)
    print("1. EXACT INTERIOR L1/L2/L3/L4 PRODUCT-LAW AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 11),
        (5, 13),
        (5, 15),
        (5, 17),
        (7, 15),
        (7, 17),
        (7, 19),
        (9, 21),
        (9, 23),
        (11, 23),
        (11, 25),
        (13, 25),
        (13, 27),
        (15, 27),
        (17, 29),
    ]

    total = 0
    failures = 0

    for kv, ellv in cases:
        _, layers = exact_layers(kv, ellv, max_r=4)

        local = 0

        for r in range(1, 5):
            if ellv - r < 0:
                continue

            for a in range(kv):
                actual = layer_coefficient(
                    layers[r],
                    a,
                    ellv,
                    r
                )

                expected = interior_P(
                    r,
                    kv,
                    a,
                    ellv
                )

                total += 1
                local += 1

                if not exact_equal(actual, expected):
                    failures += 1

        status = "PASS" if local >= 0 and all(
            exact_equal(
                layer_coefficient(
                    exact_layers(kv, ellv, max_r=4)[1][r],
                    a,
                    ellv,
                    r
                ),
                interior_P(r, kv, a, ellv)
            )
            for r in range(1, 5)
            for a in range(kv)
            if ellv - r >= 0
        ) else "FAIL"

        print(
            f"k={kv:2d} ell={ellv:2d} "
            f"interior tested={local:3d} {status}"
        )

    print()
    print(f"tested interior coefficients = {total}")
    print(f"interior formula failures = {failures}")
    print()


# =============================================================================
# EXACT BOUNDARY CORRECTIONS
#
# Delta_r(k,j,L)
# =
# exact coefficient at a=k+j
# -
# interior extrapolation at a=k+j
#
# This is intentionally computed directly from the exact kernel.
# =============================================================================

def boundary_delta(r, k, ell, j):
    a = k + j

    G, layers = exact_layers(k, ell, max_r=r)

    actual = layer_coefficient(
        layers[r],
        a,
        ell,
        r
    )

    interior = interior_P(
        r,
        k,
        a,
        ell
    )

    return sp.expand(actual - interior)


# =============================================================================
# DISPLAY EXACT BOUNDARY DATA
# =============================================================================

def boundary_data(r, j, k_values):
    print("=" * 78)
    print(f"BOUNDARY DATA r={r}, j={j}, a=k+{j}")
    print("=" * 78)

    for kv in k_values:
        values = []

        # Use several admissible ell values.
        base = kv + 4
        for offset in range(0, 6):
            ellv = base + 2 * offset

            if ellv <= kv:
                continue

            delta = sp.factor(
                boundary_delta(r, kv, ellv, j)
            )

            values.append((ellv, delta))

        print(f"k={kv}")
        for ellv, delta in values:
            print(
                f"  ell={ellv:2d} "
                f"Delta={str(delta)}"
            )
        print()


# =============================================================================
# FINITE-DIFFERENCE DEGREE TEST
#
# Exact differences in ell with step h=2.
# =============================================================================

def finite_differences(values):
    current = list(values)
    rows = [current]

    while len(current) > 1:
        nxt = [
            sp.expand(current[i + 1] - current[i])
            for i in range(len(current) - 1)
        ]
        rows.append(nxt)
        current = nxt

    return rows


def constant_difference_order(values):
    rows = finite_differences(values)

    for order, row in enumerate(rows):
        if len(row) <= 1:
            return order

        if all(exact_equal(v, row[0]) for v in row):
            return order

    return len(rows) - 1


def degree_audit(r, j, k_values):
    print("=" * 78)
    print(f"FINITE-DIFFERENCE DEGREE AUDIT r={r}, j={j}")
    print("=" * 78)

    for kv in k_values:
        base = kv + 4

        ell_values = [
            base + 2 * offset
            for offset in range(6)
        ]

        deltas = [
            boundary_delta(r, kv, ellv, j)
            for ellv in ell_values
        ]

        order = constant_difference_order(deltas)

        print(
            f"k={kv:2d} "
            f"Delta={list(map(str, deltas))} "
            f"degree={order} "
            f"expected<={j}"
        )

    print()


# =============================================================================
# SYMBOLIC INTERPOLATION IN L
#
# We interpolate Delta_r(k,j,L) for selected numeric k.
# This is diagnostic only.
# =============================================================================

def interpolate_boundary_in_L(r, j, kv):
    ellv = [
        kv + 4 + 2*i
        for i in range(max(j + 3, 5))
    ]

    values = [
        boundary_delta(r, kv, ellv_i, j)
        for ellv_i in ellv
    ]

    poly = sp.interpolate(
        [(sp.Integer(ellv[i]), values[i])
         for i in range(len(ellv))],
        L
    )

    return sp.factor(sp.expand(poly))


def interpolation_audit(r, j, k_values):
    print("=" * 78)
    print(f"INTERPOLATION AUDIT r={r}, j={j}")
    print("=" * 78)

    for kv in k_values:
        poly = interpolate_boundary_in_L(r, j, kv)

        print(f"k={kv}")
        print(f"  Delta(L) = {poly}")
        print(f"  expanded = {sp.expand(poly)}")
        print(f"  factored = {sp.factor(poly)}")
        print()

    print()


# =============================================================================
# SYMBOLIC CROSS-k FIT OF THE COEFFICIENTS
#
# For each r,j, extract polynomial coefficients in L
# and test simple dependence on k.
#
# This does NOT declare a theorem.
# =============================================================================

def coefficient_table(r, j, k_values):
    print("=" * 78)
    print(f"CROSS-k COEFFICIENT TABLE r={r}, j={j}")
    print("=" * 78)

    rows = []

    for kv in k_values:
        poly = interpolate_boundary_in_L(r, j, kv)
        coeffs = sp.Poly(
            sp.expand(poly),
            L
        ).all_coeffs()

        rows.append((kv, coeffs))

    max_len = max(len(coeffs) for _, coeffs in rows)

    for kv, coeffs in rows:
        padded = [sp.Integer(0)] * (max_len - len(coeffs)) + coeffs

        print(
            f"k={kv:2d}: "
            + "  ".join(
                f"c{i}={str(v)}"
                for i, v in enumerate(padded)
            )
        )

    print()


# =============================================================================
# FACTORIZATION AUDIT
# =============================================================================

def factor_audit(r, j, k_values):
    print("=" * 78)
    print(f"BOUNDARY FACTORIZATION AUDIT r={r}, j={j}")
    print("=" * 78)

    for kv in k_values:
        poly = interpolate_boundary_in_L(r, j, kv)

        print(f"k={kv}")
        print(f"  polynomial = {sp.factor(poly)}")

        roots = sp.solve(
            sp.Eq(poly, 0),
            L
        )

        print(
            "  roots      = "
            + str([sp.factor(root) for root in roots])
        )

        print()

    print()


# =============================================================================
# SPECIAL ATTENTION: r=4
#
# The point of the experiment is to see whether the boundary ladder
# continues one more rung:
#
#   j=0 -> degree 0
#   j=1 -> degree 1
#   j=2 -> degree 2
#   j=3 -> degree 3
#   j=4 -> degree 4
#
# without assuming it beforehand.
# =============================================================================

def r4_boundary_audit():
    k_values = [
        3, 5, 7, 9, 11, 13
    ]

    for j in range(0, 5):
        boundary_data(
            4,
            j,
            k_values
        )

        degree_audit(
            4,
            j,
            k_values
        )

        interpolation_audit(
            4,
            j,
            k_values
        )

        coefficient_table(
            4,
            j,
            k_values
        )

        factor_audit(
            4,
            j,
            k_values
        )


# =============================================================================
# CROSS-CHECK r=2 AND r=3
# =============================================================================

def established_boundary_crosscheck():
    print("=" * 78)
    print("EXACT r=2 / r=3 BOUNDARY CROSS-CHECK")
    print("=" * 78)

    k_values = [3, 5, 7, 9, 11, 13]

    for r in (2, 3):
        for j in range(0, r + 1):
            print(
                f"r={r} j={j}"
            )

            for kv in k_values:
                ellv = kv + 9

                delta = boundary_delta(
                    r,
                    kv,
                    ellv,
                    j
                )

                print(
                    f"  k={kv:2d} "
                    f"ell={ellv:2d} "
                    f"Delta={delta}"
                )

            print()

    print()


# =============================================================================
# TERMINAL SUPPORT CHECK
#
# For r=4, inspect a=k+4 and a>k+4.
# This tells us whether the next terminal cancellation begins
# exactly where the natural boundary ladder predicts.
# =============================================================================

def terminal_support_audit():
    print("=" * 78)
    print("r=4 TERMINAL SUPPORT / CANCELLATION AUDIT")
    print("=" * 78)

    cases = [
        (3, 9),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
        (13, 25),
    ]

    for kv, ellv in cases:
        G, layers = exact_layers(
            kv,
            ellv,
            max_r=4
        )

        print(f"k={kv} ell={ellv}")

        for j in range(0, 7):
            a = kv + j

            actual = layer_coefficient(
                layers[4],
                a,
                ellv,
                4
            )

            interior = interior_P(
                4,
                kv,
                a,
                ellv
            )

            delta = sp.expand(
                actual - interior
            )

            print(
                f"  j={j} a={a:2d} "
                f"actual={actual} "
                f"interior={sp.factor(interior)} "
                f"Delta={sp.factor(delta)}"
            )

        print()


# =============================================================================
# SYMBOLIC PRODUCT LAW r=4
# =============================================================================

def symbolic_r4_product_law():
    print("=" * 78)
    print("SYMBOLIC r=4 INTERIOR PRODUCT LAW")
    print("=" * 78)

    r = 4

    candidate = interior_P(
        r,
        K,
        A,
        L
    )

    expected = sp.factor(
        -sp.binomial(K + 4, A)
        * (L - A - 1)
        * (L - A - 2)
        * (L - A - 3)
        * ((K + 4) * L - K * A)
        / (24 * (K + 4))
    )

    residual = sp.factor(
        sp.expand(candidate - expected)
    )

    print("candidate =")
    print(candidate)

    print()
    print("expected =")
    print(expected)

    print()
    print("residual =")
    print(residual)

    print(
        "PASS" if residual == 0 else "FAIL"
    )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 243")
    print("EXACT L4 INTERIOR + BOUNDARY LADDER DERIVATION")
    print("=" * 78)
    print()
    print("All arithmetic is exact over QQ.")
    print("No previous experiment imported.")
    print("No filesystem access.")
    print("No floating-point arithmetic.")
    print("No L5 analysis.")
    print()

    symbolic_r4_product_law()

    interior_audit()

    established_boundary_crosscheck()

    # Focus on r=4 boundary ladder.
    r4_boundary_audit()

    terminal_support_audit()

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "This experiment derives L4 directly from the exact pq kernel."
    )
    print()
    print(
        "The main unresolved question is whether the boundary correction"
    )
    print(
        "Delta_4(k,k+j,L) has a genuine degree-j ladder for j=0,...,4."
    )
    print()
    print(
        "No boundary formula is accepted from interpolation alone."
    )
    print(
        "Every reported Delta is obtained as:"
    )
    print(
        "  exact kernel coefficient - interior extrapolation"
    )
    print()
    print(
        "If the r=4 ladder closes, the next step is to compare the"
    )
    print(
        "r=2,3,4 boundary products and search for a general formula"
    )
    print(
        "Delta_r(k,k+j,L)."
    )
    print()
    print(
        "Do not move to L5 until the r=4 boundary mechanism is understood."
    )
    print()


if __name__ == "__main__":
    main()

