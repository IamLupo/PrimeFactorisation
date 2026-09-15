import sympy as sp


# ==============================================================================
# EXPERIMENT 244
# EXACT L4 DOMAIN OF VALIDITY + d-BASED BOUNDARY LADDER
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No floating-point arithmetic
# No L5 analysis
#
# Goals:
#
#   1. Identify the exact domain of validity of the general interior
#      product law for r=4.
#
#   2. Do NOT treat the exceptional smallest ell-values as representative
#      of the general boundary polynomial.
#
#   3. Re-express all L4 boundary corrections using
#
#           d = ell - k
#
#      which is the natural variable visible in the data.
#
#   4. Determine exact finite-difference degree in d for
#
#           Delta_4(k, k+j, ell)
#
#      for j=0,1,2,3,4.
#
#   5. Search for exact cross-k formulas for the coefficients.
#
#   6. Independently verify the terminal j=4 product:
#
#         Delta_4(k,k+4,L)
#           = (L-k)(L-k-5)(L-k-6)(L-k-7)/24
#
#   7. Separate:
#
#        ordinary boundary ladder
#
#      from
#
#        minimal-ell exceptional points.
#
# ==============================================================================


p, q = sp.symbols("p q")
N, X = sp.symbols("N X")

K, A, L, D = sp.symbols("K A L D")


# ==============================================================================
# EXACT pq KERNEL
# ==============================================================================

def exact_F(k, ell):
    k = sp.Integer(k)
    ell = sp.Integer(ell)

    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ==============================================================================
# EXACT G(N,X)
# ==============================================================================

def exact_G(k, ell):
    F = exact_F(k, ell)

    sym_expr, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True
    )

    if sp.expand(remainder) != 0:
        raise AssertionError(
            f"Nonzero symmetrization remainder: {remainder}"
        )

    s1 = mapping[0][0]
    s2 = mapping[1][0]

    return sp.expand(
        sym_expr.subs(
            {
                s1: X - 1,
                s2: N,
            }
        )
    )


# ==============================================================================
# HOMOGENEOUS LAYER
# ==============================================================================

def homogeneous_layer(G, degree):
    poly = sp.Poly(
        sp.expand(G),
        N,
        X
    )

    result = sp.Integer(0)

    for (n_exp, x_exp), coeff in poly.terms():
        if n_exp + x_exp == degree:
            result += coeff * N**n_exp * X**x_exp

    return sp.expand(result)


def exact_L4(k, ell):
    G = exact_G(k, ell)
    return homogeneous_layer(
        G,
        ell - 4
    )


# ==============================================================================
# EXACT COEFFICIENT
# ==============================================================================

def coeff_L4(k, ell, a):
    layer = exact_L4(k, ell)

    x_exp = ell - 4 - a

    if x_exp < 0:
        return sp.Integer(0)

    return sp.expand(
        sp.Poly(
            layer,
            N,
            X
        ).coeff_monomial(
            N**a * X**x_exp
        )
    )


# ==============================================================================
# GENERAL INTERIOR PRODUCT LAW
#
# P_r =
# (-1)^(r+1)/r!
# * C(k+r,a)
# * product_{j=1}^{r-1}(ell-a-j)
# * ((k+r)ell-ka)/(k+r)
# ==============================================================================

def interior_P4(k, a, ell):
    k = sp.sympify(k)
    a = sp.sympify(a)
    ell = sp.sympify(ell)

    return sp.expand(
        -sp.binomial(k + 4, a)
        * (ell - a - 1)
        * (ell - a - 2)
        * (ell - a - 3)
        * (
            ((k + 4) * ell - k * a)
            / (k + 4)
        )
        / sp.Integer(24)
    )


# ==============================================================================
# EXACT EQUALITY
# ==============================================================================

def exact_equal(x, y):
    return sp.expand(x - y) == 0


# ==============================================================================
# PART 1
# DOMAIN AUDIT
#
# Search over small k and ell and find precisely where the r=4 interior
# law works for every a<k.
# ==============================================================================

def interior_domain_audit():
    print("=" * 78)
    print("1. EXACT L4 INTERIOR DOMAIN OF VALIDITY")
    print("=" * 78)

    failures = []

    for k in range(1, 14, 2):
        print(f"k={k}")

        for ell in range(k + 4, k + 17, 2):

            local_failures = []

            for a in range(k):
                actual = coeff_L4(
                    k,
                    ell,
                    a
                )

                expected = interior_P4(
                    k,
                    a,
                    ell
                )

                if not exact_equal(
                    actual,
                    expected
                ):
                    local_failures.append(
                        (
                            a,
                            actual,
                            expected,
                            sp.expand(actual - expected)
                        )
                    )

            status = (
                "PASS"
                if not local_failures
                else "FAIL"
            )

            print(
                f"  ell={ell:2d} "
                f"d={ell-k:2d} "
                f"{status}"
            )

            if local_failures:
                failures.append(
                    (k, ell, local_failures)
                )

                for (
                    a,
                    actual,
                    expected,
                    residual
                ) in local_failures:
                    print(
                        f"    a={a} "
                        f"actual={actual} "
                        f"expected={expected} "
                        f"residual={residual}"
                    )

        print()

    print(
        f"interior-domain failure cases = {len(failures)}"
    )
    print()


# ==============================================================================
# PART 2
# BOUNDARY CORRECTION
#
# Delta_4(k,j,L)
#   = exact coefficient at a=k+j
#     - interior extrapolation
#
# We immediately rewrite ell = k+d.
# ==============================================================================

def delta_L4_d(k, j, d):
    ell = sp.Integer(k) + sp.Integer(d)
    a = sp.Integer(k) + sp.Integer(j)

    actual = coeff_L4(
        k,
        ell,
        a
    )

    interior = interior_P4(
        k,
        a,
        ell
    )

    return sp.expand(
        actual - interior
    )


# ==============================================================================
# PART 3
# PRINT RAW d-SEQUENCES
#
# We deliberately begin at d=4.
# Later we also inspect d>=6 separately.
# ==============================================================================

def raw_boundary_sequences():
    print("=" * 78)
    print("2. RAW L4 BOUNDARY SEQUENCES IN d = ell-k")
    print("=" * 78)

    k_values = [3, 5, 7, 9, 11, 13]

    for k in k_values:
        print(f"k={k}")

        for j in range(5):

            rows = []

            for d in range(4, 16, 2):
                delta = delta_L4_d(
                    k,
                    j,
                    d
                )

                rows.append(
                    (
                        d,
                        delta
                    )
                )

            print(
                f"  j={j}: "
                + str(
                    [(d, str(v)) for d, v in rows]
                )
            )

        print()


# ==============================================================================
# PART 4
# FINITE DIFFERENCE IN d
# ==============================================================================

def finite_difference_rows(values):
    rows = [list(values)]
    current = list(values)

    while len(current) > 1:
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


def first_constant_difference_order(values):
    rows = finite_difference_rows(values)

    for order, row in enumerate(rows):

        if len(row) <= 1:
            return order

        if all(
            exact_equal(v, row[0])
            for v in row
        ):
            return order

    return len(rows) - 1


def degree_in_d_audit():
    print("=" * 78)
    print("3. EXACT FINITE-DIFFERENCE DEGREE IN d")
    print("=" * 78)

    k_values = [3, 5, 7, 9, 11, 13]

    for k in k_values:
        print(f"k={k}")

        for j in range(5):

            values = [
                delta_L4_d(
                    k,
                    j,
                    d
                )
                for d in range(
                    4,
                    16,
                    2
                )
            ]

            rows = finite_difference_rows(
                values
            )

            print(
                f"  j={j}"
            )

            for order, row in enumerate(rows):
                print(
                    f"    Delta^{order} = "
                    + str(
                        row
                    )
                )

                if len(row) > 1 and all(
                    exact_equal(
                        v,
                        row[0]
                    )
                    for v in row
                ):
                    print(
                        f"    first constant difference order = "
                        f"{order}"
                    )
                    break

        print()


# ==============================================================================
# PART 5
# IMPORTANT: REMOVE THE MINIMAL-d TRANSIENT
#
# For j >= 0 the first few d-values can touch the support boundary.
# We therefore repeat the finite-difference audit using only
#
#       d >= 6
#
# and compare it to the d>=4 sequence.
# ==============================================================================

def stable_degree_audit():
    print("=" * 78)
    print("4. STABLE DEGREE AUDIT USING d >= 6")
    print("=" * 78)

    k_values = [3, 5, 7, 9, 11, 13]

    for k in k_values:
        print(f"k={k}")

        for j in range(5):

            values = [
                delta_L4_d(
                    k,
                    j,
                    d
                )
                for d in range(
                    6,
                    18,
                    2
                )
            ]

            order = first_constant_difference_order(
                values
            )

            print(
                f"  j={j} "
                f"values={list(map(str, values))} "
                f"degree={order}"
            )

        print()


# ==============================================================================
# PART 6
# INTERPOLATE IN d ONLY AFTER THE DIFFERENCE TEST
#
# We fit using the stable d>=6 values.
# ==============================================================================

def stable_interpolant(k, j):
    samples = []

    for d in range(
        6,
        18,
        2
    ):
        samples.append(
            (
                sp.Integer(d),
                delta_L4_d(
                    k,
                    j,
                    d
                )
            )
        )

    return sp.factor(
        sp.interpolate(
            samples,
            D
        )
    )


def stable_interpolation_audit():
    print("=" * 78)
    print("5. STABLE d-POLYNOMIALS")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_interpolant(
                k,
                j
            )

            print(
                f"  k={k}: "
                f"{poly}"
            )

        print()


# ==============================================================================
# PART 7
# LEADING COEFFICIENTS OF THE STABLE d-POLYNOMIAL
# ==============================================================================

def stable_coefficient_audit():
    print("=" * 78)
    print("6. STABLE d-POLYNOMIAL COEFFICIENTS")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = sp.Poly(
                sp.expand(
                    stable_interpolant(
                        k,
                        j
                    )
                ),
                D
            )

            coeffs = poly.all_coeffs()

            print(
                f"  k={k}: "
                + str(
                    coeffs
                )
            )

        print()


# ==============================================================================
# PART 8
# FACTORIZATION OF STABLE d-POLYNOMIALS
# ==============================================================================

def stable_factorization_audit():
    print("=" * 78)
    print("7. STABLE BOUNDARY FACTORIZATION IN d")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_interpolant(
                k,
                j
            )

            print(
                f"  k={k}: "
                f"{sp.factor(poly)}"
            )

        print()


# ==============================================================================
# PART 9
# CROSS-k COMPARISON
#
# Check whether the stable polynomials exhibit a simple dependence
# on k, especially in their roots and leading coefficients.
# ==============================================================================

def cross_k_structure():
    print("=" * 78)
    print("8. CROSS-k STRUCTURE OF STABLE BOUNDARY POLYNOMIALS")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_interpolant(
                k,
                j
            )

            roots = sp.solve(
                sp.Eq(poly, 0),
                D
            )

            roots = [
                sp.factor(root)
                for root in roots
            ]

            print(
                f"  k={k}: "
                f"roots={roots}"
            )

        print()


# ==============================================================================
# PART 10
# TERMINAL J=4 EXACT FORMULA AUDIT
# ==============================================================================

def terminal_j4_audit():
    print("=" * 78)
    print("9. EXACT TERMINAL j=4 PRODUCT AUDIT")
    print("=" * 78)

    failures = 0

    for k in [3, 5, 7, 9, 11, 13, 15, 17]:

        for d in range(
            4,
            16,
            2
        ):

            actual = delta_L4_d(
                k,
                4,
                d
            )

            expected = sp.expand(
                (d - 5)
                * (d - 6)
                * (d - 7)
                * d
                / sp.Integer(24)
            )

            residual = sp.expand(
                actual - expected
            )

            if residual != 0:
                failures += 1

            print(
                f"k={k:2d} "
                f"d={d:2d} "
                f"actual={actual} "
                f"expected={expected} "
                f"residual={residual} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

    print()
    print(
        f"terminal j=4 failures = {failures}"
    )
    print()


# ==============================================================================
# PART 11
# TERMINAL CANCELLATION
#
# For j >= 5 the coefficient should be zero in the observed support.
# ==============================================================================

def terminal_zero_audit():
    print("=" * 78)
    print("10. TERMINAL ZERO-TAIL AUDIT")
    print("=" * 78)

    failures = 0

    for k in [3, 5, 7, 9, 11, 13]:

        ell = k + 10

        L4 = exact_L4(
            k,
            ell
        )

        for j in range(
            5,
            9
        ):

            a = k + j

            actual = sp.expand(
                sp.Poly(
                    L4,
                    N,
                    X
                ).coeff_monomial(
                    N**a
                    * X**(ell - 4 - a)
                )
            )

            expected = sp.Integer(0)

            residual = sp.expand(
                actual - expected
            )

            if residual != 0:
                failures += 1

            print(
                f"k={k:2d} "
                f"ell={ell:2d} "
                f"j={j} "
                f"actual={actual} "
                f"residual={residual} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

    print()
    print(
        f"terminal zero-tail failures = {failures}"
    )
    print()


# ==============================================================================
# PART 12
# CHECK r=4 INTERIOR LAW AFTER THE MINIMAL EDGE
#
# The k=3, ell=7 failures should disappear once ell is separated
# from the minimal layer degree.
# ==============================================================================

def edge_exception_audit():
    print("=" * 78)
    print("11. MINIMAL-ELL EDGE EXCEPTION AUDIT")
    print("=" * 78)

    for k in [3, 5, 7, 9]:

        print(
            f"k={k}"
        )

        for ell in range(
            k + 4,
            k + 13,
            2
        ):

            failures = []

            for a in range(k):

                actual = coeff_L4(
                    k,
                    ell,
                    a
                )

                expected = interior_P4(
                    k,
                    a,
                    ell
                )

                if not exact_equal(
                    actual,
                    expected
                ):
                    failures.append(
                        (
                            a,
                            sp.expand(
                                actual - expected
                            )
                        )
                    )

            print(
                f"  ell={ell:2d} "
                f"d={ell-k:2d} "
                f"failures={failures}"
            )

        print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 244")
    print("EXACT L4 DOMAIN + d-BASED BOUNDARY LADDER")
    print("=" * 78)
    print()
    print("Exact arithmetic over QQ.")
    print("Standalone main.py.")
    print("No previous experiment imported.")
    print("No filesystem access.")
    print("No L5 analysis.")
    print()

    interior_domain_audit()

    raw_boundary_sequences()

    degree_in_d_audit()

    stable_degree_audit()

    stable_interpolation_audit()

    stable_coefficient_audit()

    stable_factorization_audit()

    cross_k_structure()

    terminal_j4_audit()

    terminal_zero_audit()

    edge_exception_audit()

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "The purpose of Experiment 244 is to separate"
    )
    print(
        "genuine boundary structure from minimal-ell edge effects."
    )
    print()
    print(
        "The natural variable is d = ell-k."
    )
    print()
    print(
        "The experiment does NOT accept an interpolated polynomial"
    )
    print(
        "unless its finite-difference degree is independently visible."
    )
    print()
    print(
        "The terminal j=4 formula is tested directly for many"
    )
    print(
        "k and d values."
    )
    print()
    print(
        "Only after the L4 boundary ladder is structurally stable"
    )
    print(
        "should the project attempt a general Delta_r formula."
    )
    print()
    print(
        "No L5 analysis is performed."
    )


if __name__ == "__main__":
    main()

