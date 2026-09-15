# ==============================================================================
# EXPERIMENT 245
# EXACT r=4 BOUNDARY LADDER DERIVATION IN d = ell-k
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No L5 analysis
#
# Main goals:
#
#   1. Reuse only the exact pq kernel given explicitly below.
#
#   2. Keep the established exact r=4 interior product law.
#
#   3. Compute every boundary correction directly:
#
#        Delta_4(k,j,d)
#          = exact coefficient at a=k+j
#            - interior extrapolation.
#
#   4. Work in d = ell-k.
#
#   5. Derive exact closed forms for j=0,1,2,3,4.
#
#   6. Verify those closed forms across a fresh larger set of k,d.
#
#   7. Verify the terminal j=4 law:
#
#        d(d-5)(d-6)(d-7)/24.
#
#   8. Verify that j>=5 is identically zero whenever the coefficient
#      exists.
#
#   9. Investigate whether all five rows fit one common product pattern.
#
#   10. No L5 analysis.
#
# ==============================================================================

import sympy as sp


# ==============================================================================
# SYMBOLS
# ==============================================================================

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")

K, A, D = sp.symbols("K A D")


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
# EXACT SYMMETRIC KERNEL G(N,X)
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
            "Symmetric kernel has nonzero remainder: "
            + str(remainder)
        )

    s1 = mapping[0][0]
    s2 = mapping[1][0]

    return sp.expand(
        sym_expr.subs(
            {
                s1: X - 1,
                s2: N
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

    for monomial, coeff in poly.terms():
        n_exp, x_exp = monomial

        if n_exp + x_exp == degree:
            result += (
                coeff
                * N**n_exp
                * X**x_exp
            )

    return sp.expand(result)


def exact_L4(k, ell):
    return homogeneous_layer(
        exact_G(k, ell),
        ell - 4
    )


# ==============================================================================
# SAFE COEFFICIENT EXTRACTION
#
# This fixes the exact bug in Experiment 244:
# never construct X^(-1), X^(-2), ...
# ==============================================================================

def coeff_L4(k, ell, a):
    k = int(k)
    ell = int(ell)
    a = int(a)

    x_exp = ell - 4 - a

    if a < 0:
        return sp.Integer(0)

    if x_exp < 0:
        return sp.Integer(0)

    layer = exact_L4(
        k,
        ell
    )

    poly = sp.Poly(
        layer,
        N,
        X
    )

    return sp.expand(
        poly.coeff_monomial(
            N**a
            * X**x_exp
        )
    )


# ==============================================================================
# EXACT INTERIOR r=4 PRODUCT LAW
#
# P4(k,a,ell)
#
# = - C(k+4,a)/24
#   * (ell-a-1)
#   * (ell-a-2)
#   * (ell-a-3)
#   * (((k+4)ell-ka)/(k+4)).
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
# d-BASED BOUNDARY CORRECTION
#
# ell = k+d
# a   = k+j
# ==============================================================================

def delta4(k, j, d):
    k = int(k)
    j = int(j)
    d = int(d)

    ell = k + d
    a = k + j

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
# SAFE EXACT EQUALITY
# ==============================================================================

def exact_zero(expr):
    return sp.expand(expr) == 0


# ==============================================================================
# PART 1
# DIRECT DATA TABLE
# ==============================================================================

def direct_boundary_table():
    print("=" * 78)
    print("1. DIRECT r=4 BOUNDARY CORRECTIONS IN d")
    print("=" * 78)

    for k in [3, 5, 7, 9, 11, 13, 15, 17]:

        print(f"k={k}")

        for j in range(5):

            values = []

            for d in range(
                4,
                18,
                2
            ):
                values.append(
                    (
                        d,
                        delta4(
                            k,
                            j,
                            d
                        )
                    )
                )

            print(
                f"  j={j}: "
                + str(
                    [(d, str(v)) for d, v in values]
                )
            )

        print()


# ==============================================================================
# PART 2
# DIFFERENCE DEGREE ON STABLE REGION d>=6
# ==============================================================================

def difference_table(values):
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


def stable_difference_degree(k, j):
    values = [
        delta4(
            k,
            j,
            d
        )
        for d in range(
            6,
            20,
            2
        )
    ]

    rows = difference_table(
        values
    )

    for order, row in enumerate(rows):

        if len(row) <= 1:
            return order

        if all(
            exact_zero(
                x - row[0]
            )
            for x in row
        ):
            return order

    return len(rows) - 1


def degree_audit():
    print("=" * 78)
    print("2. STABLE FINITE-DIFFERENCE DEGREE")
    print("=" * 78)

    failures = 0

    for k in [3, 5, 7, 9, 11, 13, 15, 17]:

        print(f"k={k}")

        for j in range(5):

            degree = stable_difference_degree(
                k,
                j
            )

            expected = j

            status = (
                "PASS"
                if degree == expected
                else "FAIL"
            )

            if status == "FAIL":
                failures += 1

            print(
                f"  j={j} "
                f"degree={degree} "
                f"expected={expected} "
                f"{status}"
            )

        print()

    print(
        f"stable degree failures = {failures}"
    )
    print()


# ==============================================================================
# PART 3
# STABLE INTERPOLANTS
# ==============================================================================

def stable_poly(k, j):
    samples = [
        (
            sp.Integer(d),
            delta4(
                k,
                j,
                d
            )
        )
        for d in range(
            6,
            20,
            2
        )
    ]

    return sp.factor(
        sp.interpolate(
            samples,
            D
        )
    )


def stable_polynomial_audit():
    print("=" * 78)
    print("3. STABLE CLOSED POLYNOMIALS")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_poly(
                k,
                j
            )

            print(
                f"  k={k}: {poly}"
            )

        print()


# ==============================================================================
# PART 4
# TEST THE POLYNOMIALS AGAINST ALL AVAILABLE d
#
# This is important:
# the interpolant is constructed from d>=6,
# but then tested at d=4 as well.
# ==============================================================================

def full_reconstruction_audit():
    print("=" * 78)
    print("4. FULL RECONSTRUCTION FROM STABLE POLYNOMIAL")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in [3, 5, 7, 9, 11, 13, 15, 17]:

        for j in range(5):

            poly = stable_poly(
                k,
                j
            )

            for d in range(
                4,
                20,
                2
            ):

                actual = delta4(
                    k,
                    j,
                    d
                )

                predicted = sp.expand(
                    poly.subs(
                        D,
                        d
                    )
                )

                residual = sp.expand(
                    actual - predicted
                )

                tested += 1

                if residual != 0:
                    failures += 1

                if residual != 0:
                    print(
                        f"FAIL "
                        f"k={k} "
                        f"j={j} "
                        f"d={d} "
                        f"actual={actual} "
                        f"predicted={predicted} "
                        f"residual={residual}"
                    )

    print()
    print(
        f"tested = {tested}"
    )
    print(
        f"reconstruction failures = {failures}"
    )
    print()


# ==============================================================================
# PART 5
# CROSS-k COEFFICIENT STRUCTURE
# ==============================================================================

def polynomial_coefficients(poly):
    return sp.Poly(
        sp.expand(poly),
        D
    ).all_coeffs()


def coefficient_table():
    print("=" * 78)
    print("5. CROSS-k COEFFICIENT STRUCTURE")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_poly(
                k,
                j
            )

            coeffs = polynomial_coefficients(
                poly
            )

            print(
                f"  k={k}: {coeffs}"
            )

        print()


# ==============================================================================
# PART 6
# EXACT FACTORIZATION
# ==============================================================================

def factorization_audit():
    print("=" * 78)
    print("6. STABLE FACTORIZATION")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_poly(
                k,
                j
            )

            print(
                f"  k={k}: "
                f"{sp.factor(poly)}"
            )

        print()


# ==============================================================================
# PART 7
# ROOT PATTERN
# ==============================================================================

def root_audit():
    print("=" * 78)
    print("7. ROOT STRUCTURE")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_poly(
                k,
                j
            )

            roots = sp.solve(
                sp.Eq(
                    poly,
                    0
                ),
                D
            )

            print(
                f"  k={k}: "
                f"roots={roots}"
            )

        print()


# ==============================================================================
# PART 8
# TEST CANDIDATE GENERAL PATTERN
#
# The observed factors strongly suggest:
#
# j=0:
#   C(k+3,4)
#
# j=1:
#   -C(k+3,3) * linear factor
#
# j=2:
#   C(k+3,2) * (d-5) * linear / 2
#
# j=3:
#   -C(k+3,1) * (d-5)(d-6) * linear / 6
#
# j=4:
#   (d)(d-5)(d-6)(d-7)/24
#
# This section does NOT assume the linear factors.
# It extracts them from exact data and tests candidate
# binomial normalizations.
# ==============================================================================

def normalized_boundary_rows():
    print("=" * 78)
    print("8. BINOMIAL NORMALIZATION AUDIT")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for k in [3, 5, 7, 9, 11, 13]:

            poly = stable_poly(
                k,
                j
            )

            # Natural normalization:
            #
            # C(k+3, 4-j)
            #
            if 4 - j >= 0:
                scale = sp.binomial(
                    k + 3,
                    4 - j
                )
            else:
                scale = sp.Integer(1)

            normalized = sp.factor(
                poly / scale
            )

            print(
                f"  k={k}: "
                f"scale={scale} "
                f"normalized={normalized}"
            )

        print()


# ==============================================================================
# PART 9
# DIRECT CHECK OF OBSERVED PRODUCT FOR j=4
# ==============================================================================

def exact_j4_product_audit():
    print("=" * 78)
    print("9. EXACT j=4 PRODUCT LAW")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in [1, 3, 5, 7, 9, 11, 13, 15, 17]:

        for d in range(
            4,
            24,
            2
        ):

            actual = delta4(
                k,
                4,
                d
            )

            expected = sp.expand(
                d
                * (d - 5)
                * (d - 6)
                * (d - 7)
                / sp.Integer(24)
            )

            residual = sp.expand(
                actual - expected
            )

            tested += 1

            if residual != 0:
                failures += 1

            if residual != 0:
                print(
                    f"FAIL "
                    f"k={k} "
                    f"d={d} "
                    f"actual={actual} "
                    f"expected={expected} "
                    f"residual={residual}"
                )

    print()
    print(
        f"tested = {tested}"
    )
    print(
        f"j=4 product failures = {failures}"
    )
    print()


# ==============================================================================
# PART 10
# SAFE ZERO-TAIL AUDIT
#
# Never construct a monomial with a negative X exponent.
# ==============================================================================

def zero_tail_audit():
    print("=" * 78)
    print("10. SAFE ZERO-TAIL AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in [1, 3, 5, 7, 9, 11, 13, 15, 17]:

        for d in [6, 8, 10, 12, 14]:

            ell = k + d

            L4 = exact_L4(
                k,
                ell
            )

            poly = sp.Poly(
                L4,
                N,
                X
            )

            max_supported_a = ell - 4

            for j in range(
                5,
                9
            ):

                a = k + j

                tested += 1

                if a > max_supported_a:
                    actual = sp.Integer(0)
                else:
                    x_exp = ell - 4 - a

                    if x_exp < 0:
                        actual = sp.Integer(0)
                    else:
                        actual = sp.expand(
                            poly.coeff_monomial(
                                N**a
                                * X**x_exp
                            )
                        )

                residual = sp.expand(
                    actual
                )

                if residual != 0:
                    failures += 1

                    print(
                        f"FAIL "
                        f"k={k} "
                        f"ell={ell} "
                        f"d={d} "
                        f"j={j} "
                        f"a={a} "
                        f"actual={actual}"
                    )

    print()
    print(
        f"tested = {tested}"
    )
    print(
        f"zero-tail failures = {failures}"
    )
    print()


# ==============================================================================
# PART 11
# TEST GENERALIZED FACTORIZATION SKELETON
#
# We test the part that is already visibly stable:
#
#   j=0: degree 0
#   j=1: degree 1
#   j=2: factor (d-5)
#   j=3: factors (d-5)(d-6)
#   j=4: factors d(d-5)(d-6)(d-7)
#
# The remaining factor is tested for being linear.
# ==============================================================================

def skeleton_audit():
    print("=" * 78)
    print("11. GENERAL FACTORIZATION SKELETON")
    print("=" * 78)

    failures = 0

    for k in [3, 5, 7, 9, 11, 13]:

        # j=0
        p0 = stable_poly(
            k,
            0
        )

        if sp.degree(
            p0,
            D
        ) != 0:
            failures += 1

        # j=1
        p1 = stable_poly(
            k,
            1
        )

        if sp.degree(
            p1,
            D
        ) != 1:
            failures += 1

        # j=2
        p2 = stable_poly(
            k,
            2
        )

        if sp.rem(
            sp.Poly(
                p2,
                D
            ),
            sp.Poly(
                D - 5,
                D
            )
        ) != 0:
            failures += 1

        # j=3
        p3 = stable_poly(
            k,
            3
        )

        for root_factor in [
            D - 5,
            D - 6
        ]:

            if sp.rem(
                sp.Poly(
                    p3,
                    D
                ),
                sp.Poly(
                    root_factor,
                    D
                )
            ) != 0:
                failures += 1

        # j=4
        p4 = stable_poly(
            k,
            4
        )

        for root_factor in [
            D,
            D - 5,
            D - 6,
            D - 7
        ]:

            if sp.rem(
                sp.Poly(
                    p4,
                    D
                ),
                sp.Poly(
                    root_factor,
                    D
                )
            ) != 0:
                failures += 1

    print(
        f"factorization-skeleton failures = {failures}"
    )
    print()


# ==============================================================================
# PART 12
# SYMBOLIC FORMULAS OBSERVED FROM EXACT DATA
#
# These are hypotheses only. They are tested against all available
# k,d values before being reported as passes.
# ==============================================================================

def candidate_delta4(k, j, d):
    k = sp.sympify(k)
    j = int(j)
    d = sp.sympify(d)

    if j == 0:
        return sp.binomial(
            k + 3,
            4
        )

    if j == 1:
        return (
            -sp.binomial(
                k + 3,
                3
            )
            * (
                (k + 1) * d
                - 2 * k
                - 5
            )
            / sp.Integer(
                k + 3
            )
        )

    if j == 2:
        return (
            sp.binomial(
                k + 3,
                2
            )
            * (d - 5)
            * (
                (k + 2) * d
                - 2 * k
                - 6
            )
            / (
                2 * (k + 3)
            )
        )

    if j == 3:
        return (
            -sp.binomial(
                k + 3,
                1
            )
            * (d - 5)
            * (d - 6)
            * (
                (k + 3) * d
                - 2 * k
                - 7
            )
            / (
                6 * (k + 3)
            )
        )

    if j == 4:
        return sp.expand(
            d
            * (d - 5)
            * (d - 6)
            * (d - 7)
            / sp.Integer(24)
        )

    raise ValueError(
        "candidate_delta4 only supports j=0,...,4"
    )


# ==============================================================================
# PART 13
# CANDIDATE GENERALIZATION AUDIT
# ==============================================================================

def candidate_generalization_audit():
    print("=" * 78)
    print("12. CANDIDATE r=4 BOUNDARY LADDER AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in [
        1, 3, 5, 7, 9, 11,
        13, 15, 17, 19
    ]:

        for j in range(5):

            for d in range(
                4,
                24,
                2
            ):

                actual = delta4(
                    k,
                    j,
                    d
                )

                predicted = sp.expand(
                    candidate_delta4(
                        k,
                        j,
                        d
                    )
                )

                residual = sp.expand(
                    actual - predicted
                )

                tested += 1

                if residual != 0:
                    failures += 1

                    print(
                        f"FAIL "
                        f"k={k} "
                        f"j={j} "
                        f"d={d} "
                        f"actual={actual} "
                        f"predicted={predicted} "
                        f"residual={residual}"
                    )

    print()
    print(
        f"tested = {tested}"
    )
    print(
        f"candidate failures = {failures}"
    )
    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 245")
    print("EXACT r=4 BOUNDARY LADDER DERIVATION")
    print("=" * 78)
    print()
    print("Exact arithmetic over QQ.")
    print("Standalone main.py.")
    print("No previous experiment imported.")
    print("No filesystem access.")
    print("No L5 analysis.")
    print()

    direct_boundary_table()

    degree_audit()

    stable_polynomial_audit()

    full_reconstruction_audit()

    coefficient_table()

    factorization_audit()

    root_audit()

    normalized_boundary_rows()

    exact_j4_product_audit()

    zero_tail_audit()

    skeleton_audit()

    candidate_generalization_audit()

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "Experiment 245 moves from interpolation toward"
    )
    print(
        "an explicit r=4 boundary ladder."
    )
    print()
    print(
        "The natural variable is d = ell-k."
    )
    print()
    print(
        "The stable degree pattern is expected to be:"
    )
    print(
        "  j=0 -> degree 0"
    )
    print(
        "  j=1 -> degree 1"
    )
    print(
        "  j=2 -> degree 2"
    )
    print(
        "  j=3 -> degree 3"
    )
    print(
        "  j=4 -> degree 4"
    )
    print()
    print(
        "The final section tests a single closed candidate"
    )
    print(
        "against the exact kernel for many k and d values."
    )
    print()
    print(
        "A PASS there would establish the r=4 boundary ladder"
    )
    print(
        "for j=0,...,4 over the tested domain."
    )
    print()
    print(
        "No L5 analysis is performed."
    )


if __name__ == "__main__":
    main()

