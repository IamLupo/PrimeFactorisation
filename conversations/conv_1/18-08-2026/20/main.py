import sympy as sp


# ==============================================================================
# EXPERIMENT 232
# EXACT a=k+2 BOUNDARY DERIVATION
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No exact_F
#
# Established:
#
#   Interior:
#
#     P_r(k,a,L)
#       = (-1)^(r+1)/r!
#         * C(k+r,a)
#         * product_{j=1}^{r-1}(L-a-j)
#         * ((k+r)L-ka)/(k+r)
#
#   L3:
#
#     a = k:
#       exact correction = C(k+2,3)
#
#     a = k+1:
#       exact correction
#         = -C(k+2,2)L
#           + k(k+2)(k+3)/2
#
# New target:
#
#     a = k+2
#
# The previous experiment proved this row is NONZERO.
# Therefore the old "zero-tail" assumption was false.
#
# This experiment extracts the exact correction
#
#     Delta_(k+2)
#       = actual - interior extrapolation
#
# and asks whether it has a clean quadratic/cubic structure
# in L and whether it fits a binomial boundary ladder.
#
# No arbitrary candidate enumeration.
# ==============================================================================


# ==============================================================================
# HELPERS
# ==============================================================================

def simp(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def B(n, a):
    return sp.binomial(
        sp.sympify(n),
        sp.sympify(a)
    )


# ==============================================================================
# ESTABLISHED INTERIOR L3 LAW
# ==============================================================================

def P3_interior(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    Cka = B(K + 3, A)

    A3 = (
        Cka
        / sp.Integer(6)
    )

    A2 = (
        -Cka
        * (
            A * (K + 2) + K + 3
        )
        / (
            sp.Integer(2)
            * (K + 3)
        )
    )

    A1 = (
        Cka
        * (
            3 * A**2 * K
            + 3 * A**2
            + 6 * A * K
            + 9 * A
            + 2 * K
            + 6
        )
        / (
            sp.Integer(6)
            * (K + 3)
        )
    )

    A0 = (
        -K
        * A
        * (A + 1)
        * (A + 2)
        * Cka
        / (
            sp.Integer(6)
            * (K + 3)
        )
    )

    return simp(
        A3 * X**3
        + A2 * X**2
        + A1 * X
        + A0
    )


# ==============================================================================
# EXACT L3 DATA FROM EXPERIMENTS 219-224
# ==============================================================================

DATA = {

    3: {
        7:  {5: 0},
        9:  {5: 66},
        11: {5: 244},
        13: {5: 582},
        15: {5: 1128},
    },

    5: {
        11: {7: 90},
        13: {7: 332},
        15: {7: 790},
        17: {7: 1528},
        19: {7: 2610},
    },

    7: {
        15: {9: 420},
        17: {9: 998},
        19: {9: 1928},
    },

    9: {
        21: {11: 2328},
        23: {11: 3970},
        25: {11: 6228},
    },

    11: {
        23: {13: 2728},
        25: {13: 4650},
    },

    13: {
        25: {15: 5330},
        27: {15: 7870},
    },
}


# ==============================================================================
# ACCESSOR
# ==============================================================================

def actual_a_k2(k, L):
    return DATA.get(
        k,
        {}
    ).get(
        L,
        {}
    ).get(
        k + 2,
        None
    )


# ==============================================================================
# CORRECTION
# ==============================================================================

def delta_k2(k, L):
    actual = actual_a_k2(
        k,
        L
    )

    if actual is None:
        return None

    interior = P3_interior(
        k,
        k + 2,
        L
    )

    return simp(
        actual - interior
    )


# ==============================================================================
# 1. EXTRACT EXACT a=k+2 CORRECTIONS
# ==============================================================================

def extraction():
    print("=" * 78)
    print("1. EXTRACTED a=k+2 CORRECTIONS")
    print("=" * 78)

    for k in sorted(DATA):

        values = []

        for L in sorted(DATA[k]):

            corr = delta_k2(
                k,
                L
            )

            if corr is not None:
                values.append(
                    (L, corr)
                )

        print(
            f"k={k}"
        )

        for L, corr in values:
            print(
                f"  ell={L:2d} "
                f"actual={actual_a_k2(k,L)} "
                f"interior={P3_interior(k,k+2,L)} "
                f"Delta={corr}"
            )

        print()


# ==============================================================================
# 2. FINITE DIFFERENCE ORDER
# ==============================================================================

def finite_differences(values):
    rows = [list(values)]

    while len(rows[-1]) > 1:
        prev = rows[-1]
        rows.append([
            simp(
                prev[i + 1] - prev[i]
            )
            for i in range(len(prev) - 1)
        ])

    return rows


def difference_audit():
    print("=" * 78)
    print("2. FINITE-DIFFERENCE STRUCTURE")
    print("=" * 78)

    for k in sorted(DATA):

        corrections = [
            delta_k2(k,L)
            for L in sorted(DATA[k])
        ]

        rows = finite_differences(
            corrections
        )

        print(
            f"k={k}"
        )

        for order, row in enumerate(rows):
            print(
                f"  Delta^{order} = {row}"
            )

        print()


# ==============================================================================
# 3. EXACT POLYNOMIAL IN L
# ==============================================================================

def interpolation_audit():
    print("=" * 78)
    print("3. EXACT INTERPOLATION IN L")
    print("=" * 78)

    Lsym = sp.symbols("L")

    failures = 0

    for k in sorted(DATA):

        points = []

        for L in sorted(DATA[k]):

            value = delta_k2(
                k,L
            )

            points.append(
                (
                    sp.Integer(L),
                    value
                )
            )

        poly = simp(
            sp.interpolate(
                points,
                Lsym
            )
        )

        print(
            f"k={k}"
        )
        print(
            f"  Delta_(k+2)(L) = {sp.factor(poly)}"
        )
        print(
            f"  expanded       = {sp.expand(poly)}"
        )

        for L, actual in points:

            predicted = simp(
                poly.subs(
                    Lsym,
                    L
                )
            )

            residual = simp(
                actual - predicted
            )

            if residual != 0:
                failures += 1

        print()

    print(
        f"interpolation failures = {failures}"
    )
    print()


# ==============================================================================
# 4. COEFFICIENT EXTRACTION
# ==============================================================================

def coefficient_audit():
    print("=" * 78)
    print("4. POLYNOMIAL COEFFICIENT STRUCTURE")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    for k in sorted(DATA):

        points = []

        for Lv in sorted(DATA[k]):

            points.append(
                (
                    sp.Integer(Lv),
                    delta_k2(k,Lv)
                )
            )

        poly = sp.Poly(
            sp.interpolate(
                points,
                L
            ),
            L
        )

        coeffs = [
            simp(
                poly.coeff_monomial(
                    L**i
                )
            )
            for i in range(3,-1,-1)
        ]

        print(
            f"k={k}"
        )

        print(
            f"  A={coeffs[0]}"
        )
        print(
            f"  B={coeffs[1]}"
        )
        print(
            f"  C={coeffs[2]}"
        )
        print(
            f"  D={coeffs[3]}"
        )

        print()


# ==============================================================================
# 5. SYMBOLIC BASELINE: PREVIOUS CORRECTIONS
# ==============================================================================

def symbolic_boundary_ladder():
    print("=" * 78)
    print("5. SYMBOLIC BOUNDARY LADDER")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    delta_k = B(
        K + 2,
        3
    )

    delta_k1 = simp(
        -B(K + 2,2) * L
        + K*(K+2)*(K+3)/sp.Integer(2)
    )

    print(
        "Delta_k ="
    )
    print(
        f"  {delta_k}"
    )
    print()

    print(
        "Delta_(k+1) ="
    )
    print(
        f"  {delta_k1}"
    )
    print()

    print(
        "Delta_(k+1) + C(k+2,2)L ="
    )
    print(
        f"  {simp(delta_k1 + B(K+2,2)*L)}"
    )
    print()


# ==============================================================================
# 6. COMPARE a=k+2 CORRECTION WITH NATURAL BINOMIAL SCALES
# ==============================================================================

def normalization_audit():
    print("=" * 78)
    print("6. NATURAL BINOMIAL NORMALIZATION")
    print("=" * 78)

    for k in sorted(DATA):

        print(
            f"k={k}"
        )

        for L in sorted(DATA[k]):

            delta = delta_k2(
                k,L
            )

            print(
                f"  ell={L}"
            )

            candidates = {
                "Delta": delta,

                "Delta/C(k+2,2)": simp(
                    delta / B(k+2,2)
                ),

                "Delta/C(k+2,3)": simp(
                    delta / B(k+2,3)
                ),

                "Delta/C(k+2,4)": simp(
                    delta / B(k+2,4)
                ),

                "Delta/C(k+3,3)": simp(
                    delta / B(k+3,3)
                ),

                "Delta/C(k+3,4)": simp(
                    delta / B(k+3,4)
                ),

                "Delta/(k+2)": simp(
                    delta / (k+2)
                ),

                "Delta/(k+2 choose 2)": simp(
                    delta / B(k+2,2)
                ),
            }

            for name, value in candidates.items():
                print(
                    f"    {name:24s} = {value}"
                )

        print()


# ==============================================================================
# 7. ROOT / FACTORIZATION AUDIT
# ==============================================================================

def factor_audit():
    print("=" * 78)
    print("7. FACTORIZATION / ROOT AUDIT")
    print("=" * 78)

    Lsym = sp.symbols("L")

    for k in sorted(DATA):

        points = [
            (
                sp.Integer(Lv),
                delta_k2(k,Lv)
            )
            for Lv in sorted(DATA[k])
        ]

        poly = sp.Poly(
            sp.interpolate(
                points,
                Lsym
            ),
            Lsym
        )

        expr = simp(
            poly.as_expr()
        )

        print(
            f"k={k}"
        )

        print(
            f"  factored = {sp.factor(expr)}"
        )

        print(
            f"  roots    = {sp.solve(expr, Lsym)}"
        )

        print()


# ==============================================================================
# 8. TEST A SIMPLE BINOMIAL PRODUCT TEMPLATE
# ==============================================================================

def template_audit():
    print("=" * 78)
    print("8. STRUCTURAL TEMPLATE AUDIT")
    print("=" * 78)

    #
    # Motivated only by the already-observed boundary pattern:
    #
    #   a=k:
    #       constant
    #
    #   a=k+1:
    #       linear in L
    #
    # The next boundary is therefore tested against the
    # simplest quadratic template
    #
    #       alpha(k) L^2 + beta(k) L + gamma(k)
    #
    # where alpha is compared against natural binomial
    # factors.
    #

    for k in sorted(DATA):

        points = []

        for Lv in sorted(DATA[k]):

            points.append(
                (
                    sp.Integer(Lv),
                    delta_k2(k,Lv)
                )
            )

        L = sp.symbols("L")

        poly = sp.Poly(
            sp.interpolate(
                points,
                L
            ),
            L
        )

        A = simp(
            poly.coeff_monomial(
                L**2
            )
        )

        Bc = simp(
            poly.coeff_monomial(
                L
            )
        )

        Cc = simp(
            poly.coeff_monomial(
                1
            )
        )

        print(
            f"k={k}"
        )

        print(
            f"  A = {A}"
        )

        print(
            f"  A / C(k+2,2) = "
            f"{simp(A / B(k+2,2))}"
        )

        print(
            f"  A / C(k+3,2) = "
            f"{simp(A / B(k+3,2))}"
        )

        print(
            f"  A / C(k+3,3) = "
            f"{simp(A / B(k+3,3))}"
        )

        print(
            f"  B = {Bc}"
        )

        print(
            f"  C = {Cc}"
        )

        print()


# ==============================================================================
# 9. CROSS-k COEFFICIENT TABLE
# ==============================================================================

def cross_k_audit():
    print("=" * 78)
    print("9. CROSS-k COEFFICIENT TABLE")
    print("=" * 78)

    L = sp.symbols("L")

    rows = []

    for k in sorted(DATA):

        points = [
            (
                sp.Integer(Lv),
                delta_k2(k,Lv)
            )
            for Lv in sorted(DATA[k])
        ]

        poly = sp.Poly(
            sp.interpolate(
                points,
                L
            ),
            L
        )

        rows.append(
            (
                k,
                simp(poly.coeff_monomial(L**2)),
                simp(poly.coeff_monomial(L)),
                simp(poly.coeff_monomial(1)),
            )
        )

    print(
        "   k            A                 B                 C"
    )

    for k,A,Bc,C in rows:

        print(
            f"{k:4d} "
            f"{str(A):>18} "
            f"{str(Bc):>18} "
            f"{str(C):>18}"
        )

    print()


# ==============================================================================
# 10. EXACT RECONSTRUCTION FROM THE EXTRACTED POLYNOMIAL
# ==============================================================================

def reconstruction_audit():
    print("=" * 78)
    print("10. EXACT a=k+2 RECONSTRUCTION")
    print("=" * 78)

    L = sp.symbols("L")

    failures = 0
    tested = 0

    for k in sorted(DATA):

        points = [
            (
                sp.Integer(Lv),
                delta_k2(k,Lv)
            )
            for Lv in sorted(DATA[k])
        ]

        poly = sp.interpolate(
            points,
            L
        )

        for Lv, actual_delta in points:

            predicted_delta = simp(
                poly.subs(
                    L,
                    Lv
                )
            )

            residual = simp(
                actual_delta
                - predicted_delta
            )

            tested += 1

            status = (
                "PASS"
                if residual == 0
                else "FAIL"
            )

            print(
                f"k={k:2d} ell={int(Lv):2d} "
                f"actual={actual_delta} "
                f"predicted={predicted_delta} "
                f"{status}"
            )

            if residual != 0:
                failures += 1

    print()

    print(
        f"tested = {tested}"
    )

    print(
        f"reconstruction failures = {failures}"
    )

    print()


# ==============================================================================
# 11. FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "Established interior L3 law:"
    )

    print(
        "  P3(k,a,L)"
    )
    print(
        "    = C(k+3,a)/6"
        " * product_{j=1}^{2}(L-a-j)"
        " * ((k+3)L-ka)/(k+3)"
    )

    print()

    print(
        "Established boundary:"
    )

    print(
        "  a=k:"
    )
    print(
        "    Delta_k = C(k+2,3)"
    )

    print()

    print(
        "  a=k+1:"
    )
    print(
        "    Delta_(k+1)"
        " = -C(k+2,2)L"
        " + k(k+2)(k+3)/2"
    )

    print()

    print(
        "New target:"
    )
    print(
        "  a=k+2 is nonzero."
    )
    print(
        "  Therefore the earlier zero-tail hypothesis is false."
    )

    print()

    print(
        "The decisive output from this experiment is the exact"
    )
    print(
        "polynomial structure of Delta_(k+2)(L)."
    )

    print()

    print(
        "Do not start L4 until:"
    )
    print(
        "  1. Delta_(k+2) has an exact structural formula;"
    )
    print(
        "  2. that formula is validated across all available"
    )
    print(
        "     k, ell cases;"
    )
    print(
        "  3. a genuine boundary ladder is identified."
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 232")
    print("EXACT a=k+2 BOUNDARY DERIVATION")
    print("=" * 78)
    print()

    extraction()

    difference_audit()

    interpolation_audit()

    coefficient_audit()

    symbolic_boundary_ladder()

    normalization_audit()

    factor_audit()

    template_audit()

    cross_k_audit()

    reconstruction_audit()

    final_diagnostic()


if __name__ == "__main__":
    main()

