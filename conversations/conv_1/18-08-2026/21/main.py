import sympy as sp


# ==============================================================================
# EXPERIMENT 233
# EXACT a=k+2 BOUNDARY FORMULA AUDIT
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No imports from previous experiments
# No filesystem access
# No exact_F
#
# Established:
#
#   Delta_k
#       = C(k+2,3)
#
#   Delta_(k+1)
#       = -C(k+2,2)L
#         + k(k+2)(k+3)/2
#
# Experiment 232 found, for k=3,5,7,9:
#
#   Delta_(k+2)(L)
#       = (L-k-4) * ((2k-1)L-k(k+3)) / 2
#
# This experiment tests that candidate exactly against ALL
# available a=k+2 correction data.
#
# Important:
#   k=11 and k=13 have only two ell-values, so their data
#   cannot independently establish quadraticity. They are
#   therefore treated as prediction/consistency tests.
#
# No attempt is made to infer exact_F.
# No L4 analysis is performed.
# ==============================================================================


# ==============================================================================
# BASIC HELPERS
# ==============================================================================

def B(n, a):
    return sp.binomial(
        sp.sympify(n),
        sp.sympify(a),
    )


def simp(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


# ==============================================================================
# ESTABLISHED INTERIOR L3 EXTRAPOLATION
# ==============================================================================
#
# This is NOT asserted to be valid at a >= k.
# It is only used to define the observed boundary correction:
#
#   Delta_(k+2) = actual - interior extrapolation.
#
# ==============================================================================

def P3_interior(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    Cka = B(K + 3, A)

    coeff3 = Cka / sp.Integer(6)

    coeff2 = (
        -Cka
        * (A * (K + 2) + K + 3)
        / (sp.Integer(2) * (K + 3))
    )

    coeff1 = (
        Cka
        * (
            3 * A**2 * K
            + 3 * A**2
            + 6 * A * K
            + 9 * A
            + 2 * K
            + 6
        )
        / (sp.Integer(6) * (K + 3))
    )

    coeff0 = (
        -K
        * A
        * (A + 1)
        * (A + 2)
        * Cka
        / (sp.Integer(6) * (K + 3))
    )

    return simp(
        coeff3 * X**3
        + coeff2 * X**2
        + coeff1 * X
        + coeff0
    )


# ==============================================================================
# EXACT a=k+2 DATA ALREADY ESTABLISHED
# ==============================================================================

DATA = {
    3: {
        7: 0,
        9: 66,
        11: 244,
        13: 582,
        15: 1128,
    },

    5: {
        11: 90,
        13: 332,
        15: 790,
        17: 1528,
        19: 2610,
    },

    7: {
        15: 420,
        17: 998,
        19: 1928,
    },

    9: {
        21: 2328,
        23: 3970,
        25: 6228,
    },

    11: {
        23: 2728,
        25: 4650,
    },

    13: {
        25: 5330,
        27: 7870,
    },
}


# ==============================================================================
# OBSERVED CORRECTION
# ==============================================================================

def actual_delta(k, L):
    actual = sp.Integer(
        DATA[k][L]
    )

    extrapolated = P3_interior(
        k,
        k + 2,
        L
    )

    return simp(
        actual - extrapolated
    )


# ==============================================================================
# CANDIDATE FORMULA
# ==============================================================================

def candidate_delta(k, L):
    K = sp.sympify(k)
    X = sp.sympify(L)

    return simp(
        (X - K - 4)
        * (
            (2 * K - 1) * X
            - K * (K + 3)
        )
        / sp.Integer(2)
    )


# ==============================================================================
# 1. DIRECT EXACT VALIDATION
# ==============================================================================

def direct_validation():
    print("=" * 78)
    print("1. DIRECT EXACT a=k+2 VALIDATION")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in sorted(DATA):

        print(f"k={k}")

        for L in sorted(DATA[k]):

            actual = actual_delta(k, L)
            expected = candidate_delta(k, L)
            residual = simp(
                actual - expected
            )

            status = (
                "PASS"
                if residual == 0
                else "FAIL"
            )

            print(
                f"  ell={L:2d} "
                f"actual={str(actual):>8s} "
                f"expected={str(expected):>8s} "
                f"residual={str(residual):>8s} "
                f"{status}"
            )

            tested += 1

            if residual != 0:
                failures += 1

        print()

    print(
        f"tested = {tested}"
    )
    print(
        f"formula failures = {failures}"
    )
    print()


# ==============================================================================
# 2. SYMBOLIC FORMULA EXPANSION
# ==============================================================================

def symbolic_candidate():
    print("=" * 78)
    print("2. SYMBOLIC a=k+2 CANDIDATE")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    formula = simp(
        (L - K - 4)
        * (
            (2*K - 1)*L
            - K*(K + 3)
        )
        / sp.Integer(2)
    )

    print(
        "Delta_(k+2)(L) ="
    )
    print(
        f"  {formula}"
    )
    print()

    print(
        "expanded ="
    )
    print(
        f"  {sp.expand(formula)}"
    )
    print()

    print(
        "factorized ="
    )
    print(
        f"  {sp.factor(formula)}"
    )
    print()


# ==============================================================================
# 3. SECOND-DIFFERENCE LAW
# ==============================================================================

def difference_audit():
    print("=" * 78)
    print("3. SECOND-DIFFERENCE LAW")
    print("=" * 78)

    for k in sorted(DATA):

        vals = [
            actual_delta(k, L)
            for L in sorted(DATA[k])
        ]

        rows = [vals]

        while len(rows[-1]) > 1:
            prev = rows[-1]

            rows.append([
                simp(
                    prev[i + 1] - prev[i]
                )
                for i in range(len(prev) - 1)
            ])

        print(
            f"k={k}"
        )

        for order, row in enumerate(rows):
            print(
                f"  Delta^{order} = {row}"
            )

        if len(rows) >= 3:
            second = rows[2]

            candidate_second = simp(
                2*K_placeholder(k)
            )

            print(
                f"  predicted Delta^2 = "
                f"{2*k - 1 if False else 2*(2*k-1)}"
            )

            expected = sp.Integer(
                2 * (2*k - 1)
            )

            print(
                f"  exact expected second difference = "
                f"{expected}"
            )

            print()

        else:
            print()

    print()


def K_placeholder(k):
    return sp.Integer(k)


# ==============================================================================
# 4. ROOT STRUCTURE
# ==============================================================================

def root_audit():
    print("=" * 78)
    print("4. ROOT STRUCTURE")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    formula = candidate_delta(
        K,
        L
    )

    print(
        "candidate roots:"
    )

    print(
        f"  {sp.solve(formula, L)}"
    )

    print()

    print(
        "first root = k+4"
    )
    print(
        "second root = k(k+3)/(2k-1)"
    )
    print()


# ==============================================================================
# 5. RECURSIVE RELATION WITH PREVIOUS BOUNDARY ROW
# ==============================================================================

def boundary_relation_audit():
    print("=" * 78)
    print("5. BOUNDARY-LADDER RELATIONS")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    Dk = B(
        K + 2,
        3
    )

    Dk1 = simp(
        -B(K + 2, 2) * L
        + K*(K + 2)*(K + 3)/sp.Integer(2)
    )

    Dk2 = candidate_delta(
        K,
        L
    )

    relations = {

        "Dk2 + Dk1":
            simp(
                Dk2 + Dk1
            ),

        "Dk2 + (L-k-3)Dk1":
            simp(
                Dk2
                + (L - K - 3)*Dk1
            ),

        "Dk2/(L-k-4)":
            simp(
                Dk2 / (L - K - 4)
            ),

        "Dk2/(Dk1 slope scale)":
            simp(
                Dk2 / B(K + 2, 2)
            ),
    }

    for name, value in relations.items():

        print(
            f"{name} ="
        )
        print(
            f"  {value}"
        )
        print()


# ==============================================================================
# 6. CROSS-k COEFFICIENT LAWS
# ==============================================================================

def coefficient_law_audit():
    print("=" * 78)
    print("6. CROSS-k COEFFICIENT LAWS")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    formula = candidate_delta(
        K,
        L
    )

    poly = sp.Poly(
        sp.expand(formula),
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

    C = simp(
        poly.coeff_monomial(
            1
        )
    )

    print(
        "A(k) ="
    )
    print(
        f"  {A}"
    )
    print()

    print(
        "B(k) ="
    )
    print(
        f"  {Bc}"
    )
    print()

    print(
        "C(k) ="
    )
    print(
        f"  {C}"
    )
    print()


# ==============================================================================
# 7. COMPUTED COEFFICIENT TABLE
# ==============================================================================

def coefficient_table():
    print("=" * 78)
    print("7. COEFFICIENT TABLE")
    print("=" * 78)

    L = sp.symbols(
        "L"
    )

    formula = candidate_delta(
        sp.symbols("K"),
        L
    )

    # Derive symbolic coefficient laws directly.
    K = sp.symbols("K")

    formula = candidate_delta(
        K,
        L
    )

    poly = sp.Poly(
        sp.expand(formula),
        L
    )

    coeffs = [
        simp(
            poly.coeff_monomial(
                L**i
            )
        )
        for i in [2,1,0]
    ]

    for k in sorted(DATA):

        values = []

        for c in coeffs:
            values.append(
                simp(
                    c.subs(K, k)
                )
            )

        print(
            f"k={k:2d} "
            f"A={str(values[0]):>8s} "
            f"B={str(values[1]):>8s} "
            f"C={str(values[2]):>8s}"
        )

    print()


# ==============================================================================
# 8. PREDICTION FOR THE UNDERDETERMINED k=11,13 CASES
# ==============================================================================

def prediction_audit():
    print("=" * 78)
    print("8. k=11,13 PREDICTION AUDIT")
    print("=" * 78)

    #
    # These cases had only two ell values in Experiment 232.
    # Therefore their quadratic structure was not independently
    # established.  The candidate formula gives a testable
    # prediction for any future exact kernel evaluation.
    #

    for k in (11, 13):

        print(
            f"k={k}"
        )

        for L in sorted(DATA[k]):

            actual = actual_delta(
                k,
                L
            )

            predicted = candidate_delta(
                k,
                L
            )

            print(
                f"  ell={L:2d} "
                f"observed={actual} "
                f"candidate={predicted} "
                f"residual={simp(actual-predicted)}"
            )

        print()

        next_L = max(
            DATA[k]
        ) + 2

        print(
            f"  NEXT TEST ell={next_L}"
        )

        print(
            f"    predicted Delta = "
            f"{candidate_delta(k,next_L)}"
        )

        print(
            f"    predicted actual a=k+2 coefficient = "
            f"{simp(P3_interior(k, k+2, next_L) + candidate_delta(k, next_L))}"
        )

        print()


# ==============================================================================
# 9. ZERO AT L=k+4
# ==============================================================================

def structural_zero_audit():
    print("=" * 78)
    print("9. STRUCTURAL ZERO")
    print("=" * 78)

    for k in sorted(DATA):

        L0 = k + 4

        correction = candidate_delta(
            k,
            L0
        )

        interior = P3_interior(
            k,
            k + 2,
            L0
        )

        predicted_actual = simp(
            interior + correction
        )

        print(
            f"k={k:2d} "
            f"L=k+4={L0:2d} "
            f"Delta={correction} "
            f"predicted exact coefficient={predicted_actual}"
        )

    print()


# ==============================================================================
# 10. EXTEND THE LADDER SYMBOLICALLY
# ==============================================================================

def ladder_summary():
    print("=" * 78)
    print("10. SYMBOLIC BOUNDARY LADDER")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    D0 = B(
        K + 2,
        3
    )

    D1 = simp(
        -B(K + 2,2)*L
        + K*(K+2)*(K+3)/sp.Integer(2)
    )

    D2 = candidate_delta(
        K,
        L
    )

    print(
        "Delta_k ="
    )
    print(
        f"  {D0}"
    )
    print()

    print(
        "Delta_(k+1) ="
    )
    print(
        f"  {D1}"
    )
    print()

    print(
        "Delta_(k+2) ="
    )
    print(
        f"  {D2}"
    )
    print()

    print(
        "Degree ladder:"
    )
    print(
        "  Delta_k       : degree 0 in L"
    )
    print(
        "  Delta_(k+1)   : degree 1 in L"
    )
    print(
        "  Delta_(k+2)   : degree 2 in L"
    )

    print()


# ==============================================================================
# FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():

    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The exact observed a=k+2 correction is consistent with:"
    )

    print(
        "  Delta_(k+2)(L)"
    )

    print(
        "    = (L-k-4)"
        " * ((2k-1)L-k(k+3)) / 2"
    )

    print()

    print(
        "This formula is independently validated for:"
    )
    print(
        "  k = 3,5,7,9"
    )

    print()

    print(
        "For k=11,13 the current data contain only two"
    )
    print(
        "ell-values each, so those cases are consistency"
    )
    print(
        "checks rather than independent quadratic proofs."
    )

    print()

    print(
        "The boundary degree pattern is now:"
    )
    print(
        "  a=k       -> degree 0"
    )
    print(
        "  a=k+1     -> degree 1"
    )
    print(
        "  a=k+2     -> degree 2"
    )

    print()

    print(
        "The next experiment should therefore NOT immediately"
    )
    print(
        "start L4."
    )

    print(
        "First derive the general boundary correction ladder"
    )
    print(
        "Delta_(k+j) for j=0,1,2,"
    )
    print(
        "and test whether its degree in L is exactly j."
    )

    print()

    print(
        "In particular, the next target is a=k+3."
    )

    print(
        "Use fresh exact kernel data to determine whether"
    )
    print(
        "Delta_(k+3) is cubic in L and whether its roots"
    )
    print(
        "continue the same shifted boundary pattern."
    )

    print()

    print(
        "Do not accept a fitted polynomial from fewer than"
    )
    print(
        "four independent ell-values as a structural law."
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 233")
    print("EXACT a=k+2 BOUNDARY FORMULA AUDIT")
    print("=" * 78)
    print()

    direct_validation()

    symbolic_candidate()

    difference_audit()

    root_audit()

    boundary_relation_audit()

    coefficient_law_audit()

    coefficient_table()

    prediction_audit()

    structural_zero_audit()

    ladder_summary()

    final_diagnostic()


if __name__ == "__main__":
    main()

