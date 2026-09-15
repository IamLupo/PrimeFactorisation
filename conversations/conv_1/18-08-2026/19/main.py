import sympy as sp


# ==============================================================================
# EXPERIMENT 231
# EXACT L3 BOUNDARY CORRECTION LADDER
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No imports from previous experiments
# No filesystem access
# No exact_F
#
# Purpose:
#
#   1. Verify the established a=k correction:
#          Delta_k = C(k+2,3)
#
#   2. Extract the exact a=k+1 correction.
#
#   3. Test the candidate:
#          Delta_(k+1)
#             = -C(k+2,2)*L
#               + k(k+2)(k+3)/2
#
#   4. Verify the observed zero tail a>=k+2.
#
#   5. Compute the cancellation that would be required at a=k+2.
#
#   6. Look for a structural relation between the boundary corrections.
#
# No L4 analysis.
# No arbitrary candidate enumeration.
# No floating-point arithmetic.
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


def eq_zero(expr):
    return simp(expr) == 0


def B(n, a):
    return sp.binomial(
        sp.sympify(n),
        sp.sympify(a)
    )


# ==============================================================================
# ESTABLISHED INTERIOR L3 LAW
# ==============================================================================

def P3_interior(k, a, L):
    """
    Established exact interior L3 coefficient law:

      P_3(k,a,L)
        = A L^3 + B L^2 + C L + D

    with

      A = C(k+3,a)/6

      B = -C(k+3,a)
          * (a(k+2)+k+3)
          / (2(k+3))

      C = C(k+3,a)
          * (3a^2k + 3a^2 + 6ak + 9a + 2k + 6)
          / (6(k+3))

      D = -k a(a+1)(a+2) C(k+3,a)
          / (6(k+3))

    The formula is used only as the established interior
    extrapolation, including when evaluating boundary points.
    """

    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    choose = B(K + 3, A)

    term3 = (
        choose
        * X**3
        / sp.Integer(6)
    )

    term2 = (
        -choose
        * (
            A * (K + 2) + K + 3
        )
        * X**2
        / (
            sp.Integer(2)
            * (K + 3)
        )
    )

    term1 = (
        choose
        * (
            3 * A**2 * K
            + 3 * A**2
            + 6 * A * K
            + 9 * A
            + 2 * K
            + 6
        )
        * X
        / (
            sp.Integer(6)
            * (K + 3)
        )
    )

    term0 = (
        -K
        * A
        * (A + 1)
        * (A + 2)
        * choose
        / (
            sp.Integer(6)
            * (K + 3)
        )
    )

    return simp(
        term3 + term2 + term1 + term0
    )


# ==============================================================================
# EXACT DATA ALREADY ESTABLISHED IN EXPERIMENTS 219-224
#
# Only values explicitly present in the supplied experimental output are used.
# Missing values are simply absent and are never invented.
#
# DATA[k][ell][a] = exact L3 coefficient [N^a X^(ell-3-a)]
# ==============================================================================

DATA = {

    3: {
        7: {
            3: sp.Integer(120),
            4: sp.Integer(0),
        },

        9: {
            3: sp.Integer(510),
            4: sp.Integer(165),
            5: sp.Integer(66),
            6: sp.Integer(0),
        },

        11: {
            3: sp.Integer(1340),
            4: sp.Integer(610),
            5: sp.Integer(244),
            6: sp.Integer(0),
        },

        13: {
            3: sp.Integer(2770),
            4: sp.Integer(1455),
            5: sp.Integer(582),
            6: sp.Integer(0),
        },

        15: {
            3: sp.Integer(4960),
            4: sp.Integer(2820),
            5: sp.Integer(1128),
            6: sp.Integer(0),
        },
    },

    5: {
        11: {
            5: sp.Integer(1505),
            6: sp.Integer(315),
            7: sp.Integer(90),
            8: sp.Integer(0),
        },

        13: {
            5: sp.Integer(3906),
            6: sp.Integer(1162),
            7: sp.Integer(332),
            8: sp.Integer(0),
        },

        15: {
            5: sp.Integer(8015),
            6: sp.Integer(2765),
            7: sp.Integer(790),
            8: sp.Integer(0),
        },

        17: {
            5: sp.Integer(14280),
            6: sp.Integer(5348),
            7: sp.Integer(1528),
            8: sp.Integer(0),
        },

        19: {
            5: sp.Integer(23149),
            6: sp.Integer(9135),
            7: sp.Integer(2610),
            8: sp.Integer(0),
        },
    },

    7: {
        15: {
            7: sp.Integer(8568),
            8: sp.Integer(1890),
            9: sp.Integer(420),
            10: sp.Integer(0),
        },

        17: {
            7: sp.Integer(17508),
            8: sp.Integer(4491),
            9: sp.Integer(998),
            10: sp.Integer(0),
        },

        19: {
            7: sp.Integer(31104),
            8: sp.Integer(8676),
            9: sp.Integer(1928),
            10: sp.Integer(0),
        },
    },

    9: {
        21: {
            9: sp.Integer(57640),
            10: sp.Integer(12804),
            11: sp.Integer(2328),
            12: sp.Integer(0),
        },

        23: {
            9: sp.Integer(93115),
            10: sp.Integer(21835),
            11: sp.Integer(3970),
            12: sp.Integer(0),
        },

        25: {
            9: sp.Integer(140690),
            10: sp.Integer(34254),
            11: sp.Integer(6228),
            12: sp.Integer(0),
        },
    },

    11: {
        23: {
            11: sp.Integer(96096),
            12: sp.Integer(17732),
            13: sp.Integer(2728),
            14: sp.Integer(0),
        },

        25: {
            11: sp.Integer(155090),
            12: sp.Integer(30225),
            13: sp.Integer(4650),
            14: sp.Integer(0),
        },
    },

    13: {
        25: {
            13: sp.Integer(148680),
            14: sp.Integer(23460),
            15: sp.Integer(5330),
            16: sp.Integer(0),
        },

        27: {
            13: sp.Integer(239785),
            14: sp.Integer(39975),
            16: sp.Integer(0),
        },
    },
}


# ==============================================================================
# ACCESSOR
# ==============================================================================

def actual_coefficient(k, L, a):
    return DATA.get(
        k,
        {}
    ).get(
        L,
        {}
    ).get(
        a,
        None
    )


# ==============================================================================
# a=k CORRECTION
# ==============================================================================

def correction_a_k(k, L):
    actual = actual_coefficient(
        k,
        L,
        k
    )

    if actual is None:
        return None

    naive = P3_interior(
        k,
        k,
        L
    )

    return simp(
        actual - naive
    )


# ==============================================================================
# a=k+1 CORRECTION
# ==============================================================================

def correction_a_k1(k, L):
    actual = actual_coefficient(
        k,
        L,
        k + 1
    )

    if actual is None:
        return None

    naive = P3_interior(
        k,
        k + 1,
        L
    )

    return simp(
        actual - naive
    )


# ==============================================================================
# SUSPECTED EXACT a=k+1 CORRECTION
# ==============================================================================

def candidate_delta_k1(k, L):
    K = sp.sympify(k)
    X = sp.sympify(L)

    return simp(
        -B(K + 2, 2) * X
        +
        K * (K + 2) * (K + 3)
        / sp.Integer(2)
    )


# ==============================================================================
# REQUIRED a=k+2 CANCELLATION IF THE EXACT TAIL IS ZERO
# ==============================================================================

def required_delta_k2(k, L):
    return simp(
        -P3_interior(
            k,
            k + 2,
            L
        )
    )


# ==============================================================================
# 1. a=k EXACT CORRECTION AUDIT
# ==============================================================================

def audit_a_k():
    print("=" * 78)
    print("1. EXACT a=k CORRECTION")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in sorted(DATA):

        expected = B(
            k + 2,
            3
        )

        print(
            f"k={k} expected correction={expected}"
        )

        for L in sorted(DATA[k]):

            actual = correction_a_k(
                k,
                L
            )

            if actual is None:
                continue

            tested += 1

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
                f"actual={actual} "
                f"correction={actual} "
                f"expected={expected} "
                f"{status}"
            )

            if residual != 0:
                failures += 1

        print()

    print(
        f"a=k correction failures = {failures}"
    )
    print(
        f"a=k cases tested = {tested}"
    )
    print()


# ==============================================================================
# 2. EXTRACT a=k+1 CORRECTIONS
# ==============================================================================

def extract_a_k1():
    print("=" * 78)
    print("2. EXTRACTED a=k+1 CORRECTIONS")
    print("=" * 78)

    for k in sorted(DATA):

        values = []

        for L in sorted(DATA[k]):

            corr = correction_a_k1(
                k,
                L
            )

            if corr is not None:
                values.append(
                    (L, corr)
                )

        if not values:
            continue

        print(
            f"k={k}"
        )

        for L, corr in values:
            print(
                f"  ell={L:2d} "
                f"Delta_(k+1)={corr}"
            )

        print()


# ==============================================================================
# 3. EXACT a=k+1 CANDIDATE AUDIT
# ==============================================================================

def audit_a_k1_candidate():
    print("=" * 78)
    print("3. EXACT a=k+1 CORRECTION LAW")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in sorted(DATA):

        for L in sorted(DATA[k]):

            actual = correction_a_k1(
                k,
                L
            )

            if actual is None:
                continue

            expected = candidate_delta_k1(
                k,
                L
            )

            residual = simp(
                actual - expected
            )

            tested += 1

            print(
                f"k={k:2d} ell={L:2d} "
                f"actual={actual} "
                f"expected={expected} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

            if residual != 0:
                failures += 1

    print()

    print(
        f"tested = {tested}"
    )

    print(
        f"candidate failures = {failures}"
    )

    print()


# ==============================================================================
# 4. SYMBOLIC a=k+1 FORM
# ==============================================================================

def symbolic_a_k1():
    print("=" * 78)
    print("4. SYMBOLIC a=k+1 CORRECTION")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    delta = candidate_delta_k1(
        K,
        L
    )

    print(
        "Delta_(k+1)(L) ="
    )
    print(
        f"  {simp(delta)}"
    )
    print()

    print(
        "expanded ="
    )
    print(
        f"  {sp.expand(delta)}"
    )
    print()

    slope = simp(
        sp.diff(
            delta,
            L
        )
    )

    intercept = simp(
        delta.subs(
            L,
            0
        )
    )

    second_derivative = simp(
        sp.diff(
            delta,
            L,
            2
        )
    )

    print(
        "d/dL Delta_(k+1) ="
    )
    print(
        f"  {slope}"
    )
    print()

    print(
        "Delta_(k+1)(0) ="
    )
    print(
        f"  {intercept}"
    )
    print()

    print(
        "d^2/dL^2 Delta_(k+1) ="
    )
    print(
        f"  {second_derivative}"
    )
    print()


# ==============================================================================
# 5. FINITE-DIFFERENCE AUDIT OF a=k+1
# ==============================================================================

def finite_difference_a_k1():
    print("=" * 78)
    print("5. FINITE-DIFFERENCE STRUCTURE OF a=k+1")
    print("=" * 78)

    for k in sorted(DATA):

        points = []

        for L in sorted(DATA[k]):

            corr = correction_a_k1(
                k,
                L
            )

            if corr is not None:
                points.append(
                    corr
                )

        if len(points) < 3:
            continue

        d1 = [
            simp(
                points[i + 1]
                - points[i]
            )
            for i in range(len(points) - 1)
        ]

        d2 = [
            simp(
                d1[i + 1]
                - d1[i]
            )
            for i in range(len(d1) - 1)
        ]

        print(
            f"k={k}"
        )

        print(
            f"  Delta={points}"
        )

        print(
            f"  Delta^1={d1}"
        )

        print(
            f"  Delta^2={d2}"
        )

        expected_d2 = sp.Integer(0)

        status = all(
            value == expected_d2
            for value in d2
        )

        print(
            f"  quadratic? {'YES' if status else 'NO'}"
        )
        print()


# ==============================================================================
# 6. a=k+2 ZERO TAIL AUDIT
# ==============================================================================

def audit_zero_tail():
    print("=" * 78)
    print("6. a=k+2 ZERO-TAIL AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in sorted(DATA):

        for L in sorted(DATA[k]):

            actual = actual_coefficient(
                k,
                L,
                k + 2
            )

            if actual is None:
                continue

            tested += 1

            residual = simp(
                actual
            )

            print(
                f"k={k:2d} ell={L:2d} "
                f"a={k+2:2d} "
                f"actual={actual} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

            if residual != 0:
                failures += 1

    print()

    print(
        f"tested = {tested}"
    )

    print(
        f"zero-tail failures = {failures}"
    )

    print()


# ==============================================================================
# 7. REQUIRED a=k+2 CANCELLATION
# ==============================================================================

def symbolic_k2_cancellation():
    print("=" * 78)
    print("7. REQUIRED a=k+2 CANCELLATION")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    interior = simp(
        P3_interior(
            K,
            K + 2,
            L
        )
    )

    required = simp(
        -interior
    )

    print(
        "Interior extrapolation at a=k+2:"
    )
    print(
        f"  {interior}"
    )
    print()

    print(
        "Required correction for an exact zero:"
    )
    print(
        f"  {required}"
    )
    print()

    print(
        "factorized:"
    )
    print(
        f"  {sp.factor(required)}"
    )
    print()


# ==============================================================================
# 8. COMPARE BOUNDARY CORRECTIONS
# ==============================================================================

def boundary_comparison():
    print("=" * 78)
    print("8. BOUNDARY CORRECTION COMPARISON")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    delta_k = simp(
        B(
            K + 2,
            3
        )
    )

    delta_k1 = candidate_delta_k1(
        K,
        L
    )

    delta_k2 = simp(
        required_delta_k2(
            K,
            L
        )
    )

    print(
        "Delta_k:"
    )
    print(
        f"  {delta_k}"
    )
    print()

    print(
        "Delta_(k+1):"
    )
    print(
        f"  {delta_k1}"
    )
    print()

    print(
        "Required Delta_(k+2):"
    )
    print(
        f"  {delta_k2}"
    )
    print()

    print(
        "Delta_(k+1) + C(k+2,2)L:"
    )
    print(
        f"  {simp(delta_k1 + B(K + 2, 2) * L)}"
    )
    print()


# ==============================================================================
# 9. k-SCALING AUDIT
# ==============================================================================

def k_scaling_audit():
    print("=" * 78)
    print("9. k-SCALING AUDIT")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    delta = candidate_delta_k1(
        K,
        L
    )

    expressions = {
        "Delta_(k+1)/(k+2)": simp(
            delta / (K + 2)
        ),

        "Delta_(k+1)/C(k+2,2)": simp(
            delta / B(K + 2, 2)
        ),

        "-dDelta/dL": simp(
            -sp.diff(
                delta,
                L
            )
        ),

        "intercept": simp(
            delta.subs(
                L,
                0
            )
        ),

        "intercept / (k+2)": simp(
            delta.subs(
                L,
                0
            )
            / (K + 2)
        ),

        "intercept / C(k+2,2)": simp(
            delta.subs(
                L,
                0
            )
            / B(K + 2, 2)
        ),
    }

    for name, expression in expressions.items():

        print(
            f"{name}:"
        )

        print(
            f"  {expression}"
        )

        print()


# ==============================================================================
# 10. NUMERICAL RECONSTRUCTION USING THE CANDIDATE
# ==============================================================================

def reconstruction_audit():
    print("=" * 78)
    print("10. COMPLETE a=k+1 RECONSTRUCTION")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in sorted(DATA):

        for L in sorted(DATA[k]):

            actual = actual_coefficient(
                k,
                L,
                k + 1
            )

            if actual is None:
                continue

            interior = P3_interior(
                k,
                k + 1,
                L
            )

            correction = candidate_delta_k1(
                k,
                L
            )

            predicted = simp(
                interior + correction
            )

            residual = simp(
                actual - predicted
            )

            tested += 1

            print(
                f"k={k:2d} ell={L:2d} "
                f"actual={actual} "
                f"predicted={predicted} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
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
        "Established interior mechanism:"
    )
    print(
        "  P_r(k,a,L)"
    )
    print(
        "    = (-1)^(r+1)/r!"
    )
    print(
        "      * C(k+r,a)"
    )
    print(
        "      * product_{j=1}^{r-1}(L-a-j)"
    )
    print(
        "      * ((k+r)L-ka)/(k+r)"
    )
    print()

    print(
        "Established L3 boundary:"
    )
    print(
        "  a=k:"
    )
    print(
        "    Delta_k = C(k+2,3)"
    )
    print()

    print(
        "Working a=k+1 boundary candidate:"
    )
    print(
        "  Delta_(k+1)"
    )
    print(
        "    = -C(k+2,2)L"
    )
    print(
        "      + k(k+2)(k+3)/2"
    )
    print()

    print(
        "Observed tail:"
    )
    print(
        "  a >= k+2 -> coefficient = 0"
    )
    print()

    print(
        "This experiment therefore asks whether the a=k+1"
    )
    print(
        "correction is an exact boundary law and whether"
    )
    print(
        "the a=k+2 cancellation reveals the next member"
    )
    print(
        "of the boundary sequence."
    )
    print()

    print(
        "No L4 analysis is performed."
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 231")
    print("EXACT L3 BOUNDARY CORRECTION LADDER")
    print("=" * 78)
    print()

    audit_a_k()

    extract_a_k1()

    audit_a_k1_candidate()

    symbolic_a_k1()

    finite_difference_a_k1()

    audit_zero_tail()

    symbolic_k2_cancellation()

    boundary_comparison()

    k_scaling_audit()

    reconstruction_audit()

    final_diagnostic()


if __name__ == "__main__":
    main()