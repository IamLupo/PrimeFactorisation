import sympy as sp


# ==============================================================================
# EXPERIMENT 236
# k=13 a=k+2 BOUNDARY ANOMALY / FINITE-DIFFERENCE CONSISTENCY AUDIT
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No imports from previous experiments
# No filesystem access
# No exact_F required
#
# Purpose:
#   1. Recompute the a=k+2 correction exactly from the published L3 data.
#   2. Correctly audit step-2 finite differences.
#   3. Separate the interpolation law from the actual observed correction.
#   4. Determine exactly what additional fresh data are required.
#   5. Test whether the k=13 anomaly is compatible with:
#        - the known quadratic law,
#        - a shifted quadratic,
#        - a quadratic with the same leading coefficient,
#        - or an unresolved boundary transition.
#
# No new boundary formula is accepted.
# ==============================================================================


def simp(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def C(n, a):
    return sp.binomial(
        sp.sympify(n),
        sp.sympify(a)
    )


# ==============================================================================
# ESTABLISHED INTERIOR L3 FORMULA
# ==============================================================================

def P3_interior(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        C(K + 3, A)
        / sp.Integer(6)
        * (X - A - 1)
        * (X - A - 2)
        * ((K + 3) * X - K * A)
        / (K + 3)
    )


# ==============================================================================
# OBSERVED a=k+2 L3 COEFFICIENTS
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
# CURRENT OBSERVED CORRECTION
# ==============================================================================

def observed_delta(k, L):
    actual = sp.Integer(
        DATA[k][L]
    )

    interior = P3_interior(
        k,
        k + 2,
        L
    )

    return simp(
        actual - interior
    )


# ==============================================================================
# CURRENT QUADRATIC CANDIDATE
# ==============================================================================

def candidate_delta(k, L):
    K = sp.sympify(k)
    X = sp.sympify(L)

    return simp(
        (X - K - 4)
        * (
            (K + 2) * X
            - K * (K + 3)
        )
        / sp.Integer(2)
    )


# ==============================================================================
# 1. EXACT OBSERVED CORRECTIONS
# ==============================================================================

def observed_audit():
    print("=" * 78)
    print("1. EXACT OBSERVED a=k+2 CORRECTIONS")
    print("=" * 78)

    for k in sorted(DATA):

        print(
            f"k={k}"
        )

        for L in sorted(DATA[k]):

            actual = sp.Integer(
                DATA[k][L]
            )

            interior = P3_interior(
                k,
                k + 2,
                L
            )

            delta = observed_delta(
                k,
                L
            )

            print(
                f"  ell={L:2d} "
                f"actual={str(actual):>8s} "
                f"interior={str(interior):>8s} "
                f"Delta={str(delta):>8s}"
            )

        print()


# ==============================================================================
# 2. CORRECT FINITE-DIFFERENCE AUDIT
# ==============================================================================

def finite_difference_audit():
    print("=" * 78)
    print("2. CORRECT STEP-2 FINITE-DIFFERENCE AUDIT")
    print("=" * 78)

    print(
        "For D(L)=A L^2+B L+C with step h=2:"
    )
    print(
        "  Delta_h^2 D = 2*A*h^2 = 8A"
    )
    print()

    for k in sorted(DATA):

        points = sorted(
            DATA[k]
        )

        if len(points) < 3:
            print(
                f"k={k}: insufficient points"
            )
            continue

        values = [
            observed_delta(
                k,
                L
            )
            for L in points
        ]

        d1 = [
            simp(
                values[i + 1] - values[i]
            )
            for i in range(
                len(values) - 1
            )
        ]

        d2 = [
            simp(
                d1[i + 1] - d1[i]
            )
            for i in range(
                len(d1) - 1
            )
        ]

        print(
            f"k={k}"
        )
        print(
            f"  Delta^1 = {d1}"
        )
        print(
            f"  Delta^2 = {d2}"
        )
        print(
            f"  expected candidate = {4 * (k + 2)}"
        )

        status = all(
            x == 4 * (k + 2)
            for x in d2
        )

        print(
            f"  candidate quadratic second-difference = "
            f"{'PASS' if status else 'FAIL'}"
        )
        print()


# ==============================================================================
# 3. CANDIDATE COMPARISON
# ==============================================================================

def candidate_audit():
    print("=" * 78)
    print("3. CANDIDATE LAW AUDIT")
    print("=" * 78)

    tested = 0
    failures = 0

    for k in sorted(DATA):

        for L in sorted(DATA[k]):

            actual = observed_delta(
                k,
                L
            )

            expected = candidate_delta(
                k,
                L
            )

            residual = simp(
                actual - expected
            )

            tested += 1

            if residual != 0:
                failures += 1

            print(
                f"k={k:2d} ell={L:2d} "
                f"actual={str(actual):>8s} "
                f"candidate={str(expected):>8s} "
                f"residual={str(residual):>8s} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

        print()

    print(
        f"tested = {tested}"
    )
    print(
        f"failures = {failures}"
    )
    print()


# ==============================================================================
# 4. INTERPOLATE ONLY WHEN MATHEMATICALLY JUSTIFIED
# ==============================================================================

def interpolation_audit():
    print("=" * 78)
    print("4. INDEPENDENT QUADRATIC INTERPOLATION")
    print("=" * 78)

    L = sp.symbols(
        "L"
    )

    for k in sorted(DATA):

        points = sorted(
            DATA[k]
        )

        if len(points) < 3:

            print(
                f"k={k}: only {len(points)} point(s); "
                f"quadratic law NOT testable"
            )
            continue

        samples = [
            (
                sp.Integer(x),
                observed_delta(
                    k,
                    x
                )
            )
            for x in points
        ]

        q = simp(
            sp.interpolate(
                samples,
                L
            )
        )

        residuals = [
            simp(
                q.subs(
                    L,
                    x
                ) - y
            )
            for x, y in samples
        ]

        print(
            f"k={k}"
        )
        print(
            f"  quadratic = {sp.factor(q)}"
        )
        print(
            f"  residuals = {residuals}"
        )

        if len(samples) >= 4:

            leave_one_out = True

            for omitted in range(
                len(samples)
            ):

                reduced = (
                    samples[:omitted]
                    + samples[omitted + 1:]
                )

                q2 = simp(
                    sp.interpolate(
                        reduced,
                        L
                    )
                )

                x, y = samples[omitted]

                check = simp(
                    q2.subs(
                        L,
                        x
                    ) - y
                )

                if check != 0:
                    leave_one_out = False

            print(
                "  leave-one-out consistency = "
                + (
                    "PASS"
                    if leave_one_out
                    else "FAIL"
                )
            )

        print()


# ==============================================================================
# 5. k=13 ANOMALY DECOMPOSITION
# ==============================================================================

def k13_audit():
    print("=" * 78)
    print("5. k=13 ANOMALY DECOMPOSITION")
    print("=" * 78)

    k = 13

    for L in sorted(DATA[k]):

        actual = observed_delta(
            k,
            L
        )

        candidate = candidate_delta(
            k,
            L
        )

        residual = simp(
            actual - candidate
        )

        print(
            f"ell={L}"
        )
        print(
            f"  observed Delta = {actual}"
        )
        print(
            f"  candidate Delta = {candidate}"
        )
        print(
            f"  residual = {residual}"
        )
        print()

    print(
        "There are only two k=13 points."
    )
    print(
        "Therefore they cannot establish the degree in ell."
    )
    print()


# ==============================================================================
# 6. SAME-LEADING-COEFFICIENT TEST
# ==============================================================================

def same_leading_coefficient_test():
    print("=" * 78)
    print("6. SAME-LEADING-COEFFICIENT TEST")
    print("=" * 78)

    L = sp.symbols(
        "L"
    )

    for k in sorted(DATA):

        points = sorted(
            DATA[k]
        )

        if len(points) < 3:
            continue

        samples = [
            (
                sp.Integer(x),
                observed_delta(
                    k,
                    x
                )
            )
            for x in points
        ]

        q = sp.Poly(
            sp.expand(
                sp.interpolate(
                    samples,
                    L
                )
            ),
            L
        )

        A_actual = simp(
            q.coeff_monomial(
                L**2
            )
        )

        A_candidate = simp(
            sp.Rational(
                k + 2,
                2
            )
        )

        print(
            f"k={k} "
            f"A_actual={A_actual} "
            f"A_candidate={A_candidate} "
            f"{'PASS' if A_actual == A_candidate else 'FAIL'}"
        )

    print()


# ==============================================================================
# 7. ROOT STRUCTURE OF ACTUAL QUADRATICS
# ==============================================================================

def root_audit():
    print("=" * 78)
    print("7. ROOT STRUCTURE OF ACTUAL QUADRATICS")
    print("=" * 78)

    L = sp.symbols(
        "L"
    )

    for k in sorted(DATA):

        points = sorted(
            DATA[k]
        )

        if len(points) < 3:
            continue

        samples = [
            (
                sp.Integer(x),
                observed_delta(
                    k,
                    x
                )
            )
            for x in points
        ]

        q = simp(
            sp.interpolate(
                samples,
                L
            )
        )

        print(
            f"k={k}"
        )
        print(
            f"  polynomial = {sp.factor(q)}"
        )
        print(
            f"  roots = {sp.solve(q, L)}"
        )
        print()


# ==============================================================================
# 8. WHAT FRESH DATA ARE NEEDED?
# ==============================================================================

def fresh_data_requirements():
    print("=" * 78)
    print("8. REQUIRED FRESH DATA")
    print("=" * 78)
    print()

    print(
        "To resolve k=13, two additional exact L3 coefficients"
    )
    print(
        "are required at new ell-values."
    )
    print()

    print(
        "Recommended:"
    )
    print(
        "  k=13, ell=29, a=15"
    )
    print(
        "  k=13, ell=31, a=15"
    )
    print()
    print(
        "Better:"
    )
    print(
        "  k=13, ell=29,31,33"
    )
    print(
        "so the quadratic second difference can be tested."
    )
    print()

    print(
        "Also obtain:"
    )
    print(
        "  k=11, ell=27,29,31, a=13"
    )
    print(
        "to determine whether the transition begins at k=13"
    )
    print(
        "or whether the apparent anomaly is caused by insufficient"
    )
    print(
        "data/extraction."
    )
    print()


# ==============================================================================
# 9. FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The finite-difference calculation in Experiment 235"
    )
    print(
        "used the wrong expected value."
    )
    print(
        "For step size 2, the correct quadratic second difference is"
    )
    print(
        "  8A = 4(k+2)."
    )
    print()

    print(
        "Thus the observed k=3,5,7,9 sequences are"
    )
    print(
        "fully consistent with the candidate quadratic."
    )
    print()

    print(
        "The genuine unresolved issue is k=13."
    )
    print(
        "Only two k=13 points are available, so no quadratic"
    )
    print(
        "boundary law can yet be established there."
    )
    print()

    print(
        "Do not infer a new formula from these two points."
    )
    print(
        "Obtain fresh exact kernel data first."
    )
    print()

    print(
        "Do NOT start L4."
    )
    print(
        "Resolve the a=k+2 boundary mechanism completely."
    )
    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 236")
    print("k=13 a=k+2 BOUNDARY ANOMALY / CONSISTENCY AUDIT")
    print("=" * 78)
    print()

    observed_audit()

    finite_difference_audit()

    candidate_audit()

    interpolation_audit()

    k13_audit()

    same_leading_coefficient_test()

    root_audit()

    fresh_data_requirements()

    final_diagnostic()


if __name__ == "__main__":
    main()

