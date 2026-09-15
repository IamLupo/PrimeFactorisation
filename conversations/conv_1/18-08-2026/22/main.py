import sympy as sp


# ==============================================================================
# EXPERIMENT 235
# a=k+2 BOUNDARY-LAW STRESS TEST
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No exact_F required
#
# Purpose:
#   * validate the current a=k+2 candidate where data support it;
#   * explicitly expose the k=13 contradiction;
#   * separate interpolation from genuine structural evidence;
#   * test the finite-difference law correctly;
#   * inspect whether the boundary correction can be written in a
#     higher-level binomial form.
#
# IMPORTANT:
#   The candidate
#
#       Delta_(k+2)(L)
#         = (L-k-4)*((k+2)L-k(k+3))/2
#
# is NOT assumed to be universally true.
#
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

    return simp(
        B(K + 3, A)
        / sp.Integer(6)
        * (X - A - 1)
        * (X - A - 2)
        * ((K + 3) * X - K * A)
        / (K + 3)
    )


# ==============================================================================
# EXACT OBSERVED a=k+2 DATA
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
# CURRENT CANDIDATE
# ==============================================================================

def candidate_delta_k2(k, L):
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
# ACTUAL CORRECTION
# ==============================================================================

def actual_delta_k2(k, L):
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
# 1. DIRECT CANDIDATE AUDIT
# ==============================================================================

def direct_audit():
    print("=" * 78)
    print("1. DIRECT a=k+2 CANDIDATE AUDIT")
    print("=" * 78)

    tested = 0
    failures = 0

    for k in sorted(DATA):

        print(
            f"k={k}"
        )

        for L in sorted(DATA[k]):

            actual = actual_delta_k2(
                k,
                L
            )

            expected = candidate_delta_k2(
                k,
                L
            )

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
# 2. SYMBOLIC COEFFICIENTS
# ==============================================================================

def symbolic_candidate():
    print("=" * 78)
    print("2. SYMBOLIC CURRENT CANDIDATE")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    D = candidate_delta_k2(
        K,
        L
    )

    P = sp.Poly(
        sp.expand(D),
        L
    )

    A = simp(
        P.coeff_monomial(
            L**2
        )
    )

    Bc = simp(
        P.coeff_monomial(
            L
        )
    )

    C = simp(
        P.coeff_monomial(
            1
        )
    )

    print(
        "Delta_(k+2)(L) ="
    )
    print(
        f"  {sp.factor(D)}"
    )
    print()

    print(
        "A(k) ="
    )
    print(
        f"  {A}"
    )

    print(
        "B(k) ="
    )
    print(
        f"  {Bc}"
    )

    print(
        "C(k) ="
    )
    print(
        f"  {C}"
    )
    print()


# ==============================================================================
# 3. CORRECT FINITE-DIFFERENCE AUDIT
# ==============================================================================

def finite_difference_audit():
    print("=" * 78)
    print("3. EXACT FINITE-DIFFERENCE AUDIT")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    D = candidate_delta_k2(
        K,
        L
    )

    P = sp.Poly(
        sp.expand(D),
        L
    )

    A = simp(
        P.coeff_monomial(
            L**2
        )
    )

    symbolic_second = simp(
        sp.diff(
            D,
            L,
            2
        )
    )

    print(
        "symbolic second derivative ="
    )
    print(
        f"  {symbolic_second}"
    )

    print()

    print(
        "For ell-step h=2:"
    )
    print(
        "  Delta^2_h = 2*h^2*A = 8A"
    )
    print(
        "  equivalently Delta^2_h = 2*(k+2)"
    )
    print()

    expected_symbolic = simp(
        2 * (K + 2)
    )

    print(
        "expected step-2 second difference ="
    )
    print(
        f"  {expected_symbolic}"
    )
    print()

    for k in sorted(DATA):

        values = [
            actual_delta_k2(
                k,
                L
            )
            for L in sorted(DATA[k])
        ]

        if len(values) < 3:
            continue

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

        expected = simp(
            expected_symbolic.subs(
                K,
                k
            )
        )

        status = all(
            x == expected
            for x in d2
        )

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
            f"  expected = {expected}"
        )
        print(
            f"  {('PASS' if status else 'FAIL')}"
        )
        print()

    print()


# ==============================================================================
# 4. K=13 CONTRADICTION AUDIT
# ==============================================================================

def k13_audit():
    print("=" * 78)
    print("4. k=13 CONTRADICTION AUDIT")
    print("=" * 78)

    k = 13

    for L in sorted(DATA[k]):

        actual = actual_delta_k2(
            k,
            L
        )

        candidate = candidate_delta_k2(
            k,
            L
        )

        residual = simp(
            actual - candidate
        )

        print(
            f"ell={L:2d} "
            f"actual={actual} "
            f"candidate={candidate} "
            f"residual={residual}"
        )

    print()

    vals = [
        actual_delta_k2(
            k,
            L
        )
        for L in sorted(DATA[k])
    ]

    print(
        "actual correction values:"
    )
    print(
        f"  {vals}"
    )

    if len(vals) >= 3:

        d1 = [
            simp(
                vals[i + 1] - vals[i]
            )
            for i in range(
                len(vals) - 1
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
            "first differences:"
        )
        print(
            f"  {d1}"
        )

        print(
            "second differences:"
        )
        print(
            f"  {d2}"
        )

    print()


# ==============================================================================
# 5. INDEPENDENT QUADRATIC INTERPOLATION
# ==============================================================================

def interpolation_audit():
    print("=" * 78)
    print("5. INDEPENDENT QUADRATIC INTERPOLATION")
    print("=" * 78)

    L = sp.symbols(
        "L"
    )

    for k in sorted(DATA):

        pairs = sorted(
            DATA[k].items()
        )

        if len(pairs) < 3:
            print(
                f"k={k}: fewer than 3 points; no quadratic fit"
            )
            continue

        points = [
            (
                sp.Integer(x),
                actual_delta_k2(
                    k,
                    x
                )
            )
            for x, _ in pairs
        ]

        q = sp.interpolate(
            points,
            L
        )

        q = simp(q)

        residuals = []

        for x, y in points:

            residuals.append(
                simp(
                    q.subs(
                        L,
                        x
                    ) - y
                )
            )

        print(
            f"k={k}"
        )
        print(
            f"  interpolated = {sp.factor(q)}"
        )
        print(
            f"  residuals    = {residuals}"
        )

        if len(points) >= 4:

            # Leave-one-out structural test.
            loo_ok = True

            for omit in range(
                len(points)
            ):

                reduced = (
                    points[:omit]
                    + points[omit + 1:]
                )

                if len(reduced) < 3:
                    continue

                q2 = simp(
                    sp.interpolate(
                        reduced,
                        L
                    )
                )

                x, y = points[omit]

                if simp(
                    q2.subs(L, x) - y
                ) != 0:
                    loo_ok = False

            print(
                f"  leave-one-out cubic/quadratic consistency = "
                f"{'PASS' if loo_ok else 'FAIL'}"
            )

        print()


# ==============================================================================
# 6. COMPARE k=3..13 FITS
# ==============================================================================

def cross_k_fit():
    print("=" * 78)
    print("6. CROSS-k QUADRATIC COEFFICIENTS")
    print("=" * 78)

    L = sp.symbols(
        "L"
    )

    print(
        "   k          A             B             C"
    )

    for k in sorted(DATA):

        pairs = sorted(
            DATA[k].items()
        )

        if len(pairs) < 3:
            print(
                f"{k:4d}   insufficient data"
            )
            continue

        points = [
            (
                sp.Integer(x),
                actual_delta_k2(
                    k,
                    x
                )
            )
            for x, _ in pairs
        ]

        q = sp.Poly(
            sp.expand(
                sp.interpolate(
                    points,
                    L
                )
            ),
            L
        )

        A = simp(
            q.coeff_monomial(
                L**2
            )
        )

        Bc = simp(
            q.coeff_monomial(
                L
            )
        )

        C = simp(
            q.coeff_monomial(
                1
            )
        )

        print(
            f"{k:4d} "
            f"{str(A):>12s} "
            f"{str(Bc):>12s} "
            f"{str(C):>12s}"
        )

    print()


# ==============================================================================
# 7. STRUCTURAL ROOT AUDIT OF ACTUAL INTERPOLANTS
# ==============================================================================

def interpolation_root_audit():
    print("=" * 78)
    print("7. ROOT AUDIT OF ACTUAL INTERPOLANTS")
    print("=" * 78)

    L = sp.symbols(
        "L"
    )

    for k in sorted(DATA):

        pairs = sorted(
            DATA[k].items()
        )

        if len(pairs) < 3:
            continue

        points = [
            (
                sp.Integer(x),
                actual_delta_k2(
                    k,
                    x
                )
            )
            for x, _ in pairs
        ]

        q = simp(
            sp.interpolate(
                points,
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
            f"  roots      = {sp.solve(q, L)}"
        )
        print()


# ==============================================================================
# 8. CURRENT CANDIDATE ROOTS
# ==============================================================================

def candidate_root_audit():
    print("=" * 78)
    print("8. CURRENT CANDIDATE ROOTS")
    print("=" * 78)

    K, L = sp.symbols(
        "K L"
    )

    D = candidate_delta_k2(
        K,
        L
    )

    print(
        f"candidate = {sp.factor(D)}"
    )
    print(
        f"roots = {sp.solve(D, L)}"
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
        "The principal result is that the current"
    )
    print(
        "a=k+2 formula cannot yet be declared universal."
    )
    print()

    print(
        "It exactly reproduces k=3,5,7,9,11,"
    )
    print(
        "but fails at k=13."
    )
    print()

    print(
        "Therefore k=13 is not a numerical nuisance."
    )
    print(
        "It is evidence that the apparent boundary"
    )
    print(
        "formula changes or that the underlying"
    )
    print(
        "a=k+2 correction has additional terms."
    )
    print()

    print(
        "The next research target is:"
    )
    print()
    print(
        "  derive Delta_(k+2) directly from the exact"
    )
    print(
        "  pq kernel, rather than interpolate in ell."
    )
    print()

    print(
        "The finite-difference law must also be checked"
    )
    print(
        "against the actual data, not only against the"
    )
    print(
        "candidate polynomial."
    )
    print()

    print(
        "Do NOT start L4."
    )
    print(
        "First resolve the k=13 contradiction."
    )
    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 235")
    print("a=k+2 BOUNDARY-LAW STRESS TEST")
    print("=" * 78)
    print()

    direct_audit()

    symbolic_candidate()

    finite_difference_audit()

    k13_audit()

    interpolation_audit()

    cross_k_fit()

    interpolation_root_audit()

    candidate_root_audit()

    final_diagnostic()


if __name__ == "__main__":
    main()