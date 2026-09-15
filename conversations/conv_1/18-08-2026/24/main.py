import sympy as sp


# ==============================================================================
# EXPERIMENT 237
# L2/L3 BOUNDARY CONSISTENCY AUDIT
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No exact_F required
#
# Purpose:
#
#   The k=13 anomaly in Experiment 236 may originate in the L2 peel,
#   rather than in the genuine L3 boundary.
#
#   This experiment therefore:
#
#     1. records the established L2 boundary formulas;
#     2. records the established L3 interior formula;
#     3. reconstructs the exact numerical L2/L3 boundary data already known;
#     4. checks whether the L3 anomaly can be represented as an L2
#        boundary-propagation error;
#     5. compares the anomaly at equal d = ell-k across k;
#     6. checks whether the anomaly is concentrated at k=13;
#     7. derives the minimum fresh data required before another L3 law
#        is accepted.
#
# No new boundary formula is accepted.
# No L4 analysis is performed.
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
# ESTABLISHED INTERIOR L2 LAW
# ==============================================================================

def L2_interior(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        -C(K + 2, A)
        * X**2
        / 2
        + C(K + 2, A)
        * (
            2 * (K + 1) * A
            + K + 2
        )
        * X
        / (
            2 * (K + 2)
        )
        - K * A * (A + 1)
        * C(K + 2, A)
        / (
            2 * (K + 2)
        )
    )


# ==============================================================================
# ESTABLISHED L2 BOUNDARY CORRECTIONS
# ==============================================================================

def Delta_L2_k(k):
    K = sp.sympify(k)

    return simp(
        C(K + 2, 3)
    )


def Delta_L2_k1(k, L):
    K = sp.sympify(k)
    X = sp.sympify(L)

    return simp(
        -C(K + 2, 2) * X
        + K * (K + 2) * (K + 3)
        / 2
    )


def L2_boundary_k(k, L):
    return simp(
        L2_interior(k, k, L)
        + Delta_L2_k(k)
    )


def L2_boundary_k1(k, L):
    return simp(
        L2_interior(k, k + 1, L)
        + Delta_L2_k1(k, L)
    )


# ==============================================================================
# ESTABLISHED INTERIOR L3 LAW
# ==============================================================================

def L3_interior(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        C(K + 3, A)
        / 6
        * (X - A - 1)
        * (X - A - 2)
        * (
            (K + 3) * X
            - K * A
        )
        / (K + 3)
    )


# ==============================================================================
# ESTABLISHED L3 a=k CORRECTION
# ==============================================================================

def Delta_L3_k(k):
    K = sp.sympify(k)

    return simp(
        C(K + 2, 3)
    )


def L3_boundary_k(k, L):
    return simp(
        L3_interior(k, k, L)
        + Delta_L3_k(k)
    )


# ==============================================================================
# OBSERVED L3 a=k+2 DATA
#
# These are the exact values appearing in the previous experiments.
# ==============================================================================

L3_K2_DATA = {
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
# L3 INTERIOR EXTRAPOLATION AND OBSERVED CORRECTION
# ==============================================================================

def L3_observed_delta_k2(k, L):
    actual = sp.Integer(
        L3_K2_DATA[k][L]
    )

    interior = L3_interior(
        k,
        k + 2,
        L
    )

    return simp(
        actual - interior
    )


# ==============================================================================
# 1. L2 BOUNDARY FORMULAS
# ==============================================================================

def l2_boundary_audit():
    print("=" * 78)
    print("1. ESTABLISHED L2 BOUNDARY FORMULAS")
    print("=" * 78)
    print()

    K, L = sp.symbols(
        "K L"
    )

    print(
        "a=k:"
    )
    print(
        "  interior ="
    )
    print(
        "   ",
        sp.factor(
            L2_interior(K, K, L)
        )
    )
    print(
        "  correction =",
        Delta_L2_k(K)
    )
    print(
        "  corrected ="
    )
    print(
        "   ",
        sp.factor(
            L2_boundary_k(K, L)
        )
    )
    print()

    print(
        "a=k+1:"
    )
    print(
        "  interior ="
    )
    print(
        "   ",
        sp.factor(
            L2_interior(K, K + 1, L)
        )
    )
    print(
        "  correction ="
    )
    print(
        "   ",
        sp.factor(
            Delta_L2_k1(K, L)
        )
    )
    print(
        "  corrected ="
    )
    print(
        "   ",
        sp.factor(
            L2_boundary_k1(K, L)
        )
    )
    print()


# ==============================================================================
# 2. L3 k+2 OBSERVED CORRECTION
# ==============================================================================

def l3_k2_audit():
    print("=" * 78)
    print("2. OBSERVED L3 a=k+2 CORRECTIONS")
    print("=" * 78)
    print()

    for k in sorted(L3_K2_DATA):

        print(
            f"k={k}"
        )

        for L in sorted(
            L3_K2_DATA[k]
        ):

            actual = sp.Integer(
                L3_K2_DATA[k][L]
            )

            interior = L3_interior(
                k,
                k + 2,
                L
            )

            delta = simp(
                actual - interior
            )

            print(
                f"  ell={L:2d} "
                f"actual={str(actual):>8s} "
                f"interior={str(interior):>8s} "
                f"Delta={str(delta):>8s}"
            )

        print()


# ==============================================================================
# 3. ANOMALY RELATIVE TO THE LOW-k QUADRATIC
# ==============================================================================

def low_k_quadratic(k, L):
    K = sp.sympify(k)
    X = sp.sympify(L)

    return simp(
        (X - K - 4)
        * (
            (K + 2) * X
            - K * (K + 3)
        )
        / 2
    )


def anomaly_audit():
    print("=" * 78)
    print("3. DEVIATION FROM THE LOW-k a=k+2 LAW")
    print("=" * 78)
    print()

    for k in sorted(
        L3_K2_DATA
    ):

        print(
            f"k={k}"
        )

        for L in sorted(
            L3_K2_DATA[k]
        ):

            observed = L3_observed_delta_k2(
                k,
                L
            )

            candidate = low_k_quadratic(
                k,
                L
            )

            anomaly = simp(
                observed - candidate
            )

            print(
                f"  ell={L:2d} "
                f"observed={str(observed):>8s} "
                f"candidate={str(candidate):>8s} "
                f"anomaly={str(anomaly):>8s}"
            )

        print()


# ==============================================================================
# 4. SAME d = ell-k COMPARISON
# ==============================================================================

def fixed_d_audit():
    print("=" * 78)
    print("4. SAME d = ell-k CROSS-k AUDIT")
    print("=" * 78)
    print()

    groups = {}

    for k in sorted(
        L3_K2_DATA
    ):

        for L in sorted(
            L3_K2_DATA[k]
        ):

            d = L - k

            groups.setdefault(
                d,
                []
            ).append(
                (
                    k,
                    L,
                    L3_observed_delta_k2(
                        k,
                        L
                    )
                )
            )

    for d in sorted(groups):

        entries = groups[d]

        if len(entries) < 2:
            continue

        print(
            f"d={d}"
        )

        for k, L, value in entries:

            print(
                f"  k={k:2d} "
                f"ell={L:2d} "
                f"Delta={value}"
            )

        print()


# ==============================================================================
# 5. CHECK WHETHER k=13 IS AN ISOLATED OUTLIER
# ==============================================================================

def k13_isolation_audit():
    print("=" * 78)
    print("5. k=13 ISOLATION AUDIT")
    print("=" * 78)
    print()

    for L in sorted(
        L3_K2_DATA[13]
    ):

        observed = L3_observed_delta_k2(
            13,
            L
        )

        candidate = low_k_quadratic(
            13,
            L
        )

        anomaly = simp(
            observed - candidate
        )

        print(
            f"k=13 ell={L}"
        )
        print(
            f"  observed Delta = {observed}"
        )
        print(
            f"  low-k candidate = {candidate}"
        )
        print(
            f"  anomaly = {anomaly}"
        )
        print()


# ==============================================================================
# 6. L2 ERROR PROPAGATION TEST
#
# If an incorrect L2 boundary correction E(k,L) had been used while
# extracting L3, its homogeneous subtraction would appear with exactly
# the same degree and monomial support in L3.
#
# This section therefore tests the simplest possibilities:
#
#   E = constant
#   E = linear in L
#   E = quadratic in L
#
# against the observed k=13 L3 anomaly.
#
# This does NOT assert that any such error exists.
# ==============================================================================

def l2_error_propagation_audit():
    print("=" * 78)
    print("6. L2-ERROR PROPAGATION DIAGNOSTIC")
    print("=" * 78)
    print()

    K, L = sp.symbols(
        "K L"
    )

    k = 13

    observed_points = []

    for ell in sorted(
        L3_K2_DATA[k]
    ):

        anomaly = simp(
            L3_observed_delta_k2(
                k,
                ell
            )
            - low_k_quadratic(
                k,
                ell
            )
        )

        observed_points.append(
            (
                sp.Integer(ell),
                anomaly
            )
        )

    print(
        "k=13 anomaly sequence:"
    )
    print(
        observed_points
    )
    print()

    if len(observed_points) >= 3:

        x = sp.symbols(
            "x"
        )

        q = simp(
            sp.interpolate(
                observed_points,
                x
            )
        )

        print(
            "Interpolating anomaly:"
        )
        print(
            "  ",
            sp.factor(q)
        )
        print()

        degree = sp.Poly(
            q,
            x
        ).degree()

        print(
            "anomaly degree =",
            degree
        )

        if degree <= 2:
            print(
                "The k=13 anomaly is compatible with a quadratic "
                "propagation error."
            )
        else:
            print(
                "The k=13 anomaly is not quadratic in the available data."
            )

    else:

        print(
            "Only two anomaly points are currently available."
        )
        print(
            "No polynomial degree can be inferred."
        )

    print()


# ==============================================================================
# 7. FRESH DATA THAT SHOULD BE GENERATED
# ==============================================================================

def fresh_data_plan():
    print("=" * 78)
    print("7. REQUIRED FRESH DATA")
    print("=" * 78)
    print()

    print(
        "The next kernel computation should produce, at minimum:"
    )
    print()

    print(
        "  k=13, ell=29, a=15"
    )
    print(
        "  k=13, ell=31, a=15"
    )
    print(
        "  k=13, ell=33, a=15"
    )
    print()

    print(
        "But this time also extract the complete L2 rows:"
    )
    print(
        "  k=13, ell=29, a=13"
    )
    print(
        "  k=13, ell=29, a=14"
    )
    print(
        "  k=13, ell=29, a=15"
    )
    print()

    print(
        "Repeat for ell=31."
    )
    print()

    print(
        "This distinguishes:"
    )
    print(
        "  (A) a genuine L3 boundary transition"
    )
    print(
        "from"
    )
    print(
        "  (B) an L2 peeling error propagating into L3."
    )
    print()


# ==============================================================================
# 8. FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The low-k a=k+2 quadratic is strongly supported for"
    )
    print(
        "k=3,5,7,9."
    )
    print()

    print(
        "The k=13 discrepancy is real in the supplied data,"
    )
    print(
        "but it is not yet established as a new boundary law."
    )
    print()

    print(
        "The crucial unresolved question is:"
    )
    print(
        "  Does the anomaly already exist in the L2 boundary peel?"
    )
    print(
        "or"
    )
    print(
        "  Does it first appear in the genuine L3 coefficient?"
    )
    print()

    print(
        "Therefore the next exact-kernel run must extract L2 and L3"
    )
    print(
        "together for k=13."
    )
    print()

    print(
        "Do NOT start L4."
    )
    print(
        "Do NOT fit a new universal Delta_(k+2) law."
    )
    print(
        "First localize the k=13 anomaly."
    )
    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 237")
    print("L2/L3 BOUNDARY CONSISTENCY AUDIT")
    print("=" * 78)
    print()

    l2_boundary_audit()

    l3_k2_audit()

    anomaly_audit()

    fixed_d_audit()

    k13_isolation_audit()

    l2_error_propagation_audit()

    fresh_data_plan()

    final_diagnostic()


if __name__ == "__main__":
    main()

