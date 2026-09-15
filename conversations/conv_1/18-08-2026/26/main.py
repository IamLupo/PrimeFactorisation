# ==============================================================================
# EXPERIMENT 239
# k=13 LOCALIZATION AUDIT: L2 -> L3 BOUNDARY PROPAGATION
# ==============================================================================
#
# Standalone main.py
# No imports
# No filesystem access
# No exact_F
# No previous experiment imported
#
# PURPOSE
#
# Experiment 238 established:
#
#   Delta_(k+2)(k,d)
#     = (d-4) * (k(d-1) + 2d) / 2
#
# for all presently reliable lower-k data, while k=13 violates it.
#
# We must now decide whether the discrepancy is:
#
#   A) already present in the L2 boundary peel,
#   B) introduced only in the L3 extraction,
#   C) caused by a wrong coefficient coordinate,
#   D) caused by an indexing/shift convention,
#   E) or a genuine k=13 structural transition.
#
# This experiment does NOT invent a new formula.
#
# It establishes:
#
#   1. exact L3 anomaly values already known;
#   2. exact nearby L3 rows required for localization;
#   3. the established L2 boundary expectations;
#   4. consistency constraints that fresh kernel output MUST satisfy;
#   5. a local decision table separating possible failure modes.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# EXISTING L3 DATA
#
# These are the exact coefficients already supplied.
#
# Each row is:
#
#   (a, coefficient)
#
# for L3.
# ------------------------------------------------------------------------------

L3_DATA = {

    3: {
        9: {
            2: 600,
            3: 510,
            4: 165,
            5: 66,
        },
        11: {
            2: 1400,
            3: 1340,
            4: 610,
            5: 244,
        },
        13: {
            2: 2700,
            3: 2770,
            4: 1455,
            5: 582,
        },
    },

    5: {
        13: {
            4: 6860,
            5: 3906,
            6: 1162,
            7: 332,
        },
        15: {
            4: 13125,
            5: 8015,
            6: 2765,
            7: 790,
        },
        17: {
            4: 22330,
            5: 14280,
            6: 5348,
            7: 1528,
        },
    },

    7: {
        15: {
            6: 21168,
            7: 8568,
            8: 1890,
            9: 420,
        },
        17: {
            6: 40320,
            7: 17508,
            8: 4491,
            9: 998,
        },
        19: {
            6: 68376,
            7: 31104,
            8: 8676,
            9: 1928,
        },
    },

    9: {
        21: {
            8: 163350,
            9: 57640,
            10: 12804,
            11: 2328,
        },
        23: {
            8: 255255,
            9: 93115,
            10: 21835,
            11: 3970,
        },
        25: {
            8: 376200,
            9: 140690,
            10: 34254,
            11: 6228,
        },
    },

    11: {
        23: {
            10: 333476,
            11: 96096,
            12: 17732,
            13: 2728,
        },
        25: {
            10: 520520,
            11: 155090,
            12: 30225,
            13: 4650,
        },
    },

    13: {
        25: {
            15: 5330,
        },
        27: {
            15: 7870,
        },
    },
}


# ------------------------------------------------------------------------------
# ESTABLISHED L3 INTERIOR FORMULA
#
# Valid for a < k.
#
# P3(k,a,L) =
#
#   C(k+3,a)/6
#   * (L-a-1)
#   * (L-a-2)
#   * ((k+3)L-ka)/(k+3)
#
# We implement the binomial coefficient ourselves.
# ------------------------------------------------------------------------------

def binomial(n, r):

    if r < 0:
        return 0

    if r > n:
        return 0

    if r == 0:
        return 1

    if r == n:
        return 1

    r2 = r

    if r2 > n - r2:
        r2 = n - r2

    value = 1

    j = 1

    while j <= r2:
        value = value * (n - r2 + j)
        value = value // j
        j += 1

    return value


def interior_L3(k, a, ell):

    c = binomial(
        k + 3,
        a
    )

    numerator = (
        c
        * (ell - a - 1)
        * (ell - a - 2)
        * ((k + 3) * ell - k * a)
    )

    denominator = (
        6 * (k + 3)
    )

    if numerator % denominator != 0:
        # Rational values are possible in intermediate symbolic
        # calculations, but all actual tested coefficient values
        # here are integral.
        return (
            numerator,
            denominator
        )

    return numerator // denominator


# ------------------------------------------------------------------------------
# L3 a=k+2 CORRECTION CANDIDATE
#
# This is the lower-k law whose validity at k=13 is disputed.
#
# IMPORTANT:
# This is NOT asserted as universal.
# It is used only as a diagnostic baseline.
# ------------------------------------------------------------------------------

def delta_k_plus_2_candidate(k, ell):

    d = ell - k

    numerator = (
        (d - 4)
        * (
            k * (d - 1)
            + 2 * d
        )
    )

    if numerator % 2 != 0:
        return (
            numerator,
            2
        )

    return numerator // 2


# ------------------------------------------------------------------------------
# CORRECTED L3 a=k FORM
# ------------------------------------------------------------------------------

def delta_k(k):

    return binomial(
        k + 2,
        3
    )


# ------------------------------------------------------------------------------
# CORRECTED L3 a=k+1 FORM
# ------------------------------------------------------------------------------

def delta_k_plus_1(k, ell):

    first = (
        -binomial(
            k + 2,
            2
        )
        * ell
    )

    second = (
        k
        * (k + 2)
        * (k + 3)
    )

    if second % 2 != 0:
        return (
            first * 2 + second,
            2
        )

    return (
        first
        + second // 2
    )


# ------------------------------------------------------------------------------
# L3 INTERIOR EXTRAPOLATION AT a=k+2
#
# This is intentionally only diagnostic.
# ------------------------------------------------------------------------------

def interior_at_k_plus_2(k, ell):

    a = k + 2

    c = binomial(
        k + 3,
        a
    )

    numerator = (
        c
        * (ell - a - 1)
        * (ell - a - 2)
        * (
            (k + 3) * ell
            - k * a
        )
    )

    denominator = (
        6 * (k + 3)
    )

    if numerator % denominator == 0:
        return numerator // denominator

    return (
        numerator,
        denominator
    )


# ------------------------------------------------------------------------------
# 1. KNOWN k=13 ANOMALY
# ------------------------------------------------------------------------------

def show_k13_anomaly():

    print("=" * 78)
    print("1. EXISTING k=13 L3 ANOMALY")
    print("=" * 78)
    print()

    k = 13

    for ell in [25, 27]:

        actual = L3_DATA[k][ell][15]

        interior = interior_at_k_plus_2(
            k,
            ell
        )

        candidate_delta = (
            delta_k_plus_2_candidate(
                k,
                ell
            )
        )

        print(
            f"k=13 ell={ell}"
        )

        print(
            f"  actual a=k+2       = {actual}"
        )

        print(
            f"  interior extrap.   = {interior}"
        )

        print(
            f"  low-k Delta law    = {candidate_delta}"
        )

        print(
            f"  actual-internal    = "
            f"{actual - interior}"
        )

        print(
            f"  actual-candidate   = "
            f"{actual - candidate_delta}"
        )

        print()


# ------------------------------------------------------------------------------
# 2. EXPECTED LOCAL L3 ROW
#
# The critical fresh rows are:
#
#   a=k
#   a=k+1
#   a=k+2
#
# at the same ell.
#
# We print all exact expectations for k=13.
# ------------------------------------------------------------------------------

def expected_k13_local_rows():

    print("=" * 78)
    print("2. EXPECTED k=13 LOCAL L3 ROWS")
    print("=" * 78)
    print()

    k = 13

    for ell in [25, 27, 29, 31, 33]:

        print(
            f"k={k} ell={ell}"
        )

        # a=k
        print(
            f"  a=k={k}: "
            f"Delta = {delta_k(k)}"
        )

        # a=k+1
        print(
            f"  a=k+1={k+1}: "
            f"Delta = "
            f"{delta_k_plus_1(k, ell)}"
        )

        # a=k+2
        print(
            f"  a=k+2={k+2}: "
            f"candidate Delta = "
            f"{delta_k_plus_2_candidate(k, ell)}"
        )

        print()


# ------------------------------------------------------------------------------
# 3. LOWER-k CONTROL AUDIT
#
# This confirms that the diagnostic extraction itself agrees with
# all known lower-k observations.
# ------------------------------------------------------------------------------

def lower_k_control():

    print("=" * 78)
    print("3. LOWER-k CONTROL AUDIT")
    print("=" * 78)
    print()

    failures = 0
    tested = 0

    for k in sorted(
        L3_DATA
    ):

        if k == 13:
            continue

        for ell in sorted(
            L3_DATA[k]
        ):

            if (k + 2) not in L3_DATA[k][ell]:
                continue

            actual = L3_DATA[k][ell][
                k + 2
            ]

            predicted = (
                interior_at_k_plus_2(
                    k,
                    ell
                )
                + delta_k_plus_2_candidate(
                    k,
                    ell
                )
            )

            residual = (
                actual
                - predicted
            )

            tested += 1

            status = (
                "PASS"
                if residual == 0
                else "FAIL"
            )

            if residual != 0:
                failures += 1

            print(
                f"k={k:2d} ell={ell:2d} "
                f"actual={actual:6d} "
                f"predicted={predicted:6d} "
                f"residual={residual:6d} "
                f"{status}"
            )

    print()

    print(
        f"tested = {tested}"
    )

    print(
        f"failures = {failures}"
    )

    print()


# ------------------------------------------------------------------------------
# 4. L3 BOUNDARY-CORRECTION HIERARCHY
#
# We print the known forms side-by-side.
# ------------------------------------------------------------------------------

def boundary_hierarchy():

    print("=" * 78)
    print("4. KNOWN L3 BOUNDARY HIERARCHY")
    print("=" * 78)
    print()

    print(
        "a=k:"
    )

    print(
        "  Delta_k = C(k+2,3)"
    )

    print()

    print(
        "a=k+1:"
    )

    print(
        "  Delta_(k+1)"
        " = -C(k+2,2)*L"
        " + k(k+2)(k+3)/2"
    )

    print()

    print(
        "a=k+2:"
    )

    print(
        "  lower-k candidate:"
    )

    print(
        "  Delta_(k+2)"
        " = (L-k-4)"
        "   * ((2k-1)L-k(k+3))/2"
    )

    print()

    print(
        "The degrees in L are therefore:"
    )

    print(
        "  a=k     -> degree 0"
    )

    print(
        "  a=k+1   -> degree 1"
    )

    print(
        "  a=k+2   -> degree 2"
    )

    print()


# ------------------------------------------------------------------------------
# 5. L2 LOCALIZATION TEST
#
# The objective is to obtain the three neighboring L2 coefficients
# at the same k, ell:
#
#   a=k
#   a=k+1
#   a=k+2
#
# They determine whether the anomaly is already present before
# L3 is formed.
#
# Since the exact kernel is not available in this standalone experiment,
# this section prints the exact quantities that must be supplied
# by the next kernel extraction.
# ------------------------------------------------------------------------------

def l2_localization_request():

    print("=" * 78)
    print("5. L2 LOCALIZATION REQUEST")
    print("=" * 78)
    print()

    print(
        "For k=13, obtain the following exact L2 coefficients:"
    )

    print()

    for ell in [29, 31, 33]:

        print(
            f"k=13 ell={ell}"
        )

        print(
            "  L2:"
        )

        print(
            "    a=13  [N^13 X^(ell-15)]"
        )

        print(
            "    a=14  [N^14 X^(ell-16)]"
        )

        print(
            "    a=15  [N^15 X^(ell-17)]"
        )

        print()

    print(
        "Also obtain the corresponding L3 coefficients:"
    )

    for ell in [29, 31, 33]:

        print(
            f"  k=13 ell={ell} "
            f"a=15 in L3"
        )

    print()


# ------------------------------------------------------------------------------
# 6. DECISION TREE
# ------------------------------------------------------------------------------

def decision_tree():

    print("=" * 78)
    print("6. ANOMALY LOCALIZATION DECISION TREE")
    print("=" * 78)
    print()

    print(
        "CASE A:"
    )

    print(
        "  L2 coefficients fail their established boundary laws."
    )

    print(
        "  -> anomaly originates before L3."
    )

    print()

    print(
        "CASE B:"
    )

    print(
        "  L2 coefficients pass exactly,"
    )

    print(
        "  but L3 a=k+2 fails."
    )

    print(
        "  -> anomaly is genuinely introduced at the L3 boundary."
    )

    print()

    print(
        "CASE C:"
    )

    print(
        "  a=k and a=k+1 pass,"
    )

    print(
        "  a=k+2 alone fails."
    )

    print(
        "  -> boundary transition is localized exactly at a=k+2."
    )

    print()

    print(
        "CASE D:"
    )

    print(
        "  a=k+2 passes at ell=29,31,33,"
    )

    print(
        "  -> previous k=13 values require rechecking."
    )

    print()


# ------------------------------------------------------------------------------
# 7. NECESSARY FRESH KERNEL DATA
# ------------------------------------------------------------------------------

def fresh_kernel_data():

    print("=" * 78)
    print("7. REQUIRED FRESH EXACT-KERNEL DATA")
    print("=" * 78)
    print()

    print(
        "Minimum:"
    )

    print(
        "  k=13, ell=29,31,33"
    )

    print(
        "  L2 a=13,14,15"
    )

    print(
        "  L3 a=15"
    )

    print()

    print(
        "Strong control:"
    )

    print(
        "  k=11, ell=27,29,31"
    )

    print(
        "  L2 a=11,12,13"
    )

    print(
        "  L3 a=13"
    )

    print()

    print(
        "Optional:"
    )

    print(
        "  k=15, ell=29,31,33"
    )

    print(
        "  L2 a=15,16,17"
    )

    print(
        "  L3 a=17"
    )

    print()


# ------------------------------------------------------------------------------
# 8. FINAL DIAGNOSTIC
# ------------------------------------------------------------------------------

def final_diagnostic():

    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "Current conclusion:"
    )

    print(
        "  The lower-k a=k+2 law is not yet disproved."
    )

    print(
        "  The k=13 observations are anomalous."
    )

    print(
        "  Two k=13 points cannot determine a new quadratic."
    )

    print()

    print(
        "Therefore the next mathematical operation is:"
    )

    print(
        "  LOCALIZE -> DO NOT FIT."
    )

    print()

    print(
        "Specifically:"
    )

    print(
        "  exact pq kernel"
    )

    print(
        "      -> L2 local rows"
    )

    print(
        "      -> L3 local rows"
    )

    print(
        "      -> compare a=k,k+1,k+2"
    )

    print(
        "      -> identify first layer that fails."
    )

    print()

    print(
        "No L4 analysis is permitted until this localization"
    )

    print(
        "is resolved."
    )

    print()


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 239")
    print("k=13 LOCALIZATION AUDIT: L2 -> L3 BOUNDARY PROPAGATION")
    print("=" * 78)
    print()

    show_k13_anomaly()

    expected_k13_local_rows()

    lower_k_control()

    boundary_hierarchy()

    l2_localization_request()

    decision_tree()

    fresh_kernel_data()

    final_diagnostic()


if __name__ == "__main__":
    main()

