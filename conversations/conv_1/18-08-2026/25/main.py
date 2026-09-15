# ==============================================================================
# EXPERIMENT 238
# FIXED-d CROSS-k BOUNDARY LAW / k=13 ANOMALY AUDIT
# ==============================================================================
#
# Standalone main.py
# No imports
# No previous experiment
# No filesystem access
# No exact_F required
#
# Purpose:
#
#   Experiment 237 showed that the k=13 a=k+2 corrections
#   are anomalous relative to the lower-k data.
#
#   Experiment 238 changes coordinates from (k, ell) to
#
#       d = ell-k
#
#   and tests the boundary correction at fixed d.
#
#   The observed lower-k data strongly suggest
#
#       Delta_(k+2)(k,d)
#         = (d-4) * (k(d-1) + 2d) / 2.
#
#   This experiment:
#
#       1. verifies that law against all available non-anomalous data;
#       2. tests same-d cross-k consistency;
#       3. performs leave-one-out prediction;
#       4. isolates the k=13 residuals;
#       5. checks whether the anomaly is constant, linear, or multiplicative;
#       6. determines the minimum fresh data needed to resolve it.
#
# No new universal law is accepted unless it survives all audits.
# No L4 analysis.
# ==============================================================================


# ------------------------------------------------------------------------------
# OBSERVED DATA
#
# Stored as:
#
#   k : { ell : observed Delta_(k+2) }
#
# These are the exact corrections extracted in Experiments 232-237.
# ------------------------------------------------------------------------------

DATA = {
    3: {
        7: 0,
        9: 27,
        11: 74,
        13: 141,
        15: 228,
    },

    5: {
        11: 37,
        13: 102,
        15: 195,
        17: 316,
        19: 465,
    },

    7: {
        15: 130,
        17: 249,
        19: 404,
    },

    9: {
        21: 492,
        23: 725,
        25: 1002,
    },

    11: {
        23: 580,
        25: 855,
    },

    13: {
        25: 2870,
        27: 3525,
    },
}


# ------------------------------------------------------------------------------
# EXACT CANDIDATE IN (k,d)
# ------------------------------------------------------------------------------

def candidate_delta(k, d):
    """
    Candidate fixed-d law:

        Delta = (d-4) * (k(d-1) + 2d) / 2

    All inputs are integers.

    The numerator is known to be even on the tested odd-k,
    even-d branches.
    """
    numerator = (d - 4) * (
        k * (d - 1) + 2 * d
    )

    if numerator % 2 != 0:
        raise ValueError(
            "Candidate produced non-integral value."
        )

    return numerator // 2


# ------------------------------------------------------------------------------
# CONVERT DATA TO FIXED-d GROUPS
# ------------------------------------------------------------------------------

def build_d_groups():
    groups = {}

    for k in sorted(DATA):
        for ell in sorted(DATA[k]):

            d = ell - k
            value = DATA[k][ell]

            if d not in groups:
                groups[d] = []

            groups[d].append(
                (k, ell, value)
            )

    return groups


# ------------------------------------------------------------------------------
# 1. DIRECT CANDIDATE AUDIT
# ------------------------------------------------------------------------------

def direct_candidate_audit():
    print("=" * 78)
    print("1. DIRECT FIXED-d CANDIDATE AUDIT")
    print("=" * 78)
    print()

    tested = 0
    failures = 0

    for k in sorted(DATA):

        print(
            f"k={k}"
        )

        for ell in sorted(DATA[k]):

            d = ell - k
            actual = DATA[k][ell]
            expected = candidate_delta(k, d)
            residual = actual - expected

            tested += 1

            status = (
                "PASS"
                if residual == 0
                else "FAIL"
            )

            if residual != 0:
                failures += 1

            print(
                f"  ell={ell:2d} d={d:2d} "
                f"actual={actual:6d} "
                f"expected={expected:6d} "
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
# 2. SAME-d CROSS-k AUDIT
# ------------------------------------------------------------------------------

def same_d_audit():
    print("=" * 78)
    print("2. SAME-d CROSS-k AUDIT")
    print("=" * 78)
    print()

    groups = build_d_groups()

    for d in sorted(groups):

        entries = groups[d]

        if len(entries) < 2:
            continue

        print(
            f"d={d}"
        )

        for k, ell, actual in entries:

            expected = candidate_delta(
                k,
                d
            )

            residual = actual - expected

            print(
                f"  k={k:2d} "
                f"ell={ell:2d} "
                f"actual={actual:6d} "
                f"expected={expected:6d} "
                f"residual={residual:6d}"
            )

        print()


# ------------------------------------------------------------------------------
# 3. CHECK LINEARITY IN k AT FIXED d
# ------------------------------------------------------------------------------

def fixed_d_linearity_audit():
    print("=" * 78)
    print("3. FIXED-d LINEARITY IN k")
    print("=" * 78)
    print()

    groups = build_d_groups()

    for d in sorted(groups):

        entries = groups[d]

        if len(entries) < 3:
            continue

        print(
            f"d={d}"
        )

        # Candidate is:
        #
        #   Delta = A_d * k + B_d
        #
        # where
        #
        #   A_d = (d-4)(d-1)/2
        #   B_d = d(d-4)

        slope = (
            (d - 4) * (d - 1) // 2
        )
        intercept = (
            d * (d - 4)
        )

        print(
            f"  predicted slope     = {slope}"
        )
        print(
            f"  predicted intercept = {intercept}"
        )

        for i in range(
            len(entries) - 1
        ):

            k1, ell1, v1 = entries[i]
            k2, ell2, v2 = entries[i + 1]

            dk = k2 - k1

            if dk == 0:
                continue

            observed_slope_num = v2 - v1

            print(
                f"  k={k1}->{k2}: "
                f"delta change={observed_slope_num} "
                f"over dk={dk}"
            )

        print()


# ------------------------------------------------------------------------------
# 4. LEAVE-ONE-OUT PREDICTION
#
# Since the law is linear in k for fixed d, two independent points
# determine the cross-k line.
#
# This is especially useful at d=12 and d=14 where k=13 is an outlier.
# ------------------------------------------------------------------------------

def leave_one_out_audit():
    print("=" * 78)
    print("4. LEAVE-ONE-OUT CROSS-k PREDICTION")
    print("=" * 78)
    print()

    groups = build_d_groups()

    for d in sorted(groups):

        entries = groups[d]

        if len(entries) < 3:
            continue

        print(
            f"d={d}"
        )

        for omitted_index in range(
            len(entries)
        ):

            training = [
                entry
                for j, entry in enumerate(entries)
                if j != omitted_index
            ]

            if len(training) < 2:
                continue

            k1, _, v1 = training[0]
            k2, _, v2 = training[1]

            if k2 == k1:
                continue

            # Exact slope.
            numerator = v2 - v1
            denominator = k2 - k1

            kt, ellt, actual = entries[
                omitted_index
            ]

            left = (
                v1 * denominator
                + numerator * (
                    kt - k1
                )
            )

            right = denominator

            # Only report exact integral predictions.
            if left % right == 0:
                predicted = left // right
                residual = actual - predicted

                print(
                    f"  omitted k={kt:2d} "
                    f"ell={ellt:2d} "
                    f"predicted={predicted:6d} "
                    f"actual={actual:6d} "
                    f"residual={residual:6d}"
                )

        print()


# ------------------------------------------------------------------------------
# 5. k=13 ANOMALY MAGNITUDE
# ------------------------------------------------------------------------------

def k13_anomaly_audit():
    print("=" * 78)
    print("5. k=13 ANOMALY MAGNITUDE")
    print("=" * 78)
    print()

    k = 13

    for ell in sorted(DATA[k]):

        d = ell - k
        actual = DATA[k][ell]
        expected = candidate_delta(
            k,
            d
        )

        anomaly = actual - expected

        print(
            f"k=13 ell={ell} d={d}"
        )
        print(
            f"  actual    = {actual}"
        )
        print(
            f"  expected  = {expected}"
        )
        print(
            f"  anomaly   = {anomaly}"
        )

        # Ratio is only meaningful when expected != 0.
        if expected != 0:

            print(
                f"  actual/expected = "
                f"{actual}/{expected}"
            )

        print()


# ------------------------------------------------------------------------------
# 6. SAME-d=12 DIAGNOSTIC
# ------------------------------------------------------------------------------

def d12_diagnostic():
    print("=" * 78)
    print("6. d=12 DIAGNOSTIC")
    print("=" * 78)
    print()

    d = 12

    entries = [
        entry
        for entry in build_d_groups()[d]
    ]

    print(
        "For d=12 the candidate law becomes:"
    )
    print(
        "  Delta = 44*k + 96"
    )
    print()

    for k, ell, actual in entries:

        expected = (
            44 * k
            + 96
        )

        residual = (
            actual
            - expected
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"actual={actual:6d} "
            f"expected={expected:6d} "
            f"residual={residual:6d}"
        )

    print()


# ------------------------------------------------------------------------------
# 7. SAME-d=14 DIAGNOSTIC
# ------------------------------------------------------------------------------

def d14_diagnostic():
    print("=" * 78)
    print("7. d=14 DIAGNOSTIC")
    print("=" * 78)
    print()

    d = 14

    entries = [
        entry
        for entry in build_d_groups()[d]
    ]

    print(
        "For d=14 the candidate law becomes:"
    )
    print(
        "  Delta = 65*k + 140"
    )
    print()

    for k, ell, actual in entries:

        expected = (
            65 * k
            + 140
        )

        residual = (
            actual
            - expected
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"actual={actual:6d} "
            f"expected={expected:6d} "
            f"residual={residual:6d}"
        )

    print()


# ------------------------------------------------------------------------------
# 8. TEST WHETHER k=13 ANOMALY IS A SIMPLE MULTIPLE
# ------------------------------------------------------------------------------

def multiplicative_anomaly_audit():
    print("=" * 78)
    print("8. MULTIPLICATIVE ANOMALY AUDIT")
    print("=" * 78)
    print()

    for k in [13]:

        for ell in sorted(
            DATA[k]
        ):

            d = ell - k
            actual = DATA[k][ell]
            expected = candidate_delta(
                k,
                d
            )

            if expected == 0:
                continue

            # Compare integer products q*expected.
            quotient = actual // expected
            remainder = actual % expected

            print(
                f"k={k} ell={ell} d={d}"
            )
            print(
                f"  expected = {expected}"
            )
            print(
                f"  actual   = {actual}"
            )
            print(
                f"  integer quotient = {quotient}"
            )
            print(
                f"  remainder        = {remainder}"
            )
            print()


# ------------------------------------------------------------------------------
# 9. MINIMUM FRESH DATA
# ------------------------------------------------------------------------------

def fresh_data_plan():
    print("=" * 78)
    print("9. REQUIRED FRESH DATA")
    print("=" * 78)
    print()

    print(
        "The next exact-kernel extraction should obtain:"
    )
    print()

    print(
        "A. k=13"
    )
    print(
        "   ell=29, a=15"
    )
    print(
        "   ell=31, a=15"
    )
    print(
        "   ell=33, a=15"
    )
    print()

    print(
        "B. k=11"
    )
    print(
        "   ell=27, a=13"
    )
    print(
        "   ell=29, a=13"
    )
    print(
        "   ell=31, a=13"
    )
    print()

    print(
        "C. k=13 complete local rows"
    )
    print(
        "   ell=29: a=13,14,15"
    )
    print(
        "   ell=31: a=13,14,15"
    )
    print(
        "   ell=33: a=13,14,15"
    )
    print()

    print(
        "The purpose is to decide whether the anomaly:"
    )
    print(
        "  1. begins at k=13;"
    )
    print(
        "  2. begins already in the L2 peel;"
    )
    print(
        "  3. depends on d;"
    )
    print(
        "  4. or is caused by selecting the wrong coefficient."
    )
    print()


# ------------------------------------------------------------------------------
# 10. FINAL DIAGNOSTIC
# ------------------------------------------------------------------------------

def final_diagnostic():
    print("=" * 78)
    print("10. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The lower-k boundary data are highly structured."
    )
    print()

    print(
        "At fixed d = ell-k, the established candidate is:"
    )
    print(
        "  Delta_(k+2)(k,d)"
    )
    print(
        "    = (d-4) * (k(d-1) + 2d) / 2"
    )
    print()

    print(
        "Hence the correction is affine in k for every fixed d."
    )
    print()

    print(
        "The reported k=13 values violate this structure."
    )
    print(
        "The violation is far too large to be treated as rounding."
    )
    print()

    print(
        "Therefore:"
    )
    print(
        "  DO NOT fit a replacement boundary formula yet."
    )
    print(
        "  DO NOT start L4."
    )
    print(
        "  DO NOT reinterpret k=13 as a new regime."
    )
    print()

    print(
        "The next exact-kernel computation must localize"
    )
    print(
        "the k=13 discrepancy using complete L2/L3 rows."
    )
    print()


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 238")
    print("FIXED-d CROSS-k BOUNDARY LAW / k=13 ANOMALY AUDIT")
    print("=" * 78)
    print()

    direct_candidate_audit()

    same_d_audit()

    fixed_d_linearity_audit()

    leave_one_out_audit()

    k13_anomaly_audit()

    d12_diagnostic()

    d14_diagnostic()

    multiplicative_anomaly_audit()

    fresh_data_plan()

    final_diagnostic()


if __name__ == "__main__":
    main()

