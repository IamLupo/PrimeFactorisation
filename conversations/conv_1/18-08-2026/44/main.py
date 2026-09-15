"""
==============================================================================
EXPERIMENT 251
EXACT r=5 BOUNDARY BREAK-POINT / LOCAL LADDER FORENSICS
==============================================================================

Standalone.
Exact arithmetic over QQ.
No filesystem.
No pq kernel.
No L6.
No large symbolic expansion.

Purpose
-------
The r=5 data from Experiment 249/250 are suspicious:

  * every j has apparent degree 5 when all six d-values are used;
  * nevertheless many early values agree with the old r<=4-style law;
  * j=5 is especially revealing:
        d=10  -> -2
        d=12  -> -36
        d=14  -> -196
     which exactly matches the old candidate for several k,
     but later values suddenly diverge.

This experiment asks:

  1. Where does the degree-j law first break?
  2. Do short consecutive windows still have degree j?
  3. Which d-value is the first anomalous point?
  4. Does removing the anomalous tail restore the old law?
  5. Is the anomaly concentrated at the largest d only?
  6. For each j, what is the maximal prefix that obeys the old law?
  7. Does the discrepancy have a common onset d?
  8. Does the discrepancy itself have a simple factor in d?

Interpolation is used only locally to diagnose the observed finite data.
No interpolated polynomial is accepted as a theorem.
==============================================================================

"""

import sympy as sp


# =============================================================================
# DATA
# =============================================================================

DS = [6, 8, 10, 12, 14, 16]


DATA = {
    3: {
        0: [-441, -3409, -12663, -50184, -127127, -273273],
        1: [-105, -1995, -15560, -44352, -107415, -245245],
        2: [0, -511, -3074, -5859, -51674, -137137],
        3: [0, 0, -903, -15918, -30562, -46431],
        4: [0, 0, -83, 2310, 28144, 3270],
        5: [0, 0, -2, -36, -23996, -31236],
    },

    5: {
        0: [-2142, -18522, -75999, -219024, -589246, -1272726],
        1: [-336, -6300, -35280, -139160, -345072, -756756],
        2: [0, -1140, -9660, -33514, -68904, -285240],
        3: [0, 0, -1485, -8190, -58200, -102015],
        4: [0, 0, -105, -875, 10691, 102730],
        5: [0, 0, -2, -36, -196, -52472],
    },

    7: {
        0: [-7062, -60522, -256674, -744534, -1818124, -4071882],
        1: [-825, -15345, -85470, -285516, -797895, -1832523],
        2: [0, -2145, -18095, -73513, -151284, -348205],
        3: [0, 0, -2211, -12166, -14014, -78287],
        4: [0, 0, -127, -1057, 1323, 140908],
        5: [0, 0, -2, -36, -196, 55384],
    },

    9: {
        0: [-18447, -157157, -663949, -1912482, -4552218, -10123763],
        1: [-1716, -31746, -176176, -576576, -1305876, -3369014],
        2: [0, -3614, -30394, -122031, -176280, -157730],
        3: [0, 0, -3081, -16926, 13650, 516945],
        4: [0, 0, -149, -1239, 5782, 378750],
        5: [0, 0, -2, -36, -196, 131520],
    },

    11: {
        0: [-41223, -349713, -1473381, -4135131, -8999991, -20449143],
        1: [-3185, -58695, -324870, -1039402, -1696695, -3132129],
        2: [0, -5635, -47285, -187502, 8470, 2166395],
        3: [0, 0, -4095, -22470, 84078, 2217105],
        4: [0, 0, -171, -1421, 13458, 1048278],
        5: [0, 0, -2, -36, -196, 276462],
    },

    13: {
        0: [-82348, -696388, -2927876, -7973952, -13749736, -27084876],
        1: [-5440, -99960, -552160, -1722576, -915552, 7331896],
        2: [0, -8296, -69496, -271932, 640968, 11736868],
        3: [0, 0, -5253, -28798, 228276, 6941253],
        4: [0, 0, -193, -1603, 25641, 2501016],
        5: [0, 0, -2, -36, -196, 531860],
    },
}


# =============================================================================
# SYMBOLS
# =============================================================================

D = sp.symbols("D")


# =============================================================================
# OLD UNIVERSAL LAW
# =============================================================================

def old_universal_delta(r, j, k, d):
    """
    Important:
    `d` is allowed to be symbolic.

    This fixes the bug from the previous experiment where:
        sp.Integer(D)
    was attempted.
    """

    r = sp.Integer(r)
    j = sp.Integer(j)
    k = sp.Integer(k)

    if j == 0:
        return sp.binomial(k + r - 1, r)

    return sp.factor(
        sp.Rational((-1) ** j, sp.factorial(j))
        * sp.binomial(k + r - 1, r - j)
        * sp.prod(
            d - r - m
            for m in range(1, j)
        )
        * (
            d
            - sp.Rational((r - j) * k, k + j)
        )
    )


# =============================================================================
# FINITE DIFFERENCES
# =============================================================================

def differences(values):
    rows = [list(map(sp.Integer, values))]

    while len(rows[-1]) > 1:
        prev = rows[-1]

        rows.append([
            sp.simplify(prev[i + 1] - prev[i])
            for i in range(len(prev) - 1)
        ])

    return rows


def exact_degree(values):
    rows = differences(values)

    for order, row in enumerate(rows[:-1]):
        if len(row) <= 1:
            break

        if all(x == row[0] for x in row):
            return order

    return len(values) - 1


# =============================================================================
# LOCAL POLYNOMIAL THROUGH A WINDOW
# =============================================================================

def interpolate_window(xs, ys):
    return sp.factor(
        sp.interpolate(
            list(zip(xs, ys)),
            D
        )
    )


# =============================================================================
# CHECK A POLYNOMIAL AGAINST A POINT
# =============================================================================

def residual_at(poly, d, actual):
    return sp.factor(
        poly.subs(D, d) - actual
    )


# =============================================================================
# 1. LOCAL WINDOW DEGREE
# =============================================================================

def section_1_local_window_degree():
    print("=" * 78)
    print("1. LOCAL WINDOW DEGREE TEST")
    print("=" * 78)

    for k in sorted(DATA):
        print(f"\nk={k}")

        for j in range(6):
            values = DATA[k][j]

            print(f"  j={j}")

            for size in range(2, 7):
                degrees = []

                for start in range(0, len(DS) - size + 1):
                    ys = values[start:start + size]
                    deg = exact_degree(ys)
                    degrees.append(deg)

                print(
                    f"    window={size} degrees={degrees}"
                )

    print()


# =============================================================================
# 2. FIRST BREAK FROM OLD LAW
# =============================================================================

def section_2_old_law_breakpoints():
    print("=" * 78)
    print("2. FIRST BREAK FROM OLD UNIVERSAL LAW")
    print("=" * 78)

    for k in sorted(DATA):
        print(f"\nk={k}")

        for j in range(6):

            matches = []

            for d, actual in zip(
                DS,
                DATA[k][j]
            ):
                expected = sp.simplify(
                    old_universal_delta(
                        5,
                        j,
                        k,
                        sp.Integer(d)
                    )
                )

                matches.append(
                    actual == expected
                )

            first_failure = None

            for idx, ok in enumerate(matches):
                if not ok:
                    first_failure = DS[idx]
                    break

            print(
                f"  j={j} "
                f"matches={matches} "
                f"first_failure={first_failure}"
            )

    print()


# =============================================================================
# 3. MAXIMAL PREFIX AGREEMENT
# =============================================================================

def section_3_prefix_lengths():
    print("=" * 78)
    print("3. MAXIMAL PREFIX AGREEMENT WITH OLD LAW")
    print("=" * 78)

    for k in sorted(DATA):
        print(f"\nk={k}")

        for j in range(6):

            count = 0

            for d, actual in zip(
                DS,
                DATA[k][j]
            ):
                expected = sp.simplify(
                    old_universal_delta(
                        5,
                        j,
                        k,
                        sp.Integer(d)
                    )
                )

                if actual == expected:
                    count += 1
                else:
                    break

            print(
                f"  j={j}: "
                f"prefix_length={count} "
                f"last_matching_d="
                f"{DS[count - 1] if count else None}"
            )

    print()


# =============================================================================
# 4. LOCAL J-DEGREE FIT
# =============================================================================

def section_4_degree_j_fit():
    print("=" * 78)
    print("4. LOCAL DEGREE-j FIT TEST")
    print("=" * 78)

    for k in sorted(DATA):
        print(f"\nk={k}")

        for j in range(6):

            required = j + 1

            if required > len(DS):
                continue

            print(f"  j={j}")

            # Sliding windows with exactly j+1 points.
            for start in range(
                0,
                len(DS) - required + 1
            ):
                xs = DS[start:start + required]
                ys = DATA[k][j][start:start + required]

                poly = interpolate_window(
                    xs,
                    ys
                )

                print(
                    f"    window={xs} "
                    f"poly={poly}"
                )

    print()


# =============================================================================
# 5. J+1 POINT PREDICTION TEST
# =============================================================================

def section_5_prediction_test():
    print("=" * 78)
    print("5. DEGREE-j PREDICTION TEST")
    print("=" * 78)

    """
    A true degree-j law is determined by j+1 points.

    Take the first j+1 points and predict every later point.
    This is much more informative than simply interpolating all six points.
    """

    failures = 0

    for k in sorted(DATA):
        print(f"\nk={k}")

        for j in range(6):

            n = j + 1

            xs = DS[:n]
            ys = DATA[k][j][:n]

            poly = interpolate_window(
                xs,
                ys
            )

            print(f"  j={j}")
            print(f"    basis polynomial={poly}")

            for d, actual in zip(
                DS[n:],
                DATA[k][j][n:]
            ):
                pred = sp.factor(
                    poly.subs(D, d)
                )

                residual = sp.factor(
                    actual - pred
                )

                ok = residual == 0

                if not ok:
                    failures += 1

                print(
                    f"    d={d} "
                    f"actual={actual} "
                    f"predicted={pred} "
                    f"residual={residual} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

    print()
    print(
        "degree-j prediction failures =",
        failures
    )
    print()


# =============================================================================
# 6. DISCREPANCY SEQUENCES
# =============================================================================

def section_6_discrepancy_sequences():
    print("=" * 78)
    print("6. DISCREPANCY SEQUENCES")
    print("=" * 78)

    for k in sorted(DATA):
        print(f"\nk={k}")

        for j in range(6):

            residuals = []

            for d, actual in zip(
                DS,
                DATA[k][j]
            ):
                expected = sp.simplify(
                    old_universal_delta(
                        5,
                        j,
                        k,
                        sp.Integer(d)
                    )
                )

                residuals.append(
                    sp.Integer(actual) - expected
                )

            print(
                f"  j={j}: {residuals}"
            )

    print()


# =============================================================================
# 7. FIRST-FAILURE SHAPE
# =============================================================================

def section_7_first_failure_shape():
    print("=" * 78)
    print("7. FIRST-FAILURE RESIDUAL SHAPE")
    print("=" * 78)

    """
    For each j, find the earliest d where the old law fails.

    Then compare the residual across k.

    If residual / a simple binomial factor is independent of k,
    we may have isolated the missing term.
    """

    for j in range(6):

        print(f"\nj={j}")

        for k in sorted(DATA):

            first = None

            for d, actual in zip(
                DS,
                DATA[k][j]
            ):
                expected = sp.simplify(
                    old_universal_delta(
                        5,
                        j,
                        k,
                        sp.Integer(d)
                    )
                )

                residual = sp.Integer(actual) - expected

                if residual != 0:
                    first = (
                        d,
                        sp.factor(residual)
                    )
                    break

            print(
                f"  k={k}: first_failure={first}"
            )

    print()


# =============================================================================
# 8. COMMON ONSET AUDIT
# =============================================================================

def section_8_common_onset():
    print("=" * 78)
    print("8. COMMON ONSET d AUDIT")
    print("=" * 78)

    for j in range(6):

        onset = {}

        for k in sorted(DATA):

            for d, actual in zip(
                DS,
                DATA[k][j]
            ):

                expected = sp.simplify(
                    old_universal_delta(
                        5,
                        j,
                        k,
                        sp.Integer(d)
                    )
                )

                if actual != expected:
                    onset[k] = d
                    break

        print(
            f"j={j}: onsets={onset}"
        )

    print()


# =============================================================================
# 9. TERMINAL-TAIL AUDIT
# =============================================================================

def section_9_terminal_tail():
    print("=" * 78)
    print("9. TERMINAL TAIL AUDIT")
    print("=" * 78)

    for j in range(6):

        print(f"\nj={j}")

        for k in sorted(DATA):

            print(
                f"  k={k}:",
                end=" "
            )

            for d, actual in zip(
                DS,
                DATA[k][j]
            ):
                if d >= 10:
                    expected = sp.simplify(
                        old_universal_delta(
                            5,
                            j,
                            k,
                            sp.Integer(d)
                        )
                    )

                    residual = sp.factor(
                        sp.Integer(actual) - expected
                    )

                    print(
                        f"d={d}:r={residual}",
                        end=" | "
                    )

            print()

    print()


# =============================================================================
# 10. CROSS-k NORMALIZED DISCREPANCY
# =============================================================================

def section_10_normalized_discrepancy():
    print("=" * 78)
    print("10. CROSS-k NORMALIZED DISCREPANCY")
    print("=" * 78)

    """
    Compare the discrepancy after dividing by the natural k-dependent
    scale C(k+4, 5-j).

    We are not claiming this normalization is correct; it is a diagnostic.
    """

    for j in range(6):

        print(f"\nj={j}")

        for d in DS:

            vals = []

            for k in sorted(DATA):

                actual = sp.Integer(
                    DATA[k][j][DS.index(d)]
                )

                expected = sp.simplify(
                    old_universal_delta(
                        5,
                        j,
                        k,
                        sp.Integer(d)
                    )
                )

                residual = actual - expected

                scale = sp.binomial(
                    k + 4,
                    5 - j
                )

                vals.append(
                    sp.factor(
                        residual / scale
                    )
                )

            print(
                f"  d={d}: {vals}"
            )

    print()


# =============================================================================
# 11. FINAL DIAGNOSTIC
# =============================================================================

def section_11_final():
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
Interpretation:

A) If the first j+1 points predict the remaining points exactly,
   the old degree-j ladder survives locally.

B) If the prediction fails beginning at one particular d for many k,
   the problem is probably a support/domain transition rather than
   a new polynomial law.

C) If the onset d is common across k, that is especially strong evidence
   of a kernel-support boundary.

D) If the discrepancy first appears at the largest d only, do NOT fit
   a degree-5 polynomial to all six points. The degree-5 interpolation
   is then an artifact of one contaminated endpoint.

E) If j=5 remains exact through d=14 but fails at d=16, that strongly
   suggests that the apparent r=5 failure is caused by the extraction
   becoming invalid at the terminal coefficient range.

No r=6 and no kernel expansion are allowed here.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    section_1_local_window_degree()
    section_2_old_law_breakpoints()
    section_3_prefix_lengths()
    section_4_degree_j_fit()
    section_5_prediction_test()
    section_6_discrepancy_sequences()
    section_7_first_failure_shape()
    section_8_common_onset()
    section_9_terminal_tail()
    section_10_normalized_discrepancy()
    section_11_final()


if __name__ == "__main__":
    main()

