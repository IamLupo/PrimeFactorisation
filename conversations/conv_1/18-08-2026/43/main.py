"""
==============================================================================
EXPERIMENT 250
EXACT r=5 BOUNDARY DATA FORENSICS
==============================================================================

Purpose
-------
Investigate the r=5 boundary discrepancy WITHOUT using the pq kernel.

The previous experiments show:

  * r <= 4 closes perfectly.
  * the proposed universal law works symbolically for r <= 4.
  * r=5 interior law is correct.
  * r=5 raw/interior arithmetic is internally consistent.
  * r=5 boundary data do NOT follow the old universal law.
  * r=5 extracted boundary sequences become degree-5 objects instead
    of the expected degree-j ladder.

This experiment therefore treats the observed r=5 deltas as exact data
and asks:

  1. What polynomial is actually determined by the data?
  2. Which roots are forced exactly?
  3. Which roots fail?
  4. Do the "stable" values contain a common factor across k?
  5. Can the r=5 discrepancy be decomposed into a universal contaminant?
  6. Is there a low-degree correction common to every j?
  7. Does removing the same contaminant restore degree j?

No interpolation assumption is used to claim a law.
Interpolation is used only as an exact diagnostic of the finite data.

All arithmetic is exact over QQ.
No floating point.
No filesystem.
No kernel construction.
No r=6.
No L5.
==============================================================================
"""

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

D = sp.symbols("D")


# =============================================================================
# RAW r=5 DELTA DATA FROM EXPERIMENT 249
# =============================================================================

DATA = {
    3: {
        0: [(-441), (-3409), (-12663), (-50184), (-127127), (-273273)],
        1: [(-105), (-1995), (-15560), (-44352), (-107415), (-245245)],
        2: [0, (-511), (-3074), (-5859), (-51674), (-137137)],
        3: [0, 0, (-903), (-15918), (-30562), (-46431)],
        4: [0, 0, (-83), 2310, 28144, 3270],
        5: [0, 0, (-2), (-36), (-23996), (-31236)],
    },

    5: {
        0: [(-2142), (-18522), (-75999), (-219024), (-589246), (-1272726)],
        1: [(-336), (-6300), (-35280), (-139160), (-345072), (-756756)],
        2: [0, (-1140), (-9660), (-33514), (-68904), (-285240)],
        3: [0, 0, (-1485), (-8190), (-58200), (-102015)],
        4: [0, 0, (-105), (-875), 10691, 102730],
        5: [0, 0, (-2), (-36), (-196), (-52472)],
    },

    7: {
        0: [(-7062), (-60522), (-256674), (-744534), (-1818124), (-4071882)],
        1: [(-825), (-15345), (-85470), (-285516), (-797895), (-1832523)],
        2: [0, (-2145), (-18095), (-73513), (-151284), (-348205)],
        3: [0, 0, (-2211), (-12166), (-14014), (-78287)],
        4: [0, 0, (-127), (-1057), 1323, 140908],
        5: [0, 0, (-2), (-36), (-196), 55384],
    },

    9: {
        0: [(-18447), (-157157), (-663949), (-1912482), (-4552218), (-10123763)],
        1: [(-1716), (-31746), (-176176), (-576576), (-1305876), (-3369014)],
        2: [0, (-3614), (-30394), (-122031), (-176280), (-157730)],
        3: [0, 0, (-3081), (-16926), 13650, 516945],
        4: [0, 0, (-149), (-1239), 5782, 378750],
        5: [0, 0, (-2), (-36), (-196), 131520],
    },

    11: {
        0: [(-41223), (-349713), (-1473381), (-4135131), (-8999991), (-20449143)],
        1: [(-3185), (-58695), (-324870), (-1039402), (-1696695), (-3132129)],
        2: [0, (-5635), (-47285), (-187502), 8470, 2166395],
        3: [0, 0, (-4095), (-22470), 84078, 2217105],
        4: [0, 0, (-171), (-1421), 13458, 1048278],
        5: [0, 0, (-2), (-36), (-196), 276462],
    },

    13: {
        0: [(-82348), (-696388), (-2927876), (-7973952), (-13749736), (-27084876)],
        1: [(-5440), (-99960), (-552160), (-1722576), (-915552), 7331896],
        2: [0, (-8296), (-69496), (-271932), 640968, 11736868],
        3: [0, 0, (-5253), (-28798), 228276, 6941253],
        4: [0, 0, (-193), (-1603), 25641, 2501016],
        5: [0, 0, (-2), (-36), (-196), 531860],
    },
}


# Stable d-values corresponding to the six data points.
DS = [6, 8, 10, 12, 14, 16]


# =============================================================================
# OLD PROPOSED UNIVERSAL LAW
# =============================================================================

def old_universal_delta(r, j, k, d):
    """
    The candidate law from the r<=4 experiments.
    Included only as a comparison baseline.
    """
    d = sp.Integer(d)
    k = sp.Integer(k)

    if j == 0:
        return sp.binomial(k + r - 1, r)

    return sp.simplify(
        ((-1) ** j)
        / sp.factorial(j)
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
# EXACT INTERPOLATION
# =============================================================================

def interpolate_exact(xs, ys):
    xs = list(map(sp.Integer, xs))
    ys = list(map(sp.Integer, ys))

    return sp.factor(
        sp.interpolate(
            list(zip(xs, ys)),
            D
        )
    )


# =============================================================================
# FINITE DIFFERENCES
# =============================================================================

def finite_difference_table(values):
    rows = [list(map(sp.Integer, values))]

    while len(rows[-1]) > 1:
        prev = rows[-1]
        rows.append([
            sp.simplify(prev[i + 1] - prev[i])
            for i in range(len(prev) - 1)
        ])

    return rows


def polynomial_degree_from_differences(values):
    rows = finite_difference_table(values)

    for order, row in enumerate(rows):
        if len(row) <= 1:
            return order

        if all(x == row[0] for x in row):
            return order

    return len(values) - 1


# =============================================================================
# EXPECTED KNOWN ROOTS FOR OLD LAW
# =============================================================================

def expected_roots_old(r, j, k):
    roots = []

    if j >= 2:
        roots.extend(
            sp.Integer(r + m)
            for m in range(1, j)
        )

    moving = sp.Rational(
        (r - j) * k,
        k + j
    )

    roots.append(moving)

    return roots


# =============================================================================
# ROOT TESTING
# =============================================================================

def exact_root_residual(poly, root):
    return sp.factor(poly.subs(D, root))


# =============================================================================
# SECTION 1
# =============================================================================

def section_1_actual_degrees():
    print("=" * 78)
    print("1. ACTUAL r=5 FINITE-DIFFERENCE DEGREES")
    print("=" * 78)

    for k in DATA:
        print(f"k={k}")

        for j in range(6):
            values = DATA[k][j]
            degree = polynomial_degree_from_differences(values)

            print(
                f"  j={j} degree={degree} expected_old={j}"
            )

        print()


# =============================================================================
# SECTION 2
# =============================================================================

def section_2_exact_polynomials():
    print("=" * 78)
    print("2. EXACT INTERPOLATED POLYNOMIALS")
    print("=" * 78)

    polys = {}

    for k in DATA:
        polys[k] = {}

        print(f"\nk={k}")

        for j in range(6):
            poly = interpolate_exact(DS, DATA[k][j])
            polys[k][j] = poly

            print(
                f"  j={j}: {poly}"
            )

    return polys


# =============================================================================
# SECTION 3
# =============================================================================

def section_3_forced_root_audit(polys):
    print("=" * 78)
    print("3. FORCED-ROOT AUDIT")
    print("=" * 78)

    for k in sorted(polys):
        print(f"\nk={k}")

        for j in range(6):
            poly = polys[k][j]

            print(f"  j={j}")

            # Roots predicted by old law.
            roots = expected_roots_old(
                5, j, k
            )

            for root in roots:
                residual = exact_root_residual(
                    poly, root
                )

                print(
                    f"    old-root {root}: "
                    f"residual={residual} "
                    f"{'PASS' if residual == 0 else 'FAIL'}"
                )

            # Fixed roots that appear empirically from the first valid d points.
            for root in [6, 7, 8, 9, 10]:
                residual = exact_root_residual(
                    poly, root
                )

                if residual == 0:
                    print(
                        f"    empirical-root {root}: PASS"
                    )

    print()


# =============================================================================
# SECTION 4
# =============================================================================

def section_4_factorization(polys):
    print("=" * 78)
    print("4. EXACT FACTORIZATION")
    print("=" * 78)

    for k in sorted(polys):
        print(f"\nk={k}")

        for j in range(6):
            poly = polys[k][j]

            print(
                f"  j={j}: {sp.factor(poly)}"
            )

    print()


# =============================================================================
# SECTION 5
# =============================================================================

def section_5_old_law_difference(polys):
    print("=" * 78)
    print("5. ACTUAL POLYNOMIAL - OLD UNIVERSAL LAW")
    print("=" * 78)

    for k in sorted(polys):
        print(f"\nk={k}")

        for j in range(6):
            poly = polys[k][j]

            old = sp.expand(
                old_universal_delta(
                    5,
                    j,
                    k,
                    D
                )
            )

            correction = sp.factor(
                poly - old
            )

            print(
                f"  j={j}: correction = {correction}"
            )

    print()


# =============================================================================
# SECTION 6
# =============================================================================

def section_6_common_correction_factor(polys):
    print("=" * 78)
    print("6. COMMON-CORRECTION FACTOR AUDIT")
    print("=" * 78)

    """
    For each j, compute

        C_j(k,D) = Actual_j(k,D) - OldLaw_j(k,D)

    Then determine whether C_j has factors common across all k.

    A genuinely universal contaminant should survive the gcd across k.
    """

    for j in range(6):
        print(f"\nj={j}")

        corrections = []

        for k in sorted(polys):
            actual = polys[k][j]

            old = sp.expand(
                old_universal_delta(
                    5,
                    j,
                    k,
                    D
                )
            )

            correction = sp.Poly(
                sp.together(actual - old),
                D,
                domain="QQ"
            )

            corrections.append(
                correction
            )

        common = corrections[0]

        for poly in corrections[1:]:
            common = sp.gcd(common, poly)

        print(
            "  gcd across k =",
            sp.factor(common.as_expr())
        )

        print(
            "  degree =",
            common.degree()
        )

    print()


# =============================================================================
# SECTION 7
# =============================================================================

def section_7_shifted_basis_test(polys):
    print("=" * 78)
    print("7. SHIFTED-BASIS TEST")
    print("=" * 78)

    """
    Test whether replacing D by D+s in the old law gives the observed
    degree-5 polynomial family.

    We only scan small exact integer shifts.

    This is NOT a claim of correctness; it is a fast falsification test.
    """

    shifts = list(range(-8, 9))

    for j in range(6):
        print(f"\nj={j}")

        for shift in shifts:
            ok = True

            for k in sorted(DATA):
                old_shifted = sp.expand(
                    old_universal_delta(
                        5,
                        j,
                        k,
                        D + shift
                    )
                )

                actual = polys[k][j]

                if sp.simplify(actual - old_shifted) != 0:
                    ok = False
                    break

            if ok:
                print(
                    f"  EXACT common integer shift found: {shift}"
                )

    print()


# =============================================================================
# SECTION 8
# =============================================================================

def section_8_common_leading_terms(polys):
    print("=" * 78)
    print("8. CROSS-k LEADING-COEFFICIENT AUDIT")
    print("=" * 78)

    for j in range(6):
        print(f"\nj={j}")

        for k in sorted(polys):
            poly = sp.Poly(
                polys[k][j],
                D
            )

            print(
                f"  k={k:2d} "
                f"degree={poly.degree()} "
                f"leading={poly.LC()}"
            )

    print()


# =============================================================================
# SECTION 9
# =============================================================================

def section_9_normalized_polynomials(polys):
    print("=" * 78)
    print("9. NORMALIZED POLYNOMIALS")
    print("=" * 78)

    """
    Divide each polynomial by its leading coefficient.

    This removes the obvious k-dependent scale and reveals whether the
    remaining shape is common across k.
    """

    for j in range(6):
        print(f"\nj={j}")

        normalized = {}

        for k in sorted(polys):
            poly = sp.Poly(
                polys[k][j],
                D
            )

            lc = poly.LC()

            normalized[k] = sp.factor(
                poly.as_expr() / lc
            )

            print(
                f"  k={k}: {normalized[k]}"
            )

        first = normalized[min(normalized)]

        common_shape = True

        for k in normalized:
            if sp.simplify(
                normalized[k] - first
            ) != 0:
                common_shape = False
                break

        print(
            "  common normalized shape =",
            common_shape
        )

    print()


# =============================================================================
# SECTION 10
# =============================================================================

def section_10_low_degree_fit():
    print("=" * 78)
    print("10. LOW-DEGREE FIT TEST")
    print("=" * 78)

    """
    Determine whether subsets of the stable data agree with degree-j laws.

    This helps identify whether the later d-values are entering another
    regime rather than following one global degree-5 polynomial.
    """

    subsets = {
        "first3": DS[:3],
        "first4": DS[:4],
        "first5": DS[:5],
        "last3": DS[-3:],
        "last4": DS[-4:],
    }

    for k in sorted(DATA):
        print(f"\nk={k}")

        for j in range(6):
            print(f"  j={j}")

            values = DATA[k][j]

            for name, xs in subsets.items():
                ys = [
                    values[DS.index(x)]
                    for x in xs
                ]

                degree = polynomial_degree_from_differences(
                    ys
                )

                poly = interpolate_exact(
                    xs,
                    ys
                )

                print(
                    f"    {name:6s} "
                    f"degree={degree} "
                    f"poly={poly}"
                )

    print()


# =============================================================================
# SECTION 11
# =============================================================================

def section_11_rank_test():
    print("=" * 78)
    print("11. CROSS-k DATA MATRIX RANK TEST")
    print("=" * 78)

    """
    For each j and each d-value, form the matrix:

        row = k
        column = d

    Then inspect exact rank.

    A universal polynomial family with low-dimensional k dependence can
    exhibit unexpectedly small rank after suitable normalization.
    """

    for j in range(6):

        matrix = []

        for k in sorted(DATA):
            matrix.append([
                sp.Integer(v)
                for v in DATA[k][j]
            ])

        M = sp.Matrix(matrix)

        print(
            f"j={j} raw rank={M.rank()} "
            f"shape={M.shape}"
        )

    print()


# =============================================================================
# SECTION 12
# =============================================================================

def section_12_summary(polys):
    print("=" * 78)
    print("12. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
Interpretation rules:

1. If the actual r=5 polynomials have a common correction factor
   across k, investigate that factor first.

2. If normalized shapes coincide across k, the discrepancy is probably
   a universal basis/extraction issue rather than a genuinely new
   r=5 boundary mechanism.

3. If low-d subsets have degree j but the full six-point set jumps
   to degree 5, the boundary is likely changing regime with d.

4. If neither phenomenon occurs, the r=5 boundary mechanism genuinely
   differs from r<=4.

No conclusion is accepted from interpolation alone.

This experiment deliberately does not use the pq kernel.
Its purpose is to tell us WHAT must be explained next before touching
another expensive kernel computation.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    section_1_actual_degrees()

    polys = section_2_exact_polynomials()

    section_3_forced_root_audit(polys)
    section_4_factorization(polys)
    section_5_old_law_difference(polys)
    section_6_common_correction_factor(polys)
    section_7_shifted_basis_test(polys)
    section_8_common_leading_terms(polys)
    section_9_normalized_polynomials(polys)
    section_10_low_degree_fit()
    section_11_rank_test()
    section_12_summary(polys)


if __name__ == "__main__":
    main()