"""
==============================================================================
EXPERIMENT 252
EXACT r=5 DISCREPANCY FACTORIZATION
==============================================================================

Standalone.
Exact arithmetic over QQ.
No kernel.
No filesystem.
No r=6.
No large bivariate interpolation.

Goal
----
Study

    E_j(k,d) = actual_r5_boundary(j,k,d)
               - old_universal_boundary(j,k,d)

The previous experiment established:

    j=0  : discrepancy already nonzero at d=6
    j=1  : discrepancy starts at d=8, zero at d=6
    j=2  : discrepancy starts at d=8, zero at d=6
    j=3  : discrepancy starts at d=8, zero at d=6
    j=4  : discrepancy starts at d=10, zero at d=6,8
    j=5  : discrepancy is zero through d=12, then breaks

Hypothesis to test:
-------------------
The discrepancy has its OWN boundary/support product.

Candidate onset factors:

    j=1 : (d-6)
    j=2 : (d-6)
    j=3 : (d-6)
    j=4 : (d-6)(d-8)
    j=5 : (d-6)(d-8)(d-10)(d-12)

We do NOT assume these factors are correct.
We test them exactly.

Questions:
----------
1. Does E_j vanish at the observed onset points?
2. After dividing by the candidate d-factor, does the quotient have
   a simple dependence on k?
3. Is the quotient polynomial in k of low degree?
4. Does the quotient factor?
5. Is there a common pattern across j?
6. Does the anomalous k=3,j=5,d=14 point remain anomalous after
   normalization?
7. Are the huge degree-5 interpolants merely the consequence of
   multiplying a simple support factor by another low-degree object?

"""

import sympy as sp


# ============================================================================
# DATA
# ============================================================================

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


K = sp.symbols("K")
D = sp.symbols("D")


# ============================================================================
# OLD UNIVERSAL LAW
# ============================================================================

def old_delta(j, k, d):
    r = sp.Integer(5)
    j = sp.Integer(j)
    k = sp.Integer(k)

    if j == 0:
        return sp.binomial(k + 4, 5)

    return sp.factor(
        sp.Rational((-1) ** j, sp.factorial(j))
        * sp.binomial(k + 4, 5 - j)
        * sp.prod(
            d - 5 - m
            for m in range(1, j)
        )
        * (
            d
            - sp.Rational((5 - j) * k, k + j)
        )
    )


def discrepancy(k, j, d):
    actual = sp.Integer(
        DATA[k][j][DS.index(d)]
    )

    return sp.factor(
        actual - old_delta(j, k, sp.Integer(d))
    )


# ============================================================================
# CANDIDATE SUPPORT FACTOR
# ============================================================================

def support_factor(j, d):
    """
    Candidate factor suggested by the common zero onset.

    This is a diagnostic only.
    """

    if j <= 0:
        return sp.Integer(1)

    roots = {
        1: [6],
        2: [6],
        3: [6],
        4: [6, 8],
        5: [6, 8, 10, 12],
    }[j]

    out = sp.Integer(1)

    for root in roots:
        out *= d - root

    return sp.factor(out)


# ============================================================================
# 1. EXACT ZERO AUDIT
# ============================================================================

def section_1_zero_audit():
    print("=" * 78)
    print("1. EXACT DISCREPANCY ZERO AUDIT")
    print("=" * 78)

    for j in range(6):

        roots = []

        for d in DS:
            vals = [
                discrepancy(k, j, d)
                for k in sorted(DATA)
            ]

            if all(v == 0 for v in vals):
                roots.append(d)

        print(
            f"j={j}: common_zero_d={roots}"
        )

    print()


# ============================================================================
# 2. FACTORED DISCREPANCIES
# ============================================================================

def section_2_factored_discrepancies():
    print("=" * 78)
    print("2. RAW DISCREPANCY FACTORIZATION")
    print("=" * 78)

    for j in range(6):

        print(f"\nj={j}")

        for k in sorted(DATA):

            factors = []

            for d in DS:
                e = discrepancy(k, j, d)

                if e != 0:
                    factors.append(
                        f"d={d}:{sp.factor(e)}"
                    )

            print(
                f"  k={k}: "
                + " | ".join(factors)
            )

    print()


# ============================================================================
# 3. NORMALIZE BY CANDIDATE SUPPORT FACTOR
# ============================================================================

def section_3_normalized():
    print("=" * 78)
    print("3. DISCREPANCY / CANDIDATE SUPPORT FACTOR")
    print("=" * 78)

    for j in range(1, 6):

        print(f"\nj={j}")

        for k in sorted(DATA):

            vals = []

            for d in DS:

                sf = support_factor(
                    j,
                    sp.Integer(d)
                )

                e = discrepancy(k, j, d)

                if sf == 0:
                    if e == 0:
                        q = sp.Integer(0)
                    else:
                        q = "NONZERO/0"
                else:
                    q = sp.factor(e / sf)

                vals.append(
                    f"d={d}:{q}"
                )

            print(
                f"  k={k}: "
                + " | ".join(vals)
            )

    print()


# ============================================================================
# 4. CROSS-k POLYNOMIAL DEGREE OF NORMALIZED QUOTIENT
# ============================================================================

def section_4_cross_k_degree():
    print("=" * 78)
    print("4. CROSS-k DEGREE OF NORMALIZED DISCREPANCY")
    print("=" * 78)

    ks = sorted(DATA)

    for j in range(1, 6):

        print(f"\nj={j}")

        for d in DS:

            sf = support_factor(
                j,
                sp.Integer(d)
            )

            if sf == 0:
                continue

            vals = []

            for k in ks:
                e = discrepancy(k, j, d)

                vals.append(
                    sp.Rational(e, sf)
                )

            poly = sp.factor(
                sp.interpolate(
                    list(zip(ks, vals)),
                    K
                )
            )

            degree = sp.degree(
                sp.Poly(poly, K),
                K
            )

            print(
                f"  d={d}: "
                f"degree_in_k={degree} "
                f"poly={poly}"
            )

    print()


# ============================================================================
# 5. CROSS-k FACTORIZATION
# ============================================================================

def section_5_cross_k_factorization():
    print("=" * 78)
    print("5. CROSS-k FACTORIZATION OF NORMALIZED DISCREPANCY")
    print("=" * 78)

    ks = sorted(DATA)

    for j in range(1, 6):

        print(f"\nj={j}")

        for d in DS:

            sf = support_factor(
                j,
                sp.Integer(d)
            )

            if sf == 0:
                continue

            vals = [
                sp.Rational(
                    discrepancy(k, j, d),
                    sf
                )
                for k in ks
            ]

            poly = sp.factor(
                sp.interpolate(
                    list(zip(ks, vals)),
                    K
                )
            )

            print(
                f"  d={d}: "
                f"{sp.factor(poly)}"
            )

    print()


# ============================================================================
# 6. NORMALIZE BY OBVIOUS k-SCALES
# ============================================================================

def section_6_natural_k_normalizations():
    print("=" * 78)
    print("6. NATURAL k-NORMALIZATION AUDIT")
    print("=" * 78)

    for j in range(1, 6):

        print(f"\nj={j}")

        for d in DS:

            sf = support_factor(
                j,
                sp.Integer(d)
            )

            if sf == 0:
                continue

            values = []

            for k in sorted(DATA):

                e = discrepancy(k, j, d)

                q = sp.Rational(e, sf)

                candidates = {
                    "C(k+4,5-j)":
                        sp.binomial(k + 4, 5 - j),

                    "k+5":
                        k + 5,

                    "k+j":
                        k + j,

                    "k+1":
                        k + 1,
                }

                normalized = {
                    name: sp.factor(q / scale)
                    for name, scale in candidates.items()
                    if scale != 0
                }

                values.append(
                    (k, normalized)
                )

            print(f"  d={d}")

            for k, vals in values:
                print(
                    f"    k={k}: {vals}"
                )

    print()


# ============================================================================
# 7. SPECIAL j=5 AUDIT
# ============================================================================

def section_7_j5():
    print("=" * 78)
    print("7. SPECIAL j=5 SUPPORT AUDIT")
    print("=" * 78)

    j = 5

    roots = [6, 8, 10, 12]

    print(
        "candidate support factor =",
        sp.factor(
            sp.prod(D - r for r in roots)
        )
    )

    for k in sorted(DATA):

        print(f"\nk={k}")

        for d in DS:

            e = discrepancy(k, j, d)

            sf = support_factor(
                j,
                sp.Integer(d)
            )

            if sf != 0:
                q = sp.factor(
                    sp.Rational(e, sf)
                )
            else:
                q = "undefined"

            print(
                f"  d={d} "
                f"E={e} "
                f"quotient={q}"
            )

    print()


# ============================================================================
# 8. LOOK FOR COMMON FORMULAS AT EACH d
# ============================================================================

def section_8_d_slices():
    print("=" * 78)
    print("8. NORMALIZED d-SLICES AS POLYNOMIALS IN k")
    print("=" * 78)

    ks = sorted(DATA)

    for j in range(1, 6):

        print(f"\nj={j}")

        for d in DS:

            sf = support_factor(
                j,
                sp.Integer(d)
            )

            if sf == 0:
                continue

            vals = [
                sp.Rational(
                    discrepancy(k, j, d),
                    sf
                )
                for k in ks
            ]

            poly = sp.factor(
                sp.interpolate(
                    list(zip(ks, vals)),
                    K
                )
            )

            print(
                f"  d={d}: "
                f"q_j(K)={poly}"
            )

    print()


# ============================================================================
# 9. FINAL SUMMARY
# ============================================================================

def section_9_summary():
    print("=" * 78)
    print("9. FINAL SUMMARY")
    print("=" * 78)

    print(
        """
Interpret the output as follows:

1. If the discrepancy has common zeros exactly at the candidate
   support roots, that is evidence for a second boundary mechanism.

2. If E_j / support_factor(j,d) has a low-degree polynomial dependence
   on k, that is much stronger evidence than the old r<=4 interpolation.

3. If those normalized k-polynomials themselves factor into simple
   binomial/polynomial pieces, record them as candidate correction laws.

4. j=5 is especially important:
       E_5 = 0 at d=6,8,10,12
   for every k in the data.

   That is impossible to explain as a random terminal interpolation
   artifact. It strongly suggests genuine support truncation.

5. Do NOT yet propose a replacement universal r,j formula.
   First identify the discrepancy correction.

6. Do NOT run r=6.

7. Do NOT expand the full pq kernel.

The mathematical target is now:

    actual r=5 boundary
      =
    old r<=4-type contribution
      +
    support correction.

This experiment tries to expose that correction directly.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():
    section_1_zero_audit()
    section_2_factored_discrepancies()
    section_3_normalized()
    section_4_cross_k_degree()
    section_5_cross_k_factorization()
    section_6_natural_k_normalizations()
    section_7_j5()
    section_8_d_slices()
    section_9_summary()


if __name__ == "__main__":
    main()

