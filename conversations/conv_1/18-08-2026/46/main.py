from __future__ import annotations

import sympy as sp


# ============================================================================
# EXPERIMENT 254
# EXACT r=5 SUPPORT-CORRECTION STRUCTURE
# ============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No filesystem access
# No previous experiment imported
# No r=6
# No full pq-kernel expansion
#
# This version fixes the SymPy failure:
#
#     Poly.multiplicity(...)
#
# which does not exist in the installed SymPy version.
#
# Root multiplicities are computed explicitly by repeated exact division
# by (D - root).
#
# ============================================================================


D, K = sp.symbols("D K")


# ============================================================================
# 1. EXACT r=5 DATA
# ============================================================================

DATA = {
    3: {
        0: {6: -441, 8: -3409, 10: -12663, 12: -50184, 14: -127127, 16: -273273},
        1: {6: -105, 8: -1995, 10: -15560, 12: -44352, 14: -107415, 16: -245245},
        2: {6: 0, 8: -511, 10: -3074, 12: -5859, 14: -51674, 16: -137137},
        3: {6: 0, 8: -49, 10: -903, 12: -15918, 14: -30562, 16: -46431},
        4: {6: 0, 8: 0, 10: -83, 12: 2310, 14: 28144, 16: 3270},
        5: {6: 0, 8: 0, 10: -2, 12: -36, 14: -23996, 16: -30564},
    },

    5: {
        0: {6: -2142, 8: -18522, 10: -75999, 12: -219024, 14: -589246, 16: -1272726},
        1: {6: -336, 8: -6300, 10: -35280, 12: -139160, 14: -345072, 16: -756756},
        2: {6: 0, 8: -1140, 10: -9660, 12: -33514, 14: -68904, 16: -285240},
        3: {6: 0, 8: -81, 10: -1485, 12: -8190, 14: -58200, 16: -102015},
        4: {6: 0, 8: 0, 10: -105, 12: -875, 14: 10691, 16: 102730},
        5: {6: 0, 8: 0, 10: -2, 12: -36, 14: -196, 16: -52472},
    },

    7: {
        0: {6: -7062, 8: -60522, 10: -256674, 12: -744534, 14: -1818124, 16: -4071882},
        1: {6: -825, 8: -15345, 10: -85470, 12: -285516, 14: -797895, 16: -1832523},
        2: {6: 0, 8: -2145, 10: -18095, 12: -73513, 14: -151284, 16: -348205},
        3: {6: 0, 8: -121, 10: -2211, 12: -12166, 14: -14014, 16: -78287},
        4: {6: 0, 8: 0, 10: -127, 12: -1057, 14: 1323, 16: 140908},
        5: {6: 0, 8: 0, 10: -2, 12: -36, 14: -196, 16: 55384},
    },

    9: {
        0: {6: -18447, 8: -157157, 10: -663949, 12: -1912482, 14: -4552218, 16: -10123763},
        1: {6: -1716, 8: -31746, 10: -176176, 12: -576576, 14: -1305876, 16: -3369014},
        2: {6: 0, 8: -3614, 10: -30394, 12: -122031, 14: -176280, 16: -157730},
        3: {6: 0, 8: -169, 10: -3081, 12: -16926, 14: 13650, 16: 516945},
        4: {6: 0, 8: 0, 10: -149, 12: -1239, 14: 5782, 16: 378750},
        5: {6: 0, 8: 0, 10: -2, 12: -36, 14: -196, 16: 131520},
    },

    11: {
        0: {6: -41223, 8: -349713, 10: -1473381, 12: -4135131, 14: -8999991, 16: -20449143},
        1: {6: -3185, 8: -58695, 10: -324870, 12: -1039402, 14: -1696695, 16: -3132129},
        2: {6: 0, 8: -5635, 10: -47285, 12: -187502, 14: 8470, 16: 2166395},
        3: {6: 0, 8: -225, 10: -4095, 12: -22470, 14: 84078, 16: 2217105},
        4: {6: 0, 8: 0, 10: -171, 12: -1421, 14: 13458, 16: 1048278},
        5: {6: 0, 8: 0, 10: -2, 12: -36, 14: -196, 16: 276462},
    },

    13: {
        0: {6: -82348, 8: -696388, 10: -2927876, 12: -7973952, 14: -13749736, 16: -27084876},
        1: {6: -5440, 8: -99960, 10: -552160, 12: -1722576, 14: -915552, 16: 7331896},
        2: {6: 0, 8: -8296, 10: -69496, 12: -271932, 14: 640968, 16: 11691308},
        3: {6: 0, 8: -289, 10: -5253, 12: -28798, 14: 228276, 16: 6941253},
        4: {6: 0, 8: 0, 10: -193, 12: -1603, 14: 25641, 16: 2493246},
        5: {6: 0, 8: 0, 10: -2, 12: -36, 14: 532056, 16: 533204},
    },
}


# ============================================================================
# 2. OLD r<=4 UNIVERSAL CANDIDATE
# ============================================================================

def old_candidate(r: int, j: int, k: int, d: int | sp.Expr) -> sp.Expr:
    d = sp.sympify(d)
    k = sp.sympify(k)

    if j == 0:
        return sp.binomial(k + r - 1, r)

    expr = (
        sp.Rational((-1) ** j, sp.factorial(j))
        * sp.binomial(k + r - 1, r - j)
    )

    for m in range(1, j):
        expr *= d - r - m

    expr *= d - sp.Rational((r - j) * k, 1) / (k + j)

    return sp.factor(expr)


def discrepancy(r: int, j: int, k: int, d: int) -> sp.Expr:
    actual = sp.Integer(DATA[k][j][d])
    candidate = old_candidate(r, j, k, d)
    return sp.factor(actual - candidate)


# ============================================================================
# 3. COMMON SUPPORT ZEROS
# ============================================================================

def common_zero_set(j: int) -> list[int]:
    ds = sorted(next(iter(DATA.values()))[j])

    zeros = []

    for d in ds:
        if all(discrepancy(5, j, k, d) == 0 for k in DATA):
            zeros.append(d)

    return zeros


def support_factor_from_common_zeros(j: int) -> sp.Expr:
    result = sp.Integer(1)

    for d0 in common_zero_set(j):
        result *= D - d0

    return sp.expand(result)


# ============================================================================
# 4. INTERPOLATED DISCREPANCY POLYNOMIAL
# ============================================================================

def discrepancy_polynomial(j: int, k: int) -> sp.Expr:
    points = [
        (sp.Integer(d), discrepancy(5, j, k, d))
        for d in sorted(DATA[k][j])
    ]

    return sp.factor(sp.interpolate(points, D))


# ============================================================================
# 5. SUPPORT QUOTIENT
# ============================================================================

def quotient_in_d(j: int, k: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """
    Return:
        E(D),
        support(D),
        Q(D) where E = support * Q.

    Raises no exception if exact division fails.
    """

    poly = sp.expand(discrepancy_polynomial(j, k))
    support = sp.expand(support_factor_from_common_zeros(j))

    if support == 0:
        return sp.factor(poly), sp.Integer(0), sp.Integer(0)

    p = sp.Poly(poly, D)
    s = sp.Poly(support, D)

    quotient, remainder = sp.div(p, s)

    if not remainder.is_zero:
        q = sp.factor(
            quotient.as_expr()
            + sp.Symbol("NONZERO_REMAINDER")
        )
    else:
        q = sp.factor(quotient.as_expr())

    return sp.factor(poly), sp.factor(support), q


def quotient_is_exact(q: sp.Expr) -> bool:
    return sp.Symbol("NONZERO_REMAINDER") not in q.free_symbols


# ============================================================================
# 6. CORRECT ROOT MULTIPLICITY FUNCTION
# ============================================================================

def exact_root_multiplicity(expr: sp.Expr, root: sp.Expr) -> int:
    """
    Exact multiplicity of D=root in expr.

    This deliberately avoids Poly.multiplicity(), which is not available
    in some SymPy versions.
    """

    expr = sp.factor(sp.expand(expr))

    if expr == 0:
        return sp.oo

    multiplicity = 0
    current = sp.Poly(expr, D)
    divisor = sp.Poly(D - root, D)

    while True:
        quotient, remainder = sp.div(current, divisor)

        if remainder.is_zero:
            multiplicity += 1
            current = quotient

            if current.is_zero:
                return sp.oo
        else:
            break

    return multiplicity


# ============================================================================
# 7. BASIC AUDIT
# ============================================================================

def section_1_basic_audit() -> None:
    print("=" * 78)
    print("1. BASIC EXACT DISCREPANCY AUDIT")
    print("=" * 78)

    total = 0
    nonzero = 0

    for k in sorted(DATA):
        print(f"k={k}")

        for j in range(6):
            row = []

            for d in sorted(DATA[k][j]):
                e = discrepancy(5, j, k, d)

                total += 1

                if e != 0:
                    nonzero += 1

                row.append((d, str(e)))

            print(f"  j={j}: {row}")

        print()

    print(f"total points = {total}")
    print(f"nonzero discrepancies = {nonzero}")
    print()


# ============================================================================
# 8. COMMON SUPPORT ZEROS
# ============================================================================

def section_2_common_zeros() -> None:
    print("=" * 78)
    print("2. COMMON SUPPORT ZEROS")
    print("=" * 78)

    for j in range(6):
        print(
            f"j={j}: "
            f"common_zero_d={common_zero_set(j)}"
        )

    print()


# ============================================================================
# 9. SUPPORT FACTORIZATION
# ============================================================================

def section_3_support_factorization() -> None:
    print("=" * 78)
    print("3. SUPPORT-FACTOR EXTRACTION")
    print("=" * 78)

    for j in range(6):
        print(f"j={j}")

        support = support_factor_from_common_zeros(j)

        print(
            f"  common support factor = "
            f"{sp.factor(support)}"
        )

        for k in sorted(DATA):
            poly, sf, quotient = quotient_in_d(j, k)

            print(f"  k={k}")
            print(f"    E(D) = {poly}")
            print(f"    support = {sf}")
            print(f"    Q(D) = {quotient}")

        print()


# ============================================================================
# 10. QUOTIENT DEGREE
# ============================================================================

def section_4_quotient_degree() -> None:
    print("=" * 78)
    print("4. QUOTIENT DEGREE AFTER SUPPORT FACTORIZATION")
    print("=" * 78)

    for j in range(6):
        support = support_factor_from_common_zeros(j)
        support_degree = sp.Poly(support, D).degree()

        print(
            f"j={j} "
            f"support_degree={support_degree}"
        )

        for k in sorted(DATA):
            poly, sf, quotient = quotient_in_d(j, k)

            E_degree = sp.Poly(poly, D).degree()

            if not quotient_is_exact(quotient):
                print(
                    f"  k={k}: "
                    f"E_degree={E_degree} "
                    f"Q=NONEXACT"
                )
                continue

            Q_degree = sp.Poly(quotient, D).degree()

            print(
                f"  k={k}: "
                f"E_degree={E_degree} "
                f"Q_degree={Q_degree}"
            )

        print()


# ============================================================================
# 11. NATURAL NORMALIZATION
# ============================================================================

def natural_factor(k: int, j: int) -> sp.Expr:
    return sp.binomial(k + 4, 5 - j)


def normalized_quotient(
    j: int,
    k: int,
    d: int,
) -> sp.Expr:

    support = support_factor_from_common_zeros(j)
    support_value = sp.expand(support.subs(D, d))

    if support_value == 0:
        return sp.Integer(0)

    E = discrepancy(5, j, k, d)

    return sp.factor(
        sp.together(
            E /
            (
                natural_factor(k, j)
                * support_value
            )
        )
    )


# ============================================================================
# 12. NORMALIZED VALUES
# ============================================================================

def section_5_normalized_values() -> None:
    print("=" * 78)
    print("5. NORMALIZED CORRECTION VALUES")
    print("=" * 78)

    for j in range(1, 6):
        support = support_factor_from_common_zeros(j)

        print(f"j={j}")
        print(
            f"  support factor = "
            f"{sp.factor(support)}"
        )

        for d in sorted(next(iter(DATA.values()))[j]):
            if support.subs(D, d) == 0:
                continue

            values = []

            for k in sorted(DATA):
                values.append(
                    (
                        k,
                        str(normalized_quotient(j, k, d))
                    )
                )

            print(f"  d={d}: {values}")

        print()


# ============================================================================
# 13. CROSS-k POLYNOMIALS
# ============================================================================

def interpolate_in_k(
    values: list[tuple[int, sp.Expr]]
) -> sp.Expr:

    points = [
        (sp.Integer(k), sp.sympify(v))
        for k, v in values
    ]

    return sp.factor(
        sp.interpolate(points, K)
    )


def section_6_cross_k_polynomials() -> None:
    print("=" * 78)
    print("6. CROSS-k POLYNOMIALS OF NORMALIZED CORRECTION")
    print("=" * 78)

    for j in range(1, 6):
        support = support_factor_from_common_zeros(j)

        print(f"j={j}")

        for d in sorted(next(iter(DATA.values()))[j]):
            if support.subs(D, d) == 0:
                continue

            values = [
                (
                    k,
                    normalized_quotient(j, k, d)
                )
                for k in sorted(DATA)
            ]

            poly = interpolate_in_k(values)

            print(
                f"  d={d}: "
                f"degree_in_k={sp.Poly(poly, K).degree()}"
            )

            print(
                f"       q(K)={sp.factor(poly)}"
            )

        print()


# ============================================================================
# 14. SEPARABILITY TEST
# ============================================================================

def rank_one_test(j: int) -> None:
    support = support_factor_from_common_zeros(j)

    points = []

    for k in sorted(DATA):
        for d in sorted(DATA[k][j]):
            sv = sp.expand(support.subs(D, d))

            if sv == 0:
                continue

            v = normalized_quotient(j, k, d)

            points.append((k, d, v))

    reference = None

    for item in points:
        if item[2] != 0:
            reference = item
            break

    if reference is None:
        print(f"j={j}: all normalized values zero")
        return

    k0, d0, v0 = reference

    failures = []

    for k, d, v in points:
        lhs = sp.simplify(
            v * normalized_quotient(j, k0, d0)
        )

        rhs = sp.simplify(
            normalized_quotient(j, k, d0)
            * normalized_quotient(j, k0, d)
        )

        if sp.simplify(lhs - rhs) != 0:
            failures.append((k, d))

    if failures:
        print(
            f"j={j}: NOT separable; "
            f"failures={len(failures)}"
        )
        print(
            f"  first failures={failures[:5]}"
        )
    else:
        print(
            f"j={j}: exact rank-one separation PASS"
        )


def section_7_separability() -> None:
    print("=" * 78)
    print("7. SEPARABILITY TEST")
    print("=" * 78)

    for j in range(1, 6):
        rank_one_test(j)

    print()


# ============================================================================
# 15. j=5 SPECIAL SUPPORT
# ============================================================================

def section_8_j5_special() -> None:
    print("=" * 78)
    print("8. SPECIAL j=5 SUPPORT AUDIT")
    print("=" * 78)

    support = support_factor_from_common_zeros(5)

    print(
        f"common j=5 support factor = "
        f"{sp.factor(support)}"
    )
    print()

    for k in sorted(DATA):
        print(f"k={k}")

        for d in sorted(DATA[k][5]):
            E = discrepancy(5, 5, k, d)
            sf = support.subs(D, d)

            if sf == 0:
                print(
                    f"  d={d}: "
                    f"E={E}, "
                    f"support_factor=0"
                )
            else:
                print(
                    f"  d={d}: "
                    f"E={E}, "
                    f"E/support={sp.factor(E / sf)}"
                )

        print()


# ============================================================================
# 16. FULL CORRECTION POLYNOMIALS
# ============================================================================

def section_9_full_correction_polynomials() -> None:
    print("=" * 78)
    print("9. FULL CORRECTION POLYNOMIALS")
    print("=" * 78)

    for j in range(6):
        print(f"j={j}")

        for k in sorted(DATA):
            poly, support, quotient = quotient_in_d(j, k)

            print(f"  k={k}:")
            print(f"    E(D) = {poly}")
            print(f"    support = {support}")
            print(f"    quotient = {quotient}")

        print()


# ============================================================================
# 17. FRESH-k PROJECTION
# ============================================================================

def section_10_fresh_k_projection() -> None:
    print("=" * 78)
    print("10. FRESH-k PROJECTION OF NORMALIZED CORRECTIONS")
    print("=" * 78)

    fresh_k = [15, 17, 19, 21, 23]

    for j in range(1, 6):
        support = support_factor_from_common_zeros(j)

        print(f"j={j}")

        for d in sorted(next(iter(DATA.values()))[j]):
            if support.subs(D, d) == 0:
                continue

            values = [
                (
                    k,
                    normalized_quotient(j, k, d)
                )
                for k in sorted(DATA)
            ]

            poly_k = interpolate_in_k(values)

            projections = [
                (
                    k,
                    sp.factor(poly_k.subs(K, k))
                )
                for k in fresh_k
            ]

            print(
                f"  d={d}:"
            )
            print(
                f"    q(K)={sp.factor(poly_k)}"
            )
            print(
                f"    projections={projections}"
            )

        print()


# ============================================================================
# 18. SIMPLE COMBINATORIAL SCALE SEARCH
# ============================================================================

def candidate_scales() -> list[sp.Expr]:
    return [
        sp.Integer(1),
        K + 1,
        K + 2,
        K + 3,
        K + 4,
        K + 5,
        sp.binomial(K + 4, 1),
        sp.binomial(K + 4, 2),
        sp.binomial(K + 4, 3),
        sp.binomial(K + 4, 4),
        sp.binomial(K + 4, 5),
        sp.binomial(K + 5, 1),
        sp.binomial(K + 5, 2),
        sp.binomial(K + 5, 3),
        sp.binomial(K + 5, 4),
        sp.binomial(K + 5, 5),
    ]


def section_11_simple_scale_search() -> None:
    print("=" * 78)
    print("11. SIMPLE COMBINATORIAL SCALE SEARCH")
    print("=" * 78)

    for j in range(1, 6):
        support = support_factor_from_common_zeros(j)

        print(f"j={j}")

        for d in sorted(next(iter(DATA.values()))[j]):
            if support.subs(D, d) == 0:
                continue

            values = [
                (
                    k,
                    normalized_quotient(j, k, d)
                )
                for k in sorted(DATA)
            ]

            poly = interpolate_in_k(values)

            print(f"  d={d}")
            print(f"    q(K)={sp.factor(poly)}")

            for scale in candidate_scales():
                quotient = sp.cancel(poly / scale)

                if quotient.is_polynomial(K):
                    degree = sp.Poly(
                        quotient,
                        K
                    ).degree()

                    if degree <= 2:
                        print(
                            f"    SIMPLE SCALE: "
                            f"{scale} -> "
                            f"{sp.factor(quotient)}"
                        )

        print()


# ============================================================================
# 19. COMMON ONSET
# ============================================================================

def section_12_onsets() -> None:
    print("=" * 78)
    print("12. COMMON ONSET AUDIT")
    print("=" * 78)

    for j in range(6):
        onsets = {}

        for k in sorted(DATA):
            first = None

            for d in sorted(DATA[k][j]):
                if discrepancy(5, j, k, d) != 0:
                    first = d
                    break

            onsets[k] = first

        print(f"j={j}: onsets={onsets}")

    print()


# ============================================================================
# 20. ZERO MULTIPLICITY
# ============================================================================

def section_13_zero_multiplicity() -> None:
    print("=" * 78)
    print("13. ZERO MULTIPLICITY AUDIT")
    print("=" * 78)

    for j in range(6):
        zeros = common_zero_set(j)

        print(f"j={j}")
        print(f"  common zeros={zeros}")

        for k in sorted(DATA):
            poly = discrepancy_polynomial(j, k)

            multiplicities = {}

            for z in zeros:
                multiplicities[z] = exact_root_multiplicity(
                    poly,
                    sp.Integer(z)
                )

            print(
                f"  k={k}: "
                f"multiplicities={multiplicities}"
            )

        print()


# ============================================================================
# 21. ROOT AUDIT OF THE CORRECTION
# ============================================================================

def section_14_root_audit() -> None:
    print("=" * 78)
    print("14. CORRECTION ROOT AUDIT")
    print("=" * 78)

    for j in range(6):
        print(f"j={j}")

        for k in sorted(DATA):
            poly = discrepancy_polynomial(j, k)
            factored = sp.factor(poly)

            print(
                f"  k={k}: {factored}"
            )

        print()


# ============================================================================
# 22. INTERNAL CONSISTENCY
# ============================================================================

def section_15_consistency() -> None:
    print("=" * 78)
    print("15. INTERNAL CONSISTENCY")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in range(6):
        for k in sorted(DATA):
            poly = discrepancy_polynomial(j, k)

            for d in sorted(DATA[k][j]):
                actual = discrepancy(5, j, k, d)
                predicted = sp.expand(
                    poly.subs(D, d)
                )

                tested += 1

                if sp.simplify(
                    predicted - actual
                ) != 0:
                    failures += 1

                    print(
                        "FAIL "
                        f"k={k} j={j} d={d} "
                        f"actual={actual} "
                        f"predicted={predicted}"
                    )

    print(f"tested={tested}")
    print(
        f"interpolation consistency failures="
        f"{failures}"
    )
    print()


# ============================================================================
# 23. SUPPORT-DEGREE SUMMARY
# ============================================================================

def section_16_support_summary() -> None:
    print("=" * 78)
    print("16. SUPPORT STRUCTURE SUMMARY")
    print("=" * 78)

    for j in range(6):
        support = support_factor_from_common_zeros(j)

        print(
            f"j={j}: "
            f"zeros={common_zero_set(j)} "
            f"support={sp.factor(support)} "
            f"degree={sp.Poly(support, D).degree()}"
        )

    print()


# ============================================================================
# 24. TERMINAL INTERPRETATION
# ============================================================================

def section_17_final_diagnostic() -> None:
    print("=" * 78)
    print("17. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The old r<=4 law is treated only as a baseline."
    )
    print()

    print(
        "The r=5 discrepancy is analyzed directly through "
        "its exact common support zeros."
    )
    print()

    print(
        "Important observed support:"
    )

    for j in range(6):
        print(
            f"  j={j}: "
            f"{sp.factor(support_factor_from_common_zeros(j))}"
        )

    print()

    print(
        "The main questions are:"
    )
    print(
        "  1. Which support zeros are genuine across all k?"
    )
    print(
        "  2. What is the exact multiplicity of those zeros?"
    )
    print(
        "  3. What remains after dividing them out?"
    )
    print(
        "  4. Does the quotient have a low-degree closed form?"
    )
    print(
        "  5. Can the remaining k-dependence be normalized "
        "by a simple combinatorial factor?"
    )
    print(
        "  6. Is j=5 governed by a distinct terminal-support "
        "correction?"
    )

    print()

    print(
        "No r=6 analysis is performed."
    )
    print(
        "No full pq-kernel expansion is performed."
    )
    print(
        "No replacement universal law is asserted."
    )
    print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    section_1_basic_audit()
    section_2_common_zeros()
    section_3_support_factorization()
    section_4_quotient_degree()
    section_5_normalized_values()
    section_6_cross_k_polynomials()
    section_7_separability()
    section_8_j5_special()
    section_9_full_correction_polynomials()
    section_10_fresh_k_projection()
    section_11_simple_scale_search()
    section_12_onsets()
    section_13_zero_multiplicity()
    section_14_root_audit()
    section_15_consistency()
    section_16_support_summary()
    section_17_final_diagnostic()


if __name__ == "__main__":
    main()