from __future__ import annotations

import sympy as sp


# =============================================================================
# r=5 SUPPORT-CORRECTION EXPERIMENT
# =============================================================================
#
# IMPORTANT:
# ---------
# The DATA below are the ALREADY COMPUTED EXACT DISCREPANCIES
#
#       E_j(k,D) = actual_r5(k,j,D) - old_law(k,j,D)
#
# Therefore this script MUST NOT subtract the old law again.
#
# The previous script had that structural error:
#
#     discrepancy(...) = DATA - old_candidate(...)
#
# which double-counted the old-law subtraction and corrupted the support test.
#
# This script works directly with the supplied discrepancy values.
#
# Goals:
#   1. reconstruct exact E_j(k,D) in D;
#   2. determine common exact support roots;
#   3. determine root multiplicities;
#   4. divide out common support;
#   5. examine the quotient in shifted coordinate X=D-6;
#   6. test quotient coefficient structure in k;
#   7. test Newton/binomial-basis degree in k;
#   8. search controlled combinatorial normalizations;
#   9. test j=5 terminal support separately;
#  10. search for additional exact roots outside the measured grid;
#  11. test whether support depends only on j;
#  12. reconstruct every input point exactly.
#
# Restrictions:
#   - NO r=6
#   - NO full pq-kernel expansion
#   - NO assumed replacement universal law
#
# =============================================================================


D, K, X, T = sp.symbols("D K X T")


# =============================================================================
# 1. EXACT SUPPLIED DISCREPANCY DATA
# =============================================================================

DATA: dict[int, dict[int, dict[int, sp.Integer]]] = {
    3: {
        0: {
            6: -462,
            8: -3430,
            10: -12684,
            12: -50205,
            14: -127148,
            16: -273294,
        },
        1: {
            6: 0,
            8: -1820,
            10: -15315,
            12: -44037,
            14: -107030,
            16: -244790,
        },
        2: {
            6: 0,
            8: -728,
            10: -3648,
            12: -6930,
            14: -53382,
            16: -139622,
        },
        3: {
            6: 0,
            8: 0,
            10: -525,
            12: -14763,
            14: -28014,
            16: -41706,
        },
        4: {
            6: 0,
            8: 0,
            10: -150,
            12: 1905,
            14: 26814,
            16: 0,
        },
        5: {
            6: 0,
            8: 0,
            10: 0,
            12: 0,
            14: -23800,
            16: -29892,
        },
    },

    5: {
        0: {
            6: -2268,
            8: -18648,
            10: -76125,
            12: -219150,
            14: -589372,
            16: -1272852,
        },
        1: {
            6: 0,
            8: -5712,
            10: -34440,
            12: -138068,
            14: -343728,
            16: -755160,
        },
        2: {
            6: 0,
            8: -1632,
            10: -10980,
            12: -35998,
            14: -72888,
            16: -291060,
        },
        3: {
            6: 0,
            8: 0,
            10: -855,
            12: -6255,
            14: -53916,
            16: -94050,
        },
        4: {
            6: 0,
            8: 0,
            10: -190,
            12: -1390,
            14: 8997,
            16: 98560,
        },
        5: {
            6: 0,
            8: 0,
            10: 0,
            12: 0,
            14: 0,
            16: -51800,
        },
    },

    7: {
        0: {
            6: -7524,
            8: -60984,
            10: -257136,
            12: -744996,
            14: -1818586,
            16: -4072344,
        },
        1: {
            6: 0,
            8: -13860,
            10: -83325,
            12: -282711,
            14: -794430,
            16: -1828398,
        },
        2: {
            6: 0,
            8: -3080,
            10: -20625,
            12: -78298,
            14: -158984,
            16: -359480,
        },
        3: {
            6: 0,
            8: 0,
            10: -1265,
            12: -9251,
            14: -7546,
            16: -66242,
        },
        4: {
            6: 0,
            8: 0,
            10: -230,
            12: -1682,
            14: -735,
            16: 135838,
        },
        5: {
            6: 0,
            8: 0,
            10: 0,
            12: 0,
            14: 0,
            16: 56056,
        },
    },

    9: {
        0: {
            6: -19734,
            8: -158444,
            10: -665236,
            12: -1913769,
            14: -4553505,
            16: -10125050,
        },
        1: {
            6: 0,
            8: -28600,
            10: -171600,
            12: -570570,
            14: -1298440,
            16: -3360148,
        },
        2: {
            6: 0,
            8: -5200,
            10: -34710,
            12: -130221,
            14: -189488,
            16: -177100,
        },
        3: {
            6: 0,
            8: 0,
            10: -1755,
            12: -12831,
            14: 22750,
            16: 533910,
        },
        4: {
            6: 0,
            8: 0,
            10: -270,
            12: -1974,
            14: 3360,
            16: 372780,
        },
        5: {
            6: 0,
            8: 0,
            10: 0,
            12: 0,
            14: 0,
            16: 132192,
        },
    },

    11: {
        0: {
            6: -44226,
            8: -352716,
            10: -1476384,
            12: -4138134,
            14: -9002994,
            16: -20452146,
        },
        1: {
            6: 0,
            8: -52780,
            10: -316225,
            12: -1028027,
            14: -1682590,
            16: -3115294,
        },
        2: {
            6: 0,
            8: -8120,
            10: -54075,
            12: -200417,
            14: -12390,
            16: 2135770,
        },
        3: {
            6: 0,
            8: 0,
            10: -2325,
            12: -16995,
            14: 96258,
            16: 2239830,
        },
        4: {
            6: 0,
            8: 0,
            10: -310,
            12: -2266,
            14: 10672,
            16: 1048278,
        },
        5: {
            6: 0,
            8: 0,
            10: 0,
            12: 0,
            14: 0,
            16: 277134,
        },
    },

    13: {
        0: {
            6: -88536,
            8: -702576,
            10: -2934064,
            12: -7980140,
            14: -13755924,
            16: -27091064,
        },
        1: {
            6: 0,
            8: -89760,
            10: -537200,
            12: -1702856,
            14: -891072,
            16: 7361136,
        },
        2: {
            6: 0,
            8: -11968,
            10: -79560,
            12: -291108,
            14: 609960,
            16: 11645748,
        },
        3: {
            6: 0,
            8: 0,
            10: -2975,
            12: -21743,
            14: 243984,
            16: 6970578,
        },
        4: {
            6: 0,
            8: 0,
            10: -350,
            12: -2558,
            14: 22491,
            16: 2485476,
        },
        5: {
            6: 0,
            8: 0,
            10: 0,
            12: 0,
            14: 532252,
            16: 533876,
        },
    },
}


KS = sorted(DATA)
JS = list(range(6))
DS = sorted(DATA[KS[0]][0])


# =============================================================================
# 2. BASIC DATA VALIDATION
# =============================================================================

def validate_data() -> None:
    expected_js = set(JS)
    expected_ds = set(DS)

    for k in KS:
        actual_js = set(DATA[k])
        if actual_js != expected_js:
            raise RuntimeError(
                f"k={k}: expected j={expected_js}, got {actual_js}"
            )

        for j in JS:
            actual_ds = set(DATA[k][j])
            if actual_ds != expected_ds:
                raise RuntimeError(
                    f"k={k}, j={j}: "
                    f"expected D={expected_ds}, got {actual_ds}"
                )

    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)
    print(f"k values = {KS}")
    print(f"j values = {JS}")
    print(f"D values = {DS}")
    print(f"points   = {len(KS) * len(JS) * len(DS)}")
    print()


# =============================================================================
# 3. EXACT E_j(k,D)
# =============================================================================

def E(j: int, k: int, d: int) -> sp.Integer:
    return sp.Integer(DATA[k][j][d])


# =============================================================================
# 4. EXACT D-POLYNOMIAL
# =============================================================================

def E_poly(j: int, k: int) -> sp.Expr:
    points = [
        (sp.Integer(d), E(j, k, d))
        for d in DS
    ]

    return sp.factor(
        sp.interpolate(points, D)
    )


# =============================================================================
# 5. EXACT ROOT MULTIPLICITY
# =============================================================================

def root_multiplicity(
    poly_expr: sp.Expr,
    root: int,
) -> int:
    """
    Exact multiplicity using repeated polynomial division.

    This avoids relying on a nonexistent Poly.multiplicity(...) method.
    """

    p = sp.Poly(
        sp.expand(poly_expr),
        D,
        domain="QQ",
    )

    divisor = sp.Poly(
        D - sp.Integer(root),
        D,
        domain="QQ",
    )

    multiplicity = 0

    while True:
        quotient, remainder = sp.div(
            p,
            divisor,
        )

        if not remainder.is_zero:
            return multiplicity

        multiplicity += 1
        p = quotient


# =============================================================================
# 6. COMMON ROOTS ON THE OBSERVED GRID
# =============================================================================

def common_observed_roots(j: int) -> list[int]:
    result = []

    for d in DS:
        if all(
            E(j, k, d) == 0
            for k in KS
        ):
            result.append(d)

    return result


# =============================================================================
# 7. SUPPORT FACTOR
# =============================================================================

def support_factor(j: int) -> sp.Expr:
    result = sp.Integer(1)

    for root in common_observed_roots(j):
        result *= D - root

    return sp.factor(result)


# =============================================================================
# 8. EXACT SUPPORT DIVISION
# =============================================================================

def support_quotient(j: int, k: int) -> sp.Expr:
    p = sp.Poly(
        E_poly(j, k),
        D,
        domain="QQ",
    )

    s = sp.Poly(
        support_factor(j),
        D,
        domain="QQ",
    )

    q, r = sp.div(
        p,
        s,
    )

    if not r.is_zero:
        raise RuntimeError(
            "Support does not divide exactly: "
            f"j={j}, k={k}, remainder={sp.factor(r.as_expr())}"
        )

    return sp.factor(q.as_expr())


# =============================================================================
# 9. SHIFTED QUOTIENT
# =============================================================================

def shifted_quotient(j: int, k: int) -> sp.Expr:
    return sp.factor(
        sp.expand(
            support_quotient(j, k).subs(
                D,
                X + 6,
            )
        )
    )


# =============================================================================
# 10. FINITE-DIFFERENCE DEGREE
# =============================================================================

def difference_rows(
    values: list[sp.Expr],
) -> list[list[sp.Expr]]:

    rows = [
        [sp.factor(v) for v in values]
    ]

    while len(rows[-1]) > 1:
        previous = rows[-1]

        next_row = [
            sp.factor(
                previous[i + 1]
                - previous[i]
            )
            for i in range(
                len(previous) - 1
            )
        ]

        rows.append(next_row)

    return rows


def finite_difference_degree(
    values: list[sp.Expr],
) -> int:

    rows = difference_rows(values)

    for order, row in enumerate(rows):
        if all(v == 0 for v in row):
            return order - 1

    return len(values) - 1


# =============================================================================
# 11. NEWTON/BINOMIAL BASIS
# =============================================================================

def newton_coefficients(
    values: list[sp.Expr],
) -> list[sp.Expr]:

    rows = difference_rows(values)

    return [
        sp.factor(row[0])
        for row in rows
    ]


def newton_degree(
    values: list[sp.Expr],
) -> int:

    coeffs = newton_coefficients(values)

    last = -1

    for i, c in enumerate(coeffs):
        if c != 0:
            last = i

    return last


# =============================================================================
# 12. INTERPOLATE IN k
# =============================================================================

def interpolate_k(
    values: list[tuple[int, sp.Expr]],
) -> sp.Expr:

    return sp.factor(
        sp.interpolate(
            [
                (
                    sp.Integer(k),
                    sp.sympify(v),
                )
                for k, v in values
            ],
            K,
        )
    )


# =============================================================================
# 13. DATA-POINT SUPPORT AUDIT
# =============================================================================

def section_1_support_roots() -> None:
    print("=" * 78)
    print("1. EXACT SUPPORT ROOT STRUCTURE")
    print("=" * 78)

    for j in JS:
        roots = common_observed_roots(j)
        support = support_factor(j)

        print(f"j={j}")
        print(f"  common observed roots = {roots}")
        print(f"  support = {support}")
        print(
            f"  support degree = "
            f"{sp.Poly(support, D).degree()}"
        )

        for root in roots:
            multiplicities = {
                k: root_multiplicity(
                    E_poly(j, k),
                    root,
                )
                for k in KS
            }

            print(
                f"  D={root} multiplicities="
                f"{multiplicities}"
            )

        print()


# =============================================================================
# 14. FIRST-FAILURE / ZERO PATTERN
# =============================================================================

def section_2_zero_pattern() -> None:
    print("=" * 78)
    print("2. ZERO PATTERN ON THE OBSERVED D GRID")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for k in KS:
            zeros = [
                d
                for d in DS
                if E(j, k, d) == 0
            ]

            print(
                f"  k={k}: zeros={zeros}"
            )

        print()


# =============================================================================
# 15. EXACT E POLYNOMIALS
# =============================================================================

def section_3_exact_E_polynomials() -> None:
    print("=" * 78)
    print("3. EXACT E_j(k,D) POLYNOMIALS")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for k in KS:
            print(
                f"  k={k}: "
                f"E(D)={E_poly(j, k)}"
            )

        print()


# =============================================================================
# 16. SUPPORT-REMOVED QUOTIENTS
# =============================================================================

def section_4_support_quotients() -> None:
    print("=" * 78)
    print("4. SUPPORT-REMOVED QUOTIENTS")
    print("=" * 78)

    for j in JS:
        support = support_factor(j)

        print(f"j={j}")
        print(f"  support={support}")

        for k in KS:
            print(
                f"  k={k}: "
                f"Q(D)={support_quotient(j, k)}"
            )

        print()


# =============================================================================
# 17. SHIFT D = X + 6
# =============================================================================

def section_5_shifted_quotients() -> None:
    print("=" * 78)
    print("5. SHIFTED QUOTIENTS X = D - 6")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for k in KS:
            print(
                f"  k={k}: "
                f"Q(X)={shifted_quotient(j, k)}"
            )

        print()


# =============================================================================
# 18. COEFFICIENT ARRAYS
# =============================================================================

def coefficient_vector(
    expr: sp.Expr,
    variable,
) -> list[sp.Expr]:

    p = sp.Poly(
        sp.expand(expr),
        variable,
        domain="QQ",
    )

    degree = p.degree()

    return [
        sp.factor(p.nth(i))
        for i in range(degree + 1)
    ]


def section_6_coefficient_arrays() -> None:
    print("=" * 78)
    print("6. X-COEFFICIENT ARRAYS")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            sp.Poly(
                shifted_quotient(j, k),
                X,
            ).degree()
            for k in KS
        )

        for power in range(
            max_degree + 1
        ):
            row = []

            for k in KS:
                coeffs = coefficient_vector(
                    shifted_quotient(j, k),
                    X,
                )

                value = (
                    coeffs[power]
                    if power < len(coeffs)
                    else sp.Integer(0)
                )

                row.append(
                    (k, sp.factor(value))
                )

            print(
                f"  X^{power}: {row}"
            )

        print()


# =============================================================================
# 19. FINITE DIFFERENCE DEGREE IN k
# =============================================================================

def coefficient_row(
    j: int,
    power: int,
) -> list[sp.Expr]:

    result = []

    for k in KS:
        coeffs = coefficient_vector(
            shifted_quotient(j, k),
            X,
        )

        result.append(
            coeffs[power]
            if power < len(coeffs)
            else sp.Integer(0)
        )

    return result


def section_7_k_degree() -> None:
    print("=" * 78)
    print("7. FINITE-DIFFERENCE DEGREE IN k")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            sp.Poly(
                shifted_quotient(j, k),
                X,
            ).degree()
            for k in KS
        )

        for power in range(
            max_degree + 1
        ):
            values = coefficient_row(
                j,
                power,
            )

            deg = finite_difference_degree(
                values
            )

            print(
                f"  X^{power}: "
                f"degree_in_k={deg}"
            )

        print()


# =============================================================================
# 20. NEWTON/BINOMIAL COEFFICIENTS
# =============================================================================

def section_8_newton_basis() -> None:
    print("=" * 78)
    print("8. NEWTON / BINOMIAL BASIS IN k")
    print("=" * 78)

    k0 = KS[0]
    h = KS[1] - KS[0]

    print(
        f"k-grid: k = {k0} + {h}*n"
    )
    print()

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            sp.Poly(
                shifted_quotient(j, k),
                X,
            ).degree()
            for k in KS
        )

        for power in range(
            max_degree + 1
        ):
            values = coefficient_row(
                j,
                power,
            )

            coeffs = newton_coefficients(
                values
            )

            print(
                f"  X^{power}: "
                f"newton_coeffs={coeffs}"
            )

            print(
                f"    newton_degree="
                f"{newton_degree(values)}"
            )

        print()


# =============================================================================
# 21. CONTROLLED NORMALIZATIONS
# =============================================================================

def normalization_candidates() -> list[tuple[str, sp.Expr]]:
    candidates = [
        ("1", sp.Integer(1)),
        ("K+1", K + 1),
        ("K+2", K + 2),
        ("K+3", K + 3),
        ("K+4", K + 4),
        ("K+5", K + 5),

        ("C(K+4,1)", sp.binomial(K + 4, 1)),
        ("C(K+4,2)", sp.binomial(K + 4, 2)),
        ("C(K+4,3)", sp.binomial(K + 4, 3)),
        ("C(K+4,4)", sp.binomial(K + 4, 4)),
        ("C(K+4,5)", sp.binomial(K + 4, 5)),

        ("C(K+5,1)", sp.binomial(K + 5, 1)),
        ("C(K+5,2)", sp.binomial(K + 5, 2)),
        ("C(K+5,3)", sp.binomial(K + 5, 3)),
        ("C(K+5,4)", sp.binomial(K + 5, 4)),
        ("C(K+5,5)", sp.binomial(K + 5, 5)),
    ]

    # Add shifted binomial families around K-j and K+j.
    for shift_name, shift in [
        ("K-1", -1),
        ("K", 0),
        ("K+1", 1),
        ("K+2", 2),
        ("K+3", 3),
        ("K+4", 4),
        ("K+5", 5),
    ]:
        base = K + shift

        for r in range(1, 6):
            if r <= 5:
                candidates.append(
                    (
                        f"C({shift_name},{r})",
                        sp.binomial(base, r),
                    )
                )

    # Deduplicate.
    result = []
    seen = set()

    for name, expr in candidates:
        expr = sp.factor(
            sp.expand_func(expr)
        )

        key = sp.srepr(expr)

        if key not in seen:
            seen.add(key)
            result.append(
                (
                    name,
                    expr,
                )
            )

    return result


# =============================================================================
# 22. NORMALIZED COEFFICIENT VALUES
# =============================================================================

def normalized_row(
    j: int,
    power: int,
    scale_expr: sp.Expr,
) -> list[sp.Expr] | None:

    values = []

    for k in KS:

        coeffs = coefficient_vector(
            shifted_quotient(j, k),
            X,
        )

        value = (
            coeffs[power]
            if power < len(coeffs)
            else sp.Integer(0)
        )

        scale_value = sp.sympify(
            scale_expr
        ).subs(
            K,
            sp.Integer(k),
        )

        if scale_value == 0:
            return None

        values.append(
            sp.factor(
                sp.cancel(
                    value / scale_value
                )
            )
        )

    return values


# =============================================================================
# 23. NORMALIZATION SEARCH
# =============================================================================

def section_9_normalization_search() -> None:
    print("=" * 78)
    print("9. CONTROLLED COMBINATORIAL NORMALIZATION SEARCH")
    print("=" * 78)

    candidates = normalization_candidates()

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            sp.Poly(
                shifted_quotient(j, k),
                X,
            ).degree()
            for k in KS
        )

        for power in range(
            max_degree + 1
        ):

            raw_values = coefficient_row(
                j,
                power,
            )

            raw_degree = finite_difference_degree(
                raw_values
            )

            improvements = []

            for name, scale in candidates:

                values = normalized_row(
                    j,
                    power,
                    scale,
                )

                if values is None:
                    continue

                degree = finite_difference_degree(
                    values
                )

                if degree < raw_degree:
                    poly = interpolate_k(
                        list(
                            zip(
                                KS,
                                values,
                            )
                        )
                    )

                    improvements.append(
                        (
                            degree,
                            name,
                            sp.factor(poly),
                        )
                    )

            improvements.sort(
                key=lambda t: (
                    t[0],
                    t[1],
                )
            )

            print(
                f"  X^{power}: "
                f"raw_degree={raw_degree}"
            )

            if not improvements:
                print(
                    "    no tested normalization "
                    "reduces degree"
                )
            else:
                for degree, name, poly in improvements[:10]:
                    print(
                        f"    improvement: "
                        f"degree={degree}, "
                        f"scale={name}"
                    )
                    print(
                        f"      normalized polynomial="
                        f"{poly}"
                    )

        print()


# =============================================================================
# 24. NORMALIZATION BY COMMON BINOMIAL FACTORS
# =============================================================================

def section_10_targeted_binomial_normalization() -> None:
    print("=" * 78)
    print("10. TARGETED BINOMIAL NORMALIZATION")
    print("=" * 78)

    targets = {
        0: [
            ("C(K+4,4)", sp.binomial(K + 4, 4)),
            ("C(K+5,5)", sp.binomial(K + 5, 5)),
        ],
        1: [
            ("C(K+4,3)", sp.binomial(K + 4, 3)),
            ("C(K+4,4)", sp.binomial(K + 4, 4)),
            ("C(K+5,4)", sp.binomial(K + 5, 4)),
        ],
        2: [
            ("C(K+4,2)", sp.binomial(K + 4, 2)),
            ("C(K+4,3)", sp.binomial(K + 4, 3)),
            ("C(K+5,3)", sp.binomial(K + 5, 3)),
        ],
        3: [
            ("C(K+4,2)", sp.binomial(K + 4, 2)),
            ("C(K+5,2)", sp.binomial(K + 5, 2)),
        ],
        4: [
            ("K+3", K + 3),
            ("K+4", K + 4),
            ("C(K+4,1)", sp.binomial(K + 4, 1)),
        ],
        5: [
            ("1", sp.Integer(1)),
            ("K+5", K + 5),
        ],
    }

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            sp.Poly(
                shifted_quotient(j, k),
                X,
            ).degree()
            for k in KS
        )

        for power in range(
            max_degree + 1
        ):

            raw = coefficient_row(
                j,
                power,
            )

            raw_degree = finite_difference_degree(
                raw
            )

            print(
                f"  X^{power}: "
                f"raw_degree={raw_degree}"
            )

            for name, scale in targets[j]:

                values = normalized_row(
                    j,
                    power,
                    scale,
                )

                if values is None:
                    continue

                deg = finite_difference_degree(
                    values
                )

                poly = interpolate_k(
                    list(
                        zip(
                            KS,
                            values,
                        )
                    )
                )

                print(
                    f"    {name}: "
                    f"degree={deg}"
                )
                print(
                    f"      {poly}"
                )

        print()


# =============================================================================
# 25. J=5 TERMINAL SUPPORT
# =============================================================================

def section_11_j5_terminal() -> None:
    print("=" * 78)
    print("11. J=5 TERMINAL SUPPORT AUDIT")
    print("=" * 78)

    j = 5

    support = support_factor(j)

    print(
        f"common support={support}"
    )
    print()

    for k in KS:
        polynomial = E_poly(j, k)
        quotient = support_quotient(j, k)

        print(f"k={k}")
        print(
            f"  E(D)={polynomial}"
        )
        print(
            f"  Q(D)={quotient}"
        )

        for d in DS:
            value = E(j, k, d)

            print(
                f"    D={d}: E={value}"
            )

        print()


# =============================================================================
# 26. J=5 TEST FOR ADDITIONAL D=14 ROOT
# =============================================================================

def section_12_j5_D14() -> None:
    print("=" * 78)
    print("12. J=5 D=14 ROOT TEST")
    print("=" * 78)

    j = 5

    for k in KS:
        Q = support_quotient(j, k)

        q14 = sp.factor(
            Q.subs(
                D,
                14,
            )
        )

        m = (
            root_multiplicity(
                E_poly(j, k),
                14,
            )
            if q14 == 0
            else 0
        )

        print(
            f"k={k}: "
            f"Q(14)={q14} "
            f"E multiplicity at D=14={m}"
        )

    print()


# =============================================================================
# 27. SEARCH FOR COMMON INTEGER ROOTS OUTSIDE GRID
# =============================================================================

def section_13_extra_integer_roots() -> None:
    print("=" * 78)
    print("13. COMMON INTEGER ROOT SEARCH")
    print("=" * 78)

    search_range = list(
        range(
            -4,
            41,
            2,
        )
    )

    for j in JS:
        roots = []

        for d in search_range:

            if all(
                E_poly(j, k).subs(D, d) == 0
                for k in KS
            ):
                roots.append(d)

        print(
            f"j={j}: common integer roots "
            f"in search range={roots}"
        )

    print()


# =============================================================================
# 28. ROOT STRUCTURE BY PAIRS OF k
# =============================================================================

def section_14_pairwise_support() -> None:
    print("=" * 78)
    print("14. PAIRWISE k SUPPORT INTERSECTION")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        pairwise = {}

        for i, k1 in enumerate(KS):
            for k2 in KS[i + 1:]:

                common = []

                for d in DS:
                    if (
                        E(j, k1, d) == 0
                        and
                        E(j, k2, d) == 0
                    ):
                        common.append(d)

                pairwise[(k1, k2)] = common

        for pair, roots in pairwise.items():
            print(
                f"  k={pair}: common zeros={roots}"
            )

        print()


# =============================================================================
# 29. SUPPORT STABILITY BY OMITTING ONE k
# =============================================================================

def section_15_leave_one_k_out() -> None:
    print("=" * 78)
    print("15. LEAVE-ONE-k-OUT SUPPORT STABILITY")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for omitted_k in KS:

            retained = [
                k
                for k in KS
                if k != omitted_k
            ]

            common = [
                d
                for d in DS
                if all(
                    E(j, k, d) == 0
                    for k in retained
                )
            ]

            print(
                f"  omit k={omitted_k}: "
                f"common zeros={common}"
            )

        print()


# =============================================================================
# 30. SUPPORT-VANISHING ORDER
# =============================================================================

def section_16_vanishing_order() -> None:
    print("=" * 78)
    print("16. VANISHING-ORDER PROFILE")
    print("=" * 78)

    for j in JS:

        roots = common_observed_roots(j)

        print(
            f"j={j}: roots={roots}"
        )

        for root in roots:

            profile = []

            for k in KS:
                profile.append(
                    (
                        k,
                        root_multiplicity(
                            E_poly(j, k),
                            root,
                        ),
                    )
                )

            print(
                f"  D={root}: {profile}"
            )

        print()


# =============================================================================
# 31. ZERO-SET LATTICE
# =============================================================================

def section_17_zero_lattice() -> None:
    print("=" * 78)
    print("17. ZERO-SET LATTICE")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        matrix = []

        for k in KS:
            row = []

            for d in DS:
                row.append(
                    1
                    if E(j, k, d) == 0
                    else 0
                )

            matrix.append(row)

            print(
                f"  k={k}: {row}"
            )

        print()


# =============================================================================
# 32. QUOTIENT ROOT AUDIT
# =============================================================================

def section_18_quotient_roots() -> None:
    print("=" * 78)
    print("18. SUPPORT-QUOTIENT ROOT AUDIT")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for k in KS:

            Q = support_quotient(j, k)

            exact_integer_roots = []

            for d in range(
                -20,
                51,
                2,
            ):
                if Q.subs(D, d) == 0:
                    exact_integer_roots.append(d)

            print(
                f"  k={k}: "
                f"integer roots of Q(D)="
                f"{exact_integer_roots}"
            )

        print()


# =============================================================================
# 33. CROSS-k POLYNOMIALS FOR FIXED D
# =============================================================================

def section_19_fixed_D_k_polynomials() -> None:
    print("=" * 78)
    print("19. FIXED-D CROSS-k POLYNOMIALS")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for d in DS:

            values = [
                (
                    k,
                    E(j, k, d)
                )
                for k in KS
            ]

            Pk = interpolate_k(
                values
            )

            degree = sp.Poly(
                Pk,
                K,
            ).degree()

            print(
                f"  D={d}: "
                f"degree_in_k={degree}"
            )
            print(
                f"    E(K)={Pk}"
            )

        print()


# =============================================================================
# 34. CROSS-k POLYNOMIALS FOR NORMALIZED QUOTIENTS
# =============================================================================

def section_20_fixed_X_k_polynomials() -> None:
    print("=" * 78)
    print("20. FIXED-X CROSS-k POLYNOMIALS OF QUOTIENT")
    print("=" * 78)

    for j in JS:

        max_degree = max(
            sp.Poly(
                shifted_quotient(j, k),
                X,
            ).degree()
            for k in KS
        )

        print(
            f"j={j}"
        )

        for power in range(
            max_degree + 1
        ):

            values = coefficient_row(
                j,
                power,
            )

            Pk = interpolate_k(
                list(
                    zip(
                        KS,
                        values,
                    )
                )
            )

            degree = sp.Poly(
                Pk,
                K,
            ).degree()

            print(
                f"  X^{power}: "
                f"degree_in_k={degree}"
            )
            print(
                f"    {Pk}"
            )

        print()


# =============================================================================
# 35. LEADING-D COEFFICIENT AUDIT
# =============================================================================

def section_21_leading_D_coefficients() -> None:
    print("=" * 78)
    print("21. LEADING-D COEFFICIENT AUDIT")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for k in KS:

            P = sp.Poly(
                E_poly(j, k),
                D,
            )

            degree = P.degree()
            lead = sp.factor(
                P.LC()
            )

            print(
                f"  k={k}: "
                f"degree={degree}, "
                f"leading={lead}"
            )

        print()


# =============================================================================
# 36. SUPPORT DEGREE VERSUS j
# =============================================================================

def section_22_support_degree_pattern() -> None:
    print("=" * 78)
    print("22. SUPPORT DEGREE PATTERN")
    print("=" * 78)

    degrees = []

    for j in JS:
        deg = sp.Poly(
            support_factor(j),
            D,
        ).degree()

        degrees.append(
            (
                j,
                deg,
            )
        )

    print(
        f"support degrees={degrees}"
    )

    # Simple candidate patterns for diagnostic comparison only.
    candidates = {
        "floor(j/2)": [
            j // 2
            for j in JS
        ],
        "floor((j+1)/2)": [
            (j + 1) // 2
            for j in JS
        ],
        "min(2,j)": [
            min(2, j)
            for j in JS
        ],
    }

    actual = [
        deg
        for _, deg in degrees
    ]

    for name, values in candidates.items():
        print(
            f"  {name}: "
            f"{values} "
            f"match={values == actual}"
        )

    print()


# =============================================================================
# 37. SUPPORT ROOT PATTERN VERSUS j
# =============================================================================

def section_23_support_roots_vs_j() -> None:
    print("=" * 78)
    print("23. SUPPORT ROOT PATTERN VERSUS j")
    print("=" * 78)

    for j in JS:
        print(
            f"j={j}: "
            f"roots={common_observed_roots(j)}"
        )

    print()

    # Candidate simple ladder roots.
    ladders = {
        "even_from_6": [
            6 + 2*t
            for t in range(5)
        ],
        "6_through_6+2j": {
            j: [
                6 + 2*t
                for t in range(j + 1)
            ]
            for j in JS
        },
    }

    print(
        "candidate even ladder:"
    )
    print(
        ladders["even_from_6"]
    )

    print()

    for j in JS:
        actual = common_observed_roots(j)

        if j in ladders["6_through_6+2j"]:
            candidate = ladders[
                "6_through_6+2j"
            ][j]

            print(
                f"j={j}: "
                f"actual={actual}, "
                f"candidate={candidate}"
            )

    print()


# =============================================================================
# 38. EXACT SUPPORT RECONSTRUCTION
# =============================================================================

def section_24_support_reconstruction() -> None:
    print("=" * 78)
    print("24. EXACT SUPPORT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in JS:
        S = support_factor(j)

        for k in KS:

            Efull = E_poly(
                j,
                k,
            )

            Q = support_quotient(
                j,
                k,
            )

            residual = sp.factor(
                sp.expand(
                    Efull
                    -
                    sp.expand(S * Q)
                )
            )

            tested += 1

            if residual != 0:
                failures += 1

                print(
                    f"FAIL "
                    f"j={j} k={k}: "
                    f"residual={residual}"
                )

    print(
        f"tested={tested}"
    )
    print(
        f"failures={failures}"
    )
    print()


# =============================================================================
# 39. EXACT INTERPOLATION CHECK
# =============================================================================

def section_25_interpolation_check() -> None:
    print("=" * 78)
    print("25. EXACT INTERPOLATION CHECK")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in JS:
        for k in KS:

            P = E_poly(
                j,
                k,
            )

            for d in DS:

                actual = E(
                    j,
                    k,
                    d,
                )

                predicted = sp.factor(
                    P.subs(
                        D,
                        d,
                    )
                )

                tested += 1

                if predicted != actual:
                    failures += 1

                    print(
                        f"FAIL "
                        f"j={j} "
                        f"k={k} "
                        f"D={d} "
                        f"actual={actual} "
                        f"predicted={predicted}"
                    )

    print(
        f"tested={tested}"
    )
    print(
        f"interpolation failures={failures}"
    )
    print()


# =============================================================================
# 40. DIRECT SUPPORT DIVISIBILITY
# =============================================================================

def section_26_symbolic_divisibility() -> None:
    print("=" * 78)
    print("26. SYMBOLIC SUPPORT DIVISIBILITY")
    print("=" * 78)

    for j in JS:
        S = sp.Poly(
            support_factor(j),
            D,
            domain="QQ",
        )

        print(f"j={j}")
        print(
            f"  support={sp.factor(S.as_expr())}"
        )

        for k in KS:

            P = sp.Poly(
                E_poly(j, k),
                D,
                domain="QQ",
            )

            q, r = sp.div(
                P,
                S,
            )

            print(
                f"  k={k}: "
                f"remainder={sp.factor(r.as_expr())}"
            )

        print()


# =============================================================================
# 41. COMPACT SUMMARY
# =============================================================================

def section_27_summary() -> None:
    print("=" * 78)
    print("27. COMPACT STRUCTURAL SUMMARY")
    print("=" * 78)

    for j in JS:

        roots = common_observed_roots(j)
        S = support_factor(j)

        q_degrees = [
            sp.Poly(
                support_quotient(j, k),
                D,
            ).degree()
            for k in KS
        ]

        print(
            f"j={j}: "
            f"roots={roots} "
            f"support={sp.factor(S)} "
            f"support_degree="
            f"{sp.Poly(S,D).degree()} "
            f"quotient_degrees={q_degrees}"
        )

    print()


# =============================================================================
# 42. FINAL INTERPRETATION
# =============================================================================

def section_28_final_diagnostic() -> None:
    print("=" * 78)
    print("28. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "THE INPUT DATA ARE TREATED AS ALREADY COMPUTED DISCREPANCIES."
    )
    print()

    print(
        "Therefore:"
    )
    print(
        "    E_j(k,D) = supplied discrepancy"
    )
    print()

    print(
        "No old-law subtraction is performed anywhere in this script."
    )
    print()

    print(
        "The structural decomposition being tested is:"
    )
    print()
    print(
        "    E_j(k,D) = Support_j(D) * Q_j(k,D)"
    )
    print()

    print(
        "Support is inferred only from exact common zeros across all tested k."
    )
    print()

    print(
        "Multiplicity is computed by exact repeated polynomial division."
    )
    print()

    print(
        "The quotient is then studied in X = D - 6."
    )
    print()

    print(
        "The k-structure is tested in two independent ways:"
    )
    print(
        "    (1) ordinary polynomial degree in k;"
    )
    print(
        "    (2) Newton/binomial finite-difference degree."
    )
    print()

    print(
        "A degree-5 interpolation across six k-values is NOT treated "
        "as evidence of a genuine degree-5 law."
    )
    print()

    print(
        "Controlled combinatorial normalizations are accepted only when "
        "they reduce the exact finite-difference degree."
    )
    print()

    print(
        "j=5 is analyzed separately because its terminal behaviour may "
        "come from a distinct support boundary."
    )
    print()

    print(
        "No r=6."
    )
    print(
        "No full pq-kernel expansion."
    )
    print(
        "No replacement universal law."
    )
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    validate_data()

    section_1_support_roots()
    section_2_zero_pattern()
    section_3_exact_E_polynomials()
    section_4_support_quotients()
    section_5_shifted_quotients()
    section_6_coefficient_arrays()
    section_7_k_degree()
    section_8_newton_basis()
    section_9_normalization_search()
    section_10_targeted_binomial_normalization()
    section_11_j5_terminal()
    section_12_j5_D14()
    section_13_extra_integer_roots()
    section_14_pairwise_support()
    section_15_leave_one_k_out()
    section_16_vanishing_order()
    section_17_zero_lattice()
    section_18_quotient_roots()
    section_19_fixed_D_k_polynomials()
    section_20_fixed_X_k_polynomials()
    section_21_leading_D_coefficients()
    section_22_support_degree_pattern()
    section_23_support_roots_vs_j()
    section_24_support_reconstruction()
    section_25_interpolation_check()
    section_26_symbolic_divisibility()
    section_27_summary()
    section_28_final_diagnostic()


if __name__ == "__main__":
    main()

