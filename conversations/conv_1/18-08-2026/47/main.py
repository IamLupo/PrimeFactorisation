from __future__ import annotations

import sympy as sp


# ============================================================================
# EXPERIMENT 255
# EXACT r=5 SUPPORT CORRECTION / BINOMIAL-BASIS TEST
# ============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No filesystem access
# No previous experiment imported
# No r=6
# No full pq-kernel expansion
#
# PURPOSE
# -------
# The previous experiments established:
#
#   j=1 : common zero D=6
#   j=2 : common zero D=6
#   j=3 : common zeros D=6,8
#   j=4 : common zeros D=6,8
#   j=5 : common zeros D=6,8,10,12
#
# Every detected common root has multiplicity 1.
#
# The previous k-interpolation produced degree-5 polynomials simply because
# six k-values were interpolated. That is not evidence of a true degree-5 law.
#
# This experiment therefore does NOT blindly interpolate a new universal law.
#
# It instead:
#
#   1. reconstructs every exact discrepancy polynomial in D;
#   2. confirms common support roots and multiplicities;
#   3. removes the support exactly;
#   4. expresses the quotient in the shifted variable x = D - 6;
#   5. extracts quotient coefficient rows by k;
#   6. converts those coefficient rows to a Newton/binomial basis in k;
#   7. tests whether the coefficient rows have low finite-difference degree;
#   8. searches a controlled family of combinatorial normalizations;
#   9. tests whether the normalized quotient becomes low-degree in k;
#  10. treats j=5 separately;
#  11. reports exact candidate structures without claiming a theorem.
#
# IMPORTANT
# ---------
# All numerical data below are the supplied exact r=5 discrepancy values.
#
# ============================================================================


D, K, X = sp.symbols("D K X")


# ============================================================================
# 1. EXACT DATA
# ============================================================================

DATA = {
    3: {
        0: {6: -462, 8: -3430, 10: -12684, 12: -50205, 14: -127148, 16: -273294},
        1: {6: 0, 8: -1820, 10: -15315, 12: -44037, 14: -107030, 16: -244790},
        2: {6: 0, 8: -728, 10: -3648, 12: -6930, 14: -53382, 16: -139622},
        3: {6: 0, 8: 0, 10: -525, 12: -14763, 14: -28014, 16: -41706},
        4: {6: 0, 8: 0, 10: -150, 12: 1905, 14: 26814, 16: 0},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: -23800, 16: -29892},
    },

    5: {
        0: {6: -2268, 8: -18648, 10: -76125, 12: -219150, 14: -589372, 16: -1272852},
        1: {6: 0, 8: -5712, 10: -34440, 12: -138068, 14: -343728, 16: -755160},
        2: {6: 0, 8: -1632, 10: -10980, 12: -35998, 14: -72888, 16: -291060},
        3: {6: 0, 8: 0, 10: -855, 12: -6255, 14: -53916, 16: -94050},
        4: {6: 0, 8: 0, 10: -190, 12: -1390, 14: 8997, 16: 98560},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: -51800},
    },

    7: {
        0: {6: -7524, 8: -60984, 10: -257136, 12: -744996, 14: -1818586, 16: -4072344},
        1: {6: 0, 8: -13860, 10: -83325, 12: -282711, 14: -794430, 16: -1828398},
        2: {6: 0, 8: -3080, 10: -20625, 12: -78298, 14: -158984, 16: -359480},
        3: {6: 0, 8: 0, 10: -1265, 12: -9251, 14: -7546, 16: -66242},
        4: {6: 0, 8: 0, 10: -230, 12: -1682, 14: -735, 16: 135838},
        5: {6: 0, 8: 0, 10: 0, 12: -2, 14: -196, 16: 56056},
    },

    9: {
        0: {6: -19734, 8: -158444, 10: -665236, 12: -1913769, 14: -4553505, 16: -10125050},
        1: {6: 0, 8: -28600, 10: -171600, 12: -570570, 14: -1298440, 16: -3360148},
        2: {6: 0, 8: -5200, 10: -34710, 12: -130221, 14: -189488, 16: -177100},
        3: {6: 0, 8: 0, 10: -1755, 12: -12831, 14: 22750, 16: 533910},
        4: {6: 0, 8: 0, 10: -270, 12: -1974, 14: 3360, 16: 372780},
        5: {6: 0, 8: 0, 10: 0, 12: -2, 14: -196, 16: 132192},
    },

    11: {
        0: {6: -44226, 8: -352716, 10: -1476384, 12: -4138134, 14: -9002994, 16: -20452146},
        1: {6: 0, 8: -52780, 10: -316225, 12: -1028027, 14: -1682590, 16: -3115294},
        2: {6: 0, 8: -8120, 10: -54075, 12: -200417, 14: -12390, 16: 2135770},
        3: {6: 0, 8: 0, 10: -2325, 12: -16995, 14: 96258, 16: 2239830},
        4: {6: 0, 8: 0, 10: -310, 12: -2266, 14: 10672, 16: 1048278},
        5: {6: 0, 8: 0, 10: 0, 12: -2, 14: -196, 16: 277134},
    },

    13: {
        0: {6: -88536, 8: -702576, 10: -2934064, 12: -7980140, 14: -13755924, 16: -27091064},
        1: {6: 0, 8: -89760, 10: -537200, 12: -1702856, 14: -891072, 16: 7361136},
        2: {6: 0, 8: -11968, 10: -79560, 12: -291108, 14: 609960, 16: 11645748},
        3: {6: 0, 8: 0, 10: -2975, 12: -21743, 14: 243984, 16: 6970578},
        4: {6: 0, 8: 0, 10: -350, 12: -2558, 14: 22491, 16: 2485476},
        5: {6: 0, 8: 0, 10: 0, 12: -2, 14: 532252, 16: 533876},
    },
}


KS = sorted(DATA)
JS = range(6)
DS = sorted(DATA[KS[0]][0])


# ============================================================================
# 2. OLD UNIVERSAL LAW
# ============================================================================

def old_candidate(r: int, j: int, k: int, d) -> sp.Expr:
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

    expr *= d - sp.Rational((r - j) * k, k + j)

    return sp.factor(expr)


def discrepancy(j: int, k: int, d: int) -> sp.Expr:
    return sp.factor(
        sp.Integer(DATA[k][j][d])
        - old_candidate(5, j, k, d)
    )


# ============================================================================
# 3. INTERPOLATION IN D
# ============================================================================

def discrepancy_poly(j: int, k: int) -> sp.Expr:
    pts = [
        (sp.Integer(d), discrepancy(j, k, d))
        for d in sorted(DATA[k][j])
    ]

    return sp.factor(
        sp.interpolate(pts, D)
    )


# ============================================================================
# 4. EXACT ROOT MULTIPLICITY
# ============================================================================

def root_multiplicity(expr: sp.Expr, root) -> int | sp.oo:
    expr = sp.factor(sp.expand(expr))

    if expr == 0:
        return sp.oo

    p = sp.Poly(expr, D)
    divisor = sp.Poly(D - sp.sympify(root), D)

    m = 0

    while True:
        q, r = sp.div(p, divisor)

        if not r.is_zero:
            break

        m += 1
        p = q

        if p.is_zero:
            return sp.oo

    return m


# ============================================================================
# 5. COMMON ROOTS
# ============================================================================

def common_roots(j: int) -> list[int]:
    result = []

    for d in DS:
        if all(
            discrepancy(j, k, d) == 0
            for k in KS
        ):
            result.append(d)

    return result


def support_factor(j: int) -> sp.Expr:
    result = sp.Integer(1)

    for z in common_roots(j):
        result *= D - z

    return sp.expand(result)


# ============================================================================
# 6. EXACT SUPPORT QUOTIENT
# ============================================================================

def quotient_poly(j: int, k: int) -> sp.Expr:
    E = sp.Poly(
        discrepancy_poly(j, k),
        D
    )

    S = sp.Poly(
        support_factor(j),
        D
    )

    q, r = sp.div(E, S)

    if not r.is_zero:
        raise RuntimeError(
            f"Non-exact support division at j={j}, k={k}: "
            f"remainder={sp.factor(r.as_expr())}"
        )

    return sp.factor(q.as_expr())


# ============================================================================
# 7. SHIFT D -> X + 6
# ============================================================================

def shifted_quotient(j: int, k: int) -> sp.Expr:
    Q = quotient_poly(j, k)

    return sp.factor(
        sp.expand(
            Q.subs(D, X + 6)
        )
    )


# ============================================================================
# 8. COEFFICIENT VECTOR IN X
# ============================================================================

def x_coefficients(j: int, k: int) -> list[sp.Expr]:
    Q = sp.Poly(
        shifted_quotient(j, k),
        X
    )

    degree = Q.degree()

    return [
        sp.factor(Q.nth(i))
        for i in range(degree + 1)
    ]


# ============================================================================
# 9. FINITE DIFFERENCES IN k
# ============================================================================

def finite_difference_table(
    values: list[sp.Expr]
) -> list[list[sp.Expr]]:

    rows = [[sp.sympify(v) for v in values]]

    while len(rows[-1]) > 1:
        prev = rows[-1]

        nxt = [
            sp.factor(
                prev[i + 1] - prev[i]
            )
            for i in range(len(prev) - 1)
        ]

        rows.append(nxt)

    return rows


def finite_difference_degree(
    values: list[sp.Expr]
) -> int:

    if len(values) <= 1:
        return 0

    rows = finite_difference_table(values)

    for degree, row in enumerate(rows):
        if all(v == 0 for v in row):
            return degree - 1

    return len(values) - 1


# ============================================================================
# 10. NEWTON / BINOMIAL BASIS
# ============================================================================

def newton_binomial_coefficients(
    values: list[sp.Expr]
) -> list[sp.Expr]:
    """
    Values are assumed to lie at equally spaced k-values.

    For k = k0 + h*n, the Newton expansion is

        f(k0 + h*n)
          = sum_m c_m * binomial(n,m)

    where c_m = Delta^m f(0).
    """

    table = finite_difference_table(values)

    return [
        sp.factor(row[0])
        for row in table
    ]


def reconstruct_from_binomial_basis(
    coeffs: list[sp.Expr],
    n: sp.Expr
) -> sp.Expr:

    result = sp.Integer(0)

    for m, c in enumerate(coeffs):
        result += c * sp.binomial(n, m)

    return sp.factor(result)


# ============================================================================
# 11. DIRECT K POLYNOMIAL
# ============================================================================

def interpolate_k(
    values: list[tuple[int, sp.Expr]]
) -> sp.Expr:

    return sp.factor(
        sp.interpolate(
            [
                (sp.Integer(k), sp.sympify(v))
                for k, v in values
            ],
            K
        )
    )


# ============================================================================
# 12. SECTION 1 -- SUPPORT ROOT STRUCTURE
# ============================================================================

def section_1_support_structure() -> None:
    print("=" * 78)
    print("1. SUPPORT ROOT STRUCTURE")
    print("=" * 78)

    for j in JS:
        roots = common_roots(j)
        support = support_factor(j)

        print(f"j={j}")
        print(f"  common roots = {roots}")
        print(f"  support     = {sp.factor(support)}")
        print(
            f"  degree      = "
            f"{sp.Poly(support, D).degree()}"
        )

        for root in roots:
            mults = {
                k: root_multiplicity(
                    discrepancy_poly(j, k),
                    root
                )
                for k in KS
            }

            print(
                f"  root D={root} multiplicities={mults}"
            )

        print()


# ============================================================================
# 13. SECTION 2 -- SUPPORT QUOTIENT DEGREE
# ============================================================================

def section_2_quotient_degrees() -> None:
    print("=" * 78)
    print("2. SUPPORT-QUOTIENT DEGREE")
    print("=" * 78)

    for j in JS:
        support = support_factor(j)
        support_degree = sp.Poly(
            support,
            D
        ).degree()

        print(
            f"j={j} "
            f"support_degree={support_degree}"
        )

        for k in KS:
            E = discrepancy_poly(j, k)
            Q = quotient_poly(j, k)

            print(
                f"  k={k}: "
                f"E_degree={sp.Poly(E, D).degree()} "
                f"Q_degree={sp.Poly(Q, D).degree()}"
            )

        print()


# ============================================================================
# 14. SECTION 3 -- QUOTIENTS IN X=D-6
# ============================================================================

def section_3_shifted_quotients() -> None:
    print("=" * 78)
    print("3. SUPPORT-REMOVED QUOTIENTS IN X = D - 6")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        for k in KS:
            print(
                f"  k={k}: "
                f"Q(X)={shifted_quotient(j, k)}"
            )

        print()


# ============================================================================
# 15. SECTION 4 -- COEFFICIENT ARRAYS
# ============================================================================

def section_4_coefficient_arrays() -> None:
    print("=" * 78)
    print("4. X-COEFFICIENT ARRAYS")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            len(x_coefficients(j, k))
            for k in KS
        ) - 1

        for power in range(max_degree + 1):
            row = []

            for k in KS:
                coeffs = x_coefficients(j, k)

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


# ============================================================================
# 16. SECTION 5 -- FINITE-DIFFERENCE DEGREE IN k
# ============================================================================

def section_5_k_difference_degree() -> None:
    print("=" * 78)
    print("5. FINITE-DIFFERENCE DEGREE IN k")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            len(x_coefficients(j, k))
            for k in KS
        ) - 1

        for power in range(max_degree + 1):

            values = []

            for k in KS:
                coeffs = x_coefficients(j, k)

                values.append(
                    coeffs[power]
                    if power < len(coeffs)
                    else sp.Integer(0)
                )

            deg = finite_difference_degree(
                values
            )

            print(
                f"  X^{power}: "
                f"k-difference-degree={deg}"
            )

        print()


# ============================================================================
# 17. SECTION 6 -- NEWTON/BINOMIAL COEFFICIENTS
# ============================================================================

def section_6_newton_basis() -> None:
    print("=" * 78)
    print("6. NEWTON / BINOMIAL-BASIS COEFFICIENTS IN k")
    print("=" * 78)

    k0 = KS[0]
    step = KS[1] - KS[0]

    print(
        f"k-grid: k = {k0} + {step}*n"
    )
    print()

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            len(x_coefficients(j, k))
            for k in KS
        ) - 1

        for power in range(max_degree + 1):

            values = []

            for k in KS:
                coeffs = x_coefficients(j, k)

                values.append(
                    coeffs[power]
                    if power < len(coeffs)
                    else sp.Integer(0)
                )

            delta_coeffs = newton_binomial_coefficients(
                values
            )

            print(
                f"  X^{power}:"
            )

            print(
                f"    Delta-binomial coefficients="
                f"{delta_coeffs}"
            )

        print()


# ============================================================================
# 18. SECTION 7 -- TEST LOW-DEGREE k LAWS
# ============================================================================

def section_7_low_degree_k_tests() -> None:
    print("=" * 78)
    print("7. LOW-DEGREE k-LAW TESTS")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            len(x_coefficients(j, k))
            for k in KS
        ) - 1

        for power in range(max_degree + 1):

            values = []

            for k in KS:
                coeffs = x_coefficients(j, k)

                values.append(
                    coeffs[power]
                    if power < len(coeffs)
                    else sp.Integer(0)
                )

            deg = finite_difference_degree(
                values
            )

            poly = interpolate_k(
                list(zip(KS, values))
            )

            if deg <= 2:
                status = "LOW-DEGREE"
            else:
                status = "HIGH-DEGREE"

            print(
                f"  X^{power}: "
                f"degree={deg} "
                f"{status}"
            )

            print(
                f"    k-polynomial={poly}"
            )

        print()


# ============================================================================
# 19. SECTION 8 -- CONTROLLED NORMALIZATION SEARCH
# ============================================================================

def normalization_candidates(j: int) -> list[tuple[str, sp.Expr]]:
    """
    Controlled family.

    The purpose is not to force a fit. We only accept a normalization
    when it materially reduces the finite-difference degree for an entire
    coefficient row.
    """

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

    # Additional j-aware candidates.
    candidates.extend([
        (
            f"C(K+4,{5-j})",
            sp.binomial(K + 4, 5 - j)
        ),
        (
            f"C(K+5,{5-j})",
            sp.binomial(K + 5, 5 - j)
        ),
    ])

    # Deduplicate structurally.
    result = []
    seen = set()

    for name, expr in candidates:
        key = sp.srepr(
            sp.factor(expr)
        )

        if key not in seen:
            seen.add(key)
            result.append(
                (name, expr)
            )

    return result


def normalized_x_coefficient(
    j: int,
    power: int,
    k: int,
    scale: sp.Expr
) -> sp.Expr:

    coeffs = x_coefficients(j, k)

    if power >= len(coeffs):
        value = sp.Integer(0)
    else:
        value = coeffs[power]

    scale_value = sp.sympify(scale).subs(
        K,
        sp.Integer(k)
    )

    if scale_value == 0:
        raise ZeroDivisionError

    return sp.factor(
        sp.cancel(
            value / scale_value
        )
    )


def section_8_normalization_search() -> None:
    print("=" * 78)
    print("8. CONTROLLED COMBINATORIAL NORMALIZATION SEARCH")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            len(x_coefficients(j, k))
            for k in KS
        ) - 1

        candidates = normalization_candidates(j)

        for power in range(max_degree + 1):

            values_raw = []

            for k in KS:
                coeffs = x_coefficients(j, k)
                values_raw.append(
                    coeffs[power]
                    if power < len(coeffs)
                    else sp.Integer(0)
                )

            raw_degree = finite_difference_degree(
                values_raw
            )

            improvements = []

            for name, scale in candidates:

                try:
                    values = [
                        normalized_x_coefficient(
                            j,
                            power,
                            k,
                            scale
                        )
                        for k in KS
                    ]
                except ZeroDivisionError:
                    continue

                degree = finite_difference_degree(
                    values
                )

                if degree < raw_degree:
                    poly = interpolate_k(
                        list(
                            zip(KS, values)
                        )
                    )

                    improvements.append(
                        (
                            degree,
                            name,
                            sp.factor(poly)
                        )
                    )

            improvements.sort(
                key=lambda item: (
                    item[0],
                    item[1]
                )
            )

            print(
                f"  X^{power}: "
                f"raw_degree={raw_degree}"
            )

            if improvements:
                for degree, name, poly in improvements[:5]:
                    print(
                        f"    improvement: "
                        f"degree={degree} "
                        f"scale={name}"
                    )
                    print(
                        f"      normalized="
                        f"{poly}"
                    )
            else:
                print(
                    "    no tested normalization "
                    "reduces degree"
                )

        print()


# ============================================================================
# 20. SECTION 9 -- SPECIAL j=5 FACTORIZATION
# ============================================================================

def section_9_j5() -> None:
    print("=" * 78)
    print("9. SPECIAL j=5 TERMINAL FACTORIZATION")
    print("=" * 78)

    j = 5

    support = support_factor(j)

    print(
        f"support = {sp.factor(support)}"
    )
    print()

    for k in KS:
        Q = quotient_poly(j, k)

        print(
            f"k={k}"
        )
        print(
            f"  Q(D)={sp.factor(Q)}"
        )

        print(
            f"  Q(X)={shifted_quotient(j, k)}"
        )

        coeffs = x_coefficients(j, k)

        print(
            f"  X coefficients={coeffs}"
        )

        print()


# ============================================================================
# 21. SECTION 10 -- TERMINAL FACTOR TEST
# ============================================================================

def section_10_terminal_factor_test() -> None:
    print("=" * 78)
    print("10. TERMINAL FACTOR TEST")
    print("=" * 78)

    j = 5

    known_support = support_factor(j)

    extra_roots = [
        14,
    ]

    print(
        f"base support={sp.factor(known_support)}"
    )
    print()

    for k in KS:
        Q = quotient_poly(j, k)

        print(
            f"k={k}"
        )

        for root in extra_roots:
            value = sp.factor(
                Q.subs(D, root)
            )

            mult = root_multiplicity(
                discrepancy_poly(j, k),
                root
            )

            print(
                f"  D={root}: "
                f"Q(D)={value} "
                f"E-multiplicity={mult}"
            )

        print()


# ============================================================================
# 22. SECTION 11 -- COMMON ROOT SEARCH BEYOND OBSERVED ZERO SET
# ============================================================================

def section_11_candidate_even_root_search() -> None:
    print("=" * 78)
    print("11. CANDIDATE EVEN-ROOT SEARCH")
    print("=" * 78)

    # Search integer D candidates beyond the observed grid.
    candidates = list(range(0, 31, 2))

    for j in JS:
        print(f"j={j}")

        for d in candidates:

            if d in DS:
                continue

            values = [
                discrepancy_poly(j, k).subs(
                    D,
                    d
                )
                for k in KS
            ]

            if all(v == 0 for v in values):
                print(
                    f"  COMMON EXACT ROOT: D={d}"
                )

        print()


# ============================================================================
# 23. SECTION 12 -- COEFFICIENT CROSS-K MATRIX
# ============================================================================

def section_12_coefficient_matrix() -> None:
    print("=" * 78)
    print("12. COEFFICIENT CROSS-k MATRIX")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            len(x_coefficients(j, k))
            for k in KS
        ) - 1

        for power in range(max_degree + 1):

            row = []

            for k in KS:
                coeffs = x_coefficients(j, k)

                value = (
                    coeffs[power]
                    if power < len(coeffs)
                    else sp.Integer(0)
                )

                row.append(value)

            print(
                f"  X^{power}: "
                f"{row}"
            )

        print()


# ============================================================================
# 24. SECTION 13 -- BINOMIAL-BASIS DEGREE SUMMARY
# ============================================================================

def section_13_binomial_degree_summary() -> None:
    print("=" * 78)
    print("13. BINOMIAL-BASIS DEGREE SUMMARY")
    print("=" * 78)

    for j in JS:
        print(f"j={j}")

        max_degree = max(
            len(x_coefficients(j, k))
            for k in KS
        ) - 1

        for power in range(max_degree + 1):

            values = []

            for k in KS:
                coeffs = x_coefficients(j, k)

                values.append(
                    coeffs[power]
                    if power < len(coeffs)
                    else sp.Integer(0)
                )

            coeffs = newton_binomial_coefficients(
                values
            )

            last_nonzero = -1

            for idx, value in enumerate(coeffs):
                if value != 0:
                    last_nonzero = idx

            print(
                f"  X^{power}: "
                f"binomial-degree={last_nonzero}"
            )

            print(
                f"    coefficients={coeffs}"
            )

        print()


# ============================================================================
# 25. SECTION 14 -- EXACT RECONSTRUCTION CHECK
# ============================================================================

def section_14_reconstruction() -> None:
    print("=" * 78)
    print("14. EXACT SUPPORT-QUOTIENT RECONSTRUCTION")
    print("=" * 78)

    failures = 0
    tested = 0

    for j in JS:
        support = support_factor(j)

        for k in KS:
            E = discrepancy_poly(j, k)
            Q = quotient_poly(j, k)

            reconstructed = sp.expand(
                support * Q
            )

            residual = sp.factor(
                E - reconstructed
            )

            tested += 1

            if residual != 0:
                failures += 1

                print(
                    f"FAIL j={j} k={k}: "
                    f"residual={residual}"
                )

    print(
        f"tested={tested}"
    )
    print(
        f"support reconstruction failures={failures}"
    )
    print()


# ============================================================================
# 26. SECTION 15 -- INTERPOLATION CONSISTENCY
# ============================================================================

def section_15_interpolation_consistency() -> None:
    print("=" * 78)
    print("15. INTERPOLATION CONSISTENCY")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in JS:
        for k in KS:

            poly = discrepancy_poly(j, k)

            for d in sorted(DATA[k][j]):

                actual = discrepancy(
                    j,
                    k,
                    d
                )

                predicted = sp.factor(
                    poly.subs(D, d)
                )

                tested += 1

                if sp.factor(
                    predicted - actual
                ) != 0:

                    failures += 1

                    print(
                        f"FAIL "
                        f"j={j} k={k} d={d} "
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


# ============================================================================
# 27. SECTION 16 -- SUPPORT PATTERN SUMMARY
# ============================================================================

def section_16_summary() -> None:
    print("=" * 78)
    print("16. SUPPORT PATTERN SUMMARY")
    print("=" * 78)

    print(
        "Observed common support:"
    )

    for j in JS:
        roots = common_roots(j)

        print(
            f"  j={j}: "
            f"roots={roots} "
            f"support={sp.factor(support_factor(j))}"
        )

    print()

    print(
        "Observed support degrees:"
    )

    for j in JS:
        print(
            f"  j={j}: "
            f"{sp.Poly(support_factor(j), D).degree()}"
        )

    print()

    print(
        "Observed quotient degrees:"
    )

    for j in JS:

        row = []

        for k in KS:
            row.append(
                sp.Poly(
                    quotient_poly(j, k),
                    D
                ).degree()
            )

        print(
            f"  j={j}: {row}"
        )

    print()


# ============================================================================
# 28. SECTION 17 -- FINAL DIAGNOSTIC
# ============================================================================

def section_17_final_diagnostic() -> None:
    print("=" * 78)
    print("17. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "This experiment does NOT infer a new universal r,j law."
    )
    print()

    print(
        "It treats the observed r=5 discrepancy as:"
    )
    print()
    print(
        "    E_j(k,D) = Support_j(D) * Q_j(k,D)"
    )
    print()

    print(
        "The exact common-support structure found in the supplied data is:"
    )

    for j in JS:
        print(
            f"  j={j}: "
            f"{sp.factor(support_factor(j))}"
        )

    print()

    print(
        "Every common support root is simple in every tested k."
    )
    print()

    print(
        "The decisive next question is no longer whether the old law "
        "has the right roots."
    )
    print()

    print(
        "The decisive question is whether Q_j(k,D) admits a structured "
        "k-dependence after the support is removed."
    )
    print()

    print(
        "The experiment therefore reports both ordinary polynomial "
        "degree in k and exact Newton/binomial-basis degree."
    )
    print()

    print(
        "A low finite-difference degree is meaningful."
    )
    print(
        "A generic degree-5 interpolation across six k-values is not."
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


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    section_1_support_structure()
    section_2_quotient_degrees()
    section_3_shifted_quotients()
    section_4_coefficient_arrays()
    section_5_k_difference_degree()
    section_6_newton_basis()
    section_7_low_degree_k_tests()
    section_8_normalization_search()
    section_9_j5()
    section_10_terminal_factor_test()
    section_11_candidate_even_root_search()
    section_12_coefficient_matrix()
    section_13_binomial_degree_summary()
    section_14_reconstruction()
    section_15_interpolation_consistency()
    section_16_summary()
    section_17_final_diagnostic()


if __name__ == "__main__":
    main()

