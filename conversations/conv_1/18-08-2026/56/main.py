#!/usr/bin/env python3

"""
R=5 EXPERIMENT 56
NORMALIZED Y-COEFFICIENT RECURRENCE SEARCH

Input:
    Q_j(K,D) = sum_m a_{j,m}(K) * (D-D0(j))^m

This experiment asks whether the higher Y-coefficients can be described
structurally relative to the boundary coefficient a_{j,0}(K).

Tests:

1. Reconstruct exact coefficient functions a_{j,m}(K).
2. Form exact ratios a_{j,m}/a_{j,0}.
3. Search for low-degree rational descriptions.
4. Form successive coefficient ratios a_{j,m+1}/a_{j,m}.
5. Test the previously observed K-ladders against every coefficient.
6. Compare equal-Y-power coefficients between different j.
7. Special audit of j=5.
8. Integer-root and Newton-degree profiles.
9. Exact reconstruction checks.
10. Print symbolic fresh-K targets.

Important:
    The supplied six K values are interpolation points.
    No generic degree-5 interpolation is treated as a law.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


# ============================================================================
# SYMBOLS / GRID
# ============================================================================

K = sp.Symbol("K")

K_VALUES = [3, 5, 7, 9, 11, 13]
D_VALUES = [6, 8, 10, 12, 14, 16]
J_VALUES = [0, 1, 2, 3, 4, 5]

D0 = {
    0: 6,
    1: 8,
    2: 8,
    3: 10,
    4: 10,
    5: 14,
}

Q_DEGREE = {
    0: 5,
    1: 4,
    2: 4,
    3: 3,
    4: 3,
    5: 1,
}


# ============================================================================
# Y-COEFFICIENT DATA
#
# Each row corresponds to K in [3,5,7,9,11,13].
# Each row contains coefficients of powers of Y = D-D0(j).
# ============================================================================

Y_ROWS = {

    0: [
        [
            sp.Rational(-462),
            sp.Rational(-157489, 20),
            sp.Rational(300973, 48),
            sp.Rational(-89761, 48),
            sp.Rational(34865, 192),
            sp.Rational(-7363, 960),
        ],
        [
            sp.Rational(-2268),
            sp.Rational(221041, 10),
            sp.Rational(-1278395, 48),
            sp.Rational(737443, 96),
            sp.Rational(-201385, 192),
            sp.Rational(76393, 1920),
        ],
        [
            sp.Rational(-7524),
            sp.Rational(-53233, 4),
            sp.Rational(17461, 48),
            sp.Rational(-406949, 96),
            sp.Rational(82907, 192),
            sp.Rational(-15541, 384),
        ],
        [
            sp.Rational(-19734),
            sp.Rational(-2666789, 40),
            sp.Rational(1062167, 32),
            sp.Rational(-4147205, 192),
            sp.Rational(318293, 128),
            sp.Rational(-615341, 3840),
        ],
        [
            sp.Rational(-44226),
            sp.Rational(-455572),
            sp.Rational(776451, 2),
            sp.Rational(-4942847, 32),
            sp.Rational(158457, 8),
            sp.Rational(-125801, 128),
        ],
        [
            sp.Rational(-88536),
            sp.Rational(-51962659, 30),
            sp.Rational(38256001, 24),
            sp.Rational(-27464129, 48),
            sp.Rational(6918779, 96),
            sp.Rational(-1016379, 320),
        ],
    ],

    1: [
        [
            sp.Rational(-910),
            sp.Rational(-12477, 8),
            sp.Rational(4789, 48),
            sp.Rational(-695, 32),
            sp.Rational(-149, 96),
        ],
        [
            sp.Rational(-2856),
            sp.Rational(2135, 2),
            sp.Rational(-31549, 12),
            sp.Rational(9149, 24),
            sp.Rational(-1267, 48),
        ],
        [
            sp.Rational(-6930),
            sp.Rational(-284541, 40),
            sp.Rational(5321, 5),
            sp.Rational(-10667, 20),
            sp.Rational(1679, 80),
        ],
        [
            sp.Rational(-14300),
            sp.Rational(110737, 20),
            sp.Rational(-464629, 30),
            sp.Rational(132231, 40),
            sp.Rational(-125081, 480),
        ],
        [
            sp.Rational(-26390),
            sp.Rational(1004367, 40),
            sp.Rational(-2426879, 60),
            sp.Rational(2027753, 240),
            sp.Rational(-260533, 480),
        ],
        [
            sp.Rational(-44880),
            sp.Rational(187969, 5),
            sp.Rational(-878611, 15),
            sp.Rational(132923, 15),
            sp.Rational(-3009, 40),
        ],
    ],

    2: [
        [
            sp.Rational(-364),
            sp.Rational(-289867, 120),
            sp.Rational(284213, 160),
            sp.Rational(-384079, 960),
            sp.Rational(7569, 320),
        ],
        [
            sp.Rational(-816),
            sp.Rational(7695, 4),
            sp.Rational(-7406, 3),
            sp.Rational(29213, 48),
            sp.Rational(-289, 6),
        ],
        [
            sp.Rational(-1540),
            sp.Rational(2111),
            sp.Rational(-71963, 24),
            sp.Rational(115403, 192),
            sp.Rational(-5223, 128),
        ],
        [
            sp.Rational(-2600),
            sp.Rational(12995, 3),
            sp.Rational(-531133, 96),
            sp.Rational(99215, 96),
            sp.Rational(-21077, 384),
        ],
        [
            sp.Rational(-4060),
            sp.Rational(-10409, 24),
            sp.Rational(-26999, 24),
            sp.Rational(-83783, 96),
            sp.Rational(34685, 192),
        ],
        [
            sp.Rational(-5984),
            sp.Rational(-3330691, 60),
            sp.Rational(12171677, 240),
            sp.Rational(-1987453, 120),
            sp.Rational(1605157, 960),
        ],
    ],

    3: [
        [
            sp.Rational(-525, 8),
            sp.Rational(-5117, 10),
            sp.Rational(707, 5),
            sp.Rational(-917, 80),
        ],
        [
            sp.Rational(-855, 8),
            sp.Rational(11313, 32),
            sp.Rational(-557, 2),
            sp.Rational(4051, 128),
        ],
        [
            sp.Rational(-1265, 8),
            sp.Rational(-36267, 80),
            sp.Rational(54307, 240),
            sp.Rational(-9031, 320),
        ],
        [
            sp.Rational(-1755, 8),
            sp.Rational(156),
            sp.Rational(-7631, 24),
            sp.Rational(7735, 96),
        ],
        [
            sp.Rational(-2325, 8),
            sp.Rational(4733, 2),
            sp.Rational(-17017, 8),
            sp.Rational(1679, 4),
        ],
        [
            sp.Rational(-2975, 8),
            sp.Rational(4652033, 480),
            sp.Rational(-3780851, 480),
            sp.Rational(927163, 640),
        ],
    ],

    4: [
        [
            sp.Rational(-75, 4),
            sp.Rational(-9047, 32),
            sp.Rational(14401, 64),
            sp.Rational(-473, 16),
        ],
        [
            sp.Rational(-95, 4),
            sp.Rational(-67, 192),
            sp.Rational(-1441, 48),
            sp.Rational(8315, 768),
        ],
        [
            sp.Rational(-115, 4),
            sp.Rational(69027, 320),
            sp.Rational(-175957, 960),
            sp.Rational(124993, 3840),
        ],
        [
            sp.Rational(-135, 4),
            sp.Rational(10107, 16),
            sp.Rational(-1009, 2),
            sp.Rational(5649, 64),
        ],
        [
            sp.Rational(-155, 4),
            sp.Rational(229327, 120),
            sp.Rational(-1418347, 960),
            sp.Rational(487679, 1920),
        ],
        [
            sp.Rational(-175, 4),
            sp.Rational(4518607, 960),
            sp.Rational(-1724927, 480),
            sp.Rational(783647, 1280),
        ],
    ],

    5: [
        [
            sp.Rational(-2975, 48),
            sp.Rational(22277, 960),
        ],
        [
            sp.Rational(0),
            sp.Rational(-1295, 96),
        ],
        [
            sp.Rational(0),
            sp.Rational(7007, 480),
        ],
        [
            sp.Rational(0),
            sp.Rational(1377, 40),
        ],
        [
            sp.Rational(0),
            sp.Rational(46189, 640),
        ],
        [
            sp.Rational(133063, 96),
            sp.Rational(-88641, 160),
        ],
    ],
}


# ============================================================================
# BASIC HELPERS
# ============================================================================

def poly(expr):
    return sp.Poly(sp.cancel(expr), K, domain=sp.QQ)


def degree(expr):
    p = poly(expr)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def factor(expr):
    return sp.factor(sp.cancel(expr))


def interpolate(values):
    points = list(zip(K_VALUES, values))
    return sp.factor(sp.interpolate(points, K))


def exact_division(numer, denom):
    p = poly(numer)
    q = poly(denom)

    if p.is_zero:
        return sp.Integer(0), True

    quotient, remainder = sp.div(p, q)
    return sp.factor(quotient.as_expr()), remainder.is_zero


def rational_complexity(expr):
    num, den = sp.fraction(sp.cancel(expr))
    return degree(num), degree(den)


def newton_coefficients(values):
    levels = [list(map(sp.sympify, values))]

    while len(levels[-1]) > 1:
        prev = levels[-1]
        nxt = []

        for i in range(len(prev) - 1):
            nxt.append(
                sp.simplify(prev[i + 1] - prev[i])
            )

        levels.append(nxt)

    return [
        sp.simplify(level[0])
        for level in levels
    ]


def newton_degree(values):
    coeffs = newton_coefficients(values)
    result = -sp.oo

    for idx, coeff in enumerate(coeffs):
        if coeff != 0:
            result = idx

    return result


# ============================================================================
# EXACT RATIONAL FIT
# ============================================================================

def try_rational_fit(values, numerator_degree, denominator_degree):
    """
    Search for

        p(K) / q(K)

    with:
        deg p <= numerator_degree
        deg q <= denominator_degree

    and q monic.

    With six data points this is only exploratory.
    """

    npdeg = numerator_degree
    dqdeg = denominator_degree

    # Unknown numerator coefficients.
    a = sp.symbols(
        "a0:" + str(npdeg + 1)
    )

    # Unknown denominator coefficients except the leading monic term.
    b = sp.symbols(
        "b0:" + str(dqdeg)
    )

    numerator = sum(
        a[i] * K**i
        for i in range(npdeg + 1)
    )

    if dqdeg == 0:
        denominator = sp.Integer(1)
    else:
        denominator = (
            K**dqdeg
            + sum(
                b[i] * K**i
                for i in range(dqdeg)
            )
        )

    unknowns = list(a) + list(b)

    # Need exactly enough equations for a determined system.
    if len(unknowns) != len(K_VALUES):
        return None

    equations = []

    for k0, value in zip(K_VALUES, values):
        equation = sp.expand(
            numerator.subs(K, k0)
            - value * denominator.subs(K, k0)
        )
        equations.append(equation)

    solution = sp.solve(
        equations,
        unknowns,
        dict=True,
        simplify=False,
    )

    if not solution:
        return None

    sol = solution[0]

    result_num = sp.expand(
        numerator.subs(sol)
    )

    result_den = sp.expand(
        denominator.subs(sol)
    )

    if result_den == 0:
        return None

    result = sp.factor(
        sp.cancel(
            result_num / result_den
        )
    )

    # Exact verification.
    for k0, value in zip(K_VALUES, values):
        test = sp.simplify(
            result.subs(K, k0) - value
        )

        if test != 0:
            return None

    return result


# ============================================================================
# RECONSTRUCT a_{j,m}(K)
# ============================================================================

def coefficient_polynomials():
    result = {}

    for j in J_VALUES:
        result[j] = []

        for m in range(Q_DEGREE[j] + 1):
            values = [
                Y_ROWS[j][row][m]
                for row in range(6)
            ]

            result[j].append(
                interpolate(values)
            )

    return result


# ============================================================================
# SECTION 0
# ============================================================================

def section_0_validation():
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    total_points = 0
    coefficient_rows = 0

    for j in J_VALUES:
        rows = Y_ROWS[j]

        if len(rows) != len(K_VALUES):
            raise RuntimeError(
                f"bad row count for j={j}"
            )

        expected_width = Q_DEGREE[j] + 1

        for row in rows:
            if len(row) != expected_width:
                raise RuntimeError(
                    f"bad coefficient width for j={j}: "
                    f"{len(row)} != {expected_width}"
                )

            coefficient_rows += 1
            total_points += len(D_VALUES)

    if total_points != 216:
        raise RuntimeError(
            f"bad point count: {total_points}"
        )

    print(f"k values = {K_VALUES}")
    print(f"j values = {J_VALUES}")
    print(f"D values = {D_VALUES}")
    print(f"points   = {total_points}")
    print("grid status = OK")
    print(f"checked coefficient rows = {coefficient_rows}")
    print()


# ============================================================================
# SECTION 1
# ============================================================================

def section_1_boundary_normalization(polys):
    print("=" * 78)
    print("1. BOUNDARY-NORMALIZED COEFFICIENTS")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        a0 = polys[j][0]

        print(
            f"  a_{{{j},0}}(K)={factor(a0)}"
        )

        if a0 == 0:
            print("  boundary coefficient is ZERO")
            print()
            continue

        for m in range(1, Q_DEGREE[j] + 1):
            am = polys[j][m]

            ratio = sp.factor(
                sp.cancel(am / a0)
            )

            nd, dd = rational_complexity(ratio)

            print(
                f"  m={m}: "
                f"a_{{{j},{m}}}/a_{{{j},0}}="
                f"{ratio}"
            )
            print(
                f"    numerator_degree={nd}, "
                f"denominator_degree={dd}"
            )

        print()


# ============================================================================
# SECTION 2
# ============================================================================

def section_2_low_complexity_rational_search(polys):
    print("=" * 78)
    print("2. LOW-COMPLEXITY RATIONAL SEARCH")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for m in range(1, Q_DEGREE[j] + 1):

            raw_values = [
                Y_ROWS[j][i][m]
                for i in range(6)
            ]

            boundary_values = [
                Y_ROWS[j][i][0]
                for i in range(6)
            ]

            if any(
                value == 0
                for value in boundary_values
            ):
                print(
                    f"  m={m}: skipped because "
                    f"a_{{{j},0}} vanishes on supplied K-grid"
                )
                continue

            ratio_values = [
                sp.factor(
                    sp.cancel(
                        raw_values[i]
                        / boundary_values[i]
                    )
                )
                for i in range(6)
            ]

            hits = []

            # Only search determined low-complexity cases.
            for npdeg in range(0, 4):
                for dqdeg in range(0, 4):

                    number_unknowns = (
                        npdeg + 1 + dqdeg
                    )

                    if number_unknowns != 6:
                        continue

                    result = try_rational_fit(
                        ratio_values,
                        npdeg,
                        dqdeg,
                    )

                    if result is not None:
                        hits.append(
                            (
                                npdeg,
                                dqdeg,
                                result,
                            )
                        )

            if hits:
                hits.sort(
                    key=lambda item: (
                        item[0] + item[1],
                        item[0],
                        item[1],
                    )
                )

                npdeg, dqdeg, result = hits[0]

                print(
                    f"  m={m}: "
                    f"FOUND low-complexity fit"
                )
                print(
                    f"    deg_num={npdeg}"
                )
                print(
                    f"    deg_den={dqdeg}"
                )
                print(
                    f"    ratio={result}"
                )
            else:
                print(
                    f"  m={m}: "
                    f"no determined low-degree rational fit"
                )

        print()


# ============================================================================
# SECTION 3
# ============================================================================

def section_3_successive_m_ratios(polys):
    print("=" * 78)
    print("3. SUCCESSIVE-Y-COEFFICIENT RATIOS")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for m in range(Q_DEGREE[j]):

            current = polys[j][m]
            following = polys[j][m + 1]

            if current == 0:
                print(
                    f"  m={m}->{m+1}: "
                    "current coefficient is ZERO"
                )
                continue

            ratio = sp.factor(
                sp.cancel(
                    following / current
                )
            )

            nd, dd = rational_complexity(ratio)

            print(
                f"  a_{{{j},{m+1}}}/a_{{{j},{m}}}"
                f"={ratio}"
            )
            print(
                f"    numerator_degree={nd}, "
                f"denominator_degree={dd}"
            )

        print()


# ============================================================================
# SECTION 4
# ============================================================================

def section_4_ladder_divisibility(polys):
    print("=" * 78)
    print("4. K-LADDER DIVISIBILITY ACROSS ALL Y-COEFFICIENTS")
    print("=" * 78)

    ladders = {
        0: (
            (K + 1)
            * (K + 2)
            * (K + 3)
            * (K + 4)
        ),
        1: (
            (K + 2)
            * (K + 3)
            * (K + 4)
        ),
        2: (
            (K + 3)
            * (K + 4)
        ),
        3: (
            K + 4
        ),
        4: sp.Integer(1),
        5: (
            (K - 5)
            * (K - 7)
            * (K - 9)
            * (K - 11)
        ),
    }

    for j in J_VALUES:
        ladder = ladders[j]

        print(f"j={j}")
        print(
            f"  ladder={factor(ladder)}"
        )

        for m, am in enumerate(polys[j]):

            if ladder == 1:
                print(
                    f"  a_{{{j},{m}}}: "
                    "trivially divisible"
                )
                continue

            quotient, exact = exact_division(
                am,
                ladder,
            )

            if exact:
                print(
                    f"  a_{{{j},{m}}}: YES"
                )
                print(
                    f"    quotient={factor(quotient)}"
                )
            else:
                gcd_expr = factor(
                    sp.gcd(
                        poly(am),
                        poly(ladder),
                    ).as_expr()
                )

                print(
                    f"  a_{{{j},{m}}}: NO"
                )
                print(
                    f"    gcd_with_ladder={gcd_expr}"
                )

        print()


# ============================================================================
# SECTION 5
# ============================================================================

def section_5_cross_j_same_m(polys):
    print("=" * 78)
    print("5. CROSS-j GCD TEST AT FIXED Y-POWER")
    print("=" * 78)

    max_m = max(
        Q_DEGREE[j]
        for j in J_VALUES
    )

    for m in range(max_m + 1):
        print(f"Y-power m={m}")

        active = [
            j
            for j in J_VALUES
            if m <= Q_DEGREE[j]
        ]

        found = False

        for j1, j2 in combinations(active, 2):

            p1 = poly(polys[j1][m])
            p2 = poly(polys[j2][m])

            if p1.is_zero or p2.is_zero:
                continue

            gcd_expr = factor(
                sp.gcd(
                    p1,
                    p2,
                ).as_expr()
            )

            if gcd_expr != 1:
                found = True

                print(
                    f"  gcd("
                    f"a_{{{j1},{m}}}, "
                    f"a_{{{j2},{m}}}"
                    f")={gcd_expr}"
                )

        if not found:
            print(
                "  no nontrivial cross-j gcds"
            )

        print()


# ============================================================================
# SECTION 6
# ============================================================================

def section_6_cross_j_normalized(polys):
    print("=" * 78)
    print("6. CROSS-j NORMALIZED COEFFICIENT RATIOS")
    print("=" * 78)

    for j in range(5):
        print(
            f"j={j} -> j={j+1}"
        )

        common_m = min(
            Q_DEGREE[j],
            Q_DEGREE[j + 1],
        )

        for m in range(common_m + 1):

            left = polys[j + 1][m]
            right = polys[j][m]

            if right == 0:
                print(
                    f"  m={m}: denominator ZERO"
                )
                continue

            ratio = sp.factor(
                sp.cancel(
                    left / right
                )
            )

            nd, dd = rational_complexity(ratio)

            print(
                f"  m={m}: "
                f"a_{{{j+1},{m}}}/a_{{{j},{m}}}"
                f"={ratio}"
            )
            print(
                f"    numerator_degree={nd}, "
                f"denominator_degree={dd}"
            )

        print()


# ============================================================================
# SECTION 7
# ============================================================================

def section_7_j5_terminal_audit(polys):
    print("=" * 78)
    print("7. SPECIAL j=5 TERMINAL AUDIT")
    print("=" * 78)

    a0 = polys[5][0]
    a1 = polys[5][1]

    ladder = (
        (K - 5)
        * (K - 7)
        * (K - 9)
        * (K - 11)
    )

    print(
        f"a_{{5,0}}={factor(a0)}"
    )
    print(
        f"a_{{5,1}}={factor(a1)}"
    )
    print(
        f"terminal_ladder={factor(ladder)}"
    )

    q0, ok0 = exact_division(
        a0,
        ladder,
    )

    q1, ok1 = exact_division(
        a1,
        ladder,
    )

    print(
        f"a_{{5,0}} divisible={ok0}"
    )

    if ok0:
        print(
            f"a_{{5,0}}/ladder={factor(q0)}"
        )

    print(
        f"a_{{5,1}} divisible={ok1}"
    )

    print()
    print("terminal K-values:")

    for k0 in [5, 7, 9, 11]:
        print(
            f"  K={k0}: "
            f"a_{{5,0}}={factor(a0.subs(K, k0))}, "
            f"a_{{5,1}}={factor(a1.subs(K, k0))}"
        )

    print()


# ============================================================================
# SECTION 8
# ============================================================================

def section_8_integer_root_profile(polys):
    print("=" * 78)
    print("8. INTEGER ROOT PROFILE")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for m, am in enumerate(polys[j]):

            if am == 0:
                print(
                    f"  a_{{{j},{m}}}: ZERO"
                )
                continue

            roots = []

            p = poly(am)

            for root in sp.ground_roots(p).keys():

                if root.is_Integer:
                    roots.append(
                        int(root)
                    )

            print(
                f"  a_{{{j},{m}}}: "
                f"integer_roots={roots}"
            )

        print()


# ============================================================================
# SECTION 9
# ============================================================================

def section_9_newton_profile():
    print("=" * 78)
    print("9. NEWTON / FINITE-DIFFERENCE DEGREE PROFILE")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for m in range(Q_DEGREE[j] + 1):

            values = [
                Y_ROWS[j][i][m]
                for i in range(6)
            ]

            coeffs = newton_coefficients(
                values
            )

            nd = -sp.oo

            for r, c in enumerate(coeffs):
                if c != 0:
                    nd = r

            print(
                f"  a_{{{j},{m}}}: "
                f"newton_degree={nd}"
            )

        print()


# ============================================================================
# SECTION 10
# ============================================================================

def section_10_exact_reconstruction(polys):
    print("=" * 78)
    print("10. EXACT COEFFICIENT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = []

    for j in J_VALUES:

        for row_index, k0 in enumerate(K_VALUES):

            for m in range(Q_DEGREE[j] + 1):

                predicted = sp.simplify(
                    polys[j][m].subs(K, k0)
                )

                actual = sp.simplify(
                    Y_ROWS[j][row_index][m]
                )

                tested += 1

                if predicted != actual:
                    failures.append(
                        (
                            j,
                            k0,
                            m,
                            actual,
                            predicted,
                        )
                    )

    print(
        f"tested={tested}"
    )
    print(
        f"reconstruction failures={len(failures)}"
    )

    for failure in failures[:20]:
        print(
            f"  failure={failure}"
        )

    print()


# ============================================================================
# SECTION 11
# ============================================================================

def section_11_fresh_k_targets(polys):
    print("=" * 78)
    print("11. FRESH-K SYMBOLIC TARGETS")
    print("=" * 78)

    fresh_k = [
        -9,
        -7,
        -5,
        -3,
        -1,
        1,
        15,
        17,
        19,
    ]

    for j in J_VALUES:
        print(f"j={j}")

        for m in range(
            1,
            Q_DEGREE[j] + 1,
        ):

            a0 = polys[j][0]
            am = polys[j][m]

            if a0 == 0:
                print(
                    f"  m={m}: boundary coefficient ZERO"
                )
                continue

            ratio = sp.cancel(
                am / a0
            )

            print(
                f"  m={m}: "
                f"R_{{{j},{m}}}(K)=a_{{{j},{m}}}/a_{{{j},0}}"
            )

            for k0 in fresh_k:

                denominator = sp.simplify(
                    a0.subs(K, k0)
                )

                if denominator == 0:
                    print(
                        f"    K={k0}: undefined"
                    )
                    continue

                value = sp.factor(
                    sp.cancel(
                        ratio.subs(K, k0)
                    )
                )

                print(
                    f"    K={k0}: {value}"
                )

        print()


# ============================================================================
# SECTION 12
# ============================================================================

def section_12_compact_summary(polys):
    print("=" * 78)
    print("12. COMPACT STRUCTURAL SUMMARY")
    print("=" * 78)

    header = (
        "j | m | deg(a) | deg(num ratio) | "
        "deg(den ratio) | a_{j,m}"
    )

    print(header)
    print("-" * 125)

    for j in J_VALUES:

        for m, am in enumerate(polys[j]):

            if am == 0:
                print(
                    f"{j} | {m} | ZERO"
                )
                continue

            deg_a = degree(am)

            if m == 0:
                nd = "-"
                dd = "-"
            else:
                ratio = sp.cancel(
                    am / polys[j][0]
                )

                nd, dd = rational_complexity(
                    ratio
                )

            print(
                f"{j} | {m} | "
                f"{deg_a} | "
                f"{nd} | "
                f"{dd} | "
                f"{factor(am)}"
            )

    print()


# ============================================================================
# FINAL DIAGNOSTIC
# ============================================================================

def final_diagnostic():
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The previous rank experiment ruled out exact rank-one "
        "k-D separability."
    )
    print()

    print(
        "The present experiment decomposes the support-removed quotient as"
    )
    print()
    print(
        "    Q_j(K,D) = sum_m a_{j,m}(K) (D-D0(j))^m"
    )
    print()

    print(
        "The key question is whether the higher coefficients are "
        "structurally generated from a_{j,0}(K)."
    )
    print()

    print(
        "A low-degree exact ratio"
    )
    print(
        "    a_{j,m}(K) / a_{j,0}(K)"
    )
    print(
        "would be meaningful."
    )
    print()

    print(
        "A generic degree-5 / degree-5 rational expression obtained "
        "from six interpolation points is NOT treated as evidence."
    )
    print()

    print(
        "The same restriction applies to successive ratios"
    )
    print(
        "    a_{j,m+1}(K) / a_{j,m}(K)."
    )
    print()

    print(
        "The clean boundary ladders remain:"
    )
    print(
        "  j=0 : (K+1)(K+2)(K+3)(K+4)"
    )
    print(
        "  j=1 : (K+2)(K+3)(K+4)"
    )
    print(
        "  j=2 : (K+3)(K+4)"
    )
    print(
        "  j=3 : (K+4)"
    )
    print(
        "  j=4 : 1"
    )
    print(
        "  j=5 : (K-5)(K-7)(K-9)(K-11)"
    )
    print()

    print(
        "Fresh K values are needed to distinguish genuine rational "
        "structure from finite interpolation."
    )
    print()

    print("No r=6.")
    print("No full pq-kernel expansion.")
    print("No replacement universal r,j formula.")
    print("=" * 78)


# ============================================================================
# MAIN
# ============================================================================

def main():
    section_0_validation()

    polys = coefficient_polynomials()

    section_1_boundary_normalization(polys)
    section_2_low_complexity_rational_search(polys)
    section_3_successive_m_ratios(polys)
    section_4_ladder_divisibility(polys)
    section_5_cross_j_same_m(polys)
    section_6_cross_j_normalized(polys)
    section_7_j5_terminal_audit(polys)
    section_8_integer_root_profile(polys)
    section_9_newton_profile()
    section_10_exact_reconstruction(polys)
    section_11_fresh_k_targets(polys)
    section_12_compact_summary(polys)

    final_diagnostic()


if __name__ == "__main__":
    main()