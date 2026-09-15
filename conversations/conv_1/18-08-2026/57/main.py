#!/usr/bin/env python3

"""
R=5 EXPERIMENT 57
CROSS-j COEFFICIENT-SPACE / SHIFTED-K STRUCTURE

We have

    Q_j(K,D) = sum_m a_{j,m}(K) (D-D0(j))^m

This experiment tests whether the coefficient functions share
a low-dimensional polynomial structure across j, including after
j-dependent shifts of K and after removing the leading boundary
component a_{j,0}.

No universal law is inferred.
"""

from __future__ import annotations

import sympy as sp


# ============================================================================
# SYMBOLS / GRID
# ============================================================================

K = sp.Symbol("K")

K_VALUES = [3, 5, 7, 9, 11, 13]
J_VALUES = [0, 1, 2, 3, 4, 5]
D_VALUES = [6, 8, 10, 12, 14, 16]

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
# Rows correspond to K_VALUES in order.
# Each row contains coefficients of
#
#     Q_j(K,D) = sum_m a_{j,m}(K) * (D-D0(j))^m
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

def polynomial(expr):
    return sp.Poly(
        sp.expand(expr),
        K,
        domain=sp.QQ,
    )


def polynomial_degree(expr):
    p = polynomial(expr)

    if p.is_zero:
        return -sp.oo

    return p.degree()


def coefficient_vector(expr, width=6):
    p = polynomial(expr)

    return [
        p.as_expr().coeff(K, i)
        for i in range(width)
    ]


def interpolate_from_grid(values):
    points = [
        (K_VALUES[i], values[i])
        for i in range(len(K_VALUES))
    ]

    return sp.factor(
        sp.interpolate(
            points,
            K,
        )
    )


def reconstruct_coefficients():
    coefficients = {}

    for j in J_VALUES:
        coefficients[j] = []

        for m in range(Q_DEGREE[j] + 1):

            values = [
                Y_ROWS[j][row][m]
                for row in range(6)
            ]

            coefficients[j].append(
                interpolate_from_grid(values)
            )

    return coefficients


def polynomial_matrix(expressions):
    if not expressions:
        return sp.Matrix([])

    return sp.Matrix(
        [
            coefficient_vector(expr)
            for expr in expressions
        ]
    )


def polynomial_rank(expressions):
    if not expressions:
        return 0

    return polynomial_matrix(
        expressions
    ).rank()


def shift_k(expr, shift):
    return sp.expand(
        expr.subs(
            K,
            K + shift,
        )
    )


def exact_leading_projection(f, g):
    pf = polynomial(f)
    pg = polynomial(g)

    if pf.is_zero or pg.is_zero:
        return None

    return sp.cancel(
        pf.LC() / pg.LC()
    )


def subtract_leading_projection(f, g):
    c = exact_leading_projection(f, g)

    if c is None:
        return None, None

    residual = sp.factor(
        sp.expand(
            f - c * g
        )
    )

    return c, residual


def primitive_integer_form(expr):
    """
    Return a primitive integer coefficient vector in descending degree.
    Handles zero and single-coefficient polynomials safely.
    """
    p = polynomial(expr)

    if p.is_zero:
        return [0]

    coeffs = p.all_coeffs()

    denominators = [
        int(sp.denom(c))
        for c in coeffs
    ]

    if len(denominators) == 1:
        den_lcm = denominators[0]
    else:
        den_lcm = int(
            sp.ilcm(
                *denominators
            )
        )

    ints = [
        int(c * den_lcm)
        for c in coeffs
    ]

    nonzero = [
        abs(v)
        for v in ints
        if v != 0
    ]

    if not nonzero:
        return [0]

    if len(nonzero) == 1:
        content = nonzero[0]
    else:
        content = int(
            sp.igcd(
                *nonzero
            )
        )

    ints = [
        v // content
        for v in ints
    ]

    for v in ints:
        if v != 0:
            if v < 0:
                ints = [-x for x in ints]
            break

    return ints


# ============================================================================
# SECTION 0
# ============================================================================

def section_0():
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    coefficient_rows = 0

    for j in J_VALUES:

        expected_width = Q_DEGREE[j] + 1

        if j not in Y_ROWS:
            raise RuntimeError(
                f"missing Y_ROWS for j={j}"
            )

        if len(Y_ROWS[j]) != len(K_VALUES):
            raise RuntimeError(
                f"bad K-row count for j={j}: "
                f"{len(Y_ROWS[j])}"
            )

        for row in Y_ROWS[j]:

            if len(row) != expected_width:
                raise RuntimeError(
                    f"bad coefficient width for "
                    f"j={j}: "
                    f"{len(row)} != {expected_width}"
                )

            coefficient_rows += 1

    points = (
        len(K_VALUES)
        * len(J_VALUES)
        * len(D_VALUES)
    )

    print(
        f"k values = {K_VALUES}"
    )
    print(
        f"j values = {J_VALUES}"
    )
    print(
        f"D values = {D_VALUES}"
    )
    print(
        f"points   = {points}"
    )
    print(
        "grid status = OK"
    )
    print(
        f"checked coefficient rows = "
        f"{coefficient_rows}"
    )
    print()


# ============================================================================
# SECTION 1
# ============================================================================

def section_1_cross_j_rank(coeffs):
    print("=" * 78)
    print("1. CROSS-j POLYNOMIAL-SPACE RANK BY Y-COEFFICIENT")
    print("=" * 78)

    for m in range(6):

        active = [
            j
            for j in J_VALUES
            if m <= Q_DEGREE[j]
        ]

        family = [
            coeffs[j][m]
            for j in active
        ]

        rank = polynomial_rank(family)

        print(
            f"m={m}: "
            f"active_j={active} "
            f"rank={rank}/{len(active)}"
        )

        if rank == 1:
            print(
                "  RESULT: scalar-multiple family"
            )
        elif rank <= 2:
            print(
                "  RESULT: very small polynomial subspace"
            )
        elif rank < len(active):
            print(
                "  RESULT: nontrivial shared subspace"
            )
        else:
            print(
                "  RESULT: full cross-j rank"
            )

    print()


# ============================================================================
# SECTION 2
# ============================================================================

def section_2_shifted_rank(coeffs):
    print("=" * 78)
    print("2. SHIFTED-K CROSS-j RANK TEST")
    print("=" * 78)

    shifts = {
        "K+j": lambda j: j,
        "K-j": lambda j: -j,
        "K+(4-j)": lambda j: 4 - j,
        "K-(4-j)": lambda j: -(4 - j),
    }

    for label, shift_fn in shifts.items():

        print(
            f"SHIFT = {label}"
        )

        for m in range(6):

            active = [
                j
                for j in J_VALUES
                if m <= Q_DEGREE[j]
            ]

            family = [
                shift_k(
                    coeffs[j][m],
                    shift_fn(j),
                )
                for j in active
            ]

            rank = polynomial_rank(family)

            print(
                f"  m={m}: "
                f"rank={rank}/{len(active)}"
            )

        print()


# ============================================================================
# SECTION 3
# ============================================================================

def section_3_boundary_projection(coeffs):
    print("=" * 78)
    print("3. LEADING-TERM PROJECTION AGAINST a_{j,0}")
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        a0 = coeffs[j][0]

        for m in range(
            1,
            Q_DEGREE[j] + 1
        ):

            am = coeffs[j][m]

            c, residual = subtract_leading_projection(
                am,
                a0,
            )

            if c is None:
                print(
                    f"  m={m}: unavailable"
                )
                continue

            print(
                f"  m={m}: projection_c={sp.factor(c)}"
            )
            print(
                f"    original_degree="
                f"{polynomial_degree(am)}"
            )
            print(
                f"    residual_degree="
                f"{polynomial_degree(residual)}"
            )
            print(
                f"    residual="
                f"{sp.factor(residual)}"
            )

        print()


# ============================================================================
# SECTION 4
# ============================================================================

def section_4_residual_rank(coeffs):
    print("=" * 78)
    print("4. CROSS-j RANK AFTER LEADING-TERM REMOVAL")
    print("=" * 78)

    for m in range(1, 6):

        active = [
            j
            for j in J_VALUES
            if m <= Q_DEGREE[j]
        ]

        residuals = []

        for j in active:

            c, residual = subtract_leading_projection(
                coeffs[j][m],
                coeffs[j][0],
            )

            if residual is not None:
                residuals.append(
                    residual
                )

        rank = polynomial_rank(
            residuals
        )

        print(
            f"m={m}: residual rank="
            f"{rank}/{len(residuals)}"
        )

    print()


# ============================================================================
# SECTION 5
# ============================================================================

def section_5_primitive_forms(coeffs):
    print("=" * 78)
    print("5. PRIMITIVE INTEGER FORMS")
    print("=" * 78)

    for m in range(6):

        active = [
            j
            for j in J_VALUES
            if m <= Q_DEGREE[j]
        ]

        print(
            f"Y-power m={m}"
        )

        for j in active:

            primitive = primitive_integer_form(
                coeffs[j][m]
            )

            print(
                f"  j={j}: {primitive}"
            )

        print()


# ============================================================================
# SECTION 6
# ============================================================================

def section_6_cross_j_gcd(coeffs):
    print("=" * 78)
    print("6. CROSS-j GCD OF COEFFICIENT FUNCTIONS")
    print("=" * 78)

    for m in range(6):

        active = [
            j
            for j in J_VALUES
            if m <= Q_DEGREE[j]
        ]

        print(
            f"m={m}"
        )

        found = False

        for i in range(len(active)):

            for q in range(i + 1, len(active)):

                j1 = active[i]
                j2 = active[q]

                p1 = polynomial(
                    coeffs[j1][m]
                )

                p2 = polynomial(
                    coeffs[j2][m]
                )

                if p1.is_zero or p2.is_zero:
                    continue

                g = sp.factor(
                    sp.gcd(
                        p1,
                        p2,
                    ).as_expr()
                )

                if g != 1:
                    found = True

                    print(
                        f"  gcd("
                        f"a_{{{j1},{m}}}, "
                        f"a_{{{j2},{m}}}"
                        f")={g}"
                    )

        if not found:
            print(
                "  no nontrivial cross-j gcds"
            )

        print()


# ============================================================================
# SECTION 7
# ============================================================================

def section_7_shifted_root_audit(coeffs):
    print("=" * 78)
    print("7. INTEGER-ROOT AUDIT AFTER j-DEPENDENT K SHIFTS")
    print("=" * 78)

    shifts = {
        "K+j": lambda j: j,
        "K-j": lambda j: -j,
        "K+(4-j)": lambda j: 4 - j,
        "K-(4-j)": lambda j: -(4 - j),
    }

    for label, shift_fn in shifts.items():

        print(
            f"SHIFT = {label}"
        )

        for m in range(6):

            active = [
                j
                for j in J_VALUES
                if m <= Q_DEGREE[j]
            ]

            print(
                f"  m={m}"
            )

            for j in active:

                expr = shift_k(
                    coeffs[j][m],
                    shift_fn(j),
                )

                roots = []

                p = polynomial(expr)

                if not p.is_zero:

                    ground = sp.ground_roots(
                        p
                    )

                    for root in ground:

                        if root.is_Integer:
                            roots.append(
                                int(root)
                            )

                print(
                    f"    j={j}: "
                    f"integer_roots={roots}"
                )

        print()


# ============================================================================
# SECTION 8
# ============================================================================

def section_8_small_basis(coeffs):
    print("=" * 78)
    print("8. SMALL COMMON-BASIS SEARCH")
    print("=" * 78)

    for m in range(6):

        active = [
            j
            for j in J_VALUES
            if m <= Q_DEGREE[j]
        ]

        family = [
            coeffs[j][m]
            for j in active
        ]

        rank = polynomial_rank(
            family
        )

        print(
            f"m={m}: rank={rank}/{len(active)}"
        )

        if rank <= 3:
            print(
                "  SMALL BASIS CANDIDATE"
            )
        else:
            print(
                "  no basis of dimension <= 3"
            )

    print()


# ============================================================================
# SECTION 9
# ============================================================================

def section_9_exact_reconstruction(coeffs):
    print("=" * 78)
    print("9. EXACT COEFFICIENT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = []

    for j in J_VALUES:

        for row_index, k0 in enumerate(
            K_VALUES
        ):

            for m in range(
                Q_DEGREE[j] + 1
            ):

                predicted = sp.simplify(
                    coeffs[j][m].subs(
                        K,
                        k0,
                    )
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
        f"reconstruction failures="
        f"{len(failures)}"
    )

    for failure in failures[:20]:
        print(
            f"  {failure}"
        )

    print()


# ============================================================================
# SECTION 10
# ============================================================================

def section_10_fresh_targets(coeffs):
    print("=" * 78)
    print("10. FRESH-K SYMBOLIC TARGETS")
    print("=" * 78)

    fresh = [
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

        print(
            f"j={j}"
        )

        for m in range(
            Q_DEGREE[j] + 1
        ):

            expr = coeffs[j][m]

            print(
                f"  a_{{{j},{m}}}:"
            )

            for k0 in fresh:

                value = sp.factor(
                    expr.subs(
                        K,
                        k0,
                    )
                )

                print(
                    f"    K={k0}: {value}"
                )

        print()


# ============================================================================
# SECTION 11
# ============================================================================

def section_11_summary(coeffs):
    print("=" * 78)
    print("11. COMPACT RANK SUMMARY")
    print("=" * 78)

    print(
        "m | active_j | ordinary | K+j | K-j | K+(4-j) | K-(4-j)"
    )
    print("-" * 78)

    shifts = {
        "K+j": lambda j: j,
        "K-j": lambda j: -j,
        "K+(4-j)": lambda j: 4 - j,
        "K-(4-j)": lambda j: -(4 - j),
    }

    for m in range(6):

        active = [
            j
            for j in J_VALUES
            if m <= Q_DEGREE[j]
        ]

        ordinary = polynomial_rank(
            [
                coeffs[j][m]
                for j in active
            ]
        )

        shifted_ranks = []

        for shift_fn in shifts.values():

            rank = polynomial_rank(
                [
                    shift_k(
                        coeffs[j][m],
                        shift_fn(j),
                    )
                    for j in active
                ]
            )

            shifted_ranks.append(rank)

        print(
            f"{m} | "
            f"{active} | "
            f"{ordinary} | "
            f"{shifted_ranks[0]} | "
            f"{shifted_ranks[1]} | "
            f"{shifted_ranks[2]} | "
            f"{shifted_ranks[3]}"
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
        "The earlier experiments established:"
    )
    print(
        "  * exact support factors in D;"
    )
    print(
        "  * clean boundary factors in a_{j,0}(K);"
    )
    print(
        "  * loss of those factors for m>0;"
    )
    print(
        "  * full-rank behavior in the raw k-D coefficient matrix;"
    )
    print(
        "  * generic-looking degree-5 higher coefficient functions."
    )
    print()

    print(
        "This experiment tests whether that generic appearance is "
        "an artifact of the K origin."
    )
    print()

    print(
        "The decisive diagnostics are:"
    )
    print(
        "  1. cross-j polynomial rank by m;"
    )
    print(
        "  2. rank after K -> K+j and related shifts;"
    )
    print(
        "  3. residual rank after removing the leading "
        "K-degree contribution proportional to a_{j,0};"
    )
    print(
        "  4. nontrivial cross-j coefficient gcds."
    )
    print()

    print(
        "A genuine structural collapse would need to survive "
        "exact reconstruction and subsequently be checked on "
        "fresh K values."
    )
    print()

    print(
        "No universal R=5 law is inferred."
    )
    print(
        "No r=6."
    )
    print(
        "No full pq-kernel expansion."
    )
    print(
        "No replacement universal r,j formula."
    )

    print("=" * 78)


# ============================================================================
# MAIN
# ============================================================================

def main():

    section_0()

    coeffs = reconstruct_coefficients()

    section_1_cross_j_rank(coeffs)
    section_2_shifted_rank(coeffs)
    section_3_boundary_projection(coeffs)
    section_4_residual_rank(coeffs)
    section_5_primitive_forms(coeffs)
    section_6_cross_j_gcd(coeffs)
    section_7_shifted_root_audit(coeffs)
    section_8_small_basis(coeffs)
    section_9_exact_reconstruction(coeffs)
    section_10_fresh_targets(coeffs)
    section_11_summary(coeffs)

    final_diagnostic()


if __name__ == "__main__":
    main()