#!/usr/bin/env python3

"""
R=5 EXPERIMENT 59
============================================================

BOUNDARY-LADDER DEFECT DECOMPOSITION

Previous experiments established:

  * exact D-support factors;
  * clean boundary ladders in a_{j,0}(K);
  * higher a_{j,m}(K), m>0, generally lose those factors;
  * raw k-D separability fails;
  * simple K -> K +/- j shifts fail to reduce rank;
  * leading-term subtraction fails to reduce rank;
  * Newton/falling-factorial coordinates do not create sparsity
    in the higher coefficients.

This experiment studies the exact defect from the known
boundary ladder.

For each j define L_j(K), then divide every coefficient:

    a_{j,m}(K) = L_j(K) q_{j,m}(K) + r_{j,m}(K)

with

    deg(r_{j,m}) < deg(L_j).

We examine:

  1. quotient and remainder degrees;
  2. exact remainder factorization;
  3. remainder values at every ladder root;
  4. rank of the root-value matrices;
  5. rank of remainder coefficient spaces;
  6. rank of quotient coefficient spaces;
  7. whether m=0 is uniquely characterized by zero remainder;
  8. derivative values at ladder roots;
  9. primitive integer remainder forms;
 10. fresh-K quotient/remainder targets;
 11. exact reconstruction.

No arbitrary rational fitting is used.
No universal law is inferred from the supplied K-grid.
"""

from __future__ import annotations

import sympy as sp


# ============================================================
# SYMBOLS / GRID
# ============================================================

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


# ============================================================
# EXACT Y-COEFFICIENT ROWS
# ============================================================

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


# ============================================================
# BOUNDARY LADDERS
# ============================================================

LADDERS = {
    0: (K + 1) * (K + 2) * (K + 3) * (K + 4),
    1: (K + 2) * (K + 3) * (K + 4),
    2: (K + 3) * (K + 4),
    3: (K + 4),
    4: sp.Integer(1),
    5: (K - 5) * (K - 7) * (K - 9) * (K - 11),
}

LADDER_ROOTS = {
    0: [-1, -2, -3, -4],
    1: [-2, -3, -4],
    2: [-3, -4],
    3: [-4],
    4: [],
    5: [5, 7, 9, 11],
}

FRESH_K = [-9, -7, -5, -3, -1, 1, 15, 17, 19]


# ============================================================
# HELPERS
# ============================================================

def poly(expr):
    return sp.Poly(
        sp.expand(expr),
        K,
        domain=sp.QQ,
    )


def degree(expr):
    p = poly(expr)

    if p.is_zero:
        return -sp.oo

    return p.degree()


def degree_text(expr):
    d = degree(expr)

    if d == -sp.oo:
        return "-oo"

    return str(int(d))


def interpolate(values):
    pairs = [
        (K_VALUES[i], values[i])
        for i in range(len(K_VALUES))
    ]

    return sp.factor(
        sp.interpolate(pairs, K)
    )


def reconstruct_coefficients():
    coeffs = {}

    for j in J_VALUES:

        coeffs[j] = []

        for m in range(Q_DEGREE[j] + 1):

            values = [
                Y_ROWS[j][i][m]
                for i in range(len(K_VALUES))
            ]

            coeffs[j].append(
                interpolate(values)
            )

    return coeffs


def primitive_full_vector(values):
    vals = [
        sp.Rational(v)
        for v in values
    ]

    nonzero = [
        v for v in vals
        if v != 0
    ]

    if not nonzero:
        return [0] * len(vals)

    denominators = [
        int(v.q)
        for v in nonzero
    ]

    if len(denominators) == 1:
        common_den = denominators[0]
    else:
        common_den = int(
            sp.ilcm(*denominators)
        )

    ints = [
        int(v * common_den)
        for v in vals
    ]

    abs_nonzero = [
        abs(v)
        for v in ints
        if v != 0
    ]

    if len(abs_nonzero) == 1:
        content = abs_nonzero[0]
    else:
        content = int(
            sp.igcd(*abs_nonzero)
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


def polynomial_coeff_vector(expr):
    p = poly(expr)

    if p.is_zero:
        return [sp.Integer(0)]

    return [
        p.nth(i)
        for i in range(p.degree() + 1)
    ]


def matrix_from_polynomials(exprs):
    if not exprs:
        return sp.Matrix([])

    max_degree = 0

    for expr in exprs:
        p = poly(expr)

        if not p.is_zero:
            max_degree = max(
                max_degree,
                p.degree(),
            )

    rows = []

    for expr in exprs:
        p = poly(expr)

        rows.append(
            [
                p.nth(i)
                for i in range(max_degree + 1)
            ]
        )

    return sp.Matrix(rows)


# ============================================================
# 0. DATA VALIDATION
# ============================================================

def section_0():
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    checked = 0

    for j in J_VALUES:

        expected_width = (
            Q_DEGREE[j] + 1
        )

        if len(Y_ROWS[j]) != len(K_VALUES):
            raise RuntimeError(
                f"bad row count for j={j}: "
                f"{len(Y_ROWS[j])}"
            )

        for row in Y_ROWS[j]:

            if len(row) != expected_width:
                raise RuntimeError(
                    f"bad coefficient count "
                    f"for j={j}: {len(row)}"
                )

            checked += 1

    points = (
        len(K_VALUES)
        * len(J_VALUES)
        * len(D_VALUES)
    )

    print(f"k values = {K_VALUES}")
    print(f"j values = {J_VALUES}")
    print(f"D values = {D_VALUES}")
    print(f"points   = {points}")
    print("grid status = OK")
    print(
        f"checked coefficient rows = {checked}"
    )
    print()


# ============================================================
# 1. EXACT LADDER DIVISION
# ============================================================

def section_1_division(coeffs):
    print("=" * 78)
    print("1. EXACT LADDER DIVISION")
    print("=" * 78)

    divisions = {}

    for j in J_VALUES:

        ladder = sp.factor(
            LADDERS[j]
        )

        divisions[j] = []

        print(f"j={j}")
        print(f"  ladder={ladder}")
        print(
            f"  ladder_degree={degree_text(ladder)}"
        )

        for m in range(Q_DEGREE[j] + 1):

            a = coeffs[j][m]

            q, r = sp.div(
                poly(a),
                poly(ladder),
            )

            q = sp.factor(
                q.as_expr()
            )

            r = sp.factor(
                r.as_expr()
            )

            divisions[j].append(
                (q, r)
            )

            print(
                f"  m={m}: "
                f"deg(a)={degree_text(a)} "
                f"q_degree={degree_text(q)} "
                f"r_degree={degree_text(r)}"
            )
            print(
                f"    quotient={q}"
            )
            print(
                f"    remainder={r}"
            )

        print()

    return divisions


# ============================================================
# 2. REMAINDER VALUES AT LADDER ROOTS
# ============================================================

def section_2_remainder_root_values(
    divisions
):
    print("=" * 78)
    print(
        "2. REMAINDER VALUES AT BOUNDARY ROOTS"
    )
    print("=" * 78)

    for j in J_VALUES:

        roots = LADDER_ROOTS[j]

        print(f"j={j}")
        print(f"  roots={roots}")

        if not roots:
            print(
                "  no nontrivial ladder roots"
            )
            print()
            continue

        for m, (_, r) in enumerate(
            divisions[j]
        ):

            values = [
                sp.factor(
                    r.subs(K, root)
                )
                for root in roots
            ]

            print(
                f"  m={m}: r(root)={values}"
            )

        print()


# ============================================================
# 3. BOUNDARY-ROOT VALUE MATRICES
# ============================================================

def section_3_root_matrix(divisions):
    print("=" * 78)
    print(
        "3. BOUNDARY-ROOT VALUE MATRICES"
    )
    print("=" * 78)

    for j in J_VALUES:

        roots = LADDER_ROOTS[j]

        if not roots:
            continue

        rows = []

        for root in roots:

            row = []

            for _, r in divisions[j]:

                row.append(
                    sp.factor(
                        r.subs(K, root)
                    )
                )

            rows.append(row)

        M = sp.Matrix(rows)

        print(f"j={j}")
        print(f"  shape={M.shape}")
        print(f"  rank={M.rank()}")
        print(f"  matrix={M}")
        print()


# ============================================================
# 4. REMAINDER DEGREE PROFILE
# ============================================================

def section_4_remainder_degree_profile(
    divisions
):
    print("=" * 78)
    print(
        "4. REMAINDER / QUOTIENT DEGREE PROFILE"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(f"j={j}")

        for m, (q, r) in enumerate(
            divisions[j]
        ):

            print(
                f"  m={m}: "
                f"q_degree={degree_text(q)} "
                f"r_degree={degree_text(r)}"
            )

        print()


# ============================================================
# 5. REMAINDER COEFFICIENT SPACE RANK
# ============================================================

def section_5_remainder_coeff_rank(
    divisions
):
    print("=" * 78)
    print(
        "5. RANK OF REMAINDER COEFFICIENT MATRICES"
    )
    print("=" * 78)

    for j in J_VALUES:

        ladder_degree = degree(
            LADDERS[j]
        )

        if ladder_degree <= 0:
            print(
                f"j={j}: "
                "ladder degree 0; "
                "remainders are identically zero"
            )
            continue

        exprs = [
            r
            for _, r in divisions[j]
        ]

        M = matrix_from_polynomials(
            exprs
        )

        print(
            f"j={j}: "
            f"shape={M.shape} "
            f"rank={M.rank()}"
        )

    print()


# ============================================================
# 6. PRIMITIVE REMAINDER FORMS
# ============================================================

def section_6_remainder_primitive_forms(
    divisions
):
    print("=" * 78)
    print(
        "6. PRIMITIVE INTEGER FORMS OF REMAINDERS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(f"j={j}")

        for m, (_, r) in enumerate(
            divisions[j]
        ):

            coeffs = (
                polynomial_coeff_vector(r)
            )

            primitive = (
                primitive_full_vector(
                    coeffs
                )
            )

            print(
                f"  m={m}:"
            )
            print(
                f"    remainder={r}"
            )
            print(
                f"    coeffs={coeffs}"
            )
            print(
                f"    primitive={primitive}"
            )

        print()


# ============================================================
# 7. M=0 EXCEPTION AUDIT
# ============================================================

def section_7_m0_exception_audit(
    divisions
):
    print("=" * 78)
    print("7. m=0 EXCEPTION AUDIT")
    print("=" * 78)

    for j in J_VALUES:

        print(f"j={j}")

        for m, (_, r) in enumerate(
            divisions[j]
        ):

            is_zero = (
                sp.expand(r) == 0
            )

            print(
                f"  m={m}: "
                f"remainder_zero={is_zero}"
            )

        print()


# ============================================================
# 8. ROOT-VALUE RATIOS
# ============================================================

def section_8_root_value_ratios(
    divisions
):
    print("=" * 78)
    print(
        "8. ROOT-VALUE RATIOS ACROSS m"
    )
    print("=" * 78)

    for j in J_VALUES:

        roots = LADDER_ROOTS[j]

        if not roots:
            continue

        print(f"j={j}")

        for root in roots:

            values = [
                sp.factor(
                    r.subs(K, root)
                )
                for _, r in divisions[j]
            ]

            print(
                f"  K={root}: values={values}"
            )

            base = None

            for value in values:
                if value != 0:
                    base = value
                    break

            if base is None:
                continue

            ratios = []

            for value in values:

                if value == 0:
                    ratios.append("zero")
                else:
                    ratios.append(
                        sp.factor(
                            value / base
                        )
                    )

            print(
                "    normalized_to_first_nonzero="
                f"{ratios}"
            )

        print()


# ============================================================
# 9. CROSS-J REMAINDER SPACE RANK
# ============================================================

def section_9_cross_j_remainder_rank(
    divisions
):
    print("=" * 78)
    print(
        "9. CROSS-j REMAINDER SPACE RANK"
    )
    print("=" * 78)

    for m in range(6):

        rows = []
        active = []

        for j in J_VALUES:

            if m > Q_DEGREE[j]:
                continue

            r = divisions[j][m][1]

            if r == 0:
                continue

            rows.append(r)
            active.append(j)

        if not rows:
            print(
                f"m={m}: no nonzero remainders"
            )
            continue

        M = matrix_from_polynomials(
            rows
        )

        print(
            f"m={m}: "
            f"active_j={active} "
            f"shape={M.shape} "
            f"rank={M.rank()}/{len(active)}"
        )

    print()


# ============================================================
# 10. CROSS-J QUOTIENT SPACE RANK
# ============================================================

def section_10_quotient_rank(
    divisions
):
    print("=" * 78)
    print(
        "10. CROSS-j QUOTIENT SPACE RANK"
    )
    print("=" * 78)

    for m in range(6):

        rows = []
        active = []

        for j in J_VALUES:

            if m > Q_DEGREE[j]:
                continue

            q = divisions[j][m][0]

            if q == 0:
                continue

            rows.append(q)
            active.append(j)

        if not rows:
            print(
                f"m={m}: no nonzero quotients"
            )
            continue

        M = matrix_from_polynomials(
            rows
        )

        print(
            f"m={m}: "
            f"active_j={active} "
            f"shape={M.shape} "
            f"rank={M.rank()}/{len(active)}"
        )

    print()


# ============================================================
# 11. DIRECT VALUE / DERIVATIVE AUDIT
# ============================================================

def section_11_derivative_root_audit(
    coeffs
):
    print("=" * 78)
    print(
        "11. DIRECT VALUE / DERIVATIVE AUDIT AT LADDER ROOTS"
    )
    print("=" * 78)

    for j in J_VALUES:

        roots = LADDER_ROOTS[j]

        if not roots:
            continue

        print(f"j={j}")

        for root in roots:

            print(
                f"  K={root}"
            )

            for m in range(
                Q_DEGREE[j] + 1
            ):

                a = coeffs[j][m]

                value = sp.factor(
                    a.subs(K, root)
                )

                derivative = sp.factor(
                    sp.diff(a, K).subs(
                        K, root
                    )
                )

                print(
                    f"    m={m}: "
                    f"value={value}, "
                    f"first_derivative={derivative}"
                )

        print()


# ============================================================
# 12. EXACT DIVISION RECONSTRUCTION
# ============================================================

def section_12_division_reconstruction(
    divisions,
    coeffs
):
    print("=" * 78)
    print(
        "12. EXACT DIVISION RECONSTRUCTION"
    )
    print("=" * 78)

    tested = 0
    failures = []

    for j in J_VALUES:

        ladder = LADDERS[j]

        for m in range(
            Q_DEGREE[j] + 1
        ):

            q, r = divisions[j][m]

            recovered = sp.factor(
                sp.expand(
                    ladder * q + r
                )
            )

            expected = sp.factor(
                coeffs[j][m]
            )

            tested += 1

            if sp.expand(
                recovered - expected
            ) != 0:

                failures.append(
                    (
                        j,
                        m,
                        expected,
                        recovered,
                    )
                )

    print(f"tested={tested}")
    print(
        "division reconstruction failures="
        f"{len(failures)}"
    )

    for item in failures[:20]:
        print(f"  {item}")

    print()


# ============================================================
# 13. FRESH-K TARGETS
# ============================================================

def section_13_fresh_k_targets(
    divisions
):
    print("=" * 78)
    print(
        "13. FRESH-K QUOTIENT / REMAINDER TARGETS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(f"j={j}")

        for m, (q, r) in enumerate(
            divisions[j]
        ):

            print(f"  m={m}:")

            for k0 in FRESH_K:

                q_value = sp.factor(
                    q.subs(K, k0)
                )

                r_value = sp.factor(
                    r.subs(K, k0)
                )

                print(
                    f"    K={k0}: "
                    f"q={q_value}, "
                    f"r={r_value}"
                )

        print()


# ============================================================
# 14. NORMALIZED DEFECT
# ============================================================

def section_14_normalized_defect(
    divisions
):
    print("=" * 78)
    print(
        "14. NORMALIZED BOUNDARY DEFECT"
    )
    print("=" * 78)

    for j in J_VALUES:

        ladder = LADDERS[j]
        ladder_degree = degree(ladder)

        print(f"j={j}")
        print(f"  ladder={ladder}")
        print(
            f"  ladder_degree="
            f"{degree_text(ladder)}"
        )

        for m, (q, r) in enumerate(
            divisions[j]
        ):

            print(
                f"  m={m}:"
            )
            print(
                f"    quotient={q}"
            )
            print(
                f"    remainder={r}"
            )
            print(
                f"    q_degree={degree_text(q)}"
            )
            print(
                f"    r_degree={degree_text(r)}"
            )

            if (
                ladder_degree > 0
                and q != 0
            ):

                qpoly = poly(q)

                print(
                    f"    quotient_leading="
                    f"{qpoly.LC()}"
                )

        print()


# ============================================================
# 15. COMPACT STRUCTURAL SUMMARY
# ============================================================

def section_15_structural_summary(
    divisions
):
    print("=" * 78)
    print(
        "15. COMPACT STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    print(
        "j | ladder_degree | m | deg(a) | "
        "deg(q) | deg(r) | remainder_zero"
    )
    print("-" * 78)

    for j in J_VALUES:

        ladder_degree = degree(
            LADDERS[j]
        )

        for m, (q, r) in enumerate(
            divisions[j]
        ):

            deg_a = (
                degree(
                    q * LADDERS[j] + r
                )
            )

            print(
                f"{j} | "
                f"{int(ladder_degree):13d} | "
                f"{m} | "
                f"{degree_text(deg_a):>6} | "
                f"{degree_text(q):>6} | "
                f"{degree_text(r):>6} | "
                f"{sp.expand(r) == 0}"
            )

    print()


# ============================================================
# FINAL DIAGNOSTIC
# ============================================================

def final_diagnostic():
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The Newton experiment produced a sharp split:"
    )
    print(
        "  * a_{j,0} is sparse in a boundary-adapted basis;"
    )
    print(
        "  * higher a_{j,m}, m>0, remain fully supported."
    )
    print()

    print(
        "This experiment therefore studies the exact defect"
    )
    print(
        "from the known boundary ladder:"
    )
    print()
    print(
        "    a_{j,m} = L_j q_{j,m} + r_{j,m}"
    )
    print()

    print(
        "The central questions are:"
    )
    print(
        "  1. Does r_{j,m} have unusually small degree?"
    )
    print(
        "  2. Are the root-value matrices low-rank?"
    )
    print(
        "  3. Is m=0 uniquely characterized by r=0?"
    )
    print(
        "  4. Do the quotients q_{j,m} simplify?"
    )
    print(
        "  5. Do derivative values at the ladder roots"
    )
    print(
        "     exhibit a smaller structure?"
    )
    print(
        "  6. Does anything survive at fresh K values?"
    )
    print()

    print(
        "This is an exact divisibility/defect experiment."
    )
    print(
        "No arbitrary rational fit is treated as evidence."
    )
    print()

    print(
        "No universal R=5 law is inferred."
    )
    print("No r=6.")
    print("No full pq-kernel expansion.")
    print("No replacement universal r,j formula.")

    print("=" * 78)


# ============================================================
# MAIN
# ============================================================

def main():

    section_0()

    coeffs = reconstruct_coefficients()

    divisions = section_1_division(
        coeffs
    )

    section_2_remainder_root_values(
        divisions
    )

    section_3_root_matrix(
        divisions
    )

    section_4_remainder_degree_profile(
        divisions
    )

    section_5_remainder_coeff_rank(
        divisions
    )

    section_6_remainder_primitive_forms(
        divisions
    )

    section_7_m0_exception_audit(
        divisions
    )

    section_8_root_value_ratios(
        divisions
    )

    section_9_cross_j_remainder_rank(
        divisions
    )

    section_10_quotient_rank(
        divisions
    )

    section_11_derivative_root_audit(
        coeffs
    )

    section_12_division_reconstruction(
        divisions,
        coeffs
    )

    section_13_fresh_k_targets(
        divisions
    )

    section_14_normalized_defect(
        divisions
    )

    section_15_structural_summary(
        divisions
    )

    final_diagnostic()


if __name__ == "__main__":
    main()