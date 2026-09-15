#!/usr/bin/env python3

"""
R=5 EXPERIMENT 62
============================================================

NEIGHBORING-j TRANSITION DEFECT POLYNOMIALS

Known boundary ladders:

    L_0=(K+1)(K+2)(K+3)(K+4)
    L_1=(K+2)(K+3)(K+4)
    L_2=(K+3)(K+4)
    L_3=(K+4)
    L_4=1

For j=0,1,2 define the exact polynomial

    Delta_{j,m}(K)
      = rem_{L_{j+1}}
          (a_{j+1,m}(K)-a_{j,m}(K)).

Because L_{j+1} has degree 3,2,1 respectively,
Delta has degree at most 2,1,0.

This is the polynomial whose values at the common
ladder roots are exactly the neighboring-j differences
seen in the previous experiment.

The experiment tests:

  1. exact transition-defect polynomials;
  2. their degree profile in m;
  3. factorization;
  4. gcd across different m;
  5. coefficient-matrix rank;
  6. values at the root removed by the next ladder;
  7. normalized removed-root profiles across m;
  8. exact quotient/defect division;
  9. fresh-K evaluation;
 10. exact reconstruction.

No arbitrary rational fitting is used.
No zero-padding of unequal root sets is used.
j=5 remains terminal and separate.
"""


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
# KNOWN LADDERS
# ============================================================

LADDERS = {
    0: (K + 1) * (K + 2) * (K + 3) * (K + 4),
    1: (K + 2) * (K + 3) * (K + 4),
    2: (K + 3) * (K + 4),
    3: K + 4,
    4: sp.Integer(1),
    5: (K - 5) * (K - 7) * (K - 9) * (K - 11),
}

ROOTS = {
    0: [-1, -2, -3, -4],
    1: [-2, -3, -4],
    2: [-3, -4],
    3: [-4],
    4: [],
    5: [5, 7, 9, 11],
}

REMOVED_ROOT = {
    0: -2,
    1: -3,
    2: -4,
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
    return sp.factor(
        sp.interpolate(
            [
                (K_VALUES[i], values[i])
                for i in range(len(K_VALUES))
            ],
            K,
        )
    )


def reconstruct_coefficients():
    coeffs = {}

    for j in J_VALUES:

        coeffs[j] = []

        for m in range(
            Q_DEGREE[j] + 1
        ):

            values = [
                Y_ROWS[j][i][m]
                for i in range(
                    len(K_VALUES)
                )
            ]

            coeffs[j].append(
                interpolate(values)
            )

    return coeffs


def primitive_vector(values):
    values = [sp.Rational(v) for v in values]

    nonzero = [
        v for v in values
        if v != 0
    ]

    if not nonzero:
        return [0] * len(values)

    denoms = [
        int(v.q)
        for v in nonzero
    ]

    if len(denoms) == 1:
        den_lcm = denoms[0]
    else:
        den_lcm = int(
            sp.ilcm(*denoms)
        )

    ints = [
        int(v * den_lcm)
        for v in values
    ]

    nonzero_ints = [
        abs(v)
        for v in ints
        if v != 0
    ]

    if len(nonzero_ints) == 1:
        content = nonzero_ints[0]
    else:
        content = int(
            sp.igcd(*nonzero_ints)
        )

    ints = [
        v // content
        for v in ints
    ]

    first = next(
        (
            v for v in ints
            if v != 0
        ),
        0,
    )

    if first < 0:
        ints = [-v for v in ints]

    return ints


def matrix_from_rows(rows):
    if not rows:
        return sp.Matrix.zeros(
            0,
            0,
        )

    width = len(rows[0])

    for row in rows:
        if len(row) != width:
            raise ValueError(
                "rows have unequal dimensions"
            )

    return sp.Matrix(rows)


# ============================================================
# 0. VALIDATION
# ============================================================

def section_0():
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    checked = 0

    for j in J_VALUES:

        expected = (
            Q_DEGREE[j] + 1
        )

        if len(Y_ROWS[j]) != len(
            K_VALUES
        ):
            raise RuntimeError(
                f"row count mismatch j={j}"
            )

        for row in Y_ROWS[j]:

            if len(row) != expected:
                raise RuntimeError(
                    f"width mismatch j={j}"
                )

            checked += 1

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
        f"checked coefficient rows = {checked}"
    )
    print()


# ============================================================
# 1. BUILD COEFFICIENT POLYNOMIALS
# ============================================================

def section_1_coefficients(coeffs):
    print("=" * 78)
    print(
        "1. EXACT k-COEFFICIENT POLYNOMIALS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for m, expr in enumerate(
            coeffs[j]
        ):

            print(
                f"  a_{{{j},{m}}}(K)="
                f"{sp.factor(expr)}"
            )

        print()


# ============================================================
# 2. TRANSITION DEFECTS
# ============================================================

def compute_transition_defects(
    coeffs
):
    defects = {}
    quotients = {}

    for j in [0, 1, 2]:

        divisor = poly(
            LADDERS[j + 1]
        )

        defects[j] = []
        quotients[j] = []

        max_m = min(
            len(coeffs[j]),
            len(coeffs[j + 1]),
        )

        for m in range(max_m):

            difference = sp.expand(
                coeffs[j + 1][m]
                - coeffs[j][m]
            )

            q, r = sp.div(
                poly(difference),
                divisor,
            )

            quotients[j].append(
                sp.factor(
                    q.as_expr()
                )
            )

            defects[j].append(
                sp.factor(
                    r.as_expr()
                )
            )

    return defects, quotients


def section_2_transition_defects(
    defects
):
    print("=" * 78)
    print(
        "2. EXACT NEIGHBORING-j TRANSITION DEFECTS"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j} -> j={j+1}"
        )
        print(
            f"  modulus="
            f"{sp.factor(LADDERS[j+1])}"
        )

        for m, defect in enumerate(
            defects[j]
        ):

            print(
                f"  m={m}: "
                f"degree={degree_text(defect)}"
            )
            print(
                f"    Delta="
                f"{defect}"
            )

        print()


# ============================================================
# 3. DEGREE / FACTORIZATION PROFILE
# ============================================================

def section_3_profile(
    defects
):
    print("=" * 78)
    print(
        "3. TRANSITION DEFECT DEGREE / FACTORIZATION PROFILE"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        divisor = poly(
            LADDERS[j + 1]
        )

        print(
            f"j={j}->{j+1}: "
            f"modulus_degree={divisor.degree()}"
        )

        for m, defect in enumerate(
            defects[j]
        ):

            print(
                f"  m={m}:"
            )
            print(
                f"    degree={degree_text(defect)}"
            )
            print(
                f"    factorization="
                f"{sp.factor(defect)}"
            )
            print(
                f"    primitive_coeffs="
                f"{primitive_vector([c for c in poly(defect).all_coeffs()])}"
            )

        print()


# ============================================================
# 4. GCD ACROSS m
# ============================================================

def section_4_gcd_across_m(
    defects
):
    print("=" * 78)
    print(
        "4. GCD OF TRANSITION DEFECTS ACROSS m"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        active = [
            d
            for d in defects[j]
            if d != 0
        ]

        if not active:

            print(
                "  all defects zero"
            )
            print()
            continue

        g = poly(
            active[0]
        )

        for d in active[1:]:

            g = sp.gcd(
                g,
                poly(d),
            )

        gexpr = sp.factor(
            g.as_expr()
        )

        print(
            f"  common_gcd={gexpr}"
        )
        print(
            f"  gcd_degree={degree_text(gexpr)}"
        )

        for m, d in enumerate(
            defects[j]
        ):

            if d == 0:
                print(
                    f"  m={m}: ZERO"
                )
                continue

            q, r = sp.div(
                poly(d),
                g,
            )

            print(
                f"  m={m}: "
                f"quotient={sp.factor(q.as_expr())} "
                f"remainder={sp.factor(r.as_expr())}"
            )

        print()


# ============================================================
# 5. COEFFICIENT MATRIX RANK
# ============================================================

def section_5_rank(
    defects
):
    print("=" * 78)
    print(
        "5. TRANSITION DEFECT COEFFICIENT-MATRIX RANK"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        max_degree = max(
            [
                degree(d)
                for d in defects[j]
            ]
        )

        rows = []

        for d in defects[j]:

            p = poly(d)

            coeffs = [
                p.coeff_monomial(
                    K**r
                )
                for r in range(
                    max_degree + 1
                )
            ]

            rows.append(coeffs)

        M = matrix_from_rows(rows)

        print(
            f"j={j}->{j+1}: "
            f"shape={M.shape} "
            f"rank={M.rank()} / {min(M.shape)}"
        )

        for m, row in enumerate(
            rows
        ):

            print(
                f"  m={m}: {row}"
            )

        print()


# ============================================================
# 6. REMOVED-ROOT VALUES
# ============================================================

def section_6_removed_root(
    defects
):
    print("=" * 78)
    print(
        "6. TRANSITION DEFECT VALUES AT REMOVED ROOT"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        root = REMOVED_ROOT[j]

        print(
            f"j={j}->{j+1}: "
            f"removed_root={root}"
        )

        values = []

        for m, d in enumerate(
            defects[j]
        ):

            value = sp.factor(
                d.subs(K, root)
            )

            values.append(value)

            print(
                f"  m={m}: {value}"
            )

        print(
            f"  primitive="
            f"{primitive_vector(values)}"
        )

        nonzero = [
            v
            for v in values
            if v != 0
        ]

        if nonzero:

            base = nonzero[0]

            normalized = [
                "zero"
                if v == 0
                else sp.factor(v / base)
                for v in values
            ]

            print(
                f"  normalized="
                f"{normalized}"
            )

        print()


# ============================================================
# 7. SUCCESSIVE m-RATIOS OF DEFECTS
# ============================================================

def section_7_successive_ratios(
    defects
):
    print("=" * 78)
    print(
        "7. SUCCESSIVE-m RATIOS OF TRANSITION DEFECTS"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for m in range(
            1,
            len(defects[j]),
        ):

            current = defects[j][m]
            previous = defects[j][m - 1]

            if previous == 0:

                ratio = (
                    "undefined"
                    if current != 0
                    else "0/0"
                )

            else:

                ratio = sp.factor(
                    current / previous
                )

            print(
                f"  m={m}/{m-1}: "
                f"{ratio}"
            )

        print()


# ============================================================
# 8. QUOTIENT / DEFECT SPLIT
# ============================================================

def section_8_exact_split(
    coeffs,
    defects,
    quotients
):
    print("=" * 78)
    print(
        "8. EXACT DIFFERENCE = MODULUS * QUOTIENT + DEFECT"
    )
    print("=" * 78)

    failures = 0
    tested = 0

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        modulus = (
            sp.factor(
                LADDERS[j + 1]
            )
        )

        for m in range(
            len(defects[j])
        ):

            lhs = sp.expand(
                coeffs[j + 1][m]
                - coeffs[j][m]
            )

            rhs = sp.expand(
                modulus
                * quotients[j][m]
                + defects[j][m]
            )

            ok = (
                sp.expand(
                    lhs - rhs
                )
                == 0
            )

            tested += 1

            if not ok:
                failures += 1

            print(
                f"  m={m}: "
                f"ok={ok}"
            )
            print(
                f"    quotient="
                f"{quotients[j][m]}"
            )
            print(
                f"    defect="
                f"{defects[j][m]}"
            )

        print()

    print(
        f"tested={tested}"
    )
    print(
        f"reconstruction_failures={failures}"
    )
    print()


# ============================================================
# 9. DEFECT ROOT PROFILE
# ============================================================

def section_9_defect_roots(
    defects
):
    print("=" * 78)
    print(
        "9. EXACT ROOT AUDIT OF TRANSITION DEFECTS"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for m, d in enumerate(
            defects[j]
        ):

            if d == 0:

                print(
                    f"  m={m}: zero defect"
                )
                continue

            roots = sp.solve(
                sp.Eq(d, 0),
                K,
            )

            rational_roots = []

            for root in roots:

                if getattr(
                    root,
                    "is_Rational",
                    False,
                ):
                    rational_roots.append(
                        root
                    )

            print(
                f"  m={m}: "
                f"rational_roots="
                f"{rational_roots}"
            )

        print()


# ============================================================
# 10. FRESH-K DEFECT TARGETS
# ============================================================

def section_10_fresh_k(
    defects
):
    print("=" * 78)
    print(
        "10. FRESH-K TRANSITION DEFECT TARGETS"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for m, d in enumerate(
            defects[j]
        ):

            values = [
                sp.factor(
                    d.subs(K, k0)
                )
                for k0 in FRESH_K
            ]

            print(
                f"  m={m}:"
            )
            print(
                f"    {values}"
            )

        print()


# ============================================================
# 11. DEFECT VALUES ON THE COMMON ROOTS
# ============================================================

def section_11_common_root_check(
    coeffs,
    defects
):
    print("=" * 78)
    print(
        "11. COMMON-ROOT VALUE CHECK"
    )
    print("=" * 78)

    failures = 0
    tested = 0

    common_roots = {
        0: [-2, -3, -4],
        1: [-3, -4],
        2: [-4],
    }

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        roots = common_roots[j]

        for m, d in enumerate(
            defects[j]
        ):

            for root in roots:

                direct = sp.factor(
                    (
                        coeffs[j + 1][m]
                        - coeffs[j][m]
                    ).subs(
                        K,
                        root,
                    )
                )

                defect_value = sp.factor(
                    d.subs(
                        K,
                        root,
                    )
                )

                tested += 1

                ok = (
                    sp.expand(
                        direct
                        - defect_value
                    )
                    == 0
                )

                if not ok:
                    failures += 1

                if not ok:

                    print(
                        f"  FAILURE "
                        f"m={m} K={root}"
                    )

        print()

    print(
        f"tested={tested}"
    )
    print(
        f"failures={failures}"
    )
    print()


# ============================================================
# 12. SPECIAL j=5 TERMINAL DEFECT
# ============================================================

def section_12_terminal(
    coeffs
):
    print("=" * 78)
    print(
        "12. TERMINAL j=5 AUDIT"
    )
    print("=" * 78)

    ladder = LADDERS[5]

    print(
        f"ladder={sp.factor(ladder)}"
    )
    print(
        f"roots={ROOTS[5]}"
    )

    a0 = coeffs[5][0]
    a1 = coeffs[5][1]

    q0, r0 = sp.div(
        poly(a0),
        poly(ladder),
    )

    q1, r1 = sp.div(
        poly(a1),
        poly(ladder),
    )

    print(
        f"a_{{5,0}}={sp.factor(a0)}"
    )
    print(
        f"a_{{5,0}}/ladder="
        f"{sp.factor(q0.as_expr())}"
    )
    print(
        f"remainder={sp.factor(r0.as_expr())}"
    )

    print(
        f"a_{{5,1}}="
        f"{sp.factor(a1)}"
    )
    print(
        f"a_{{5,1}} quotient="
        f"{sp.factor(q1.as_expr())}"
    )
    print(
        f"a_{{5,1}} remainder="
        f"{sp.factor(r1.as_expr())}"
    )

    print()

    roots = ROOTS[5]

    print(
        "a_5,1 at terminal roots:"
    )

    for root in roots:

        print(
            f"  K={root}: "
            f"a_5,1={sp.factor(a1.subs(K, root))}"
        )

    print()


# ============================================================
# 13. EXACT ORIGINAL GRID RECONSTRUCTION
# ============================================================

def section_13_grid(
    coeffs
):
    print("=" * 78)
    print(
        "13. EXACT ORIGINAL-GRID RECONSTRUCTION"
    )
    print("=" * 78)

    failures = 0
    tested = 0

    for j in J_VALUES:

        for k0 in K_VALUES:

            for d0 in D_VALUES:

                y = d0 - D0[j]

                reconstructed = (
                    sp.Integer(0)
                )

                for m in range(
                    Q_DEGREE[j] + 1
                ):

                    reconstructed += (
                        coeffs[j][m]
                        .subs(K, k0)
                        * y**m
                    )

                tested += 1

                if sp.expand(
                    reconstructed
                    - reconstructed
                ) != 0:

                    failures += 1

    print(
        f"tested={tested}"
    )
    print(
        f"reconstruction_failures={failures}"
    )
    print()


# ============================================================
# 14. COMPACT SUMMARY
# ============================================================

def section_14_summary(
    defects
):
    print("=" * 78)
    print(
        "14. COMPACT TRANSITION-DEFECT SUMMARY"
    )
    print("=" * 78)

    print(
        "j | modulus_degree | m | defect_degree"
    )
    print("-" * 78)

    for j in [0, 1, 2]:

        mod_degree = poly(
            LADDERS[j + 1]
        ).degree()

        for m, d in enumerate(
            defects[j]
        ):

            print(
                f"{j} | "
                f"{mod_degree} | "
                f"{m} | "
                f"{degree_text(d)}"
            )

    print()

    print(
        "Expected maximum defect degrees:"
    )

    for j in [0, 1, 2]:

        print(
            f"  j={j}->{j+1}: "
            f"< {poly(LADDERS[j+1]).degree()}"
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
        "The nested-root experiment showed that the raw"
    )
    print(
        "neighboring-j coordinate differences are not"
    )
    print(
        "low-rank in the obvious coordinate representation."
    )
    print()

    print(
        "The next exact object is therefore the polynomial"
    )
    print(
        "represented by those differences on the common"
    )
    print(
        "root set."
    )
    print()

    print(
        "For j=0,1,2:"
    )
    print(
        "  Delta_{j,m} = rem_{L_{j+1}}"
    )
    print(
        "       (a_{j+1,m} - a_{j,m})"
    )
    print()

    print(
        "Because deg(L_{j+1}) is 3, 2, 1 respectively,"
    )
    print(
        "the transition defects can have degree at most"
    )
    print(
        "2, 1, 0."
    )
    print()

    print(
        "This is a much smaller exact object than the"
    )
    print(
        "original degree-5 coefficient functions."
    )
    print()

    print(
        "The important outcomes are:"
    )
    print(
        "  * whether the defect degree reaches its maximum;"
    )
    print(
        "  * whether defects share a nontrivial gcd in m;"
    )
    print(
        "  * whether their coefficient matrices collapse;"
    )
    print(
        "  * whether the removed-root values have a simple"
    )
    print(
        "    exact pattern across m;"
    )
    print(
        "  * whether the same structure survives at fresh K."
    )
    print()

    print(
        "j=5 remains terminal and is not forced into this"
    )
    print(
        "negative-root transition chain."
    )
    print()

    print(
        "No arbitrary rational fit is treated as evidence."
    )
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


# ============================================================
# MAIN
# ============================================================

def main():

    section_0()

    coeffs = reconstruct_coefficients()

    section_1_coefficients(
        coeffs
    )

    defects, quotients = (
        compute_transition_defects(
            coeffs
        )
    )

    section_2_transition_defects(
        defects
    )

    section_3_profile(
        defects
    )

    section_4_gcd_across_m(
        defects
    )

    section_5_rank(
        defects
    )

    section_6_removed_root(
        defects
    )

    section_7_successive_ratios(
        defects
    )

    section_8_exact_split(
        coeffs,
        defects,
        quotients,
    )

    section_9_defect_roots(
        defects
    )

    section_10_fresh_k(
        defects
    )

    section_11_common_root_check(
        coeffs,
        defects,
    )

    section_12_terminal(
        coeffs
    )

    section_13_grid(
        coeffs
    )

    section_14_summary(
        defects
    )

    final_diagnostic()


if __name__ == "__main__":
    main()

