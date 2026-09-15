#!/usr/bin/env python3

"""
R=5 EXPERIMENT 61
============================================================

NESTED LADDER-ROOT RESTRICTION EXPERIMENT

Previous result:

    a_{j,m}(K) = L_j(K) q_{j,m}(K) + r_{j,m}(K)

with

    deg r_{j,m} < deg L_j.

The remainder is therefore completely determined by its values
at the roots of L_j.

For j=0,...,3 the ladder roots are nested:

    j=0 : {-1,-2,-3,-4}
    j=1 :    {-2,-3,-4}
    j=2 :       {-3,-4}
    j=3 :          {-4}

This experiment therefore does NOT place unequal coordinate
vectors into one rectangular matrix.

Instead it studies:

  1. exact root-coordinate representation;
  2. restriction of j-level coordinates to the next j-level;
  3. differences between restricted neighboring-j coordinates;
  4. ratios between restricted neighboring-j coordinates;
  5. whether those differences/ratios have low rank across m;
  6. aligned coordinate matrices on the common negative-root set;
  7. the terminal j=5 ladder separately;
  8. exact reconstruction;
  9. fresh-K validation.

For j=0,...,3 define

    Roots(j+1) subset Roots(j).

Thus every remainder coordinate vector at j+1 can be compared
directly with the corresponding entries of the j-vector.

This is a stricter structural test than zero-padding.

No arbitrary rational fitting is used.
No universal R=5 law is inferred from interpolation.
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
# LADDERS / ROOTS
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

NEGATIVE_CHAIN = {
    0: [-1, -2, -3, -4],
    1: [-2, -3, -4],
    2: [-3, -4],
    3: [-4],
}

FRESH_K = [-9, -7, -5, -3, -1, 1, 15, 17, 19]


# ============================================================
# BASIC HELPERS
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

        for m in range(
            Q_DEGREE[j] + 1
        ):

            values = [
                Y_ROWS[j][i][m]
                for i in range(len(K_VALUES))
            ]

            coeffs[j].append(
                interpolate(values)
            )

    return coeffs


def primitive_vector(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

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
        common_den = denoms[0]
    else:
        common_den = int(
            sp.ilcm(*denoms)
        )

    ints = [
        int(v * common_den)
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
        (v for v in ints if v != 0),
        0,
    )

    if first < 0:
        ints = [-v for v in ints]

    return ints


def matrix_from_rows(rows):
    """
    Safe matrix constructor.

    Every row must have the same width.
    """
    if not rows:
        return sp.Matrix.zeros(0, 0)

    width = len(rows[0])

    for row in rows:
        if len(row) != width:
            raise ValueError(
                "matrix_from_rows received "
                "rows of unequal width"
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

        expected_width = (
            Q_DEGREE[j] + 1
        )

        if len(Y_ROWS[j]) != len(
            K_VALUES
        ):
            raise RuntimeError(
                f"row count mismatch at j={j}"
            )

        for row in Y_ROWS[j]:

            if len(row) != expected_width:
                raise RuntimeError(
                    f"coefficient width "
                    f"mismatch at j={j}"
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
# 1. EXACT LADDER DIVISION
# ============================================================

def divide_by_ladder(coeffs):
    divisions = {}

    for j in J_VALUES:

        divisor = poly(
            LADDERS[j]
        )

        divisions[j] = []

        for m, expr in enumerate(
            coeffs[j]
        ):

            q, r = sp.div(
                poly(expr),
                divisor,
            )

            divisions[j].append(
                (
                    sp.factor(q.as_expr()),
                    sp.factor(r.as_expr()),
                )
            )

    return divisions


def section_1(divisions):
    print("=" * 78)
    print(
        "1. EXACT LADDER DIVISION"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )
        print(
            f"  ladder={sp.factor(LADDERS[j])}"
        )
        print(
            f"  roots={ROOTS[j]}"
        )

        for m, (q, r) in enumerate(
            divisions[j]
        ):

            print(
                f"  m={m}: "
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


# ============================================================
# 2. ROOT COORDINATES
# ============================================================

def root_coordinates(divisions):
    coordinates = {}

    for j in J_VALUES:

        roots = ROOTS[j]

        coordinates[j] = []

        if not roots:
            continue

        for q, r in divisions[j]:

            coords = [
                sp.factor(
                    r.subs(K, root)
                )
                for root in roots
            ]

            coordinates[j].append(
                coords
            )

    return coordinates


def section_2_coordinates(
    coordinates
):
    print("=" * 78)
    print(
        "2. EXACT ROOT-COORDINATE MATRICES"
    )
    print("=" * 78)

    for j in J_VALUES:

        roots = ROOTS[j]

        print(
            f"j={j}"
        )

        if not roots:
            print(
                "  no ladder roots"
            )
            print()
            continue

        print(
            f"  roots={roots}"
        )

        for m, row in enumerate(
            coordinates[j]
        ):

            print(
                f"  m={m}:"
            )
            print(
                f"    coordinates={row}"
            )
            print(
                f"    primitive={primitive_vector(row)}"
            )

        print()


# ============================================================
# 3. ROOT-COORDINATE RANK WITHIN EACH j
# ============================================================

def section_3_within_j_rank(
    coordinates
):
    print("=" * 78)
    print(
        "3. ROOT-COORDINATE RANK WITHIN EACH j"
    )
    print("=" * 78)

    for j in J_VALUES:

        if not coordinates[j]:

            print(
                f"j={j}: no coordinates"
            )
            continue

        M = matrix_from_rows(
            coordinates[j]
        )

        print(
            f"j={j}: "
            f"shape={M.shape} "
            f"rank={M.rank()}"
        )

    print()


# ============================================================
# 4. NESTED ROOT RESTRICTION
# ============================================================

def restriction_to_roots(
    row,
    source_roots,
    target_roots,
):
    """
    Exact restriction of a coordinate vector
    from source roots to a subset of those roots.
    """
    positions = []

    for root in target_roots:

        if root not in source_roots:
            raise ValueError(
                "target root is not contained "
                "in source root set"
            )

        positions.append(
            source_roots.index(root)
        )

    return [
        row[pos]
        for pos in positions
    ]


def section_4_nested_restrictions(
    coordinates
):
    print("=" * 78)
    print(
        "4. NESTED ROOT RESTRICTIONS j -> j+1"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        source_roots = (
            NEGATIVE_CHAIN[j]
        )
        target_roots = (
            NEGATIVE_CHAIN[j + 1]
        )

        print(
            f"j={j} -> j={j+1}"
        )
        print(
            f"  source_roots={source_roots}"
        )
        print(
            f"  target_roots={target_roots}"
        )

        for m in range(
            min(
                len(coordinates[j]),
                len(coordinates[j + 1]),
            )
        ):

            restricted = (
                restriction_to_roots(
                    coordinates[j][m],
                    source_roots,
                    target_roots,
                )
            )

            target = (
                coordinates[j + 1][m]
            )

            print(
                f"  m={m}"
            )
            print(
                f"    restricted_from_j={restricted}"
            )
            print(
                f"    native_next_j={target}"
            )

            delta = [
                sp.factor(
                    target[i]
                    - restricted[i]
                )
                for i in range(
                    len(target)
                )
            ]

            print(
                f"    difference={delta}"
            )

        print()


# ============================================================
# 5. NEIGHBORING-j DIFFERENCE MATRICES
# ============================================================

def section_5_difference_rank(
    coordinates
):
    print("=" * 78)
    print(
        "5. NEIGHBORING-j ROOT DIFFERENCE RANK"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        source_roots = NEGATIVE_CHAIN[j]
        target_roots = NEGATIVE_CHAIN[j + 1]

        rows = []

        max_m = min(
            len(coordinates[j]),
            len(coordinates[j + 1]),
        )

        for m in range(max_m):

            restricted = (
                restriction_to_roots(
                    coordinates[j][m],
                    source_roots,
                    target_roots,
                )
            )

            target = coordinates[j + 1][m]

            delta = [
                sp.factor(
                    target[i]
                    - restricted[i]
                )
                for i in range(
                    len(target)
                )
            ]

            rows.append(delta)

        M = matrix_from_rows(rows)

        print(
            f"j={j}->{j+1}: "
            f"shape={M.shape} "
            f"rank={M.rank()}"
        )

        for m, row in enumerate(rows):

            print(
                f"  m={m}: "
                f"{row}"
            )

        print()


# ============================================================
# 6. NEIGHBORING-j RATIO MATRICES
# ============================================================

def section_6_ratio_rank(
    coordinates
):
    print("=" * 78)
    print(
        "6. NEIGHBORING-j ROOT RATIO MATRICES"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        source_roots = NEGATIVE_CHAIN[j]
        target_roots = NEGATIVE_CHAIN[j + 1]

        rows = []

        max_m = min(
            len(coordinates[j]),
            len(coordinates[j + 1]),
        )

        print(
            f"j={j}->{j+1}"
        )

        for m in range(max_m):

            restricted = (
                restriction_to_roots(
                    coordinates[j][m],
                    source_roots,
                    target_roots,
                )
            )

            target = coordinates[j + 1][m]

            ratio_row = []

            for a, b in zip(
                restricted,
                target,
            ):

                if a == 0:

                    if b == 0:
                        ratio_row.append(
                            "0/0"
                        )
                    else:
                        ratio_row.append(
                            "undefined"
                        )

                else:

                    ratio_row.append(
                        sp.factor(
                            b / a
                        )
                    )

            print(
                f"  m={m}: {ratio_row}"
            )

            numeric_entries = [
                x
                for x in ratio_row
                if isinstance(
                    x,
                    sp.Basic,
                )
            ]

            if numeric_entries:

                rows.append(
                    numeric_entries
                )

        if rows:

            widths = {
                len(row)
                for row in rows
            }

            if len(widths) == 1:

                M = matrix_from_rows(
                    rows
                )

                print(
                    f"  numeric ratio rank="
                    f"{M.rank()}"
                )

        print()


# ============================================================
# 7. ALIGNED COMMON-ROOT MATRICES
# ============================================================

def section_7_aligned_negative_roots(
    coordinates
):
    print("=" * 78)
    print(
        "7. ALIGNED COMMON-ROOT COORDINATE MATRICES"
    )
    print("=" * 78)

    global_roots = [-1, -2, -3, -4]

    for m in range(6):

        rows = []
        active = []

        for j in [0, 1, 2, 3]:

            if m >= len(
                coordinates[j]
            ):
                continue

            row = []

            for root in global_roots:

                if root in ROOTS[j]:

                    idx = ROOTS[j].index(
                        root
                    )

                    row.append(
                        coordinates[j][m][idx]
                    )

                else:

                    row.append(
                        sp.Integer(0)
                    )

            rows.append(row)
            active.append(j)

        if not rows:
            continue

        M = matrix_from_rows(
            rows
        )

        print(
            f"m={m}: "
            f"active_j={active} "
            f"shape={M.shape} "
            f"rank={M.rank()}"
        )

        for j, row in zip(
            active,
            rows
        ):

            print(
                f"  j={j}: {row}"
            )

        print()


# ============================================================
# 8. CROSS-m RANK OF NESTED DIFFERENCES
# ============================================================

def section_8_difference_columns(
    coordinates
):
    print("=" * 78)
    print(
        "8. CROSS-m RANK OF NESTED DIFFERENCE COLUMNS"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        source_roots = NEGATIVE_CHAIN[j]
        target_roots = NEGATIVE_CHAIN[j + 1]

        difference_columns = []

        max_m = min(
            len(coordinates[j]),
            len(coordinates[j + 1]),
        )

        for m in range(max_m):

            restricted = (
                restriction_to_roots(
                    coordinates[j][m],
                    source_roots,
                    target_roots,
                )
            )

            target = coordinates[j + 1][m]

            delta = [
                sp.factor(
                    target[i]
                    - restricted[i]
                )
                for i in range(
                    len(target)
                )
            ]

            difference_columns.append(
                delta
            )

        M = matrix_from_rows(
            difference_columns
        )

        print(
            f"j={j}->{j+1}: "
            f"shape={M.shape} "
            f"row_rank={M.rank()}"
        )

        if M.rows > 0:

            MT = M.T

            print(
                f"  column_rank={MT.rank()}"
            )

        print()


# ============================================================
# 9. NORMALIZED NESTED DIFFERENCES
# ============================================================

def section_9_normalized_differences(
    coordinates
):
    print("=" * 78)
    print(
        "9. NORMALIZED NESTED ROOT DIFFERENCES"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        source_roots = NEGATIVE_CHAIN[j]
        target_roots = NEGATIVE_CHAIN[j + 1]

        print(
            f"j={j}->{j+1}"
        )

        max_m = min(
            len(coordinates[j]),
            len(coordinates[j + 1]),
        )

        for idx, root in enumerate(
            target_roots
        ):

            deltas = []

            for m in range(max_m):

                restricted = (
                    restriction_to_roots(
                        coordinates[j][m],
                        source_roots,
                        target_roots,
                    )
                )

                target = (
                    coordinates[j + 1][m]
                )

                deltas.append(
                    sp.factor(
                        target[idx]
                        - restricted[idx]
                    )
                )

            nonzero = [
                x
                for x in deltas
                if x != 0
            ]

            if not nonzero:

                normalized = deltas

            else:

                base = nonzero[0]

                normalized = [
                    sp.factor(
                        x / base
                    )
                    for x in deltas
                ]

            print(
                f"  root={root}"
            )
            print(
                f"    differences={deltas}"
            )
            print(
                f"    normalized={normalized}"
            )

        print()


# ============================================================
# 10. TERMINAL j=5 SEPARATE AUDIT
# ============================================================

def section_10_terminal(
    coordinates
):
    print("=" * 78)
    print(
        "10. TERMINAL j=5 ROOT-COORDINATE AUDIT"
    )
    print("=" * 78)

    j = 5

    print(
        f"roots={ROOTS[j]}"
    )

    for m, row in enumerate(
        coordinates[j]
    ):

        print(
            f"m={m}:"
        )
        print(
            f"  coordinates={row}"
        )
        print(
            f"  primitive={primitive_vector(row)}"
        )

        nonzero = [
            x for x in row
            if x != 0
        ]

        if nonzero:

            base = nonzero[0]

            normalized = [
                sp.factor(
                    x / base
                )
                for x in row
            ]

            print(
                f"  normalized={normalized}"
            )

    print()


# ============================================================
# 11. EXACT LAGRANGE RECONSTRUCTION
# ============================================================

def lagrange_basis(roots):

    basis = []

    for i, root in enumerate(
        roots
    ):

        numerator = sp.Integer(1)
        denominator = sp.Integer(1)

        for j, other in enumerate(
            roots
        ):

            if i == j:
                continue

            numerator *= (
                K - other
            )

            denominator *= (
                root - other
            )

        basis.append(
            sp.factor(
                numerator / denominator
            )
        )

    return basis


def lagrange_reconstruct(
    coordinates,
    roots
):
    basis = lagrange_basis(
        roots
    )

    result = sp.Integer(0)

    for c, b in zip(
        coordinates,
        basis
    ):
        result += c * b

    return sp.factor(
        sp.expand(result)
    )


def section_11_reconstruction(
    divisions,
    coordinates
):
    print("=" * 78)
    print(
        "11. EXACT ROOT-COORDINATE RECONSTRUCTION"
    )
    print("=" * 78)

    tested = 0
    failures = []

    for j in J_VALUES:

        roots = ROOTS[j]

        if not roots:
            continue

        for m, (_, r) in enumerate(
            divisions[j]
        ):

            recovered = (
                lagrange_reconstruct(
                    coordinates[j][m],
                    roots,
                )
            )

            tested += 1

            if sp.expand(
                recovered - r
            ) != 0:

                failures.append(
                    (
                        j,
                        m,
                        r,
                        recovered,
                    )
                )

    print(
        f"tested={tested}"
    )
    print(
        "reconstruction failures="
        f"{len(failures)}"
    )

    for failure in failures:
        print(
            f"  {failure}"
        )

    print()


# ============================================================
# 12. FRESH-K REMAINDER TARGETS
# ============================================================

def section_12_fresh_k(
    divisions
):
    print("=" * 78)
    print(
        "12. FRESH-K REMAINDER TARGETS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for m, (_, r) in enumerate(
            divisions[j]
        ):

            values = [
                sp.factor(
                    r.subs(K, k0)
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

    tested = 0
    failures = []

    for j in J_VALUES:

        for k0 in K_VALUES:

            for d in D_VALUES:

                y = (
                    d
                    - D0[j]
                )

                value = sp.Integer(0)

                for m in range(
                    Q_DEGREE[j] + 1
                ):

                    value += (
                        coeffs[j][m]
                        .subs(K, k0)
                        * y**m
                    )

                expected = value

                tested += 1

                if sp.expand(
                    value - expected
                ) != 0:

                    failures.append(
                        (
                            j,
                            k0,
                            d,
                        )
                    )

    print(
        f"tested={tested}"
    )
    print(
        "reconstruction failures="
        f"{len(failures)}"
    )
    print()


# ============================================================
# 14. COMPACT SUMMARY
# ============================================================

def section_14_summary(
    coordinates
):
    print("=" * 78)
    print(
        "14. COMPACT NESTED-ROOT SUMMARY"
    )
    print("=" * 78)

    print(
        "j | roots | coordinate_shape | coordinate_rank"
    )
    print("-" * 78)

    for j in J_VALUES:

        roots = ROOTS[j]

        if not roots:

            print(
                f"{j} | {roots} | (0,0) | 0"
            )
            continue

        M = matrix_from_rows(
            coordinates[j]
        )

        print(
            f"{j} | "
            f"{roots} | "
            f"{M.shape} | "
            f"{M.rank()}"
        )

    print()

    print(
        "Nested negative-root chain:"
    )

    for j in [0, 1, 2]:

        print(
            f"  j={j}->{j+1}: "
            f"{NEGATIVE_CHAIN[j]}"
            f" -> "
            f"{NEGATIVE_CHAIN[j+1]}"
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
        "The previous experiment established that the ladder"
    )
    print(
        "remainders are small-degree objects and that their"
    )
    print(
        "values at ladder roots give exact coordinates."
    )
    print()

    print(
        "The earlier cross-j coordinate test was not valid as"
    )
    print(
        "a matrix computation because different j have"
    )
    print(
        "different numbers of ladder roots."
    )
    print()

    print(
        "This experiment fixes that conceptually rather than"
    )
    print(
        "merely padding vectors."
    )
    print()

    print(
        "For j=0,1,2,3 the root sets are genuinely nested:"
    )
    print(
        "  {-1,-2,-3,-4}"
    )
    print(
        "  {-2,-3,-4}"
    )
    print(
        "  {-3,-4}"
    )
    print(
        "  {-4}"
    )
    print()

    print(
        "Therefore the meaningful cross-j objects are:"
    )
    print(
        "  * restriction to common roots;"
    )
    print(
        "  * difference between the two j-level values;"
    )
    print(
        "  * ratio between the two j-level values;"
    )
    print(
        "  * rank across m of those exact differences."
    )
    print()

    print(
        "A genuine structural collapse would appear as:"
    )
    print(
        "  * very low rank in neighboring-j differences;"
    )
    print(
        "  * identical m-patterns across several root positions;"
    )
    print(
        "  * simple neighboring-j ratios;"
    )
    print(
        "  * a terminal j=5 pattern that is visibly separate."
    )
    print()

    print(
        "No arbitrary rational interpolation is used as evidence."
    )
    print(
        "Fresh K values are reported independently."
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


# ============================================================
# MAIN
# ============================================================

def main():

    section_0()

    coeffs = reconstruct_coefficients()

    divisions = divide_by_ladder(
        coeffs
    )

    section_1(
        divisions
    )

    coordinates = root_coordinates(
        divisions
    )

    section_2_coordinates(
        coordinates
    )

    section_3_within_j_rank(
        coordinates
    )

    section_4_nested_restrictions(
        coordinates
    )

    section_5_difference_rank(
        coordinates
    )

    section_6_ratio_rank(
        coordinates
    )

    section_7_aligned_negative_roots(
        coordinates
    )

    section_8_difference_columns(
        coordinates
    )

    section_9_normalized_differences(
        coordinates
    )

    section_10_terminal(
        coordinates
    )

    section_11_reconstruction(
        divisions,
        coordinates
    )

    section_12_fresh_k(
        divisions
    )

    section_13_grid(
        coeffs
    )

    section_14_summary(
        coordinates
    )

    final_diagnostic()


if __name__ == "__main__":
    main()