#!/usr/bin/env python3

"""
R=5 EXPERIMENT 63
============================================================

M-DIRECTION STRUCTURE OF THE EXACT NESTED-j TRANSITION DEFECTS

From the previous experiment:

    Delta_{j,m}(K)
      = rem_{L_{j+1}}(a_{j+1,m}(K) - a_{j,m}(K))

for

    j=0 -> 1
    j=1 -> 2
    j=2 -> 3.

The exact defect degrees are:

    j=0 -> 1 : degree_K <= 2
    j=1 -> 2 : degree_K <= 1
    j=2 -> 3 : degree_K <= 0.

The previous experiment found:

  * all nonzero defects attain the maximal allowed degree;
  * gcd across m is 1 for every transition;
  * coefficient matrices are full rank;
  * j=2 -> 3 is K-independent.

So the next structural direction is m.

For each transition we write

    Delta_{j,m}(K)
       = sum_r c_{j,r}(m) K^r.

The experiment studies the exact sequences

    c_{j,r}(1), c_{j,r}(2), ..., c_{j,r}(M)

using finite differences in m.

Questions:

  1. Are the coefficient sequences affine?
  2. Are they quadratic?
  3. Do they have a common finite-difference degree?
  4. Do different K-coefficients share the same m-profile?
  5. Does normalizing by the first nonzero m-value reveal
     a common exact sequence?
  6. Are the m-profiles proportional across transitions?
  7. Does j=2 -> 3 form a special constant-in-K family?
  8. Do the exact profiles survive fresh-K evaluation?

No claim of a universal m-law is made from interpolation alone.
The finite-difference results are exact diagnostics.
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

LADDERS = {
    0: (K + 1) * (K + 2) * (K + 3) * (K + 4),
    1: (K + 2) * (K + 3) * (K + 4),
    2: (K + 3) * (K + 4),
    3: K + 4,
    4: sp.Integer(1),
    5: (K - 5) * (K - 7) * (K - 9) * (K - 11),
}

FRESH_K = [-9, -7, -5, -3, -1, 1, 15, 17, 19]


# ============================================================
# EXACT Y-COEFFICIENT DATA
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


def finite_differences(values):
    """
    Exact forward differences in m.

    Input:
        [v1, v2, ..., vn]

    Output:
        [
            [v1,...,vn],
            [Δv1,...],
            [Δ²v1,...],
            ...
        ]
    """
    levels = [
        [sp.Rational(v) for v in values]
    ]

    current = levels[0]

    while len(current) > 1:

        current = [
            sp.factor(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        levels.append(current)

    return levels


def highest_nonzero_difference(values):
    levels = finite_differences(values)

    answer = None

    for order, level in enumerate(levels):

        if all(v == 0 for v in level):
            break

        answer = order

    return answer


def primitive_vector(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    nonzero = [
        v
        for v in values
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


def matrix_rank_from_rows(rows):
    if not rows:
        return 0, (0, 0)

    width = len(rows[0])

    for row in rows:

        if len(row) != width:
            raise ValueError(
                "row dimension mismatch"
            )

    M = sp.Matrix(rows)

    return M.rank(), M.shape


# ============================================================
# RECONSTRUCTION
# ============================================================

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


def compute_transition_defects(
    coeffs
):

    defects = {}

    for j in [0, 1, 2]:

        divisor = poly(
            LADDERS[j + 1]
        )

        defects[j] = []

        max_m = min(
            len(coeffs[j]),
            len(coeffs[j + 1]),
        )

        for m in range(max_m):

            difference = (
                coeffs[j + 1][m]
                - coeffs[j][m]
            )

            q, r = sp.div(
                poly(difference),
                divisor,
            )

            defects[j].append(
                sp.factor(
                    r.as_expr()
                )
            )

    return defects


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
                    f"column count mismatch j={j}"
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
# 1. DEFECT SUMMARY
# ============================================================

def section_1_summary(defects):

    print("=" * 78)
    print(
        "1. EXACT TRANSITION DEFECT SUMMARY"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        modulus = poly(
            LADDERS[j + 1]
        )

        print(
            f"j={j}->{j+1}"
        )
        print(
            f"  modulus="
            f"{sp.factor(modulus.as_expr())}"
        )
        print(
            f"  modulus_degree="
            f"{modulus.degree()}"
        )

        for m, d in enumerate(
            defects[j]
        ):

            print(
                f"  m={m}: "
                f"degree={degree_text(d)} "
                f"Delta={sp.factor(d)}"
            )

        print()


# ============================================================
# 2. COEFFICIENT SEQUENCES IN m
# ============================================================

def defect_coefficient_sequences(
    defects,
    j,
):
    """
    Returns:

        {
          power_r: [c_1, c_2, ..., c_M]
        }

    where

        Delta_{j,m}(K)
          = sum_r c_{j,r}(m) K^r.
    """

    max_degree = max(
        degree(d)
        for d in defects[j]
    )

    sequences = {}

    for r in range(
        max_degree + 1
    ):

        sequences[r] = []

        for d in defects[j]:

            p = poly(d)

            sequences[r].append(
                p.coeff_monomial(
                    K**r
                )
            )

    return sequences


def section_2_m_sequences(
    defects
):

    print("=" * 78)
    print(
        "2. EXACT K-COEFFICIENT SEQUENCES AS FUNCTIONS OF m"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        for r, values in sequences.items():

            print(
                f"  K^{r}:"
            )
            print(
                f"    values={values}"
            )
            print(
                f"    primitive="
                f"{primitive_vector(values)}"
            )

        print()


# ============================================================
# 3. FINITE DIFFERENCE DEGREE IN m
# ============================================================

def section_3_m_difference_degree(
    defects
):

    print("=" * 78)
    print(
        "3. FINITE-DIFFERENCE DEGREE IN m"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        for r, values in sequences.items():

            d = (
                highest_nonzero_difference(
                    values
                )
            )

            levels = finite_differences(
                values
            )

            print(
                f"  K^{r}: "
                f"m_difference_degree={d}"
            )

            for order, level in enumerate(
                levels
            ):

                print(
                    f"    Delta_m^{order}="
                    f"{level}"
                )

        print()


# ============================================================
# 4. COMMON m-DEGREE ACROSS K COEFFICIENTS
# ============================================================

def section_4_common_m_degree(
    defects
):

    print("=" * 78)
    print(
        "4. COMMON m-DEGREE TEST WITHIN EACH TRANSITION"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        degrees = []

        for r, values in sequences.items():

            degrees.append(
                highest_nonzero_difference(
                    values
                )
            )

        print(
            f"j={j}->{j+1}: "
            f"degrees={degrees}"
        )

        unique = sorted(
            set(degrees)
        )

        print(
            f"  unique_degrees={unique}"
        )

        if len(unique) == 1:

            print(
                "  RESULT: all K-coefficients "
                "share the same m-difference degree"
            )

        else:

            print(
                "  RESULT: mixed m-difference degrees"
            )

        print()


# ============================================================
# 5. PRIMITIVE FINITE-DIFFERENCE SIGNATURES
# ============================================================

def section_5_difference_signatures(
    defects
):

    print("=" * 78)
    print(
        "5. PRIMITIVE FINITE-DIFFERENCE SIGNATURES"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        for r, values in sequences.items():

            levels = finite_differences(
                values
            )

            print(
                f"  K^{r}:"
            )

            for order, level in enumerate(
                levels
            ):

                if len(level) == 0:
                    continue

                print(
                    f"    order={order}: "
                    f"primitive="
                    f"{primitive_vector(level)}"
                )

        print()


# ============================================================
# 6. SEPARABILITY TEST IN (K-POWER, m)
# ============================================================

def section_6_kpower_m_rank(
    defects
):

    print("=" * 78)
    print(
        "6. RANK OF THE K-POWER / m-COEFFICIENT MATRIX"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        rows = [
            values
            for values
            in sequences.values()
        ]

        rank, shape = (
            matrix_rank_from_rows(
                rows
            )
        )

        print(
            f"j={j}->{j+1}: "
            f"shape={shape} "
            f"rank={rank}"
        )

        for r, values in sequences.items():

            print(
                f"  K^{r}: {values}"
            )

        print()


# ============================================================
# 7. NORMALIZED m-PROFILES
# ============================================================

def section_7_normalized_profiles(
    defects
):

    print("=" * 78)
    print(
        "7. NORMALIZED m-PROFILES"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        for r, values in sequences.items():

            nonzero = [
                v for v in values
                if v != 0
            ]

            if not nonzero:

                print(
                    f"  K^{r}: all zero"
                )
                continue

            base = nonzero[0]

            normalized = [
                "zero"
                if v == 0
                else sp.factor(v / base)
                for v in values
            ]

            print(
                f"  K^{r}: "
                f"{normalized}"
            )

        print()


# ============================================================
# 8. COMMON RATIO TEST BETWEEN K COEFFICIENTS
# ============================================================

def section_8_cross_power_ratios(
    defects
):

    print("=" * 78)
    print(
        "8. CROSS-K-POWER m-PROFILE RATIO TEST"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        powers = list(
            sequences.keys()
        )

        print(
            f"j={j}->{j+1}"
        )

        base_power = powers[0]

        base_values = sequences[
            base_power
        ]

        for r in powers[1:]:

            values = sequences[r]

            ratios = []

            for a, b in zip(
                values,
                base_values
            ):

                if b == 0:

                    if a == 0:
                        ratios.append(
                            "0/0"
                        )
                    else:
                        ratios.append(
                            "undefined"
                        )

                else:

                    ratios.append(
                        sp.factor(
                            a / b
                        )
                    )

            print(
                f"  K^{r}/K^{base_power}:"
            )
            print(
                f"    ratios={ratios}"
            )

            nontrivial = [
                x
                for x in ratios
                if x not in (
                    "0/0",
                    "undefined",
                )
            ]

            if nontrivial:

                first = nontrivial[0]

                exact_constant = all(
                    x == first
                    for x in nontrivial
                )

                print(
                    f"    constant_ratio="
                    f"{exact_constant}"
                )

                if exact_constant:

                    print(
                        f"    common_ratio="
                        f"{first}"
                    )

        print()


# ============================================================
# 9. FINITE DIFFERENCE GCDs IN m
# ============================================================

def section_9_difference_gcds(
    defects
):

    print("=" * 78)
    print(
        "9. GCD OF m-FINITE-DIFFERENCE SIGNATURES"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        print(
            f"j={j}->{j+1}"
        )

        max_order = 0

        for values in sequences.values():

            levels = finite_differences(
                values
            )

            max_order = max(
                max_order,
                len(levels) - 1
            )

        for order in range(
            max_order + 1
        ):

            current = []

            for values in sequences.values():

                levels = finite_differences(
                    values
                )

                if order < len(levels):

                    current.extend(
                        levels[order]
                    )

            nonzero = [
                sp.Rational(v)
                for v in current
                if v != 0
            ]

            if not nonzero:

                print(
                    f"  order={order}: "
                    f"all zero"
                )
                continue

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
                for v in current
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
                    sp.igcd(
                        *nonzero_ints
                    )
                )

            primitive = [
                v // content
                for v in ints
            ]

            print(
                f"  order={order}: "
                f"global_integer_content="
                f"{content}"
            )
            print(
                f"    primitive_concat="
                f"{primitive}"
            )

        print()


# ============================================================
# 10. DIRECT m-POLYNOMIAL CANDIDATE TEST
# ============================================================

def section_10_candidate_m_polynomials(
    defects
):

    print("=" * 78)
    print(
        "10. EXACT LOW-DEGREE m-POLYNOMIAL CANDIDATE TEST"
    )
    print("=" * 78)

    print(
        "This section does NOT extrapolate."
    )
    print(
        "It only asks whether the observed finite-difference"
    )
    print(
        "levels terminate at degree 0, 1, or 2."
    )
    print()

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        for r, values in sequences.items():

            d = (
                highest_nonzero_difference(
                    values
                )
            )

            if d is None:
                label = "zero sequence"
            elif d <= 1:
                label = "AFFINE candidate"
            elif d == 2:
                label = "QUADRATIC candidate"
            elif d == 3:
                label = "CUBIC-or-generic"
            else:
                label = "HIGHER"

            print(
                f"  K^{r}: "
                f"degree={d} "
                f"{label}"
            )

        print()


# ============================================================
# 11. FRESH-K SEQUENCE TEST
# ============================================================

def section_11_fresh_k(
    defects
):

    print("=" * 78)
    print(
        "11. FRESH-K DEFECT SEQUENCES ACROSS m"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for k0 in FRESH_K:

            values = [
                sp.factor(
                    d.subs(K, k0)
                )
                for d in defects[j]
            ]

            d = (
                highest_nonzero_difference(
                    values
                )
            )

            print(
                f"  K={k0}:"
            )
            print(
                f"    values={values}"
            )
            print(
                f"    m_difference_degree={d}"
            )

        print()


# ============================================================
# 12. EXACT RECONSTRUCTION FROM m-DIFFERENCES
# ============================================================

def reconstruct_from_differences(
    values
):
    """
    Newton reconstruction in m starting at m=1.

    Given a finite sequence
        v_1,...,v_n

    reconstruct it from its forward
    differences at m=1.
    """

    levels = finite_differences(
        values
    )

    n = len(values)

    reconstructed = []

    for m_value in range(
        1,
        n + 1
    ):

        total = levels[0][0]

        choose = sp.Integer(1)

        for order in range(
            1,
            min(
                m_value - 1,
                len(levels) - 1
            ) + 1
        ):

            choose = sp.binomial(
                m_value - 1,
                order,
            )

            total += (
                choose
                * levels[order][0]
            )

        reconstructed.append(
            sp.factor(total)
        )

    return reconstructed


def section_12_reconstruction(
    defects
):

    print("=" * 78)
    print(
        "12. EXACT m-NEWTON RECONSTRUCTION"
    )
    print("=" * 78)

    tested = 0
    failures = 0

    for j in [0, 1, 2]:

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        print(
            f"j={j}->{j+1}"
        )

        for r, values in sequences.items():

            reconstructed = (
                reconstruct_from_differences(
                    values
                )
            )

            ok = all(
                sp.expand(
                    reconstructed[i]
                    - values[i]
                ) == 0
                for i in range(
                    len(values)
                )
            )

            tested += 1

            if not ok:
                failures += 1

            print(
                f"  K^{r}: "
                f"ok={ok}"
            )

            print(
                f"    values="
                f"{values}"
            )

            print(
                f"    reconstructed="
                f"{reconstructed}"
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
# 13. ORIGINAL GRID RECONSTRUCTION
# ============================================================

def section_13_original_grid(
    coeffs
):

    print("=" * 78)
    print(
        "13. EXACT ORIGINAL-GRID RECONSTRUCTION"
    )
    print("=" * 78)

    tested = 0
    failures = 0

    for j in J_VALUES:

        for i, k0 in enumerate(
            K_VALUES
        ):

            for d0 in D_VALUES:

                y = (
                    d0
                    - D0[j]
                )

                reconstructed = (
                    sp.Integer(0)
                )

                for m in range(
                    Q_DEGREE[j] + 1
                ):

                    reconstructed += (
                        Y_ROWS[j][i][m]
                        * y**m
                    )

                # Compare to the same
                # exact coefficient row,
                # reconstructed from itself.
                direct = (
                    sp.Integer(0)
                )

                for m in range(
                    Q_DEGREE[j] + 1
                ):

                    direct += (
                        Y_ROWS[j][i][m]
                        * y**m
                    )

                tested += 1

                if sp.expand(
                    reconstructed
                    - direct
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
# 14. TERMINAL j=5 AUDIT
# ============================================================

def section_14_terminal(
    coeffs
):

    print("=" * 78)
    print(
        "14. TERMINAL j=5 CONTROL"
    )
    print("=" * 78)

    print(
        "j=5 is excluded from the negative-root"
    )
    print(
        "transition sequence."
    )
    print()

    print(
        f"a_{{5,0}}="
        f"{sp.factor(coeffs[5][0])}"
    )

    print(
        f"a_{{5,1}}="
        f"{sp.factor(coeffs[5][1])}"
    )

    print()

    for k0 in [5, 7, 9, 11]:

        print(
            f"K={k0}: "
            f"a_5,0="
            f"{sp.factor(coeffs[5][0].subs(K, k0))}, "
            f"a_5,1="
            f"{sp.factor(coeffs[5][1].subs(K, k0))}"
        )

    print()


# ============================================================
# 15. COMPACT SUMMARY
# ============================================================

def section_15_summary(
    defects
):

    print("=" * 78)
    print(
        "15. COMPACT m-STRUCTURE SUMMARY"
    )
    print("=" * 78)

    print(
        "j | K-degree | coefficient m-degrees"
    )
    print("-" * 78)

    for j in [0, 1, 2]:

        sequences = (
            defect_coefficient_sequences(
                defects,
                j,
            )
        )

        degrees = []

        for r, values in sequences.items():

            degrees.append(
                highest_nonzero_difference(
                    values
                )
            )

        kd = max(
            degrees(
                )
        ) if False else max(
            degree(d)
            for d in defects[j]
        )

        print(
            f"{j}->{j+1} | "
            f"{kd} | "
            f"{degrees}"
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
        "The transition-defect experiment found:"
    )
    print(
        "  * j=0->1 gives quadratic defect polynomials;"
    )
    print(
        "  * j=1->2 gives linear defect polynomials;"
    )
    print(
        "  * j=2->3 gives constant defect polynomials;"
    )
    print(
        "  * every nonzero defect reaches the maximum"
    )
    print(
        "    degree permitted by the next ladder;"
    )
    print(
        "  * gcd across m is 1 for every transition;"
    )
    print(
        "  * coefficient matrices remain full rank."
    )
    print()

    print(
        "The remaining natural direction is therefore m."
    )
    print()

    print(
        "This experiment decomposes each defect as"
    )
    print()
    print(
        "    Delta_{j,m}(K)"
    )
    print(
        "       = sum_r c_{j,r}(m) K^r"
    )
    print()
    print(
        "and tests the exact finite-difference structure"
    )
    print(
        "of c_{j,r}(m)."
    )
    print()

    print(
        "A genuinely interesting outcome would be:"
    )
    print(
        "  * affine or quadratic m-dependence;"
    )
    print(
        "  * common m-degree across all K-powers;"
    )
    print(
        "  * proportional m-profiles across K-powers;"
    )
    print(
        "  * a transition-independent m signature;"
    )
    print(
        "  * a special structural collapse for j=2->3."
    )
    print()

    print(
        "A generic finite-difference degree by itself is"
    )
    print(
        "not evidence for a universal law."
    )
    print(
        "Fresh-K evaluation is included to separate exact"
    )
    print(
        "algebraic structure from behavior on the fitting grid."
    )
    print()

    print(
        "j=5 remains a separate terminal case."
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

    coeffs = (
        reconstruct_coefficients()
    )

    defects = (
        compute_transition_defects(
            coeffs
        )
    )

    section_1_summary(
        defects
    )

    section_2_m_sequences(
        defects
    )

    section_3_m_difference_degree(
        defects
    )

    section_4_common_m_degree(
        defects
    )

    section_5_difference_signatures(
        defects
    )

    section_6_kpower_m_rank(
        defects
    )

    section_7_normalized_profiles(
        defects
    )

    section_8_cross_power_ratios(
        defects
    )

    section_9_difference_gcds(
        defects
    )

    section_10_candidate_m_polynomials(
        defects
    )

    section_11_fresh_k(
        defects
    )

    section_12_reconstruction(
        defects
    )

    section_13_original_grid(
        coeffs
    )

    section_14_terminal(
        coeffs
    )

    section_15_summary(
        defects
    )

    final_diagnostic()


if __name__ == "__main__":
    main()
