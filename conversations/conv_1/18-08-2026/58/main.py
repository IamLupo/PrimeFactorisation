#!/usr/bin/env python3

"""
R=5 NEXT EXPERIMENT
============================================================
NEWTON / FALLING-FACTORIAL BASIS STRUCTURE IN K

Goal
----

Previous experiments established:

  * exact D-support factors;
  * clean boundary K-ladders in a_{j,0}(K);
  * higher a_{j,m}(K), m>0, are generically degree 5;
  * simple K -> K +/- j shifts do not reduce cross-j rank;
  * removing the leading K^5 contribution does not produce
    a low-rank collapse.

The next structural question is:

    Are the degree-5 coefficient functions sparse when written
    in a Newton / falling-factorial basis adapted to the boundary
    roots?

For a chosen integer anchor s, use

    B_r(K;s) = prod_{t=0}^{r-1} (K-s-t)

with B_0 = 1.

Equivalently,

    B_r(K;s) = r! * binomial(K-s, r)

for integer r >= 0.

Because the clean a_{j,0} factors are consecutive linear factors,
this basis is more structurally meaningful than ordinary powers.

For every a_{j,m} we compute:

  1. exact Newton coefficients;
  2. number of nonzero Newton coefficients;
  3. highest occupied Newton index;
  4. zero-pattern across r;
  5. comparison of Newton coefficient rows across j;
  6. exact finite differences in the Newton basis;
  7. boundary-adapted anchors;
  8. reconstruction;
  9. fresh-K symbolic targets.

No universal law is inferred from the existing six K-values.
"""

from __future__ import annotations

import sympy as sp

K = sp.Symbol("K")


# ============================================================
# GRID
# ============================================================

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
# BOUNDARY INFORMATION
# ============================================================

LADDERS = {
    0: (K + 1) * (K + 2) * (K + 3) * (K + 4),
    1: (K + 2) * (K + 3) * (K + 4),
    2: (K + 3) * (K + 4),
    3: (K + 4),
    4: sp.Integer(1),
    5: (K - 5) * (K - 7) * (K - 9) * (K - 11),
}


# Candidate anchors for Newton basis.
ANCHORS = {
    "K": lambda j: 0,
    "K-j": lambda j: j,
    "K+j": lambda j: -j,
    "K+(4-j)": lambda j: j - 4,
    "K-(4-j)": lambda j: 4 - j,
    "boundary_left": {
        0: -1,
        1: -2,
        2: -3,
        3: -4,
        4: -5,
        5: 5,
    },
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

        for m in range(Q_DEGREE[j] + 1):
            values = [
                Y_ROWS[j][i][m]
                for i in range(len(K_VALUES))
            ]

            coeffs[j].append(
                interpolate(values)
            )

    return coeffs


def basis_term(anchor, r):
    """
    B_r(K; anchor) =
        product_{t=0}^{r-1}(K-anchor-t)

    B_0 = 1.
    """
    if r == 0:
        return sp.Integer(1)

    out = sp.Integer(1)

    for t in range(r):
        out *= (
            K - anchor - t
        )

    return sp.expand(out)


def newton_coefficients(expr, anchor):
    """
    Expand a polynomial exactly in the basis

        B_r(K;anchor)
          = product_{t=0}^{r-1}(K-anchor-t).

    Coefficients are obtained by repeated forward differences
    evaluated at K=anchor:

        c_r = Delta^r f(anchor) / r!
    """
    p = poly(expr)

    if p.is_zero:
        return [sp.Integer(0)]

    d = p.degree()

    values = [
        sp.simplify(
            p.as_expr().subs(
                K,
                anchor + n,
            )
        )
        for n in range(d + 1)
    ]

    coeffs = []

    current = values

    for r in range(d + 1):

        coeffs.append(
            sp.simplify(
                current[0]
                / sp.factorial(r)
            )
        )

        if len(current) == 1:
            break

        current = [
            sp.simplify(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

    return coeffs


def reconstruct_from_newton(coeffs, anchor):
    out = sp.Integer(0)

    for r, c in enumerate(coeffs):
        out += c * basis_term(
            anchor,
            r,
        )

    return sp.factor(
        sp.expand(out)
    )


def primitive_vector(values):
    """
    Primitive integer normalization of a rational vector.
    """
    values = [
        sp.Rational(v)
        for v in values
        if v != 0
    ]

    if not values:
        return [0]

    denoms = [
        int(v.q)
        for v in values
    ]

    if len(denoms) == 1:
        lcm_den = denoms[0]
    else:
        lcm_den = int(
            sp.ilcm(
                *denoms
            )
        )

    ints = [
        int(v * lcm_den)
        for v in values
    ]

    abs_nonzero = [
        abs(x)
        for x in ints
        if x != 0
    ]

    if len(abs_nonzero) == 1:
        g = abs_nonzero[0]
    else:
        g = int(
            sp.igcd(
                *abs_nonzero
            )
        )

    ints = [
        x // g
        for x in ints
    ]

    # Restore omitted zero entries is handled by caller.
    return ints


def primitive_full_vector(values):
    """
    Preserve zeros while normalizing the entire vector.
    """
    vals = [
        sp.Rational(v)
        for v in values
    ]

    nonzero = [
        v
        for v in vals
        if v != 0
    ]

    if not nonzero:
        return [0] * len(vals)

    denoms = [
        int(v.q)
        for v in nonzero
    ]

    if len(denoms) == 1:
        lcm_den = denoms[0]
    else:
        lcm_den = int(
            sp.ilcm(
                *denoms
            )
        )

    ints = [
        int(v * lcm_den)
        for v in vals
    ]

    nz = [
        abs(x)
        for x in ints
        if x != 0
    ]

    if len(nz) == 1:
        g = nz[0]
    else:
        g = int(
            sp.igcd(
                *nz
            )
        )

    ints = [
        x // g
        for x in ints
    ]

    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-y for y in ints]
            break

    return ints


def coefficient_matrix(rows):
    if not rows:
        return sp.Matrix([])

    width = max(
        len(row)
        for row in rows
    )

    padded = []

    for row in rows:
        row = list(row)

        if len(row) < width:
            row += [
                sp.Integer(0)
            ] * (
                width - len(row)
            )

        padded.append(row)

    return sp.Matrix(padded)


# ============================================================
# SECTION 0
# ============================================================

def section_0():
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    checked = 0

    for j in J_VALUES:

        width = (
            Q_DEGREE[j] + 1
        )

        if len(Y_ROWS[j]) != 6:
            raise RuntimeError(
                f"bad row count for j={j}"
            )

        for row in Y_ROWS[j]:

            if len(row) != width:
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
# SECTION 1
# ============================================================

def section_1_newton_profiles(coeffs):
    print("=" * 78)
    print(
        "1. NEWTON COEFFICIENT PROFILES "
        "IN ORDINARY K BASIS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for m in range(
            Q_DEGREE[j] + 1
        ):

            expr = coeffs[j][m]

            c = newton_coefficients(
                expr,
                0,
            )

            nz = [
                r
                for r, x in enumerate(c)
                if x != 0
            ]

            highest = (
                max(nz)
                if nz
                else -sp.oo
            )

            print(
                f"  a_{{{j},{m}}}:"
            )
            print(
                f"    coeffs={c}"
            )
            print(
                f"    nonzero_indices={nz}"
            )
            print(
                f"    support_size={len(nz)}"
            )
            print(
                f"    highest_index={highest}"
            )

        print()


# ============================================================
# SECTION 2
# ============================================================

def section_2_anchor_comparison(coeffs):
    print("=" * 78)
    print(
        "2. NEWTON SPARSITY ACROSS K-ANCHORS"
    )
    print("=" * 78)

    for name, rule in ANCHORS.items():

        if isinstance(rule, dict):
            continue

        print(
            f"ANCHOR = {name}"
        )

        for j in J_VALUES:

            anchor = rule(j)

            total = 0
            max_support = 0
            supports = []

            for m in range(
                Q_DEGREE[j] + 1
            ):

                c = newton_coefficients(
                    coeffs[j][m],
                    anchor,
                )

                nz = sum(
                    x != 0
                    for x in c
                )

                total += nz
                max_support = max(
                    max_support,
                    nz,
                )

                supports.append(nz)

            print(
                f"  j={j}: "
                f"anchor={anchor} "
                f"support_sizes={supports} "
                f"total_nonzero={total}"
            )

        print()


# ============================================================
# SECTION 3
# ============================================================

def section_3_boundary_anchor(coeffs):
    print("=" * 78)
    print(
        "3. BOUNDARY-ADAPTED NEWTON BASIS"
    )
    print("=" * 78)

    rule = ANCHORS["boundary_left"]

    for j in J_VALUES:

        anchor = rule[j]

        print(
            f"j={j}"
        )
        print(
            f"  anchor={anchor}"
        )

        for m in range(
            Q_DEGREE[j] + 1
        ):

            c = newton_coefficients(
                coeffs[j][m],
                anchor,
            )

            print(
                f"  m={m}:"
            )
            print(
                f"    coeffs={c}"
            )
            print(
                f"    primitive="
                f"{primitive_full_vector(c)}"
            )

        print()


# ============================================================
# SECTION 4
# ============================================================

def section_4_zero_pattern(coeffs):
    print("=" * 78)
    print(
        "4. EXACT NEWTON ZERO-PATTERN MATRIX"
    )
    print("=" * 78)

    anchor_rule = ANCHORS[
        "boundary_left"
    ]

    for j in J_VALUES:

        anchor = anchor_rule[j]

        print(
            f"j={j}, anchor={anchor}"
        )

        for m in range(
            Q_DEGREE[j] + 1
        ):

            c = newton_coefficients(
                coeffs[j][m],
                anchor,
            )

            pattern = [
                1 if x != 0 else 0
                for x in c
            ]

            print(
                f"  m={m}: {pattern}"
            )

        print()


# ============================================================
# SECTION 5
# ============================================================

def section_5_cross_j_newton_rank(coeffs):
    print("=" * 78)
    print(
        "5. CROSS-j RANK OF NEWTON COEFFICIENT MATRICES"
    )
    print("=" * 78)

    anchor_rule = ANCHORS[
        "boundary_left"
    ]

    for m in range(6):

        rows = []
        active = []

        for j in J_VALUES:

            if m > Q_DEGREE[j]:
                continue

            anchor = anchor_rule[j]

            c = newton_coefficients(
                coeffs[j][m],
                anchor,
            )

            rows.append(c)
            active.append(j)

        if not rows:
            continue

        M = coefficient_matrix(
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
# SECTION 6
# ============================================================

def section_6_newton_coefficient_ratios(coeffs):
    print("=" * 78)
    print(
        "6. SUCCESSIVE NEWTON-COEFFICIENT RATIOS"
    )
    print("=" * 78)

    anchor_rule = ANCHORS[
        "boundary_left"
    ]

    for j in J_VALUES:

        anchor = anchor_rule[j]

        print(
            f"j={j}, anchor={anchor}"
        )

        for m in range(
            Q_DEGREE[j] + 1
        ):

            c = newton_coefficients(
                coeffs[j][m],
                anchor,
            )

            ratios = []

            for r in range(
                len(c) - 1
            ):

                if c[r] == 0:
                    ratios.append(
                        "undefined"
                    )
                else:
                    ratios.append(
                        sp.factor(
                            c[r + 1]
                            / c[r]
                        )
                    )

            print(
                f"  m={m}:"
            )
            print(
                f"    ratios={ratios}"
            )

        print()


# ============================================================
# SECTION 7
# ============================================================

def section_7_cross_j_newton_gcd(coeffs):
    print("=" * 78)
    print(
        "7. CROSS-j GCD OF NEWTON-BASIS COEFFICIENT POLYNOMIALS"
    )
    print("=" * 78)

    anchor_rule = ANCHORS[
        "boundary_left"
    ]

    max_r = 5

    for r in range(
        max_r + 1
    ):

        print(
            f"Newton index r={r}"
        )

        active = []

        exprs = []

        for j in J_VALUES:

            anchor = anchor_rule[j]

            for_m = Q_DEGREE[j]

            for m in range(
                for_m + 1
            ):

                c = newton_coefficients(
                    coeffs[j][m],
                    anchor,
                )

                if r >= len(c):
                    continue

                value = c[r]

                if value == 0:
                    continue

                exprs.append(
                    sp.Integer(value)
                )
                active.append(
                    (j, m)
                )

        if not exprs:
            print(
                "  no active coefficients"
            )
            print()
            continue

        print(
            f"  active_count={len(exprs)}"
        )

        # Constants have no nontrivial polynomial gcd.
        # We instead inspect rational ratios against the first
        # active coefficient.
        base = exprs[0]

        ratios = [
            sp.factor(
                e / base
            )
            for e in exprs[1:]
        ]

        print(
            f"  first={active[0]} value={base}"
        )

        for pair, ratio in zip(
            active[1:],
            ratios,
        ):

            if ratio != 1:
                print(
                    f"    {pair}: ratio={ratio}"
                )

        print()


# ============================================================
# SECTION 8
# ============================================================

def section_8_newton_reconstruction(coeffs):
    print("=" * 78)
    print(
        "8. EXACT NEWTON RECONSTRUCTION"
    )
    print("=" * 78)

    anchor_rule = ANCHORS[
        "boundary_left"
    ]

    tested = 0
    failures = []

    for j in J_VALUES:

        anchor = anchor_rule[j]

        for m in range(
            Q_DEGREE[j] + 1
        ):

            expr = coeffs[j][m]

            c = newton_coefficients(
                expr,
                anchor,
            )

            recovered = (
                reconstruct_from_newton(
                    c,
                    anchor,
                )
            )

            if sp.expand(
                recovered - expr
            ) != 0:

                failures.append(
                    (
                        j,
                        m,
                        expr,
                        recovered,
                    )
                )

            tested += 1

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


# ============================================================
# SECTION 9
# ============================================================

def section_9_fresh_k_newton(coeffs):
    print("=" * 78)
    print(
        "9. FRESH-K NEWTON-BASIS TARGETS"
    )
    print("=" * 78)

    anchor_rule = ANCHORS[
        "boundary_left"
    ]

    for j in J_VALUES:

        anchor = anchor_rule[j]

        print(
            f"j={j}, anchor={anchor}"
        )

        for m in range(
            Q_DEGREE[j] + 1
        ):

            c = newton_coefficients(
                coeffs[j][m],
                anchor,
            )

            print(
                f"  m={m}"
            )

            for k0 in FRESH_K:

                value = sp.factor(
                    sum(
                        c[r]
                        * sp.prod(
                            K - anchor - t
                            for t in range(r)
                        )
                        for r in range(
                            len(c)
                        )
                    ).subs(
                        K,
                        k0,
                    )
                )

                print(
                    f"    K={k0}: {value}"
                )

        print()


# ============================================================
# SECTION 10
# ============================================================

def section_10_original_grid_reconstruction(
    coeffs
):
    print("=" * 78)
    print(
        "10. EXACT ORIGINAL-GRID RECONSTRUCTION"
    )
    print("=" * 78)

    tested = 0
    failures = []

    for j in J_VALUES:

        for i, k0 in enumerate(
            K_VALUES
        ):

            for m in range(
                Q_DEGREE[j] + 1
            ):

                expected = sp.Rational(
                    Y_ROWS[j][i][m]
                )

                predicted = sp.factor(
                    coeffs[j][m].subs(
                        K,
                        k0,
                    )
                )

                tested += 1

                if predicted != expected:
                    failures.append(
                        (
                            j,
                            k0,
                            m,
                            expected,
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

    for item in failures[:20]:
        print(
            f"  {item}"
        )

    print()


# ============================================================
# SECTION 11
# ============================================================

def section_11_compact_summary(coeffs):
    print("=" * 78)
    print(
        "11. COMPACT NEWTON STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    anchor_rule = ANCHORS[
        "boundary_left"
    ]

    print(
        "j | m | anchor | "
        "K-degree | Newton support | "
        "highest r"
    )
    print("-" * 78)

    for j in J_VALUES:

        anchor = anchor_rule[j]

        for m in range(
            Q_DEGREE[j] + 1
        ):

            expr = coeffs[j][m]

            c = newton_coefficients(
                expr,
                anchor,
            )

            nz = [
                r
                for r, value in enumerate(c)
                if value != 0
            ]

            highest = (
                max(nz)
                if nz
                else -sp.oo
            )

            print(
                f"{j} | {m} | {anchor} | "
                f"{degree(expr)} | "
                f"{len(nz)} | "
                f"{highest}"
            )

    print()


# ============================================================
# FINAL
# ============================================================

def final_diagnostic():
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The previous experiments ruled out:"
    )
    print(
        "  * simple k-D separability;"
    )
    print(
        "  * simple j-dependent K shifts;"
    )
    print(
        "  * leading-term subtraction as a low-rank collapse;"
    )
    print(
        "  * shared higher-coefficient polynomial factors."
    )
    print()

    print(
        "The remaining clean structure is concentrated in"
    )
    print(
        "the consecutive K-factors of a_{j,0}."
    )
    print()

    print(
        "This experiment therefore changes from ordinary powers"
    )
    print(
        "of K to an exact Newton/falling-factorial basis."
    )
    print()

    print(
        "The important question is NOT whether interpolation can"
    )
    print(
        "produce another degree-5 representation."
    )
    print(
        "The important question is whether the Newton basis"
    )
    print(
        "produces exact sparsity, short support, repeated coefficient"
    )
    print(
        "patterns, or low-rank structure across j and m."
    )
    print()

    print(
        "A genuine discovery would look like:"
    )
    print(
        "  * exact zeros in large portions of the Newton arrays;"
    )
    print(
        "  * a systematic dependence of occupied r on j or m;"
    )
    print(
        "  * exact repeated Newton coefficients;"
    )
    print(
        "  * a low-rank coefficient matrix surviving exact"
    )
    print(
        "    reconstruction;"
    )
    print(
        "  * confirmation at fresh K values."
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

    section_1_newton_profiles(
        coeffs
    )

    section_2_anchor_comparison(
        coeffs
    )

    section_3_boundary_anchor(
        coeffs
    )

    section_4_zero_pattern(
        coeffs
    )

    section_5_cross_j_newton_rank(
        coeffs
    )

    section_6_newton_coefficient_ratios(
        coeffs
    )

    section_7_cross_j_newton_gcd(
        coeffs
    )

    section_8_newton_reconstruction(
        coeffs
    )

    section_9_fresh_k_newton(
        coeffs
    )

    section_10_original_grid_reconstruction(
        coeffs
    )

    section_11_compact_summary(
        coeffs
    )

    final_diagnostic()


if __name__ == "__main__":
    main()
