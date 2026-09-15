#!/usr/bin/env python3
"""
R=5 NEXT EXPERIMENT: EXACT k-COEFFICIENT FUNCTIONS

Goal
----
We already know that the support-removed Q_j(k,D) families are not
low-rank/separable in k and D.

The next question is therefore:

    Q_j(K,D) = sum_m a_{j,m}(K) * (D - D0(j))^m

What exact structure do the coefficient functions a_{j,m}(K) have?

This experiment:

1. validates the supplied coefficient matrices;
2. reconstructs each a_{j,m}(K) from the six K values;
3. factors the exact interpolating polynomial over Q;
4. finds exact rational roots;
5. computes ordinary degree in K;
6. computes finite-difference degree on n=(K-3)/2;
7. rewrites in the Newton/binomial basis;
8. searches for integer ladders K+c and 2K+c;
9. specifically audits j=5 and K=5,7,9,11;
10. checks exact reconstruction against all supplied coefficients.

No r=6.
No full pq-kernel expansion.
No replacement universal law.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, List, Tuple

import sympy as sp


# ============================================================================
# GLOBAL GRID
# ============================================================================

K_VALUES = [3, 5, 7, 9, 11, 13]
N_VALUES = list(range(len(K_VALUES)))  # n = 0,...,5 where K = 3 + 2n

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

# Number of powers of Y retained in Q_j:
# j=0 -> degree 5
# j=1 -> degree 4
# j=2 -> degree 4
# j=3 -> degree 3
# j=4 -> degree 3
# j=5 -> degree 1
Q_DEGREE = {
    0: 5,
    1: 4,
    2: 4,
    3: 3,
    4: 3,
    5: 1,
}


# ============================================================================
# EXACT Y-COEFFICIENT MATRICES
#
# Rows are K = 3,5,7,9,11,13.
# Columns are powers Y^0, Y^1, ..., Y^d.
#
# These are copied from the successful Experiment 54 output.
# ============================================================================

Y_ROWS: Dict[int, List[List[sp.Rational]]] = {

    # ------------------------------------------------------------------------
    # j = 0
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # j = 1
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # j = 2
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # j = 3
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # j = 4
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # j = 5
    # ------------------------------------------------------------------------
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
# SUPPORT POLYNOMIALS
# ============================================================================

SUPPORT = {
    0: sp.Integer(1),
    1: sp.Symbol("D") - 6,
    2: sp.Symbol("D") - 6,
    3: (sp.Symbol("D") - 8) * (sp.Symbol("D") - 6),
    4: (sp.Symbol("D") - 8) * (sp.Symbol("D") - 6),
    5: (sp.Symbol("D") - 12)
       * (sp.Symbol("D") - 10)
       * (sp.Symbol("D") - 8)
       * (sp.Symbol("D") - 6),
}


# ============================================================================
# BASIC HELPERS
# ============================================================================

K = sp.Symbol("K")
Y = sp.Symbol("Y")


def factor_exact(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(expr))


def rational_roots(expr: sp.Expr) -> List[sp.Expr]:
    expr = sp.Poly(expr, K, domain=sp.QQ)
    if expr.is_zero:
        return []
    return sp.solve(expr.as_expr(), K)


def integer_roots(expr: sp.Expr) -> List[int]:
    roots = []
    p = sp.Poly(expr, K, domain=sp.QQ)
    if p.is_zero:
        return roots

    for r in sp.solve(p.as_expr(), K):
        if r.is_Rational and r.q == 1:
            roots.append(int(r))
    return roots


def degree_or_minus_inf(expr: sp.Expr) -> int:
    p = sp.Poly(expr, K, domain=sp.QQ)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def interpolate_k(values: List[sp.Expr]) -> sp.Expr:
    """
    Interpolate in K using the six points:
        K = 3,5,7,9,11,13
    """
    points = list(zip(K_VALUES, values))
    return sp.factor(sp.interpolate(points, K))


def finite_difference_sequence(values: List[sp.Expr]) -> List[List[sp.Expr]]:
    """
    Forward differences on n=0,...,5.
    """
    levels = [list(map(sp.sympify, values))]

    while len(levels[-1]) > 1:
        old = levels[-1]
        new = [sp.simplify(old[i + 1] - old[i]) for i in range(len(old) - 1)]
        levels.append(new)

    return levels


def finite_difference_degree(values: List[sp.Expr]) -> int:
    """
    Exact finite-difference degree in n=(K-3)/2.

    If a function is polynomial of degree d on six equally spaced n values,
    its d-th differences are constant/nonzero and its (d+1)-st difference is 0.
    """
    levels = finite_difference_sequence(values)

    for d in range(len(levels) - 1):
        if len(levels[d]) == 1:
            continue

        if all(x == levels[d][0] for x in levels[d]):
            # Constant sequence.
            return d

    return len(values) - 1


def newton_binomial_coefficients(values: List[sp.Expr]) -> List[sp.Expr]:
    """
    Newton coefficients for n=0,1,...,5:

        f(n) = c0*C(n,0) + c1*C(n,1) + ... + c5*C(n,5)

    where c_r = Delta^r f(0).
    """
    levels = finite_difference_sequence(values)
    coeffs = []

    for r, level in enumerate(levels):
        coeffs.append(sp.simplify(level[0]))

    return coeffs


def newton_degree(values: List[sp.Expr]) -> int:
    coeffs = newton_binomial_coefficients(values)

    degree = -sp.oo
    for r, c in enumerate(coeffs):
        if c != 0:
            degree = r

    return degree


def primitive_integer_polynomial(expr: sp.Expr) -> Tuple[sp.Poly, sp.Rational]:
    """
    Return primitive integer polynomial p and rational scale s such that

        expr = s * p(K).

    Robust for:
      - zero;
      - constants;
      - one nonzero coefficient;
      - arbitrary rational coefficients.
    """
    p = sp.Poly(expr, K, domain=sp.QQ)

    if p.is_zero:
        return sp.Poly(0, K, domain=sp.ZZ), sp.Rational(0)

    coeffs = list(p.all_coeffs())

    denoms = [sp.denom(c) for c in coeffs]
    den_lcm = sp.ilcm(*denoms) if len(denoms) >= 2 else int(denoms[0])

    ints = [int(c * den_lcm) for c in coeffs]

    nonzero_ints = [abs(v) for v in ints if v != 0]
    if not nonzero_ints:
        return sp.Poly(0, K, domain=sp.ZZ), sp.Rational(0)

    content = sp.igcd(*nonzero_ints) if len(nonzero_ints) >= 2 else nonzero_ints[0]

    ints = [v // content for v in ints]

    # Canonical positive leading coefficient.
    if ints[0] < 0:
        ints = [-v for v in ints]

    primitive = sp.Poly.from_list(ints, gens=K, domain=sp.ZZ)
    scale = sp.simplify(sp.LC(p) / primitive.LC())

    return primitive, scale


def candidate_linear_factors(expr: sp.Expr) -> List[Tuple[sp.Expr, bool]]:
    """
    Test common integer-shift and small even-linear factors.
    """
    candidates = []

    for c in range(-15, 16):
        f = K + c
        q, r = sp.div(sp.Poly(expr, K, domain=sp.QQ),
                      sp.Poly(f, K, domain=sp.QQ))
        candidates.append((f, r.is_zero))

    for c in range(-15, 16, 2):
        f = 2 * K + c
        q, r = sp.div(sp.Poly(expr, K, domain=sp.QQ),
                      sp.Poly(f, K, domain=sp.QQ))
        candidates.append((f, r.is_zero))

    return candidates


def exact_divides(expr: sp.Expr, divisor: sp.Expr) -> bool:
    p = sp.Poly(expr, K, domain=sp.QQ)
    q = sp.Poly(divisor, K, domain=sp.QQ)

    if p.is_zero:
        return True

    _, rem = sp.div(p, q)
    return rem.is_zero


# ============================================================================
# SECTION 0
# ============================================================================

def validate_data() -> None:
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    assert K_VALUES == [3, 5, 7, 9, 11, 13]
    assert J_VALUES == [0, 1, 2, 3, 4, 5]
    assert D_VALUES == [6, 8, 10, 12, 14, 16]

    expected_rows = 6

    for j in J_VALUES:
        rows = Y_ROWS[j]
        assert len(rows) == expected_rows
        assert len(rows[0]) == Q_DEGREE[j] + 1

        for row in rows:
            assert len(row) == Q_DEGREE[j] + 1

    print(f"k values = {K_VALUES}")
    print(f"j values = {J_VALUES}")
    print(f"D values = {D_VALUES}")
    print(f"points   = {len(K_VALUES) * len(J_VALUES) * len(D_VALUES)}")
    print("grid status = OK")
    print(
        f"checked coefficient rows = "
        f"{sum(len(Y_ROWS[j]) for j in J_VALUES)}"
    )
    print()


# ============================================================================
# SECTION 1
# ============================================================================

def section_1_coefficient_polynomials() -> Dict[int, List[sp.Expr]]:
    print("=" * 78)
    print("1. EXACT k-COEFFICIENT FUNCTIONS a_{j,m}(K)")
    print("=" * 78)

    coeff_polys: Dict[int, List[sp.Expr]] = {}

    for j in J_VALUES:
        rows = Y_ROWS[j]
        d = Q_DEGREE[j]
        polys = []

        print(f"j={j}")
        print(f"  D0={D0[j]}")
        print(f"  Q degree in Y={d}")

        for m in range(d + 1):
            values = [rows[i][m] for i in range(len(K_VALUES))]
            poly = interpolate_k(values)
            polys.append(poly)

            print(f"  a_{j},{m}(K):")
            print(f"    degree={degree_or_minus_inf(poly)}")
            print(f"    polynomial={sp.expand(poly)}")
            print(f"    factorization={factor_exact(poly)}")

        coeff_polys[j] = polys
        print()

    return coeff_polys


# ============================================================================
# SECTION 2
# ============================================================================

def section_2_root_audit(coeff_polys: Dict[int, List[sp.Expr]]) -> None:
    print("=" * 78)
    print("2. EXACT ROOT AUDIT OF EACH k-COEFFICIENT")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for m, poly in enumerate(coeff_polys[j]):
            print(f"  a_{j},{m}(K):")

            if sp.expand(poly) == 0:
                print("    ZERO POLYNOMIAL")
                continue

            print(f"    rational_roots={rational_roots(poly)}")
            print(f"    integer_roots={integer_roots(poly)}")

        print()


# ============================================================================
# SECTION 3
# ============================================================================

def section_3_linear_factor_search(coeff_polys: Dict[int, List[sp.Expr]]) -> None:
    print("=" * 78)
    print("3. EXACT LINEAR FACTOR SEARCH IN k")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for m, poly in enumerate(coeff_polys[j]):
            print(f"  a_{j},{m}(K):")

            if sp.expand(poly) == 0:
                print("    ZERO POLYNOMIAL")
                continue

            found = []

            for factor, ok in candidate_linear_factors(poly):
                if ok:
                    quotient = sp.factor(sp.cancel(poly / factor))
                    found.append((factor, quotient))

            if not found:
                print("    no tested linear factors")
            else:
                for factor, quotient in found:
                    print(f"    factor={factor}")
                    print(f"      quotient={quotient}")

        print()


# ============================================================================
# SECTION 4
# ============================================================================

def section_4_finite_difference_audit(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("4. FINITE-DIFFERENCE / NEWTON DEGREE OF a_{j,m}(K)")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        rows = Y_ROWS[j]

        for m in range(Q_DEGREE[j] + 1):
            values = [rows[i][m] for i in range(6)]

            ordinary = degree_or_minus_inf(coeff_polys[j][m])
            fd_degree = finite_difference_degree(values)
            newton_deg = newton_degree(values)
            newton_coeffs = newton_binomial_coefficients(values)

            print(f"  a_{j},{m}(K):")
            print(f"    ordinary_degree={ordinary}")
            print(f"    finite_difference_degree={fd_degree}")
            print(f"    newton_degree={newton_deg}")
            print(f"    newton_coeffs={newton_coeffs}")

        print()


# ============================================================================
# SECTION 5
# ============================================================================

def section_5_primitive_forms(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("5. PRIMITIVE INTEGER FORMS OF a_{j,m}(K)")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for m, poly in enumerate(coeff_polys[j]):
            if sp.expand(poly) == 0:
                print(f"  a_{j},{m}(K): ZERO")
                continue

            primitive, scale = primitive_integer_polynomial(poly)

            print(f"  a_{j},{m}(K)")
            print(f"    primitive={primitive.as_expr()}")
            print(f"    scale={scale}")

        print()


# ============================================================================
# SECTION 6
# ============================================================================

def section_6_special_j5(coeff_polys: Dict[int, List[sp.Expr]]) -> None:
    print("=" * 78)
    print("6. SPECIAL j=5 COEFFICIENT STRUCTURE")
    print("=" * 78)

    a0 = coeff_polys[5][0]
    a1 = coeff_polys[5][1]

    terminal_factor = (
        (K - 5)
        * (K - 7)
        * (K - 9)
        * (K - 11)
    )

    print("j=5")
    print(f"  Q_5(K,D) = a_0(K) + (D-14)*a_1(K)")
    print()

    print("  a_0(K):")
    print(f"    polynomial={sp.expand(a0)}")
    print(f"    factorization={sp.factor(a0)}")
    print(f"    integer_roots={integer_roots(a0)}")
    print()

    print("  a_1(K):")
    print(f"    polynomial={sp.expand(a1)}")
    print(f"    factorization={sp.factor(a1)}")
    print(f"    integer_roots={integer_roots(a1)}")
    print()

    print("  terminal ladder:")
    print(f"    L(K)={terminal_factor}")
    print(f"    a_0 divisible by L={exact_divides(a0, terminal_factor)}")
    print(f"    a_1 divisible by L={exact_divides(a1, terminal_factor)}")

    if exact_divides(a0, terminal_factor):
        print(
            f"    a_0 / L = "
            f"{sp.factor(sp.cancel(a0 / terminal_factor))}"
        )

    if exact_divides(a1, terminal_factor):
        print(
            f"    a_1 / L = "
            f"{sp.factor(sp.cancel(a1 / terminal_factor))}"
        )

    print()
    print("  values at terminal roots:")

    for k0 in [5, 7, 9, 11]:
        v0 = sp.simplify(a0.subs(K, k0))
        v1 = sp.simplify(a1.subs(K, k0))

        print(
            f"    K={k0}: "
            f"a_0={v0}, "
            f"a_1={v1}"
        )

    print()


# ============================================================================
# SECTION 7
# ============================================================================

def section_7_middle_j5_proportionality(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("7. j=5 MIDDLE-K PROPORTIONALITY AUDIT")
    print("=" * 78)

    a0 = coeff_polys[5][0]
    a1 = coeff_polys[5][1]

    print(
        "For K=5,7,9,11 the Y^0 coefficient vanishes, so only a_1(K) "
        "remains."
    )
    print()

    middle = [5, 7, 9, 11]
    base_k = middle[0]

    base = sp.simplify(a1.subs(K, base_k))

    print(f"base K={base_k}")
    print(f"  a_1={base}")
    print()

    for k0 in middle[1:]:
        current = sp.simplify(a1.subs(K, k0))

        if base == 0:
            print(f"K={k0}: base coefficient is zero; ratio undefined")
        else:
            ratio = sp.factor(sp.cancel(current / base))
            print(
                f"K={k0}: "
                f"a_1={current}, "
                f"ratio_to_K{base_k}={ratio}"
            )

    print()


# ============================================================================
# SECTION 8
# ============================================================================

def section_8_coefficient_reconstruction(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("8. EXACT COEFFICIENT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = []

    for j in J_VALUES:
        rows = Y_ROWS[j]

        for i, k0 in enumerate(K_VALUES):
            for m in range(Q_DEGREE[j] + 1):
                predicted = sp.simplify(
                    coeff_polys[j][m].subs(K, k0)
                )
                actual = sp.simplify(rows[i][m])

                tested += 1

                if predicted != actual:
                    failures.append(
                        (j, k0, m, actual, predicted)
                    )

    print(f"tested={tested}")
    print(f"coefficient reconstruction failures={len(failures)}")

    if failures:
        for item in failures[:20]:
            print(f"  failure={item}")

    print()


# ============================================================================
# SECTION 9
# ============================================================================

def section_9_cross_coefficient_gcds(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("9. GCD / FACTOR OVERLAP BETWEEN DIFFERENT Y-COEFFICIENTS")
    print("=" * 78)

    for j in J_VALUES:
        polys = coeff_polys[j]

        print(f"j={j}")

        for m1, m2 in itertools.combinations(range(len(polys)), 2):
            p1 = sp.Poly(polys[m1], K, domain=sp.QQ)
            p2 = sp.Poly(polys[m2], K, domain=sp.QQ)

            if p1.is_zero or p2.is_zero:
                continue

            g = sp.factor(sp.gcd(p1, p2).as_expr())

            if g != 1:
                print(
                    f"  gcd(a_{j},{m1}, a_{j},{m2}) = {g}"
                )

        print()


# ============================================================================
# SECTION 10
# ============================================================================

def section_10_shifted_k_search(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("10. SHIFTED-K FACTOR SEARCH")
    print("=" * 78)

    shifts = list(range(-12, 13))

    for j in J_VALUES:
        print(f"j={j}")

        for m, poly in enumerate(coeff_polys[j]):
            if sp.expand(poly) == 0:
                print(f"  a_{j},{m}: ZERO")
                continue

            hits = []

            for s in shifts:
                f = K + s
                if exact_divides(poly, f):
                    hits.append(f)

            if hits:
                print(f"  a_{j},{m}: K-shift factors={hits}")
            else:
                print(f"  a_{j},{m}: no integer-shift roots")

        print()


# ============================================================================
# SECTION 11
# ============================================================================

def section_11_structural_summary(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("11. COMPACT COEFFICIENT-STRUCTURE SUMMARY")
    print("=" * 78)

    print(
        "j | Y-degree | coefficient | deg_K | fd_degree | factorization"
    )
    print("-" * 110)

    for j in J_VALUES:
        for m, poly in enumerate(coeff_polys[j]):
            if sp.expand(poly) == 0:
                print(
                    f"{j:1d} | "
                    f"{Q_DEGREE[j]:7d} | "
                    f"a_{j},{m:<1d} | "
                    f"ZERO"
                )
                continue

            deg_k = degree_or_minus_inf(poly)
            fd = finite_difference_degree(
                [Y_ROWS[j][i][m] for i in range(6)]
            )
            fac = factor_exact(poly)

            print(
                f"{j:1d} | "
                f"{Q_DEGREE[j]:7d} | "
                f"a_{j},{m:<1d} | "
                f"{str(deg_k):>5} | "
                f"{fd:>9} | "
                f"{fac}"
            )

    print()


# ============================================================================
# SECTION 12
# ============================================================================

def section_12_fresh_evaluation_targets(
    coeff_polys: Dict[int, List[sp.Expr]]
) -> None:
    print("=" * 78)
    print("12. FRESH-K EVALUATION TARGETS")
    print("=" * 78)

    fresh_values = [-7, -5, -3, -1, 1, 15, 17, 19]

    print(
        "These are NOT used for fitting. They are only reported as "
        "out-of-sample symbolic targets."
    )
    print()

    for j in J_VALUES:
        print(f"j={j}")

        for m, poly in enumerate(coeff_polys[j]):
            if sp.expand(poly) == 0:
                print(f"  a_{j},{m}: ZERO")
                continue

            print(f"  a_{j},{m}:")
            for k0 in fresh_values:
                value = sp.factor(
                    sp.cancel(poly.subs(K, k0))
                )
                print(f"    K={k0}: {value}")

        print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    validate_data()

    coeff_polys = section_1_coefficient_polynomials()

    section_2_root_audit(coeff_polys)
    section_3_linear_factor_search(coeff_polys)
    section_4_finite_difference_audit(coeff_polys)
    section_5_primitive_forms(coeff_polys)
    section_6_special_j5(coeff_polys)
    section_7_middle_j5_proportionality(coeff_polys)
    section_8_coefficient_reconstruction(coeff_polys)
    section_9_cross_coefficient_gcds(coeff_polys)
    section_10_shifted_k_search(coeff_polys)
    section_11_structural_summary(coeff_polys)
    section_12_fresh_evaluation_targets(coeff_polys)

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print(
        "The rank experiment showed that Q_j(k,D) is not low-rank "
        "separable in k and D."
    )
    print()
    print(
        "This experiment therefore resolves Q_j into exact coefficient "
        "functions:"
    )
    print()
    print(
        "    Q_j(K,D) = sum_m a_{j,m}(K) (D-D0(j))^m"
    )
    print()
    print(
        "Each a_{j,m}(K) is reconstructed exactly from the six supplied "
        "K-values."
    )
    print()
    print(
        "The next structural test is factorization and finite-difference "
        "degree of the individual coefficient functions."
    )
    print()
    print(
        "j=5 is separately audited for the terminal ladder "
        "(K-5)(K-7)(K-9)(K-11)."
    )
    print()
    print("No r=6.")
    print("No full pq-kernel expansion.")
    print("No replacement universal r,j law.")
    print("=" * 78)


if __name__ == "__main__":
    main()

