from __future__ import annotations

import math
from itertools import combinations

import sympy as sp


# =============================================================================
# R=5 NEXT EXPERIMENT
# CROSS-D BOUNDARY FACTORIZATION / FACTOR PERSISTENCE
# =============================================================================
#
# INPUT:
#   The discrepancies E_j(k,D) are taken exactly as supplied.
#
# GOAL:
#
#   The previous experiment discovered exact onset factors:
#
#   j=0:
#       P(K) = -(K+1)(K+2)(K+3)(K+4)(2K+5)/20
#
#   j=1:
#       P(K) = -2(K+2)(K+3)(K+4)(2K+7)/3
#
#   j=2:
#       P(K) = -4(K+3)(K+4)(2K+7)/3
#
#   j=3:
#       P(K) = -5(K+4)(2K+9)
#
#   j=4:
#       P(K) = -10(2K+9)
#
#   j=5:
#       7(K-11)(K-9)(K-7)(K-5)(19859K-68077)/960
#
# This script now asks the stronger question:
#
#   Are those k-factors merely first-boundary phenomena,
#   or do they divide the correction at later D values too?
#
# Specifically:
#
#   1. Test the integer-shift ladder across EVERY D.
#   2. Test the onset linear factor across EVERY D.
#   3. Test the combined onset factor across EVERY D.
#   4. Record exact quotients when divisibility succeeds.
#   5. Measure the first D where each factor disappears.
#   6. Test whether factor persistence is controlled by support depth.
#   7. Test pairwise gcds between fixed-D slices.
#   8. Test whether the factors form a nested boundary filtration.
#
# IMPORTANT:
#   - No r=6.
#   - No full pq-kernel expansion.
#   - No replacement universal law.
#   - No generic degree-5 fit is treated as mathematical evidence.
#
# =============================================================================


# =============================================================================
# SYMBOLS
# =============================================================================

D, K = sp.symbols("D K")


# =============================================================================
# EXACT DATA
# =============================================================================

DATA: dict[int, dict[int, dict[int, int]]] = {
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
            16: 1041408,
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


# =============================================================================
# GRID
# =============================================================================

KS = sorted(DATA)
JS = list(range(6))
DS = sorted(
    {
        d
        for k in DATA
        for j in DATA[k]
        for d in DATA[k][j]
    }
)


# =============================================================================
# BASIC HELPERS
# =============================================================================

def E(k: int, j: int, d: int) -> sp.Integer:
    return sp.Integer(DATA[k][j][d])


def factor(expr: sp.Expr) -> sp.Expr:
    return sp.factor(
        sp.cancel(
            sp.expand(expr)
        )
    )


def degree_in(
    expr: sp.Expr,
    variable: sp.Symbol,
):
    expr = sp.expand(expr)

    if expr == 0:
        return -sp.oo

    return sp.Poly(
        expr,
        variable,
        domain="QQ",
    ).degree()


def exact_roots(
    expr: sp.Expr,
):
    if expr == 0:
        return []

    poly = sp.Poly(
        expr,
        K,
        domain="QQ",
    )

    return sorted(
        list(
            sp.roots(poly.as_expr(), K).keys()
        ),
        key=sp.default_sort_key,
    )


def integer_roots(
    expr: sp.Expr,
    low: int = -50,
    high: int = 50,
):
    return [
        n
        for n in range(low, high + 1)
        if sp.expand(
            expr.subs(K, n)
        ) == 0
    ]


# =============================================================================
# VALIDATION
# =============================================================================

def validate_data():

    expected_js = set(JS)
    expected_ds = set(DS)

    problems = []

    for k in KS:

        actual_js = set(
            DATA[k]
        )

        missing_js = sorted(
            expected_js - actual_js
        )

        extra_js = sorted(
            actual_js - expected_js
        )

        if missing_js:
            problems.append(
                f"k={k}: missing j={missing_js}"
            )

        if extra_js:
            problems.append(
                f"k={k}: extra j={extra_js}"
            )

        for j in JS:

            if j not in DATA[k]:
                continue

            actual_ds = set(
                DATA[k][j]
            )

            missing_ds = sorted(
                expected_ds - actual_ds
            )

            extra_ds = sorted(
                actual_ds - expected_ds
            )

            if missing_ds:
                problems.append(
                    f"k={k}, j={j}: "
                    f"missing D={missing_ds}"
                )

            if extra_ds:
                problems.append(
                    f"k={k}, j={j}: "
                    f"extra D={extra_ds}"
                )

    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)
    print(f"k values = {KS}")
    print(f"j values = {JS}")
    print(f"D values = {DS}")
    print(
        f"points   = "
        f"{len(KS) * len(JS) * len(DS)}"
    )

    if problems:

        print()
        print(
            "GRID PROBLEMS:"
        )

        for p in problems:
            print(
                f"  {p}"
            )

        raise RuntimeError(
            "Incomplete or inconsistent input grid."
        )

    print(
        "grid status = OK"
    )
    print()


# =============================================================================
# POLYNOMIAL INTERPOLATION
# =============================================================================

def interpolate_D(
    k: int,
    j: int,
) -> sp.Expr:

    points = [
        (
            sp.Integer(d),
            E(k,j,d),
        )
        for d in DS
    ]

    return factor(
        sp.interpolate(
            points,
            D,
        )
    )


def interpolate_K(
    j: int,
    d: int,
) -> sp.Expr:

    points = [
        (
            sp.Integer(k),
            E(k,j,d),
        )
        for k in KS
    ]

    return factor(
        sp.interpolate(
            points,
            K,
        )
    )


# =============================================================================
# SUPPORT
# =============================================================================

def common_support_roots(
    j: int,
):
    return [
        d
        for d in DS
        if all(
            E(k,j,d) == 0
            for k in KS
        )
    ]


def support_factor(
    j: int,
):
    result = sp.Integer(1)

    for root in common_support_roots(j):
        result *= D-root

    return factor(result)


def first_nonzero_D(
    j: int,
):
    for d in DS:
        if any(
            E(k,j,d) != 0
            for k in KS
        ):
            return d

    return None


# =============================================================================
# DIVISIBILITY
# =============================================================================

def divide_K(
    polynomial_expr: sp.Expr,
    divisor_expr: sp.Expr,
):
    P = sp.Poly(
        sp.expand(polynomial_expr),
        K,
        domain="QQ",
    )

    Q = sp.Poly(
        sp.expand(divisor_expr),
        K,
        domain="QQ",
    )

    quotient, remainder = sp.div(
        P,
        Q,
    )

    return (
        quotient.as_expr(),
        remainder.as_expr(),
        remainder.is_zero,
    )


def divides_K(
    polynomial_expr: sp.Expr,
    divisor_expr: sp.Expr,
) -> bool:

    _, _, ok = divide_K(
        polynomial_expr,
        divisor_expr,
    )

    return ok


# =============================================================================
# FACTOR TEMPLATES
# =============================================================================

def positive_ladder(
    j: int,
):
    """
    j=0 -> (K+1)(K+2)(K+3)(K+4)
    j=1 -> (K+2)(K+3)(K+4)
    j=2 -> (K+3)(K+4)
    j=3 -> (K+4)
    j=4 -> 1
    """

    result = sp.Integer(1)

    for m in range(
        j+1,
        5,
    ):
        result *= K+m

    return sp.expand(result)


def negative_ladder():
    """
    j=5 terminal ladder:
        (K-5)(K-7)(K-9)(K-11)
    """

    return sp.expand(
        (K-5)
        * (K-7)
        * (K-9)
        * (K-11)
    )


def onset_linear(
    j: int,
    d: int,
):
    """
    Boundary linear factor suggested by the
    j=0,...,4 onset slices:

        2K + D - 1

    """

    return sp.expand(
        2*K + d - 1
    )


def onset_divisor(
    j: int,
    d: int,
):
    """
    For j=0,...,4:

        ladder_j * (2K+D-1)

    For j=5:

        negative terminal ladder

    """

    if j <= 4:

        return factor(
            positive_ladder(j)
            * onset_linear(j,d)
        )

    return factor(
        negative_ladder()
    )


# =============================================================================
# SECTION 1
# =============================================================================

def section_1_factor_persistence():

    print("=" * 78)
    print("1. POSITIVE K-LADDER PERSISTENCE ACROSS ALL D")
    print("=" * 78)

    for j in range(5):

        ladder = positive_ladder(j)

        print(
            f"j={j}"
        )
        print(
            f"  ladder={factor(ladder)}"
        )

        successes = []
        failures = []

        for d in DS:

            P = interpolate_K(j,d)

            ok = divides_K(
                P,
                ladder,
            )

            if ok:
                successes.append(d)
            else:
                failures.append(d)

        print(
            f"  divides at D={successes}"
        )
        print(
            f"  fails at D={failures}"
        )

        print()


# =============================================================================
# SECTION 2
# =============================================================================

def section_2_linear_factor_persistence():

    print("=" * 78)
    print("2. LINEAR FACTOR PERSISTENCE ACROSS ALL D")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for d in DS:

            P = interpolate_K(j,d)
            L = onset_linear(j,d)

            ok = divides_K(
                P,
                L,
            )

            print(
                f"  D={d}: "
                f"factor={L}, "
                f"divides={ok}"
            )

            if ok:

                Q, _, _ = divide_K(
                    P,
                    L,
                )

                print(
                    f"       quotient="
                    f"{factor(Q)}"
                )

        print()


# =============================================================================
# SECTION 3
# =============================================================================

def section_3_combined_onset_factor():

    print("=" * 78)
    print("3. COMBINED ONSET FACTOR PERSISTENCE")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for d in DS:

            P = interpolate_K(j,d)

            divisor = onset_divisor(
                j,
                d,
            )

            ok = divides_K(
                P,
                divisor,
            )

            print(
                f"  D={d}: "
                f"divisor={factor(divisor)}, "
                f"divides={ok}"
            )

            if ok:

                Q, _, _ = divide_K(
                    P,
                    divisor,
                )

                print(
                    f"       quotient="
                    f"{factor(Q)}"
                )

        print()


# =============================================================================
# SECTION 4
# =============================================================================

def section_4_first_failure_map():

    print("=" * 78)
    print("4. FIRST FACTOR-FAILURE MAP")
    print("=" * 78)

    for j in range(5):

        ladder = positive_ladder(j)

        ladder_successes = [
            d
            for d in DS
            if divides_K(
                interpolate_K(j,d),
                ladder,
            )
        ]

        linear_successes = [
            d
            for d in DS
            if divides_K(
                interpolate_K(j,d),
                onset_linear(j,d),
            )
        ]

        combined_successes = [
            d
            for d in DS
            if divides_K(
                interpolate_K(j,d),
                onset_divisor(j,d),
            )
        ]

        print(
            f"j={j}"
        )
        print(
            f"  ladder successes={ladder_successes}"
        )
        print(
            f"  linear successes={linear_successes}"
        )
        print(
            f"  combined successes="
            f"{combined_successes}"
        )

        if ladder_successes:
            print(
                f"  ladder first failure after onset="
                f"{next((d for d in DS if d > min(ladder_successes) and d not in ladder_successes), None)}"
            )

        print()


# =============================================================================
# SECTION 5
# =============================================================================

def section_5_factor_multiplicity_in_k():

    print("=" * 78)
    print("5. EXACT K-ROOT MULTIPLICITY PROFILE")
    print("=" * 78)

    for j in range(5):

        print(
            f"j={j}"
        )

        for d in DS:

            P = interpolate_K(
                j,
                d,
            )

            roots = exact_roots(P)

            if roots:
                print(
                    f"  D={d}: roots={roots}"
                )

        print()


# =============================================================================
# SECTION 6
# =============================================================================

def section_6_integer_root_persistence():

    print("=" * 78)
    print("6. INTEGER K-ROOT PERSISTENCE")
    print("=" * 78)

    candidate_roots = list(
        range(-15,16)
    )

    for j in JS:

        print(
            f"j={j}"
        )

        for r in candidate_roots:

            hit_D = []

            for d in DS:

                P = interpolate_K(
                    j,
                    d,
                )

                if (
                    sp.expand(
                        P.subs(K,r)
                    ) == 0
                ):
                    hit_D.append(d)

            if len(hit_D) >= 2:

                print(
                    f"  K={r}: "
                    f"zero at D={hit_D}"
                )

        print()


# =============================================================================
# SECTION 7
# =============================================================================

def section_7_gcd_of_fixed_D_slices():

    print("=" * 78)
    print("7. GCD OF FIXED-D CROSS-k SLICES")
    print("=" * 78)

    for j in JS:

        print(
            f"j={j}"
        )

        for d1, d2 in combinations(
            DS,
            2,
        ):

            P1 = sp.Poly(
                interpolate_K(j,d1),
                K,
                domain="QQ",
            )

            P2 = sp.Poly(
                interpolate_K(j,d2),
                K,
                domain="QQ",
            )

            g = factor(
                sp.gcd(
                    P1,
                    P2,
                ).as_expr()
            )

            if g != 1:

                print(
                    f"  D=({d1},{d2}): "
                    f"gcd={g}"
                )

        print()


# =============================================================================
# SECTION 8
# =============================================================================

def section_8_pairwise_j_gcd_at_fixed_D():

    print("=" * 78)
    print("8. PAIRWISE j GCD AT FIXED D")
    print("=" * 78)

    for d in DS:

        print(
            f"D={d}"
        )

        for j1, j2 in combinations(
            JS,
            2,
        ):

            P1 = sp.Poly(
                interpolate_K(j1,d),
                K,
                domain="QQ",
            )

            P2 = sp.Poly(
                interpolate_K(j2,d),
                K,
                domain="QQ",
            )

            if P1.is_zero or P2.is_zero:
                continue

            g = factor(
                sp.gcd(
                    P1,
                    P2,
                ).as_expr()
            )

            if g != 1:

                print(
                    f"  (j={j1},j={j2}): "
                    f"gcd={g}"
                )

        print()


# =============================================================================
# SECTION 9
# =============================================================================

def section_9_onset_factor_constants():

    print("=" * 78)
    print("9. ONSET CONSTANTS RECONSTRUCTION")
    print("=" * 78)

    constants = {
        0: sp.Rational(-1,20),
        1: sp.Rational(-2,3),
        2: sp.Rational(-4,3),
        3: sp.Rational(-5,1),
        4: sp.Rational(-10,1),
    }

    for j in range(5):

        d0 = next(
            d
            for d in DS
            if any(
                E(k,j,d) != 0
                for k in KS
            )
        )

        P = interpolate_K(
            j,
            d0,
        )

        expected = factor(
            constants[j]
            * positive_ladder(j)
            * onset_linear(j,d0)
        )

        print(
            f"j={j}"
        )
        print(
            f"  D0={d0}"
        )
        print(
            f"  C={constants[j]}"
        )
        print(
            f"  actual={factor(P)}"
        )
        print(
            f"  expected={expected}"
        )
        print(
            f"  exact match="
            f"{sp.expand(P-expected)==0}"
        )
        print()


# =============================================================================
# SECTION 10
# =============================================================================

def section_10_boundary_depth_test():

    print("=" * 78)
    print("10. BOUNDARY DEPTH TEST")
    print("=" * 78)

    print(
        "For each j, compare:"
    )
    print(
        "  support roots"
    )
    print(
        "  first nonzero D"
    )
    print(
        "  ladder persistence"
    )
    print(
        "  linear persistence"
    )
    print()

    for j in range(5):

        roots = common_support_roots(j)

        d0 = first_nonzero_D(j)

        ladder = positive_ladder(j)

        ladder_D = [
            d
            for d in DS
            if divides_K(
                interpolate_K(j,d),
                ladder,
            )
        ]

        linear_D = [
            d
            for d in DS
            if divides_K(
                interpolate_K(j,d),
                onset_linear(j,d),
            )
        ]

        print(
            f"j={j}"
        )
        print(
            f"  support_roots={roots}"
        )
        print(
            f"  D0={d0}"
        )
        print(
            f"  ladder_D={ladder_D}"
        )
        print(
            f"  linear_D={linear_D}"
        )
        print()


# =============================================================================
# SECTION 11
# =============================================================================

def section_11_boundary_ratio_test():

    print("=" * 78)
    print("11. BOUNDARY QUOTIENT RATIO TEST")
    print("=" * 78)

    for j in range(5):

        ladder = positive_ladder(j)

        print(
            f"j={j}"
        )

        for d in DS:

            P = interpolate_K(j,d)

            if not divides_K(
                P,
                ladder,
            ):
                continue

            Q, _, _ = divide_K(
                P,
                ladder,
            )

            print(
                f"  D={d}: "
                f"Q(K)={factor(Q)}"
            )

        print()


# =============================================================================
# SECTION 12
# =============================================================================

def section_12_terminal_j5_ladder():

    print("=" * 78)
    print("12. TERMINAL j=5 INTEGER LADDER")
    print("=" * 78)

    ladder = negative_ladder()

    print(
        f"terminal ladder="
        f"{factor(ladder)}"
    )
    print()

    for d in DS:

        P = interpolate_K(
            5,
            d,
        )

        ok = divides_K(
            P,
            ladder,
        )

        print(
            f"D={d}: "
            f"divides={ok}"
        )

        if ok:

            Q, _, _ = divide_K(
                P,
                ladder,
            )

            print(
                f"  quotient="
                f"{factor(Q)}"
            )

    print()


# =============================================================================
# SECTION 13
# =============================================================================

def section_13_terminal_j5_root_tracking():

    print("=" * 78)
    print("13. TERMINAL j=5 ROOT TRACKING")
    print("=" * 78)

    for d in DS:

        P = interpolate_K(
            5,
            d,
        )

        roots = integer_roots(
            P,
            -30,
            30,
        )

        print(
            f"D={d}: "
            f"integer roots={roots}"
        )

    print()


# =============================================================================
# SECTION 14
# =============================================================================

def section_14_special_D14_D16_j5():

    print("=" * 78)
    print("14. j=5 D=14 / D=16 COMPARISON")
    print("=" * 78)

    P14 = interpolate_K(5,14)
    P16 = interpolate_K(5,16)

    print(
        "D=14:"
    )
    print(
        f"  P14={factor(P14)}"
    )
    print(
        f"  roots={exact_roots(P14)}"
    )
    print(
        f"  integer roots="
        f"{integer_roots(P14)}"
    )

    print()

    print(
        "D=16:"
    )
    print(
        f"  P16={factor(P16)}"
    )
    print(
        f"  roots={exact_roots(P16)}"
    )
    print(
        f"  integer roots="
        f"{integer_roots(P16)}"
    )

    print()

    g = factor(
        sp.gcd(
            sp.Poly(P14,K,domain="QQ"),
            sp.Poly(P16,K,domain="QQ"),
        ).as_expr()
    )

    print(
        f"gcd(P14,P16)={g}"
    )

    print()


# =============================================================================
# SECTION 15
# =============================================================================

def section_15_fixed_D_degree_table():

    print("=" * 78)
    print("15. FIXED-D CROSS-k DEGREE TABLE")
    print("=" * 78)

    for j in JS:

        print(
            f"j={j}"
        )

        for d in DS:

            P = interpolate_K(
                j,
                d,
            )

            print(
                f"  D={d}: "
                f"degree={degree_in(P,K)}"
            )

        print()


# =============================================================================
# SECTION 16
# =============================================================================

def section_16_support_vs_boundary_factor():

    print("=" * 78)
    print("16. SUPPORT VS BOUNDARY-FACTOR TABLE")
    print("=" * 78)

    print(
        "j | support | D0 | ladder | "
        "linear at D0 | combined at D0"
    )
    print(
        "-" * 78
    )

    for j in JS:

        S = support_factor(j)
        d0 = first_nonzero_D(j)

        if j <= 4:

            ladder = positive_ladder(j)
            linear = onset_linear(j,d0)
            combined = onset_divisor(j,d0)

        else:

            ladder = negative_ladder()
            linear = onset_linear(j,d0)
            combined = ladder

        P = interpolate_K(
            j,
            d0,
        )

        print(
            f"{j} | "
            f"{factor(S)} | "
            f"{d0} | "
            f"{factor(ladder)} | "
            f"{factor(linear)} | "
            f"{factor(combined)}"
        )

        print(
            f"    actual onset="
            f"{factor(P)}"
        )

    print()


# =============================================================================
# SECTION 17
# =============================================================================

def section_17_factor_survival_matrix():

    print("=" * 78)
    print("17. FACTOR SURVIVAL MATRIX")
    print("=" * 78)

    for j in range(5):

        ladder = positive_ladder(j)

        print(
            f"j={j}"
        )

        print(
            "  D : "
            + " ".join(
                f"{d:>4}"
                for d in DS
            )
        )

        print(
            "  L : "
            + " ".join(
                f"{1 if divides_K(interpolate_K(j,d),ladder) else 0:>4}"
                for d in DS
            )
        )

        print(
            "  A : "
            + " ".join(
                f"{1 if divides_K(interpolate_K(j,d),onset_linear(j,d)) else 0:>4}"
                for d in DS
            )
        )

        print()


# =============================================================================
# SECTION 18
# =============================================================================

def section_18_exact_reconstruction():

    print("=" * 78)
    print("18. EXACT INPUT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = 0

    for k in KS:

        for j in JS:

            P = interpolate_D(
                k,
                j,
            )

            for d in DS:

                tested += 1

                predicted = sp.expand(
                    P.subs(D,d)
                )

                actual = E(k,j,d)

                if predicted != actual:

                    failures += 1

                    print(
                        "FAIL:",
                        f"k={k}",
                        f"j={j}",
                        f"D={d}",
                        f"actual={actual}",
                        f"predicted={predicted}",
                    )

    print(
        f"tested={tested}"
    )
    print(
        f"failures={failures}"
    )
    print()


# =============================================================================
# SECTION 19
# =============================================================================

def section_19_support_interpolation_consistency():

    print("=" * 78)
    print("19. SUPPORT-QUOTIENT CONSISTENCY")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in JS:

        S = support_factor(j)

        print(
            f"j={j}: support={S}"
        )

        for k in KS:

            P_D = sp.Poly(
                interpolate_D(k,j),
                D,
                domain="QQ",
            )

            S_D = sp.Poly(
                S,
                D,
                domain="QQ",
            )

            Q, R = sp.div(
                P_D,
                S_D,
            )

            tested += 1

            if not R.is_zero:

                failures += 1

                print(
                    f"  FAIL k={k}: "
                    f"remainder={R.as_expr()}"
                )

    print(
        f"tested={tested}"
    )
    print(
        f"failures={failures}"
    )
    print()


# =============================================================================
# SECTION 20
# =============================================================================

def section_20_final_structural_summary():

    print("=" * 78)
    print("20. FINAL STRUCTURAL SUMMARY")
    print("=" * 78)

    for j in JS:

        roots = common_support_roots(j)
        support = support_factor(j)
        d0 = first_nonzero_D(j)

        if j <= 4:

            onset = interpolate_K(
                j,
                d0,
            )

            ladder = positive_ladder(j)
            linear = onset_linear(j,d0)

            ladder_success = [
                d
                for d in DS
                if divides_K(
                    interpolate_K(j,d),
                    ladder,
                )
            ]

            linear_success = [
                d
                for d in DS
                if divides_K(
                    interpolate_K(j,d),
                    linear if j <= 4 else 1,
                )
            ]

        else:

            onset = interpolate_K(
                j,
                d0,
            )

            ladder = negative_ladder()

            ladder_success = [
                d
                for d in DS
                if divides_K(
                    interpolate_K(j,d),
                    ladder,
                )
            ]

            linear_success = []

        print(
            f"j={j}"
        )
        print(
            f"  support_roots={roots}"
        )
        print(
            f"  support={factor(support)}"
        )
        print(
            f"  D0={d0}"
        )
        print(
            f"  onset={factor(onset)}"
        )

        print(
            f"  ladder="
            f"{factor(ladder)}"
        )

        print(
            f"  ladder_survives_D="
            f"{ladder_success}"
        )

        if j <= 4:

            print(
                f"  linear_survives_D="
                f"{linear_success}"
            )

        print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    validate_data()

    section_1_factor_persistence()
    section_2_linear_factor_persistence()
    section_3_combined_onset_factor()
    section_4_first_failure_map()
    section_5_factor_multiplicity_in_k()
    section_6_integer_root_persistence()
    section_7_gcd_of_fixed_D_slices()
    section_8_pairwise_j_gcd_at_fixed_D()
    section_9_onset_factor_constants()
    section_10_boundary_depth_test()
    section_11_boundary_ratio_test()
    section_12_terminal_j5_ladder()
    section_13_terminal_j5_root_tracking()
    section_14_special_D14_D16_j5()
    section_15_fixed_D_degree_table()
    section_16_support_vs_boundary_factor()
    section_17_factor_survival_matrix()
    section_18_exact_reconstruction()
    section_19_support_interpolation_consistency()
    section_20_final_structural_summary()


if __name__ == "__main__":
    main()

