from __future__ import annotations

import sympy as sp
from itertools import combinations


# =============================================================================
# R=5 NEXT EXPERIMENT:
# BOUNDARY-FACTOR DEPTH + FACTOR TRANSITION MAP
# =============================================================================
#
# This script uses the supplied discrepancies directly.
#
# It does NOT:
#   - subtract an old law;
#   - analyze r=6;
#   - expand a full pq-kernel;
#   - assert a replacement universal formula;
#   - treat a generic degree-5 interpolation as a discovered law.
#
# The previous experiment found:
#
# j=0:
#   D0=6
#   onset = -(K+1)(K+2)(K+3)(K+4)(2K+5)/20
#
# j=1:
#   D0=8
#   onset = -2(K+2)(K+3)(K+4)(2K+7)/3
#
# j=2:
#   D0=8
#   onset = -4(K+3)(K+4)(2K+7)/3
#
# j=3:
#   D0=10
#   onset = -5(K+4)(2K+9)
#
# j=4:
#   D0=10
#   onset = -10(2K+9)
#
# j=5:
#   D0=14
#   onset = 7(K-11)(K-9)(K-7)(K-5)(19859K-68077)/960
#
# The new question is:
#
#   HOW DOES THE EXACT K-FACTOR STRUCTURE CHANGE AS D MOVES
#   AWAY FROM THE FIRST-NONZERO BOUNDARY?
#
# We therefore compute:
#
#   1. exact factorization for every fixed-D slice;
#   2. exact multiplicity of every candidate linear factor;
#   3. factor-survival depth;
#   4. first D at which each onset factor disappears;
#   5. whether disappearance occurs simultaneously;
#   6. whether support depth predicts factor depth;
#   7. exact gcds of neighboring D slices;
#   8. exact transition ratios after removing common factors;
#   9. a boundary transition table;
#  10. structural tests against simple "factor survives up to D0" rules;
#  11. exact reconstruction checks.
#
# Special handling:
#
#   A zero polynomial has every factor, algebraically speaking, but that is
#   NOT counted as evidence of factor persistence. We report zero slices
#   separately and only count factor divisibility for nonzero slices.
#
# =============================================================================


# =============================================================================
# SYMBOLS
# =============================================================================

D, K = sp.symbols("D K")


# =============================================================================
# EXACT INPUT DATA
# =============================================================================

DATA: dict[int, dict[int, dict[int, int]]] = {
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
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: 56056},
    },

    9: {
        0: {6: -19734, 8: -158444, 10: -665236, 12: -1913769, 14: -4553505, 16: -10125050},
        1: {6: 0, 8: -28600, 10: -171600, 12: -570570, 14: -1298440, 16: -3360148},
        2: {6: 0, 8: -5200, 10: -34710, 12: -130221, 14: -189488, 16: -177100},
        3: {6: 0, 8: 0, 10: -1755, 12: -12831, 14: 22750, 16: 533910},
        4: {6: 0, 8: 0, 10: -270, 12: -1974, 14: 3360, 16: 372780},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: 132192},
    },

    11: {
        0: {6: -44226, 8: -352716, 10: -1476384, 12: -4138134, 14: -9002994, 16: -20452146},
        1: {6: 0, 8: -52780, 10: -316225, 12: -1028027, 14: -1682590, 16: -3115294},
        2: {6: 0, 8: -8120, 10: -54075, 12: -200417, 14: -12390, 16: 2135770},
        3: {6: 0, 8: 0, 10: -2325, 12: -16995, 14: 96258, 16: 2239830},
        4: {6: 0, 8: 0, 10: -310, 12: -2266, 14: 10672, 16: 1041408},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: 277134},
    },

    13: {
        0: {6: -88536, 8: -702576, 10: -2934064, 12: -7980140, 14: -13755924, 16: -27091064},
        1: {6: 0, 8: -89760, 10: -537200, 12: -1702856, 14: -891072, 16: 7361136},
        2: {6: 0, 8: -11968, 10: -79560, 12: -291108, 14: 609960, 16: 11645748},
        3: {6: 0, 8: 0, 10: -2975, 12: -21743, 14: 243984, 16: 6970578},
        4: {6: 0, 8: 0, 10: -350, 12: -2558, 14: 22491, 16: 2485476},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 532252, 16: 533876},
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
        for k in KS
        for j in JS
        for d in DATA[k][j]
    }
)


# =============================================================================
# BASIC HELPERS
# =============================================================================

def E(k: int, j: int, d: int) -> sp.Integer:
    return sp.Integer(DATA[k][j][d])


def factor(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.expand(expr)))


def degree_in(expr: sp.Expr, variable: sp.Symbol):
    expr = sp.expand(expr)

    if expr == 0:
        return -sp.oo

    return sp.Poly(
        expr,
        variable,
        domain="QQ",
    ).degree()


def is_zero_poly(expr: sp.Expr) -> bool:
    return sp.expand(expr) == 0


def interpolate_D(k: int, j: int) -> sp.Expr:
    pts = [
        (sp.Integer(d), E(k, j, d))
        for d in DS
    ]

    return factor(
        sp.interpolate(
            pts,
            D,
        )
    )


def interpolate_K(j: int, d: int) -> sp.Expr:
    pts = [
        (sp.Integer(k), E(k, j, d))
        for k in KS
    ]

    return factor(
        sp.interpolate(
            pts,
            K,
        )
    )


def divide_K(
    expr: sp.Expr,
    divisor: sp.Expr,
):
    P = sp.Poly(
        sp.expand(expr),
        K,
        domain="QQ",
    )

    Q = sp.Poly(
        sp.expand(divisor),
        K,
        domain="QQ",
    )

    quotient, remainder = sp.div(P, Q)

    return (
        factor(quotient.as_expr()),
        factor(remainder.as_expr()),
        remainder.is_zero,
    )


def nonzero_divides(
    expr: sp.Expr,
    divisor: sp.Expr,
) -> bool:
    """
    A zero polynomial is NOT counted as evidence that a factor persists.
    """
    if is_zero_poly(expr):
        return False

    _, _, ok = divide_K(expr, divisor)
    return ok


def exact_roots(expr: sp.Expr):
    if is_zero_poly(expr):
        return []

    p = sp.Poly(
        sp.expand(expr),
        K,
        domain="QQ",
    )

    return sorted(
        list(
            sp.roots(
                p.as_expr(),
                K,
            ).keys()
        ),
        key=sp.default_sort_key,
    )


def integer_roots(
    expr: sp.Expr,
    low: int = -30,
    high: int = 30,
):
    if is_zero_poly(expr):
        return []

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

        if set(DATA[k]) != expected_js:
            problems.append(
                f"k={k}: "
                f"j-grid mismatch: "
                f"{sorted(DATA[k])}"
            )

        for j in JS:

            actual = set(
                DATA[k][j]
            )

            if actual != expected_ds:
                problems.append(
                    f"k={k}, j={j}: "
                    f"bad D-grid: "
                    f"{sorted(actual)}"
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
        print("GRID PROBLEMS:")
        for p in problems:
            print(f"  {p}")

        raise RuntimeError(
            "Input grid is incomplete or inconsistent."
        )

    print("grid status = OK")
    print()


# =============================================================================
# SUPPORT
# =============================================================================

def common_support_roots(j: int):
    return [
        d
        for d in DS
        if all(
            E(k, j, d) == 0
            for k in KS
        )
    ]


def support_factor(j: int):
    S = sp.Integer(1)

    for root in common_support_roots(j):
        S *= D-root

    return factor(S)


def first_nonzero_D(j: int):
    for d in DS:
        if any(
            E(k, j, d) != 0
            for k in KS
        ):
            return d

    return None


# =============================================================================
# ONSET FACTORS
# =============================================================================

def positive_ladder(j: int):
    """
    j=0: (K+1)(K+2)(K+3)(K+4)
    j=1: (K+2)(K+3)(K+4)
    j=2: (K+3)(K+4)
    j=3: (K+4)
    j=4: 1
    """
    result = sp.Integer(1)

    for m in range(j + 1, 5):
        result *= K + m

    return sp.expand(result)


def terminal_ladder_j5():
    return sp.expand(
        (K - 5)
        * (K - 7)
        * (K - 9)
        * (K - 11)
    )


def onset_linear(d: int):
    return sp.expand(
        2*K + d - 1
    )


def onset_constant(j: int):
    return {
        0: sp.Rational(-1, 20),
        1: sp.Rational(-2, 3),
        2: sp.Rational(-4, 3),
        3: sp.Rational(-5, 1),
        4: sp.Rational(-10, 1),
    }[j]


def onset_template(j: int):
    d0 = first_nonzero_D(j)

    if j <= 4:
        return factor(
            onset_constant(j)
            * positive_ladder(j)
            * onset_linear(d0)
        )

    return factor(
        sp.Rational(7, 960)
        * terminal_ladder_j5()
        * (19859*K - 68077)
    )


# =============================================================================
# FACTOR MULTIPLICITIES
# =============================================================================

def linear_factor_multiplicity(
    expr: sp.Expr,
    linear_factor: sp.Expr,
):
    """
    Exact multiplicity of a linear factor.

    Zero polynomial gets None, because its multiplicity is not
    meaningful for this structural experiment.
    """

    if is_zero_poly(expr):
        return None

    P = sp.Poly(
        sp.expand(expr),
        K,
        domain="QQ",
    )

    Q = sp.Poly(
        sp.expand(linear_factor),
        K,
        domain="QQ",
    )

    multiplicity = 0
    current = P

    while True:

        quotient, remainder = sp.div(
            current,
            Q,
        )

        if not remainder.is_zero:
            break

        multiplicity += 1
        current = quotient

    return multiplicity


# =============================================================================
# SECTION 1
# =============================================================================

def section_1_exact_fixed_D_factoring():

    print("=" * 78)
    print("1. EXACT FIXED-D CROSS-k FACTORIZATION")
    print("=" * 78)

    for j in JS:

        print(f"j={j}")

        for d in DS:

            P = interpolate_K(j, d)

            if is_zero_poly(P):
                print(
                    f"  D={d}: ZERO POLYNOMIAL"
                )
                continue

            print(
                f"  D={d}: "
                f"degree={degree_in(P,K)}"
            )
            print(
                f"       P(K)={factor(P)}"
            )
            print(
                f"       roots={exact_roots(P)}"
            )

        print()


# =============================================================================
# SECTION 2
# =============================================================================

def section_2_factor_multiplicity_matrix():

    print("=" * 78)
    print("2. EXACT CANDIDATE-FACTOR MULTIPLICITY MATRIX")
    print("=" * 78)

    for j in range(5):

        print(f"j={j}")

        ladder = positive_ladder(j)

        ladder_factors = [
            K + m
            for m in range(
                j + 1,
                5,
            )
        ]

        for d in DS:

            P = interpolate_K(
                j,
                d,
            )

            if is_zero_poly(P):
                print(
                    f"  D={d}: ZERO"
                )
                continue

            pieces = []

            for f in ladder_factors:

                mult = linear_factor_multiplicity(
                    P,
                    f,
                )

                pieces.append(
                    f"{f}: {mult}"
                )

            L = onset_linear(d)

            L_mult = linear_factor_multiplicity(
                P,
                L,
            )

            pieces.append(
                f"{L}: {L_mult}"
            )

            print(
                f"  D={d}: "
                + ", ".join(pieces)
            )

        print()


# =============================================================================
# SECTION 3
# =============================================================================

def section_3_factor_survival_depth():

    print("=" * 78)
    print("3. FACTOR SURVIVAL DEPTH")
    print("=" * 78)

    for j in range(5):

        print(f"j={j}")

        factors = [
            K + m
            for m in range(
                j + 1,
                5,
            )
        ]

        for f in factors:

            surviving_D = []

            for d in DS:

                P = interpolate_K(
                    j,
                    d,
                )

                if nonzero_divides(
                    P,
                    f,
                ):
                    surviving_D.append(d)

            print(
                f"  factor={f}"
            )
            print(
                f"    nonzero-survives-at-D="
                f"{surviving_D}"
            )

            if surviving_D:
                print(
                    f"    max surviving D="
                    f"{max(surviving_D)}"
                )

                print(
                    f"    last consecutive surviving D="
                    f"{surviving_D[-1]}"
                )

        print()


# =============================================================================
# SECTION 4
# =============================================================================

def section_4_linear_factor_depth():

    print("=" * 78)
    print("4. LINEAR BOUNDARY FACTOR DEPTH")
    print("=" * 78)

    for j in range(5):

        print(f"j={j}")

        survival = []

        for d in DS:

            P = interpolate_K(
                j,
                d,
            )

            L = onset_linear(d)

            ok = nonzero_divides(
                P,
                L,
            )

            if ok:
                survival.append(d)

        print(
            f"  factor = 2K + D - 1"
        )
        print(
            f"  nonzero survival D="
            f"{survival}"
        )

        print()


# =============================================================================
# SECTION 5
# =============================================================================

def section_5_zero_vs_factor_separation():

    print("=" * 78)
    print("5. ZERO-SLICE / FACTOR-SLICE SEPARATION")
    print("=" * 78)

    for j in JS:

        zero_D = [
            d
            for d in DS
            if is_zero_poly(
                interpolate_K(j,d)
            )
        ]

        nonzero_D = [
            d
            for d in DS
            if not is_zero_poly(
                interpolate_K(j,d)
            )
        ]

        print(
            f"j={j}"
        )
        print(
            f"  zero slices={zero_D}"
        )
        print(
            f"  nonzero slices={nonzero_D}"
        )
        print()

    print(
        "IMPORTANT:"
    )
    print(
        "A factor is only counted as surviving when the "
        "fixed-D polynomial is nonzero."
    )
    print()


# =============================================================================
# SECTION 6
# =============================================================================

def section_6_factor_transition_points():

    print("=" * 78)
    print("6. FIRST FACTOR TRANSITION POINTS")
    print("=" * 78)

    for j in range(5):

        print(f"j={j}")

        factors = [
            K + m
            for m in range(
                j + 1,
                5,
            )
        ]

        for f in factors:

            statuses = []

            for d in DS:

                P = interpolate_K(
                    j,
                    d,
                )

                if is_zero_poly(P):
                    status = "ZERO"
                elif nonzero_divides(P,f):
                    status = "YES"
                else:
                    status = "NO"

                statuses.append(
                    (d,status)
                )

            transitions = []

            for idx in range(
                1,
                len(statuses),
            ):

                previous = statuses[idx-1][1]
                current = statuses[idx][1]

                if previous != current:
                    transitions.append(
                        (
                            statuses[idx-1][0],
                            previous,
                            statuses[idx][0],
                            current,
                        )
                    )

            print(
                f"  factor={f}"
            )

            print(
                f"    statuses={statuses}"
            )

            print(
                f"    transitions={transitions}"
            )

        print()


# =============================================================================
# SECTION 7
# =============================================================================

def section_7_onset_linear_transitions():

    print("=" * 78)
    print("7. 2K+D-1 TRANSITION TABLE")
    print("=" * 78)

    for j in range(5):

        print(f"j={j}")

        for d in DS:

            P = interpolate_K(
                j,
                d,
            )

            L = onset_linear(d)

            if is_zero_poly(P):

                print(
                    f"  D={d}: ZERO"
                )

                continue

            quotient, remainder, ok = divide_K(
                P,
                L,
            )

            if ok:

                print(
                    f"  D={d}: YES, "
                    f"quotient={quotient}"
                )

            else:

                print(
                    f"  D={d}: NO"
                )

        print()


# =============================================================================
# SECTION 8
# =============================================================================

def section_8_gcd_neighboring_D():

    print("=" * 78)
    print("8. GCD OF NEIGHBORING NONZERO D-SLICES")
    print("=" * 78)

    for j in JS:

        print(f"j={j}")

        for d1, d2 in zip(
            DS,
            DS[1:],
        ):

            P1 = interpolate_K(j,d1)
            P2 = interpolate_K(j,d2)

            if is_zero_poly(P1) or is_zero_poly(P2):
                print(
                    f"  D=({d1},{d2}): ZERO involved"
                )
                continue

            p1 = sp.Poly(
                sp.expand(P1),
                K,
                domain="QQ",
            )

            p2 = sp.Poly(
                sp.expand(P2),
                K,
                domain="QQ",
            )

            g = factor(
                sp.gcd(
                    p1,
                    p2,
                ).as_expr()
            )

            print(
                f"  D=({d1},{d2}): "
                f"gcd={g}"
            )

        print()


# =============================================================================
# SECTION 9
# =============================================================================

def section_9_gcd_removed_transition():

    print("=" * 78)
    print("9. GCD-REMOVED NEIGHBORING D-SLICES")
    print("=" * 78)

    for j in JS:

        print(f"j={j}")

        for d1, d2 in zip(
            DS,
            DS[1:],
        ):

            P1 = interpolate_K(j,d1)
            P2 = interpolate_K(j,d2)

            if is_zero_poly(P1) or is_zero_poly(P2):
                continue

            p1 = sp.Poly(
                sp.expand(P1),
                K,
                domain="QQ",
            )

            p2 = sp.Poly(
                sp.expand(P2),
                K,
                domain="QQ",
            )

            g = sp.gcd(
                p1,
                p2,
            )

            if g.as_expr() == 1:
                print(
                    f"  D=({d1},{d2}): gcd=1"
                )
                continue

            q1 = factor(
                sp.div(
                    p1,
                    g,
                )[0].as_expr()
            )

            q2 = factor(
                sp.div(
                    p2,
                    g,
                )[0].as_expr()
            )

            print(
                f"  D=({d1},{d2})"
            )
            print(
                f"    gcd={factor(g.as_expr())}"
            )
            print(
                f"    reduced D={d1}: {q1}"
            )
            print(
                f"    reduced D={d2}: {q2}"
            )

        print()


# =============================================================================
# SECTION 10
# =============================================================================

def section_10_support_depth_relation():

    print("=" * 78)
    print("10. SUPPORT DEGREE VS FACTOR DEPTH")
    print("=" * 78)

    print(
        "j | support roots | support degree | "
        "D0 | first nonzero factor depth"
    )
    print(
        "-" * 78
    )

    for j in JS:

        roots = common_support_roots(j)
        support = support_factor(j)
        sdeg = degree_in(
            support,
            D,
        )
        d0 = first_nonzero_D(j)

        if j <= 4:

            ladder = positive_ladder(j)

            surviving_nonzero = [
                d
                for d in DS
                if nonzero_divides(
                    interpolate_K(j,d),
                    ladder,
                )
            ]

        else:

            ladder = terminal_ladder_j5()

            surviving_nonzero = [
                d
                for d in DS
                if nonzero_divides(
                    interpolate_K(j,d),
                    ladder,
                )
            ]

        print(
            f"{j} | "
            f"{roots} | "
            f"{sdeg} | "
            f"{d0} | "
            f"{surviving_nonzero}"
        )

    print()


# =============================================================================
# SECTION 11
# =============================================================================

def section_11_simple_boundary_rule_test():

    print("=" * 78)
    print("11. SIMPLE BOUNDARY RULE TEST")
    print("=" * 78)

    print(
        "Test rule for j=0,...,4:"
    )
    print(
        "  A. fixed-D slice is zero before D0;"
    )
    print(
        "  B. onset K-ladder survives at D0;"
    )
    print(
        "  C. 2K+D-1 survives at D0;"
    )
    print(
        "  D. after D0, each factor may disappear."
    )
    print()

    for j in range(5):

        d0 = first_nonzero_D(j)

        print(
            f"j={j}: D0={d0}"
        )

        before = [
            d
            for d in DS
            if d < d0
            and is_zero_poly(
                interpolate_K(j,d)
            )
        ]

        print(
            f"  all earlier D zero="
            f"{before == [d for d in DS if d < d0]}"
        )

        ladder = positive_ladder(j)

        print(
            f"  ladder at D0="
            f"{nonzero_divides(
                interpolate_K(j,d0),
                ladder
            )}"
        )

        linear = onset_linear(d0)

        print(
            f"  linear at D0="
            f"{nonzero_divides(
                interpolate_K(j,d0),
                linear
            )}"
        )

        print()


# =============================================================================
# SECTION 12
# =============================================================================

def section_12_j5_terminal_transition():

    print("=" * 78)
    print("12. TERMINAL j=5 FACTOR TRANSITION")
    print("=" * 78)

    ladder = terminal_ladder_j5()

    print(
        f"terminal ladder={factor(ladder)}"
    )
    print()

    for d in DS:

        P = interpolate_K(
            5,
            d,
        )

        if is_zero_poly(P):

            print(
                f"D={d}: ZERO"
            )

            continue

        if nonzero_divides(
            P,
            ladder,
        ):

            Q, _, _ = divide_K(
                P,
                ladder,
            )

            print(
                f"D={d}: ladder survives"
            )
            print(
                f"  quotient={factor(Q)}"
            )

        else:

            print(
                f"D={d}: ladder disappears"
            )
            print(
                f"  P(K)={factor(P)}"
            )

    print()


# =============================================================================
# SECTION 13
# =============================================================================

def section_13_j5_root_transition():

    print("=" * 78)
    print("13. j=5 INTEGER ROOT TRANSITION")
    print("=" * 78)

    terminal_roots = [
        5,
        7,
        9,
        11,
    ]

    print(
        "Candidate terminal roots="
        f"{terminal_roots}"
    )
    print()

    for d in DS:

        P = interpolate_K(
            5,
            d,
        )

        if is_zero_poly(P):

            print(
                f"D={d}: ZERO"
            )

            continue

        statuses = []

        for root in terminal_roots:

            statuses.append(
                (
                    root,
                    sp.expand(
                        P.subs(
                            K,
                            root,
                        )
                    ) == 0,
                )
            )

        print(
            f"D={d}: "
            f"{statuses}"
        )

    print()


# =============================================================================
# SECTION 14
# =============================================================================

def section_14_pairwise_j_factor_overlap():

    print("=" * 78)
    print("14. PAIRWISE j FACTOR OVERLAP")
    print("=" * 78)

    for d in DS:

        print(
            f"D={d}"
        )

        for j1, j2 in combinations(
            JS,
            2,
        ):

            P1 = interpolate_K(
                j1,
                d,
            )

            P2 = interpolate_K(
                j2,
                d,
            )

            if (
                is_zero_poly(P1)
                or is_zero_poly(P2)
            ):
                continue

            p1 = sp.Poly(
                sp.expand(P1),
                K,
                domain="QQ",
            )

            p2 = sp.Poly(
                sp.expand(P2),
                K,
                domain="QQ",
            )

            g = factor(
                sp.gcd(
                    p1,
                    p2,
                ).as_expr()
            )

            if g != 1:

                print(
                    f"  (j={j1},j={j2}): "
                    f"gcd={g}"
                )

        print()


# =============================================================================
# SECTION 15
# =============================================================================

def section_15_onset_reconstruction():

    print("=" * 78)
    print("15. EXACT ONSET RECONSTRUCTION")
    print("=" * 78)

    for j in range(5):

        d0 = first_nonzero_D(j)

        actual = interpolate_K(
            j,
            d0,
        )

        expected = onset_template(j)

        print(
            f"j={j}"
        )
        print(
            f"  D0={d0}"
        )
        print(
            f"  actual={factor(actual)}"
        )
        print(
            f"  expected={factor(expected)}"
        )
        print(
            f"  match="
            f"{sp.expand(actual-expected)==0}"
        )
        print()


# =============================================================================
# SECTION 16
# =============================================================================

def section_16_factor_depth_table():

    print("=" * 78)
    print("16. COMPACT FACTOR-DEPTH TABLE")
    print("=" * 78)

    for j in range(5):

        d0 = first_nonzero_D(j)

        print(
            f"j={j}, D0={d0}"
        )

        if j <= 4:

            factors = [
                K + m
                for m in range(
                    j + 1,
                    5,
                )
            ]

            factors.append(
                onset_linear(d0)
            )

        else:

            factors = [
                K - 5,
                K - 7,
                K - 9,
                K - 11,
            ]

        for f in factors:

            survival = []

            for d in DS:

                P = interpolate_K(
                    j,
                    d,
                )

                if is_zero_poly(P):
                    state = "ZERO"
                elif nonzero_divides(P,f):
                    state = "YES"
                else:
                    state = "NO"

                survival.append(
                    (d,state)
                )

            print(
                f"  factor={f}"
            )
            print(
                f"    {survival}"
            )

        print()


# =============================================================================
# SECTION 17
# =============================================================================

def section_17_degree_and_factor_transition():

    print("=" * 78)
    print("17. DEGREE DROP VS FACTOR TRANSITION")
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

            if is_zero_poly(P):

                print(
                    f"  D={d}: ZERO"
                )

                continue

            deg = degree_in(
                P,
                K,
            )

            roots = exact_roots(P)

            print(
                f"  D={d}: "
                f"degree={deg}, "
                f"roots={roots}"
            )

        print()


# =============================================================================
# SECTION 18
# =============================================================================

def section_18_exact_interpolation_check():

    print("=" * 78)
    print("18. EXACT INTERPOLATION CHECK")
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
                    P.subs(
                        D,
                        d,
                    )
                )

                actual = E(
                    k,
                    j,
                    d,
                )

                if predicted != actual:

                    failures += 1

                    print(
                        "FAIL:"
                        f" k={k}"
                        f" j={j}"
                        f" D={d}"
                        f" actual={actual}"
                        f" predicted={predicted}"
                    )

    print(
        f"tested={tested}"
    )
    print(
        f"interpolation failures="
        f"{failures}"
    )
    print()


# =============================================================================
# SECTION 19
# =============================================================================

def section_19_support_divisibility_check():

    print("=" * 78)
    print("19. SYMBOLIC SUPPORT DIVISIBILITY CHECK")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in JS:

        S = support_factor(j)

        print(
            f"j={j}"
        )
        print(
            f"  support={factor(S)}"
        )

        for k in KS:

            P = sp.Poly(
                interpolate_D(
                    k,
                    j,
                ),
                D,
                domain="QQ",
            )

            SD = sp.Poly(
                S,
                D,
                domain="QQ",
            )

            _, remainder = sp.div(
                P,
                SD,
            )

            tested += 1

            if not remainder.is_zero:

                failures += 1

                print(
                    f"  FAIL k={k}: "
                    f"remainder="
                    f"{factor(remainder.as_expr())}"
                )

    print(
        f"tested={tested}"
    )
    print(
        f"support divisibility failures="
        f"{failures}"
    )
    print()


# =============================================================================
# SECTION 20
# =============================================================================

def section_20_final_transition_summary():

    print("=" * 78)
    print("20. FINAL TRANSITION SUMMARY")
    print("=" * 78)

    print(
        "j | support | D0 | nonzero ladder D | "
        "nonzero linear D"
    )
    print(
        "-" * 78
    )

    for j in range(5):

        S = support_factor(j)
        d0 = first_nonzero_D(j)

        ladder = positive_ladder(j)

        ladder_D = [
            d
            for d in DS
            if nonzero_divides(
                interpolate_K(j,d),
                ladder,
            )
        ]

        linear_D = [
            d
            for d in DS
            if nonzero_divides(
                interpolate_K(j,d),
                onset_linear(d),
            )
        ]

        print(
            f"{j} | "
            f"{factor(S)} | "
            f"{d0} | "
            f"{ladder_D} | "
            f"{linear_D}"
        )

    print()

    print(
        "j=5 terminal ladder:"
    )

    terminal = terminal_ladder_j5()

    terminal_D = [
        d
        for d in DS
        if nonzero_divides(
            interpolate_K(5,d),
            terminal,
        )
    ]

    print(
        f"  ladder={factor(terminal)}"
    )
    print(
        f"  nonzero ladder D={terminal_D}"
    )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    validate_data()

    section_1_exact_fixed_D_factoring()
    section_2_factor_multiplicity_matrix()
    section_3_factor_survival_depth()
    section_4_linear_factor_depth()
    section_5_zero_vs_factor_separation()
    section_6_factor_transition_points()
    section_7_onset_linear_transitions()
    section_8_gcd_neighboring_D()
    section_9_gcd_removed_transition()
    section_10_support_depth_relation()
    section_11_simple_boundary_rule_test()
    section_12_j5_terminal_transition()
    section_13_j5_root_transition()
    section_14_pairwise_j_factor_overlap()
    section_15_onset_reconstruction()
    section_16_factor_depth_table()
    section_17_degree_and_factor_transition()
    section_18_exact_interpolation_check()
    section_19_support_divisibility_check()
    section_20_final_transition_summary()


if __name__ == "__main__":
    main()

