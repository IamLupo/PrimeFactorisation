from __future__ import annotations

import math
from typing import Optional

import sympy as sp


# =============================================================================
# R=5 CROSS-k BOUNDARY FACTORIZATION EXPERIMENT
# =============================================================================
#
# IMPORTANT:
#   - The supplied discrepancies are treated as the exact input data.
#   - No old-law subtraction is performed.
#   - No r=6.
#   - No full pq-kernel expansion.
#   - No replacement universal r,j law is asserted.
#
# Main object:
#
#     E_j(k,D)
#
# The experiment investigates the first nonzero D slice for each j and its
# exact cross-k factorization.
#
# =============================================================================


# =============================================================================
# SYMBOLS
# =============================================================================

D, K = sp.symbols("D K")


# =============================================================================
# EXACT INPUT DATA
# =============================================================================
#
# NOTE:
# k=13, j=5, D=14 IS INCLUDED HERE:
#
#     D=14 -> 532252
#
# This was the missing entry responsible for the validation error:
#
#     bad D-grid for k=13, j=5: [6, 8, 10, 12, 16]
#
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
            16: 1048278,
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

            # IMPORTANT: DO NOT REMOVE THIS ENTRY.
            14: 532252,

            16: 533876,
        },
    },
}


# =============================================================================
# OPTIONAL FRESH DATA
# =============================================================================
#
# Add genuinely new k-values here later.
#
# Example:
#
# FRESH_DATA = {
#     15: {
#         0: {6: ...},
#         1: {8: ...},
#         ...
#     }
# }
#
# These points are never used to determine the original onset templates.
# They are used only for out-of-sample validation.
# =============================================================================

FRESH_DATA: dict[int, dict[int, dict[int, int]]] = {}


# =============================================================================
# GRID DEFINITIONS
# =============================================================================

KS = sorted(DATA.keys())
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
# HELPERS
# =============================================================================

def E(k: int, j: int, d: int) -> sp.Integer:
    return sp.Integer(DATA[k][j][d])


def fresh_E(k: int, j: int, d: int) -> sp.Integer:
    return sp.Integer(FRESH_DATA[k][j][d])


def degree_in(expr: sp.Expr, variable) -> int:
    expr = sp.expand(expr)

    if expr == 0:
        return -math.inf

    return int(
        sp.Poly(
            expr,
            variable,
            domain="QQ",
        ).degree()
    )


def factor(expr: sp.Expr) -> sp.Expr:
    return sp.factor(
        sp.cancel(
            sp.expand(expr)
        )
    )


def exact_roots(expr: sp.Expr) -> list[sp.Expr]:
    if expr == 0:
        return []

    return sorted(
        list(
            sp.roots(
                sp.Poly(
                    expr,
                    K,
                    domain="QQ",
                ).as_expr(),
                K,
            ).keys()
        ),
        key=sp.default_sort_key,
    )


def rational_roots(expr: sp.Expr) -> list[sp.Rational]:
    return [
        sp.Rational(r)
        for r in exact_roots(expr)
        if r.is_Rational
    ]


# =============================================================================
# ROBUST DATA VALIDATION
# =============================================================================
#
# Instead of using a single set comparison and immediately raising on the
# first error, this reports EVERY missing and extra D-entry.
# =============================================================================

def validate_data() -> None:

    expected_js = set(JS)
    expected_ds = set(DS)

    problems = []

    for k in KS:

        actual_js = set(DATA[k])

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
            "DATA GRID PROBLEMS:"
        )

        for problem in problems:
            print(
                f"  {problem}"
            )

        raise RuntimeError(
            "Input data grid is incomplete or inconsistent. "
            "See the complete list above."
        )

    print(
        "grid status = OK"
    )

    if FRESH_DATA:
        print(
            f"fresh k values = "
            f"{sorted(FRESH_DATA)}"
        )
    else:
        print(
            "fresh k values = []"
        )

    print()


# =============================================================================
# INTERPOLATION IN D
# =============================================================================

def interpolate_D(
    j: int,
    k: int,
) -> sp.Expr:

    points = [
        (
            sp.Integer(d),
            E(k, j, d),
        )
        for d in DS
    ]

    return factor(
        sp.interpolate(
            points,
            D,
        )
    )


# =============================================================================
# INTERPOLATION IN K AT FIXED D
# =============================================================================

def interpolate_K_at_D(
    j: int,
    d: int,
) -> sp.Expr:

    points = [
        (
            sp.Integer(k),
            E(k, j, d),
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
# COMMON SUPPORT ROOTS
# =============================================================================

def common_support_roots(
    j: int,
) -> list[int]:

    return [
        d
        for d in DS
        if all(
            E(k, j, d) == 0
            for k in KS
        )
    ]


def support_factor(
    j: int,
) -> sp.Expr:

    result = sp.Integer(1)

    for root in common_support_roots(j):
        result *= (
            D - root
        )

    return factor(result)


# =============================================================================
# FIRST NONZERO D
# =============================================================================

def first_nonzero_D(
    j: int,
) -> Optional[int]:

    for d in DS:

        if any(
            E(k, j, d) != 0
            for k in KS
        ):
            return d

    return None


# =============================================================================
# SUPPORT-REMOVED QUOTIENT
# =============================================================================

def support_quotient(
    j: int,
    k: int,
) -> sp.Expr:

    P = sp.Poly(
        interpolate_D(j, k),
        D,
        domain="QQ",
    )

    S = sp.Poly(
        support_factor(j),
        D,
        domain="QQ",
    )

    Q, R = sp.div(
        P,
        S,
    )

    if not R.is_zero:
        raise RuntimeError(
            f"support division failed: "
            f"k={k}, j={j}, "
            f"remainder={R.as_expr()}"
        )

    return factor(
        Q.as_expr()
    )


# =============================================================================
# ONSET TEMPLATE FOR j=0,...,4
# =============================================================================
#
# Candidate:
#
#   P_j(K)
#     =
#   C_j *
#   product_{m=j+1}^4 (K+m) *
#   (2K + D0 - 1)
#
# =============================================================================

def onset_base(
    j: int,
    d0: int,
) -> sp.Expr:

    result = sp.Integer(1)

    for m in range(
        j + 1,
        5,
    ):
        result *= (
            K + m
        )

    result *= (
        2*K + d0 - 1
    )

    return sp.expand(result)


def extract_onset_constant(
    j: int,
    d0: int,
    P: sp.Expr,
) -> Optional[sp.Rational]:

    base = sp.Poly(
        onset_base(j, d0),
        K,
        domain="QQ",
    )

    p = sp.Poly(
        P,
        K,
        domain="QQ",
    )

    if p.degree() != base.degree():
        return None

    return sp.Rational(
        p.LC(),
        base.LC(),
    )


def onset_candidate(
    j: int,
    d0: int,
    constant: sp.Rational,
) -> sp.Expr:

    return factor(
        constant *
        onset_base(j, d0)
    )


def exact_onset_test(
    j: int,
):

    d0 = first_nonzero_D(j)

    if d0 is None:
        return (
            False,
            None,
            None,
            None,
            None,
        )

    actual = interpolate_K_at_D(
        j,
        d0,
    )

    base = onset_base(
        j,
        d0,
    )

    constant = extract_onset_constant(
        j,
        d0,
        actual,
    )

    if constant is None:
        return (
            False,
            d0,
            actual,
            base,
            None,
        )

    candidate = onset_candidate(
        j,
        d0,
        constant,
    )

    ok = (
        sp.expand(
            actual - candidate
        ) == 0
    )

    return (
        ok,
        d0,
        actual,
        base,
        constant,
    )


# =============================================================================
# SECTION 1
# =============================================================================

def section_1_exact_support() -> None:

    print("=" * 78)
    print("1. EXACT SUPPORT ROOT STRUCTURE")
    print("=" * 78)

    for j in JS:

        roots = common_support_roots(j)

        print(
            f"j={j}"
        )
        print(
            f"  common observed roots = "
            f"{roots}"
        )
        print(
            f"  support = "
            f"{support_factor(j)}"
        )
        print(
            f"  support degree = "
            f"{degree_in(support_factor(j),D)}"
        )

        for root in roots:

            multiplicities = {}

            for k in KS:

                P = sp.Poly(
                    interpolate_D(j,k),
                    D,
                    domain="QQ",
                )

                multiplicity = 0
                Q = P

                root_poly = sp.Poly(
                    D-root,
                    D,
                    domain="QQ",
                )

                while True:

                    quotient, remainder = sp.div(
                        Q,
                        root_poly,
                    )

                    if remainder.is_zero:

                        multiplicity += 1
                        Q = quotient

                    else:
                        break

                multiplicities[k] = multiplicity

            print(
                f"  D={root} multiplicities="
                f"{multiplicities}"
            )

        print()


# =============================================================================
# SECTION 2
# =============================================================================

def section_2_cross_k_degree_map() -> None:

    print("=" * 78)
    print("2. CROSS-k DEGREE DROP MAP")
    print("=" * 78)

    for j in JS:

        print(
            f"j={j}"
        )

        drops = []

        for d in DS:

            P = interpolate_K_at_D(
                j,
                d,
            )

            deg = degree_in(
                P,
                K,
            )

            print(
                f"  D={d}: degree={deg}"
            )

            if deg != 5 and P != 0:
                drops.append(
                    (d, deg)
                )

        print(
            f"  degree-drop slices="
            f"{drops}"
        )

        print()


# =============================================================================
# SECTION 3
# =============================================================================

def section_3_exact_cross_k_factorization() -> None:

    print("=" * 78)
    print("3. EXACT FACTORIZATION OF CROSS-k D-SLICES")
    print("=" * 78)

    for j in JS:

        print(
            f"j={j}"
        )

        for d in DS:

            P = interpolate_K_at_D(
                j,
                d,
            )

            deg = degree_in(
                P,
                K,
            )

            print(
                f"  D={d}:"
            )
            print(
                f"    degree={deg}"
            )
            print(
                f"    E(K)={factor(P)}"
            )
            print(
                f"    roots={exact_roots(P)}"
            )

        print()


# =============================================================================
# SECTION 4
# =============================================================================

def section_4_distinguished_slices() -> None:

    print("=" * 78)
    print("4. DISTINGUISHED FIRST-NONZERO D SLICES")
    print("=" * 78)

    for j in JS:

        d0 = first_nonzero_D(j)

        print(
            f"j={j}: first_nonzero_D={d0}"
        )

        if d0 is not None:

            P = interpolate_K_at_D(
                j,
                d0,
            )

            print(
                f"  degree={degree_in(P,K)}"
            )
            print(
                f"  factor={factor(P)}"
            )
            print(
                f"  roots={exact_roots(P)}"
            )

        print()


# =============================================================================
# SECTION 5
# =============================================================================

def section_5_onset_template_tests() -> None:

    print("=" * 78)
    print("5. EXACT ONSET TEMPLATE TESTS")
    print("=" * 78)

    print(
        "Candidate for j=0,...,4:"
    )
    print(
        "  P_j(K) = C_j * "
        "product_{m=j+1}^{4}(K+m) * "
        "(2K + D0 - 1)"
    )
    print()

    for j in range(5):

        (
            ok,
            d0,
            actual,
            base,
            constant,
        ) = exact_onset_test(j)

        print(
            f"j={j}"
        )
        print(
            f"  D0={d0}"
        )
        print(
            f"  actual degree="
            f"{degree_in(actual,K)}"
        )
        print(
            f"  expected degree="
            f"{degree_in(base,K)}"
        )
        print(
            f"  C_j={constant}"
        )
        print(
            f"  base={factor(base)}"
        )
        print(
            f"  candidate="
            f"{factor(onset_candidate(j,d0,constant))}"
            if constant is not None
            else "  candidate=<unavailable>"
        )
        print(
            f"  actual={factor(actual)}"
        )
        print(
            f"  EXACT MATCH={ok}"
        )
        print()


# =============================================================================
# SECTION 6
# =============================================================================

def section_6_onset_constants() -> None:

    print("=" * 78)
    print("6. ONSET CONSTANTS")
    print("=" * 78)

    constants = {}

    for j in range(5):

        d0 = first_nonzero_D(j)

        actual = interpolate_K_at_D(
            j,
            d0,
        )

        c = extract_onset_constant(
            j,
            d0,
            actual,
        )

        constants[j] = c

        print(
            f"j={j}: C_j={c}"
        )

    print()

    print(
        "successive ratios:"
    )

    for j in range(1,5):

        ratio = factor(
            constants[j] /
            constants[j-1]
        )

        print(
            f"  C_{j}/C_{j-1}={ratio}"
        )

    print()


# =============================================================================
# SECTION 7
# =============================================================================
#
# Search simple normalizations for C_j.
#
# =============================================================================

def section_7_constant_normalization_search() -> None:

    print("=" * 78)
    print("7. SIMPLE C_j NORMALIZATION SEARCH")
    print("=" * 78)

    constants = {}

    for j in range(5):

        d0 = first_nonzero_D(j)
        P = interpolate_K_at_D(j,d0)

        constants[j] = extract_onset_constant(
            j,
            d0,
            P,
        )

    candidates = [
        (
            "(-1)^j",
            lambda j: sp.Integer(-1) ** j,
        ),
        (
            "2^j",
            lambda j: sp.Integer(2) ** j,
        ),
        (
            "(-2)^j",
            lambda j: sp.Integer(-2) ** j,
        ),
        (
            "j+1",
            lambda j: sp.Integer(j+1),
        ),
        (
            "j!",
            lambda j: sp.factorial(j),
        ),
        (
            "(j+1)!",
            lambda j: sp.factorial(j+1),
        ),
        (
            "C(4,j)",
            lambda j: sp.binomial(4,j),
        ),
        (
            "C(5,j)",
            lambda j: sp.binomial(5,j),
        ),
        (
            "C(4,j+1)",
            lambda j: sp.binomial(4,j+1),
        ),
        (
            "C(5,j+1)",
            lambda j: sp.binomial(5,j+1),
        ),
    ]

    for name, fn in candidates:

        values = [
            factor(
                constants[j] /
                fn(j)
            )
            for j in range(5)
        ]

        constant = all(
            v == values[0]
            for v in values
        )

        print(
            f"{name}: "
            f"normalized={values} "
            f"constant={constant}"
        )

    print()


# =============================================================================
# SECTION 8
# =============================================================================

def section_8_linear_factor() -> None:

    print("=" * 78)
    print("8. ONSET LINEAR FACTOR TEST")
    print("=" * 78)

    for j in JS:

        d0 = first_nonzero_D(j)

        P = sp.Poly(
            interpolate_K_at_D(j,d0),
            K,
            domain="QQ",
        )

        linear = 2*K + d0 - 1

        quotient, remainder = sp.div(
            P,
            sp.Poly(
                linear,
                K,
                domain="QQ",
            ),
        )

        print(
            f"j={j}, D0={d0}"
        )
        print(
            f"  linear factor={linear}"
        )
        print(
            f"  divides={remainder.is_zero}"
        )

        if remainder.is_zero:

            print(
                f"  quotient="
                f"{factor(quotient.as_expr())}"
            )

        print()


# =============================================================================
# SECTION 9
# =============================================================================

def section_9_consecutive_factor_ladder() -> None:

    print("=" * 78)
    print("9. CONSECUTIVE K-FACTOR LADDER")
    print("=" * 78)

    for j in range(5):

        d0 = first_nonzero_D(j)

        P = sp.Poly(
            interpolate_K_at_D(j,d0),
            K,
            domain="QQ",
        )

        ladder = sp.Integer(1)

        for m in range(j+1,5):
            ladder *= K+m

        quotient, remainder = sp.div(
            P,
            sp.Poly(
                ladder,
                K,
                domain="QQ",
            ),
        )

        print(
            f"j={j}:"
        )
        print(
            f"  ladder={factor(ladder)}"
        )
        print(
            f"  divides={remainder.is_zero}"
        )

        if remainder.is_zero:
            print(
                f"  quotient="
                f"{factor(quotient.as_expr())}"
            )

        print()


# =============================================================================
# SECTION 10
# =============================================================================

def section_10_j5_terminal_structure() -> None:

    print("=" * 78)
    print("10. J=5 TERMINAL STRUCTURE")
    print("=" * 78)

    j = 5
    d0 = first_nonzero_D(j)

    P = interpolate_K_at_D(
        j,
        d0,
    )

    print(
        f"j=5: D0={d0}"
    )
    print(
        f"degree={degree_in(P,K)}"
    )
    print(
        f"P(K)={factor(P)}"
    )
    print(
        f"roots={exact_roots(P)}"
    )
    print(
        f"rational roots={rational_roots(P)}"
    )

    integer_roots = [
        r
        for r in range(-20,21)
        if sp.expand(
            P.subs(K,r)
        ) == 0
    ]

    print(
        f"integer roots="
        f"{integer_roots}"
    )

    if integer_roots:

        ladder = sp.prod(
            K-r
            for r in integer_roots
        )

        quotient, remainder = sp.div(
            sp.Poly(P,K,domain="QQ"),
            sp.Poly(ladder,K,domain="QQ"),
        )

        print(
            f"integer-root product="
            f"{factor(ladder)}"
        )
        print(
            f"divides={remainder.is_zero}"
        )

        if remainder.is_zero:
            print(
                f"remaining factor="
                f"{factor(quotient.as_expr())}"
            )

    print()


# =============================================================================
# SECTION 11
# =============================================================================

def section_11_j5_D14_D16() -> None:

    print("=" * 78)
    print("11. J=5 D=14 / D=16 SIDE-BY-SIDE")
    print("=" * 78)

    P14 = interpolate_K_at_D(
        5,
        14,
    )

    P16 = interpolate_K_at_D(
        5,
        16,
    )

    print(
        f"D=14:"
    )
    print(
        f"  degree={degree_in(P14,K)}"
    )
    print(
        f"  P(K)={factor(P14)}"
    )
    print(
        f"  roots={exact_roots(P14)}"
    )

    print()

    print(
        f"D=16:"
    )
    print(
        f"  degree={degree_in(P16,K)}"
    )
    print(
        f"  P(K)={factor(P16)}"
    )
    print(
        f"  roots={exact_roots(P16)}"
    )

    gcd = sp.factor(
        sp.gcd(
            sp.Poly(P14,K),
            sp.Poly(P16,K),
        ).as_expr()
    )

    print()

    print(
        f"gcd={gcd}"
    )
    print()


# =============================================================================
# SECTION 12
# =============================================================================

def section_12_onset_gcd() -> None:

    print("=" * 78)
    print("12. GCD BETWEEN SUCCESSIVE ONSET SLICES")
    print("=" * 78)

    onset = {}

    for j in JS:

        d0 = first_nonzero_D(j)

        onset[j] = interpolate_K_at_D(
            j,
            d0,
        )

    for j in range(5):

        g = sp.factor(
            sp.gcd(
                sp.Poly(
                    onset[j],
                    K,
                ),
                sp.Poly(
                    onset[j+1],
                    K,
                ),
            ).as_expr()
        )

        print(
            f"j={j}, j+1={j+1}: gcd={g}"
        )

    print()


# =============================================================================
# SECTION 13
# =============================================================================

def section_13_successive_ratios() -> None:

    print("=" * 78)
    print("13. SUCCESSIVE ONSET RATIOS")
    print("=" * 78)

    onset = {}

    for j in JS:

        d0 = first_nonzero_D(j)

        onset[j] = interpolate_K_at_D(
            j,
            d0,
        )

    for j in range(5):

        ratio = factor(
            sp.cancel(
                onset[j+1] /
                onset[j]
            )
        )

        print(
            f"P_{j+1}/P_{j}="
            f"{ratio}"
        )

    print()


# =============================================================================
# SECTION 14
# =============================================================================

def section_14_onset_D_pattern() -> None:

    print("=" * 78)
    print("14. FIRST-NONZERO D PATTERN")
    print("=" * 78)

    values = []

    for j in JS:

        d0 = first_nonzero_D(j)

        values.append(
            (j,d0)
        )

        print(
            f"j={j}: D0={d0}"
        )

    print()

    print(
        "increments:"
    )

    for i in range(1,len(values)):

        a = values[i-1]
        b = values[i]

        if a[1] is not None and b[1] is not None:

            print(
                f"  {b[0]}-{a[0]}: "
                f"{b[1]-a[1]}"
            )

    print()


# =============================================================================
# SECTION 15
# =============================================================================

def section_15_support_and_onset_degrees() -> None:

    print("=" * 78)
    print("15. SUPPORT DEGREE / ONSET DEGREE")
    print("=" * 78)

    for j in JS:

        S = support_factor(j)

        support_deg = degree_in(
            S,
            D,
        )

        d0 = first_nonzero_D(j)

        P = interpolate_K_at_D(
            j,
            d0,
        )

        onset_deg = degree_in(
            P,
            K,
        )

        print(
            f"j={j}: "
            f"support_degree={support_deg}, "
            f"onset_degree={onset_deg}"
        )

    print()


# =============================================================================
# SECTION 16
# =============================================================================

def section_16_zero_lattice() -> None:

    print("=" * 78)
    print("16. ZERO LATTICE")
    print("=" * 78)

    for j in JS:

        print(
            f"j={j}"
        )

        for k in KS:

            zeros = [
                d
                for d in DS
                if E(k,j,d) == 0
            ]

            print(
                f"  k={k}: "
                f"{zeros}"
            )

        print()


# =============================================================================
# SECTION 17
# =============================================================================

def section_17_support_reconstruction() -> None:

    print("=" * 78)
    print("17. EXACT SUPPORT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in JS:

        S = support_factor(j)

        for k in KS:

            P = interpolate_D(
                j,
                k,
            )

            Q = support_quotient(
                j,
                k,
            )

            reconstructed = factor(
                S * Q
            )

            tested += 1

            if sp.expand(
                reconstructed - P
            ) != 0:

                failures += 1

                print(
                    "FAIL:",
                    f"k={k}",
                    f"j={j}",
                    f"original={P}",
                    f"reconstructed={reconstructed}",
                )

    print(
        f"tested={tested}"
    )
    print(
        f"failures={failures}"
    )
    print()


# =============================================================================
# SECTION 18
# =============================================================================

def section_18_interpolation_check() -> None:

    print("=" * 78)
    print("18. EXACT INTERPOLATION CHECK")
    print("=" * 78)

    tested = 0
    failures = 0

    for k in KS:

        for j in JS:

            P = interpolate_D(
                j,
                k,
            )

            for d in DS:

                tested += 1

                actual = E(
                    k,
                    j,
                    d,
                )

                predicted = sp.expand(
                    P.subs(
                        D,
                        d,
                    )
                )

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
        f"interpolation failures={failures}"
    )
    print()


# =============================================================================
# SECTION 19
# =============================================================================
#
# Fresh-k validation of the j<=4 onset template.
#
# =============================================================================

def candidate_onset(
    j: int,
) -> Optional[sp.Expr]:

    if j > 4:
        return None

    d0 = first_nonzero_D(j)

    P = interpolate_K_at_D(
        j,
        d0,
    )

    c = extract_onset_constant(
        j,
        d0,
        P,
    )

    if c is None:
        return None

    candidate = onset_candidate(
        j,
        d0,
        c,
    )

    if sp.expand(
        candidate - P
    ) != 0:
        return None

    return candidate


def section_19_fresh_validation() -> None:

    print("=" * 78)
    print("19. FRESH-k OUT-OF-SAMPLE VALIDATION")
    print("=" * 78)

    if not FRESH_DATA:

        print(
            "No fresh data supplied."
        )
        print()

        return

    fresh_ks = sorted(
        FRESH_DATA
    )

    print(
        f"fresh k values={fresh_ks}"
    )
    print()

    tested = 0
    failures = 0

    for j in range(5):

        candidate = candidate_onset(
            j
        )

        if candidate is None:
            print(
                f"j={j}: "
                "no exact onset template"
            )
            continue

        d0 = first_nonzero_D(j)

        print(
            f"j={j}, D0={d0}"
        )
        print(
            f"  candidate="
            f"{factor(candidate)}"
        )

        for k in fresh_ks:

            if j not in FRESH_DATA[k]:
                continue

            if d0 not in FRESH_DATA[k][j]:
                continue

            actual = fresh_E(
                k,
                j,
                d0,
            )

            predicted = factor(
                candidate.subs(
                    K,
                    k,
                )
            )

            tested += 1

            ok = (
                actual == predicted
            )

            if not ok:
                failures += 1

            print(
                f"  k={k}: "
                f"actual={actual}, "
                f"predicted={predicted}, "
                f"match={ok}"
            )

        print()

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

def section_20_compact_table() -> None:

    print("=" * 78)
    print("20. KEY STRUCTURAL TABLE")
    print("=" * 78)

    print(
        "j | support | support_deg | D0 | "
        "onset_degree | onset_factor"
    )
    print(
        "-" * 78
    )

    for j in JS:

        S = support_factor(j)
        support_deg = degree_in(
            S,
            D,
        )

        d0 = first_nonzero_D(j)

        if d0 is None:

            print(
                f"{j} | {S} | {support_deg} | "
                f"-- | -- | none"
            )

            continue

        P = interpolate_K_at_D(
            j,
            d0,
        )

        print(
            f"{j} | "
            f"{factor(S)} | "
            f"{support_deg} | "
            f"{d0} | "
            f"{degree_in(P,K)} | "
            f"{factor(P)}"
        )

    print()


# =============================================================================
# SECTION 21
# =============================================================================

def section_21_final_diagnostic() -> None:

    print("=" * 78)
    print("21. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The input discrepancies are used directly."
    )
    print()

    print(
        "First nonzero D by j:"
    )

    for j in JS:

        d0 = first_nonzero_D(j)

        print(
            f"  j={j}: D0={d0}"
        )

    print()

    print(
        "For j=0,...,4 the onset slices are tested against:"
    )
    print()
    print(
        "  P_j(K) = C_j"
        " * product_{m=j+1}^{4}(K+m)"
        " * (2K + D0 - 1)"
    )
    print()

    for j in range(5):

        ok, d0, actual, base, c = (
            exact_onset_test(j)
        )

        print(
            f"  j={j}: "
            f"match={ok}, "
            f"C_j={c}"
        )

    print()

    P5 = interpolate_K_at_D(
        5,
        first_nonzero_D(5),
    )

    print(
        "j=5 is kept separate:"
    )
    print(
        f"  D0={first_nonzero_D(5)}"
    )
    print(
        f"  degree={degree_in(P5,K)}"
    )
    print(
        f"  factor={factor(P5)}"
    )
    print(
        f"  roots={exact_roots(P5)}"
    )

    print()

    print(
        "No generic degree-5 interpolation is treated as a law."
    )
    print(
        "Fresh k-values are required for genuine out-of-sample"
        " confirmation."
    )
    print()

    print(
        "No r=6."
    )
    print(
        "No full pq-kernel expansion."
    )
    print(
        "No replacement universal r,j formula."
    )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    validate_data()

    section_1_exact_support()
    section_2_cross_k_degree_map()
    section_3_exact_cross_k_factorization()
    section_4_distinguished_slices()
    section_5_onset_template_tests()
    section_6_onset_constants()
    section_7_constant_normalization_search()
    section_8_linear_factor()
    section_9_consecutive_factor_ladder()
    section_10_j5_terminal_structure()
    section_11_j5_D14_D16()
    section_12_onset_gcd()
    section_13_successive_ratios()
    section_14_onset_D_pattern()
    section_15_support_and_onset_degrees()
    section_16_zero_lattice()
    section_17_support_reconstruction()
    section_18_interpolation_check()
    section_19_fresh_validation()
    section_20_compact_table()
    section_21_final_diagnostic()


if __name__ == "__main__":
    main()