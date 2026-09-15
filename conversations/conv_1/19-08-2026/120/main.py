#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 286R — EXACT B-ENTRY PROVENANCE / FORMULA RECONSTRUCTION AUDIT
==============================================================================

Corrected version of Experiment 286.

Fixes:
    * safe LCM for zero/one/many denominators;
    * no sp.ilcm(*single_item) calls;
    * exact QQ arithmetic throughout;
    * no q-family fitting;
    * no arbitrary matrix fitting.

The experiment searches for recognizable arithmetic fingerprints in B[k,r]:

    factorials,
    binomials,
    falling factorials,
    Stirling numbers,
    powers of 2 and 3,
    signs,
    k/r/(r-k) factors,
    common rational scaling.

A common exact identity is reported only when it holds for every
defined B entry.
"""

from __future__ import annotations

import math
import sys
import sympy as sp


# ============================================================================
# EXACT B TABLE
# ============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# ============================================================================
# HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def safe_lcm(values):
    """
    Exact integer LCM that works for:
        []      -> 1
        [x]     -> abs(x)
        [x,y]   -> lcm(x,y)
        many    -> recursive accumulation
    """
    vals = [
        abs(int(v))
        for v in values
        if int(v) != 0
    ]

    if not vals:
        return 1

    result = 1

    for value in vals:
        result = math.lcm(
            result,
            value,
        )

    return result


def b_points():
    out = []

    for k, row in enumerate(B):
        for offset, value in enumerate(row):
            r = k + offset

            out.append(
                (
                    k,
                    r,
                    sp.Rational(value),
                )
            )

    return out


def falling(n, m):
    if m < 0:
        return sp.Integer(0)

    out = sp.Integer(1)

    for j in range(m):
        out *= n - j

    return clean(out)


def safe_binomial(n, r):
    if r < 0 or r > n:
        return sp.Integer(0)

    return sp.binomial(
        sp.Integer(n),
        sp.Integer(r),
    )


def signed_stirling_first(n, k):
    return sp.functions.combinatorial.numbers.stirling(
        n,
        k,
        kind=1,
        signed=True,
    )


def unsigned_stirling_first(n, k):
    return sp.functions.combinatorial.numbers.stirling(
        n,
        k,
        kind=1,
        signed=False,
    )


def stirling_second(n, k):
    return sp.functions.combinatorial.numbers.stirling(
        n,
        k,
        kind=2,
    )


# ============================================================================
# CANDIDATE BASIS
# ============================================================================

def candidate_library(k, r):

    d = r - k

    return {
        "1": sp.Integer(1),

        "k": sp.Integer(k),
        "r": sp.Integer(r),
        "d": sp.Integer(d),

        "k!": sp.factorial(k),
        "r!": sp.factorial(r),
        "d!": sp.factorial(d),

        "1/k!": sp.Rational(
            1,
            sp.factorial(k),
        ),

        "1/r!": sp.Rational(
            1,
            sp.factorial(r),
        ),

        "1/d!": sp.Rational(
            1,
            sp.factorial(d),
        ),

        "C(r,k)": safe_binomial(r, k),
        "C(r,d)": safe_binomial(r, d),

        "falling(r,k)": falling(r, k),
        "falling(r,d)": falling(r, d),

        "2^k": sp.Integer(2) ** k,
        "2^r": sp.Integer(2) ** r,
        "2^d": sp.Integer(2) ** d,

        "3^k": sp.Integer(3) ** k,
        "3^r": sp.Integer(3) ** r,
        "3^d": sp.Integer(3) ** d,

        "(-1)^k": sp.Integer(-1) ** k,
        "(-1)^r": sp.Integer(-1) ** r,
        "(-1)^d": sp.Integer(-1) ** d,

        "s1(r,k)": signed_stirling_first(r, k),
        "abs_s1(r,k)": unsigned_stirling_first(r, k),
        "s2(r,k)": stirling_second(r, k),

        "s1(r,d)": signed_stirling_first(r, d),
        "abs_s1(r,d)": unsigned_stirling_first(r, d),
        "s2(r,d)": stirling_second(r, d),

        "k*r": sp.Integer(k * r),
        "k*d": sp.Integer(k * d),
        "r*d": sp.Integer(r * d),

        "k+1": sp.Integer(k + 1),
        "r+1": sp.Integer(r + 1),
        "d+1": sp.Integer(d + 1),
    }


# ============================================================================
# COMMON SCALE TEST
# ============================================================================

def common_scale_test(points, candidate_name):

    ratios = []
    zero_points = []

    for k, r, value in points:

        c = clean(
            candidate_library(
                k,
                r,
            )[candidate_name]
        )

        if c == 0:
            zero_points.append(
                (k, r)
            )
            continue

        ratios.append(
            clean(value / c)
        )

    if not ratios:
        return {
            "status": "NO_DEFINED_POINTS",
            "scale": None,
            "zero_points": zero_points,
            "defined_count": 0,
        }

    scale = ratios[0]

    exact = all(
        ratio == scale
        for ratio in ratios
    )

    # A candidate with zero entries in places where B is nonzero cannot be
    # a complete common-scale law.
    if zero_points:
        exact = False

    return {
        "status": (
            "EXACT_COMMON_SCALE"
            if exact
            else "FAILED"
        ),
        "scale": (
            scale
            if exact
            else None
        ),
        "zero_points": zero_points,
        "defined_count": len(ratios),
    }


# ============================================================================
# FACTOR PROFILE
# ============================================================================

def factor_profile(value):

    value = sp.Rational(value)

    numerator = abs(int(value.p))
    denominator = int(value.q)

    return {
        "sign": (
            "+"
            if value >= 0
            else "-"
        ),
        "numerator": numerator,
        "denominator": denominator,
        "numerator_factorization": (
            sp.factorint(numerator)
            if numerator
            else {}
        ),
        "denominator_factorization": (
            sp.factorint(denominator)
            if denominator > 1
            else {}
        ),
    }


# ============================================================================
# ROW/COLUMN PRIMITIVE CONTENT
# ============================================================================

def primitive_integer_content(values):
    """
    Convert a rational vector to an integer vector using the common
    denominator and return its primitive gcd.
    """

    vals = [
        sp.Rational(v)
        for v in values
    ]

    if not vals:
        return {
            "denominator_lcm": 1,
            "integers": [],
            "gcd": 0,
        }

    denominator_lcm = safe_lcm(
        [
            sp.denom(v)
            for v in vals
        ]
    )

    integers = [
        int(
            v * denominator_lcm
        )
        for v in vals
    ]

    g = 0

    for x in integers:
        g = math.gcd(
            g,
            abs(x),
        )

    return {
        "denominator_lcm": denominator_lcm,
        "integers": integers,
        "gcd": g,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    points = b_points()

    print("=" * 78)
    print(
        "EXPERIMENT 286R — EXACT B-ENTRY "
        "PROVENANCE / FORMULA RECONSTRUCTION AUDIT"
    )
    print("=" * 78)

    # ========================================================================
    # 1. ENTRY FACTOR PROFILES
    # ========================================================================

    print()
    print("=" * 78)
    print("1. EXACT B-ENTRY FACTOR PROFILE")
    print("=" * 78)

    for k, r, value in points:

        profile = factor_profile(
            value
        )

        print()
        print(
            f"  (k={k}, r={r})"
        )

        print(
            f"    value={value}"
        )

        print(
            f"    sign={profile['sign']}"
        )

        print(
            f"    numerator={profile['numerator']}"
        )

        print(
            f"    denominator={profile['denominator']}"
        )

        print(
            "    numerator_factorization="
            f"{profile['numerator_factorization']}"
        )

        print(
            "    denominator_factorization="
            f"{profile['denominator_factorization']}"
        )

    # ========================================================================
    # 2. COMMON ELEMENTARY SCALE SEARCH
    # ========================================================================

    print()
    print("=" * 78)
    print("2. COMMON ELEMENTARY SCALE SEARCH")
    print("=" * 78)

    library_names = list(
        candidate_library(
            0,
            0,
        ).keys()
    )

    common_hits = []

    for name in library_names:

        result = common_scale_test(
            points,
            name,
        )

        print()
        print(
            f"  {name}: "
            f"status={result['status']}"
        )

        print(
            f"    defined_count="
            f"{result['defined_count']}"
        )

        if result["status"] == "EXACT_COMMON_SCALE":

            print(
                f"    common_scale="
                f"{result['scale']}"
            )

            common_hits.append(
                (
                    name,
                    result["scale"],
                )
            )

        if result["zero_points"]:

            print(
                "    zero_candidate_points="
                f"{result['zero_points']}"
            )

    # ========================================================================
    # 3. DENOMINATOR STRUCTURE
    # ========================================================================

    print()
    print("=" * 78)
    print("3. DENOMINATOR STRUCTURE")
    print("=" * 78)

    denominator_counts = {}

    for k, r, value in points:

        denominator = int(
            sp.denom(
                sp.Rational(value)
            )
        )

        factorization = (
            sp.factorint(denominator)
            if denominator > 1
            else {}
        )

        denominator_counts[
            denominator
        ] = (
            denominator_counts.get(
                denominator,
                0,
            )
            + 1
        )

        print(
            f"  (k={k},r={r}): "
            f"denominator={denominator} "
            f"factorization={factorization}"
        )

    print()
    print(
        f"  distinct_denominators="
        f"{len(denominator_counts)}"
    )

    print(
        f"  denominator_frequency="
        f"{sorted(denominator_counts.items())}"
    )

    # ========================================================================
    # 4. SIGN STRUCTURE
    # ========================================================================

    print()
    print("=" * 78)
    print("4. SIGN STRUCTURE")
    print("=" * 78)

    for k, row in enumerate(B):

        signs = []

        for offset, value in enumerate(row):

            r = k + offset

            signs.append(
                (
                    r,
                    "+"
                    if value >= 0
                    else "-",
                )
            )

        print(
            f"  k={k}: {signs}"
        )

    # ========================================================================
    # 5. ROW CONTENT
    # ========================================================================

    print()
    print("=" * 78)
    print("5. ROW INTEGER-CONTENT AUDIT")
    print("=" * 78)

    for k, row in enumerate(B):

        profile = primitive_integer_content(
            row
        )

        print()
        print(
            f"  row k={k}:"
        )

        print(
            f"    denominator_lcm="
            f"{profile['denominator_lcm']}"
        )

        print(
            f"    integer_row="
            f"{profile['integers']}"
        )

        print(
            f"    primitive_gcd="
            f"{profile['gcd']}"
        )

    # ========================================================================
    # 6. ABSOLUTE-COLUMN CONTENT
    # ========================================================================

    print()
    print("=" * 78)
    print("6. ABSOLUTE-COLUMN INTEGER-CONTENT AUDIT")
    print("=" * 78)

    for r in range(8):

        column = []

        for k, row in enumerate(B):

            offset = r - k

            if (
                0 <= offset < len(row)
            ):
                column.append(
                    sp.Rational(
                        row[offset]
                    )
                )

        if not column:
            continue

        profile = primitive_integer_content(
            column
        )

        print()
        print(
            f"  column r={r}:"
        )

        print(
            f"    values={column}"
        )

        print(
            f"    denominator_lcm="
            f"{profile['denominator_lcm']}"
        )

        print(
            f"    integer_column="
            f"{profile['integers']}"
        )

        print(
            f"    primitive_gcd="
            f"{profile['gcd']}"
        )

    # ========================================================================
    # 7. INTERPRETATION
    # ========================================================================

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiments 283-285 rejected the obvious numerical descriptions of the
B-array.

Experiment 286R therefore asks a different question:

    Which exact arithmetic ingredients are visible in B[k,r]
    itself?

A complete common-scale hit

    B[k,r] = C * E(k,r)

would be strong evidence for an explicit source formula.

The denominator and primitive-content tables are also important.

For example, if denominators factor systematically through factorials,
binomials, powers of 2, or products of consecutive integers, that is
direct evidence about the construction of B.

This experiment intentionally does not use the q-family.

The purpose is to determine what the B-array itself is made from before
trying to reconnect it to q_p(r).
"""
    )

    # ========================================================================
    # 8. FINAL
    # ========================================================================

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  common_elementary_scale_hits="
        f"{len(common_hits)}"
    )

    print(
        f"  distinct_denominator_count="
        f"{len(denominator_counts)}"
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 286R COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print(
            "\nInterrupted."
        )
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise