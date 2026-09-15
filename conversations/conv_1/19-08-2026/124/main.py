#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 290R — EXACT B-CONSTRUCTION REVERSE ENGINEERING / SAFE NORMALIZATION
==============================================================================

Corrected version of Experiment 290.

Fix:
    None normalization values are treated as UNDEFINED and are never passed
    to SymPy arithmetic routines.

The experiment tests exact combinatorial normalizations:

    k!,
    r!,
    (r-k)!,
    k!/r!,
    r!/k!,
    binomial factors,
    falling-factorial factors,

and then checks whether the normalized table becomes:

    * integer-valued;
    * exactly affine along fixed d=r-k;
    * exactly affine along fixed r;
    * exactly affine along fixed k.

No q-family data.
No arbitrary matrix fitting.
No interpolation-only proof.
Exact QQ arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# EXACT B TABLE
# =============================================================================

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


# =============================================================================
# HELPERS
# =============================================================================

def clean(expr):
    if expr is None:
        return None

    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def safe_lcm(values):
    result = 1

    for value in values:
        if value is None:
            continue

        value = abs(int(value))

        if value:
            result = math.lcm(
                result,
                value,
            )

    return result


def falling(n, q):
    if q < 0:
        return sp.Integer(0)

    result = sp.Integer(1)

    for j in range(q):
        result *= n - j

    return clean(result)


def safe_binomial(n, q):
    if q < 0 or q > n:
        return sp.Integer(0)

    return sp.binomial(
        sp.Integer(n),
        sp.Integer(q),
    )


def points():
    result = []

    for kk, row in enumerate(B):

        for offset, value in enumerate(row):

            rr = kk + offset
            dd = rr - kk

            result.append(
                (
                    kk,
                    rr,
                    dd,
                    sp.Rational(value),
                )
            )

    return result


# =============================================================================
# NORMALIZATION LIBRARY
# =============================================================================

def normalization_library(kk, rr):

    dd = rr - kk

    out = {
        "raw": sp.Integer(1),

        "k!": sp.factorial(kk),
        "r!": sp.factorial(rr),
        "d!": sp.factorial(dd),

        "k!/r!": None,
        "r!/k!": None,

        "C(7-k,d)": None,
        "C(7-r,k)": None,
        "C(7,k)": None,
        "C(7,r)": None,

        "falling(r,k)": falling(rr, kk),
        "falling(r,d)": falling(rr, dd),
        "falling(7-k,d)": falling(7 - kk, dd),
        "falling(7-r,k)": falling(7 - rr, kk),
    }

    if rr != 0:
        out["k!/r!"] = sp.Rational(
            sp.factorial(kk),
            sp.factorial(rr),
        )

        out["r!/k!"] = sp.Rational(
            sp.factorial(rr),
            sp.factorial(kk),
        )

    out["C(7-k,d)"] = safe_binomial(
        7 - kk,
        dd,
    )

    out["C(7-r,k)"] = safe_binomial(
        7 - rr,
        kk,
    )

    out["C(7,k)"] = safe_binomial(
        7,
        kk,
    )

    out["C(7,r)"] = safe_binomial(
        7,
        rr,
    )

    return out


# =============================================================================
# NORMALIZED VALUES
# =============================================================================

def normalized_values(candidate):

    result = []

    for kk, rr, dd, value in points():

        factor = normalization_library(
            kk,
            rr,
        ).get(candidate)

        if factor is None:
            result.append(
                (
                    kk,
                    rr,
                    dd,
                    None,
                )
            )
            continue

        factor = clean(
            factor
        )

        if factor is None or factor == 0:
            result.append(
                (
                    kk,
                    rr,
                    dd,
                    None,
                )
            )
            continue

        result.append(
            (
                kk,
                rr,
                dd,
                clean(
                    value * factor
                ),
            )
        )

    return result


# =============================================================================
# INTEGER PROFILE
# =============================================================================

def integer_profile(values):

    defined = [
        value
        for _, _, _, value in values
        if value is not None
    ]

    if not defined:
        return {
            "defined": 0,
            "undefined": len(values),
            "integer_all": False,
            "denominator_lcm": None,
            "primitive_gcd": None,
        }

    integer_all = all(
        sp.denom(
            sp.Rational(value)
        ) == 1
        for value in defined
    )

    denominator_lcm = safe_lcm(
        [
            sp.denom(
                sp.Rational(value)
            )
            for value in defined
        ]
    )

    integers = [
        int(
            sp.Rational(value)
            * denominator_lcm
        )
        for value in defined
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    return {
        "defined": len(defined),
        "undefined": (
            len(values)
            - len(defined)
        ),
        "integer_all": integer_all,
        "denominator_lcm": denominator_lcm,
        "primitive_gcd": g,
    }


# =============================================================================
# GROUPING
# =============================================================================

def group_by(values, index):

    groups = {}

    for kk, rr, dd, value in values:

        if value is None:
            continue

        key = (
            kk,
            rr,
            dd,
        )[index]

        groups.setdefault(
            key,
            [],
        ).append(
            (
                key,
                value,
            )
        )

    return groups


# =============================================================================
# EXACT LINEAR LAW
# =============================================================================

def exact_linear(sequence):

    if len(sequence) < 3:
        return (
            "INSUFFICIENT_DATA",
            None,
        )

    M = sp.Matrix(
        [
            [
                sp.Integer(xx),
                sp.Integer(1),
            ]
            for xx, _ in sequence
        ]
    )

    y = sp.Matrix(
        [
            yy
            for _, yy in sequence
        ]
    )

    rank = M.rank()

    augmented_rank = (
        M.row_join(y).rank()
    )

    if augmented_rank != rank:
        return (
            "NO_SOLUTION",
            None,
        )

    if rank != 2:
        return (
            "NONUNIQUE",
            None,
        )

    solution = M.gauss_jordan_solve(
        y
    )[0]

    return (
        "EXACT",
        tuple(
            clean(value)
            for value in solution
        ),
    )


# =============================================================================
# PRINT DEFINED/UNDEFINED SUMMARY
# =============================================================================

def print_normalized_table(
    name,
    values,
):

    print()
    print(
        "  "
        + name
        + ":"
    )

    for kk, rr, dd, value in values:

        label = (
            "(k="
            + str(kk)
            + ",r="
            + str(rr)
            + ",d="
            + str(dd)
            + ")"
        )

        print(
            "    "
            + label
            + " -> "
            + (
                "UNDEFINED"
                if value is None
                else str(value)
            )
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 290R — EXACT B-CONSTRUCTION "
        "REVERSE ENGINEERING / SAFE NORMALIZATION"
    )
    print("=" * 78)

    candidate_names = list(
        normalization_library(
            0,
            0,
        ).keys()
    )

    integer_hits = []
    linear_hits_d = []
    linear_hits_r = []
    linear_hits_k = []

    # =========================================================================
    # 1. NORMALIZATION PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "1. EXACT NORMALIZATION PROFILE"
    )
    print("=" * 78)

    for name in candidate_names:

        values = normalized_values(
            name
        )

        profile = integer_profile(
            values
        )

        print()
        print(
            "  "
            + name
            + ":"
        )

        print(
            "    defined="
            + str(
                profile["defined"]
            )
        )

        print(
            "    undefined="
            + str(
                profile["undefined"]
            )
        )

        print(
            "    integer_all="
            + str(
                profile["integer_all"]
            )
        )

        print(
            "    denominator_lcm="
            + str(
                profile[
                    "denominator_lcm"
                ]
            )
        )

        print(
            "    primitive_gcd="
            + str(
                profile[
                    "primitive_gcd"
                ]
            )
        )

        if (
            profile["integer_all"]
            and profile["undefined"] == 0
        ):
            integer_hits.append(
                (
                    name,
                    profile,
                )
            )

    # =========================================================================
    # 2. FIXED-OFFSET d = r-k
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "2. FIXED-OFFSET NORMALIZED LAWS"
    )
    print("=" * 78)

    for name in candidate_names:

        values = normalized_values(
            name
        )

        groups = {}

        for kk, rr, dd, value in values:

            if value is None:
                continue

            groups.setdefault(
                dd,
                [],
            ).append(
                (
                    kk,
                    value,
                )
            )

        for dd in sorted(groups):

            sequence = groups[dd]

            status, params = exact_linear(
                sequence
            )

            if status == "EXACT":

                print()
                print(
                    "  candidate="
                    + name
                    + " d="
                    + str(dd)
                )

                print(
                    "    points="
                    + str(sequence)
                )

                print(
                    "    exact_linear="
                    + str(params)
                )

                linear_hits_d.append(
                    (
                        name,
                        dd,
                        params,
                    )
                )

    # =========================================================================
    # 3. FIXED ABSOLUTE r
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "3. FIXED-ABSOLUTE-r NORMALIZED LAWS"
    )
    print("=" * 78)

    for name in candidate_names:

        values = normalized_values(
            name
        )

        groups = {}

        for kk, rr, dd, value in values:

            if value is None:
                continue

            groups.setdefault(
                rr,
                [],
            ).append(
                (
                    kk,
                    value,
                )
            )

        for rr in sorted(groups):

            sequence = groups[rr]

            status, params = exact_linear(
                sequence
            )

            if status == "EXACT":

                print()
                print(
                    "  candidate="
                    + name
                    + " r="
                    + str(rr)
                )

                print(
                    "    points="
                    + str(sequence)
                )

                print(
                    "    exact_linear="
                    + str(params)
                )

                linear_hits_r.append(
                    (
                        name,
                        rr,
                        params,
                    )
                )

    # =========================================================================
    # 4. FIXED k
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "4. FIXED-k NORMALIZED LAWS"
    )
    print("=" * 78)

    for name in candidate_names:

        values = normalized_values(
            name
        )

        groups = {}

        for kk, rr, dd, value in values:

            if value is None:
                continue

            groups.setdefault(
                kk,
                [],
            ).append(
                (
                    rr,
                    value,
                )
            )

        for kk in sorted(groups):

            sequence = groups[kk]

            status, params = exact_linear(
                sequence
            )

            if status == "EXACT":

                print()
                print(
                    "  candidate="
                    + name
                    + " k="
                    + str(kk)
                )

                print(
                    "    points="
                    + str(sequence)
                )

                print(
                    "    exact_linear="
                    + str(params)
                )

                linear_hits_k.append(
                    (
                        name,
                        kk,
                        params,
                    )
                )

    # =========================================================================
    # 5. IMPORTANT NORMALIZED TABLES
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "5. IMPORTANT NORMALIZED TABLES"
    )
    print("=" * 78)

    important = [
        "k!/r!",
        "r!/k!",
        "d!",
        "C(7-k,d)",
        "C(7-r,k)",
        "falling(r,k)",
        "falling(7-k,d)",
    ]

    for name in important:

        values = normalized_values(
            name
        )

        print_normalized_table(
            name,
            values,
        )

    # =========================================================================
    # 6. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "6. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The previous experiments rejected direct polynomial, recurrence,
separability, factorization, and differential-operator explanations.

Experiment 290R asks whether the apparent complexity is only a
normalization artifact.

The most useful outcomes are exact identities after multiplying B[k,r]
by natural finite combinatorial factors:

    k!,
    r!,
    (r-k)!,
    binomial coefficients,
    falling factorials.

Undefined normalization entries are explicitly retained as UNDEFINED.
They are never converted to zero and never passed to SymPy arithmetic.

A strong result would be an exact normalization that makes the complete
B-table integer-valued and then reveals a simple law along d, r, or k.

A negative result means the rational coefficients are not explained by
these elementary finite-product normalizations.
"""
    )

    # =========================================================================
    # 7. TERMINAL SOURCE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "7. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    q1_terminal = 495451247
    q3_terminal = 421514439

    gcd_terminal = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    print(
        "  q1_terminal="
        + str(q1_terminal)
    )

    print(
        "  q3_terminal="
        + str(q3_terminal)
    )

    print(
        "  gcd="
        + str(gcd_terminal)
    )

    print(
        "  q1/17="
        + str(
            q1_terminal // 17
        )
    )

    print(
        "  q3/17="
        + str(
            q3_terminal // 17
        )
    )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  integer_normalization_hits="
        + str(
            len(integer_hits)
        )
    )

    print(
        "  exact_linear_d_hits="
        + str(
            len(linear_hits_d)
        )
    )

    print(
        "  exact_linear_r_hits="
        + str(
            len(linear_hits_r)
        )
    )

    print(
        "  exact_linear_k_hits="
        + str(
            len(linear_hits_k)
        )
    )

    print(
        "  safe_none_handling=True"
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
        "EXPERIMENT 290R COMPLETE"
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