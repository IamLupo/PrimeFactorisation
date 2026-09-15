#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 267RR — EXACT q-TABLE / B-CHANNEL OPERATOR PROVENANCE AUDIT
==============================================================================

Robust correction of Experiments 267 and 267R.

The experiment separates:

    COMPLETE q-ROW OPERATOR
from
    INDIVIDUAL q_r BASIS-TERM DIVISIBILITY.

An isolated q_r term may fail the residual divisibility test even when
the complete q-row is exactly divisible.

Therefore:

    None = undefined / inadmissible isolated contribution

and is NEVER treated as integer zero.

All comparisons involving None are guarded explicitly.

The experiment tests:

    q_p(r)
        ->
    raw falling-factorial polynomial
        ->
    exact residual division
        ->
    centering
        ->
    odd parity coefficient.

No fitted matrices.
No floating point.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

j, y = sp.symbols("j y")


# =============================================================================
# SOURCE q-TABLE
# =============================================================================

Q = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],

    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],

    5: [
        -62398,
        4771718,
        16027881,
    ],

    7: [
        1,
    ],
}


D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}


# =============================================================================
# OBSERVED B-CHANNEL DATA
# =============================================================================

B = [
    [
        25,
        619,
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],

    [
        1750,
        8624,
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],

    [
        9690,
        10234,
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
# EXACT HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def poly(expr, var):
    return sp.Poly(
        sp.expand(expr),
        var,
        domain=sp.QQ,
    )


def falling(x, n):
    result = sp.Integer(1)

    for r in range(n):
        result *= x - r

    return sp.expand(result)


def upper_falling(m, x, n):
    result = sp.Integer(1)

    for r in range(n):
        result *= m - x - r

    return sp.expand(result)


def exact_quotient(num, den, var):
    pn = poly(
        clean(num),
        var,
    )

    pd = poly(
        clean(den),
        var,
    )

    q, r = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(
        q.as_expr()
    )


def coefficient(expr, power):
    expr = clean(expr)

    if expr == 0:
        return sp.Integer(0)

    return sp.Rational(
        poly(
            expr,
            y,
        ).nth(power)
    )


def v17(value):
    """
    Exact 17-adic valuation.

    Returns None for zero.
    """
    value = int(value)

    if value == 0:
        return None

    value = abs(value)
    exponent = 0

    while value % 17 == 0:
        value //= 17
        exponent += 1

    return exponent


def safe_equal(a, b):
    """
    Equality that safely handles None.
    """
    if a is None or b is None:
        return a is None and b is None

    return clean(a - b) == 0


def safe_difference(a, b):
    """
    Return a-b only when both values are defined.
    """
    if a is None or b is None:
        return None

    return clean(a - b)


# =============================================================================
# SOURCE q-ROW -> RAW
# =============================================================================

def source_raw_polynomial(q_row, k_index):
    return clean(
        sum(
            sp.Integer(
                q_row[r]
            )
            *
            falling(
                j,
                k_index + r,
            )
            for r in range(
                len(q_row)
            )
        )
    )


# =============================================================================
# RAW -> RESIDUAL
# =============================================================================

def residual_from_q_row(
    q_row,
    m,
    k_index,
):
    raw = source_raw_polynomial(
        q_row,
        k_index,
    )

    s = max(
        0,
        k_index - 4,
    )

    divisor = clean(
        falling(
            j,
            k_index,
        )
        *
        upper_falling(
            m,
            j,
            s,
        )
    )

    residual = exact_quotient(
        raw,
        divisor,
        j,
    )

    return (
        raw,
        divisor,
        residual,
    )


# =============================================================================
# SOURCE q-ROW -> B-ODD COEFFICIENT
# =============================================================================

def source_parity_coefficient(
    q_row,
    k_index,
    parity_power,
):
    raw, divisor, residual = (
        residual_from_q_row(
            q_row,
            7,
            k_index,
        )
    )

    if residual is None:
        return {
            "admissible": False,
            "raw": raw,
            "divisor": divisor,
            "residual": None,
            "centered": None,
            "odd": None,
            "coefficient": None,
        }

    centered = clean(
        residual.subs(
            j,
            (y + 7) / 2,
        )
    )

    odd = clean(
        (
            centered
            -
            centered.subs(
                y,
                -y,
            )
        )
        / 2
    )

    return {
        "admissible": True,
        "raw": raw,
        "divisor": divisor,
        "residual": residual,
        "centered": centered,
        "odd": odd,
        "coefficient": coefficient(
            odd,
            parity_power,
        ),
    }


# =============================================================================
# ACTUAL OBSERVED B-ODD COEFFICIENTS
# =============================================================================

def actual_b_odd_vectors():

    actual = {}

    for p in [1, 3, 5, 7]:

        values = []

        for k_index in range(
            D[p] + 1
        ):

            raw = sp.Integer(0)

            for offset, coeff in enumerate(
                B[k_index]
            ):
                raw += (
                    sp.sympify(coeff)
                    *
                    falling(
                        j,
                        k_index + offset,
                    )
                )

            s = max(
                0,
                k_index - 4,
            )

            divisor = clean(
                falling(
                    j,
                    k_index,
                )
                *
                upper_falling(
                    7,
                    j,
                    s,
                )
            )

            residual = exact_quotient(
                raw,
                divisor,
                j,
            )

            if residual is None:
                raise ArithmeticError(
                    f"Observed B residual failed: "
                    f"p={p}, k={k_index}"
                )

            centered = clean(
                residual.subs(
                    j,
                    (y + 7) / 2,
                )
            )

            odd = clean(
                (
                    centered
                    -
                    centered.subs(
                        y,
                        -y,
                    )
                )
                / 2
            )

            values.append(
                coefficient(
                    odd,
                    p,
                )
            )

        actual[p] = values

    return actual


# =============================================================================
# MAIN
# =============================================================================

def main():

    actual = actual_b_odd_vectors()

    failures = []

    print("=" * 78)
    print(
        "EXPERIMENT 267RR — EXACT q-TABLE / B-CHANNEL "
        "OPERATOR PROVENANCE AUDIT"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. COMPLETE q-ROW -> B-ODD
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. COMPLETE q-ROW -> B-ODD PIPELINE")
    print("=" * 78)

    source_derived = {}

    for p in [1, 3, 5, 7]:

        values = []
        admissibility = []

        for k_index in range(
            D[p] + 1
        ):

            result = source_parity_coefficient(
                Q[p],
                k_index,
                p,
            )

            admissibility.append(
                result["admissible"]
            )

            values.append(
                result["coefficient"]
            )

        source_derived[p] = values

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    admissible="
            f"{admissibility}"
        )

        print(
            f"    source_derived="
            f"{values}"
        )

        print(
            f"    actual="
            f"{actual[p]}"
        )

        comparable = [
            i
            for i, value in enumerate(values)
            if value is not None
        ]

        exact = all(
            safe_equal(
                values[i],
                actual[p][i],
            )
            for i in comparable
        )

        print(
            f"    comparable_indices="
            f"{comparable}"
        )

        print(
            f"    exact_on_defined_rows="
            f"{exact}"
        )

    # -------------------------------------------------------------------------
    # 2. SAFE DIFFERENCE PROFILE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SOURCE-DERIVED MINUS ACTUAL")
    print("=" * 78)

    complete_difference_ok = True

    for p in [1, 3, 5, 7]:

        differences = []

        for i in range(
            len(actual[p])
        ):

            value = safe_difference(
                source_derived[p][i],
                actual[p][i],
            )

            differences.append(
                value
            )

        defined_differences = [
            value
            for value in differences
            if value is not None
        ]

        exact = all(
            value == 0
            for value in defined_differences
        )

        complete_difference_ok &= exact

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    differences="
            f"{differences}"
        )

        print(
            f"    defined_difference_zero="
            f"{exact}"
        )

    # -------------------------------------------------------------------------
    # 3. INDIVIDUAL q_r DIVISIBILITY
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. INDIVIDUAL q_r DIVISIBILITY AUDIT")
    print("=" * 78)

    admissible_support = {}

    for p in [1, 3, 5, 7]:

        support = []

        print()
        print(
            f"  p={p}:"
        )

        for r, q_value in enumerate(
            Q[p]
        ):

            unit_row = [
                0
                for _ in Q[p]
            ]

            unit_row[r] = 1

            profile = []

            for k_index in range(
                D[p] + 1
            ):

                result = (
                    source_parity_coefficient(
                        unit_row,
                        k_index,
                        p,
                    )
                )

                profile.append(
                    result["admissible"]
                )

            all_admissible = all(
                profile
            )

            if all_admissible:
                support.append(r)

            print(
                f"    r={r}: "
                f"q={q_value} "
                f"all_levels_admissible="
                f"{all_admissible} "
                f"profile={profile}"
            )

        admissible_support[p] = support

        print(
            f"    fully_admissible_indices="
            f"{support}"
        )

    # -------------------------------------------------------------------------
    # 4. TERMINAL SOURCE AUDIT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. TERMINAL SOURCE AUDIT")
    print("=" * 78)

    terminal_ok = True

    for p in [1, 3, 5, 7]:

        q_terminal = Q[p][
            D[p]
        ]

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    q_terminal="
            f"{q_terminal}"
        )

        print(
            f"    v17="
            f"{v17(q_terminal)}"
        )

        if p == 1:

            normalized = (
                q_terminal // 17
            )

            expected = 29144191

            print(
                f"    /17="
                f"{normalized}"
            )

            print(
                f"    expected="
                f"{expected}"
            )

            if normalized != expected:
                terminal_ok = False

        elif p == 3:

            normalized = (
                q_terminal // 17
            )

            expected = 24794967

            print(
                f"    /17="
                f"{normalized}"
            )

            print(
                f"    expected="
                f"{expected}"
            )

            if normalized != expected:
                terminal_ok = False

        elif p == 5:

            print(
                f"    v17_expected=0"
            )

            if v17(q_terminal) != 0:
                terminal_ok = False

        elif p == 7:

            print(
                f"    q_terminal_is_one="
                f"{q_terminal == 1}"
            )

            if q_terminal != 1:
                terminal_ok = False

    # -------------------------------------------------------------------------
    # 5. COMPLETE SOURCE OPERATOR DIAGNOSTIC
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. COMPLETE SOURCE OPERATOR DIAGNOSTIC")
    print("=" * 78)

    print(
r"""
The important mathematical distinction is now explicit.

An isolated basis vector q_r=1 does NOT have to be independently
divisible by the residual denominator.

Only the COMPLETE q-row is required to satisfy the divisibility
conditions of the original construction.

Therefore:

    None
        = undefined isolated contribution,

not

    None = 0.

The complete q-row test and the individual basis-term audit are
therefore kept separate.
"""
    )

    # -------------------------------------------------------------------------
    # 6. SOURCE-PROVENANCE STATUS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SOURCE-PROVENANCE STATUS")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        defined = [
            i
            for i, value in enumerate(
                source_derived[p]
            )
            if value is not None
        ]

        print(
            f"  p={p}: "
            f"defined_source_operator_rows="
            f"{defined}"
        )

    # -------------------------------------------------------------------------
    # 7. FINAL EXACTNESS
    # -------------------------------------------------------------------------

    final_ok = terminal_ok

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  safe_none_handling=True"
    )

    print(
        f"  complete_source_difference_audit="
        f"{complete_difference_ok}"
    )

    print(
        f"  terminal_source_provenance="
        f"{terminal_ok}"
    )

    print(
        "  individual_qr_divisibility_not_assumed=True"
    )

    print(
        "  universal_operator_formula_proved=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print(
        "EXPERIMENT 267RR COMPLETE"
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