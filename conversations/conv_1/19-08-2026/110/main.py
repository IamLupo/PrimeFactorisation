#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 276 — EXACT COMPLETE-ROW RESIDUAL / COLLECTIVE-OPERATOR AUDIT
==============================================================================

This is the corrected version of Experiment 275.

Critical SymPy fix:

    WRONG:
        sp.Rational(y + m, 2)

    RIGHT:
        (y + m) / 2

The residual pipeline is treated collectively:

    complete B-row
        ->
    polynomial P_k(j)
        ->
    exact division by row divisor
        ->
    centered residual
        ->
    odd projection
        ->
    observed C[k].

Individual B-basis terms are tested separately only as diagnostics.

An isolated basis term is allowed to fail exact division and is recorded
as UNDEFINED_ISOLATED. It is never converted to zero.

This experiment answers:

    1. Does the COMPLETE B-row reproduce the observed B-odd coefficients?
    2. Which individual B-basis terms are independently divisible?
    3. Is the residual extraction intrinsically collective?
"""

from __future__ import annotations

import sys
import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

j, y = sp.symbols("j y")


# ============================================================================
# SOURCE q-TABLES
# ============================================================================

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


# ============================================================================
# EXACT B-CHANNEL ROWS FROM EXPERIMENT 104
# ============================================================================

B_ROWS = [
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
        sp.Rational(51337, 30),
        sp.Rational(-52447, 180),
        sp.Rational(-14333, 210),
        sp.Rational(710501, 13440),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# ============================================================================
# OBSERVED ODD COEFFICIENTS [y^1, y^3, y^5, ...]
# ============================================================================

OBSERVED_C = [
    [
        sp.Rational(-584531, 35840),
        sp.Rational(-2908483, 1920),
        sp.Rational(-31233169, 13440),
        sp.Rational(-1446167, 1344),
        sp.Rational(-22259149, 40320),
        sp.Rational(-301, 240),
    ],
    [
        sp.Rational(-59257, 46080),
        sp.Rational(186547, 1440),
        sp.Rational(367433, 1680),
        sp.Rational(126549, 448),
        sp.Rational(-162139, 40320),
    ],
    [
        sp.Rational(4457, 46080),
        sp.Rational(-16819, 5760),
        sp.Rational(-5769, 896),
    ],
    [
        sp.Rational(-421, 322560),
    ],
    [
        # k=4
        sp.Rational(-22259149, 40320),
        sp.Rational(-162139, 40320),
        sp.Rational(0),
    ],
    [
        # k=5
        sp.Rational(-301, 240),
    ],
]


# ============================================================================
# HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
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


def as_poly(expr, variable):
    return sp.Poly(
        sp.expand(expr),
        variable,
        domain=sp.QQ,
    )


def exact_quotient(
    numerator,
    denominator,
):
    """
    Exact polynomial division over QQ.

    Returns:
        (True, quotient)
        (False, None)
    """

    numerator = clean(numerator)
    denominator = clean(denominator)

    if denominator == 0:
        return False, None

    quotient, remainder = sp.div(
        as_poly(numerator, j),
        as_poly(denominator, j),
        domain=sp.QQ,
    )

    if not remainder.is_zero:
        return False, None

    return True, clean(
        quotient.as_expr()
    )


def centered(
    expr,
    m=7,
):
    """
    j = (y+m)/2.

    IMPORTANT:
        use division, not sp.Rational(y+m, 2).
    """

    if expr is None:
        return None

    substitution = (
        y + sp.Integer(m)
    ) / sp.Integer(2)

    return clean(
        expr.subs(
            j,
            substitution,
        )
    )


def odd_part(expr):
    if expr is None:
        return None

    return clean(
        (
            expr
            - expr.subs(y, -y)
        )
        / sp.Integer(2)
    )


def coefficient_vector_odd(expr):
    """
    Return coefficients in ascending odd powers:

        [y^1, y^3, y^5, ...]
    """

    if expr is None:
        return None

    expr = clean(expr)

    if expr == 0:
        return []

    p = as_poly(
        expr,
        y,
    )

    degree = int(
        p.degree()
    )

    return [
        sp.Rational(
            p.nth(power)
        )
        for power in range(
            1,
            degree + 1,
            2,
        )
    ]


def row_polynomial(row):
    result = sp.Integer(0)

    for r, coeff in enumerate(row):
        result += (
            sp.sympify(coeff)
            * falling(
                j,
                r,
            )
        )

    return clean(result)


def row_divisor(k):
    """
    Exact residual divisor used in the Experiment-104 B channel.
    """

    return clean(
        falling(
            j,
            k,
        )
        * upper_falling(
            7,
            j,
            max(
                0,
                k - 4,
            ),
        )
    )


def complete_pipeline(
    row,
    k,
):
    P = row_polynomial(row)

    divisor = row_divisor(k)

    exact, residual = exact_quotient(
        P,
        divisor,
    )

    if not exact:
        return {
            "division_exact": False,
            "P": P,
            "divisor": divisor,
            "residual": None,
            "centered": None,
            "odd": None,
            "coefficients": None,
        }

    F = centered(
        residual,
        7,
    )

    O = odd_part(
        F
    )

    coeffs = coefficient_vector_odd(
        O
    )

    return {
        "division_exact": True,
        "P": P,
        "divisor": divisor,
        "residual": residual,
        "centered": F,
        "odd": O,
        "coefficients": coeffs,
    }


def isolated_pipeline(
    coeff,
    r,
    k,
):
    P = clean(
        sp.sympify(coeff)
        * falling(
            j,
            r,
        )
    )

    divisor = row_divisor(k)

    exact, residual = exact_quotient(
        P,
        divisor,
    )

    if not exact:
        return {
            "status": "UNDEFINED_ISOLATED",
            "P": P,
            "divisor": divisor,
            "residual": None,
            "centered": None,
            "odd": None,
            "coefficients": None,
        }

    F = centered(
        residual,
        7,
    )

    O = odd_part(
        F
    )

    return {
        "status": "DEFINED_ISOLATED",
        "P": P,
        "divisor": divisor,
        "residual": residual,
        "centered": F,
        "odd": O,
        "coefficients": coefficient_vector_odd(
            O
        ),
    }


def vector_equal(a, b):
    if a is None or b is None:
        return False

    if len(a) != len(b):
        return False

    return all(
        clean(x - y) == 0
        for x, y in zip(a, b)
    )


def pad_vector(v, width):
    if v is None:
        return None

    return list(v) + [
        sp.Integer(0)
        for _ in range(
            max(
                0,
                width - len(v),
            )
        )
    ]


def v17(n):
    """
    Exact 17-adic valuation for integer n.

    Returns None for zero.
    """

    n = int(n)

    if n == 0:
        return None

    n = abs(n)
    e = 0

    while n % 17 == 0:
        n //= 17
        e += 1

    return e


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 276 — EXACT COMPLETE-ROW RESIDUAL / "
        "COLLECTIVE-OPERATOR AUDIT"
    )
    print("=" * 78)

    complete_ok = True
    isolated_undefined = False

    complete_results = {}
    isolated_results = {}

    # ------------------------------------------------------------------------
    # 1. COMPLETE ROW PIPELINE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. COMPLETE B-ROW PIPELINE")
    print("=" * 78)

    for k, row in enumerate(
        B_ROWS
    ):

        result = complete_pipeline(
            row,
            k,
        )

        complete_results[k] = result

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    divisor={result['divisor']}"
        )

        print(
            f"    division_exact="
            f"{result['division_exact']}"
        )

        if not result[
            "division_exact"
        ]:
            complete_ok = False

            print(
                "    FATAL: complete row "
                "failed exact division."
            )

            continue

        derived = result[
            "coefficients"
        ]

        print(
            f"    residual="
            f"{result['residual']}"
        )

        print(
            f"    centered="
            f"{result['centered']}"
        )

        print(
            f"    odd_coefficients="
            f"{derived}"
        )

        if k < len(
            OBSERVED_C
        ):
            observed = OBSERVED_C[k]

            width = max(
                len(derived),
                len(observed),
            )

            a = pad_vector(
                derived,
                width,
            )

            b = pad_vector(
                observed,
                width,
            )

            exact = vector_equal(
                a,
                b,
            )

            print(
                f"    observed="
                f"{observed}"
            )

            print(
                f"    observed_match="
                f"{exact}"
            )

            complete_ok &= exact

    # ------------------------------------------------------------------------
    # 2. ISOLATED BASIS AUDIT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ISOLATED B-BASIS DIVISIBILITY AUDIT")
    print("=" * 78)

    for k, row in enumerate(
        B_ROWS
    ):

        isolated_results[k] = []

        print()
        print(
            f"  k={k}:"
        )

        defined = []
        undefined = []

        for r, coeff in enumerate(
            row
        ):

            result = isolated_pipeline(
                coeff,
                r,
                k,
            )

            isolated_results[k].append(
                result
            )

            if result[
                "status"
            ] == "DEFINED_ISOLATED":

                defined.append(r)

            else:

                undefined.append(r)
                isolated_undefined = True

            print(
                f"    r={r}: "
                f"{result['status']}"
            )

            if (
                result[
                    "status"
                ]
                == "DEFINED_ISOLATED"
            ):
                print(
                    f"      coefficients="
                    f"{result['coefficients']}"
                )

        print(
            f"    defined_indices={defined}"
        )

        print(
            f"    undefined_indices={undefined}"
        )

    # ------------------------------------------------------------------------
    # 3. COLLECTIVE-DIVISIBILITY PROFILE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. COLLECTIVE-DIVISIBILITY PROFILE")
    print("=" * 78)

    for k, row in enumerate(
        B_ROWS
    ):

        total_P = row_polynomial(
            row
        )

        divisor = row_divisor(
            k
        )

        full_exact, _ = exact_quotient(
            total_P,
            divisor,
        )

        individual_exact = 0

        for r, coeff in enumerate(
            row
        ):

            result = isolated_results[
                k
            ][r]

            if result[
                "status"
            ] == "DEFINED_ISOLATED":
                individual_exact += 1

        print()
        print(
            f"  k={k}: "
            f"complete_exact={full_exact} "
            f"individual_exact="
            f"{individual_exact}/{len(row)}"
        )

    # ------------------------------------------------------------------------
    # 4. TERMINAL SOURCE REFERENCE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. TERMINAL PROJECTIVE SOURCE")
    print("=" * 78)

    q1_terminal = Q[1][-1]
    q3_terminal = Q[3][-1]

    print(
        f"  q1_terminal={q1_terminal}"
    )

    print(
        f"  q3_terminal={q3_terminal}"
    )

    print(
        f"  gcd="
        f"{sp.gcd(
            q1_terminal,
            q3_terminal,
        )}"
    )

    print(
        f"  q1_v17={v17(q1_terminal)}"
    )

    print(
        f"  q3_v17={v17(q3_terminal)}"
    )

    print(
        f"  q1/17="
        f"{q1_terminal // 17}"
    )

    print(
        f"  q3/17="
        f"{q3_terminal // 17}"
    )

    # ------------------------------------------------------------------------
    # 5. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
The key distinction is:

    complete B-row divisibility
        versus
    isolated basis-term divisibility.

The complete row is the legitimate object of the residual construction.

An individual term

    B[k][r] * j_(r)

can fail exact division even when

    sum_r B[k][r] j_(r)

is exactly divisible.

Such a term is therefore marked

    UNDEFINED_ISOLATED

and is NOT treated as zero.

If all complete rows reconstruct their observed B-odd coefficients,
while many isolated terms are undefined, the residual operator is
collective.

That would explain why the earlier attempts to map individual q_p(r)
values directly to B-channel coefficients failed.

The unresolved source map is then precisely:

    q_p(r)
       ->
    complete B[k][r] row.

That is the level at which the original construction must next be
recovered.
"""
    )

    # ------------------------------------------------------------------------
    # 6. FINAL
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  complete_B_pipeline_exact="
        f"{complete_ok}"
    )

    print(
        f"  isolated_undefined_terms_present="
        f"{isolated_undefined}"
    )

    print(
        "  isolated_undefined_not_treated_as_zero=True"
    )

    print(
        "  residual_operator_collective_candidate="
        f"{complete_ok and isolated_undefined}"
    )

    print(
        "  q_to_B_operator_still_unresolved=True"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  failures="
        f"{0 if complete_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS="
        f"{complete_ok}"
    )

    print()
    print(
        "EXPERIMENT 276 COMPLETE"
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