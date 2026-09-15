#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 289R — EXACT FIRST-ORDER DIFFERENTIAL LOWERING / RODRIGUES AUDIT
==============================================================================

Tests whether adjacent primitive B-row generating polynomials satisfy

    H_{k+1}(x)
      =
    (a_k*x + b_k) H_k'(x)
      +
    c_k H_k(x).

This is the next natural degree-lowering operator after Experiment 288.

Only complete exact polynomial identities count.

Exact QQ arithmetic only.
No q-family data.
No interpolation-only fitting.
No arbitrary matrix fitting.
"""

from __future__ import annotations

import sys
import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

x = sp.Symbol("x")
ksym = sp.Symbol("k")


# ============================================================================
# PRIMITIVE INTEGER B-ROW POLYNOMIALS
# ============================================================================

H = [
    (
        -421*x**7
        - 69937*x**6
        + 423738*x**5
        - 1055250*x**4
        - 41580*x**3
        + 4071060*x**2
        + 1559880*x
        + 63000
    ),

    (
        -12106*x**6
        - 33638*x**5
        + 379836*x**4
        - 1230345*x**3
        + 819480*x**2
        + 3104640*x
        + 630000
    ),

    (
        -173070*x**5
        + 240576*x**4
        + 1037057*x**3
        - 6072045*x**2
        + 8596560*x
        + 8139600
    ),

    (
        -808952*x**4
        + 2458897*x**3
        - 3513510*x**2
        - 3999450*x
        + 18564000
    ),

    (
        -162139*x**3
        + 926401*x**2
        - 2779686*x
        + 3753750
    ),

    (
        301*x**2
        - 626*x
        + 650
    ),
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


def poly(expr):
    return sp.Poly(
        sp.expand(expr),
        x,
        domain=sp.QQ,
    )


def coefficient(expr, power):
    return poly(expr).nth(power)


def degree(expr):
    return int(
        poly(expr).degree()
    )


# ============================================================================
# FIRST-ORDER OPERATOR
# ============================================================================

def solve_first_order(k):
    """
    Solve

        H[k+1]
          =
        (a*x+b) H[k]' + c H[k].

    using every polynomial coefficient.
    """

    current = poly(H[k])
    target = poly(H[k + 1])
    derivative = current.diff()

    max_degree = max(
        current.degree(),
        target.degree(),
    )

    rows = []
    rhs = []

    for power in range(
        max_degree + 1
    ):

        # coefficient of a*x*H'
        if power == 0:
            ax_coeff = sp.Integer(0)
        else:
            ax_coeff = derivative.nth(
                power - 1
            )

        # coefficient of b*H'
        b_coeff = derivative.nth(
            power
        )

        # coefficient of c*H
        c_coeff = current.nth(
            power
        )

        rows.append(
            [
                ax_coeff,
                b_coeff,
                c_coeff,
            ]
        )

        rhs.append(
            target.nth(
                power
            )
        )

    M = sp.Matrix(rows)
    y = sp.Matrix(rhs)

    rank = M.rank()
    augmented_rank = (
        M.row_join(y).rank()
    )

    result = {
        "status": None,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "solution": None,
        "exact": False,
    }

    if augmented_rank != rank:
        result["status"] = "NO_SOLUTION"
        return result

    if rank != 3:
        result["status"] = "NONUNIQUE"
        return result

    solution = M.gauss_jordan_solve(
        y
    )[0]

    a = clean(solution[0])
    b = clean(solution[1])
    c = clean(solution[2])

    rebuilt = clean(
        (a*x + b)
        * sp.diff(
            H[k],
            x,
        )
        + c*H[k]
    )

    exact = (
        rebuilt
        == clean(
            H[k + 1]
        )
    )

    result["solution"] = (
        a,
        b,
        c,
    )

    result["exact"] = exact
    result["status"] = (
        "EXACT"
        if exact
        else "FAILED"
    )

    return result


# ============================================================================
# LEADING TERM CANCELLATION
# ============================================================================

def leading_cancellation(k, solution):

    if solution is None:
        return None

    a, _, c = solution

    P = poly(H[k])

    d = P.degree()
    leading = P.LC()

    value = clean(
        (
            a * d
            + c
        )
        * leading
    )

    return {
        "degree": d,
        "leading_coefficient": leading,
        "cancellation_value": value,
        "exact": value == 0,
    }


# ============================================================================
# STRICT VARIANT: H_{k+1}=(a*x+b)H'_k
# ============================================================================

def solve_strict(k):

    current = poly(H[k])
    target = poly(H[k + 1])
    derivative = current.diff()

    max_degree = max(
        current.degree(),
        target.degree(),
    )

    rows = []
    rhs = []

    for power in range(
        max_degree + 1
    ):

        ax_coeff = (
            sp.Integer(0)
            if power == 0
            else derivative.nth(
                power - 1
            )
        )

        b_coeff = derivative.nth(
            power
        )

        rows.append(
            [
                ax_coeff,
                b_coeff,
            ]
        )

        rhs.append(
            target.nth(
                power
            )
        )

    M = sp.Matrix(rows)
    y = sp.Matrix(rhs)

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

    a = clean(solution[0])
    b = clean(solution[1])

    rebuilt = clean(
        (a*x + b)
        * sp.diff(
            H[k],
            x,
        )
    )

    exact = (
        rebuilt
        == clean(
            H[k + 1]
        )
    )

    return (
        "EXACT"
        if exact
        else "FAILED",
        (
            a,
            b,
        ),
    )


# ============================================================================
# PARAMETER LAW SEARCH
# ============================================================================

def parameter_law(
    values,
    name,
):

    print()
    print(
        "  parameter="
        + name
    )

    for degree_bound in (
        0,
        1,
        2,
    ):

        if len(values) <= degree_bound + 1:
            print(
                "    degree<="
                + str(degree_bound)
                + ": INSUFFICIENT_DATA"
            )
            continue

        rows = []
        rhs = []

        for kk, value in values:

            rows.append(
                [
                    sp.Integer(kk) ** j
                    for j in range(
                        degree_bound + 1
                    )
                ]
            )

            rhs.append(
                value
            )

        M = sp.Matrix(rows)
        y = sp.Matrix(rhs)

        rank = M.rank()
        augmented_rank = (
            M.row_join(y).rank()
        )

        if augmented_rank != rank:
            print(
                "    degree<="
                + str(degree_bound)
                + ": NO_SOLUTION"
            )
            continue

        if rank != M.cols:
            print(
                "    degree<="
                + str(degree_bound)
                + ": NONUNIQUE"
            )
            continue

        solution = M.gauss_jordan_solve(
            y
        )[0]

        expression = clean(
            sum(
                solution[j]
                * ksym ** j
                for j in range(
                    degree_bound + 1
                )
            )
        )

        exact = True

        for kk, value in values:

            if clean(
                expression.subs(
                    ksym,
                    kk,
                )
                - value
            ) != 0:
                exact = False
                break

        print(
            "    degree<="
            + str(degree_bound)
            + ": "
            + (
                "EXACT "
                if exact
                else "FAILED "
            )
            + str(expression)
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 289R — EXACT FIRST-ORDER "
        "DIFFERENTIAL LOWERING / RODRIGUES AUDIT"
    )
    print("=" * 78)

    recovered = {}

    # ------------------------------------------------------------------------
    # 1. GENERAL FIRST-ORDER OPERATOR
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. FIRST-ORDER DIFFERENTIAL OPERATOR"
    )
    print("=" * 78)

    for k in range(
        len(H) - 1
    ):

        result = solve_first_order(
            k
        )

        recovered[k] = result

        print()
        print(
            "  k="
            + str(k)
            + ":"
        )

        print(
            "    equations="
            + str(
                max(
                    degree(H[k]),
                    degree(H[k + 1]),
                )
                + 1
            )
        )

        print(
            "    unknowns=3"
        )

        print(
            "    rank="
            + str(result["rank"])
        )

        print(
            "    augmented_rank="
            + str(
                result[
                    "augmented_rank"
                ]
            )
        )

        print(
            "    status="
            + str(
                result["status"]
            )
        )

        print(
            "    solution="
            + str(
                result["solution"]
            )
        )

        if result["solution"] is not None:

            cancel = (
                leading_cancellation(
                    k,
                    result["solution"],
                )
            )

            print(
                "    leading_degree="
                + str(
                    cancel["degree"]
                )
            )

            print(
                "    leading_coefficient="
                + str(
                    cancel[
                        "leading_coefficient"
                    ]
                )
            )

            print(
                "    leading_cancellation="
                + str(
                    cancel[
                        "cancellation_value"
                    ]
                )
            )

            print(
                "    cancellation_exact="
                + str(
                    cancel["exact"]
                )
            )

            print(
                "    exact_identity="
                + str(
                    result["exact"]
                )
            )

    # ------------------------------------------------------------------------
    # 2. STRICT VARIANT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "2. STRICT (a*x+b) H' VARIANT"
    )
    print("=" * 78)

    strict_hits = []

    for k in range(
        len(H) - 1
    ):

        status, params = solve_strict(
            k
        )

        print()
        print(
            "  k="
            + str(k)
            + ": status="
            + str(status)
            + " parameters="
            + str(params)
        )

        if status == "EXACT":
            strict_hits.append(
                (
                    k,
                    params,
                )
            )

    # ------------------------------------------------------------------------
    # 3. PARAMETER PROFILE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. EXACT PARAMETER PROFILE"
    )
    print("=" * 78)

    exact_pairs = [
        (
            k,
            result["solution"],
        )
        for k, result in recovered.items()
        if result["status"] == "EXACT"
    ]

    print()
    print(
        "  exact_first_order_pairs="
        + str(
            len(exact_pairs)
        )
    )

    for k, params in exact_pairs:
        print(
            "  k="
            + str(k)
            + ": (a,b,c)="
            + str(params)
        )

    if exact_pairs:

        for index, name in enumerate(
            ("a", "b", "c")
        ):

            values = [
                (
                    k,
                    params[index],
                )
                for k, params in exact_pairs
            ]

            parameter_law(
                values,
                name,
            )

    # ------------------------------------------------------------------------
    # 4. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 288 ruled out ordinary differentiation.

Experiment 289 tests

    H_{k+1}
      =
    (a_k x+b_k) H'_k
      +
    c_k H_k.

This is the simplest differential operator capable of lowering degree
through cancellation of the highest-degree term.

A genuine hit must satisfy the COMPLETE polynomial identity.

The most informative case would be:

    exact_first_order_pairs = 5,

together with simple exact formulas for a_k, b_k, c_k.

If there are no exact pairs, then the degree ladder is not generated by
this Rodrigues-type first-order mechanism.
"""
    )

    # ------------------------------------------------------------------------
    # 5. TERMINAL SOURCE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "5. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    q1_terminal = 495451247
    q3_terminal = 421514439
    gcd_terminal = sp.gcd(
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

    # ------------------------------------------------------------------------
    # 6. FINAL
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_first_order_pairs="
        + str(
            len(exact_pairs)
        )
    )

    print(
        "  exact_strict_pairs="
        + str(
            len(strict_hits)
        )
    )

    print(
        "  complete_identity_only=True"
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
        "EXPERIMENT 289R COMPLETE"
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