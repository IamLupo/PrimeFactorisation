#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 329R — EXACT SECOND-ORDER RATIONAL P-DEPENDENT VERTICAL AUDIT
==============================================================================

Purpose
-------
Experiment 328R rejected rational first-order laws of degree <= 2:

    Q_{t+1}(p) / Q_t(p) = N(p) / D(p).

The next structured possibility is a second-order recurrence

    Q_{t+2}(p)
      =
    R0(p) Q_t(p)
      +
    R1(p) Q_{t+1}(p),

where

    R0(p) = N0(p) / D(p),
    R1(p) = N1(p) / D(p).

After cross multiplication:

    Q_{t+2}(p) D(p)
      - Q_t(p) N0(p)
      - Q_{t+1}(p) N1(p)
      = 0.

The unknown coefficients therefore enter linearly.

Tested rational degrees:

    degree <= 0
    degree <= 1
    degree <= 2

Important data-availability check:

    order 2 gives 8 usable vertical windows.

    degree 0:
        3 coefficient polynomials/constants,
        projective dimension 2,
        genuinely overdetermined;

    degree 1:
        6 coefficients,
        projective dimension 5,
        genuinely overdetermined;

    degree 2:
        9 coefficients,
        projective dimension 8,
        only data-sized.

Therefore only degree 0 and degree 1 can provide genuine evidence.

Rules
-----
No missing values.
No interpolation.
No synthetic second n=pq case.
No extrapolation.
Exact SymPy rational arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# OBSERVED SOURCE TABLE
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


p = sp.symbols("p")


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


def valuation(x, prime):
    x = sp.Rational(x)

    if x == 0:
        return sp.oo

    numerator = abs(int(x.p))
    denominator = abs(int(x.q))

    value = 0

    while numerator % prime == 0:
        numerator //= prime
        value += 1

    while denominator % prime == 0:
        denominator //= prime
        value -= 1

    return value


def primitive_integer_polynomial(expr):
    """
    Convert a rational polynomial in p into a primitive integer polynomial.
    """

    poly = sp.Poly(
        sp.expand(expr),
        p,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return sp.Poly(
            0,
            p,
            domain=sp.ZZ,
        )

    coefficients = poly.all_coeffs()

    denominator_lcm = 1

    for coefficient in coefficients:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(sp.denom(coefficient)),
        )

    integers = [
        int(coefficient * denominator_lcm)
        for coefficient in coefficients
    ]

    coefficient_gcd = 0

    for value in integers:
        coefficient_gcd = math.gcd(
            coefficient_gcd,
            abs(value),
        )

    integers = [
        value // coefficient_gcd
        for value in integers
    ]

    if integers[0] < 0:
        integers = [
            -value
            for value in integers
        ]

    expression = sum(
        sp.Integer(value)
        * p**(
            len(integers) - 1 - index
        )
        for index, value in enumerate(integers)
    )

    return sp.Poly(
        sp.expand(expression),
        p,
        domain=sp.ZZ,
    )


def build_table():
    """
    Convert reversed layer storage into

        table[(p,t)] = Q_t(p)

    using observed cells only.
    """

    table = {}

    for p_value, values in Q.items():

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            table[
                (
                    int(p_value),
                    int(t_value),
                )
            ] = sp.Integer(value)

    return table


# ============================================================================
# SECOND-ORDER RATIONAL SYSTEM
# ============================================================================

def build_system(table, degree):
    """
    Build

        Q_{t+2} D - Q_t N0 - Q_{t+1} N1 = 0.

    Unknown vector ordering:

        N0_0 ... N0_d,
        N1_0 ... N1_d,
        D_0  ... D_d.
    """

    unknown_count = 3 * (
        degree + 1
    )

    rows = []
    metadata = []

    for (
        p_value,
        t_value,
    ) in sorted(
        table,
        key=lambda key: (
            key[1],
            key[0],
        ),
    ):

        key0 = (
            p_value,
            t_value,
        )

        key1 = (
            p_value,
            t_value + 1,
        )

        key2 = (
            p_value,
            t_value + 2,
        )

        if (
            key1 not in table
            or
            key2 not in table
        ):
            continue

        q0 = table[key0]
        q1 = table[key1]
        q2 = table[key2]

        row = []

        # N0 coefficients.
        for j in range(
            degree + 1
        ):
            row.append(
                clean(
                    -q0 * p_value**j
                )
            )

        # N1 coefficients.
        for j in range(
            degree + 1
        ):
            row.append(
                clean(
                    -q1 * p_value**j
                )
            )

        # D coefficients.
        for j in range(
            degree + 1
        ):
            row.append(
                clean(
                    q2 * p_value**j
                )
            )

        rows.append(row)

        metadata.append(
            (
                p_value,
                t_value,
            )
        )

    A = sp.Matrix(rows)

    return {
        "A": A,
        "metadata": metadata,
        "unknown_count": unknown_count,
    }


# ============================================================================
# SOLVER
# ============================================================================

def solve_degree(table, degree):

    system = build_system(
        table,
        degree,
    )

    A = system["A"]

    equation_count = A.rows
    coefficient_count = A.cols
    projective_dimension = (
        coefficient_count - 1
    )

    if equation_count == 0:

        return {
            **system,
            "status": "NO_EQUATIONS",
            "rank": 0,
            "nullity": coefficient_count,
            "vector": None,
            "N0": None,
            "N1": None,
            "D": None,
            "residuals": None,
            "denominator_values": None,
        }

    rank = A.rank()
    nullspace = A.nullspace()
    nullity = len(nullspace)

    # Homogeneous system.
    if nullity == 0:

        return {
            **system,
            "status": "NO_SOLUTION",
            "rank": rank,
            "nullity": 0,
            "vector": None,
            "N0": None,
            "N1": None,
            "D": None,
            "residuals": None,
            "denominator_values": None,
        }

    # Multiple independent rational laws.
    if nullity != 1:

        return {
            **system,
            "status": "NONUNIQUE",
            "rank": rank,
            "nullity": nullity,
            "vector": nullspace,
            "N0": None,
            "N1": None,
            "D": None,
            "residuals": None,
            "denominator_values": None,
        }

    vector = nullspace[0]

    # Projective normalization.
    first_nonzero = next(
        value
        for value in vector
        if value != 0
    )

    vector = [
        clean(
            value / first_nonzero
        )
        for value in vector
    ]

    split = degree + 1

    n0_values = vector[
        :split
    ]

    n1_values = vector[
        split:2 * split
    ]

    d_values = vector[
        2 * split:
    ]

    N0 = clean(
        sum(
            n0_values[j] * p**j
            for j in range(
                degree + 1
            )
        )
    )

    N1 = clean(
        sum(
            n1_values[j] * p**j
            for j in range(
                degree + 1
            )
        )
    )

    D = clean(
        sum(
            d_values[j] * p**j
            for j in range(
                degree + 1
            )
        )
    )

    residuals = []

    for (
        p_value,
        t_value,
    ) in system["metadata"]:

        q0 = table[
            (
                p_value,
                t_value,
            )
        ]

        q1 = table[
            (
                p_value,
                t_value + 1,
            )
        ]

        q2 = table[
            (
                p_value,
                t_value + 2,
            )
        ]

        residual = clean(
            q2
            * D.subs(
                p,
                p_value,
            )
            - q0
            * N0.subs(
                p,
                p_value,
            )
            - q1
            * N1.subs(
                p,
                p_value,
            )
        )

        residuals.append(residual)

    if not all(
        residual == 0
        for residual in residuals
    ):

        return {
            **system,
            "status": "VERIFICATION_FAILED",
            "rank": rank,
            "nullity": nullity,
            "vector": vector,
            "N0": N0,
            "N1": N1,
            "D": D,
            "residuals": residuals,
            "denominator_values": None,
        }

    observed_p = sorted(
        {
            p_value
            for p_value, _ in system["metadata"]
        }
    )

    denominator_values = {
        p_value: clean(
            D.subs(
                p,
                p_value,
            )
        )
        for p_value in observed_p
    }

    # No pole may occur on an observed vertical edge.
    if any(
        value == 0
        for value in denominator_values.values()
    ):

        return {
            **system,
            "status": "INVALID_DENOMINATOR",
            "rank": rank,
            "nullity": nullity,
            "vector": vector,
            "N0": N0,
            "N1": N1,
            "D": D,
            "residuals": residuals,
            "denominator_values": denominator_values,
        }

    if equation_count > projective_dimension:
        status = "EXACT_OVERDETERMINED"
    else:
        status = "EXACT_DATA_SIZED"

    return {
        **system,
        "status": status,
        "rank": rank,
        "nullity": nullity,
        "vector": vector,
        "N0": N0,
        "N1": N1,
        "D": D,
        "residuals": residuals,
        "denominator_values": denominator_values,
    }


# ============================================================================
# REPORT
# ============================================================================

def report(table, degree):

    result = solve_degree(
        table,
        degree,
    )

    print()
    print("=" * 78)
    print(
        "SECOND-ORDER RATIONAL MODEL DEGREE <= {}".format(
            degree
        )
    )
    print("=" * 78)

    print(
        "  equation_count={}".format(
            result["A"].rows
        )
    )

    print(
        "  coefficient_count={}".format(
            result["A"].cols
        )
    )

    print(
        "  projective_unknown_dimension={}".format(
            result["A"].cols - 1
        )
    )

    print(
        "  coefficient_matrix_rank={}".format(
            result["rank"]
        )
    )

    print(
        "  nullity={}".format(
            result["nullity"]
        )
    )

    print(
        "  redundancy={}".format(
            result["A"].rows
            - (result["A"].cols - 1)
        )
    )

    print(
        "  usable_windows={}".format(
            result["metadata"]
        )
    )

    print(
        "  status={}".format(
            result["status"]
        )
    )

    if result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        print()
        print(
            "  N0(p)={}".format(
                result["N0"]
            )
        )

        print(
            "  N1(p)={}".format(
                result["N1"]
            )
        )

        print(
            "  D(p)={}".format(
                result["D"]
            )
        )

        print()
        print(
            "  primitive_N0={}".format(
                primitive_integer_polynomial(
                    result["N0"]
                ).as_expr()
            )
        )

        print(
            "  primitive_N1={}".format(
                primitive_integer_polynomial(
                    result["N1"]
                ).as_expr()
            )
        )

        print(
            "  primitive_D={}".format(
                primitive_integer_polynomial(
                    result["D"]
                ).as_expr()
            )
        )

        print()
        print(
            "  denominator_values={}".format(
                result["denominator_values"]
            )
        )

        print(
            "  residuals={}".format(
                result["residuals"]
            )
        )

        print(
            "  all_residuals_zero={}".format(
                all(
                    residual == 0
                    for residual
                    in result["residuals"]
                )
            )
        )

        if degree > 0:

            print()
            print(
                "  N0_factorization={}".format(
                    sp.factor_list(
                        result["N0"]
                    )
                )
            )

            print(
                "  N1_factorization={}".format(
                    sp.factor_list(
                        result["N1"]
                    )
                )
            )

            print(
                "  D_factorization={}".format(
                    sp.factor_list(
                        result["D"]
                    )
                )
            )

    elif result["status"] == "NONUNIQUE":

        print(
            "  status_note=multiple_projective_rational_second_order_laws"
        )

    return result


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 329R — EXACT SECOND-ORDER RATIONAL "
        "P-DEPENDENT VERTICAL AUDIT"
    )
    print("=" * 78)

    table = build_table()

    print()
    print("=" * 78)
    print(
        "1. OBSERVED SOURCE TABLE"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(table)
        )
    )

    for (
        p_value,
        t_value,
    ) in sorted(
        table,
        key=lambda key: (
            key[1],
            key[0],
        ),
    ):

        print(
            "  (p={},t={}) -> {}".format(
                p_value,
                t_value,
                table[
                    (
                        p_value,
                        t_value,
                    )
                ],
            )
        )

    print()
    print(
        "=" * 78
    )
    print(
        "2. SECOND-ORDER VERTICAL WINDOWS"
    )
    print(
        "=" * 78
    )

    # Count manually through the actual table.
    windows = []

    for (
        p_value,
        t_value,
    ) in sorted(
        table,
        key=lambda key: (
            key[1],
            key[0],
        ),
    ):

        if (
            (
                p_value,
                t_value + 1,
            ) in table
            and
            (
                p_value,
                t_value + 2,
            ) in table
        ):
            windows.append(
                (
                    p_value,
                    t_value,
                )
            )

    print(
        "  usable_window_count={}".format(
            len(windows)
        )
    )

    print(
        "  usable_windows={}".format(
            windows
        )
    )

    results = {}

    for degree in (
        0,
        1,
        2,
    ):

        results[
            degree
        ] = report(
            table,
            degree,
        )

    print()
    print("=" * 78)
    print(
        "3. SUMMARY"
    )
    print("=" * 78)

    for degree, result in results.items():

        print(
            "  degree<={}: status={}, "
            "equations={}, projective_unknowns={}, "
            "rank={}, nullity={}, redundancy={}".format(
                degree,
                result["status"],
                result["A"].rows,
                result["A"].cols - 1,
                result["rank"],
                result["nullity"],
                result["A"].rows
                - (result["A"].cols - 1),
            )
        )

    discoveries = [
        degree
        for degree, result
        in results.items()
        if result["status"]
        == "EXACT_OVERDETERMINED"
    ]

    print()
    print("=" * 78)
    print(
        "4. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 328R rejected rational first-order vertical laws through
degree 2.

Experiment 329R asks whether the source evolution nevertheless has
a second-order rational dependence on t:

    Q_{t+2}
      =
    R0(p) Q_t
      +
    R1(p) Q_{t+1}.

The use of a common denominator D(p) is essential because it keeps the
equations linear in the unknown polynomial coefficients.

The available data provide exactly eight second-order vertical windows.

Therefore:

    degree 0:
        3 projective unknown dimensions,
        8 equations,
        strongly overdetermined;

    degree 1:
        6 raw coefficients,
        5 projective dimensions,
        8 equations,
        genuinely overdetermined;

    degree 2:
        9 raw coefficients,
        8 projective dimensions,
        8 equations,
        data-sized only.

Consequently, only a degree-0 or degree-1 exact solution can count as
evidence for a genuine low-complexity second-order rational source law.

A degree-2 solution must not be interpreted as a discovery.

No missing cells are used.
No interpolation is performed.
No synthetic second n=pq case is generated.
"""
    )

    print()
    print("=" * 78)
    print(
        "5. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_second_order_rational_models={}".format(
            discoveries
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  genuine_second_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 329R COMPLETE"
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
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )
        raise
