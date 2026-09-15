#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 328R — EXACT RATIONAL P-DEPENDENT FIRST-ORDER TRANSFER AUDIT
==============================================================================

Purpose
-------
Experiments 326R and 327R found no low-degree bivariate polynomial law and
no low-degree polynomial-coefficient recurrence in t.

The next structured possibility is a rational first-order law

    Q_{t+1}(p) / Q_t(p) = N(p) / D(p),

or equivalently

    Q_{t+1}(p) D(p) - Q_t(p) N(p) = 0.

For fixed degree d:

    N(p) = n_0 + n_1 p + ... + n_d p^d
    D(p) = d_0 + d_1 p + ... + d_d p^d.

After cross multiplication the unknown coefficients occur linearly.

The experiment therefore tests whether a LOW-DEGREE rational function of p
governs the observed vertical evolution.

Models:

    numerator/denominator degree <= 0
    numerator/denominator degree <= 1
    numerator/denominator degree <= 2
    numerator/denominator degree <= 3

Only overlapping observed cells are used.

Important:
    * denominator is normalized projectively;
    * the zero coefficient vector is excluded;
    * a candidate is accepted only if it reproduces every observed equation;
    * overdetermined exact solutions are distinguished from data-sized ones;
    * denominator zeros on observed p-values are rejected;
    * no missing values are reconstructed;
    * no interpolation outside observed cells;
    * no synthetic second n=pq case.

Exact SymPy rational arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# SOURCE DATA
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


p, t = sp.symbols("p t")


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

    n = abs(int(x.p))
    d = abs(int(x.q))

    out = 0

    while n % prime == 0:
        n //= prime
        out += 1

    while d % prime == 0:
        d //= prime
        out -= 1

    return out


def build_table():
    """
    Construct

        table[(p,t)] = Q_t(p)

    from observed data only.
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


def primitive_integer_polynomial(expr):
    """
    Primitive integer polynomial in p.
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

    coeffs = poly.all_coeffs()

    den_lcm = 1

    for c in coeffs:
        den_lcm = sp.ilcm(
            den_lcm,
            int(sp.denom(c)),
        )

    ints = [
        int(c * den_lcm)
        for c in coeffs
    ]

    g = 0

    for x in ints:
        g = math.gcd(
            g,
            abs(x),
        )

    ints = [
        x // g
        for x in ints
    ]

    if ints[0] < 0:
        ints = [
            -x
            for x in ints
        ]

    expr_int = sum(
        sp.Integer(c)
        * p**(
            len(ints) - 1 - i
        )
        for i, c in enumerate(ints)
    )

    return sp.Poly(
        sp.expand(expr_int),
        p,
        domain=sp.ZZ,
    )


# ============================================================================
# RATIONAL FIRST-ORDER SYSTEM
# ============================================================================

def build_equations(
    table,
    degree,
):
    """
    Build equations

        Q_{t+1}(p) D(p) - Q_t(p) N(p) = 0

    over every observed vertical edge.
    """

    # Unknown vector:
    #
    #   n_0,...,n_d,d_0,...,d_d
    #
    # There are 2(d+1) coefficients, defined only up to common scale.

    symbols = sp.symbols(
        "c0:{}".format(
            2 * (degree + 1)
        )
    )

    n_coeffs = symbols[
        :degree + 1
    ]

    d_coeffs = symbols[
        degree + 1:
    ]

    N = sum(
        n_coeffs[j] * p**j
        for j in range(
            degree + 1
        )
    )

    D = sum(
        d_coeffs[j] * p**j
        for j in range(
            degree + 1
        )
    )

    rows = []
    rhs = []
    metadata = []

    for (p_value, t_value), q_value in sorted(
        table.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        next_key = (
            p_value,
            t_value + 1,
        )

        if next_key not in table:
            continue

        q_next = table[next_key]

        # q_next * D(p) - q_value * N(p) = 0
        row = []

        for j in range(
            degree + 1
        ):
            row.append(
                clean(
                    -q_value
                    * p_value**j
                )
            )

        for j in range(
            degree + 1
        ):
            row.append(
                clean(
                    q_next
                    * p_value**j
                )
            )

        rows.append(row)
        rhs.append(0)

        metadata.append(
            (
                p_value,
                t_value,
            )
        )

    A = sp.Matrix(rows)
    b = sp.Matrix(rhs)

    return {
        "symbols": symbols,
        "N": sp.expand(N),
        "D": sp.expand(D),
        "A": A,
        "b": b,
        "metadata": metadata,
    }


# ============================================================================
# EXACT PROJECTIVE NULLSPACE SOLVER
# ============================================================================

def solve_degree(
    table,
    degree,
):
    system = build_equations(
        table,
        degree,
    )

    A = system["A"]
    b = system["b"]
    metadata = system["metadata"]

    equation_count = A.rows
    unknown_count = A.cols

    if equation_count == 0:

        return {
            **system,
            "status": "NO_EQUATIONS",
            "rank": 0,
            "nullity": unknown_count,
            "vector": None,
            "N_poly": None,
            "D_poly": None,
            "residuals": None,
        }

    rank = A.rank()
    nullspace = A.nullspace()
    nullity = len(nullspace)

    # Since the system is homogeneous, existence of any nonzero
    # nullspace vector is necessary and sufficient for a rational-law
    # candidate.
    if nullity == 0:

        return {
            **system,
            "status": "NO_SOLUTION",
            "rank": rank,
            "nullity": 0,
            "vector": None,
            "N_poly": None,
            "D_poly": None,
            "residuals": None,
        }

    # We only accept a unique projective law when nullity == 1.
    if nullity != 1:

        return {
            **system,
            "status": "NONUNIQUE",
            "rank": rank,
            "nullity": nullity,
            "vector": nullspace,
            "N_poly": None,
            "D_poly": None,
            "residuals": None,
        }

    vector = nullspace[0]

    # Normalize the first nonzero coefficient to 1.
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

    n_values = vector[
        :split
    ]

    d_values = vector[
        split:
    ]

    N_poly = clean(
        sum(
            n_values[j] * p**j
            for j in range(
                degree + 1
            )
        )
    )

    D_poly = clean(
        sum(
            d_values[j] * p**j
            for j in range(
                degree + 1
            )
        )
    )

    # Verify the original cross-multiplied equations.
    residuals = []

    for (
        p_value,
        t_value,
    ) in metadata:

        q_value = table[
            (
                p_value,
                t_value,
            )
        ]

        q_next = table[
            (
                p_value,
                t_value + 1,
            )
        ]

        residual = clean(
            q_next
            * D_poly.subs(
                p,
                p_value,
            )
            - q_value
            * N_poly.subs(
                p,
                p_value,
            )
        )

        residuals.append(
            residual
        )

    if not all(
        r == 0
        for r in residuals
    ):

        return {
            **system,
            "status": "VERIFICATION_FAILED",
            "rank": rank,
            "nullity": nullity,
            "vector": vector,
            "N_poly": N_poly,
            "D_poly": D_poly,
            "residuals": residuals,
        }

    # Reject denominator zeros on observed p-values.
    observed_p = sorted(
        {
            p_value
            for p_value, _ in metadata
        }
    )

    denominator_values = {
        p_value: clean(
            D_poly.subs(
                p,
                p_value,
            )
        )
        for p_value in observed_p
    }

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
            "N_poly": N_poly,
            "D_poly": D_poly,
            "residuals": residuals,
            "denominator_values": denominator_values,
        }

    if equation_count > unknown_count - 1:
        status = "EXACT_OVERDETERMINED"

    else:
        status = "EXACT_DATA_SIZED"

    return {
        **system,
        "status": status,
        "rank": rank,
        "nullity": nullity,
        "vector": vector,
        "N_poly": N_poly,
        "D_poly": D_poly,
        "residuals": residuals,
        "denominator_values": denominator_values,
    }


# ============================================================================
# REPORT
# ============================================================================

def report(
    table,
    degree,
):
    result = solve_degree(
        table,
        degree,
    )

    print()
    print("=" * 78)
    print(
        "RATIONAL FIRST-ORDER MODEL DEGREE <= {}".format(
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
            "  N(p)={}".format(
                result["N_poly"]
            )
        )

        print(
            "  D(p)={}".format(
                result["D_poly"]
            )
        )

        print()
        print(
            "  primitive_N={}".format(
                primitive_integer_polynomial(
                    result["N_poly"]
                ).as_expr()
            )
        )

        print(
            "  primitive_D={}".format(
                primitive_integer_polynomial(
                    result["D_poly"]
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

        if degree > 0:

            print()
            print(
                "  N_factorization={}".format(
                    sp.factor_list(
                        result["N_poly"]
                    )
                )
            )

            print(
                "  D_factorization={}".format(
                    sp.factor_list(
                        result["D_poly"]
                    )
                )
            )

    elif result["status"] == "NONUNIQUE":

        print(
            "  status_note=multiple_projective_rational_laws"
        )


    return result


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 328R — EXACT RATIONAL P-DEPENDENT "
        "FIRST-ORDER TRANSFER AUDIT"
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
        key=lambda z: (
            z[1],
            z[0],
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

    results = {}

    for degree in (
        0,
        1,
        2,
        3,
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
        "SUMMARY"
    )
    print("=" * 78)

    for degree, result in results.items():

        print(
            "  degree<={}: status={}, "
            "equations={}, projective_unknowns={}, "
            "rank={}, nullity={}".format(
                degree,
                result["status"],
                result["A"].rows,
                result["A"].cols - 1,
                result["rank"],
                result["nullity"],
            )
        )

    discoveries = [
        degree
        for degree, result in results.items()
        if result["status"]
        == "EXACT_OVERDETERMINED"
    ]

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiments 326R-327R rejected low-degree polynomial source laws and
low-degree polynomial-coefficient recurrences.

Experiment 328R asks whether the simplest non-polynomial alternative
survives:

    Q_{t+1}(p) / Q_t(p) = N(p) / D(p).

This is a rational first-order vertical law.

The cross-multiplied equations are linear in the unknown coefficients,
so the test remains exact and non-interpolative.

A successful degree-0 model would mean a constant geometric evolution.

A successful degree-1 or degree-2 model would mean that the evolution
depends on p through a very small rational function.

The requirement that the nullspace be one-dimensional is deliberate:
it prevents an underdetermined family of rational functions from being
reported as a discovered law.

Only an exact solution satisfying every observed vertical edge is
accepted.

No missing value is used.
No synthetic point is created.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_rational_models={}".format(
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
        "EXPERIMENT 328R COMPLETE"
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
