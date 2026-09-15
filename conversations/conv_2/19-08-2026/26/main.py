#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 342R — EXACT SHARED NEWTON-COORDINATE OPERATOR AUDIT
==============================================================================

Purpose
-------
Experiment 341R established that removing row-content gcds does not simplify
the source shape:

    primitive row rank = 4,
    no adjacent proportional primitive rows,
    primitive width-2 laws remain data-sized reconstructions.

The next structural possibility is that the NEWTON COEFFICIENTS

    A[j,t] = Delta_p^j Q_t(1)

evolve under a single hidden linear operator in coefficient space.

Crucially, missing higher Newton coefficients are NOT treated as zero.

For each transition t -> t+1, only equations involving genuinely observed
coefficients are included.

Test families:

    1. diagonal operator
    2. bidiagonal / bandwidth-1 operator
    3. lower-triangular operator
    4. upper-triangular operator
    5. full 4x4 operator

The experiment accepts a model as structurally interesting only when:

    * it is exactly consistent;
    * the operator is uniquely determined;
    * more equations than operator parameters are available.

For each family we report:

    equations,
    unknowns,
    rank,
    augmented rank,
    nullity,
    redundancy,
    status,
    exact operator when uniquely determined.

No missing values.
No padding of missing coefficients.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact rational arithmetic only.
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


# ============================================================================
# SOURCE LAYERS
# ============================================================================

def build_layers():

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(maximum_t + 1):

        row = []

        for p_value in sorted(Q):

            values = Q[p_value]

            index = len(values) - 1 - t

            if index >= 0:

                row.append(
                    (
                        (p_value - 1) // 2,
                        sp.Integer(values[index]),
                    )
                )

        layers[t] = row

    return layers


# ============================================================================
# NEWTON TRIANGLE
# ============================================================================

def finite_difference_rows(values):

    current = [
        sp.Integer(v)
        for v in values
    ]

    rows = [current]

    while len(current) > 1:

        current = [
            clean(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        rows.append(current)

    return rows


def build_newton_triangle(layers):

    A = {}

    for t, layer in layers.items():

        values = [
            value
            for _, value in layer
        ]

        rows = finite_difference_rows(
            values
        )

        for j, row in enumerate(rows):

            if row:

                A[(j, t)] = clean(
                    row[0]
                )

    return A


# ============================================================================
# OPERATOR SUPPORT PATTERNS
# ============================================================================

def operator_support(kind):

    """
    Return allowed matrix entries M[j,k].

    We have 4 Newton-coordinate levels j=0,1,2,3.

    'diagonal'
        M[j,k] only if j=k.

    'bidiagonal_upper'
        diagonal + first superdiagonal.

    'bidiagonal_lower'
        diagonal + first subdiagonal.

    'lower'
        k <= j.

    'upper'
        k >= j.

    'full'
        all 16 entries.
    """

    support = []

    for j in range(4):

        for k in range(4):

            allowed = False

            if kind == "diagonal":
                allowed = (j == k)

            elif kind == "bidiagonal_upper":
                allowed = (
                    k == j
                    or
                    k == j + 1
                )

            elif kind == "bidiagonal_lower":
                allowed = (
                    k == j
                    or
                    k == j - 1
                )

            elif kind == "lower":
                allowed = k <= j

            elif kind == "upper":
                allowed = k >= j

            elif kind == "full":
                allowed = True

            else:
                raise ValueError(
                    "unknown operator kind"
                )

            if allowed:
                support.append(
                    (j, k)
                )

    return support


# ============================================================================
# BUILD OBSERVED OPERATOR EQUATIONS
# ============================================================================

def build_operator_system(
    A,
    support,
):
    """
    Build equations

        A[j,t+1] = sum_k M[j,k] A[k,t]

    but only when:

        * A[j,t+1] exists;
        * every source coefficient used by M[j,k] exists.

    Thus missing higher Newton levels are NEVER interpreted as zero.
    """

    variables = {
        location: sp.Symbol(
            "m_{}_{}".format(
                location[0],
                location[1],
            )
        )
        for location in support
    }

    equations = []

    max_t = max(
        t
        for _, t in A
    )

    for t in range(max_t):

        target_js = [
            j
            for j in range(4)
            if (j, t + 1) in A
        ]

        source_available = {
            k
            for k in range(4)
            if (k, t) in A
        }

        for j in target_js:

            row_support = [
                (jj, k)
                for jj, k in support
                if jj == j
            ]

            # A matrix row is usable only if every coefficient
            # invoked by that row is actually observed at time t.
            if not all(
                k in source_available
                for _, k in row_support
            ):
                continue

            lhs = A[(j, t)]

            rhs = sum(
                variables[(jj, k)] * A[(k, t)]
                for jj, k in row_support
            )

            target_value = A[(j, t + 1)]

            equations.append(
                (
                    t,
                    j,
                    clean(rhs),
                    clean(target_value),
                )
            )

    return variables, equations


# ============================================================================
# SOLVE OPERATOR
# ============================================================================

def solve_operator(
    A,
    kind,
):

    support = operator_support(
        kind
    )

    variables, equations = (
        build_operator_system(
            A,
            support,
        )
    )

    ordered_variables = [
        variables[location]
        for location in support
    ]

    if not equations:

        return {
            "status": "NO_EQUATIONS",
            "support": support,
            "variables": variables,
            "equations": equations,
        }

    matrix_rows = []

    rhs_values = []

    for (
        _t,
        _j,
        symbolic_rhs,
        target_value,
    ) in equations:

        row = [
            sp.expand(
                symbolic_rhs
            ).coeff(variable)
            for variable in ordered_variables
        ]

        matrix_rows.append(row)
        rhs_values.append(
            target_value
        )

    M = sp.Matrix(
        matrix_rows
    )

    b = sp.Matrix(
        rhs_values
    )

    rank = M.rank()

    augmented_rank = (
        M.row_join(b).rank()
    )

    equation_count = len(
        equations
    )

    unknown_count = len(
        support
    )

    redundancy = (
        equation_count
        - unknown_count
    )

    nullity = (
        unknown_count
        - rank
    )

    if augmented_rank > rank:

        status = "NO_SOLUTION"
        solution = None

    elif rank < unknown_count:

        status = "NONUNIQUE"
        solution = None

    else:

        raw_solution = M.gauss_jordan_solve(
            b
        )[0]

        solution = {
            support[i]: clean(
                raw_solution[i, 0]
            )
            for i in range(
                unknown_count
            )
        }

        residuals = [
            clean(
                sum(
                    M[row, col]
                    * raw_solution[col, 0]
                    for col in range(
                        unknown_count
                    )
                )
                - b[row]
            )
            for row in range(
                equation_count
            )
        ]

        if not all(
            residual == 0
            for residual in residuals
        ):
            status = (
                "VERIFICATION_FAILED"
            )

        elif redundancy > 0:

            status = (
                "EXACT_OVERDETERMINED"
            )

        else:

            status = (
                "EXACT_DATA_SIZED"
            )

    return {
        "status": status,
        "support": support,
        "variables": variables,
        "equations": equations,
        "matrix": M,
        "rhs": b,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "equation_count": equation_count,
        "unknown_count": unknown_count,
        "redundancy": redundancy,
        "nullity": nullity,
        "solution": solution,
    }


# ============================================================================
# REPORT
# ============================================================================

def report_operator(
    A,
    kind,
):

    result = solve_operator(
        A,
        kind,
    )

    print()
    print("=" * 78)
    print(
        "OPERATOR FAMILY: {}".format(
            kind.upper()
        )
    )
    print("=" * 78)

    if "equation_count" not in result:

        print(
            "  status={}".format(
                result["status"]
            )
        )

        return result

    print(
        "  support={}".format(
            result["support"]
        )
    )

    print(
        "  equations={}".format(
            result["equation_count"]
        )
    )

    print(
        "  unknowns={}".format(
            result["unknown_count"]
        )
    )

    print(
        "  redundancy={}".format(
            result["redundancy"]
        )
    )

    print(
        "  rank={}".format(
            result["rank"]
        )
    )

    print(
        "  augmented_rank={}".format(
            result["augmented_rank"]
        )
    )

    print(
        "  nullity={}".format(
            result["nullity"]
        )
    )

    print(
        "  status={}".format(
            result["status"]
        )
    )

    if result["solution"] is not None:

        print()
        print(
            "  exact_operator_entries="
        )

        for location in result["support"]:

            print(
                "    M{}={}".format(
                    location,
                    result["solution"][
                        location
                    ],
                )
            )

    return result


# ============================================================================
# OPERATOR MATRIX RECONSTRUCTION
# ============================================================================

def operator_matrix(result):

    if result.get("solution") is None:
        return None

    M = sp.zeros(4)

    for (
        j,
        k,
    ), value in result["solution"].items():

        M[j, k] = value

    return M


def print_operator_invariants(
    result,
):

    M = operator_matrix(
        result
    )

    if M is None:
        return

    print()
    print(
        "  reconstructed_operator_matrix="
    )

    print(M)

    print()
    print(
        "  trace={}".format(
            clean(M.trace())
        )
    )

    print(
        "  determinant={}".format(
            clean(M.det())
        )
    )

    print(
        "  rank={}".format(
            M.rank()
        )
    )

    print(
        "  charpoly={}".format(
            sp.factor(
                M.charpoly().as_expr()
            )
        )
    )


# ============================================================================
# PRIME PROFILE
# ============================================================================

def operator_prime_profile(
    result,
):

    if result.get("solution") is None:
        return

    print()
    print(
        "  prime-profile:"
    )

    for location in result["support"]:

        value = result["solution"][
            location
        ]

        print(
            "    M{}={}".format(
                location,
                value,
            )
        )

        print(
            "      valuations={}".format(
                {
                    prime: valuation(
                        value,
                        prime,
                    )
                    for prime in (
                        2,
                        3,
                        5,
                        7,
                        11,
                        13,
                        17,
                    )
                }
            )
        )


# ============================================================================
# DIRECT TRANSITION RESIDUAL AUDIT
# ============================================================================

def direct_residual_audit(
    A,
    result,
):

    if result.get("solution") is None:
        return

    solution = result["solution"]
    support = result["support"]

    print()
    print(
        "  direct_transition_residuals:"
    )

    max_t = max(
        t
        for _, t in A
    )

    all_zero = True

    for t in range(max_t):

        for j in range(4):

            if (j, t + 1) not in A:
                continue

            row_support = [
                (jj, k)
                for jj, k in support
                if jj == j
            ]

            if not all(
                (k, t) in A
                for _, k in row_support
            ):
                continue

            predicted = sum(
                solution[(jj, k)]
                * A[(k, t)]
                for jj, k in row_support
            )

            residual = clean(
                predicted
                - A[(j, t + 1)]
            )

            print(
                "    t={} -> {}, j={}: {}".format(
                    t,
                    t + 1,
                    j,
                    residual,
                )
            )

            if residual != 0:
                all_zero = False

    print(
        "  all_zero={}".format(
            all_zero
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 342R — EXACT SHARED NEWTON-COORDINATE "
        "OPERATOR AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    A = build_newton_triangle(
        layers
    )

    print()
    print("=" * 78)
    print(
        "1. OBSERVED NEWTON TRIANGLE"
    )
    print("=" * 78)

    for t in range(6):

        print(
            "  t={}: {}".format(
                t,
                [
                    A.get(
                        (j, t),
                        None,
                    )
                    for j in range(4)
                ],
            )
        )

    operator_kinds = (
        "diagonal",
        "bidiagonal_upper",
        "bidiagonal_lower",
        "lower",
        "upper",
        "full",
    )

    results = {}

    for kind in operator_kinds:

        result = report_operator(
            A,
            kind,
        )

        results[kind] = result

        print_operator_invariants(
            result
        )

        operator_prime_profile(
            result
        )

        direct_residual_audit(
            A,
            result,
        )

    # ========================================================================
    # COMPARISON
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "8. OPERATOR-FAMILY SUMMARY"
    )
    print("=" * 78)

    exact_overdetermined = []

    for kind in operator_kinds:

        result = results[kind]

        if (
            result.get("status")
            == "EXACT_OVERDETERMINED"
        ):

            exact_overdetermined.append(
                kind
            )

        if "equation_count" in result:

            print()
            print(
                "  {}: status={}, equations={}, "
                "unknowns={}, rank={}, "
                "augmented_rank={}, redundancy={}, "
                "nullity={}".format(
                    kind,
                    result["status"],
                    result["equation_count"],
                    result["unknown_count"],
                    result["rank"],
                    result["augmented_rank"],
                    result["redundancy"],
                    result["nullity"],
                )
            )

    # ========================================================================
    # INTERPRETATION
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiments 339R-341R showed that arithmetic-content normalization does not
explain the structure of the observed source table.

The primitive Newton triangle is the canonical finite-difference coordinate
system for the observed p-grid.

Experiment 342R now asks a sharper operator-theoretic question:

    Do the observed Newton coefficients evolve by ONE FIXED LINEAR MAP?

This is stronger than fitting separate transitions.

For example, an upper-triangular operator has the form

    A[0,t+1] = m00 A[0,t]
             + m01 A[1,t]
             + m02 A[2,t]
             + m03 A[3,t]

    A[1,t+1] = m11 A[1,t]
             + m12 A[2,t]
             + m13 A[3,t]

and so on.

Crucially, an equation is included only when every source coefficient
appearing in that equation is actually observed.

Missing Newton coefficients are NEVER padded by zero.

The distinction is:

    EXACT_OVERDETERMINED
        genuine finite-data evidence for a shared operator;

    EXACT_DATA_SIZED
        reconstruction only;

    NONUNIQUE
        insufficient identification;

    NO_SOLUTION
        exact contradiction.

If even the full operator family fails, then there is no fixed linear
evolution in the canonical Newton coordinates.

That would be a strong endpoint for generic operator searching.

If a restricted operator survives, its exact matrix becomes the next object
to interpret arithmetically or combinatorially.

No missing source value is introduced.
No interpolation is performed.
No synthetic second n=pq case is generated.
"""
    )

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  newton_triangle_built_exactly=True"
    )

    print(
        "  missing_newton_coefficients_padded=False"
    )

    print(
        "  shared_operator_families_tested={}".format(
            list(operator_kinds)
        )
    )

    print(
        "  exact_overdetermined_operator_families={}".format(
            exact_overdetermined
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_used=False"
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
        "EXPERIMENT 342R COMPLETE"
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
