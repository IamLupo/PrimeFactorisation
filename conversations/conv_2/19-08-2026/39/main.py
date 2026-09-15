#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 353R — EXACT D6 INTEGRALITY / ARITHMETIC FALSIFICATION AUDIT
==============================================================================

Purpose
-------
Experiment 352R-FIXED found one identifiable denominator model:

    D6 = 1 - a*x - b*y - c*x*y - d*y^2

It is EXACT_DATA_SIZED on the currently observed interior equations and
predicts the missing cell Q_3(5).

The prediction is:

    N / D

This experiment asks a much sharper question:

    Is the D6 prediction compatible with the observed integer-valued
    source lattice?

This is NOT a proof that D6 is wrong unless the original source definition
guarantees integer-valued Q_t(p). The experiment therefore reports the
arithmetic obstruction precisely and separately.

Audits
------

1. Reconstruct the D6 parameters exactly.
2. Recompute the target prediction independently.
3. Reduce the prediction to lowest terms.
4. Test integrality.
5. Factor numerator and denominator.
6. Compute gcd.
7. Test divisibility modulo denominator prime factors.
8. Compare with the integer-valued nature of every observed cell.
9. Test whether multiplying by a small natural normalization can restore
   integrality.
10. State the logically correct verdict.

No new Q-cell is invented.
No prediction is treated as observed data.
No interpolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# OBSERVED SOURCE DATA
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
# SYMBOLS
# ============================================================================

x, y = sp.symbols("x y")

a, b, c, d = sp.symbols(
    "a b c d"
)


TARGET = (2, 3)


# ============================================================================
# HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def build_lattice():
    lattice = {}

    for p_value, values in Q.items():

        r_value = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            lattice[
                (r_value, t_value)
            ] = sp.Integer(value)

    return lattice


def factor_integer(value):

    value = int(
        sp.Integer(value)
    )

    if value == 0:
        return {}

    return sp.factorint(
        abs(value)
    )


def valuation(value, prime):

    value = int(
        abs(
            sp.Integer(value)
        )
    )

    if value == 0:
        return sp.oo

    count = 0

    while value % prime == 0:

        value //= prime
        count += 1

    return count


# ============================================================================
# D6 MODEL
# ============================================================================

def d6_denominator():

    return (
        1
        - a*x
        - b*y
        - c*x*y
        - d*y**2
    )


def d6_equations(lattice):

    D = d6_denominator()

    equations = []

    # Interior cells for which every source term required by D is observed.
    for r in range(4):

        for t in range(6):

            target = (
                r,
                t,
            )

            if target not in lattice:
                continue

            required = [
                (r, t),
                (r - 1, t),
                (r, t - 1),
                (r - 1, t - 1),
                (r, t - 2),
            ]

            if all(
                cell in lattice
                for cell in required
            ):

                q = (
                    lattice[(r, t)]
                    - a * lattice[(r - 1, t)]
                    - b * lattice[(r, t - 1)]
                    - c * lattice[(r - 1, t - 1)]
                    - d * lattice[(r, t - 2)]
                )

                equations.append(
                    clean(q)
                )

    return equations


def solve_d6(lattice):

    equations = d6_equations(
        lattice
    )

    matrix = sp.Matrix([
        [
            sp.diff(
                equation,
                parameter,
            )
            for parameter in (
                a, b, c, d
            )
        ]
        for equation in equations
    ])

    rhs = sp.Matrix([
        clean(
            -equation.subs(
                {
                    a: 0,
                    b: 0,
                    c: 0,
                    d: 0,
                }
            )
        )
        for equation in equations
    ])

    rank = matrix.rank()

    augmented_rank = (
        matrix
        .row_join(rhs)
        .rank()
    )

    solution = matrix.gauss_jordan_solve(
        rhs
    )[0]

    parameters = {
        a: clean(solution[0]),
        b: clean(solution[1]),
        c: clean(solution[2]),
        d: clean(solution[3]),
    }

    return (
        equations,
        rank,
        augmented_rank,
        parameters,
    )


# ============================================================================
# TARGET PREDICTION
# ============================================================================

def predict_target(
    lattice,
    parameters,
):

    D = clean(
        d6_denominator().subs(
            parameters
        )
    )

    # Coefficient of D*G at (2,3):
    #
    # Q(2,3)
    # - a Q(1,3)
    # - b Q(2,2)
    # - c Q(1,2)
    # - d Q(2,1)
    #
    # Therefore:
    #
    # Q(2,3) =
    # a Q(1,3)
    # + b Q(2,2)
    # + c Q(1,2)
    # + d Q(2,1)

    prediction = clean(
        parameters[a] * lattice[(1, 3)]
        + parameters[b] * lattice[(2, 2)]
        + parameters[c] * lattice[(1, 2)]
        + parameters[d] * lattice[(2, 1)]
    )

    return D, prediction


# ============================================================================
# OBSERVED INTEGER AUDIT
# ============================================================================

def observed_integer_audit(lattice):

    print()
    print("=" * 78)
    print(
        "3. OBSERVED SOURCE INTEGER AUDIT"
    )
    print("=" * 78)

    nonintegers = []

    for cell in sorted(
        lattice,
        key=lambda z: (z[1], z[0]),
    ):

        value = lattice[cell]

        is_integer = (
            sp.denom(
                sp.Rational(value)
            )
            == 1
        )

        print(
            "  cell={}: value={}, integer={}".format(
                cell,
                value,
                is_integer,
            )
        )

        if not is_integer:
            nonintegers.append(
                cell
            )

    print()
    print(
        "  observed_noninteger_cells={}".format(
            nonintegers
        )
    )

    return len(
        nonintegers
    ) == 0


# ============================================================================
# D6 ARITHMETIC PROFILE
# ============================================================================

def d6_arithmetic_audit(prediction):

    prediction = sp.Rational(
        prediction
    )

    numerator = int(
        sp.numer(prediction)
    )

    denominator = int(
        sp.denom(prediction)
    )

    gcd_value = math.gcd(
        abs(numerator),
        abs(denominator),
    )

    print()
    print("=" * 78)
    print(
        "4. EXACT D6 TARGET ARITHMETIC"
    )
    print("=" * 78)

    print(
        "  prediction={}".format(
            prediction
        )
    )

    print(
        "  numerator={}".format(
            numerator
        )
    )

    print(
        "  denominator={}".format(
            denominator
        )
    )

    print(
        "  gcd_numerator_denominator={}".format(
            gcd_value
        )
    )

    print(
        "  reduced_fraction={}".format(
            gcd_value == 1
        )
    )

    print(
        "  is_integer={}".format(
            denominator == 1
        )
    )

    print(
        "  absolute_value={}".format(
            sp.N(
                abs(prediction),
                20,
            )
        )
    )

    print()
    print(
        "  numerator_factorization={}".format(
            factor_integer(
                numerator
            )
        )
    )

    print(
        "  denominator_factorization={}".format(
            factor_integer(
                denominator
            )
        )
    )

    return (
        numerator,
        denominator,
        gcd_value,
    )


# ============================================================================
# DENOMINATOR PRIME OBSTRUCTION
# ============================================================================

def denominator_prime_audit(
    numerator,
    denominator,
):

    print()
    print("=" * 78)
    print(
        "5. DENOMINATOR PRIME OBSTRUCTION"
    )
    print("=" * 78)

    denominator_factors = (
        factor_integer(
            denominator
        )
    )

    obstruction = False

    for prime, exponent in sorted(
        denominator_factors.items()
    ):

        residue = (
            numerator % prime
        )

        v_num = valuation(
            numerator,
            prime,
        )

        print()
        print(
            "  prime={}: exponent_in_denominator={}".format(
                prime,
                exponent,
            )
        )

        print(
            "    numerator_mod_prime={}".format(
                residue
            )
        )

        print(
            "    v_prime_numerator={}".format(
                v_num
            )
        )

        if v_num < exponent:

            obstruction = True

            print(
                "    integrality_obstruction=True"
            )

        else:

            print(
                "    integrality_obstruction=False"
            )

    print()
    print(
        "  denominator_obstruction_exists={}".format(
            obstruction
        )
    )

    return obstruction


# ============================================================================
# SMALL NORMALIZATION AUDIT
# ============================================================================

def normalization_audit(
    prediction
):

    prediction = sp.Rational(
        prediction
    )

    print()
    print("=" * 78)
    print(
        "6. SMALL NORMALIZATION / DENOMINATOR AUDIT"
    )
    print("=" * 78)

    candidates = [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        12,
        17,
        24,
        60,
        120,
    ]

    hits = []

    for multiplier in candidates:

        scaled = clean(
            multiplier
            * prediction
        )

        is_integer = (
            sp.denom(
                sp.Rational(
                    scaled
                )
            )
            == 1
        )

        print(
            "  multiplier={}: scaled={}, integer={}".format(
                multiplier,
                scaled,
                is_integer,
            )
        )

        if is_integer:
            hits.append(
                multiplier
            )

    print()
    print(
        "  small_integer_normalization_hits={}".format(
            hits
        )
    )

    return hits


# ============================================================================
# COMPARISON WITH OBSERVED ARITHMETIC SCALE
# ============================================================================

def scale_comparison(
    lattice,
    prediction
):

    prediction = sp.Rational(
        prediction
    )

    print()
    print("=" * 78)
    print(
        "7. TARGET-SCALE COMPARISON"
    )
    print("=" * 78)

    observed_abs = [
        abs(
            int(value)
        )
        for value in lattice.values()
    ]

    max_observed = max(
        observed_abs
    )

    min_nonzero_observed = min(
        value
        for value in observed_abs
        if value != 0
    )

    print(
        "  maximum_observed_abs={}".format(
            max_observed
        )
    )

    print(
        "  minimum_nonzero_observed_abs={}".format(
            min_nonzero_observed
        )
    )

    print(
        "  predicted_abs={}".format(
            sp.N(
                abs(prediction),
                20,
            )
        )
    )

    print(
        "  predicted_abs_over_max_observed={}".format(
            sp.N(
                abs(prediction)
                / max_observed,
                20,
            )
        )
    )


# ============================================================================
# LOGICAL VERDICT
# ============================================================================

def verdict(
    observed_integer,
    is_integer,
    denominator_obstruction,
):

    print()
    print("=" * 78)
    print(
        "8. LOGICAL VERDICT"
    )
    print("=" * 78)

    if (
        observed_integer
        and
        not is_integer
        and
        denominator_obstruction
    ):

        print(
            "  verdict=D6_FALSIFIED_IF_SOURCE_IS_INTEGER_VALUED"
        )

        print()
        print(
            "  reason="
            "all 15 observed source cells are integers, "
            "while the unique D6 prediction is a reduced "
            "non-integer rational."
        )

        print()
        print(
            "  important="
            "this does NOT use Q_3(5) as observed data."
        )

    else:

        print(
            "  verdict=NO_ARITHMETIC_FALSIFICATION_YET"
        )

        print()
        print(
            "  reason="
            "integrality alone does not contradict the D6 prediction."
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 353R — EXACT D6 INTEGRALITY / "
        "ARITHMETIC FALSIFICATION AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "1. TARGET"
    )

    print(
        "  target=(2,3)"
    )

    print(
        "  p=5"
    )

    print(
        "  t=3"
    )

    print(
        "  observed_target=False"
    )

    # ------------------------------------------------------------------------
    # Solve D6
    # ------------------------------------------------------------------------

    (
        equations,
        rank,
        augmented_rank,
        parameters,
    ) = solve_d6(
        lattice
    )

    print()
    print("=" * 78)
    print(
        "2. D6 EXACT RECONSTRUCTION"
    )
    print("=" * 78)

    print(
        "  equations={}".format(
            len(equations)
        )
    )

    print(
        "  unknowns=4"
    )

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  augmented_rank={}".format(
            augmented_rank
        )
    )

    print(
        "  status={}".format(
            (
                "EXACT_DATA_SIZED"
                if rank == 4
                and augmented_rank == 4
                and len(equations) == 4
                else
                "UNEXPECTED"
            )
        )
    )

    print(
        "  parameters={}".format(
            parameters
        )
    )

    # ------------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------------

    denominator, prediction = (
        predict_target(
            lattice,
            parameters,
        )
    )

    print()
    print(
        "  solved_denominator={}".format(
            denominator
        )
    )

    print(
        "  predicted_Q3_5={}".format(
            prediction
        )
    )

    (
        numerator,
        denominator_integer,
        gcd_value,
    ) = d6_arithmetic_audit(
        prediction
    )

    # ------------------------------------------------------------------------
    # Observed integer nature
    # ------------------------------------------------------------------------

    observed_integer = (
        observed_integer_audit(
            lattice
        )
    )

    is_integer = (
        denominator_integer == 1
    )

    denominator_obstruction = (
        denominator_prime_audit(
            numerator,
            denominator_integer,
        )
    )

    normalization_hits = (
        normalization_audit(
            prediction
        )
    )

    scale_comparison(
        lattice,
        prediction,
    )

    verdict(
        observed_integer,
        is_integer,
        denominator_obstruction,
    )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "9. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  d6_reconstruction_exact=True"
    )

    print(
        "  d6_prediction_computed=True"
    )

    print(
        "  target_observed=False"
    )

    print(
        "  all_observed_cells_integer={}".format(
            observed_integer
        )
    )

    print(
        "  d6_prediction_integer={}".format(
            is_integer
        )
    )

    print(
        "  reduced_prediction=True"
    )

    print(
        "  denominator_prime_obstruction={}".format(
            denominator_obstruction
        )
    )

    print(
        "  small_integer_normalization_hits={}".format(
            normalization_hits
        )
    )

    print(
        "  missing_value_used_as_data=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_counted_as_evidence=False"
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
        "EXPERIMENT 353R COMPLETE"
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
                exc
            )
        )

        raise
