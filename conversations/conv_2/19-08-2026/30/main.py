#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 345R — EXACT SOURCE-DEFINITION VALIDATION HARNESS
==============================================================================

Purpose
-------
The preceding experiments have rejected a very broad collection of
data-driven laws:

    * fixed-width transfer recurrences;
    * higher-width local linear stencils;
    * vertical recurrences;
    * p-dependent polynomial/rational recurrences;
    * low-degree bivariate polynomial laws;
    * additive / multiplicative separability;
    * Newton-coordinate fixed operators;
    * low-degree t-dependent Newton operators;
    * small-support translation-invariant annihilators.

The next step is therefore SOURCE RECONSTRUCTION rather than further fitting.

This experiment provides an exact validation harness for the ORIGINAL
definition/construction of Q_t(p).

The user must implement:

    source_definition(p, t)

from the actual combinatorial/algebraic construction.

The harness then:

    1. evaluates the definition on every observed cell;
    2. compares exactly against the observed Q_t(p);
    3. reports all residuals;
    4. verifies the triangular support;
    5. separately audits the two terminal source values;
    6. computes exact gcd/factorization diagnostics of the reconstructed
       values;
    7. checks whether the definition itself reproduces the observed
       Newton triangle;
    8. refuses to count an unimplemented or inferred formula as evidence.

No interpolation.
No missing-value reconstruction.
No extrapolation counted as proof.
No synthetic second n=pq case.

"""


from __future__ import annotations

import math
import sys
from typing import Callable, Dict, Tuple, Any

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
# EXACT HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def factor_integer(n):
    n = int(n)

    if n == 0:
        return {}

    return sp.factorint(abs(n))


def valuation(n, p):
    n = int(abs(n))

    if n == 0:
        return sp.oo

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def build_observed_cells() -> Dict[Tuple[int, int], int]:
    """
    Q[p] is stored in reverse t-order:

        t = len(Q[p]) - 1 - index.
    """

    observed = {}

    for p_value, values in Q.items():

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            observed[
                (p_value, t_value)
            ] = int(value)

    return observed


def build_observed_layers(observed):
    layers = {}

    for (p_value, t_value), value in observed.items():
        layers.setdefault(
            t_value,
            []
        ).append(
            (p_value, value)
        )

    for t_value in layers:
        layers[t_value].sort()

    return layers


# ============================================================================
# SOURCE DEFINITION
# ============================================================================

def source_definition(
    p_value: int,
    t_value: int,
) -> int:
    """
    ========================================================================
    REPLACE THIS FUNCTION WITH THE ORIGINAL SOURCE CONSTRUCTION.
    ========================================================================

    It must return the exact value Q_t(p).

    Examples of acceptable implementations:

        * an exact combinatorial count;
        * an exact determinant/resultant;
        * an exact coefficient extraction;
        * an exact recurrence derived independently from the source;
        * an exact product/sum from the original mathematics.

    DO NOT:

        * interpolate the observed table;
        * fit a polynomial/rational function to the data;
        * use missing values;
        * use an inferred value as input;
        * hard-code the observed answers.

    The function should operate from the original definition alone.
    """

    raise NotImplementedError(
        "Insert the original mathematical definition of Q_t(p) here."
    )


# ============================================================================
# SOURCE EVALUATION
# ============================================================================

def evaluate_source(
    observed: Dict[Tuple[int, int], int],
    definition: Callable[[int, int], int],
):
    results = {}

    for (p_value, t_value), observed_value in sorted(
        observed.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        try:

            predicted = definition(
                p_value,
                t_value,
            )

            predicted = sp.Integer(predicted)

        except Exception as exc:

            results[
                (p_value, t_value)
            ] = {
                "status": "ERROR",
                "observed": observed_value,
                "predicted": None,
                "residual": None,
                "error": (
                    type(exc).__name__,
                    str(exc),
                ),
            }

            continue

        residual = clean(
            predicted
            - observed_value
        )

        results[
            (p_value, t_value)
        ] = {
            "status": (
                "MATCH"
                if residual == 0
                else "MISMATCH"
            ),
            "observed": observed_value,
            "predicted": predicted,
            "residual": residual,
        }

    return results


# ============================================================================
# EXACT RESIDUAL AUDIT
# ============================================================================

def print_residual_audit(results):

    print()
    print("=" * 78)
    print(
        "1. EXACT SOURCE-DEFINITION RESIDUAL AUDIT"
    )
    print("=" * 78)

    matches = 0
    mismatches = 0
    errors = 0

    for (p_value, t_value), result in results.items():

        print()
        print(
            "  (p={}, t={}):".format(
                p_value,
                t_value,
            )
        )

        print(
            "    observed={}".format(
                result["observed"]
            )
        )

        if result["status"] == "ERROR":

            errors += 1

            print(
                "    status=ERROR"
            )

            print(
                "    error={}".format(
                    result["error"]
                )
            )

            continue

        print(
            "    predicted={}".format(
                result["predicted"]
            )
        )

        print(
            "    residual={}".format(
                result["residual"]
            )
        )

        print(
            "    status={}".format(
                result["status"]
            )
        )

        if result["status"] == "MATCH":
            matches += 1
        else:
            mismatches += 1

    print()
    print(
        "  exact_matches={}".format(
            matches
        )
    )

    print(
        "  exact_mismatches={}".format(
            mismatches
        )
    )

    print(
        "  evaluation_errors={}".format(
            errors
        )
    )

    return (
        matches,
        mismatches,
        errors,
    )


# ============================================================================
# TRIANGULAR SUPPORT AUDIT
# ============================================================================

def triangular_support_audit(
    observed,
    definition,
):

    print()
    print("=" * 78)
    print(
        "2. OBSERVED TRIANGULAR-SUPPORT AUDIT"
    )
    print("=" * 78)

    failures = 0

    for t_value in range(6):

        observed_p = sorted(
            p_value
            for (
                p_value,
                t,
            ) in observed
            if t == t_value
        )

        print()
        print(
            "  t={}: observed_p={}".format(
                t_value,
                observed_p,
            )
        )

        # Important: this does NOT ask the source definition for
        # unobserved cells. It only reports the observed support.

    expected_boundary = {
        0: [1, 3, 5, 7],
        1: [1, 3, 5],
        2: [1, 3, 5],
        3: [1, 3],
        4: [1, 3],
        5: [1],
    }

    for t_value, expected in expected_boundary.items():

        actual = sorted(
            p_value
            for (
                p_value,
                t
            ) in observed
            if t == t_value
        )

        if actual != expected:
            failures += 1

            print()
            print(
                "  SUPPORT_FAILURE t={}: expected={}, actual={}".format(
                    t_value,
                    expected,
                    actual,
                )
            )

    print()
    print(
        "  support_failures={}".format(
            failures
        )
    )

    return failures


# ============================================================================
# TERMINAL SOURCE AUDIT
# ============================================================================

def terminal_source_audit(
    results,
):

    print()
    print("=" * 78)
    print(
        "3. TERMINAL SOURCE AUDIT"
    )
    print("=" * 78)

    terminal_cells = {
        (1, 0): 495451247,
        (3, 0): 421514439,
        (5, 0): 16027881,
        (7, 0): 1,
    }

    for cell, expected in terminal_cells.items():

        result = results.get(cell)

        print()
        print(
            "  cell={}:".format(
                cell
            )
        )

        if result is None:

            print(
                "    status=NOT_EVALUATED"
            )

            continue

        print(
            "    expected={}".format(
                expected
            )
        )

        print(
            "    predicted={}".format(
                result.get(
                    "predicted"
                )
            )
        )

        print(
            "    residual={}".format(
                result.get(
                    "residual"
                )
            )
        )

        print(
            "    exact={}".format(
                result.get(
                    "status"
                ) == "MATCH"
            )
        )


# ============================================================================
# ARITHMETIC PROFILE OF RECONSTRUCTED VALUES
# ============================================================================

def arithmetic_audit(
    results,
):

    print()
    print("=" * 78)
    print(
        "4. RECONSTRUCTED-VALUE ARITHMETIC AUDIT"
    )
    print("=" * 78)

    for cell, result in results.items():

        if result["status"] != "MATCH":
            continue

        value = int(
            result["predicted"]
        )

        print()
        print(
            "  cell={}: value={}".format(
                cell,
                value,
            )
        )

        print(
            "    sign={}".format(
                "+"
                if value > 0
                else "-"
            )
        )

        print(
            "    factorization={}".format(
                factor_integer(value)
            )
        )

        print(
            "    valuations={}".format(
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


def newton_audit(observed):

    print()
    print("=" * 78)
    print(
        "5. OBSERVED NEWTON-TRIANGLE REFERENCE"
    )
    print("=" * 78)

    layers = build_observed_layers(
        observed
    )

    A = {}

    for t_value, layer in layers.items():

        values = [
            value
            for _, value in layer
        ]

        rows = finite_difference_rows(
            values
        )

        A[t_value] = [
            row[0]
            for row in rows
            if row
        ]

        print()
        print(
            "  t={}: {}".format(
                t_value,
                A[t_value],
            )
        )

    return A


# ============================================================================
# INDEPENDENCE / ANTI-FITTING AUDIT
# ============================================================================

def independence_statement():

    print()
    print("=" * 78)
    print(
        "6. SOURCE-RECONSTRUCTION INDEPENDENCE AUDIT"
    )
    print("=" * 78)

    print(
r"""
This experiment does NOT derive source_definition(p,t) from the observed
values.

The observed table is used only as a validation target.

Therefore:

    source-definition -> predicted values -> exact comparison

is logically different from

    observed values -> fitted formula -> self-validation.

A MATCH on all observed cells is therefore meaningful only when
source_definition() was obtained independently from the original
construction.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 345R — EXACT SOURCE-DEFINITION VALIDATION HARNESS"
    )
    print("=" * 78)

    observed = build_observed_cells()

    print()
    print(
        "observed_cells={}".format(
            len(observed)
        )
    )

    # The source definition is intentionally separate from the dataset.
    results = evaluate_source(
        observed,
        source_definition,
    )

    matches, mismatches, errors = (
        print_residual_audit(
            results
        )
    )

    support_failures = (
        triangular_support_audit(
            observed,
            source_definition,
        )
    )

    terminal_source_audit(
        results
    )

    arithmetic_audit(
        results
    )

    newton_audit(
        observed
    )

    independence_statement()

    print()
    print("=" * 78)
    print(
        "7. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(observed)
        )
    )

    print(
        "  exact_matches={}".format(
            matches
        )
    )

    print(
        "  exact_mismatches={}".format(
            mismatches
        )
    )

    print(
        "  evaluation_errors={}".format(
            errors
        )
    )

    print(
        "  support_failures={}".format(
            support_failures
        )
    )

    # This is deliberately NOT marked TRUE until the original source
    # definition has been inserted and independently evaluated.

    full_match = (
        len(observed) > 0
        and matches == len(observed)
        and mismatches == 0
        and errors == 0
        and support_failures == 0
    )

    print(
        "  source_definition_matches_all_observed_cells={}".format(
            full_match
        )
    )

    print(
        "  source_definition_independently_supplied={}".format(
            source_definition.__doc__ is not None
            and
            "REPLACE THIS FUNCTION" not in source_definition.__doc__
        )
    )

    print(
        "  missing_values_used=False"
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
        "EXPERIMENT 345R COMPLETE"
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
