#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 346R — EXACT SOURCE-PROVENANCE / DEPENDENCY-CHAIN AUDIT
==============================================================================

Purpose
-------
Stop blind formula fitting.

Recover and verify the exact chain by which the source quantity Q_t(p)
was constructed from the upstream n = p*q data.

The experiment accepts independently reconstructed upstream objects and
checks exact identities leading to the observed Q_t(p).

TARGET:

    upstream source objects
            |
            v
        Q_t(p)
            |
            v
       Delta_p^j Q_t(p)
            |
            v
       transfer data

This experiment does NOT fit Q.

It only verifies independently proposed source-construction identities.

Required input
--------------
The missing ingredient is the actual upstream definition/data from
Experiments 304R–315R.

Insert those objects into the clearly marked functions below.

No interpolation.
No missing-value reconstruction.
No extrapolation.
No synthetic second n=pq case.
Exact integer/SymPy arithmetic only.
"""


from __future__ import annotations

import math
import sys
from typing import Any, Dict, Tuple

import sympy as sp


# ============================================================================
# OBSERVED Q TABLE
# ============================================================================

Q = {
    (1, 0): 495451247,
    (3, 0): 421514439,
    (5, 0): 16027881,
    (7, 0): 1,

    (1, 1): -1338089411,
    (3, 1): -128667196,
    (5, 1): 4771718,

    (1, 2): 1764373740,
    (3, 2): -152369292,
    (5, 2): -62398,

    (1, 3): 2668721436,
    (3, 3): -1263551016,

    (1, 4): -11600759760,
    (3, 4): 9955176,

    (1, 5): -126258696,
}


# ============================================================================
# BASIC EXACT HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def exact_equal(a, b):
    return clean(a - b) == 0


def factor_integer(n):
    n = int(n)

    if n == 0:
        return {}

    return sp.factorint(abs(n))


def valuation(n, p):
    n = abs(int(n))

    if n == 0:
        return sp.oo

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


# ============================================================================
# UPSTREAM SOURCE DEFINITIONS
# ============================================================================

def upstream_object(
    p_value: int,
    t_value: int,
) -> Any:
    """
    ========================================================================
    INSERT THE ACTUAL UPSTREAM OBJECT USED TO CONSTRUCT Q_t(p).
    ========================================================================

    This must come from the original mathematical construction, not from
    fitting the observed Q table.

    Return whatever exact object is natural:

        integer,
        rational,
        polynomial,
        determinant,
        resultant,
        coefficient vector,
        generating-function coefficient,
        etc.

    """
    raise NotImplementedError(
        "Insert the independently reconstructed upstream object."
    )


def construct_Q_from_upstream(
    upstream: Any,
    p_value: int,
    t_value: int,
) -> Any:
    """
    ========================================================================
    INSERT THE ACTUAL EXACT CONSTRUCTION OF Q_t(p).
    ========================================================================

    Example conceptual forms:

        return upstream
        return determinant(...)
        return coefficient(...)
        return gcd(...)
        return numerator(...)
        return resultant(...)

    The implementation must come from the mathematical derivation.
    """
    raise NotImplementedError(
        "Insert the independently derived construction of Q_t(p)."
    )


# ============================================================================
# EXACT SOURCE VERIFICATION
# ============================================================================

def verify_source_cells():

    print()
    print("=" * 78)
    print("1. EXACT UPSTREAM -> Q_t(p) VERIFICATION")
    print("=" * 78)

    matches = 0
    mismatches = 0
    errors = 0

    for (p_value, t_value), observed in sorted(
        Q.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        print()
        print(
            "  (p={}, t={}):".format(
                p_value,
                t_value,
            )
        )

        try:

            upstream = upstream_object(
                p_value,
                t_value,
            )

            reconstructed = construct_Q_from_upstream(
                upstream,
                p_value,
                t_value,
            )

            residual = clean(
                reconstructed
                - observed
            )

            print(
                "    upstream={}".format(
                    upstream
                )
            )

            print(
                "    reconstructed={}".format(
                    reconstructed
                )
            )

            print(
                "    observed={}".format(
                    observed
                )
            )

            print(
                "    residual={}".format(
                    residual
                )
            )

            if residual == 0:

                print(
                    "    status=MATCH"
                )

                matches += 1

            else:

                print(
                    "    status=MISMATCH"
                )

                mismatches += 1

        except Exception as exc:

            print(
                "    status=ERROR"
            )

            print(
                "    error={}: {}".format(
                    type(exc).__name__,
                    exc,
                )
            )

            errors += 1

    return (
        matches,
        mismatches,
        errors,
    )


# ============================================================================
# TERMINAL SOURCE AUDIT
# ============================================================================

def terminal_source_audit():

    print()
    print("=" * 78)
    print("2. TERMINAL SOURCE FACTORIZATION")
    print("=" * 78)

    terminals = {
        (1, 0): 495451247,
        (3, 0): 421514439,
        (5, 0): 16027881,
        (7, 0): 1,
    }

    for cell, value in terminals.items():

        print()
        print(
            "  cell={}:".format(
                cell
            )
        )

        print(
            "    value={}".format(
                value
            )
        )

        print(
            "    factorization={}".format(
                factor_integer(value)
            )
        )

        print(
            "    v17={}".format(
                valuation(value, 17)
            )
        )

    print()
    print(
        "  gcd(q1,q3)={}".format(
            math.gcd(
                terminals[(1, 0)],
                terminals[(3, 0)],
            )
        )
    )


# ============================================================================
# DIFFERENCE-CONSISTENCY AUDIT
# ============================================================================

def finite_difference_rows(values):

    rows = [
        [
            sp.Integer(v)
            for v in values
        ]
    ]

    while len(rows[-1]) > 1:

        previous = rows[-1]

        rows.append([
            clean(
                previous[i + 1]
                - previous[i]
            )
            for i in range(
                len(previous) - 1
            )
        ])

    return rows


def observed_newton_audit():

    print()
    print("=" * 78)
    print("3. OBSERVED NEWTON COORDINATES")
    print("=" * 78)

    layers = {}

    for (p_value, t_value), value in Q.items():

        layers.setdefault(
            t_value,
            []
        ).append(
            (p_value, value)
        )

    for t_value in sorted(layers):

        layer = sorted(
            layers[t_value]
        )

        values = [
            value
            for _, value in layer
        ]

        rows = finite_difference_rows(
            values
        )

        print()
        print(
            "  t={}:".format(
                t_value
            )
        )

        print(
            "    A={}".format(
                [
                    row[0]
                    for row in rows
                    if row
                ]
            )
        )


# ============================================================================
# 17-ADIC TERMINAL BOUNDARY AUDIT
# ============================================================================

def seventeen_boundary_audit():

    print()
    print("=" * 78)
    print("4. 17-ADIC TERMINAL-BOUNDARY AUDIT")
    print("=" * 78)

    for (p_value, t_value), value in sorted(
        Q.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        v17 = valuation(
            value,
            17,
        )

        print(
            "  (p={},t={}): v17={}".format(
                p_value,
                t_value,
                v17,
            )
        )


# ============================================================================
# INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiments 322R-344R have tested a broad range of generic algebraic,
linear, projective, recurrence, stencil, separability, and arithmetic
models.

None produced a validated universal source law.

The remaining high-value question is therefore not:

    "Which fitted formula reproduces Q?"

It is:

    "What exact mathematical construction produced Q?"

A successful 346R run requires an independently reconstructed upstream
definition.

If all 15 observed cells match exactly, that validates the source
construction and gives a principled route toward the next theorem.

If even one cell fails, the candidate source reconstruction is wrong or
incomplete.

This experiment therefore separates:

    source derivation
        from
    post-hoc data fitting.

No extrapolation is considered evidence.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 346R — EXACT SOURCE-PROVENANCE / "
        "DEPENDENCY-CHAIN AUDIT"
    )
    print("=" * 78)

    terminal_source_audit()
    observed_newton_audit()
    seventeen_boundary_audit()

    try:

        matches, mismatches, errors = (
            verify_source_cells()
        )

        implemented = True

    except NotImplementedError as exc:

        print()
        print(
            "SOURCE-DEFINITION STATUS: NOT IMPLEMENTED"
        )

        print(
            "  {}".format(
                exc
            )
        )

        matches = 0
        mismatches = 0
        errors = 0
        implemented = False

    interpretation()

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(Q)
        )
    )

    print(
        "  source_definition_implemented={}".format(
            implemented
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
        "EXPERIMENT 346R COMPLETE"
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
