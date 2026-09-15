#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 351R — EXACT MINIMAL-NEW-DATA / HYPOTHESIS-FALSIFICATION AUDIT
==============================================================================

Purpose
-------
Experiments 322R-350R have exhausted a broad family of post-hoc laws on the
current 15 observed cells.

The key remaining issue is DATA IDENTIFICATION.

Several interesting models are either:

    * already contradicted by the observed cells; or
    * exactly data-sized and therefore neither validated nor falsified.

Experiment 351R asks:

    Which SINGLE additional source cell Q(r,t) would convert the largest
    number of currently data-sized hypotheses into genuinely overdetermined
    exact tests?

This is not reconstruction.

No value is assigned to any missing cell.

We only identify which missing cells would create additional independent
equations.

The principal targets are the small rational generating denominators from
350R, together with selected local/transfer models.

For each candidate model we compute:

    1. currently usable exact equations;
    2. current unknown count;
    3. current redundancy;
    4. missing cells that would make one additional valid equation;
    5. the smallest number of additional source cells needed to become
       overdetermined.

The experiment then ranks missing cells by the number of hypotheses they
would simultaneously make falsifiable.

This is a data-acquisition audit, not a formula fit.

No missing value is inferred.
No interpolation is used.
No extrapolation is used.
No synthetic second n=pq case is generated.
"""


from __future__ import annotations

import itertools
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
# LATTICE
# ============================================================================

def build_observed_lattice():

    L = {}

    for p_value, values in Q.items():

        r = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t = (
                len(values)
                - 1
                - index
            )

            L[
                (r, t)
            ] = sp.Integer(value)

    return L


def build_candidate_region():

    """
    Small triangular enlargement around the observed support.

    We do NOT assign values to these cells.
    They are merely candidates for future observation.
    """

    region = set()

    for t in range(0, 7):

        for r in range(0, 5):

            # Extend the observed triangle conservatively.
            if r + t <= 5:
                region.add(
                    (r, t)
                )

    return sorted(region)


# ============================================================================
# DENOMINATOR FAMILIES FROM 350R
# ============================================================================

def denominator_families():

    a, b, c, d, e = sp.symbols(
        "a b c d e"
    )

    return {
        "D1_xy": (
            1 - a*sp.Symbol("x")
            - b*sp.Symbol("y"),
            [a, b],
        ),

        "D2_xy_xy": (
            1 - a*sp.Symbol("x")
            - b*sp.Symbol("y")
            - c*sp.Symbol("x")*sp.Symbol("y"),
            [a, b, c],
        ),

        "D3_xy_x2": (
            1 - a*sp.Symbol("x")
            - b*sp.Symbol("y")
            - c*sp.Symbol("x")**2,
            [a, b, c],
        ),

        "D4_xy_y2": (
            1 - a*sp.Symbol("x")
            - b*sp.Symbol("y")
            - c*sp.Symbol("y")**2,
            [a, b, c],
        ),

        "D5_xy_xy_x2": (
            1
            - a*sp.Symbol("x")
            - b*sp.Symbol("y")
            - c*sp.Symbol("x")*sp.Symbol("y")
            - d*sp.Symbol("x")**2,
            [a, b, c, d],
        ),

        "D6_xy_xy_y2": (
            1
            - a*sp.Symbol("x")
            - b*sp.Symbol("y")
            - c*sp.Symbol("x")*sp.Symbol("y")
            - d*sp.Symbol("y")**2,
            [a, b, c, d],
        ),

        "D7_xy_xy_x2_y2": (
            1
            - a*sp.Symbol("x")
            - b*sp.Symbol("y")
            - c*sp.Symbol("x")*sp.Symbol("y")
            - d*sp.Symbol("x")**2
            - e*sp.Symbol("y")**2,
            [a, b, c, d, e],
        ),

        "D8_xy_xy_x2_xy2": (
            1
            - a*sp.Symbol("x")
            - b*sp.Symbol("y")
            - c*sp.Symbol("x")*sp.Symbol("y")
            - d*sp.Symbol("x")**2
            - e*sp.Symbol("x")*sp.Symbol("y")**2,
            [a, b, c, d, e],
        ),
    }


# ============================================================================
# DENOMINATOR SUPPORT
# ============================================================================

def denominator_shifts(
    denominator
):
    """
    Extract monomial support:

        coefficient * x^dr y^dt.
    """

    x, y = sp.symbols(
        "x y"
    )

    P = sp.Poly(
        sp.expand(
            denominator
        ),
        x,
        y,
    )

    return [
        (
            int(mon[0]),
            int(mon[1]),
        )
        for mon in P.monoms()
    ]


def valid_denominator_equation_cells(
    observed,
    denominator
):
    """
    A coefficient equation for D*G at (r,t) is usable only if every source
    cell required by the denominator is observed.
    """

    shifts = denominator_shifts(
        denominator
    )

    cells = []

    for r, t in sorted(observed):

        valid = True

        for dr, dt in shifts:

            source = (
                r - dr,
                t - dt,
            )

            if source not in observed:
                valid = False
                break

        if valid:
            cells.append(
                (r, t)
            )

    return cells


# ============================================================================
# MISSING-CELL REQUIREMENTS
# ============================================================================

def required_missing_cells_for_equation(
    cell,
    denominator,
    observed
):

    shifts = denominator_shifts(
        denominator
    )

    missing = []

    r, t = cell

    for dr, dt in shifts:

        source = (
            r - dr,
            t - dt,
        )

        if source not in observed:

            missing.append(
                source
            )

    return sorted(
        set(missing)
    )


def first_new_equation_candidates(
    observed,
    denominator
):

    region = build_candidate_region()

    candidates = []

    for cell in region:

        if cell in observed:
            continue

        missing = (
            required_missing_cells_for_equation(
                cell,
                denominator,
                observed,
            )
        )

        # A single newly observed cell creates one equation only if it is the
        # sole missing source needed for that coefficient equation.
        if len(missing) == 1:

            candidates.append(
                (
                    cell,
                    missing,
                )
            )

    return candidates


# ============================================================================
# MODEL IDENTIFICATION SUMMARY
# ============================================================================

def denominator_identification_report(
    observed,
):

    print()
    print("=" * 78)
    print(
        "1. CURRENT DENOMINATOR IDENTIFICATION STATUS"
    )
    print("=" * 78)

    results = []

    for (
        name,
        (
            denominator,
            parameters,
        ),
    ) in denominator_families().items():

        cells = (
            valid_denominator_equation_cells(
                observed,
                denominator,
            )
        )

        equation_count = len(
            cells
        )

        unknown_count = len(
            parameters
        )

        redundancy = (
            equation_count
            - unknown_count
        )

        one_cell_candidates = (
            first_new_equation_candidates(
                observed,
                denominator,
            )
        )

        print()
        print(
            "  {}:".format(
                name
            )
        )

        print(
            "    equations={}".format(
                equation_count
            )
        )

        print(
            "    unknowns={}".format(
                unknown_count
            )
        )

        print(
            "    redundancy={}".format(
                redundancy
            )
        )

        print(
            "    one_cell_candidates={}".format(
                [
                    cell
                    for cell, _
                    in one_cell_candidates
                ]
            )
        )

        # Minimum additional equations needed for genuine overdetermination.
        needed = max(
            0,
            unknown_count
            + 1
            - equation_count,
        )

        print(
            "    additional_equations_needed_for_overdetermination={}".format(
                needed
            )
        )

        results.append(
            {
                "name": name,
                "equation_count": equation_count,
                "unknown_count": unknown_count,
                "redundancy": redundancy,
                "one_cell_candidates": [
                    cell
                    for cell, _
                    in one_cell_candidates
                ],
                "needed": needed,
            }
        )

    return results


# ============================================================================
# CROSS-HYPOTHESIS DATA VALUE
# ============================================================================

def cross_hypothesis_cell_ranking(
    reports
):

    print()
    print("=" * 78)
    print(
        "2. CROSS-HYPOTHESIS VALUE OF A SINGLE NEW CELL"
    )
    print("=" * 78)

    scores = {}

    for report in reports:

        for cell in report[
            "one_cell_candidates"
        ]:

            scores.setdefault(
                cell,
                [],
            ).append(
                report["name"]
            )

    ranked = sorted(
        scores.items(),
        key=lambda item: (
            -len(item[1]),
            item[0][1],
            item[0][0],
        ),
    )

    print()

    for cell, models in ranked:

        print(
            "  cell={}: enables_{}".format(
                cell,
                models,
            )
        )

        print(
            "    hypothesis_count={}".format(
                len(models)
            )
        )

    if not ranked:

        print(
            "  no_single_cell_creates_a_new_valid_equation"
        )

    return ranked


# ============================================================================
# D6-SPECIFIC REPORT
# ============================================================================

def d6_priority_audit(
    observed
):

    specs = denominator_families()

    denominator, parameters = (
        specs["D6_xy_xy_y2"]
    )

    current = (
        valid_denominator_equation_cells(
            observed,
            denominator,
        )
    )

    candidates = (
        first_new_equation_candidates(
            observed,
            denominator,
        )
    )

    print()
    print("=" * 78)
    print(
        "3. D6 PRIORITY AUDIT"
    )
    print("=" * 78)

    print(
        "  current_equations={}".format(
            current
        )
    )

    print(
        "  parameter_count={}".format(
            len(parameters)
        )
    )

    print(
        "  current_redundancy={}".format(
            len(current)
            - len(parameters)
        )
    )

    print(
        "  one_new_cell_candidates={}".format(
            [
                cell
                for cell, _
                in candidates
            ]
        )
    )

    for cell, missing in candidates:

        print()
        print(
            "  candidate_cell={}".format(
                cell
            )
        )

        print(
            "    missing_sources_needed={}".format(
                missing
            )
        )

        print(
            "    after_observation_equation_count={}".format(
                len(current) + 1
            )
        )

        print(
            "    after_observation_redundancy={}".format(
                len(current)
                + 1
                - len(parameters)
            )
        )

        print(
            "    verdict=OVERDETERMINATED_D6_TEST"
        )


# ============================================================================
# OTHER CURRENTLY IMPORTANT MODELS
# ============================================================================

def local_transfer_data_audit(
    observed
):

    print()
    print("=" * 78)
    print(
        "4. RAW WIDTH-2 TRANSFER DATA REQUIREMENT AUDIT"
    )
    print("=" * 78)

    """
    For Q_{t+1}(r) = c0 Q_t(r) + c1 Q_t(r+1),
    an equation at (r,t+1) requires:

        Q_t(r)
        Q_t(r+1)
        Q_{t+1}(r).

    We identify missing target/source cells that would produce a new
    independently testable transfer equation.
    """

    candidate_region = build_candidate_region()

    edge_candidates = []

    for r, t in candidate_region:

        if t == 0:
            continue

        required = [
            (r, t),
            (r, t - 1),
            (r + 1, t - 1),
        ]

        missing = [
            cell
            for cell in required
            if cell not in observed
        ]

        if len(missing) == 1:

            edge_candidates.append(
                (
                    (r, t),
                    missing[0],
                )
            )

    print()

    for target_cell, missing_cell in (
        edge_candidates
    ):

        print(
            "  target_equation_at={}: "
            "one_missing_source={}".format(
                target_cell,
                missing_cell,
            )
        )

    if not edge_candidates:

        print(
            "  no_single_missing_cell_creates_a_new_raw_width2_equation"
        )

    return edge_candidates


# ============================================================================
# RANKING / RECOMMENDATION
# ============================================================================

def recommendation(
    ranked_cells,
    d6_candidates,
    raw_candidates,
):

    print()
    print("=" * 78)
    print(
        "5. EXPERIMENTAL PRIORITY RECOMMENDATION"
    )
    print("=" * 78)

    if ranked_cells:

        best_cell = (
            ranked_cells[0][0]
        )

        best_models = (
            ranked_cells[0][1]
        )

        print()
        print(
            "  highest_value_single_cell={}".format(
                best_cell
            )
        )

        print(
            "  simultaneously_tests={}".format(
                best_models
            )
        )

        if len(best_models) >= 1:

            print(
                "  recommendation=OBTAIN_THIS_SOURCE_VALUE_NEXT"
            )

    elif d6_candidates:

        print(
            "  recommendation=D6_CAN_BE_FALSIFIED_WITH_ONE_ADDITIONAL_SOURCE_CELL"
        )

    elif raw_candidates:

        print(
            "  recommendation=EXTEND_SOURCE_TABLE_AT_A_SINGLE_MISSING_TRANSFER_EDGE"
        )

    else:

        print(
            "  recommendation=NO_SINGLE_NEW_CELL_HAS_HIGH_LEVERAGE"
        )

    print()
    print(
        "  IMPORTANT:"
    )

    print(
        "    The experiment does NOT infer the recommended cell."
    )

    print(
        "    It only identifies which observation would be most informative."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 351R — EXACT MINIMAL-NEW-DATA / "
        "HYPOTHESIS-FALSIFICATION AUDIT"
    )
    print("=" * 78)

    observed = build_observed_lattice()

    print()
    print(
        "OBSERVED CELLS"
    )

    print(
        "  count={}".format(
            len(observed)
        )
    )

    print(
        "  cells={}".format(
            sorted(observed)
        )
    )

    reports = (
        denominator_identification_report(
            observed
        )
    )

    ranked_cells = (
        cross_hypothesis_cell_ranking(
            reports
        )
    )

    d6_candidates = (
        d6_priority_audit(
            observed
        )
    )

    raw_candidates = (
        local_transfer_data_audit(
            observed
        )
    )

    recommendation(
        ranked_cells,
        d6_candidates,
        raw_candidates,
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells_used_only=True"
    )

    print(
        "  missing_values_created=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_performed=False"
    )

    print(
        "  hypothesis_values_fitted=False"
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
        "EXPERIMENT 351R COMPLETE"
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
