#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 366R-FIXED — EXACT BASE-LATTICE RECONSTRUCTION OF TRANSFORMED
                         INTEGER ANNIHILATORS / SUPPORT-SAFE Z AUDIT
==============================================================================

Purpose
-------
Reconstruct the finite-difference discoveries from Experiment 363R directly
from the ORIGINAL triangular source lattice.

For

    H_ab = Delta_r^a Delta_t^b Q

the transform is expanded directly as

    H_ab(r,t)
      = sum_i sum_j
          (-1)^(a-i+b-j)
          C(a,i) C(b,j)
          Q(r+i,t+j).

No transformed lattice is recursively constructed.

A transformed value exists only when EVERY underlying original source cell
needed by that difference exists.

The experiment then:

    1. preserves Z = Q_3(5) symbolically;
    2. identifies support-safe transformed windows;
    3. computes nullspaces from fully observed transformed windows;
    4. evaluates symbolic windows directly from base data;
    5. extracts Z from the ACTUAL symbolic equation;
    6. verifies the resulting Z exactly;
    7. compares discoveries across transforms;
    8. compares them with the previously reported 363R values.

No missing value is inserted numerically.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict

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

Z = sp.Symbol("Z")


# ============================================================================
# TRANSFORMS
# ============================================================================

TRANSFORMS = {
    "H00": (0, 0),
    "H10": (1, 0),
    "H01": (0, 1),
    "H11": (1, 1),
    "H20": (2, 0),
    "H02": (0, 2),
    "H21": (2, 1),
    "H12": (1, 2),
    "H22": (2, 2),
}


STENCILS = {
    "horizontal_width3": [
        (0, 0),
        (1, 0),
        (2, 0),
    ],
    "vertical_width3": [
        (0, 0),
        (0, 1),
        (0, 2),
    ],
    "rectangle_2x2": [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
    ],
    "transport_3": [
        (0, 0),
        (0, 1),
        (1, 1),
    ],
    "lower_transport_3": [
        (0, 0),
        (1, 0),
        (1, 1),
    ],
    "diamond_5": [
        (0, 0),
        (1, 0),
        (0, 1),
        (0, 2),
        (1, 1),
    ],
}


DISCOVERY_KEYS = [
    ("H10", "transport_3"),
    ("H10", "lower_transport_3"),
    ("H01", "horizontal_width3"),
    ("H01", "diamond_5"),
    ("H11", "vertical_width3"),
    ("H02", "transport_3"),
    ("H02", "lower_transport_3"),
]


# ============================================================================
# 363R REPORTED VALUES — COMPARISON ONLY
# ============================================================================

REPORTED_363R = {
    ("H10", "transport_3"):
        sp.Rational(
            -3191330642348143851622391184,
            1320984232593184369,
        ),

    ("H10", "lower_transport_3"):
        sp.Rational(
            -392992110572009658827397946,
            240269201882469729,
        ),

    ("H01", "horizontal_width3"):
        sp.Rational(
            -50926698968037761526399730,
            1750377005640251053,
        ),

    ("H01", "diamond_5"):
        sp.Rational(
            732126390015674197979598697500986324215505626,
            5342832968590334980318521932388509869,
        ),

    ("H11", "vertical_width3"):
        sp.Rational(
            -42367633024449172779805183330,
            12359557018829527669,
        ),

    ("H02", "transport_3"):
        sp.Rational(
            -1636649764753763519072419014,
            79729788761536106653,
        ),

    ("H02", "lower_transport_3"):
        sp.Rational(
            -42490000284788818260262826,
            4210540774600727807,
        ),
}


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


def integer_gcd(values):
    result = 0

    for value in values:
        result = math.gcd(
            result,
            abs(int(value)),
        )

    return result


def primitive_integer_vector(values):
    integers = [
        int(sp.Integer(value))
        for value in values
    ]

    gcd_value = integer_gcd(
        integers
    )

    if gcd_value != 0:
        integers = [
            value // gcd_value
            for value in integers
        ]

    for value in integers:
        if value != 0:
            if value < 0:
                integers = [
                    -v
                    for v in integers
                ]
            break

    return integers


def rational_is_integer(value):
    rational = sp.Rational(value)
    return rational.q == 1


def build_base_lattice(
    include_Z=False,
):
    lattice = {}

    for p_value, values in Q.items():

        r_value = (
            p_value - 1
        ) // 2

        for index, value in enumerate(
            values
        ):

            t_value = (
                len(values)
                - 1
                - index
            )

            lattice[
                (r_value, t_value)
            ] = sp.Integer(value)

    if include_Z:
        lattice[
            (2, 3)
        ] = Z

    return lattice


# ============================================================================
# DIRECT FINITE-DIFFERENCE EXPANSION
# ============================================================================

def difference_terms(
    r,
    t,
    order_r,
    order_t,
):
    terms = []

    for i in range(
        order_r + 1
    ):

        for j in range(
            order_t + 1
        ):

            coefficient = (
                (-1) ** (
                    order_r - i
                    + order_t - j
                )
                * math.comb(
                    order_r,
                    i,
                )
                * math.comb(
                    order_t,
                    j,
                )
            )

            terms.append(
                (
                    coefficient,
                    (
                        r + i,
                        t + j,
                    ),
                )
            )

    return terms


def transformed_value(
    lattice,
    r,
    t,
    order_r,
    order_t,
):
    terms = difference_terms(
        r,
        t,
        order_r,
        order_t,
    )

    for _, cell in terms:
        if cell not in lattice:
            return None

    return clean(
        sum(
            coefficient * lattice[cell]
            for coefficient, cell in terms
        )
    )


def transformed_support(
    lattice,
    order_r,
    order_t,
):
    if not lattice:
        return set()

    max_r = max(
        r
        for r, _ in lattice
    )

    max_t = max(
        t
        for _, t in lattice
    )

    support = set()

    for r in range(
        max_r + 1
    ):

        for t in range(
            max_t + 1
        ):

            value = transformed_value(
                lattice,
                r,
                t,
                order_r,
                order_t,
            )

            if value is not None:
                support.add(
                    (r, t)
                )

    return support


def transformed_windows(
    support,
    offsets,
):
    if not support:
        return []

    max_r = max(
        r
        for r, _ in support
    )

    max_t = max(
        t
        for _, t in support
    )

    max_dr = max(
        dr
        for dr, _ in offsets
    )

    max_dt = max(
        dt
        for _, dt in offsets
    )

    windows = []

    for r0 in range(
        max_r - max_dr + 1
    ):

        for t0 in range(
            max_t - max_dt + 1
        ):

            valid = all(
                (
                    r0 + dr,
                    t0 + dt,
                ) in support
                for dr, dt in offsets
            )

            if valid:
                windows.append(
                    (r0, t0)
                )

    return windows


# ============================================================================
# OBSERVED / SYMBOLIC VALUES
# ============================================================================

def observed_transformed_value(
    observed_base,
    r,
    t,
    order_r,
    order_t,
):
    value = transformed_value(
        observed_base,
        r,
        t,
        order_r,
        order_t,
    )

    if value is None:
        return None

    return sp.Integer(value)


def symbolic_transformed_value(
    symbolic_base,
    r,
    t,
    order_r,
    order_t,
):
    return transformed_value(
        symbolic_base,
        r,
        t,
        order_r,
        order_t,
    )


def window_values(
    symbolic_base,
    window,
    offsets,
    order_r,
    order_t,
):
    r0, t0 = window

    values = []

    for dr, dt in offsets:

        value = symbolic_transformed_value(
            symbolic_base,
            r0 + dr,
            t0 + dt,
            order_r,
            order_t,
        )

        values.append(
            clean(value)
        )

    return values


# ============================================================================
# OBSERVED NULLSPACE
# ============================================================================

def observed_matrix(
    observed_base,
    windows,
    offsets,
    order_r,
    order_t,
):
    rows = []
    usable = []

    for window in windows:

        values = []

        valid = True

        for dr, dt in offsets:

            value = observed_transformed_value(
                observed_base,
                window[0] + dr,
                window[1] + dt,
                order_r,
                order_t,
            )

            if value is None:
                valid = False
                break

            values.append(
                value
            )

        if valid:
            rows.append(values)
            usable.append(window)

    if not rows:
        return (
            sp.zeros(0, len(offsets)),
            usable,
        )

    return (
        sp.Matrix(rows),
        usable,
    )


def primitive_nullspace_vectors(
    matrix,
):
    if matrix.rows == 0:
        return []

    vectors = matrix.nullspace()

    results = []

    for vector in vectors:

        denominator_lcm = 1

        for value in vector:

            denominator_lcm = sp.ilcm(
                denominator_lcm,
                int(
                    sp.denom(
                        sp.Rational(
                            value
                        )
                    )
                ),
            )

        integers = [
            int(
                sp.Rational(value)
                * denominator_lcm
            )
            for value in vector
        ]

        results.append(
            primitive_integer_vector(
                integers
            )
        )

    return results


# ============================================================================
# FORCE Z FROM ACTUAL EQUATION
# ============================================================================

def force_Z_from_equation(
    values,
    relation,
):
    equation = clean(
        sum(
            sp.Integer(
                relation[index]
            )
            * values[index]
            for index in range(
                len(values)
            )
        )
    )

    expanded = sp.expand(
        equation
    )

    coefficient_Z = clean(
        expanded.coeff(
            Z
        )
    )

    constant_part = clean(
        expanded.subs(
            Z,
            0,
        )
    )

    if coefficient_Z == 0:
        return None

    forced_Z = clean(
        -constant_part
        / coefficient_Z
    )

    residual = clean(
        equation.subs(
            Z,
            forced_Z,
        )
    )

    return {
        "equation": equation,
        "coefficient_Z": coefficient_Z,
        "constant_part": constant_part,
        "forced_Z": forced_Z,
        "residual": residual,
    }


# ============================================================================
# ONE DISCOVERY RECONSTRUCTION
# ============================================================================

def reconstruct_discoveries(
    transform_name,
    stencil_name,
):
    order_r, order_t = TRANSFORMS[
        transform_name
    ]

    offsets = STENCILS[
        stencil_name
    ]

    observed_base = build_base_lattice(
        include_Z=False
    )

    symbolic_base = build_base_lattice(
        include_Z=True
    )

    observed_support = transformed_support(
        observed_base,
        order_r,
        order_t,
    )

    symbolic_support = transformed_support(
        symbolic_base,
        order_r,
        order_t,
    )

    all_windows = transformed_windows(
        symbolic_support,
        offsets,
    )

    matrix, fully_observed_windows = (
        observed_matrix(
            observed_base,
            all_windows,
            offsets,
            order_r,
            order_t,
        )
    )

    symbolic_windows = [
        window
        for window in all_windows
        if window not in fully_observed_windows
    ]

    nullspace = primitive_nullspace_vectors(
        matrix
    )

    discoveries = []

    for relation in nullspace:

        for window in symbolic_windows:

            values = window_values(
                symbolic_base,
                window,
                offsets,
                order_r,
                order_t,
            )

            forced = force_Z_from_equation(
                values,
                relation,
            )

            if forced is None:
                continue

            discoveries.append(
                {
                    "window": window,
                    "relation": relation,
                    "values": values,
                    "equation": forced[
                        "equation"
                    ],
                    "coefficient_Z": forced[
                        "coefficient_Z"
                    ],
                    "constant_part": forced[
                        "constant_part"
                    ],
                    "forced_Z": forced[
                        "forced_Z"
                    ],
                    "residual": forced[
                        "residual"
                    ],
                }
            )

    return {
        "transform": transform_name,
        "stencil": stencil_name,
        "order_r": order_r,
        "order_t": order_t,
        "observed_support": observed_support,
        "symbolic_support": symbolic_support,
        "all_windows": all_windows,
        "fully_observed_windows": (
            fully_observed_windows
        ),
        "symbolic_windows": symbolic_windows,
        "nullspace_dimension": len(
            nullspace
        ),
        "nullspace_relations": nullspace,
        "discoveries": discoveries,
    }


# ============================================================================
# RUN AUDIT
# ============================================================================

def run_audit():
    results = {}

    print()
    print("=" * 78)
    print(
        "1. SUPPORT-SAFE DIRECT RECONSTRUCTION"
    )
    print("=" * 78)

    for key in DISCOVERY_KEYS:

        result = reconstruct_discoveries(
            key[0],
            key[1],
        )

        results[key] = result

        print()
        print(
            "  {} / {}".format(
                key[0],
                key[1],
            )
        )

        print(
            "    transform_orders=({}, {})".format(
                result["order_r"],
                result["order_t"],
            )
        )

        print(
            "    observed_support={}".format(
                sorted(
                    result[
                        "observed_support"
                    ]
                )
            )
        )

        print(
            "    symbolic_support={}".format(
                sorted(
                    result[
                        "symbolic_support"
                    ]
                )
            )
        )

        print(
            "    all_windows={}".format(
                result[
                    "all_windows"
                ]
            )
        )

        print(
            "    fully_observed_windows={}".format(
                result[
                    "fully_observed_windows"
                ]
            )
        )

        print(
            "    symbolic_windows={}".format(
                result[
                    "symbolic_windows"
                ]
            )
        )

        print(
            "    nullspace_dimension={}".format(
                result[
                    "nullspace_dimension"
                ]
            )
        )

        print(
            "    nullspace_relations={}".format(
                result[
                    "nullspace_relations"
                ]
            )
        )

        print(
            "    discovery_count={}".format(
                len(
                    result[
                        "discoveries"
                    ]
                )
            )
        )

        for index, discovery in enumerate(
            result["discoveries"]
        ):

            print()
            print(
                "    discovery_{}:".format(
                    index
                )
            )

            print(
                "      window={}".format(
                    discovery[
                        "window"
                    ]
                )
            )

            print(
                "      relation={}".format(
                    discovery[
                        "relation"
                    ]
                )
            )

            print(
                "      transformed_values={}".format(
                    discovery[
                        "values"
                    ]
                )
            )

            print(
                "      equation={}".format(
                    discovery[
                        "equation"
                    ]
                )
            )

            print(
                "      coefficient_of_Z={}".format(
                    discovery[
                        "coefficient_Z"
                    ]
                )
            )

            print(
                "      forced_Z={}".format(
                    discovery[
                        "forced_Z"
                    ]
                )
            )

            print(
                "      residual={}".format(
                    discovery[
                        "residual"
                    ]
                )
            )

    return results


# ============================================================================
# CROSS-DISCOVERY CONSISTENCY
# ============================================================================

def collect_discoveries(
    results,
):
    records = []

    for key in DISCOVERY_KEYS:

        for discovery in results[
            key
        ]["discoveries"]:

            records.append(
                (
                    key,
                    discovery,
                )
            )

    return records


def consistency_audit(
    records,
):
    print()
    print("=" * 78)
    print(
        "2. CROSS-DISCOVERY Z CONSISTENCY"
    )
    print("=" * 78)

    for key, discovery in records:

        print()
        print(
            "  {} / {}:".format(
                key[0],
                key[1],
            )
        )

        print(
            "    window={}".format(
                discovery[
                    "window"
                ]
            )
        )

        print(
            "    forced_Z={}".format(
                discovery[
                    "forced_Z"
                ]
            )
        )

    equal_pairs = []

    for i in range(
        len(records)
    ):

        for j in range(
            i + 1,
            len(records),
        ):

            zi = sp.Rational(
                records[i][1]["forced_Z"]
            )

            zj = sp.Rational(
                records[j][1]["forced_Z"]
            )

            if zi == zj:
                equal_pairs.append(
                    (
                        records[i][0],
                        records[j][0],
                        zi,
                    )
                )

    distinct_Z = len({
        sp.Rational(
            discovery["forced_Z"]
        )
        for _, discovery in records
    })

    print()
    print(
        "  discovery_count={}".format(
            len(records)
        )
    )

    print(
        "  distinct_Z_count={}".format(
            distinct_Z
        )
    )

    print(
        "  equal_pairs={}".format(
            equal_pairs
        )
    )

    return equal_pairs, distinct_Z


# ============================================================================
# RESIDUAL AUDIT
# ============================================================================

def residual_audit(
    records,
):
    print()
    print("=" * 78)
    print(
        "3. EXACT RESIDUAL VERIFICATION"
    )
    print("=" * 78)

    failures = 0

    for key, discovery in records:

        residual = clean(
            discovery[
                "residual"
            ]
        )

        print(
            "  {} / {} window={}: residual={}".format(
                key[0],
                key[1],
                discovery[
                    "window"
                ],
                residual,
            )
        )

        if residual != 0:
            failures += 1

    print()
    print(
        "  residual_failures={}".format(
            failures
        )
    )

    return failures


# ============================================================================
# 363R COMPARISON
# ============================================================================

def compare_363r(
    records,
):
    print()
    print("=" * 78)
    print(
        "4. COMPARISON WITH 363R"
    )
    print("=" * 78)

    mismatches = 0

    for key, discovery in records:

        reconstructed = sp.Rational(
            discovery[
                "forced_Z"
            ]
        )

        reported = REPORTED_363R.get(
            key
        )

        if reported is None:
            continue

        matches = (
            reconstructed
            == reported
        )

        print()
        print(
            "  {} / {}:".format(
                key[0],
                key[1],
            )
        )

        print(
            "    reconstructed={}".format(
                reconstructed
            )
        )

        print(
            "    reported_363R={}".format(
                reported
            )
        )

        print(
            "    matches={}".format(
                matches
            )
        )

        if not matches:
            mismatches += 1

    print()
    print(
        "  reported_value_mismatches={}".format(
            mismatches
        )
    )

    return mismatches


# ============================================================================
# INTEGER AUDIT
# ============================================================================

def integer_audit(
    records,
):
    print()
    print("=" * 78)
    print(
        "5. INTEGER-Z COMPATIBILITY"
    )
    print("=" * 78)

    integer_count = 0

    for key, discovery in records:

        forced_Z = discovery[
            "forced_Z"
        ]

        compatible = rational_is_integer(
            forced_Z
        )

        print(
            "  {} / {}: Z_integer={}".format(
                key[0],
                key[1],
                compatible,
            )
        )

        if compatible:
            integer_count += 1

    print()
    print(
        "  integer_compatible_count={}".format(
            integer_count
        )
    )

    return integer_count


# ============================================================================
# TRANSFORM SUMMARY
# ============================================================================

def transform_summary(
    records,
):
    print()
    print("=" * 78)
    print(
        "6. TRANSFORM-LEVEL CONSISTENCY"
    )
    print("=" * 78)

    grouped = defaultdict(list)

    for key, discovery in records:

        grouped[
            key[0]
        ].append(
            sp.Rational(
                discovery[
                    "forced_Z"
                ]
            )
        )

    inconsistent = []

    for transform_name in sorted(
        grouped
    ):

        values = grouped[
            transform_name
        ]

        distinct = len(
            set(values)
        )

        internally_consistent = (
            distinct == 1
        )

        print()
        print(
            "  {}:".format(
                transform_name
            )
        )

        print(
            "    forced_Z_values={}".format(
                values
            )
        )

        print(
            "    distinct_count={}".format(
                distinct
            )
        )

        print(
            "    internally_consistent={}".format(
                internally_consistent
            )
        )

        if not internally_consistent:
            inconsistent.append(
                transform_name
            )

    return inconsistent


# ============================================================================
# FINAL VERDICT
# ============================================================================

def final_verdict(
    records,
    residual_failures,
    reported_mismatches,
):
    if residual_failures != 0:
        verdict = (
            "IMPLEMENTATION_ERROR_OR_INVALID_WINDOW"
        )

    elif not records:
        verdict = (
            "NO_SUPPORT_SAFE_Z_FORCING_DISCOVERIES"
        )

    elif reported_mismatches != 0:
        verdict = (
            "363R_REPORTED_VALUES_NOT_REPRODUCED"
        )

    else:

        distinct_Z = len({
            sp.Rational(
                discovery[
                    "forced_Z"
                ]
            )
            for _, discovery in records
        })

        if distinct_Z == 1:
            verdict = (
                "COMMON_Z_RECONSTRUCTED"
            )
        else:
            verdict = (
                "CROSS_DISCOVERY_INCONSISTENCY"
            )

    print()
    print("=" * 78)
    print(
        "7. STRUCTURAL VERDICT"
    )
    print("=" * 78)

    print(
        "  reconstructed_discoveries={}".format(
            len(records)
        )
    )

    print(
        "  residual_failures={}".format(
            residual_failures
        )
    )

    print(
        "  reported_363R_mismatches={}".format(
            reported_mismatches
        )
    )

    print(
        "  verdict={}".format(
            verdict
        )
    )

    return verdict


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 366R-FIXED — EXACT BASE-LATTICE RECONSTRUCTION "
        "OF TRANSFORMED INTEGER ANNIHILATORS / SUPPORT-SAFE Z AUDIT"
    )
    print("=" * 78)

    results = run_audit()

    records = collect_discoveries(
        results
    )

    equal_pairs, distinct_Z = (
        consistency_audit(
            records
        )
    )

    residual_failures = residual_audit(
        records
    )

    reported_mismatches = compare_363r(
        records
    )

    integer_count = integer_audit(
        records
    )

    transform_inconsistencies = (
        transform_summary(
            records
        )
    )

    verdict = final_verdict(
        records,
        residual_failures,
        reported_mismatches,
    )

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  direct_base_lattice_transform_expansion=True"
    )

    print(
        "  recursive_transformed_lattice_used=False"
    )

    print(
        "  original_support_checked_for_every_difference_term=True"
    )

    print(
        "  symbolic_Z_preserved=True"
    )

    print(
        "  actual_coefficient_of_Z_computed=True"
    )

    print(
        "  stored_relation_position_assumption=False"
    )

    print(
        "  reconstructed_discoveries={}".format(
            len(records)
        )
    )

    print(
        "  distinct_reconstructed_Z={}".format(
            distinct_Z
        )
    )

    print(
        "  equal_discovery_pairs={}".format(
            len(equal_pairs)
        )
    )

    print(
        "  residual_failures={}".format(
            residual_failures
        )
    )

    print(
        "  reported_363R_mismatches={}".format(
            reported_mismatches
        )
    )

    print(
        "  integer_compatible_discoveries={}".format(
            integer_count
        )
    )

    print(
        "  transform_internal_inconsistencies={}".format(
            transform_inconsistencies
        )
    )

    print(
        "  missing_value_inserted_numerically=False"
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

    # Execution/basic consistency check.
    basic_ok = (
        residual_failures == 0
    )

    print(
        "  failures={}".format(
            residual_failures
        )
    )

    print(
        "  ALL BASIC CHECKS PASS={}".format(
            basic_ok
        )
    )

    print()
    print(
        "EXPERIMENT 366R-FIXED COMPLETE"
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