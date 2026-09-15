#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 365R — EXACT RECONSTRUCTION OF TRANSFORMED INTEGER ANNIHILATORS /
                  Z-FORCING PROVENANCE AUDIT
==============================================================================

Purpose
-------
Correct the provenance issue exposed by 364R.

364R compared the previously reported forced-Z values correctly, but its
symbolic residual cross-check incorrectly assumed that the final entry of
every primitive relation multiplies Z.

365R does NOT reuse that assumption.

For every transformed lattice and every discovered stencil:

    1. rebuild the transformed lattice symbolically from the observed source;
    2. identify the actual symbolic-Z-containing window;
    3. build the exact integer homogeneous equation in its true support order;
    4. derive Z directly from that equation;
    5. verify the derived Z by exact substitution;
    6. compare the independently reconstructed Z values across discoveries.

The experiment therefore separates:

    discovery compatibility
        from
    bookkeeping/order errors in stored relation vectors.

No missing value is inserted numerically.
Z remains symbolic throughout calibration.
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
# SOURCE
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
# TRANSFORM DEFINITIONS
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
    g = 0

    for value in values:
        g = math.gcd(
            g,
            abs(int(value)),
        )

    return g


def primitive_vector(values):
    values = [
        int(sp.Integer(v))
        for v in values
    ]

    g = integer_gcd(values)

    if g:
        values = [
            value // g
            for value in values
        ]

    for value in values:
        if value != 0:

            if value < 0:
                values = [
                    -v
                    for v in values
                ]

            break

    return values


def rational_is_integer(value):
    value = sp.Rational(value)
    return value.q == 1


def signed_rational(value):
    return sp.factor(
        sp.Rational(value)
    )


# ============================================================================
# BASE LATTICE
# ============================================================================

def build_base_lattice():
    lattice = {}

    for p_value, values in Q.items():

        r_value = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t_value = (
                len(values) - 1 - index
            )

            lattice[
                (r_value, t_value)
            ] = sp.Integer(value)

    return lattice


# ============================================================================
# SYMBOLIC BASE LATTICE
# ============================================================================

def build_symbolic_base_lattice():
    lattice = build_base_lattice()

    # This is the ONLY symbolic missing source value used.
    lattice[(2, 3)] = Z

    return lattice


# ============================================================================
# FINITE DIFFERENCES
# ============================================================================

def delta_r(lattice):
    result = {}

    for r, t in lattice:

        source = (
            r + 1,
            t,
        )

        if source not in lattice:
            continue

        result[
            (r, t)
        ] = clean(
            lattice[source]
            - lattice[(r, t)]
        )

    return result


def delta_t(lattice):
    result = {}

    for r, t in lattice:

        source = (
            r,
            t + 1,
        )

        if source not in lattice:
            continue

        result[
            (r, t)
        ] = clean(
            lattice[source]
            - lattice[(r, t)]
        )

    return result


def apply_delta_r(
    lattice,
    order,
):
    result = lattice

    for _ in range(order):
        result = delta_r(result)

    return result


def apply_delta_t(
    lattice,
    order,
):
    result = lattice

    for _ in range(order):
        result = delta_t(result)

    return result


def transformed_lattice(
    base,
    r_order,
    t_order,
):
    result = base

    result = apply_delta_r(
        result,
        r_order,
    )

    result = apply_delta_t(
        result,
        t_order,
    )

    return result


# ============================================================================
# STENCIL WINDOWS
# ============================================================================

def all_windows(
    lattice,
    offsets,
):
    """
    Return all translations for which every stencil cell exists.
    """

    if not lattice:
        return []

    min_r = min(
        r for r, _ in lattice
    )
    max_r = max(
        r for r, _ in lattice
    )
    min_t = min(
        t for _, t in lattice
    )
    max_t = max(
        t for _, t in lattice
    )

    max_dr = max(
        dr for dr, _ in offsets
    )
    max_dt = max(
        dt for _, dt in offsets
    )

    windows = []

    for r0 in range(
        min_r,
        max_r - max_dr + 1,
    ):

        for t0 in range(
            min_t,
            max_t - max_dt + 1,
        ):

            cells = [
                (
                    r0 + dr,
                    t0 + dt,
                )
                for dr, dt in offsets
            ]

            if all(
                cell in lattice
                for cell in cells
            ):
                windows.append(
                    (
                        r0,
                        t0,
                    )
                )

    return windows


# ============================================================================
# INTEGER ANNIHILATOR
# ============================================================================

def symbolic_window_coefficients(
    lattice,
    window,
    offsets,
):
    r0, t0 = window

    values = [
        clean(
            lattice[
                (
                    r0 + dr,
                    t0 + dt,
                )
            ]
        )
        for dr, dt in offsets
    ]

    return values


def exact_integer_nullspace_for_observed_windows(
    symbolic_lattice,
    observed_lattice,
    offsets,
):
    """
    Build the homogeneous matrix from windows containing ONLY observed
    transformed entries.

    Since Z is symbolic, any window containing Z is excluded from the
    calibration nullspace.

    The resulting nullspace therefore describes integer relations that
    are already present before using the missing cell.
    """

    observed_windows = []
    symbolic_windows = []

    all_symbolic_windows = all_windows(
        symbolic_lattice,
        offsets,
    )

    for window in all_symbolic_windows:

        values = symbolic_window_coefficients(
            symbolic_lattice,
            window,
            offsets,
        )

        if any(
            sp.sympify(value).has(Z)
            for value in values
        ):
            symbolic_windows.append(
                window
            )
        else:
            observed_windows.append(
                window
            )

    if not observed_windows:
        return {
            "observed_windows": [],
            "symbolic_windows": symbolic_windows,
            "nullspace": [],
        }

    matrix = sp.Matrix([
        [
            sp.Integer(
                observed_lattice[
                    (
                        window[0] + offsets[j][0],
                        window[1] + offsets[j][1],
                    )
                ]
            )
            for j in range(len(offsets))
        ]
        for window in observed_windows
    ])

    nullspace = matrix.nullspace()

    return {
        "observed_windows": observed_windows,
        "symbolic_windows": symbolic_windows,
        "nullspace": nullspace,
    }


# ============================================================================
# SYMBOLIC FORCED-Z FROM ACTUAL WINDOW ORDER
# ============================================================================

def forced_z_from_window(
    lattice,
    window,
    offsets,
    relation,
):
    """
    Given one exact relation and one window, identify which actual stencil
    entry contains Z.

    No positional assumption is made.
    """

    r0, t0 = window

    values = [
        clean(
            lattice[
                (
                    r0 + dr,
                    t0 + dt,
                )
            ]
        )
        for dr, dt in offsets
    ]

    symbolic_positions = [
        index
        for index, value in enumerate(values)
        if sp.sympify(value).has(Z)
    ]

    if len(symbolic_positions) != 1:
        return None

    z_position = symbolic_positions[0]

    coefficient_z = sp.Integer(
        relation[z_position]
    )

    if coefficient_z == 0:
        return None

    constant = sp.Integer(0)

    for index, coefficient in enumerate(
        relation
    ):

        if index == z_position:
            continue

        constant += (
            sp.Integer(coefficient)
            * values[index]
        )

    forced = clean(
        -constant / coefficient_z
    )

    symbolic_residual = clean(
        sum(
            sp.Integer(relation[i])
            * values[i]
            for i in range(len(values))
        )
    )

    return {
        "z_position": z_position,
        "values": values,
        "forced_Z": forced,
        "symbolic_residual": symbolic_residual,
        "coefficient_Z": coefficient_z,
        "constant_part": clean(constant),
    }


# ============================================================================
# DISCOVERY RECONSTRUCTION
# ============================================================================

def reconstruct_discovery(
    transform_name,
    stencil_name,
):
    r_order, t_order = TRANSFORMS[
        transform_name
    ]

    offsets = STENCILS[
        stencil_name
    ]

    base_observed = build_base_lattice()
    base_symbolic = (
        build_symbolic_base_lattice()
    )

    observed_transformed = (
        transformed_lattice(
            base_observed,
            r_order,
            t_order,
        )
    )

    symbolic_transformed = (
        transformed_lattice(
            base_symbolic,
            r_order,
            t_order,
        )
    )

    result = (
        exact_integer_nullspace_for_observed_windows(
            symbolic_transformed,
            observed_transformed,
            offsets,
        )
    )

    candidate_records = []

    for vector in result["nullspace"]:

        denominators = [
            sp.denom(
                sp.Rational(v)
            )
            for v in vector
        ]

        lcm = 1

        for d in denominators:
            lcm = sp.ilcm(
                lcm,
                int(d),
            )

        integer_values = [
            int(
                sp.Rational(v)
                * lcm
            )
            for v in vector
        ]

        relation = primitive_vector(
            integer_values
        )

        for window in result[
            "symbolic_windows"
        ]:

            forced = forced_z_from_window(
                symbolic_transformed,
                window,
                offsets,
                relation,
            )

            if forced is not None:

                candidate_records.append(
                    {
                        "relation": relation,
                        "window": window,
                        **forced,
                    }
                )

    return {
        "transform": transform_name,
        "stencil": stencil_name,
        "observed_windows": result[
            "observed_windows"
        ],
        "symbolic_windows": result[
            "symbolic_windows"
        ],
        "nullspace": result[
            "nullspace"
        ],
        "candidate_records": candidate_records,
    }


# ============================================================================
# 1. RECONSTRUCT ALL SEVEN DISCOVERIES
# ============================================================================

def reconstruction_audit():

    print()
    print("=" * 78)
    print(
        "1. DIRECT RECONSTRUCTION OF 363R DISCOVERIES"
    )
    print("=" * 78)

    results = {}

    for key in DISCOVERY_KEYS:

        transform_name, stencil_name = key

        result = reconstruct_discovery(
            transform_name,
            stencil_name,
        )

        results[key] = result

        print()
        print(
            "  {} / {}".format(
                transform_name,
                stencil_name,
            )
        )

        print(
            "    observed_windows={}".format(
                result[
                    "observed_windows"
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
                len(
                    result["nullspace"]
                )
            )
        )

        print(
            "    symbolic_candidates={}".format(
                len(
                    result[
                        "candidate_records"
                    ]
                )
            )
        )

        for index, candidate in enumerate(
            result["candidate_records"]
        ):

            print()
            print(
                "    candidate_{}:".format(
                    index
                )
            )

            print(
                "      window={}".format(
                    candidate["window"]
                )
            )

            print(
                "      z_position={}".format(
                    candidate["z_position"]
                )
            )

            print(
                "      relation={}".format(
                    candidate["relation"]
                )
            )

            print(
                "      forced_Z={}".format(
                    candidate["forced_Z"]
                )
            )

    return results


# ============================================================================
# 2. CANONICAL CANDIDATE SELECTION
# ============================================================================

def select_canonical_candidates(
    results
):
    """
    Select exactly one symbolic-Z candidate per discovery.

    The preferred candidate is the first one in deterministic order.
    This is NOT a mathematical preference; it merely gives a stable record.

    All candidates are retained for consistency checks.
    """

    selected = {}

    for key, result in results.items():

        candidates = result[
            "candidate_records"
        ]

        if not candidates:
            continue

        selected[key] = candidates[0]

    print()
    print("=" * 78)
    print(
        "2. CANONICAL RECONSTRUCTED Z VALUES"
    )
    print("=" * 78)

    for key in DISCOVERY_KEYS:

        if key not in selected:
            print()
            print(
                "  {} / {}: NO_RECONSTRUCTED_CANDIDATE".format(
                    key[0],
                    key[1],
                )
            )
            continue

        candidate = selected[key]

        print()
        print(
            "  {} / {}".format(
                key[0],
                key[1],
            )
        )

        print(
            "    window={}".format(
                candidate["window"]
            )
        )

        print(
            "    Z={}".format(
                candidate["forced_Z"]
            )
        )

        print(
            "    integer={}".format(
                rational_is_integer(
                    candidate["forced_Z"]
                )
            )
        )

    return selected


# ============================================================================
# 3. EXACT CONSISTENCY MATRIX
# ============================================================================

def consistency_matrix(selected):

    print()
    print("=" * 78)
    print(
        "3. EXACT CROSS-DISCOVERY CONSISTENCY"
    )
    print("=" * 78)

    keys = [
        key
        for key in DISCOVERY_KEYS
        if key in selected
    ]

    for key_i in keys:

        row = []

        for key_j in keys:

            zi = sp.Rational(
                selected[key_i][
                    "forced_Z"
                ]
            )

            zj = sp.Rational(
                selected[key_j][
                    "forced_Z"
                ]
            )

            row.append(
                zi == zj
            )

        print(
            "  {}: {}".format(
                " / ".join(key_i),
                row,
            )
        )

    equal_pairs = []

    for i in range(
        len(keys)
    ):

        for j in range(
            i + 1,
            len(keys),
        ):

            zi = sp.Rational(
                selected[keys[i]][
                    "forced_Z"
                ]
            )

            zj = sp.Rational(
                selected[keys[j]][
                    "forced_Z"
                ]
            )

            if zi == zj:

                equal_pairs.append(
                    (
                        keys[i],
                        keys[j],
                        zi,
                    )
                )

    print()
    print(
        "  equal_pairs={}".format(
            equal_pairs
        )
    )

    print(
        "  distinct_Z_count={}".format(
            len({
                sp.Rational(
                    selected[key]["forced_Z"]
                )
                for key in keys
            })
        )
    )

    return keys, equal_pairs


# ============================================================================
# 4. SYMBOLIC RESIDUAL VERIFICATION
# ============================================================================

def residual_verification(
    results,
    selected,
):

    print()
    print("=" * 78)
    print(
        "4. EXACT SYMBOLIC RESIDUAL VERIFICATION"
    )
    print("=" * 78)

    failures = 0

    for key, candidate in selected.items():

        result = results[key]

        # Rebuild the transformed symbolic lattice.
        r_order, t_order = TRANSFORMS[
            key[0]
        ]

        offsets = STENCILS[
            key[1]
        ]

        lattice = transformed_lattice(
            build_symbolic_base_lattice(),
            r_order,
            t_order,
        )

        window = candidate[
            "window"
        ]

        relation = candidate[
            "relation"
        ]

        values = [
            clean(
                lattice[
                    (
                        window[0] + dr,
                        window[1] + dt,
                    )
                ]
            )
            for dr, dt in offsets
        ]

        residual = clean(
            sum(
                sp.Integer(
                    relation[i]
                ) * values[i]
                for i in range(
                    len(values)
                )
            )
        )

        evaluated = clean(
            residual.subs(
                Z,
                candidate["forced_Z"],
            )
        )

        print()
        print(
            "  {} / {}".format(
                key[0],
                key[1],
            )
        )

        print(
            "    symbolic_residual={}".format(
                residual
            )
        )

        print(
            "    residual_at_forced_Z={}".format(
                evaluated
            )
        )

        if evaluated != 0:
            failures += 1

    print()
    print(
        "  residual_failures={}".format(
            failures
        )
    )

    return failures


# ============================================================================
# 5. COMPARE WITH 363R REPORTED VALUES
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


def compare_reported_values(
    selected
):

    print()
    print("=" * 78)
    print(
        "5. 363R REPORTED-VALUE COMPARISON"
    )
    print("=" * 78)

    mismatches = []

    for key, candidate in selected.items():

        reconstructed = sp.Rational(
            candidate["forced_Z"]
        )

        reported = REPORTED_363R[
            key
        ]

        matches = (
            reconstructed
            == reported
        )

        print()
        print(
            "  {} / {}".format(
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
            mismatches.append(
                key
            )

    print()
    print(
        "  reported_value_mismatches={}".format(
            len(mismatches)
        )
    )

    return mismatches


# ============================================================================
# 6. INTEGER COMPATIBILITY
# ============================================================================

def integer_compatibility(
    selected
):

    print()
    print("=" * 78)
    print(
        "6. EXACT INTEGER-Z COMPATIBILITY"
    )
    print("=" * 78)

    integer_keys = []

    for key, candidate in selected.items():

        value = sp.Rational(
            candidate["forced_Z"]
        )

        if rational_is_integer(
            value
        ):
            integer_keys.append(
                key
            )

        print()
        print(
            "  {} / {}: integer={}".format(
                key[0],
                key[1],
                rational_is_integer(
                    value
                ),
            )
        )

    print()
    print(
        "  integer_compatible_count={}".format(
            len(integer_keys)
        )
    )

    return integer_keys


# ============================================================================
# 7. FINAL VERDICT
# ============================================================================

def final_verdict(
    selected,
    residual_failures,
    reported_mismatches,
):
    distinct = len({
        sp.Rational(
            candidate["forced_Z"]
        )
        for candidate in selected.values()
    })

    count = len(selected)

    if residual_failures:
        verdict = (
            "RECONSTRUCTION_BUG_OR_INVALID_DISCOVERY"
        )

    elif reported_mismatches:
        verdict = (
            "363R_REPORTED_VALUES_DO_NOT_MATCH_RECONSTRUCTION"
        )

    elif count == 0:
        verdict = (
            "NO_RECONSTRUCTED_DISCOVERIES"
        )

    elif distinct == 1:
        verdict = (
            "COMMON_RECONSTRUCTED_Z"
        )

    elif distinct == count:
        verdict = (
            "CROSS_TRANSFORM_INCONSISTENCY_CONFIRMED"
        )

    else:
        verdict = (
            "PARTIAL_CROSS_TRANSFORM_COMPATIBILITY"
        )

    print()
    print("=" * 78)
    print(
        "7. FINAL STRUCTURAL VERDICT"
    )
    print("=" * 78)

    print(
        "  reconstructed_discoveries={}".format(
            count
        )
    )

    print(
        "  distinct_reconstructed_Z={}".format(
            distinct
        )
    )

    print(
        "  residual_failures={}".format(
            residual_failures
        )
    )

    print(
        "  reported_value_mismatches={}".format(
            len(reported_mismatches)
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
        "EXPERIMENT 365R — EXACT RECONSTRUCTION OF TRANSFORMED "
        "INTEGER ANNIHILATORS / Z-FORCING PROVENANCE AUDIT"
    )
    print("=" * 78)

    results = reconstruction_audit()

    selected = select_canonical_candidates(
        results
    )

    consistency_matrix(
        selected
    )

    residual_failures = (
        residual_verification(
            results,
            selected,
        )
    )

    reported_mismatches = (
        compare_reported_values(
            selected
        )
    )

    integer_keys = (
        integer_compatibility(
            selected
        )
    )

    verdict = final_verdict(
        selected,
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
        "  source_values_observed_only=True"
    )

    print(
        "  symbolic_Z_preserved=True"
    )

    print(
        "  transformed_lattices_rebuilt_directly=True"
    )

    print(
        "  Z_position_identified_from_actual_window=True"
    )

    print(
        "  stored_relation_order_assumption=False"
    )

    print(
        "  symbolic_residuals_rebuilt=True"
    )

    print(
        "  symbolic_residual_failures={}".format(
            residual_failures
        )
    )

    print(
        "  reported_363R_mismatches={}".format(
            len(reported_mismatches)
        )
    )

    print(
        "  integer_compatible_reconstructions={}".format(
            len(integer_keys)
        )
    )

    print(
        "  final_verdict={}".format(
            verdict
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
        "EXPERIMENT 365R COMPLETE"
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
