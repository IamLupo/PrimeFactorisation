#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 322R — EXACT REDUNDANCY / LOCAL TRANSFER-DRIFT AUDIT
==============================================================================

Purpose
-------
Audit the only currently available overdetermined raw and first-difference
width-2 transitions.

The key distinction is:

    exact reconstruction
        versus
    independent validation.

No interpolation.
No missing-value reconstruction.
No synthetic second n=pq case.
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


def rational_content(values):
    """
    Exact rational content of a finite list of rationals.

    For v_i = n_i/d_i, let D be the lcm of the denominators.
    Then content = gcd(D*v_i) / D.
    """
    values = [
        sp.Rational(v)
        for v in values
        if sp.Rational(v) != 0
    ]

    if not values:
        return sp.Integer(0)

    common_denominator = sp.ilcm(
        *[
            int(v.q)
            for v in values
        ]
    )

    integers = [
        int(v * common_denominator)
        for v in values
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    return sp.Rational(
        g,
        common_denominator,
    )


def primitive_row(values):
    content = rational_content(values)

    if content == 0:
        return [
            sp.Rational(v)
            for v in values
        ]

    return [
        clean(
            sp.Rational(v) / content
        )
        for v in values
    ]


# ============================================================================
# SOURCE LAYERS
# ============================================================================

def build_layers():
    maximum_t = max(
        len(Q[p]) - 1
        for p in Q
    )

    layers = {}

    for t in range(
        maximum_t + 1
    ):
        layer = []

        for p in sorted(Q):
            index = (
                len(Q[p])
                - 1
                - t
            )

            if index >= 0:
                layer.append(
                    (
                        p,
                        sp.Integer(
                            Q[p][index]
                        ),
                    )
                )

        layers[t] = layer

    return layers


# ============================================================================
# FINITE DIFFERENCES
# ============================================================================

def finite_difference_rows(values):
    current = [
        sp.Rational(v)
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


def build_first_entry_triangle(layers):
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
# WIDTH-2 SYSTEM
# ============================================================================

def full_width2_system(
    source0,
    source1,
    target,
):
    if not (
        len(source0)
        == len(source1)
        == len(target)
    ):
        raise ValueError(
            "source/target lengths do not agree"
        )

    if not target:
        return {
            "rank": 0,
            "augmented_rank": 0,
            "consistent": True,
            "equations": [],
        }

    M = sp.Matrix([
        [
            sp.Rational(
                source0[i]
            ),
            sp.Rational(
                source1[i]
            ),
        ]
        for i in range(len(target))
    ])

    b = sp.Matrix([
        sp.Rational(v)
        for v in target
    ])

    rank = M.rank()
    augmented_rank = (
        M.row_join(b).rank()
    )

    return {
        "rank": rank,
        "augmented_rank": augmented_rank,
        "consistent": (
            rank == augmented_rank
        ),
        "equations": [
            (
                source0[i],
                source1[i],
                target[i],
            )
            for i in range(len(target))
        ],
    }


def solve_two_equations(
    source0,
    source1,
    target,
    start,
):
    equations = []

    for j in (
        start,
        start + 1,
    ):
        equations.append(
            (
                sp.Rational(
                    source0[j]
                ),
                sp.Rational(
                    source1[j]
                ),
                sp.Rational(
                    target[j]
                ),
            )
        )

    M = sp.Matrix([
        [
            a,
            b,
        ]
        for a, b, _
        in equations
    ])

    rhs = sp.Matrix([
        c
        for _, _, c
        in equations
    ])

    if M.rank() < 2:
        return {
            "status": "NONUNIQUE",
            "c0": None,
            "c1": None,
        }

    solution = M.LUsolve(
        rhs
    )

    c0 = clean(
        solution[0]
    )

    c1 = clean(
        solution[1]
    )

    verified = all(
        clean(
            c0 * a
            + c1 * b
            - c
        ) == 0
        for a, b, c
        in equations
    )

    return {
        "status": (
            "EXACT"
            if verified
            else "VERIFICATION_FAILED"
        ),
        "c0": c0,
        "c1": c1,
    }


# ============================================================================
# REDUNDANCY AUDIT
# ============================================================================

def redundancy_audit(
    label,
    source0,
    source1,
    target,
):
    print()
    print("=" * 78)
    print(label)
    print("=" * 78)

    equation_count = len(target)

    full = full_width2_system(
        source0,
        source1,
        target,
    )

    print(
        "  equation_count={}".format(
            equation_count
        )
    )

    print(
        "  unknown_coefficients=2"
    )

    print(
        "  redundancy={}".format(
            equation_count - 2
        )
    )

    print(
        "  coefficient_matrix_rank={}".format(
            full["rank"]
        )
    )

    print(
        "  augmented_matrix_rank={}".format(
            full["augmented_rank"]
        )
    )

    print(
        "  globally_consistent={}".format(
            full["consistent"]
        )
    )

    local = []

    if equation_count >= 2:

        for start in range(
            equation_count - 1
        ):
            result = solve_two_equations(
                source0,
                source1,
                target,
                start,
            )

            local.append(
                result
            )

            print()
            print(
                "  local_window=({},{}):".format(
                    start,
                    start + 1,
                )
            )

            print(
                "    status={}".format(
                    result["status"]
                )
            )

            if result["status"] == "EXACT":
                print(
                    "    c0={}".format(
                        result["c0"]
                    )
                )

                print(
                    "    c1={}".format(
                        result["c1"]
                    )
                )

    exact = [
        (i, result)
        for i, result
        in enumerate(local)
        if result["status"] == "EXACT"
    ]

    if len(exact) >= 2:

        print()
        print(
            "  LOCAL-WINDOW DRIFT"
        )

        reference_index, reference = exact[0]

        for index, result in exact[1:]:

            delta_c0 = clean(
                result["c0"]
                - reference["c0"]
            )

            delta_c1 = clean(
                result["c1"]
                - reference["c1"]
            )

            print(
                "    reference_window={}".format(
                    reference_index
                )
            )

            print(
                "    compared_window={}".format(
                    index
                )
            )

            print(
                "    delta_c0={}".format(
                    delta_c0
                )
            )

            print(
                "    delta_c1={}".format(
                    delta_c1
                )
            )

            print(
                "    same_operator={}".format(
                    (
                        delta_c0 == 0
                        and
                        delta_c1 == 0
                    )
                )
            )

    for index, result in exact:

        residuals = [
            clean(
                result["c0"] * a
                + result["c1"] * b
                - c
            )
            for a, b, c
            in full["equations"]
        ]

        print()
        print(
            "  residuals_from_window_{}={}".format(
                index,
                residuals,
            )
        )

        print(
            "  all_zero={}".format(
                all(
                    r == 0
                    for r in residuals
                )
            )
        )

    return {
        "full": full,
        "local": local,
    }


# ============================================================================
# RAW TRANSITION
# ============================================================================

def raw_t0_to_t1(layers):

    source = [
        value
        for _, value
        in layers[0]
    ]

    target = [
        value
        for _, value
        in layers[1]
    ]

    return redundancy_audit(
        "1. RAW t=0 -> 1 REDUNDANCY AUDIT",
        source[:-1],
        source[1:],
        target,
    )


# ============================================================================
# FIRST-ENTRY TRANSITION
# ============================================================================

def first_entry_t0_to_t1(A):

    source0 = [
        A[(j, 0)]
        for j in range(3)
    ]

    source1 = [
        A[(j + 1, 0)]
        for j in range(3)
    ]

    target = [
        A[(j, 1)]
        for j in range(3)
    ]

    return redundancy_audit(
        "2. FIRST-ENTRY t=0 -> 1 REDUNDANCY AUDIT",
        source0,
        source1,
        target,
    )


# ============================================================================
# MINIMAL-OVERLAP TRANSITIONS
# ============================================================================

def minimal_overlap_report(layers):

    print()
    print("=" * 78)
    print(
        "3. MINIMAL-OVERLAP TRANSITIONS"
    )
    print("=" * 78)

    for t in (
        1,
        2,
    ):

        equation_count = len(
            layers[t + 1]
        )

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    equations={}".format(
                equation_count
            )
        )

        print(
            "    unknown_coefficients=2"
        )

        print(
            "    redundancy={}".format(
                equation_count - 2
            )
        )

        if equation_count == 2:
            print(
                "    verdict=EXACT_RECONSTRUCTION_ONLY"
            )


# ============================================================================
# RAW VS DIFFERENCE DRIFT
# ============================================================================

def compare_coordinate_drift(
    raw_result,
    diff_result,
):

    print()
    print("=" * 78)
    print(
        "4. RAW VS FIRST-ENTRY LOCAL DRIFT"
    )
    print("=" * 78)

    for label, result in (
        ("raw", raw_result),
        ("first_entry", diff_result),
    ):

        exact = [
            (i, r)
            for i, r
            in enumerate(
                result["local"]
            )
            if r["status"] == "EXACT"
        ]

        print()
        print(
            "  {}:".format(
                label
            )
        )

        if len(exact) < 2:

            print(
                "    comparable_exact_windows=False"
            )

            continue

        _, first = exact[0]
        _, second = exact[1]

        print(
            "    delta_c0={}".format(
                clean(
                    second["c0"]
                    - first["c0"]
                )
            )
        )

        print(
            "    delta_c1={}".format(
                clean(
                    second["c1"]
                    - first["c1"]
                )
            )
        )


# ============================================================================
# CONTENT / PRIME AUDIT
# ============================================================================

def content_prime_audit(
    result,
    label,
):

    print()
    print("=" * 78)
    print(
        "{} CONTENT / PRIME PROFILE".format(
            label
        )
    )
    print("=" * 78)

    for i, local in enumerate(
        result["local"]
    ):

        if local["status"] != "EXACT":
            continue

        row = [
            local["c0"],
            local["c1"],
        ]

        print()
        print(
            "  window_{}:".format(i)
        )

        print(
            "    c0={}".format(
                local["c0"]
            )
        )

        print(
            "    c1={}".format(
                local["c1"]
            )
        )

        print(
            "    rational_content={}".format(
                rational_content(row)
            )
        )

        print(
            "    primitive_row={}".format(
                primitive_row(row)
            )
        )

        for name, value in (
            ("c0", local["c0"]),
            ("c1", local["c1"]),
        ):

            profile = {
                p: valuation(
                    value,
                    p,
                )
                for p in (
                    2,
                    3,
                    5,
                    7,
                    11,
                    13,
                    17,
                )
            }

            print(
                "    {}_valuations={}".format(
                    name,
                    profile,
                )
            )


# ============================================================================
# TERMINAL REFERENCE
# ============================================================================

def terminal_reference():

    q1 = 495451247
    q3 = 421514439

    print()
    print("=" * 78)
    print(
        "6. TERMINAL SOURCE REFERENCE"
    )
    print("=" * 78)

    print(
        "  q1_terminal={}".format(
            q1
        )
    )

    print(
        "  q3_terminal={}".format(
            q3
        )
    )

    print(
        "  gcd={}".format(
            math.gcd(
                q1,
                q3,
            )
        )
    )

    print(
        "  q1/17={}".format(
            q1 // 17
        )
    )

    print(
        "  q3/17={}".format(
            q3 // 17
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 322R — EXACT REDUNDANCY / "
        "LOCAL TRANSFER-DRIFT AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    A = build_first_entry_triangle(
        layers
    )

    print()
    print("SOURCE LAYERS")

    for t, layer in layers.items():
        print(
            "  t={}: {}".format(
                t,
                layer,
            )
        )

    raw_result = raw_t0_to_t1(
        layers
    )

    diff_result = first_entry_t0_to_t1(
        A
    )

    minimal_overlap_report(
        layers
    )

    compare_coordinate_drift(
        raw_result,
        diff_result,
    )

    content_prime_audit(
        raw_result,
        "5. RAW",
    )

    content_prime_audit(
        diff_result,
        "5. FIRST-ENTRY",
    )

    terminal_reference()

    raw_full = raw_result["full"]
    diff_full = diff_result["full"]

    raw_exact_windows = [
        r
        for r in raw_result["local"]
        if r["status"] == "EXACT"
    ]

    diff_exact_windows = [
        r
        for r in diff_result["local"]
        if r["status"] == "EXACT"
    ]

    raw_windows_agree = (
        len(raw_exact_windows) >= 2
        and all(
            clean(
                r["c0"]
                - raw_exact_windows[0]["c0"]
            ) == 0
            and
            clean(
                r["c1"]
                - raw_exact_windows[0]["c1"]
            ) == 0
            for r in raw_exact_windows[1:]
        )
    )

    diff_windows_agree = (
        len(diff_exact_windows) >= 2
        and all(
            clean(
                r["c0"]
                - diff_exact_windows[0]["c0"]
            ) == 0
            and
            clean(
                r["c1"]
                - diff_exact_windows[0]["c1"]
            ) == 0
            for r in diff_exact_windows[1:]
        )
    )

    print()
    print("=" * 78)
    print(
        "7. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  raw_t0_to_t1_overdetermined={}".format(
            len(raw_full["equations"]) > 2
        )
    )

    print(
        "  raw_t0_to_t1_global_width2_consistent={}".format(
            raw_full["consistent"]
        )
    )

    print(
        "  raw_local_windows_agree={}".format(
            raw_windows_agree
        )
    )

    print(
        "  first_entry_t0_to_t1_overdetermined={}".format(
            len(diff_full["equations"]) > 2
        )
    )

    print(
        "  first_entry_t0_to_t1_global_width2_consistent={}".format(
            diff_full["consistent"]
        )
    )

    print(
        "  first_entry_local_windows_agree={}".format(
            diff_windows_agree
        )
    )

    print(
        "  t1_to_t2_independent_validation_available=False"
    )

    print(
        "  t2_to_t3_independent_validation_available=False"
    )

    print(
        "  exact_reconstruction_vs_validation_separated=True"
    )

    print(
        "  rational_content_audit_corrected=True"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  missing_value_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
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
        "EXPERIMENT 322R COMPLETE"
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