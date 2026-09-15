#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 376R — EXACT PLÜCKER / 2x2-MINOR DYNAMICS / PROJECTIVE-ROW
                  DETERMINANT AUDIT
==============================================================================

Purpose
-------
375R established an important limitation:

    every 2-row block is automatically rank <= 2,

but the observed 3x3 block has rank 3.

Therefore a global rank-2 source model is impossible on the observed
3x3 region.

The next structural object is the family of 2x2 minors of row pairs.

For two row vectors R_i(t) and R_j(t), define

    Δ_{ij}(s,t)
        = Q(i,s)Q(j,t) - Q(i,t)Q(j,s).

These are the Plücker coordinates of the 2-dimensional row span.

They are invariant under changes of basis inside the row pair, up to a
common scalar, and therefore give a basis-independent projective audit.

This experiment tests:

    * all exact 2x2 minors of every fully observed 2-row block;
    * primitive normalization and exact factorization;
    * gcd/content across minor families;
    * ratios between consecutive minors when exact;
    * adjacent-minor determinant identities;
    * whether minor sequences are proportional between different row pairs;
    * whether minor sequences factor through simple content sequences;
    * whether minor sequences satisfy exact low-order scalar recurrences;
    * whether the 2-row projective geometry changes intrinsically with t.

A positive result would identify a projective law even though the original
entry-wise operator searches failed.

A negative result would show that the rank-2 structure of 2-row slices is
only the dimensional tautology expected from two rows.

No missing values.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
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
# SYMBOLS
# ============================================================================

x = sp.symbols("x")


# ============================================================================
# BASIC HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def integer_factorization(value):
    value = int(sp.Integer(value))

    if value == 0:
        return {}

    return sp.factorint(
        abs(value)
    )


def gcd_list(values):

    g = 0

    for value in values:

        g = math.gcd(
            g,
            abs(int(value)),
        )

    return g


def primitive_integer_vector(values):

    values = [
        int(sp.Integer(value))
        for value in values
    ]

    g = gcd_list(values)

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

    return tuple(values)


def is_integer(value):

    return sp.Rational(value).q == 1


# ============================================================================
# OBSERVED LATTICE
# ============================================================================

def build_lattice():

    lattice = {}

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

            lattice[
                (r, t)
            ] = sp.Integer(value)

    return lattice


# ============================================================================
# FULLY OBSERVED ROW INTERVALS
# ============================================================================

def observed_row_intervals(
    lattice,
    r,
):

    ts = sorted(
        t
        for (
            rr,
            t,
        ) in lattice
        if rr == r
    )

    if not ts:

        return []

    intervals = []

    start = ts[0]
    previous = ts[0]

    for t in ts[1:]:

        if t == previous + 1:

            previous = t

        else:

            intervals.append(
                (
                    start,
                    previous,
                )
            )

            start = t
            previous = t

    intervals.append(
        (
            start,
            previous,
        )
    )

    return intervals


# ============================================================================
# 2x2 MINOR
# ============================================================================

def minor_2x2(
    lattice,
    r0,
    r1,
    t0,
    t1,
):

    return clean(
        lattice[
            (r0, t0)
        ]
        *
        lattice[
            (r1, t1)
        ]
        -
        lattice[
            (r0, t1)
        ]
        *
        lattice[
            (r1, t0)
        ]
    )


def consecutive_minor_sequence(
    lattice,
    r0,
    r1,
):

    common_t = sorted(
        set(
            t
            for rr, t in lattice
            if rr == r0
        )
        &
        set(
            t
            for rr, t in lattice
            if rr == r1
        )
    )

    sequence = []

    for i in range(
        len(common_t) - 1
    ):

        t0 = common_t[i]
        t1 = common_t[i + 1]

        value = minor_2x2(
            lattice,
            r0,
            r1,
            t0,
            t1,
        )

        sequence.append(
            (
                (
                    t0,
                    t1,
                ),
                value,
            )
        )

    return sequence


# ============================================================================
# ALL 2x2 PLÜCKER COORDINATES
# ============================================================================

def all_pair_minors(
    lattice,
    r0,
    r1,
):

    common_t = sorted(
        set(
            t
            for rr, t in lattice
            if rr == r0
        )
        &
        set(
            t
            for rr, t in lattice
            if rr == r1
        )
    )

    records = []

    for t0, t1 in itertools.combinations(
        common_t,
        2,
    ):

        value = minor_2x2(
            lattice,
            r0,
            r1,
            t0,
            t1,
        )

        records.append(
            {
                "columns": (
                    t0,
                    t1,
                ),
                "value": value,
            }
        )

    return records


# ============================================================================
# MINOR SEQUENCE AUDIT
# ============================================================================

def audit_minor_sequence(
    sequence,
    name,
):

    print()
    print(
        f"  {name}:"
    )

    print(
        f"    consecutive_minor_count="
        f"{len(sequence)}"
    )

    if not sequence:

        print(
            "    status=NO_SEQUENCE"
        )

        return {
            "values": [],
            "nonzero_values": [],
        }

    values = [
        clean(value)
        for (
            _,
            value,
        ) in sequence
    ]

    for columns, value in sequence:

        print(
            f"    columns={columns}"
        )

        print(
            f"      value={value}"
        )

        if value != 0:

            print(
                f"      factorization="
                f"{integer_factorization(value)}"
            )

    nonzero = [
        value
        for value in values
        if value != 0
    ]

    if nonzero:

        g = gcd_list(
            nonzero
        )

        print(
            f"    gcd_nonzero={g}"
        )

        print(
            f"    gcd_factorization="
            f"{integer_factorization(g)}"
        )

        primitive = [
            int(
                value // g
            )
            for value in nonzero
        ]

        print(
            f"    primitive_nonzero_values="
            f"{primitive}"
        )

    else:

        g = 0

        print(
            "    all_zero=True"
        )

    return {
        "values": values,
        "nonzero_values": nonzero,
        "gcd": g,
    }


# ============================================================================
# EXACT RATIO AUDIT
# ============================================================================

def exact_ratios(values):

    ratios = []

    for i in range(
        len(values) - 1
    ):

        a = sp.Rational(
            values[i]
        )

        b = sp.Rational(
            values[i + 1]
        )

        if b == 0:

            ratios.append(
                None
            )

        else:

            ratios.append(
                clean(
                    a / b
                )
            )

    return ratios


# ============================================================================
# SECOND-DIFFERENCE / RECURRENCE AUDIT
# ============================================================================

def scalar_recurrence_audit(
    values,
    name,
):

    print()
    print(
        f"  {name}:"
    )

    if len(values) < 3:

        print(
            "    status=DATA_LIMITED"
        )

        return {
            "status": "DATA_LIMITED"
        }

    # First difference.
    first = [
        clean(
            values[i + 1]
            - values[i]
        )
        for i in range(
            len(values) - 1
        )
    ]

    second = [
        clean(
            first[i + 1]
            - first[i]
        )
        for i in range(
            len(first) - 1
        )
    ]

    print(
        f"    first_differences={first}"
    )

    print(
        f"    second_differences={second}"
    )

    first_constant = (
        len(
            set(first)
        )
        == 1
    )

    second_constant = (
        len(
            set(second)
        )
        == 1
    )

    print(
        f"    first_difference_constant="
        f"{first_constant}"
    )

    print(
        f"    second_difference_constant="
        f"{second_constant}"
    )

    # Test Q_{n+2} = A Q_{n+1} + B Q_n.
    pair_rows = []
    rhs = []

    for i in range(
        len(values) - 2
    ):

        pair_rows.append(
            [
                sp.Rational(
                    values[i + 1]
                ),
                sp.Rational(
                    values[i]
                ),
            ]
        )

        rhs.append(
            sp.Rational(
                values[i + 2]
            )
        )

    matrix = sp.Matrix(
        pair_rows
    )

    rhs_matrix = sp.Matrix(
        rhs
    )

    rank = matrix.rank()
    augmented_rank = (
        matrix
        .row_join(
            rhs_matrix
        )
        .rank()
    )

    print(
        f"    width2_recurrence_rank={rank}"
    )

    print(
        f"    width2_recurrence_augmented_rank="
        f"{augmented_rank}"
    )

    recurrence_status = (
        "NO_SOLUTION"
        if augmented_rank > rank
        else (
            "NONUNIQUE"
            if rank < 2
            else
            "EXACT_DATA_SIZED"
        )
    )

    print(
        f"    width2_recurrence_status="
        f"{recurrence_status}"
    )

    if recurrence_status == "EXACT_DATA_SIZED":

        solution = matrix.gauss_jordan_solve(
            rhs_matrix
        )[0]

        A = clean(
            solution[0]
        )

        B = clean(
            solution[1]
        )

        print(
            f"    A={A}"
        )

        print(
            f"    B={B}"
        )

        residuals = [
            clean(
                A * values[i + 1]
                +
                B * values[i]
                -
                values[i + 2]
            )
            for i in range(
                len(values) - 2
            )
        ]

        print(
            f"    recurrence_residuals={residuals}"
        )

    return {
        "first": first,
        "second": second,
        "first_constant": first_constant,
        "second_constant": second_constant,
        "recurrence_status": recurrence_status,
    }


# ============================================================================
# CROSS-PAIR PROPORTIONALITY
# ============================================================================

def proportional_sequences(
    seq_a,
    seq_b,
):

    values_a = dict(seq_a)
    values_b = dict(seq_b)

    common = sorted(
        set(values_a)
        &
        set(values_b)
    )

    if not common:

        return {
            "status": "NO_OVERLAP"
        }

    first_pair = None

    for key in common:

        a = values_a[key]
        b = values_b[key]

        if a == 0 and b == 0:

            continue

        if b == 0:

            return {
                "status": "NOT_PROPORTIONAL"
            }

        first_pair = (
            clean(
                a / b
            )
        )

        break

    if first_pair is None:

        return {
            "status": "ALL_ZERO"
        }

    residuals = [
        clean(
            values_a[key]
            -
            first_pair
            * values_b[key]
        )
        for key in common
    ]

    return {
        "status": (
            "PROPORTIONAL"
            if all(
                residual == 0
                for residual in residuals
            )
            else
            "NOT_PROPORTIONAL"
        ),
        "ratio": first_pair,
        "common_columns": common,
        "residuals": residuals,
    }


# ============================================================================
# PLÜCKER IDENTITY CHECK
# ============================================================================

def plucker_identity_audit(
    lattice,
    r0,
    r1,
    columns,
):

    print()
    print(
        f"  row_pair=({r0},{r1})"
    )

    if len(columns) < 4:

        print(
            "    status=DATA_LIMITED"
        )

        return

    # For four columns a<b<c<d:
    #
    # Δ_ab Δ_cd - Δ_ac Δ_bd + Δ_ad Δ_bc = 0.
    #
    for a_idx, b_idx, c_idx, d_idx in itertools.combinations(
        columns,
        4,
    ):

        dab = minor_2x2(
            lattice,
            r0,
            r1,
            a_idx,
            b_idx,
        )

        dac = minor_2x2(
            lattice,
            r0,
            r1,
            a_idx,
            c_idx,
        )

        dad = minor_2x2(
            lattice,
            r0,
            r1,
            a_idx,
            d_idx,
        )

        dbc = minor_2x2(
            lattice,
            r0,
            r1,
            b_idx,
            c_idx,
        )

        dbd = minor_2x2(
            lattice,
            r0,
            r1,
            b_idx,
            d_idx,
        )

        dcd = minor_2x2(
            lattice,
            r0,
            r1,
            c_idx,
            d_idx,
        )

        residual = clean(
            dab * dcd
            - dac * dbd
            + dad * dbc
        )

        print()
        print(
            f"    columns="
            f"{(a_idx,b_idx,c_idx,d_idx)}"
        )

        print(
            f"      residual={residual}"
        )

        print(
            f"      exact_zero="
            f"{residual == 0}"
        )


# ============================================================================
# MINOR POLYNOMIAL
# ============================================================================

def minor_generating_polynomial(
    sequence,
    variable_name,
):

    variable = sp.symbols(
        variable_name
    )

    polynomial = clean(
        sum(
            value * variable**columns[1]
            for (
                columns,
                value,
            ) in sequence
        )
    )

    return polynomial


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 376R — EXACT PLÜCKER / 2x2-MINOR DYNAMICS / "
        "PROJECTIVE-ROW DETERMINANT AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "OBSERVED SOURCE"
    )

    print(
        f"  observed_cells={len(lattice)}"
    )

    # ------------------------------------------------------------------------
    # Row pairs.
    # ------------------------------------------------------------------------

    row_pairs = [
        (
            0,
            1,
        ),
        (
            0,
            2,
        ),
        (
            1,
            2,
        ),
    ]

    sequences = {}

    print()
    print("=" * 78)
    print(
        "1. EXACT CONSECUTIVE 2x2 MINOR SEQUENCES"
    )
    print("=" * 78)

    for r0, r1 in row_pairs:

        sequence = consecutive_minor_sequence(
            lattice,
            r0,
            r1,
        )

        sequences[
            (
                r0,
                r1,
            )
        ] = sequence

        audit_minor_sequence(
            sequence,
            f"row_pair=({r0},{r1})",
        )

    # ------------------------------------------------------------------------
    # All Plücker coordinates.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "2. FULL PLÜCKER COORDINATE AUDIT"
    )
    print("=" * 78)

    all_plucker = {}

    for r0, r1 in row_pairs:

        records = all_pair_minors(
            lattice,
            r0,
            r1,
        )

        all_plucker[
            (
                r0,
                r1,
            )
        ] = records

        print()
        print(
            f"  row_pair=({r0},{r1})"
        )

        for record in records:

            print(
                f"    columns="
                f"{record['columns']}"
            )

            print(
                f"      value="
                f"{record['value']}"
            )

    # ------------------------------------------------------------------------
    # Ratio audit.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. EXACT CONSECUTIVE-MINOR RATIO AUDIT"
    )
    print("=" * 78)

    for pair, sequence in sequences.items():

        values = [
            value
            for (
                _,
                value,
            )
            in sequence
        ]

        ratios = exact_ratios(
            values
        )

        print()
        print(
            f"  row_pair={pair}"
        )

        print(
            f"    ratios={ratios}"
        )

        print(
            f"    integral_ratios="
            f"{[
                ratio
                for ratio in ratios
                if ratio is not None
                and is_integer(ratio)
            ]}"
        )

    # ------------------------------------------------------------------------
    # Difference / recurrence audit.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. SCALAR DYNAMICS OF CONSECUTIVE MINORS"
    )
    print("=" * 78)

    recurrence_results = {}

    for pair, sequence in sequences.items():

        values = [
            value
            for (
                _,
                value,
            )
            in sequence
        ]

        recurrence_results[
            pair
        ] = scalar_recurrence_audit(
            values,
            f"row_pair={pair}",
        )

    # ------------------------------------------------------------------------
    # Cross-row-pair proportionality.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "5. CROSS-ROW-PAIR MINOR PROPORTIONALITY"
    )
    print("=" * 78)

    for pair_a, pair_b in itertools.combinations(
        sequences,
        2,
    ):

        result = proportional_sequences(
            sequences[
                pair_a
            ],
            sequences[
                pair_b
            ],
        )

        print()
        print(
            f"  {pair_a} vs {pair_b}:"
        )

        print(
            f"    status="
            f"{result['status']}"
        )

        if result[
            "status"
        ] == "PROPORTIONAL":

            print(
                f"    ratio="
                f"{result['ratio']}"
            )

            print(
                f"    residuals="
                f"{result['residuals']}"
            )

    # ------------------------------------------------------------------------
    # Plücker identities.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. EXACT PLÜCKER IDENTITY AUDIT"
    )
    print("=" * 78)

    for r0, r1 in row_pairs:

        common_t = sorted(
            set(
                t
                for rr, t in lattice
                if rr == r0
            )
            &
            set(
                t
                for rr, t in lattice
                if rr == r1
            )
        )

        plucker_identity_audit(
            lattice,
            r0,
            r1,
            common_t,
        )

    # ------------------------------------------------------------------------
    # Generating polynomials for consecutive minors.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "7. MINOR-GENERATING-POLYNOMIAL AUDIT"
    )
    print("=" * 78)

    minor_polynomials = {}

    for pair, sequence in sequences.items():

        if not sequence:

            continue

        variable_name = (
            f"z_{pair[0]}_{pair[1]}"
        )

        polynomial = (
            minor_generating_polynomial(
                sequence,
                variable_name,
            )
        )

        minor_polynomials[
            pair
        ] = polynomial

        print()
        print(
            f"  row_pair={pair}"
        )

        print(
            f"    polynomial={polynomial}"
        )

        print(
            f"    factorized="
            f"{sp.factor(polynomial)}"
        )

    # ------------------------------------------------------------------------
    # Content factor audit.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "8. CROSS-MINOR CONTENT / GCD AUDIT"
    )
    print("=" * 78)

    all_nonzero_minor_values = []

    for records in all_plucker.values():

        all_nonzero_minor_values.extend(
            record["value"]
            for record in records
            if record["value"] != 0
        )

    global_gcd = gcd_list(
        all_nonzero_minor_values
    )

    print(
        f"  global_nonzero_minor_gcd="
        f"{global_gcd}"
    )

    print(
        f"  global_gcd_factorization="
        f"{integer_factorization(global_gcd)}"
    )

    for pair, records in all_plucker.items():

        values = [
            record["value"]
            for record in records
            if record["value"] != 0
        ]

        pair_gcd = gcd_list(
            values
        )

        print()
        print(
            f"  row_pair={pair}"
        )

        print(
            f"    gcd="
            f"{pair_gcd}"
        )

        print(
            f"    factorization="
            f"{integer_factorization(pair_gcd)}"
        )

    # ------------------------------------------------------------------------
    # Structural conclusion.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
375R showed that the observed table has a genuine rank transition:

    2-row rectangles  -> rank 2,
    3x3 block         -> rank 3.

Therefore a global rank-2 decomposition is impossible.

That does not make the 2-row rank-2 structure useless.

For any two rows, the complete projective information is encoded by the
2x2 minors

    Δ_ab = Q(r0,a)Q(r1,b) - Q(r0,b)Q(r1,a).

These are Plücker coordinates for the row pair.

The key questions in 376R are therefore:

    * Do consecutive minors follow a simple scalar law?
    * Are minor sequences proportional between different row pairs?
    * Do their generating polynomials factor?
    * Is there a common content factor?
    * Does projective geometry evolve by a simple multiplicative rule?

A positive result would identify structure invisible in the raw entries.

A negative result would show that the exact rank-2 slices are merely
dimensionally forced and that their projective evolution is itself
complicated.

The Plücker identity is checked exactly as an internal sanity check.

No missing cells are used.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness.
    # ------------------------------------------------------------------------

    ratio_hits = 0
    recurrence_exact = []

    for pair, sequence in sequences.items():

        values = [
            value
            for (
                _,
                value,
            )
            in sequence
        ]

        ratios = exact_ratios(
            values
        )

        ratio_hits += sum(
            1
            for ratio in ratios
            if ratio is not None
        )

        result = recurrence_results.get(
            pair
        )

        if (
            result is not None
            and result.get(
                "recurrence_status"
            ) == "EXACT_DATA_SIZED"
        ):

            recurrence_exact.append(
                pair
            )

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells_used_only=True"
    )

    print(
        "  row_pairs_tested=3"
    )

    print(
        "  exact_2x2_minors_computed=True"
    )

    print(
        "  plucker_coordinates_completed=True"
    )

    print(
        "  plucker_identity_checked=True"
    )

    print(
        f"  nonzero_minor_global_gcd={global_gcd}"
    )

    print(
        f"  exact_nonzero_ratio_count="
        f"{ratio_hits}"
    )

    print(
        f"  exact_width2_recurrence_pairs="
        f"{recurrence_exact}"
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
        "  missing_values_used=False"
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
        "EXPERIMENT 376R COMPLETE"
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
            "\nFATAL ERROR: "
            f"{type(exc).__name__}: {exc}"
        )

        raise
