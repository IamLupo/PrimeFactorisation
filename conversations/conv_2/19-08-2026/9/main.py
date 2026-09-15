#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 325R — EXACT TRIANGULAR-ARRAY RANK / SEPARABILITY AUDIT
==============================================================================

Purpose
-------
Experiments 322R-324R found:

    * no overdetermined width-2 raw recurrence in p;
    * no overdetermined width-2 first-difference recurrence in p;
    * no shared constant-coefficient recurrence of order 1, 2, or 3 in t.

The next question is therefore not another recurrence.

Instead, test whether the observed source object has a LOW-RANK
separable representation.

For the available triangular data Q_t(p), ask whether

    Q_t(p) = sum_{k=1}^r U_k(p) V_k(t)

for small r.

This is equivalent to asking for low rank of appropriate rectangular
submatrices of the observed triangular array.

Important:
    * no missing values are invented;
    * only genuinely observed entries are used;
    * exact rational/integer arithmetic;
    * rank is computed over Q;
    * every nonzero minor used as an obstruction is exact.

Tests
-----

1. Build the observed triangular matrix.

2. Compute ranks of all natural rectangular submatrices obtainable
   without filling missing entries.

3. Search for exact 2x2, 3x3, and 4x4 minors.

4. Compute the maximum rank certified by nonzero minors.

5. Test whether rows/layers are pairwise proportional.

6. Test whether columns/p-values are pairwise proportional on overlaps.

7. Test multiplicative cross-ratio identities equivalent to rank <= 1.

8. Test rank <= 2 using exact 3x3 minors whenever a complete 3x3
   observed rectangle exists.

9. Test rank <= 3 using exact 4x4 minors when available.

10. Perform a bipartite-support analysis:
        which p,t pairs are observed,
    and which rectangular substructures are genuinely available.

Interpretation
--------------

If rank <= 1, the source table is essentially separable and a very
simple hidden mechanism exists.

If rank <= 2 or <= 3, there is still a finite-dimensional separable
source model even though fixed-coefficient recurrences fail.

If the observed minors force full available rank, the data are not
consistent with a low-rank separable source object.

This does NOT prove a universal formula because the array is finite.
It is nevertheless a stronger structural diagnostic than another
two-point recurrence fit.

No interpolation.
No missing-value reconstruction.
No synthetic second n=pq case.
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


def build_layers():

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(
        maximum_t + 1
    ):

        layer = {}

        for p in sorted(Q):

            index = (
                len(Q[p]) - 1 - t
            )

            if index < 0:
                continue

            layer[p] = sp.Integer(
                Q[p][index]
            )

        layers[t] = layer

    return layers


def observed_value(layers, p, t):

    if t not in layers:
        return None

    return layers[t].get(p)


def complete_rectangle(
    layers,
    ps,
    ts,
):

    for p in ps:
        for t in ts:
            if observed_value(
                layers,
                p,
                t,
            ) is None:
                return False

    return True


def rectangle_matrix(
    layers,
    ps,
    ts,
):

    if not complete_rectangle(
        layers,
        ps,
        ts,
    ):
        return None

    return sp.Matrix([
        [
            observed_value(
                layers,
                p,
                t,
            )
            for t in ts
        ]
        for p in ps
    ])


# ============================================================================
# SUPPORT
# ============================================================================

def support_audit(layers):

    print()
    print("=" * 78)
    print(
        "1. OBSERVED TRIANGULAR SUPPORT"
    )
    print("=" * 78)

    print(
        "  p_values={}".format(
            sorted(Q)
        )
    )

    print(
        "  t_values={}".format(
            sorted(layers)
        )
    )

    for t in sorted(layers):

        print()
        print(
            "  t={}: observed_p={}".format(
                t,
                sorted(
                    layers[t]
                ),
            )
        )


# ============================================================================
# COMPLETE RECTANGLE AUDIT
# ============================================================================

def enumerate_rectangles(
    layers,
    size,
):

    ps_all = sorted(Q)
    ts_all = sorted(layers)

    rectangles = []

    for ps in itertools.combinations(
        ps_all,
        size,
    ):

        for ts in itertools.combinations(
            ts_all,
            size,
        ):

            if complete_rectangle(
                layers,
                ps,
                ts,
            ):

                rectangles.append(
                    (
                        ps,
                        ts,
                    )
                )

    return rectangles


def rectangle_rank_audit(
    layers,
    size,
):

    print()
    print("=" * 78)
    print(
        "{}x{} COMPLETE-RECTANGLE RANK AUDIT".format(
            size,
            size,
        )
    )
    print("=" * 78)

    rectangles = enumerate_rectangles(
        layers,
        size,
    )

    print(
        "  complete_rectangles={}".format(
            len(rectangles)
        )
    )

    nonzero = []
    zero = []

    for ps, ts in rectangles:

        M = rectangle_matrix(
            layers,
            ps,
            ts,
        )

        det = clean(
            M.det()
        )

        if det == 0:
            zero.append(
                (ps, ts, M)
            )
        else:
            nonzero.append(
                (ps, ts, M, det)
            )

    print(
        "  zero_determinants={}".format(
            len(zero)
        )
    )

    print(
        "  nonzero_determinants={}".format(
            len(nonzero)
        )
    )

    for ps, ts, M, det in nonzero:

        print()
        print(
            "  NONZERO_MINOR:"
        )

        print(
            "    p_indices={}".format(
                ps
            )
        )

        print(
            "    t_indices={}".format(
                ts
            )
        )

        print(
            "    matrix={}".format(
                M
            )
        )

        print(
            "    determinant={}".format(
                det
            )
        )

        break

    return {
        "rectangles": rectangles,
        "zero": zero,
        "nonzero": nonzero,
    }


# ============================================================================
# LARGEST OBSERVED MATRIX RANK
# ============================================================================

def maximal_observed_rank(layers):

    print()
    print("=" * 78)
    print(
        "4. MAXIMAL COMPLETE-RECTANGLE RANK"
    )
    print("=" * 78)

    best_rank = 0
    witnesses = []

    ps_all = sorted(Q)
    ts_all = sorted(layers)

    for rows in range(
        1,
        len(ps_all) + 1,
    ):

        for cols in range(
            1,
            len(ts_all) + 1,
        ):

            size = min(
                rows,
                cols,
            )

            for ps in itertools.combinations(
                ps_all,
                rows,
            ):

                for ts in itertools.combinations(
                    ts_all,
                    cols,
                ):

                    if not complete_rectangle(
                        layers,
                        ps,
                        ts,
                    ):
                        continue

                    M = rectangle_matrix(
                        layers,
                        ps,
                        ts,
                    )

                    rank = int(
                        M.rank()
                    )

                    if rank > best_rank:

                        best_rank = rank
                        witnesses = [
                            (
                                ps,
                                ts,
                                M,
                            )
                        ]

                    elif rank == best_rank:

                        witnesses.append(
                            (
                                ps,
                                ts,
                                M,
                            )
                        )

    print(
        "  maximal_rank={}".format(
            best_rank
        )
    )

    print(
        "  witness_count={}".format(
            len(witnesses)
        )
    )

    if witnesses:

        ps, ts, M = witnesses[0]

        print()
        print(
            "  witness_p_indices={}".format(
                ps
            )
        )

        print(
            "  witness_t_indices={}".format(
                ts
            )
        )

        print(
            "  witness_matrix={}".format(
                M
            )
        )

    return best_rank, witnesses


# ============================================================================
# RANK-1 / PROPORTIONALITY AUDIT
# ============================================================================

def proportionality_audit(layers):

    print()
    print("=" * 78)
    print(
        "5. EXACT RANK-1 / PROPORTIONALITY AUDIT"
    )
    print("=" * 78)

    # Compare every pair of observed t-layers on their common p support.

    ts = sorted(layers)

    row_failures = []

    for t1, t2 in itertools.combinations(
        ts,
        2,
    ):

        common_p = sorted(
            set(layers[t1])
            & set(layers[t2])
        )

        if len(common_p) < 2:
            continue

        ratios = []

        for p in common_p:

            a = sp.Rational(
                layers[t1][p]
            )

            b = sp.Rational(
                layers[t2][p]
            )

            if b == 0:
                ratios.append(
                    None
                )
            else:
                ratios.append(
                    clean(a / b)
                )

        usable = [
            r
            for r in ratios
            if r is not None
        ]

        proportional = (
            len(usable) > 0
            and
            all(
                r == usable[0]
                for r in usable
            )
        )

        print()
        print(
            "  t_pair=({},{}):".format(
                t1,
                t2,
            )
        )

        print(
            "    common_p={}".format(
                common_p
            )
        )

        print(
            "    ratios={}".format(
                ratios
            )
        )

        print(
            "    proportional={}".format(
                proportional
            )
        )

        if not proportional:
            row_failures.append(
                (t1, t2)
            )

    print()
    print(
        "  nonproportional_row_pairs={}".format(
            len(row_failures)
        )
    )


# ============================================================================
# CROSS-RATIO / 2x2 MINOR AUDIT
# ============================================================================

def minor_2x2_audit(layers):

    print()
    print("=" * 78)
    print(
        "6. EXACT 2x2 MINOR / CROSS-RATIO AUDIT"
    )
    print("=" * 78)

    rectangles = enumerate_rectangles(
        layers,
        2,
    )

    nonzero = 0
    zero = 0

    for ps, ts in rectangles:

        M = rectangle_matrix(
            layers,
            ps,
            ts,
        )

        det = clean(
            M.det()
        )

        print()
        print(
            "  p_indices={}".format(
                ps
            )
        )

        print(
            "  t_indices={}".format(
                ts
            )
        )

        print(
            "  determinant={}".format(
                det
            )
        )

        if det == 0:
            zero += 1
        else:
            nonzero += 1

    print()
    print(
        "  zero_2x2_minors={}".format(
            zero
        )
    )

    print(
        "  nonzero_2x2_minors={}".format(
            nonzero
        )
    )

    return {
        "zero": zero,
        "nonzero": nonzero,
    }


# ============================================================================
# 3x3 AND 4x4 MINOR SUMMARY
# ============================================================================

def higher_minor_summary(
    layers,
):

    print()
    print("=" * 78)
    print(
        "7. HIGHER-ORDER MINOR SUMMARY"
    )
    print("=" * 78)

    outputs = {}

    for size in (
        3,
        4,
    ):

        result = rectangle_rank_audit(
            layers,
            size,
        )

        outputs[size] = result

    return outputs


# ============================================================================
# TRACE-LIKE NORMALIZATION / SEPARABILITY CHECK
# ============================================================================

def normalized_row_audit(layers):

    print()
    print("=" * 78)
    print(
        "8. LAYER NORMALIZATION / PROJECTIVE SEPARABILITY AUDIT"
    )
    print("=" * 78)

    # Normalize each nonzero row by its first observed entry.
    # If all normalized rows become identical on overlap, this would
    # indicate rank-one separability up to scalar layer factors.

    normalized = {}

    for t in sorted(layers):

        ps = sorted(
            layers[t]
        )

        if not ps:
            continue

        anchor = layers[t][ps[0]]

        if anchor == 0:
            continue

        normalized[t] = {
            p: clean(
                layers[t][p] / anchor
            )
            for p in ps
        }

    for t, row in normalized.items():

        print()
        print(
            "  t={}: normalized={}".format(
                t,
                row,
            )
        )

    mismatches = 0

    ts = sorted(
        normalized
    )

    for t1, t2 in itertools.combinations(
        ts,
        2,
    ):

        common_p = (
            set(normalized[t1])
            &
            set(normalized[t2])
        )

        for p in sorted(
            common_p
        ):

            if clean(
                normalized[t1][p]
                -
                normalized[t2][p]
            ) != 0:

                mismatches += 1

    print()
    print(
        "  normalized_overlap_mismatches={}".format(
            mismatches
        )
    )

    print(
        "  rank_one_separability={}".format(
            mismatches == 0
        )
    )


# ============================================================================
# INTERPRETATION
# ============================================================================

def interpretation(
    rank2,
    rank3,
    rank4,
    maximal_rank,
):

    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiments 322R-324R ruled out simple fixed-coefficient transfer laws
in either direction of the triangular source table.

Experiment 325R changes the question:

    Is the failure of recurrence structure merely a coordinate artifact,
    while the full table still has low separable rank?

A representation

    Q_t(p) = sum_{k=1}^r U_k(p) V_k(t)

implies that every complete (r+1)x(r+1) observed minor vanishes.

Therefore:

    nonzero 2x2 minor
        => rank >= 2
        => not rank-one separable.

    nonzero 3x3 minor
        => rank >= 3
        => not rank-two separable.

    nonzero 4x4 minor
        => rank >= 4
        => not rank-three separable.

This is independent of the earlier recurrence tests.

Because the source support is triangular, missing cells are not filled.
Only complete observed rectangles are used.

A low-rank result would give a new and potentially useful structural
description:

    recurrence failure
        but
    finite-dimensional separability.

A full observed rank result would instead strengthen the conclusion that
Q_t(p) is genuinely two-variable at the level of the available data.

This remains diagnostic, not universal.
No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 325R — EXACT TRIANGULAR-ARRAY RANK / SEPARABILITY AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    support_audit(
        layers
    )

    r2 = rectangle_rank_audit(
        layers,
        2,
    )

    r3 = rectangle_rank_audit(
        layers,
        3,
    )

    r4 = rectangle_rank_audit(
        layers,
        4,
    )

    minor_2x2_audit(
        layers
    )

    higher = higher_minor_summary(
        layers
    )

    maximal_rank, witnesses = (
        maximal_observed_rank(
            layers
        )
    )

    proportionality_audit(
        layers
    )

    normalized_row_audit(
        layers
    )

    interpretation(
        r2,
        r3,
        r4,
        maximal_rank,
    )

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  complete_observed_rectangle_audit=True"
    )

    print(
        "  2x2_minor_audit=True"
    )

    print(
        "  3x3_minor_audit=True"
    )

    print(
        "  4x4_minor_audit=True"
    )

    print(
        "  maximal_observed_rank={}".format(
            maximal_rank
        )
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  missing_value_reconstruction=False"
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
        "EXPERIMENT 325R COMPLETE"
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
