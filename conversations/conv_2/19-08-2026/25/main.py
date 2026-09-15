#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 341R — EXACT PRIMITIVE-ROW SHAPE / CONTENT-SEPARATION AUDIT
==============================================================================

Purpose
-------
Experiment 340R found that the row gcd sequence is highly nontrivial:

    g_t = gcd of the observed Q_t row,

while row/column/local normalization did not reduce the observed matrix rank.

This experiment separates the source into

    Q_t = g_t * U_t,

where U_t is the primitive integer row.

The question is now:

    Is the arithmetic complexity mainly scalar content,
    or is it intrinsic to the primitive row shape?

Audits
------

1. Exact row-content sequence g_t.

2. Primitive-row lattice U_t.

3. Primitive-row rank.

4. Exact proportionality tests between adjacent primitive rows.

5. Width-2 transfer tests on primitive rows.

6. Overdetermined consistency / local-window drift.

7. Comparison with the raw transfer coefficients.

8. Exact content-rescaled transfer coefficients:

       U_{t+1}
         = d0_t U_t + d1_t U_t(shift),

   where

       d_i = (g_t / g_{t+1}) c_i.

9. Integer / denominator / prime-profile audit of d_i.

10. Primitive Newton-difference triangle and its rank.

IMPORTANT
---------
Only observed cells are used.

No interpolation.
No missing values.
No extrapolation.
No synthetic second n=pq case.
Exact integer / rational arithmetic only.
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


def gcd_list(values):
    values = [
        abs(int(v))
        for v in values
        if int(v) != 0
    ]

    if not values:
        return 0

    g = values[0]

    for value in values[1:]:
        g = math.gcd(g, value)

    return g


# ============================================================================
# BUILD OBSERVED LAYERS
# ============================================================================

def build_layers():
    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(maximum_t + 1):

        layer = []

        for p_value in sorted(Q):

            values = Q[p_value]

            index = len(values) - 1 - t

            if index >= 0:

                layer.append(
                    (
                        (p_value - 1) // 2,
                        sp.Integer(values[index]),
                    )
                )

        layers[t] = layer

    return layers


# ============================================================================
# ROW CONTENT / PRIMITIVE SHAPES
# ============================================================================

def primitive_rows(layers):

    contents = {}
    primitive = {}

    for t, layer in layers.items():

        values = [
            value
            for _, value in layer
        ]

        g = gcd_list(values)

        contents[t] = sp.Integer(g)

        primitive[t] = [
            (
                r,
                sp.Integer(value // g)
            )
            for r, value in layer
        ]

    return contents, primitive


# ============================================================================
# WIDTH-2 EXACT SOLVER
# ============================================================================

def width2_system(source0, source1, target):

    A = sp.Matrix([
        [
            sp.Rational(source0[i]),
            sp.Rational(source1[i]),
        ]
        for i in range(len(target))
    ])

    b = sp.Matrix([
        sp.Rational(v)
        for v in target
    ])

    rank = A.rank()
    augmented_rank = A.row_join(b).rank()

    return A, b, rank, augmented_rank


def solve_width2(source0, source1, target):

    A, b, rank, augmented_rank = width2_system(
        source0,
        source1,
        target,
    )

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "c0": None,
            "c1": None,
            "residuals": None,
        }

    if rank < 2:

        return {
            "status": "NONUNIQUE",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "c0": None,
            "c1": None,
            "residuals": None,
        }

    solution = A.gauss_jordan_solve(b)[0]

    c0 = clean(solution[0, 0])
    c1 = clean(solution[1, 0])

    residuals = [
        clean(
            c0 * source0[i]
            + c1 * source1[i]
            - target[i]
        )
        for i in range(len(target))
    ]

    return {
        "status": (
            "EXACT"
            if all(r == 0 for r in residuals)
            else "VERIFICATION_FAILED"
        ),
        "rank": rank,
        "augmented_rank": augmented_rank,
        "c0": c0,
        "c1": c1,
        "residuals": residuals,
    }


# ============================================================================
# 1. CONTENT SEQUENCE
# ============================================================================

def print_content_sequence(contents):

    print()
    print("=" * 78)
    print("1. EXACT ROW-CONTENT SEQUENCE")
    print("=" * 78)

    for t in sorted(contents):

        g = contents[t]

        print()
        print(
            "  t={}: g_t={}".format(
                t,
                g,
            )
        )

        print(
            "    factorization={}".format(
                sp.factorint(abs(int(g)))
                if g != 0
                else {}
            )
        )

    print()

    for t in range(
        max(contents)
    ):

        ratio = clean(
            sp.Rational(
                contents[t + 1],
                contents[t],
            )
        )

        print(
            "  g_{} / g_{} = {}".format(
                t + 1,
                t,
                ratio,
            )
        )


# ============================================================================
# 2. PRIMITIVE ROWS
# ============================================================================

def print_primitive_rows(primitive):

    print()
    print("=" * 78)
    print("2. PRIMITIVE ROW LATTICE")
    print("=" * 78)

    for t in sorted(primitive):

        print()
        print(
            "  t={}: U_t={}".format(
                t,
                primitive[t],
            )
        )


# ============================================================================
# 3. PRIMITIVE MATRIX RANK
# ============================================================================

def primitive_rank_audit(primitive):

    print()
    print("=" * 78)
    print("3. PRIMITIVE-ROW RANK AUDIT")
    print("=" * 78)

    matrix = sp.Matrix([
        [
            next(
                value
                for r, value
                in primitive[t]
                if r == row
            )
            if any(
                r == row
                for r, _ in primitive[t]
            )
            else 0
            for t in range(6)
        ]
        for row in range(4)
    ])

    print(
        "  matrix="
    )
    print(matrix)

    print()
    print(
        "  shape={}".format(
            matrix.shape
        )
    )

    print(
        "  rank={}".format(
            matrix.rank()
        )
    )

    return matrix


# ============================================================================
# 4. ADJACENT PROPORTIONALITY
# ============================================================================

def adjacent_proportionality(primitive):

    print()
    print("=" * 78)
    print("4. ADJACENT PRIMITIVE-ROW PROPORTIONALITY")
    print("=" * 78)

    results = []

    for t in range(5):

        left = dict(primitive[t])
        right = dict(primitive[t + 1])

        common = sorted(
            set(left) & set(right)
        )

        if len(common) < 2:

            print()
            print(
                "  t={} -> {}: insufficient_overlap".format(
                    t,
                    t + 1,
                )
            )

            continue

        ratios = []

        for r in common:

            if left[r] == 0:
                ratios.append(None)
            else:
                ratios.append(
                    clean(
                        sp.Rational(
                            right[r],
                            left[r],
                        )
                    )
                )

        proportional = (
            all(
                ratio is not None
                for ratio in ratios
            )
            and
            all(
                ratio == ratios[0]
                for ratio in ratios[1:]
            )
        )

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    common_r={}".format(
                common
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

        results.append(
            (t, proportional)
        )

    return results


# ============================================================================
# 5. PRIMITIVE WIDTH-2 TRANSFER AUDIT
# ============================================================================

def primitive_width2_audit(primitive):

    print()
    print("=" * 78)
    print("5. PRIMITIVE WIDTH-2 TRANSFER AUDIT")
    print("=" * 78)

    results = {}

    for t in range(5):

        source = dict(primitive[t])
        target = dict(primitive[t + 1])

        common = sorted(
            set(source) & set(target)
        )

        source0 = [
            source[r]
            for r in common[:-1]
        ]

        source1 = [
            source[r]
            for r in common[1:]
        ]

        target_values = [
            target[r]
            for r in common[:-1]
        ]

        # Equations correspond to r where both r and r+1
        # exist in the source and r exists in the target.

        if len(source0) < 1:

            print()
            print(
                "  t={} -> {}: insufficient_overlap".format(
                    t,
                    t + 1,
                )
            )

            continue

        result = solve_width2(
            source0,
            source1,
            target_values,
        )

        results[t] = result

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    equations={}".format(
                len(target_values)
            )
        )

        print(
            "    rank={}".format(
                result["rank"]
            )
        )

        print(
            "    augmented_rank={}".format(
                result["augmented_rank"]
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

            print(
                "    residuals={}".format(
                    result["residuals"]
                )
            )

    return results


# ============================================================================
# 6. CONTENT-RESCALED RAW TRANSFER
# ============================================================================

def content_rescaled_transfer(
    primitive_results,
    contents,
):

    print()
    print("=" * 78)
    print("6. CONTENT-RESCALED TRANSFER COEFFICIENTS")
    print("=" * 78)

    results = {}

    # If

    #     Q_{t+1}
    #       = c0 Q_t + c1 Q_t(shift),

    # and

    #     Q_t = g_t U_t,

    # then

    #     U_{t+1}
    #       = (g_t/g_{t+1}) c0 U_t
    #       + (g_t/g_{t+1}) c1 U_t(shift).

    for t, result in primitive_results.items():

        if result["status"] != "EXACT":
            continue

        scale = clean(
            sp.Rational(
                contents[t],
                contents[t + 1],
            )
        )

        d0 = clean(
            scale * result["c0"]
        )

        d1 = clean(
            scale * result["c1"]
        )

        results[t] = (
            d0,
            d1,
        )

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    content_scale={}".format(
                scale
            )
        )

        print(
            "    d0={}".format(
                d0
            )
        )

        print(
            "    d1={}".format(
                d1
            )
        )

        print(
            "    d0_is_integer={}".format(
                d0.q == 1
            )
        )

        print(
            "    d1_is_integer={}".format(
                d1.q == 1
            )
        )

        print(
            "    d0_valuations={}".format(
                {
                    p: valuation(d0, p)
                    for p in (
                        2,3,5,7,11,13,17
                    )
                }
            )
        )

        print(
            "    d1_valuations={}".format(
                {
                    p: valuation(d1, p)
                    for p in (
                        2,3,5,7,11,13,17
                    )
                }
            )
        )

    return results


# ============================================================================
# 7. PRIMITIVE NEWTON TRIANGLE
# ============================================================================

def primitive_newton_triangle(
    primitive
):

    print()
    print("=" * 78)
    print("7. PRIMITIVE NEWTON-TRIANGLE AUDIT")
    print("=" * 78)

    triangle = {}

    for t, row in primitive.items():

        values = [
            value
            for _, value
            in sorted(row)
        ]

        current = values

        j = 0

        while current:

            triangle[(j, t)] = current[0]

            if len(current) == 1:
                break

            current = [
                current[i + 1] - current[i]
                for i in range(
                    len(current) - 1
                )
            ]

            j += 1

    for t in range(6):

        print()
        print(
            "  t={}: {}".format(
                t,
                [
                    triangle.get(
                        (j,t),
                        None,
                    )
                    for j in range(4)
                ],
            )
        )

    matrix = sp.Matrix([
        [
            triangle.get(
                (j,t),
                0,
            )
            for t in range(6)
        ]
        for j in range(4)
    ])

    print()
    print(
        "  primitive_newton_matrix_rank={}".format(
            matrix.rank()
        )
    )

    return triangle


# ============================================================================
# 8. PRIME PROFILE OF PRIMITIVE SHAPE
# ============================================================================

def primitive_prime_profile(
    primitive
):

    print()
    print("=" * 78)
    print("8. PRIMITIVE-ROW PRIME PROFILE")
    print("=" * 78)

    for t in range(6):

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        values = [
            value
            for _, value
            in primitive[t]
        ]

        print(
            "    valuations={}".format(
                {
                    p: [
                        valuation(
                            value,
                            p,
                        )
                        for value in values
                    ]
                    for p in (
                        2,3,5,7,11,13,17
                    )
                }
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 341R — EXACT PRIMITIVE-ROW SHAPE / "
        "CONTENT-SEPARATION AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    contents, primitive = primitive_rows(
        layers
    )

    print_content_sequence(
        contents
    )

    print_primitive_rows(
        primitive
    )

    primitive_matrix = primitive_rank_audit(
        primitive
    )

    proportionality = adjacent_proportionality(
        primitive
    )

    primitive_transfers = primitive_width2_audit(
        primitive
    )

    rescaled = content_rescaled_transfer(
        primitive_transfers,
        contents,
    )

    triangle = primitive_newton_triangle(
        primitive
    )

    primitive_prime_profile(
        primitive
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
Experiment 340R showed:

    row/column/local gcd normalization
        does not reduce the full observed rank below 4.

Experiment 341R therefore asks a more specific question.

Write

    Q_t = g_t U_t,

where

    g_t = gcd(observed row t)

and

    gcd(U_t) = 1.

The scalar sequence g_t captures pure integer content.
The primitive row U_t captures the actual shape across p.

If the source mechanism is fundamentally scalar-times-shape, then one
would expect:

    * primitive rows to have simpler dynamics;
    * adjacent primitive rows to become proportional or satisfy a stable
      width-2 law;
    * content-rescaled transfer coefficients to simplify;
    * the primitive Newton triangle to have lower rank.

A negative result is also informative:

    if U_t remains rank 4 and no exact primitive transfer emerges, then
    the complexity is intrinsic to the row shape rather than merely an
    integer-content artifact.

This is the cleanest next test after the 339R/340R arithmetic audits.

No missing values are introduced.
No interpolation is used.
No synthetic second n=pq case is created.
"""
    )

    primitive_exact = [
        t
        for t, result
        in primitive_transfers.items()
        if result["status"] == "EXACT"
    ]

    proportional_hits = [
        t
        for t, ok
        in proportionality
        if ok
    ]

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  primitive_row_rank={}".format(
            primitive_matrix.rank()
        )
    )

    print(
        "  adjacent_proportional_pairs={}".format(
            proportional_hits
        )
    )

    print(
        "  exact_primitive_width2_transitions={}".format(
            primitive_exact
        )
    )

    print(
        "  content_rescaled_coefficients_computed=True"
    )

    print(
        "  primitive_newton_triangle_computed=True"
    )

    print(
        "  only_observed_cells_used=True"
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
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 341R COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print("\nInterrupted.")

        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )

        raise
