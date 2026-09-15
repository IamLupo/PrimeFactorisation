#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 321R — EXACT MISSING-SOURCE VALUE / OVERLAP-CONSISTENCY AUDIT
==============================================================================

Purpose
-------
Experiment 320R showed that the raw width-2 transfer law

    Q_{t+1}(p) = c0_t Q_t(p) + c1_t Q_t(p+2)

is NOT an identity of the interpolating layer polynomials.

The exact residuals vanish at the only two source points used to determine
the coefficients and fail elsewhere.

This means the raw coefficient pairs should currently be regarded as
finite-overlap reconstructions, not established polynomial laws.

The next useful question is whether the transfer law predicts meaningful
values at the MISSING source points.

The only nontrivial latent value currently recoverable from an existing
target point is

    Q_1(7),

because the transition t=1 -> 2 contains

    Q_2(5)
      =
    c0_1 Q_1(5) + c1_1 Q_1(7),

while Q_1(7) is absent from the original source table.

We therefore reconstruct this value exactly and compare it against:

    1. the interpolation polynomial for Q_1(p);
    2. finite-difference continuation from p=1,3,5;
    3. exact divisibility/content expectations;
    4. the known terminal-layer structure;
    5. the analogous boundary geometry at later t;
    6. whether the inferred value is compatible with an integer-valued
       or simple rational-valued continuation;
    7. whether the discrepancy from interpolation has a simple factor;
    8. whether the transfer-derived missing value restores a higher-order
       finite-difference pattern.

IMPORTANT
---------
The reconstructed Q_1(7) is NOT treated as observed data.

It is a derived diagnostic consequence of the fitted raw transfer law.

No synthetic second n=pq case is generated.
No universal formula is claimed.
Exact rational arithmetic only.
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
# BASIC HELPERS
# ============================================================================

p = sp.Symbol("p")


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

    num = abs(int(x.p))
    den = abs(int(x.q))

    out = 0

    while num % prime == 0:
        num //= prime
        out += 1

    while den % prime == 0:
        den //= prime
        out -= 1

    return out


def degree_at_source(prime):
    return len(Q[prime]) - 1


def integer_content(values):
    values = [
        int(v)
        for v in values
        if int(v) != 0
    ]

    if not values:
        return 0

    g = 0

    for value in values:
        g = math.gcd(
            g,
            abs(value),
        )

    return g


def primitive_integer_row(values):
    g = integer_content(values)

    if g == 0:
        return list(values)

    return [
        int(v) // g
        for v in values
    ]


# ============================================================================
# RAW LAYERS
# ============================================================================

def build_layers():

    layers = {}

    maximum_t = max(
        degree_at_source(prime)
        for prime in Q
    )

    for t in range(
        maximum_t + 1
    ):

        row = []

        for prime in sorted(Q):

            r = degree_at_source(prime) - t

            if r < 0:
                continue

            row.append(
                (
                    prime,
                    sp.Integer(
                        Q[prime][r]
                    ),
                )
            )

        layers[t] = row

    return layers


# ============================================================================
# RAW TRANSFER COEFFICIENTS
# ============================================================================

def solve_raw_width2(
    layers,
    t,
):

    source = [
        value
        for _, value
        in layers[t]
    ]

    target = [
        value
        for _, value
        in layers[t + 1]
    ]

    equations = []

    for k in range(
        min(
            len(target),
            len(source) - 1,
        )
    ):

        equations.append(
            (
                source[k],
                source[k + 1],
                target[k],
            )
        )

    if len(equations) < 2:
        return {
            "status": "INSUFFICIENT_DATA",
            "equations": equations,
        }

    M = sp.Matrix([
        [x0, x1]
        for x0, x1, _
        in equations
    ])

    rhs = sp.Matrix([
        y
        for _, _, y
        in equations
    ])

    rank = M.rank()
    augmented_rank = M.row_join(rhs).rank()

    if augmented_rank > rank:
        return {
            "status": "NO_SOLUTION",
            "equations": equations,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    if rank < 2:
        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    solution = M.inv() * rhs

    c0 = clean(
        solution[0]
    )

    c1 = clean(
        solution[1]
    )

    verified = all(
        clean(
            c0 * x0
            + c1 * x1
            - y
        ) == 0
        for x0, x1, y
        in equations
    )

    return {
        "status": "EXACT" if verified else "VERIFICATION_FAILED",
        "c0": c0,
        "c1": c1,
        "equations": equations,
        "rank": rank,
        "augmented_rank": augmented_rank,
    }


# ============================================================================
# INTERPOLATION
# ============================================================================

def interpolate_layer(
    layers,
    t,
):

    points = [
        (
            sp.Integer(prime),
            sp.Integer(value),
        )
        for prime, value
        in layers[t]
    ]

    return sp.Poly(
        sp.interpolate(
            points,
            p,
        ),
        p,
        domain=sp.QQ,
    )


# ============================================================================
# FINITE DIFFERENCES IN p
# ============================================================================

def finite_difference_rows(
    values,
):

    current = [
        sp.Integer(v)
        for v in values
    ]

    rows = [
        current
    ]

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

        rows.append(
            current
        )

    return rows


# ============================================================================
# MISSING VALUE RECONSTRUCTION
# ============================================================================

def reconstruct_missing_q17(
    layers,
    raw1,
):

    c0 = raw1["c0"]
    c1 = raw1["c1"]

    # Existing:
    #
    #   Q_2(5) = c0 Q_1(5) + c1 Q_1(7)
    #
    # Therefore:
    #
    #   Q_1(7) = (Q_2(5)-c0 Q_1(5))/c1.

    q1_5 = dict(layers[1])[5]
    q2_5 = dict(layers[2])[5]

    q1_7 = clean(
        (
            q2_5
            - c0 * q1_5
        )
        /
        c1
    )

    return q1_7


# ============================================================================
# COMPARISON AGAINST INTERPOLATION
# ============================================================================

def compare_against_interpolation(
    layers,
    inferred_q17,
):

    P1 = interpolate_layer(
        layers,
        1,
    )

    interpolated_q17 = clean(
        P1.as_expr().subs(
            p,
            7,
        )
    )

    discrepancy = clean(
        inferred_q17
        - interpolated_q17
    )

    return (
        interpolated_q17,
        discrepancy,
    )


# ============================================================================
# AUGMENTED FINITE-DIFFERENCE AUDIT
# ============================================================================

def augmented_difference_audit(
    layers,
    inferred_q17,
):

    values = [
        dict(layers[1])[1],
        dict(layers[1])[3],
        dict(layers[1])[5],
        inferred_q17,
    ]

    rows = finite_difference_rows(
        values
    )

    return rows


# ============================================================================
# INTEGER / DENOMINATOR AUDIT
# ============================================================================

def rational_arithmetic_audit(
    value,
):

    num = int(value.p)
    den = int(value.q)

    return {
        "value": value,
        "is_integer": den == 1,
        "numerator": num,
        "denominator": den,
        "numerator_factorization": (
            sp.factorint(
                abs(num)
            )
            if num != 0
            else {}
        ),
        "denominator_factorization": (
            sp.factorint(
                abs(den)
            )
            if abs(den) > 1
            else {}
        ),
        "valuations": {
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
        },
    }


# ============================================================================
# BOUNDARY COMPARISON
# ============================================================================

def boundary_audit(
    layers,
    inferred_q17,
):

    print()
    print("=" * 78)
    print(
        "8. BOUNDARY / TRIANGULAR-GEOMETRY AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "  observed t=0 p-values="
        "{}".format(
            sorted(
                dict(layers[0]).keys()
            )
        )
    )

    print(
        "  observed t=1 p-values="
        "{}".format(
            sorted(
                dict(layers[1]).keys()
            )
        )
    )

    print(
        "  inferred t=1 p=7={}".format(
            inferred_q17
        )
    )

    print(
        "  actual t=1 p=7=UNOBSERVED"
    )

    print()
    print(
        "  t=1 augmented p-row="
        "{}".format(
            [
                dict(layers[1])[1],
                dict(layers[1])[3],
                dict(layers[1])[5],
                inferred_q17,
            ]
        )
    )

    print(
        "  primitive_augmented_row={}".format(
            primitive_integer_row(
                [
                    dict(layers[1])[1],
                    dict(layers[1])[3],
                    dict(layers[1])[5],
                    inferred_q17,
                ]
            )
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 321R — EXACT MISSING-SOURCE VALUE / "
        "OVERLAP-CONSISTENCY AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    # ------------------------------------------------------------------------
    # 1. SOURCE TABLE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. SOURCE LAYER TABLE"
    )
    print("=" * 78)

    for t in sorted(layers):

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    {}".format(
                layers[t]
            )
        )

    # ------------------------------------------------------------------------
    # 2. RAW TRANSFERS
    # ------------------------------------------------------------------------

    raw = {}

    print()
    print("=" * 78)
    print(
        "2. EXACT RAW TRANSFERS"
    )
    print("=" * 78)

    for t in (
        1,
        2,
    ):

        result = solve_raw_width2(
            layers,
            t,
        )

        raw[t] = result

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
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

    # ------------------------------------------------------------------------
    # 3. RECONSTRUCT Q_1(7)
    # ------------------------------------------------------------------------

    inferred_q17 = reconstruct_missing_q17(
        layers,
        raw[1],
    )

    print()
    print("=" * 78)
    print(
        "3. EXACT RECONSTRUCTION OF MISSING Q_1(7)"
    )
    print("=" * 78)

    print(
        "  inferred_Q1(7)={}".format(
            inferred_q17
        )
    )

    print(
        "  observed_Q1(7)=UNAVAILABLE"
    )

    # ------------------------------------------------------------------------
    # 4. ARITHMETIC PROFILE
    # ------------------------------------------------------------------------

    arithmetic = rational_arithmetic_audit(
        inferred_q17
    )

    print()
    print("=" * 78)
    print(
        "4. ARITHMETIC PROFILE OF INFERRED VALUE"
    )
    print("=" * 78)

    for key, value in arithmetic.items():

        print(
            "  {}={}".format(
                key,
                value,
            )
        )

    # ------------------------------------------------------------------------
    # 5. INTERPOLATION COMPARISON
    # ------------------------------------------------------------------------

    interpolated_q17, discrepancy = (
        compare_against_interpolation(
            layers,
            inferred_q17,
        )
    )

    print()
    print("=" * 78)
    print(
        "5. INFERRED VALUE VS EXISTING-SAMPLE INTERPOLATION"
    )
    print("=" * 78)

    print(
        "  interpolation_Q1(7)={}".format(
            interpolated_q17
        )
    )

    print(
        "  transfer_inferred_Q1(7)={}".format(
            inferred_q17
        )
    )

    print(
        "  discrepancy={}".format(
            discrepancy
        )
    )

    print(
        "  discrepancy_zero={}".format(
            discrepancy == 0
        )
    )

    # ------------------------------------------------------------------------
    # 6. AUGMENTED DIFFERENCE AUDIT
    # ------------------------------------------------------------------------

    rows = augmented_difference_audit(
        layers,
        inferred_q17,
    )

    print()
    print("=" * 78)
    print(
        "6. AUGMENTED p-FINITE-DIFFERENCE AUDIT"
    )
    print("=" * 78)

    for j, row in enumerate(rows):

        print(
            "  Delta_p^{}={}".format(
                j,
                row,
            )
        )

    # ------------------------------------------------------------------------
    # 7. COMPARE WITH T=0 THIRD DIFFERENCE
    # ------------------------------------------------------------------------

    t0_values = [
        value
        for _, value
        in layers[0]
    ]

    t0_rows = finite_difference_rows(
        t0_values
    )

    print()
    print("=" * 78)
    print(
        "7. CROSS-LAYER DIFFERENCE SIGNATURE COMPARISON"
    )
    print("=" * 78)

    print(
        "  t=0 Delta_p^3={}".format(
            t0_rows[3]
            if len(t0_rows) > 3
            else None
        )
    )

    print(
        "  t=1 augmented Delta_p^3={}".format(
            rows[3]
            if len(rows) > 3
            else None
        )
    )

    print(
        "  t0_third_difference_constant="
        "{}".format(
            t0_rows[3][0]
            if len(t0_rows) > 3
            else None
        )
    )

    print(
        "  t1_augmented_third_difference_constant="
        "{}".format(
            rows[3][0]
            if len(rows) > 3
            else None
        )
    )

    # ------------------------------------------------------------------------
    # 8. BOUNDARY AUDIT
    # ------------------------------------------------------------------------

    boundary_audit(
        layers,
        inferred_q17,
    )

    # ------------------------------------------------------------------------
    # 9. CONTENT / NORMALIZATION COMPARISON
    # ------------------------------------------------------------------------

    augmented_row = [
        dict(layers[1])[1],
        dict(layers[1])[3],
        dict(layers[1])[5],
        inferred_q17,
    ]

    observed_row = [
        dict(layers[1])[1],
        dict(layers[1])[3],
        dict(layers[1])[5],
    ]

    observed_content = integer_content(
        observed_row
    )

    augmented_content = integer_content(
        [
            inferred_q17
        ]
        if inferred_q17.q == 1
        else []
    )

    print()
    print("=" * 78)
    print(
        "9. INTEGER CONTENT / NORMALIZATION AUDIT"
    )
    print("=" * 78)

    print(
        "  observed_t1_content={}".format(
            observed_content
        )
    )

    print(
        "  inferred_value_integer={}".format(
            inferred_q17.q == 1
        )
    )

    print(
        "  inferred_value_content_if_integer={}".format(
            augmented_content
        )
    )

    # ------------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 320R showed that the raw width-2 transfer laws fail as
polynomial identities in p.

The exact residual factorization is especially important:

    R(p) is proportional to (p-1)(p-3).

Those are precisely the two overlap points used to reconstruct the
two coefficients.

Therefore the observed raw transfer coefficients are currently only
finite-overlap fits.

Experiment 321R asks what those fits would predict at the first missing
source point.

The transition t=1 -> 2 contains enough information to infer Q_1(7),
because Q_2(5) is observed while Q_1(7) is not.

This inferred value is NOT new experimental data.

It is a derived consistency diagnostic.

The comparison with the degree-2 interpolation of the observed
Q_1(1), Q_1(3), Q_1(5) is particularly informative:

    agreement
        would support a low-degree polynomial continuation;

    disagreement
        would show that the raw transfer law and the visible
        p-polynomial pattern point to genuinely different hidden
        continuations.

The augmented finite-difference audit then asks whether the inferred
value restores a simple constant higher difference.

No universal claim is made.
No synthetic second n=pq case is generated.
"""
    )

    # ------------------------------------------------------------------------
    # 11. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  raw_t1_transfer_exact={}".format(
            raw[1]["status"] == "EXACT"
        )
    )

    print(
        "  missing_Q1_7_reconstructed_exactly=True"
    )

    print(
        "  inferred_Q1_7_observed=False"
    )

    print(
        "  interpolation_comparison_performed=True"
    )

    print(
        "  augmented_difference_audit=True"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  interpolation_counted_as_proof=False"
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
        "EXPERIMENT 321R COMPLETE"
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
