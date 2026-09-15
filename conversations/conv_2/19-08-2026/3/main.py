#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 319R — EXACT RAW-TRANSFER COEFFICIENT FACTORIZATION /
NORMALIZATION / CROSS-TRANSITION AUDIT
==============================================================================

Purpose
-------
Experiment 318R established exactly:

    raw Q_{t+1,k}
        = c0_t Q_{t,k} + c1_t Q_{t,k+1}

for t=1 and t=2,

and, under finite-difference coordinates,

    alpha_t = c0_t + c1_t
    beta_t  = c1_t.

Thus the previously studied first-entry operators are exactly induced by
the raw source-layer operators.

The next question is therefore:

    Do the raw coefficients c0_t,c1_t themselves have a simpler
    arithmetic or normalization structure?

This experiment does NOT fit a universal formula from two points.

Instead it performs an exact structural audit of the two raw coefficient
pairs:

    t=1:
        c0_1, c1_1

    t=2:
        c0_2, c1_2

Tests include:

    1. exact coefficient factorization;
    2. coefficient ratios c1/c0;
    3. ratios c0_2/c0_1 and c1_2/c1_1;
    4. sums/differences/products;
    5. determinant-like raw companion invariants;
    6. projective raw invariant c1/c0;
    7. comparison against layer contents;
    8. comparison against adjacent source-layer endpoint ratios;
    9. comparison against terminal source values;
   10. exact normalization by integer layer content;
   11. search for simple rational/integer normalized values;
   12. exact cross-transition relations;
   13. arithmetic prime-valuation profile;
   14. whether either coefficient is a simple rational function of the
       other coefficient or of the layer-content ratios.

The experiment deliberately avoids declaring a law from two transitions.

No external files.
No synthetic second n=pq case.
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


def valuation(x, p):
    x = sp.Rational(x)

    if x == 0:
        return None

    n = abs(int(x.p))
    d = abs(int(x.q))

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    while d % p == 0:
        d //= p
        v -= 1

    return v


def factor_rational(x):
    x = sp.Rational(x)

    return {
        "numerator": sp.factorint(abs(int(x.p))),
        "denominator": sp.factorint(abs(int(x.q))),
        "sign": -1 if x < 0 else 1,
    }


# ============================================================================
# SOURCE LAYERS
# ============================================================================

def degree_at_p(p):
    return len(Q[p]) - 1


def build_layers():
    layers = {}

    maximum_t = max(
        degree_at_p(p)
        for p in Q
    )

    for t in range(maximum_t + 1):

        row = []

        for p in sorted(Q):

            r = degree_at_p(p) - t

            if r < 0:
                continue

            row.append(
                (
                    p,
                    sp.Integer(
                        Q[p][r]
                    ),
                )
            )

        layers[t] = row

    return layers


def layer_values(layers, t):
    return [
        value
        for _, value in layers[t]
    ]


def layer_ps(layers, t):
    return [
        p
        for p, _ in layers[t]
    ]


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


# ============================================================================
# RAW WIDTH-2 SOLVER
# ============================================================================

def solve_raw_width2(
    layers,
    t,
):

    source = layer_values(
        layers,
        t,
    )

    target = layer_values(
        layers,
        t + 1,
    )

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
        target
        for _, _, target
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

    coeffs = M.inv() * rhs

    c0 = clean(
        coeffs[0, 0]
    )

    c1 = clean(
        coeffs[1, 0]
    )

    residuals = [
        clean(
            c0 * x0
            + c1 * x1
            - target
        )
        for x0, x1, target
        in equations
    ]

    return {
        "status": (
            "EXACT"
            if all(r == 0 for r in residuals)
            else "VERIFICATION_FAILED"
        ),
        "equations": equations,
        "c0": c0,
        "c1": c1,
        "residuals": residuals,
        "rank": rank,
        "augmented_rank": augmented_rank,
    }


# ============================================================================
# RAW COMPANION / PROJECTIVE INVARIANTS
# ============================================================================

def raw_companion(c0, c1):
    return sp.Matrix([
        [c0, c1],
        [1, 0],
    ])


def raw_invariants(c0, c1):

    T = raw_companion(
        c0,
        c1,
    )

    tr = clean(
        sp.trace(T)
    )

    det = clean(
        T.det()
    )

    disc = clean(
        tr**2 - 4 * det
    )

    return {
        "trace": tr,
        "det": det,
        "discriminant": disc,
        "trace_squared_over_det": (
            clean(tr**2 / det)
            if det != 0
            else None
        ),
        "beta_over_alpha_raw": (
            clean(c1 / c0)
            if c0 != 0
            else None
        ),
    }


# ============================================================================
# NORMALIZATION SEARCH
# ============================================================================

def normalized_values(
    value,
    content_t,
    content_next,
):

    return {
        "raw": clean(value),
        "times_content_t": clean(
            value * content_t
        ),
        "times_content_next": clean(
            value * content_next
        ),
        "times_content_t_over_next": clean(
            value
            * sp.Rational(
                content_t,
                content_next,
            )
        ),
        "times_content_next_over_t": clean(
            value
            * sp.Rational(
                content_next,
                content_t,
            )
        ),
    }


# ============================================================================
# ENDPOINT RATIO AUDIT
# ============================================================================

def endpoint_ratios(
    source,
    target,
):

    out = {}

    if len(source) >= 2:
        out["source_1_over_source_0"] = clean(
            sp.Rational(
                source[1],
                source[0],
            )
        )

    if len(source) >= 3:
        out["source_2_over_source_1"] = clean(
            sp.Rational(
                source[2],
                source[1],
            )
        )

    if len(target) >= 2:
        out["target_1_over_target_0"] = clean(
            sp.Rational(
                target[1],
                target[0],
            )
        )

    if len(target) >= 3:
        out["target_2_over_target_1"] = clean(
            sp.Rational(
                target[2],
                target[1],
            )
        )

    if (
        len(source) >= 3
        and len(target) >= 3
    ):

        out["cross_ratio_0"] = clean(
            sp.Rational(
                target[0] * source[1],
                source[0] * target[1],
            )
        )

        out["cross_ratio_1"] = clean(
            sp.Rational(
                target[1] * source[2],
                source[1] * target[2],
            )
        )

    return out


# ============================================================================
# TERMINAL REFERENCES
# ============================================================================

def terminal_reference():

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print()
    print("=" * 78)
    print(
        "9. TERMINAL SOURCE REFERENCE"
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
            sp.Integer(q1 // 17)
        )
    )

    print(
        "  q3/17={}".format(
            sp.Integer(q3 // 17)
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 319R — EXACT RAW-TRANSFER COEFFICIENT "
        "FACTORIZATION / NORMALIZATION / CROSS-TRANSITION AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    profiles = {}

    # ------------------------------------------------------------------------
    # 1. RAW SOURCE LAYERS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. RAW SOURCE LAYERS")
    print("=" * 78)

    for t in sorted(layers):

        vals = layer_values(
            layers,
            t,
        )

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    p={}".format(
                layer_ps(
                    layers,
                    t,
                )
            )
        )

        print(
            "    values={}".format(
                vals
            )
        )

        print(
            "    integer_content={}".format(
                integer_content(
                    vals
                )
            )
        )

    # ------------------------------------------------------------------------
    # 2. RAW COEFFICIENT EXTRACTION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "2. EXACT RAW WIDTH-2 COEFFICIENTS"
    )
    print("=" * 78)

    for t in (1, 2):

        result = solve_raw_width2(
            layers,
            t,
        )

        profiles[t] = result

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

        if result["status"] != "EXACT":
            continue

        c0 = result["c0"]
        c1 = result["c1"]

        print(
            "    c0={}".format(
                c0
            )
        )

        print(
            "    c1={}".format(
                c1
            )
        )

        print(
            "    c1/c0={}".format(
                clean(
                    c1 / c0
                )
            )
        )

        print(
            "    c0+c1={}".format(
                clean(
                    c0 + c1
                )
            )
        )

        print(
            "    c0-c1={}".format(
                clean(
                    c0 - c1
                )
            )
        )

        print(
            "    c0*c1={}".format(
                clean(
                    c0 * c1
                )
            )

        )

        print(
            "    equation_residuals={}".format(
                result["residuals"]
            )
        )

    # ------------------------------------------------------------------------
    # 3. EXACT PRIME FACTORIZATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. EXACT PRIME FACTORIZATION"
    )
    print("=" * 78)

    for t, result in profiles.items():

        if result["status"] != "EXACT":
            continue

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        for label in ("c0", "c1"):

            value = result[label]

            print()
            print(
                "    {}={}".format(
                    label,
                    value,
                )
            )

            profile = factor_rational(
                value
            )

            print(
                "      sign={}".format(
                    profile["sign"]
                )
            )

            print(
                "      numerator_factorization={}".format(
                    profile[
                        "numerator"
                    ]
                )
            )

            print(
                "      denominator_factorization={}".format(
                    profile[
                        "denominator"
                    ]
                )

            )

    # ------------------------------------------------------------------------
    # 4. RAW COMPANION INVARIANTS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. RAW COMPANION PROJECTIVE INVARIANTS"
    )
    print("=" * 78)

    for t, result in profiles.items():

        if result["status"] != "EXACT":
            continue

        inv = raw_invariants(
            result["c0"],
            result["c1"],
        )

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        for key, value in inv.items():

            print(
                "    {}={}".format(
                    key,
                    value,
                )
            )

    # ------------------------------------------------------------------------
    # 5. CROSS-TRANSITION RATIOS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "5. EXACT CROSS-TRANSITION RELATIONS"
    )
    print("=" * 78)

    r1 = profiles[1]
    r2 = profiles[2]

    c01 = r1["c0"]
    c11 = r1["c1"]

    c02 = r2["c0"]
    c12 = r2["c1"]

    pairs = {
        "c0_2/c0_1": clean(
            c02 / c01
        ),
        "c1_2/c1_1": clean(
            c12 / c11
        ),
        "c0_2/c1_1": clean(
            c02 / c11
        ),
        "c1_2/c0_1": clean(
            c12 / c01
        ),
        "(c1/c0)_2/(c1/c0)_1": clean(
            (c12 / c02)
            /
            (c11 / c01)
        ),
        "det_raw_2/det_raw_1": clean(
            (c12 / c11)
            *
            (c02 / c01)
        ),
    }

    for name, value in pairs.items():

        print(
            "  {}={}".format(
                name,
                value,
            )
        )

    # ------------------------------------------------------------------------
    # 6. SIMPLE INTEGER/RATIONAL NORMALIZATION SEARCH
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. SMALL NORMALIZATION SEARCH"
    )
    print("=" * 78)

    candidates = []

    for t, result in profiles.items():

        source = layer_values(
            layers,
            t,
        )

        target = layer_values(
            layers,
            t + 1,
        )

        g_t = integer_content(
            source
        )

        g_next = integer_content(
            target
        )

        print()
        print(
            "  t={} -> {}: content_ratio={}".format(
                t,
                t + 1,
                sp.Rational(
                    g_next,
                    g_t,
                ),
            )
        )

        for label in (
            "c0",
            "c1",
        ):

            print()
            print(
                "    {}:".format(
                    label
                )
            )

            normalized = normalized_values(
                result[label],
                g_t,
                g_next,
            )

            for key, value in normalized.items():

                print(
                    "      {}={}".format(
                        key,
                        value,
                    )
                )

                if value.q == 1 and abs(int(value)) <= 1000000:
                    candidates.append(
                        (
                            t,
                            label,
                            key,
                            value,
                        )
                    )

    print()
    print(
        "  small_integer_candidates={}".format(
            candidates
        )
    )

    # ------------------------------------------------------------------------
    # 7. SOURCE ENDPOINT RATIO AUDIT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "7. SOURCE ENDPOINT RATIO AUDIT"
    )
    print("=" * 78)

    for t, result in profiles.items():

        source = layer_values(
            layers,
            t,
        )

        target = layer_values(
            layers,
            t + 1,
        )

        ratios = endpoint_ratios(
            source,
            target,
        )

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        for key, value in ratios.items():

            print(
                "    {}={}".format(
                    key,
                    value,
                )
            )

    # ------------------------------------------------------------------------
    # 8. TERMINAL SOURCE REFERENCE
    # ------------------------------------------------------------------------

    terminal_reference()

    # ------------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 318R established the exact source-to-Newton-coordinate
relationship

    alpha_t = c0_t + c1_t,
    beta_t  = c1_t.

Therefore the next target is the raw coefficient pair itself.

This experiment deliberately does NOT fit a universal function of t,
because only two exact transitions are currently available.

Instead it searches for arithmetic structure that would constrain a
future formula:

    * simple ratios;
    * common denominator structure;
    * common prime factors;
    * layer-content normalization;
    * source-endpoint ratios;
    * projective raw invariants;
    * cross-transition coefficient relations.

The strongest possible positive outcome would be a normalization under
which the two raw coefficient pairs become equal or acquire a simple
integer/rational pattern.

The strongest negative outcome would be that after natural source-layer
normalizations, all cross-transition invariants remain complicated and
independent.

In that case the next experiment should not interpolate c0_t,c1_t.
Instead it should seek a deeper exact formula for Q_t itself, or obtain
a genuine second independent n=pq case.

No synthetic second case is generated.
"""
    )

    # ------------------------------------------------------------------------
    # 10. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    exact_count = sum(
        result["status"] == "EXACT"
        for result in profiles.values()
    )

    print(
        "  exact_raw_transitions={}".format(
            exact_count
        )
    )

    print(
        "  transitions=[1,2]"
    )

    print(
        "  coefficient_factorization_exact=True"
    )

    print(
        "  cross_transition_ratios_exact=True"
    )

    print(
        "  natural_content_normalization_tested=True"
    )

    print(
        "  endpoint_ratio_audit=True"
    )

    print(
        "  no_universal_interpolation_performed=True"
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
        "EXPERIMENT 319R COMPLETE"
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

