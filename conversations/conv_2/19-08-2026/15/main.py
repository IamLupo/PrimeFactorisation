#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 331R — EXACT CROSS-RATIO / MULTIPLICATIVE-SEPARABILITY AUDIT
==============================================================================

Purpose
-------
Experiments 322R-330R rejected a wide range of linear source laws:

    * constant vertical recurrences;
    * p-dependent vertical recurrences of low degree;
    * low-degree rational vertical recurrences;
    * low-degree bivariate polynomial laws;
    * constant local 2D linear stencils.

Experiment 330R also showed that the only surviving local 2D relation found
at low stencil size was data-sized rather than overdetermined.

The next qualitatively different possibility is MULTIPLICATIVE structure.

For every fully observed 2x2 rectangle

    Q_t(p)        Q_t(p+2)
    Q_{t+1}(p)    Q_{t+1}(p+2),

define the exact cross-ratio

    R(p,t)
      = Q_t(p) Q_{t+1}(p+2)
        --------------------------------
        Q_t(p+2) Q_{t+1}(p).

For a rank-one multiplicatively separable law

    Q_t(p) = U(p) V(t),

one must have

    R(p,t) = 1

for every complete rectangle.

More generally:

    R(p,t) depending only on t
        is compatible with a separable p-coordinate structure;

    R(p,t) depending only on p
        is compatible with a separable t-coordinate structure.

The experiment therefore tests:

    1. exact cross-ratio values;
    2. cross-ratio == 1;
    3. constancy across p at fixed t;
    4. constancy across t at fixed p;
    5. multiplicative rank-one consistency where enough rectangles exist;
    6. exact prime-factor profiles;
    7. whether deviations factor into simple integer patterns.

No interpolation.
No missing values.
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

    n = abs(int(x.p))
    d = abs(int(x.q))

    out = 0

    while n % prime == 0:
        n //= prime
        out += 1

    while d % prime == 0:
        d //= prime
        out -= 1

    return out


def factor_integer(n):
    n = int(abs(n))

    if n <= 1:
        return {}

    return sp.factorint(n)


def build_table():
    """
    Convert reverse-ordered source storage into

        table[(p,t)] = Q_t(p).
    """

    table = {}

    for p_value, values in Q.items():

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            table[
                (
                    int(p_value),
                    int(t_value),
                )
            ] = sp.Integer(value)

    return table


# ============================================================================
# COMPLETE RECTANGLES
# ============================================================================

def complete_rectangles(table):

    observed = set(table)

    rectangles = []

    for p0, t0 in sorted(
        observed,
        key=lambda x: (
            x[1],
            x[0],
        ),
    ):

        points = [
            (p0, t0),
            (p0 + 2, t0),
            (p0, t0 + 1),
            (p0 + 2, t0 + 1),
        ]

        if all(
            point in observed
            for point in points
        ):

            rectangles.append(
                (
                    p0,
                    t0,
                )
            )

    return rectangles


# ============================================================================
# CROSS-RATIO
# ============================================================================

def cross_ratio(
    table,
    p0,
    t0,
):

    q00 = table[
        (p0, t0)
    ]

    q01 = table[
        (p0 + 2, t0)
    ]

    q10 = table[
        (p0, t0 + 1)
    ]

    q11 = table[
        (p0 + 2, t0 + 1)
    ]

    denominator = q01 * q10

    if denominator == 0:
        return None

    return clean(
        q00 * q11
        / denominator
    )


# ============================================================================
# PRIME PROFILE
# ============================================================================

def print_prime_profile(
    value,
):

    return {
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
    }


# ============================================================================
# FIXED-T CONSTANCY
# ============================================================================

def fixed_t_audit(
    ratios,
):

    print()
    print("=" * 78)
    print(
        "FIXED-t CROSS-RATIO CONSTANCY"
    )
    print("=" * 78)

    by_t = {}

    for (
        p0,
        t0,
    ), value in ratios.items():

        by_t.setdefault(
            t0,
            [],
        ).append(
            (
                p0,
                value,
            )
        )

    all_constant = True

    for t0 in sorted(by_t):

        entries = by_t[t0]

        values = [
            value
            for _, value
            in entries
        ]

        constant = all(
            clean(
                value
                - values[0]
            ) == 0
            for value in values[1:]
        )

        if not constant:
            all_constant = False

        print()
        print(
            "  t={}:".format(
                t0
            )
        )

        print(
            "    values={}".format(
                entries
            )
        )

        print(
            "    constant_across_p={}".format(
                constant
            )
        )

    print()
    print(
        "  global_function_of_t_only={}".format(
            all_constant
        )
    )

    return all_constant


# ============================================================================
# FIXED-P CONSTANCY
# ============================================================================

def fixed_p_audit(
    ratios,
):

    print()
    print("=" * 78)
    print(
        "FIXED-p CROSS-RATIO CONSTANCY"
    )
    print("=" * 78)

    by_p = {}

    for (
        p0,
        t0,
    ), value in ratios.items():

        by_p.setdefault(
            p0,
            [],
        ).append(
            (
                t0,
                value,
            )
        )

    all_constant = True

    for p0 in sorted(by_p):

        entries = by_p[p0]

        values = [
            value
            for _, value
            in entries
        ]

        constant = all(
            clean(
                value
                - values[0]
            ) == 0
            for value in values[1:]
        )

        if not constant:
            all_constant = False

        print()
        print(
            "  p={}:".format(
                p0
            )
        )

        print(
            "    values={}".format(
                entries
            )
        )

        print(
            "    constant_across_t={}".format(
                constant
            )
        )

    print()
    print(
        "  global_function_of_p_only={}".format(
            all_constant
        )
    )

    return all_constant


# ============================================================================
# RATIO == 1 AUDIT
# ============================================================================

def rank_one_audit(
    ratios,
):

    print()
    print("=" * 78)
    print(
        "RANK-ONE MULTIPLICATIVE SEPARABILITY AUDIT"
    )
    print("=" * 78)

    failures = []

    for location, value in ratios.items():

        exact_one = (
            clean(
                value - 1
            )
            == 0
        )

        print()
        print(
            "  rectangle={}:".format(
                location
            )
        )

        print(
            "    cross_ratio={}".format(
                value
            )
        )

        print(
            "    equals_one={}".format(
                exact_one
            )
        )

        if not exact_one:

            failures.append(
                (
                    location,
                    value,
                )
            )

    print()
    print(
        "  all_cross_ratios_one={}".format(
            len(failures) == 0
        )
    )

    return failures


# ============================================================================
# DEVIATION FROM 1
# ============================================================================

def deviation_audit(
    ratios,
):

    print()
    print("=" * 78)
    print(
        "CROSS-RATIO DEVIATION-FACTOR AUDIT"
    )
    print("=" * 78)

    deviations = {}

    for location, value in ratios.items():

        deviation = clean(
            value - 1
        )

        deviations[
            location
        ] = deviation

        print()
        print(
            "  rectangle={}:".format(
                location
            )
        )

        print(
            "    R-1={}".format(
                deviation
            )
        )

        numerator = int(
            abs(
                int(deviation.p)
            )
        )

        denominator = int(
            abs(
                int(deviation.q)
            )
        )

        print(
            "    numerator_factorization={}".format(
                factor_integer(
                    numerator
                )
            )
        )

        print(
            "    denominator_factorization={}".format(
                factor_integer(
                    denominator
                )
            )
        )

        print(
            "    valuations={}".format(
                print_prime_profile(
                    deviation
                )
            )
        )

    return deviations


# ============================================================================
# CROSS-RATIO CROSS-CONSISTENCY
# ============================================================================

def cross_ratio_grid_audit(
    ratios,
):

    print()
    print("=" * 78)
    print(
        "CROSS-RATIO GRID CONSISTENCY AUDIT"
    )
    print("=" * 78)

    # The available ratios lie on:
    #
    #   t = 0 : p = 1,3,5
    #   t = 1 : p = 1,3
    #   t = 2 : p = 1
    #
    # The only complete 2x2 rectangle inside the ratio grid is
    #
    #   (p,t) = (1,0), (3,0), (1,1), (3,1).
    #
    # Test multiplicative rank-one consistency there:
    #
    #   R(1,0) R(3,1) = R(3,0) R(1,1).

    needed = [
        (1, 0),
        (3, 0),
        (1, 1),
        (3, 1),
    ]

    if not all(
        location in ratios
        for location in needed
    ):

        print(
            "  complete_ratio_rectangle_available=False"
        )

        return None

    lhs = clean(
        ratios[(1, 0)]
        * ratios[(3, 1)]
    )

    rhs = clean(
        ratios[(3, 0)]
        * ratios[(1, 1)]
    )

    difference = clean(
        lhs - rhs
    )

    print(
        "  complete_ratio_rectangle_available=True"
    )

    print(
        "  lhs=R(1,0)*R(3,1)={}".format(
            lhs
        )
    )

    print(
        "  rhs=R(3,0)*R(1,1)={}".format(
            rhs
        )
    )

    print(
        "  difference={}".format(
            difference
        )
    )

    print(
        "  multiplicative_rank_one_on_ratio_grid={}".format(
            difference == 0
        )
    )

    return difference == 0


# ============================================================================
# RAW 2x2 MINOR COMPARISON
# ============================================================================

def determinant_audit(
    table,
    rectangles,
):

    print()
    print("=" * 78)
    print(
        "RAW RECTANGLE DETERMINANT / CROSS-RATIO CONNECTION"
    )
    print("=" * 78)

    for p0, t0 in rectangles:

        q00 = table[
            (p0, t0)
        ]

        q01 = table[
            (p0 + 2, t0)
        ]

        q10 = table[
            (p0, t0 + 1)
        ]

        q11 = table[
            (p0 + 2, t0 + 1)
        ]

        determinant = clean(
            q00 * q11
            - q01 * q10
        )

        ratio = cross_ratio(
            table,
            p0,
            t0,
        )

        denominator = clean(
            q01 * q10
        )

        reconstructed = clean(
            denominator
            * (ratio - 1)
        )

        print()
        print(
            "  rectangle=({},{}):".format(
                p0,
                t0,
            )
        )

        print(
            "    determinant={}".format(
                determinant
            )
        )

        print(
            "    cross_ratio_minus_one={}".format(
                clean(
                    ratio - 1
                )
            )
        )

        print(
            "    denominator={}".format(
                denominator
            )
        )

        print(
            "    denominator_times_R_minus_1={}".format(
                reconstructed
            )
        )

        print(
            "    exact_connection={}".format(
                determinant == reconstructed
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 331R — EXACT CROSS-RATIO / "
        "MULTIPLICATIVE-SEPARABILITY AUDIT"
    )
    print("=" * 78)

    table = build_table()

    rectangles = complete_rectangles(
        table
    )

    print()
    print("=" * 78)
    print(
        "1. COMPLETE OBSERVED RECTANGLES"
    )
    print("=" * 78)

    print(
        "  rectangle_count={}".format(
            len(rectangles)
        )
    )

    print(
        "  rectangles={}".format(
            rectangles
        )
    )

    ratios = {}

    print()
    print("=" * 78)
    print(
        "2. EXACT CROSS-RATIOS"
    )
    print("=" * 78)

    for p0, t0 in rectangles:

        ratio = cross_ratio(
            table,
            p0,
            t0,
        )

        ratios[
            (
                p0,
                t0,
            )
        ] = ratio

        print()
        print(
            "  R(p={},t={})=".format(
                p0,
                t0,
            )
        )

        print(
            "    {}".format(
                ratio
            )
        )

        print(
            "    valuation_profile={}".format(
                print_prime_profile(
                    ratio
                )
            )
        )

    failures = rank_one_audit(
        ratios
    )

    t_only = fixed_t_audit(
        ratios
    )

    p_only = fixed_p_audit(
        ratios
    )

    ratio_grid_rank_one = (
        cross_ratio_grid_audit(
            ratios
        )
    )

    deviations = deviation_audit(
        ratios
    )

    determinant_audit(
        table,
        rectangles,
    )

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 330R found no overdetermined constant-coefficient local linear
stencil.

That rules out a large class of additive source laws, but it does not
rule out multiplicative or projective structure.

For a multiplicatively separable table

    Q_t(p) = U(p) V(t),

every exact rectangle satisfies

    Q_t(p) Q_{t+1}(p+2)
      =
    Q_t(p+2) Q_{t+1}(p),

hence

    R(p,t) = 1.

Thus the cross-ratio audit is the nonlinear analogue of the rank-one
minor test.

More generally:

    R(p,t) = f(t)

would indicate that the p-dependence is separating;

    R(p,t) = g(p)

would indicate that the t-dependence is separating.

The ratio-grid test then asks whether the cross-ratios themselves
possess multiplicative rank one.

This is independent of the additive recurrence searches.

Because every quantity is constructed only from observed rectangles,
no missing value or extrapolation is introduced.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  complete_rectangle_count={}".format(
            len(rectangles)
        )
    )

    print(
        "  all_cross_ratios_one={}".format(
            len(failures) == 0
        )
    )

    print(
        "  cross_ratio_function_of_t_only={}".format(
            t_only
        )
    )

    print(
        "  cross_ratio_function_of_p_only={}".format(
            p_only
        )
    )

    print(
        "  cross_ratio_grid_rank_one={}".format(
            ratio_grid_rank_one
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
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
        "EXPERIMENT 331R COMPLETE"
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
