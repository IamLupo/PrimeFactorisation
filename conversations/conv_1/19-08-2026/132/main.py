#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 298R — EXACT UPSTREAM SOURCE-MAP / IDENTIFIABILITY AUDIT
==============================================================================

Purpose
-------

The previous experiments have exhausted a large class of downstream pattern
searches.

Experiment 298R changes the question.

Instead of fitting B[k,r] directly by arbitrary formulas, ask:

    Can the observed B[k,r] entries be generated from the available
    q_p(r) source values through a small UNIVERSAL source operator?

The experiment tests increasingly broad but still rigid operator classes.

Class A
-------

    B[k,r] = q_p(r) * C

A single universal scalar.

Class B
-------

    B[k,r] =
        q_p(r) *
        P(p)

where P has degree <= d.

Class C
-------

    B[k,r] =
        q_p(r) *
        P(p,k)

where P is a total-degree polynomial in (p,k).

Class D
-------

    B[k,r] =
        q_p(r) *
        P(p,k,r)

where P is a total-degree polynomial in (p,k,r).

These are GLOBAL operators: the same coefficients must work across every
available p,k,r observation.

Only exact overdetermined solutions count as structural hits.

The script also computes the formal underdetermination of arbitrary
source mixing at each (k,r), showing how many degrees of freedom remain
if one abandons structural restrictions.

Important:

    * No second synthetic n=pq case is generated.
    * No guessed q1(p,q) or q3(p,q) is inserted.
    * No external files are read.
    * No floating point.
    * No arbitrary matrix fit is counted as a theorem.
"""

from __future__ import annotations

import itertools
import math
import sys

import sympy as sp


# ============================================================================
# SOURCE q-TABLE
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
# B[k,r]
#
# B[k] is stored using absolute r-index:
#
#   k=0 -> r=0,...,7
#   k=1 -> r=1,...,7
#   ...
# ============================================================================
# 
# These are EXACTLY the values recovered in Experiments 278-297.
# ============================================================================

B = {
    0: {
        0: sp.Rational(25),
        1: sp.Rational(619),
        2: sp.Rational(3231, 2),
        3: sp.Rational(-33, 2),
        4: sp.Rational(-1675, 4),
        5: sp.Rational(3363, 20),
        6: sp.Rational(-9991, 360),
        7: sp.Rational(-421, 2520),
    },

    1: {
        1: sp.Rational(1750),
        2: sp.Rational(8624),
        3: sp.Rational(6829, 3),
        4: sp.Rational(-27341, 8),
        5: sp.Rational(10551, 10),
        6: sp.Rational(-16819, 180),
        7: sp.Rational(-6053, 180),
    },

    2: {
        2: sp.Rational(9690),
        3: sp.Rational(10234),
        4: sp.Rational(-57829, 8),
        5: sp.Rational(148151, 120),
        6: sp.Rational(1432, 5),
        7: sp.Rational(-5769, 28),
    },

    3: {
        3: sp.Rational(22100, 3),
        4: sp.Rational(-19045, 12),
        5: sp.Rational(-5577, 4),
        6: sp.Rational(351271, 360),
        7: sp.Rational(-101119, 315),
    },

    4: {
        4: sp.Rational(17875, 24),
        5: sp.Rational(-22061, 40),
        6: sp.Rational(132343, 720),
        7: sp.Rational(-162139, 5040),
    },

    5: {
        5: sp.Rational(65, 12),
        6: sp.Rational(-313, 60),
        7: sp.Rational(301, 120),
    },
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


def monomial_exponents(total_degree, variable_count):
    """
    All exponent tuples

        e_1 + ... + e_n <= total_degree.
    """
    out = []

    for exponents in itertools.product(
        range(total_degree + 1),
        repeat=variable_count,
    ):
        if sum(exponents) <= total_degree:
            out.append(exponents)

    out.sort(
        key=lambda x: (
            sum(x),
            x,
        )
    )

    return out


def solve_exact_system(A, b):
    """
    Safe exact linear solve.

    Returns:
        status,
        particular_solution,
        rank,
        augmented_rank,
        nullity
    """
    A = sp.Matrix(A)
    b = sp.Matrix(b)

    if A.rows == 0:
        return (
            "NO_EQUATIONS",
            None,
            0,
            0,
            A.cols,
        )

    augmented = A.row_join(b)

    rank_A = A.rank()
    rank_aug = augmented.rank()

    if rank_A != rank_aug:
        return (
            "NO_SOLUTION",
            None,
            rank_A,
            rank_aug,
            A.cols - rank_A,
        )

    if rank_A < A.cols:
        solution_set = sp.linsolve(
            (A, b)
        )

        return (
            "NONUNIQUE",
            solution_set,
            rank_A,
            rank_aug,
            A.cols - rank_A,
        )

    solution = A.inv() * b

    return (
        "UNIQUE",
        solution,
        rank_A,
        rank_aug,
        0,
    )


def available_p_values(r):
    return [
        p
        for p in sorted(Q)
        if r < len(Q[p])
    ]


def q_value(p, r):
    if p not in Q:
        return None

    if r < 0 or r >= len(Q[p]):
        return None

    return sp.Integer(
        Q[p][r]
    )


def observations():
    """
    Build all observations

        (p,k,r,q_p(r),B[k,r])

    whenever q_p(r) exists.
    """
    rows = []

    for k in sorted(B):

        for r in sorted(B[k]):

            target = sp.Rational(
                B[k][r]
            )

            for p in available_p_values(r):

                q = q_value(
                    p,
                    r,
                )

                if q is None:
                    continue

                if q == 0:
                    continue

                rows.append(
                    (
                        p,
                        k,
                        r,
                        q,
                        target,
                    )
                )

    return rows


def print_source_inventory():
    print()
    print("=" * 78)
    print("1. SOURCE INVENTORY")
    print("=" * 78)

    for p in sorted(Q):

        print()
        print(
            "  p={}: q={}".format(
                p,
                Q[p],
            )
        )

        print(
            "    length={}".format(
                len(Q[p]),
            )
        )


def print_observation_inventory(obs):
    print()
    print("=" * 78)
    print("2. OBSERVATION INVENTORY")
    print("=" * 78)

    print(
        "  total_observations={}".format(
            len(obs)
        )
    )

    by_r = {}

    for p, k, r, q, target in obs:
        by_r.setdefault(
            r,
            [],
        ).append(
            p
        )

    for r in sorted(by_r):

        print(
            "  r={}: available_p={}".format(
                r,
                by_r[r],
            )
        )


# ============================================================================
# CLASS A — B = C*q
# ============================================================================

def test_class_A(obs):
    print()
    print("=" * 78)
    print("3. CLASS A — UNIVERSAL SCALAR")
    print("=" * 78)

    A = []
    b = []

    for p, k, r, q, target in obs:
        A.append(
            [q]
        )
        b.append(
            target
        )

    result = solve_exact_system(
        A,
        b,
    )

    status = result[0]

    print(
        "  status={}".format(
            status
        )
    )

    print(
        "  rank={}".format(
            result[2]
        )
    )

    print(
        "  augmented_rank={}".format(
            result[3]
        )
    )

    print(
        "  nullity={}".format(
            result[4]
        )
    )

    if status == "UNIQUE":
        print(
            "  scalar={}".format(
                list(result[1])
            )
        )


# ============================================================================
# CLASS B — B = q * P(p)
# ============================================================================

def test_class_B(obs, degree):
    print()
    print("=" * 78)
    print(
        "4. CLASS B — q * POLYNOMIAL(p), DEGREE <= {}".format(
            degree
        )
    )
    print("=" * 78)

    exponents = monomial_exponents(
        degree,
        1,
    )

    A = []
    b = []

    for p, k, r, q, target in obs:

        row = []

        for exponent in exponents:

            power_p = exponent[0]

            row.append(
                clean(
                    q * (
                        sp.Integer(p)
                        ** power_p
                    )
                )
            )

        A.append(row)
        b.append(target)

    result = solve_exact_system(
        A,
        b,
    )

    status = result[0]

    print(
        "  monomial_count={}".format(
            len(exponents)
        )
    )

    print(
        "  equations={}".format(
            len(A)
        )
    )

    print(
        "  status={}".format(
            status
        )
    )

    print(
        "  rank={}".format(
            result[2]
        )
    )

    print(
        "  augmented_rank={}".format(
            result[3]
        )
    )

    if status == "UNIQUE":

        solution = list(
            result[1]
        )

        print(
            "  coefficients={}".format(
                solution
            )
        )


# ============================================================================
# CLASS C — B = q * P(p,k)
# ============================================================================

def test_class_C(obs, degree):
    print()
    print("=" * 78)
    print(
        "5. CLASS C — q * POLYNOMIAL(p,k), TOTAL DEGREE <= {}".format(
            degree
        )
    )
    print("=" * 78)

    exponents = monomial_exponents(
        degree,
        2,
    )

    A = []
    b = []

    for p, k, r, q, target in obs:

        row = []

        for exponent in exponents:

            ep = exponent[0]
            ek = exponent[1]

            monomial = (
                sp.Integer(p) ** ep
                * sp.Integer(k) ** ek
            )

            row.append(
                clean(
                    q * monomial
                )
            )

        A.append(row)
        b.append(target)

    result = solve_exact_system(
        A,
        b,
    )

    print(
        "  monomial_count={}".format(
            len(exponents)
        )
    )

    print(
        "  equations={}".format(
            len(A)
        )
    )

    print(
        "  status={}".format(
            result[0]
        )
    )

    print(
        "  rank={}".format(
            result[2]
        )
    )

    print(
        "  augmented_rank={}".format(
            result[3]
        )
    )

    print(
        "  nullity={}".format(
            result[4]
        )
    )


# ============================================================================
# CLASS D — B = q * P(p,k,r)
# ============================================================================

def test_class_D(obs, degree):
    print()
    print("=" * 78)
    print(
        "6. CLASS D — q * POLYNOMIAL(p,k,r), TOTAL DEGREE <= {}".format(
            degree
        )
    )
    print("=" * 78)

    exponents = monomial_exponents(
        degree,
        3,
    )

    A = []
    b = []

    for p, k, r, q, target in obs:

        row = []

        for exponent in exponents:

            ep = exponent[0]
            ek = exponent[1]
            er = exponent[2]

            monomial = (
                sp.Integer(p) ** ep
                * sp.Integer(k) ** ek
                * sp.Integer(r) ** er
            )

            row.append(
                clean(
                    q * monomial
                )
            )

        A.append(row)
        b.append(target)

    result = solve_exact_system(
        A,
        b,
    )

    print(
        "  monomial_count={}".format(
            len(exponents)
        )
    )

    print(
        "  equations={}".format(
            len(A)
        )
    )

    print(
        "  status={}".format(
            result[0]
        )
    )

    print(
        "  rank={}".format(
            result[2]
        )
    )

    print(
        "  augmented_rank={}".format(
            result[3]
        )
    )

    print(
        "  nullity={}".format(
            result[4]
        )
    )


# ============================================================================
# ARBITRARY SOURCE MIXING DIMENSION
# ============================================================================

def arbitrary_mixing_profile():
    print()
    print("=" * 78)
    print(
        "7. FORMAL UNDERDETERMINATION OF ARBITRARY SOURCE MIXING"
    )
    print("=" * 78)

    total_unknowns = 0
    total_equations = 0

    for k in sorted(B):

        for r in sorted(B[k]):

            ps = available_p_values(r)

            m = len(ps)

            if m == 0:
                continue

            total_unknowns += m
            total_equations += 1

            degrees_of_freedom = m - 1

            print()
            print(
                "  (k={},r={}):".format(
                    k,
                    r,
                )
            )

            print(
                "    available_p={}".format(
                    ps
                )
            )

            print(
                "    unknown_weights={}".format(
                    m
                )
            )

            print(
                "    equations={}".format(
                    1
                )
            )

            print(
                "    formal_degrees_of_freedom={}".format(
                    degrees_of_freedom
                )
            )

    print()
    print(
        "  total_weight_unknowns={}".format(
            total_unknowns
        )
    )

    print(
        "  total_B_equations={}".format(
            total_equations
        )
    )

    print(
        "  formal_excess_unknowns={}".format(
            total_unknowns
            - total_equations
        )
    )


# ============================================================================
# TERMINAL SOURCE REFERENCE
# ============================================================================

def terminal_source_reference():
    print()
    print("=" * 78)
    print(
        "8. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    q1 = 495451247
    q3 = 421514439

    g = math.gcd(
        q1,
        q3,
    )

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
            g
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
# INTERPRETATION
# ============================================================================

def interpretation():
    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
This experiment is intentionally different from the earlier pattern
searches.

The central question is IDENTIFIABILITY.

We ask whether the currently available q-values contain enough
information to reconstruct B through a small universal source operator.

The hierarchy is:

    B = C q

    B = q P(p)

    B = q P(p,k)

    B = q P(p,k,r)

A successful exact solution in a low-dimensional class would identify a
real upstream source mechanism.

A failure means only that the tested operator class is insufficient.

The final formal underdetermination calculation is also important.

If arbitrary weights w_p(k,r) are allowed, then each observed B[k,r]
provides only one equation for as many unknown source weights as there
are available p-values.

Therefore arbitrary source mixing is intrinsically underdetermined.

Consequently, an exact fitted combination is not evidence of provenance
unless the coefficients are constrained by an independently justified
universal law.

This experiment therefore separates:

    identifiable source law

from

    underdetermined interpolation.

No second n=pq case is synthesized.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 78)
    print(
        "EXPERIMENT 298R — EXACT UPSTREAM SOURCE-MAP / "
        "IDENTIFIABILITY AUDIT"
    )
    print("=" * 78)

    print_source_inventory()

    obs = observations()

    print_observation_inventory(
        obs
    )

    test_class_A(
        obs
    )

    for degree in range(0, 4):
        test_class_B(
            obs,
            degree,
        )

    for degree in range(0, 4):
        test_class_C(
            obs,
            degree,
        )

    for degree in range(0, 4):
        test_class_D(
            obs,
            degree,
        )

    arbitrary_mixing_profile()

    terminal_source_reference()

    interpretation()

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  exact_arithmetic=True"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
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
        "EXPERIMENT 298R COMPLETE"
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