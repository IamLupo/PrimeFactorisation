#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 305R — EXACT CROSS-LAYER TRANSFER-COEFFICIENT LAW AUDIT
==============================================================================

Purpose
-------
Experiment 304RR found two exact width-2 cross-layer transitions:

    A[j,t+1] = alpha_t A[j,t] + beta_t A[j+1,t]

for t=1 and t=2, while t=0 had no unique solution.

This experiment asks the next structural question:

    Are the surviving alpha_t, beta_t coefficients themselves governed by
    a simple exact law in t?

Tests
-----
1. Reconstruct the terminal-distance layers Q_t.
2. Reconstruct the finite-difference triangle.
3. Extract the first-entry triangle A[j,t].
4. Solve width-2 transfer coefficients independently for every transition.
5. Test:
       alpha_t, beta_t = constants
       alpha_t, beta_t = affine in t
       alpha_t, beta_t = quadratic in t
       alpha_t / beta_t = simple rational/polynomial relation
6. Perform leave-one-transition-out prediction whenever enough transitions
   are available.
7. Test a common transfer law directly on all overdetermined equations.
8. Keep insufficient-data cases separate from exact failures.

No external files.
No synthetic second n=pq case.
Exact rational arithmetic throughout.
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

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def finite_differences(values):
    current = [
        sp.Rational(v)
        for v in values
    ]

    rows = [current]

    while len(current) > 1:
        current = [
            clean(
                current[i + 1] - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]
        rows.append(current)

    return rows


def D(p):
    return len(Q[p]) - 1


def valuation(value, prime):
    value = sp.Rational(value)

    if value == 0:
        return None

    num = abs(int(value.p))
    den = abs(int(value.q))

    v = 0

    while num % prime == 0:
        num //= prime
        v += 1

    while den % prime == 0:
        den //= prime
        v -= 1

    return v


def primitive_integer_list(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    lcm_den = 1

    for value in values:
        lcm_den = math.lcm(
            lcm_den,
            int(value.q),
        )

    integers = [
        int(value * lcm_den)
        for value in values
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    if g:
        integers = [
            value // g
            for value in integers
        ]

    for value in integers:
        if value != 0:
            if value < 0:
                integers = [
                    -x
                    for x in integers
                ]
            break

    return integers


# ============================================================================
# TERMINAL-DISTANCE LAYERS
# ============================================================================

def build_layers():

    layers = {}

    max_t = max(
        D(p)
        for p in Q
    )

    for t in range(
        max_t + 1
    ):

        row = []

        for p in sorted(Q):

            r = D(p) - t

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


# ============================================================================
# DIFFERENCE TRIANGLE
# ============================================================================

def build_difference_triangle(layers):

    triangle = {}

    for t in sorted(layers):

        values = [
            value
            for _, value in layers[t]
        ]

        differences = finite_differences(
            values
        )

        for j, row in enumerate(
            differences
        ):

            triangle[
                (j, t)
            ] = row

    return triangle


def build_A(triangle):

    A = {}

    for (j, t), row in triangle.items():

        if row:
            A[
                (j, t)
            ] = clean(
                row[0]
            )

    return A


# ============================================================================
# WIDTH-2 TRANSFER SOLVER
# ============================================================================

def solve_width2(A, t):

    equations = []

    max_j = max(
        j
        for j, tt in A
        if tt == t
    )

    for j in range(
        max_j + 1
    ):

        lhs = A.get(
            (j, t + 1)
        )

        x0 = A.get(
            (j, t)
        )

        x1 = A.get(
            (j + 1, t)
        )

        if (
            lhs is None
            or x0 is None
            or x1 is None
        ):
            continue

        equations.append(
            (
                sp.Rational(x0),
                sp.Rational(x1),
                sp.Rational(lhs),
            )
        )

    if len(equations) < 2:
        return {
            "status": "INSUFFICIENT_DATA",
            "equations": equations,
            "alpha": None,
            "beta": None,
        }

    alpha, beta = sp.symbols(
        "alpha beta"
    )

    matrix = sp.Matrix([
        [x0, x1]
        for x0, x1, _ in equations
    ])

    rhs = sp.Matrix([
        lhs
        for _, _, lhs in equations
    ])

    rank = matrix.rank()
    augmented_rank = matrix.row_join(
        rhs
    ).rank()

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "equations": equations,
            "alpha": None,
            "beta": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    if rank < 2:

        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "alpha": None,
            "beta": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    solution = sp.solve(
        [
            sp.Eq(
                alpha * x0
                + beta * x1,
                lhs,
            )
            for x0, x1, lhs in equations
        ],
        [alpha, beta],
        dict=True,
    )

    if len(solution) != 1:

        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "alpha": None,
            "beta": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    sol = solution[0]

    if (
        alpha not in sol
        or beta not in sol
    ):

        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "alpha": None,
            "beta": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    a = clean(sol[alpha])
    b = clean(sol[beta])

    exact = all(
        clean(
            a * x0
            + b * x1
            - lhs
        ) == 0
        for x0, x1, lhs in equations
    )

    if not exact:

        return {
            "status": "VERIFICATION_FAILED",
            "equations": equations,
            "alpha": a,
            "beta": b,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    return {
        "status": "EXACT",
        "equations": equations,
        "alpha": a,
        "beta": b,
        "rank": rank,
        "augmented_rank": augmented_rank,
    }


# ============================================================================
# DISPLAY SOURCE
# ============================================================================

def print_source_layers(layers):

    print()
    print("=" * 78)
    print(
        "1. TERMINAL-DISTANCE SOURCE LAYERS"
    )
    print("=" * 78)

    for t in sorted(layers):

        print()
        print(
            "  t={}:".format(t)
        )

        print(
            "    available_p={}".format(
                [
                    p
                    for p, _ in layers[t]
                ]
            )
        )

        print(
            "    Q_t={}".format(
                [
                    value
                    for _, value
                    in layers[t]
                ]
            )
        )


# ============================================================================
# DIFFERENCE TRIANGLE DISPLAY
# ============================================================================

def print_difference_triangle(
    triangle
):

    print()
    print("=" * 78)
    print(
        "2. FIRST-ENTRY DIFFERENCE TRIANGLE"
    )
    print("=" * 78)

    for j in sorted(
        set(
            jj
            for jj, _ in triangle
        )
    ):

        entries = []

        for t in sorted(
            tt
            for jj, tt in triangle
            if jj == j
        ):

            entries.append(
                (
                    t,
                    triangle[
                        (j, t)
                    ][0],
                )
            )

        print()
        print(
            "  j={}: {}".format(
                j,
                entries,
            )
        )


# ============================================================================
# TRANSFER COEFFICIENT PROFILE
# ============================================================================

def transfer_profile(A):

    print()
    print("=" * 78)
    print(
        "3. EXACT WIDTH-2 TRANSFER COEFFICIENTS"
    )
    print("=" * 78)

    profile = {}

    max_t = max(
        tt
        for _, tt in A
    )

    for t in range(max_t):

        result = solve_width2(
            A,
            t,
        )

        profile[t] = result

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
                "    alpha={}".format(
                    result["alpha"]
                )
            )

            print(
                "    beta={}".format(
                    result["beta"]
                )
            )

            print(
                "    equations={}".format(
                    len(
                        result["equations"]
                    )
                )
            )

    return profile


# ============================================================================
# COMMON-CONSTANT TRANSFER TEST
# ============================================================================

def common_constant_transfer(A):

    print()
    print("=" * 78)
    print(
        "4. COMMON CONSTANT WIDTH-2 TRANSFER"
    )
    print("=" * 78)

    alpha, beta = sp.symbols(
        "alpha beta"
    )

    equations = []

    for (j, t), lhs in A.items():

        target = A.get(
            (j, t + 1)
        )

        right0 = A.get(
            (j, t)
        )

        right1 = A.get(
            (j + 1, t)
        )

        if (
            target is None
            or right0 is None
            or right1 is None
        ):
            continue

        equations.append(
            sp.Eq(
                alpha * right0
                + beta * right1,
                target,
            )
        )

    if len(equations) < 2:

        print(
            "  status=INSUFFICIENT_DATA"
        )

        return None

    matrix = []
    rhs = []

    for eq in equations:

        expr = sp.expand(
            eq.lhs - eq.rhs
        )

        coeff_a = expr.coeff(
            alpha
        )

        coeff_b = expr.coeff(
            beta
        )

        constant = -(
            expr.subs(
                {
                    alpha: 0,
                    beta: 0,
                }
            )
        )

        matrix.append(
            [
                coeff_a,
                coeff_b,
            ]
        )

        rhs.append(
            constant
        )

    M = sp.Matrix(matrix)
    b = sp.Matrix(rhs)

    rank = M.rank()
    augmented_rank = M.row_join(
        b
    ).rank()

    print(
        "  equations={}".format(
            len(equations)
        )
    )

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  augmented_rank={}".format(
            augmented_rank
        )
    )

    if augmented_rank > rank:

        print(
            "  status=NO_SOLUTION"
        )

        return None

    if rank < 2:

        print(
            "  status=NONUNIQUE"
        )

        return None

    solution = sp.solve(
        equations,
        [alpha, beta],
        dict=True,
    )

    if len(solution) != 1:

        print(
            "  status=NONUNIQUE"
        )

        return None

    sol = solution[0]

    if (
        alpha not in sol
        or beta not in sol
    ):

        print(
            "  status=NONUNIQUE"
        )

        return None

    a = clean(sol[alpha])
    b = clean(sol[beta])

    verified = all(
        clean(
            a * A[(j, t)]
            + b * A[(j + 1, t)]
            - A[(j, t + 1)]
        ) == 0
        for (j, t) in A
        if (
            (j, t + 1) in A
            and
            (j + 1, t) in A
        )
    )

    if verified:

        print(
            "  status=EXACT"
        )

        print(
            "  alpha={}".format(
                a
            )
        )

        print(
            "  beta={}".format(
                b
            )
        )

        return (
            a,
            b,
        )

    print(
        "  status=VERIFICATION_FAILED"
    )

    return None


# ============================================================================
# PARAMETER LAW FIT
# ============================================================================

def exact_polynomial_fit(
    points,
    degree,
):

    if len(points) < degree + 1:
        return None

    tvar = sp.Symbol(
        "t"
    )

    polynomial = sp.interpolate(
        [
            (
                sp.Integer(t),
                sp.Rational(value),
            )
            for t, value in points
        ],
        tvar,
    )

    polynomial = clean(
        polynomial
    )

    actual_degree = sp.Poly(
        polynomial,
        tvar,
        domain=sp.QQ,
    ).degree()

    if actual_degree > degree:
        return None

    return polynomial


def parameter_law_audit(profile):

    print()
    print("=" * 78)
    print(
        "5. EXACT alpha_t / beta_t LAW SEARCH"
    )
    print("=" * 78)

    exact_points = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    print()
    print(
        "  exact_transition_points={}".format(
            [
                t
                for t, _, _
                in exact_points
            ]
        )
    )

    if not exact_points:

        print(
            "  no exact width-2 transitions available"
        )

        return {
            "alpha_hits": [],
            "beta_hits": [],
        }

    alpha_points = [
        (
            t,
            alpha,
        )
        for t, alpha, _
        in exact_points
    ]

    beta_points = [
        (
            t,
            beta,
        )
        for t, _, beta
        in exact_points
    ]

    alpha_hits = []
    beta_hits = []

    for degree in (
        0,
        1,
        2,
    ):

        alpha_law = exact_polynomial_fit(
            alpha_points,
            degree,
        )

        beta_law = exact_polynomial_fit(
            beta_points,
            degree,
        )

        print()
        print(
            "  degree<={}:".format(
                degree
            )
        )

        print(
            "    alpha_law={}".format(
                alpha_law
            )
        )

        print(
            "    beta_law={}".format(
                beta_law
            )
        )

        if alpha_law is not None:
            alpha_hits.append(
                (
                    degree,
                    alpha_law,
                )
            )

        if beta_law is not None:
            beta_hits.append(
                (
                    degree,
                    beta_law,
                )
            )

    return {
        "alpha_hits": alpha_hits,
        "beta_hits": beta_hits,
    }


# ============================================================================
# LEAVE-ONE-TRANSITION-OUT
# ============================================================================

def leave_one_transition_out(
    profile
):

    print()
    print("=" * 78)
    print(
        "6. LEAVE-ONE-TRANSITION-OUT PARAMETER PREDICTION"
    )
    print("=" * 78)

    points = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    if len(points) < 3:

        print(
            "  status=INSUFFICIENT_EXACT_TRANSITIONS"
        )

        return 0, 0

    tested = 0
    exact = 0

    tvar = sp.Symbol(
        "t"
    )

    for omitted_index in range(
        len(points)
    ):

        training = [
            point
            for i, point
            in enumerate(points)
            if i != omitted_index
        ]

        omitted_t, actual_a, actual_b = (
            points[
                omitted_index
            ]
        )

        for degree in (
            0,
            1,
        ):

            if len(training) < degree + 1:
                continue

            a_training = [
                (
                    t,
                    alpha,
                )
                for t, alpha, _
                in training
            ]

            b_training = [
                (
                    t,
                    beta,
                )
                for t, _, beta
                in training
            ]

            a_law = exact_polynomial_fit(
                a_training,
                degree,
            )

            b_law = exact_polynomial_fit(
                b_training,
                degree,
            )

            if (
                a_law is None
                or b_law is None
            ):
                continue

            predicted_a = clean(
                a_law.subs(
                    tvar,
                    omitted_t,
                )
            )

            predicted_b = clean(
                b_law.subs(
                    tvar,
                    omitted_t,
                )
            )

            good = (
                predicted_a
                == actual_a
                and
                predicted_b
                == actual_b
            )

            tested += 1

            if good:
                exact += 1

            print(
                "  omitted_t={} degree<={}: "
                "predicted_alpha={} actual_alpha={} "
                "predicted_beta={} actual_beta={} "
                "exact={}".format(
                    omitted_t,
                    degree,
                    predicted_a,
                    actual_a,
                    predicted_b,
                    actual_b,
                    good,
                )
            )

    print()
    print(
        "  tested={}".format(
            tested
        )
    )

    print(
        "  exact={}".format(
            exact
        )
    )

    return (
        tested,
        exact,
    )


# ============================================================================
# RATIO AUDIT
# ============================================================================

def ratio_audit(profile):

    print()
    print("=" * 78)
    print(
        "7. EXACT alpha/beta AND CROSS-TRANSITION RATIO AUDIT"
    )
    print("=" * 78)

    exact_points = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    if not exact_points:

        print(
            "  status=NO_EXACT_TRANSITIONS"
        )

        return

    for t, alpha, beta in exact_points:

        print()
        print(
            "  t={}:".format(t)
        )

        if beta != 0:

            print(
                "    alpha/beta={}".format(
                    clean(
                        alpha / beta
                    )
                )
            )

        else:

            print(
                "    alpha/beta=UNDEFINED"
            )

        print(
            "    alpha_sign={}".format(
                "+" if alpha > 0 else "-"
            )
        )

        print(
            "    beta_sign={}".format(
                "+" if beta > 0 else "-"
            )
        )

        for prime in (
            2,
            3,
            5,
            7,
            17,
        ):

            print(
                "    v_{}(alpha)={} "
                "v_{}(beta)={}".format(
                    prime,
                    valuation(
                        alpha,
                        prime,
                    ),
                    prime,
                    valuation(
                        beta,
                        prime,
                    ),
                )
            )


# ============================================================================
# FINAL SOURCE REFERENCE
# ============================================================================

def print_terminal_reference():

    q1_terminal = Q[1][-1]
    q3_terminal = Q[3][-1]

    gcd_value = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    print()
    print("=" * 78)
    print(
        "8. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    print(
        "  q1_terminal={}".format(
            q1_terminal
        )
    )

    print(
        "  q3_terminal={}".format(
            q3_terminal
        )
    )

    print(
        "  gcd={}".format(
            gcd_value
        )
    )

    print(
        "  q1/17={}".format(
            q1_terminal // 17
            if q1_terminal % 17 == 0
            else None
        )
    )

    print(
        "  q3/17={}".format(
            q3_terminal // 17
            if q3_terminal % 17 == 0
            else None
        )
    )


# ============================================================================
# STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 304RR established that the terminal-distance coordinate is
compatible with exact mixed finite differences and that the first-entry
difference triangle has full rank 4.

It also exposed two exact width-2 transitions, but those transitions
were fitted independently.

Experiment 305R therefore asks the next identifiability question:

    Are the transfer coefficients alpha_t and beta_t themselves
    controlled by a simple law?

The hierarchy is:

    common constants
        ->
    affine functions of t
        ->
    quadratic functions of t
        ->
    prediction on an omitted transition.

A polynomial passing through every surviving transition is only a fit.

The decisive evidence is leave-one-transition-out prediction.

Because the present data contain only a small number of transitions
with enough rows to determine both coefficients, the experiment
explicitly reports insufficient-data cases rather than converting them
into failures or successes.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 305R — EXACT CROSS-LAYER "
        "TRANSFER-COEFFICIENT LAW AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    triangle = build_difference_triangle(
        layers
    )

    A = build_A(
        triangle
    )

    print_source_layers(
        layers
    )

    print_difference_triangle(
        triangle
    )

    profile = transfer_profile(
        A
    )

    common_constant = (
        common_constant_transfer(
            A
        )
    )

    law_profile = parameter_law_audit(
        profile
    )

    (
        loo_tested,
        loo_exact,
    ) = leave_one_transition_out(
        profile
    )

    ratio_audit(
        profile
    )

    print_terminal_reference()

    interpretation()

    exact_transitions = [
        t
        for t, result
        in profile.items()
        if result["status"] == "EXACT"
    ]

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  terminal_layer_count={}".format(
            len(layers)
        )
    )

    print(
        "  exact_width2_transitions={}".format(
            len(exact_transitions)
        )
    )

    print(
        "  exact_transition_indices={}".format(
            exact_transitions
        )
    )

    print(
        "  common_constant_width2_exact={}".format(
            common_constant is not None
        )
    )

    print(
        "  alpha_polynomial_law_candidates={}".format(
            len(
                law_profile[
                    "alpha_hits"
                ]
            )
        )
    )

    print(
        "  beta_polynomial_law_candidates={}".format(
            len(
                law_profile[
                    "beta_hits"
                ]
            )
        )
    )

    print(
        "  leave_one_transition_tested={}".format(
            loo_tested
        )
    )

    print(
        "  leave_one_transition_exact={}".format(
            loo_exact
        )
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
        "EXPERIMENT 305R COMPLETE"
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

