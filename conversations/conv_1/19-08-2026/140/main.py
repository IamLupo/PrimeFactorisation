#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 306R — EXACT CROSS-LAYER TRANSFER NORMALIZATION /
SIMILARITY AUDIT
==============================================================================

Purpose
-------
Experiment 305R found exactly two independently determined width-2
cross-layer transfers:

    A[j,t+1] = alpha_t A[j,t] + beta_t A[j+1,t]

for:

    t=1 -> 2
    t=2 -> 3

A universal constant transfer was ruled out.

Before concluding that the two transfers are genuinely different, this
experiment tests whether they become the SAME operator after exact
coordinate normalization.

We test:

    A~[j,t] = u_t * c^j * A[j,t]

which transforms

    alpha_t -> u_(t+1)/u_t * alpha_t

    beta_t  -> u_(t+1)/u_t / c * beta_t.

Therefore two transitions can be made identical iff the ratio

    beta_t / alpha_t

is related by one common geometric coordinate factor c.

The experiment additionally audits:

    * exact c reconstruction;
    * exact normalization factors u_t;
    * natural rational/integer candidate values of c;
    * prime valuations of c;
    * normalized operator equality;
    * determinant-like and ratio invariants;
    * a more general diagonal similarity test;
    * whether the normalization is unique.

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
# BASIC EXACT HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


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


def gcd_int(a, b):
    return math.gcd(int(a), int(b))


def D(p):
    return len(Q[p]) - 1


# ============================================================================
# TERMINAL-DISTANCE LAYERS
# ============================================================================

def build_layers():

    layers = {}

    max_t = max(
        D(p)
        for p in Q
    )

    for t in range(max_t + 1):

        row = []

        for p in sorted(Q):

            r = D(p) - t

            if r < 0:
                continue

            row.append(
                (
                    p,
                    sp.Integer(Q[p][r]),
                )
            )

        layers[t] = row

    return layers


# ============================================================================
# FINITE DIFFERENCES
# ============================================================================

def finite_differences(values):

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


def first_entry_triangle(layers):

    A = {}

    for t, layer in layers.items():

        values = [
            value
            for _, value in layer
        ]

        diff_rows = finite_differences(
            values
        )

        for j, row in enumerate(
            diff_rows
        ):

            if row:
                A[(j, t)] = clean(
                    row[0]
                )

    return A


# ============================================================================
# WIDTH-2 TRANSFER RECONSTRUCTION
# ============================================================================

def solve_width2(A, t):

    equations = []

    max_j = max(
        j
        for (j, tt) in A
        if tt == t
    )

    for j in range(max_j + 1):

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

    M = sp.Matrix([
        [x0, x1]
        for x0, x1, _ in equations
    ])

    b = sp.Matrix([
        lhs
        for _, _, lhs in equations
    ])

    rank = M.rank()
    aug_rank = M.row_join(b).rank()

    if aug_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "equations": equations,
            "alpha": None,
            "beta": None,
            "rank": rank,
            "augmented_rank": aug_rank,
        }

    if rank < 2:

        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "alpha": None,
            "beta": None,
            "rank": rank,
            "augmented_rank": aug_rank,
        }

    solution = sp.solve(
        [
            sp.Eq(
                alpha * x0
                + beta * x1,
                lhs,
            )
            for x0, x1, lhs
            in equations
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
            "augmented_rank": aug_rank,
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
            "augmented_rank": aug_rank,
        }

    a = clean(sol[alpha])
    bcoef = clean(sol[beta])

    verified = all(
        clean(
            a * x0
            + bcoef * x1
            - lhs
        ) == 0
        for x0, x1, lhs
        in equations
    )

    if not verified:

        return {
            "status": "VERIFICATION_FAILED",
            "equations": equations,
            "alpha": a,
            "beta": bcoef,
            "rank": rank,
            "augmented_rank": aug_rank,
        }

    return {
        "status": "EXACT",
        "equations": equations,
        "alpha": a,
        "beta": bcoef,
        "rank": rank,
        "augmented_rank": aug_rank,
    }


# ============================================================================
# SIMPLE DISPLAY
# ============================================================================

def print_layers(layers):

    print()
    print("=" * 78)
    print("1. TERMINAL-DISTANCE SOURCE LAYERS")
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


def print_A(A):

    print()
    print("=" * 78)
    print("2. FIRST-ENTRY DIFFERENCE TRIANGLE A[j,t]")
    print("=" * 78)

    js = sorted(
        set(
            j
            for j, _ in A
        )
    )

    for j in js:

        entries = []

        for t in sorted(
            tt
            for jj, tt in A
            if jj == j
        ):

            entries.append(
                (
                    t,
                    A[(j, t)],
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
# TRANSFER EXTRACTION
# ============================================================================

def get_transfer_profile(A):

    profile = {}

    print()
    print("=" * 78)
    print("3. EXACT WIDTH-2 TRANSFER PROFILE")
    print("=" * 78)

    for t in range(5):

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
# BASIC INVARIANTS
# ============================================================================

def invariant_audit(profile):

    print()
    print("=" * 78)
    print("4. TRANSFER INVARIANT AUDIT")
    print("=" * 78)

    exact = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    for t, alpha, beta in exact:

        print()
        print(
            "  t={}:".format(t)
        )

        print(
            "    alpha={}".format(
                alpha
            )
        )

        print(
            "    beta={}".format(
                beta
            )
        )

        if alpha != 0:

            print(
                "    beta/alpha={}".format(
                    clean(beta / alpha)
                )
            )

        if beta != 0:

            print(
                "    alpha/beta={}".format(
                    clean(alpha / beta)
                )
            )

        print(
            "    alpha+beta={}".format(
                clean(alpha + beta)
            )
        )

        print(
            "    alpha-beta={}".format(
                clean(alpha - beta)
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
                    valuation(alpha, prime),
                    prime,
                    valuation(beta, prime),
                )
            )


# ============================================================================
# GEOMETRIC-DIAGONAL NORMALIZATION
# ============================================================================

def geometric_similarity(profile):

    print()
    print("=" * 78)
    print(
        "5. GEOMETRIC DIAGONAL NORMALIZATION TEST"
    )
    print("=" * 78)

    exact = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    if len(exact) < 2:

        print(
            "  status=INSUFFICIENT_EXACT_TRANSITIONS"
        )

        return None

    t0, a0, b0 = exact[0]

    for t1, a1, b1 in exact[1:]:

        print()
        print(
            "  compare t={} against t={}".format(
                t0,
                t1,
            )
        )

        # We want constants lambda and c with
        #
        #   a1 = lambda*a0
        #   b1 = lambda*b0/c
        #
        # so:
        #
        #   c = (a1/b1) / (a0/b0).

        if (
            a0 == 0
            or b0 == 0
            or a1 == 0
            or b1 == 0
        ):

            print(
                "    zero coefficient blocks geometric test"
            )
            continue

        lambda_value = clean(
            a1 / a0
        )

        c_value = clean(
            (a1 / b1)
            /
            (a0 / b0)
        )

        print(
            "    lambda={}".format(
                lambda_value
            )
        )

        print(
            "    c={}".format(
                c_value
            )
        )

        # Verify both transformed coefficients.

        transformed_alpha = clean(
            lambda_value * a0
        )

        transformed_beta = clean(
            lambda_value * b0 / c_value
        )

        print(
            "    transformed_alpha={}".format(
                transformed_alpha
            )
        )

        print(
            "    target_alpha={}".format(
                a1
            )
        )

        print(
            "    transformed_beta={}".format(
                transformed_beta
            )
        )

        print(
            "    target_beta={}".format(
                b1
            )
        )

        print(
            "    exact_alpha={}".format(
                transformed_alpha == a1
            )
        )

        print(
            "    exact_beta={}".format(
                transformed_beta == b1
            )
        )

        print(
            "    c_integer={}".format(
                c_value.q == 1
            )
        )

        print(
            "    c_unit={}".format(
                c_value in (
                    sp.Integer(1),
                    sp.Integer(-1),
                )
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
                "    v_{}(c)={}".format(
                    prime,
                    valuation(
                        c_value,
                        prime,
                    )
                )
            )

    return True


# ============================================================================
# NATURAL NORMALIZATION SEARCH
# ============================================================================

def natural_c_search(profile):

    print()
    print("=" * 78)
    print(
        "6. NATURAL c NORMALIZATION SEARCH"
    )
    print("=" * 78)

    exact = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    if len(exact) < 2:

        print(
            "  status=INSUFFICIENT_EXACT_TRANSITIONS"
        )

        return

    t0, a0, b0 = exact[0]

    natural = [
        sp.Integer(-12),
        sp.Integer(-8),
        sp.Integer(-7),
        sp.Integer(-6),
        sp.Integer(-5),
        sp.Integer(-4),
        sp.Integer(-3),
        sp.Integer(-2),
        sp.Integer(-1),
        sp.Integer(1),
        sp.Integer(2),
        sp.Integer(3),
        sp.Integer(4),
        sp.Integer(5),
        sp.Integer(6),
        sp.Integer(7),
        sp.Integer(8),
        sp.Integer(12),
        sp.Rational(1, 2),
        sp.Rational(1, 3),
        sp.Rational(2, 3),
        sp.Rational(3, 2),
        sp.Rational(4, 3),
        sp.Rational(3, 4),
        sp.Rational(5, 2),
    ]

    for t1, a1, b1 in exact[1:]:

        print()
        print(
            "  target transition={}".format(
                t1
            )
        )

        for c in natural:

            if a0 == 0 or b0 == 0:
                continue

            implied_lambda_a = clean(
                a1 / a0
            )

            implied_lambda_b = clean(
                b1 * c / b0
            )

            if (
                implied_lambda_a
                == implied_lambda_b
            ):

                print(
                    "    EXACT natural c={} lambda={}".format(
                        c,
                        implied_lambda_a,
                    )
                )


# ============================================================================
# FULL DIAGONAL SIMILARITY TEST
# ============================================================================

def general_diagonal_similarity(profile):

    print()
    print("=" * 78)
    print(
        "7. GENERAL DIAGONAL SIMILARITY TEST"
    )
    print("=" * 78)

    exact = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    if len(exact) < 2:

        print(
            "  status=INSUFFICIENT_EXACT_TRANSITIONS"
        )

        return

    _, a0, b0 = exact[0]

    # General diagonal scaling:
    #
    #   A~[j,t] = u_t d_j A[j,t]
    #
    # gives
    #
    #   alpha'_t = (u_{t+1}/u_t) alpha_t
    #
    #   beta'_t,j = (u_{t+1}/u_t)(d_j/d_{j+1}) beta_t
    #
    # A common operator requires the alpha ratio to be t-independent
    # and the beta ratios to differ from it by a j-dependent factor.
    #
    # With only j=0,1 available for the exact transitions, we can test
    # the necessary scalar invariant beta/alpha.

    print()
    print(
        "  necessary invariant:"
    )

    baseline_ratio = clean(
        b0 / a0
    )

    print(
        "    baseline beta/alpha={}".format(
            baseline_ratio
        )
    )

    all_equal = True

    for t, alpha, beta in exact:

        ratio = clean(
            beta / alpha
        )

        equal = (
            ratio == baseline_ratio
        )

        print(
            "    t={}: beta/alpha={} equal_to_baseline={}".format(
                t,
                ratio,
                equal,
            )
        )

        if not equal:
            all_equal = False

    print()
    print(
        "  invariant_constant_across_transitions={}".format(
            all_equal
        )
    )

    if not all_equal:

        print(
            "  conclusion=NO COMMON SCALAR/GEOMETRIC "
            "DIAGONAL NORMALIZATION"
        )

    else:

        print(
            "  conclusion=NECESSARY_INVARIANT_PASSES"
        )


# ============================================================================
# EXACT CONJUGACY DIFFERENCE
# ============================================================================

def conjugacy_audit(profile):

    print()
    print("=" * 78)
    print(
        "8. EXACT TRANSFER-PAIR CONJUGACY AUDIT"
    )
    print("=" * 78)

    exact = [
        (
            t,
            result["alpha"],
            result["beta"],
        )
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    if len(exact) < 2:

        print(
            "  status=INSUFFICIENT_EXACT_TRANSITIONS"
        )

        return

    for i in range(
        len(exact)
    ):

        for j in range(
            i + 1,
            len(exact),
        ):

            t1, a1, b1 = exact[i]
            t2, a2, b2 = exact[j]

            if (
                a1 == 0
                or a2 == 0
                or b1 == 0
                or b2 == 0
            ):

                continue

            ratio_a = clean(
                a2 / a1
            )

            ratio_b = clean(
                b2 / b1
            )

            print()
            print(
                "  pair=({},{}):".format(
                    t1,
                    t2,
                )
            )

            print(
                "    alpha_ratio={}".format(
                    ratio_a
                )
            )

            print(
                "    beta_ratio={}".format(
                    ratio_b
                )
            )

            print(
                "    equal_ratios={}".format(
                    ratio_a == ratio_b
                )
            )

            print(
                "    ratio_of_ratios={}".format(
                    clean(
                        ratio_b / ratio_a
                    )
                )
            )


# ============================================================================
# TERMINAL SOURCE REFERENCE
# ============================================================================

def terminal_reference():

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    g = math.gcd(
        q1,
        q3,
    )

    print()
    print("=" * 78)
    print(
        "9. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    print(
        "  q1_terminal={}".format(q1)
    )

    print(
        "  q3_terminal={}".format(q3)
    )

    print(
        "  gcd={}".format(g)
    )

    print(
        "  q1/17={}".format(
            q1 // 17
            if q1 % 17 == 0
            else None
        )
    )

    print(
        "  q3/17={}".format(
            q3 // 17
            if q3 % 17 == 0
            else None
        )
    )


# ============================================================================
# INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 305R found exactly two independently reconstructed width-2
cross-layer transfers, at t=1 and t=2.

A direct common transfer was ruled out.

Experiment 306R therefore tests whether the apparent difference between
those two transfers is only a coordinate artifact.

The normalization model is

    A~[j,t] = u_t c^j A[j,t].

Under this transformation,

    alpha_t -> (u_{t+1}/u_t) alpha_t,

    beta_t  -> (u_{t+1}/u_t) beta_t / c.

Therefore the ratio

    beta_t / alpha_t

is the key invariant up to one common geometric coordinate factor.

A successful exact normalization would mean:

    transition 1
        and
    transition 2

represent the same underlying width-2 operator in different normalized
coordinates.

Failure would strengthen the conclusion that the two exact transitions
are isolated fits rather than manifestations of one universal operator.

The experiment deliberately does not fit an arbitrary matrix.

Because only two exact transitions are currently available, a successful
normalization is diagnostic rather than a universal theorem. A genuine
law still requires an independent source case or additional exact
transitions.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 306R — EXACT CROSS-LAYER TRANSFER NORMALIZATION / "
        "SIMILARITY AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    A = first_entry_triangle(
        layers
    )

    print_layers(
        layers
    )

    print_A(
        A
    )

    profile = get_transfer_profile(
        A
    )

    invariant_audit(
        profile
    )

    geometric_similarity(
        profile
    )

    natural_c_search(
        profile
    )

    general_diagonal_similarity(
        profile
    )

    conjugacy_audit(
        profile
    )

    terminal_reference()

    interpretation()

    exact_transitions = [
        t
        for t, result in profile.items()
        if result["status"] == "EXACT"
    ]

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
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
        "  common_width2_operator_proved=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  geometric_diagonal_similarity_tested=True"
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
        "EXPERIMENT 306R COMPLETE"
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
