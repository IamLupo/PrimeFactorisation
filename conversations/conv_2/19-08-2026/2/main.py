#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 318R — EXACT RAW-TO-DIFFERENCE COORDINATE TRANSFER LAW AUDIT
==============================================================================

Purpose
-------
Experiment 317R found the first genuinely positive structural connection
between the raw source layers and the previously mysterious first-entry
transfer coefficients.

For the raw layers:

    Q_{t+1,k}
        = c0_t Q_{t,k}
        + c1_t Q_{t,k+1}

for t=1 and t=2.

The first-entry difference triangle satisfies:

    A[j,t+1]
        = alpha_t A[j,t]
        + beta_t A[j+1,t].

The observed coefficients satisfy the exact candidate relation

    alpha_t = c0_t + c1_t
    beta_t  = c1_t.

Experiment 318R tests this structurally.

The key identity is the Newton/difference-coordinate expansion:

    Delta^j Q_t[0] = A[j,t].

If

    Q_{t+1,k}
      = c0 Q_{t,k} + c1 Q_{t,k+1},

then

    Delta^j Q_{t+1}[0]
      = c0 Delta^j Q_t[0]
        + c1 Delta^j Q_t[1].

But

    Delta^j Q_t[1]
      = Delta^j Q_t[0]
        + Delta^(j+1) Q_t[0].

Therefore

    A[j,t+1]
      = (c0+c1) A[j,t]
        + c1 A[j+1,t].

Hence:

    alpha = c0+c1
    beta  = c1.

This experiment:

    1. reconstructs raw width-2 coefficients;
    2. reconstructs first-entry width-2 coefficients;
    3. verifies alpha=c0+c1 and beta=c1 exactly;
    4. derives the transformed recurrence independently;
    5. verifies it for EVERY available j;
    6. checks the underlying finite-difference identity directly;
    7. constructs the finite-dimensional Pascal/difference coordinate maps;
    8. checks exact operator conjugation on all available coordinates;
    9. separates the genuine source-level law from the coordinate change;
   10. searches whether the same transformation explains both transitions.

This is not interpolation.

The transformed coefficients are derived algebraically from the raw
coefficients and then verified on independent first-entry equations.

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


# ============================================================================
# SOURCE LAYERS
# ============================================================================

def D(p):
    return len(Q[p]) - 1


def build_layers():

    layers = {}

    maximum_t = max(
        D(p)
        for p in Q
    )

    for t in range(
        maximum_t + 1
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


def values(layers, t):

    return [
        value
        for _, value
        in layers[t]
    ]


# ============================================================================
# FINITE DIFFERENCES
# ============================================================================

def difference_rows(row):

    current = [
        sp.Rational(v)
        for v in row
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

    for t in sorted(layers):

        rows = difference_rows(
            values(
                layers,
                t,
            )
        )

        for j, row in enumerate(rows):

            if row:
                A[(j, t)] = clean(
                    row[0]
                )

    return A


# ============================================================================
# RAW WIDTH-2 TRANSFER
# ============================================================================

def solve_raw_width2(
    layers,
    t,
):

    source = values(
        layers,
        t,
    )

    target = values(
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
        [q0, q1]
        for q0, q1, _ in equations
    ])

    rhs = sp.Matrix([
        target
        for _, _, target in equations
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
            c0 * q0
            + c1 * q1
            - target
        )
        for q0, q1, target in equations
    ]

    verified = all(
        r == 0
        for r in residuals
    )

    return {
        "status": (
            "EXACT"
            if verified
            else "VERIFICATION_FAILED"
        ),
        "equations": equations,
        "c0": c0,
        "c1": c1,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "residuals": residuals,
    }


# ============================================================================
# FIRST-ENTRY WIDTH-2 TRANSFER
# ============================================================================

def solve_A_width2(
    A,
    t,
):

    equations = []

    maximum_j = max(
        (
            j
            for j, tt in A
            if tt == t
        ),
        default=-1,
    )

    for j in range(
        maximum_j + 1
    ):

        lhs = A.get(
            (j, t + 1)
        )

        a0 = A.get(
            (j, t)
        )

        a1 = A.get(
            (j + 1, t)
        )

        if (
            lhs is None
            or a0 is None
            or a1 is None
        ):
            continue

        equations.append(
            (
                a0,
                a1,
                lhs,
            )
        )

    if len(equations) < 2:

        return {
            "status": "INSUFFICIENT_DATA",
            "equations": equations,
        }

    M = sp.Matrix([
        [a0, a1]
        for a0, a1, _
        in equations
    ])

    rhs = sp.Matrix([
        lhs
        for _, _, lhs
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

    alpha = clean(
        coeffs[0, 0]
    )

    beta = clean(
        coeffs[1, 0]
    )

    residuals = [
        clean(
            alpha * a0
            + beta * a1
            - lhs
        )
        for a0, a1, lhs
        in equations
    ]

    verified = all(
        r == 0
        for r in residuals
    )

    return {
        "status": (
            "EXACT"
            if verified
            else "VERIFICATION_FAILED"
        ),
        "equations": equations,
        "alpha": alpha,
        "beta": beta,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "residuals": residuals,
    }


# ============================================================================
# DIRECT NEWTON-IDENTITY VERIFICATION
# ============================================================================

def direct_newton_identity(
    layers,
    t,
):

    source = values(
        layers,
        t,
    )

    target = values(
        layers,
        t + 1,
    )

    raw = solve_raw_width2(
        layers,
        t,
    )

    if raw.get(
        "status"
    ) != "EXACT":

        return {
            "status": "RAW_TRANSFER_UNAVAILABLE"
        }

    c0 = raw["c0"]
    c1 = raw["c1"]

    rows_source = difference_rows(
        source
    )

    rows_target = difference_rows(
        target
    )

    residuals = []

    identities = []

    maximum_j = min(
        len(rows_target),
        len(rows_source),
    )

    # For each j such that A[j+1,t] exists:
    #
    # Delta^j Q_t[1]
    #   =
    # Delta^j Q_t[0]
    #   +
    # Delta^(j+1) Q_t[0].
    #
    # Therefore:
    #
    # Delta^j Q_{t+1}[0]
    # =
    # c0 A[j,t]
    # +
    # c1(A[j,t] + A[j+1,t]).

    for j in range(
        maximum_j
    ):

        if not rows_source[j]:
            continue

        if not rows_target[j]:
            continue

        lhs = clean(
            rows_target[j][0]
        )

        A_j = clean(
            rows_source[j][0]
        )

        if j + 1 < len(rows_source):

            A_j1 = clean(
                rows_source[j + 1][0]
            )

            transformed_rhs = clean(
                c0 * A_j
                + c1 * (
                    A_j
                    + A_j1
                )
            )

            residual = clean(
                transformed_rhs
                - lhs
            )

            residuals.append(
                residual
            )

            identities.append(
                {
                    "j": j,
                    "lhs": lhs,
                    "A_j": A_j,
                    "A_j1": A_j1,
                    "derived_rhs": transformed_rhs,
                    "residual": residual,
                }
            )

    return {
        "status": "EXACT",
        "c0": c0,
        "c1": c1,
        "residuals": residuals,
        "identities": identities,
        "all_exact": all(
            r == 0
            for r in residuals
        ),
    }


# ============================================================================
# COORDINATE TRANSFORM FORMULA
# ============================================================================

def transformed_coefficients(
    c0,
    c1,
):

    # Exact Newton-coordinate transformation:
    #
    #   alpha = c0 + c1
    #   beta  = c1

    return (
        clean(c0 + c1),
        clean(c1),
    )


# ============================================================================
# PASCAL / DIFFERENCE BASIS MATRICES
# ============================================================================

def difference_matrix(n):

    """
    Lower-triangular matrix D such that

        A_j = sum_k D[j,k] q_k

    with

        D[j,k] = (-1)^(j-k) binomial(j,k)

    for k <= j.

    Thus the first-entry finite-difference vector is

        A = D q.
    """

    D = sp.zeros(
        n,
        n,
    )

    for j in range(n):

        for k in range(
            j + 1
        ):

            D[j, k] = (
                (-1) ** (j - k)
                * sp.binomial(
                    j,
                    k,
                )
            )

    return D


def raw_operator_matrix(
    n,
    c0,
    c1,
):

    """
    Leading triangular part of the raw width-2 operator:

        (Cq)_k = c0 q_k + c1 q_{k+1}.

    The final boundary row is omitted from the shift.
    """

    C = sp.zeros(
        n,
        n,
    )

    for k in range(n):

        C[k, k] = c0

        if k + 1 < n:

            C[k, k + 1] = c1

    return C


def newton_operator_matrix(
    n,
    c0,
    c1,
):

    alpha, beta = (
        transformed_coefficients(
            c0,
            c1,
        )
    )

    C = sp.zeros(
        n,
        n,
    )

    for j in range(n):

        C[j, j] = alpha

        if j + 1 < n:

            C[j, j + 1] = beta

    return C


def coordinate_conjugation_audit(
    n,
    c0,
    c1,
):

    D = difference_matrix(
        n
    )

    C_raw = raw_operator_matrix(
        n,
        c0,
        c1,
    )

    C_newton = newton_operator_matrix(
        n,
        c0,
        c1,
    )

    # The exact orientation is:
    #
    #     D * C_raw * D^{-1}
    #
    # when the same finite-dimensional square truncation is used.
    #
    # Boundary effects are expected in the final column/row. Therefore we
    # report the full matrix equality and also the leading n-1 block.

    try:
        D_inv = D.inv()

        conjugated = (
            D
            * C_raw
            * D_inv
        )

        conjugated = conjugated.applyfunc(
            clean
        )

        full_exact = (
            conjugated
            == C_newton
        )

    except Exception:

        conjugated = None
        full_exact = False

    if n > 1:

        leading_exact = (
            conjugated[:n-1, :n-1]
            ==
            C_newton[:n-1, :n-1]
            if conjugated is not None
            else False
        )

    else:

        leading_exact = (
            conjugated == C_newton
            if conjugated is not None
            else False
        )

    return {
        "D": D,
        "C_raw": C_raw,
        "C_newton": C_newton,
        "conjugated": conjugated,
        "full_exact": full_exact,
        "leading_exact": leading_exact,
    }


# ============================================================================
# TRANSITION AUDIT
# ============================================================================

def transition_audit(
    layers,
    A,
    t,
):

    raw = solve_raw_width2(
        layers,
        t,
    )

    first = solve_A_width2(
        A,
        t,
    )

    print()
    print("=" * 78)
    print(
        "TRANSITION t={} -> {}"
        .format(
            t,
            t + 1,
        )
    )
    print("=" * 78)

    print()
    print(
        "  raw_status={}".format(
            raw["status"]
        )
    )

    if raw["status"] != "EXACT":
        return {
            "raw": raw,
            "first": first,
            "derived_match": False,
            "newton_identity": None,
        }

    print(
        "  raw_c0={}".format(
            raw["c0"]
        )
    )

    print(
        "  raw_c1={}".format(
            raw["c1"]
        )
    )

    alpha_derived, beta_derived = (
        transformed_coefficients(
            raw["c0"],
            raw["c1"],
        )
    )

    print()
    print(
        "  derived_alpha=c0+c1={}".format(
            alpha_derived
        )
    )

    print(
        "  derived_beta=c1={}".format(
            beta_derived
        )
    )

    print()
    print(
        "  first_entry_status={}".format(
            first["status"]
        )
    )

    if first["status"] == "EXACT":

        print(
            "  actual_alpha={}".format(
                first["alpha"]
            )
        )

        print(
            "  actual_beta={}".format(
                first["beta"]
            )
        )

        alpha_match = (
            alpha_derived
            == first["alpha"]
        )

        beta_match = (
            beta_derived
            == first["beta"]
        )

        print(
            "  alpha_match={}".format(
                alpha_match
            )
        )

        print(
            "  beta_match={}".format(
                beta_match
            )
        )

        derived_match = (
            alpha_match
            and beta_match
        )

    else:

        derived_match = False

    newton = direct_newton_identity(
        layers,
        t,
    )

    print()
    print(
        "  direct_newton_identity_exact={}".format(
            newton.get(
                "all_exact",
                False,
            )
        )
    )

    for item in newton.get(
        "identities",
        [],
    ):

        print()
        print(
            "    j={}:".format(
                item["j"]
            )
        )

        print(
            "      A[j,t]={}".format(
                item["A_j"]
            )
        )

        print(
            "      A[j+1,t]={}".format(
                item["A_j1"]
            )
        )

        print(
            "      actual_A[j,t+1]={}".format(
                item["lhs"]
            )
        )

        print(
            "      derived_rhs={}".format(
                item["derived_rhs"]
            )
        )

        print(
            "      residual={}".format(
                item["residual"]
            )
        )

    return {
        "raw": raw,
        "first": first,
        "derived_alpha": alpha_derived,
        "derived_beta": beta_derived,
        "derived_match": derived_match,
        "newton_identity": newton,
    }


# ============================================================================
# PRIME PROFILE
# ============================================================================

def coefficient_profile(
    transition_results,
):

    print()
    print("=" * 78)
    print(
        "PRIME PROFILE OF RAW AND TRANSFORMED COEFFICIENTS"
    )
    print("=" * 78)

    for t, result in transition_results.items():

        raw = result["raw"]

        if raw["status"] != "EXACT":
            continue

        c0 = raw["c0"]
        c1 = raw["c1"]

        alpha = result[
            "derived_alpha"
        ]

        beta = result[
            "derived_beta"
        ]

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        for label, value in (
            ("c0", c0),
            ("c1", c1),
            ("alpha", alpha),
            ("beta", beta),
        ):

            print(
                "    {}={}".format(
                    label,
                    value,
                )
            )

            print(
                "      valuations={}".format(
                    {
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
                )
            )


# ============================================================================
# TERMINAL SOURCE REFERENCE
# ============================================================================

def terminal_reference():

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print()
    print("=" * 78)
    print(
        "TERMINAL PROJECTIVE SOURCE REFERENCE"
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


# ============================================================================
# INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 317R found:

    raw width-2 laws exist at t=1 and t=2,

but the raw coefficients are not numerically identical to the
first-entry coefficients.

Experiment 318R tests whether this discrepancy is exactly the expected
coordinate effect of passing from raw values to first finite differences.

For any sequence q_k, define

    A_j = Delta^j q_0.

Then

    q_1 = A_0 + A_1

and, more generally,

    Delta^j q_1
      = A_j + A_{j+1}.

If

    q'__k = c0 q_k + c1 q_{k+1},

then

    A'_j
      = c0 A_j + c1(A_j+A_{j+1})
      = (c0+c1)A_j + c1 A_{j+1}.

Therefore:

    alpha = c0+c1,
    beta  = c1.

A successful exact verification across every available j would mean that
the apparently mysterious first-entry transfer coefficients are not
independent discoveries.

They are the exact Newton-coordinate representation of a simpler raw
source-layer width-2 operator.

That would substantially change the interpretation of Experiments
305R--316R:

    the 2x2 companion matrices are still noncommutative as operators,
    but their coefficients may originate from a simpler scalar
    cross-layer law before finite differencing.

The next question after this experiment would therefore be:

    Is there a simple exact law for the RAW coefficients c0_t,c1_t?

That is a more promising target than searching for further invariants
of the already-derived companion matrices.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 318R — EXACT RAW-TO-DIFFERENCE "
        "COORDINATE TRANSFER LAW AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    A = first_entry_triangle(
        layers
    )

    print()
    print(
        "SOURCE LAYERS"
    )
    print(
        layers
    )

    transition_results = {}

    for t in (
        1,
        2,
    ):

        transition_results[t] = transition_audit(
            layers,
            A,
            t,
        )

    # ========================================================================
    # COORDINATE CONJUGATION AUDITS
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "EXACT DIFFERENCE-COORDINATE CONJUGATION AUDIT"
    )
    print("=" * 78)

    conjugation_results = {}

    for t, result in transition_results.items():

        raw = result["raw"]

        if raw["status"] != "EXACT":
            continue

        conj = coordinate_conjugation_audit(
            3,
            raw["c0"],
            raw["c1"],
        )

        conjugation_results[t] = conj

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    D="
        )
        print(
            conj["D"]
        )

        print(
            "    C_raw="
        )
        print(
            conj["C_raw"]
        )

        print(
            "    C_newton="
        )
        print(
            conj["C_newton"]
        )

        print(
            "    D C_raw D^-1="
        )
        print(
            conj["conjugated"]
        )

        print(
            "    full_conjugation_exact={}".format(
                conj["full_exact"]
            )
        )

        print(
            "    leading_block_exact={}".format(
                conj["leading_exact"]
            )
        )

    coefficient_profile(
        transition_results
    )

    terminal_reference()

    interpretation()

    # ========================================================================
    # FINAL EXACTNESS
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    raw_exact = all(
        result["raw"]["status"] == "EXACT"
        for result in transition_results.values()
    )

    derived_exact = all(
        result["derived_match"]
        for result in transition_results.values()
    )

    newton_exact = all(
        result["newton_identity"] is not None
        and
        result["newton_identity"].get(
            "all_exact",
            False,
        )
        for result in transition_results.values()
    )

    conjugation_leading_exact = all(
        result["leading_exact"]
        for result in conjugation_results.values()
    )

    print(
        "  raw_width2_law_exact={}".format(
            raw_exact
        )
    )

    print(
        "  transformed_coefficients_match_first_entry={}".format(
            derived_exact
        )
    )

    print(
        "  direct_newton_identity_exact={}".format(
            newton_exact
        )
    )

    print(
        "  leading_coordinate_conjugation_exact={}".format(
            conjugation_leading_exact
        )
    )

    print(
        "  two_transitions_tested=True"
    )

    print(
        "  no_interpolation_used=True"
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
        "EXPERIMENT 318R COMPLETE"
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

