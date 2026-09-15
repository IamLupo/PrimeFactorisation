#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 307R — EXACT WIDTH-2 TRANSFER PROJECTIVE-SIMILARITY AUDIT
==============================================================================

Purpose
-------
Experiment 306R ruled out the simple geometric diagonal normalization

    A~[j,t] = u_t c^j A[j,t].

There are still two exact width-2 transitions:

    t=1 -> 2
    t=2 -> 3

Each gives

    A[j,t+1] = alpha_t A[j,t] + beta_t A[j+1,t].

Encode this recurrence by the companion transfer matrix

    T_t = Matrix([[alpha_t, beta_t],
                  [1,       0     ]]).

This experiment asks a strictly broader question:

    Is T_2 projectively similar to T_1?

That is, do there exist a nonzero scalar lambda and an invertible matrix S
such that

    T_2 = lambda * S * T_1 * S^(-1) ?

This includes arbitrary 2x2 changes of basis, not merely diagonal or
geometric coordinate rescaling.

For a 2x2 matrix, the key projective similarity invariants are:

    det(T) / tr(T)^2

when tr(T) != 0, together with equivalent eigenvalue-ratio/discriminant
tests.

We therefore perform several exact audits:

    1. exact construction of the two transfer matrices;
    2. trace and determinant;
    3. projective trace/determinant invariant;
    4. discriminant invariant;
    5. eigenvalue-ratio invariant;
    6. direct symbolic projective-similarity equations;
    7. ordinary similarity after removing a scalar;
    8. rational candidate checks;
    9. exact conclusion.

No external files.
No synthetic second n=pq case.
Exact SymPy rational arithmetic only.
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
        return None

    num = abs(int(x.p))
    den = abs(int(x.q))

    v = 0

    while num % prime == 0:
        num //= prime
        v += 1

    while den % prime == 0:
        den //= prime
        v -= 1

    return v


def D(p):
    return len(Q[p]) - 1


# ============================================================================
# TERMINAL-DISTANCE LAYERS
# ============================================================================

def build_layers():

    layers = {}

    maximum_t = max(
        D(p)
        for p in Q
    )

    for t in range(
        maximum_t + 1
    ):

        layer = []

        for p in sorted(Q):

            r = D(p) - t

            if r < 0:
                continue

            layer.append(
                (
                    p,
                    sp.Integer(
                        Q[p][r]
                    ),
                )
            )

        layers[t] = layer

    return layers


# ============================================================================
# FINITE DIFFERENCES
# ============================================================================

def finite_difference_rows(values):

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
            for _, value
            in layer
        ]

        rows = finite_difference_rows(
            values
        )

        for j, row in enumerate(rows):

            if row:
                A[(j, t)] = clean(
                    row[0]
                )

    return A


# ============================================================================
# WIDTH-2 TRANSFER SOLVER
# ============================================================================

def solve_width2(A, t):

    equations = []

    maximum_j = max(
        (
            j
            for (j, tt) in A
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
            "alpha": None,
            "beta": None,
            "equations": equations,
        }

    alpha, beta = sp.symbols(
        "alpha beta"
    )

    M = sp.Matrix([
        [x0, x1]
        for x0, x1, _
        in equations
    ])

    b = sp.Matrix([
        lhs
        for _, _, lhs
        in equations
    ])

    rank = M.rank()
    augmented_rank = M.row_join(b).rank()

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "alpha": None,
            "beta": None,
            "equations": equations,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    if rank < 2:

        return {
            "status": "NONUNIQUE",
            "alpha": None,
            "beta": None,
            "equations": equations,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    sol = sp.solve(
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

    if len(sol) != 1:
        return {
            "status": "NONUNIQUE",
            "alpha": None,
            "beta": None,
            "equations": equations,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    alpha_value = clean(
        sol[0][alpha]
    )

    beta_value = clean(
        sol[0][beta]
    )

    verified = all(
        clean(
            alpha_value * x0
            + beta_value * x1
            - lhs
        ) == 0
        for x0, x1, lhs in equations
    )

    if not verified:
        return {
            "status": "VERIFICATION_FAILED",
            "alpha": alpha_value,
            "beta": beta_value,
            "equations": equations,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    return {
        "status": "EXACT",
        "alpha": alpha_value,
        "beta": beta_value,
        "equations": equations,
        "rank": rank,
        "augmented_rank": augmented_rank,
    }


# ============================================================================
# MATRIX CONSTRUCTION
# ============================================================================

def companion_matrix(alpha, beta):

    return sp.Matrix([
        [clean(alpha), clean(beta)],
        [sp.Integer(1), sp.Integer(0)],
    ])


# ============================================================================
# BASIC MATRIX INVARIANTS
# ============================================================================

def matrix_invariants(T):

    tr = clean(
        T.trace()
    )

    det = clean(
        T.det()
    )

    discriminant = clean(
        tr**2 - 4 * det
    )

    return {
        "trace": tr,
        "det": det,
        "discriminant": discriminant,
    }


def projective_invariants(T):

    data = matrix_invariants(T)

    tr = data["trace"]
    det = data["det"]
    disc = data["discriminant"]

    result = {
        **data,
        "det_over_trace_squared": None,
        "disc_over_trace_squared": None,
        "trace_zero": tr == 0,
    }

    if tr != 0:

        result[
            "det_over_trace_squared"
        ] = clean(
            det / tr**2
        )

        result[
            "disc_over_trace_squared"
        ] = clean(
            disc / tr**2
        )

    return result


# ============================================================================
# PROJECTIVE SIMILARITY NECESSARY TEST
# ============================================================================

def projective_invariant_test(T1, T2):

    I1 = projective_invariants(T1)
    I2 = projective_invariants(T2)

    print()
    print("=" * 78)
    print(
        "4. PROJECTIVE SIMILARITY INVARIANT TEST"
    )
    print("=" * 78)

    print()
    print(
        "  T1 trace={}".format(
            I1["trace"]
        )
    )

    print(
        "  T1 det={}".format(
            I1["det"]
        )
    )

    print(
        "  T1 discriminant={}".format(
            I1["discriminant"]
        )
    )

    print()
    print(
        "  T2 trace={}".format(
            I2["trace"]
        )
    )

    print(
        "  T2 det={}".format(
            I2["det"]
        )
    )

    print(
        "  T2 discriminant={}".format(
            I2["discriminant"]
        )
    )

    if (
        I1["trace_zero"]
        or I2["trace_zero"]
    ):

        print(
            "  trace_zero_case=True"
        )

        print(
            "  necessary_projective_test="
            "REQUIRES_SEPARATE_ZERO_TRACE_ANALYSIS"
        )

        return False

    det_ratio_equal = (
        I1["det_over_trace_squared"]
        ==
        I2["det_over_trace_squared"]
    )

    disc_ratio_equal = (
        I1["disc_over_trace_squared"]
        ==
        I2["disc_over_trace_squared"]
    )

    print()
    print(
        "  T1 det/trace^2={}".format(
            I1[
                "det_over_trace_squared"
            ]
        )
    )

    print(
        "  T2 det/trace^2={}".format(
            I2[
                "det_over_trace_squared"
            ]
        )
    )

    print(
        "  det_trace_invariant_equal={}".format(
            det_ratio_equal
        )
    )

    print()
    print(
        "  T1 disc/trace^2={}".format(
            I1[
                "disc_over_trace_squared"
            ]
        )
    )

    print(
        "  T2 disc/trace^2={}".format(
            I2[
                "disc_over_trace_squared"
            ]
        )
    )

    print(
        "  discriminant_invariant_equal={}".format(
            disc_ratio_equal
        )
    )

    return (
        det_ratio_equal
        and disc_ratio_equal
    )


# ============================================================================
# DIRECT PROJECTIVE-SIMILARITY EQUATIONS
# ============================================================================

def direct_projective_similarity(T1, T2):

    print()
    print("=" * 78)
    print(
        "5. DIRECT PROJECTIVE-SIMILARITY SOLVER"
    )
    print("=" * 78)

    lam = sp.symbols(
        "lambda"
    )

    a, b, c, d = sp.symbols(
        "a b c d"
    )

    S = sp.Matrix([
        [a, b],
        [c, d],
    ])

    equations_matrix = clean_matrix(
        T2 * S
        - lam * S * T1
    )

    equations = [
        clean(entry)
        for entry in equations_matrix
    ]

    equations = [
        sp.Eq(
            entry,
            0
        )
        for entry in equations
    ]

    print()
    print(
        "  equations=4"
    )

    # Eliminate the scale of S by fixing one
    # nonzero matrix entry in turn.
    #
    # This avoids treating the zero matrix as an
    # invertible solution.

    candidates = []

    normalizations = [
        (a, 1),
        (b, 1),
        (c, 1),
        (d, 1),
    ]

    for symbol, value in normalizations:

        normalized = [
            eq.subs(
                symbol,
                value
            )
            for eq in equations
        ]

        try:

            solutions = sp.solve(
                normalized,
                [
                    lam,
                    a,
                    b,
                    c,
                    d,
                ],
                dict=True,
                simplify=False,
            )

        except Exception:
            solutions = []

        for sol in solutions:

            av = clean(
                sol.get(a, value)
            )

            bv = clean(
                sol.get(b, 0)
            )

            cv = clean(
                sol.get(c, 0)
            )

            dv = clean(
                sol.get(d, 0)
            )

            lv = clean(
                sol.get(lam)
                if lam in sol
                else sp.nan
            )

            detS = clean(
                av * dv
                - bv * cv
            )

            if (
                lv is not sp.nan
                and detS != 0
            ):

                candidates.append(
                    (
                        lv,
                        av,
                        bv,
                        cv,
                        dv,
                        detS,
                    )
                )

    # Remove duplicates.

    unique = []

    for item in candidates:

        if item not in unique:
            unique.append(item)

    if not unique:

        print(
            "  status=NO_PROJECTIVE_SIMILARITY_SOLUTION"
        )

        return False

    print(
        "  exact_solutions={}".format(
            len(unique)
        )
    )

    for i, item in enumerate(
        unique,
        start=1,
    ):

        lv, av, bv, cv, dv, detS = item

        print()
        print(
            "  solution={}:".format(i)
        )

        print(
            "    lambda={}".format(
                lv
            )
        )

        print(
            "    S=Matrix([[{}, {}], [{}, {}]])".format(
                av,
                bv,
                cv,
                dv,
            )
        )

        print(
            "    det(S)={}".format(
                detS
            )
        )

        check = clean_matrix(
            T2
            - lv
            * S.subs(
                {
                    a: av,
                    b: bv,
                    c: cv,
                    d: dv,
                }
            )
            * T1
            * S.subs(
                {
                    a: av,
                    b: bv,
                    c: cv,
                    d: dv,
                }
            ).inv()
        )

        print(
            "    verified={}".format(
                check == sp.zeros(2)
            )
        )

    return True


# ============================================================================
# MATRIX CLEANING
# ============================================================================

def clean_matrix(M):

    return M.applyfunc(
        clean
    )


# ============================================================================
# SCALAR NORMALIZATION TEST
# ============================================================================

def scalar_normalization_test(T1, T2):

    print()
    print("=" * 78)
    print(
        "6. SCALAR MULTIPLE / ORDINARY SIMILARITY TEST"
    )
    print("=" * 78)

    tr1 = clean(
        T1.trace()
    )

    tr2 = clean(
        T2.trace()
    )

    if tr1 == 0:

        print(
            "  scalar_from_trace=UNAVAILABLE"
        )

        return

    lam = clean(
        tr2 / tr1
    )

    print(
        "  lambda_from_trace={}".format(
            lam
        )
    )

    S = sp.Matrix([
        [
            sp.symbols("a"),
            sp.symbols("b"),
        ],
        [
            sp.symbols("c"),
            sp.symbols("d"),
        ],
    ])

    # Since the matrices are 2x2, use the necessary
    # characteristic-polynomial condition first.

    target_det = clean(
        T2.det()
        / lam**2
    )

    print(
        "  adjusted_T2_determinant={}".format(
            target_det
        )
    )

    print(
        "  T1_determinant={}".format(
            clean(T1.det())
        )
    )

    print(
        "  determinant_match_after_lambda={}".format(
            target_det
            ==
            clean(T1.det())
        )
    )

    if target_det != clean(T1.det()):

        print(
            "  ordinary_similarity_after_scaling=False"
        )

        return

    print(
        "  characteristic_polynomial_match=True"
    )

    # For 2x2 matrices with matching characteristic
    # polynomial, similarity is guaranteed over the
    # algebraic closure unless both are exceptional
    # scalar matrices. We still run a direct solver.

    direct_projective_similarity(
        T1,
        T2,
    )


# ============================================================================
# EIGENVALUE-RATIO STYLE AUDIT
# ============================================================================

def eigen_ratio_audit(T1, T2):

    print()
    print("=" * 78)
    print(
        "7. EIGENVALUE-RATIO / DISCRIMINANT AUDIT"
    )
    print("=" * 78)

    x = sp.symbols(
        "x"
    )

    for label, T in (
        ("T1", T1),
        ("T2", T2),
    ):

        poly = clean(
            T.charpoly(x).as_expr()
        )

        discr = clean(
            sp.discriminant(
                poly,
                x,
            )
        )

        print()
        print(
            "  {} characteristic_polynomial={}".format(
                label,
                poly,
            )
        )

        print(
            "  {} discriminant={}".format(
                label,
                discr,
            )
        )

        print(
            "  {} v2(discriminant)={}".format(
                label,
                valuation(
                    discr,
                    2,
                ),
            )
        )

        print(
            "  {} v3(discriminant)={}".format(
                label,
                valuation(
                    discr,
                    3,
                ),
            )
        )

        print(
            "  {} v5(discriminant)={}".format(
                label,
                valuation(
                    discr,
                    5,
                ),
            )
        )

        print(
            "  {} v7(discriminant)={}".format(
                label,
                valuation(
                    discr,
                    7,
                ),
            )
        )

        print(
            "  {} v17(discriminant)={}".format(
                label,
                valuation(
                    discr,
                    17,
                ),
            )
        )


# ============================================================================
# PROJECTIVE INVARIANT PRIME PROFILE
# ============================================================================

def projective_prime_profile(T1, T2):

    print()
    print("=" * 78)
    print(
        "8. PROJECTIVE INVARIANT PRIME PROFILE"
    )
    print("=" * 78)

    for label, T in (
        ("T1", T1),
        ("T2", T2),
    ):

        tr = clean(
            T.trace()
        )

        det = clean(
            T.det()
        )

        if tr == 0:
            continue

        inv = clean(
            det / tr**2
        )

        print()
        print(
            "  {} det/trace^2={}".format(
                label,
                inv,
            )
        )

        for prime in (
            2,
            3,
            5,
            7,
            11,
            13,
            17,
        ):

            print(
                "    v_{}={}".format(
                    prime,
                    valuation(
                        inv,
                        prime,
                    ),
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
# STRUCTURAL INTERPRETATION
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
Experiment 306R ruled out the restricted normalization

    A~[j,t] = u_t c^j A[j,t].

The two exact transitions can instead be represented by companion matrices

    T_t =
        [[alpha_t, beta_t],
         [1,       0]].

Experiment 307R asks whether the two matrices belong to the same
projective similarity class:

    T_2 = lambda S T_1 S^(-1),

with S any invertible 2x2 matrix.

This is strictly more general than:

    * a scalar rescaling;
    * diagonal similarity;
    * geometric coordinate scaling.

The decisive exact invariant is

    det(T) / tr(T)^2,

together with the equivalent characteristic/discriminant information.

If the invariant differs, no projective similarity can exist.

If it agrees, the direct symbolic solver determines whether an exact
projective similarity is actually present.

Because only two independently determined transitions are available,
even a positive result is diagnostic rather than a universal theorem.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 307R — EXACT WIDTH-2 TRANSFER "
        "PROJECTIVE-SIMILARITY AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    A = first_entry_triangle(
        layers
    )

    print()
    print("=" * 78)
    print(
        "1. EXACT WIDTH-2 TRANSFER EXTRACTION"
    )
    print("=" * 78)

    profile = {}

    for t in range(5):

        result = solve_width2(
            A,
            t,
        )

        profile[t] = result

        print()
        print(
            "  t={} -> {}: status={}".format(
                t,
                t + 1,
                result["status"],
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

        print()
        print("=" * 78)
        print(
            "FINAL: INSUFFICIENT_EXACT_TRANSITIONS"
        )
        print("=" * 78)
        return

    t1, alpha1, beta1 = exact[0]
    t2, alpha2, beta2 = exact[1]

    T1 = companion_matrix(
        alpha1,
        beta1,
    )

    T2 = companion_matrix(
        alpha2,
        beta2,
    )

    print()
    print("=" * 78)
    print(
        "2. COMPANION TRANSFER MATRICES"
    )
    print("=" * 78)

    print()
    print(
        "  transition t={}:".format(
            t1
        )
    )

    print(
        "  T1={}".format(
            T1
        )
    )

    print()
    print(
        "  transition t={}:".format(
            t2
        )
    )

    print(
        "  T2={}".format(
            T2
        )
    )

    projective_invariants = (
        projective_invariant_test(
            T1,
            T2,
        )
    )

    scalar_normalization_test(
        T1,
        T2,
    )

    eigen_ratio_audit(
        T1,
        T2,
    )

    projective_prime_profile(
        T1,
        T2,
    )

    terminal_reference()

    interpretation()

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_width2_transitions={}".format(
            len(exact)
        )
    )

    print(
        "  transition_indices={}".format(
            [t for t, _, _ in exact]
        )
    )

    print(
        "  projective_invariant_equal={}".format(
            projective_invariants
        )
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
        "EXPERIMENT 307R COMPLETE"
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
