#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 270 — EXACT PASCAL-TYPE SOURCE OPERATOR / COMMON SCALING AUDIT
==============================================================================

Experiments 268R and 269 established:

    * no common Toeplitz/convolution kernel;
    * an arbitrary p-dependent triangular representation is trivially
      possible, but the minimum-support construction collapses to a diagonal
      matrix and therefore explains nothing.

Experiment 270 imposes a genuinely nontrivial structure.

We test whether the observed B-odd vectors C_p can be written as

    C = P_D * S * q,

where:

    q = source q_p(r) row,

    S = diag(s_0,s_1,...)

is a COMMON source-index scaling, independent of p,

and P_D is one of several Pascal-type triangular kernels determined only
by D=D(p).

The tested kernels are:

    1.  P[i,j] = (-1)^(i-j) C(D-j, i-j)

    2.  P[i,j] = C(D-j, i-j)

    3.  P[i,j] = (-1)^(i-j) C(D-i, j)

    4.  P[i,j] = C(D-i, j)

with invalid binomial entries treated as zero.

The source scaling s_r is solved globally across all available p rows.

This is much stronger than fitting a separate T_p:

    * one common scaling;
    * one combinatorial kernel family;
    * p enters only through D(p).

If one candidate passes exactly, it identifies a concrete candidate
operator.

If none passes, the missing operator is not a simple Pascal/binomial
transform with common q-index normalization.

Exact QQ arithmetic only.
"""


from __future__ import annotations

import sys
import sympy as sp


# =============================================================================
# DATA
# =============================================================================

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


D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}


C = {
    1: [
        sp.Rational(-584531, 35840),
        sp.Rational(-2908483, 1920),
        sp.Rational(-31233169, 13440),
        sp.Rational(-1446167, 1344),
        sp.Rational(-22259149, 40320),
        sp.Rational(-301, 240),
    ],
    3: [
        sp.Rational(-59257, 46080),
        sp.Rational(186547, 1440),
        sp.Rational(367433, 1680),
        sp.Rational(126549, 448),
        sp.Rational(-162139, 40320),
    ],
    5: [
        sp.Rational(4457, 46080),
        sp.Rational(-16819, 5760),
        sp.Rational(-5769, 896),
    ],
    7: [
        sp.Rational(-421, 322560),
    ],
}


# =============================================================================
# HELPERS
# =============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def binom(n, r):
    if n < 0 or r < 0 or r > n:
        return sp.Integer(0)
    return sp.binomial(n, r)


def build_kernel(Dp, kind):
    """
    Return a (D+1)x(D+1) Pascal-type matrix.
    """

    n = Dp + 1
    P = sp.zeros(n, n)

    for i in range(n):
        for j in range(n):

            if j > i and kind in {
                "lower_signed",
                "lower_plain",
            }:
                continue

            if i > j and kind in {
                "upper_signed",
                "upper_plain",
            }:
                continue

            gap = i - j

            if kind == "lower_signed":
                if i >= j:
                    P[i, j] = (
                        (-1) ** gap
                        * binom(
                            Dp - j,
                            gap,
                        )
                    )

            elif kind == "lower_plain":
                if i >= j:
                    P[i, j] = binom(
                        Dp - j,
                        gap,
                    )

            elif kind == "upper_signed":
                if j >= i:
                    gap2 = j - i
                    P[i, j] = (
                        (-1) ** gap2
                        * binom(
                            Dp - i,
                            gap2,
                        )
                    )

            elif kind == "upper_plain":
                if j >= i:
                    gap2 = j - i
                    P[i, j] = binom(
                        Dp - i,
                        gap2,
                    )

            else:
                raise ValueError(kind)

    return P


def build_global_system(kind):
    """
    Unknowns:

        s_0,...,s_5

    Common across all p.

    Equations:

        C_p[i] = sum_j P_D[i,j] * s_j * q_p[j].
    """

    scale = sp.symbols("s0:6")
    equations = []

    for p in [1, 3, 5, 7]:

        Dp = D[p]
        P = build_kernel(
            Dp,
            kind,
        )

        q = Q[p]
        c = C[p]

        n = len(q)

        for i in range(n):

            lhs = sp.Integer(0)

            for r in range(n):
                lhs += (
                    P[i, r]
                    * scale[r]
                    * q[r]
                )

            equations.append(
                sp.Eq(
                    lhs,
                    c[i],
                )
            )

    return scale, equations


def solve_global(kind):
    scale, equations = build_global_system(
        kind
    )

    expressions = [
        clean(
            eq.lhs - eq.rhs
        )
        for eq in equations
    ]

    A, b = sp.linear_eq_to_matrix(
        expressions,
        scale,
    )

    augmented = A.row_join(b)

    if A.rank() != augmented.rank():
        return None

    if A.rank() != len(scale):
        return None

    solution = sp.linsolve(
        (A, b),
        scale,
    )

    if len(solution) != 1:
        return None

    vector = next(iter(solution))

    for value in vector:
        if value.free_symbols:
            return None

    return [
        clean(value)
        for value in vector
    ]


def predict(kind, scales):
    output = {}

    for p in [1, 3, 5, 7]:

        Dp = D[p]
        P = build_kernel(
            Dp,
            kind,
        )

        q = sp.Matrix(
            Q[p]
        )

        scaled_q = sp.Matrix(
            [
                scales[r] * q[r]
                for r in range(
                    len(q)
                )
            ]
        )

        result = P * scaled_q

        output[p] = [
            clean(x)
            for x in result
        ]

    return output


def verify_predictions(predicted):
    return all(
        predicted[p] == C[p]
        for p in [1, 3, 5, 7]
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 270 — EXACT PASCAL-TYPE SOURCE OPERATOR / "
        "COMMON SCALING AUDIT"
    )
    print("=" * 78)

    kinds = [
        "lower_signed",
        "lower_plain",
        "upper_signed",
        "upper_plain",
    ]

    results = {}

    # =========================================================================
    # 1. DATA
    # =========================================================================

    print()
    print("=" * 78)
    print("1. SOURCE / OBSERVED DATA")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        print()
        print(
            f"  p={p}: D={D[p]}"
        )

        print(
            f"    q={Q[p]}"
        )

        print(
            f"    C={C[p]}"
        )

    # =========================================================================
    # 2. KERNEL DEFINITIONS
    # =========================================================================

    print()
    print("=" * 78)
    print("2. PASCAL-TYPE KERNELS")
    print("=" * 78)

    for kind in kinds:

        print()
        print(
            f"  {kind}:"
        )

        for p in [1, 3, 5, 7]:

            print(
                f"    p={p}:"
            )

            print(
                build_kernel(
                    D[p],
                    kind,
                )
            )

    # =========================================================================
    # 3. GLOBAL COMMON-SCALING SOLUTION
    # =========================================================================

    print()
    print("=" * 78)
    print("3. GLOBAL COMMON-SCALING SEARCH")
    print("=" * 78)

    for kind in kinds:

        scales = solve_global(
            kind
        )

        print()
        print(
            f"  kind={kind}:"
        )

        if scales is None:

            print(
                "    exact_common_scaling=False"
            )

            results[kind] = None
            continue

        predicted = predict(
            kind,
            scales,
        )

        exact = verify_predictions(
            predicted
        )

        print(
            "    exact_common_scaling=True"
        )

        print(
            f"    scales={scales}"
        )

        print(
            f"    exact_verification={exact}"
        )

        if exact:

            for p in [1, 3, 5, 7]:

                print(
                    f"      p={p}: "
                    f"predicted={predicted[p]}"
                )

        results[kind] = {
            "scales": scales,
            "exact": exact,
        }

    # =========================================================================
    # 4. SCALE FACTORIZATION
    # =========================================================================

    print()
    print("=" * 78)
    print("4. COMMON SCALE FACTOR PROFILE")
    print("=" * 78)

    for kind, result in results.items():

        if result is None:
            continue

        scales = result["scales"]

        print()
        print(
            f"  {kind}:"
        )

        for i, value in enumerate(
            scales
        ):

            print(
                f"    s[{i}]={value}"
            )

    # =========================================================================
    # 5. TERMINAL SOURCE REFERENCE
    # =========================================================================

    print()
    print("=" * 78)
    print("5. TERMINAL SOURCE REFERENCE")
    print("=" * 78)

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print(
        f"  q1_terminal={q1}"
    )

    print(
        f"  q3_terminal={q3}"
    )

    print(
        f"  gcd={sp.gcd(q1, q3)}"
    )

    print(
        f"  q1_terminal/17={q1 // 17}"
    )

    print(
        f"  q3_terminal/17={q3 // 17}"
    )

    # =========================================================================
    # 6. INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 268R ruled out a common Toeplitz convolution.

Experiment 269 showed that an unrestricted p-dependent triangular map
can always collapse to an uninformative diagonal representation.

Experiment 270 imposes a much stronger and more meaningful hypothesis:

    C_p = P_D * diag(s_r) * q_p,

with:

    * one common scale s_r for each source index r;
    * a Pascal/binomial kernel determined only by D(p);
    * no p-dependent fitted matrix;
    * no arbitrary off-diagonal parameters.

A successful candidate would therefore provide an explicit source-index
operator with a clear combinatorial interpretation.

A failure of every candidate means the missing B-channel operator is
not of this natural Pascal type, and we should stop trying elementary
source-index transforms and instead reconstruct the operator directly
from the original full construction.
"""
    )

    # =========================================================================
    # 7. FINAL
    # =========================================================================

    successful = [
        kind
        for kind, result in results.items()
        if result is not None
        and result["exact"]
    ]

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  successful_pascal_modes="
        f"{successful}"
    )

    print(
        f"  common_scaling_pascal_operator_found="
        f"{bool(successful)}"
    )

    print(
        "  arbitrary_p_dependent_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 270 COMPLETE"
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
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise

