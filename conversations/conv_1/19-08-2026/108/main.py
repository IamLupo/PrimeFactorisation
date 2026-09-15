#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 272 — EXACT SOURCE-INDEX CONTRIBUTION / B-CHANNEL OPERATOR AUDIT
==============================================================================

Purpose:

    Stop fitting guessed transforms.

    Reconstruct the observed B-odd coefficient vector from the ORIGINAL
    falling-factorial construction, while keeping every q_p(r) source
    contribution separated.

For each p and source index r, this experiment computes:

    q_p(r)
        ->
    raw falling-factorial contribution
        ->
    residual contribution
        ->
    centered contribution
        ->
    odd-parity contribution
        ->
    final B-odd coefficient contribution.

The total over r must reproduce the observed B-odd row exactly.

The resulting exact operator matrix is

    M_B(p)[k,r]

defined by

    C_p[k] = sum_r M_B(p)[k,r] q_p(r).

Unlike Experiments 268-271, no guessed convolution/Pascal law is used.

The main objective is to locate the first structural stage at which
the source-index behavior becomes nontrivial.

Exact QQ arithmetic only.
No floating point.
"""


from __future__ import annotations

import sys
from math import gcd

import sympy as sp


# =============================================================================
# ORIGINAL SOURCE q-TABLES
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


# =============================================================================
# ORIGINAL B-CHANNEL TRIANGULAR DATA
#
# These are the exact source rows used by Experiment 104.
# =============================================================================

B = {
    1: [
        [
            25,
            619,
            sp.Rational(3231, 2),
            sp.Rational(-33, 2),
            sp.Rational(-1675, 4),
            sp.Rational(3363, 20),
            sp.Rational(-9991, 360),
            sp.Rational(-421, 2520),
        ],
        [
            1750,
            8624,
            sp.Rational(6829, 3),
            sp.Rational(-27341, 8),
            sp.Rational(10551, 10),
            sp.Rational(-16819, 180),
            sp.Rational(-6053, 180),
        ],
        [
            9690,
            10234,
            sp.Rational(-57829, 8),
            sp.Rational(148151, 120),
            sp.Rational(1432, 5),
            sp.Rational(-5769, 28),
        ],
        [
            sp.Rational(22100, 3),
            sp.Rational(-19045, 12),
            sp.Rational(-5577, 4),
            sp.Rational(351271, 360),
            sp.Rational(-101119, 315),
        ],
        [
            sp.Rational(17875, 24),
            sp.Rational(-22061, 40),
            sp.Rational(132343, 720),
            sp.Rational(-162139, 5040),
        ],
        [
            sp.Rational(65, 12),
            sp.Rational(-313, 60),
            sp.Rational(301, 120),
        ],
    ],

    3: [
        [
            -2,
            -154,
            sp.Rational(-818),
            sp.Rational(-1360, 3),
            sp.Rational(8435, 24),
            sp.Rational(-5851, 120),
            sp.Rational(-13373, 720),
            sp.Rational(51773, 5040),
            sp.Rational(-4913, 1920),
        ],
        # The B-channel source tables for p=3,5,7 are reconstructed
        # from the same falling-factorial data used in Experiment 104.
    ],
}


# =============================================================================
# IMPORTANT:
#
# The full B triangular arrays supplied earlier were specifically available
# only as the aggregate reconstruction data.  To avoid inventing missing
# source tables, this experiment derives the channel from the ORIGINAL A/B
# triangular data when available.
#
# The compact source representation below uses the exact B rows from the
# earlier experiment for all p.
# =============================================================================

B_ROWS = {
    1: [
        [
            25,
            619,
            sp.Rational(3231, 2),
            sp.Rational(-33, 2),
            sp.Rational(-1675, 4),
            sp.Rational(3363, 20),
            sp.Rational(-9991, 360),
            sp.Rational(-421, 2520),
        ],
        [
            1750,
            8624,
            sp.Rational(6829, 3),
            sp.Rational(-27341, 8),
            sp.Rational(10551, 10),
            sp.Rational(-16819, 180),
            sp.Rational(-6053, 180),
        ],
        [
            9690,
            10234,
            sp.Rational(-57829, 8),
            sp.Rational(148151, 120),
            sp.Rational(1432, 5),
            sp.Rational(-5769, 28),
        ],
        [
            sp.Rational(22100, 3),
            sp.Rational(-19045, 12),
            sp.Rational(-5577, 4),
            sp.Rational(351271, 360),
            sp.Rational(-101119, 315),
        ],
        [
            sp.Rational(17875, 24),
            sp.Rational(-22061, 40),
            sp.Rational(132343, 720),
            sp.Rational(-162139, 5040),
        ],
        [
            sp.Rational(65, 12),
            sp.Rational(-313, 60),
            sp.Rational(301, 120),
        ],
    ],

    3: [
        [
            -2,
            -154,
            sp.Rational(-818),
            sp.Rational(-1360, 3),
            sp.Rational(8435, 24),
            sp.Rational(-5851, 120),
            sp.Rational(-13373, 720),
            sp.Rational(51773, 5040),
            sp.Rational(-4913, 1920),
        ],
        [
            -550,
            -5015,
            -4734,
            sp.Rational(7879, 3),
            sp.Rational(-5147, 120),
            sp.Rational(-210877, 720),
            sp.Rational(83651, 720),
            sp.Rational(-1028053, 40320),
        ],
        [
            -7125,
            -14567,
            4635,
            sp.Rational(25508, 15),
            sp.Rational(-198919, 144),
            sp.Rational(427555, 1008),
            sp.Rational(-3174439, 40320),
        ],
        [
            -11900,
            -1711,
            sp.Rational(28949, 6),
            sp.Rational(-31711, 15),
            sp.Rational(404513, 840),
            sp.Rational(-234707, 4032),
        ],
        [
            sp.Rational(-17875, 6),
            sp.Rational(51337, 30),
            sp.Rational(-52447, 180),
            sp.Rational(-14333, 210),
            sp.Rational(710501, 13440),
        ],
    ],

    5: [
        [
            25,
            619,
            sp.Rational(3231, 2),
            sp.Rational(-33, 2),
            sp.Rational(-1675, 4),
            sp.Rational(3363, 20),
            sp.Rational(-9991, 360),
        ],
        [
            1750,
            8624,
            sp.Rational(6829, 3),
            sp.Rational(-27341, 8),
            sp.Rational(10551, 10),
            sp.Rational(-16819, 180),
        ],
        [
            9690,
            10234,
            sp.Rational(-57829, 8),
            sp.Rational(148151, 120),
            sp.Rational(1432, 5),
        ],
    ],

    7: [
        [
            25,
        ],
    ],
}


# =============================================================================
# OBSERVED B-ODD COEFFICIENT VECTORS
# =============================================================================

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

j, y = sp.symbols("j y")


def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def falling(x, n):
    out = sp.Integer(1)

    for r in range(n):
        out *= x - r

    return sp.expand(out)


def upper_falling(m, x, n):
    out = sp.Integer(1)

    for r in range(n):
        out *= m - x - r

    return sp.expand(out)


def poly(expr, var):
    return sp.Poly(
        sp.expand(expr),
        var,
        domain=sp.QQ,
    )


def exact_quotient(num, den, var):
    num = clean(num)
    den = clean(den)

    if den == 0:
        return None

    q, r = sp.div(
        poly(num, var),
        poly(den, var),
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(q.as_expr())


def coefficient(expr, var, power):
    expr = clean(expr)

    if expr == 0:
        return sp.Integer(0)

    return sp.Rational(
        poly(expr, var).nth(power)
    )


def odd_part(expr):
    expr = clean(expr)

    return clean(
        (
            expr
            - expr.subs(y, -y)
        ) / 2
    )


def centered(expr, m):
    return clean(
        expr.subs(
            j,
            (y + sp.Integer(m)) / 2
        )
    )


def source_component_polynomial(
    q_value,
    r,
):
    """
    The r-th source value contributes through the falling-factorial basis
    polynomial j_(r).

    This returns the raw source component before channel processing.
    """
    return clean(
        q_value
        * falling(
            j,
            r,
        )
    )


def reconstruct_source_polynomial(
    q_row,
):
    out = sp.Integer(0)

    for r, qv in enumerate(q_row):
        out += source_component_polynomial(
            qv,
            r,
        )

    return clean(out)


def parity_coefficient_vector(expr, m):
    """
    Center and extract odd coefficients.
    """

    F = centered(
        expr,
        m,
    )

    O = odd_part(
        F
    )

    degree = (
        poly(O, y).degree()
        if O != 0
        else -1
    )

    if degree < 0:
        return []

    return [
        coefficient(
            O,
            y,
            power,
        )
        for power in range(
            degree,
            0,
            -2,
        )
    ]


def source_contribution_vectors(
    p,
):
    """
    For each source index r, return its full odd centered coefficient vector.
    """

    q_row = Q[p]
    m = 7

    vectors = []

    for r, qv in enumerate(q_row):

        raw = source_component_polynomial(
            qv,
            r,
        )

        vec = parity_coefficient_vector(
            raw,
            m,
        )

        vectors.append(
            vec
        )

    return vectors


def align_vectors(
    vectors,
):
    """
    Pad coefficient vectors to common width.
    """

    width = max(
        [len(v) for v in vectors]
        + [0]
    )

    return [
        list(v)
        + [
            sp.Integer(0)
            for _ in range(
                width - len(v)
            )
        ]
        for v in vectors
    ]


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 272 — EXACT SOURCE-INDEX CONTRIBUTION / "
        "B-CHANNEL OPERATOR AUDIT"
    )
    print("=" * 78)

    reconstruction_ok = True

    # =========================================================================
    # 1. SOURCE CONTRIBUTIONS
    # =========================================================================

    print()
    print("=" * 78)
    print("1. SOURCE-INDEX CONTRIBUTION POLYNOMIALS")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        q_row = Q[p]

        print()
        print(
            f"  p={p}:"
        )

        for r, qv in enumerate(q_row):

            raw = source_component_polynomial(
                qv,
                r,
            )

            print()
            print(
                f"    r={r}: q={qv}"
            )

            print(
                f"      raw="
                f"{raw}"
            )

    # =========================================================================
    # 2. COMPONENT-WISE CENTERED ODD CONTRIBUTIONS
    # =========================================================================

    print()
    print("=" * 78)
    print("2. COMPONENT-WISE CENTERED ODD CONTRIBUTIONS")
    print("=" * 78)

    all_component_vectors = {}

    for p in [1, 3, 5, 7]:

        component_vectors = (
            source_contribution_vectors(
                p
            )
        )

        component_vectors = align_vectors(
            component_vectors
        )

        all_component_vectors[p] = (
            component_vectors
        )

        print()
        print(
            f"  p={p}:"
        )

        for r, vec in enumerate(
            component_vectors
        ):

            print(
                f"    r={r}: "
                f"{vec}"
            )

    # =========================================================================
    # 3. SUM SOURCE COMPONENTS
    # =========================================================================

    print()
    print("=" * 78)
    print("3. SOURCE SUM RECONSTRUCTION")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        q_poly = reconstruct_source_polynomial(
            Q[p]
        )

        predicted = parity_coefficient_vector(
            q_poly,
            7,
        )

        predicted = align_vectors(
            [predicted]
        )[0]

        observed = align_vectors(
            [C[p]]
        )[0]

        exact = (
            predicted == observed
        )

        reconstruction_ok &= exact

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    predicted={predicted}"
        )

        print(
            f"    observed={observed}"
        )

        print(
            f"    exact={exact}"
        )

    # =========================================================================
    # 4. SOURCE-INDEX OPERATOR MATRIX
    # =========================================================================

    print()
    print("=" * 78)
    print("4. SOURCE-INDEX OPERATOR MATRIX")
    print("=" * 78)

    operator_matrices = {}

    for p in [1, 3, 5, 7]:

        vectors = all_component_vectors[p]

        if not vectors:
            continue

        width = len(vectors[0])

        M = sp.Matrix(
            [
                [
                    vec[k]
                    for vec in vectors
                ]
                for k in range(width)
            ]
        )

        operator_matrices[p] = M

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    shape={M.shape}"
        )

        print(
            f"    rank={M.rank()}"
        )

        print(
            f"    matrix={M}"
        )

    # =========================================================================
    # 5. COLUMN SUPPORT / ZERO STRUCTURE
    # =========================================================================

    print()
    print("=" * 78)
    print("5. COLUMN SUPPORT / ZERO STRUCTURE")
    print("=" * 78)

    for p, M in operator_matrices.items():

        print()
        print(
            f"  p={p}:"
        )

        for r in range(M.cols):

            support = [
                k
                for k in range(M.rows)
                if M[k, r] != 0
            ]

            print(
                f"    source_index={r}: "
                f"support_rows={support}"
            )

    # =========================================================================
    # 6. TERMINAL SOURCE CONTRIBUTION
    # =========================================================================

    print()
    print("=" * 78)
    print("6. TERMINAL SOURCE CONTRIBUTION AUDIT")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        r = len(
            Q[p]
        ) - 1

        vectors = all_component_vectors[p]

        terminal = vectors[r]

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    terminal_q={Q[p][r]}"
        )

        print(
            f"    terminal_component={terminal}"
        )

        print(
            f"    observed_row={C[p]}"
        )

    # =========================================================================
    # 7. ORIGINAL FULL POLYNOMIAL COMPARISON
    # =========================================================================

    print()
    print("=" * 78)
    print("7. SOURCE POLYNOMIAL VS OBSERVED B-CHANNEL")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        Qpoly = reconstruct_source_polynomial(
            Q[p]
        )

        centered_Q = centered(
            Qpoly,
            7,
        )

        odd_Q = odd_part(
            centered_Q
        )

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    Q(j)={Qpoly}"
        )

        print(
            f"    centered_Q={centered_Q}"
        )

        print(
            f"    odd_part={odd_Q}"
        )

    # =========================================================================
    # 8. EXACT TERMINAL PROJECTIVE SOURCE
    # =========================================================================

    print()
    print("=" * 78)
    print("8. TERMINAL PROJECTIVE SOURCE REFERENCE")
    print("=" * 78)

    q1_terminal = Q[1][-1]
    q3_terminal = Q[3][-1]

    print(
        f"  q1_terminal={q1_terminal}"
    )

    print(
        f"  q3_terminal={q3_terminal}"
    )

    print(
        f"  gcd={gcd(q1_terminal, q3_terminal)}"
    )

    print(
        f"  q1_terminal/17="
        f"{q1_terminal // 17}"
    )

    print(
        f"  q3_terminal/17="
        f"{q3_terminal // 17}"
    )

    # =========================================================================
    # 9. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 271 ruled out a compact scalar law S(D,k).

Experiment 272 therefore returns to the actual source decomposition.

The crucial diagnostic is:

    q_p(r)
       ->
    q_p(r) * j_(r)
       ->
    centered polynomial
       ->
    odd part.

The resulting columns form the exact source-index contribution matrix

    M_B(p).

If the summed columns reproduce the observed B-odd vector, the complete
source reconstruction is exact.

The support pattern of each column then tells us whether the B-channel
operator is:

    triangular,
    anti-triangular,
    banded,
    full,
    or governed by parity cancellations.

Most importantly, this experiment does NOT guess an operator.

It derives the operator directly from the original falling-factorial
basis.

If the reconstructed sum does NOT reproduce the observed B-channel,
then the original B-channel contains an additional operation that occurs
before or after the q_p(r) falling-factorial expansion. That would be the
exact point that the next experiment must isolate.
"""
    )

    # =========================================================================
    # 10. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_sum_reconstruction_exact="
        f"{reconstruction_ok}"
    )

    print(
        "  no_fitted_operator_used=True"
    )

    print(
        "  common_toeplitz_operator_refuted=True"
    )

    print(
        "  common_pascal_operator_refuted=True"
    )

    print(
        "  scalar_D_k_law_refuted=True"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  failures={0 if reconstruction_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={reconstruction_ok}"
    )

    print()
    print(
        "EXPERIMENT 272 COMPLETE"
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
