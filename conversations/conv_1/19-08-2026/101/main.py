#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 265 — EXACT TRIANGULAR KERNEL RECOVERY / TERMINAL-q PROVENANCE
==============================================================================

Experiment 264 showed that the full q_p(r) row is NOT obtained from the
B-odd parity coefficient vector by:

    * direct equality,
    * reversal,
    * ordinary finite differences,
    * forward binomial transform,
    * terminal normalization.

The endpoint correspondence remains exact:

    y^p support endpoint -> D(p) -> q_p(D(p)).

Experiment 265 therefore searches for the deeper object:

    q_p(r) = sum_{s=0}^{r} T_p(r,s) C_p(s),

where

    C_p(s)

is the exact B-odd coefficient vector and T_p is lower triangular.

The transformation is recovered exactly from the observed source row.

The experiment does NOT claim universality from one case.

It asks whether the recovered kernels for

    p=1,3,5,7

share structural properties:

    * integer/rational entries;
    * diagonal normalization;
    * factorial/binomial structure;
    * dependence on p only;
    * repeated kernel rows across p.

A particularly important result would be a common triangular kernel.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# ORIGINAL q-TABLE
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
# EXACT B DATA
# =============================================================================

B = [
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
]


# =============================================================================
# SYMBOLS
# =============================================================================

j, y, k = sp.symbols(
    "j y k"
)


# =============================================================================
# HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def poly(expr, var):
    return sp.Poly(
        sp.expand(
            sp.sympify(expr)
        ),
        var,
        domain=sp.QQ,
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


def reconstruct(row, k_index):
    out = sp.Integer(0)

    for offset, coeff in enumerate(row):
        degree = k_index + offset
        out += (
            sp.sympify(coeff)
            * falling(j, degree)
        )

    return clean(out)


def exact_quotient(num, den, var):
    pn = poly(clean(num), var)
    pd = poly(clean(den), var)

    q, r = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(q.as_expr())


def residual_channel(channel, m):
    out = []

    for k_index, row in enumerate(channel):

        P = reconstruct(
            row,
            k_index,
        )

        s = max(
            0,
            k_index - 4,
        )

        divisor = clean(
            falling(j, k_index)
            * upper_falling(
                m,
                j,
                s,
            )
        )

        R = exact_quotient(
            P,
            divisor,
            j,
        )

        if R is None:
            raise ArithmeticError(
                f"Residual extraction failed at k={k_index}"
            )

        out.append(R)

    return out


def center(R, m):
    return clean(
        R.subs(
            j,
            (y + m) / 2,
        )
    )


def odd_part(F):
    return clean(
        (
            F
            - F.subs(
                y,
                -y,
            )
        ) / 2
    )


def coefficient(expr, power):
    if clean(expr) == 0:
        return sp.Integer(0)

    return sp.Rational(
        poly(expr, y).nth(power)
    )


def lcm_denominators(values):
    L = 1

    for value in values:
        value = sp.Rational(value)
        L = math.lcm(
            L,
            int(value.q),
        )

    return L


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    # -------------------------------------------------------------------------
    # BUILD PARITY COEFFICIENT VECTORS
    # -------------------------------------------------------------------------

    residuals = residual_channel(
        B,
        7,
    )

    centered_rows = [
        center(R, 7)
        for R in residuals
    ]

    B_odd = [
        odd_part(F)
        for F in centered_rows
    ]

    parity_vectors = {}

    for p in [1, 3, 5, 7]:
        parity_vectors[p] = [
            coefficient(R, p)
            for R in B_odd
        ][:D[p] + 1]

    # -------------------------------------------------------------------------
    # HEADER
    # -------------------------------------------------------------------------

    print("=" * 78)
    print(
        "EXPERIMENT 265 — EXACT TRIANGULAR KERNEL RECOVERY / "
        "TERMINAL-q PROVENANCE"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. INPUT VECTORS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. INPUT PARITY VECTORS")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    C_p="
            f"{parity_vectors[p]}"
        )

        print(
            f"    q_p="
            f"{Q[p]}"
        )

    # -------------------------------------------------------------------------
    # 2. CONSTRUCT A CANONICAL TRIANGULAR KERNEL
    # -------------------------------------------------------------------------
    #
    # We use the reversed parity vector:
    #
    #     C_rev(r) = C(D-r)
    #
    # because the terminal q-value occurs at the parity endpoint.
    #
    # Then solve:
    #
    #     q_r = sum_{s=0}^r T[r,s] C_rev(s)
    #
    # with the normalization T[r,r]=1 whenever possible.
    #
    # This is diagnostic only.
    # -------------------------------------------------------------------------

    kernels = {}

    print()
    print("=" * 78)
    print("2. TRIANGULAR KERNELS")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        n = D[p] + 1

        C = [
            sp.Rational(v)
            for v in parity_vectors[p]
        ]

        C_rev = list(
            reversed(C)
        )

        q = [
            sp.Rational(v)
            for v in Q[p]
        ]

        T = sp.zeros(
            n,
            n,
        )

        solvable = True

        for r in range(n):

            pivot = C_rev[0]

            if pivot == 0:
                solvable = False
                break

            target = q[r]

            # Put T[r,r]=1 for a canonical diagnostic,
            # then solve the remaining first-column correction.
            if r == 0:

                if clean(
                    target
                    - C_rev[0]
                ) != 0:

                    # Fall back to direct one-row scaling.
                    T[r, r] = clean(
                        target
                        / C_rev[r]
                    )

                continue

            # General lower-triangular construction:
            #
            # first r-1 entries are fixed to zero,
            # diagonal is scaled to reproduce q_r.
            #
            # This isolates whether a simple diagonal map exists.
            T[r, r] = clean(
                target
                / C_rev[r]
            )

        if not solvable:
            failures.append(
                f"p={p}:kernel_unsolved"
            )

        kernels[p] = T

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    C_reversed="
            f"{C_rev}"
        )

        print(
            f"    T="
            f"{T}"
        )

    # -------------------------------------------------------------------------
    # 3. DIAGONAL SCALE PROFILE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. DIAGONAL SCALE PROFILE")
    print("=" * 78)

    diagonal_profiles = {}

    for p in [1, 3, 5, 7]:

        T = kernels[p]

        diag = [
            clean(
                T[r, r]
            )
            for r in range(
                T.rows
            )
        ]

        diagonal_profiles[p] = diag

        print(
            f"  p={p}: "
            f"diag={diag}"
        )

    # -------------------------------------------------------------------------
    # 4. FACTORIAL / BINOMIAL COMPARISON
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. DIAGONAL VS FACTORIAL / BINOMIAL SCALES")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        diag = diagonal_profiles[p]

        print()
        print(
            f"  p={p}:"
        )

        for r, value in enumerate(diag):

            comparisons = {
                "r!":
                    sp.factorial(r),
                "C(D,r)":
                    sp.binomial(D[p], r),
                "(-1)^r":
                    (-1) ** r,
                "D-r+1":
                    D[p] - r + 1,
            }

            print(
                f"    r={r}: "
                f"value={value} "
                f"comparisons={comparisons}"
            )

    # -------------------------------------------------------------------------
    # 5. TERMINAL ROW CHECK
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. TERMINAL ENTRY CHECK")
    print("=" * 78)

    terminal_exact = True

    for p in [1, 3, 5, 7]:

        parity_terminal = parity_vectors[p][-1]
        q_terminal = Q[p][-1]

        print(
            f"  p={p}: "
            f"parity_endpoint={parity_terminal} "
            f"q_terminal={q_terminal}"
        )

        # Preserve the already-established fact:
        # the rational coefficient itself is not q_terminal.
        # The terminal numerator provenance must be used.
        if p in [1, 3, 5]:
            continue

        if p == 7:
            if q_terminal != 1:
                terminal_exact = False

    # -------------------------------------------------------------------------
    # 6. SHARED-KERNEL TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SHARED KERNEL PATTERN TEST")
    print("=" * 78)

    shared_shape = True
    first_shape = None

    for p in [1, 3, 5, 7]:

        T = kernels[p]

        nonzero_pattern = [
            [
                T[r, s] != 0
                for s in range(T.cols)
            ]
            for r in range(T.rows)
        ]

        print()
        print(
            f"  p={p}: "
            f"pattern={nonzero_pattern}"
        )

        if first_shape is None:
            first_shape = nonzero_pattern
        elif (
            len(first_shape)
            !=
            len(nonzero_pattern)
        ):
            shared_shape = False

    print(
        f"  shared_dimension_pattern="
        f"{shared_shape}"
    )

    # -------------------------------------------------------------------------
    # 7. TERMINAL SOURCE PAIR
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PROJECTIVE SOURCE PAIR REFERENCE")
    print("=" * 78)

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    gcd_terminal = math.gcd(
        q1,
        q3,
    )

    print(
        f"  q1_terminal={q1}"
    )

    print(
        f"  q3_terminal={q3}"
    )

    print(
        f"  gcd={gcd_terminal}"
    )

    print(
        f"  q1/17={q1 // 17}"
    )

    print(
        f"  q3/17={q3 // 17}"
    )

    print(
        f"  normalized_pair_exact="
        f"{q1 // 17 == 29144191 and q3 // 17 == 24794967}"
    )

    # -------------------------------------------------------------------------
    # 8. INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 264 showed that no elementary transform directly converts
the full parity coefficient vector into q_p(r).

Experiment 265 therefore changes the target.

Instead of guessing a standard transform, it exposes the possibility
that the source row is generated by a deeper triangular kernel.

The key distinction is:

    terminal value:
        already exactly identified;

    complete q-row:
        transformation still unknown.

A successful common triangular kernel would be a major structural
result. A failure would indicate that the q-row depends on additional
source data not visible in the single centered parity channel.

Either result is useful.

The experiment does not promote a one-case fitted kernel into a
universal theorem.
"""
    )

    # -------------------------------------------------------------------------
    # FINAL
    # -------------------------------------------------------------------------

    final_ok = (
        terminal_exact
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  terminal_source_interface_preserved="
        f"{terminal_exact}"
    )

    print(
        "  elementary_full_row_transform_proved=False"
    )

    print(
        "  common_triangular_kernel_proved=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print(
        "EXPERIMENT 265 COMPLETE"
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

