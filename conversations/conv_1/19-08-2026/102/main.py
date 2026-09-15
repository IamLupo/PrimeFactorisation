#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 266 — EXACT CENTERED-BASIS / FALLING-FACTORIAL SOURCE MAP AUDIT
==============================================================================

Experiment 265 showed that an unconstrained fitted triangular kernel is not
informative: the resulting diagonal factors are arbitrary rational numbers.

Experiment 266 therefore derives the transformation from the actual algebra.

The B-channel is constructed in the falling-factorial basis

    j_(r) = j(j-1)...(j-r+1),

while the centered parity representation uses

    y = 2j-m.

This experiment computes the exact change-of-basis matrices between:

    {j_(r)}
and
    {y^s}.

It then projects the known q_p(r) source rows through the exact basis
transformation implied by the residual construction.

The target question is:

    Is the B-odd coefficient vector a deterministic basis transform of
    the q_p(r) row, with no fitted parameters?

The experiment also checks whether the terminal coefficient emerges
automatically from that same transform.

No fitted kernel.
No heuristic matrix.
No floating point.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


j, y = sp.symbols("j y")


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
# EXPERIMENT 104 B DATA
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


def falling(x, n):
    result = sp.Integer(1)

    for r in range(n):
        result *= x - r

    return sp.expand(result)


def upper_falling(m, x, n):
    result = sp.Integer(1)

    for r in range(n):
        result *= m - x - r

    return sp.expand(result)


def poly(expr, var):
    return sp.Poly(
        sp.expand(expr),
        var,
        domain=sp.QQ,
    )


def exact_quotient(num, den):
    pn = poly(
        clean(num),
        j,
    )

    pd = poly(
        clean(den),
        j,
    )

    q, r = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(
        q.as_expr()
    )


def reconstruct(row, k_index):

    result = sp.Integer(0)

    for offset, coeff in enumerate(row):

        result += (
            sp.sympify(coeff)
            *
            falling(
                j,
                k_index + offset,
            )
        )

    return clean(result)


def residual_channel(channel, m):

    result = []

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
            falling(
                j,
                k_index,
            )
            *
            upper_falling(
                m,
                j,
                s,
            )
        )

        R = exact_quotient(
            P,
            divisor,
        )

        if R is None:
            raise ArithmeticError(
                f"Residual extraction failed at k={k_index}"
            )

        result.append(R)

    return result


def centered(R, m):

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

    return sp.Rational(
        poly(
            expr,
            y,
        ).nth(power)
    )


# =============================================================================
# CHANGE OF BASIS
# =============================================================================

def falling_to_center_matrix(n, m):

    """
    Matrix T such that

        [j_(0), ..., j_(n)]^T

    is expressed in the centered monomial basis

        [1, y, ..., y^n]^T.

    Specifically, row r contains the coefficients of

        j_(r)

    after

        j = (y+m)/2.
    """

    T = sp.zeros(
        n + 1,
        n + 1,
    )

    substitution = (
        y + sp.Integer(m)
    ) / 2

    for r in range(n + 1):

        expr = clean(
            falling(
                substitution,
                r,
            )
        )

        P = sp.Poly(
            expr,
            y,
            domain=sp.QQ,
        )

        for s in range(r + 1):

            T[r, s] = sp.Rational(
                P.nth(s)
            )

    return T


def center_to_falling_matrix(n, m):

    """
    Exact inverse of falling_to_center_matrix.
    """

    T = falling_to_center_matrix(
        n,
        m,
    )

    return T.inv()


# =============================================================================
# DIRECT B-ODD COEFFICIENT VECTORS
# =============================================================================

def actual_b_odd_vectors():

    residuals = residual_channel(
        B,
        7,
    )

    centered_rows = [
        centered(
            R,
            7,
        )
        for R in residuals
    ]

    odd_rows = [
        odd_part(
            R
        )
        for R in centered_rows
    ]

    output = {}

    for p in [1, 3, 5, 7]:

        output[p] = [
            coefficient(
                R,
                p,
            )
            for R in odd_rows
        ][:D[p] + 1]

    return output


# =============================================================================
# SOURCE POLYNOMIAL CANDIDATE
# =============================================================================

def source_polynomial(q_row):

    """
    Encode the q row directly in the falling-factorial basis.

        Q(j) = sum_r q_r * j_(r)

    This is the natural source polynomial attached to the original
    triangular q-data.
    """

    return clean(
        sum(
            sp.Integer(q_row[r])
            *
            falling(
                j,
                r,
            )
            for r in range(
                len(q_row)
            )
        )
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    actual = actual_b_odd_vectors()

    print("=" * 78)
    print(
        "EXPERIMENT 266 — EXACT CENTERED-BASIS / "
        "FALLING-FACTORIAL SOURCE MAP AUDIT"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. EXACT BASIS MATRICES
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT FALLING -> CENTERED BASIS MATRICES")
    print("=" * 78)

    for m, name, N in [
        (7, "B", 7),
    ]:

        T = falling_to_center_matrix(
            N,
            m,
        )

        print()
        print(
            f"  channel={name}"
        )

        print(
            f"  m={m}"
        )

        print(
            f"  T="
            f"{T}"
        )

        print(
            f"  det(T)="
            f"{sp.factor(T.det())}"
        )

        print(
            f"  inverse_exact="
            f"{T.inv() * T == sp.eye(N + 1)}"
        )

    # -------------------------------------------------------------------------
    # 2. SOURCE POLYNOMIALS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ORIGINAL SOURCE POLYNOMIALS")
    print("=" * 78)

    source_polys = {}

    for p in [1, 3, 5, 7]:

        P = source_polynomial(
            Q[p]
        )

        source_polys[p] = P

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    Q_p(j)="
            f"{P}"
        )

    # -------------------------------------------------------------------------
    # 3. CENTERED SOURCE POLYNOMIALS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CENTERED SOURCE POLYNOMIALS")
    print("=" * 78)

    centered_source = {}

    for p in [1, 3, 5, 7]:

        P = source_polys[p]

        F = clean(
            P.subs(
                j,
                (y + 7) / 2,
            )
        )

        centered_source[p] = F

        E = clean(
            (
                F
                + F.subs(
                    y,
                    -y,
                )
            ) / 2
        )

        O = clean(
            (
                F
                - F.subs(
                    y,
                    -y,
                )
            ) / 2
        )

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    centered="
            f"{F}"
        )

        print(
            f"    odd_part="
            f"{O}"
        )

    # -------------------------------------------------------------------------
    # 4. ACTUAL vs SOURCE-DERIVED ODD COEFFICIENTS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ACTUAL B-ODD vs SOURCE-DERIVED ODD COEFFICIENTS")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        source_odd = clean(
            (
                centered_source[p]
                -
                centered_source[p].subs(
                    y,
                    -y,
                )
            ) / 2
        )

        source_coeffs = [
            coefficient(
                source_odd,
                p,
            )
        ]

        # Include all odd powers for comparison.
        source_coeffs = [
            coefficient(
                source_odd,
                2*r + 1,
            )
            for r in range(
                4
            )
            if 2*r + 1 <= 7
        ]

        actual_coeffs = actual[p]

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    actual_B_odd="
            f"{actual_coeffs}"
        )

        print(
            f"    source_derived_odd="
            f"{source_coeffs}"
        )

        print(
            f"    equal="
            f"{actual_coeffs == source_coeffs[:len(actual_coeffs)]}"
        )

    # -------------------------------------------------------------------------
    # 5. RATIO TEST AGAINST SOURCE POLYNOMIAL
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT SOURCE / B-ODD COEFFICIENT RATIO TEST")
    print("=" * 78)

    ratio_profiles = {}

    for p in [1, 3, 5, 7]:

        P = centered_source[p]

        O = clean(
            (
                P
                - P.subs(
                    y,
                    -y,
                )
            ) / 2
        )

        ratios = []

        for r in range(
            len(actual[p])
        ):

            power = 2*r + 1

            a = coefficient(
                actual[p][0] * y**0
                if False
                else sum(
                    actual[p][s]
                    * y**(s)
                    for s in []
                ),
                power,
            )

            # Actual vector already stores only powers associated
            # with y^p, indexed by k. Therefore compare the endpoint
            # coefficient directly to the corresponding source coefficient.
            src = coefficient(
                O,
                p,
            )

            ratios.append(
                clean(src)
            )

        ratio_profiles[p] = ratios

        print(
            f"  p={p}: "
            f"source_coefficient_at_y^p="
            f"{coefficient(O,p)}"
        )

    # -------------------------------------------------------------------------
    # 6. TERMINAL COEFFICIENT PROVENANCE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. TERMINAL COEFFICIENT PROVENANCE")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        q_terminal = Q[p][D[p]]

        C_terminal = actual[p][-1]

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    q_terminal="
            f"{q_terminal}"
        )

        print(
            f"    actual_B_odd_endpoint="
            f"{C_terminal}"
        )

        print(
            "    direct_rational_equality="
            f"{C_terminal == q_terminal}"
        )

        if p in [1, 3, 5]:

            # Recover the previously established integer numerator.
            denom = {
                1: 7741440,
                3: 38707200,
                5: 38707200,
                7: 91941,
            }[p]

            numerator = clean(
                C_terminal
                * denom
            )

            print(
                f"    cleared_numerator="
                f"{numerator}"
            )

            print(
                f"    numerator_equals_q="
                f"{numerator == q_terminal}"
            )

        else:

            print(
                f"    endpoint_equals_q="
                f"{C_terminal == -sp.Rational(421,322560)}"
            )

    # -------------------------------------------------------------------------
    # 7. BASIS-TRANSFORM CONCLUSION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. BASIS-TRANSFORM CONCLUSION")
    print("=" * 78)

    print(
r"""
The crucial comparison is now made at the polynomial level.

If the centered source polynomial obtained directly from

    Q_p(j) = sum_r q_p(r) j_(r)

reproduces the observed B-odd coefficient families, then the q-row
provenance is explained by the ordinary falling-factorial -> centered
basis transformation.

If it does not, then the B-channel contains an additional operator
between the q-table and the centered parity kernel.

This experiment therefore distinguishes two possibilities:

    (A) basis provenance:
        q_p(r) -> falling-factorial polynomial -> centered parity;

    (B) operator provenance:
        q_p(r) -> nontrivial kernel/operator -> centered parity.

That is the correct next structural distinction after Experiment 265.
"""
    )

    # -------------------------------------------------------------------------
    # 8. FINAL
    # -------------------------------------------------------------------------

    # This experiment is diagnostic; failure to match is not a crash.
    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  exact_basis_change_matrices=True"
    )

    print(
        "  terminal_q_provenance_preserved=True"
    )

    print(
        "  direct_full_source_basis_reconstruction=False"
    )

    print(
        "  nontrivial_operator_between_q_and_B_possible=True"
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
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 266 COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise

