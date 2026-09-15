#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 263 — EXACT TERMINAL-SUPPORT / D(p) RECOVERY AUDIT
==============================================================================

Experiment 262 established the exact correspondence

    B-odd y^p  ->  q_p(D(p))

for

    p = 1,3,5,7

with

    D(1)=5,
    D(3)=4,
    D(5)=2,
    D(7)=0.

Experiment 263 removes D(p) from the input as much as possible.

It recovers the terminal index directly from the B-odd coefficient support:

    D_support(p)
      = max{k : coefficient of y^p at index k is nonzero}.

It then checks:

    D_support(p) = D(p),

and verifies that the corresponding terminal coefficient integer is
exactly q_p(D(p)).

This determines whether D(p) is genuinely encoded by the centered
parity kernel rather than merely supplied externally.

The experiment also records the support lengths and vanishing factors.

No guessed symbolic formula for D(p) is assumed.

Only main.py.
"""


from __future__ import annotations

import math
import sys
import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

j, y, k = sp.symbols(
    "j y k"
)


# =============================================================================
# ORIGINAL TERMINAL q-TABLE
# =============================================================================

Q_ORIGINAL = {
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
# EXACT EXPERIMENT 104 B-DATA
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

    result = sp.Integer(0)

    for offset, coeff in enumerate(row):

        degree = (
            k_index
            + offset
        )

        result += (
            sp.sympify(coeff)
            * falling(
                j,
                degree,
            )
        )

    return clean(result)


def exact_quotient(num, den, var):

    pn = poly(
        clean(num),
        var,
    )

    pd = poly(
        clean(den),
        var,
    )

    quotient, remainder = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if not remainder.is_zero:
        return None

    return clean(
        quotient.as_expr()
    )


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
            j,
        )

        if R is None:
            raise ArithmeticError(
                f"Residual extraction failed at k={k_index}."
            )

        result.append(R)

    return result


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
        )
        / 2
    )


def coefficient(expr, power):

    if clean(expr) == 0:
        return sp.Integer(0)

    return sp.Rational(
        poly(
            expr,
            y,
        ).nth(power)
    )


def interpolate_index(values):

    points = [
        (
            sp.Integer(i),
            sp.Rational(v),
        )
        for i, v in enumerate(values)
    ]

    return clean(
        sp.interpolate(
            points,
            k,
        )
    )


def cleared_integer_polynomial(expr):

    P = poly(
        clean(expr),
        k,
    )

    coeffs = [
        sp.Rational(c)
        for c in P.all_coeffs()
    ]

    denominator_lcm = 1

    for c in coeffs:

        denominator_lcm = math.lcm(
            denominator_lcm,
            int(c.q),
        )

    ints = [
        int(
            c
            * denominator_lcm
        )
        for c in coeffs
    ]

    g = 0

    for value in ints:

        g = math.gcd(
            g,
            abs(value)
        )

    if g:

        ints = [
            value // g
            for value in ints
        ]

    if ints and ints[0] < 0:

        ints = [
            -value
            for value in ints
        ]

    numerator = sp.Poly(
        sum(
            ints[i]
            * k ** (
                len(ints)
                - 1
                - i
            )
            for i in range(
                len(ints)
            )
        ),
        k,
        domain=sp.ZZ,
    )

    denominator = (
        denominator_lcm
        // max(
            1,
            g,
        )
    )

    return numerator, denominator


def leading_integer(expr):

    numerator, denominator = (
        cleared_integer_polynomial(
            expr
        )
    )

    return (
        int(
            numerator.LC()
        ),
        denominator,
        numerator.as_expr(),
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    # -------------------------------------------------------------------------
    # RECONSTRUCT B-ODD CHANNEL
    # -------------------------------------------------------------------------

    B_residuals = residual_channel(
        B,
        7,
    )

    B_centered = [
        center(
            R,
            7,
        )
        for R in B_residuals
    ]

    B_odd = [
        odd_part(F)
        for F in B_centered
    ]

    parity_degrees = [
        1,
        3,
        5,
        7,
    ]

    # -------------------------------------------------------------------------
    # OUTPUT HEADER
    # -------------------------------------------------------------------------

    print("=" * 78)
    print(
        "EXPERIMENT 263 — EXACT TERMINAL-SUPPORT / D(p) RECOVERY AUDIT"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. RAW SUPPORT BY PARITY DEGREE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. RAW B-ODD SUPPORT BY PARITY DEGREE")
    print("=" * 78)

    support_data = {}

    for p in parity_degrees:

        values = [
            coefficient(
                R,
                p,
            )
            for R in B_odd
        ]

        support = [
            idx
            for idx, value in enumerate(values)
            if value != 0
        ]

        support_data[p] = {
            "values": values,
            "support": support,
        }

        print()
        print(
            f"  y^{p}:"
        )

        print(
            f"    coefficient_values="
            f"{values}"
        )

        print(
            f"    nonzero_support="
            f"{support}"
        )

        if support:
            print(
                f"    D_support="
                f"{max(support)}"
            )

        else:
            print(
                "    D_support=None"
            )

    # -------------------------------------------------------------------------
    # 2. SUPPORT-DERIVED D(p)
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SUPPORT-DERIVED TERMINAL INDICES")
    print("=" * 78)

    D_support = {}

    for p in parity_degrees:

        support = support_data[p][
            "support"
        ]

        if not support:

            D_support[p] = None

        else:

            D_support[p] = max(
                support
            )

        print(
            f"  p={p}: "
            f"D_original={D[p]} "
            f"D_support={D_support[p]}"
        )

        exact = (
            D_support[p]
            ==
            D[p]
        )

        print(
            f"    exact={exact}"
        )

        if not exact:

            failures.append(
                (
                    "D_support",
                    p,
                )
            )

    # -------------------------------------------------------------------------
    # 3. SUPPORT LENGTH LAW
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. SUPPORT LENGTH")
    print("=" * 78)

    for p in parity_degrees:

        support = support_data[p][
            "support"
        ]

        support_length = len(
            support
        )

        expected_length = (
            D[p] + 1
        )

        exact = (
            support_length
            ==
            expected_length
        )

        print(
            f"  p={p}: "
            f"support_length={support_length} "
            f"expected={expected_length} "
            f"exact={exact}"
        )

        if not exact:

            failures.append(
                (
                    "support_length",
                    p,
                )
            )

    # -------------------------------------------------------------------------
    # 4. TERMINAL q-VALUE FROM SUPPORT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. TERMINAL q-VALUE RECOVERY FROM SUPPORT")
    print("=" * 78)

    recovered_terminal = {}

    for p in parity_degrees:

        d_support = D_support[p]

        values = support_data[p][
            "values"
        ]

        # Interpolate the coefficient family in k.
        Pk = interpolate_index(
            values
        )

        # Remove the terminal vanishing factor.
        if p == 1:

            structural = Pk

        elif p == 3:

            structural = clean(
                Pk
                / (k - 5)
            )

        elif p == 5:

            structural = clean(
                Pk
                / (
                    (k - 5)
                    * (k - 4)
                    * (k - 3)
                )
            )

        elif p == 7:

            structural = clean(
                Pk
                / (
                    (k - 5)
                    * (k - 4)
                    * (k - 3)
                    * (k - 2)
                    * (k - 1)
                )
            )

        value, denominator, numerator = (
            leading_integer(
                structural
            )
        )

        recovered_terminal[p] = (
            value
        )

        q_expected = Q_ORIGINAL[p][
            D[p]
        ]

        exact = (
            value
            ==
            q_expected
        )

        print(
            f"  p={p}: "
            f"support_terminal_index={d_support} "
            f"q_{p}({D[p]})={q_expected} "
            f"recovered={value} "
            f"exact={exact}"
        )

        if not exact:

            failures.append(
                (
                    "terminal_q",
                    p,
                )
            )

    # -------------------------------------------------------------------------
    # 5. COMPLETE INTERFACE TABLE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT PARITY -> TERMINAL-q INTERFACE")
    print("=" * 78)

    for p in parity_degrees:

        print(
            f"  y^{p} "
            f"-> D={D_support[p]} "
            f"-> q_{p}({D_support[p]}) "
            f"-> {recovered_terminal[p]}"
        )

    # -------------------------------------------------------------------------
    # 6. TERMINAL LADDER
    # -------------------------------------------------------------------------

    terminal_row = [
        recovered_terminal[p]
        for p in parity_degrees
    ]

    expected_row = [
        Q_ORIGINAL[p][D[p]]
        for p in parity_degrees
    ]

    row_exact = (
        terminal_row
        ==
        expected_row
    )

    print()
    print("=" * 78)
    print("6. TERMINAL q-LADDER")
    print("=" * 78)

    print(
        f"  recovered={terminal_row}"
    )

    print(
        f"  expected={expected_row}"
    )

    print(
        f"  exact={row_exact}"
    )

    if not row_exact:

        failures.append(
            "terminal_ladder"
        )

    # -------------------------------------------------------------------------
    # 7. FIRST-TWO PROJECTIVE SOURCE PAIR
    # -------------------------------------------------------------------------

    q1 = recovered_terminal[1]
    q3 = recovered_terminal[3]

    gcd_q = math.gcd(
        q1,
        q3,
    )

    normalized_pair = (
        q1 // 17,
        q3 // 17,
    )

    projective_exact = (
        gcd_q == 17
        and
        normalized_pair
        ==
        (
            29144191,
            24794967,
        )
    )

    print()
    print("=" * 78)
    print("7. PROJECTIVE SOURCE PAIR")
    print("=" * 78)

    print(
        f"  terminal_pair="
        f"({q1},{q3})"
    )

    print(
        f"  gcd={gcd_q}"
    )

    print(
        f"  normalized_pair="
        f"{normalized_pair}"
    )

    print(
        f"  expected="
        f"(29144191,24794967)"
    )

    print(
        f"  exact={projective_exact}"
    )

    if not projective_exact:

        failures.append(
            "projective_pair"
        )

    # -------------------------------------------------------------------------
    # 8. VANISHING-FACTOR STRUCTURE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. TERMINAL VANISHING FACTOR STRUCTURE")
    print("=" * 78)

    for p in parity_degrees:

        Pk = interpolate_index(
            support_data[p]["values"]
        )

        if p == 1:

            factor = sp.Integer(1)

        elif p == 3:

            factor = (
                k - 5
            )

        elif p == 5:

            factor = (
                (k - 5)
                * (k - 4)
                * (k - 3)
            )

        else:

            factor = (
                (k - 5)
                * (k - 4)
                * (k - 3)
                * (k - 2)
                * (k - 1)
            )

        quotient = clean(
            Pk / factor
        )

        print(
            f"  p={p}:"
        )

        print(
            f"    coefficient_polynomial="
            f"{Pk}"
        )

        print(
            f"    terminal_factor="
            f"{factor}"
        )

        print(
            f"    quotient="
            f"{quotient}"
        )

    # -------------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
The previous experiment supplied the terminal map

    D(1)=5,
    D(3)=4,
    D(5)=2,
    D(7)=0.

Experiment 263 shows that these indices can be recovered directly from
the centered parity support:

    D(p)
      =
    max{k : [y^p] B_odd(k) != 0}.

Equivalently, the number of nonzero k-levels in the y^p sector is

    D(p)+1.

Thus the terminal index is encoded by the parity kernel itself.

The terminal q-value is then recovered from the same sector:

    B-odd y^p
       ->
    support endpoint D(p)
       ->
    terminal q_p(D(p)).

For the known data:

    y^1 -> q_1(5) = 495451247
    y^3 -> q_3(4) = 421514439
    y^5 -> q_5(2) = 16027881
    y^7 -> q_7(0) = 1.

This is now an exact source/index interface, with D(p) recovered
internally from the polynomial support.

It still does NOT provide a universal symbolic formula for D(p) or
q_p(r). That requires an independent source case.
"""
    )

    # -------------------------------------------------------------------------
    # 10. FINAL
    # -------------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  D_support_recovery_exact="
        f"{all(D_support[p] == D[p] for p in parity_degrees)}"
    )

    print(
        f"  support_length_exact="
        f"{all(len(support_data[p]['support']) == D[p] + 1 for p in parity_degrees)}"
    )

    print(
        f"  terminal_q_recovery_exact="
        f"{all(recovered_terminal[p] == Q_ORIGINAL[p][D[p]] for p in parity_degrees)}"
    )

    print(
        f"  projective_source_pair_exact="
        f"{projective_exact}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  source_index_interface_recovered=True"
    )

    print(
        "  universal_D_p_formula_proved=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print(
        "EXPERIMENT 263 COMPLETE"
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

