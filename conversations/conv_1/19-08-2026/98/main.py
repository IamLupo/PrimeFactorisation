#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 262 — EXACT PARITY-DEGREE / TERMINAL-q INDEX CORRESPONDENCE AUDIT
==============================================================================

Experiment 261 recovered the exact B-odd terminal integer ladder:

    y^1 -> 495451247
    y^3 -> 421514439
    y^5 -> 16027881
    y^7 -> 1

The original q-table has terminal indices

    D(1)=5,
    D(3)=4,
    D(5)=2,
    D(7)=0.

Experiment 262 tests the exact correspondence

    y^p  ->  q_p(D(p)).

For the known source data this should give

    y^1 -> q_1(5) = 495451247
    y^3 -> q_3(4) = 421514439
    y^5 -> q_5(2) = 16027881
    y^7 -> q_7(0) = 1.

This is a provenance audit, not a guessed symbolic formula.

It also verifies that the recovered values occur literally inside the
original q-table and checks the full terminal row

    [q_1(5), q_3(4), q_5(2), q_7(0)].

Only main.py is used.
No floating point.
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
# ORIGINAL q-TABLE FROM EXPERIMENT 184
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

        degree = k_index + offset

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
            c * denominator_lcm
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
        int(numerator.LC()),
        denominator,
        numerator.as_expr(),
    )


def integer_v17(x):

    x = abs(int(x))

    if x == 0:
        return None

    value = 0

    while x % 17 == 0:

        x //= 17
        value += 1

    return value


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    # -------------------------------------------------------------------------
    # 1. RECONSTRUCT B-ODD PARITY SOURCE INTEGERS
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

    recovered = {}

    for power in parity_degrees:

        values = [
            coefficient(
                R,
                power,
            )
            for R in B_odd
        ]

        Pk = interpolate_index(
            values
        )

        # Remove the exact terminal factor predicted by the
        # vanishing pattern of the coefficient family.
        if power == 1:

            structural = Pk

        elif power == 3:

            structural = clean(
                Pk
                / (
                    k - 5
                )
            )

        elif power == 5:

            structural = clean(
                Pk
                / (
                    (k - 5)
                    * (k - 4)
                    * (k - 3)
                )
            )

        elif power == 7:

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

        recovered[power] = {
            "value": value,
            "denominator": denominator,
            "numerator": numerator,
            "polynomial": structural,
        }

    # -------------------------------------------------------------------------
    # OUTPUT
    # -------------------------------------------------------------------------

    print("=" * 78)
    print(
        "EXPERIMENT 262 — EXACT PARITY-DEGREE / "
        "TERMINAL-q INDEX CORRESPONDENCE AUDIT"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. TERMINAL INDEX MAP
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. TERMINAL q-INDEX MAP")
    print("=" * 78)

    print(
        f"  D={D}"
    )

    for p in parity_degrees:

        print(
            f"  p={p}: "
            f"parity_degree=y^{p} "
            f"D(p)={D[p]}"
        )

    # -------------------------------------------------------------------------
    # 2. PARITY-RECOVERED VALUES
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PARITY-RECOVERED TERMINAL VALUES")
    print("=" * 78)

    for p in parity_degrees:

        print()
        print(
            f"  y^{p}:"
        )

        print(
            f"    recovered_integer="
            f"{recovered[p]['value']}"
        )

        print(
            f"    denominator="
            f"{recovered[p]['denominator']}"
        )

        print(
            f"    structural_polynomial="
            f"{recovered[p]['polynomial']}"
        )

    # -------------------------------------------------------------------------
    # 3. DIRECT q-TABLE LOOKUP
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. DIRECT ORIGINAL q-TABLE LOOKUP")
    print("=" * 78)

    lookup_exact = True

    for p in parity_degrees:

        r = D[p]

        q_value = Q_ORIGINAL[p][r]

        recovered_value = recovered[p][
            "value"
        ]

        exact = (
            q_value
            ==
            recovered_value
        )

        lookup_exact &= exact

        print(
            f"  p={p}: "
            f"q_{p}({r})={q_value} "
            f"parity_recovered={recovered_value} "
            f"exact={exact}"
        )

        if not exact:

            failures.append(
                (
                    "q_lookup",
                    p,
                    r,
                )
            )

    # -------------------------------------------------------------------------
    # 4. FULL TERMINAL LADDER
    # -------------------------------------------------------------------------

    terminal_q_row = [
        Q_ORIGINAL[p][D[p]]
        for p in parity_degrees
    ]

    parity_row = [
        recovered[p]["value"]
        for p in parity_degrees
    ]

    print()
    print("=" * 78)
    print("4. COMPLETE TERMINAL q-LADDER")
    print("=" * 78)

    print(
        f"  terminal_indices="
        f"{[(p, D[p]) for p in parity_degrees]}"
    )

    print(
        f"  q_terminal_row="
        f"{terminal_q_row}"
    )

    print(
        f"  parity_recovered_row="
        f"{parity_row}"
    )

    row_exact = (
        terminal_q_row
        ==
        parity_row
    )

    print(
        f"  complete_row_exact="
        f"{row_exact}"
    )

    if not row_exact:

        failures.append(
            "complete_terminal_row"
        )

    # -------------------------------------------------------------------------
    # 5. FACTOR / 17 PROFILE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. TERMINAL q-LADDER FACTOR PROFILE")
    print("=" * 78)

    for p in parity_degrees:

        value = recovered[p][
            "value"
        ]

        print(
            f"  q_{p}({D[p]})={value}"
        )

        print(
            f"    factorization="
            f"{sp.factorint(abs(value)) if value else {}}"
        )

        print(
            f"    v17={integer_v17(value)}"
        )

    # -------------------------------------------------------------------------
    # 6. PROJECTIVE NORMALIZATION OF THE FIRST TWO TERMINAL VALUES
    # -------------------------------------------------------------------------

    q1_terminal = Q_ORIGINAL[1][5]
    q3_terminal = Q_ORIGINAL[3][4]

    gcd_terminal = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    normalized_q1 = (
        q1_terminal // 17
    )

    normalized_q3 = (
        q3_terminal // 17
    )

    normalized_exact = (
        gcd_terminal == 17
        and
        normalized_q1 == 29144191
        and
        normalized_q3 == 24794967
    )

    print()
    print("=" * 78)
    print("6. PROJECTIVE SOURCE NORMALIZATION")
    print("=" * 78)

    print(
        f"  q_1(5)={q1_terminal}"
    )

    print(
        f"  q_3(4)={q3_terminal}"
    )

    print(
        f"  gcd={gcd_terminal}"
    )

    print(
        f"  q_1(5)/17={normalized_q1}"
    )

    print(
        f"  q_3(4)/17={normalized_q3}"
    )

    print(
        f"  normalized_pair_exact="
        f"{normalized_exact}"
    )

    if not normalized_exact:

        failures.append(
            "projective_normalization"
        )

    # -------------------------------------------------------------------------
    # 7. INDEX/PARITY CORRESPONDENCE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PARITY DEGREE / TERMINAL INDEX CORRESPONDENCE")
    print("=" * 78)

    correspondence = True

    for p in parity_degrees:

        recovered_value = recovered[p][
            "value"
        ]

        expected_value = Q_ORIGINAL[p][
            D[p]
        ]

        exact = (
            recovered_value
            ==
            expected_value
        )

        correspondence &= exact

        print(
            f"  y^{p}"
            f" -> "
            f"q_{p}({D[p]})"
            f" -> "
            f"{expected_value}"
            f" "
            f"exact={exact}"
        )

    print(
        f"  correspondence_exact="
        f"{correspondence}"
    )

    if not correspondence:

        failures.append(
            "parity_index_correspondence"
        )

    # -------------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
The exact four-term source ladder is now identified as

    y^1 -> q_1(D(1)) = q_1(5)
    y^3 -> q_3(D(3)) = q_3(4)
    y^5 -> q_5(D(5)) = q_5(2)
    y^7 -> q_7(D(7)) = q_7(0).

For the known source table this gives

    [495451247,
     421514439,
     16027881,
     1].

Thus the parity degree is not merely a descriptive label. It exactly
matches the p-index of the terminal q-family, while the terminal
vanishing factor determines the r-index D(p).

In compact form, the observed correspondence is

    B-odd y^p
       ->
    q_p(D(p)).

This is a substantially stronger provenance statement than simply
observing two matching integers.

It does NOT yet prove a universal symbolic formula for q_p(r) or for
D(p). It establishes the exact source/index interface for the known
kernel.
"""
    )

    # -------------------------------------------------------------------------
    # 9. FUTURE FAMILY INTERFACE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FUTURE n=pq SOURCE INTERFACE")
    print("=" * 78)

    print(
        """
The next family-level question is now precise.

For another genuine n=pq construction, determine whether the
corresponding centered parity kernel again produces

    B-odd y^p
       =
    q_p(D(p))

for the same terminal-index map.

A second genuine source case should therefore record, for every
terminal odd p,

    p,
    D(p),
    q_p(D(p)).

Only then should we attempt a symbolic p,q formula.
"""
    )

    # -------------------------------------------------------------------------
    # FINAL
    # -------------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  complete_terminal_row_exact="
        f"{row_exact}"
    )

    print(
        f"  parity_index_correspondence_exact="
        f"{correspondence}"
    )

    print(
        f"  projective_normalization_exact="
        f"{normalized_exact}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  source_index_interface_identified=True"
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
        "EXPERIMENT 262 COMPLETE"
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

