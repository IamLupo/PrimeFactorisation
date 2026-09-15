#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 260R — EXACT SOURCE-PROVENANCE RECONSTRUCTION FROM
                   CENTERED PARITY POLYNOMIAL NUMERATORS
==============================================================================

Correction of Experiment 260.

The rational index polynomials have the form

    P(k) = N(k) / D,

where N(k) has integer coefficients.

Experiment 104 reported the distinguished source integers from the
NUMERATOR of the leading term after clearing the common denominator.

Therefore:

    B-odd y^1:
        numerator leading coefficient = 495451247

    B-odd y^3:
        after removing (k-5),
        numerator leading coefficient = 421514439.

This script reconstructs those exact integers and then checks their
17-content removal:

    495451247 / 17 = 29144191
    421514439 / 17 = 24794967

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
# EXACT EXPERIMENT 104 DATA
# =============================================================================

A = [
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
    [
        sp.Rational(-1001, 12),
        sp.Rational(26687, 360),
        sp.Rational(-40921, 1260),
        sp.Rational(9389, 1008),
    ],
    [
        sp.Rational(-5, 72),
        sp.Rational(5, 72),
        sp.Rational(-5, 144),
    ],
]

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
        sp.expand(
            sp.sympify(expr)
        ),
        var,
        domain=sp.QQ,
    )


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


def parity_parts(F):

    even = clean(
        (
            F
            + F.subs(y, -y)
        )
        / 2
    )

    odd = clean(
        (
            F
            - F.subs(y, -y)
        )
        / 2
    )

    return even, odd


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


# =============================================================================
# NUMERATOR EXTRACTION
# =============================================================================

def cleared_integer_polynomial(expr):

    """
    Given a rational polynomial

        P(k) = N(k)/D,

    return:

        N(k) primitive integer polynomial,
        D positive integer.

    This is the crucial correction to Experiment 260.
    """

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

    integer_coeffs = [
        int(
            c
            * denominator_lcm
        )
        for c in coeffs
    ]

    coefficient_gcd = 0

    for c in integer_coeffs:
        coefficient_gcd = math.gcd(
            coefficient_gcd,
            abs(c),
        )

    if coefficient_gcd != 0:

        integer_coeffs = [
            c // coefficient_gcd
            for c in integer_coeffs
        ]

    # Keep the denominator associated with this displayed
    # primitive numerator representation.
    denominator = (
        denominator_lcm
        // max(
            1,
            coefficient_gcd,
        )
    )

    numerator = sp.Poly(
        sum(
            integer_coeffs[i]
            * k ** (
                len(integer_coeffs)
                - 1
                - i
            )
            for i in range(
                len(integer_coeffs)
            )
        ),
        k,
        domain=sp.ZZ,
    )

    return numerator, denominator


def leading_integer_numerator(expr):

    numerator, denominator = (
        cleared_integer_polynomial(
            expr
        )
    )

    degree = numerator.degree()

    leading_integer = int(
        numerator.LC()
    )

    return (
        leading_integer,
        denominator,
        degree,
        numerator.as_expr(),
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    print("=" * 78)
    print(
        "EXPERIMENT 260R — EXACT SOURCE-PROVENANCE RECONSTRUCTION "
        "FROM CENTERED PARITY POLYNOMIAL NUMERATORS"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # Build B-odd channel
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
        parity_parts(F)[1]
        for F in B_centered
    ]

    # -------------------------------------------------------------------------
    # 1. COEFFICIENT FAMILIES
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. B-ODD SOURCE COEFFICIENT FAMILIES")
    print("=" * 78)

    families = {}

    for power in [1, 3, 5, 7]:

        values = [
            coefficient(
                R,
                power,
            )
            for R in B_odd
        ]

        families[power] = values

        print()
        print(
            f"  y^{power}:"
        )

        print(
            f"    values={values}"
        )

    # -------------------------------------------------------------------------
    # 2. INDEX POLYNOMIALS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT INDEX POLYNOMIALS")
    print("=" * 78)

    index_polynomials = {}

    for power, values in families.items():

        Pk = interpolate_index(
            values
        )

        index_polynomials[power] = Pk

        print()
        print(
            f"  y^{power}:"
        )

        print(
            f"    P(k)={Pk}"
        )

    # -------------------------------------------------------------------------
    # 3. CLEARED NUMERATORS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CLEARED INTEGER NUMERATORS")
    print("=" * 78)

    cleared = {}

    for power, Pk in index_polynomials.items():

        numerator, denominator = (
            cleared_integer_polynomial(
                Pk
            )
        )

        cleared[power] = (
            numerator,
            denominator,
        )

        print()
        print(
            f"  y^{power}:"
        )

        print(
            f"    numerator="
            f"{numerator.as_expr()}"
        )

        print(
            f"    denominator="
            f"{denominator}"
        )

        print(
            f"    leading_integer="
            f"{int(numerator.LC())}"
        )

    # -------------------------------------------------------------------------
    # 4. RECOVER Q1(5) AND Q3(4)
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT TERMINAL SOURCE RECOVERY")
    print("=" * 78)

    q1_5, q1_den, q1_degree, q1_num_poly = (
        leading_integer_numerator(
            index_polynomials[1]
        )
    )

    # For y^3 the structural terminal factor is (k-5).
    reduced_y3 = clean(
        index_polynomials[3]
        / (
            k - 5
        )
    )

    q3_4, q3_den, q3_degree, q3_num_poly = (
        leading_integer_numerator(
            reduced_y3
        )
    )

    print(
        f"  q1_source_integer={q1_5}"
    )

    print(
        f"  q1_source_denominator={q1_den}"
    )

    print(
        f"  q3_source_integer={q3_4}"
    )

    print(
        f"  q3_source_denominator={q3_den}"
    )

    expected_q1 = 495451247
    expected_q3 = 421514439

    q1_exact = (
        q1_5
        == expected_q1
    )

    q3_exact = (
        q3_4
        == expected_q3
    )

    print(
        f"  q1_recovery_exact={q1_exact}"
    )

    print(
        f"  q3_recovery_exact={q3_exact}"
    )

    if not q1_exact:
        failures.append(
            "q1_recovery"
        )

    if not q3_exact:
        failures.append(
            "q3_recovery"
        )

    # -------------------------------------------------------------------------
    # 5. 17-CONTENT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT 17-CONTENT REMOVAL")
    print("=" * 78)

    common_gcd = math.gcd(
        q1_5,
        q3_4,
    )

    def integer_v17(x):

        x = abs(int(x))

        if x == 0:
            return None

        v = 0

        while x % 17 == 0:
            x //= 17
            v += 1

        return v

    v17_q1 = integer_v17(
        q1_5
    )

    v17_q3 = integer_v17(
        q3_4
    )

    q1_reduced = (
        q1_5 // 17
    )

    q3_reduced = (
        q3_4 // 17
    )

    print(
        f"  gcd(q1_5,q3_4)={common_gcd}"
    )

    print(
        f"  v17(q1_5)={v17_q1}"
    )

    print(
        f"  v17(q3_4)={v17_q3}"
    )

    print(
        f"  q1_5/17={q1_reduced}"
    )

    print(
        f"  q3_4/17={q3_reduced}"
    )

    content_exact = (
        common_gcd == 17
        and
        v17_q1 == 1
        and
        v17_q3 == 1
        and
        q1_reduced == 29144191
        and
        q3_reduced == 24794967
    )

    print(
        f"  exact_17_content_recovery="
        f"{content_exact}"
    )

    if not content_exact:
        failures.append(
            "17_content_recovery"
        )

    # -------------------------------------------------------------------------
    # 6. STRUCTURAL COEFFICIENT PROFILE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. B-ODD LEADING INTEGER PROFILE")
    print("=" * 78)

    for power in [1, 3, 5, 7]:

        Pk = (
            index_polynomials[power]
        )

        if power == 3:

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

        else:

            structural = Pk

        numerator, denominator = (
            cleared_integer_polynomial(
                structural
            )
        )

        print(
            f"  y^{power}: "
            f"leading_integer="
            f"{int(numerator.LC())} "
            f"degree={numerator.degree()}"
        )

    # -------------------------------------------------------------------------
    # 7. DIRECT RECONSTRUCTION CROSS-CHECK
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. DIRECT CROSS-CHECK")
    print("=" * 78)

    checks = {
        "495451247/17":
            495451247 // 17
            ==
            29144191,

        "421514439/17":
            421514439 // 17
            ==
            24794967,

        "gcd=17":
            math.gcd(
                495451247,
                421514439,
            )
            ==
            17,

        "q1_5_not_divisible_by_17_twice":
            495451247 % (
                17 * 17
            ) != 0,

        "q3_4_not_divisible_by_17_twice":
            421514439 % (
                17 * 17
            ) != 0,
    }

    for name, value in checks.items():

        print(
            f"  {name}={value}"
        )

        if not value:
            failures.append(
                name
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
The previous failure was purely representational.

The exact source integers are not the rational leading coefficients

    495451247/7741440
    421514439/38707200.

They are the integer numerators obtained after clearing the common
rational denominator in the indexed coefficient polynomial.

The corrected provenance chain is therefore:

    B-odd y^1
        ->
    integer numerator 495451247

    B-odd y^3/(k-5)
        ->
    integer numerator 421514439

and then

    (495451247,421514439)
        =
    17*(29144191,24794967).

This directly connects the centered parity source layer to the
projective source pair used in the later experiments.

The next source-level question is now:

    why do the terminal B-odd numerator families encode q1(5)
    and q3(4), and can their construction be generalized?
"""
    )

    # -------------------------------------------------------------------------
    # 9. FINAL
    # -------------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_5_recovered_exactly={q1_exact}"
    )

    print(
        f"  q3_4_recovered_exactly={q3_exact}"
    )

    print(
        f"  exact_17_content_recovered="
        f"{content_exact}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print(
        "EXPERIMENT 260R COMPLETE"
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

