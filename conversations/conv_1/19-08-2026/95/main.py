#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 260 — EXACT SOURCE-PROVENANCE RECONSTRUCTION FROM
                  CENTERED PARITY POLYNOMIALS
==============================================================================

Goal:

Recover the two distinguished source integers

    q1(5) = 495451247
    q3(4) = 421514439

directly from the exact centered parity data of Experiment 104.

Then verify:

    gcd(q1(5), q3(4)) = 17

and therefore

    q1(5)/17 = 29144191
    q3(4)/17 = 24794967.

The experiment also determines whether these numbers occupy a
structurally distinguished position in the B-odd coefficient family.

No floating point.
No numerical fitting.
No other Python files.
Only main.py.
"""

from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

j, y, k = sp.symbols("j y k")


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
# EXACT HELPERS
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
    num = clean(num)
    den = clean(den)

    pn = poly(num, var)
    pd = poly(den, var)

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
    out = sp.Integer(0)

    for offset, coeff in enumerate(row):
        r = k_index + offset
        out += (
            sp.sympify(coeff)
            * falling(j, r)
        )

    return clean(out)


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
                f"Residual extraction failed "
                f"at k={k_index}."
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
        (F + F.subs(y, -y)) / 2
    )

    odd = clean(
        (F - F.subs(y, -y)) / 2
    )

    return even, odd


def coefficient(expr, power):
    return sp.Rational(
        poly(
            expr,
            y,
        ).nth(power)
    ) if expr != 0 else sp.Integer(0)


def interpolate_index(values):
    pts = [
        (
            sp.Integer(i),
            sp.Rational(v),
        )
        for i, v in enumerate(values)
    ]

    return clean(
        sp.interpolate(
            pts,
            k,
        )
    )


def primitive_integer_signature(values):
    vals = [
        sp.Rational(v)
        for v in values
    ]

    if not vals:
        return []

    denominator_lcm = 1

    for v in vals:
        denominator_lcm = math.lcm(
            denominator_lcm,
            int(v.q),
        )

    ints = [
        int(v * denominator_lcm)
        for v in vals
    ]

    g = 0

    for z in ints:
        g = math.gcd(
            g,
            abs(z),
        )

    if g:
        ints = [
            z // g
            for z in ints
        ]

    for z in ints:
        if z != 0:
            if z < 0:
                ints = [
                    -q
                    for q in ints
                ]
            break

    return ints


# =============================================================================
# SOURCE-PROVENANCE EXTRACTION
# =============================================================================

def build_centered_parity_data():

    A_residuals = residual_channel(
        A,
        8,
    )

    B_residuals = residual_channel(
        B,
        7,
    )

    A_centered = [
        center(R, 8)
        for R in A_residuals
    ]

    B_centered = [
        center(R, 7)
        for R in B_residuals
    ]

    A_parity = [
        parity_parts(F)
        for F in A_centered
    ]

    B_parity = [
        parity_parts(F)
        for F in B_centered
    ]

    return {
        "A_residuals": A_residuals,
        "B_residuals": B_residuals,
        "A_even": [
            pair[0]
            for pair in A_parity
        ],
        "A_odd": [
            pair[1]
            for pair in A_parity
        ],
        "B_even": [
            pair[0]
            for pair in B_parity
        ],
        "B_odd": [
            pair[1]
            for pair in B_parity
        ],
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    data = build_centered_parity_data()

    B_odd = data["B_odd"]

    # -------------------------------------------------------------------------
    # 1. EXTRACT B-ODD y^1 AND y^3 COEFFICIENT FAMILIES
    # -------------------------------------------------------------------------

    print("=" * 78)
    print(
        "EXPERIMENT 260 — EXACT SOURCE-PROVENANCE RECONSTRUCTION "
        "FROM CENTERED PARITY POLYNOMIALS"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print("1. B-ODD SOURCE COEFFICIENT FAMILIES")
    print("=" * 78)

    powers = [1, 3, 5, 7]

    coefficient_families = {}

    for power in powers:

        values = [
            coefficient(
                R,
                power,
            )
            for R in B_odd
        ]

        coefficient_families[power] = values

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

    for power, values in coefficient_families.items():

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

        print(
            f"    degree={poly(Pk, k).degree()}"
        )

    # -------------------------------------------------------------------------
    # 3. RECOVER DISTINGUISHED NUMBERS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. TERMINAL SOURCE RECOVERY")
    print("=" * 78)

    P_y1 = index_polynomials[1]
    P_y3 = index_polynomials[3]

    lc_y1 = sp.LC(
        poly(
            P_y1,
            k,
        )
    )

    # P_y3 has an explicit (k-5) factor.
    factored_y3 = sp.factor(
        P_y3
    )

    inner_y3 = sp.cancel(
        P_y3 / (k - 5)
    )

    lc_y3 = sp.LC(
        poly(
            inner_y3,
            k,
        )
    )

    print(
        f"  B-odd y^1 leading_coefficient="
        f"{lc_y1}"
    )

    print(
        f"  B-odd y^3 polynomial="
        f"{factored_y3}"
    )

    print(
        f"  B-odd y^3/(k-5) leading_coefficient="
        f"{lc_y3}"
    )

    q1_terminal = int(
        lc_y1
    )

    q3_terminal = int(
        lc_y3
    )

    expected_q1 = 495451247
    expected_q3 = 421514439

    q1_exact = (
        q1_terminal
        ==
        expected_q1
    )

    q3_exact = (
        q3_terminal
        ==
        expected_q3
    )

    print(
        f"  recovered_q1_5="
        f"{q1_terminal}"
    )

    print(
        f"  recovered_q3_4="
        f"{q3_terminal}"
    )

    print(
        f"  q1_recovery_exact="
        f"{q1_exact}"
    )

    print(
        f"  q3_recovery_exact="
        f"{q3_exact}"
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
    # 4. 17-ADIC CONTENT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT 17-CONTENT REMOVAL")
    print("=" * 78)

    g = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    print(
        f"  gcd(q1_5,q3_4)={g}"
    )

    v17_q1 = 0
    x = abs(q1_terminal)

    while x % 17 == 0:
        v17_q1 += 1
        x //= 17

    v17_q3 = 0
    x = abs(q3_terminal)

    while x % 17 == 0:
        v17_q3 += 1
        x //= 17

    print(
        f"  v17(q1_5)={v17_q1}"
    )

    print(
        f"  v17(q3_4)={v17_q3}"
    )

    q1_reduced = (
        q1_terminal // 17
    )

    q3_reduced = (
        q3_terminal // 17
    )

    print(
        f"  q1_5/17={q1_reduced}"
    )

    print(
        f"  q3_4/17={q3_reduced}"
    )

    print(
        f"  reduced_q1_expected=29144191"
    )

    print(
        f"  reduced_q3_expected=24794967"
    )

    content_exact = (
        g == 17
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
    # 5. CHECK THE SOURCE NUMBERS ARE NOT RANDOM COEFFICIENTS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SOURCE-COEFFICIENT STRUCTURAL CHECK")
    print("=" * 78)

    print(
        """
The recovered integers occupy two different layers of the same
B-odd coefficient family:

    y^1:
        leading coefficient = q1(5)

    y^3:
        after removal of its terminal zero factor (k-5),
        leading coefficient = q3(4).

This is a provenance relation inside the centered parity kernel,
rather than an externally inserted pair of integers.
"""
    )

    structural_checks = {
        "q1_from_B_odd_y1":
            q1_exact,

        "q3_from_B_odd_y3":
            q3_exact,

        "17_content_common":
            g == 17,

        "q1_17_primitive":
            q1_reduced % 17 != 0,

        "q3_17_primitive":
            q3_reduced % 17 != 0,
    }

    for name, value in structural_checks.items():

        print(
            f"  {name}={value}"
        )

        if not value:
            failures.append(
                name
            )

    # -------------------------------------------------------------------------
    # 6. FULL B-ODD LEADING-COEFFICIENT PROFILE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. B-ODD LEADING-COEFFICIENT PROFILE")
    print("=" * 78)

    for power in powers:

        Pk = index_polynomials[power]

        degree = int(
            poly(
                Pk,
                k,
            ).degree()
        )

        if degree < 0:
            continue

        lc = sp.LC(
            poly(
                Pk,
                k,
            )
        )

        print(
            f"  y^{power}: "
            f"degree={degree} "
            f"leading_coefficient={lc}"
        )

    # -------------------------------------------------------------------------
    # 7. COMPARE WITH EXPERIMENT 184 PROVENANCE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. CROSS-CHECK AGAINST THE 17-NORMALIZED SOURCE PAIR")
    print("=" * 78)

    print(
        f"  original_q1_5={q1_terminal}"
    )

    print(
        f"  original_q3_4={q3_terminal}"
    )

    print(
        f"  normalized_q1_5={q1_reduced}"
    )

    print(
        f"  normalized_q3_4={q3_reduced}"
    )

    print(
        f"  normalized_pair="
        f"({q1_reduced},{q3_reduced})"
    )

    print(
        f"  expected_projective_pair="
        f"(29144191,24794967)"
    )

    normalized_pair_exact = (
        (
            q1_reduced,
            q3_reduced,
        )
        ==
        (
            29144191,
            24794967,
        )
    )

    print(
        f"  normalized_pair_exact="
        f"{normalized_pair_exact}"
    )

    if not normalized_pair_exact:
        failures.append(
            "normalized_pair"
        )

    # -------------------------------------------------------------------------
    # 8. WHAT THIS ACTUALLY PROVES
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        r"""
This experiment establishes a direct provenance chain inside the exact
centered parity kernel:

    B-odd y^1 coefficient
        ->
    q1(5)=495451247

    B-odd y^3 coefficient
        ->
    q3(4)=421514439

followed by the exact common-content removal

    (q1(5),q3(4))/17
        =
    (29144191,24794967).

Therefore the projective source pair used in the later 7-adic experiments
is not an arbitrary pair of large integers. It is inherited from the
earlier centered parity coefficient structure.

The next unresolved question is now sharper:

    why do these particular B-odd coefficient leading terms equal
    q1(5) and q3(4), and what is the general q_p(r) construction?

That is the source-level problem to attack next.
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
        f"  q1_5_recovered_exactly="
        f"{q1_exact}"
    )

    print(
        f"  q3_4_recovered_exactly="
        f"{q3_exact}"
    )

    print(
        f"  exact_17_content_recovered="
        f"{content_exact}"
    )

    print(
        f"  normalized_pair_exact="
        f"{normalized_pair_exact}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 260 COMPLETE")


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

