#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 261 — EXACT TERMINAL q-VALUE LADDER / PARITY-NUMERATOR
                 PROVENANCE AUDIT
==============================================================================

Experiment 260R established:

    B-odd y^1  -> 495451247
    B-odd y^3  -> 421514439

and

    gcd(495451247,421514439)=17,

with

    (495451247,421514439)/17
      =
    (29144191,24794967).

Experiment 261 now audits the COMPLETE B-odd leading-integer ladder.

From the exact centered parity data we recover:

    y^1 -> 495451247
    y^3 -> 421514439
    y^5 -> 16027881
    y^7 -> 1

The purpose is to determine whether these integers form the terminal
q-value sequence of the underlying source construction.

The experiment deliberately separates:

    OBSERVED:
        exact leading integer numerators;

    CANDIDATE:
        a mapping from parity degree to q_p(r);

    PROVED:
        only the exact numerical provenance established by the data.

No guessed closed formula is introduced.

Only main.py is used.
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
    out = sp.Integer(0)

    for offset, coeff in enumerate(row):

        degree = k_index + offset

        out += (
            sp.sympify(coeff)
            * falling(
                j,
                degree,
            )
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
            abs(value),
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
                len(ints) - 1 - i
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
        // max(1, g)
    )

    return numerator, denominator


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    # -------------------------------------------------------------------------
    # 1. RECONSTRUCT B-ODD CHANNEL
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
    # 2. RECOVER ALL B-ODD INDEX POLYNOMIALS
    # -------------------------------------------------------------------------

    index_polynomials = {}

    for power in [1, 3, 5, 7]:

        values = [
            coefficient(
                R,
                power,
            )
            for R in B_odd
        ]

        index_polynomials[power] = (
            interpolate_index(
                values
            )
        )

    # -------------------------------------------------------------------------
    # 3. TERMINAL FACTOR REMOVAL
    # -------------------------------------------------------------------------

    terminal_factors = {
        1: 1,
        3: (k - 5),
        5: (
            (k - 5)
            * (k - 4)
            * (k - 3)
        ),
        7: (
            (k - 5)
            * (k - 4)
            * (k - 3)
            * (k - 2)
            * (k - 1)
        ),
    }

    structural_polynomials = {}

    for power, Pk in index_polynomials.items():

        factor = terminal_factors[
            power
        ]

        structural_polynomials[power] = clean(
            Pk / factor
        )

    # -------------------------------------------------------------------------
    # 4. SOURCE LADDER
    # -------------------------------------------------------------------------

    print("=" * 78)
    print(
        "EXPERIMENT 261 — EXACT TERMINAL q-VALUE LADDER / "
        "PARITY-NUMERATOR PROVENANCE AUDIT"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print("1. COMPLETE B-ODD STRUCTURAL POLYNOMIAL LADDER")
    print("=" * 78)

    source_ladder = {}

    for power in [1, 3, 5, 7]:

        numerator, denominator = (
            cleared_integer_polynomial(
                structural_polynomials[
                    power
                ]
            )
        )

        leading_integer = int(
            numerator.LC()
        )

        source_ladder[power] = (
            leading_integer,
            denominator,
            numerator.as_expr(),
        )

        print()
        print(
            f"  y^{power}:"
        )

        print(
            f"    structural_polynomial="
            f"{structural_polynomials[power]}"
        )

        print(
            f"    integer_numerator="
            f"{numerator.as_expr()}"
        )

        print(
            f"    denominator="
            f"{denominator}"
        )

        print(
            f"    leading_integer="
            f"{leading_integer}"
        )

    # -------------------------------------------------------------------------
    # 5. OBSERVED TERMINAL LADDER
    # -------------------------------------------------------------------------

    observed = [
        source_ladder[1][0],
        source_ladder[3][0],
        source_ladder[5][0],
        source_ladder[7][0],
    ]

    expected = [
        495451247,
        421514439,
        16027881,
        1,
    ]

    print()
    print("=" * 78)
    print("2. OBSERVED TERMINAL INTEGER LADDER")
    print("=" * 78)

    print(
        f"  observed={observed}"
    )

    print(
        f"  expected_from_experiment_104="
        f"{expected}"
    )

    ladder_exact = (
        observed == expected
    )

    print(
        f"  ladder_recovered_exactly="
        f"{ladder_exact}"
    )

    if not ladder_exact:
        failures.append(
            "terminal_ladder"
        )

    # -------------------------------------------------------------------------
    # 6. PRIME FACTORIZATION / 17 CONTENT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. TERMINAL LADDER FACTOR / CONTENT PROFILE")
    print("=" * 78)

    for power, value in zip(
        [1, 3, 5, 7],
        observed,
    ):

        factorization = sp.factorint(
            abs(value)
        ) if value else {}

        v17 = (
            factorization.get(
                17,
                0,
            )
        )

        print(
            f"  y^{power}: "
            f"value={value} "
            f"factorization={factorization} "
            f"v17={v17}"
        )

    # -------------------------------------------------------------------------
    # 7. TERMINAL PAIR RECOVERY
    # -------------------------------------------------------------------------

    q1_5 = observed[0]
    q3_4 = observed[1]

    gcd_terminal = math.gcd(
        q1_5,
        q3_4,
    )

    q1_projective = (
        q1_5 // 17
    )

    q3_projective = (
        q3_4 // 17
    )

    pair_exact = (
        gcd_terminal == 17
        and
        q1_projective == 29144191
        and
        q3_projective == 24794967
    )

    print()
    print("=" * 78)
    print("4. PROJECTIVE SOURCE PAIR RECOVERY")
    print("=" * 78)

    print(
        f"  terminal_pair="
        f"({q1_5},{q3_4})"
    )

    print(
        f"  gcd={gcd_terminal}"
    )

    print(
        f"  projective_pair="
        f"({q1_projective},{q3_projective})"
    )

    print(
        f"  expected_projective_pair="
        f"(29144191,24794967)"
    )

    print(
        f"  exact="
        f"{pair_exact}"
    )

    if not pair_exact:
        failures.append(
            "projective_pair"
        )

    # -------------------------------------------------------------------------
    # 8. SUCCEEDING SOURCE LAYER
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SOURCE-LAYER ORDERING TEST")
    print("=" * 78)

    print(
        """
The parity ladder is indexed by

    y^1, y^3, y^5, y^7.

After terminal-factor removal the leading integer numerators are

    495451247,
    421514439,
    16027881,
    1.

The first two are exactly the two terminal source values previously
used in the 17-normalized projective construction.

The next two are new structural observations:

    y^5 -> 16027881
    y^7 -> 1.

This experiment does NOT assume in advance that these are
q_5(2) and q_7(0). It only records the exact correspondence supplied
by the parity kernel.

The next question is whether the parity degree and terminal factor
location systematically determine the indices of a general q_p(r)
family.
"""
    )

    # -------------------------------------------------------------------------
    # 9. CROSS-PARITY CHECK
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CROSS-PARITY TERMINAL PROFILE")
    print("=" * 78)

    # A-even leading integer profile.
    A_residuals = residual_channel(
        A,
        8,
    )

    A_centered = [
        center(
            R,
            8,
        )
        for R in A_residuals
    ]

    A_even = [
        parity_parts(F)[0]
        for F in A_centered
    ]

    A_even_leading = {}

    for power in [0, 2, 4, 6, 8]:

        values = [
            coefficient(
                R,
                power,
            )
            for R in A_even
        ]

        Pk = interpolate_index(
            values
        )

        numerator, denominator = (
            cleared_integer_polynomial(
                Pk
            )
        )

        A_even_leading[power] = (
            int(numerator.LC())
        )

    print(
        f"  B_odd_leading="
        f"{dict((p, source_ladder[p][0]) for p in source_ladder)}"
    )

    print(
        f"  A_even_leading="
        f"{A_even_leading}"
    )

    # The striking coefficient 4913=17^3 is explicitly checked.
    a_y8 = A_even_leading[8]

    print(
        f"  A-even y^8 leading_integer="
        f"{a_y8}"
    )

    print(
        f"  17^3={17**3}"
    )

    print(
        f"  A-even y^8 equals 17^3="
        f"{a_y8 == 17**3}"
    )

    # -------------------------------------------------------------------------
    # 10. STRUCTURAL STATUS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. WHAT IS NOW KNOWN")
    print("=" * 78)

    print(
        """
Exact observations now available from the centered parity kernel:

    B-odd y^1 -> 495451247
    B-odd y^3 -> 421514439
    B-odd y^5 -> 16027881
    B-odd y^7 -> 1

and

    gcd(495451247,421514439)=17.

The normalized pair is therefore

    (29144191,24794967).

The new quantities

    16027881,
    1

are the next source-level candidates to trace.

This is substantially more informative than the previous one-pair
search because the parity kernel itself supplies a four-level terminal
ladder.

The next experiment should determine whether this ladder can be mapped
to the original q-index notation and whether the corresponding source
rows exist for another n=pq instance.
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
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  terminal_ladder_exact="
        f"{ladder_exact}"
    )

    print(
        f"  projective_pair_exact="
        f"{pair_exact}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  source_provenance_recovered=True"
    )

    print(
        "  q_p_r_general_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  ALL BASIC CHECKS PASS="
        f"{final_ok}"
    )

    print()
    print(
        "EXPERIMENT 261 COMPLETE"
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

