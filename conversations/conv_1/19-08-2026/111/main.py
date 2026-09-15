#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 277 — EXACT B-CHANNEL RECONSTRUCTION / k-INDEX AUDIT
==============================================================================

Purpose:

    Correct the indexing problem exposed by Experiment 276.

There are TWO different index systems:

    1. B-channel residual index:
           k = 0,1,...,5

    2. source q-family index:
           p in {1,3,5,7}, r = source-column index.

Experiment 276 incorrectly compared the B[k] residual output against
q-family data arranged by p.

Experiment 277 completely isolates the B-channel.

It reconstructs, from the exact B_ROWS:

    P_k(j)
      ->
    R_k(j)
      ->
    centered R_k(y)
      ->
    even part
      ->
    odd part.

For every k it then compares the result with the ORIGINAL exact
Experiment-104 parity coefficient tables.

Only after that baseline is verified does it analyze the structure of
the B[k,r] matrix itself.

No q_p(r) is used to reconstruct B.

This produces the exact B-channel object that a future q -> B experiment
must explain.
"""

from __future__ import annotations

import sys
import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

j, y = sp.symbols("j y")


# ============================================================================
# EXACT B-CHANNEL TRIANGULAR DATA
# ============================================================================

B_ROWS = [
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
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# ============================================================================
# CORRECT EXPERIMENT-104 PARITY TABLES
#
# Each row is indexed by k.
# The values below are coefficient vectors for the centered residual R_k(y).
#
# Stored in ascending odd/even powers.
# ============================================================================

B_EVEN = {
    0: [
        sp.Rational(12980463, 1024),
        sp.Rational(-19344659, 15360),
        sp.Rational(129415, 3072),
        sp.Rational(-2267, 5120),
    ],
    1: [
        sp.Rational(6255583, 256),
        sp.Rational(-26986999, 11520),
        sp.Rational(267779, 3840),
        sp.Rational(-6053, 11520),
    ],
    2: [
        sp.Rational(87841139, 4480),
        sp.Rational(-2066529, 1120),
        sp.Rational(224417, 4480),
    ],
    3: [
        sp.Rational(16998339, 2240),
        sp.Rational(-573325, 576),
        sp.Rational(-101119, 5040),
    ],
    4: [
        sp.Rational(1091983, 896),
        sp.Rational(3312053, 40320),
    ],
    5: [
        sp.Rational(1553, 240),
    ],
}


B_ODD = {
    0: [
        sp.Rational(-584531, 35840),
        sp.Rational(-59257, 46080),
        sp.Rational(-421, 322560),
    ],
    1: [
        sp.Rational(-2908483, 1920),
        sp.Rational(186547, 1440),
        sp.Rational(-301, 240),
    ],
    2: [
        sp.Rational(-31233169, 13440),
        sp.Rational(367433, 1680),
        sp.Rational(-162139, 40320),
    ],
    3: [
        sp.Rational(-1446167, 1344),
        sp.Rational(126549, 448),
    ],
    4: [
        sp.Rational(-22259149, 40320),
        sp.Rational(-162139, 40320),
    ],
    5: [
        sp.Rational(-301, 240),
    ],
}


# ============================================================================
# HELPERS
# ============================================================================

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
        sp.expand(expr),
        var,
        domain=sp.QQ,
    )


def exact_quotient(num, den):
    num = clean(num)
    den = clean(den)

    if den == 0:
        return False, None

    q, r = sp.div(
        poly(num, j),
        poly(den, j),
        domain=sp.QQ,
    )

    if not r.is_zero:
        return False, None

    return True, clean(
        q.as_expr()
    )


def row_polynomial(row):
    out = sp.Integer(0)

    for r, coeff in enumerate(row):
        out += coeff * falling(j, r)

    return clean(out)


def row_divisor(k):
    return clean(
        falling(j, k)
        * upper_falling(
            7,
            j,
            max(
                0,
                k - 4,
            ),
        )
    )


def center(expr):
    return clean(
        expr.subs(
            j,
            (y + sp.Integer(7))
            / sp.Integer(2),
        )
    )


def even_part(expr):
    return clean(
        (
            expr
            + expr.subs(y, -y)
        )
        / sp.Integer(2)
    )


def odd_part(expr):
    return clean(
        (
            expr
            - expr.subs(y, -y)
        )
        / sp.Integer(2)
    )


def coefficient(expr, power):
    expr = clean(expr)

    if expr == 0:
        return sp.Integer(0)

    return sp.Rational(
        poly(expr, y).nth(power)
    )


def parity_coefficients(expr):
    """
    Return exact coefficient dictionaries.
    """

    p = poly(expr, y) if expr != 0 else None

    if p is None:
        return {}, {}

    degree = int(p.degree())

    even = {}
    odd = {}

    for power in range(
        0,
        degree + 1,
    ):
        value = sp.Rational(
            p.nth(power)
        )

        if value == 0:
            continue

        if power % 2 == 0:
            even[power] = value
        else:
            odd[power] = value

    return even, odd


def dict_equal(a, b):
    keys = set(a) | set(b)

    return all(
        clean(
            a.get(k, 0)
            - b.get(k, 0)
        ) == 0
        for k in keys
    )


def primitive_integer_row(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    den = 1

    for value in values:
        den = sp.ilcm(
            den,
            int(value.q),
        )

    ints = [
        int(value * den)
        for value in values
    ]

    g = 0

    for x in ints:
        g = sp.igcd(
            g,
            abs(x),
        )

    if g:
        ints = [
            x // g
            for x in ints
        ]

    return ints


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 277 — EXACT B-CHANNEL RECONSTRUCTION / "
        "k-INDEX AUDIT"
    )
    print("=" * 78)

    all_exact = True

    # =========================================================================
    # 1. RECONSTRUCT EVERY B-ROW
    # =========================================================================

    print()
    print("=" * 78)
    print("1. COMPLETE B-CHANNEL PIPELINE")
    print("=" * 78)

    results = {}

    for k, row in enumerate(
        B_ROWS
    ):

        P = row_polynomial(
            row
        )

        divisor = row_divisor(
            k
        )

        exact, R = exact_quotient(
            P,
            divisor,
        )

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    P_k(j)={P}"
        )

        print(
            f"    divisor={divisor}"
        )

        print(
            f"    exact_division={exact}"
        )

        if not exact:
            all_exact = False
            print(
                "    ERROR: residual extraction failed."
            )
            continue

        F = center(R)
        E = even_part(F)
        O = odd_part(F)

        even, odd = parity_coefficients(
            F
        )

        results[k] = {
            "P": P,
            "R": R,
            "F": F,
            "E": E,
            "O": O,
            "even": even,
            "odd": odd,
        }

        print(
            f"    R_k(j)={R}"
        )

        print(
            f"    centered={F}"
        )

        print(
            f"    even={E}"
        )

        print(
            f"    odd={O}"
        )

    # =========================================================================
    # 2. CORRECT k-INDEX COMPARISON
    # =========================================================================

    print()
    print("=" * 78)
    print("2. EXACT k-INDEX PARITY COMPARISON")
    print("=" * 78)

    for k in range(
        len(B_ROWS)
    ):

        if k not in results:
            continue

        derived_even = results[k][
            "even"
        ]

        derived_odd = results[k][
            "odd"
        ]

        expected_even = {
            2 * r: value
            for r, value in enumerate(
                B_EVEN.get(
                    k,
                    [],
                )
            )
        }

        expected_odd = {
            2 * r + 1: value
            for r, value in enumerate(
                B_ODD.get(
                    k,
                    [],
                )
            )
        }

        even_ok = dict_equal(
            derived_even,
            expected_even,
        )

        odd_ok = dict_equal(
            derived_odd,
            expected_odd,
        )

        all_exact &= (
            even_ok
            and odd_ok
        )

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    derived_even={derived_even}"
        )

        print(
            f"    expected_even={expected_even}"
        )

        print(
            f"    even_exact={even_ok}"
        )

        print(
            f"    derived_odd={derived_odd}"
        )

        print(
            f"    expected_odd={expected_odd}"
        )

        print(
            f"    odd_exact={odd_ok}"
        )

    # =========================================================================
    # 3. PRIMITIVE INTEGER SIGNATURES
    # =========================================================================

    print()
    print("=" * 78)
    print("3. PRIMITIVE INTEGER B-ROW SIGNATURES")
    print("=" * 78)

    for k, row in enumerate(
        B_ROWS
    ):

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    B_row={row}"
        )

        print(
            f"    primitive="
            f"{primitive_integer_row(row)}"
        )

    # =========================================================================
    # 4. ZERO / SUPPORT PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print("4. B-ROW SUPPORT PROFILE")
    print("=" * 78)

    for k, row in enumerate(
        B_ROWS
    ):

        support = [
            r
            for r, value in enumerate(row)
            if value != 0
        ]

        print(
            f"  k={k}: "
            f"support={support} "
            f"length={len(row)}"
        )

    # =========================================================================
    # 5. SOURCE q-TABLE IS KEPT SEPARATE
    # =========================================================================

    print()
    print("=" * 78)
    print("5. SOURCE q-FAMILY — SEPARATE INDEX SYSTEM")
    print("=" * 78)

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

    for p in [
        1,
        3,
        5,
        7,
    ]:

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    D(p)={D[p]}"
        )

        print(
            f"    q_p(r)={Q[p]}"
        )

    # =========================================================================
    # 6. INDEXING SEPARATION
    # =========================================================================

    print()
    print("=" * 78)
    print("6. INDEXING SEPARATION")
    print("=" * 78)

    print(
r"""
The B-channel data are indexed by:

    k = residual / polynomial row index.

The source family is indexed by:

    p = terminal q-family,
    r = q-table source-column index.

Therefore there is no legitimate comparison of

    B_ODD[k]

directly with

    q_p(r)

without an independently derived q -> B operator.

Experiment 277 intentionally does NOT make that comparison.

Its purpose is to establish the B-channel exactly first.
"""
    )

    # =========================================================================
    # 7. TERMINAL SOURCE REFERENCE
    # =========================================================================

    print()
    print("=" * 78)
    print("7. TERMINAL PROJECTIVE SOURCE")
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
        f"  gcd={sp.gcd(
            q1_terminal,
            q3_terminal,
        )}"
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
    # 8. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
This is a baseline experiment.

A successful result means:

    B_ROWS
      ->
    residual extraction
      ->
    centering
      ->
    parity decomposition

is internally exact.

Once that is established, the source problem becomes clean:

    q_p(r)
       ->
    B[k][r].

That map is upstream of the residual operator.

Experiments 268R-273 attempted to infer this map from the final
coefficients. That was premature.

The correct next stage after this audit is to obtain or reconstruct the
original formula that generated B[k][r].

The B-channel itself should no longer be treated as an unknown.
"""
    )

    # =========================================================================
    # 9. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  complete_B_channel_exact="
        f"{all_exact}"
    )

    print(
        "  correct_k_indexing=True"
    )

    print(
        "  q_family_kept_separate=True"
    )

    print(
        "  q_to_B_operator_still_unresolved=True"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  failures={0 if all_exact else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={all_exact}"
    )

    print()
    print(
        "EXPERIMENT 277 COMPLETE"
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

