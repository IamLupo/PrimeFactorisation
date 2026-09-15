#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 279 — EXACT ORIGINAL B-CHANNEL / SHIFTED FALLING-FACTORIAL AUDIT
==============================================================================

Experiment 278 exposed the remaining data-integrity issue.

The original Experiment-104 reconstruction uses:

    r = k + offset

not:

    r = offset.

Thus a B-row indexed by k is:

    P_k(j)
      =
    sum_offset B[k][offset] * j_(k+offset).

This shifted falling-factorial convention is essential because the residual
divisor contains j_(k):

    j_(k) = j(j-1)...(j-k+1).

Experiment 279 restores that convention exactly.

For each k:

    B[k][0]        -> j_(k)
    B[k][1]        -> j_(k+1)
    ...
    B[k][m]        -> j_(k+m)

Then:

    P_k(j)
       ->
    exact residual division
       ->
    center at y = 2j-7
       ->
    even/odd decomposition.

This experiment contains NO q_p(r) provenance inference.

It is strictly a reconstruction of Experiment 104.
"""


from __future__ import annotations

import sys
import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

j, y = sp.symbols("j y")


# =============================================================================
# ORIGINAL B DATA
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
# EXPECTED EXPERIMENT-104 PARITY DATA
# =============================================================================

EXPECTED_EVEN = {
    0: {
        0: sp.Rational(12980463, 1024),
        2: sp.Rational(-19344659, 15360),
        4: sp.Rational(129415, 3072),
        6: sp.Rational(-2267, 5120),
    },
    1: {
        0: sp.Rational(6255583, 256),
        2: sp.Rational(-26986999, 11520),
        4: sp.Rational(267779, 3840),
        6: sp.Rational(-6053, 11520),
    },
    2: {
        0: sp.Rational(87841139, 4480),
        2: sp.Rational(-2066529, 1120),
        4: sp.Rational(224417, 4480),
    },
    3: {
        0: sp.Rational(16998339, 2240),
        2: sp.Rational(-573325, 576),
        4: sp.Rational(3312053, 40320),
    },
    4: {
        0: sp.Rational(1091983, 896),
        2: sp.Rational(3312053, 40320),
    },
    5: {
        0: sp.Rational(1553, 240),
    },
}


EXPECTED_ODD = {
    0: {
        1: sp.Rational(-584531, 35840),
        3: sp.Rational(-59257, 46080),
        5: sp.Rational(4457, 46080),
        7: sp.Rational(-421, 322560),
    },
    1: {
        1: sp.Rational(-2908483, 1920),
        3: sp.Rational(186547, 1440),
        5: sp.Rational(-16819, 5760),
    },
    2: {
        1: sp.Rational(-31233169, 13440),
        3: sp.Rational(367433, 1680),
        5: sp.Rational(-5769, 896),
    },
    3: {
        1: sp.Rational(-1446167, 1344),
        3: sp.Rational(126549, 448),
    },
    4: {
        1: sp.Rational(-22259149, 40320),
        3: sp.Rational(-162139, 40320),
    },
    5: {
        1: sp.Rational(-301, 240),
    },
}


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


def exact_poly(expr, var):
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
        exact_poly(num, j),
        exact_poly(den, j),
        domain=sp.QQ,
    )

    if not r.is_zero:
        return False, None

    return True, clean(
        q.as_expr()
    )


def shifted_row_polynomial(
    row,
    k,
):
    """
    IMPORTANT:

        r = k + offset

    This is the convention from Experiment 104.
    """

    out = sp.Integer(0)

    for offset, coeff in enumerate(row):

        r = k + offset

        out += (
            sp.sympify(coeff)
            * falling(j, r)
        )

    return clean(out)


def divisor_for_k(k):
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
            (
                y
                + sp.Integer(7)
            )
            / sp.Integer(2),
        )
    )


def even_part(expr):
    return clean(
        (
            expr
            + expr.subs(
                y,
                -y,
            )
        )
        / sp.Integer(2)
    )


def odd_part(expr):
    return clean(
        (
            expr
            - expr.subs(
                y,
                -y,
            )
        )
        / sp.Integer(2)
    )


def coefficient_dictionary(expr):
    expr = clean(expr)

    if expr == 0:
        return {}

    P = sp.Poly(
        sp.expand(expr),
        y,
        domain=sp.QQ,
    )

    out = {}

    for power in range(
        0,
        int(P.degree()) + 1,
    ):

        value = sp.Rational(
            P.nth(power)
        )

        if value != 0:
            out[power] = value

    return out


def dict_equal(a, b):
    powers = set(a) | set(b)

    return all(
        clean(
            a.get(power, 0)
            - b.get(power, 0)
        ) == 0
        for power in powers
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 279 — EXACT ORIGINAL B-CHANNEL / "
        "SHIFTED FALLING-FACTORIAL AUDIT"
    )
    print("=" * 78)

    all_exact = True

    # =========================================================================
    # 1. SHIFTED FALLING-FACTORIAL RECONSTRUCTION
    # =========================================================================

    print()
    print("=" * 78)
    print("1. SHIFTED FALLING-FACTORIAL RECONSTRUCTION")
    print("=" * 78)

    results = {}

    for k, row in enumerate(B):

        P = shifted_row_polynomial(
            row,
            k,
        )

        divisor = divisor_for_k(
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
            f"    shifted_indices="
            f"{list(range(k, k + len(row)))}"
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
                "    ERROR: shifted reconstruction "
                "still fails exact division."
            )
            continue

        F = center(R)
        E = even_part(F)
        O = odd_part(F)

        edict = coefficient_dictionary(E)
        odict = coefficient_dictionary(O)

        results[k] = {
            "P": P,
            "R": R,
            "F": F,
            "E": E,
            "O": O,
            "even": edict,
            "odd": odict,
        }

        print(
            f"    residual={R}"
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
    # 2. EXACT PARITY COMPARISON
    # =========================================================================

    print()
    print("=" * 78)
    print("2. EXACT PARITY COMPARISON")
    print("=" * 78)

    for k in range(
        len(B)
    ):

        if k not in results:
            all_exact = False
            continue

        even_derived = results[k][
            "even"
        ]

        odd_derived = results[k][
            "odd"
        ]

        even_expected = EXPECTED_EVEN.get(
            k,
            {},
        )

        odd_expected = EXPECTED_ODD.get(
            k,
            {},
        )

        even_ok = dict_equal(
            even_derived,
            even_expected,
        )

        odd_ok = dict_equal(
            odd_derived,
            odd_expected,
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
            f"    even_derived={even_derived}"
        )

        print(
            f"    even_expected={even_expected}"
        )

        print(
            f"    even_exact={even_ok}"
        )

        print(
            f"    odd_derived={odd_derived}"
        )

        print(
            f"    odd_expected={odd_expected}"
        )

        print(
            f"    odd_exact={odd_ok}"
        )

    # =========================================================================
    # 3. SHIFTED INDEX TABLE
    # =========================================================================

    print()
    print("=" * 78)
    print("3. SHIFTED SOURCE-INDEX PROFILE")
    print("=" * 78)

    for k, row in enumerate(B):

        indices = list(
            range(
                k,
                k + len(row),
            )
        )

        print(
            f"  k={k}: "
            f"source_indices={indices}"
        )

    # =========================================================================
    # 4. WHY EXPERIMENT 278 FAILED
    # =========================================================================

    print()
    print("=" * 78)
    print("4. PREVIOUS-ERROR DIAGNOSIS")
    print("=" * 78)

    print(
r"""
The failure of Experiment 278 is completely explained by the missing
index shift.

For k=1 the first coefficient is attached to

    j_(1),

not

    j_(0).

For k=2 the first coefficient is attached to

    j_(2),

not

    j_(0).

In general:

    B[k][offset]
        multiplies
    j_(k+offset).

Therefore the factor j_(k) is built into every B-row from its first
basis term onward.

The apparent lack of divisibility in Experiment 278 was therefore an
artifact of using the wrong falling-factorial basis indices.
"""
    )

    # =========================================================================
    # 5. TERMINAL SOURCE
    # =========================================================================

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

    print()
    print("=" * 78)
    print("5. TERMINAL PROJECTIVE SOURCE")
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
        f"  q1/17="
        f"{q1_terminal // 17}"
    )

    print(
        f"  q3/17="
        f"{q3_terminal // 17}"
    )

    # =========================================================================
    # 6. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
If Experiment 279 passes, the B-channel baseline is finally secure.

The exact construction is:

    B[k][offset]
        ->
    j_(k+offset)
        ->
    P_k(j)
        ->
    divide by the known residual factor
        ->
    center at y = 2j-7
        ->
    even/odd parts.

Only after this passes should the source-provenance problem be attacked.

The remaining unknown is then genuinely upstream:

    q_p(r)
        ->
    B[k][offset].

This is the correct target for Experiment 280.
"""
    )

    # =========================================================================
    # 7. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  shifted_B_pipeline_exact="
        f"{all_exact}"
    )

    print(
        "  shifted_falling_factorial_convention_verified="
        f"{all_exact}"
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
        "EXPERIMENT 279 COMPLETE"
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

