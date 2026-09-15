#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 278 — EXACT ORIGINAL B-CHANNEL RECONSTRUCTION / DATA-INTEGRITY AUDIT
==============================================================================

Purpose:

    Repair the transcription error exposed by Experiment 277.

This experiment uses the ORIGINAL exact B data from Experiment 104 and
does nothing speculative.

For every k:

    B[k,r]
       ->
    P_k(j) = sum_r B[k,r] j_(r)
       ->
    exact residual division
       ->
    centered polynomial
       ->
    even/odd decomposition.

The derived coefficient dictionaries are compared against the exact
Experiment-104 coefficient families.

Important:

    k is the B-channel row index.

    It is NOT the same thing as the source-family index p.

This experiment must pass before any q_p(r) -> B[k,r] provenance analysis
is considered meaningful.

Exact QQ arithmetic only.
"""


from __future__ import annotations

import sys
import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

j, y = sp.symbols("j y")


# =============================================================================
# EXACT ORIGINAL B DATA FROM EXPERIMENT 104
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
# EXACT ORIGINAL A/B PARITY COEFFICIENTS FROM EXPERIMENT 104
#
# Stored by k.
#
# Each dictionary is:
#
#     power -> coefficient
#
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
        4: sp.Rational(-101119, 315),
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
        7: sp.Rational(-421, 0) if False else sp.Integer(0),
    },
    2: {
        1: sp.Rational(-31233169, 13440),
        3: sp.Rational(367433, 1680),
        5: sp.Rational(-5769, 896),
        7: sp.Integer(0),
    },
    3: {
        1: sp.Rational(-1446167, 1344),
        3: sp.Rational(126549, 448),
        5: sp.Integer(0),
        7: sp.Integer(0),
    },
    4: {
        1: sp.Rational(-22259149, 40320),
        3: sp.Rational(-162139, 40320),
        5: sp.Integer(0),
        7: sp.Integer(0),
    },
    5: {
        1: sp.Rational(-301, 240),
        3: sp.Integer(0),
        5: sp.Integer(0),
        7: sp.Integer(0),
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


def exact_poly(expr):
    return sp.Poly(
        sp.expand(expr),
        j,
        domain=sp.QQ,
    )


def exact_quotient(num, den):
    num = clean(num)
    den = clean(den)

    if den == 0:
        return False, None

    q, r = sp.div(
        exact_poly(num),
        exact_poly(den),
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
    # Correct symbolic substitution:
    # j = (y + 7)/2
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
    if clean(expr) == 0:
        return {}

    p = sp.Poly(
        sp.expand(expr),
        y,
        domain=sp.QQ,
    )

    result = {}

    for power in range(
        0,
        int(p.degree()) + 1,
    ):
        value = sp.Rational(
            p.nth(power)
        )

        if value != 0:
            result[power] = value

    return result


def dict_equal(a, b):
    keys = set(a) | set(b)

    return all(
        clean(
            a.get(power, 0)
            - b.get(power, 0)
        ) == 0
        for power in keys
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 278 — EXACT ORIGINAL B-CHANNEL "
        "RECONSTRUCTION / DATA-INTEGRITY AUDIT"
    )
    print("=" * 78)

    all_exact = True

    # =========================================================================
    # 1. RAW B ROWS
    # =========================================================================

    print()
    print("=" * 78)
    print("1. EXACT B-ROW DATA")
    print("=" * 78)

    for k, row in enumerate(B):
        print()
        print(
            f"  k={k}:"
        )
        print(
            f"    row={row}"
        )
        print(
            f"    primitive_length={len(row)}"
        )

    # =========================================================================
    # 2. COMPLETE PIPELINE
    # =========================================================================

    print()
    print("=" * 78)
    print("2. COMPLETE RESIDUAL PIPELINE")
    print("=" * 78)

    derived_even = {}
    derived_odd = {}

    for k, row in enumerate(B):

        P = row_polynomial(
            row
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
                "    ERROR: complete row "
                "is not exactly divisible."
            )

            continue

        F = center(R)
        E = even_part(F)
        O = odd_part(F)

        even_dict = coefficient_dictionary(
            E
        )

        odd_dict = coefficient_dictionary(
            O
        )

        derived_even[k] = even_dict
        derived_odd[k] = odd_dict

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
    # 3. EXACT PARITY COMPARISON
    # =========================================================================

    print()
    print("=" * 78)
    print("3. EXACT PARITY COMPARISON")
    print("=" * 78)

    for k in range(
        len(B)
    ):

        even_ok = dict_equal(
            derived_even.get(k, {}),
            EXPECTED_EVEN.get(k, {}),
        )

        odd_ok = dict_equal(
            derived_odd.get(k, {}),
            EXPECTED_ODD.get(k, {}),
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
            f"    derived_even="
            f"{derived_even.get(k, {})}"
        )

        print(
            f"    expected_even="
            f"{EXPECTED_EVEN.get(k, {})}"
        )

        print(
            f"    even_exact={even_ok}"
        )

        print(
            f"    derived_odd="
            f"{derived_odd.get(k, {})}"
        )

        print(
            f"    expected_odd="
            f"{EXPECTED_ODD.get(k, {})}"
        )

        print(
            f"    odd_exact={odd_ok}"
        )

    # =========================================================================
    # 4. COEFFICIENT SUPPORT
    # =========================================================================

    print()
    print("=" * 78)
    print("4. EXACT PARITY SUPPORT")
    print("=" * 78)

    for k in range(
        len(B)
    ):

        even_support = sorted(
            derived_even.get(k, {}).keys()
        )

        odd_support = sorted(
            derived_odd.get(k, {}).keys()
        )

        print()
        print(
            f"  k={k}: "
            f"even_support={even_support} "
            f"odd_support={odd_support}"
        )

    # =========================================================================
    # 5. DISTINGUISH DATA-INTEGRITY FROM MATHEMATICS
    # =========================================================================

    print()
    print("=" * 78)
    print("5. DATA-INTEGRITY RESULT")
    print("=" * 78)

    print(
r"""
This experiment deliberately contains no q_p(r) data.

Its only question is:

    Are the original Experiment-104 B rows,
    the residual divisor,
    the centering map,
    and the parity tables mutually consistent?

If this passes, the B-channel is established as an exact object.

If it fails, NO q -> B provenance experiment should be trusted yet.

The earlier Experiment 277 failure is therefore treated as a
data-integrity problem first, not as evidence against the mathematics.
"""
    )

    # =========================================================================
    # 6. TERMINAL SOURCE REFERENCE
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
    print("6. TERMINAL PROJECTIVE SOURCE REFERENCE")
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
    # 7. NEXT FORK
    # =========================================================================

    print()
    print("=" * 78)
    print("7. NEXT STRUCTURAL FORK")
    print("=" * 78)

    if all_exact:

        print(
r"""
BASELINE PASSES.

The B-channel is internally exact.

The next experiment should therefore reconstruct the upstream map:

    q_p(r) -> B[k][r].

Do NOT compare q_p(r) directly against the final parity coefficients.

The correct target is the full B coefficient table.
"""
        )

    else:

        print(
r"""
BASELINE FAILS.

Do not proceed to q -> B inference.

First repair the exact B-row source data, divisor convention, or parity
tables until this baseline reproduces Experiment 104 exactly.
"""
        )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  original_B_pipeline_exact="
        f"{all_exact}"
    )

    print(
        "  q_family_not_used_in_B_reconstruction=True"
    )

    print(
        "  k_indexing_explicit=True"
    )

    print(
        "  q_to_B_operator_not_yet_tested=True"
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
        "EXPERIMENT 278 COMPLETE"
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

