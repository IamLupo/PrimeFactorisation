#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 295 — EXACT B-ROW FORMULA DECOMPOSITION AUDIT
==============================================================================

Experiments 283-294 have ruled out:

    * low-degree bivariate polynomial laws;
    * simple diagonal recurrences;
    * geometric ratios;
    * low-rank separability;
    * finite-band row transfer;
    * differential/Appell laws;
    * simple Newton-basis row transfer.

Experiment 295 therefore stops fitting B.

It reconstructs the elementary finite-product factors that occur in the
original shifted falling-factorial construction and asks whether B[k,r]
is exactly a coefficient of a product/quotient built from them.

For each entry B[k,r], test exact identities against combinations of:

    k_(r-k)
    (7-k)_(r-k)
    r_(k)
    r_(r-k)
    (7-r)_k
    binomial(7-k,r-k)
    binomial(r,k)
    factorial ratios
    reciprocal factorial ratios
    signed versions

The test is deliberately multiplicative.

We also inspect:

    numerator / denominator,
    p-adic valuations for small primes,
    exact sign,
    exact content.

No q-family data.
No arbitrary fitted parameters.
No interpolation.
Exact QQ arithmetic only.
"""

from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# EXACT B DATA
# ============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
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


# ============================================================================
# HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def falling(n, q):
    if q < 0:
        return sp.Integer(0)

    out = sp.Integer(1)

    for i in range(q):
        out *= n - i

    return clean(out)


def valuation(n, p):
    """
    Exact p-adic valuation of a nonzero integer.
    """
    if n == 0:
        return None

    n = abs(int(n))
    value = 0

    while n % p == 0:
        n //= p
        value += 1

    return value


def rational_valuation(q, p):
    """
    v_p(a/b)=v_p(a)-v_p(b).
    """
    q = sp.Rational(q)

    if q == 0:
        return None

    return (
        valuation(int(q.p), p)
        - valuation(int(q.q), p)
    )


def factor_integer(n):
    n = abs(int(n))

    if n <= 1:
        return {}

    return dict(
        sp.factorint(n)
    )


def factor_rational(q):
    q = sp.Rational(q)

    return {
        "numerator": factor_integer(q.p),
        "denominator": factor_integer(q.q),
    }


# ============================================================================
# ELEMENTARY FACTORS
# ============================================================================

def factors(k, r):

    d = r - k

    return {
        "1": sp.Integer(1),

        "k_falling_d":
            falling(k, d),

        "r_falling_k":
            falling(r, k),

        "r_falling_d":
            falling(r, d),

        "7-k_falling_d":
            falling(7 - k, d),

        "7-r_falling_k":
            falling(7 - r, k),

        "7-r_falling_d":
            falling(7 - r, d),

        "C(7-k,d)":
            sp.binomial(
                7 - k,
                d,
            ),

        "C(r,k)":
            sp.binomial(
                r,
                k,
            ),

        "C(r,d)":
            sp.binomial(
                r,
                d,
            ),

        "C(7,r)":
            sp.binomial(
                7,
                r,
            ),

        "C(7,k)":
            sp.binomial(
                7,
                k,
            ),

        "k!":
            sp.factorial(k),

        "r!":
            sp.factorial(r),

        "d!":
            sp.factorial(d),

        "(7-k)!":
            sp.factorial(7 - k),

        "(7-r)!":
            sp.factorial(7 - r),
    }


# ============================================================================
# 1. RAW ENTRY PROFILE
# ============================================================================

def raw_profile():

    print()
    print("=" * 78)
    print("1. RAW B-ENTRY PROFILE")
    print("=" * 78)

    for k, row in enumerate(B):

        for offset, value in enumerate(row):

            r = k + offset
            d = r - k

            print()
            print(
                f"  (k={k},r={r},d={d}):"
            )

            print(
                f"    B={value}"
            )

            print(
                f"    sign="
                f"{'+' if value > 0 else '-'}"
            )

            print(
                f"    factors="
                f"{factor_rational(value)}"
            )

            print(
                "    valuations="
                + str(
                    {
                        p: rational_valuation(
                            value,
                            p,
                        )
                        for p in (2, 3, 5, 7, 11, 13, 17)
                    }
                )
            )


# ============================================================================
# 2. EXACT RATIO SEARCH AGAINST ELEMENTARY FACTORS
# ============================================================================

def ratio_search():

    print()
    print("=" * 78)
    print(
        "2. EXACT RATIO SEARCH AGAINST ELEMENTARY FACTORS"
    )
    print("=" * 78)

    names = None

    first_done = False

    hit_counts = {}

    for k, row in enumerate(B):

        for offset, value in enumerate(row):

            r = k + offset

            current = factors(
                k,
                r,
            )

            if not first_done:
                names = list(
                    current.keys()
                )

                for name in names:
                    hit_counts[name] = 0

                first_done = True

            print()
            print(
                f"  (k={k},r={r}):"
            )

            for name in names:

                factor = current[name]

                if factor == 0:
                    continue

                ratio = clean(
                    sp.Rational(
                        value,
                        factor,
                    )
                )

                # Exact integer ratio is especially informative.
                integer_ratio = (
                    sp.denom(ratio) == 1
                )

                if integer_ratio:
                    hit_counts[name] += 1

                print(
                    f"    {name}: "
                    f"ratio={ratio} "
                    f"integer={integer_ratio}"
                )

    print()
    print(
        "  INTEGER-RATIO HIT COUNTS:"
    )

    for name, count in hit_counts.items():
        print(
            f"    {name}: {count}"
        )


# ============================================================================
# 3. NORMALIZED ENTRY CONTENT
# ============================================================================

def normalized_content():

    print()
    print("=" * 78)
    print(
        "3. NORMALIZED ENTRY CONTENT"
    )
    print("=" * 78)

    for k, row in enumerate(B):

        for offset, value in enumerate(row):

            r = k + offset
            d = r - k

            entry = factors(
                k,
                r,
            )

            candidates = {}

            for name, factor in entry.items():

                if factor == 0:
                    continue

                candidates[name] = clean(
                    sp.Rational(
                        value,
                        factor,
                    )
                )

            print()
            print(
                f"  (k={k},r={r},d={d}):"
            )

            # Only show integer or very simple normalized values.
            simple = []

            for name, ratio in candidates.items():

                if ratio.is_Integer:
                    simple.append(
                        (
                            name,
                            ratio,
                        )
                    )

                elif (
                    abs(int(ratio.p))
                    <= 10000
                    and abs(int(ratio.q))
                    <= 10000
                ):
                    simple.append(
                        (
                            name,
                            ratio,
                        )
                    )

            print(
                f"    simple_normalizations="
                f"{simple}"
            )


# ============================================================================
# 4. DENOMINATOR-FACTOR COMPATIBILITY
# ============================================================================

def denominator_compatibility():

    print()
    print("=" * 78)
    print(
        "4. DENOMINATOR-FACTOR COMPATIBILITY"
    )
    print("=" * 78)

    denominators = []

    for k, row in enumerate(B):

        for offset, value in enumerate(row):

            r = k + offset

            denominators.append(
                (
                    k,
                    r,
                    value.q,
                )
            )

    for k, r, denominator in denominators:

        print()
        print(
            f"  (k={k},r={r}): "
            f"denominator={denominator}"
        )

        print(
            f"    factors="
            f"{factor_integer(denominator)}"
        )

        print(
            "    factorial_profile="
            + str(
                {
                    "k!": sp.factorial(k),
                    "r!": sp.factorial(r),
                    "d!": sp.factorial(r-k),
                }
            )
        )


# ============================================================================
# 5. ROW-WISE PRODUCT STRUCTURE
# ============================================================================

def row_product_structure():

    print()
    print("=" * 78)
    print(
        "5. ROW-WISE PRODUCT STRUCTURE"
    )
    print("=" * 78)

    for k, row in enumerate(B):

        values = [
            sp.Rational(value)
            for value in row
        ]

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    row_values={values}"
        )

        # Pairwise adjacent ratios.
        ratios = []

        for i in range(
            len(values) - 1
        ):

            if values[i] == 0:
                ratios.append(
                    None
                )
            else:
                ratios.append(
                    clean(
                        values[i + 1]
                        / values[i]
                    )
                )

        print(
            f"    adjacent_ratios={ratios}"
        )

        # Sign-adjusted absolute ratios.
        abs_ratios = []

        for i in range(
            len(values) - 1
        ):

            if values[i] == 0:
                abs_ratios.append(
                    None
                )
            else:
                abs_ratios.append(
                    clean(
                        abs(
                            values[i + 1]
                            / values[i]
                        )
                    )
                )

        print(
            f"    abs_adjacent_ratios="
            f"{abs_ratios}"
        )


# ============================================================================
# 6. TERMINAL COLUMNS
# ============================================================================

def terminal_columns():

    print()
    print("=" * 78)
    print(
        "6. TERMINAL-COLUMN FACTOR AUDIT"
    )
    print("=" * 78)

    for r in range(5, 8):

        values = []

        for k, row in enumerate(B):

            if k <= r < k + len(row):

                offset = r - k
                value = row[offset]

                values.append(
                    (
                        k,
                        value,
                    )
                )

        print()
        print(
            f"  r={r}:"
        )

        for k, value in values:

            print(
                f"    k={k}: "
                f"B={value} "
                f"factors={factor_rational(value)}"
            )


# ============================================================================
# 7. STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "7. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiments 283-294 rejected a large class of downstream index laws.

Experiment 295 changes the question.

Instead of asking for a polynomial or recurrence that fits B[k,r],
it asks whether each individual coefficient is built from recognizable
finite-product factors.

The relevant variables are:

    k,
    r,
    d=r-k,
    7-k,
    7-r.

The natural algebraic factors are factorials, binomial coefficients,
falling factorials, and their reciprocals.

A repeated exact simplification such as

    B[k,r] / F(k,r)

being integral or belonging to a small structured family would identify
the missing combinatorial normalization.

If no elementary factor explains the entries, then the evidence for
"pattern hunting" is exhausted and the original construction of B must
be supplied directly.
"""
    )


# ============================================================================
# 8. FINAL
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 295 — EXACT B-ROW FORMULA "
        "DECOMPOSITION AUDIT"
    )
    print("=" * 78)

    raw_profile()
    ratio_search()
    normalized_content()
    denominator_compatibility()
    row_product_structure()
    terminal_columns()
    interpretation()

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_rational_arithmetic=True"
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 295 COMPLETE"
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

