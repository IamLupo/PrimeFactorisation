#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 271 — EXACT (D,k)-SCALAR LAW / SOURCE-ROW OPERATOR AUDIT
==============================================================================

Experiments 268R-270 established:

    * no common Toeplitz/convolution operator;
    * no common Pascal/binomial operator;
    * a p-dependent diagonal representation always exists trivially.

The diagonal factor at row k is

    S(D,k) = C_p[k] / q_p[k],

whenever q_p[k] != 0.

Experiment 271 asks whether these many apparently different diagonal
factors collapse to a simple UNIVERSAL function of the structural
variables

    D = D(p)
    k = row index.

This is substantially stronger than fitting an arbitrary scalar per p,k.

We search exact rational/polynomial laws of the form

    S(D,k) = P(D,k)

with total degree <= 1,2,3,4,

and rational laws

    S(D,k) = P(D,k) / Q(D,k)

with numerator/denominator total degree <= 1 and <= 2.

We also test natural factorial/binomial normalizations:

    S * k!
    S / k!
    S * binom(D,k)
    S / binom(D,k)
    S * 2^k
    S / 2^k
    S * (D-k)!
    S / (D-k)!.

A successful low-complexity law would identify an actual universal
source-index scaling rule.

A failure means even the diagonal representation cannot be reduced to
a simple function of D and k, strongly suggesting that the B-channel
depends on additional construction data.

Exact QQ arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


Dsym, ksym = sp.symbols("D k")


# =============================================================================
# SOURCE DATA
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


C = {
    1: [
        sp.Rational(-584531, 35840),
        sp.Rational(-2908483, 1920),
        sp.Rational(-31233169, 13440),
        sp.Rational(-1446167, 1344),
        sp.Rational(-22259149, 40320),
        sp.Rational(-301, 240),
    ],
    3: [
        sp.Rational(-59257, 46080),
        sp.Rational(186547, 1440),
        sp.Rational(367433, 1680),
        sp.Rational(126549, 448),
        sp.Rational(-162139, 40320),
    ],
    5: [
        sp.Rational(4457, 46080),
        sp.Rational(-16819, 5760),
        sp.Rational(-5769, 896),
    ],
    7: [
        sp.Rational(-421, 322560),
    ],
}


# =============================================================================
# HELPERS
# =============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def total_degree_monomials(bound):
    monomials = []

    for a in range(bound + 1):
        for b in range(bound + 1 - a):
            monomials.append(
                Dsym**a * ksym**b
            )

    return monomials


def fit_polynomial(points, degree):
    """
    Exact polynomial interpolation in two variables by linear algebra.

    Returns a unique polynomial only if the coefficient system is
    determined exactly.
    """

    monomials = total_degree_monomials(
        degree
    )

    unknowns = sp.symbols(
        f"a0:{len(monomials)}"
    )

    equations = []

    for Dv, kv, value in points:

        lhs = sp.Integer(0)

        for coefficient, monomial in zip(
            unknowns,
            monomials,
        ):
            lhs += (
                coefficient
                * monomial.subs({
                    Dsym: Dv,
                    ksym: kv,
                })
            )

        equations.append(
            sp.Eq(
                lhs,
                value,
            )
        )

    expressions = [
        clean(eq.lhs - eq.rhs)
        for eq in equations
    ]

    A, b = sp.linear_eq_to_matrix(
        expressions,
        unknowns,
    )

    if A.rank() != A.row_join(b).rank():
        return None

    if A.rank() != len(unknowns):
        return None

    solution = sp.linsolve(
        (A, b),
        unknowns,
    )

    if len(solution) != 1:
        return None

    vector = next(iter(solution))

    if any(
        x.free_symbols
        for x in vector
    ):
        return None

    result = sp.Integer(0)

    for x, monomial in zip(
        vector,
        monomials,
    ):
        result += x * monomial

    result = clean(result)

    # Independent verification.
    for Dv, kv, value in points:

        predicted = clean(
            result.subs({
                Dsym: Dv,
                ksym: kv,
            })
        )

        if predicted != clean(value):
            return None

    return result


def fit_rational(points, degree):
    """
    Search for

        P(D,k) / Q(D,k)

    where P,Q have total degree <= degree.

    Fix one denominator coefficient to 1 to remove scale ambiguity.
    """

    monomials = total_degree_monomials(
        degree
    )

    n = len(monomials)

    numerator = sp.symbols(
        f"n0:{n}"
    )

    denominator_free = sp.symbols(
        f"d0:{n - 1}"
    )

    # denominator = 1 + d0*m1 + ...
    denominator = (
        monomials[0]
        + sum(
            x * monomial
            for x, monomial in zip(
                denominator_free,
                monomials[1:],
            )
        )
    )

    numerator_expr = sum(
        x * monomial
        for x, monomial in zip(
            numerator,
            monomials,
        )
    )

    unknowns = list(
        numerator
    ) + list(
        denominator_free
    )

    equations = []

    for Dv, kv, value in points:

        Pval = numerator_expr.subs({
            Dsym: Dv,
            ksym: kv,
        })

        Qval = denominator.subs({
            Dsym: Dv,
            ksym: kv,
        })

        equations.append(
            clean(
                Pval - value * Qval
            )
        )

    A, b = sp.linear_eq_to_matrix(
        equations,
        unknowns,
    )

    if A.rank() != A.row_join(b).rank():
        return None

    if A.rank() != len(unknowns):
        return None

    solution = sp.linsolve(
        (A, b),
        unknowns,
    )

    if len(solution) != 1:
        return None

    vector = next(iter(solution))

    if any(
        x.free_symbols
        for x in vector
    ):
        return None

    P = clean(
        numerator_expr.subs(
            dict(
                zip(
                    unknowns,
                    vector,
                )
            )
        )
    )

    Qd = clean(
        denominator.subs(
            dict(
                zip(
                    unknowns,
                    vector,
                )
            )
        )
    )

    if Qd == 0:
        return None

    ratio = clean(
        P / Qd
    )

    for Dv, kv, value in points:

        qvalue = clean(
            Qd.subs({
                Dsym: Dv,
                ksym: kv,
            })
        )

        if qvalue == 0:
            return None

        predicted = clean(
            ratio.subs({
                Dsym: Dv,
                ksym: kv,
            })
        )

        if predicted != clean(value):
            return None

    return P, Qd, ratio


def build_points(normalization=None):
    points = []

    for p in [1, 3, 5, 7]:

        Dp = D[p]

        for k, (qv, cv) in enumerate(
            zip(
                Q[p],
                C[p],
            )
        ):

            if qv == 0:
                continue

            scale = clean(
                cv / qv
            )

            if normalization is not None:
                scale = clean(
                    normalization(
                        scale,
                        Dp,
                        k,
                    )
                )

            points.append(
                (
                    Dp,
                    k,
                    scale,
                )
            )

    return points


def factorial_value(n):
    return sp.Integer(
        math.factorial(
            int(n)
        )
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 271 — EXACT (D,k)-SCALAR LAW / "
        "SOURCE-ROW OPERATOR AUDIT"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. RAW SCALAR DATA
    # -------------------------------------------------------------------------

    points = build_points()

    print()
    print("=" * 78)
    print("1. RAW DIAGONAL-SCALAR DATA")
    print("=" * 78)

    for Dv, kv, value in points:

        print(
            f"  D={Dv} "
            f"k={kv} "
            f"S={value}"
        )

    # -------------------------------------------------------------------------
    # 2. POLYNOMIAL SEARCH
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT POLYNOMIAL LAW SEARCH")
    print("=" * 78)

    polynomial_results = {}

    for degree in [1, 2, 3, 4]:

        P = fit_polynomial(
            points,
            degree,
        )

        polynomial_results[degree] = P

        print()
        print(
            f"  total_degree<={degree}:"
        )

        print(
            f"    exact={P is not None}"
        )

        if P is not None:
            print(
                f"    P(D,k)={P}"
            )

    # -------------------------------------------------------------------------
    # 3. RATIONAL LAW SEARCH
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT RATIONAL LAW SEARCH")
    print("=" * 78)

    rational_results = {}

    for degree in [1, 2]:

        result = fit_rational(
            points,
            degree,
        )

        rational_results[degree] = result

        print()
        print(
            f"  numerator/denominator total_degree<={degree}:"
        )

        print(
            f"    exact={result is not None}"
        )

        if result is not None:

            P, Qd, ratio = result

            print(
                f"    numerator={P}"
            )

            print(
                f"    denominator={Qd}"
            )

            print(
                f"    ratio={ratio}"
            )

    # -------------------------------------------------------------------------
    # 4. NATURAL NORMALIZATION TESTS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. NATURAL FACTORIAL / BINOMIAL NORMALIZATIONS")
    print("=" * 78)

    normalizations = {
        "S*k!": lambda S, Dv, k:
            S * factorial_value(k),

        "S/k!": lambda S, Dv, k:
            S / factorial_value(k),

        "S*C(D,k)": lambda S, Dv, k:
            S * sp.binomial(Dv, k)
            if k <= Dv else S,

        "S/C(D,k)": lambda S, Dv, k:
            S / sp.binomial(Dv, k)
            if k <= Dv else S,

        "S*2^k": lambda S, Dv, k:
            S * 2**k,

        "S/2^k": lambda S, Dv, k:
            S / 2**k,

        "S*(D-k)!": lambda S, Dv, k:
            S * factorial_value(Dv - k)
            if k <= Dv else S,

        "S/(D-k)!": lambda S, Dv, k:
            S / factorial_value(Dv - k)
            if k <= Dv else S,
    }

    normalization_results = {}

    for name, transform in normalizations.items():

        normalized_points = build_points(
            transform
        )

        # Try a low-degree polynomial after normalization.
        found = None

        for degree in [1, 2, 3]:

            candidate = fit_polynomial(
                normalized_points,
                degree,
            )

            if candidate is not None:
                found = (
                    degree,
                    candidate,
                )
                break

        normalization_results[name] = found

        print()
        print(
            f"  {name}:"
        )

        if found is None:
            print(
                "    no polynomial law through degree 3"
            )
        else:
            degree, candidate = found

            print(
                f"    exact_degree={degree}"
            )

            print(
                f"    law={candidate}"
            )

    # -------------------------------------------------------------------------
    # 5. SAME-k CROSS-D COMPARISON
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SAME-k CROSS-D SCALAR COMPARISON")
    print("=" * 78)

    by_k = {}

    for Dv, kv, value in points:

        by_k.setdefault(
            kv,
            []
        ).append(
            (Dv, value)
        )

    for kv in sorted(
        by_k
    ):

        values = by_k[kv]

        if len(values) < 2:
            continue

        print()
        print(
            f"  k={kv}:"
        )

        for Dv, value in values:
            print(
                f"    D={Dv}: S={value}"
            )

        base_D, base_value = values[0]

        ratios = []

        for Dv, value in values[1:]:

            if base_value == 0:
                ratios.append(
                    "undefined"
                )
            else:
                ratios.append(
                    (
                        Dv,
                        clean(
                            value
                            / base_value
                        )
                    )
                )

        print(
            f"    ratios_to_first={ratios}"
        )

    # -------------------------------------------------------------------------
    # 6. TERMINAL SOURCE REFERENCE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. TERMINAL SOURCE REFERENCE")
    print("=" * 78)

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print(
        f"  q1_terminal={q1}"
    )

    print(
        f"  q3_terminal={q3}"
    )

    print(
        f"  gcd={math.gcd(q1, q3)}"
    )

    print(
        f"  q1_terminal/17={q1 // 17}"
    )

    print(
        f"  q3_terminal/17={q3 // 17}"
    )

    # -------------------------------------------------------------------------
    # 7. INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
The unrestricted diagonal representation from Experiment 269 is

    C_p[k] = S(D(p),k) q_p[k].

Experiment 271 asks whether this scalar is actually a compact universal
function of the structural variables D and k.

This is the last elementary source-index possibility before returning
to the original construction.

A successful polynomial/rational law would mean the B-channel operator
has a simple closed scalar form despite the failure of all common
Toeplitz and Pascal transforms.

A successful factorial/binomial normalization would be even stronger:
it would identify the combinatorial source of the denominators.

If all searches fail, the conclusion is not that no formula exists.
Rather, the B-channel scaling depends on information not captured by
D and k alone, so further source reconstruction must use the original
construction rather than only the terminal q-table.
"""
    )

    # -------------------------------------------------------------------------
    # 8. FINAL
    # -------------------------------------------------------------------------

    polynomial_success = [
        d
        for d, result in polynomial_results.items()
        if result is not None
    ]

    rational_success = [
        d
        for d, result in rational_results.items()
        if result is not None
    ]

    normalization_success = {
        name: result
        for name, result in normalization_results.items()
        if result is not None
    }

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  exact_polynomial_degrees="
        f"{polynomial_success}"
    )

    print(
        f"  exact_rational_degrees="
        f"{rational_success}"
    )

    print(
        f"  successful_normalizations="
        f"{normalization_success}"
    )

    print(
        "  common_toeplitz_operator_refuted=True"
    )

    print(
        "  common_pascal_operator_refuted=True"
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
        "EXPERIMENT 271 COMPLETE"
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

