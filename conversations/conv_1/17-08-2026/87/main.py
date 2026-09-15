#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 151
DIRECT CLOSED-FORM LAURENT KERNEL MODES

GOAL
----
Experiment 150 discovered:

    C_0(t) = 2*binom(ell,k)*t^(2k)
           = 2*binom(ell,k)*N^k,

while every informative nonzero Laurent mode retains hidden u-dependence.

This experiment derives the GENERAL coefficient

    C_j(t) = [u^j] F_(k,ell)(t*u, t/u)

DIRECTLY from the four finite binomial terms.

No quotient is constructed.

For

    F_(k,ell)(p,q)
      = p^k(1+q)^ell
      + q^k(1+p)^ell
      - p^ell(1+q)^k
      - q^ell(1+p)^k,

after

    p=t*u, q=t*u^(-1),

derive an exact finite binomial formula for every Laurent mode.

TARGET
------
Find and verify a formula of the form

    C_j(t)
      = sum of a small number of binomial terms
        times powers of t.

Then test:

    1. symmetry:
           C_j = C_(-j)

    2. zero-mode:
           C_0 = 2*binom(ell,k)*N^k

    3. exact closed formulas for j=1,2,3

    4. separation of t-power from the integer/binomial coefficient.

This is the input needed for Experiment 152, where

    F = (1+t(u+u^(-1))) Q

will be converted into an exact recurrence for the Laurent modes of Q.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO C/D TENSOR
NO POLYNOMIAL FITTING
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
t, u = sp.symbols("t u")


# ============================================================================
# Original finite kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Direct Laurent kernel
# ============================================================================

def direct_laurent(k, ell):
    F = kernel_F(k, ell)

    return sp.expand(
        F.subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


# ============================================================================
# Exact Laurent coefficient from expansion
# ============================================================================

def expanded_laurent_coeff(expr, j):
    """
    Extract [u^j] from the explicit Laurent polynomial.
    """

    expr = sp.expand(expr)

    support = []

    for term in sp.Add.make_args(expr):

        exponent = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not exponent.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent: {term}"
            )

        support.append(
            int(exponent)
        )

    if not support:
        return sp.Integer(0)

    shift = max(
        0,
        -min(support)
    )

    shifted = sp.expand(
        expr * u**shift
    )

    poly = sp.Poly(
        shifted,
        u,
        domain=sp.QQ.frac_field(t)
    )

    target = j + shift

    if target < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(
            u**target
        )
    )


# ============================================================================
# Direct finite-binomial coefficient formula
# ============================================================================

def binomial_term(
    outer_power,
    inner_power,
    target_j,
    sign
):
    """
    Compute the Laurent contribution of

        (t u)^outer_power
        (1 + t u^(-1))^inner_power

    to u^target_j.

    Expansion:

        t^outer * sum_r binom(inner,r)
        t^r u^(outer-r).

    Therefore:

        outer-r = j
        r = outer-j.

    Contribution:

        sign * binom(inner, outer-j) * t^(outer + outer-j)
              if the binomial index is valid.
    """

    r = outer_power - target_j

    if r < 0 or r > inner_power:
        return sp.Integer(0)

    return sp.expand(
        sign
        * sp.binomial(
            inner_power,
            r
        )
        * t**(
            outer_power + r
        )
    )


def predicted_mode(k, ell, j):
    """
    Direct four-term Laurent coefficient.

    Term 1:
        p^k (1+q)^ell

    Term 2:
        q^k (1+p)^ell

    Term 3:
        -p^ell (1+q)^k

    Term 4:
        -q^ell (1+p)^k

    The second and fourth terms are easiest to write by symmetry
    as contributions at -j.
    """

    # ------------------------------------------------------------------------
    # Term 1
    # ------------------------------------------------------------------------

    c1 = binomial_term(
        outer_power=k,
        inner_power=ell,
        target_j=j,
        sign=1
    )

    # ------------------------------------------------------------------------
    # Term 2:
    #
    # (t/u)^k (1+t*u)^ell
    #
    # Expansion exponent:
    #
    #   -k + r = j
    #
    # so:
    #
    #   r = j+k.
    # ------------------------------------------------------------------------

    r2 = j + k

    if 0 <= r2 <= ell:

        c2 = sp.expand(
            sp.binomial(
                ell,
                r2
            )
            * t**(
                k + r2
            )
        )

    else:
        c2 = sp.Integer(0)

    # ------------------------------------------------------------------------
    # Term 3:
    #
    # - (t*u)^ell (1+t*u^-1)^k
    # ------------------------------------------------------------------------

    c3 = binomial_term(
        outer_power=ell,
        inner_power=k,
        target_j=j,
        sign=-1
    )

    # ------------------------------------------------------------------------
    # Term 4:
    #
    # - (t/u)^ell (1+t*u)^k
    #
    # -ell + r = j
    # r = j+ell.
    # ------------------------------------------------------------------------

    r4 = j + ell

    if 0 <= r4 <= k:

        c4 = sp.expand(
            -sp.binomial(
                k,
                r4
            )
            * t**(
                ell + r4
            )
        )

    else:
        c4 = sp.Integer(0)

    return sp.expand(
        c1 + c2 + c3 + c4
    )


# ============================================================================
# Simplify coefficient by parity / common t-power
# ============================================================================

def factor_t_power(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return {
            "power": None,
            "coefficient": sp.Integer(0)
        }

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ
    )

    powers = [
        int(monom[0])
        for monom in poly.monoms()
    ]

    common = min(powers)

    reduced = sp.expand(
        expr / t**common
    )

    return {
        "power": common,
        "coefficient": sp.factor(reduced)
    }


# ============================================================================
# Exact zero-mode formula
# ============================================================================

def expected_zero_mode(k, ell):
    return sp.expand(
        2
        * sp.binomial(ell, k)
        * t**(2*k)
    )


# ============================================================================
# Exact mode formula diagnostics
# ============================================================================

def verify_formula(k, ell, j):
    R = direct_laurent(
        k,
        ell
    )

    exact = expanded_laurent_coeff(
        R,
        j
    )

    predicted = predicted_mode(
        k,
        ell,
        j
    )

    residual = sp.factor(
        sp.expand(
            exact - predicted
        )
    )

    return exact, predicted, residual


# ============================================================================
# Closed-form mode table
# ============================================================================

def print_mode_table(k, ell, max_j):
    print(
        f"({k},{ell})"
    )

    R = direct_laurent(
        k,
        ell
    )

    support = []

    for term in sp.Add.make_args(
        sp.expand(R)
    ):

        exponent = int(
            sp.sympify(
                term.as_powers_dict().get(u, 0)
            )
        )

        support.append(
            exponent
        )

    support = sorted(
        set(support)
    )

    print(
        "  support =",
        support
    )

    for j in range(
        0,
        min(max_j, max(abs(v) for v in support)) + 1
    ):

        exact = expanded_laurent_coeff(
            R,
            j
        )

        if exact == 0:
            continue

        data = factor_t_power(
            exact
        )

        print(
            f"  C_{j}(t) = {sp.factor(exact)}"
        )

        print(
            f"      common t-power = {data['power']}"
        )

        print(
            f"      reduced coefficient = "
            f"{data['coefficient']}"
        )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 151")
    print("DIRECT CLOSED-FORM LAURENT KERNEL MODES")
    print("=" * 78)
    print()

    PAIRS = [
        (1, 3),
        (1, 5),
        (1, 7),
        (1, 9),
        (3, 5),
        (3, 7),
        (3, 9),
        (5, 7),
        (5, 9),
    ]

    TEST_J = [
        0,
        1,
        2,
        3,
    ]

    total_formula_failures = 0
    total_symmetry_failures = 0
    total_zero_failures = 0

    # ------------------------------------------------------------------------
    # 1. Exact binomial formula verification
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. EXACT FINITE-BINOMIAL MODE FORMULA")
    print("=" * 78)

    for k, ell in PAIRS:

        print()
        print(
            f"({k},{ell})"
        )

        for j in TEST_J:

            exact, predicted, residual = verify_formula(
                k,
                ell,
                j
            )

            status = (
                "PASS"
                if residual == 0
                else "FAIL"
            )

            print(
                f"  j={j}: {status}"
            )

            if residual != 0:

                total_formula_failures += 1

                print(
                    "    exact    =",
                    sp.factor(exact)
                )

                print(
                    "    predicted=",
                    sp.factor(predicted)
                )

                print(
                    "    residual =",
                    residual
                )

        # --------------------------------------------------------------------
        # Symmetry
        # --------------------------------------------------------------------

        R = direct_laurent(
            k,
            ell
        )

        support = []

        for term in sp.Add.make_args(
            sp.expand(R)
        ):

            power = int(
                sp.sympify(
                    term.as_powers_dict().get(u, 0)
                )
            )

            support.append(
                power
            )

        symmetry_failed = False

        for j in sorted(
            set(support)
        ):

            cj = expanded_laurent_coeff(
                R,
                j
            )

            cmj = expanded_laurent_coeff(
                R,
                -j
            )

            if sp.expand(
                cj - cmj
            ) != 0:

                symmetry_failed = True
                total_symmetry_failures += 1

                print(
                    f"  symmetry FAIL at j={j}"
                )

        if not symmetry_failed:

            print(
                "  Laurent symmetry = PASS"
            )

        # --------------------------------------------------------------------
        # Zero mode
        # --------------------------------------------------------------------

        c0 = expanded_laurent_coeff(
            R,
            0
        )

        expected0 = expected_zero_mode(
            k,
            ell
        )

        zero_residual = sp.factor(
            sp.expand(
                c0 - expected0
            )
        )

        print(
            "  zero mode:",
            sp.factor(c0)
        )

        print(
            "  expected:",
            sp.factor(expected0)
        )

        print(
            "  zero-mode status =",
            "PASS"
            if zero_residual == 0
            else "FAIL"
        )

        if zero_residual != 0:
            total_zero_failures += 1

    # ------------------------------------------------------------------------
    # 2. Explicit low-mode examples
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. LOW NONZERO MODE STRUCTURE")
    print("=" * 78)

    for k, ell in PAIRS:

        print_mode_table(
            k,
            ell,
            max_j=3
        )

        print()

    # ------------------------------------------------------------------------
    # 3. General coefficient shape
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. GENERAL COEFFICIENT SHAPE")
    print("=" * 78)

    print(
        "For the first kernel term,"
    )

    print(
        "  p^k(1+q)^ell"
    )

    print(
        "the coefficient of u^j is"
    )

    print(
        "  binom(ell, k-j) * t^(2k-j),"
    )

    print(
        "whenever 0 <= k-j <= ell."
    )

    print()

    print(
        "For the symmetric term"
    )

    print(
        "  q^k(1+p)^ell,"
    )

    print(
        "the coefficient is"
    )

    print(
        "  binom(ell, k+j) * t^(2k+j),"
    )

    print(
        "whenever 0 <= k+j <= ell."
    )

    print()

    print(
        "The ell-outer terms contribute the corresponding"
    )

    print(
        "negative binomial corrections."
    )

    # ------------------------------------------------------------------------
    # 4. N-only classification of C0
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ZERO-MODE CLOSED FORM")
    print("=" * 78)

    print(
        "Verified identity:"
    )

    print(
        "  C_0(t) = 2*binom(ell,k)*t^(2k)"
    )

    print(
        "           = 2*binom(ell,k)*N^k"
    )

    # ------------------------------------------------------------------------
    # 5. Final diagnostic
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "formula failures =",
        total_formula_failures
    )

    print(
        "symmetry failures =",
        total_symmetry_failures
    )

    print(
        "zero-mode failures =",
        total_zero_failures
    )

    if (
        total_formula_failures == 0
        and total_symmetry_failures == 0
        and total_zero_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The nonzero Laurent coefficients are now given"
        )

        print(
            "directly by finite binomial expressions."
        )

        print()
        print(
            "The zero mode is exactly:"
        )

        print(
            "    C_0 = 2*binom(ell,k)*N^k."
        )

        print()
        print(
            "The next step is Experiment 152:"
        )

        print(
            "derive the exact quotient-mode recurrence from"
        )

        print(
            "    F = (1+t(u+u^-1)) Q."
        )

    else:

        print(
            "STATUS = FAIL"
        )

    print("=" * 78)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print()
        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)

