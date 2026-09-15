#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 158R
UNIVERSAL GREEN-FUNCTION WEIGHTS -- FIXED LAURENT-t HANDLING

FIX
---
The previous script failed on Laurent coefficients such as

    1/t,
    -1/t^2,

because SymPy Poly(expr, t) expects ordinary polynomial expressions.

This version explicitly separates the minimum t-exponent using
as_powers_dict(), then converts the remaining even t-polynomial to N=t^2.

TARGET
------
Verify the universal recurrence

    W_0 = 1
    W_1 = -1
    W_(r+1) = -W_r - N W_(r-1)

and its closed finite-binomial form

    W_r(N)
      = sum_{m=0}^{floor(r/2)}
          (-1)^(r-m) binom(r-m,m) N^m.

Then verify that the exact antisymmetric quotient source weights satisfy

    source_weight_r
      = W_r(N) / t^(r+1)

with r = m - 2k,

and derive

    D_(k,ell)(N)
      = (q_(2k-1)-q_(2k+1)) / (t N^(k-1))

directly from the finite binomial kernel.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
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
N = sp.symbols("N")


# ============================================================================
# Exact kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Exact quotient
# ============================================================================

def quotient_Q_pq(k, ell):
    F = kernel_F(k, ell)

    PF = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )

    PD = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    Q, R = sp.div(PF, PD)

    rem = sp.expand(R.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = {sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Laurent helpers
# ============================================================================

def laurent_F(k, ell):
    return sp.expand(
        kernel_F(k, ell).subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


def laurent_Q(k, ell):
    return sp.expand(
        quotient_Q_pq(k, ell).subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


def laurent_support(expr):
    support = []

    for term in sp.Add.make_args(sp.expand(expr)):
        power = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent in u: {term}"
            )

        support.append(int(power))

    return sorted(set(support))


def laurent_coeff(expr, j):
    support = laurent_support(expr)

    if not support:
        return sp.Integer(0)

    shift = max(0, -min(support))

    poly = sp.Poly(
        sp.expand(expr * u**shift),
        u,
        domain=sp.QQ.frac_field(t),
    )

    target = j + shift

    if target < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(u**target)
    )


# ============================================================================
# Robust Laurent-t decomposition
# ============================================================================

def split_t_laurent(expr):
    """
    Return

        expr = t^base * reduced,

    where reduced has no negative t powers.

    This works for Laurent polynomials such as

        1/t,
        (2*N-1)/t^3,
        t^4 + 2*t^2 + 1.
    """

    expr = sp.expand(expr)

    if expr == 0:
        return 0, sp.Integer(0)

    min_power = None

    for term in sp.Add.make_args(expr):
        power = sp.sympify(
            term.as_powers_dict().get(t, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"non-integer t exponent: {term}"
            )

        power = int(power)

        if min_power is None or power < min_power:
            min_power = power

    reduced = sp.expand(
        expr / t**min_power
    )

    return min_power, sp.expand(reduced)


def even_t_to_N(expr):
    """
    Convert an ordinary even polynomial in t to a polynomial in N=t^2.

    Returns None if odd t-powers remain.
    """

    expr = sp.expand(expr)

    if expr == 0:
        return sp.Integer(0)

    # First ensure this is an ordinary polynomial in t.
    for term in sp.Add.make_args(expr):
        power = sp.sympify(
            term.as_powers_dict().get(t, 0)
        )

        if not power.is_Integer:
            return None

        if int(power) < 0:
            return None

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ,
    )

    result = sp.Integer(0)

    for (power,), coeff in poly.terms():
        power = int(power)

        if power % 2 != 0:
            return None

        result += coeff * N**(power // 2)

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Universal Green-function recurrence
# ============================================================================

def W_recursive(max_r):
    W = {
        0: sp.Integer(1),
        1: sp.Integer(-1),
    }

    for r in range(1, max_r):
        W[r + 1] = sp.factor(
            sp.expand(
                -W[r] - N*W[r - 1]
            )
        )

    return W


def W_closed(r):
    return sp.factor(
        sp.expand(
            sum(
                (-1)**(r-m)
                * sp.binomial(r-m, m)
                * N**m
                for m in range(
                    0,
                    r//2 + 1
                )
            )
        )
    )


# ============================================================================
# Source representation of q_j
# ============================================================================

def source_weights_for_q(k, ell, target_j):
    """
    Express q_target_j as a finite exact linear combination of C_m.

    Boundary:
        q_(ell-1) = C_ell/t

    Recurrence:
        q_(j-1)
          = C_j/t - q_j/t - q_(j+1).
    """

    source = {}

    source[ell - 1] = {
        ell: sp.cancel(1/t)
    }

    for j in range(
        ell - 1,
        target_j,
        -1
    ):
        current = source.get(
            j,
            {}
        )

        next_one = source.get(
            j + 1,
            {}
        )

        new_map = {
            j: sp.cancel(1/t)
        }

        for m, coeff in current.items():
            new_map[m] = sp.expand(
                new_map.get(m, 0)
                - coeff/t
            )

        for m, coeff in next_one.items():
            new_map[m] = sp.expand(
                new_map.get(m, 0)
                - coeff
            )

        source[j - 1] = {
            m: sp.factor(c)
            for m, c in new_map.items()
            if sp.expand(c) != 0
        }

    return {
        m: sp.factor(c)
        for m, c in source.get(
            target_j,
            {}
        ).items()
        if sp.expand(c) != 0
    }


def antisymmetric_source_weights(k, ell):
    left = source_weights_for_q(
        k,
        ell,
        2*k - 1
    )

    right = source_weights_for_q(
        k,
        ell,
        2*k + 1
    )

    combined = {}

    for m, coeff in left.items():
        combined[m] = sp.expand(
            combined.get(m, 0) + coeff
        )

    for m, coeff in right.items():
        combined[m] = sp.expand(
            combined.get(m, 0) - coeff
        )

    return {
        m: sp.factor(c)
        for m, c in combined.items()
        if sp.expand(c) != 0
    }


# ============================================================================
# Normalize source weight against W_r
# ============================================================================

def normalized_source_weight(k, m, coeff):
    """
    For r = m - 2k, test the conjecture

        coeff = W_r(N) / t^(r+1).

    Equivalent normalized quantity:

        coeff * t^(r+1) = W_r(N).
    """

    r = m - 2*k

    if r < 0:
        return {
            "r": r,
            "valid": False,
            "reason": "negative Green index",
        }

    normalized = sp.expand(
        coeff * t**(r + 1)
    )

    Nform = even_t_to_N(
        normalized
    )

    if Nform is None:
        return {
            "r": r,
            "valid": False,
            "reason": "normalized source is not N-only",
        }

    expected = W_closed(r)

    residual = sp.factor(
        sp.expand(
            Nform - expected
        )
    )

    return {
        "r": r,
        "valid": residual == 0,
        "normalized": Nform,
        "expected": expected,
        "residual": residual,
    }


# ============================================================================
# Closed kernel coefficient
# ============================================================================

def kernel_mode_closed(k, ell, j):
    # p^k (1+q)^ell
    r1 = k - j

    if 0 <= r1 <= ell:
        c1 = (
            sp.binomial(ell, r1)
            * t**(k + r1)
        )
    else:
        c1 = sp.Integer(0)

    # q^k (1+p)^ell
    r2 = j + k

    if 0 <= r2 <= ell:
        c2 = (
            sp.binomial(ell, r2)
            * t**(k + r2)
        )
    else:
        c2 = sp.Integer(0)

    # -p^ell (1+q)^k
    r3 = ell - j

    if 0 <= r3 <= k:
        c3 = (
            -sp.binomial(k, r3)
            * t**(ell + r3)
        )
    else:
        c3 = sp.Integer(0)

    # -q^ell (1+p)^k
    r4 = ell + j

    if 0 <= r4 <= k:
        c4 = (
            -sp.binomial(k, r4)
            * t**(ell + r4)
        )
    else:
        c4 = sp.Integer(0)

    return sp.expand(
        c1 + c2 + c3 + c4
    )


# ============================================================================
# Recover D from finite source
# ============================================================================

def D_from_source(k, ell):
    weights = antisymmetric_source_weights(
        k,
        ell
    )

    raw = sp.Integer(0)

    for m, weight in weights.items():
        raw += (
            weight
            * kernel_mode_closed(
                k,
                ell,
                m
            )
        )

    raw = sp.cancel(
        raw
        / (
            t
            * t**(2*(k-1))
        )
    )

    base, reduced = split_t_laurent(
        raw
    )

    # D should have no residual t-power.
    if base != 0:
        raise ArithmeticError(
            f"D source retained t-power {base} "
            f"for ({k},{ell})"
        )

    Nform = even_t_to_N(
        reduced
    )

    if Nform is None:
        raise ArithmeticError(
            f"D source retained odd t-powers "
            f"for ({k},{ell})"
        )

    return sp.factor(
        Nform
    )


# ============================================================================
# Direct D from quotient
# ============================================================================

def D_direct(k, ell):
    Q = laurent_Q(k, ell)

    qm = laurent_coeff(
        Q,
        2*k - 1
    )

    qp = laurent_coeff(
        Q,
        2*k + 1
    )

    raw = sp.cancel(
        (qm - qp)
        / (
            t
            * t**(2*(k-1))
        )
    )

    base, reduced = split_t_laurent(
        raw
    )

    if base != 0:
        raise ArithmeticError(
            f"D direct retained t-power {base} "
            f"for ({k},{ell})"
        )

    Nform = even_t_to_N(
        reduced
    )

    if Nform is None:
        raise ArithmeticError(
            f"D direct retained odd t-powers "
            f"for ({k},{ell})"
        )

    return sp.factor(
        Nform
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 158R")
    print("UNIVERSAL GREEN-FUNCTION WEIGHTS -- FIXED")
    print("=" * 78)
    print()

    TRAIN = [
        (1,3),
        (1,5),
        (1,7),
        (1,9),
        (1,11),
        (1,13),
        (3,7),
        (3,9),
        (3,11),
        (3,13),
        (5,11),
        (5,13),
    ]

    FORWARD = [
        (1,15),
        (1,17),
        (3,15),
        (3,17),
        (5,15),
        (5,17),
    ]

    # ------------------------------------------------------------------------
    # 1. Green recurrence
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. UNIVERSAL GREEN POLYNOMIAL")
    print("=" * 78)

    W_MAX = 20
    W = W_recursive(W_MAX)

    green_failures = 0

    for r in range(W_MAX + 1):

        rec = W[r]
        closed = W_closed(r)

        residual = sp.factor(
            sp.expand(rec - closed)
        )

        print(
            f"r={r:2d} "
            f"W_r={rec} "
            f"closed={closed} "
            f"residual={residual}"
        )

        if residual != 0:
            green_failures += 1

    # ------------------------------------------------------------------------
    # 2. Source-weight matching
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SOURCE-WEIGHT MATCHING")
    print("=" * 78)

    source_failures = 0

    for k, ell in TRAIN:

        if ell < 2*k + 1:
            continue

        weights = antisymmetric_source_weights(
            k,
            ell
        )

        print()
        print(
            f"({k},{ell})"
        )

        for m in sorted(weights):

            result = normalized_source_weight(
                k,
                m,
                weights[m]
            )

            r = result["r"]

            if not result["valid"]:

                source_failures += 1

                print(
                    f"  C_{m}: FAIL "
                    f"r={r} "
                    f"reason={result.get('reason')}"
                )

                if "normalized" in result:
                    print(
                        "      normalized =",
                        result["normalized"]
                    )

                    print(
                        "      expected =",
                        result["expected"]
                    )

                    print(
                        "      residual =",
                        result["residual"]
                    )

                continue

            print(
                f"  C_{m}: "
                f"r={r} "
                f"normalized={result['normalized']} "
                f"W_r={result['expected']} "
                f"PASS"
            )

    # ------------------------------------------------------------------------
    # 3. D from source
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. D FROM UNIVERSAL GREEN SOURCE")
    print("=" * 78)

    D_failures = 0
    boundary_failures = 0

    for k, ell in TRAIN:

        if ell < 2*k + 1:
            continue

        direct = D_direct(
            k,
            ell
        )

        source = D_from_source(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - source
            )
        )

        boundary = sp.factor(
            direct.subs(N, 0)
        )

        print()
        print(
            f"({k},{ell})"
        )

        print(
            "  D_direct =",
            direct
        )

        print(
            "  D_source =",
            source
        )

        print(
            "  residual =",
            residual
        )

        print(
            "  D(0) =",
            boundary
        )

        if residual != 0:
            D_failures += 1

        if boundary != 1:
            boundary_failures += 1

    # ------------------------------------------------------------------------
    # 4. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FORWARD ELL HOLDOUT")
    print("=" * 78)

    holdout_failures = 0

    for k, ell in FORWARD:

        direct = D_direct(
            k,
            ell
        )

        source = D_from_source(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - source
            )
        )

        boundary = sp.factor(
            direct.subs(N, 0)
        )

        ok = (
            residual == 0
            and boundary == 1
        )

        print(
            f"({k},{ell}) "
            f"status={'PASS' if ok else 'FAIL'}"
        )

        print(
            "  D(N) =",
            direct
        )

        if not ok:
            holdout_failures += 1

    # ------------------------------------------------------------------------
    # 5. Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "Green closed-form failures =",
        green_failures
    )

    print(
        "source-weight failures =",
        source_failures
    )

    print(
        "D source residual failures =",
        D_failures
    )

    print(
        "D(0) failures =",
        boundary_failures
    )

    print(
        "forward failures =",
        holdout_failures
    )

    if (
        green_failures == 0
        and source_failures == 0
        and D_failures == 0
        and boundary_failures == 0
        and holdout_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The universal Green-function source weights"
        )

        print(
            "are exactly"
        )

        print(
            "    W_0=1,"
        )

        print(
            "    W_1=-1,"
        )

        print(
            "    W_(r+1)=-W_r-N W_(r-1),"
        )

        print(
            "with the finite closed form"
        )

        print(
            "    W_r="
        )

        print(
            "      sum_m (-1)^(r-m)"
        )

        print(
            "          binom(r-m,m) N^m."
        )

        print()
        print(
            "The same universal kernel generates D_(k,ell)."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "compress the D convolution to a single finite"
        )

        print(
            "binomial/hypergeometric sum."
        )

    else:

        print()
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