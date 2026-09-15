#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 156
ANTISYMMETRIC THREE-MODE RECURRENCE
EXACT FORMULA FOR R(N) - T(N)

From Experiment 155:

    q_(2k-1) = t R(N)
    q_(2k)   = N^k P(N)
    q_(2k+1) = t T(N)

and

    C_j = q_j + t q_(j-1) + t q_(j+1).

The central equation gives R+T.

This experiment derives information on R-T by combining the two neighboring
recurrence equations at j=2k-1 and j=2k+1.

The goal is to determine whether R-T has a direct closed form, preferably
in terms of finite binomial coefficients or already-known edge polynomials.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


p, q = sp.symbols("p q")
t, u = sp.symbols("t u")
N = sp.symbols("N")


# ============================================================================
# Kernel and quotient
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


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
    out = []

    for term in sp.Add.make_args(sp.expand(expr)):
        power = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent: {term}"
            )

        out.append(int(power))

    return sorted(set(out))


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
# Convert even t-polynomials to N
# ============================================================================

def even_t_to_N(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return sp.Integer(0)

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ,
    )

    result = sp.Integer(0)

    for (power,), coeff in poly.terms():
        power = int(power)

        if power % 2:
            return None

        result += coeff * N**(power // 2)

    return sp.factor(sp.expand(result))


# ============================================================================
# Extract the three central modes
# ============================================================================

def three_modes(k, ell):
    Q = laurent_Q(k, ell)
    j = 2*k

    qm = sp.factor(laurent_coeff(Q, j-1))
    q0 = sp.factor(laurent_coeff(Q, j))
    qp = sp.factor(laurent_coeff(Q, j+1))

    return qm, q0, qp


def normalize_three(k, qm, q0, qp):
    R = even_t_to_N(
        sp.cancel(qm / t)
    )

    P = even_t_to_N(
        sp.cancel(q0 / t**(2*k))
    )

    T = even_t_to_N(
        sp.cancel(qp / t)
    )

    return R, P, T


# ============================================================================
# Kernel coefficients at neighboring indices
# ============================================================================

def kernel_C(k, ell, j):
    return sp.factor(
        laurent_coeff(
            laurent_F(k, ell),
            j
        )
    )


def kernel_C_N(k, ell, j):
    C = kernel_C(k, ell, j)

    # Determine whether the natural t-power is odd/even.
    poly = sp.Poly(
        sp.expand(C),
        t,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return sp.Integer(0)

    powers = [int(m[0]) for m in poly.monoms()]
    base = min(powers)

    normalized = sp.expand(C / t**base)

    Nform = even_t_to_N(normalized)

    return base, Nform


# ============================================================================
# Exact antisymmetric recurrence identities
# ============================================================================

def antisymmetric_identity(k, ell):
    """
    At j0 = 2k-1:

        C_(j0) = q_(j0) + t q_(j0-1) + t q_(j0+1)

                  = tR + t q_(2k-2) + t*(N^k P)

    At j1 = 2k+1:

        C_(j1) = t*(N^k P) + t q_(2k)?? 
                 + neighboring terms.

    We compute these identities directly and then form differences.
    """

    Q = laurent_Q(k, ell)

    j0 = 2*k - 1
    j1 = 2*k + 1

    C0 = laurent_coeff(
        laurent_F(k, ell),
        j0
    )

    C1 = laurent_coeff(
        laurent_F(k, ell),
        j1
    )

    # All quotient modes required for the two equations.
    q_m2 = laurent_coeff(Q, 2*k-2)
    q_m1 = laurent_coeff(Q, 2*k-1)
    q_0  = laurent_coeff(Q, 2*k)
    q_p1 = laurent_coeff(Q, 2*k+1)
    q_p2 = laurent_coeff(Q, 2*k+2)

    eq_left = sp.expand(
        C0
        - q_m1
        - t*q_m2
        - t*q_0
    )

    eq_right = sp.expand(
        C1
        - q_p1
        - t*q_0
        - t*q_p2
    )

    return (
        sp.factor(eq_left),
        sp.factor(eq_right)
    )


# ============================================================================
# Direct R-T extraction
# ============================================================================

def direct_R_minus_T(k, ell):
    qm, q0, qp = three_modes(k, ell)

    R, P, T = normalize_three(
        k,
        qm,
        q0,
        qp
    )

    if R is None or T is None:
        return None

    return sp.factor(
        sp.expand(R - T)
    )


# ============================================================================
# Divisibility / factor probes
# ============================================================================

def factor_profile(expr):
    if expr == 0:
        return {
            "factor": sp.Integer(0),
            "constant": sp.Integer(0)
        }

    return {
        "factor": sp.factor(expr),
        "at_N0": sp.factor(expr.subs(N, 0))
    }


# ============================================================================
# Search for relation to central P
# ============================================================================

def relation_to_P(RminusT, P):
    """
    Report simple exact combinations already suggested by the data:

        R-T
        (R-T)/P
        (R-T)/(P+1)
        (R-T)/(2P+...)
    """

    candidates = {
        "R-T": RminusT,
        "(R-T)/P":
            sp.cancel(RminusT / P)
            if P != 0 else None,
        "(R-T)/(P+1)":
            sp.cancel(RminusT / (P+1)),
        "R-T + P":
            sp.expand(RminusT + P),
        "R-T - P":
            sp.expand(RminusT - P),
    }

    return {
        key: (
            None
            if value is None
            else sp.factor(value)
        )
        for key, value in candidates.items()
    }


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 156")
    print("ANTISYMMETRIC THREE-MODE RECURRENCE")
    print("EXACT FORMULA FOR R(N)-T(N)")
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

    failures = 0

    # ------------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. ANTISYMMETRIC THREE-MODE CERTIFICATE")
    print("=" * 78)

    records = []

    for k, ell in TRAIN:

        qm, q0, qp = three_modes(k, ell)

        R, P, T = normalize_three(
            k,
            qm,
            q0,
            qp
        )

        if R is None or P is None or T is None:
            print(
                f"({k},{ell}) STRUCTURE FAIL"
            )
            failures += 1
            continue

        diff = sp.factor(
            sp.expand(R - T)
        )

        sm = sp.factor(
            sp.expand(R + T)
        )

        # Verify the two neighboring recurrence equations exactly.
        left, right = antisymmetric_identity(
            k,
            ell
        )

        if left != 0 or right != 0:
            failures += 1

        print()
        print(
            f"({k},{ell})"
        )

        print(
            "  R(N) =",
            R
        )

        print(
            "  P(N) =",
            P
        )

        print(
            "  T(N) =",
            T
        )

        print(
            "  R+T =",
            sm
        )

        print(
            "  R-T =",
            diff
        )

        print(
            "  neighboring recurrence =",
            "PASS"
            if left == 0 and right == 0
            else "FAIL"
        )

        print(
            "  factor profile =",
            factor_profile(diff)
        )

        relations = relation_to_P(
            diff,
            P
        )

        print(
            "  relation diagnostics:"
        )

        for name, value in relations.items():
            print(
                f"    {name} = {value}"
            )

        records.append(
            {
                "k": k,
                "ell": ell,
                "R": R,
                "P": P,
                "T": T,
                "sum": sm,
                "diff": diff,
            }
        )

    # ------------------------------------------------------------------------
    # Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FORWARD ELL HOLDOUT")
    print("=" * 78)

    forward_failures = 0

    for k, ell in FORWARD:

        qm, q0, qp = three_modes(
            k,
            ell
        )

        R, P, T = normalize_three(
            k,
            qm,
            q0,
            qp
        )

        if R is None or P is None or T is None:
            forward_failures += 1
            print(
                f"({k},{ell}) STRUCTURE FAIL"
            )
            continue

        left, right = antisymmetric_identity(
            k,
            ell
        )

        ok = (
            left == 0
            and right == 0
        )

        diff = sp.factor(
            sp.expand(R - T)
        )

        print()
        print(
            f"({k},{ell}) status="
            + ("PASS" if ok else "FAIL")
        )

        print(
            "  R-T =",
            diff
        )

        if not ok:
            forward_failures += 1

    # ------------------------------------------------------------------------
    # 3. Look for a universal normalization
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. R-T NORMALIZATION TEST")
    print("=" * 78)

    for rec in records:

        k = rec["k"]
        ell = rec["ell"]
        diff = rec["diff"]

        # Expected universal powers suggested by the observed data.
        #
        # The natural first normalization is N^(k-1).
        #
        # We only report it; we do not claim a law from fitting.
        #
        base = N**max(k-1, 0)

        normalized = sp.factor(
            sp.cancel(diff / base)
        )

        print()
        print(
            f"({k},{ell})"
        )

        print(
            "  (R-T)/N^(k-1) =",
            normalized
        )

        print(
            "  value at N=0 =",
            sp.factor(
                normalized.subs(N, 0)
            )
        )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "training antisymmetric recurrence failures =",
        failures
    )

    print(
        "forward failures =",
        forward_failures
    )

    if (
        failures == 0
        and forward_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The neighboring Laurent equations provide an exact"
        )

        print(
            "algebraic system for R and T."
        )

        print()
        print(
            "R+T is already determined by the central equation."
        )

        print(
            "The remaining local invariant is R-T."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "derive a direct finite-binomial formula for R-T."
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

