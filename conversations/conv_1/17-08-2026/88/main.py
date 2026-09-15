#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 152
EXACT QUOTIENT-MODE RECURRENCE FROM THE LAURENT KERNEL

GOAL
----
Experiment 151 gave the exact Laurent coefficients

    F_(k,ell)(t,u) = sum_j C_j(t) u^j

directly from the finite binomial kernel.

Since

    F_(k,ell) = (p+q+1) Q_(k,ell),

the substitution

    p = t u,
    q = t u^(-1)

gives

    F(t,u)
      = (1 + t(u+u^(-1))) Q(t,u).

Writing

    Q(t,u) = sum_j q_j(t) u^j,

we obtain the exact recurrence

    C_j = q_j + t q_(j-1) + t q_(j+1).

This experiment:

    1. derives q_j recursively from the finite Laurent support of F;
    2. verifies the result against the independently constructed quotient;
    3. proves the Laurent reconstruction identity exactly;
    4. extracts q_0, q_1, q_2, q_3;
    5. factors their t-powers and rewrites even powers as N=t^2;
    6. tests whether the first quotient modes have a simple N-envelope
       plus unavoidable sqrt(N) / branch dependence.

IMPORTANT
---------
We do NOT construct the C/D tensor.

We do NOT construct the full Newton tensor.

We do NOT search factor pairs.

The quotient-mode recurrence is derived only from

    F = (p+q+1)Q.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
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
# Exact finite kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Exact quotient in p
# ============================================================================

def quotient_Q_pq(k, ell):
    F = kernel_F(k, ell)

    PF = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )

    divisor = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    Q, R = sp.div(
        PF,
        divisor
    )

    rem = sp.expand(
        R.as_expr()
    )

    if rem != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = "
            f"{sp.factor(rem)}"
        )

    return sp.expand(
        Q.as_expr()
    )


# ============================================================================
# Direct Laurent substitutions
# ============================================================================

def laurent_F(k, ell):
    F = kernel_F(k, ell)

    return sp.expand(
        F.subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


def laurent_Q_direct(k, ell):
    Q = quotient_Q_pq(k, ell)

    return sp.expand(
        Q.subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


# ============================================================================
# Laurent support / coefficient
# ============================================================================

def laurent_support(expr):
    expr = sp.expand(expr)

    support = []

    for term in sp.Add.make_args(expr):

        power = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent: {term}"
            )

        support.append(
            int(power)
        )

    return sorted(
        set(support)
    )


def laurent_coeff(expr, j):
    expr = sp.expand(expr)

    support = laurent_support(expr)

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
# Explicit quotient-mode recurrence
# ============================================================================

def recover_q_modes_from_F(k, ell):
    """
    Recover the entire quotient Laurent sequence q_j from the finite
    kernel coefficients C_j.

    We use the highest-support boundary.

    If

        C_j = q_j + t q_(j-1) + t q_(j+1),

    and q_j has finite support [qmin, qmax], then

        qmax = Fmax - 1.

    At j=Fmax:

        C_Fmax = t q_(Fmax-1).

    This gives the highest q mode.

    We then descend:

        q_(j-1)
          = (C_j - q_j - t q_(j+1))/t.
    """

    F = laurent_F(
        k,
        ell
    )

    F_support = laurent_support(
        F
    )

    if not F_support:
        return {}

    fmin = min(F_support)
    fmax = max(F_support)

    qmax = fmax - 1
    qmin = fmin + 1

    q_modes = {}

    # Highest boundary.
    C_top = laurent_coeff(
        F,
        fmax
    )

    q_modes[qmax] = sp.cancel(
        C_top / t
    )

    # Descend from fmax-1 down to fmin+1.
    #
    # At each j:
    #
    # C_j = q_j + t q_(j-1) + t q_(j+1)
    #
    # therefore
    #
    # q_(j-1) = (C_j - q_j - t q_(j+1))/t.
    #
    for j in range(
        fmax - 1,
        fmin,
        -1
    ):

        q_j = q_modes.get(
            j,
            sp.Integer(0)
        )

        q_j_plus_1 = q_modes.get(
            j + 1,
            sp.Integer(0)
        )

        C_j = laurent_coeff(
            F,
            j
        )

        q_prev = sp.cancel(
            (
                C_j
                - q_j
                - t*q_j_plus_1
            ) / t
        )

        q_modes[j - 1] = sp.expand(
            q_prev
        )

    # Restrict to expected finite support.
    q_modes = {
        j: sp.factor(expr)
        for j, expr in sorted(
            q_modes.items()
        )
        if expr != 0
    }

    return q_modes


# ============================================================================
# Exact recurrence residual
# ============================================================================

def recurrence_residuals(
    k,
    ell,
    q_modes
):
    F = laurent_F(
        k,
        ell
    )

    F_support = laurent_support(
        F
    )

    if not F_support:
        return []

    lo = min(F_support)
    hi = max(F_support)

    failures = []

    for j in range(
        lo,
        hi + 1
    ):

        qj = q_modes.get(
            j,
            sp.Integer(0)
        )

        qm = q_modes.get(
            j - 1,
            sp.Integer(0)
        )

        qp = q_modes.get(
            j + 1,
            sp.Integer(0)
        )

        Cj = laurent_coeff(
            F,
            j
        )

        residual = sp.factor(
            sp.expand(
                Cj
                - qj
                - t*qm
                - t*qp
            )
        )

        if residual != 0:
            failures.append(
                (
                    j,
                    residual
                )
            )

    return failures


# ============================================================================
# Laurent reconstruction from q-modes
# ============================================================================

def reconstruct_F_from_q(k, ell, q_modes):
    Q = sp.Integer(0)

    for j, coeff in q_modes.items():
        Q += coeff*u**j

    return sp.expand(
        (
            1
            + t*(u + u**-1)
        )
        * Q
    )


# ============================================================================
# Difference from direct quotient
# ============================================================================

def compare_with_direct_quotient(
    k,
    ell,
    q_modes
):
    direct = laurent_Q_direct(
        k,
        ell
    )

    reconstructed = sp.Integer(0)

    for j, coeff in q_modes.items():
        reconstructed += coeff*u**j

    residual = sp.factor(
        sp.expand(
            direct - reconstructed
        )
    )

    return residual


# ============================================================================
# Factor t-power
# ============================================================================

def factor_common_t_power(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return {
            "power": None,
            "reduced": sp.Integer(0)
        }

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ
    )

    powers = [
        int(mon[0])
        for mon in poly.monoms()
    ]

    common = min(
        powers
    )

    reduced = sp.factor(
        sp.expand(
            expr / t**common
        )
    )

    return {
        "power": common,
        "reduced": reduced,
    }


# ============================================================================
# Rewrite even polynomial in t as polynomial in N
# ============================================================================

def rewrite_even_t_as_N(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return sp.Integer(0)

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ
    )

    result = sp.Integer(0)

    for (power,), coeff in poly.terms():

        power = int(power)

        if power % 2 != 0:
            return None

        result += (
            coeff
            * N**(power // 2)
        )

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Print low modes
# ============================================================================

def print_low_modes(
    k,
    ell,
    q_modes,
    max_j=3
):
    print(
        f"({k},{ell})"
    )

    for j in range(
        0,
        max_j + 1
    ):

        expr = q_modes.get(
            j,
            sp.Integer(0)
        )

        if expr == 0:
            continue

        shape = factor_common_t_power(
            expr
        )

        print(
            f"  q_{j}(t) = "
            f"{sp.factor(expr)}"
        )

        print(
            f"      common t-power = "
            f"{shape['power']}"
        )

        print(
            f"      reduced = "
            f"{shape['reduced']}"
        )

        n_form = rewrite_even_t_as_N(
            shape["reduced"]
        )

        if n_form is not None:

            print(
                f"      reduced N-form = "
                f"{n_form}"
            )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 152")
    print("EXACT QUOTIENT-MODE RECURRENCE FROM THE LAURENT KERNEL")
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

    total_recurrence_failures = 0
    total_direct_failures = 0
    total_reconstruction_failures = 0

    # ------------------------------------------------------------------------
    # 1. Recover quotient Laurent modes.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. QUOTIENT-MODE RECOVERY")
    print("=" * 78)

    all_modes = {}

    for k, ell in PAIRS:

        print()
        print(
            f"({k},{ell})"
        )

        q_modes = recover_q_modes_from_F(
            k,
            ell
        )

        all_modes[
            (k, ell)
        ] = q_modes

        print(
            "  q-support =",
            sorted(q_modes)
        )

        recurrence_failures = recurrence_residuals(
            k,
            ell,
            q_modes
        )

        print(
            "  recurrence failures =",
            len(recurrence_failures)
        )

        if recurrence_failures:

            total_recurrence_failures += len(
                recurrence_failures
            )

            for j, residual in recurrence_failures[:5]:

                print(
                    f"    j={j}: "
                    f"{residual}"
                )

        reconstruction = reconstruct_F_from_q(
            k,
            ell,
            q_modes
        )

        direct_F = laurent_F(
            k,
            ell
        )

        reconstruction_residual = sp.factor(
            sp.expand(
                reconstruction - direct_F
            )
        )

        print(
            "  direct F reconstruction =",
            "PASS"
            if reconstruction_residual == 0
            else "FAIL"
        )

        if reconstruction_residual != 0:

            total_reconstruction_failures += 1

            print(
                "    residual =",
                reconstruction_residual
            )

        direct_quotient_residual = (
            compare_with_direct_quotient(
                k,
                ell,
                q_modes
            )
        )

        print(
            "  independent quotient comparison =",
            "PASS"
            if direct_quotient_residual == 0
            else "FAIL"
        )

        if direct_quotient_residual != 0:

            total_direct_failures += 1

            print(
                "    residual =",
                direct_quotient_residual
            )

    # ------------------------------------------------------------------------
    # 2. Low-mode structure.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. LOW QUOTIENT-MODE STRUCTURE")
    print("=" * 78)

    for k, ell in PAIRS:

        print_low_modes(
            k,
            ell,
            all_modes[
                (k, ell)
            ],
            max_j=3
        )

        print()

    # ------------------------------------------------------------------------
    # 3. Branch information classification.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. MODE INFORMATION CLASSIFICATION")
    print("=" * 78)

    zero_mode_count = 0
    pure_N_reduced_count = 0
    nonzero_mode_count = 0

    for k, ell in PAIRS:

        q_modes = all_modes[
            (k, ell)
        ]

        for j, coeff in q_modes.items():

            if coeff == 0:
                continue

            if j == 0:
                zero_mode_count += 1

                shape = factor_common_t_power(
                    coeff
                )

                n_form = rewrite_even_t_as_N(
                    shape["reduced"]
                )

                if n_form is not None:
                    pure_N_reduced_count += 1

            else:
                nonzero_mode_count += 1

    print(
        "nonzero q-modes =",
        nonzero_mode_count
    )

    print(
        "zero q-modes checked =",
        zero_mode_count
    )

    print(
        "zero/reduced-N forms =",
        pure_N_reduced_count
    )

    # ------------------------------------------------------------------------
    # 4. Candidate q0 recurrence.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. q0 N-ONLY TEST")
    print("=" * 78)

    for k, ell in PAIRS:

        q0 = all_modes[
            (k, ell)
        ].get(
            0,
            sp.Integer(0)
        )

        print(
            f"({k},{ell}) q_0(t) = "
            f"{sp.factor(q0)}"
        )

        shape = factor_common_t_power(
            q0
        )

        n_form = rewrite_even_t_as_N(
            shape["reduced"]
        )

        print(
            "    common t-power =",
            shape["power"]
        )

        if n_form is not None:

            print(
                "    reduced N-form =",
                n_form
            )

        else:

            print(
                "    reduced N-form = NONE"
            )

    # ------------------------------------------------------------------------
    # 5. Final diagnostic.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "quotient recurrence failures =",
        total_recurrence_failures
    )

    print(
        "F reconstruction failures =",
        total_reconstruction_failures
    )

    print(
        "direct quotient comparison failures =",
        total_direct_failures
    )

    if (
        total_recurrence_failures == 0
        and total_reconstruction_failures == 0
        and total_direct_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The quotient Laurent coefficients are recovered"
        )

        print(
            "exactly from the finite kernel through"
        )

        print()
        print(
            "    C_j = q_j + t q_(j-1) + t q_(j+1)."
        )

        print()
        print(
            "The recovered modes agree with the independently"
        )

        print(
            "constructed quotient."
        )

        print()
        print(
            "The next target is now to derive a closed form for"
        )

        print(
            "the first informative quotient modes q_1,q_2,q_3."
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

