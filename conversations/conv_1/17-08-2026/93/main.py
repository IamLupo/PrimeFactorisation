#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 157
EXACT FINITE-BINOMIAL FORMULA FOR THE ANTISYMMETRIC MODE

TARGET
------
Experiment 156 established

    q_(2k-1) = t R(N)
    q_(2k)   = N^k P(N)
    q_(2k+1) = t T(N)

and therefore

    R(N) - T(N)
        = N^(k-1) D_(k,ell)(N),

with

    D_(k,ell)(0) = 1

on all tested training and forward pairs.

This experiment derives D directly from the finite Laurent kernel.

The quotient recurrence is

    C_j = q_j + t q_(j-1) + t q_(j+1).

Because Q has finite Laurent support, the recurrence can be solved from
the upper boundary. We derive the antisymmetric combination

    q_(2k-1) - q_(2k+1)

as an EXACT finite linear combination of kernel coefficients C_j.

Then we substitute the exact finite-binomial formula for C_j from
Experiment 151.

The result is normalized as

    D_(k,ell)(N)
      =
    (q_(2k-1)-q_(2k+1)) / (t N^(k-1)).

The experiment checks:

    1. exact finite-source representation;
    2. exact equality with direct quotient extraction;
    3. exact N-only reduction;
    4. D(0)=1;
    5. forward-ell holdout;
    6. coefficient sequence of D;
    7. whether the coefficient source is a finite binomial sum.

IMPORTANT
---------
No polynomial interpolation is used.

The coefficients C_j are computed directly from binomial coefficients.

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
# Laurent expressions
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

    for term in sp.Add.make_args(
        sp.expand(expr)
    ):
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
# Direct finite-binomial Laurent coefficient
# ============================================================================

def binomial_contribution(
    outer_power,
    inner_power,
    j,
    sign,
):
    """
    Contribution from

        sign * (t u)^outer (1+t u^-1)^inner.

    u^j requires

        outer-r = j,

    hence

        r = outer-j.
    """

    r = outer_power - j

    if r < 0 or r > inner_power:
        return sp.Integer(0)

    return sp.expand(
        sign
        * sp.binomial(inner_power, r)
        * t**(outer_power + r)
    )


def kernel_mode_closed(k, ell, j):
    """
    Exact [u^j] coefficient of F_(k,ell)(tu,t/u).
    """

    # p^k (1+q)^ell
    c1 = binomial_contribution(
        k,
        ell,
        j,
        +1,
    )

    # q^k (1+p)^ell
    r2 = j + k

    if 0 <= r2 <= ell:
        c2 = sp.expand(
            sp.binomial(ell, r2)
            * t**(k + r2)
        )
    else:
        c2 = sp.Integer(0)

    # -p^ell (1+q)^k
    c3 = binomial_contribution(
        ell,
        k,
        j,
        -1,
    )

    # -q^ell (1+p)^k
    r4 = j + ell

    if 0 <= r4 <= k:
        c4 = sp.expand(
            -sp.binomial(k, r4)
            * t**(ell + r4)
        )
    else:
        c4 = sp.Integer(0)

    return sp.expand(
        c1 + c2 + c3 + c4
    )


# ============================================================================
# Solve quotient recurrence symbolically from upper boundary
# ============================================================================

def recover_q_from_closed_C(k, ell):
    """
    Recover q_j using ONLY the closed finite-binomial C_j.

    Q support:

        -(ell-1) <= j <= ell-1.

    Start with

        C_ell = t q_(ell-1),

    then descend.
    """

    fmax = ell

    q_modes = {}

    q_modes[ell - 1] = sp.cancel(
        kernel_mode_closed(
            k,
            ell,
            ell
        ) / t
    )

    for j in range(
        ell - 1,
        -ell,
        -1
    ):
        Cj = kernel_mode_closed(
            k,
            ell,
            j
        )

        qj = q_modes.get(
            j,
            sp.Integer(0)
        )

        qjp = q_modes.get(
            j + 1,
            sp.Integer(0)
        )

        q_modes[j - 1] = sp.expand(
            sp.cancel(
                (
                    Cj
                    - qj
                    - t*qjp
                ) / t
            )
        )

    return {
        j: sp.factor(expr)
        for j, expr in q_modes.items()
        if expr != 0
    }


# ============================================================================
# Direct quotient mode
# ============================================================================

def direct_q_mode(k, ell, j):
    return sp.factor(
        laurent_coeff(
            laurent_Q(
                k,
                ell
            ),
            j
        )
    )


# ============================================================================
# Convert even polynomial in t to N
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

        if power % 2 != 0:
            return None

        result += coeff * N**(power // 2)

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Direct D
# ============================================================================

def D_direct(k, ell):
    qm = direct_q_mode(
        k,
        ell,
        2*k - 1
    )

    qp = direct_q_mode(
        k,
        ell,
        2*k + 1
    )

    raw = sp.expand(
        sp.cancel(
            (qm - qp)
            / (
                t
                * N**(k - 1)
            )
        )
    )

    # Replace N by t^2 before checking parity.
    raw_t = sp.expand(
        raw.subs(
            N,
            t**2
        )
    )

    return even_t_to_N(
        raw_t
    )


# ============================================================================
# D from finite-binomial C-source
# ============================================================================

def D_from_closed_C(k, ell):
    q_modes = recover_q_from_closed_C(
        k,
        ell
    )

    qm = q_modes.get(
        2*k - 1,
        sp.Integer(0)
    )

    qp = q_modes.get(
        2*k + 1,
        sp.Integer(0)
    )

    raw_t = sp.expand(
        sp.cancel(
            (qm - qp)
            / (
                t
                * t**(2*(k-1))
            )
        )
    )

    return even_t_to_N(
        raw_t
    )


# ============================================================================
# Explicit finite source weights
# ============================================================================

def derive_source_weights(
    k,
    ell
):
    """
    Derive the linear combination of C_j coefficients that produces
    q_(2k-1) - q_(2k+1).

    This is done symbolically by propagating abstract source symbols.

    q_j is linear in the C_m. We represent it as a dictionary

        j -> {m: coefficient}.

    This gives an exact finite source formula independent of the
    explicit kernel coefficients.
    """

    source = {}

    # q_(ell-1) = C_ell / t
    source[ell - 1] = {
        ell: sp.cancel(1/t)
    }

    for j in range(
        ell - 1,
        -ell,
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

        new_map = {}

        # q_(j-1) =
        #   C_j/t
        #   - q_j/t
        #   - q_(j+1)

        new_map[j] = sp.cancel(
            1/t
        )

        for m, coeff in current.items():
            new_map[m] = sp.expand(
                new_map.get(
                    m,
                    0
                )
                - coeff/t
            )

        for m, coeff in next_one.items():
            new_map[m] = sp.expand(
                new_map.get(
                    m,
                    0
                )
                - coeff
            )

        source[j - 1] = {
            m: sp.factor(c)
            for m, c in new_map.items()
            if c != 0
        }

    a = source.get(
        2*k - 1,
        {}
    )

    b = source.get(
        2*k + 1,
        {}
    )

    combined = {}

    for m, coeff in a.items():
        combined[m] = sp.expand(
            combined.get(m, 0)
            + coeff
        )

    for m, coeff in b.items():
        combined[m] = sp.expand(
            combined.get(m, 0)
            - coeff
        )

    return {
        m: sp.factor(c)
        for m, c in combined.items()
        if sp.expand(c) != 0
    }


# ============================================================================
# Coefficient list
# ============================================================================

def coefficient_list(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return []

    return [
        sp.factor(
            poly.coeff_monomial(N**r)
        )
        for r in range(
            poly.degree() + 1
        )
    ]


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 157")
    print("EXACT FINITE-BINOMIAL FORMULA FOR D_(k,ell)")
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

    total_source_failures = 0
    total_N_failures = 0
    total_boundary_failures = 0
    total_holdout_failures = 0

    # ------------------------------------------------------------------------
    # 1. Exact D source derivation
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. FINITE-BINOMIAL SOURCE DERIVATION")
    print("=" * 78)

    for k, ell in TRAIN:

        if ell < 2*k + 1:
            continue

        print()
        print(
            f"({k},{ell})"
        )

        direct = D_direct(
            k,
            ell
        )

        source_based = D_from_closed_C(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - source_based
            )
        )

        print(
            "  D_direct(N) =",
            direct
        )

        print(
            "  D_source(N) =",
            source_based
        )

        print(
            "  source residual =",
            residual
        )

        if residual != 0:
            total_source_failures += 1

        if direct is None:
            total_N_failures += 1

        else:
            at_zero = sp.factor(
                direct.subs(N, 0)
            )

            print(
                "  D(0) =",
                at_zero
            )

            if at_zero != 1:
                total_boundary_failures += 1

            print(
                "  coefficients =",
                coefficient_list(direct)
            )

    # ------------------------------------------------------------------------
    # 2. Display abstract finite source weights
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ABSTRACT SOURCE WEIGHTS")
    print("=" * 78)

    # Use one representative pair because the weights depend only on
    # support geometry and k,ell.
    for k, ell in [
        (1, 9),
        (3, 11),
        (5, 13),
    ]:

        weights = derive_source_weights(
            k,
            ell
        )

        print()
        print(
            f"({k},{ell})"
        )

        print(
            "  q_(2k-1)-q_(2k+1) ="
        )

        for m in sorted(
            weights
        ):
            print(
                f"    C_{m}: "
                f"{weights[m]}"
            )

    # ------------------------------------------------------------------------
    # 3. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    for k, ell in FORWARD:

        direct = D_direct(
            k,
            ell
        )

        source_based = D_from_closed_C(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - source_based
            )
        )

        ok = (
            residual == 0
            and direct is not None
            and sp.expand(
                direct.subs(N, 0) - 1
            ) == 0
        )

        print(
            f"({k},{ell}) "
            f"status={'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            total_holdout_failures += 1

        print(
            "  D(N) =",
            direct
        )

    # ------------------------------------------------------------------------
    # 4. Leading/constant diagnostic
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. D NORMALIZATION DIAGNOSTIC")
    print("=" * 78)

    for k, ell in TRAIN:

        if ell < 2*k + 1:
            continue

        D = D_direct(
            k,
            ell
        )

        if D is None:
            continue

        poly = sp.Poly(
            D,
            N,
            domain=sp.QQ,
        )

        print(
            f"({k},{ell}) "
            f"degree={poly.degree()} "
            f"constant={sp.factor(poly.coeff_monomial(1))} "
            f"leading={sp.factor(poly.LC())}"
        )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "finite-source residual failures =",
        total_source_failures
    )

    print(
        "N-only conversion failures =",
        total_N_failures
    )

    print(
        "D(0)!=1 failures =",
        total_boundary_failures
    )

    print(
        "forward holdout failures =",
        total_holdout_failures
    )

    if (
        total_source_failures == 0
        and total_N_failures == 0
        and total_boundary_failures == 0
        and total_holdout_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The antisymmetric local invariant"
        )

        print(
            "    D_(k,ell)(N)"
        )

        print(
            "has been recovered directly from the finite"
        )

        print(
            "binomial Laurent source."
        )

        print()
        print(
            "In particular:"
        )

        print(
            "    R-T = N^(k-1) D_(k,ell)(N),"
        )

        print(
            "    D_(k,ell)(0) = 1."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "simplify the finite source weights to a closed"
        )

        print(
            "binomial/hypergeometric expression for D_(k,ell)."
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

