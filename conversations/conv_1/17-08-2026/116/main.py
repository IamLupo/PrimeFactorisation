#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 180
GENERAL ODD-k C-BRANCH: SYMBOLIC + EXACT VANDERMONDE COLLAPSE
==============================================================================

Fixes Experiment 179:

The previous script used the numerical helper Cbin() during a symbolic
calculation. That helper contains Python comparisons such as

    r < 0

which fail when r is symbolic.

This version uses:

    Cbin_num(n,r)     for numerical exact evaluation
    Cbin_sym(n,r)     for symbolic derivation

The C-law under test is

    C(k,ell,s)
      = (-1)^s [
          binomial(ell-s-1, s-k+1)
          + binomial(ell-s-2, s-k)
        ]

The symbolic derivation keeps binomial() symbolic and performs all
support-sensitive work only in the numerical verification stage.
==============================================================================

"""

import sympy as sp
import sys


# ============================================================================
# NUMERICAL BINOMIAL
# ============================================================================

def Cbin_num(n, r):
    """
    Exact binomial with zero outside the standard support.
    Both n and r are assumed numerical integers.
    """
    n = int(n)
    r = int(r)

    if r < 0 or n < 0 or r > n:
        return sp.Integer(0)

    return sp.binomial(n, r)


# ============================================================================
# SYMBOLIC BINOMIAL
# ============================================================================

def Cbin_sym(n, r):
    """
    Pure symbolic binomial.

    IMPORTANT:
    Do not perform Python inequalities on symbolic arguments.
    """
    return sp.binomial(n, r)


# ============================================================================
# V KERNEL -- NUMERICAL
# ============================================================================

def V_num(r, m):

    if r < 0 or m < 0 or 2*m > r:
        return sp.Integer(0)

    return sp.expand(
        (-1)**(r-m) *
        (
            Cbin_num(r-m, m)
            +
            Cbin_num(r-m-1, m-1)
        )
    )


# ============================================================================
# EXACT C TERM
# ============================================================================

def C_term(k, ell, s, j):

    if j < ell-k or j > ell:
        return sp.Integer(0)

    a = j - (ell-k)

    r = j - 2*k

    m_num = (
        2*s
        - (2*ell-j)
        + r
        + 2
    )

    if m_num % 2:
        return sp.Integer(0)

    m = m_num // 2

    return sp.expand(
        -Cbin_num(k, a) * V_num(r, m)
    )


def C_exact(k, ell, s):

    return sp.factor(
        sum(
            C_term(k, ell, s, j)
            for j in range(
                ell-k,
                ell+1
            )
        )
    )


# ============================================================================
# SYMBOLIC VANDERMONDE DERIVATION
# ============================================================================

def symbolic_C_derivation():

    k, ell, s, a = sp.symbols(
        "k ell s a",
        integer=True,
        nonnegative=True
    )

    print()
    print("="*78)
    print("1. SYMBOLIC C-COLLAPSE")
    print("="*78)

    print()
    print("Set")
    print("  r = ell - 3*k + a")
    print("  m = s - 2*k + a + 1")

    r_minus_m = sp.expand(
        (ell - 3*k + a)
        - (s - 2*k + a + 1)
    )

    print()
    print("Then")
    print("  r-m =", r_minus_m)

    expected_q = ell-k-s-1

    if sp.simplify(
        r_minus_m-expected_q
    ) != 0:
        print("FAIL: r-m simplification")
        return 1

    q = expected_q
    t = s-2*k+1

    print()
    print("Define")
    print("  q =", q)
    print("  t =", t)

    # ------------------------------------------------------------------------
    # FIRST V PIECE
    #
    # sum_a C(k,a) C(q,t+a)
    #
    # Reindex b=k-a:
    #
    # C(q,t+k-b)
    #
    # Vandermonde:
    #
    # sum_b C(k,b) C(q,t+k-b)
    #   = C(q+k,t+k)
    # ------------------------------------------------------------------------

    first_sum = sp.Function(
        "SumFirst"
    )

    first_closed = Cbin_sym(
        q+k,
        t+k
    )

    first_target = Cbin_sym(
        ell-s-1,
        s-k+1
    )

    print()
    print("First V piece:")
    print("  sum_a C(k,a) C(q,t+a)")
    print("  -> C(q+k,t+k)")

    first_residual = sp.simplify(
        first_closed-first_target
    )

    print(
        "  symbolic residual =",
        first_residual
    )

    # ------------------------------------------------------------------------
    # SECOND V PIECE
    #
    # sum_a C(k,a) C(q-1,t+a-1)
    #
    # Reindex b=k-a:
    #
    # C(q-1,t+k-1-b)
    #
    # Vandermonde:
    #
    # = C(q+k-1,t+k-1)
    # ------------------------------------------------------------------------

    second_closed = Cbin_sym(
        q+k-1,
        t+k-1
    )

    second_target = Cbin_sym(
        ell-s-2,
        s-k
    )

    print()
    print("Second V piece:")
    print("  sum_a C(k,a) C(q-1,t+a-1)")
    print("  -> C(q+k-1,t+k-1)")

    second_residual = sp.simplify(
        second_closed-second_target
    )

    print(
        "  symbolic residual =",
        second_residual
    )

    # The symbolic comparison above should be identically zero after
    # substitution. SymPy can preserve equivalent binomial forms, so also
    # compare their expanded arguments directly.

    arg_residual_1 = (
        sp.expand(q+k)
        - sp.expand(ell-s-1)
    )

    idx_residual_1 = (
        sp.expand(t+k)
        - sp.expand(s-k+1)
    )

    arg_residual_2 = (
        sp.expand(q+k-1)
        - sp.expand(ell-s-2)
    )

    idx_residual_2 = (
        sp.expand(t+k-1)
        - sp.expand(s-k)
    )

    print()
    print("Argument checks:")
    print("  first upper residual  =", arg_residual_1)
    print("  first lower residual  =", idx_residual_1)
    print("  second upper residual =", arg_residual_2)
    print("  second lower residual =", idx_residual_2)

    failures = 0

    if arg_residual_1 != 0:
        failures += 1

    if idx_residual_1 != 0:
        failures += 1

    if arg_residual_2 != 0:
        failures += 1

    if idx_residual_2 != 0:
        failures += 1

    print()
    if failures == 0:
        print("SYMBOLIC COLLAPSE = PASS")
    else:
        print(
            "SYMBOLIC COLLAPSE = FAIL:",
            failures
        )

    return failures


# ============================================================================
# CLOSED C FORMULA
# ============================================================================

def C_closed(k, ell, s):

    return sp.factor(
        (-1)**s *
        (
            Cbin_num(
                ell-s-1,
                s-k+1
            )
            +
            Cbin_num(
                ell-s-2,
                s-k
            )
        )
    )


# ============================================================================
# C VERIFICATION
# ============================================================================

def test_C():

    failures = 0

    print()
    print("="*78)
    print("2. EXACT C VERIFICATION")
    print("="*78)

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+2,3),
            52,
            2
        ):

            max_s = (ell-1)//2

            for s in range(
                max_s+1
            ):

                exact = C_exact(
                    k, ell, s
                )

                closed = C_closed(
                    k, ell, s
                )

                residual = sp.factor(
                    exact-closed
                )

                if residual != 0:

                    failures += 1

                    print(
                        "FAIL:",
                        f"k={k}, ell={ell}, s={s}",
                        f"exact={exact}",
                        f"closed={closed}",
                        f"residual={residual}"
                    )

    print()
    print(
        "C failures =",
        failures
    )

    return failures


# ============================================================================
# REPRESENTATIVE CHECKS
# ============================================================================

def representatives():

    print()
    print("="*78)
    print("3. REPRESENTATIVE C VALUES")
    print("="*78)

    samples = [
        (1,17,0),
        (1,17,1),
        (1,17,2),
        (1,17,5),
        (3,17,2),
        (3,17,3),
        (3,17,6),
        (5,25,4),
        (5,25,5),
        (5,25,8),
        (7,31,6),
        (7,31,7),
        (9,35,8),
        (9,35,9),
    ]

    failures = 0

    for k, ell, s in samples:

        exact = C_exact(
            k, ell, s
        )

        closed = C_closed(
            k, ell, s
        )

        residual = sp.factor(
            exact-closed
        )

        print()
        print(
            f"(k,ell,s)=({k},{ell},{s})"
        )
        print(
            f"  exact  = {exact}"
        )
        print(
            f"  closed = {closed}"
        )
        print(
            f"  residual = {residual}"
        )

        if residual != 0:
            failures += 1

    print()
    print(
        "representative failures =",
        failures
    )

    return failures


# ============================================================================
# B FORMULA
# ============================================================================

def B_term(k, ell, s, j):

    m = s-k
    r = j-2*k

    if m < 0:
        return sp.Integer(0)

    if r < 0:
        return sp.Integer(0)

    return sp.expand(
        Cbin_num(ell,j+k)
        * V_num(r,m)
    )


def B_exact(k, ell, s):

    return sp.factor(
        sum(
            B_term(k,ell,s,j)
            for j in range(ell+1)
        )
    )


# ============================================================================
# ESTABLISHED PIECEWISE B LAW
# ============================================================================

def B_closed(k, ell, s):

    # s < k: vanishes
    if s < k:
        return sp.Integer(0)

    # edge s = k:
    # established exact theorem
    if s == k:
        return Cbin_num(
            ell-1,
            3*k-1
        )

    # For s >= k+1, use the already verified shifted
    # two-term Vandermonde representation.
    #
    # This is intentionally kept separate from the C experiment.

    L = ell-2*s-k

    if L < 0:
        return sp.Integer(0)

    pref = (
        (-1)**(s-k)
        * Cbin_num(
            ell,
            2*s+k
        )
    )

    # First 2F1 collapse
    first = sp.Rational(1,1)

    for n in range(L+1):
        first *= sp.Rational(
            s+k+1+n,
            2*s+k+1+n
        )

    # Second 2F1 collapse
    second = sp.Rational(1,1)

    for n in range(L+1):
        second *= sp.Rational(
            s+k+2+n,
            2*s+k+1+n
        )

    return sp.factor(
        pref*(first+second)
    )


# ============================================================================
# COMPLETE D
# ============================================================================

def D_exact(k,ell,s):

    return sp.factor(
        B_exact(k,ell,s)
        +
        C_exact(k,ell,s)
    )


def D_closed(k,ell,s):

    return sp.factor(
        B_closed(k,ell,s)
        +
        C_closed(k,ell,s)
    )


# ============================================================================
# COMPLETE RECONSTRUCTION
# ============================================================================

def test_D():

    failures = 0

    print()
    print("="*78)
    print("4. COMPLETE B+C RECONSTRUCTION")
    print("="*78)

    for k in [1,3,5,7,9]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+2,3),
            42,
            2
        ):

            max_s = (ell-1)//2

            for s in range(
                max_s+1
            ):

                exact = D_exact(
                    k,ell,s
                )

                closed = D_closed(
                    k,ell,s
                )

                residual = sp.factor(
                    exact-closed
                )

                if residual != 0:

                    failures += 1

                    print(
                        "FAIL:",
                        f"k={k}, ell={ell}, s={s}",
                        f"exact={exact}",
                        f"closed={closed}",
                        f"residual={residual}"
                    )

    print()
    print(
        "D failures =",
        failures
    )

    return failures


# ============================================================================
# FRESH HOLDOUT
# ============================================================================

def test_holdout():

    failures = 0

    print()
    print("="*78)
    print("5. FRESH HOLDOUT")
    print("="*78)

    for k in [1,3,5,7,9]:

        for ell in [43,45,47,49]:

            status = "PASS"

            for s in range(
                (ell-1)//2 + 1
            ):

                exact = D_exact(
                    k,ell,s
                )

                closed = D_closed(
                    k,ell,s
                )

                if sp.factor(
                    exact-closed
                ) != 0:

                    failures += 1
                    status = "FAIL"

                    print(
                        "FAIL:",
                        f"k={k}, ell={ell}, s={s}",
                        f"exact={exact}",
                        f"closed={closed}"
                    )

                    break

            print(
                f"(k={k},ell={ell}) "
                f"status={status}"
            )

    return failures


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print("KAPPA EXPERIMENT 180")
    print("GENERAL ODD-k C-BRANCH: EXACT VANDERMONDE COLLAPSE")
    print("="*78)

    symbolic_failures = (
        symbolic_C_derivation()
    )

    c_failures = test_C()

    representative_failures = (
        representatives()
    )

    d_failures = test_D()

    holdout_failures = (
        test_holdout()
    )

    print()
    print("="*78)
    print("FINAL DIAGNOSTIC")
    print("="*78)

    print(
        "symbolic failures       =",
        symbolic_failures
    )

    print(
        "C failures              =",
        c_failures
    )

    print(
        "representative failures =",
        representative_failures
    )

    print(
        "D failures              =",
        d_failures
    )

    print(
        "holdout failures        =",
        holdout_failures
    )

    print()

    if (
        symbolic_failures == 0
        and c_failures == 0
        and representative_failures == 0
        and d_failures == 0
        and holdout_failures == 0
    ):

        print("STATUS = PASS")
        print()
        print(
            "The general C Vandermonde formula is"
            " exact and the B+C reconstruction passes."
        )

    else:

        print("STATUS = PARTIAL/FAIL")
        print()
        print(
            "The symbolic C derivation is separated"
            " correctly from numerical support logic."
        )
        print(
            "Inspect the first failing layer before"
            " attempting a broader odd-k theorem."
        )

    print("="*78)


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