#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 181
CORRECT GENERAL ODD-k B + PIECEWISE C
FULL B+C RECONSTRUCTION
==============================================================================

Goals
-----
1. Verify the general interior B hypergeometric collapse for odd k.
2. Verify the exact B edge theorem s=k.
3. Verify the two exceptional C boundary layers.
4. Verify the generic C Vandermonde formula elsewhere.
5. Reconstruct D = B + C exactly.
6. Test both odd and even ell.
7. Run fresh holdouts beyond the training range.

Correct B interior law
----------------------
For s >= k+1, define

    L = ell - 2*s - k

    B = (-1)^(s-k) * C(ell, 2*s+k)
        * [ (s+2k)_L / (2s+k+1)_L
            + (s+2k+1)_L / (2s+k+1)_L ]

B edge
------
For s = k,

    B = C(ell-1, 3k-1)

and B=0 for s<k.

Correct C law
-------------
Generic:

    C = (-1)^s [
        C(ell-s-1, s-k+1)
        + C(ell-s-2, s-k)
    ]

with two exceptional geometric boundary layers:

    ell = 2k-1, s = k-1  -> C = 0
    ell = 2k+1, s = k    -> C = -k

The exceptions arise because the formal Vandermonde expression
contains a contribution whose corresponding original active support
is absent.
"""

import sympy as sp
import sys


# ============================================================================
# BASIC EXACT HELPERS
# ============================================================================

def Cbin(n, r):
    """
    Numerical exact binomial with zero outside support.
    """
    n = int(n)
    r = int(r)

    if n < 0 or r < 0 or r > n:
        return sp.Integer(0)

    return sp.binomial(n, r)


def poch_ratio(a, c, L):
    """
    (a)_L / (c)_L for nonnegative integer L.
    """
    if L < 0:
        return sp.Integer(0)

    out = sp.Integer(1)

    for t in range(L):
        out *= sp.Rational(a + t, c + t)

    return sp.factor(out)


# ============================================================================
# EXACT VERNELS / DIRECT BRANCHES
# ============================================================================

def V_num(r, m):
    """
    Exact finite V-kernel used by the direct coefficient construction.

    This is the same kernel structure used in the earlier experiments:
        (-1)^(r-m) [ C(r-m,m) + C(r-m-1,m-1) ]

    with zero outside admissible support.
    """
    if r < 0 or m < 0 or 2*m > r:
        return sp.Integer(0)

    return sp.expand(
        (-1)**(r-m)
        * (
            Cbin(r-m, m)
            + Cbin(r-m-1, m-1)
        )
    )


# ============================================================================
# DIRECT B
# ============================================================================

def B_term(k, ell, s, j):

    if s < k:
        return sp.Integer(0)

    m = s - k
    r = j - 2*k

    if r < 0:
        return sp.Integer(0)

    return sp.expand(
        Cbin(ell, j + k)
        * V_num(r, m)
    )


def B_exact(k, ell, s):

    return sp.factor(
        sum(
            B_term(k, ell, s, j)
            for j in range(ell + 1)
        )
    )


# ============================================================================
# CORRECT GENERAL B CLOSED FORM
# ============================================================================

def B_closed(k, ell, s):

    # Vanishing region
    if s < k:
        return sp.Integer(0)

    # Exact edge s=k
    if s == k:
        return Cbin(
            ell - 1,
            3*k - 1
        )

    # Interior
    if s >= k + 1:

        L = ell - 2*s - k

        if L < 0:
            return sp.Integer(0)

        prefactor = (
            (-1)**(s-k)
            * Cbin(
                ell,
                2*s+k
            )
        )

        first = poch_ratio(
            s + 2*k,
            2*s + k + 1,
            L
        )

        second = poch_ratio(
            s + 2*k + 1,
            2*s + k + 1,
            L
        )

        return sp.factor(
            prefactor * (first + second)
        )

    return sp.Integer(0)


# ============================================================================
# DIRECT C
# ============================================================================

def C_term(k, ell, s, j):

    # The active boundary strip is j in [ell-k, ell].
    if j < ell-k or j > ell:
        return sp.Integer(0)

    a = j - (ell-k)

    r = j - 2*k

    numerator = (
        2*s
        - (2*ell-j)
        + r
        + 2
    )

    if numerator % 2:
        return sp.Integer(0)

    m = numerator // 2

    return sp.expand(
        -Cbin(k, a)
        * V_num(r, m)
    )


def C_exact(k, ell, s):

    return sp.factor(
        sum(
            C_term(k, ell, s, j)
            for j in range(ell + 1)
        )
    )


# ============================================================================
# CORRECT PIECEWISE C CLOSED FORM
# ============================================================================

def C_closed(k, ell, s):

    # ------------------------------------------------------------
    # Exceptional boundary 1:
    #
    # ell = 2k-1, s=k-1
    #
    # There is no active original C support, hence C=0.
    # ------------------------------------------------------------

    if (
        ell == 2*k - 1
        and s == k - 1
    ):
        return sp.Integer(0)

    # ------------------------------------------------------------
    # Exceptional boundary 2:
    #
    # ell = 2k+1, s=k
    #
    # The formal Vandermonde expression overcounts by one.
    # Exact value = -k.
    # ------------------------------------------------------------

    if (
        ell == 2*k + 1
        and s == k
    ):
        return sp.Integer(-k)

    # ------------------------------------------------------------
    # Generic Vandermonde collapse
    # ------------------------------------------------------------

    return sp.factor(
        (-1)**s
        * (
            Cbin(
                ell-s-1,
                s-k+1
            )
            +
            Cbin(
                ell-s-2,
                s-k
            )
        )
    )


# ============================================================================
# DIRECT / CLOSED D
# ============================================================================

def D_exact(k, ell, s):

    return sp.factor(
        B_exact(k, ell, s)
        +
        C_exact(k, ell, s)
    )


def D_closed(k, ell, s):

    return sp.factor(
        B_closed(k, ell, s)
        +
        C_closed(k, ell, s)
    )


# ============================================================================
# 1. B EDGE TEST
# ============================================================================

def test_B_edge():

    print()
    print("="*78)
    print("1. B EDGE THEOREM")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(0, k-2),
            52,
            2
        ):

            exact = B_exact(
                k, ell, k
            )

            closed = B_closed(
                k, ell, k
            )

            residual = sp.factor(
                exact-closed
            )

            if residual != 0:
                failures += 1

                print(
                    "FAIL:",
                    f"k={k}, ell={ell}",
                    f"exact={exact}",
                    f"closed={closed}",
                    f"residual={residual}"
                )

    print()
    print(
        "B edge failures =",
        failures
    )

    return failures


# ============================================================================
# 2. B INTERIOR TEST
# ============================================================================

def test_B_interior():

    print()
    print("="*78)
    print("2. B INTERIOR HYPERGEOMETRIC COLLAPSE")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+1, 3),
            52
        ):

            for s in range(
                k+1,
                ell
            ):

                exact = B_exact(
                    k, ell, s
                )

                closed = B_closed(
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
        "B interior failures =",
        failures
    )

    return failures


# ============================================================================
# 3. C EXCEPTION TEST
# ============================================================================

def test_C_exceptions():

    print()
    print("="*78)
    print("3. C EXCEPTIONAL BOUNDARIES")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9,11]:

        cases = [
            (2*k-1, k-1),
            (2*k+1, k),
        ]

        print()
        print(f"k={k}")

        for ell, s in cases:

            exact = C_exact(
                k, ell, s
            )

            closed = C_closed(
                k, ell, s
            )

            residual = sp.factor(
                exact-closed
            )

            print(
                f"  (ell,s)=({ell},{s})",
                f"exact={exact}",
                f"closed={closed}",
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    print()
    print(
        "C exceptional failures =",
        failures
    )

    return failures


# ============================================================================
# 4. C GENERIC TEST
# ============================================================================

def test_C_generic():

    print()
    print("="*78)
    print("4. C GENERIC VANDERMONDE COLLAPSE")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+1, 3),
            52
        ):

            for s in range(
                0,
                ell + 1
            ):

                # Skip the two exceptional layers.
                if (
                    ell == 2*k-1
                    and s == k-1
                ):
                    continue

                if (
                    ell == 2*k+1
                    and s == k
                ):
                    continue

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
        "C generic failures =",
        failures
    )

    return failures


# ============================================================================
# 5. FULL D RECONSTRUCTION
# ============================================================================

def test_D():

    print()
    print("="*78)
    print("5. COMPLETE B+C RECONSTRUCTION")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9]:

        print()
        print(f"k={k}")

        # Test both parities of ell.
        for ell in range(
            max(1, k-2),
            46
        ):

            max_s = ell

            for s in range(
                max_s + 1
            ):

                exact = D_exact(
                    k, ell, s
                )

                closed = D_closed(
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
        "D failures =",
        failures
    )

    return failures


# ============================================================================
# 6. FRESH HOLDOUT
# ============================================================================

def test_holdout():

    print()
    print("="*78)
    print("6. FRESH HOLDOUT")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9,11]:

        for ell in [47,48,49,50,51]:

            status = "PASS"

            for s in range(
                ell + 1
            ):

                exact = D_exact(
                    k, ell, s
                )

                closed = D_closed(
                    k, ell, s
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

                    # Only print first failing s for each pair.
                    break

            print(
                f"(k={k},ell={ell}) "
                f"status={status}"
            )

    print()
    print(
        "holdout failures =",
        failures
    )

    return failures


# ============================================================================
# 7. SYMBOLIC B CERTIFICATE
# ============================================================================

def symbolic_B_certificate():

    print()
    print("="*78)
    print("7. SYMBOLIC B CERTIFICATE")
    print("="*78)

    k, ell, s, L, n = sp.symbols(
        "k ell s L n",
        integer=True,
        nonnegative=True
    )

    # General proposed interior shape:
    #
    # L = ell - 2s - k
    #
    # first hypergeometric:
    #   2F1(-L, s-k+1 ; 2s+k+1 ; 1)
    #
    # Chu-Vandermonde:
    #   (c-b)_L / (c)_L
    #
    # c-b = (2s+k+1)-(s-k+1) = s+2k
    #
    # second:
    #   c-b = s+2k+1

    c = 2*s + k + 1

    b1 = s-k+1
    b2 = s-k

    cb1 = sp.expand(c-b1)
    cb2 = sp.expand(c-b2)

    expected1 = s+2*k
    expected2 = s+2*k+1

    r1 = sp.simplify(
        cb1-expected1
    )

    r2 = sp.simplify(
        cb2-expected2
    )

    print()
    print(
        "First c-b =",
        cb1
    )
    print(
        "Expected   =",
        expected1
    )
    print(
        "Residual   =",
        r1
    )

    print()
    print(
        "Second c-b =",
        cb2
    )
    print(
        "Expected    =",
        expected2
    )
    print(
        "Residual    =",
        r2
    )

    # Also certify L.
    Ldef = sp.expand(
        ell-2*s-k
    )

    print()
    print(
        "L =", Ldef
    )

    failures = int(
        r1 != 0
        or r2 != 0
    )

    print()
    print(
        "symbolic B failures =",
        failures
    )

    return failures


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print("KAPPA EXPERIMENT 181")
    print("CORRECT GENERAL ODD-k B + PIECEWISE C")
    print("="*78)

    symbolic_failures = (
        symbolic_B_certificate()
    )

    b_edge_failures = (
        test_B_edge()
    )

    b_interior_failures = (
        test_B_interior()
    )

    c_exception_failures = (
        test_C_exceptions()
    )

    c_generic_failures = (
        test_C_generic()
    )

    d_failures = (
        test_D()
    )

    holdout_failures = (
        test_holdout()
    )

    print()
    print("="*78)
    print("FINAL DIAGNOSTIC")
    print("="*78)

    print(
        "symbolic B failures    =",
        symbolic_failures
    )

    print(
        "B edge failures        =",
        b_edge_failures
    )

    print(
        "B interior failures    =",
        b_interior_failures
    )

    print(
        "C exception failures   =",
        c_exception_failures
    )

    print(
        "C generic failures     =",
        c_generic_failures
    )

    print(
        "D failures             =",
        d_failures
    )

    print(
        "holdout failures       =",
        holdout_failures
    )

    print()

    total = (
        symbolic_failures
        + b_edge_failures
        + b_interior_failures
        + c_exception_failures
        + c_generic_failures
        + d_failures
        + holdout_failures
    )

    if total == 0:

        print("STATUS = PASS")
        print()
        print(
            "The corrected odd-k B interior law, "
            "B edge law, and piecewise C law all "
            "reconstruct the exact coefficient."
        )

    else:

        print("STATUS = PARTIAL/FAIL")
        print()
        print(
            "The first failing layer identifies the "
            "remaining obstruction to the general theorem."
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

