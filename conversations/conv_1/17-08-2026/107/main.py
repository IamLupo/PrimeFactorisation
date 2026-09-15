#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 170R
PIECEWISE-EXACT COMPLETE k=1 COEFFICIENT LAW
==============================================================================

FIXES EXPERIMENT 169
--------------------
1. s=0 is handled explicitly:
       B=0
       C=1
       D=1

2. s=1 B has ONLY the first Chu-Vandermonde component:
       B = C(ell-1,2)

   The second hypergeometric component is absent.

3. s>=2 retains the two-term closed B formula.

GOAL
----
Verify the complete coefficient law:

    D_(1,ell,0) = 1

    D_(1,ell,1)
       = C(ell-1,2) - (ell-1)

    for s>=2,

    D_(1,ell,s)
       = (-1)^(s-1)
         [ C(ell-s,s+1)
           + C(ell-s+1,s+2) ]
       + (-1)^s
         [ C(ell-s-1,s)
           + C(ell-s-2,s-1) ]

Then reconstruct the complete polynomial exactly.

NO FLOATS
NO FITTING
NO INTERPOLATION
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Exact V coefficient
# ============================================================================

def v_coeff(r, m):
    if r < 0 or m < 0 or m > r // 2:
        return sp.Integer(0)

    first = sp.binomial(r - m, m)

    second = (
        sp.Integer(0)
        if m == 0
        else sp.binomial(r - m - 1, m - 1)
    )

    return sp.expand(
        (-1) ** (r - m) * (first + second)
    )


# ============================================================================
# Exact B branch
# ============================================================================

def B_term(ell, s, n):
    m = s - 1

    if m < 0:
        return sp.Integer(0)

    if m > n // 2:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (n - m)
        * sp.binomial(
            ell,
            n + 3
        )
        * (
            sp.binomial(
                n - m,
                m
            )
            +
            (
                sp.Integer(0)
                if m == 0
                else sp.binomial(
                    n - m - 1,
                    m - 1
                )
            )
        )
    )


def B_support(ell, s):
    if s == 0:
        return []

    lo = 2 * (s - 1)
    hi = ell - 3

    if lo > hi:
        return []

    return list(
        range(
            lo,
            hi + 1
        )
    )


def B_exact(ell, s):
    if s == 0:
        return sp.Integer(0)

    return sp.factor(
        sum(
            B_term(
                ell,
                s,
                n
            )
            for n in B_support(
                ell,
                s
            )
        )
    )


# ============================================================================
# Exact C branch
# ============================================================================

def C_term(ell, s, j):
    if s == 0:
        # Explicit boundary value
        return (
            sp.Integer(1)
            if j == ell
            else sp.Integer(0)
        )

    if j not in (
        ell - 1,
        ell
    ):
        return sp.Integer(0)

    r = j - 2

    if r < 0:
        return sp.Integer(0)

    e = 2 * ell - j

    numerator = e - r - 2

    if numerator % 2:
        return sp.Integer(0)

    m = s - numerator // 2

    if m < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        -v_coeff(
            r,
            m
        )
    )


def C_exact(ell, s):
    return sp.factor(
        sum(
            C_term(
                ell,
                s,
                j
            )
            for j in (
                ell - 1,
                ell
            )
        )
    )


# ============================================================================
# Correct piecewise closed forms
# ============================================================================

def B_closed(ell, s):

    # IMPORTANT: s=0
    if s == 0:
        return sp.Integer(0)

    # IMPORTANT: s=1
    if s == 1:
        return sp.binomial(
            ell - 1,
            2
        )

    # s>=2
    return sp.factor(
        (-1) ** (s - 1)
        * (
            sp.binomial(
                ell - s,
                s + 1
            )
            +
            sp.binomial(
                ell - s + 1,
                s + 2
            )
        )
    )


def C_closed(ell, s):

    # IMPORTANT: s=0
    if s == 0:
        return sp.Integer(1)

    # s>=1
    return sp.factor(
        (-1) ** s
        * (
            sp.binomial(
                ell - s - 1,
                s
            )
            +
            sp.binomial(
                ell - s - 2,
                s - 1
            )
        )
    )


def D_closed(ell, s):
    return sp.factor(
        B_closed(
            ell,
            s
        )
        +
        C_closed(
            ell,
            s
        )
    )


def D_exact(ell, s):
    return sp.factor(
        B_exact(
            ell,
            s
        )
        +
        C_exact(
            ell,
            s
        )
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 170R")
    print("PIECEWISE-EXACT COMPLETE k=1 COEFFICIENT LAW")
    print("=" * 78)
    print()

    TRAIN = [
        9, 11, 13, 15, 17,
        19, 21, 23
    ]

    HOLDOUT = [
        25, 27, 29, 31, 33
    ]

    failures = 0

    # ------------------------------------------------------------------------
    # 1. B certificate
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. B-BRANCH")
    print("=" * 78)

    for ell in TRAIN:

        print()
        print(
            f"ell={ell}"
        )

        max_s = (
            ell - 1
        ) // 2

        for s in range(
            0,
            max_s + 1
        ):

            exact = B_exact(
                ell,
                s
            )

            closed = B_closed(
                ell,
                s
            )

            residual = sp.factor(
                exact - closed
            )

            print(
                f"  s={s}: "
                f"exact={exact}, "
                f"closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 2. C certificate
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. C-BOUNDARY")
    print("=" * 78)

    for ell in TRAIN:

        print()
        print(
            f"ell={ell}"
        )

        max_s = (
            ell - 1
        ) // 2

        for s in range(
            0,
            max_s + 1
        ):

            exact = C_exact(
                ell,
                s
            )

            closed = C_closed(
                ell,
                s
            )

            residual = sp.factor(
                exact - closed
            )

            print(
                f"  s={s}: "
                f"exact={exact}, "
                f"closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 3. Complete D certificate
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. COMPLETE D COEFFICIENT")
    print("=" * 78)

    for ell in TRAIN:

        print()
        print(
            f"ell={ell}"
        )

        max_s = (
            ell - 1
        ) // 2

        for s in range(
            0,
            max_s + 1
        ):

            exact = D_exact(
                ell,
                s
            )

            closed = D_closed(
                ell,
                s
            )

            residual = sp.factor(
                exact - closed
            )

            print(
                f"  s={s}: "
                f"D={exact}, "
                f"closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 4. Full polynomial
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FULL k=1 POLYNOMIAL")
    print("=" * 78)

    N = sp.symbols(
        "N"
    )

    for ell in TRAIN:

        degree = (
            ell - 3
        ) // 2

        direct_poly = sp.factor(
            sum(
                D_exact(
                    ell,
                    s
                ) * N**s
                for s in range(
                    degree + 1
                )
            )
        )

        closed_poly = sp.factor(
            sum(
                D_closed(
                    ell,
                    s
                ) * N**s
                for s in range(
                    degree + 1
                )
            )
        )

        residual = sp.factor(
            direct_poly
            - closed_poly
        )

        print()
        print(
            f"ell={ell}"
        )

        print(
            "  D_direct =",
            direct_poly
        )

        print(
            "  D_closed =",
            closed_poly
        )

        print(
            "  residual =",
            residual
        )

        if residual != 0:
            failures += 1

    # ------------------------------------------------------------------------
    # 5. Holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FORWARD ELL HOLDOUT")
    print("=" * 78)

    holdout_failures = 0

    for ell in HOLDOUT:

        degree = (
            ell - 3
        ) // 2

        ok = True

        for s in range(
            0,
            degree + 1
        ):

            residual = sp.factor(
                D_exact(
                    ell,
                    s
                )
                -
                D_closed(
                    ell,
                    s
                )
            )

            if residual != 0:
                ok = False

        print(
            f"ell={ell}: "
            f"status={'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            holdout_failures += 1
            failures += 1

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "coefficient failures =",
        failures
    )

    print(
        "holdout failures =",
        holdout_failures
    )

    if failures == 0:

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The complete k=1 coefficient law is exact."
        )

        print()
        print(
            "The boundary cases s=0 and s=1 are now"
        )

        print(
            "handled separately, while s>=2 retains"
        )

        print(
            "the Chu-Vandermonde-derived formula."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "compress the resulting coefficient formula"
        )

        print(
            "and then generalize the B/C derivation from"
        )

        print(
            "k=1 to k=3 and arbitrary odd k."
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

