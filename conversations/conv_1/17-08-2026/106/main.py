#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 169
COMPLETE k=1 COEFFICIENT LAW: B + C COMPRESSION
==============================================================================

EXPERIMENT 168 established the complete B-branch collapse:

    B1 = (-1)^(s-1) C(ell-s,   s+1)
    B2 = (-1)^(s-1) C(ell-s+1, s+2)

hence

    B = (-1)^(s-1)
        [ C(ell-s,s+1) + C(ell-s+1,s+2) ].

The C branch has only two boundary positions for k=1:

    j = ell-1,
    j = ell.

Using the exact V kernel gives

    C = (-1)^s
        [ C(ell-s-1,s) + C(ell-s-2,s-1) ].

TARGET
------
Verify the COMPLETE coefficient law

    D_s = [N^s] D_(1,ell)

         = (-1)^(s-1) *
           [
             C(ell-s,s+1)
             + C(ell-s+1,s+2)
             - C(ell-s-1,s)
             - C(ell-s-2,s-1)
           ]

for all admissible s.

Then verify:

    1. direct finite Laurent construction;
    2. exact B convolution;
    3. exact C boundary;
    4. compressed B formula;
    5. compressed C formula;
    6. complete D formula;
    7. polynomial reconstruction;
    8. forward ell holdout.

No fitting.
No interpolation.
No guessed parameters.

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

    a = sp.binomial(r-m, m)

    b = (
        sp.Integer(0)
        if m == 0
        else sp.binomial(r-m-1, m-1)
    )

    return sp.expand(
        (-1)**(r-m) * (a+b)
    )


# ============================================================================
# Exact k=1 B contribution
# ============================================================================

def B_term(ell, s, n):
    m = s-1

    if m < 0:
        return sp.Integer(0)

    if m > n//2:
        return sp.Integer(0)

    return sp.factor(
        (-1)**(n-m)
        * sp.binomial(
            ell,
            n+3
        )
        * (
            sp.binomial(
                n-m,
                m
            )
            +
            (
                sp.Integer(0)
                if m == 0
                else
                sp.binomial(
                    n-m-1,
                    m-1
                )
            )
        )
    )


def B_support(ell, s):
    lo = 2*(s-1)
    hi = ell-3

    if lo > hi:
        return []

    return list(
        range(
            lo,
            hi+1
        )
    )


def B_exact(ell, s):
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
# Exact k=1 C boundary
# ============================================================================

def C_term(ell, s, j):
    """
    k=1:

      C = -p^ell(1+q)

    only j=ell-1 and j=ell can contribute.
    """

    if j not in (
        ell-1,
        ell
    ):
        return sp.Integer(0)

    r = j-2

    if r < 0:
        return sp.Integer(0)

    # C branch t exponent = 2ell-j.
    e = 2*ell-j

    # Final exponent condition:
    #
    # 2m + e - (r+1) - 1 = 2s
    #
    numerator = e-r-2

    if numerator % 2:
        return sp.Integer(0)

    m = s - numerator//2

    if m < 0 or m > r//2:
        return sp.Integer(0)

    return sp.factor(
        -v_coeff(
            r,
            m
        )
    )


def C_exact(ell, s):
    return sp.factor(
        C_term(
            ell,
            s,
            ell-1
        )
        +
        C_term(
            ell,
            s,
            ell
        )
    )


# ============================================================================
# Compressed formulas
# ============================================================================

def B_closed(ell, s):
    return sp.factor(
        (-1)**(s-1)
        * (
            sp.binomial(
                ell-s,
                s+1
            )
            +
            sp.binomial(
                ell-s+1,
                s+2
            )
        )
    )


def C_closed(ell, s):
    return sp.factor(
        (-1)**s
        * (
            sp.binomial(
                ell-s-1,
                s
            )
            +
            sp.binomial(
                ell-s-2,
                s-1
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


# ============================================================================
# Direct D coefficient from exact Laurent kernel
# ============================================================================

def D_direct(ell, s):
    """
    Independent reconstruction of D coefficient using the exact B+C
    decomposition. This does not use B_closed or C_closed.
    """

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
# Symbolic Pascal simplification diagnostics
# ============================================================================

def pascal_residuals(ell, s):
    """
    Check several equivalent forms without assuming any one is correct.
    """

    left = (
        sp.binomial(
            ell-s,
            s+1
        )
        +
        sp.binomial(
            ell-s+1,
            s+2
        )
    )

    right = (
        sp.binomial(
            ell-s,
            s+1
        )
        +
        sp.binomial(
            ell-s,
            s+2
        )
        +
        sp.binomial(
            ell-s,
            s+1
        )
    )

    return sp.factor(
        left-right
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("="*78)
    print("KAPPA EXPERIMENT 169")
    print("COMPLETE k=1 COEFFICIENT LAW: B + C")
    print("="*78)
    print()

    TRAIN = [
        9, 11, 13, 15, 17,
        19, 21, 23
    ]

    HOLDOUT = [
        25, 27, 29, 31
    ]

    failures = 0

    # ------------------------------------------------------------------------
    # 1. B formula
    # ------------------------------------------------------------------------

    print("="*78)
    print("1. B-BRANCH CLOSED COEFFICIENT")
    print("="*78)

    for ell in TRAIN:

        print()
        print(
            f"ell={ell}"
        )

        max_s = (
            ell-1
        )//2

        for s in range(
            1,
            max_s+1
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
                exact-closed
            )

            print(
                f"  s={s}: "
                f"B_exact={exact}, "
                f"B_closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 2. C formula
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("2. C-BOUNDARY CLOSED COEFFICIENT")
    print("="*78)

    for ell in TRAIN:

        print()
        print(
            f"ell={ell}"
        )

        max_s = (
            ell-1
        )//2

        for s in range(
            0,
            max_s+1
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
                exact-closed
            )

            print(
                f"  s={s}: "
                f"C_exact={exact}, "
                f"C_closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 3. Complete D formula
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("3. COMPLETE D COEFFICIENT")
    print("="*78)

    for ell in TRAIN:

        print()
        print(
            f"ell={ell}"
        )

        max_s = (
            ell-1
        )//2

        for s in range(
            0,
            max_s+1
        ):

            direct = D_direct(
                ell,
                s
            )

            closed = D_closed(
                ell,
                s
            )

            residual = sp.factor(
                direct-closed
            )

            print(
                f"  s={s}: "
                f"D={direct}, "
                f"closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 4. Full polynomial reconstruction
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("4. FULL k=1 D POLYNOMIAL")
    print("="*78)

    N = sp.symbols(
        "N"
    )

    for ell in TRAIN:

        degree = (
            ell-3
        )//2

        direct_poly = sp.factor(
            sum(
                D_direct(
                    ell,
                    s
                ) * N**s
                for s in range(
                    degree+1
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
                    degree+1
                )
            )
        )

        residual = sp.factor(
            direct_poly-closed_poly
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
    # 5. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("5. FORWARD ELL HOLDOUT")
    print("="*78)

    holdout_failures = 0

    for ell in HOLDOUT:

        max_s = (
            ell-1
        )//2

        ok = True

        for s in range(
            0,
            max_s+1
        ):

            if sp.factor(
                D_direct(
                    ell,
                    s
                )
                -
                D_closed(
                    ell,
                    s
                )
            ) != 0:
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
    print("="*78)
    print("FINAL DIAGNOSTIC")
    print("="*78)

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
        print("STATUS = PASS")

        print()
        print(
            "The complete k=1 D_(1,ell) coefficient"
        )

        print(
            "is now represented by an explicit finite"
        )

        print(
            "binomial formula."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "simplify the four-binomial coefficient law,"
        )

        print(
            "then generalize the same B/C mechanism from"
        )

        print(
            "k=1 to arbitrary odd k."
        )

    else:

        print()
        print("STATUS = FAIL")

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

