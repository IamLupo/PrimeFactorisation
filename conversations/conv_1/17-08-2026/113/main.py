#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 176
GENERAL ODD-k B+C COEFFICIENT LAW
BINOMIAL COLLAPSE OF BOTH BRANCHES
==============================================================================

GOAL
----
Move from the separately verified k=1 and k=3 cases to a general odd-k law.

TEST VALUES:
    k in {1,3,5,7}
    ell odd
    ell > k

The proposed exact formulas are

C_(k,ell,s)
    = (-1)^s [
          C(ell-s-1, s-k+1)
        + C(ell-s-2, s-k)
      ]

and

B_(k,ell,s) = 0,                         s < k

B_(k,ell,k)
    = C(ell-1, 2k-1)

B_(k,ell,s)
    = (-1)^(s-k) [
          C(ell+k-s-1, s+2k-1)
        + C(ell+k-s,   s+2k)
      ],
                                             s >= k+1.

The experiment checks these formulas directly against the finite
Laurent/Green-kernel construction.

It also tests the structural consequence

    D_(k,ell,s) = 0       for s < k-1
    D_(k,ell,k-1) = 1

for every tested odd k.

No interpolation.
No numerical fitting.
Exact SymPy arithmetic only.
==============================================================================

"""

import sys
import sympy as sp


# ============================================================================
# Universal Green kernel
# ============================================================================

def V(r, m):
    """
    Universal antisymmetric Green polynomial coefficient.
    """

    if r < 0 or m < 0 or m > r // 2:
        return sp.Integer(0)

    a = sp.binomial(r - m, m)

    b = (
        sp.Integer(0)
        if m == 0
        else sp.binomial(r - m - 1, m - 1)
    )

    return sp.expand(
        (-1) ** (r - m) * (a + b)
    )


# ============================================================================
# General odd-k B branch
# ============================================================================

def B_term(k, ell, s, j):

    m = s - k
    r = j - 2 * k

    if m < 0:
        return sp.Integer(0)

    if r < 0:
        return sp.Integer(0)

    if m > r // 2:
        return sp.Integer(0)

    top = k + j

    if top < 0 or top > ell:
        return sp.Integer(0)

    return sp.factor(
        sp.binomial(
            ell,
            top
        )
        * V(r, m)
    )


def B_exact(k, ell, s):

    return sp.factor(
        sum(
            B_term(
                k,
                ell,
                s,
                j
            )
            for j in range(
                0,
                ell + 1
            )
        )
    )


# ============================================================================
# General odd-k C boundary
# ============================================================================

def C_term(k, ell, s, j):

    if j < ell - k:
        return sp.Integer(0)

    if j > ell:
        return sp.Integer(0)

    a = ell - j

    if a < 0 or a > k:
        return sp.Integer(0)

    r = j - 2 * k

    if r < 0:
        return sp.Integer(0)

    e = 2 * ell - j

    # General t-degree condition:
    #
    # e + 2m - (r+1) - 1 = 2s
    #
    num = e - r - 2

    if num % 2 != 0:
        return sp.Integer(0)

    m = s - num // 2

    if m < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        -sp.binomial(
            k,
            a
        )
        * V(
            r,
            m
        )
    )


def C_exact(k, ell, s):

    return sp.factor(
        sum(
            C_term(
                k,
                ell,
                s,
                j
            )
            for j in range(
                ell - k,
                ell + 1
            )
        )
    )


# ============================================================================
# Proposed general B formula
# ============================================================================

def B_closed(k, ell, s):

    if s < k:
        return sp.Integer(0)

    # Special central edge s=k:
    if s == k:
        return sp.binomial(
            ell - 1,
            2 * k - 1
        )

    return sp.factor(
        (-1) ** (s - k)
        * (
            sp.binomial(
                ell + k - s - 1,
                s + 2 * k - 1
            )
            +
            sp.binomial(
                ell + k - s,
                s + 2 * k
            )
        )
    )


# ============================================================================
# Proposed general C formula
# ============================================================================

def C_closed(k, ell, s):

    return sp.factor(
        (-1) ** s
        * (
            sp.binomial(
                ell - s - 1,
                s - k + 1
            )
            +
            sp.binomial(
                ell - s - 2,
                s - k
            )
        )
    )


# ============================================================================
# Complete coefficient
# ============================================================================

def D_exact(k, ell, s):

    return sp.factor(
        B_exact(
            k,
            ell,
            s
        )
        +
        C_exact(
            k,
            ell,
            s
        )
    )


def D_closed(k, ell, s):

    return sp.factor(
        B_closed(
            k,
            ell,
            s
        )
        +
        C_closed(
            k,
            ell,
            s
        )
    )


# ============================================================================
# Exact finite-kernel support
# ============================================================================

def B_support(k, ell, s):

    return [
        j
        for j in range(
            0,
            ell + 1
        )
        if B_term(
            k,
            ell,
            s,
            j
        ) != 0
    ]


def C_support(k, ell, s):

    return [
        j
        for j in range(
            ell - k,
            ell + 1
        )
        if C_term(
            k,
            ell,
            s,
            j
        ) != 0
    ]


# ============================================================================
# Leading coefficient prediction
# ============================================================================

def predicted_leading(k, ell):

    """
    This is NOT assumed as a theorem.

    For k=3 the observed law was

        (-1)^((ell-9)/2) ((ell-3)/2)^2.

    For general k we therefore test whether the natural pattern is

        (-1)^((ell-(2k+3))/2)
        * binom((ell-k-1)/2, k-1)

    or another explicit binomial expression.

    We print the empirical leading coefficient rather than hard-code
    a theorem for k>3.
    """

    return None


# ============================================================================
# Main experiment
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 176")
    print("GENERAL ODD-k B+C COEFFICIENT LAW")
    print("BINOMIAL COLLAPSE OF BOTH BRANCHES")
    print("=" * 78)

    K_VALUES = [
        1,
        3,
        5,
        7
    ]

    failures_B = 0
    failures_C = 0
    failures_D = 0
    failures_low = 0
    failures_holdout = 0

    # ------------------------------------------------------------------------
    # 1. B formula
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. GENERAL B-FORMULA")
    print("=" * 78)

    for k in K_VALUES:

        print()
        print(f"k={k}")

        TRAIN_ELL = list(
            range(
                k + 2 if (k + 2) % 2 == 1 else k + 3,
                32,
                2
            )
        )

        for ell in TRAIN_ELL:

            max_s = (
                ell - 1
            ) // 2

            print()
            print(f"  ell={ell}")

            for s in range(
                max_s + 1
            ):

                exact = B_exact(
                    k,
                    ell,
                    s
                )

                closed = B_closed(
                    k,
                    ell,
                    s
                )

                residual = sp.factor(
                    exact - closed
                )

                print(
                    f"    s={s}: "
                    f"exact={exact}, "
                    f"closed={closed}, "
                    f"residual={residual}"
                )

                if residual != 0:
                    failures_B += 1

    # ------------------------------------------------------------------------
    # 2. C formula
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. GENERAL C-FORMULA")
    print("=" * 78)

    for k in K_VALUES:

        print()
        print(f"k={k}")

        TRAIN_ELL = list(
            range(
                k + 2 if (k + 2) % 2 == 1 else k + 3,
                32,
                2
            )
        )

        for ell in TRAIN_ELL:

            max_s = (
                ell - 1
            ) // 2

            print()
            print(f"  ell={ell}")

            for s in range(
                max_s + 1
            ):

                exact = C_exact(
                    k,
                    ell,
                    s
                )

                closed = C_closed(
                    k,
                    ell,
                    s
                )

                residual = sp.factor(
                    exact - closed
                )

                print(
                    f"    s={s}: "
                    f"support={C_support(k,ell,s)}, "
                    f"exact={exact}, "
                    f"closed={closed}, "
                    f"residual={residual}"
                )

                if residual != 0:
                    failures_C += 1

    # ------------------------------------------------------------------------
    # 3. Complete D formula
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. GENERAL COMPLETE D-FORMULA")
    print("=" * 78)

    for k in K_VALUES:

        print()
        print(f"k={k}")

        TRAIN_ELL = list(
            range(
                k + 2 if (k + 2) % 2 == 1 else k + 3,
                32,
                2
            )
        )

        for ell in TRAIN_ELL:

            max_s = (
                ell - 1
            ) // 2

            print()
            print(
                f"  ell={ell}"
            )

            for s in range(
                max_s + 1
            ):

                exact = D_exact(
                    k,
                    ell,
                    s
                )

                closed = D_closed(
                    k,
                    ell,
                    s
                )

                residual = sp.factor(
                    exact - closed
                )

                if residual != 0:
                    failures_D += 1

                print(
                    f"    s={s}: "
                    f"D={exact}, "
                    f"closed={closed}, "
                    f"residual={residual}"
                )

    # ------------------------------------------------------------------------
    # 4. Universal low-order structure
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. UNIVERSAL LOW-ORDER LAW")
    print("=" * 78)

    for k in K_VALUES:

        print()
        print(f"k={k}")

        ell = k + 6
        if ell % 2 == 0:
            ell += 1

        for s in range(
            0,
            k + 1
        ):

            actual = D_closed(
                k,
                ell,
                s
            )

            if s < k - 1:
                expected = sp.Integer(0)

            elif s == k - 1:
                expected = sp.Integer(1)

            else:
                expected = None

            if expected is not None:

                residual = sp.factor(
                    actual - expected
                )

                print(
                    f"  s={s}: "
                    f"actual={actual}, "
                    f"expected={expected}, "
                    f"residual={residual}"
                )

                if residual != 0:
                    failures_low += 1

    # ------------------------------------------------------------------------
    # 5. Compare known k=1 and k=3 formulas
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. k=1 / k=3 CROSS-CHECK")
    print("=" * 78)

    for k in [1, 3]:

        ell = 17

        max_s = (
            ell - 1
        ) // 2

        print()
        print(
            f"k={k}, ell={ell}"
        )

        for s in range(
            max_s + 1
        ):

            print(
                f"  s={s}: "
                f"B={B_closed(k,ell,s)}, "
                f"C={C_closed(k,ell,s)}, "
                f"D={D_closed(k,ell,s)}"
            )

    # ------------------------------------------------------------------------
    # 6. Fresh forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FRESH ODD-k HOLDOUT")
    print("=" * 78)

    HOLDOUTS = [
        (1, 33),
        (1, 35),
        (3, 33),
        (3, 35),
        (5, 33),
        (5, 35),
        (7, 33),
        (7, 35),
    ]

    for k, ell in HOLDOUTS:

        if ell <= k:
            continue

        max_s = (
            ell - 1
        ) // 2

        ok = True

        for s in range(
            max_s + 1
        ):

            rB = sp.factor(
                B_exact(
                    k,
                    ell,
                    s
                )
                -
                B_closed(
                    k,
                    ell,
                    s
                )
            )

            rC = sp.factor(
                C_exact(
                    k,
                    ell,
                    s
                )
                -
                C_closed(
                    k,
                    ell,
                    s
                )
            )

            rD = sp.factor(
                D_exact(
                    k,
                    ell,
                    s
                )
                -
                D_closed(
                    k,
                    ell,
                    s
                )
            )

            if rB != 0 or rC != 0 or rD != 0:

                ok = False
                failures_holdout += 1

                print(
                    f"  FAIL: k={k}, ell={ell}, s={s}"
                )

        print(
            f"  (k={k}, ell={ell}) "
            f"status={'PASS' if ok else 'FAIL'}"
        )

    # ------------------------------------------------------------------------
    # 7. Polynomial leading coefficients
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. LEADING COEFFICIENT DATA")
    print("=" * 78)

    N = sp.symbols(
        "N"
    )

    for k in K_VALUES:

        print()
        print(
            f"k={k}"
        )

        for ell in [15, 17, 19, 21, 23, 25]:

            if ell <= k:
                continue

            degree = (
                ell - 1
            ) // 2

            poly = sp.expand(
                sum(
                    D_closed(
                        k,
                        ell,
                        s
                    ) * N**s
                    for s in range(
                        degree + 1
                    )
                )
            )

            quotient = sp.cancel(
                poly / N**(k - 1)
            )

            P = sp.Poly(
                quotient,
                N
            )

            leading = sp.LC(
                P
            )

            print(
                f"  ell={ell}: "
                f"leading={leading}"
            )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "B formula failures =",
        failures_B
    )

    print(
        "C formula failures =",
        failures_C
    )

    print(
        "D formula failures =",
        failures_D
    )

    print(
        "low-order failures =",
        failures_low
    )

    print(
        "holdout failures =",
        failures_holdout
    )

    if (
        failures_B == 0
        and failures_C == 0
        and failures_D == 0
        and failures_low == 0
        and failures_holdout == 0
    ):

        print()
        print("STATUS = PASS")

        print()
        print(
            "The B and C coefficient laws now generalize"
        )

        print(
            "from k=1 and k=3 to all tested odd k."
        )

        print()
        print(
            "The universal structural consequence is:"
        )

        print(
            "    D_(k,ell,s)=0 for s<k-1"
        )

        print(
            "    D_(k,ell,k-1)=1"
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "derive a single closed expression for"
        )

        print(
            "D_(k,ell,s)=B_(k,ell,s)+C_(k,ell,s),"
        )

        print(
            "then analyze the resulting polynomial"
        )

        print(
            "D_(k,ell)(N)=N^(k-1) P_(k,ell)(N)"
        )

        print(
            "for a general odd k."
        )

    else:

        print()
        print("STATUS = PARTIAL")

        print(
            "At least one proposed generalization"
        )

        print(
            "failed exact verification."
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

