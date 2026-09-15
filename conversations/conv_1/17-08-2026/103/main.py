#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 166
CORRECT PCHAMMER HYPERGEOMETRIC B-BRANCH + C BOUNDARY CERTIFICATE

TARGET
------
Correct the sign normalization exposed by Experiment 165R.

For k=1:

s=1:
    B_n / B_0 = (-L)_n / (4)_n,
    L = ell-3.

s=2:
    B_n / B_0 =
        (-L)_n (3)_n / ((2)_n (6)_n),
    L = ell-5.

The experiment verifies these identities term-by-term and then isolates
the short C boundary contribution.

It also derives the corresponding finite hypergeometric sums.

NO POLYNOMIAL FITTING
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp

n = sp.symbols("n", integer=True, nonnegative=True)


# ============================================================================
# UNIVERSAL V COEFFICIENT
# ============================================================================

def v_coeff(r, m):
    if r < 0 or m < 0 or m > r // 2:
        return sp.Integer(0)

    a = sp.binomial(r - m, m)

    b = (
        sp.Integer(0)
        if m == 0
        else sp.binomial(r - m - 1, m - 1)
    )

    return sp.expand(
        (-1)**(r-m) * (a+b)
    )


# ============================================================================
# BRANCH B
# ============================================================================

def B_term(k, ell, s, j):
    """
    Exact branch-B contribution to [N^s]D.
    """

    lo = -k
    hi = ell-k

    if not (lo <= j <= hi):
        return sp.Integer(0)

    r = j - 2*k

    if r < 0:
        return sp.Integer(0)

    # B has final N degree k+m.
    m = s-k

    if m < 0 or m > r//2:
        return sp.Integer(0)

    return sp.factor(
        sp.binomial(ell, k+j)
        * v_coeff(r, m)
    )


def B_support(k, ell, s):
    return [
        j
        for j in range(
            -k,
            ell-k+1
        )
        if B_term(
            k,
            ell,
            s,
            j
        ) != 0
    ]


def B_total(k, ell, s):
    return sp.factor(
        sum(
            B_term(k, ell, s, j)
            for j in B_support(
                k,
                ell,
                s
            )
        )
    )


# ============================================================================
# BRANCH C
# ============================================================================

def C_term(k, ell, s, j):
    lo = ell-k
    hi = ell

    if not (lo <= j <= hi):
        return sp.Integer(0)

    r = j-2*k

    if r < 0:
        return sp.Integer(0)

    e = 2*ell-j

    residual = e-r-2*k

    if residual % 2:
        return sp.Integer(0)

    m = s-residual//2

    if m < 0 or m > r//2:
        return sp.Integer(0)

    return sp.factor(
        -sp.binomial(k, ell-j)
        * v_coeff(r,m)
    )


def C_support(k, ell, s):
    return [
        j
        for j in range(
            ell-k,
            ell+1
        )
        if C_term(
            k,
            ell,
            s,
            j
        ) != 0
    ]


def C_total(k, ell, s):
    return sp.factor(
        sum(
            C_term(k, ell, s, j)
            for j in C_support(
                k,
                ell,
                s
            )
        )
    )


# ============================================================================
# K=1, S=1 PCHAMMER FORM
# ============================================================================

def B_k1_s1_pochhammer(ell):
    """
    j = 2+n.

    Exact ratio:

        B_(n+1)/B_n = -(ell-3-n)/(n+1)

    and therefore

        B_n = B_0 * (-ell+3)_n / (4)_n.
    """

    k = 1
    s = 1

    support = B_support(
        k,
        ell,
        s
    )

    if not support:
        return None

    j0 = support[0]
    B0 = B_term(
        k,
        ell,
        s,
        j0
    )

    L = ell-3

    candidate = sp.factor(
        B0
        * sp.rf(-L,n)
        / sp.rf(1+3,n)
    )

    return (
        j0,
        support,
        candidate
    )


# ============================================================================
# K=1, S=2 PCHAMMER FORM
# ============================================================================

def B_k1_s2_pochhammer(ell):
    """
    j = 4+n.

    Exact ratio:

       -(ell-5-n)(n+3) / ((n+2)(n+6))

    so

       B_n =
       B_0 *
       (-ell+5)_n (3)_n
       --------------------------------
       (2)_n (6)_n.
    """

    k = 1
    s = 2

    support = B_support(
        k,
        ell,
        s
    )

    if not support:
        return None

    j0 = support[0]

    B0 = B_term(
        k,
        ell,
        s,
        j0
    )

    L = ell-5

    candidate = sp.factor(
        B0
        * sp.rf(-L,n)
        * sp.rf(3,n)
        /
        (
            sp.rf(2,n)
            * sp.rf(6,n)
        )
    )

    return (
        j0,
        support,
        candidate
    )


# ============================================================================
# VERIFY C BOUNDARY
# ============================================================================

def verify_boundary(k, ell, s):
    values = []

    for j in C_support(
        k,
        ell,
        s
    ):
        values.append(
            (
                j,
                C_term(
                    k,
                    ell,
                    s,
                    j
                )
            )
        )

    return values


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print("KAPPA EXPERIMENT 166")
    print("CORRECT PCHAMMER HYPERGEOMETRIC B-BRANCH")
    print("="*78)
    print()

    failures = 0

    TEST_ELL = [9,11,13,15,17]

    # ------------------------------------------------------------------------
    # 1. s=1
    # ------------------------------------------------------------------------

    print("="*78)
    print("1. k=1, s=1 PCHAMMER CERTIFICATE")
    print("="*78)

    for ell in TEST_ELL:

        result = B_k1_s1_pochhammer(
            ell
        )

        if result is None:
            continue

        j0, support, candidate = result

        print()
        print(
            f"(1,{ell})"
        )

        print(
            "  support =",
            support
        )

        print(
            "  candidate =",
            candidate
        )

        for jj in support:

            nn = jj-j0

            predicted = sp.factor(
                candidate.subs(
                    n,
                    nn
                )
            )

            actual = B_term(
                1,
                ell,
                1,
                jj
            )

            residual = sp.factor(
                predicted-actual
            )

            print(
                f"    j={jj}, n={nn}, "
                f"actual={actual}, "
                f"predicted={predicted}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

        # Sum check
        candidate_sum = sp.factor(
            sum(
                candidate.subs(n,nn)
                for nn in range(
                    len(support)
                )
            )
        )

        exact_sum = B_total(
            1,
            ell,
            1
        )

        print(
            "  candidate sum =",
            candidate_sum
        )

        print(
            "  exact B sum =",
            exact_sum
        )

        if sp.factor(
            candidate_sum-exact_sum
        ) != 0:
            failures += 1

    # ------------------------------------------------------------------------
    # 2. s=2
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("2. k=1, s=2 PCHAMMER CERTIFICATE")
    print("="*78)

    for ell in TEST_ELL:

        result = B_k1_s2_pochhammer(
            ell
        )

        if result is None:
            continue

        j0, support, candidate = result

        print()
        print(
            f"(1,{ell})"
        )

        print(
            "  support =",
            support
        )

        print(
            "  candidate =",
            candidate
        )

        for jj in support:

            nn = jj-j0

            predicted = sp.factor(
                candidate.subs(
                    n,
                    nn
                )
            )

            actual = B_term(
                1,
                ell,
                2,
                jj
            )

            residual = sp.factor(
                predicted-actual
            )

            print(
                f"    j={jj}, n={nn}, "
                f"actual={actual}, "
                f"predicted={predicted}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

        candidate_sum = sp.factor(
            sum(
                candidate.subs(n,nn)
                for nn in range(
                    len(support)
                )
            )
        )

        exact_sum = B_total(
            1,
            ell,
            2
        )

        print(
            "  candidate sum =",
            candidate_sum
        )

        print(
            "  exact B sum =",
            exact_sum
        )

        if sp.factor(
            candidate_sum-exact_sum
        ) != 0:
            failures += 1

    # ------------------------------------------------------------------------
    # 3. Boundary C
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("3. C BOUNDARY CONTRIBUTIONS")
    print("="*78)

    for ell in TEST_ELL:

        print()
        print(
            f"(1,{ell})"
        )

        for s in [0,1,2]:

            values = verify_boundary(
                1,
                ell,
                s
            )

            total = C_total(
                1,
                ell,
                s
            )

            print(
                f"  s={s}: "
                f"values={values}, "
                f"total={total}"
            )

    # ------------------------------------------------------------------------
    # 4. Simple B+C reconstruction for k=1
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("4. k=1 B+C RECONSTRUCTION")
    print("="*78)

    for ell in TEST_ELL:

        print()
        print(
            f"(1,{ell})"
        )

        for s in range(
            0,
            ell-1
        ):

            total = sp.factor(
                B_total(
                    1,
                    ell,
                    s
                )
                +
                C_total(
                    1,
                    ell,
                    s
                )
            )

            print(
                f"  [N^{s}] D = {total}"
            )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("FINAL DIAGNOSTIC")
    print("="*78)

    print(
        "Pochhammer verification failures =",
        failures
    )

    if failures == 0:

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The corrected Pochhammer normalization matches"
        )

        print(
            "the exact B branch term-by-term."
        )

        print()
        print(
            "The next target is to evaluate the finite"
        )

        print(
            "hypergeometric sums using Vandermonde /"
        )

        print(
            "Chu-Vandermonde identities."
        )

    else:

        print()
        print(
            "STATUS = FAIL"
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

