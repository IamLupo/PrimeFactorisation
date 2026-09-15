#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 178
EXACT ODD-k EDGE THEOREM
B_(k,ell,k) = C(ell-1,3k-1)
C_(k,ell,k-1) = 1
C_(k,ell,k)   = -(ell-k)

AND COMPLETE B+C RECONSTRUCTION
==============================================================================

OBJECTIVE
---------
Experiment 177 isolated the problem:

  * generic interior B formula: exact
  * proposed B edge hypergeometric specialization: wrong
  * C boundary generic law: exact except for truncated-support cases

The numerical data reveal a much simpler exact edge structure.

CONJECTURED ACTIVE-RANGE LAWS
-----------------------------

For odd k >= 1 and ell >= 3k:

    B(k,ell,k) = binomial(ell-1, 3k-1)

For odd k >= 3 and ell >= 2k+1:

    C(k,ell,k-1) = 1

For odd k >= 1 and ell >= 2k+1:

    C(k,ell,k) = -(ell-k)

Hence

    D(k,ell,k)
      = binomial(ell-1, 3k-1) - (ell-k)

for ell >= 3k.

At s=k-1:

    D(k,ell,k-1) = 1

when that layer is active.

For k=1 the first layer is s=0, so separately:

    D(1,ell,0) = 1.

This experiment performs an exact symbolic verification and also
tests whether the B edge identity can be obtained DIRECTLY from the
finite V-sum rather than assumed from numerical interpolation.

IMPORTANT
---------
No floating point.
No guessed hypergeometric expression at the edge.
Boundary support is respected explicitly.
==============================================================================

"""

import sys
import sympy as sp


# ============================================================================
# V KERNEL
# ============================================================================

def V(r, m):

    if r < 0:
        return sp.Integer(0)

    if m < 0 or 2*m > r:
        return sp.Integer(0)

    return sp.expand(
        (-1)**(r-m) *
        (
            sp.binomial(r-m, m)
            +
            (
                sp.binomial(r-m-1, m-1)
                if m >= 1
                else 0
            )
        )
    )


# ============================================================================
# EXACT B SOURCE
# ============================================================================

def B_term(k, ell, s, j):

    m = s-k
    r = j-2*k

    if m < 0 or r < 0 or 2*m > r:
        return sp.Integer(0)

    top = j+k

    if top < 0 or top > ell:
        return sp.Integer(0)

    return sp.expand(
        sp.binomial(ell, top) * V(r,m)
    )


def B_exact(k, ell, s):

    return sp.factor(
        sum(
            B_term(k,ell,s,j)
            for j in range(ell+1)
        )
    )


# ============================================================================
# EXACT C SOURCE
# ============================================================================

def C_term(k, ell, s, j):

    if not (ell-k <= j <= ell):
        return sp.Integer(0)

    a = ell-j

    if a < 0 or a > k:
        return sp.Integer(0)

    r = j-2*k

    if r < 0:
        return sp.Integer(0)

    # exponent condition
    e = 2*ell-j

    # e + 2m -(r+1)-1 = 2s
    m_num = 2*s - e + r + 2

    if m_num % 2:
        return sp.Integer(0)

    m = m_num//2

    if m < 0 or 2*m > r:
        return sp.Integer(0)

    return sp.expand(
        -sp.binomial(k,a)*V(r,m)
    )


def C_exact(k, ell, s):

    return sp.factor(
        sum(
            C_term(k,ell,s,j)
            for j in range(ell-k,ell+1)
        )
    )


# ============================================================================
# INTERIOR FORMULAS ALREADY VERIFIED
# ============================================================================

def B_interior(k, ell, s):

    if s < k:
        return sp.Integer(0)

    if s == k:
        raise ValueError("edge layer")

    return sp.factor(
        (-1)**(s-k) *
        (
            sp.binomial(
                ell+k-s-1,
                s+2*k-1
            )
            +
            sp.binomial(
                ell+k-s,
                s+2*k
            )
        )
    )


def C_generic(k, ell, s):

    return sp.factor(
        (-1)**s *
        (
            sp.binomial(
                ell-s-1,
                s-k+1
            )
            +
            sp.binomial(
                ell-s-2,
                s-k
            )
        )
    )


# ============================================================================
# EXACT B EDGE CONJECTURE
# ============================================================================

def B_edge_closed(k, ell):

    return sp.binomial(
        ell-1,
        3*k-1
    )


# ============================================================================
# EXACT C EDGE CONJECTURES
# ============================================================================

def C_k_minus_1_closed(k, ell):

    return sp.Integer(1)


def C_k_closed(k, ell):

    return -(ell-k)


# ============================================================================
# COMPLETE EDGE D
# ============================================================================

def D_edge_closed(k, ell):

    return sp.factor(
        B_edge_closed(k,ell)
        +
        C_k_closed(k,ell)
    )


# ============================================================================
# DIRECT SYMBOLIC B EDGE EXPANSION
# ============================================================================

def B_edge_direct_terms(k, ell):

    """
    At s=k, m=0, so

        V(r,0) = (-1)^r.

    Therefore

        B(k,ell,k)
        = sum_j (-1)^(j-2k) C(ell,j+k).

    Reindex q=j+k:

        q runs over [3k, ell+k]

    but the binomial support truncates to q <= ell.

    Hence the actual active range is

        q = 3k,...,ell.

    The experiment prints the direct alternating sum and compares it
    against C(ell-1,3k-1).
    """

    if ell < 3*k:
        return []

    return [
        (
            q,
            sp.expand(
                (-1)**(q-k) *
                sp.binomial(ell,q)
            )
        )
        for q in range(3*k,ell+1)
    ]


# ============================================================================
# SYMBOLIC PASCAL COLLAPSE CHECK
# ============================================================================

def alternating_binomial_sum(k, ell):

    terms = B_edge_direct_terms(k,ell)

    return sp.factor(
        sum(
            term
            for _,term in terms
        )
    )


# ============================================================================
# TEST B EDGE
# ============================================================================

def test_B_edge():

    failures = 0

    print()
    print("="*78)
    print("1. B EDGE THEOREM")
    print("="*78)

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+1,3),
            52,
            2
        ):

            exact = B_exact(
                k,ell,k
            )

            closed = B_edge_closed(
                k,ell
            )

            # The binomial is automatically zero in the inactive range.
            residual = sp.factor(
                exact-closed
            )

            print(
                f"  ell={ell}: "
                f"exact={exact}, "
                f"closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

            # Direct alternating representation only where active.
            if ell >= 3*k:

                direct = alternating_binomial_sum(
                    k,ell
                )

                direct_residual = sp.factor(
                    direct-closed
                )

                if direct_residual != 0:

                    failures += 1

                    print(
                        "    DIRECT FAIL:",
                        f"direct={direct},",
                        f"closed={closed},",
                        f"residual={direct_residual}"
                    )

    print(
        "B edge failures =",
        failures
    )

    return failures


# ============================================================================
# TEST C EDGE
# ============================================================================

def test_C_edges():

    failures = 0

    print()
    print("="*78)
    print("2. C EDGE THEOREMS")
    print("="*78)

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+1,3),
            52,
            2
        ):

            # k-1 layer
            s1 = k-1

            exact1 = C_exact(
                k,ell,s1
            )

            closed1 = C_k_minus_1_closed(
                k,ell
            ) if ell >= 2*k+1 else sp.Integer(0)

            residual1 = sp.factor(
                exact1-closed1
            )

            # k layer
            s2 = k

            exact2 = C_exact(
                k,ell,s2
            )

            closed2 = (
                C_k_closed(k,ell)
                if ell >= 2*k+1
                else sp.Integer(0)
            )

            residual2 = sp.factor(
                exact2-closed2
            )

            print(
                f"  ell={ell}: "
                f"C[k-1]={exact1} vs {closed1}, "
                f"res={residual1}; "
                f"C[k]={exact2} vs {closed2}, "
                f"res={residual2}"
            )

            if residual1 != 0:
                failures += 1

            if residual2 != 0:
                failures += 1

    print(
        "C edge failures =",
        failures
    )

    return failures


# ============================================================================
# TEST INTERIOR B/C
# ============================================================================

def test_interior():

    failures = 0

    print()
    print("="*78)
    print("3. INTERIOR FORMULAS")
    print("="*78)

    for k in [1,3,5,7,9]:

        for ell in range(
            k+3,
            42,
            2
        ):

            max_s = (ell-1)//2

            for s in range(
                k+1,
                max_s+1
            ):

                b = B_exact(
                    k,ell,s
                )

                b_closed = B_interior(
                    k,ell,s
                )

                c = C_exact(
                    k,ell,s
                )

                c_closed = C_generic(
                    k,ell,s
                )

                rb = sp.factor(
                    b-b_closed
                )

                rc = sp.factor(
                    c-c_closed
                )

                if rb != 0:

                    failures += 1

                    print(
                        "B FAIL:",
                        k,ell,s,
                        rb
                    )

                if rc != 0:

                    failures += 1

                    print(
                        "C FAIL:",
                        k,ell,s,
                        rc
                    )

    print(
        "interior failures =",
        failures
    )

    return failures


# ============================================================================
# COMPLETE D LAW
# ============================================================================

def D_closed(k,ell,s):

    if s < k-1:
        return sp.Integer(0)

    if k == 1 and s == 0:
        return sp.Integer(1)

    if s == k-1:

        # active boundary
        if ell < 2*k+1:
            return sp.Integer(0)

        return sp.Integer(1)

    if s == k:

        if ell < 3*k:
            # Must use exact edge support.
            return sp.factor(
                B_exact(k,ell,s)
                +
                C_exact(k,ell,s)
            )

        return D_edge_closed(
            k,ell
        )

    return sp.factor(
        B_interior(k,ell,s)
        +
        C_generic(k,ell,s)
    )


def test_D():

    failures = 0

    print()
    print("="*78)
    print("4. COMPLETE D RECONSTRUCTION")
    print("="*78)

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+1,3),
            46,
            2
        ):

            max_s = (ell-1)//2

            for s in range(
                max_s+1
            ):

                exact = sp.factor(
                    B_exact(
                        k,ell,s
                    )
                    +
                    C_exact(
                        k,ell,s
                    )
                )

                closed = sp.factor(
                    D_closed(
                        k,ell,s
                    )
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

    print(
        "D failures =",
        failures
    )

    return failures


# ============================================================================
# EDGE D TABLE
# ============================================================================

def print_edge_D():

    print()
    print("="*78)
    print("5. EDGE D FORMULA")
    print("="*78)

    for k in [1,3,5,7,9]:

        print()
        print(f"k={k}")

        for ell in range(
            3*k,
            3*k+18,
            2
        ):

            value = D_edge_closed(
                k,ell
            )

            exact = sp.factor(
                B_exact(k,ell,k)
                +
                C_exact(k,ell,k)
            )

            print(
                f"  ell={ell}: "
                f"D_k={exact}, "
                f"closed={value}, "
                f"residual={sp.factor(exact-value)}"
            )


# ============================================================================
# POLYNOMIAL / LEADING DATA
# ============================================================================

def test_leading():

    print()
    print("="*78)
    print("6. EDGE LEADING DATA")
    print("="*78)

    for k in [1,3,5,7,9]:

        print()
        print(f"k={k}")

        for ell in range(
            max(3*k,15),
            42,
            2
        ):

            # D coefficient at s=k is the boundary coefficient.
            value = D_edge_closed(
                k,ell
            )

            print(
                f"  ell={ell}: "
                f"D_k={value}"
            )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print("KAPPA EXPERIMENT 178")
    print("EXACT ODD-k EDGE THEOREM")
    print("="*78)

    b_fail = test_B_edge()
    c_fail = test_C_edges()
    i_fail = test_interior()
    d_fail = test_D()

    print_edge_D()
    test_leading()

    print()
    print("="*78)
    print("FINAL DIAGNOSTIC")
    print("="*78)

    print(
        "B edge failures      =",
        b_fail
    )

    print(
        "C edge failures      =",
        c_fail
    )

    print(
        "interior failures    =",
        i_fail
    )

    print(
        "D reconstruction     =",
        d_fail
    )

    if (
        b_fail == 0
        and c_fail == 0
        and i_fail == 0
        and d_fail == 0
    ):

        print()
        print("STATUS = PASS")

        print()
        print(
            "The odd-k boundary structure is exact."
        )

        print()
        print(
            "In particular, the active edge satisfies:"
        )
        print(
            "  B(k,ell,k) = C(ell-1,3k-1)"
        )
        print(
            "  C(k,ell,k-1) = 1"
        )
        print(
            "  C(k,ell,k) = -(ell-k)"
        )

        print()
        print(
            "Therefore:"
        )
        print(
            "  D(k,ell,k)"
            " = C(ell-1,3k-1) - (ell-k)"
        )

        print()
        print(
            "NEXT TARGET:"
        )
        print(
            "derive a uniform symbolic formula for"
        )
        print(
            "all s >= k-1 and then compare the"
        )
        print(
            "result across k=1,3,5,7,9,11."
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

