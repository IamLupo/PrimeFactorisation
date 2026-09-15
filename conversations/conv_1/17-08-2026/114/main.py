#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 177
EXCEPTIONAL B/C BOUNDARY LAYERS
GENERAL ODD-k EDGE FORMULAS
==============================================================================

PURPOSE
-------
Experiment 176 shows that the generic B/C formulas are essentially correct
away from the first active coefficient layer.

The failures are concentrated at:

    B: s = k
    C: smallest admissible s near k-1 / k

This experiment derives those edge layers directly from the finite kernel
and searches for exact binomial closed forms.

We do NOT alter the verified interior formulas.

For B, use the exact Chu-Vandermonde representation for s=k:

    B_(k,ell,k)
      = C(ell, 2k+3)
        * (k+6)_L / (2k+4)_L,

    L = ell - 2k - 3,

with the convention that the second branch is absent.

For C, explicitly enumerate the active boundary terms and compare them
against candidate binomial forms.

The main goals are:

  1. identify the exact B_(k,ell,k) simplification;
  2. identify the exact C_(k,ell,k-1) and C_(k,ell,k) formulas;
  3. verify that the generic formulas remain valid for s >= k+1;
  4. test fresh k=1,3,5,7,9 holdouts;
  5. reconstruct the complete D coefficient law.

Exact SymPy arithmetic only.
==============================================================================

"""

import sys
import sympy as sp

from math import comb


# ============================================================================
# UNIVERSAL GREEN POLYNOMIAL
# ============================================================================

def V(r, m):
    if r < 0 or m < 0 or 2*m > r:
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

    m = s - k
    r = j - 2*k

    if m < 0 or r < 0 or 2*m > r:
        return sp.Integer(0)

    top = j + k

    if top < 0 or top > ell:
        return sp.Integer(0)

    return sp.expand(
        sp.binomial(ell, top) * V(r, m)
    )


def B_exact(k, ell, s):

    return sp.factor(
        sum(
            B_term(k, ell, s, j)
            for j in range(ell + 1)
        )
    )


# ============================================================================
# EXACT C SOURCE
# ============================================================================

def C_term(k, ell, s, j):

    if j < ell-k or j > ell:
        return sp.Integer(0)

    a = ell-j

    if a < 0 or a > k:
        return sp.Integer(0)

    r = j - 2*k

    if r < 0:
        return sp.Integer(0)

    e = 2*ell-j

    # t-power condition
    #
    # e + 2m - (r+1) - 1 = 2s
    #
    m_num = 2*s - e + r + 2

    if m_num % 2 != 0:
        return sp.Integer(0)

    m = m_num // 2

    if m < 0 or 2*m > r:
        return sp.Integer(0)

    return sp.expand(
        -sp.binomial(k, a) * V(r, m)
    )


def C_exact(k, ell, s):

    return sp.factor(
        sum(
            C_term(k, ell, s, j)
            for j in range(ell-k, ell+1)
        )
    )


# ============================================================================
# GENERIC INTERIOR FORMULAS
# ============================================================================

def B_generic(k, ell, s):

    if s < k:
        return sp.Integer(0)

    # Deliberately exclude s=k.
    if s == k:
        raise ValueError("B_generic is for s >= k+1")

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
# B EDGE: DIRECT CHU-VANDERMONDE FORM
# ============================================================================

def B_edge_hyper(k, ell):

    L = ell - 2*k - 3

    if L < 0:
        return sp.Integer(0)

    # First terminating 2F1 only.
    #
    # 2F1(-L, k-2 ; 2k+4 ; 1)
    # = (k+6)_L / (2k+4)_L

    return sp.factor(
        sp.binomial(
            ell,
            2*k+3
        )
        *
        sp.rf(k+6, L)
        /
        sp.rf(2*k+4, L)
    )


# ============================================================================
# B EDGE: SEARCH FOR SIMPLE BINOMIAL EQUIVALENT
# ============================================================================

def find_binomial_equivalent(k, ell):

    exact = B_exact(
        k,
        ell,
        k
    )

    candidates = []

    # Search small affine binomial indices.
    #
    # C(ell-a, b*ell+c) is too broad, so focus on
    # affine-in-ell indices with coefficient 1/2.
    for a in range(0, 10):
        for c_num in range(-15, 16):

            num = ell + c_num

            if num % 2 != 0:
                continue

            r = num // 2

            if r < 0 or r > ell-a:
                continue

            cand = sp.binomial(
                ell-a,
                r
            )

            if sp.simplify(
                cand - exact
            ) == 0:

                candidates.append(
                    ("half-index", a, c_num, cand)
                )

    # Also search fixed small lower indices.
    for a in range(0, 10):
        for r in range(0, 20):

            if r > ell-a:
                continue

            cand = sp.binomial(
                ell-a,
                r
            )

            if sp.simplify(
                cand - exact
            ) == 0:

                candidates.append(
                    ("fixed-index", a, r, cand)
                )

    return candidates


# ============================================================================
# C EDGE SUPPORT
# ============================================================================

def C_active(k, ell, s):

    return [
        (
            j,
            sp.factor(
                C_term(
                    k,
                    ell,
                    s,
                    j
                )
            )
        )
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


# ============================================================================
# C EDGE CANDIDATE SEARCH
# ============================================================================

def search_C_candidate(k, ell, s):

    exact = C_exact(
        k,
        ell,
        s
    )

    candidates = []

    # Search:
    #   +/- C(ell-a, r)
    #   +/- C(ell-a,r1) +/- C(ell-a,r2)

    one_terms = []

    for sign in [1, -1]:

        for a in range(0, 10):

            for r in range(0, 20):

                if r > ell-a:
                    continue

                cand = sign * sp.binomial(
                    ell-a,
                    r
                )

                if sp.simplify(
                    cand - exact
                ) == 0:

                    one_terms.append(
                        (sign, a, r, cand)
                    )

    candidates.extend(
        ("one-term", x)
        for x in one_terms
    )

    for a in range(0, 8):

        for r1 in range(0, 15):

            for r2 in range(0, 15):

                if r1 > ell-a or r2 > ell-a:
                    continue

                cand = (
                    sp.binomial(
                        ell-a,
                        r1
                    )
                    +
                    sp.binomial(
                        ell-a,
                        r2
                    )
                )

                for sign in [1, -1]:

                    test = sign*cand

                    if sp.simplify(
                        test - exact
                    ) == 0:

                        candidates.append(
                            (
                                "two-term",
                                sign,
                                a,
                                r1,
                                r2,
                                test
                            )
                        )

    return candidates


# ============================================================================
# DERIVE A RATIONAL FORM FROM CHU EDGE
# ============================================================================

def B_edge_simplified(k, ell):

    expr = B_edge_hyper(
        k,
        ell
    )

    return sp.factor(
        sp.cancel(expr)
    )


# ============================================================================
# COMPLETE D
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


def B_fixed(k, ell, s):

    if s < k:
        return sp.Integer(0)

    if s == k:
        return B_edge_hyper(
            k,
            ell
        )

    return B_generic(
        k,
        ell,
        s
    )


def C_fixed(k, ell, s):

    return C_generic(
        k,
        ell,
        s
    )


def D_fixed(k, ell, s):

    return sp.factor(
        B_fixed(
            k,
            ell,
            s
        )
        +
        C_fixed(
            k,
            ell,
            s
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print("KAPPA EXPERIMENT 177")
    print("EXCEPTIONAL B/C BOUNDARY LAYERS")
    print("GENERAL ODD-k EDGE FORMULAS")
    print("="*78)

    K_VALUES = [
        1, 3, 5, 7, 9
    ]

    failures_B_edge = 0
    failures_B_interior = 0
    failures_C_edge = 0
    failures_D = 0

    # ------------------------------------------------------------------------
    # 1. B EDGE
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("1. B EDGE s=k")
    print("="*78)

    for k in K_VALUES:

        print()
        print(f"k={k}")

        for ell in range(
            k+2,
            42,
            2
        ):

            exact = B_exact(
                k,
                ell,
                k
            )

            hyper = B_edge_hyper(
                k,
                ell
            )

            residual = sp.factor(
                exact - hyper
            )

            print(
                f"  ell={ell}: "
                f"exact={exact}, "
                f"hyper={hyper}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures_B_edge += 1

            simple = find_binomial_equivalent(
                k,
                ell
            )

            if simple:
                print(
                    "    simple binomial candidate:",
                    simple[:3]
                )

    # ------------------------------------------------------------------------
    # 2. INTERIOR B
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("2. INTERIOR B s>=k+1")
    print("="*78)

    for k in K_VALUES:

        for ell in range(
            k+2,
            36,
            2
        ):

            max_s = (
                ell-1
            )//2

            for s in range(
                k+1,
                max_s+1
            ):

                exact = B_exact(
                    k,
                    ell,
                    s
                )

                generic = B_generic(
                    k,
                    ell,
                    s
                )

                residual = sp.factor(
                    exact-generic
                )

                if residual != 0:

                    failures_B_interior += 1

                    print(
                        "FAIL:",
                        f"k={k}, ell={ell}, s={s},",
                        f"exact={exact},",
                        f"generic={generic},",
                        f"residual={residual}"
                    )

    print(
        "interior B failures =",
        failures_B_interior
    )

    # ------------------------------------------------------------------------
    # 3. C EDGE LAYERS
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("3. C EDGE LAYERS")
    print("="*78)

    for k in K_VALUES:

        print()
        print(f"k={k}")

        for ell in range(
            k+2,
            36,
            2
        ):

            for s in [
                k-1,
                k
            ]:

                if s < 0:
                    continue

                exact = C_exact(
                    k,
                    ell,
                    s
                )

                generic = C_generic(
                    k,
                    ell,
                    s
                )

                residual = sp.factor(
                    exact-generic
                )

                print()
                print(
                    f"  ell={ell}, s={s}"
                )
                print(
                    "    active =",
                    C_active(
                        k,
                        ell,
                        s
                    )
                )
                print(
                    "    exact =",
                    exact
                )
                print(
                    "    generic =",
                    generic
                )
                print(
                    "    residual =",
                    residual
                )

                if residual != 0:

                    failures_C_edge += 1

                    candidates = search_C_candidate(
                        k,
                        ell,
                        s
                    )

                    print(
                        "    candidates =",
                        candidates[:8]
                    )

    # ------------------------------------------------------------------------
    # 4. SEARCH FOR EDGE PATTERN BY k
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("4. EDGE TABLE")
    print("="*78)

    for k in K_VALUES:

        print()
        print(f"k={k}")

        for ell in range(
            k+4,
            32,
            2
        ):

            b = B_exact(
                k,
                ell,
                k
            )

            c1 = C_exact(
                k,
                ell,
                k-1
            )

            c2 = C_exact(
                k,
                ell,
                k
            )

            print(
                f"  ell={ell}: "
                f"B_k={b}, "
                f"C_(k-1)={c1}, "
                f"C_k={c2}"
            )

    # ------------------------------------------------------------------------
    # 5. COMPLETE D WITH EDGE CORRECTION
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("5. COMPLETE D RECONSTRUCTION")
    print("="*78)

    for k in K_VALUES:

        for ell in range(
            k+4,
            36,
            2
        ):

            max_s = (
                ell-1
            )//2

            for s in range(
                max_s+1
            ):

                exact = D_exact(
                    k,
                    ell,
                    s
                )

                closed = D_fixed(
                    k,
                    ell,
                    s
                )

                residual = sp.factor(
                    exact-closed
                )

                if residual != 0:

                    failures_D += 1

                    print(
                        "FAIL:",
                        f"k={k}, ell={ell}, s={s},",
                        f"exact={exact},",
                        f"closed={closed},",
                        f"residual={residual}"
                    )

    print(
        "complete D failures =",
        failures_D
    )

    # ------------------------------------------------------------------------
    # 6. LOW-ORDER STRUCTURE
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("6. LOW-ORDER STRUCTURE")
    print("="*78)

    for k in K_VALUES:

        ell = k+10

        if ell % 2 == 0:
            ell += 1

        print(
            f"k={k}, ell={ell}"
        )

        for s in range(
            max(0,k-1),
            k+2
        ):

            exact = D_exact(
                k,
                ell,
                s
            )

            print(
                f"  s={s}: D={exact}"
            )

    # ------------------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("FINAL DIAGNOSTIC")
    print("="*78)

    print(
        "B edge failures =",
        failures_B_edge
    )

    print(
        "B interior failures =",
        failures_B_interior
    )

    print(
        "C edge discrepancies =",
        failures_C_edge
    )

    print(
        "D reconstruction failures =",
        failures_D
    )

    if (
        failures_B_edge == 0
        and failures_B_interior == 0
        and failures_D == 0
    ):

        print()
        print(
            "STATUS = PASS/PARTIAL"
        )

        print()
        print(
            "The generic interior B law is exact,"
        )
        print(
            "and the exceptional s=k layer is exactly"
        )
        print(
            "represented by the Chu-Vandermonde edge form."
        )

        print()
        print(
            "The remaining problem is the C boundary edge."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "derive symbolic formulas for C_(k,ell,k-1)"
        )
        print(
            "and C_(k,ell,k), then combine them with"
        )
        print(
            "the B edge to obtain a genuinely uniform"
        )
        print(
            "odd-k coefficient theorem."
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

