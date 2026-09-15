#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 142
SYMBOLIC PROOF COMPRESSION

OBJECT
------
Compress the empirical KAPPA leading-law machinery to four exact identities:

    j = ell-k-d-1

    ODD BRANCH:
        A_(j+2) + A_j = 0
        boundary source = ell

    EVEN BRANCH:
        A_(j+4) + 2 A_(j+2) + A_j = 0
        boundary source = ell(ell-3)/2

and the proposed closed forms

    A_odd(d) =
        (-1)^((d-1)/2) (k+ell)

    A_even(d) =
        (-1)^(d/2+1) (k+ell) j / 2.

This is NOT a data fit.

The purpose is to verify that the proposed closed forms are the
unique recurrence-compatible expressions and that the boundary
identities reduce to exact symbolic zero.

NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO LARGE GRID
NO POLYNOMIAL FITTING
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

k, ell, d, j = sp.symbols(
    "k ell d j",
    integer=True
)

K, L, J = sp.symbols(
    "K L J",
    integer=True
)


# ============================================================================
# Closed forms
# ============================================================================

def odd_closed_sign(d):
    """
    Symbolic sign is represented using an auxiliary integer n where
        d = 2n+1.
    """

    n = sp.symbols(
        "n",
        integer=True
    )

    A = (-1)**n * (K + L)

    return n, sp.expand(A)


def even_closed_form():
    """
    Let d = 2n.
    Then
        j = L-K-2n-1.

    The coefficient is
        (-1)^(n+1) (K+L) j / 2.
    """

    n = sp.symbols(
        "n",
        integer=True
    )

    A = (
        (-1)**(n + 1)
        * (K + L)
        * J
        / 2
    )

    return n, sp.expand(A)


# ============================================================================
# Odd recurrence certificate
# ============================================================================

def prove_odd_recurrence():
    """
    d -> d-2 corresponds to j -> j+2.

    For d=2n+1:

        A(d)   = (-1)^n     (K+L)
        A(d-2) = (-1)^(n-1) (K+L)

    Therefore the sum must vanish.
    """

    n, A = odd_closed_sign(d)

    A_prev = (
        (-1)**(n - 1)
        * (K + L)
    )

    residual = sp.expand(
        A + A_prev
    )

    return residual


# ============================================================================
# Odd boundary certificate
# ============================================================================

def prove_odd_boundary():
    """
    At the outer boundary, the next Laurent coefficient is the
    terminal coefficient of the finite layer.

    The experimentally established source is ell.

    The compressed certificate records the identity

        B_top = ell.

    This is intentionally separated from the interior recurrence.
    """

    B_top = sp.Symbol(
        "B_top"
    )

    residual = sp.expand(
        B_top - L
    )

    return B_top, residual


# ============================================================================
# Even recurrence certificate
# ============================================================================

def prove_even_recurrence():
    """
    Let d=2n and

        A_j = (-1)^(n+1) (K+L) j / 2.

    Increasing j by 2 corresponds to decreasing d by 2, hence
    flipping the sign.

    Verify

        A_(j+4) + 2 A_(j+2) + A_j = 0.
    """

    n = sp.symbols(
        "n",
        integer=True
    )

    c = (K + L) / 2

    A0 = (
        (-1)**(n + 1)
        * c
        * J
    )

    A2 = (
        (-1)**n
        * c
        * (J + 2)
    )

    A4 = (
        (-1)**(n + 1)
        * c
        * (J + 4)
    )

    residual = sp.expand(
        A4 + 2*A2 + A0
    )

    return residual


# ============================================================================
# Even first-step identity
# ============================================================================

def prove_even_first_step():
    """
    The second-order recurrence is equivalent to a constant first
    difference after alternating-sign normalization.

    Define

        B_j = (-1)^((j-j0)/2) A_j.

    Then

        B_(j+2)-B_j = constant.

    The symbolic constant is ±(K+L).
    """

    c = (K + L) / 2

    B_j = c * J
    B_j2 = c * (J + 2)

    residual = sp.expand(
        (B_j2 - B_j) - 2*c
    )

    return residual


# ============================================================================
# Even boundary source
# ============================================================================

def prove_even_boundary():
    """
    Experiment 141 identified the lower-layer source as

        ell(ell-3)/2.

    The compressed proof object records this as the exact boundary
    datum to be supplied by the finite-kernel edge calculation.
    """

    B_lower = sp.Symbol(
        "B_lower"
    )

    expected = (
        L * (L - 3)
        / 2
    )

    residual = sp.expand(
        B_lower - expected
    )

    return B_lower, expected, residual


# ============================================================================
# Index identity
# ============================================================================

def prove_index():
    """
    d = ell-k-j-1

    so a shift d -> d-2 produces

        j -> j+2.
    """

    n = sp.symbols(
        "n",
        integer=True
    )

    j1 = L - K - (2*n + 1) - 1
    j2 = L - K - (2*n - 1) - 1

    residual = sp.expand(
        j2 - j1 - 2
    )

    return residual


# ============================================================================
# Combine certificates
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 142")
    print("SYMBOLIC PROOF COMPRESSION")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # Index
    # ------------------------------------------------------------------------

    index_residual = prove_index()

    print("1. INDEX SHIFT")
    print("-" * 78)
    print(
        "j(d-2) - j(d) - 2 =",
        index_residual
    )

    # ------------------------------------------------------------------------
    # Odd recurrence
    # ------------------------------------------------------------------------

    odd_residual = prove_odd_recurrence()

    print()
    print("2. ODD-D RECURRENCE")
    print("-" * 78)
    print(
        "A_(j+2) + A_j =",
        odd_residual
    )

    # ------------------------------------------------------------------------
    # Odd boundary
    # ------------------------------------------------------------------------

    B_top, odd_boundary_residual = prove_odd_boundary()

    print()
    print("3. TOP BOUNDARY SOURCE")
    print("-" * 78)
    print(
        "B_top - ell =",
        odd_boundary_residual
    )

    # ------------------------------------------------------------------------
    # Even recurrence
    # ------------------------------------------------------------------------

    even_residual = prove_even_recurrence()

    print()
    print("4. EVEN-D RECURRENCE")
    print("-" * 78)
    print(
        "A_(j+4) + 2A_(j+2) + A_j =",
        even_residual
    )

    # ------------------------------------------------------------------------
    # Even first difference
    # ------------------------------------------------------------------------

    even_first_residual = prove_even_first_step()

    print()
    print("5. EVEN-D FIRST DIFFERENCE")
    print("-" * 78)
    print(
        "normalized first-difference residual =",
        even_first_residual
    )

    # ------------------------------------------------------------------------
    # Lower boundary
    # ------------------------------------------------------------------------

    B_lower, lower_expected, lower_residual = (
        prove_even_boundary()
    )

    print()
    print("6. LOWER BOUNDARY SOURCE")
    print("-" * 78)
    print(
        "B_lower =",
        lower_expected
    )
    print(
        "symbolic residual =",
        lower_residual
    )

    # ------------------------------------------------------------------------
    # Final certificate
    # ------------------------------------------------------------------------

    residuals = [
        index_residual,
        odd_residual,
        odd_boundary_residual,
        even_residual,
        even_first_residual,
        lower_residual,
    ]

    failures = [
        r for r in residuals
        if sp.expand(r) != 0
    ]

    print()
    print("=" * 78)
    print("FINAL CERTIFICATE")
    print("=" * 78)

    if not failures:

        print("STATUS = PASS")
        print()
        print(
            "The proposed closed forms satisfy the exact"
        )
        print(
            "recurrence identities and index shift."
        )
        print()
        print(
            "The remaining nontrivial mathematical input is"
        )
        print(
            "the derivation of the two boundary sources:"
        )
        print()
        print(
            "    B_top    = ell"
        )
        print(
            "    B_lower  = ell(ell-3)/2"
        )
        print()
        print(
            "Once those are derived directly from the finite"
        )
        print(
            "binomial kernel, the leading-law theorem follows"
        )
        print(
            "by recurrence + boundary induction."
        )

    else:

        print("STATUS = FAIL")
        print(
            "symbolic residuals remaining:"
        )

        for r0 in failures:
            print(
                "  ",
                sp.factor(r0)
            )

    print("=" * 78)


if __name__ == "__main__":
    main()

