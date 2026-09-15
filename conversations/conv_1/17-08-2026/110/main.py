#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 173
k=3 COMPLETE B+C COEFFICIENT COMPRESSION
==============================================================================

GOAL
----
Experiment 172 proved the exact k=3 B-branch Chu-Vandermonde collapse.

Experiment 171 already gives the exact C-boundary contribution.

This experiment combines them and searches for an explicit closed formula
for

    D_(3,ell,s) = B_(3,ell,s) + C_(3,ell,s).

The experiment is deliberately conservative:

  * exact integer arithmetic only;
  * no interpolation;
  * no polynomial fitting;
  * no guessing from numerical values;
  * every proposed reduction is verified symbolically/numerically against
    the original finite kernel.

We test:

  1. exact B closed form;
  2. exact C finite boundary;
  3. simplified C candidate families;
  4. exact D coefficients;
  5. leading-coefficient law
         lead(P_(3,ell)) = (ell-6)^2
     after factoring the forced N^2;
  6. forward ell holdout.

The main purpose is to determine whether D_(3,ell,s) admits the same compact
binomial structure as the k=1 family, or whether the k=3 case naturally
requires a piecewise boundary correction.

==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Universal V kernel
# ============================================================================

def V(r, m):
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
# Exact k=3 B branch
# ============================================================================

def B_term(ell, s, j):
    k = 3
    m = s - k
    r = j - 2 * k

    if m < 0 or r < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        sp.binomial(ell, k + j) * V(r, m)
    )


def B_support(ell, s):
    k = 3
    m = s - k

    if m < 0:
        return []

    j_min = 2 * k + 2 * m
    j_max = ell - k

    if j_min > j_max:
        return []

    return list(range(j_min, j_max + 1))


def B_exact(ell, s):
    return sp.factor(
        sum(
            B_term(ell, s, j)
            for j in B_support(ell, s)
        )
    )


# ============================================================================
# Exact k=3 C boundary
# ============================================================================

def C_term(ell, s, j):
    k = 3

    if not (ell - k <= j <= ell):
        return sp.Integer(0)

    r = j - 2 * k

    if r < 0:
        return sp.Integer(0)

    # t-degree condition
    # e + 2m - (r+1) - 1 = 2s
    e = 2 * ell - j
    numerator = e - r - 2

    if numerator % 2:
        return sp.Integer(0)

    m = s - numerator // 2

    if m < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        -sp.binomial(k, ell - j) * V(r, m)
    )


def C_support(ell, s):
    return [
        j
        for j in range(ell - 3, ell + 1)
        if C_term(ell, s, j) != 0
    ]


def C_exact(ell, s):
    return sp.factor(
        sum(
            C_term(ell, s, j)
            for j in C_support(ell, s)
        )
    )


# ============================================================================
# Exact complete coefficient
# ============================================================================

def D_exact(ell, s):
    return sp.factor(
        B_exact(ell, s) + C_exact(ell, s)
    )


# ============================================================================
# B closed form from Experiment 172
# ============================================================================

def B_closed(ell, s):
    if s < 3:
        return sp.Integer(0)

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    first = (
        (-1) ** (s - 3)
        * sp.binomial(ell, 2 * s + 3)
        * sp.rf(s + 6, L)
        / sp.rf(2 * s + 4, L)
    )

    second = sp.Integer(0)

    if s >= 4:
        second = (
            (-1) ** (s - 3)
            * sp.binomial(ell, 2 * s + 3)
            * sp.rf(s + 7, L)
            / sp.rf(2 * s + 4, L)
        )

    return sp.factor(first + second)


# ============================================================================
# Candidate C simplification
# ============================================================================
#
# We do NOT assume one universal C formula.
#
# First derive the exact boundary representation in binomial notation.
# Then test several structurally motivated compressions.
#
# The raw C branch has at most four boundary j-values:
#
#     j = ell-3, ell-2, ell-1, ell.
#
# For each s we print the nonzero contributions and attempt to compress
# them through Pascal identities.
# ============================================================================

def C_boundary_terms(ell, s):
    out = []

    for j in range(ell - 3, ell + 1):
        value = C_term(ell, s, j)

        if value != 0:
            out.append((j, sp.factor(value)))

    return out


# ============================================================================
# Symbolic boundary simplification
# ============================================================================

def simplify_C_by_ell(ell, s):
    """
    Keep ell numeric but expose the exact dependence on s.
    This is a diagnostic rather than a fitted formula.
    """

    exact = C_exact(ell, s)

    return sp.factor(exact)


# ============================================================================
# Candidate tests
# ============================================================================

def candidate_C_family_A(ell, s):
    """
    Natural k-shift candidate:
        (-1)^s C(ell-3, s) + correction.
    """

    return sp.factor(
        (-1) ** s * sp.binomial(ell - 3, s)
    )


def candidate_C_family_B(ell, s):
    """
    Neighboring Pascal combination.
    """

    return sp.factor(
        (-1) ** s * (
            sp.binomial(ell - 3, s)
            + sp.binomial(ell - 4, s - 1)
        )
    )


def candidate_C_family_C(ell, s):
    """
    Binomial-square style candidate suggested by the four boundary modes.
    """

    return sp.factor(
        (-1) ** s * (
            sp.binomial(ell - 3, s)
            - sp.binomial(ell - 5, s - 2)
        )
    )


# ============================================================================
# D closed form using exact C, pending simplification
# ============================================================================

def D_from_B_and_C(ell, s):
    return sp.factor(
        B_closed(ell, s)
        + C_exact(ell, s)
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 173")
    print("k=3 COMPLETE B+C COEFFICIENT COMPRESSION")
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
    candidate_A_failures = 0
    candidate_B_failures = 0
    candidate_C_failures = 0

    # ------------------------------------------------------------------------
    # 1. B closed certificate
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. B-CLOSED VERIFICATION")
    print("=" * 78)

    for ell in TRAIN:

        max_s = (ell - 1) // 2

        print()
        print(f"ell={ell}")

        for s in range(0, max_s + 1):

            exact = B_exact(ell, s)
            closed = B_closed(ell, s)
            residual = sp.factor(exact - closed)

            print(
                f"  s={s}: "
                f"exact={exact}, "
                f"closed={closed}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 2. C boundary anatomy
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. C-BOUNDARY ANATOMY")
    print("=" * 78)

    for ell in TRAIN:

        max_s = (ell - 1) // 2

        print()
        print(f"ell={ell}")

        for s in range(0, max_s + 1):

            terms = C_boundary_terms(
                ell,
                s
            )

            total = C_exact(
                ell,
                s
            )

            print(
                f"  s={s}: "
                f"terms={terms}, "
                f"C={total}"
            )

    # ------------------------------------------------------------------------
    # 3. Candidate C compression search
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. STRUCTURAL C-COMPRESSION TEST")
    print("=" * 78)

    for ell in TRAIN:

        max_s = (ell - 1) // 2

        print()
        print(f"ell={ell}")

        for s in range(0, max_s + 1):

            exact = C_exact(
                ell,
                s
            )

            A = candidate_C_family_A(
                ell,
                s
            )

            B = candidate_C_family_B(
                ell,
                s
            )

            C = candidate_C_family_C(
                ell,
                s
            )

            ra = sp.factor(
                exact - A
            )

            rb = sp.factor(
                exact - B
            )

            rc = sp.factor(
                exact - C
            )

            if ra != 0:
                candidate_A_failures += 1

            if rb != 0:
                candidate_B_failures += 1

            if rc != 0:
                candidate_C_failures += 1

            print(
                f"  s={s}: "
                f"C={exact} | "
                f"Ares={ra} | "
                f"Bres={rb} | "
                f"Cres={rc}"
            )

    # ------------------------------------------------------------------------
    # 4. Complete D coefficients
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. COMPLETE D COEFFICIENT")
    print("=" * 78)

    for ell in TRAIN:

        max_s = (ell - 1) // 2

        print()
        print(f"ell={ell}")

        for s in range(
            0,
            max_s + 1
        ):

            exact = D_exact(
                ell,
                s
            )

            via_closed_B = D_from_B_and_C(
                ell,
                s
            )

            residual = sp.factor(
                exact - via_closed_B
            )

            print(
                f"  s={s}: "
                f"D={exact}, "
                f"Bclosed+C={via_closed_B}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 5. Leading coefficient law
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. LEADING COEFFICIENT TEST")
    print("=" * 78)

    N = sp.symbols("N")

    leading_failures = 0

    for ell in TRAIN + HOLDOUT:

        max_s = (ell - 1) // 2

        poly = sp.expand(
            sum(
                D_exact(
                    ell,
                    s
                ) * N ** s
                for s in range(
                    max_s + 1
                )
            )
        )

        reduced = sp.expand(
            sp.cancel(
                poly / N**2
            )
        )

        degree = sp.degree(
            reduced,
            N
        )

        leading = sp.LC(
            sp.Poly(
                reduced,
                N
            )
        )

        predicted = sp.Integer(
            (ell - 6) ** 2
        )

        residual = sp.factor(
            leading - predicted
        )

        print(
            f"ell={ell}: "
            f"degree(P)={degree}, "
            f"leading={leading}, "
            f"predicted=(ell-6)^2={predicted}, "
            f"residual={residual}"
        )

        if residual != 0:
            leading_failures += 1

    # ------------------------------------------------------------------------
    # 6. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FORWARD ELL HOLDOUT")
    print("=" * 78)

    holdout_failures = 0

    for ell in HOLDOUT:

        max_s = (ell - 1) // 2

        ok = True

        for s in range(
            0,
            max_s + 1
        ):

            residual = sp.factor(
                D_exact(
                    ell,
                    s
                )
                -
                D_from_B_and_C(
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
    # 7. Representative full polynomials
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. REPRESENTATIVE k=3 POLYNOMIALS")
    print("=" * 78)

    for ell in [
        11, 13, 15, 17, 21
    ]:

        max_s = (ell - 1) // 2

        poly = sp.factor(
            sum(
                D_exact(
                    ell,
                    s
                ) * N**s
                for s in range(
                    max_s + 1
                )
            )
        )

        print(
            f"ell={ell}: D_(3,{ell})(N) = {poly}"
        )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "B closed failures =",
        failures
    )

    print(
        "C candidate A failures =",
        candidate_A_failures
    )

    print(
        "C candidate B failures =",
        candidate_B_failures
    )

    print(
        "C candidate C failures =",
        candidate_C_failures
    )

    print(
        "leading-law failures =",
        leading_failures
    )

    print(
        "holdout failures =",
        holdout_failures
    )

    if failures == 0 and holdout_failures == 0:

        print()
        print("STATUS = PASS")

        print()
        print(
            "The collapsed k=3 B branch reproduces the exact"
        )
        print(
            "finite B kernel, and B+C reproduces every exact"
        )
        print(
            "k=3 coefficient."
        )

        print()
        if leading_failures == 0:
            print(
                "The leading-coefficient conjecture"
            )
            print(
                "    lead(P_(3,ell)) = (ell-6)^2"
            )
            print(
                "also survives the tested forward ell."
            )
        else:
            print(
                "The leading-coefficient conjecture is NOT"
            )
            print(
                "yet proved by this dataset."
            )

        print()
        print("NEXT TARGET:")
        print(
            "derive the exact closed C-boundary formula for k=3"
        )
        print(
            "and simplify B+C into a single coefficient formula."
        )
        print(
            "After that, compare the resulting k=1 and k=3"
        )
        print(
            "laws to isolate the general odd-k pattern."
        )

    else:

        print()
        print("STATUS = FAIL")

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

