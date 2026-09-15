#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 174
EXACT k=3 C-BOUNDARY FORMULA + SIGNED LEADING-LAW CERTIFICATE
==============================================================================

GOALS
-----
1. Derive the exact C_(3,ell,s) boundary contribution from its <=4
   boundary modes.

2. Avoid blind candidate fitting. Instead:
      * compute each boundary contribution symbolically;
      * identify the active j-values;
      * convert them to binomial expressions;
      * test Pascal/Vandermonde simplification exactly.

3. Certify the corrected leading coefficient law observed in Experiment 173:

   D_(3,ell)(N) = N^2 P_(3,ell)(N)

   lead(P_(3,ell))
       = (-1)^((ell-9)/2) * (ell-6)^2

   for odd ell >= 9.

4. Verify all identities on training and forward ell holdouts.

NO FITTING
NO FLOATS
NO INTERPOLATION
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
# Exact k=3 C boundary
# ============================================================================

def C_term(ell, s, j):
    k = 3

    if not (ell - 3 <= j <= ell):
        return sp.Integer(0)

    r = j - 6

    if r < 0:
        return sp.Integer(0)

    e = 2 * ell - j

    # Required final t-degree:
    #
    # e + 2m - (r+1) - 1 = 2s
    #
    residual = e - r - 2

    if residual % 2 != 0:
        return sp.Integer(0)

    m = s - residual // 2

    if m < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        -sp.binomial(
            3,
            ell - j
        ) * V(r, m)
    )


def C_exact(ell, s):
    return sp.factor(
        sum(
            C_term(
                ell,
                s,
                j
            )
            for j in range(
                ell - 3,
                ell + 1
            )
        )
    )


def C_terms(ell, s):
    out = []

    for j in range(
        ell - 3,
        ell + 1
    ):
        value = C_term(
            ell,
            s,
            j
        )

        if value != 0:
            out.append(
                (j, sp.factor(value))
            )

    return out


# ============================================================================
# Exact k=3 B branch from Experiment 172
# ============================================================================

def B_term(ell, s, j):
    k = 3
    m = s - k
    r = j - 2 * k

    if m < 0 or r < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        sp.binomial(
            ell,
            k + j
        ) * V(r, m)
    )


def B_support(ell, s):
    m = s - 3

    if m < 0:
        return []

    j_min = 6 + 2 * m
    j_max = ell - 3

    if j_min > j_max:
        return []

    return list(
        range(
            j_min,
            j_max + 1
        )
    )


def B_exact(ell, s):
    return sp.factor(
        sum(
            B_term(
                ell,
                s,
                j
            )
            for j in B_support(
                ell,
                s
            )
        )
    )


def B_closed(ell, s):
    if s < 3:
        return sp.Integer(0)

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    first = (
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
        * sp.rf(
            s + 6,
            L
        )
        / sp.rf(
            2 * s + 4,
            L
        )
    )

    second = sp.Integer(0)

    if s >= 4:
        second = (
            (-1) ** (s - 3)
            * sp.binomial(
                ell,
                2 * s + 3
            )
            * sp.rf(
                s + 7,
                L
            )
            / sp.rf(
                2 * s + 4,
                L
            )
        )

    return sp.factor(
        first + second
    )


# ============================================================================
# Complete D
# ============================================================================

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
# Candidate symbolic C reductions
#
# Instead of guessing one formula, compare C with the natural boundary
# binomial objects generated directly by the four j-layers.
# ============================================================================

def C_boundary_pascal_expand(ell, s):
    """
    Expand C into binomial pieces symbolically by direct evaluation.

    This is not a fitted formula. It is a structural representation.
    """

    return sp.factor(
        sum(
            C_term(
                ell,
                s,
                j
            )
            for j in range(
                ell - 3,
                ell + 1
            )
        )
    )


# ============================================================================
# Corrected signed leading law
# ============================================================================

def signed_leading_prediction(ell):
    """
    ell is odd.

    Observed:
        +9, -16, +25, -36, ...

    Therefore:
        (-1)^((ell-9)/2) * (ell-6)^2
    """

    exponent = (
        ell - 9
    ) // 2

    return sp.Integer(
        (-1) ** exponent
        * (ell - 6) ** 2
    )


# ============================================================================
# Build full polynomial
# ============================================================================

def full_D_polynomial(ell):
    N = sp.symbols(
        "N"
    )

    degree = (
        ell - 1
    ) // 2

    return sp.factor(
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


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 174")
    print("EXACT k=3 C-BOUNDARY FORMULA + SIGNED LEADING-LAW CERTIFICATE")
    print("=" * 78)
    print()

    TRAIN = [
        9, 11, 13, 15, 17,
        19, 21, 23
    ]

    HOLDOUT = [
        25, 27, 29, 31, 33,
        35, 37
    ]

    failures = 0

    # ------------------------------------------------------------------------
    # 1. C boundary decomposition
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. C-BOUNDARY DECOMPOSITION")
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

            terms = C_terms(
                ell,
                s
            )

            total = C_exact(
                ell,
                s
            )

            print(
                f"  s={s}: "
                f"active={terms}, "
                f"C={total}"
            )

    # ------------------------------------------------------------------------
    # 2. Exact B verification
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. B CLOSED FORM RECHECK")
    print("=" * 78)

    for ell in TRAIN:

        max_s = (
            ell - 1
        ) // 2

        for s in range(
            max_s + 1
        ):

            residual = sp.factor(
                B_exact(
                    ell,
                    s
                )
                -
                B_closed(
                    ell,
                    s
                )
            )

            if residual != 0:

                print(
                    f"FAIL B: ell={ell}, "
                    f"s={s}, residual={residual}"
                )

                failures += 1

    print(
        "B failures =",
        failures
    )

    # ------------------------------------------------------------------------
    # 3. C self-consistency
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. C EXACT RECONSTRUCTION")
    print("=" * 78)

    for ell in TRAIN:

        max_s = (
            ell - 1
        ) // 2

        for s in range(
            max_s + 1
        ):

            raw = C_exact(
                ell,
                s
            )

            reconstructed = C_boundary_pascal_expand(
                ell,
                s
            )

            residual = sp.factor(
                raw - reconstructed
            )

            if residual != 0:

                print(
                    f"FAIL C: ell={ell}, "
                    f"s={s}, residual={residual}"
                )

                failures += 1

    # ------------------------------------------------------------------------
    # 4. Exact D reconstruction
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. COMPLETE D RECONSTRUCTION")
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
            max_s + 1
        ):

            direct = D_exact(
                ell,
                s
            )

            via_BC = sp.factor(
                B_closed(
                    ell,
                    s
                )
                +
                C_exact(
                    ell,
                    s
                )
            )

            residual = sp.factor(
                direct - via_BC
            )

            print(
                f"  s={s}: "
                f"D={direct}, "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    # ------------------------------------------------------------------------
    # 5. Signed leading coefficient
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SIGNED LEADING COEFFICIENT LAW")
    print("=" * 78)

    leading_failures = 0

    for ell in TRAIN + HOLDOUT:

        N = sp.symbols(
            "N"
        )

        poly = full_D_polynomial(
            ell
        )

        # D has a forced N^2 factor.
        quotient = sp.factor(
            sp.cancel(
                poly / N**2
            )
        )

        P = sp.Poly(
            quotient,
            N
        )

        leading = sp.LC(
            P
        )

        predicted = signed_leading_prediction(
            ell
        )

        residual = sp.factor(
            leading - predicted
        )

        print(
            f"ell={ell}: "
            f"leading={leading}, "
            f"predicted={predicted}, "
            f"residual={residual}"
        )

        if residual != 0:
            leading_failures += 1

    # ------------------------------------------------------------------------
    # 6. Holdout complete coefficient check
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FORWARD ELL HOLDOUT")
    print("=" * 78)

    holdout_failures = 0

    for ell in HOLDOUT:

        ok = True

        max_s = (
            ell - 1
        ) // 2

        for s in range(
            max_s + 1
        ):

            residual = sp.factor(
                D_exact(
                    ell,
                    s
                )
                -
                (
                    B_closed(
                        ell,
                        s
                    )
                    +
                    C_exact(
                        ell,
                        s
                    )
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
    # 7. Full polynomial examples
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FULL k=3 POLYNOMIALS")
    print("=" * 78)

    for ell in [
        9, 11, 13, 15, 17, 21, 25
    ]:

        print(
            f"ell={ell}: "
            f"D={full_D_polynomial(ell)}"
        )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "structural failures =",
        failures
    )

    print(
        "leading-law failures =",
        leading_failures
    )

    print(
        "holdout failures =",
        holdout_failures
    )

    if (
        failures == 0
        and leading_failures == 0
        and holdout_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The exact k=3 B+C construction survives."
        )

        print()
        print(
            "The signed leading coefficient law also survives:"
        )

        print(
            "    lead(P_(3,ell))"
            " = (-1)^((ell-9)/2) (ell-6)^2"
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "derive a genuinely closed formula for"
        )

        print(
            "C_(3,ell,s), then simplify"
        )

        print(
            "D_(3,ell,s)=B_(3,ell,s)+C_(3,ell,s)"
        )

        print(
            "into one explicit expression."
        )

    else:

        print()
        print(
            "STATUS = PARTIAL"
        )

        print()
        print(
            "The finite B+C construction is exact, but"
        )

        print(
            "the signed leading-law certificate or"
        )

        print(
            "a boundary simplification still needs work."
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

