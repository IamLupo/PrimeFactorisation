#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 168
GENERAL k=1 B-BRANCH: CHU-VANDERMONDE COLLAPSE
==============================================================================

GOAL
----
Turn the exact finite B-branch convolution from Experiment 167 into a
closed symbolic formula for arbitrary coefficient layer s.

For k=1:

    n = 2s-2+r

and

    L = ell-2s-1.

The exact B summand is

    (-1)^(s-1+r)
    C(ell, 2s+1+r)
    [
        C(s-1+r, s-1)
        +
        C(s-2+r, s-2)
    ]

with the second term absent when s=1.

After factoring the first binomial:

    B1 =
      (-1)^(s-1) C(ell,2s+1)
      sum_r
        (-L)_r (s)_r
        / ((2s+2)_r r!)

and, for s>=2,

    B2 =
      (-1)^(s-1) C(ell,2s+1)
      sum_r
        (-L)_r (s-1)_r
        / ((2s+2)_r r!).

These are terminating 2F1(1) sums.

Chu-Vandermonde:

    2F1(-L,a;c;1)
       = (c-a)_L / (c)_L.

Therefore the predicted closed form is

    B =
      (-1)^(s-1) C(ell,2s+1)
      [
        (s+2)_L/(2s+2)_L
        +
        (s+3)_L/(2s+2)_L
      ]

for s>=2,

and only the first term for s=1.

THIS EXPERIMENT VERIFIES:
    1. exact finite sum = symbolic two-2F1 decomposition;
    2. each 2F1 sum = Chu-Vandermonde value;
    3. combined closed formula = exact B coefficient;
    4. factorized/rational equivalent forms;
    5. forward ell holdout.

NO FITTING
NO POLYNOMIAL INTERPOLATION
NO CSV
NO SKLEARN
NO FACTOR-PAIR SEARCH
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp

# Symbols
ell = sp.symbols("ell", integer=True, positive=True)
s = sp.symbols("s", integer=True, positive=True)
r = sp.symbols("r", integer=True, nonnegative=True)


# ============================================================================
# UNIVERSAL V COEFFICIENT
# ============================================================================

def v_coeff(R, m):
    if R < 0 or m < 0 or m > R // 2:
        return sp.Integer(0)

    first = sp.binomial(R - m, m)

    second = (
        sp.Integer(0)
        if m == 0
        else sp.binomial(R - m - 1, m - 1)
    )

    return sp.expand(
        (-1) ** (R - m) * (first + second)
    )


# ============================================================================
# EXACT k=1 B TERM
# ============================================================================

def B_term(ell0, s0, n0):
    """
    Experiment-167 exact branch-B term.

    k=1
    m=s-1
    r=n
    """

    m0 = s0 - 1

    if m0 < 0:
        return sp.Integer(0)

    R = n0

    if m0 > R // 2:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (R - m0)
        * sp.binomial(
            ell0,
            R + 3
        )
        * (
            sp.binomial(
                R - m0,
                m0
            )
            +
            (
                sp.Integer(0)
                if m0 == 0
                else sp.binomial(
                    R - m0 - 1,
                    m0 - 1
                )
            )
        )
    )


# ============================================================================
# SUPPORT
# ============================================================================

def B_support(ell0, s0):
    m0 = s0 - 1

    if m0 < 0:
        return []

    lo = 2 * m0
    hi = ell0 - 3

    if lo > hi:
        return []

    return list(range(lo, hi + 1))


# ============================================================================
# EXACT FINITE B TOTAL
# ============================================================================

def B_exact(ell0, s0):
    return sp.factor(
        sum(
            B_term(
                ell0,
                s0,
                n0
            )
            for n0 in B_support(
                ell0,
                s0
            )
        )
    )


# ============================================================================
# SHIFTED FINITE SUM
# ============================================================================

def shifted_first_sum(ell0, s0):
    """
    First V piece.

    r = 0,...,L
    L = ell-2s-1
    """

    L = ell0 - 2*s0 - 1

    if L < 0:
        return sp.Integer(0)

    prefactor = (
        (-1) ** (s0 - 1)
        * sp.binomial(
            ell0,
            2*s0 + 1
        )
    )

    total = sp.Integer(0)

    for rr in range(L + 1):

        total += (
            sp.rf(-L, rr)
            * sp.rf(s0, rr)
            /
            (
                sp.rf(
                    2*s0 + 2,
                    rr
                )
                * sp.factorial(rr)
            )
        )

    return sp.factor(
        prefactor * total
    )


def shifted_second_sum(ell0, s0):
    """
    Second V piece.

    Empty when s=1.
    """

    if s0 < 2:
        return sp.Integer(0)

    L = ell0 - 2*s0 - 1

    if L < 0:
        return sp.Integer(0)

    prefactor = (
        (-1) ** (s0 - 1)
        * sp.binomial(
            ell0,
            2*s0 + 1
        )
    )

    total = sp.Integer(0)

    for rr in range(L + 1):

        total += (
            sp.rf(-L, rr)
            * sp.rf(s0 - 1, rr)
            /
            (
                sp.rf(
                    2*s0 + 2,
                    rr
                )
                * sp.factorial(rr)
            )
        )

    return sp.factor(
        prefactor * total
    )


# ============================================================================
# CHU-VANDERMONDE CLOSED PIECES
# ============================================================================

def chu_first(ell0, s0):
    L = ell0 - 2*s0 - 1

    if L < 0:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (s0 - 1)
        * sp.binomial(
            ell0,
            2*s0 + 1
        )
        * sp.rf(
            s0 + 2,
            L
        )
        /
        sp.rf(
            2*s0 + 2,
            L
        )
    )


def chu_second(ell0, s0):
    if s0 < 2:
        return sp.Integer(0)

    L = ell0 - 2*s0 - 1

    if L < 0:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (s0 - 1)
        * sp.binomial(
            ell0,
            2*s0 + 1
        )
        * sp.rf(
            s0 + 3,
            L
        )
        /
        sp.rf(
            2*s0 + 2,
            L
        )
    )


def chu_closed(ell0, s0):
    return sp.factor(
        chu_first(
            ell0,
            s0
        )
        +
        chu_second(
            ell0,
            s0
        )
    )


# ============================================================================
# DIRECT SYMBOLIC 2F1 CHECK
# ============================================================================

def hyper_first(ell0, s0):
    L = ell0 - 2*s0 - 1

    if L < 0:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (s0 - 1)
        * sp.binomial(
            ell0,
            2*s0 + 1
        )
        * sp.hyper(
            [
                -L,
                s0
            ],
            [
                2*s0 + 2
            ],
            1
        )
    )


def hyper_second(ell0, s0):
    if s0 < 2:
        return sp.Integer(0)

    L = ell0 - 2*s0 - 1

    if L < 0:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (s0 - 1)
        * sp.binomial(
            ell0,
            2*s0 + 1
        )
        * sp.hyper(
            [
                -L,
                s0 - 1
            ],
            [
                2*s0 + 2
            ],
            1
        )
    )


# ============================================================================
# ALL NONNEGATIVE s UP TO SUPPORT
# ============================================================================

def check_pair(ell0, max_s=None):

    if max_s is None:
        max_s = (ell0 - 1) // 2

    failures = []

    for s0 in range(
        1,
        max_s + 1
    ):

        exact = B_exact(
            ell0,
            s0
        )

        finite1 = shifted_first_sum(
            ell0,
            s0
        )

        finite2 = shifted_second_sum(
            ell0,
            s0
        )

        finite_total = sp.factor(
            finite1 + finite2
        )

        chu1 = chu_first(
            ell0,
            s0
        )

        chu2 = chu_second(
            ell0,
            s0
        )

        chu_total = chu_closed(
            ell0,
            s0
        )

        hyper_total = sp.factor(
            hyper_first(
                ell0,
                s0
            )
            +
            hyper_second(
                ell0,
                s0
            )
        )

        checks = [
            (
                "finite-vs-exact",
                finite_total - exact
            ),
            (
                "chu-first",
                finite1 - chu1
            ),
            (
                "chu-second",
                finite2 - chu2
            ),
            (
                "chu-total",
                chu_total - exact
            ),
            (
                "hyper-total",
                hyper_total - exact
            )
        ]

        for label, residual in checks:

            residual = sp.factor(
                sp.simplify(
                    residual
                )
            )

            if residual != 0:
                failures.append(
                    (
                        s0,
                        label,
                        residual
                    )
                )

    return failures


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 168")
    print("GENERAL k=1 B-BRANCH: CHU-VANDERMONDE COLLAPSE")
    print("=" * 78)
    print()

    TRAIN = [
        9,
        11,
        13,
        15,
        17,
        19,
        21
    ]

    HOLDOUT = [
        23,
        25,
        27
    ]

    total_failures = 0

    # ------------------------------------------------------------------------
    # 1. Formula statement
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. SYMBOLIC FORM")
    print("=" * 78)

    print()
    print(
        "L = ell - 2*s - 1"
    )

    print(
        "B1 = (-1)^(s-1) C(ell,2s+1)"
    )

    print(
        "     * 2F1(-L,s;2s+2;1)"
    )

    print()
    print(
        "B2 = (-1)^(s-1) C(ell,2s+1)"
    )

    print(
        "     * 2F1(-L,s-1;2s+2;1)"
    )

    print()
    print(
        "B = B1 + B2"
    )

    # ------------------------------------------------------------------------
    # 2. Training
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. TRAINING")
    print("=" * 78)

    for ell0 in TRAIN:

        failures = check_pair(
            ell0
        )

        print()
        print(
            f"ell={ell0}: "
            f"{'PASS' if not failures else 'FAIL'}"
        )

        if failures:

            for item in failures:
                print(
                    "  ",
                    item
                )

            total_failures += len(
                failures
            )

        else:

            for s0 in range(
                1,
                (ell0-1)//2 + 1
            ):

                exact = B_exact(
                    ell0,
                    s0
                )

                closed = chu_closed(
                    ell0,
                    s0
                )

                print(
                    f"  s={s0}: "
                    f"B={exact}, "
                    f"closed={closed}"
                )

    # ------------------------------------------------------------------------
    # 3. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    holdout_failures = 0

    for ell0 in HOLDOUT:

        failures = check_pair(
            ell0
        )

        print(
            f"ell={ell0}: "
            f"{'PASS' if not failures else 'FAIL'}"
        )

        if failures:
            holdout_failures += len(
                failures
            )

            for item in failures:
                print(
                    "  ",
                    item
                )

    # ------------------------------------------------------------------------
    # 4. Factorial simplification
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FACTORIAL FORM")
    print("=" * 78)

    Lsym = sp.symbols(
        "L",
        integer=True,
        nonnegative=True
    )

    ssym = sp.symbols(
        "s",
        integer=True,
        positive=True
    )

    # Formal ratios:
    ratio1 = sp.simplify(
        sp.rf(ssym+2, Lsym)
        /
        sp.rf(2*ssym+2, Lsym)
    )

    ratio2 = sp.simplify(
        sp.rf(ssym+3, Lsym)
        /
        sp.rf(2*ssym+2, Lsym)
    )

    print(
        "(s+2)_L/(2s+2)_L =",
        ratio1
    )

    print(
        "(s+3)_L/(2s+2)_L =",
        ratio2
    )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "training failures =",
        total_failures
    )

    print(
        "holdout failures =",
        holdout_failures
    )

    if (
        total_failures == 0
        and holdout_failures == 0
    ):

        print()
        print("STATUS = PASS")

        print()
        print(
            "The complete k=1 B coefficient layer is"
        )

        print(
            "represented by two terminating 2F1(1) sums."
        )

        print()
        print(
            "Both sums collapse exactly by"
        )

        print(
            "Chu-Vandermonde."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "combine the closed B_(1,ell,s) formula with"
        )

        print(
            "the exact C-boundary coefficient and derive"
        )

        print(
            "a closed formula for [N^s] D_(1,ell)."
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

