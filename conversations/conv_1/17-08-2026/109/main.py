#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 172
k=3 B-BRANCH: EXACT CHU-VANDERMONDE COLLAPSE
==============================================================================

GOAL
----
Experiment 171 established the exact k=3 B branch and the vanishing law

    B_(3,ell,s) = 0,     s < 3.

Set

    m = s - 3,
    r = j - 6.

Then

    B_(3,ell,s)
      = sum_r C(ell,r+9) V_r(m).

For s>=3, shift

    r = 2m + n = 2s-6+n,

so that

    C(ell,r+9) = C(ell,2s+3+n).

The universal Green coefficient becomes

    V_(2m+n)(m)
      = (-1)^(m+n)
        [ C(n+m,m) + C(n+m-1,m-1) ].

Therefore, for s>=4,

    B = (-1)^(s-3) C(ell,2s+3)
        * [
          2F1(-L, s-2; 2s+4; 1)
          +
          2F1(-L, s-3; 2s+4; 1)
        ],

where

    L = ell - 2s - 3.

For s=3, m=0 and the second piece is absent:

    B_(3,ell,3)
      = C(ell,9) * 2F1(-L,1;10;1),
      L=ell-9.

Chu-Vandermonde predicts

    2F1(-L,a;c;1)
      = (c-a)_L / (c)_L.

This experiment checks:

  1. exact shifted finite summand;
  2. exact Pochhammer term representation;
  3. term-by-term agreement;
  4. finite-sum equality;
  5. Chu-Vandermonde equality;
  6. full B closed form;
  7. forward ell holdout.

NO FITTING
NO INTERPOLATION
NO CSV
NO SKLEARN
NO FACTOR-PAIR SEARCH
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Universal V kernel
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
        (-1) ** (r - m) * (a + b)
    )


# ============================================================================
# Exact k=3 B term from the original finite kernel
# ============================================================================

def B_term_exact(ell, s, j):

    k = 3
    m = s - k
    r = j - 2 * k

    if m < 0:
        return sp.Integer(0)

    if r < 0:
        return sp.Integer(0)

    if m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        sp.binomial(
            ell,
            k + j
        )
        * v_coeff(
            r,
            m
        )
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

    return list(
        range(
            j_min,
            j_max + 1
        )
    )


def B_exact(ell, s):

    return sp.factor(
        sum(
            B_term_exact(
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


# ============================================================================
# Shifted exact terms
# ============================================================================

def shifted_B_term(ell, s, n):

    m = s - 3
    r = 2 * m + n
    j = r + 6

    return sp.factor(
        B_term_exact(
            ell,
            s,
            j
        )
    )


# ============================================================================
# Pochhammer representation of individual shifted terms
# ============================================================================

def pochhammer_first_term(ell, s, n):

    if s < 3:
        return sp.Integer(0)

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    prefactor = (
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
    )

    return sp.factor(
        prefactor
        * sp.rf(
            -L,
            n
        )
        * sp.rf(
            s - 2,
            n
        )
        /
        (
            sp.rf(
                2 * s + 4,
                n
            )
            * sp.factorial(n)
        )
    )


def pochhammer_second_term(ell, s, n):

    if s < 4:
        return sp.Integer(0)

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    prefactor = (
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
    )

    return sp.factor(
        prefactor
        * sp.rf(
            -L,
            n
        )
        * sp.rf(
            s - 3,
            n
        )
        /
        (
            sp.rf(
                2 * s + 4,
                n
            )
            * sp.factorial(n)
        )
    )


def pochhammer_shifted_term(ell, s, n):

    return sp.factor(
        pochhammer_first_term(
            ell,
            s,
            n
        )
        +
        pochhammer_second_term(
            ell,
            s,
            n
        )
    )


# ============================================================================
# Finite hypergeometric sums
# ============================================================================

def finite_hyper_first(ell, s):

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    total = sp.Integer(0)

    for n in range(
        L + 1
    ):

        total += (
            sp.rf(
                -L,
                n
            )
            * sp.rf(
                s - 2,
                n
            )
            /
            (
                sp.rf(
                    2 * s + 4,
                    n
                )
                * sp.factorial(n)
            )
        )

    return sp.factor(
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
        * total
    )


def finite_hyper_second(ell, s):

    if s < 4:
        return sp.Integer(0)

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    total = sp.Integer(0)

    for n in range(
        L + 1
    ):

        total += (
            sp.rf(
                -L,
                n
            )
            * sp.rf(
                s - 3,
                n
            )
            /
            (
                sp.rf(
                    2 * s + 4,
                    n
                )
                * sp.factorial(n)
            )
        )

    return sp.factor(
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
        * total
    )


# ============================================================================
# Chu-Vandermonde values
# ============================================================================

def chu_first(ell, s):

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
        * sp.rf(
            s + 6,
            L
        )
        /
        sp.rf(
            2 * s + 4,
            L
        )
    )


def chu_second(ell, s):

    if s < 4:
        return sp.Integer(0)

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    return sp.factor(
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
        * sp.rf(
            s + 7,
            L
        )
        /
        sp.rf(
            2 * s + 4,
            L
        )
    )


def B_closed(ell, s):

    if s < 3:
        return sp.Integer(0)

    return sp.factor(
        chu_first(
            ell,
            s
        )
        +
        chu_second(
            ell,
            s
        )
    )


# ============================================================================
# Independent explicit hypergeometric form
# ============================================================================

def hyper_closed(ell, s):

    if s < 3:
        return sp.Integer(0)

    L = ell - 2 * s - 3

    if L < 0:
        return sp.Integer(0)

    first = sp.hyper(
        [
            -L,
            s - 2
        ],
        [
            2 * s + 4
        ],
        1
    )

    second = (
        sp.Integer(0)
        if s < 4
        else sp.hyper(
            [
                -L,
                s - 3
            ],
            [
                2 * s + 4
            ],
            1
        )
    )

    return sp.factor(
        (-1) ** (s - 3)
        * sp.binomial(
            ell,
            2 * s + 3
        )
        * (
            first
            +
            second
        )
    )


# ============================================================================
# Test one ell
# ============================================================================

def check_ell(ell):

    failures = []

    max_s = (
        ell - 1
    ) // 2

    # s < 3
    for s in range(
        0,
        3
    ):

        exact = B_exact(
            ell,
            s
        )

        closed = B_closed(
            ell,
            s
        )

        if sp.factor(
            exact - closed
        ) != 0:

            failures.append(
                (
                    s,
                    "vanishing",
                    sp.factor(
                        exact - closed
                    )
                )
            )

    # s >= 3
    for s in range(
        3,
        max_s + 1
    ):

        support = B_support(
            ell,
            s
        )

        if not support:
            continue

        L = ell - 2 * s - 3

        # ---------------------------------------------------------------
        # Term-by-term Pochhammer certificate
        # ---------------------------------------------------------------

        for n in range(
            L + 1
        ):

            actual = shifted_B_term(
                ell,
                s,
                n
            )

            predicted = pochhammer_shifted_term(
                ell,
                s,
                n
            )

            residual = sp.factor(
                actual - predicted
            )

            if residual != 0:
                failures.append(
                    (
                        s,
                        f"term n={n}",
                        residual
                    )
                )

        # ---------------------------------------------------------------
        # Finite sum certificate
        # ---------------------------------------------------------------

        exact = B_exact(
            ell,
            s
        )

        finite = sp.factor(
            finite_hyper_first(
                ell,
                s
            )
            +
            finite_hyper_second(
                ell,
                s
            )
        )

        chu = B_closed(
            ell,
            s
        )

        hyper = hyper_closed(
            ell,
            s
        )

        tests = [
            (
                "finite",
                finite - exact
            ),
            (
                "chu",
                chu - exact
            ),
            (
                "hyper",
                hyper - exact
            )
        ]

        for label, residual in tests:

            residual = sp.factor(
                sp.simplify(
                    residual
                )
            )

            if residual != 0:
                failures.append(
                    (
                        s,
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
    print("KAPPA EXPERIMENT 172")
    print("k=3 B-BRANCH: EXACT CHU-VANDERMONDE COLLAPSE")
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
        27,
        29,
        31
    ]

    total_failures = 0

    # ------------------------------------------------------------------------
    # 1. Symbolic formula
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. SYMBOLIC k=3 FORM")
    print("=" * 78)

    print()
    print(
        "m = s - 3"
    )

    print(
        "L = ell - 2*s - 3"
    )

    print()
    print(
        "First term:"
    )

    print(
        "(-1)^(s-3) C(ell,2s+3)"
        " * 2F1(-L,s-2;2s+4;1)"
    )

    print()
    print(
        "Second term (s>=4):"
    )

    print(
        "(-1)^(s-3) C(ell,2s+3)"
        " * 2F1(-L,s-3;2s+4;1)"
    )

    print()
    print(
        "Chu values:"
    )

    print(
        "2F1(-L,s-2;2s+4;1)"
        " = (s+6)_L/(2s+4)_L"
    )

    print(
        "2F1(-L,s-3;2s+4;1)"
        " = (s+7)_L/(2s+4)_L"
    )

    # ------------------------------------------------------------------------
    # 2. Training
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. TRAINING")
    print("=" * 78)

    for ell in TRAIN:

        failures = check_ell(
            ell
        )

        print()
        print(
            f"ell={ell}: "
            f"{'PASS' if not failures else 'FAIL'}"
        )

        if failures:

            for item in failures[:20]:
                print(
                    " ",
                    item
                )

            total_failures += len(
                failures
            )

        else:

            max_s = (
                ell - 1
            ) // 2

            for s in range(
                3,
                max_s + 1
            ):

                print(
                    f"  s={s}: "
                    f"B_exact={B_exact(ell,s)}, "
                    f"B_closed={B_closed(ell,s)}"
                )

    # ------------------------------------------------------------------------
    # 3. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    holdout_failures = 0

    for ell in HOLDOUT:

        failures = check_ell(
            ell
        )

        print(
            f"ell={ell}: "
            f"status={'PASS' if not failures else 'FAIL'}"
        )

        if failures:

            holdout_failures += 1
            total_failures += len(
                failures
            )

            for item in failures[:10]:
                print(
                    " ",
                    item
                )

    # ------------------------------------------------------------------------
    # 4. Representative shifted terms
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. REPRESENTATIVE PCHAMMER TERMS")
    print("=" * 78)

    for ell, s in [
        (13, 3),
        (15, 4),
        (17, 5),
        (21, 6)
    ]:

        if ell < 2*s + 3:
            continue

        L = ell - 2*s - 3

        print()
        print(
            f"(3,{ell}), s={s}, L={L}"
        )

        for n in range(
            L + 1
        ):

            actual = shifted_B_term(
                ell,
                s,
                n
            )

            predicted = pochhammer_shifted_term(
                ell,
                s,
                n
            )

            print(
                f"  n={n}: "
                f"actual={actual}, "
                f"predicted={predicted}"
            )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "total failures =",
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
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The k=3 B branch is an exact sum of"
        )

        print(
            "two terminating 2F1(1) series."
        )

        print()
        print(
            "Both collapse exactly by"
        )

        print(
            "Chu-Vandermonde."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "combine the k=3 closed B formula with"
        )

        print(
            "the exact C boundary and search for a"
        )

        print(
            "compressed coefficient law in ell and s."
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

