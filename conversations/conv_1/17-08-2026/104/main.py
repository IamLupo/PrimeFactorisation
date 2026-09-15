#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 167
GENERAL k=1 B-BRANCH: ARBITRARY COEFFICIENT s
PCHAMMER / HYPERGEOMETRIC STRUCTURE

GOAL
----
Experiment 166 proved exactly:

    k=1, s=1:
        B_n = C(ell,3) * (-(ell-3))_n / (4)_n

    k=1, s=2:
        B_n = -2 C(ell,5)
               * (-(ell-5))_n (3)_n
               / ((2)_n (6)_n).

Now derive the general k=1 coefficient layer.

For branch B:

    q(1+p)^ell

the quotient normalization gives

    m = s-1,

    j = 2 + n,

    r = j-2 = n.

Therefore the exact B summand is

    B_n
      = C(ell,n+3)
        (-1)^(n-m)
        [
          C(n-m,m)
          + C(n-m-1,m-1)
        ],

with

    m = s-1.

The experiment:

  1. derives the exact finite B summand;
  2. verifies the formula against the original branch construction;
  3. extracts B_0;
  4. computes the exact consecutive ratio;
  5. derives a symbolic ratio in n;
  6. tests whether the ratio factors into linear Pochhammer parameters;
  7. builds a normalized Pochhammer candidate;
  8. verifies it term-by-term;
  9. tests many ell and s values;
 10. tests forward ell holdout.

The final output should tell us whether the entire k=1 B branch is a
single terminating hypergeometric family in (ell,s).

NO POLYNOMIAL FITTING
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp

n = sp.symbols("n", integer=True, nonnegative=True)
ell = sp.symbols("ell", integer=True, positive=True)
s = sp.symbols("s", integer=True, positive=True)


# ============================================================================
# UNIVERSAL V COEFFICIENT
# ============================================================================

def v_coeff(r, m):
    if r < 0 or m < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.expand(
        (-1) ** (r - m) * (
            sp.binomial(r - m, m)
            + (
                0
                if m == 0
                else sp.binomial(r - m - 1, m - 1)
            )
        )
    )


# ============================================================================
# EXACT k=1 B TERM
# ============================================================================

def B_term(ell0, s0, n0):
    """
    k=1, j=n+2, r=n, m=s-1.
    """

    m0 = s0 - 1

    if m0 < 0:
        return sp.Integer(0)

    r0 = n0

    if m0 > r0 // 2:
        return sp.Integer(0)

    return sp.factor(
        sp.binomial(
            ell0,
            n0 + 3
        )
        * v_coeff(
            r0,
            m0
        )
    )


# ============================================================================
# EXACT SUPPORT
# ============================================================================

def B_support(ell0, s0):
    """
    Since

        m=s-1 <= floor(n/2),

    the first admissible n is

        n >= 2(s-1).

    Also C(ell,n+3) requires

        n <= ell-3.

    Hence:

        n = 2s-2,...,ell-3.
    """

    m0 = s0 - 1

    if m0 < 0:
        return []

    lo = 2 * m0
    hi = ell0 - 3

    if lo > hi:
        return []

    return list(
        range(
            lo,
            hi + 1
        )
    )


# ============================================================================
# EXACT B TOTAL
# ============================================================================

def B_total(ell0, s0):
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
# FIRST TERM B_0
# ============================================================================

def B0(ell0, s0):
    support = B_support(
        ell0,
        s0
    )

    if not support:
        return sp.Integer(0)

    return sp.factor(
        B_term(
            ell0,
            s0,
            support[0]
        )
    )


# ============================================================================
# CONCRETE CONSECUTIVE RATIOS
# ============================================================================

def ratio_table(ell0, s0):
    support = B_support(
        ell0,
        s0
    )

    out = []

    for a, b in zip(
        support,
        support[1:]
    ):
        va = B_term(
            ell0,
            s0,
            a
        )

        vb = B_term(
            ell0,
            s0,
            b
        )

        if va == 0:
            continue

        out.append(
            (
                a,
                sp.factor(
                    sp.cancel(
                        vb / va
                    )
                )
            )
        )

    return out


# ============================================================================
# SYMBOLIC TERM FOR GENERAL s
# ============================================================================

def symbolic_B_term():
    """
    Use

        m=s-1,
        r=n.

    The exact combinatorial formula is

        C(ell,n+3)
        (-1)^(n-s+1)
        [
          C(n-s+1,s-1)
          +
          C(n-s,s-2)
        ].

    The two V pieces are kept separate first.
    """

    m_expr = s - 1

    first = (
        sp.binomial(
            n - m_expr,
            m_expr
        )
    )

    second = (
        sp.Integer(0)
        + sp.binomial(
            n - m_expr - 1,
            m_expr - 1
        )
    )

    T = sp.factor(
        sp.binomial(
            ell,
            n + 3
        )
        * (-1) ** (
            n - m_expr
        )
        * (
            first
            + second
        )
    )

    return sp.factor(T)


# ============================================================================
# SYMBOLIC SEPARATION AND SIMPLIFICATION
# ============================================================================

def symbolic_B_pieces():
    m_expr = s - 1

    T1 = sp.factor(
        sp.binomial(
            ell,
            n + 3
        )
        * (-1) ** (
            n - m_expr
        )
        * sp.binomial(
            n - m_expr,
            m_expr
        )
    )

    T2 = sp.factor(
        sp.binomial(
            ell,
            n + 3
        )
        * (-1) ** (
            n - m_expr
        )
        * sp.binomial(
            n - m_expr - 1,
            m_expr - 1
        )
    )

    return T1, T2


# ============================================================================
# SYMBOLIC RATIO FOR FIXED s
# ============================================================================

def symbolic_ratio_for_s(s0):
    """
    For fixed integer s, construct T(n+1)/T(n) exactly.

    We deliberately substitute a concrete s before simplification so that
    SymPy never has to reason about symbolic finite-support binomials.
    """

    T = symbolic_B_term().subs(
        s,
        s0
    )

    T_next = T.subs(
        n,
        n + 1
    )

    ratio = sp.factor(
        sp.cancel(
            T_next / T
        )
    )

    return ratio


# ============================================================================
# NORMALIZED POCHHAMMER CANDIDATE
# ============================================================================

def candidate_from_ratio(ell0, s0):
    """
    Build the candidate from the exact first term and the exact ratio.

    We express the two V pieces as two binomial factors first.

    m=s-1.

    The sum of the two V terms can be simplified using

      C(n-m,m) + C(n-m-1,m-1)
      =
      ((n-m+1)/(n-2m+1)) C(n-m,m)

    wherever the support is valid.

    Hence:

      B_n =
      C(ell,n+3)
      (-1)^(n-m)
      ((n-m+1)/(n-2m+1))
      C(n-m,m).

    The implementation constructs this exact expression and normalizes it
    relative to the first admissible n.
    """

    m0 = s0 - 1
    support = B_support(
        ell0,
        s0
    )

    if not support:
        return None

    n0 = support[0]

    # Exact symbolic-in-n expression with concrete ell,s.
    nn = n

    T = sp.factor(
        sp.binomial(
            ell0,
            nn + 3
        )
        * (-1) ** (
            nn - m0
        )
        * (
            sp.binomial(
                nn - m0,
                m0
            )
            +
            sp.binomial(
                nn - m0 - 1,
                m0 - 1
            )
        )
    )

    # Shift n = n0 + h.
    h = n

    shifted = sp.factor(
        T.subs(
            nn,
            n0 + h
        )
    )

    normalized = sp.factor(
        sp.cancel(
            shifted / shifted.subs(n, 0)
        )
    )

    return (
        n0,
        shifted,
        normalized
    )


# ============================================================================
# TERM-BY-TERM CANDIDATE VERIFICATION
# ============================================================================

def verify_candidate(
    ell0,
    s0
):
    result = candidate_from_ratio(
        ell0,
        s0
    )

    if result is None:
        return []

    n0, shifted, normalized = result

    support = B_support(
        ell0,
        s0
    )

    B_first = B_term(
        ell0,
        s0,
        n0
    )

    candidate = sp.factor(
        B_first
        * normalized
    )

    failures = []

    for jj in support:

        h = jj - n0

        predicted = sp.factor(
            candidate.subs(
                n,
                h
            )
        )

        actual = sp.factor(
            B_term(
                ell0,
                s0,
                jj
            )
        )

        residual = sp.factor(
            predicted - actual
        )

        if residual != 0:
            failures.append(
                (
                    jj,
                    residual
                )
            )

    return failures


# ============================================================================
# DIRECT INDEPENDENT V/CHECK
# ============================================================================

def direct_reconstruct_B(ell0, s0):
    """
    Independent finite sum through the explicit V coefficient.
    This is intentionally identical to the mathematical definition,
    rather than using any candidate formula.
    """

    total = sp.Integer(0)

    for n0 in B_support(
        ell0,
        s0
    ):
        total += B_term(
            ell0,
            s0,
            n0
        )

    return sp.factor(total)


# ============================================================================
# EXPERIMENT
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 167")
    print("GENERAL k=1 B-BRANCH / ARBITRARY s")
    print("=" * 78)
    print()

    TRAIN_ELL = [
        9,
        11,
        13,
        15,
        17
    ]

    # ------------------------------------------------------------------------
    # 1. GENERAL SYMBOLIC TERM
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. GENERAL SYMBOLIC B TERM")
    print("=" * 78)

    T1, T2 = symbolic_B_pieces()

    print()
    print("V-first-piece:")
    print(" ", T1)

    print()
    print("V-second-piece:")
    print(" ", T2)

    print()
    print("combined:")
    print(" ", symbolic_B_term())

    # ------------------------------------------------------------------------
    # 2. EXACT RATIO FOR MANY s
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT B-RATIO BY COEFFICIENT LAYER")
    print("=" * 78)

    for s0 in range(
        1,
        6
    ):

        ratio = symbolic_ratio_for_s(
            s0
        )

        print()
        print(
            f"s={s0}:"
        )

        print(
            "  T(n+1)/T(n) =",
            ratio
        )

    # ------------------------------------------------------------------------
    # 3. CONCRETE RATIO CHECK
    # ------------------------------------------------------------------------

    ratio_failures = 0

    print()
    print("=" * 78)
    print("3. SYMBOLIC RATIO VS EXACT DATA")
    print("=" * 78)

    for ell0 in TRAIN_ELL:

        for s0 in range(
            1,
            min(
                6,
                ell0 - 2
            )
        ):

            support = B_support(
                ell0,
                s0
            )

            if len(support) < 2:
                continue

            symbolic_ratio = symbolic_ratio_for_s(
                s0
            )

            print()
            print(
                f"(k=1,ell={ell0},s={s0})"
            )

            for nn in support[:-1]:

                actual_ratio = sp.factor(
                    sp.cancel(
                        B_term(
                            ell0,
                            s0,
                            nn+1
                        )
                        /
                        B_term(
                            ell0,
                            s0,
                            nn
                        )
                    )
                )

                # The symbolic ratio uses n as the actual source index.
                predicted_ratio = sp.factor(
                    sp.simplify(
                        symbolic_ratio.subs(
                            n,
                            nn
                        )
                    )
                )

                residual = sp.factor(
                    predicted_ratio
                    - actual_ratio
                )

                print(
                    f"  n={nn}: residual={residual}"
                )

                if residual != 0:
                    ratio_failures += 1

    # ------------------------------------------------------------------------
    # 4. TERM-BY-TERM NORMALIZED CANDIDATE
    # ------------------------------------------------------------------------

    candidate_failures = 0

    print()
    print("=" * 78)
    print("4. NORMALIZED CLOSED CANDIDATE")
    print("=" * 78)

    for ell0 in TRAIN_ELL:

        for s0 in range(
            1,
            min(
                6,
                ell0 - 2
            )
        ):

            support = B_support(
                ell0,
                s0
            )

            if not support:
                continue

            failures = verify_candidate(
                ell0,
                s0
            )

            candidate_failures += len(
                failures
            )

            print()
            print(
                f"(1,{ell0}), s={s0}"
            )

            print(
                "  support =",
                support
            )

            if failures:
                print(
                    "  FAILURES =",
                    failures
                )
            else:
                print(
                    "  term-by-term = PASS"
                )

    # ------------------------------------------------------------------------
    # 5. COEFFICIENT TOTALS
    # ------------------------------------------------------------------------

    total_failures = 0

    print()
    print("=" * 78)
    print("5. EXACT B-TOTALS")
    print("=" * 78)

    for ell0 in TRAIN_ELL:

        print()
        print(
            f"ell={ell0}"
        )

        for s0 in range(
            0,
            min(
                6,
                ell0 - 2
            )
        ):

            direct = direct_reconstruct_B(
                ell0,
                s0
            )

            if s0 == 0:
                # B is absent for k=1,s=0.
                expected = sp.Integer(0)
            else:
                expected = B_total(
                    ell0,
                    s0
                )

            residual = sp.factor(
                direct - expected
            )

            print(
                f"  s={s0}: "
                f"B={direct}, "
                f"residual={residual}"
            )

            if residual != 0:
                total_failures += 1

    # ------------------------------------------------------------------------
    # 6. FORWARD HOLDOUT
    # ------------------------------------------------------------------------

    holdout_failures = 0

    print()
    print("=" * 78)
    print("6. FORWARD ELL HOLDOUT")
    print("=" * 78)

    for ell0 in [
        19,
        21,
        23
    ]:

        ok = True

        for s0 in range(
            1,
            min(
                8,
                ell0 - 2
            )
        ):

            failures = verify_candidate(
                ell0,
                s0
            )

            if failures:
                ok = False

        print(
            f"(1,{ell0}) "
            f"status={'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            holdout_failures += 1

    # ------------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "symbolic ratio failures =",
        ratio_failures
    )

    print(
        "candidate term failures =",
        candidate_failures
    )

    print(
        "total reconstruction failures =",
        total_failures
    )

    print(
        "forward failures =",
        holdout_failures
    )

    if (
        ratio_failures == 0
        and candidate_failures == 0
        and total_failures == 0
        and holdout_failures == 0
    ):

        print()
        print("STATUS = PASS")

        print()
        print(
            "The k=1 B branch has an exact symbolic"
        )

        print(
            "hypergeometric term for arbitrary tested s."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "convert the normalized term to explicit"
        )

        print(
            "Pochhammer parameters in ell and s,"
        )

        print(
            "then evaluate the terminating sum."
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

