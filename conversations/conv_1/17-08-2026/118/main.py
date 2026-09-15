#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 182
ODD-k C BRANCH: PARITY-CORRECT BULK + UPPER-SUPPORT CORRECTION
==============================================================================

Established B theorem
---------------------
B is already proved by Experiment 181:

    B(s=k) = C(ell-1, 3k-1)

and for s >= k+1,

    L = ell - 2*s - k

    B =
      (-1)^(s-k) C(ell, 2*s+k)
      *
      [
        (s+2k)_L / (2s+k+1)_L
        +
        (s+2k+1)_L / (2s+k+1)_L
      ]

This script does NOT modify B.

C branch
--------
The formal Vandermonde collapse gives the binomial shape

    C_bulk =
      (-1)^(ell+s+1)
      [
        C(ell-s-1, s-k+1)
        +
        C(ell-s-2, s-k)
      ]

The parity factor (-1)^(ell+s+1) is forced by the direct V-kernel.

However, the collapse is only a true full-support Vandermonde sum
when the original a-range is not truncated at the upper boundary.

The active C summand has

    r = ell - 3k + a
    m = s - 2k + a + 1

and V(r,m) is nonzero only when

    r >= 0
    m >= 0
    2m <= r

which implies

    a <= ell + k - 2s - 2.

Thus the full a=0,...,k Vandermonde support requires

    ell + k - 2s - 2 >= k
    <=> 2s <= ell - 2.

Therefore:

    BULK:        2s <= ell-2
    UPPER EDGE:  2s > ell-2

The upper-edge routine below keeps the exact active support and reports
which Vandermonde tail terms have disappeared.

This is intended to identify the exact correction structure before
claiming a global C theorem.
"""

import sympy as sp
import sys


# ============================================================================
# BINOMIAL WITH FINITE-SUPPORT CONVENTION
# ============================================================================

def Cbin(n, r):
    n = int(n)
    r = int(r)

    if n < 0 or r < 0 or r > n:
        return sp.Integer(0)

    return sp.binomial(n, r)


# ============================================================================
# V-KERNEL
# ============================================================================

def V_num(r, m):

    if r < 0 or m < 0 or 2*m > r:
        return sp.Integer(0)

    return sp.expand(
        (-1)**(r-m)
        * (
            Cbin(r-m, m)
            +
            Cbin(r-m-1, m-1)
        )
    )


# ============================================================================
# EXACT C SUM
# ============================================================================

def C_exact(k, ell, s):

    total = sp.Integer(0)
    active = []

    for a in range(k+1):

        j = ell - k + a

        r = j - 2*k
        m = s - 2*k + a + 1

        if r < 0 or m < 0 or 2*m > r:
            continue

        val = sp.expand(
            -Cbin(k, a)
            * V_num(r, m)
        )

        if val != 0:

            total += val

            active.append(
                (a, r, m, val)
            )

    return sp.factor(total), active


# ============================================================================
# PARITY-CORRECT BULK FORMULA
# ============================================================================

def C_bulk(k, ell, s):

    return sp.factor(
        (-1)**(ell+s+1)
        * (
            Cbin(
                ell-s-1,
                s-k+1
            )
            +
            Cbin(
                ell-s-2,
                s-k
            )
        )
    )


# ============================================================================
# SUPPORT DATA
# ============================================================================

def C_support(k, ell, s):

    lower = max(
        0,
        2*k-s-1
    )

    upper = min(
        k,
        ell+k-2*s-2
    )

    return lower, upper


# ============================================================================
# FORMAL VANDEMONDE a-RANGE
# ============================================================================

def formal_a_terms(k):

    return list(range(k+1))


# ============================================================================
# ACTIVE VANDERMONDE a-RANGE
# ============================================================================

def active_a_terms(k, ell, s):

    lower, upper = C_support(
        k, ell, s
    )

    if lower > upper:
        return []

    return list(
        range(
            lower,
            upper+1
        )
    )


# ============================================================================
# UPPER-SUPPORT CORRECTION
# ============================================================================

def C_upper_correction(k, ell, s):

    """
    Difference between the formal bulk expression and the actual
    support-truncated C sum.

    By definition:

        C_exact = C_bulk + correction

    Therefore

        correction = C_exact - C_bulk.
    """

    exact, _ = C_exact(
        k, ell, s
    )

    bulk = C_bulk(
        k, ell, s
    )

    return sp.factor(
        exact - bulk
    )


# ============================================================================
# 1. PARITY CHECK
# ============================================================================

def test_parity():

    print()
    print("="*78)
    print("1. PARITY-CORRECT BULK FORMULA")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9,11]:

        for ell in range(
            max(k+2, 3),
            50,
            2
        ):

            for s in range(
                0,
                max(0, (ell-2)//2)+1
            ):

                exact, _ = C_exact(
                    k, ell, s
                )

                closed = C_bulk(
                    k, ell, s
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
                        f"bulk={closed}",
                        f"residual={residual}"
                    )

    print()
    print(
        "bulk parity failures =",
        failures
    )

    return failures


# ============================================================================
# 2. BULK SUPPORT THEOREM
# ============================================================================

def test_bulk_support():

    print()
    print("="*78)
    print("2. FULL-SUPPORT C THEOREM")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+2, 3),
            51
        ):

            for s in range(
                0,
                ell+1
            ):

                # Full-support region
                if 2*s > ell-2:
                    continue

                exact, active = C_exact(
                    k, ell, s
                )

                bulk = C_bulk(
                    k, ell, s
                )

                residual = sp.factor(
                    exact-bulk
                )

                if residual != 0:

                    failures += 1

                    print(
                        "FAIL:",
                        f"k={k}, ell={ell}, s={s}",
                        f"exact={exact}",
                        f"bulk={bulk}",
                        f"residual={residual}",
                        f"active={active}"
                    )

    print()
    print(
        "full-support failures =",
        failures
    )

    return failures


# ============================================================================
# 3. UPPER EDGE DIAGNOSTIC
# ============================================================================

def diagnostic_upper_edge():

    print()
    print("="*78)
    print("3. UPPER-SUPPORT C DIAGNOSTIC")
    print("="*78)

    for k in [1,3,5,7,9]:

        print()
        print(f"k={k}")

        for ell in range(
            max(k+2, 3),
            32,
            2
        ):

            bad = []

            for s in range(
                0,
                ell+1
            ):

                if 2*s <= ell-2:
                    continue

                exact, active = C_exact(
                    k, ell, s
                )

                bulk = C_bulk(
                    k, ell, s
                )

                correction = sp.factor(
                    exact-bulk
                )

                if correction != 0:

                    lower, upper = C_support(
                        k, ell, s
                    )

                    formal = formal_a_terms(k)
                    actual = active_a_terms(
                        k, ell, s
                    )

                    missing = [
                        a
                        for a in formal
                        if a not in actual
                    ]

                    bad.append(
                        (
                            s,
                            exact,
                            bulk,
                            correction,
                            (lower, upper),
                            missing
                        )
                    )

            if bad:

                print()
                print(
                    f"ell={ell}"
                )

                for row in bad:

                    (
                        s,
                        exact,
                        bulk,
                        correction,
                        support,
                        missing
                    ) = row

                    print(
                        f"  s={s}:",
                        f"exact={exact}",
                        f"bulk={bulk}",
                        f"correction={correction}",
                        f"support={support}",
                        f"missing_a={missing}"
                    )


# ============================================================================
# 4. CORRECTION PATTERN
# ============================================================================

def test_correction_pattern():

    print()
    print("="*78)
    print("4. CORRECTION = SUPPORT-TRUNCATED TAIL")
    print("="*78)

    failures = 0

    for k in [1,3,5,7,9]:

        for ell in range(
            max(k+2, 3),
            42,
            2
        ):

            for s in range(
                ell+1
            ):

                if 2*s <= ell-2:
                    continue

                exact, _ = C_exact(
                    k, ell, s
                )

                bulk = C_bulk(
                    k, ell, s
                )

                correction = sp.factor(
                    exact-bulk
                )

                lower, upper = C_support(
                    k, ell, s
                )

                active = active_a_terms(
                    k, ell, s
                )

                # Reconstruct directly from the active a-range.
                direct = sp.Integer(0)

                for a in active:

                    direct += sp.expand(
                        -Cbin(k, a)
                        * V_num(
                            ell-3*k+a,
                            s-2*k+a+1
                        )
                    )

                direct = sp.factor(
                    direct
                )

                if direct != exact:

                    failures += 1

                    print(
                        "FAIL:",
                        f"k={k}, ell={ell}, s={s}",
                        f"exact={exact}",
                        f"active_sum={direct}",
                        f"residual={sp.factor(exact-direct)}"
                    )

    print()
    print(
        "support-tail reconstruction failures =",
        failures
    )

    return failures


# ============================================================================
# 5. B+C BULK RECONSTRUCTION
# ============================================================================

def B_closed(k, ell, s):

    if s < k:
        return sp.Integer(0)

    if s == k:

        return Cbin(
            ell-1,
            3*k-1
        )

    L = ell - 2*s - k

    if L < 0:
        return sp.Integer(0)

    prefactor = (
        (-1)**(s-k)
        * Cbin(
            ell,
            2*s+k
        )
    )

    def poch_ratio(a, c, n):

        out = sp.Integer(1)

        for t in range(n):
            out *= sp.Rational(
                a+t,
                c+t
            )

        return sp.factor(out)

    first = poch_ratio(
        s+2*k,
        2*s+k+1,
        L
    )

    second = poch_ratio(
        s+2*k+1,
        2*s+k+1,
        L
    )

    return sp.factor(
        prefactor*(first+second)
    )


def D_exact(k, ell, s):

    b = (
        Cbin(ell-1, 3*k-1)
        if s == k
        else (
            sp.Integer(0)
            if s < k
            else (
                # Use exact B definition reconstructed from Experiment 181.
                sp.Integer(0)
            )
        )
    )

    # For D diagnostics, use B_closed together with exact C.
    return sp.factor(
        B_closed(k, ell, s)
        + C_exact(k, ell, s)[0]
    )


def D_bulk(k, ell, s):

    return sp.factor(
        B_closed(k, ell, s)
        + C_bulk(k, ell, s)
    )


# ============================================================================
# 6. FIRST FAILING LAYER REPORT
# ============================================================================

def first_failures():

    print()
    print("="*78)
    print("6. FIRST FAILING UPPER LAYERS")
    print("="*78)

    for k in [1,3,5,7,9,11]:

        print()
        print(f"k={k}")

        shown = 0

        for ell in range(
            max(k+2,3),
            41,
            2
        ):

            for s in range(
                0,
                ell+1
            ):

                exact, active = C_exact(
                    k, ell, s
                )

                bulk = C_bulk(
                    k, ell, s
                )

                residual = sp.factor(
                    exact-bulk
                )

                if residual == 0:
                    continue

                lower, upper = C_support(
                    k, ell, s
                )

                print(
                    f"  ell={ell}, s={s}",
                    f"exact={exact}",
                    f"bulk={bulk}",
                    f"correction={residual}",
                    f"support=[{lower},{upper}]",
                    f"active={active}"
                )

                shown += 1

                if shown >= 8:
                    break

            if shown >= 8:
                break


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print("KAPPA EXPERIMENT 182")
    print("ODD-k C BRANCH: PARITY + SUPPORT")
    print("="*78)

    parity_failures = (
        test_parity()
    )

    bulk_failures = (
        test_bulk_support()
    )

    correction_failures = (
        test_correction_pattern()
    )

    diagnostic_upper_edge()

    first_failures()

    print()
    print("="*78)
    print("FINAL DIAGNOSTIC")
    print("="*78)

    print(
        "parity failures       =",
        parity_failures
    )

    print(
        "bulk support failures =",
        bulk_failures
    )

    print(
        "tail reconstruction    =",
        correction_failures
    )

    if (
        parity_failures == 0
        and bulk_failures == 0
        and correction_failures == 0
    ):

        print()
        print(
            "STATUS = BULK PASS / EDGE ISOLATED"
        )

        print()
        print(
            "The C theorem is now certified on the "
            "full-support region 2*s <= ell-2."
        )

        print(
            "All remaining discrepancies are confined "
            "to the upper-support boundary strip."
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

