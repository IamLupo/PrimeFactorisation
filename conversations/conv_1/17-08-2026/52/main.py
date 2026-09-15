#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 119R
QUADRATIC TRACE/BRANCH DECOMPOSITION
CORRECT SYMMETRIC DETECTOR CONSTRUCTION
Q(N,S) = C(N,X) + (2S-1) D(N,X)
X = S(S-1)
DISCRIMINANT BRIDGE: (2S-1)^2 = 4X + 1

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================

119R fixes the detector-construction bug in 119.

The detector is constructed directly from

    H_(a,b) = p^a q^b + p^b q^a

using

    H_(a,b) = N^min(a,b) P_|a-b|

and

    P_0 = 2
    P_1 = S
    P_j = S P_(j-1) - N P_(j-2).

This gives an exact symmetric polynomial in N,S.

The experiment then tests the quadratic trace decomposition

    Q(N,S)
      = C(N,X) + (2S-1) D(N,X),

where

    X = S(S-1).

Equivalently,

    Q(N,S) + Q(N,1-S) = 2 C(N,X)

and

    Q(N,S) - Q(N,1-S) = 2(2S-1) D(N,X).

Finally,

    (Q-C)^2 = (4X+1) D^2.

No N-only factoring claim is made.
"""

from __future__ import annotations

import sys
import time

import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

N, S, X, Z = sp.symbols("N S X Z")


# ============================================================================
# BASIC HELPERS
# ============================================================================

def canon(expr):
    return sp.factor(sp.expand(expr))


def poly_S(expr):
    return sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N),
    )


def poly_X(expr):
    return sp.Poly(
        sp.expand(expr),
        X,
        domain=sp.QQ.frac_field(N),
    )


def degree_X(expr):
    p = poly_X(expr)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def degree_N(expr):
    p = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ.frac_field(X),
    )
    if p.is_zero:
        return -sp.oo
    return p.degree()


# ============================================================================
# NEWTON POWER SUMS
# ============================================================================

def power_sums(max_j: int):
    """
    P_j = p^j + q^j

    P_0 = 2
    P_1 = S
    P_j = S P_(j-1) - N P_(j-2)
    """

    P = {
        0: sp.Integer(2),
        1: S,
    }

    for j in range(2, max_j + 1):
        P[j] = canon(
            S * P[j - 1] - N * P[j - 2]
        )

    return P


# ============================================================================
# SYMMETRIC MONOMIAL BLOCK
# ============================================================================

def H(a: int, b: int, P):
    """
    H(a,b) = p^a q^b + p^b q^a

    If a=b:
        2 N^a

    Otherwise:
        N^min(a,b) P_|a-b|
    """

    if a == b:
        return sp.Integer(2) * N**a

    m = min(a, b)
    d = abs(a - b)

    return canon(
        N**m * P[d]
    )


# ============================================================================
# CORRECT DIRECT DETECTOR
# ============================================================================

def detector_F(k: int, ell: int) -> sp.Expr:
    """
    Construct

        F_(k,l) =
          p^k(1+q)^l + q^k(1+p)^l
          - p^l(1+q)^k - q^l(1+p)^k

    directly in the symmetric basis.
    """

    max_j = max(k, ell)

    P = power_sums(max_j)

    F_left = sp.Integer(0)
    F_right = sp.Integer(0)

    for j in range(ell + 1):
        F_left += sp.binomial(ell, j) * H(k, j, P)

    for j in range(k + 1):
        F_right += sp.binomial(k, j) * H(ell, j, P)

    return canon(
        F_left - F_right
    )


# ============================================================================
# EXACT QUOTIENT
# ============================================================================

def quotient_polynomial(k: int, ell: int) -> sp.Expr:
    """
    Q = F/(S+1).
    """

    F = detector_F(k, ell)

    P = poly_S(F)

    divisor = sp.Poly(
        S + 1,
        S,
        domain=sp.QQ.frac_field(N),
    )

    Q, rem = sp.div(
        P,
        divisor,
    )

    rem_expr = canon(rem.as_expr())

    if rem_expr != 0:
        raise ArithmeticError(
            f"Nonzero remainder for ({k},{ell}): "
            f"{rem_expr}"
        )

    return canon(Q.as_expr())


# ============================================================================
# DIRECT QUOTIENT VALIDATION
# ============================================================================

DETECTORS = [
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
    (1, 9),
    (3, 9),
    (5, 9),
    (7, 9),
    (1, 11),
    (3, 11),
    (5, 11),
]


def direct_validation(Q_cache):
    print("\n" + "=" * 78)
    print("1. DIRECT QUOTIENT VALIDATION")
    print("-" * 78)

    failures = 0

    for pair, Q in Q_cache.items():
        k, ell = pair

        F = detector_F(k, ell)

        residual = canon(
            F - (S + 1) * Q
        )

        ok = residual == 0

        print(
            f"Q_({k},{ell}): exact={ok}"
        )

        if not ok:
            failures += 1
            print(
                f"  residual = {residual}"
            )

    print(
        f"\nquotient validation failures = "
        f"{failures}/{len(Q_cache)}"
    )

    if failures:
        raise ArithmeticError(
            "Direct quotient validation failed."
        )

    print("STATUS = PASS")


# ============================================================================
# INVOLUTION TEST
# ============================================================================

def involution_test(Q_cache):
    print("\n" + "=" * 78)
    print("2. INVOLUTION TEST")
    print("-" * 78)

    failures = 0

    for pair, Q in Q_cache.items():
        Q_reflected = canon(
            Q.subs(S, 1 - S)
        )

        diff = canon(
            Q_reflected - Q
        )

        if diff == 0:
            print(
                f"Q_{pair}: invariant=True"
            )
            continue

        # Any nonzero odd part must contain (2S-1).
        P = poly_S(diff)

        divisor = sp.Poly(
            2 * S - 1,
            S,
            domain=sp.QQ.frac_field(N),
        )

        quotient, rem = sp.div(
            P,
            divisor,
        )

        rem_expr = canon(
            rem.as_expr()
        )

        ok = rem_expr == 0

        print(
            f"Q_{pair}: invariant=False "
            f"odd-factor={ok}"
        )

        if not ok:
            failures += 1
            print(
                f"  remainder={rem_expr}"
            )

    print(
        f"\ninvolution failures = "
        f"{failures}/{len(Q_cache)}"
    )

    if failures:
        raise ArithmeticError(
            "Involution structure failed."
        )

    print("STATUS = PASS")


# ============================================================================
# REDUCE POLYNOMIAL MODULO S^2-S-X
# ============================================================================

def reduce_quadratic(expr):
    """
    In the quotient ring

        S^2 = S + X

    every polynomial becomes

        A(N,X) + S B(N,X).

    Return A,B.
    """

    P = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N, X),
    )

    if P.is_zero:
        return sp.Integer(0), sp.Integer(0)

    # Powers represented as
    #
    # S^j = A_j + S B_j.
    #
    # j=0
    A = sp.Integer(1)
    B = sp.Integer(0)

    total_A = sp.Integer(0)
    total_B = sp.Integer(0)

    max_deg = P.degree()

    for j in range(max_deg + 1):

        coeff = canon(
            P.coeff_monomial(S**j)
        )

        total_A += coeff * A
        total_B += coeff * B

        total_A = canon(total_A)
        total_B = canon(total_B)

        # Multiply A + S B by S:
        #
        # S(A + SB)
        # = SA + S^2 B
        # = SA + (S+X)B
        # = XB + S(A+B)

        new_A = canon(
            X * B
        )

        new_B = canon(
            A + B
        )

        A, B = new_A, new_B

    return (
        canon(total_A),
        canon(total_B),
    )


# ============================================================================
# EXACT C + (2S-1)D DECOMPOSITION
# ============================================================================

def branch_decomposition(Q):
    """
    Compute

        Q = C(N,X) + (2S-1)D(N,X).

    """

    Q_reflected = canon(
        Q.subs(S, 1 - S)
    )

    C_raw = canon(
        (Q + Q_reflected) / 2
    )

    odd_raw = canon(
        (Q - Q_reflected) / 2
    )

    P = poly_S(odd_raw)

    divisor = sp.Poly(
        2 * S - 1,
        S,
        domain=sp.QQ.frac_field(N),
    )

    D_S, rem = sp.div(
        P,
        divisor,
    )

    rem_expr = canon(
        rem.as_expr()
    )

    if rem_expr != 0:
        raise ArithmeticError(
            "Odd part is not divisible by (2S-1): "
            f"{rem_expr}"
        )

    D_raw = canon(
        D_S.as_expr()
    )

    # C must be X-only.
    C_A, C_B = reduce_quadratic(
        C_raw
    )

    # D must be X-only.
    D_A, D_B = reduce_quadratic(
        D_raw
    )

    if canon(C_B) != 0:
        raise ArithmeticError(
            "C(N,S) is not X-only.\n"
            f"S coefficient = {C_B}"
        )

    if canon(D_B) != 0:
        raise ArithmeticError(
            "D(N,S) is not X-only.\n"
            f"S coefficient = {D_B}"
        )

    C = canon(C_A)
    D = canon(D_A)

    return C, D


# ============================================================================
# DECOMPOSITION VALIDATION
# ============================================================================

def decomposition_test(Q_cache):
    print("\n" + "=" * 78)
    print("3. EXACT QUADRATIC BRANCH DECOMPOSITION")
    print("-" * 78)

    failures = 0
    result = {}

    for pair, Q in Q_cache.items():

        try:
            C, D = branch_decomposition(Q)
        except Exception as exc:
            failures += 1
            print(
                f"Q_{pair}: ERROR: {exc}"
            )
            continue

        Xs = S * (S - 1)

        reconstructed = canon(
            C.subs(X, Xs)
            + (2 * S - 1)
            * D.subs(X, Xs)
        )

        residual = canon(
            Q - reconstructed
        )

        ok = residual == 0

        print(
            f"Q_{pair}: exact={ok}"
        )

        if not ok:
            failures += 1
            print(
                f"  residual={residual}"
            )
            continue

        result[pair] = (C, D)

        print(
            f"  C(N,X) = {C}"
        )
        print(
            f"  D(N,X) = {D}"
        )

    print(
        f"\ndecomposition failures = "
        f"{failures}/{len(Q_cache)}"
    )

    if failures:
        raise ArithmeticError(
            "Quadratic branch decomposition failed."
        )

    print("STATUS = PASS")

    return result


# ============================================================================
# ELIMINATED RELATION
# ============================================================================

def elimination_test(Q_cache, decomposition):
    print("\n" + "=" * 78)
    print("4. ELIMINATED QUADRATIC RELATION")
    print("-" * 78)

    failures = 0

    for pair, Q in Q_cache.items():

        C, D = decomposition[pair]

        Xs = S * (S - 1)

        C_S = canon(
            C.subs(X, Xs)
        )

        D_S = canon(
            D.subs(X, Xs)
        )

        residual = canon(
            (Q - C_S)**2
            - (4 * Xs + 1)
            * D_S**2
        )

        ok = residual == 0

        print(
            f"Q_{pair}: "
            f"(Q-C)^2-(4X+1)D^2=0 -> {ok}"
        )

        if not ok:
            failures += 1
            print(
                f"  residual={residual}"
            )

    print(
        f"\nelimination failures = "
        f"{failures}/{len(Q_cache)}"
    )

    if failures:
        raise ArithmeticError(
            "Eliminated quadratic relation failed."
        )

    print("STATUS = PASS")


# ============================================================================
# DISCRIMINANT BRIDGE
# ============================================================================

def discriminant_bridge():
    print("\n" + "=" * 78)
    print("5. DISCRIMINANT BRIDGE")
    print("-" * 78)

    residual = canon(
        (2 * S - 1)**2
        - (
            4 * S * (S - 1)
            + 1
        )
    )

    print(
        "(2S-1)^2 - (4S(S-1)+1) =",
        residual,
    )

    if residual != 0:
        raise ArithmeticError(
            "Discriminant bridge failed."
        )

    print("STATUS = PASS")


# ============================================================================
# COMPLEXITY PROFILE
# ============================================================================

def complexity_profile(decomposition):
    print("\n" + "=" * 78)
    print("6. C/D COMPLEXITY PROFILE")
    print("-" * 78)

    for pair, (C, D) in decomposition.items():

        c_poly = sp.Poly(
            sp.expand(C),
            X,
            N,
            domain=sp.QQ,
        )

        d_poly = sp.Poly(
            sp.expand(D),
            X,
            N,
            domain=sp.QQ,
        )

        print(f"Q_{pair}")
        print(
            f"  C: deg_X={degree_X(C)} "
            f"deg_N={degree_N(C)} "
            f"terms={len(c_poly.terms())}"
        )
        print(
            f"  D: deg_X={degree_X(D)} "
            f"deg_N={degree_N(D)} "
            f"terms={len(d_poly.terms())}"
        )


# ============================================================================
# CROSS-DETECTOR STRUCTURE
# ============================================================================

def cross_detector_test(decomposition):
    print("\n" + "=" * 78)
    print("7. CROSS-DETECTOR C/D STRUCTURE")
    print("-" * 78)

    # Test whether the same C or D appears up to simple N-only factors.
    #
    # We do not guess aggressively. We simply report exact polynomial gcds
    # between the C parts and between the D parts.

    pairs = list(decomposition.keys())

    for i in range(len(pairs)):
        p1 = pairs[i]
        C1, D1 = decomposition[p1]

        for j in range(i + 1, len(pairs)):
            p2 = pairs[j]
            C2, D2 = decomposition[p2]

            gC = canon(
                sp.gcd(
                    sp.Poly(
                        sp.expand(C1),
                        X,
                        domain=sp.QQ.frac_field(N),
                    ),
                    sp.Poly(
                        sp.expand(C2),
                        X,
                        domain=sp.QQ.frac_field(N),
                    ),
                ).as_expr()
            )

            gD = canon(
                sp.gcd(
                    sp.Poly(
                        sp.expand(D1),
                        X,
                        domain=sp.QQ.frac_field(N),
                    ),
                    sp.Poly(
                        sp.expand(D2),
                        X,
                        domain=sp.QQ.frac_field(N),
                    ),
                ).as_expr()
            )

            if gC != 1 or gD != 1:
                print(
                    f"{p1} vs {p2}: "
                    f"gcd_C={gC} gcd_D={gD}"
                )


# ============================================================================
# SIMPLE SYMBOLIC SAMPLE
# ============================================================================

def sample_identities(Q_cache, decomposition):
    print("\n" + "=" * 78)
    print("8. SAMPLE IDENTITIES")
    print("-" * 78)

    for pair in [
        (1, 3),
        (1, 5),
        (3, 5),
        (1, 7),
    ]:
        if pair not in Q_cache:
            continue

        Q = Q_cache[pair]
        C, D = decomposition[pair]

        print(f"\nQ_{pair}")
        print(
            "  Q(N,S) =",
            Q
        )
        print(
            "  C(N,X) =",
            C
        )
        print(
            "  D(N,X) =",
            D
        )
        print(
            "  eliminated relation:"
        )
        print(
            "  (Z-C)^2 - (4X+1)D^2 = 0"
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 119R")
    print("QUADRATIC TRACE/BRANCH DECOMPOSITION")
    print("CORRECT SYMMETRIC DETECTOR CONSTRUCTION")
    print("Q(N,S) = C(N,X) + (2S-1)D(N,X)")
    print("X = S(S-1)")
    print("DISCRIMINANT BRIDGE: (2S-1)^2 = 4X + 1")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. BUILD QUOTIENTS
    # ----------------------------------------------------------------------

    print("\n1. DIRECT DETECTOR VALIDATION")
    print("-" * 78)

    Q_cache = {}

    for k, ell in DETECTORS:
        Q = quotient_polynomial(
            k,
            ell,
        )

        Q_cache[(k, ell)] = Q

        print(
            f"Q_({k},{ell}) = {Q}"
        )

    direct_validation(
        Q_cache
    )

    # ----------------------------------------------------------------------
    # 2
    # ----------------------------------------------------------------------

    involution_test(
        Q_cache
    )

    # ----------------------------------------------------------------------
    # 3
    # ----------------------------------------------------------------------

    decomposition = decomposition_test(
        Q_cache
    )

    # ----------------------------------------------------------------------
    # 4
    # ----------------------------------------------------------------------

    elimination_test(
        Q_cache,
        decomposition,
    )

    # ----------------------------------------------------------------------
    # 5
    # ----------------------------------------------------------------------

    discriminant_bridge()

    # ----------------------------------------------------------------------
    # 6
    # ----------------------------------------------------------------------

    complexity_profile(
        decomposition
    )

    # ----------------------------------------------------------------------
    # 7
    # ----------------------------------------------------------------------

    cross_detector_test(
        decomposition
    )

    # ----------------------------------------------------------------------
    # 8
    # ----------------------------------------------------------------------

    sample_identities(
        Q_cache,
        decomposition,
    )

    # ----------------------------------------------------------------------
    # FINAL
    # ----------------------------------------------------------------------

    invariant_count = 0

    for Q in Q_cache.values():
        if canon(
            Q.subs(S, 1 - S) - Q
        ) == 0:
            invariant_count += 1

    print("\n" + "=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        f"invariant detectors = "
        f"{invariant_count}/{len(Q_cache)}"
    )

    print(
        "\nExact target:"
    )

    print(
        "    Q(N,S) = C(N,X) + (2S-1)D(N,X)"
    )

    print(
        "    X = S(S-1)"
    )

    print(
        "    (2S-1)^2 = 4X+1"
    )

    print(
        "\nTherefore:"
    )

    print(
        "    (Q-C)^2 = (4X+1)D^2"
    )

    print(
        "\nThe key question is whether the detector's hidden-trace "
        "dependence really reduces to one quadratic branch variable."
    )

    elapsed = time.perf_counter() - t0

    print(
        f"\ntotal runtime = {elapsed:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 119R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(
            f"\nFATAL: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise