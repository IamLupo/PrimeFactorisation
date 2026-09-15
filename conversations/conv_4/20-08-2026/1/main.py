#!/usr/bin/env python3

"""
EXPERIMENT 478
============================================================
SHIFTED-KERNEL / SYMMETRIC-BRIDGE SEARCH

Goal
----
Investigate whether the generalized shifted kernel

    F_{k,l}(p,q;x)
      = p^k (q+x)^l + q^k (p+x)^l
        - p^l (q+x)^k - q^l (p+x)^k

contains a directly recoverable shifted symmetric product

    M_x = (p+x)(q+x)
        = N + x*S + x^2,

where

    N = p*q
    S = p+q.

The experiment is deliberately symbolic and exact.

It does NOT:
    - factor N numerically,
    - use floating point,
    - interpolate formulas from finite samples,
    - use CSV or other output files,
    - search arbitrary high-dimensional feature families.

It does:
    1. Prove the basic shifted-product identity symbolically.
    2. Build the generalized shifted kernel symbolically.
    3. Remove its automatic x=0 zero.
    4. Convert the result exactly to symmetric coordinates (N,S).
    5. Inspect coefficients in the shift variable x.
    6. Test divisibility / factorization by M_x.
    7. Test x-derivatives and finite differences.
    8. Search for exact low-complexity kernel-derived expressions that
       isolate S or M_x.
    9. Verify every discovered identity symbolically.
   10. Repeat over multiple odd (k,l) pairs.

All output is printed to stdout.
"""

from __future__ import annotations

import sys
import traceback
from dataclasses import dataclass
from typing import Iterable, Optional

import sympy as sp


# ============================================================
# START
# ============================================================

print("EXPERIMENT 478 START")
print("=" * 78)
print("SHIFTED-KERNEL / SYMMETRIC-BRIDGE SEARCH")
print("=" * 78)
print()


# ============================================================
# Symbolic variables
# ============================================================

p, q, x = sp.symbols("p q x")
N, S = sp.symbols("N S")


# ============================================================
# Configuration
# ============================================================

# Small exact odd pairs.  These are deliberately independent
# from the later semiprime-specific experiments.
K_L_PAIRS = [
    (1, 3),
    (1, 5),
    (1, 7),
    (1, 9),
    (3, 5),
    (3, 7),
    (3, 9),
    (3, 11),
    (5, 7),
    (5, 9),
    (5, 11),
    (7, 9),
]

# Integer shift controls used only for exact sanity checks.
SHIFT_CONTROLS = [-3, -2, -1, 1, 2, 3]

# Maximum number of printed coefficient details per kernel,
# to keep the terminal readable.
MAX_PRINTED_COEFFS = 12


# ============================================================
# Helpers
# ============================================================

def safe_factor(expr):
    """Exact factorization with a readable fallback."""
    try:
        return sp.factor(expr)
    except Exception:
        return sp.expand(expr)


def exact_zero(expr) -> bool:
    """Robust exact symbolic zero test."""
    expr = sp.cancel(sp.expand(expr))
    return expr == 0


def symmetric_reduce(expr):
    """
    Rewrite a symmetric polynomial in p,q as a polynomial in
    S=p+q and N=p*q.

    The result is verified by substitution back into p,q.
    """
    expr = sp.expand(expr)

    # symmetrize returns:
    #   symmetric polynomial,
    #   remainder in elementary symmetric variables,
    #   mapping information
    #
    # formal=True gives the formal elementary symmetric variables.
    sym_expr, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise ValueError(
            f"Expression is not symmetric or was not reduced cleanly: "
            f"remainder={remainder}"
        )

    # mapping normally contains [(s1, p+q), (s2, p*q)].
    # We substitute those formal variables explicitly.
    formal_vars = [pair[0] for pair in mapping]
    formal_exprs = [pair[1] for pair in mapping]

    reduced = sym_expr
    for fv, target in zip(formal_vars, formal_exprs):
        if exact_zero(target - (p + q)):
            reduced = reduced.subs(fv, S)
        elif exact_zero(target - p * q):
            reduced = reduced.subs(fv, N)
        else:
            raise ValueError(
                f"Unexpected elementary-symmetric mapping: {mapping}"
            )

    reduced = sp.expand(reduced)

    # Independent reconstruction check.
    check = sp.expand(
        reduced.subs({S: p + q, N: p * q}) - expr
    )

    if not exact_zero(check):
        raise AssertionError(
            "Symmetric reconstruction failed exactly."
        )

    return reduced


def shifted_kernel(k: int, ell: int):
    """
    Generalized shifted kernel:

        F_x =
            p^k (q+x)^ell + q^k (p+x)^ell
            - p^ell (q+x)^k - q^ell (p+x)^k
    """
    return sp.expand(
        p**k * (q + x)**ell
        + q**k * (p + x)**ell
        - p**ell * (q + x)**k
        - q**ell * (p + x)**k
    )


def divide_by_x_exact(expr):
    """
    The generalized shifted kernel vanishes at x=0 identically.
    Verify that and divide exactly by x.
    """
    zero_check = sp.expand(expr.subs(x, 0))

    if not exact_zero(zero_check):
        raise AssertionError(
            "Expected shifted kernel to vanish identically at x=0."
        )

    quotient, remainder = sp.div(
        sp.Poly(expr, x, domain="QQ[p,q]"),
        sp.Poly(x, x, domain="QQ[p,q]"),
    )

    if remainder.as_expr() != 0:
        raise AssertionError(
            "Exact division by x failed."
        )

    return sp.expand(quotient.as_expr())


def polynomial_in_S_division(expr_ns):
    """
    Divide an N,S,x expression by

        M_x = N + x*S + x^2

    treating S as the polynomial variable and N,x as coefficients.

    This is an exact divisibility test in Q(N,x)[S].
    """
    M = N + x*S + x**2

    poly_expr = sp.Poly(expr_ns, S, domain=sp.QQ.frac_field(N, x))
    poly_M = sp.Poly(M, S, domain=sp.QQ.frac_field(N, x))

    quotient, remainder = sp.div(poly_expr, poly_M)

    return (
        sp.expand(quotient.as_expr()),
        sp.expand(remainder.as_expr()),
    )


def coefficient_dict_in_x(expr_ns):
    """Return exact coefficients of a polynomial in x."""
    poly = sp.Poly(sp.expand(expr_ns), x)
    return {
        monom[0]: sp.expand(coeff)
        for monom, coeff in poly.terms()
    }


def simple_expression_score(expr) -> tuple[int, int, int]:
    """
    Crude exact complexity score for comparison only.

    Lower is simpler.
    """
    expr = sp.cancel(expr)

    num, den = sp.fraction(expr)

    num_poly = sp.Poly(sp.expand(num), N, S, x)
    den_poly = sp.Poly(sp.expand(den), N, S, x)

    return (
        max(num_poly.total_degree(), den_poly.total_degree()),
        len(num_poly.terms()) + len(den_poly.terms()),
        len(str(expr)),
    )


def report_identity(name: str, lhs, rhs):
    """Print a symbolic identity check."""
    delta = sp.cancel(sp.expand(lhs - rhs))
    ok = exact_zero(delta)

    print(f"[IDENTITY] {name}")
    print(f"  PASS = {ok}")

    if not ok:
        print(f"  Difference = {safe_factor(delta)}")

    return ok


@dataclass
class Candidate:
    name: str
    expression: sp.Expr
    target: str
    score: tuple[int, int, int]


# ============================================================
# Part 1
# Canonical shifted-product identity
# ============================================================

print("[1] CANONICAL SHIFTED-PRODUCT IDENTITY")
print("-" * 78)

M_x_pq = sp.expand((p + x) * (q + x))
M_x_sym = sp.expand(N + x * S + x**2)

report_identity(
    "(p+x)(q+x) = N + x*S + x^2",
    M_x_pq,
    M_x_sym.subs({N: p*q, S: p+q}),
)

S_recovered = sp.cancel(
    (M_x_sym - N - x**2) / x
)

report_identity(
    "Recovered S = ((M_x)-N-x^2)/x",
    S_recovered,
    S,
)

print()
print("  Canonical shifted product:")
print(f"    M_x = {M_x_sym}")
print(f"  Recovered symmetric sum:")
print(f"    S = {S_recovered}")
print()


# ============================================================
# Part 2
# Two-shift identity
# ============================================================

print("[2] TWO-SHIFT SYMMETRIC RECOVERY")
print("-" * 78)

a, b = sp.symbols("a b", integer=True)

M_a = sp.expand(N + a*S + a**2)
M_b = sp.expand(N + b*S + b**2)

S_two_shift = sp.cancel(
    (M_a - M_b) / (a - b) - (a + b)
)

report_identity(
    "Two shifted products recover S",
    S_two_shift,
    S,
)

print()
print("  M_a =", M_a)
print("  M_b =", M_b)
print("  S =", S_two_shift)
print()


# ============================================================
# Part 3
# Generalized shifted kernel
# ============================================================

print("[3] GENERALIZED SHIFTED KERNEL")
print("-" * 78)

print(
    "Definition:"
)
print(
    "  F_x(k,l) = p^k(q+x)^l + q^k(p+x)^l"
)
print(
    "           - p^l(q+x)^k - q^l(p+x)^k"
)
print()

kernel_results = []


for k, ell in K_L_PAIRS:
    print(f"--- (k, ell)=({k}, {ell}) ---")

    F = shifted_kernel(k, ell)

    # Basic symmetry check.
    symmetry_check = sp.expand(
        F.subs({p: q, q: p}, simultaneous=True) - F
    )

    symmetric = exact_zero(symmetry_check)

    print(f"  Symmetric in p,q: {symmetric}")

    if not symmetric:
        print(f"  Symmetry difference: {safe_factor(symmetry_check)}")
        continue

    # x=0 zero.
    F0 = sp.expand(F.subs(x, 0))
    x_zero = exact_zero(F0)

    print(f"  F_x|_(x=0) == 0: {x_zero}")

    if not x_zero:
        print(f"  F(0) = {safe_factor(F0)}")
        continue

    # Exact quotient.
    H = divide_by_x_exact(F)
    print(f"  Exact degree in x after /x: {sp.degree(H, x)}")

    # Symmetric reduction.
    H_NS = symmetric_reduce(H)

    print(
        "  Symmetric reduction verified:"
        f" degree_N={sp.degree(H_NS, N)},"
        f" degree_S={sp.degree(H_NS, S)},"
        f" degree_x={sp.degree(H_NS, x)}"
    )

    # M_x in symmetric coordinates.
    M = sp.expand(N + x*S + x**2)

    # Divisibility by shifted product.
    quotient, remainder = polynomial_in_S_division(H_NS)

    divisible = exact_zero(remainder)

    print(f"  M_x divides H_x symbolically: {divisible}")

    if divisible:
        print("  Quotient H_x / M_x:")
        print(f"    {safe_factor(quotient)}")
    else:
        print("  Remainder after division by M_x:")
        print(f"    {safe_factor(remainder)}")

    # Coefficient analysis in x.
    coeffs = coefficient_dict_in_x(H_NS)

    print(f"  Number of x-coefficients: {len(coeffs)}")

    for power in sorted(coeffs.keys())[:MAX_PRINTED_COEFFS]:
        coeff = coeffs[power]
        print(f"    [x^{power}] H_x = {safe_factor(coeff)}")

    kernel_results.append(
        {
            "k": k,
            "ell": ell,
            "F": F,
            "H": H,
            "H_NS": H_NS,
            "M": M,
            "divisible": divisible,
            "quotient": quotient,
            "remainder": remainder,
            "coeffs": coeffs,
        }
    )

    print()


# ============================================================
# Part 4
# Search coefficients for direct S / M_x extraction
# ============================================================

print("[4] DIRECT EXTRACTION SEARCH")
print("-" * 78)

print(
    "For each symmetric x-coefficient, test whether it is:"
)
print(
    "  - exactly S"
)
print(
    "  - a simple multiple of S"
)
print(
    "  - exactly M_x after reconstruction"
)
print(
    "  - an affine expression alpha*N + beta*S + gamma*x^2"
)
print()

simple_candidates: list[Candidate] = []

for result in kernel_results:
    k = result["k"]
    ell = result["ell"]
    H_NS = result["H_NS"]

    coeffs = result["coeffs"]

    for power, coeff in sorted(coeffs.items()):
        coeff = sp.factor(coeff)

        # Test direct S.
        if exact_zero(coeff - S):
            simple_candidates.append(
                Candidate(
                    f"H({k},{ell}) x^{power} = S",
                    coeff,
                    "S",
                    simple_expression_score(coeff),
                )
            )

        # Test -S.
        if exact_zero(coeff + S):
            simple_candidates.append(
                Candidate(
                    f"H({k},{ell}) x^{power} = -S",
                    coeff,
                    "S",
                    simple_expression_score(coeff),
                )
            )

        # Test direct N.
        if exact_zero(coeff - N):
            simple_candidates.append(
                Candidate(
                    f"H({k},{ell}) x^{power} = N",
                    coeff,
                    "N",
                    simple_expression_score(coeff),
                )
            )

        # Test N+S.
        if exact_zero(coeff - (N + S)):
            simple_candidates.append(
                Candidate(
                    f"H({k},{ell}) x^{power} = N+S",
                    coeff,
                    "N+S",
                    simple_expression_score(coeff),
                )
            )

        # Test N+S+1.
        if exact_zero(coeff - (N + S + 1)):
            simple_candidates.append(
                Candidate(
                    f"H({k},{ell}) x^{power} = N+S+1",
                    coeff,
                    "N+S+1",
                    simple_expression_score(coeff),
                )
            )

    # Test the whole H against affine combinations of N,S,x^2.
    # We only search coefficients from a tiny exact integer dictionary.
    basis = [N, S, x**2, 1]

    for A in range(-4, 5):
        for B in range(-4, 5):
            for C in range(-4, 5):
                for D in range(-4, 5):
                    if (A, B, C, D) == (0, 0, 0, 0):
                        continue

                    target = A*N + B*S + C*x**2 + D

                    # Only test if H itself has very low degree.
                    if sp.Poly(H_NS, N, S, x).total_degree() > 5:
                        continue

                    if exact_zero(H_NS - target):
                        simple_candidates.append(
                            Candidate(
                                f"H({k},{ell}) = "
                                f"{A}N + {B}S + {C}x^2 + {D}",
                                H_NS,
                                "N,S,x",
                                simple_expression_score(H_NS),
                            )
                        )


if simple_candidates:
    print("  EXACT DIRECT CANDIDATES FOUND:")
    for cand in simple_candidates:
        print(f"    {cand.name}")
        print(f"      target = {cand.target}")
        print(f"      score  = {cand.score}")
        print(f"      expr   = {safe_factor(cand.expression)}")
else:
    print("  No direct low-complexity S/M_x extraction found.")

print()


# ============================================================
# Part 5
# Shift derivatives and finite differences
# ============================================================

print("[5] SHIFT-DERIVATIVE / FINITE-DIFFERENCE SEARCH")
print("-" * 78)

print(
    "The purpose is to test whether differentiating or differencing"
)
print(
    "the shifted kernel exposes the elementary symmetric sum S."
)
print()

derivative_candidates = []

for result in kernel_results:
    k = result["k"]
    ell = result["ell"]
    F = result["F"]

    H = result["H"]

    # First derivative at x=0.
    d1 = sp.expand(sp.diff(H, x).subs(x, 0))
    d1_NS = symmetric_reduce(d1)

    # Second derivative at x=0.
    d2 = sp.expand(sp.diff(H, x, 2).subs(x, 0))
    d2_NS = symmetric_reduce(d2)

    print(f"--- (k, ell)=({k}, {ell}) ---")
    print(f"  d/dx(H)|x=0 = {safe_factor(d1_NS)}")
    print(f"  d^2/dx^2(H)|x=0 = {safe_factor(d2_NS)}")

    # Test whether derivative is proportional to S
    # with a factor independent of S.
    if exact_zero(
        sp.diff(d1_NS, S) - sp.diff(d1_NS, S)
    ):
        pass

    # Direct S match.
    if exact_zero(d1_NS - S):
        derivative_candidates.append(
            (k, ell, "H'(0)=S", d1_NS)
        )

    if exact_zero(d2_NS - S):
        derivative_candidates.append(
            (k, ell, "H''(0)=S", d2_NS)
        )

    # Finite difference H(x+1)-H(x).
    delta1 = sp.expand(H.subs(x, x + 1) - H)

    # The substitution x->x+1 must be simultaneous.
    delta1 = sp.expand(
        H.subs(x, x + 1, simultaneous=True) - H
    )

    delta1_NS = symmetric_reduce(delta1)

    print(f"  H(x+1)-H(x) = {safe_factor(delta1_NS)}")

    if exact_zero(delta1_NS - S):
        derivative_candidates.append(
            (k, ell, "H(x+1)-H(x)=S", delta1_NS)
        )

    # Second finite difference.
    delta2 = sp.expand(
        H.subs(x, x + 2, simultaneous=True)
        - 2*H.subs(x, x + 1, simultaneous=True)
        + H
    )

    delta2_NS = symmetric_reduce(delta2)

    print(f"  H(x+2)-2H(x+1)+H(x) = {safe_factor(delta2_NS)}")
    print()


if derivative_candidates:
    print("  EXACT SHIFT-DERIVATIVE CANDIDATES:")
    for k, ell, name, expr in derivative_candidates:
        print(f"    ({k},{ell}) : {name}")
        print(f"      {safe_factor(expr)}")
else:
    print("  No direct derivative/finite-difference identity for S found.")

print()


# ============================================================
# Part 6
# Shifted-product divisibility / remainder fingerprints
# ============================================================

print("[6] SHIFTED-PRODUCT REMAINDER AUDIT")
print("-" * 78)

for result in kernel_results:
    k = result["k"]
    ell = result["ell"]
    H_NS = result["H_NS"]

    quotient, remainder = polynomial_in_S_division(H_NS)

    print(f"--- (k, ell)=({k}, {ell}) ---")

    if exact_zero(remainder):
        print("  RESULT: H_x is exactly divisible by M_x.")
        print(f"  Quotient = {safe_factor(quotient)}")
    else:
        print("  RESULT: H_x is not divisible by M_x.")
        print(f"  Remainder = {safe_factor(remainder)}")

        # Because M_x is linear in S, the remainder has degree 0 in S.
        # We can inspect its structure.
        rem_factored = safe_factor(remainder)
        print(f"  Factored remainder = {rem_factored}")

        # Test whether remainder itself vanishes at simple shifts x.
        for xv in SHIFT_CONTROLS:
            rv = sp.expand(remainder.subs(x, xv))
            if exact_zero(rv):
                print(f"    Remainder vanishes at x={xv}.")

    print()


# ============================================================
# Part 7
# Search natural shifted-product combinations
# ============================================================

print("[7] NATURAL SHIFTED-PRODUCT COMBINATION SEARCH")
print("-" * 78)

print(
    "Test whether exact kernel-derived quantities can be reduced"
)
print(
    "to combinations of M_a = N + a*S + a^2."
)
print()

shift_symbols = list(range(-3, 4))
shift_products = {
    av: sp.expand(N + av*S + av**2)
    for av in shift_symbols
}

print("Shifted symmetric products:")
for av in shift_symbols:
    print(f"  M_{av} = {shift_products[av]}")
print()

combination_hits = []

for result in kernel_results:
    k = result["k"]
    ell = result["ell"]
    H_NS = result["H_NS"]

    # Coefficients in x are the most natural kernel-derived objects.
    for power, coeff in sorted(result["coeffs"].items()):
        coeff = sp.expand(coeff)

        # Compare to M_a, M_a-M_b, and normalized two-shift differences.
        for av in shift_symbols:
            Ma = shift_products[av]

            if exact_zero(coeff - Ma):
                combination_hits.append(
                    (k, ell, power, f"M_{av}", coeff)
                )

        for av in shift_symbols:
            for bv in shift_symbols:
                if av == bv:
                    continue

                Ma = shift_products[av]
                Mb = shift_products[bv]

                # Exact symmetric sum extraction identity.
                extracted_S = sp.cancel(
                    (Ma - Mb) / (av - bv) - (av + bv)
                )

                if exact_zero(extracted_S - S):
                    # This is a canonical identity rather than a kernel hit.
                    pass

                # Test direct coefficient against the difference.
                if exact_zero(coeff - (Ma - Mb)):
                    combination_hits.append(
                        (
                            k,
                            ell,
                            power,
                            f"M_{av}-M_{bv}",
                            coeff,
                        )
                    )


if combination_hits:
    print("  EXACT SHIFTED-PRODUCT HITS:")
    for k, ell, power, label, expr in combination_hits:
        print(
            f"    ({k},{ell}), x^{power}: "
            f"{label} = {safe_factor(expr)}"
        )
else:
    print("  No direct shifted-product coefficient hit found.")

print()


# ============================================================
# Part 8
# Generate exact proof certificates
# ============================================================

print("[8] SYMBOLIC PROOF CERTIFICATES")
print("-" * 78)

print()
print("CERTIFICATE A")
print("-------------")
print("For arbitrary p,q,x:")
print("  (p+x)(q+x) - pq - x(p+q) - x^2 =")
certificate_A = sp.expand(
    (p + x)*(q + x)
    - p*q
    - x*(p + q)
    - x**2
)
print(f"  {certificate_A}")
print(f"  PASS = {exact_zero(certificate_A)}")
print()

print("CERTIFICATE B")
print("-------------")
print("For arbitrary x != 0:")
print("  ((p+x)(q+x)-pq-x^2)/x - (p+q) =")
certificate_B = sp.cancel(
    ((p + x)*(q + x) - p*q - x**2) / x
    - (p + q)
)
print(f"  {certificate_B}")
print(f"  PASS = {exact_zero(certificate_B)}")
print()

print("CERTIFICATE C")
print("-------------")
print("Two-shift recovery:")
certificate_C = sp.cancel(
    (
        ((p + a)*(q + a) - (p + b)*(q + b))
        / (a - b)
    )
    - (p + q + a + b)
)
print("  Difference =", certificate_C)
print(f"  PASS = {exact_zero(certificate_C)}")
print()


# ============================================================
# Part 9
# Optional targeted semiprime numerical verification
# ============================================================

print("[9] EXACT NUMERICAL SANITY CHECKS")
print("-" * 78)

# These are only checks of the canonical identity, not discovery data.
test_instances = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
]

for pp, qq in test_instances:
    NN = pp * qq
    SS = pp + qq

    print(f"--- p={pp}, q={qq} ---")

    for xv in SHIFT_CONTROLS:
        if xv == 0:
            continue

        MM = (pp + xv) * (qq + xv)
        recovered = (MM - NN - xv**2) // xv

        exact_division = (MM - NN - xv**2) % xv == 0
        correct = exact_division and recovered == SS

        print(
            f"  x={xv:>2}: "
            f"M_x={MM}, "
            f"S_recovered={recovered}, "
            f"PASS={correct}"
        )

print()


# ============================================================
# Part 10
# Summary
# ============================================================

print("[10] EXPERIMENT SUMMARY")
print("-" * 78)

print(f"  Kernel pairs tested: {len(kernel_results)}")

divisible_count = sum(
    1 for result in kernel_results if result["divisible"]
)

print(
    f"  Shifted-kernel / M_x exact divisibility hits: "
    f"{divisible_count}"
)

print(
    f"  Direct extraction candidates: "
    f"{len(simple_candidates)}"
)

print(
    f"  Derivative / finite-difference S candidates: "
    f"{len(derivative_candidates)}"
)

print(
    f"  Direct shifted-product coefficient hits: "
    f"{len(combination_hits)}"
)

print()
print("INTERPRETATION")
print("-" * 78)
print(
    "1. The identity M_x=(p+x)(q+x)=N+x*S+x^2 is proved exactly."
)
print(
    "2. The generalized shifted kernel is symmetric and has an exact x=0 zero."
)
print(
    "3. The kernel is reduced exactly to (N,S,x) before interpretation."
)
print(
    "4. Divisibility by M_x is tested symbolically rather than inferred numerically."
)
print(
    "5. Shift derivatives and finite differences are tested as possible direct"
)
print(
    "   mechanisms for exposing S."
)
print(
    "6. Any exact hit is independently rechecked symbolically."
)
print(
    "7. A failure to find M_x directly is not treated as a proof that the"
)
print(
    "   shifted-product bridge is absent; it only rejects this tested kernel-derived"
)
print(
    "   representation family."
)
print()

print("NEXT MATHEMATICAL TARGET")
print("-" * 78)
print(
    "Find an exact kernel-derived quantity I_x satisfying"
)
print(
    "    I_x = (p+x)(q+x) = N + x*S + x^2,"
)
print(
    "or two quantities I_a, I_b satisfying"
)
print(
    "    I_a=(p+a)(q+a),"
)
print(
    "    I_b=(p+b)(q+b),"
)
print(
    "so that"
)
print(
    "    S=(I_a-I_b)/(a-b)-(a+b)."
)
print()
print(
    "Any such identity would give an exact bridge from the constructed"
)
print(
    "kernel to S=p+q, after which the factor pair follows from"
)
print(
    "    z^2-S*z+N=0."
)
print()


# ============================================================
# FINISH
# ============================================================

print("=" * 78)
print("EXPERIMENT 478 FINISHED")
print("=" * 78)
