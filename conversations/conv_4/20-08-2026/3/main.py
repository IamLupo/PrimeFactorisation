#!/usr/bin/env python3

"""
EXPERIMENT 480
==============================================================================
SHIFTED-KERNEL COEFFICIENT TRANSFER

Goal
----

Experiment 479 established the exact theorem

    [x^(ell-1)] H_(k,ell)(p,q;x) = p^k + q^k,

where

    H_(k,ell)(p,q;x)
      =
      [p^k(q+x)^ell + q^k(p+x)^ell
       - p^ell(q+x)^k - q^ell(p+x)^k] / x.

For k=1,

    [x^(ell-1)] H_(1,ell) = p+q = S.

The unresolved question is no longer whether the shifted kernel contains S.
It does.

The new question is:

    Can the top shifted-kernel coefficient be reconstructed from the
    already-existing unshifted kernel / N,S homogeneous structure?

The original KAPPA kernel is exactly the x=1 specialization:

    F_(k,ell)(p,q) = F_(k,ell)(p,q;1).

This experiment therefore studies the complete x-coefficient family

    C_r(k,ell) = [x^r] H_(k,ell)(p,q;x)

and asks whether the coefficient C_(ell-1)=p^k+q^k can be recovered from
natural structural transforms of the x=1 specialization.

Main tasks
----------

1. Derive the exact closed coefficient formula for [x^r] H_(k,ell).
2. Verify the formula symbolically against direct expansion.
3. Prove the top coefficient theorem again from the coefficient formula.
4. Express the entire coefficient family in symmetric variables N=pq,
   S=p+q.
5. Specialize x=1 and compare the resulting coefficients with the original
   unshifted KAPPA kernel F_(k,ell).
6. Search for exact linear differential / finite-difference reconstruction
   operators which recover the top coefficient from the x=1 specialization.
7. Search structured binomial transforms and reversal transforms.
8. Explicitly test whether the top coefficient can be generated from the
   lower x-coefficients using universal coefficients depending only on ell.
9. Investigate the normalized remainder factors from Experiment 478 after
   substituting y=x^2/N.
10. Produce exact proof certificates for every successful identity.

No floating point.
No CSV.
No external files.
No numerical interpolation as evidence.
No factorization of N.

All output is printed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import sympy as sp


# ============================================================================
# START
# ============================================================================

print("EXPERIMENT 480 START")
print("=" * 78)
print("SHIFTED-KERNEL COEFFICIENT TRANSFER")
print("=" * 78)
print()


# ============================================================================
# Symbols
# ============================================================================

p, q, x = sp.symbols("p q x")
N, S = sp.symbols("N S")
y = sp.symbols("y")


# ============================================================================
# Configuration
# ============================================================================

PAIRS = [
    (1, 2),
    (1, 3),
    (1, 4),
    (1, 5),
    (1, 6),
    (1, 7),
    (1, 8),
    (1, 9),
    (1, 10),
    (1, 11),
    (2, 3),
    (2, 5),
    (2, 7),
    (3, 4),
    (3, 5),
    (3, 7),
    (3, 9),
    (3, 11),
    (4, 5),
    (5, 7),
    (5, 9),
    (5, 11),
    (6, 7),
    (7, 9),
    (7, 11),
]


# ============================================================================
# Helpers
# ============================================================================

def exact_zero(expr) -> bool:
    return sp.cancel(sp.expand(expr)) == 0


def shifted_kernel(k: int, ell: int):
    return sp.expand(
        p**k * (q + x)**ell
        + q**k * (p + x)**ell
        - p**ell * (q + x)**k
        - q**ell * (p + x)**k
    )


def shifted_H(k: int, ell: int):
    """
    Exact H = F/x.
    """
    F = shifted_kernel(k, ell)

    if not exact_zero(F.subs(x, 0)):
        raise AssertionError(
            f"F({k},{ell};0) != 0"
        )

    quotient, remainder = sp.div(
        sp.Poly(F, x, domain=sp.QQ.frac_field(p, q)),
        sp.Poly(x, x, domain=sp.QQ.frac_field(p, q)),
    )

    if remainder.as_expr() != 0:
        raise AssertionError(
            f"Division by x failed for ({k},{ell})"
        )

    return sp.expand(quotient.as_expr())


def symmetric_reduce(expr):
    """
    Reduce a symmetric p,q polynomial to N=pq and S=p+q.
    """
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Symmetrization remainder is not zero: {remainder}"
        )

    result = reduced

    for formal_var, concrete in mapping:
        if exact_zero(concrete - (p + q)):
            result = result.subs(formal_var, S)
        elif exact_zero(concrete - p*q):
            result = result.subs(formal_var, N)
        else:
            raise AssertionError(
                f"Unexpected elementary symmetric mapping: {mapping}"
            )

    result = sp.expand(result)

    check = sp.expand(
        result.subs({
            S: p + q,
            N: p*q,
        }) - expr
    )

    if not exact_zero(check):
        raise AssertionError(
            "Symmetric reduction reconstruction failed."
        )

    return result


def coefficient_formula(k: int, ell: int, r: int):
    """
    Exact coefficient formula.

    Since

        H = F/x,

    [x^r] H = [x^(r+1)] F.

    For the first two terms:
        [x^(r+1)] p^k(q+x)^ell
            = p^k * binom(ell,r+1) * q^(ell-r-1)

        [x^(r+1)] q^k(p+x)^ell
            = q^k * binom(ell,r+1) * p^(ell-r-1)

    For the two negative terms:
        -p^ell(q+x)^k
        -q^ell(p+x)^k

    they contribute only when r+1 <= k.

    Therefore

      C_r =
        binom(ell,r+1)
        [p^k q^(ell-r-1) + q^k p^(ell-r-1)]
        -
        1_{r+1 <= k}
        binom(k,r+1)
        [p^ell q^(k-r-1) + q^ell p^(k-r-1)].

    This formula is exact and does not rely on expansion.
    """
    j = r + 1

    positive = (
        sp.binomial(ell, j)
        * (
            p**k * q**(ell-j)
            + q**k * p**(ell-j)
        )
    )

    if j <= k:
        negative = (
            sp.binomial(k, j)
            * (
                p**ell * q**(k-j)
                + q**ell * p**(k-j)
            )
        )
    else:
        negative = sp.Integer(0)

    return sp.expand(positive - negative)


def coefficient_from_H(k: int, ell: int, r: int):
    H = shifted_H(k, ell)

    return sp.expand(
        sp.Poly(H, x).coeff_monomial(x**r)
    )


def homogeneous_degree(expr, variables):
    try:
        return sp.Poly(expr, *variables).total_degree()
    except Exception:
        return None


# ============================================================================
# [1] Experiment theorem
# ============================================================================

print("[1] EXACT COEFFICIENT FORMULA")
print("-" * 78)

print(
    "For r >= 0 define C_r = [x^r] H_(k,l)."
)
print()
print(
    "Since H = F/x:"
)
print(
    "  C_r = [x^(r+1)] F."
)
print()
print(
    "The exact coefficient formula is:"
)
print(
    "  C_r = binom(l,r+1)"
)
print(
    "        * (p^k q^(l-r-1) + q^k p^(l-r-1))"
)
print(
    "        - 1_(r+1<=k) binom(k,r+1)"
)
print(
    "        * (p^l q^(k-r-1) + q^l p^(k-r-1))"
)
print()


# ============================================================================
# [2] Verify coefficient formula against full symbolic expansion
# ============================================================================

print("[2] FULL COEFFICIENT-FORMULA AUDIT")
print("-" * 78)

coefficient_failures = []

for k, ell in PAIRS:
    H = shifted_H(k, ell)
    H_poly = sp.Poly(H, x)

    case_pass = True

    for r in range(ell):
        actual = sp.expand(
            H_poly.coeff_monomial(x**r)
        )
        predicted = coefficient_formula(k, ell, r)

        if not exact_zero(actual - predicted):
            case_pass = False
            coefficient_failures.append(
                (k, ell, r)
            )

    print(
        f"  (k,l)=({k:2d},{ell:2d}) "
        f"PASS={case_pass}"
    )

print(
    f"  TOTAL COEFFICIENT FAILURES = "
    f"{len(coefficient_failures)}"
)

if coefficient_failures:
    print(
        f"  Failed locations = {coefficient_failures}"
    )

print()


# ============================================================================
# [3] Top coefficient theorem from the general coefficient formula
# ============================================================================

print("[3] TOP COEFFICIENT THEOREM FROM GENERAL FORMULA")
print("-" * 78)

top_failures = []

for k, ell in PAIRS:
    top = coefficient_formula(k, ell, ell - 1)
    expected = p**k + q**k

    ok = exact_zero(top - expected)

    print(
        f"  (k,l)=({k:2d},{ell:2d}) "
        f"top={sp.factor(top)} "
        f"PASS={ok}"
    )

    if not ok:
        top_failures.append((k, ell))

print(
    f"  TOP THEOREM FAILURES = {len(top_failures)}"
)
print()


# ============================================================================
# [4] Exact x-support / coefficient boundary
# ============================================================================

print("[4] X-SUPPORT STRUCTURE")
print("-" * 78)

support_failures = []

for k, ell in PAIRS:
    H = shifted_H(k, ell)
    degree = sp.degree(H, x)

    expected_degree = ell - 1

    ok_degree = degree == expected_degree

    top = sp.Poly(H, x).coeff_monomial(
        x**(ell - 1)
    )

    above = sp.Poly(H, x).coeff_monomial(
        x**ell
    )

    ok_above = exact_zero(above)

    ok = ok_degree and ok_above

    print(
        f"  (k,l)=({k:2d},{ell:2d}) "
        f"degree={degree} expected={expected_degree} "
        f"PASS={ok}"
    )

    if not ok:
        support_failures.append((k, ell))

print(
    f"  SUPPORT FAILURES = {len(support_failures)}"
)
print()


# ============================================================================
# [5] Symmetric reduction of every coefficient
# ============================================================================

print("[5] SYMMETRIC COEFFICIENT REPRESENTATION")
print("-" * 78)

for k, ell in PAIRS:
    print(f"--- (k,l)=({k},{ell}) ---")

    for r in range(ell):
        coeff = coefficient_formula(k, ell, r)
        coeff_NS = symmetric_reduce(coeff)

        # Only print selected structurally important coefficients.
        if r in (0, 1, ell - 2, ell - 1):
            print(
                f"  [x^{r}] = "
                f"{sp.factor(coeff_NS)}"
            )

    print()


# ============================================================================
# [6] Top coefficient = Newton power sum
# ============================================================================

print("[6] NEWTON POWER-SUM IDENTIFICATION")
print("-" * 78)

print(
    "Define P_k = p^k + q^k."
)
print(
    "Then C_(ell-1) = P_k."
)
print()

for k in range(1, 9):
    value = coefficient_formula(k, k + 1, k)

    expected = p**k + q**k

    ok = exact_zero(value - expected)

    print(
        f"  k={k}: "
        f"C_k={sp.factor(value)} "
        f"PASS={ok}"
    )

print()


# ============================================================================
# [7] Derive Newton recurrence entirely from coefficient data
# ============================================================================

print("[7] POWER-SUM RECURRENCE FROM COEFFICIENT DATA")
print("-" * 78)

P = {0: sp.Integer(2)}

# P_1 is the top coefficient for any ell > 1.
P[1] = coefficient_formula(1, 2, 1)

recurrence_failures = []

for k in range(2, 13):
    P[k] = coefficient_formula(k, k + 1, k)

    lhs = symmetric_reduce(P[k])

    rhs = sp.expand(
        S * symmetric_reduce(P[k - 1])
        - N * symmetric_reduce(P[k - 2])
    )

    ok = exact_zero(lhs - rhs)

    print(
        f"  k={k:2d}: "
        f"P_k={sp.factor(lhs)} "
        f"PASS={ok}"
    )

    if not ok:
        recurrence_failures.append(k)

print(
    f"  RECURRENCE FAILURES = {recurrence_failures}"
)
print()


# ============================================================================
# [8] x=1 specialization = original unshifted kernel
# ============================================================================

print("[8] X=1 SPECIALIZATION / ORIGINAL KAPPA KERNEL")
print("-" * 78)

print(
    "The generalized shifted kernel at x=1 is:"
)
print(
    "  F_shifted(p,q;1) ="
)
print(
    "    p^k(1+q)^l + q^k(1+p)^l"
)
print(
    "    - p^l(1+q)^k - q^l(1+p)^k"
)
print(
)
print(
    "This is the original unshifted KAPPA kernel."
)
print()

specialization_failures = []

for k, ell in PAIRS:
    shifted_at_1 = sp.expand(
        shifted_kernel(k, ell).subs(x, 1)
    )

    original = sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )

    ok = exact_zero(shifted_at_1 - original)

    print(
        f"  (k,l)=({k:2d},{ell:2d}) "
        f"PASS={ok}"
    )

    if not ok:
        specialization_failures.append((k, ell))

print(
    f"  SPECIALIZATION FAILURES = "
    f"{len(specialization_failures)}"
)
print()


# ============================================================================
# [9] Can the x-top coefficient be reconstructed from x=1 alone?
# ============================================================================

print("[9] CAN TOP COEFFICIENT BE RECOVERED FROM x=1?")
print("-" * 78)

print(
    "A single evaluation at x=1 cannot generally determine every x-coefficient."
)
print(
    "We therefore test structured transforms rather than arbitrary interpolation."
)
print()

# Candidate operators:
#
#   finite differences in ell
#   factorial-normalized combinations
#   endpoint coefficient extraction from the x=1 polynomial in N,S
#
# The experiment deliberately does NOT fit arbitrary coefficients.

transform_results = []

for k, ell in PAIRS:
    H = shifted_H(k, ell)

    top = sp.expand(
        sp.Poly(H, x).coeff_monomial(x**(ell - 1))
    )

    # x=1 value of H.
    H1 = sp.expand(H.subs(x, 1))

    # Derivative relation:
    top_from_derivative = sp.expand(
        sp.diff(H, x, ell - 1).subs(x, 0)
        / sp.factorial(ell - 1)
    )

    # Exact check.
    deriv_ok = exact_zero(
        top_from_derivative - top
    )

    # Re-expand the x=1 value in N,S.
    H1_NS = symmetric_reduce(H1)

    # Record complexity only; do not interpret it as a discovery.
    transform_results.append(
        {
            "k": k,
            "ell": ell,
            "top": top,
            "H1_NS": H1_NS,
            "derivative_ok": deriv_ok,
        }
    )

    print(
        f"  (k,l)=({k},{ell}) "
        f"derivative extraction PASS={deriv_ok}"
    )

print()


# ============================================================================
# [10] Structured x=1 coefficient transforms
# ============================================================================

print("[10] STRUCTURED x=1 TRANSFORM SEARCH")
print("-" * 78)

print(
    "For H(x)=sum_r C_r x^r, H(1)=sum_r C_r."
)
print(
    "We test exact finite-difference/binomial identities that could isolate"
)
print(
    "the terminal coefficient without arbitrary interpolation."
)
print()

structured_hits = []

for k, ell in PAIRS:
    H = shifted_H(k, ell)

    coefficients = [
        sp.expand(
            sp.Poly(H, x).coeff_monomial(x**r)
        )
        for r in range(ell)
    ]

    top = coefficients[-1]

    # Test alternating binomial sums corresponding to endpoint extraction.
    # These are exact identities for finite polynomial coefficient vectors,
    # but only retained as structural candidates if they reduce to top.
    for m in range(0, min(ell - 1, 5) + 1):

        # Finite-difference functional
        functional = sp.expand(
            sum(
                (-1)**(m-r)
                * sp.binomial(m, r)
                * coefficients[r]
                for r in range(m + 1)
            )
        )

        if exact_zero(functional - top):
            structured_hits.append(
                (k, ell, "binomial-prefix", m)
            )

    # Reverse binomial functional.
    for m in range(0, min(ell - 1, 5) + 1):

        functional = sp.expand(
            sum(
                (-1)**(r)
                * sp.binomial(m, r)
                * coefficients[ell - 1 - r]
                for r in range(m + 1)
            )
        )

        if exact_zero(functional - top):
            structured_hits.append(
                (k, ell, "reverse-binomial", m)
            )

print(
    f"  STRUCTURED TOP-COEFFICIENT HITS = "
    f"{len(structured_hits)}"
)

if structured_hits:
    for hit in structured_hits:
        print(
            f"    {hit}"
        )
else:
    print(
        "  No low-order structured x=1 transform exactly reconstructs the top coefficient."
    )

print()


# ============================================================================
# [11] Remainder factors normalized by y=x^2/N
# ============================================================================

print("[11] NORMALIZED M_x REMAINDER STRUCTURE")
print("-" * 78)

print(
    "Experiment 478 produced remainders after division by"
)
print(
    "  M_x = N + S*x + x^2."
)
print()
print(
    "These remainders contain factors depending only on x^2/N."
)
print(
    "We now normalize them using y=x^2/N."
)
print()

for k, ell in PAIRS[:12]:

    # Build H in symmetric coordinates.
    H = shifted_H(k, ell)
    H_NS = sp.expand(
        symmetric_reduce(H)
    )

    M = N + S*x + x**2

    # Polynomial division in S.
    poly_H = sp.Poly(
        H_NS,
        S,
        domain=sp.QQ.frac_field(N, x),
    )

    poly_M = sp.Poly(
        M,
        S,
        domain=sp.QQ.frac_field(N, x),
    )

    quotient, remainder = sp.div(
        poly_H,
        poly_M,
    )

    remainder = sp.factor(
        remainder.as_expr()
    )

    print(f"--- (k,l)=({k},{ell}) ---")
    print(
        f"  remainder = {remainder}"
    )

    # Normalize x^2 = N*y.
    # The remainder may contain x denominators. We still perform the
    # substitution exactly after clearing the obvious x/N scale.
    if remainder != 0:
        rem = sp.factor(remainder)

        # Replace x^2 by N*y by introducing a temporary substitution
        # after collecting in even powers. For this audit we only inspect
        # the factorized numerator.
        num, den = sp.fraction(rem)

        num_y = sp.expand(
            num.subs(x**2, N*y)
        )

        # SymPy may not recursively replace powers > 2, so use Poly in x
        # and reconstruct powers of y exactly.
        poly_num = sp.Poly(num, x)

        if all(
            monom[0] % 2 == 0
            for monom, _ in poly_num.terms()
        ):
            normalized = sp.Integer(0)

            for monom, coeff in poly_num.terms():
                power_x = monom[0]
                normalized += coeff * (
                    N*y
                )**(power_x // 2)

            print(
                "  numerator as polynomial in y:"
            )
            print(
                f"    {sp.factor(normalized)}"
            )
        else:
            print(
                "  numerator contains odd x-powers; "
                "no pure y=x^2/N normalization."
            )

    print()


# ============================================================================
# [12] Search whether M_x appears through a coefficient combination
# ============================================================================

print("[12] LINEAR COEFFICIENT-COMBINATION SEARCH")
print("-" * 78)

print(
    "The direct coefficient search found S as the terminal coefficient."
)
print(
    "Now test whether simple neighboring coefficient combinations produce"
)
print(
    "M_x = N + S*x + x^2."
)
print()

combination_hits = []

for k, ell in PAIRS:
    H = shifted_H(k, ell)
    coeffs = [
        sp.expand(
            sp.Poly(H, x).coeff_monomial(x**r)
        )
        for r in range(ell)
    ]

    # Test C_top + x*C_(top-1)
    if ell >= 2:
        C_top = coeffs[-1]
        C_prev = coeffs[-2]

        candidate_1 = sp.expand(
            C_top * x + C_prev
        )

        candidate_2 = sp.expand(
            C_top + C_prev * x
        )

        target = N + S*x + x**2

        if exact_zero(candidate_1 - target):
            combination_hits.append(
                (k, ell, "x*C_top + C_prev")
            )

        if exact_zero(candidate_2 - target):
            combination_hits.append(
                (k, ell, "C_top + x*C_prev")
            )

    # Test C_top*(x+1), which naturally generates S*x+S.
    C_top = coeffs[-1]

    if exact_zero(
        sp.expand(C_top * (x + 1))
        - S * (x + 1)
    ):
        print(
            f"  (k,l)=({k},{ell}) has C_top*(x+1)=S*(x+1)"
        )

print(
    f"  DIRECT M_x COMBINATION HITS = "
    f"{len(combination_hits)}"
)

if combination_hits:
    for hit in combination_hits:
        print(f"    {hit}")

print()


# ============================================================================
# [13] Exact relation to the original kernel's top symmetric coordinate
# ============================================================================

print("[13] ORIGINAL KAPPA KERNEL / TOP SHIFT COEFFICIENT")
print("-" * 78)

print(
    "At x=1:"
)
print(
    "  F_(k,l) = F_shifted(k,l;x=1)."
)
print()
print(
    "The top shifted coefficient is nevertheless:"
)
print(
    "  C_(l-1) = p^k + q^k."
)
print()
print(
    "The experiment now records the exact symmetric form of the"
)
print(
    "x=1 specialization and its relation to P_k."
)
print()

for k, ell in [
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
]:
    H = shifted_H(k, ell)

    H1 = sp.expand(H.subs(x, 1))
    H1_NS = symmetric_reduce(H1)

    top = symmetric_reduce(
        sp.Poly(H, x).coeff_monomial(
            x**(ell - 1)
        )
    )

    print(f"--- (k,l)=({k},{ell}) ---")
    print(
        f"  H(1) = {sp.factor(H1_NS)}"
    )
    print(
        f"  top coefficient = {sp.factor(top)}"
    )

print()


# ============================================================================
# [14] Independent symbolic proof certificate
# ============================================================================

print("[14] PROOF CERTIFICATE")
print("-" * 78)

print(
    "For symbolic k<ell, the coefficient of x^ell in the numerator F is:"
)
print(
    "  p^k + q^k."
)
print()
print(
    "Reason:"
)
print(
    "  p^k(q+x)^ell contributes p^k*x^ell."
)
print(
    "  q^k(p+x)^ell contributes q^k*x^ell."
)
print(
    "  -p^ell(q+x)^k and -q^ell(p+x)^k have degree <= k < ell."
)
print()
print(
    "Therefore:"
)
print(
    "  [x^ell]F = p^k+q^k."
)
print(
    "Since F=x*H:"
)
print(
    "  [x^(ell-1)]H = p^k+q^k."
)
print()

# Computational certificate for representative symbolic pairs.
proof_certificate_failures = []

for k, ell in [
    (1, 3),
    (1, 7),
    (3, 5),
    (5, 9),
    (7, 11),
]:
    F = shifted_kernel(k, ell)

    coefficient = sp.expand(
        sp.Poly(F, x).coeff_monomial(x**ell)
    )

    expected = p**k + q**k

    ok = exact_zero(
        coefficient - expected
    )

    print(
        f"  coefficient proof check ({k},{ell}) "
        f"PASS={ok}"
    )

    if not ok:
        proof_certificate_failures.append(
            (k, ell)
        )

print(
    f"  PROOF CERTIFICATE FAILURES = "
    f"{len(proof_certificate_failures)}"
)
print()


# ============================================================================
# [15] Research conclusion
# ============================================================================

print("[15] RESEARCH CONCLUSION")
print("-" * 78)

all_exact = (
    len(coefficient_failures) == 0
    and len(top_failures) == 0
    and len(support_failures) == 0
    and len(recurrence_failures) == 0
    and len(specialization_failures) == 0
    and len(proof_certificate_failures) == 0
)

print(
    f"  General coefficient formula exact: "
    f"{len(coefficient_failures) == 0}"
)

print(
    f"  Top coefficient theorem exact: "
    f"{len(top_failures) == 0}"
)

print(
    f"  Degree/support theorem exact: "
    f"{len(support_failures) == 0}"
)

print(
    f"  Newton recurrence exact: "
    f"{len(recurrence_failures) == 0}"
)

print(
    f"  x=1 specialization exact: "
    f"{len(specialization_failures) == 0}"
)

print(
    f"  Proof certificates exact: "
    f"{len(proof_certificate_failures) == 0}"
)

print(
    f"  OVERALL SYMBOLIC AUDIT = {all_exact}"
)

print()
print(
    "ESTABLISHED:"
)
print(
    "  C_r(k,l) has an exact closed coefficient formula."
)
print(
    "  C_(l-1)(k,l) = p^k + q^k."
)
print(
    "  The top coefficient is independent of ell."
)
print(
    "  The top coefficients form the Newton power-sum sequence."
)
print(
    "  P_k = S*P_(k-1) - N*P_(k-2)."
)
print(
    "  k=1 gives P_1=S=p+q."
)

print()
print(
    "UNRESOLVED:"
)
print(
    "  Whether the top shifted coefficient can be reconstructed from"
)
print(
    "  the original x=1 / N-side construction without already knowing p,q."
)

print()
print(
    "NEXT TARGET:"
)
print(
    "  Identify an exact transfer operator"
)
print(
    "      original KAPPA kernel"
)
print(
    "            -> top shifted coefficient"
)
print(
    "            -> S=p+q"
)
print(
    "            -> factor recovery."
)
print()


# ============================================================================
# FINISH
# ============================================================================

print("=" * 78)
print("EXPERIMENT 480 FINISHED")
print("=" * 78)
