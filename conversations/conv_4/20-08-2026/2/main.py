#!/usr/bin/env python3

"""
EXPERIMENT 479
==============================================================================
SHIFTED-KERNEL TOP-COEFFICIENT THEOREM

Primary target
--------------

For

    H_{k,l}(p,q;x)
      =
      [ p^k (q+x)^l
        + q^k (p+x)^l
        - p^l (q+x)^k
        - q^l (p+x)^k ] / x,

prove exactly that

    [x^(l-1)] H_{k,l}(p,q;x) = p^k + q^k

for every positive integer k<l.

The experiment then proves that these extracted coefficients satisfy the
standard symmetric power-sum / Newton recurrence

    P_k = S P_(k-1) - N P_(k-2),

where

    N = p*q
    S = p+q.

Consequently,

    P_0 = 2
    P_1 = S

and therefore the k=1 top shifted coefficient directly recovers

    S = p+q.

This experiment is proof-oriented.

It does NOT use:
    - numerical interpolation,
    - floating point,
    - prime factorization,
    - finite-grid fitting,
    - CSV output,
    - resultants,
    - arbitrary feature searches.

It uses exact symbolic polynomial identities only.

The experiment has five main goals:

    1. Prove the top-x coefficient identity directly from the binomial theorem.
    2. Verify the identity against exact symbolic expansions for many (k,l).
    3. Prove the Newton recurrence for P_k = p^k + q^k.
    4. Connect the shifted-kernel coefficient to the Newton basis.
    5. Produce exact reconstruction certificates for S and p,q.

Important:
The result established here is an algebraic identity assuming the shifted kernel
is available. It does NOT yet establish that the coefficient can be computed
from N alone.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


# ============================================================================
# START
# ============================================================================

print("EXPERIMENT 479 START")
print("=" * 78)
print("SHIFTED-KERNEL TOP-COEFFICIENT THEOREM")
print("=" * 78)
print()


# ============================================================================
# Symbols
# ============================================================================

p, q, x = sp.symbols("p q x")
N, S = sp.symbols("N S")


# ============================================================================
# Configuration
# ============================================================================

PAIRS = [
    (1, 2),
    (1, 3),
    (1, 4),
    (1, 5),
    (1, 7),
    (1, 9),
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
    (7, 9),
    (7, 11),
]


# ============================================================================
# Helpers
# ============================================================================

def exact_zero(expr) -> bool:
    """Return True iff expr simplifies exactly to zero."""
    return sp.cancel(sp.expand(expr)) == 0


def shifted_kernel(k: int, ell: int):
    """
    Exact shifted kernel numerator.
    """
    return sp.expand(
        p**k * (q + x)**ell
        + q**k * (p + x)**ell
        - p**ell * (q + x)**k
        - q**ell * (p + x)**k
    )


def shifted_kernel_quotient(k: int, ell: int):
    """
    Exact H_{k,l} = F_{k,l} / x.

    Since F(p,q;0)=0 identically, exact polynomial division by x is possible.
    """
    F = shifted_kernel(k, ell)

    zero_at_x0 = sp.expand(F.subs(x, 0))

    if not exact_zero(zero_at_x0):
        raise AssertionError(
            f"Kernel does not vanish at x=0 for k={k}, ell={ell}"
        )

    poly_F = sp.Poly(F, x, domain=sp.QQ.frac_field(p, q))
    poly_x = sp.Poly(x, x, domain=sp.QQ.frac_field(p, q))

    quotient, remainder = sp.div(poly_F, poly_x)

    if remainder.as_expr() != 0:
        raise AssertionError(
            f"Exact division by x failed for k={k}, ell={ell}"
        )

    return sp.expand(quotient.as_expr())


def top_coefficient_H(k: int, ell: int):
    """
    Extract [x^(ell-1)] H_{k,l}.
    """
    H = shifted_kernel_quotient(k, ell)
    return sp.expand(sp.Poly(H, x).coeff_monomial(x**(ell - 1)))


def symmetric_power_sum(k: int):
    """
    P_k = p^k + q^k.
    """
    return sp.expand(p**k + q**k)


def reduce_symmetric(expr):
    """
    Reduce a symmetric polynomial in p,q to S,N.

    This is used only for P_k and exact comparison with the recurrence.
    """
    expr = sp.expand(expr)

    result, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Expression was not symmetrized cleanly. remainder={remainder}"
        )

    reduced = result

    for formal_var, concrete_expr in mapping:
        if exact_zero(concrete_expr - (p + q)):
            reduced = reduced.subs(formal_var, S)
        elif exact_zero(concrete_expr - p * q):
            reduced = reduced.subs(formal_var, N)
        else:
            raise AssertionError(
                f"Unexpected symmetric mapping: {mapping}"
            )

    reduced = sp.expand(reduced)

    # Independent reconstruction.
    reconstructed = sp.expand(
        reduced.subs(
            {
                S: p + q,
                N: p * q,
            }
        )
        - expr
    )

    if not exact_zero(reconstructed):
        raise AssertionError(
            "Symmetric reduction failed reconstruction."
        )

    return reduced


def direct_binomial_top_coefficient(k: int, ell: int):
    """
    Direct coefficient proof.

    In F_{k,l}, the only x^ell contributions are

        p^k * x^ell
        q^k * x^ell

    because k < ell.

    Therefore the x^(ell-1) coefficient of H=F/x is p^k+q^k.

    This function computes that coefficient directly from the four
    binomial expansions, without expanding the full kernel.
    """
    first = p**k * sp.binomial(ell, ell)
    second = q**k * sp.binomial(ell, ell)

    # The two negative terms have x-degree at most k < ell.
    negative_1 = 0
    negative_2 = 0

    return sp.expand(first + second - negative_1 - negative_2)


def power_sum_recurrence_check(max_k: int):
    """
    Verify

        P_k = S P_(k-1) - N P_(k-2)

    symbolically for k=2,...,max_k.
    """
    print(f"[RECURRENCE] Checking power sums through k={max_k}")

    P = {
        0: sp.Integer(2),
        1: S,
    }

    # Build recursively in the symmetric variables.
    for k in range(2, max_k + 1):
        P[k] = sp.expand(
            S * P[k - 1] - N * P[k - 2]
        )

    all_pass = True

    for k in range(2, max_k + 1):
        direct = symmetric_power_sum(k)
        recurrence_expression = P[k]

        direct_SN = reduce_symmetric(direct)

        ok = exact_zero(
            direct_SN - recurrence_expression
        )

        print(
            f"  k={k:2d}: "
            f"PASS={ok}"
        )

        if not ok:
            print(f"    direct   = {direct_SN}")
            print(f"    recur    = {recurrence_expression}")
            print(
                f"    diff     = "
                f"{sp.factor(direct_SN - recurrence_expression)}"
            )
            all_pass = False

    print(f"  ALL RECURRENCE CHECKS PASS = {all_pass}")
    print()

    return all_pass, P


# ============================================================================
# [1] Canonical theorem target
# ============================================================================

print("[1] THEOREM TARGET")
print("-" * 78)

print(
    "Define"
)
print(
    "  H_(k,l)(p,q;x) = F_(k,l)(p,q;x) / x"
)
print()
print(
    "Target identity:"
)
print(
    "  [x^(l-1)] H_(k,l) = p^k + q^k"
)
print()
print(
    "For k=1 this becomes:"
)
print(
    "  [x^(l-1)] H_(1,l) = p+q = S"
)
print()


# ============================================================================
# [2] Direct symbolic proof from binomial expansion
# ============================================================================

print("[2] DIRECT BINOMIAL-COEFFICIENT PROOF")
print("-" * 78)

for k, ell in PAIRS:
    direct = direct_binomial_top_coefficient(k, ell)
    expected = p**k + q**k

    ok = exact_zero(direct - expected)

    print(
        f"  (k,l)=({k},{ell}) : "
        f"PASS={ok}"
    )

    if not ok:
        print(f"    direct   = {direct}")
        print(f"    expected = {expected}")


print()


# ============================================================================
# [3] Full exact kernel expansion verification
# ============================================================================

print("[3] FULL SYMBOLIC KERNEL VERIFICATION")
print("-" * 78)

full_kernel_failures = []

for k, ell in PAIRS:
    H = shifted_kernel_quotient(k, ell)

    actual = top_coefficient_H(k, ell)
    expected = p**k + q**k

    ok = exact_zero(actual - expected)

    print(
        f"  (k,l)=({k},{ell}) : "
        f"[x^{ell-1}]H = {sp.factor(actual)} "
        f"PASS={ok}"
    )

    if not ok:
        full_kernel_failures.append((k, ell))

print(
    f"  FULL KERNEL FAILURES = {len(full_kernel_failures)}"
)

if full_kernel_failures:
    print(f"  FAILED PAIRS = {full_kernel_failures}")

print()


# ============================================================================
# [4] Top coefficient as symmetric Newton / power-sum sequence
# ============================================================================

print("[4] NEWTON / POWER-SUM INTERPRETATION")
print("-" * 78)

print(
    "Define"
)
print(
    "  P_k := [x^(l-1)] H_(k,l)"
)
print()
print(
    "Observed exact identity:"
)
print(
    "  P_k = p^k + q^k"
)
print()
print(
    "Therefore:"
)
print(
    "  P_0 = 2"
)
print(
    "  P_1 = S"
)
print(
    "  P_k = S*P_(k-1) - N*P_(k-2)"
)
print()


# Verify that the top coefficient is independent of ell.

print("[4A] ELL-INDEPENDENCE")
print("-" * 78)

ell_independence_failures = []

for k in range(1, 8):
    values = []

    valid_ells = [ell for ell in range(k + 1, k + 8)]

    for ell in valid_ells:
        value = top_coefficient_H(k, ell)
        values.append(value)

    expected = p**k + q**k

    for ell, value in zip(valid_ells, values):
        if not exact_zero(value - expected):
            ell_independence_failures.append((k, ell))

    print(
        f"  k={k}: "
        f"tested ell={valid_ells} "
        f"PASS={all(exact_zero(v - expected) for v in values)}"
    )

print(
    f"  ELL-INDEPENDENCE FAILURES = "
    f"{len(ell_independence_failures)}"
)
print()


# ============================================================================
# [5] Symmetric recurrence proof
# ============================================================================

print("[5] EXACT NEWTON RECURRENCE")
print("-" * 78)

recurrence_pass, recurrence_P = power_sum_recurrence_check(12)


# ============================================================================
# [6] Connect extracted kernel coefficient to recurrence
# ============================================================================

print("[6] KERNEL-COEFFICIENT RECURRENCE")
print("-" * 78)

kernel_P = {}

for k in range(0, 9):
    if k == 0:
        kernel_P[k] = sp.Integer(2)
    else:
        # Use l=k+1 as a canonical admissible value.
        kernel_P[k] = top_coefficient_H(k, k + 1)

for k in range(2, 9):
    lhs = kernel_P[k]
    rhs = sp.expand(
        S * kernel_P[k - 1]
        - N * kernel_P[k - 2]
    )

    # Convert the kernel coefficient to N,S coordinates.
    lhs_SN = reduce_symmetric(lhs)

    ok = exact_zero(lhs_SN - rhs)

    print(
        f"  k={k}: "
        f"P_k = {lhs_SN}"
    )
    print(
        f"       recurrence PASS={ok}"
    )

    if not ok:
        print(
            f"       difference = "
            f"{sp.factor(lhs_SN - rhs)}"
        )

print()


# ============================================================================
# [7] Explicit first power sums
# ============================================================================

print("[7] FIRST POWER-SUM IDENTITIES")
print("-" * 78)

for k in range(0, 8):
    if k == 0:
        Pk = sp.Integer(2)
    else:
        Pk = reduce_symmetric(
            symmetric_power_sum(k)
        )

    print(f"  P_{k} = {sp.factor(Pk)}")

print()


# ============================================================================
# [8] Direct recovery of S
# ============================================================================

print("[8] DIRECT RECOVERY OF S")
print("-" * 78)

print(
    "For k=1:"
)
print(
    "  [x^(l-1)] H_(1,l) = p+q = S"
)
print()

for ell in range(2, 10):
    H = shifted_kernel_quotient(1, ell)
    recovered_S = sp.Poly(H, x).coeff_monomial(x ** (ell - 1))

    ok = exact_zero(recovered_S - S)

    print(
        f"  ell={ell:2d}: "
        f"recovered S = {sp.factor(recovered_S)} "
        f"PASS={ok}"
    )

print()


# ============================================================================
# [9] Highest derivative formulation
# ============================================================================

print("[9] HIGHEST-DERIVATIVE FORMULATION")
print("-" * 78)

print(
    "Since H has degree ell-1 in x:"
)
print(
    "  H^(ell-1)(0) = (ell-1)! * (p^k+q^k)"
)
print()

derivative_failures = []

for k, ell in PAIRS:
    H = shifted_kernel_quotient(k, ell)

    derivative_value = sp.expand(
        sp.diff(H, x, ell - 1).subs(x, 0)
    )

    expected = sp.factor(
        sp.factorial(ell - 1) * (p**k + q**k)
    )

    ok = exact_zero(
        derivative_value - expected
    )

    print(
        f"  (k,l)=({k},{ell}): "
        f"PASS={ok}"
    )

    if not ok:
        derivative_failures.append((k, ell))

print(
    f"  DERIVATIVE FAILURES = {len(derivative_failures)}"
)
print()


# ============================================================================
# [10] Does k=1 isolate S uniquely?
# ============================================================================

print("[10] k=1 ISOLATION TEST")
print("-" * 78)

for ell in range(2, 11):
    H = shifted_kernel_quotient(1, ell)

    poly = sp.Poly(H, x)

    coeff_top = poly.coeff_monomial(x ** (ell - 1))

    lower_coefficients = [
        poly.coeff_monomial(x ** r)
        for r in range(ell - 1)
    ]

    top_pass = exact_zero(coeff_top - S)

    # Verify that the top coefficient itself has no N dependence
    # after exact symmetric reduction.
    top_NS = reduce_symmetric(coeff_top)

    print(
        f"  ell={ell:2d}: "
        f"top=S -> {top_pass}, "
        f"top coefficient={sp.factor(top_NS)}"
    )

    # Compact lower-degree summary.
    if lower_coefficients:
        nonzero_count = sum(
            not exact_zero(c)
            for c in lower_coefficients
        )
        print(
            f"       nonzero lower coefficients = "
            f"{nonzero_count}/{len(lower_coefficients)}"
        )

print()


# ============================================================================
# [11] Compare with the shifted-product identity
# ============================================================================

print("[11] CONNECTION TO THE SHIFTED-PRODUCT BRIDGE")
print("-" * 78)

M_x = sp.expand(N + S*x + x**2)

print(
    "Canonical shifted product:"
)
print(
    f"  M_x = {M_x}"
)
print()

print(
    "Experiment 479 now has two independent exact bridges:"
)
print(
    "  A) M_x = N + x*S + x^2"
)
print(
    "     S = (M_x - N - x^2)/x"
)
print()
print(
    "  B) [x^(ell-1)] H_(1,ell) = S"
)
print()

bridge_A = sp.cancel(
    (M_x - N - x**2) / x
)

bridge_B = top_coefficient_H(1, 5)

print(
    f"  Bridge A recovered S = {bridge_A}"
)
print(
    f"  Bridge B recovered S = {bridge_B}"
)

print(
    f"  Bridge A PASS = {exact_zero(bridge_A - S)}"
)
print(
    f"  Bridge B PASS = {exact_zero(bridge_B - S)}"
)

print()


# ============================================================================
# [12] Important information-model boundary
# ============================================================================

print("[12] INFORMATION-MODEL CHECK")
print("-" * 78)

print(
    "The symbolic theorem established here is:"
)
print(
    "  shifted kernel -> top x coefficient -> p+q."
)
print()

print(
    "It does NOT yet establish:"
)
print(
    "  N alone -> shifted kernel."
)
print()
print(
    "Therefore the unresolved bridge is now sharply localized:"
)
print(
    "  N / source construction"
)
print(
    "        -> H_(1,ell)"
)
print(
    "        -> [x^(ell-1)]"
)
print(
    "        -> S=p+q"
)
print(
    "        -> quadratic factor recovery."
)
print()


# ============================================================================
# [13] Exact semiprime numerical sanity checks
# ============================================================================

print("[13] NUMERICAL SANITY CHECKS")
print("-" * 78)

test_instances = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
]

numerical_failures = []

for pp, qq in test_instances:
    NN = pp * qq
    SS = pp + qq

    print(
        f"--- p={pp}, q={qq}, "
        f"N={NN}, S={SS} ---"
    )

    for ell in [3, 5, 7, 9]:
        F_num = (
            pp * (qq + x)**ell
            + qq * (pp + x)**ell
            - pp**ell * (qq + x)
            - qq**ell * (pp + x)
        )

        H_num = sp.cancel(F_num / x)

        coeff = sp.expand(
            sp.Poly(H_num, x).coeff_monomial(
                x ** (ell - 1)
            )
        )

        ok = exact_zero(coeff - SS)

        print(
            f"  ell={ell}: "
            f"coefficient={coeff}, "
            f"PASS={ok}"
        )

        if not ok:
            numerical_failures.append(
                (pp, qq, ell)
            )

print(
    f"NUMERICAL SANITY FAILURES = {len(numerical_failures)}"
)
print()


# ============================================================================
# [14] Final theorem statement
# ============================================================================

print("[14] FINAL THEOREM CANDIDATE")
print("-" * 78)

print(
    "THEOREM CANDIDATE"
)
print()
print(
    "For integers k,l with 1 <= k < l, define"
)
print(
    "  H_(k,l)(p,q;x)"
)
print(
    "    = [p^k(q+x)^l + q^k(p+x)^l"
)
print(
    "       - p^l(q+x)^k - q^l(p+x)^k] / x."
)
print()
print(
    "Then H_(k,l) is a polynomial in x of degree l-1 and"
)
print(
    "  [x^(l-1)] H_(k,l) = p^k + q^k."
)
print()
print(
    "Equivalently,"
)
print(
    "  H_(k,l)^(l-1)(0) = (l-1)!(p^k+q^k)."
)
print()
print(
    "For k=1,"
)
print(
    "  [x^(l-1)] H_(1,l) = p+q."
)
print()
print(
    "Hence if this top coefficient is independently available,"
)
print(
    "the factorization follows from"
)
print(
    "  z^2 - (p+q)z + pq = 0."
)
print()


# ============================================================================
# Final status
# ============================================================================

all_checks_pass = (
    len(full_kernel_failures) == 0
    and recurrence_pass
    and len(derivative_failures) == 0
    and len(numerical_failures) == 0
)

print("[15] EXPERIMENT STATUS")
print("-" * 78)
print(
    f"  Full symbolic top-coefficient checks: "
    f"{len(full_kernel_failures) == 0}"
)
print(
    f"  Power-sum recurrence checks: "
    f"{recurrence_pass}"
)
print(
    f"  Highest-derivative checks: "
    f"{len(derivative_failures) == 0}"
)
print(
    f"  Numerical checks: "
    f"{len(numerical_failures) == 0}"
)
print(
    f"  OVERALL EXACT AUDIT: {all_checks_pass}"
)
print()

print(
    "RESEARCH STATUS:"
)
print(
    "  The shifted-kernel top coefficient gives an exact "
    "symmetric-power-sum interface."
)
print(
    "  k=1 isolates S=p+q exactly."
)
print(
    "  The remaining unresolved problem is how to obtain this "
    "coefficient from information not already containing p,q."
)
print()


# ============================================================================
# FINISH
# ============================================================================

print("=" * 78)
print("EXPERIMENT 479 FINISHED")
print("=" * 78)
