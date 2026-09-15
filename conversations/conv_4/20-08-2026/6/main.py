#!/usr/bin/env python3

"""
EXPERIMENT 483
==============================================================================
DIRECT TWO-VALUE RECOVERY OF S FROM THE ORIGINAL KAPPA KERNEL
==============================================================================

Main target:

    F_l =
        p(1+q)^l
        + q(1+p)^l
        - p^l(1+q)
        - q^l(1+p)

For k=1, the experiment tests the exact identities

    F_2 = 6N - S^2 + S

    F_3 = (S+1)F_2

and therefore, whenever F_2 != 0,

    S = F_3/F_2 - 1.

The resulting factorization bridge is

    F_2,F_3 -> S -> z^2-Sz+N -> p,q.

This experiment does not establish

    N -> F_2,F_3.

No floating point.
No factorization of N.
No interpolation.
No CSV.
"""

from __future__ import annotations

import math
import sympy as sp


# ============================================================================
# START
# ============================================================================

print("EXPERIMENT 483 START")
print("=" * 78)
print("DIRECT TWO-VALUE RECOVERY OF S FROM THE ORIGINAL KAPPA KERNEL")
print("=" * 78)
print()


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ============================================================================
# Test instances
# ============================================================================

TEST_PAIRS = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]


# ============================================================================
# Helpers
# ============================================================================

def exact_zero(expr) -> bool:
    return sp.cancel(sp.expand(expr)) == 0


def symmetric_reduce(expr):
    """
    Exact p,q -> S,N reduction.
    """
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Symmetrization failed: {remainder}"
        )

    result = reduced

    for formal_var, concrete in mapping:
        if exact_zero(concrete - (p + q)):
            result = result.subs(formal_var, S)
        elif exact_zero(concrete - p * q):
            result = result.subs(formal_var, N)
        else:
            raise AssertionError(
                f"Unexpected symmetric mapping: {mapping}"
            )

    result = sp.expand(result)

    check = sp.expand(
        result.subs({
            S: p + q,
            N: p * q,
        }) - expr
    )

    if not exact_zero(check):
        raise AssertionError(
            "Symmetric reconstruction failed."
        )

    return result


def F(ell: int):
    """
    Original k=1 KAPPA kernel.
    """
    return sp.expand(
        p * (1 + q) ** ell
        + q * (1 + p) ** ell
        - p ** ell * (1 + q)
        - q ** ell * (1 + p)
    )


# ============================================================================
# [1] Direct symbolic values
# ============================================================================

print("[1] DIRECT SYMBOLIC VALUES")
print("-" * 78)

F2_pq = sp.factor(F(2))
F3_pq = sp.factor(F(3))

print(f"  F_2(p,q) = {F2_pq}")
print(f"  F_3(p,q) = {F3_pq}")
print()

F2_NS = symmetric_reduce(F2_pq)
F3_NS = symmetric_reduce(F3_pq)

print(f"  F_2(N,S) = {sp.factor(F2_NS)}")
print(f"  F_3(N,S) = {sp.factor(F3_NS)}")
print()


# ============================================================================
# [2] Exact two-value identity
# ============================================================================

print("[2] DIRECT TWO-VALUE IDENTITY")
print("-" * 78)

identity_difference = sp.factor(
    F3_NS - (S + 1) * F2_NS
)

identity_pass = exact_zero(
    identity_difference
)

print("Target:")
print("  F_3 = (S+1) F_2")
print()

print(
    f"  F_3 - (S+1)F_2 = {identity_difference}"
)

print(
    f"  SYMBOLIC IDENTITY PASS = {identity_pass}"
)
print()


# ============================================================================
# [3] Direct recovery formula
# ============================================================================

print("[3] EXACT RECOVERY FORMULA")
print("-" * 78)

recovered_S = sp.cancel(
    F3_NS / F2_NS - 1
)

ratio_identity = sp.factor(
    recovered_S - S
)

ratio_pass = exact_zero(
    ratio_identity
)

print("Whenever F_2 != 0:")
print()
print("  S = F_3/F_2 - 1")
print()

print(
    f"  F_3/F_2 - 1 = {recovered_S}"
)

print(
    f"  (F_3/F_2 - 1) - S = {ratio_identity}"
)

print(
    f"  RATIO RECOVERY PASS = {ratio_pass}"
)
print()


# ============================================================================
# [4] F2 in p,q
# ============================================================================

print("[4] EXCEPTIONAL LOCUS F_2 = 0")
print("-" * 78)

F2_expanded_pq = sp.expand(F2_pq)

expected_F2 = (
    4 * p * q
    - p ** 2
    - q ** 2
    + p
    + q
)

exception_check = exact_zero(
    F2_expanded_pq - expected_F2
)

print(
    f"  F_2 = {F2_expanded_pq}"
)
print()

print(
    "Equivalent:"
)
print(
    "  F_2 = 4pq - p^2 - q^2 + p + q"
)

print(
    f"  EXCEPTIONAL-LOCUS FORM PASS = "
    f"{exception_check}"
)
print()


# ============================================================================
# [5] Exceptional locus as quadratic in q
# ============================================================================

print("[5] EXCEPTIONAL LOCUS AS A DIOPHANTINE EQUATION")
print("-" * 78)

exception_quadratic = sp.expand(
    q ** 2
    - (4 * p + 1) * q
    + p ** 2
    - p
)

disc_q = sp.factor(
    sp.discriminant(
        exception_quadratic,
        q,
    )
)

print(
    "Setting F_2 = 0 gives:"
)
print(
    "  q^2 - (4p+1)q + p^2 - p = 0"
)
print()

print(
    "Discriminant in q:"
)
print(
    f"  Delta_q = {disc_q}"
)
print()

print(
    "Therefore integer exceptional points require:"
)
print(
    "  12*p^2 + 12*p + 1"
)
print(
    "to be a perfect square."
)
print()


# ============================================================================
# [6] Small positive integer exceptional-point search
# ============================================================================

print("[6] SMALL INTEGER EXCEPTIONAL-POINT SEARCH")
print("-" * 78)

INTEGER_LIMIT = 500

integer_exceptions = []

for pp in range(1, INTEGER_LIMIT + 1):

    disc = 12 * pp ** 2 + 12 * pp + 1

    # Correct exact integer square root:
    root = math.isqrt(disc)

    if root * root == disc:

        numerator_plus = 4 * pp + 1 + root
        numerator_minus = 4 * pp + 1 - root

        candidates = []

        if numerator_plus % 2 == 0:
            candidates.append(
                numerator_plus // 2
            )

        if numerator_minus % 2 == 0:
            candidates.append(
                numerator_minus // 2
            )

        for qq in candidates:
            if qq >= 1:

                value = (
                    4 * pp * qq
                    - pp ** 2
                    - qq ** 2
                    + pp
                    + qq
                )

                if value == 0:
                    integer_exceptions.append(
                        (pp, qq, root)
                    )

print(
    f"  Positive integer p tested through "
    f"{INTEGER_LIMIT}"
)

print(
    f"  Exceptional integer points found = "
    f"{len(integer_exceptions)}"
)

for pp, qq, root in integer_exceptions[:30]:
    print(
        f"    p={pp}, q={qq}, sqrt(Delta)={root}"
    )

if len(integer_exceptions) > 30:
    print(
        f"    ... {len(integer_exceptions) - 30} more"
    )

print()


# ============================================================================
# [7] Prime-pair exceptional search
# ============================================================================

print("[7] PRIME-PAIR EXCEPTIONAL SEARCH")
print("-" * 78)

prime_limit = 2000
primes = list(
    sp.primerange(
        2,
        prime_limit + 1,
    )
)

prime_exceptions = []

for i, pp in enumerate(primes):

    for qq in primes[i + 1:]:

        value = (
            4 * pp * qq
            - pp ** 2
            - qq ** 2
            + pp
            + qq
        )

        if value == 0:
            prime_exceptions.append(
                (pp, qq)
            )

print(
    f"  Prime limit = {prime_limit}"
)

print(
    f"  Prime-pair exceptional cases = "
    f"{len(prime_exceptions)}"
)

if prime_exceptions:
    for pp, qq in prime_exceptions:
        print(
            f"    ({pp}, {qq})"
        )
else:
    print(
        "  No prime-pair exceptional cases found."
    )

print()


# ============================================================================
# [8] Nearby ratio structure
# ============================================================================

print("[8] LOW-INDEX RATIO STRUCTURE")
print("-" * 78)

values = {
    ell: symmetric_reduce(
        F(ell)
    )
    for ell in range(2, 9)
}

for ell in range(2, 8):

    ratio = sp.factor(
        sp.cancel(
            values[ell + 1]
            / values[ell]
        )
    )

    print(
        f"  F_{ell+1}/F_{ell} = {ratio}"
    )

print()


# ============================================================================
# [9] Primary ratio vs next ratios
# ============================================================================

print("[9] RATIO SIMPLICITY AUDIT")
print("-" * 78)

ratio_23 = sp.factor(
    sp.cancel(
        values[3] / values[2]
    )
)

ratio_34 = sp.factor(
    sp.cancel(
        values[4] / values[3]
    )
)

print(
    f"  F_3/F_2 = {ratio_23}"
)

print(
    "  Expected = S+1"
)

print()

print(
    f"  F_4/F_3 = {ratio_34}"
)

print()

ratio_23_pass = exact_zero(
    ratio_23 - (S + 1)
)

print(
    f"  F_3/F_2 = S+1 PASS = {ratio_23_pass}"
)

print()


# ============================================================================
# [10] Direct algebraic proof
# ============================================================================

print("[10] DIRECT ALGEBRAIC PROOF CERTIFICATE")
print("-" * 78)

proof_difference = sp.factor(
    F(3) - (p + q + 1) * F(2)
)

proof_pass = exact_zero(
    proof_difference
)

print(
    "  F_3 - (p+q+1)F_2 = "
    f"{proof_difference}"
)

print(
    f"  PROOF CERTIFICATE PASS = {proof_pass}"
)

print()


# ============================================================================
# [11] Numerical S recovery
# ============================================================================

print("[11] NUMERICAL S RECOVERY")
print("-" * 78)

numerical_failures = []

for pp, qq in TEST_PAIRS:

    NN = pp * qq
    SS = pp + qq

    F2_num = (
        pp * (qq + 1) ** 2
        + qq * (pp + 1) ** 2
        - pp ** 2 * (qq + 1)
        - qq ** 2 * (pp + 1)
    )

    F3_num = (
        pp * (qq + 1) ** 3
        + qq * (pp + 1) ** 3
        - pp ** 3 * (qq + 1)
        - qq ** 3 * (pp + 1)
    )

    print(
        f"--- p={pp}, q={qq} ---"
    )

    print(
        f"  N = {NN}"
    )

    print(
        f"  true S = {SS}"
    )

    print(
        f"  F_2 = {F2_num}"
    )

    print(
        f"  F_3 = {F3_num}"
    )

    if F2_num == 0:

        print(
            "  F_2 = 0 -> primary ratio unavailable"
        )

        numerical_failures.append(
            (
                pp,
                qq,
                "F2=0",
            )
        )

        print()
        continue

    recovered = sp.Rational(
        F3_num,
        F2_num,
    ) - 1

    recovered = sp.simplify(
        recovered
    )

    ok = (
        recovered == SS
    )

    print(
        f"  recovered S = {recovered}"
    )

    print(
        f"  PASS = {ok}"
    )

    if not ok:
        numerical_failures.append(
            (
                pp,
                qq,
                "incorrect S",
            )
        )

    print()


print(
    f"  NUMERICAL S RECOVERY FAILURES = "
    f"{len(numerical_failures)}"
)

for failure in numerical_failures:
    print(
        f"    FAILURE: {failure}"
    )

print()


# ============================================================================
# [12] Factor recovery
# ============================================================================

print("[12] FACTOR RECOVERY AFTER TWO-VALUE BRIDGE")
print("-" * 78)

factor_failures = []

for pp, qq in TEST_PAIRS:

    NN = pp * qq

    F2_num = (
        pp * (qq + 1) ** 2
        + qq * (pp + 1) ** 2
        - pp ** 2 * (qq + 1)
        - qq ** 2 * (pp + 1)
    )

    F3_num = (
        pp * (qq + 1) ** 3
        + qq * (pp + 1) ** 3
        - pp ** 3 * (qq + 1)
        - qq ** 3 * (pp + 1)
    )

    if F2_num == 0:
        continue

    recovered_S = sp.Rational(
        F3_num,
        F2_num,
    ) - 1

    discriminant = sp.expand(
        recovered_S ** 2 - 4 * NN
    )

    # Exact square-root check.
    discriminant_int = int(discriminant)

    if discriminant_int < 0:
        factor_failures.append(
            (
                pp,
                qq,
                "negative discriminant",
            )
        )
        continue

    root_disc = math.isqrt(
        discriminant_int
    )

    if root_disc * root_disc != discriminant_int:

        print(
            f"  ({pp},{qq}) -> discriminant is not a square"
        )

        factor_failures.append(
            (
                pp,
                qq,
                "non-square discriminant",
            )
        )

        continue

    roots = [
        (int(recovered_S) - root_disc) // 2,
        (int(recovered_S) + root_disc) // 2,
    ]

    roots_sorted = sorted(
        roots
    )

    expected = sorted(
        [pp, qq]
    )

    ok = (
        roots_sorted == expected
        and int(recovered_S) % 2 == root_disc % 2
    )

    print(
        f"  ({pp},{qq}) -> "
        f"recovered roots {roots_sorted} "
        f"PASS={ok}"
    )

    if not ok:
        factor_failures.append(
            (
                pp,
                qq,
                "wrong recovered roots",
            )
        )

print(
    f"  FACTOR RECOVERY FAILURES = "
    f"{len(factor_failures)}"
)

print()


# ============================================================================
# [13] Information-model distinction
# ============================================================================

print("[13] INFORMATION-MODEL DISTINCTION")
print("-" * 78)

print(
    "Established:"
)

print(
    "  F_2,F_3 -> S exactly"
)

print(
    "  S,N -> p,q exactly"
)

print()

print(
    "Therefore:"
)

print(
    "  F_2,F_3,N -> p,q"
)

print()

print(
    "Still unresolved:"
)

print(
    "  N -> F_2,F_3"
)

print()

print(
    "The two-value bridge is substantially simpler than the"
)
print(
    "order-four recurrence."
)

print()


# ============================================================================
# [14] Final theorem candidate
# ============================================================================

print("[14] FINAL THEOREM CANDIDATE")
print("-" * 78)

print(
    "For the k=1 original KAPPA kernel:"
)

print()

print(
    "  F_l = p(1+q)^l + q(1+p)^l"
)

print(
    "        - p^l(1+q) - q^l(1+p)"
)

print()

print(
    "The first two nonzero values satisfy:"
)

print(
    "  F_2 = 6N - S^2 + S"
)

print(
    "  F_3 = (S+1)F_2"
)

print()

print(
    "Hence, whenever F_2 != 0:"
)

print(
    "  S = F_3/F_2 - 1"
)

print()

print(
    "and therefore:"
)

print(
    "  z^2 - Sz + N = 0"
)

print(
    "recovers p and q."
)

print()


# ============================================================================
# [15] Overall status
# ============================================================================

theorem_pass = (
    identity_pass
    and ratio_pass
    and exception_check
    and ratio_23_pass
    and proof_pass
    and len(numerical_failures) == 0
    and len(factor_failures) == 0
)

print("[15] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  F3=(S+1)F2 exact: {identity_pass}"
)

print(
    f"  S=F3/F2-1 exact: {ratio_pass}"
)

print(
    f"  F2 exceptional-locus formula: {exception_check}"
)

print(
    f"  Ratio simplification F3/F2=S+1: "
    f"{ratio_23_pass}"
)

print(
    f"  Direct proof certificate: {proof_pass}"
)

print(
    f"  Numerical S recovery: "
    f"{len(numerical_failures) == 0}"
)

print(
    f"  Numerical factor recovery: "
    f"{len(factor_failures) == 0}"
)

print(
    f"  OVERALL EXACT AUDIT = {theorem_pass}"
)

print()

print(
    "MAIN RESEARCH RESULT:"
)

print(
    "  The original KAPPA kernel exposes S using only TWO"
)

print(
    "  exact kernel values:"
)

print()

print(
    "      F_2, F_3"
)

print(
    "          -> S = F_3/F_2 - 1"
)

print(
    "          -> z^2-S*z+N"
)

print(
    "          -> p,q"
)

print()

print(
    "REMAINING CENTRAL QUESTION:"
)

print(
    "  Can F_2 and F_3 be generated from N alone"
)

print(
    "  by an exact source-free construction?"
)

print()


# ============================================================================
# FINISH
# ============================================================================

print("=" * 78)
print("EXPERIMENT 483 FINISHED")
print("=" * 78)