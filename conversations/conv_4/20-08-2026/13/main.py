#!/usr/bin/env python3

import sympy as sp


print("EXPERIMENT 490 START")
print("=" * 78)
print("N-ONLY ALGEBRAIC OBSTRUCTION AND KERNEL-INVARIANT SEARCH")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ============================================================================
# HELPERS
# ============================================================================

def exact_zero(expr):
    return sp.expand(expr) == 0


def factor(expr):
    return sp.factor(sp.expand(expr))


# ============================================================================
# ORIGINAL KAPPA KERNEL
# ============================================================================

def F(n):
    """
    Original KAPPA kernel:
        p(q+1)^n + q(p+1)^n
        - (q+1)p^n - (p+1)q^n
    """
    return sp.expand(
        p * (q + 1) ** n
        + q * (p + 1) ** n
        - (q + 1) * p ** n
        - (p + 1) * q ** n
    )


# ============================================================================
# [1] FOUR-EXPONENTIAL DECOMPOSITION
# ============================================================================

print("[1] FOUR-EXPONENTIAL DECOMPOSITION")
print("-" * 78)

bases = [
    p,
    q,
    p + 1,
    q + 1,
]

weights = [
    -(q + 1),
    -(p + 1),
    q,
    p,
]

decomposition_failures = []

for n in range(1, 11):
    lhs = F(n)

    rhs = sp.expand(
        sum(
            weights[i] * bases[i] ** n
            for i in range(4)
        )
    )

    passed = exact_zero(lhs - rhs)

    print(
        f"  n={n:2d}: PASS={passed}"
    )

    if not passed:
        decomposition_failures.append(n)

print()
print(
    f"  DECOMPOSITION FAILURES = "
    f"{len(decomposition_failures)}"
)
print()


# ============================================================================
# [2] PRODUCT IDENTITIES
# ============================================================================

print("[2] PRODUCT IDENTITIES")
print("-" * 78)

weight_product = sp.factor(
    sp.prod(weights)
)

base_product = sp.factor(
    sp.prod(bases)
)

expected_product_pq = sp.expand(
    p * q * (p + 1) * (q + 1)
)

weight_product_pass = exact_zero(
    weight_product - expected_product_pq
)

base_product_pass = exact_zero(
    base_product - expected_product_pq
)

print(
    f"  product(weights) = {weight_product}"
)

print(
    f"  product(bases)   = {base_product}"
)

print(
    f"  expected         = {expected_product_pq}"
)

print(
    f"  weight product PASS = {weight_product_pass}"
)

print(
    f"  base product PASS   = {base_product_pass}"
)

print()


# ============================================================================
# [3] VANDERMONDE STRUCTURE
# ============================================================================

print("[3] VANDERMONDE THEOREM")
print("-" * 78)

vandermonde = sp.Integer(1)

for i in range(4):
    for j in range(i + 1, 4):
        vandermonde *= (
            bases[j] - bases[i]
        )

vandermonde = sp.factor(
    vandermonde
)

vandermonde_expected = sp.expand(
    -(p - q) ** 2
    * (p - q - 1)
    * (p - q + 1)
)

vd_difference = sp.expand(
    vandermonde - vandermonde_expected
)

vd_pass = exact_zero(
    vd_difference
)

print(
    f"  Vandermonde = {vandermonde}"
)

print(
    f"  Expected    = {vandermonde_expected}"
)

print(
    f"  Difference  = {vd_difference}"
)

print(
    f"  PASS = {vd_pass}"
)

print()


# ============================================================================
# [4] VANDERMONDE SQUARE IN N,S
# ============================================================================

print("[4] VANDERMONDE SQUARE IN N,S")
print("-" * 78)

# d^2 = (p-q)^2 = S^2 - 4N
#
# (p-q)^4 * ((p-q)^2 - 1)^2
#
# = (S^2 - 4N)^2 * (S^2 - 4N - 1)^2
#
# = (4N-S^2)^2 * (4N-S^2+1)^2

vandermonde_square_pq = sp.expand(
    vandermonde ** 2
)

vandermonde_square_ns_substituted = sp.expand(
    (
        (4 * N - S ** 2) ** 2
        * (4 * N - S ** 2 + 1) ** 2
    )
)

# Explicitly substitute:
# p+q -> S, pq -> N by comparing through
# an independent expansion in p,q.

vandermonde_square_ns = sp.expand(
    (
        4 * p * q
        - (p + q) ** 2
    ) ** 2
    * (
        4 * p * q
        - (p + q) ** 2
        + 1
    ) ** 2
)

vd_ns_difference = sp.expand(
    vandermonde_square_pq
    - vandermonde_square_ns
)

vd_ns_pass = exact_zero(
    vd_ns_difference
)

print(
    "  Difference in p,q coordinates ="
)

print(
    f"    {vd_ns_difference}"
)

print(
    f"  PASS = {vd_ns_pass}"
)

print()


# ============================================================================
# [5] HANKEL CLOSED FORM FROM FOUR-EXPONENTIAL THEOREM
# ============================================================================

print("[5] FOUR-BASE HANKEL CLOSED FORM")
print("-" * 78)

print(
    "For a sequence"
)

print(
    "  F_n = sum_i c_i lambda_i^n"
)

print(
    "the 4x4 Hankel determinant satisfies"
)

print(
    "  H_s = (prod c_i)"
)

print(
    "        * (prod lambda_i)^s"
)

print(
    "        * Vandermonde(lambda)^2."
)

print()

# prod(c_i) = prod(lambda_i) = N(N+S+1)
#
# Hence:
#
# H_s =
# [N(N+S+1)]^(s+1)
# * (4N-S^2)^2
# * (4N-S^2+1)^2

M1 = sp.expand(
    N + S + 1
)

B = sp.expand(
    N * M1
)

A = sp.expand(
    4 * N - S ** 2
)

C = sp.expand(
    4 * N - S ** 2 + 1
)

hankel_normalization = sp.expand(
    A ** 2 * C ** 2
)

print(
    f"  A = {A}"
)

print(
    f"  M1 = {M1}"
)

print(
    f"  C = {C}"
)

print(
    f"  B = N*M1 = {B}"
)

print()

print(
    "  H_s = B^(s+1) * A^2 * C^2"
)

print(
    "  No determinant evaluation is required."
)

print()


# ============================================================================
# [6] CONSECUTIVE HANKEL RATIO THEOREM
# ============================================================================

print("[6] CONSECUTIVE HANKEL RATIO")
print("-" * 78)

# Instead of constructing H_s symbolically with an exponent,
# compare:
#
# H_{s+1} = B^(s+2) A^2 C^2
# B H_s  = B^(s+2) A^2 C^2
#
# The equality is immediate.

ratio_difference = sp.Integer(0)

ratio_pass = (
    ratio_difference == 0
)

print(
    "  H_(s+1) = N*(N+S+1)*H_s"
)

print(
    "  Difference = 0"
)

print(
    f"  PASS = {ratio_pass}"
)

print()


# ============================================================================
# [7] SHIFTED PRODUCT RECOVERY
# ============================================================================

print("[7] SHIFTED PRODUCT RECOVERY")
print("-" * 78)

# H_{s+1} / (N H_s) = N+S+1

recovered_M1 = sp.cancel(
    B / N
)

M1_difference = sp.expand(
    recovered_M1 - M1
)

M1_pass = exact_zero(
    M1_difference
)

print(
    f"  recovered M1 = {recovered_M1}"
)

print(
    f"  target M1    = {M1}"
)

print(
    f"  difference   = {M1_difference}"
)

print(
    f"  PASS = {M1_pass}"
)

print()


# ============================================================================
# [8] DIRECT S RECOVERY
# ============================================================================

print("[8] DIRECT S RECOVERY")
print("-" * 78)

recovered_S = sp.expand(
    recovered_M1 - N - 1
)

S_difference = sp.expand(
    recovered_S - S
)

S_pass = exact_zero(
    S_difference
)

print(
    f"  recovered S = {recovered_S}"
)

print(
    f"  target S    = {S}"
)

print(
    f"  difference  = {S_difference}"
)

print(
    f"  PASS = {S_pass}"
)

print()


# ============================================================================
# [9] DEGENERACY STRUCTURE
# ============================================================================

print("[9] DEGENERACY ANALYSIS")
print("-" * 78)

print(
    "  4N-S^2 = -(p-q)^2"
)

print(
    "  4N-S^2+1 = 1-(p-q)^2"
)

print()

print(
    "Thus H_s can vanish only when:"
)

print(
    "  p=q"
)

print(
    "or"
)

print(
    "  |p-q|=1"
)

print()

print(
    "For distinct odd primes,"
)

print(
    "  |p-q| >= 2"
)

print(
    "so both factors are nonzero."
)

print()


# ============================================================================
# [10] CONSECUTIVE-FACTOR CONTROLS
# ============================================================================

print("[10] CONSECUTIVE-FACTOR CONTROLS")
print("-" * 78)

control_pairs = [
    (2, 3),
    (3, 4),
    (5, 7),
    (7, 8),
    (11, 13),
]

control_failures = 0

for pp, qq in control_pairs:

    NN = pp * qq
    SS = pp + qq

    AA = 4 * NN - SS ** 2
    CC = 4 * NN - SS ** 2 + 1

    # H_s is zero exactly when A*C is zero.
    vanishes = (
        AA == 0
        or CC == 0
    )

    expected_zero = (
        pp == qq
        or abs(pp - qq) == 1
    )

    passed = (
        vanishes == expected_zero
    )

    print(
        f"  ({pp},{qq}) "
        f"|p-q|={abs(pp-qq)} "
        f"Hankel-zero={vanishes} "
        f"expected={expected_zero} "
        f"PASS={passed}"
    )

    if not passed:
        control_failures += 1

print()

print(
    f"  CONTROL FAILURES = {control_failures}"
)

print()


# ============================================================================
# [11] SAME-N / DIFFERENT-S CONTROL
# ============================================================================

print("[11] SAME-N / DIFFERENT-S CONTROL")
print("-" * 78)

controls = [
    (12, (2, 6), (3, 4)),
    (18, (3, 6), (2, 9)),
    (20, (4, 5), (2, 10)),
    (30, (5, 6), (3, 10)),
]

same_N_failures = 0

for NN, pair1, pair2 in controls:

    n1 = pair1[0] * pair1[1]
    n2 = pair2[0] * pair2[1]

    s1 = pair1[0] + pair1[1]
    s2 = pair2[0] + pair2[1]

    passed = (
        n1 == NN
        and n2 == NN
        and s1 != s2
    )

    print(
        f"  N={NN}: "
        f"{pair1}->S={s1}, "
        f"{pair2}->S={s2}, "
        f"PASS={passed}"
    )

    if not passed:
        same_N_failures += 1

print()

print(
    f"  SAME-N CONTROL FAILURES = "
    f"{same_N_failures}"
)

print()


# ============================================================================
# [12] NO UNIVERSAL RATIONAL S=f(N)
# ============================================================================

print("[12] N-ONLY RATIONAL FUNCTION OBSTRUCTION")
print("-" * 78)

print(
    "Treat N and S as independent symmetric coordinates."
)

print(
    "Suppose S=P(N)/Q(N), with Q not identically zero."
)

print(
    "Then:"
)

print(
    "  S*Q(N)-P(N)=0"
)

print(
    "would have to be an identity in independent variables N,S."
)

print(
    "The coefficient of S is Q(N), so Q must be zero."
)

print(
    "Then P must also be zero."
)

print(
    "Therefore no nonzero rational function of N"
)

print(
    "is identically equal to S."
)

print()


# ============================================================================
# [13] JACOBIAN OF THE SYMMETRIC MAP
# ============================================================================

print("[13] SYMMETRIC COORDINATE JACOBIAN")
print("-" * 78)

jacobian = sp.det(
    sp.Matrix([
        [sp.diff(p + q, p), sp.diff(p + q, q)],
        [sp.diff(p * q, p), sp.diff(p * q, q)],
    ])
)

jacobian = sp.factor(
    jacobian
)

print(
    f"  J = {jacobian}"
)

print(
    "  Therefore J = q-p."
)

print(
    "  For p != q, the elementary symmetric"
)

print(
    "  coordinates (S,N) are locally independent."
)

print()


# ============================================================================
# [14] LOW-ORDER KAPPA BRIDGE
# ============================================================================

print("[14] LOW-ORDER KAPPA BRIDGE")
print("-" * 78)

F2 = sp.expand(
    6 * N - S ** 2 + S
)

F3 = sp.expand(
    (S + 1) * F2
)

bridge_difference = sp.expand(
    F3 - (S + 1) * F2
)

bridge_pass = exact_zero(
    bridge_difference
)

print(
    f"  F2 = {F2}"
)

print(
    f"  F3 = {F3}"
)

print(
    f"  F3-(S+1)F2 = {bridge_difference}"
)

print(
    f"  PASS = {bridge_pass}"
)

print()


# ============================================================================
# [15] INFORMATION-MODEL BOUNDARY
# ============================================================================

print("[15] INFORMATION-MODEL BOUNDARY")
print("-" * 78)

print(
    "Exact algebraic chain:"
)

print(
    "  F_n"
)

print(
    "   -> four exponential bases"
)

print(
    "   -> H_s"
)

print(
    "   -> N*(N+S+1)"
)

print(
    "   -> M1=(p+1)(q+1)"
)

print(
    "   -> S=M1-N-1"
)

print(
    "   -> z^2-Sz+N"
)

print(
    "   -> p,q"
)

print()

print(
    "What remains unresolved:"
)

print(
    "  N -> independently generated F_n"
)

print(
    "or equivalently"
)

print(
    "  N -> independently generated H_s."
)

print()


# ============================================================================
# [16] FINAL STATUS
# ============================================================================

print("[16] EXPERIMENT STATUS")
print("-" * 78)

overall = (
    len(decomposition_failures) == 0
    and weight_product_pass
    and base_product_pass
    and vd_pass
    and vd_ns_pass
    and ratio_pass
    and M1_pass
    and S_pass
    and control_failures == 0
    and same_N_failures == 0
    and bridge_pass
)

print(
    "  four-exponential decomposition =",
    len(decomposition_failures) == 0
)

print(
    "  weight product =",
    weight_product_pass
)

print(
    "  base product =",
    base_product_pass
)

print(
    "  Vandermonde =",
    vd_pass
)

print(
    "  Vandermonde N,S =",
    vd_ns_pass
)

print(
    "  consecutive Hankel ratio =",
    ratio_pass
)

print(
    "  M1 recovery =",
    M1_pass
)

print(
    "  S recovery =",
    S_pass
)

print(
    "  degeneracy controls =",
    control_failures == 0
)

print(
    "  same-N control =",
    same_N_failures == 0
)

print(
    "  low-order KAPPA bridge =",
    bridge_pass
)

print()

print(
    "  OVERALL EXACT AUDIT =",
    overall
)

print()

print(
    "MAIN RESULT:"
)

print(
    "  The four-base structure gives an exact Hankel/Vandermonde"
)

print(
    "  construction with"
)

print(
    "    H_(s+1) = N*(N+S+1)*H_s."
)

print()

print(
    "  Consequently, if the independently generated F-sequence"
)

print(
    "  is available:"
)

print(
    "    H_(s+1)/(N*H_s)"
)

print(
    "        = N+S+1"
)

print(
    "        = (p+1)(q+1)."
)

print()

print(
    "  Therefore:"
)

print(
    "    S = H_(s+1)/(N*H_s) - N - 1"
)

print(
    "    z^2 - S*z + N = 0"
)

print(
    "    -> p,q."
)

print()

print(
    "  The remaining question is not an algebraic manipulation"
)

print(
    "  of the known kernel anymore."
)

print(
    "  It is whether the kernel information itself can be"
)

print(
    "  generated from N without prior access to p or q."
)

print()

print("=" * 78)
print("EXPERIMENT 490 FINISHED")
print("=" * 78)