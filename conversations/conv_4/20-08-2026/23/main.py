#!/usr/bin/env python3

import sympy as sp

print("EXPERIMENT 500 START")
print("=" * 78)
print("SHIFTED-PRODUCT DISCRIMINANT AND GAP-SQUARE BRIDGE")
print("=" * 78)
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")
x, a, z = sp.symbols("x a z")


def fact(expr):
    return sp.factor(sp.expand(expr))


def zero(expr):
    return sp.expand(sp.cancel(expr)) == 0


# ============================================================================
# CORE DEFINITIONS
# ============================================================================

N_pq = sp.expand(p * q)
S_pq = sp.expand(p + q)

delta_NS = sp.expand(S**2 - 4*N)
delta_pq = sp.expand((p - q)**2)


def M_pq(shift):
    return sp.expand((p + shift) * (q + shift))


def M_NS(shift):
    return sp.expand(N + S*shift + shift**2)


# ============================================================================
# 1. BASIC SHIFTED PRODUCT
# ============================================================================

print("[1] SHIFTED-PRODUCT POLYNOMIAL")
print("-" * 78)

Mx_difference = M_pq(x) - M_NS(x)

print("  M_x = (p+x)(q+x)")
print("      = N + S*x + x^2")
print()
print("  Difference =", fact(Mx_difference))
print("  PASS =", zero(Mx_difference))
print()


# ============================================================================
# 2. DISCRIMINANT
# ============================================================================

print("[2] DISCRIMINANT OF THE SHIFTED-PRODUCT POLYNOMIAL")
print("-" * 78)

Mx_poly = sp.Poly(
    x**2 + S*x + N,
    x,
)

poly_discriminant = sp.discriminant(
    Mx_poly.as_expr(),
    x,
)

print("  polynomial = x^2 + S*x + N")
print("  discriminant =", fact(poly_discriminant))
print("  target       =", fact(delta_NS))
print(
    "  difference   =",
    fact(poly_discriminant - delta_NS),
)
print(
    "  PASS =",
    zero(poly_discriminant - delta_NS),
)
print()


# ============================================================================
# 3. GAP-SQUARE IDENTITY
# ============================================================================

print("[3] DISCRIMINANT = FACTOR GAP SQUARE")
print("-" * 78)

delta_gap_difference = sp.expand(
    delta_NS.subs(
        {
            N: N_pq,
            S: S_pq,
        }
    ) - delta_pq
)

print(
    "  S^2 - 4N - (p-q)^2 =",
    fact(delta_gap_difference),
)
print(
    "  PASS =",
    zero(delta_gap_difference),
)
print()


# ============================================================================
# 4. ROOT INTERPRETATION
# ============================================================================

print("[4] ROOT INTERPRETATION")
print("-" * 78)

root_poly = sp.expand(
    x**2 + S_pq*x + N_pq
)

root_p_check = sp.expand(
    root_poly.subs(x, -p)
)

root_q_check = sp.expand(
    root_poly.subs(x, -q)
)

print("  M_x = x^2 + S*x + N")
print("  roots = -p, -q")
print()

print(
    "  root x=-p:",
    fact(root_p_check),
    "PASS=",
    zero(root_p_check),
)

print(
    "  root x=-q:",
    fact(root_q_check),
    "PASS=",
    zero(root_q_check),
)

print()


# ============================================================================
# 5. OPPOSITE SHIFTS
# ============================================================================

print("[5] OPPOSITE SHIFT PRODUCTS")
print("-" * 78)

M_plus = M_NS(x)
M_minus = M_NS(-x)

sum_opposite = sp.expand(M_plus + M_minus)
difference_opposite = sp.expand(M_plus - M_minus)
product_opposite = sp.expand(M_plus * M_minus)

print("  M(+x) =", fact(M_plus))
print("  M(-x) =", fact(M_minus))
print()
print("  M(+x)+M(-x) =", fact(sum_opposite))
print("  M(+x)-M(-x) =", fact(difference_opposite))
print(
    "  (M(+x)-M(-x))/(2x) =",
    fact(sp.cancel(difference_opposite / (2*x))),
)
print()


# ============================================================================
# 6. OPPOSITE-SHIFT PRODUCT IDENTITY
# ============================================================================

print("[6] OPPOSITE-SHIFT PRODUCT IDENTITY")
print("-" * 78)

opposite_target = sp.expand(
    (N + x**2)**2 - S**2*x**2
)

opposite_difference = sp.expand(
    product_opposite - opposite_target
)

print("  M(+x)M(-x) =", fact(product_opposite))
print()
print("  target = (N+x^2)^2 - S^2*x^2")
print("  difference =", fact(opposite_difference))
print("  PASS =", zero(opposite_difference))
print()


# ============================================================================
# 7. DIRECT S^2 RECOVERY
# ============================================================================

print("[7] DIRECT S^2 RECOVERY")
print("-" * 78)

s2_recovered = sp.cancel(
    (
        (N + x**2)**2
        - product_opposite
    )
    / x**2
)

s2_difference = sp.expand(
    s2_recovered - S**2
)

print(
    "  recovered S^2 =",
    fact(s2_recovered),
)
print(
    "  difference from S^2 =",
    fact(s2_difference),
)
print(
    "  PASS =",
    zero(s2_difference),
)
print()


# ============================================================================
# 8. DIRECT DELTA RECOVERY
# ============================================================================

print("[8] DIRECT DELTA RECOVERY")
print("-" * 78)

delta_from_opposite = sp.cancel(
    s2_recovered - 4*N
)

delta_recovery_difference = sp.expand(
    delta_from_opposite - delta_NS
)

print(
    "  recovered Delta =",
    fact(delta_from_opposite),
)
print(
    "  target Delta =",
    fact(delta_NS),
)
print(
    "  difference =",
    fact(delta_recovery_difference),
)
print(
    "  PASS =",
    zero(delta_recovery_difference),
)
print()


# ============================================================================
# 9. GENERAL ±a SHIFTS
# ============================================================================

print("[9] GENERALIZED ±a SHIFT")
print("-" * 78)

M_a = M_NS(a)
M_minus_a = M_NS(-a)

a_product = sp.expand(
    M_a * M_minus_a
)

s2_from_a = sp.cancel(
    (
        (N + a**2)**2
        - a_product
    )
    / a**2
)

delta_from_a = sp.cancel(
    s2_from_a - 4*N
)

delta_a_difference = sp.expand(
    delta_from_a - delta_NS
)

print(
    "  M(+a)M(-a) =",
    fact(a_product),
)
print(
    "  recovered S^2 =",
    fact(s2_from_a),
)
print(
    "  recovered Delta =",
    fact(delta_from_a),
)
print(
    "  Delta difference =",
    fact(delta_a_difference),
)
print(
    "  PASS =",
    zero(delta_a_difference),
)
print()


# ============================================================================
# 10. SPECIAL CASE +1 / -1
# ============================================================================

print("[10] SPECIAL SHIFTS +1 AND -1")
print("-" * 78)

M1 = M_NS(1)
Mm1 = M_NS(-1)

print("  M(+1) =", fact(M1))
print("  M(-1) =", fact(Mm1))
print()

M1_sum = sp.expand(M1 + Mm1)
M1_difference = sp.expand(M1 - Mm1)
M1_product = sp.expand(M1 * Mm1)

print(
    "  M(+1)+M(-1) =",
    fact(M1_sum),
)
print(
    "  M(+1)-M(-1) =",
    fact(M1_difference),
)
print(
    "  M(+1)M(-1) =",
    fact(M1_product),
)
print()


# ============================================================================
# 11. DELTA FROM M(+1)M(-1)
# ============================================================================

print("[11] DELTA FROM M(+1)M(-1)")
print("-" * 78)

delta_from_special = sp.expand(
    (N + 1)**2 - M1_product
)

delta_special_difference = sp.expand(
    delta_from_special - delta_NS
)

print(
    "  (N+1)^2 - M(+1)M(-1) =",
    fact(delta_from_special),
)
print(
    "  target Delta =",
    fact(delta_NS),
)
print(
    "  difference =",
    fact(delta_special_difference),
)
print(
    "  PASS =",
    zero(delta_special_difference),
)
print()


# ============================================================================
# 12. S FROM OPPOSITE SHIFT PRODUCTS
# ============================================================================

print("[12] DIRECT S RECOVERY FROM M(+1), M(-1)")
print("-" * 78)

s_from_special = sp.cancel(
    (M1 - Mm1) / 2
)

s_special_difference = sp.expand(
    s_from_special - S
)

print(
    "  recovered S =",
    fact(s_from_special),
)
print(
    "  difference =",
    fact(s_special_difference),
)
print(
    "  PASS =",
    zero(s_special_difference),
)
print()


# ============================================================================
# 13. LOW-DEGREE INVARIANT SEARCH
# ============================================================================

print("[13] LOW-DEGREE OPPOSITE-SHIFT INVARIANT SEARCH")
print("-" * 78)

feature_map = {
    "sum": sum_opposite,
    "difference": difference_opposite,
    "product": product_opposite,
    "square_sum": sp.expand(M_plus**2 + M_minus**2),
    "square_difference": sp.expand(M_plus**2 - M_minus**2),
}

for name, expr in feature_map.items():
    print(f"  {name}:")
    print("    =", fact(expr))
    print("    contains S =", expr.has(S))
    print()


# ============================================================================
# 14. DISCRIMINANT-LIKE COMBINATIONS
# ============================================================================

print("[14] DISCRIMINANT-LIKE COMBINATIONS")
print("-" * 78)

disc_candidates = {
    "((N+x^2)^2-M(+x)M(-x))/x^2":
        sp.cancel(
            (
                (N + x**2)**2
                - product_opposite
            )
            / x**2
        ),

    "((M(+x)-M(-x))/(2x))^2":
        sp.cancel(
            (difference_opposite / (2*x))**2
        ),

    "((M(+x)-M(-x))/(2x))^2 - 4N":
        sp.cancel(
            (difference_opposite / (2*x))**2
            - 4*N
        ),

    "S^2 - 4N":
        delta_NS,
}

for name, expr in disc_candidates.items():
    print(f"  {name} =")
    print("   ", fact(expr))

print()


# ============================================================================
# 15. RELATION TO H2(x)
# ============================================================================

print("[15] CONNECTION TO H2(x)")
print("-" * 78)

H2 = sp.expand(
    6*N - S**2 + S*x
)

Mx = M_NS(x)

H2_candidate = sp.expand(
    Mx + N - delta_NS - x**2
)

H2_difference = sp.expand(
    H2 - H2_candidate
)

print(
    "  H2(x) =",
    fact(H2),
)
print(
    "  M_x + N - Delta - x^2 =",
    fact(H2_candidate),
)
print(
    "  difference =",
    fact(H2_difference),
)
print(
    "  PASS =",
    zero(H2_difference),
)
print()


# ============================================================================
# 16. COMPLETE ALGEBRAIC RECOVERY CHAIN
# ============================================================================

print("[16] COMPLETE RECOVERY CHAIN")
print("-" * 78)

recovered_S_chain = sp.cancel(
    (M1 - Mm1) / 2
)

recovered_delta_chain = sp.expand(
    (N + 1)**2 - M1*Mm1
)

quadratic = sp.expand(
    z**2
    - recovered_S_chain*z
    + N
)

quadratic_target = sp.expand(
    z**2 - S*z + N
)

quadratic_difference = sp.expand(
    quadratic - quadratic_target
)

print(
    "  recovered S =",
    fact(recovered_S_chain),
)
print(
    "  recovered Delta =",
    fact(recovered_delta_chain),
)
print(
    "  recovered quadratic =",
    fact(quadratic),
)
print(
    "  quadratic difference =",
    fact(quadratic_difference),
)
print(
    "  S PASS =",
    zero(recovered_S_chain - S),
)
print(
    "  Delta PASS =",
    zero(recovered_delta_chain - delta_NS),
)
print(
    "  quadratic PASS =",
    zero(quadratic_difference),
)
print()


# ============================================================================
# 17. NUMERICAL AUDIT
# ============================================================================

print("[17] NUMERICAL PRIME-PAIR AUDIT")
print("-" * 78)

instances = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]

numeric_failures = 0

for pp, qq in instances:
    nn = pp * qq
    ss = pp + qq
    true_delta = (pp - qq)**2

    for xx in [1, 2, 3, -1, -2]:

        mp_num = (pp + xx) * (qq + xx)
        mm_num = (pp - xx) * (qq - xx)

        numerator = (
            (nn + xx**2)**2
            - mp_num * mm_num
        )

        if xx == 0:
            numeric_failures += 1
            recovered_s2_num = None
            recovered_delta_num = None
            passed = False

        else:
            if numerator % (xx**2) != 0:
                numeric_failures += 1
                recovered_s2_num = None
                recovered_delta_num = None
                passed = False
            else:
                recovered_s2_num = numerator // (xx**2)
                recovered_delta_num = recovered_s2_num - 4*nn

                passed = (
                    recovered_s2_num == ss**2
                    and recovered_delta_num == true_delta
                )

                if not passed:
                    numeric_failures += 1

        print(
            f"  ({pp},{qq}), x={xx}: "
            f"S^2={recovered_s2_num}, "
            f"Delta={recovered_delta_num}, "
            f"PASS={passed}"
        )

print()
print("  NUMERICAL FAILURES =", numeric_failures)
print()


# ============================================================================
# 18. PRIME-PAIR ROOT RECOVERY
# ============================================================================

print("[18] PRIME-PAIR ROOT RECOVERY FROM M(+1), M(-1)")
print("-" * 78)

root_failures = 0

for pp, qq in instances:

    nn = pp * qq

    m1_num = (pp + 1) * (qq + 1)
    mm1_num = (pp - 1) * (qq - 1)

    recovered_s_num = (m1_num - mm1_num) // 2

    disc_num = (
        recovered_s_num**2
        - 4*nn
    )

    root_discriminant = sp.sqrt(disc_num)

    if not root_discriminant.is_Integer:
        root_failures += 1
        print(
            f"  ({pp},{qq}) discriminant not square: FAIL"
        )
        continue

    d_num = int(root_discriminant)

    r1 = (recovered_s_num - d_num) // 2
    r2 = (recovered_s_num + d_num) // 2

    passed = {r1, r2} == {pp, qq}

    if not passed:
        root_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"recovered roots={{{r1},{r2}}} "
        f"PASS={passed}"
    )

print()
print("  ROOT RECOVERY FAILURES =", root_failures)
print()


# ============================================================================
# 19. INFORMATION-MODEL BOUNDARY
# ============================================================================

print("[19] INFORMATION-MODEL BOUNDARY")
print("-" * 78)

print("The downstream structure is now exact:")
print()
print("  M_x")
print("    -> S")
print("    -> z^2-Sz+N")
print("    -> p,q")
print()

print("or, with opposite shifts:")
print()
print("  M(+x), M(-x)")
print("    -> S^2")
print("    -> Delta=(p-q)^2")
print("    -> S")
print("    -> p,q")
print()

print("For x=1:")
print("  S = (M(+1)-M(-1))/2")
print("  Delta = (N+1)^2-M(+1)M(-1)")
print()

print("Therefore the remaining research problem is upstream:")
print()
print("  N-only construction")
print("      -> M(+1), M(-1), or Delta")
print("      -> S")
print("      -> p,q")
print()


# ============================================================================
# 20. FINAL STATUS
# ============================================================================

disc_pass = zero(
    poly_discriminant - delta_NS
)

gap_pass = zero(
    delta_NS.subs(
        {
            N: N_pq,
            S: S_pq,
        }
    )
    - delta_pq
)

opposite_product_pass = zero(
    opposite_difference
)

s2_pass = zero(
    s2_difference
)

delta_pass = zero(
    delta_recovery_difference
)

general_delta_pass = zero(
    delta_a_difference
)

special_delta_pass = zero(
    delta_special_difference
)

h2_pass = zero(
    H2_difference
)

chain_s_pass = zero(
    recovered_S_chain - S
)

chain_delta_pass = zero(
    recovered_delta_chain - delta_NS
)

chain_quad_pass = zero(
    quadratic_difference
)

overall = (
    disc_pass
    and gap_pass
    and opposite_product_pass
    and s2_pass
    and delta_pass
    and general_delta_pass
    and special_delta_pass
    and h2_pass
    and chain_s_pass
    and chain_delta_pass
    and chain_quad_pass
    and numeric_failures == 0
    and root_failures == 0
)

print("[20] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  polynomial discriminant = Delta:",
    disc_pass,
)

print(
    "  Delta = (p-q)^2:",
    gap_pass,
)

print(
    "  opposite-shift product identity:",
    opposite_product_pass,
)

print(
    "  S^2 recovery:",
    s2_pass,
)

print(
    "  Delta recovery:",
    delta_pass,
)

print(
    "  generalized ±a Delta recovery:",
    general_delta_pass,
)

print(
    "  M(+1)M(-1) Delta identity:",
    special_delta_pass,
)

print(
    "  H2 = M_x + N - Delta - x^2:",
    h2_pass,
)

print(
    "  recovery-chain S:",
    chain_s_pass,
)

print(
    "  recovery-chain Delta:",
    chain_delta_pass,
)

print(
    "  recovery-chain quadratic:",
    chain_quad_pass,
)

print(
    "  numerical audit:",
    numeric_failures == 0,
)

print(
    "  root recovery audit:",
    root_failures == 0,
)

print(
    "  OVERALL EXACT AUDIT =",
    overall,
)

print()

print("MAIN RESULT")
print("-" * 78)

print("The shifted-product family is exactly")
print()
print("  M_x = (p+x)(q+x)")
print("      = x^2 + S*x + N.")
print()

print("Its discriminant as a polynomial in x is")
print()
print("  Delta = S^2 - 4N = (p-q)^2.")
print()

print("Opposite shifts satisfy")
print()
print("  M(+x)M(-x)")
print("    = (N+x^2)^2 - S^2*x^2.")
print()

print("Hence")
print()
print("  S^2")
print("    = ((N+x^2)^2-M(+x)M(-x))/x^2.")
print()

print("and")
print()
print("  Delta")
print("    = ((N+x^2)^2-M(+x)M(-x))/x^2 - 4N.")
print()

print("For x=1:")
print()
print("  Delta = (N+1)^2 - M(+1)M(-1).")
print()

print("The unresolved bridge is therefore sharply localized:")
print()
print("  N-only construction")
print("      -> Delta or shifted product(s)")
print("      -> S")
print("      -> z^2-Sz+N")
print("      -> p,q.")
print()

print("NEXT TARGET")
print("-" * 78)

print(
    "Search the existing N-only homogeneous-layer construction"
)
print(
    "for an exact object carrying the same algebraic signature as:"
)
print()
print("  (p-q)^2")
print("  (p+1)(q+1)")
print("  (p-1)(q-1)")
print()
print(
    "or their products/symmetric combinations."
)
print()
print(
    "Do not spend another experiment expanding the downstream"
)
print(
    "shifted kernel. That side is now algebraically exhausted."
)

print()
print("=" * 78)
print("EXPERIMENT 500 FINISHED")
print("=" * 78)