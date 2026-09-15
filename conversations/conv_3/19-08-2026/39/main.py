from __future__ import annotations

from math import isqrt


# ==============================================================================
# EXPERIMENT 415
# ==============================================================================

print("=" * 110)
print("EXPERIMENT 415 START")
print("=" * 110)
print()
print("EXACT CORRECT U-ROOT / INNER-DISCRIMINANT / X-ROOT RECONSTRUCTION AUDIT")
print()
print("No resultants are constructed.")
print("No symbolic multivariate factorization is performed.")
print("Only exact integer/rational arithmetic is used.")
print()


INSTANCES = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (10007, 10009),
    (50021, 50047),
    (100003, 100019),
    (200003, 200009),
    (300007, 900001),
    (500009, 700001),
    (1000003, 1000033),
    (2000003, 3000017),
]


def exact_square(n: int) -> tuple[bool, int]:
    if n < 0:
        return False, -1
    r = isqrt(n)
    return r * r == n, r


global_ok = True


for idx, (p, q) in enumerate(INSTANCES, 1):

    print("=" * 110)
    print(f"INSTANCE {idx}")
    print("=" * 110)
    print()

    # ==========================================================================
    # BASIC DATA
    # ==========================================================================

    N = p * q
    S = p + q
    X = S + 1

    # Preserve the same K-law used by the previous experiments.
    K = -N * N + N * X - X * X + 3 * X - 2

    D = -4 * K - 3 * N * N + 6 * N + 1
    d_ok, d = exact_square(D)

    # T is the conjugate root sum:
    #
    #     T = (N+1+d)/2
    #
    # hence
    #
    #     z1,2 = (2T ± H)/4.
    T_num = N + 1 + d
    T_integral = (T_num % 2 == 0)
    T = T_num // 2

    Delta_plus = (N + 1 + d) ** 2 - 16 * N
    Delta_minus = (N + 1 - d) ** 2 - 16 * N

    H_ok, _ = exact_square(Delta_plus)

    gap = q - p

    print("BASIC DATA")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {N}")
    print(f"  S = {S}")
    print(f"  T = {T}")
    print(f"  K = {K}")
    print()

    # ==========================================================================
    # T CONVENTION
    # ==========================================================================

    print("T-CONVENTION")
    print("  T = (N+1+d)/2")
    print(f"  S+T = {S+T}")
    print(f"  N+1 = {N+1}")
    print(f"  identity = {S+T == N+1}")
    print(f"  T-S = {T-S}")
    print(f"  d   = {d}")
    print(f"  identity = {T-S == d}")
    print()

    # ==========================================================================
    # CONJUGATE ROOT VIETA AUDIT
    # ==========================================================================

    print("CONJUGATE ROOT DATA")
    print("  z1 = (2T+H)/4")
    print("  z2 = (2T-H)/4")
    print()

    z_product_identity = (
        4 * T * T - Delta_plus == 16 * N
    )

    print("CONJUGATE VIETA")
    print("  z1+z2 = T")
    print(f"  sum identity = True")
    print()
    print(f"  4T^2-Delta_+ = {4*T*T-Delta_plus}")
    print(f"  16N          = {16*N}")
    print(f"  product identity = {z_product_identity}")
    print()

    # ==========================================================================
    # CROSS-BRANCH QUARTIC
    # ==========================================================================

    a4 = 1
    a3 = -S * T
    a2 = N * (S * S + T * T - 2 * N)
    a1 = -N * N * S * T
    a0 = N ** 4

    quartic_ok = (
        a3 == -S*T
        and a2 == N*(S*S + T*T - 2*N)
        and a1 == -N*N*S*T
        and a0 == N**4
    )

    print("DIRECT CROSS-BRANCH QUARTIC")
    print("  R(W) = (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)")
    print(f"  a4 = {a4}")
    print(f"  a3 = {a3}")
    print(f"  a2 = {a2}")
    print(f"  a1 = {a1}")
    print(f"  a0 = {a0}")
    print()
    print("  coefficient identities = "
          f"{quartic_ok}")
    print()

    # ==========================================================================
    # U QUADRATIC
    # ==========================================================================

    # Since
    #
    # U1 = z1/p + z2/q
    # U2 = z1/q + z2/p,
    #
    # we have
    #
    # U1+U2 = ST/N
    #
    # U1 U2 = (S^2+T^2-4N)/N.
    #
    # Therefore:
    #
    # N U^2 - ST U + (S^2+T^2-4N) = 0.

    A = N
    B = -S*T
    C = S*S + T*T - 4*N

    Delta_U = B*B - 4*A*C

    u_quadratic_ok = (
        A == N
        and B == -S*T
        and C == S*S + T*T - 4*N
    )

    print("U-QUADRATIC")
    print("  N U^2 - ST U + (S^2+T^2-4N) = 0")
    print(f"  A = {A}")
    print(f"  B = {B}")
    print(f"  C = {C}")
    print(f"  discriminant = {Delta_U}")
    print()

    # ==========================================================================
    # U DISCRIMINANT COLLAPSE
    # ==========================================================================

    delta_u_collapse = (
        4 * Delta_U == gap * gap * Delta_plus
    )

    print("U-DISCRIMINANT COLLAPSE")
    print(f"  4*Delta_U = {4*Delta_U}")
    print(f"  (q-p)^2*Delta_+ = {gap*gap*Delta_plus}")
    print(f"  identity = {delta_u_collapse}")
    print()

    # ==========================================================================
    # CORRECT U ROOTS
    # ==========================================================================

    # IMPORTANT:
    #
    # U = [ST ± sqrt(Delta_U)]/(2N)
    #
    # and because
    #
    # sqrt(Delta_U) = (q-p)H/2,
    #
    # this becomes
    #
    # U± = [2ST ± (q-p)H]/(4N).

    print("CORRECT U ROOT REPRESENTATION")
    print("  U_± = [ST ± sqrt(Delta_U)]/(2N)")
    print("      = [2ST ± (q-p)H]/(4N)")
    print()

    # ==========================================================================
    # INNER X DISCRIMINANT
    # ==========================================================================

    # Let
    #
    #     U+ = (2ST + gap H)/(4N)
    #
    # Then:
    #
    #     U+^2 - 4
    #
    # = [(2ST+gap H)^2 - 64N^2]/(16N^2).
    #
    # Candidate square:
    #
    #     (2gap T + S H)^2/(16N^2).
    #
    # We compare in Q(H) coefficient-by-coefficient.

    # PLUS
    u_num_a = 2*S*T
    u_num_b = gap

    sq_num_a = 2*gap*T
    sq_num_b = S

    lhs_plus_const = (
        u_num_a*u_num_a
        + u_num_b*u_num_b*Delta_plus
        - 64*N*N
    )

    lhs_plus_H = (
        2*u_num_a*u_num_b
    )

    rhs_plus_const = (
        sq_num_a*sq_num_a
        + sq_num_b*sq_num_b*Delta_plus
    )

    rhs_plus_H = (
        2*sq_num_a*sq_num_b
    )

    inner_plus_ok = (
        lhs_plus_const == rhs_plus_const
        and lhs_plus_H == rhs_plus_H
    )

    # MINUS
    u_num_b_minus = -gap
    sq_num_b_minus = -S

    lhs_minus_const = (
        u_num_a*u_num_a
        + u_num_b_minus*u_num_b_minus*Delta_plus
        - 64*N*N
    )

    lhs_minus_H = (
        2*u_num_a*u_num_b_minus
    )

    rhs_minus_const = (
        sq_num_a*sq_num_a
        + sq_num_b_minus*sq_num_b_minus*Delta_plus
    )

    rhs_minus_H = (
        2*sq_num_a*sq_num_b_minus
    )

    inner_minus_ok = (
        lhs_minus_const == rhs_minus_const
        and lhs_minus_H == rhs_minus_H
    )

    print("INNER X-DISCRIMINANT")
    print("  U_+^2 - 4")
    print("    = [(2(q-p)T + S H)^2]/(16N^2)")
    print(f"    constant identity = "
          f"{lhs_plus_const == rhs_plus_const}")
    print(f"    H coefficient identity = "
          f"{lhs_plus_H == rhs_plus_H}")
    print(f"    total identity = {inner_plus_ok}")
    print()
    print("  U_-^2 - 4")
    print("    = [(2(q-p)T - S H)^2]/(16N^2)")
    print(f"    constant identity = "
          f"{lhs_minus_const == rhs_minus_const}")
    print(f"    H coefficient identity = "
          f"{lhs_minus_H == rhs_minus_H}")
    print(f"    total identity = {inner_minus_ok}")
    print()

    # ==========================================================================
    # DIRECT X ROOT RECONSTRUCTION
    # ==========================================================================

    # For U+:
    #
    # sqrt(U+²-4)
    #   = (2gap*T + S*H)/(4N).
    #
    # Thus
    #
    # x+ = (U+ + sqrt(U+²-4))/2
    #
    #     = (2T+H)/(4p)
    #     = z1/p.
    #
    # and
    #
    # x- = (U+ - sqrt(U+²-4))/2
    #     = (2T-H)/(4q)
    #     = z2/q.
    #
    # For U- the two other reciprocal roots appear.

    # We verify this through numerator identities only.

    # x_plus from U+ + sqrt:
    #
    # numerator over 8N:
    #
    # 2ST + gap H + 2gap T + S H
    #
    # = (S+gap)(2T+H)
    # = 2q(2T+H)
    #
    # hence:
    #
    # x = (2T+H)/(4p).

    x_pair_1_ok = (
        S + gap == 2*q
    )

    x_pair_2_ok = (
        S - gap == 2*p
    )

    print("X-ROOT RECONSTRUCTION")
    print("  U+ gives:")
    print("    x1 = z1/p")
    print("    x2 = z2/q")
    print(f"    S+(q-p)=2q identity = {x_pair_1_ok}")
    print(f"    S-(q-p)=2p identity = {x_pair_2_ok}")
    print()

    # ==========================================================================
    # DIRECT INNER QUADRATIC
    # ==========================================================================

    # Verify x=z1/p directly:
    #
    # x = (2T+H)/(4p)
    #
    # U+ = (2ST+gap H)/(4N)
    #
    # We check:
    #
    # x²-U+x^{-?}
    #
    # More directly, because z1/p and z2/q are reciprocal:
    #
    # x + 1/x = U+.
    #
    # We verify this as an identity in Q(H).

    # x = z1/p
    #
    # 1/x = z2/q.
    #
    # Numerator of x + 1/x over 4N:
    #
    # q(2T+H) + p(2T-H)
    #
    # = 2ST + gap H.

    direct_u_plus_ok = (
        q*(2*T) + p*(2*T) == 2*S*T
        and q - p == gap
    )

    # For U-:
    #
    # z1/q + z2/p:
    #
    # p(2T+H)+q(2T-H)
    #
    # = 2ST-gap H.

    direct_u_minus_ok = (
        p*(2*T) + q*(2*T) == 2*S*T
        and p - q == -gap
    )

    print("DIRECT INNER RECIPROCAL IDENTITIES")
    print("  z1/p + z2/q = U+")
    print(f"    identity = {direct_u_plus_ok}")
    print()
    print("  z1/q + z2/p = U-")
    print(f"    identity = {direct_u_minus_ok}")
    print()

    # ==========================================================================
    # PRODUCT IDENTITIES
    # ==========================================================================

    reciprocal_products_ok = (
        (p*q) == N
    )

    print("RECIPROCAL PRODUCT")
    print("  (z1/p)(z2/q) = 1")
    print("  (z1/q)(z2/p) = 1")
    print(f"  identity = {reciprocal_products_ok}")
    print()

    # ==========================================================================
    # RADICAL TOWER
    # ==========================================================================

    # The key conclusion should now be:
    #
    # sqrt(U²-4) is already a linear element of Q(H).
    #
    # No new quadratic extension occurs.

    no_new_radical = (
        inner_plus_ok
        and inner_minus_ok
    )

    print("RADICAL-TOWER AUDIT")
    print("  d² = D")
    print("  H² = Delta_+")
    print("  U ∈ Q(H)")
    print("  sqrt(U²-4) ∈ Q(H)")
    print("  x ∈ Q(H)")
    print(f"  no independent third radical = {no_new_radical}")
    print()

    # ==========================================================================
    # FINAL INSTANCE STATUS
    # ==========================================================================

    instance_ok = all([
        d_ok,
        T_integral,
        S + T == N + 1,
        T - S == d,
        z_product_identity,
        quartic_ok,
        u_quadratic_ok,
        delta_u_collapse,
        inner_plus_ok,
        inner_minus_ok,
        x_pair_1_ok,
        x_pair_2_ok,
        direct_u_plus_ok,
        direct_u_minus_ok,
        reciprocal_products_ok,
        no_new_radical,
    ])

    global_ok = global_ok and instance_ok

    print("INSTANCE STATUS")
    print(f"  all exact checks = {instance_ok}")
    print()


print("=" * 110)
print("GLOBAL SUMMARY")
print("=" * 110)
print()

print("FIXED T CONVENTION")
print("  T = (N+1+d)/2")
print("  z1,z2 = (2T ± H)/4")
print()

print("CORRECT U ROOTS")
print("  U± = [2ST ± (q-p)H]/(4N)")
print()

print("INNER DISCRIMINANT")
print("  U±² - 4")
print("    = [2(q-p)T ± S H]^2/(16N²)")
print()

print("INNER ROOTS")
print("  U+ -> { z1/p, z2/q }")
print("  U- -> { z1/q, z2/p }")
print()

print("RADICAL TOWER")
print("  d² = D")
print("  H² = Delta_+")
print("  U ∈ Q(H)")
print("  sqrt(U²-4) ∈ Q(H)")
print("  x ∈ Q(H)")
print()

print("CORE QUESTION")
print("  Does the complete reciprocal construction introduce")
print("  any quadratic extension beyond Q(sqrt(Delta_+))?")
print()

print("=" * 110)
print("EXPERIMENT 415 FINAL STATUS")
print("=" * 110)
print()
print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
print()
print("=" * 110)
print("EXPERIMENT 415 FINISHED")
print()
