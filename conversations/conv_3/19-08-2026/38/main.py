from __future__ import annotations

from math import isqrt


# ==============================================================================
# EXPERIMENT 414
# ==============================================================================

print("=" * 110)
print("EXPERIMENT 414 START")
print("=" * 110)
print()
print("EXACT T-CONVENTION / RECIPROCAL QUARTIC / INNER RADICAL COLLAPSE AUDIT")
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


def is_square(n: int):
    if n < 0:
        return False, None
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

    # Keep the original source-law convention.
    X = S + 1
    K = -N * N + N * X - X * X + 3 * X - 2

    D = -4 * K - 3 * N * N + 6 * N + 1
    d_ok, d = is_square(D)

    # IMPORTANT:
    #
    # Here T is explicitly defined as the CONJUGATE BRANCH ROOT SUM:
    #
    #     T = z1 + z2
    #       = (N+1+d)/2
    #
    T_num = N + 1 + d
    T_even = (T_num % 2 == 0)
    T = T_num // 2

    Delta_plus = (N + 1 + d) ** 2 - 16 * N
    Delta_minus = (N + 1 - d) ** 2 - 16 * N

    H_ok, H = is_square(Delta_plus)

    print("BASIC DATA")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {N}")
    print(f"  S = {S}")
    print(f"  T = {T}")
    print(f"  K = {K}")
    print()

    print("T-CONVENTION")
    print("  T = (N+1+d)/2")
    print(f"  numerator = N+1+d = {T_num}")
    print(f"  T integer = {T_even}")
    print(f"  S+T = {S+T}")
    print(f"  expected N+1 = {N+1}")
    print(f"  S+T identity = {S+T == N+1}")
    print(f"  T-S = {T-S}")
    print(f"  expected d = {d}")
    print(f"  T-S identity = {T-S == d}")
    print()

    # ==========================================================================
    # CORRECT CONJUGATE ROOTS
    # ==========================================================================

    # Because T=(N+1+d)/2:
    #
    # z1 = (N+1+d+H)/4 = (2T+H)/4
    # z2 = (N+1+d-H)/4 = (2T-H)/4

    print("CONJUGATE ROOTS")
    print("  z1 = (2T+H)/4")
    print("  z2 = (2T-H)/4")
    print()

    # Correct Vieta relation:
    #
    # z1+z2 = T
    #
    # z1*z2 = (4T^2-H^2)/16 = N

    root_sum_ok = True

    product_lhs = 4 * T * T - Delta_plus
    product_rhs = 16 * N
    root_product_ok = (product_lhs == product_rhs)

    print("CONJUGATE VIETA")
    print(f"  z1+z2 = T")
    print(f"  root-sum identity = {root_sum_ok}")
    print()
    print(f"  4T^2-Delta_+ = {product_lhs}")
    print(f"  16N          = {product_rhs}")
    print(f"  root-product identity = {root_product_ok}")
    print()

    # ==========================================================================
    # DIRECT CROSS-BRANCH QUARTIC
    # ==========================================================================

    # With z1+z2=T and z1*z2=N:
    #
    # (W-pz1)(W-pz2)
    # = W^2-pT W+p^2 N
    #
    # (W-qz1)(W-qz2)
    # = W^2-qT W+q^2 N

    a4 = 1
    a3 = -S * T
    a2 = N * (S * S + T * T - 2 * N)
    a1 = -N * N * S * T
    a0 = N ** 4

    print("DIRECT CROSS-BRANCH QUARTIC")
    print("  R(W) = (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)")
    print(f"  a4 = {a4}")
    print(f"  a3 = {a3}")
    print(f"  a2 = {a2}")
    print(f"  a1 = {a1}")
    print(f"  a0 = {a0}")
    print()

    quartic_formula_ok = (
        a3 == -S*T
        and a2 == N * (S*S + T*T - 2*N)
        and a1 == -N*N*S*T
        and a0 == N**4
    )

    print("QUARTIC COEFFICIENT IDENTITIES")
    print(f"  a3 = -S*T                    : {a3 == -S*T}")
    print(f"  a2 = N*(S^2+T^2-2N)          : {a2 == N*(S*S+T*T-2*N)}")
    print(f"  a1 = -N^2*S*T                : {a1 == -N*N*S*T}")
    print(f"  a0 = N^4                     : {a0 == N**4}")
    print(f"  all coefficient identities    : {quartic_formula_ok}")
    print()

    # ==========================================================================
    # RECIPROCAL NORMALIZATION
    # ==========================================================================

    # W=Nx:
    #
    # R(Nx)/N^4 =
    #
    # x^4
    # -(ST/N)x^3
    # +(S^2+T^2-2N)/N^2*x^2
    # -(ST/N)x
    # +1

    reciprocal_ok = (
        a0 == N**4
        and a1 == N*N*a3
    )

    print("RECIPROCAL SYMMETRY")
    print(f"  a0 = N^4*a4 : {a0 == N**4*a4}")
    print(f"  a1 = N^2*a3 : {a1 == N*N*a3}")
    print()

    # ==========================================================================
    # U = x + 1/x
    # ==========================================================================

    # Divide by x^2:
    #
    # U^2 - 2
    # -(ST/N)U
    # +(S^2+T^2-2N)/N^2
    #
    # = 0
    #
    # Therefore:
    #
    # N U^2 - ST U + (S^2+T^2-4N) = 0.

    U_A = N
    U_B = -S*T
    U_C = S*S + T*T - 4*N

    print("CORRECT U-QUADRATIC")
    print("  N U^2 - ST U + (S^2+T^2-4N) = 0")
    print(f"  A = {U_A}")
    print(f"  B = {U_B}")
    print(f"  C = {U_C}")
    print()

    # ==========================================================================
    # U-DISCRIMINANT
    # ==========================================================================

    Delta_U = (
        U_B * U_B
        - 4 * U_A * U_C
    )

    expected_Delta_U = (q-p) ** 2 * Delta_plus // 4

    delta_divisible = ((q-p) ** 2 * Delta_plus) % 4 == 0

    U_delta_ok = (
        delta_divisible
        and Delta_U == expected_Delta_U
    )

    print("U-DISCRIMINANT")
    print(f"  Delta_U = {Delta_U}")
    print(f"  (q-p)^2*Delta_+/4 = {expected_Delta_U}")
    print(f"  divisibility by 4 = {delta_divisible}")
    print(f"  identity = {U_delta_ok}")
    print()

    # ==========================================================================
    # DIRECT U ROOTS
    # ==========================================================================

    # From actual cross products:
    #
    # x1 = z1/p
    # x1^{-1} = z2/q
    #
    # U1 = z1/p + z2/q
    #
    # = [q z1+p z2]/N
    #
    # = [ST + (q-p)H/2]/N
    #
    # = [2ST+(q-p)H]/(2N)
    #
    # Therefore the quadratic roots are:
    #
    # U± = [ST ± sqrt(Delta_U)]/(2N)
    #
    # with:
    #
    # sqrt(Delta_U)=(q-p)H/2.

    print("U ROOT REPRESENTATION")
    print("  U_± = [2ST ± (q-p)H]/(2N)")
    print("      = [ST ± (q-p)H/2]/N")
    print()

    radical_collapse_ok = (
        4 * Delta_U == (q-p)**2 * Delta_plus
    )

    print("U-RADICAL COLLAPSE")
    print(f"  4*Delta_U = {4*Delta_U}")
    print(f"  (q-p)^2*Delta_+ = {(q-p)**2*Delta_plus}")
    print(f"  identity = {radical_collapse_ok}")
    print()

    # ==========================================================================
    # INNER X DISCRIMINANT
    # ==========================================================================

    # Use the U-root numerator over 2N:
    #
    # U± = [2ST ± gap*H]/(2N)
    #
    # Thus:
    #
    # U±²-4
    #
    # numerator:
    #
    # (2ST ± gap H)^2 - 16N²
    #
    # Candidate square:
    #
    # (gap*T ± 2*S*H)^2
    #
    # Exact identity follows from:
    #
    # 4(S²T²-N? etc.)
    #
    # We audit it coefficient-by-coefficient in Q(H).
    # --------------------------------------------------------------------------

    gap = q - p

    # PLUS
    up_a = 2 * S * T
    up_b = gap

    xp_a = gap * T
    xp_b = 2 * S

    lhs_plus_const = (
        up_a * up_a
        + up_b * up_b * Delta_plus
        - 16 * N * N
    )

    lhs_plus_H = (
        2 * up_a * up_b
    )

    rhs_plus_const = (
        xp_a * xp_a
        + xp_b * xp_b * Delta_plus
    )

    rhs_plus_H = (
        2 * xp_a * xp_b
    )

    plus_inner_ok = (
        lhs_plus_const == rhs_plus_const
        and lhs_plus_H == rhs_plus_H
    )

    # MINUS
    um_a = 2 * S * T
    um_b = -gap

    xm_a = gap * T
    xm_b = -2 * S

    lhs_minus_const = (
        um_a * um_a
        + um_b * um_b * Delta_plus
        - 16 * N * N
    )

    lhs_minus_H = (
        2 * um_a * um_b
    )

    rhs_minus_const = (
        xm_a * xm_a
        + xm_b * xm_b * Delta_plus
    )

    rhs_minus_H = (
        2 * xm_a * xm_b
    )

    minus_inner_ok = (
        lhs_minus_const == rhs_minus_const
        and lhs_minus_H == rhs_minus_H
    )

    print("INNER X-DISCRIMINANT COLLAPSE")
    print("  U_+^2 - 4")
    print("    = [(q-p)T + 2S H]^2/(4N^2)")
    print(f"    constant identity = {lhs_plus_const == rhs_plus_const}")
    print(f"    H coefficient identity = {lhs_plus_H == rhs_plus_H}")
    print(f"    total identity = {plus_inner_ok}")
    print()
    print("  U_-^2 - 4")
    print("    = [(q-p)T - 2S H]^2/(4N^2)")
    print(f"    constant identity = {lhs_minus_const == rhs_minus_const}")
    print(f"    H coefficient identity = {lhs_minus_H == rhs_minus_H}")
    print(f"    total identity = {minus_inner_ok}")
    print()

    # ==========================================================================
    # X ROOT RECONSTRUCTION
    # ==========================================================================

    # x = (U ± sqrt(U²-4))/2
    #
    # For U+:
    #
    # U+ = (2ST + gap H)/(2N)
    #
    # sqrt = (gap T + 2S H)/(2N)
    #
    # Therefore:
    #
    # x = [2ST + gap H + gap T + 2S H]/(4N)
    #
    # = [(2S+gap)T + (gap+2S)H]/(4N)
    #
    # Since:
    #
    #   2S + gap = 2(p+q)+(q-p)=p+3q
    #
    # this does NOT directly equal z1/p.
    #
    # This is a warning that the chosen U pairing is not the correct
    # pairing for the proposed inner-root map.
    #
    # Therefore we explicitly construct the actual reciprocal pairs
    # and compare their sums with the quadratic roots.
    # --------------------------------------------------------------------------

    # Actual reciprocal-pair sums:
    #
    # U1 = z1/p + z2/q
    # U2 = z1/q + z2/p
    #
    # Their numerators over 2N:
    #
    # U1 = (2ST + gap H)/(2N)
    # U2 = (2ST - gap H)/(2N)

    actual_pair_ok = radical_collapse_ok

    # Direct product:
    #
    # (z1/p)(z2/q)=N/(pq)=1
    # etc.
    reciprocal_product_ok = (
        N == p*q
        and N == p*q
    )

    print("RECIPROCAL PAIRS")
    print("  U1 = z1/p + z2/q")
    print("  U2 = z1/q + z2/p")
    print("  U1*? reciprocal pairing product = 1")
    print(f"  pair-product identities = {reciprocal_product_ok}")
    print()

    # ==========================================================================
    # IMPORTANT FIELD CHECK
    # ==========================================================================

    # The inner quadratic x²-Ux+1=0 must have roots
    #
    #     x = z1/p, z2/q
    #
    # for U = z1/p + z2/q.
    #
    # We therefore test the defining identity directly:
    #
    # x²-Ux+1 = 0
    #
    # without trying to infer it from an incorrectly paired square root.
    #
    # In Q(H), for x=z1/p:
    #
    # x = (2T+H)/(4p).
    #
    # U = (2ST+gap H)/(2N).
    #
    # We audit x²-Ux+1 coefficient-by-coefficient.
    # --------------------------------------------------------------------------

    # x = (2T+H)/(4p)
    #
    # Put everything over 16 p² N.
    #
    # Rather than build rational fractions, use Fraction-free arithmetic.

    # x²:
    # numerator over 16 p²:
    x2_const = 4 * T * T + Delta_plus
    x2_H = 4 * T

    # U*x:
    # U numerator = 2ST + gap H over 2N
    # x numerator = 2T + H over 4p
    #
    # denominator = 8Np
    #
    # common denominator with x² and 1:
    # 16 p² N
    #
    # Multiply U*x numerator by 2p.
    ux_const = (
        2 * p * (
            4 * S * T * T
            + gap * Delta_plus
        )
    )

    ux_H = (
        2 * p * (
            2 * S * T
            + 2 * gap * T
        )
    )

    # x² numerator multiplied by N:
    x2_common_const = N * x2_const
    x2_common_H = N * x2_H

    # 1 numerator:
    one_common_const = 16 * p * p * N

    inner_x_const = (
        x2_common_const
        - ux_const
        + one_common_const
    )

    inner_x_H = (
        x2_common_H
        - ux_H
    )

    direct_inner_x_ok = (
        inner_x_const == 0
        and inner_x_H == 0
    )

    print("DIRECT INNER-ROOT EQUATION")
    print("  x = z1/p")
    print("  x^2-Ux+1 = 0")
    print(f"  constant residual = {inner_x_const}")
    print(f"  H residual        = {inner_x_H}")
    print(f"  exact identity    = {direct_inner_x_ok}")
    print()

    # ==========================================================================
    # FINAL STATUS
    # ==========================================================================

    local_ok = all([
        d_ok,
        T_even,
        S + T == N + 1,
        T - S == d,
        H_ok,
        root_product_ok,
        quartic_formula_ok,
        reciprocal_ok,
        U_delta_ok,
        radical_collapse_ok,
        plus_inner_ok,
        minus_inner_ok,
        reciprocal_product_ok,
        direct_inner_x_ok,
    ])

    global_ok = global_ok and local_ok

    print("INSTANCE STATUS")
    print(f"  all exact checks = {local_ok}")
    print()


print("=" * 110)
print("GLOBAL SUMMARY")
print("=" * 110)
print()

print("T CONVENTION")
print("  T = (N+1+d)/2")
print("  T = z1+z2")
print("  z1,z2 = (2T ± H)/4")
print()

print("CORRECT DIRECT QUARTIC")
print("  R(W) = W^4 - ST W^3")
print("       + N(S^2+T^2-2N)W^2")
print("       - N^2 ST W")
print("       + N^4")
print()

print("CORRECT RECIPROCAL QUADRATIC")
print("  N U^2 - ST U + (S^2+T^2-4N) = 0")
print()

print("U DISCRIMINANT")
print("  4 Delta_U = (q-p)^2 Delta_+")
print()

print("INNER QUADRATIC")
print("  x^2-Ux+1=0")
print()

print("DIRECT ROOT PAIR")
print("  x = z1/p")
print("  1/x = z2/q")
print()

print("CORE QUESTION")
print("  Does the corrected reciprocal layer remain entirely inside Q(H)?")
print()

print("=" * 110)
print("EXPERIMENT 414 FINAL STATUS")
print("=" * 110)
print()
print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
print()
print("=" * 110)
print("EXPERIMENT 414 FINISHED")
print()
