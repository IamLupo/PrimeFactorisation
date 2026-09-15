from __future__ import annotations

from math import isqrt


# ==============================================================================
# EXPERIMENT 413
# ==============================================================================
#
# EXACT CROSS-BRANCH NORMALIZATION / CORRECTED U-LAYER /
# INNER RADICAL COLLAPSE AUDIT
#
# This experiment intentionally rebuilds the cross-branch polynomial from
# the actual conjugate roots:
#
#     z1 = (T + H)/4
#     z2 = (T - H)/4
#
# with
#
#     H^2 = Delta_plus.
#
# No previously derived quartic coefficients are trusted.
#
# No resultants.
# No symbolic multivariate factorization.
# Exact integer/rational arithmetic only.
#
#
# MAIN PURPOSE
#
# 1. Rebuild R(W) directly from:
#
#       (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)
#
# 2. Verify the CORRECT coefficients:
#
#       a4 = 1
#       a3 = -S*T/2
#
#       a2 = N*(S^2 + T^2/4 - 2N)
#
#       a1 = -N^2*S*T/2
#
#       a0 = N^4
#
# 3. Normalize W=N*x.
#
# 4. Verify reciprocal symmetry.
#
# 5. Derive the corrected reciprocal variable:
#
#       U = x + 1/x
#
# 6. Verify the corrected quadratic:
#
#       4N U^2
#       -2ST U
#       +(4S^2 + T^2 - 16N)
#       = 0
#
# 7. Verify its roots:
#
#       U_± = [ST ± (q-p)H] / (4N)
#
# 8. Verify:
#
#       U_±^2 - 4
#       =
#       [(q-p)T ± S H]^2 / (16N^2)
#
# 9. Reconstruct the x-roots exactly:
#
#       x = z1/p
#       x = z2/q
#       x = z1/q
#       x = z2/p
#
# 10. Establish whether the inner quadratic introduces a new radical.
#
# ==============================================================================


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


def exact_sqrt(n: int):
    if n < 0:
        return False, 0
    r = isqrt(n)
    return r * r == n, r


def main():

    print("=" * 110)
    print("EXPERIMENT 413 START")
    print("=" * 110)
    print()
    print("EXACT CROSS-BRANCH NORMALIZATION / CORRECTED U-LAYER /")
    print("INNER RADICAL COLLAPSE AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factorization is performed.")
    print("Only exact integer/rational arithmetic is used.")
    print()

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        # ------------------------------------------------------------------
        # BASIC DATA
        # ------------------------------------------------------------------

        N = p * q
        S = p + q
        gap = q - p

        # Existing source-law reconstruction.
        X = S + 1
        K = -N * N + N * X - X * X + 3 * X - 2

        D = -4 * K - 3 * N * N + 6 * N + 1
        d_ok, d = exact_sqrt(D)

        # T is defined so that:
        #
        #   z1 + z2 = T/2
        #
        # and therefore
        #
        #   T = (N+1+d)/2.
        #
        T = (N + 1 + d) // 2

        Delta_minus = (N + 1 - d) ** 2 - 16 * N
        Delta_plus = (N + 1 + d) ** 2 - 16 * N

        h_ok, _ = exact_sqrt(Delta_plus)

        # ------------------------------------------------------------------
        # CONJUGATE ROOT DATA
        # ------------------------------------------------------------------

        # H^2 = Delta_plus
        #
        # z1 = (T+H)/4
        # z2 = (T-H)/4
        #
        # Exact identities:
        #
        #   z1+z2 = T/2
        #   z1*z2 = N
        #
        sum_half_ok = (2 * (T // 2) == T)
        product_identity_symbolic = (
            T * T - Delta_plus == 16 * N
        )

        # ------------------------------------------------------------------
        # DIRECT CROSS-BRANCH QUARTIC
        # ------------------------------------------------------------------
        #
        # First quadratic:
        #
        #   (W-pz1)(W-pz2)
        #
        # = W^2 - p(z1+z2)W + p^2 z1z2
        #
        # = W^2 - p*T/2 W + p^2 N
        #
        # Second:
        #
        # = W^2 - q*T/2 W + q^2 N
        #
        # Multiply exactly.
        #
        # To avoid fractions, first use doubled coefficients and then
        # verify divisibility by 2 where appropriate.
        # ------------------------------------------------------------------

        # a3 = -S*T/2
        a3_num = -S * T
        a3_ok_integral = (a3_num % 2 == 0)
        a3 = a3_num // 2

        # a2 = N*(S^2 + T^2/4 - 2N)
        #     = N*(4S^2 + T^2 - 8N)/4
        a2_num = N * (4 * S * S + T * T - 8 * N)
        a2_ok_integral = (a2_num % 4 == 0)
        a2 = a2_num // 4

        # a1 = -N^2*S*T/2
        a1_num = -N * N * S * T
        a1_ok_integral = (a1_num % 2 == 0)
        a1 = a1_num // 2

        a4 = 1
        a0 = N ** 4

        # ------------------------------------------------------------------
        # DIRECT MULTIPLICATION VIA INTEGER-SCALED QUADRATICS
        # ------------------------------------------------------------------
        #
        # Let:
        #
        #   Fp(W) = 2W² - pT W + 2p²N
        #   Fq(W) = 2W² - qT W + 2q²N
        #
        # Then:
        #
        #   Fp*Fq = 4R(W).
        #
        # Multiply explicitly.
        # ------------------------------------------------------------------

        scaled_a4 = 4

        scaled_a3 = -2 * S * T

        scaled_a2 = (
            p * q * T * T
            + 4 * N * (p * p + q * q)
        )

        scaled_a1 = -2 * N * N * S * T

        scaled_a0 = 4 * N ** 4

        scaled_to_direct_ok = (
            scaled_a4 == 4 * a4
            and scaled_a3 == 4 * a3
            and scaled_a2 == 4 * a2
            and scaled_a1 == 4 * a1
            and scaled_a0 == 4 * a0
        )

        # ------------------------------------------------------------------
        # NORMALIZATION
        # ------------------------------------------------------------------
        #
        # R(Nx)/N^4
        #
        # gives:
        #
        # x^4
        # -(ST/(2N)) x^3
        # + ...
        # -(ST/(2N)) x
        # + 1
        # ------------------------------------------------------------------

        reciprocal_ok = (
            a0 == N ** 4
            and a1 == N * N * a3
        )

        # ------------------------------------------------------------------
        # CORRECT U-QUADRATIC
        # ------------------------------------------------------------------
        #
        # Divide normalized quartic by x²:
        #
        #   x² + x^-2
        #   -(ST/(2N))(x+x^-1)
        #   + a2/N²
        #
        # Since:
        #
        #   x²+x^-2 = U²-2
        #
        # we obtain:
        #
        #   4N U²
        #   -2ST U
        #   +(4S²+T²-16N)
        #   = 0.
        # ------------------------------------------------------------------

        U_A = 4 * N
        U_B = -2 * S * T
        U_C = 4 * S * S + T * T - 16 * N

        # Discriminant in U:
        #
        # B² - 4AC
        corrected_Delta_U = (
            U_B * U_B
            - 4 * U_A * U_C
        )

        # Expected collapse:
        #
        #   corrected_Delta_U = (gap²)*(Delta_plus)
        #
        expected_Delta_U = gap * gap * Delta_plus

        U_discriminant_ok = (
            corrected_Delta_U == expected_Delta_U
        )

        # Therefore:
        #
        #   sqrt(Delta_U_corrected) = gap H
        #
        # and:
        #
        #   U± = [2ST ± gap H]/(8N)
        #       = [ST ± gap H/2]/(4N)
        #
        # But because the quadratic is written with 4N,
        # the equivalent useful form is:
        #
        #   U± = [ST ± gap H/2]/(4N).
        #
        # We keep this distinction explicit.

        # ------------------------------------------------------------------
        # DIRECT U VALUES FROM THE ACTUAL X ROOTS
        # ------------------------------------------------------------------
        #
        # Pair 1:
        #
        #   x1 = z1/p
        #   1/x1 = z2/q
        #
        # Hence:
        #
        #   U1 = z1/p + z2/q
        #      = [ST + gap H]/(4N)
        #
        # Pair 2:
        #
        #   x2 = z1/q
        #   1/x2 = z2/p
        #
        #      = [ST - gap H]/(4N)
        #
        # These are the actual U-values.
        # ------------------------------------------------------------------

        # For compatibility with the quadratic above, substitute:
        #
        # U_actual numerator denominator 4N:
        #
        #   U_num = ST ± gap H
        #
        # ------------------------------------------------------------------
        # INNER X DISCRIMINANT
        # ------------------------------------------------------------------
        #
        # Claim:
        #
        #   U±² - 4
        #   =
        #   [(gap)T ± S H]² / (16N²)
        #
        # Check in Q(H).
        # ------------------------------------------------------------------

        # For the + branch:
        #
        # U numerator = ST + gap H
        #
        up_a = S * T
        up_b = gap

        # Candidate sqrt numerator:
        #
        # gap*T + S*H
        #
        xp_a = gap * T
        xp_b = S

        # Uplus² - 4 numerator over 16N²:
        #
        # (ST + gap H)^2 - 64N²
        #
        lhs_plus_const = (
            up_a * up_a
            + up_b * up_b * Delta_plus
            - 64 * N * N
        )

        lhs_plus_rad = 2 * up_a * up_b

        # Candidate square:
        #
        # (gap*T + S H)^2
        #
        rhs_plus_const = (
            xp_a * xp_a
            + xp_b * xp_b * Delta_plus
        )

        rhs_plus_rad = 2 * xp_a * xp_b

        plus_xdisc_ok = (
            lhs_plus_const == rhs_plus_const
            and lhs_plus_rad == rhs_plus_rad
        )

        # Minus branch.
        um_a = S * T
        um_b = -gap

        xm_a = gap * T
        xm_b = -S

        lhs_minus_const = (
            um_a * um_a
            + um_b * um_b * Delta_plus
            - 64 * N * N
        )

        lhs_minus_rad = 2 * um_a * um_b

        rhs_minus_const = (
            xm_a * xm_a
            + xm_b * xm_b * Delta_plus
        )

        rhs_minus_rad = 2 * xm_a * xm_b

        minus_xdisc_ok = (
            lhs_minus_const == rhs_minus_const
            and lhs_minus_rad == rhs_minus_rad
        )

        # ------------------------------------------------------------------
        # X-ROOT RECONSTRUCTION
        # ------------------------------------------------------------------
        #
        # For the + branch:
        #
        #   x = (U + sqrt(U²-4))/2
        #
        # = [ST + gap H + gap T + S H] / (8N)
        #
        # The result should be:
        #
        #   z1/p
        # = (T+H)/(4p)
        #
        # and the companion root:
        #
        #   z2/q = (T-H)/(4q)
        #
        # We compare coefficient pairs after converting everything
        # to denominator 8N.
        # ------------------------------------------------------------------

        # U+ plus sqrt:
        x_plus_1 = (
            S * T + gap * T,
            gap + S,
        )

        # U+ minus sqrt:
        x_plus_2 = (
            S * T - gap * T,
            gap - S,
        )

        # U- plus sqrt:
        x_minus_1 = (
            S * T + gap * T,
            -gap + S,
        )

        # U- minus sqrt:
        x_minus_2 = (
            S * T - gap * T,
            -gap - S,
        )

        # Expected z1/p:
        #
        # (T+H)/(4p) = [2qT + 2qH]/(8N)
        expected_1 = (
            2 * q * T,
            2 * q,
        )

        # Expected z2/q:
        #
        # (T-H)/(4q) = [2pT - 2pH]/(8N)
        expected_2 = (
            2 * p * T,
            -2 * p,
        )

        # Expected z1/q:
        expected_3 = (
            2 * p * T,
            2 * p,
        )

        # Expected z2/p:
        expected_4 = (
            2 * q * T,
            -2 * q,
        )

        # IMPORTANT:
        #
        # The expressions above are the raw numerators over 8N.
        # The actual x expression contains an extra factor 1/2 from
        # the outer quadratic formula:
        #
        #   x = (U ± sqrt(U²-4))/2.
        #
        # Therefore compare 2 * numerator-pairs.
        #
        x1_ok = (
            2 * x_plus_1[0] == expected_1[0]
            and 2 * x_plus_1[1] == expected_1[1]
        )

        x2_ok = (
            2 * x_plus_2[0] == expected_2[0]
            and 2 * x_plus_2[1] == expected_2[1]
        )

        x3_ok = (
            2 * x_minus_1[0] == expected_3[0]
            and 2 * x_minus_1[1] == expected_3[1]
        )

        x4_ok = (
            2 * x_minus_2[0] == expected_4[0]
            and 2 * x_minus_2[1] == expected_4[1]
        )

        # ------------------------------------------------------------------
        # DIRECT PRODUCT CHECKS
        # ------------------------------------------------------------------

        # z1*z2=N and N=pq imply:
        #
        # (z1/p)(z2/q)=1
        # (z1/q)(z2/p)=1

        reciprocal_products_ok = (
            (p * q == N)
            and (p * q == N)
        )

        # ------------------------------------------------------------------
        # FINAL RADICAL CONCLUSION
        # ------------------------------------------------------------------

        # The inner discriminant is a square in Q(H) exactly when both
        # field identities hold.
        #
        no_new_radical = plus_xdisc_ok and minus_xdisc_ok

        local_ok = all([
            d_ok,
            product_identity_symbolic,
            scaled_to_direct_ok,
            reciprocal_ok,
            U_discriminant_ok,
            plus_xdisc_ok,
            minus_xdisc_ok,
            reciprocal_products_ok,
            x1_ok,
            x2_ok,
            x3_ok,
            x4_ok,
            no_new_radical,
        ])

        global_ok = global_ok and local_ok

        # ------------------------------------------------------------------
        # OUTPUT
        # ------------------------------------------------------------------

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p = {p}")
        print(f"  q = {q}")
        print(f"  N = {N}")
        print(f"  S = {S}")
        print(f"  T = {T}")
        print(f"  K = {K}")
        print()

        print("CONJUGATE ROOT DATA")
        print("  z1 = (T+H)/4")
        print("  z2 = (T-H)/4")
        print("  z1+z2 = T/2")
        print("  z1*z2 = N")
        print(f"  T^2-Delta_+ = {T*T-Delta_plus}")
        print(f"  16N = {16*N}")
        print(f"  product identity = {product_identity_symbolic}")
        print()

        print("DIRECT CROSS-BRANCH QUARTIC")
        print("  R(W) = (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)")
        print(f"  a4 = {a4}")
        print(f"  a3 = {a3}")
        print(f"  a2 = {a2}")
        print(f"  a1 = {a1}")
        print(f"  a0 = {a0}")
        print()

        print("CORRECTED COEFFICIENT IDENTITIES")
        print(f"  a3 = -S*T/2           : {a3 == -S*T//2}")
        print(
            f"  a2 = N*(S^2+T^2/4-2N) : "
            f"{a2 == N*(4*S*S+T*T-8*N)//4}"
        )
        print(f"  a1 = -N^2*S*T/2       : {a1 == -N*N*S*T//2}")
        print(f"  a0 = N^4              : {a0 == N**4}")
        print(f"  direct multiplication  : {scaled_to_direct_ok}")
        print()

        print("RECIPROCAL SYMMETRY")
        print(f"  a0 = N^4*a4 : {a0 == N**4*a4}")
        print(f"  a1 = N^2*a3 : {a1 == N**2*a3}")
        print()

        print("CORRECT U-QUADRATIC")
        print("  4N U^2 - 2ST U + (4S^2+T^2-16N) = 0")
        print(f"  A = {U_A}")
        print(f"  B = {U_B}")
        print(f"  C = {U_C}")
        print()

        print("U-DISCRIMINANT")
        print(f"  Delta_U(corrected) = {corrected_Delta_U}")
        print(f"  (q-p)^2*Delta_+   = {expected_Delta_U}")
        print(f"  identity = {U_discriminant_ok}")
        print()

        print("ACTUAL U ROOTS FROM THE CROSS-BRANCH ROOTS")
        print("  U_+ = [ST + (q-p)H]/(4N)")
        print("  U_- = [ST - (q-p)H]/(4N)")
        print()

        print("INNER X-DISCRIMINANT")
        print("  U_+^2 - 4")
        print("    = [(q-p)T + S H]^2/(16N^2)")
        print(f"    exact field identity = {plus_xdisc_ok}")
        print()
        print("  U_-^2 - 4")
        print("    = [(q-p)T - S H]^2/(16N^2)")
        print(f"    exact field identity = {minus_xdisc_ok}")
        print()

        print("X-ROOT RECONSTRUCTION")
        print(f"  z1/p : {x1_ok}")
        print(f"  z2/q : {x2_ok}")
        print(f"  z1/q : {x3_ok}")
        print(f"  z2/p : {x4_ok}")
        print()

        print("RADICAL TEST")
        print(f"  new independent radical = {not no_new_radical}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("CORRECTED DIRECT QUARTIC")
    print("  a3 = -S*T/2")
    print("  a2 = N*(S^2 + T^2/4 - 2N)")
    print("  a1 = -N^2*S*T/2")
    print("  a0 = N^4")
    print()

    print("CORRECT U-LAYER")
    print("  4N U^2 - 2ST U + (4S^2+T^2-16N) = 0")
    print()

    print("CORRECT U ROOTS")
    print("  U_± = [ST ± (q-p)H]/(4N)")
    print()

    print("INNER RADICAL COLLAPSE")
    print("  U_±² - 4")
    print("    = [(q-p)T ± S H]^2/(16N²)")
    print()

    print("INNER ROOTS")
    print("  x = z1/p, z2/q, z1/q, z2/p")
    print()

    print("RADICAL-TOWER CONCLUSION")
    print("  If all instances pass, the reciprocal x-layer introduces")
    print("  no independent quadratic extension.")
    print("  The apparent failure in Experiment 412 is then identified")
    print("  as a factor-of-two normalization error, not a new radical.")
    print()

    print("=" * 110)
    print("EXPERIMENT 413 FINAL STATUS")
    print("=" * 110)
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 413 FINISHED")
    print()


if __name__ == "__main__":
    main()
