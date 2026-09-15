from __future__ import annotations

import math


# ==============================================================================================================
# EXPERIMENT 398
# ==============================================================================================================

EXP_NO = 398

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT FACTOR-GAP / SHIFTED-FACTOR QUADRATIC RECONSTRUCTION")
print()
print("Objective:")
print("  1. Start from the compressed variable:")
print("       d = sqrt(-4K - 3N^2 + 6N + 1)")
print("  2. Recover:")
print("       S = p+q = (N+1-d)/2")
print("  3. Introduce the factor-gap variable:")
print("       g = q-p")
print("  4. Verify:")
print("       g^2 = S^2 - 4N")
print("  5. Verify the equivalent compressed identity:")
print("       (N+1-d)^2 - 16N = 4g^2")
print("  6. Use d+3 = (p-2)(q-2).")
print("  7. Construct the shifted factor quadratic:")
print("       Y^2 - (S-4)Y + (d+3) = 0")
print("     whose roots should be p-2 and q-2.")
print("  8. Recover p and q exactly from N and d.")
print("  9. Determine whether the smaller shifted branch is always p-2.")
print(" 10. Measure whether the factor-gap representation is simpler")
print("     than the original N,K quadratic.")
print()
print("No resultants are constructed.")
print("No symbolic multivariate factorization is performed.")
print("Only exact integer arithmetic is used.")
print()


# ==============================================================================================================
# TEST INSTANCES
# ==============================================================================================================

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


# ==============================================================================================================
# HELPERS
# ==============================================================================================================

def digits(x: int) -> int:
    x = abs(int(x))
    return 1 if x == 0 else len(str(x))


def bitlen(x: int) -> int:
    x = abs(int(x))
    return 1 if x == 0 else x.bit_length()


def exact_sqrt(n: int) -> tuple[int, bool]:
    if n < 0:
        return 0, False

    r = math.isqrt(n)
    return r, r * r == n


def gcd_abs(a: int, b: int) -> int:
    return math.gcd(abs(a), abs(b))


# ==============================================================================================================
# MAIN
# ==============================================================================================================

def main() -> None:

    global_checks = {
        "d_square": True,
        "d_factor": True,
        "S_reconstruction": True,
        "gap_square": True,
        "compressed_gap": True,
        "shifted_product": True,
        "shifted_quadratic": True,
        "shifted_discriminant": True,
        "factor_recovery": True,
        "factor_vieta": True,
        "branch_order": True,
        "all_exact": True,
    }

    rows = []

    for idx, (p, q) in enumerate(INSTANCES, start=1):

        # ------------------------------------------------------------------------------------------
        # BASIC VALUES
        # ------------------------------------------------------------------------------------------

        N = p * q
        S = p + q
        X = S + 1

        K = -N * N + N * X - X * X + 3 * X - 2

        # Known compressed variable.
        d_true = N - 2 * S + 1

        # Recover d only from N,K.
        D = -4 * K - 3 * N * N + 6 * N + 1
        d, d_square_ok = exact_sqrt(D)

        # ------------------------------------------------------------------------------------------
        # S RECONSTRUCTION
        # ------------------------------------------------------------------------------------------

        S_numerator = N + 1 - d
        S_integral = (S_numerator % 2 == 0)

        S_rec = S_numerator // 2 if S_integral else None
        S_ok = S_integral and S_rec == S

        # ------------------------------------------------------------------------------------------
        # FACTOR GAP
        # ------------------------------------------------------------------------------------------

        g_true = q - p

        gap_discriminant = S * S - 4 * N
        g_rec, gap_square_ok = exact_sqrt(gap_discriminant)

        gap_ok = (
            gap_square_ok
            and g_rec == g_true
        )

        # Equivalent compressed identity:
        #
        #   (N+1-d)^2 - 16N = 4(q-p)^2
        #
        compressed_gap_lhs = (N + 1 - d) ** 2 - 16 * N
        compressed_gap_rhs = 4 * g_true * g_true

        compressed_gap_ok = (
            compressed_gap_lhs == compressed_gap_rhs
            and compressed_gap_lhs == 4 * g_rec * g_rec
        )

        # ------------------------------------------------------------------------------------------
        # SHIFTED FACTOR PRODUCT
        # ------------------------------------------------------------------------------------------

        shifted_product = d + 3
        shifted_product_true = (p - 2) * (q - 2)

        shifted_product_ok = shifted_product == shifted_product_true

        # ------------------------------------------------------------------------------------------
        # SHIFTED FACTOR QUADRATIC
        # ------------------------------------------------------------------------------------------
        #
        # Let:
        #
        #   y1 = p-2
        #   y2 = q-2
        #
        # Then:
        #
        #   y1+y2 = p+q-4 = S-4
        #   y1*y2 = (p-2)(q-2) = d+3
        #
        # Therefore:
        #
        #   Y^2 - (S-4)Y + (d+3) = 0
        #

        shifted_sum = S - 4
        shifted_constant = d + 3

        shifted_root_1 = p - 2
        shifted_root_2 = q - 2

        shifted_poly_root_1 = (
            shifted_root_1 * shifted_root_1
            - shifted_sum * shifted_root_1
            + shifted_constant
        )

        shifted_poly_root_2 = (
            shifted_root_2 * shifted_root_2
            - shifted_sum * shifted_root_2
            + shifted_constant
        )

        shifted_quadratic_ok = (
            shifted_poly_root_1 == 0
            and shifted_poly_root_2 == 0
        )

        # ------------------------------------------------------------------------------------------
        # SHIFTED QUADRATIC VIETA
        # ------------------------------------------------------------------------------------------

        shifted_vieta_sum_ok = (
            shifted_root_1 + shifted_root_2 == shifted_sum
        )

        shifted_vieta_product_ok = (
            shifted_root_1 * shifted_root_2 == shifted_constant
        )

        shifted_vieta_ok = (
            shifted_vieta_sum_ok
            and shifted_vieta_product_ok
        )

        # ------------------------------------------------------------------------------------------
        # SHIFTED DISCRIMINANT
        # ------------------------------------------------------------------------------------------

        shifted_discriminant = (
            shifted_sum * shifted_sum
            - 4 * shifted_constant
        )

        shifted_gap_ok = (
            shifted_discriminant == g_true * g_true
        )

        shifted_gap_root, shifted_gap_square_ok = exact_sqrt(
            shifted_discriminant
        )

        shifted_gap_recovery_ok = (
            shifted_gap_square_ok
            and shifted_gap_root == g_true
        )

        # ------------------------------------------------------------------------------------------
        # RECOVER p,q FROM N,d
        # ------------------------------------------------------------------------------------------
        #
        # We know:
        #
        #   S = p+q
        #   g = q-p
        #
        # hence:
        #
        #   p = (S-g)/2
        #   q = (S+g)/2
        #

        pq_integral = (
            S_integral
            and ((S_rec - g_rec) % 2 == 0)
            and ((S_rec + g_rec) % 2 == 0)
        )

        p_rec = (S_rec - g_rec) // 2 if pq_integral else None
        q_rec = (S_rec + g_rec) // 2 if pq_integral else None

        factor_recovery_ok = (
            pq_integral
            and p_rec == p
            and q_rec == q
        )

        # ------------------------------------------------------------------------------------------
        # RECOVER p,q DIRECTLY FROM SHIFTED ROOTS
        # ------------------------------------------------------------------------------------------

        shifted_gap_for_formula = shifted_gap_root

        shifted_minus = (
            shifted_sum - shifted_gap_for_formula
        ) // 2

        shifted_plus = (
            shifted_sum + shifted_gap_for_formula
        ) // 2

        shifted_minus_integral = (
            (shifted_sum - shifted_gap_for_formula) % 2 == 0
        )

        shifted_plus_integral = (
            (shifted_sum + shifted_gap_for_formula) % 2 == 0
        )

        shifted_branch_ok = (
            shifted_minus_integral
            and shifted_plus_integral
            and shifted_minus == p - 2
            and shifted_plus == q - 2
        )

        # Recover original factors.
        p_from_shift = shifted_minus + 2
        q_from_shift = shifted_plus + 2

        shifted_factor_recovery_ok = (
            shifted_branch_ok
            and p_from_shift * q_from_shift == N
        )

        # ------------------------------------------------------------------------------------------
        # ORDER / SIZE CHECKS
        # ------------------------------------------------------------------------------------------

        branch_order_ok = (
            p < q
            and g_true > 0
            and shifted_minus < shifted_plus
            and shifted_minus == p - 2
        )

        # Some extra gcd data for the new gap.
        gcd_signatures = {
            "gcd(g,N)": gcd_abs(g_true, N),
            "gcd(g,S)": gcd_abs(g_true, S),
            "gcd(g,p)": gcd_abs(g_true, p),
            "gcd(g,q)": gcd_abs(g_true, q),
            "gcd(g,d)": gcd_abs(g_true, d),
        }

        # ------------------------------------------------------------------------------------------
        # INSTANCE STATUS
        # ------------------------------------------------------------------------------------------

        instance_ok = all([
            d_square_ok,
            d == d_true,
            S_ok,
            shifted_product_ok,
            gap_ok,
            compressed_gap_ok,
            shifted_quadratic_ok,
            shifted_vieta_ok,
            shifted_gap_ok,
            shifted_gap_recovery_ok,
            factor_recovery_ok,
            shifted_factor_recovery_ok,
            branch_order_ok,
        ])

        # ------------------------------------------------------------------------------------------
        # PRINT
        # ------------------------------------------------------------------------------------------

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p      = {p}")
        print(f"  q      = {q}")
        print(f"  N      = {N}")
        print(f"  S      = {S}")
        print(f"  X      = {X}")
        print(f"  K      = {K}")
        print()

        print("COMPRESSED d RECOVERY")
        print("  D = -4K - 3N^2 + 6N + 1")
        print(f"  D = {D}")
        print(f"  D = d^2               : {d_square_ok}")
        print(f"  recovered d           = {d}")
        print(f"  true d                = {d_true}")
        print(f"  d recovery            : {d == d_true}")
        print()

        print("S RECOVERY")
        print("  S = (N+1-d)/2")
        print(f"  numerator = {S_numerator}")
        print(f"  parity / integrality = {S_integral}")
        print(f"  reconstructed S      = {S_rec}")
        print(f"  S recovery            = {S_ok}")
        print()

        print("FACTOR GAP")
        print("  g = q-p")
        print(f"  g_true = {g_true}")
        print("  g^2 = S^2 - 4N")
        print(f"  gap discriminant = {gap_discriminant}")
        print(f"  exact square       = {gap_square_ok}")
        print(f"  recovered g        = {g_rec}")
        print(f"  gap recovery       = {gap_ok}")
        print()

        print("COMPRESSED GAP IDENTITY")
        print("  (N+1-d)^2 - 16N = 4(q-p)^2")
        print(f"  lhs = {compressed_gap_lhs}")
        print(f"  rhs = {compressed_gap_rhs}")
        print(f"  exact identity = {compressed_gap_ok}")
        print()

        print("SHIFTED PRODUCT")
        print("  d+3 = (p-2)(q-2)")
        print(f"  d+3 = {shifted_product}")
        print(f"  (p-2)(q-2) = {shifted_product_true}")
        print(f"  exact identity = {shifted_product_ok}")
        print()

        print("SHIFTED FACTOR QUADRATIC")
        print("  Y^2 - (S-4)Y + (d+3) = 0")
        print(f"  degree = 2")
        print(f"  leading coefficient = 1")
        print(f"  coefficient of Y = {-shifted_sum}")
        print(f"  constant = {shifted_constant}")
        print()

        print("SHIFTED ROOTS")
        print(f"  y1 = p-2 = {p - 2}")
        print(f"  y2 = q-2 = {q - 2}")
        print(f"  polynomial(y1) = {shifted_poly_root_1}")
        print(f"  polynomial(y2) = {shifted_poly_root_2}")
        print(f"  both exact roots = {shifted_quadratic_ok}")
        print()

        print("SHIFTED VIETA")
        print(f"  y1+y2 = S-4 : {shifted_vieta_sum_ok}")
        print(f"  y1*y2 = d+3 : {shifted_vieta_product_ok}")
        print(f"  complete Vieta = {shifted_vieta_ok}")
        print()

        print("SHIFTED DISCRIMINANT")
        print("  Delta_shift = (S-4)^2 - 4(d+3)")
        print(f"  Delta_shift = {shifted_discriminant}")
        print(f"  Delta_shift = (q-p)^2 : {shifted_gap_ok}")
        print(f"  exact sqrt = {shifted_gap_root}")
        print(f"  gap sqrt recovery = {shifted_gap_recovery_ok}")
        print()

        print("DIRECT FACTOR RECOVERY FROM N,d")
        print("  p = (S-g)/2")
        print("  q = (S+g)/2")
        print(f"  reconstructed p = {p_rec}")
        print(f"  reconstructed q = {q_rec}")
        print(f"  exact factor recovery = {factor_recovery_ok}")
        print()

        print("SHIFTED QUADRATIC FORMULA")
        print("  Y_MINUS = (S-4-sqrt(Delta_shift))/2")
        print("  Y_PLUS  = (S-4+sqrt(Delta_shift))/2")
        print(f"  Y_MINUS = {shifted_minus}")
        print(f"  Y_PLUS  = {shifted_plus}")
        print(f"  Y_MINUS == p-2 : {shifted_minus == p - 2}")
        print(f"  Y_PLUS  == q-2 : {shifted_plus == q - 2}")
        print(f"  shifted branch recovery = {shifted_branch_ok}")
        print()

        print("ORIGINAL FACTORS FROM SHIFTED ROOTS")
        print("  p = Y_MINUS + 2")
        print("  q = Y_PLUS + 2")
        print(f"  p recovered = {p_from_shift}")
        print(f"  q recovered = {q_from_shift}")
        print(f"  p*q == N = {p_from_shift * q_from_shift == N}")
        print(f"  shifted factor recovery = {shifted_factor_recovery_ok}")
        print()

        print("GAP GCD SIGNATURE")
        for name, value in gcd_signatures.items():
            print(f"  {name:<12} = {value}")
        print()

        print("SIZE DATA")
        print(f"  N digits          = {digits(N)}")
        print(f"  K digits          = {digits(K)}")
        print(f"  d digits          = {digits(d)}")
        print(f"  S digits          = {digits(S)}")
        print(f"  gap g digits      = {digits(g_true)}")
        print(f"  d+3 digits        = {digits(d + 3)}")
        print(f"  gap^2 digits      = {digits(g_true * g_true)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {instance_ok}")
        print()

        # ------------------------------------------------------------------------------------------
        # GLOBAL UPDATE
        # ------------------------------------------------------------------------------------------

        global_checks["d_square"] &= d_square_ok and (d == d_true)
        global_checks["d_factor"] &= (d == d_true)
        global_checks["S_reconstruction"] &= S_ok
        global_checks["gap_square"] &= gap_ok
        global_checks["compressed_gap"] &= compressed_gap_ok
        global_checks["shifted_product"] &= shifted_product_ok
        global_checks["shifted_quadratic"] &= shifted_quadratic_ok
        global_checks["shifted_discriminant"] &= (
            shifted_gap_ok and shifted_gap_recovery_ok
        )
        global_checks["factor_recovery"] &= (
            factor_recovery_ok and shifted_factor_recovery_ok
        )
        global_checks["factor_vieta"] &= shifted_vieta_ok
        global_checks["branch_order"] &= branch_order_ok
        global_checks["all_exact"] &= instance_ok

        rows.append({
            "instance": idx,
            "N_digits": digits(N),
            "K_digits": digits(K),
            "d_digits": digits(d),
            "S_digits": digits(S),
            "gap_digits": digits(g_true),
            "d_bitlen": bitlen(d),
            "gap_bitlen": bitlen(g_true),
            "instance_ok": instance_ok,
        })

    # ==========================================================================================================
    # GLOBAL SUMMARY
    # ==========================================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("COMPRESSED d CHAIN")
    print(f"  all d-square recoveries       = {global_checks['d_square']}")
    print(f"  all S reconstructions         = {global_checks['S_reconstruction']}")
    print()

    print("FACTOR GAP")
    print(f"  all gap-square identities     = {global_checks['gap_square']}")
    print(f"  all compressed gap identities = {global_checks['compressed_gap']}")
    print()

    print("SHIFTED FACTOR STRUCTURE")
    print(f"  all d+3 product identities    = {global_checks['shifted_product']}")
    print(f"  all shifted quadratics        = {global_checks['shifted_quadratic']}")
    print(f"  all shifted discriminants     = {global_checks['shifted_discriminant']}")
    print(f"  all shifted Vieta identities  = {global_checks['factor_vieta']}")
    print()

    print("FACTOR RECOVERY")
    print(f"  all direct p,q recoveries     = {global_checks['factor_recovery']}")
    print(f"  all branch-order checks       = {global_checks['branch_order']}")
    print()

    print("SIZE RANGES")
    print(
        f"  N digits       = "
        f"{min(r['N_digits'] for r in rows)} .. "
        f"{max(r['N_digits'] for r in rows)}"
    )
    print(
        f"  K digits       = "
        f"{min(r['K_digits'] for r in rows)} .. "
        f"{max(r['K_digits'] for r in rows)}"
    )
    print(
        f"  d digits       = "
        f"{min(r['d_digits'] for r in rows)} .. "
        f"{max(r['d_digits'] for r in rows)}"
    )
    print(
        f"  S digits       = "
        f"{min(r['S_digits'] for r in rows)} .. "
        f"{max(r['S_digits'] for r in rows)}"
    )
    print(
        f"  gap digits     = "
        f"{min(r['gap_digits'] for r in rows)} .. "
        f"{max(r['gap_digits'] for r in rows)}"
    )
    print()

    print("INSTANCE STATUS")
    for row in rows:
        print(
            f"  instance {row['instance']:>2}: "
            f"{row['instance_ok']}"
        )

    print()

    all_pass = all(row["instance_ok"] for row in rows)

    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINAL STATUS")
    print("=" * 110)
    print()
    print("The new reconstruction chain is:")
    print()
    print("  (N,K)")
    print("     -> D = -4K - 3N^2 + 6N + 1")
    print("     -> d = sqrt(D)")
    print("     -> S = (N+1-d)/2")
    print("     -> g = sqrt(S^2 - 4N)")
    print("     -> p = (S-g)/2")
    print("     -> q = (S+g)/2")
    print()
    print("Equivalently, the shifted factors satisfy:")
    print()
    print("  (p-2) + (q-2) = S-4")
    print("  (p-2)(q-2) = d+3")
    print()
    print("and therefore:")
    print()
    print("  Y^2 - (S-4)Y + (d+3) = 0")
    print()
    print("with roots:")
    print()
    print("  Y = p-2, q-2")
    print()
    print("The key new square identity is:")
    print()
    print("  (N+1-d)^2 - 16N = 4(q-p)^2")
    print()
    print("Thus the experiment tests whether the compressed d")
    print("representation exposes an exact factor-gap reconstruction.")
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {all_pass}")
    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
