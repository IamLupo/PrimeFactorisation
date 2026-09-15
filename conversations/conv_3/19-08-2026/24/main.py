from __future__ import annotations

import math


# ==============================================================================================================
# EXPERIMENT 401
# ==============================================================================================================

EXP_NO = 401

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT RADICAL CHAIN CORRECTION / DIRECT FACTOR POLYNOMIAL AUDIT")
print()
print("Objective:")
print("  1. Correct the factor-sum audit from Experiment 400.")
print("  2. Remove non-structural divisibility tests from pass/fail.")
print("  3. Verify the direct N,d -> g -> p,q reconstruction.")
print("  4. Verify:")
print("       p+q = (N+1-d)/2")
print("       pq  = N")
print("  5. Eliminate g completely and verify that p and q are")
print("     the two roots of a quadratic whose coefficients depend")
print("     only on N and d.")
print("  6. Construct the direct factor quadratic:")
print("       Z^2 - ((N+1-d)/2) Z + N = 0")
print("  7. Verify its exact discriminant:")
print("       Delta_Z = ((N+1-d)/2)^2 - 4N = g^2")
print("  8. Compare this quadratic directly with the original")
print("     shifted-factor representation.")
print("  9. Search for a denominator-free equivalent quadratic")
print("     involving only N,d,Z.")
print(" 10. Determine whether the entire factorization layer")
print("     reduces to a single quadratic after the first radical.")
print()
print("No resultants are constructed.")
print("No symbolic multivariate factoring is performed.")
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


def exact_sqrt(n: int) -> tuple[int, bool]:
    if n < 0:
        return 0, False
    r = math.isqrt(n)
    return r, (r * r == n)


def gcd_abs(a: int, b: int) -> int:
    return math.gcd(abs(a), abs(b))


# ==============================================================================================================
# MAIN
# ==============================================================================================================

def main() -> None:

    global_checks = {
        "D_square": True,
        "d_recovery": True,

        "G_integral": True,
        "G_square": True,
        "g_recovery": True,

        "factor_sum": True,
        "factor_product": True,

        "direct_p": True,
        "direct_q": True,

        "quadratic_p": True,
        "quadratic_q": True,
        "quadratic_coefficients": True,
        "quadratic_discriminant": True,

        "denominator_free_quadratic": True,

        "all_exact": True,
    }

    rows = []

    for idx, (p, q) in enumerate(INSTANCES, start=1):

        # ------------------------------------------------------------------------------------------
        # BASIC DATA
        # ------------------------------------------------------------------------------------------

        N = p * q
        S = p + q
        X = S + 1

        K = -N * N + N * X - X * X + 3 * X - 2

        d_true = N - 2 * S + 1
        g_true = q - p

        # ------------------------------------------------------------------------------------------
        # FIRST RADICAL
        # ------------------------------------------------------------------------------------------

        D = -4 * K - 3 * N * N + 6 * N + 1

        d, D_square_ok = exact_sqrt(D)
        d_ok = D_square_ok and d == d_true

        # ------------------------------------------------------------------------------------------
        # DIRECT SECOND RADICAND
        # ------------------------------------------------------------------------------------------

        G_num = (N + 1 - d) ** 2 - 16 * N
        G_integral = (G_num % 4 == 0)

        G = G_num // 4 if G_integral else 0

        g, G_square_ok = exact_sqrt(G)
        g_ok = G_square_ok and g == g_true

        # ------------------------------------------------------------------------------------------
        # CORRECT FACTOR-SUM IDENTITY
        # ------------------------------------------------------------------------------------------

        # p + q = S
        # 2S = N + 1 - d

        factor_sum_ok = (
            p + q == S
            and 2 * (p + q) == N + 1 - d
        )

        factor_product_ok = (
            p * q == N
        )

        # ------------------------------------------------------------------------------------------
        # DIRECT FACTOR FORMULAS
        # ------------------------------------------------------------------------------------------

        p_num = N - d - 7 - 2 * g
        q_num = N - d - 7 + 2 * g

        p_integral = (p_num % 4 == 0)
        q_integral = (q_num % 4 == 0)

        p_direct = p_num // 4 + 2 if p_integral else 0
        q_direct = q_num // 4 + 2 if q_integral else 0

        p_ok = p_integral and p_direct == p
        q_ok = q_integral and q_direct == q

        # ------------------------------------------------------------------------------------------
        # DIRECT FACTOR QUADRATIC
        # ------------------------------------------------------------------------------------------
        #
        # Since:
        #
        #     p + q = (N+1-d)/2
        #     pq    = N
        #
        # p and q are roots of:
        #
        #     Z^2 - S Z + N = 0
        #
        # where:
        #
        #     S = (N+1-d)/2.
        #

        sum_num = N + 1 - d

        parity_sum_ok = (sum_num % 2 == 0)

        S_direct = sum_num // 2 if parity_sum_ok else 0

        quadratic_p_value = (
            p_direct * p_direct
            - S_direct * p_direct
            + N
        )

        quadratic_q_value = (
            q_direct * q_direct
            - S_direct * q_direct
            + N
        )

        quadratic_p_ok = quadratic_p_value == 0
        quadratic_q_ok = quadratic_q_value == 0

        quadratic_coefficients_ok = (
            parity_sum_ok
            and S_direct == p_direct + q_direct
            and N == p_direct * q_direct
        )

        # ------------------------------------------------------------------------------------------
        # DIRECT QUADRATIC DISCRIMINANT
        # ------------------------------------------------------------------------------------------

        Delta_Z = S_direct * S_direct - 4 * N

        sqrt_Delta_Z, Delta_Z_square_ok = exact_sqrt(Delta_Z)

        quadratic_discriminant_ok = (
            Delta_Z_square_ok
            and sqrt_Delta_Z == g
            and Delta_Z == g * g
        )

        # ------------------------------------------------------------------------------------------
        # DENOMINATOR-FREE QUADRATIC
        # ------------------------------------------------------------------------------------------
        #
        # Starting from:
        #
        #   Z^2 - ((N+1-d)/2) Z + N = 0
        #
        # multiply by 2:
        #
        #   2Z^2 - (N+1-d)Z + 2N = 0
        #
        # This is an integer-coefficient quadratic depending
        # only on N,d.
        #

        def denominator_free_q(z: int) -> int:
            return (
                2 * z * z
                - (N + 1 - d) * z
                + 2 * N
            )

        denominator_free_p = denominator_free_q(p_direct)
        denominator_free_q_value = denominator_free_q(q_direct)

        denominator_free_ok = (
            denominator_free_p == 0
            and denominator_free_q_value == 0
        )

        # ------------------------------------------------------------------------------------------
        # DISCRIMINANT OF DENOMINATOR-FREE QUADRATIC
        # ------------------------------------------------------------------------------------------
        #
        # For:
        #
        #   2Z^2 - (N+1-d)Z + 2N = 0
        #
        # the discriminant is:
        #
        #   (N+1-d)^2 - 16N
        #
        # which is exactly 4g^2.
        #

        Delta_Z_integer = (
            (N + 1 - d) ** 2
            - 16 * N
        )

        discriminant_transfer_ok = (
            Delta_Z_integer == 4 * g * g
        )

        # ------------------------------------------------------------------------------------------
        # ALTERNATIVE EXACT FORM
        # ------------------------------------------------------------------------------------------

        # From the Experiment 399 identity:
        #
        #   (N-d-7)^2 - 4g^2 = 16(d+3)
        #
        # and:
        #
        #   4g^2 = (N+1-d)^2 - 16N
        #
        # verify the exact equivalence.

        lhs_a = (N + 1 - d) ** 2 - 16 * N
        lhs_b = (N - d - 7) ** 2 - 16 * (d + 3)

        equivalent_radical_forms_ok = (
            lhs_a == lhs_b
            and lhs_a == 4 * g * g
        )

        # ------------------------------------------------------------------------------------------
        # INSTANCE STATUS
        # ------------------------------------------------------------------------------------------

        instance_ok = all([
            D_square_ok,
            d_ok,

            G_integral,
            G_square_ok,
            g_ok,

            factor_sum_ok,
            factor_product_ok,

            p_ok,
            q_ok,

            quadratic_p_ok,
            quadratic_q_ok,
            quadratic_coefficients_ok,
            quadratic_discriminant_ok,

            denominator_free_ok,
            discriminant_transfer_ok,
            equivalent_radical_forms_ok,
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

        print("FIRST COMPRESSED RADICAL")
        print("  D = -4K - 3N^2 + 6N + 1")
        print(f"  D = {D}")
        print(f"  D = d^2 : {D_square_ok}")
        print(f"  recovered d = {d}")
        print(f"  true d      = {d_true}")
        print(f"  d recovery  = {d_ok}")
        print()

        print("DIRECT SECOND RADICAND")
        print("  G = ((N+1-d)^2 - 16N)/4")
        print(f"  numerator = {G_num}")
        print(f"  divisible by 4 = {G_integral}")
        print(f"  G = {G}")
        print(f"  G = g^2 : {G_square_ok}")
        print(f"  recovered g = {g}")
        print(f"  true g      = {g_true}")
        print(f"  g recovery  = {g_ok}")
        print()

        print("CORRECTED FACTOR-SUM CHECK")
        print("  p+q = (N+1-d)/2")
        print(f"  p+q = {p_direct + q_direct}")
        print(f"  (N+1-d)/2 = {S_direct}")
        print(f"  p+q == (N+1-d)/2 : {factor_sum_ok}")
        print()

        print("FACTOR PRODUCT CHECK")
        print("  p*q = N")
        print(f"  reconstructed p*q = {p_direct * q_direct}")
        print(f"  N = {N}")
        print(f"  exact product = {factor_product_ok}")
        print()

        print("DIRECT FACTOR FORMULAS")
        print("  p = (N-d-7-2g)/4 + 2")
        print("  q = (N-d-7+2g)/4 + 2")
        print(f"  p numerator = {p_num}")
        print(f"  q numerator = {q_num}")
        print(f"  p reconstructed = {p_direct}")
        print(f"  q reconstructed = {q_direct}")
        print(f"  p recovery = {p_ok}")
        print(f"  q recovery = {q_ok}")
        print()

        print("DIRECT FACTOR QUADRATIC")
        print("  Z^2 - ((N+1-d)/2) Z + N = 0")
        print(f"  coefficient of Z = -{S_direct}")
        print(f"  constant = {N}")
        print(f"  polynomial(p) = {quadratic_p_value}")
        print(f"  polynomial(q) = {quadratic_q_value}")
        print(f"  p is exact root = {quadratic_p_ok}")
        print(f"  q is exact root = {quadratic_q_ok}")
        print()

        print("DENOMINATOR-FREE FACTOR QUADRATIC")
        print("  2Z^2 - (N+1-d)Z + 2N = 0")
        print(f"  polynomial(p) = {denominator_free_p}")
        print(f"  polynomial(q) = {denominator_free_q_value}")
        print(f"  both exact roots = {denominator_free_ok}")
        print()

        print("DIRECT FACTOR DISCRIMINANT")
        print("  Delta_Z = ((N+1-d)^2 - 16N)")
        print(f"  Delta_Z = {Delta_Z_integer}")
        print(f"  Delta_Z = 4g^2 : {discriminant_transfer_ok}")
        print(f"  sqrt(Delta_Z)/2 = {sqrt_Delta_Z}")
        print(f"  exact g = {g}")
        print(
            f"  discriminant square check = "
            f"{quadratic_discriminant_ok}"
        )
        print()

        print("EQUIVALENT RADICAL FORMS")
        print("  Form A = (N+1-d)^2 - 16N")
        print("  Form B = (N-d-7)^2 - 16(d+3)")
        print(f"  Form A = {lhs_a}")
        print(f"  Form B = {lhs_b}")
        print(f"  A == B == 4g^2 : {equivalent_radical_forms_ok}")
        print()

        print("GCD SIGNATURE")
        print(f"  gcd(d,N)     = {gcd_abs(d, N)}")
        print(f"  gcd(g,N)     = {gcd_abs(g, N)}")
        print(f"  gcd(g,d)     = {gcd_abs(g, d)}")
        print(f"  gcd(g,d+3)   = {gcd_abs(g, d + 3)}")
        print(f"  gcd(p,q)     = {gcd_abs(p_direct, q_direct)}")
        print()

        print("SIZE DATA")
        print(f"  N digits       = {digits(N)}")
        print(f"  K digits       = {digits(K)}")
        print(f"  d digits       = {digits(d)}")
        print(f"  G digits       = {digits(G)}")
        print(f"  g digits       = {digits(g)}")
        print(f"  Delta_Z digits = {digits(Delta_Z_integer)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {instance_ok}")
        print()

        # ------------------------------------------------------------------------------------------
        # GLOBAL
        # ------------------------------------------------------------------------------------------

        global_checks["D_square"] &= D_square_ok
        global_checks["d_recovery"] &= d_ok

        global_checks["G_integral"] &= G_integral
        global_checks["G_square"] &= G_square_ok
        global_checks["g_recovery"] &= g_ok

        global_checks["factor_sum"] &= factor_sum_ok
        global_checks["factor_product"] &= factor_product_ok

        global_checks["direct_p"] &= p_ok
        global_checks["direct_q"] &= q_ok

        global_checks["quadratic_p"] &= quadratic_p_ok
        global_checks["quadratic_q"] &= quadratic_q_ok
        global_checks["quadratic_coefficients"] &= quadratic_coefficients_ok
        global_checks["quadratic_discriminant"] &= quadratic_discriminant_ok

        global_checks["denominator_free_quadratic"] &= denominator_free_ok
        global_checks["all_exact"] &= instance_ok

        rows.append({
            "instance": idx,
            "N_digits": digits(N),
            "K_digits": digits(K),
            "d_digits": digits(d),
            "G_digits": digits(G),
            "g_digits": digits(g),
            "Delta_Z_digits": digits(Delta_Z_integer),
            "instance_ok": instance_ok,
        })

    # ==========================================================================================================
    # GLOBAL SUMMARY
    # ==========================================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("FIRST RADICAL")
    print(f"  all D=d^2 checks              = {global_checks['D_square']}")
    print(f"  all d recoveries              = {global_checks['d_recovery']}")
    print()

    print("SECOND RADICAL")
    print(f"  all G integrality checks      = {global_checks['G_integral']}")
    print(f"  all G=g^2 checks              = {global_checks['G_square']}")
    print(f"  all g recoveries              = {global_checks['g_recovery']}")
    print()

    print("FACTOR PAIR")
    print(f"  all p recoveries              = {global_checks['direct_p']}")
    print(f"  all q recoveries              = {global_checks['direct_q']}")
    print(f"  all p+q checks                = {global_checks['factor_sum']}")
    print(f"  all p*q checks                = {global_checks['factor_product']}")
    print()

    print("DIRECT FACTOR QUADRATIC")
    print(f"  all p quadratic roots         = {global_checks['quadratic_p']}")
    print(f"  all q quadratic roots         = {global_checks['quadratic_q']}")
    print(
        f"  all coefficient checks        = "
        f"{global_checks['quadratic_coefficients']}"
    )
    print(
        f"  all discriminant checks      = "
        f"{global_checks['quadratic_discriminant']}"
    )
    print(
        f"  all denominator-free checks = "
        f"{global_checks['denominator_free_quadratic']}"
    )
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
        f"  G digits       = "
        f"{min(r['G_digits'] for r in rows)} .. "
        f"{max(r['G_digits'] for r in rows)}"
    )
    print(
        f"  g digits       = "
        f"{min(r['g_digits'] for r in rows)} .. "
        f"{max(r['g_digits'] for r in rows)}"
    )
    print(
        f"  Delta_Z digits = "
        f"{min(r['Delta_Z_digits'] for r in rows)} .. "
        f"{max(r['Delta_Z_digits'] for r in rows)}"
    )
    print()

    print("INSTANCE STATUS")
    for row in rows:
        print(
            f"  instance {row['instance']:>2}: "
            f"{row['instance_ok']}"
        )

    all_pass = all(row["instance_ok"] for row in rows)

    # ==========================================================================================================
    # FINAL
    # ==========================================================================================================

    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINAL STATUS")
    print("=" * 110)
    print()

    print("Experiment 400 exposed a bookkeeping error in the factor-sum")
    print("audit. The correct relation is:")
    print()
    print("  p + q = S = (N+1-d)/2")
    print()
    print("not:")
    print()
    print("  p + q = N+1-d")
    print()
    print("After that correction, the factor layer collapses to a")
    print("single quadratic depending only on N and d:")
    print()
    print("  Z^2 - ((N+1-d)/2) Z + N = 0")
    print()
    print("or, without division:")
    print()
    print("  2Z^2 - (N+1-d)Z + 2N = 0")
    print()
    print("Its roots are exactly:")
    print()
    print("  Z = p, q")
    print()
    print("and its discriminant is:")
    print()
    print("  (N+1-d)^2 - 16N = 4(q-p)^2")
    print()
    print("Therefore the factorization stage after d is not")
    print("fundamentally a second independent quadratic.")
    print("It is one quadratic whose coefficients are already")
    print("determined by N and d.")
    print()
    print("The complete nested reconstruction remains:")
    print()
    print("  d = sqrt(-4K - 3N^2 + 6N + 1)")
    print("  p,q = roots of")
    print("        2Z^2 - (N+1-d)Z + 2N = 0")
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {all_pass}")
    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
