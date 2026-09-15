from __future__ import annotations

import math


# ==============================================================================================================
# EXPERIMENT 400
# ==============================================================================================================

EXP_NO = 400

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT DIRECT N,d -> FACTOR-PAIR RECONSTRUCTION")
print()
print("Objective:")
print("  1. Start only from N and K.")
print("  2. Recover:")
print("       D = -4K - 3N^2 + 6N + 1")
print("       d = sqrt(D)")
print("  3. Eliminate S and g as explicit reconstruction stages.")
print("  4. Form the direct radicand:")
print("       G = ((N+1-d)^2 - 16N) / 4")
print("     so that:")
print("       G = (q-p)^2")
print("  5. Recover the factor gap directly from N,d:")
print("       g = sqrt(G)")
print("  6. Recover p,q directly from N,d,g.")
print("  7. Verify the equivalent one-line expressions:")
print("       p = (N-d-7-2g)/4 + 2")
print("       q = (N-d-7+2g)/4 + 2")
print("  8. Test whether the complete chain can be represented by")
print("     nested exact square roots with no S-variable.")
print("  9. Search for additional simplifications of the direct")
print("     factor radicand G.")
print(" 10. Audit whether G has a simpler factor form involving d")
print("     and low-degree expressions in N.")
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
        "direct_G_integral": True,
        "G_square": True,
        "g_recovery": True,
        "direct_p": True,
        "direct_q": True,
        "factor_product": True,
        "factor_sum": True,
        "one_line_p": True,
        "one_line_q": True,
        "nested_chain": True,
        "alternate_radical": True,
        "all_exact": True,
    }

    rows = []

    for idx, (p, q) in enumerate(INSTANCES, start=1):

        # ------------------------------------------------------------------------------------------
        # BASIC CONSTRUCTION
        # ------------------------------------------------------------------------------------------

        N = p * q
        S = p + q
        X = S + 1

        K = -N * N + N * X - X * X + 3 * X - 2

        d_true = N - 2 * S + 1
        g_true = q - p

        # ------------------------------------------------------------------------------------------
        # FIRST RADICAL: d
        # ------------------------------------------------------------------------------------------

        D = -4 * K - 3 * N * N + 6 * N + 1

        d, D_square_ok = exact_sqrt(D)
        d_ok = D_square_ok and d == d_true

        # ------------------------------------------------------------------------------------------
        # DIRECT SECOND RADICAND
        # ------------------------------------------------------------------------------------------
        #
        # From:
        #
        #   S = (N+1-d)/2
        #
        # and:
        #
        #   g^2 = S^2 - 4N
        #
        # we obtain:
        #
        #   4g^2 = (N+1-d)^2 - 16N
        #
        # Therefore define:
        #
        #   G = ((N+1-d)^2 - 16N)/4.
        #
        # This is constructed directly from N,d without first assigning S.
        #

        G_num = (N + 1 - d) ** 2 - 16 * N
        G_integral = (G_num % 4 == 0)

        G = G_num // 4 if G_integral else None

        if G is not None:
            g, G_square_ok = exact_sqrt(G)
        else:
            g = 0
            G_square_ok = False

        g_ok = G_square_ok and g == g_true

        # ------------------------------------------------------------------------------------------
        # DIRECT FACTOR FORMULAS FROM N,d,g
        # ------------------------------------------------------------------------------------------

        p_num = N - d - 7 - 2 * g
        q_num = N - d - 7 + 2 * g

        p_integral = (p_num % 4 == 0)
        q_integral = (q_num % 4 == 0)

        p_direct = p_num // 4 + 2 if p_integral else None
        q_direct = q_num // 4 + 2 if q_integral else None

        p_ok = p_integral and p_direct == p
        q_ok = q_integral and q_direct == q

        # ------------------------------------------------------------------------------------------
        # FACTOR SUM / PRODUCT
        # ------------------------------------------------------------------------------------------

        factor_sum_ok = (
            p_direct is not None
            and q_direct is not None
            and p_direct + q_direct == N + 1 - d
        )

        factor_product_ok = (
            p_direct is not None
            and q_direct is not None
            and p_direct * q_direct == N
        )

        # ------------------------------------------------------------------------------------------
        # ONE-LINE RADICAL EXPRESSIONS
        # ------------------------------------------------------------------------------------------
        #
        # g = sqrt(((N+1-d)^2 - 16N)/4)
        #
        # Since the numerator is divisible by 4 for all valid instances,
        # the direct factor pair becomes:
        #
        # p = (N-d-7 - 2*sqrt(G))/4 + 2
        # q = (N-d-7 + 2*sqrt(G))/4 + 2
        #

        one_line_p = p_direct
        one_line_q = q_direct

        one_line_p_ok = (
            p_integral
            and one_line_p == p
        )

        one_line_q_ok = (
            q_integral
            and one_line_q == q
        )

        # ------------------------------------------------------------------------------------------
        # NESTED CHAIN CHECK
        # ------------------------------------------------------------------------------------------

        nested_chain_ok = all([
            D_square_ok,
            d_ok,
            G_integral,
            G_square_ok,
            g_ok,
            p_integral,
            q_integral,
            p_ok,
            q_ok,
        ])

        # ------------------------------------------------------------------------------------------
        # ALTERNATIVE DIRECT RADICAND FORMS
        # ------------------------------------------------------------------------------------------
        #
        # Let:
        #
        #   A = N-d-7.
        #
        # We already know:
        #
        #   A^2 - 16(d+3) = 4g^2.
        #
        # Therefore:
        #
        #   G = (A^2 - 16(d+3))/4.
        #
        # This should equal the previous direct G exactly.
        #
        # Also:
        #
        #   A = N-d-7
        #
        # and:
        #
        #   (N+1-d)^2 - 16N = A^2 - 16(d+3).
        #

        A = N - d - 7

        G_alt_num = A * A - 16 * (d + 3)
        G_alt_integral = (G_alt_num % 4 == 0)
        G_alt = G_alt_num // 4 if G_alt_integral else None

        alternate_radical_ok = (
            G_integral
            and G_alt_integral
            and G == G_alt
            and G == g * g
        )

        # ------------------------------------------------------------------------------------------
        # EXTRA LOW-COMPLEXITY TESTS
        # ------------------------------------------------------------------------------------------

        expressions = {
            "N-d-7": A,
            "N-d-3": N - d - 3,
            "N-d-1": N - d - 1,
            "N+1-d": N + 1 - d,
            "N+3-d": N + 3 - d,
            "A^2": A * A,
            "16(d+3)": 16 * (d + 3),
            "A^2-16(d+3)": G_alt_num,
        }

        gcd_signature = {
            "gcd(d,N)": gcd_abs(d, N),
            "gcd(G,N)": gcd_abs(G, N) if G is not None else 0,
            "gcd(G,d)": gcd_abs(G, d) if G is not None else 0,
            "gcd(g,N)": gcd_abs(g, N),
            "gcd(g,d)": gcd_abs(g, d),
            "gcd(g,d+3)": gcd_abs(g, d + 3),
        }

        divisibility_signature = {
            "4 | N+1-d": ((N + 1 - d) % 4 == 0),
            "4 | N-d-7-2g": p_integral,
            "4 | N-d-7+2g": q_integral,
            "G | (d+3)^2": (
                G is not None
                and G != 0
                and ((d + 3) ** 2) % G == 0
            ),
        }

        instance_ok = all([
            D_square_ok,
            d_ok,
            G_integral,
            G_square_ok,
            g_ok,
            p_integral,
            q_integral,
            p_ok,
            q_ok,
            factor_sum_ok,
            factor_product_ok,
            one_line_p_ok,
            one_line_q_ok,
            nested_chain_ok,
            alternate_radical_ok,
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
        print(f"  D = d^2 = {d*d} : {D_square_ok}")
        print(f"  recovered d = {d}")
        print(f"  true d      = {d_true}")
        print(f"  d recovery  = {d_ok}")
        print()

        print("DIRECT SECOND RADICAND")
        print("  G = ((N+1-d)^2 - 16N) / 4")
        print(f"  numerator = {G_num}")
        print(f"  divisible by 4 = {G_integral}")
        print(f"  G = {G}")
        print(f"  G = g^2 : {G_square_ok}")
        print(f"  recovered g = {g}")
        print(f"  true g      = {g_true}")
        print(f"  g recovery  = {g_ok}")
        print()

        print("DIRECT N,d FACTOR FORMULAS")
        print("  p = (N-d-7-2g)/4 + 2")
        print("  q = (N-d-7+2g)/4 + 2")
        print(f"  p numerator = {p_num}")
        print(f"  q numerator = {q_num}")
        print(f"  p numerator divisible by 4 = {p_integral}")
        print(f"  q numerator divisible by 4 = {q_integral}")
        print(f"  reconstructed p = {p_direct}")
        print(f"  reconstructed q = {q_direct}")
        print(f"  p recovery = {p_ok}")
        print(f"  q recovery = {q_ok}")
        print()

        print("FACTOR PAIR INVARIANTS")
        print("  reconstructed p+q = N+1-d")
        print(
            f"  {p_direct + q_direct if p_direct is not None and q_direct is not None else None}"
            f" == {N + 1 - d} : {factor_sum_ok}"
        )
        print("  reconstructed p*q = N")
        print(
            f"  {p_direct * q_direct if p_direct is not None and q_direct is not None else None}"
            f" == {N} : {factor_product_ok}"
        )
        print()

        print("ALTERNATE DIRECT RADICAND")
        print("  A = N-d-7")
        print(f"  A = {A}")
        print("  G = (A^2 - 16(d+3))/4")
        print(f"  A^2 = {A*A}")
        print(f"  16(d+3) = {16*(d+3)}")
        print(f"  A^2 - 16(d+3) = {G_alt_num}")
        print(f"  G_alt = {G_alt}")
        print(f"  G_alt == G == g^2 : {alternate_radical_ok}")
        print()

        print("ONE-LINE N,K -> FACTOR RECONSTRUCTION")
        print("  d = sqrt(-4K - 3N^2 + 6N + 1)")
        print(
            "  g = sqrt((((N+1-d)^2 - 16N)/4))"
        )
        print(
            "  p = (N-d-7-2g)/4 + 2"
        )
        print(
            "  q = (N-d-7+2g)/4 + 2"
        )
        print(f"  one-line p = {one_line_p}")
        print(f"  one-line q = {one_line_q}")
        print(f"  one-line p exact = {one_line_p_ok}")
        print(f"  one-line q exact = {one_line_q_ok}")
        print(f"  complete nested chain = {nested_chain_ok}")
        print()

        print("LOW-COMPLEXITY EXPRESSIONS")
        for name, value in expressions.items():
            print(
                f"  {name:<18} digits={digits(value):>4}"
            )
        print()

        print("GCD SIGNATURE")
        for name, value in gcd_signature.items():
            print(f"  {name:<18} = {value}")
        print()

        print("DIVISIBILITY SIGNATURE")
        for name, value in divisibility_signature.items():
            print(f"  {name:<28} = {value}")
        print()

        print("SIZE DATA")
        print(f"  N digits          = {digits(N)}")
        print(f"  K digits          = {digits(K)}")
        print(f"  d digits          = {digits(d)}")
        print(f"  G digits          = {digits(G)}")
        print(f"  g digits          = {digits(g)}")
        print(f"  A digits          = {digits(A)}")
        print(f"  p digits          = {digits(p)}")
        print(f"  q digits          = {digits(q)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {instance_ok}")
        print()

        # ------------------------------------------------------------------------------------------
        # GLOBAL
        # ------------------------------------------------------------------------------------------

        global_checks["D_square"] &= D_square_ok
        global_checks["d_recovery"] &= d_ok
        global_checks["direct_G_integral"] &= G_integral
        global_checks["G_square"] &= G_square_ok
        global_checks["g_recovery"] &= g_ok
        global_checks["direct_p"] &= p_ok
        global_checks["direct_q"] &= q_ok
        global_checks["factor_product"] &= factor_product_ok
        global_checks["factor_sum"] &= factor_sum_ok
        global_checks["one_line_p"] &= one_line_p_ok
        global_checks["one_line_q"] &= one_line_q_ok
        global_checks["nested_chain"] &= nested_chain_ok
        global_checks["alternate_radical"] &= alternate_radical_ok
        global_checks["all_exact"] &= instance_ok

        rows.append({
            "instance": idx,
            "N_digits": digits(N),
            "K_digits": digits(K),
            "d_digits": digits(d),
            "G_digits": digits(G),
            "g_digits": digits(g),
            "A_digits": digits(A),
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
    print(f"  all D=d^2 checks            = {global_checks['D_square']}")
    print(f"  all d recoveries            = {global_checks['d_recovery']}")
    print()

    print("DIRECT SECOND RADICAND")
    print(f"  all G integrality checks    = {global_checks['direct_G_integral']}")
    print(f"  all G=g^2 checks            = {global_checks['G_square']}")
    print(f"  all g recoveries            = {global_checks['g_recovery']}")
    print()

    print("DIRECT FACTOR RECONSTRUCTION")
    print(f"  all p recoveries            = {global_checks['direct_p']}")
    print(f"  all q recoveries            = {global_checks['direct_q']}")
    print(f"  all p+q checks              = {global_checks['factor_sum']}")
    print(f"  all p*q checks              = {global_checks['factor_product']}")
    print()

    print("NESTED CLOSED CHAIN")
    print(f"  all one-line p formulas     = {global_checks['one_line_p']}")
    print(f"  all one-line q formulas     = {global_checks['one_line_q']}")
    print(f"  all nested chains           = {global_checks['nested_chain']}")
    print(f"  all alternate G forms      = {global_checks['alternate_radical']}")
    print()

    print("SIZE RANGES")
    print(
        f"  N digits = "
        f"{min(r['N_digits'] for r in rows)} .. "
        f"{max(r['N_digits'] for r in rows)}"
    )
    print(
        f"  K digits = "
        f"{min(r['K_digits'] for r in rows)} .. "
        f"{max(r['K_digits'] for r in rows)}"
    )
    print(
        f"  d digits = "
        f"{min(r['d_digits'] for r in rows)} .. "
        f"{max(r['d_digits'] for r in rows)}"
    )
    print(
        f"  G digits = "
        f"{min(r['G_digits'] for r in rows)} .. "
        f"{max(r['G_digits'] for r in rows)}"
    )
    print(
        f"  g digits = "
        f"{min(r['g_digits'] for r in rows)} .. "
        f"{max(r['g_digits'] for r in rows)}"
    )
    print(
        f"  A digits = "
        f"{min(r['A_digits'] for r in rows)} .. "
        f"{max(r['A_digits'] for r in rows)}"
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

    print("The direct reconstruction tested is:")
    print()
    print("  d = sqrt(-4K - 3N^2 + 6N + 1)")
    print()
    print("then, without explicitly introducing S:")
    print()
    print("  G = ((N+1-d)^2 - 16N)/4")
    print("  g = sqrt(G)")
    print()
    print("and finally:")
    print()
    print("  p = (N-d-7-2g)/4 + 2")
    print("  q = (N-d-7+2g)/4 + 2")
    print()
    print("The second radicand is exactly:")
    print()
    print("  G = (q-p)^2")
    print()
    print("and equivalently:")
    print()
    print("  G = ((N-d-7)^2 - 16(d+3))/4")
    print()
    print("Thus the entire tested factor-recovery chain can be")
    print("written as two nested exact square roots followed by")
    print("two affine expressions.")
    print()
    print("This experiment checks whether that direct N,K")
    print("representation is algebraically self-consistent")
    print("without treating S or g as independent inputs.")
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {all_pass}")
    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
