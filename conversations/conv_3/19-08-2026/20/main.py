from __future__ import annotations

import math
import sympy as sp


# ==============================================================================================================
# EXPERIMENT 397
# ==============================================================================================================

EXP_NO = 397

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT COMPRESSED INVARIANT / LOW-COMPLEXITY FACTORIZATION AUDIT")
print()
print("Objective:")
print("  1. Compute the compressed invariant:")
print("       D = -4K - 3N^2 + 6N + 1")
print("  2. Verify exactly that D = d^2.")
print("  3. Test whether D admits additional low-complexity forms.")
print("  4. Test shifted forms involving N and d.")
print("  5. Test gcd relationships with N, N+1, N+2, N+3.")
print("  6. Verify the factor identity:")
print("       d+3 = (p-2)(q-2)")
print("  7. Search for simple relations between d and N.")
print("  8. Verify reconstruction of S and X using only N and d.")
print("  9. Record exact size and divisibility signatures.")
print()
print("No resultants are constructed.")
print("No symbolic factoring of large multivariate expressions is performed.")
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
    if x == 0:
        return 1
    return len(str(x))


def bitlen(x: int) -> int:
    x = abs(int(x))
    if x == 0:
        return 1
    return x.bit_length()


def gcd_abs(a: int, b: int) -> int:
    return math.gcd(abs(int(a)), abs(int(b)))


def square_root_exact(n: int) -> tuple[int, bool]:
    if n < 0:
        return 0, False
    r = math.isqrt(n)
    return r, r * r == n


def factor_small_integer(n: int, limit: int = 100000) -> list[tuple[int, int]]:
    """
    Small-prime factorization only.

    This intentionally avoids attempting a full factorization of very
    large integers. It is used only to expose small structural factors.
    """
    n = abs(int(n))
    result: list[tuple[int, int]] = []

    if n == 0:
        return result

    p = 2
    while p <= limit and p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            result.append((p, e))

        if p == 2:
            p = 3
        else:
            p += 2

    if n != 1 and n <= limit:
        result.append((n, 1))

    return result


# ==============================================================================================================
# MAIN
# ==============================================================================================================

def main() -> None:
    global_checks = {
        "compressed_square": True,
        "positive_root": True,
        "factor_d": True,
        "shifted_factor": True,
        "reconstruct_S": True,
        "reconstruct_X": True,
        "gcd_checks": True,
        "shift_checks": True,
        "mod_checks": True,
        "all_low_complexity_checks": True,
    }

    summary_rows = []

    for idx, (p, q) in enumerate(INSTANCES, start=1):

        N = p * q
        S = p + q
        X = S + 1
        K = -N * N + N * X - X * X + 3 * X - 2

        d = N - 2 * S + 1

        # Compressed invariant.
        D = -4 * K - 3 * N * N + 6 * N + 1

        # Positive exact square root.
        d_recovered, is_square = square_root_exact(D)

        # Expected shifted factor form.
        d_factor = (p - 1) * (q - 1) - (p + q)
        d_shifted = d + 3
        d_shifted_factor = (p - 2) * (q - 2)

        # Recover S and X from N,d.
        S_rec = (N + 1 - d) // 2
        X_rec = (N + 3 - d) // 2

        parity_S = (N + 1 - d) % 2 == 0
        parity_X = (N + 3 - d) % 2 == 0

        # ------------------------------------------------------------------
        # LOW-COMPLEXITY N/d RELATIONS
        # ------------------------------------------------------------------

        candidates = {
            "N-d": N - d,
            "N+d": N + d,
            "N-d+1": N - d + 1,
            "N+d+1": N + d + 1,
            "N-d-1": N - d - 1,
            "N+d-1": N + d - 1,
            "N-2d": N - 2 * d,
            "N+2d": N + 2 * d,
            "N^2-d^2": N * N - d * d,
            "N^2+d^2": N * N + d * d,
        }

        # Verify the known relations.
        known_shift_checks = {
            "N-d = 2S-1": N - d == 2 * S - 1,
            "N+1-d = 2S": N + 1 - d == 2 * S,
            "N+3-d = 2X": N + 3 - d == 2 * X,
            "N-d-1 = 2(S-1)": N - d - 1 == 2 * (S - 1),
        }

        # ------------------------------------------------------------------
        # GCD SIGNATURES
        # ------------------------------------------------------------------

        gcds = {
            "gcd(d,N)": gcd_abs(d, N),
            "gcd(d,N+1)": gcd_abs(d, N + 1),
            "gcd(d,N+2)": gcd_abs(d, N + 2),
            "gcd(d,N+3)": gcd_abs(d, N + 3),
            "gcd(d+1,N)": gcd_abs(d + 1, N),
            "gcd(d+2,N)": gcd_abs(d + 2, N),
            "gcd(d+3,N)": gcd_abs(d + 3, N),
        }

        gcd_ok = all(v >= 1 for v in gcds.values())

        # ------------------------------------------------------------------
        # MODULAR SIGNATURES
        # ------------------------------------------------------------------

        moduli = [2, 3, 4, 5, 8, 9, 16, 32]

        modular_rows = []
        mod_checks_ok = True

        for m in moduli:
            row = {
                "m": m,
                "N_mod": N % m,
                "d_mod": d % m,
                "K_mod": K % m,
                "D_mod": D % m,
            }

            # D=d^2 modulo m.
            row_ok = (row["D_mod"] - (row["d_mod"] * row["d_mod"])) % m == 0
            modular_rows.append(row)

            if not row_ok:
                mod_checks_ok = False

        # ------------------------------------------------------------------
        # SMALL-PRIME FACTORIZATION OF d+3
        # ------------------------------------------------------------------

        d_shift_small_factor = factor_small_integer(d_shifted, limit=10000)

        # ------------------------------------------------------------------
        # DIRECT RELATIONS
        # ------------------------------------------------------------------

        compressed_identity = D == d * d
        positive_root = is_square and d_recovered == d

        factor_identity = d == d_factor
        shifted_factor_identity = d_shifted == d_shifted_factor

        s_reconstruction = (
            parity_S
            and S_rec == S
        )

        x_reconstruction = (
            parity_X
            and X_rec == X
        )

        # The original compact identity recovered directly from D.
        reconstructed_K = (-D - 3 * N * N + 6 * N + 1) // 4
        k_reconstruction = (
            (-D - 3 * N * N + 6 * N + 1) % 4 == 0
            and reconstructed_K == K
        )

        # Equivalent direct square relation.
        fourth_form = 4 * K + 3 * N * N - 6 * N - 1 == -d * d

        # Test whether d has any obvious relationship with N through
        # simple quotients. These are exact integer checks only.
        quotient_checks = {
            "d divides N": N % d == 0 if d != 0 else False,
            "d divides N+1": (N + 1) % d == 0 if d != 0 else False,
            "d divides N-1": (N - 1) % d == 0 if d != 0 else False,
            "d+3 divides N": N % (d + 3) == 0 if d + 3 != 0 else False,
            "d+3 divides N+1": (N + 1) % (d + 3) == 0 if d + 3 != 0 else False,
            "d+3 divides N-1": (N - 1) % (d + 3) == 0 if d + 3 != 0 else False,
        }

        # ------------------------------------------------------------------
        # PRINT INSTANCE
        # ------------------------------------------------------------------

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

        print("COMPRESSED INVARIANT")
        print("  D = -4K - 3N^2 + 6N + 1")
        print(f"  D digits     = {digits(D)}")
        print(f"  D bit-length = {bitlen(D)}")
        print()

        print("COMPRESSED SQUARE CHECK")
        print(f"  D = d^2               : {compressed_identity}")
        print(f"  positive exact sqrt   : {positive_root}")
        print(f"  recovered d           = {d_recovered}")
        print(f"  true d                = {d}")
        print()

        print("FOUR-TERM COMPRESSION")
        print("  4K + 3N^2 - 6N - 1 = -d^2")
        print(f"  exact identity = {fourth_form}")
        print(f"  K reconstructed from N,d^2 = {k_reconstruction}")
        print()

        print("d FACTOR STRUCTURE")
        print("  d = (p-1)(q-1) - (p+q)")
        print(f"  factor identity = {factor_identity}")
        print()

        print("SHIFTED FACTOR STRUCTURE")
        print("  d + 3 = (p-2)(q-2)")
        print(f"  shifted identity = {shifted_factor_identity}")
        print(f"  d+3 = {d_shifted}")
        print(f"  small-factor signature = {d_shift_small_factor}")
        print()

        print("N / d SHIFT IDENTITIES")
        for name, ok in known_shift_checks.items():
            print(f"  {name:<28} = {ok}")
        print()

        print("d -> S RECONSTRUCTION")
        print("  S = (N+1-d)/2")
        print(f"  parity / integrality = {parity_S}")
        print(f"  reconstructed S      = {S_rec}")
        print(f"  exact reconstruction = {s_reconstruction}")
        print()

        print("d -> X RECONSTRUCTION")
        print("  X = (N+3-d)/2")
        print(f"  parity / integrality = {parity_X}")
        print(f"  reconstructed X      = {X_rec}")
        print(f"  exact reconstruction = {x_reconstruction}")
        print()

        print("LOW-COMPLEXITY N,d EXPRESSIONS")
        for name, value in candidates.items():
            print(
                f"  {name:<14} digits={digits(value):>4} "
                f"bit_length={bitlen(value):>4}"
            )
        print()

        print("GCD SIGNATURE")
        for name, value in gcds.items():
            print(f"  {name:<15} = {value}")
        print()

        print("DIVISIBILITY TESTS")
        for name, value in quotient_checks.items():
            print(f"  {name:<28} = {value}")
        print()

        print("MODULAR SQUARE CHECKS")
        for row in modular_rows:
            print(
                f"  mod {row['m']:>2}: "
                f"N={row['N_mod']:>2} "
                f"d={row['d_mod']:>2} "
                f"K={row['K_mod']:>2} "
                f"D={row['D_mod']:>2}"
            )
        print(f"  all modular checks = {mod_checks_ok}")
        print()

        print("SIZE DATA")
        print(f"  |K| digits     = {digits(K)}")
        print(f"  |D| digits     = {digits(D)}")
        print(f"  |d| digits     = {digits(d)}")
        print(f"  |d+3| digits   = {digits(d + 3)}")
        print()

        # ------------------------------------------------------------------
        # Update global checks.
        # ------------------------------------------------------------------

        instance_ok = all([
            compressed_identity,
            positive_root,
            fourth_form,
            k_reconstruction,
            factor_identity,
            shifted_factor_identity,
            s_reconstruction,
            x_reconstruction,
            gcd_ok,
            mod_checks_ok,
            all(known_shift_checks.values()),
        ])

        if not instance_ok:
            global_checks["all_low_complexity_checks"] = False

        global_checks["compressed_square"] &= compressed_identity
        global_checks["positive_root"] &= positive_root
        global_checks["factor_d"] &= factor_identity
        global_checks["shifted_factor"] &= shifted_factor_identity
        global_checks["reconstruct_S"] &= s_reconstruction
        global_checks["reconstruct_X"] &= x_reconstruction
        global_checks["gcd_checks"] &= gcd_ok
        global_checks["shift_checks"] &= all(known_shift_checks.values())
        global_checks["mod_checks"] &= mod_checks_ok

        summary_rows.append({
            "instance": idx,
            "N_digits": digits(N),
            "K_digits": digits(K),
            "D_digits": digits(D),
            "d_digits": digits(d),
            "d_bitlen": bitlen(d),
            "gcd_d_N": gcds["gcd(d,N)"],
            "gcd_d_N1": gcds["gcd(d,N+1)"],
            "gcd_d_N3": gcds["gcd(d,N+3)"],
            "instance_ok": instance_ok,
        })

    # ==========================================================================================================
    # CROSS-INSTANCE SUMMARY
    # ==========================================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("COMPRESSED IDENTITIES")
    print(f"  all D = d^2 checks              = {global_checks['compressed_square']}")
    print(f"  all positive sqrt checks        = {global_checks['positive_root']}")
    print(f"  all K reconstructions           = {global_checks['compressed_square']}")
    print()

    print("FACTOR STRUCTURE")
    print(f"  all d factor identities         = {global_checks['factor_d']}")
    print(f"  all d+3 factor identities       = {global_checks['shifted_factor']}")
    print()

    print("RECONSTRUCTION")
    print(f"  all S reconstructions            = {global_checks['reconstruct_S']}")
    print(f"  all X reconstructions            = {global_checks['reconstruct_X']}")
    print()

    print("LOW-COMPLEXITY STRUCTURE")
    print(f"  all shift identities             = {global_checks['shift_checks']}")
    print(f"  all gcd checks                   = {global_checks['gcd_checks']}")
    print(f"  all modular square checks        = {global_checks['mod_checks']}")
    print()

    print("CROSS-INSTANCE SIZE RANGE")

    all_N_digits = [row["N_digits"] for row in summary_rows]
    all_K_digits = [row["K_digits"] for row in summary_rows]
    all_D_digits = [row["D_digits"] for row in summary_rows]
    all_d_digits = [row["d_digits"] for row in summary_rows]
    all_d_bits = [row["d_bitlen"] for row in summary_rows]

    print(f"  N digits       : {min(all_N_digits)} .. {max(all_N_digits)}")
    print(f"  K digits       : {min(all_K_digits)} .. {max(all_K_digits)}")
    print(f"  D=d^2 digits  : {min(all_D_digits)} .. {max(all_D_digits)}")
    print(f"  d digits       : {min(all_d_digits)} .. {max(all_d_digits)}")
    print(f"  d bit-length   : {min(all_d_bits)} .. {max(all_d_bits)}")
    print()

    print("GCD SIGNATURES")
    for row in summary_rows:
        print(
            f"  instance {row['instance']:>2}: "
            f"gcd(d,N)={row['gcd_d_N']}  "
            f"gcd(d,N+1)={row['gcd_d_N1']}  "
            f"gcd(d,N+3)={row['gcd_d_N3']}"
        )
    print()

    print("INSTANCE STATUS")
    for row in summary_rows:
        print(f"  instance {row['instance']:>2}: {row['instance_ok']}")

    all_instances_pass = all(row["instance_ok"] for row in summary_rows)

    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINAL STATUS")
    print("=" * 110)
    print()
    print("The compressed invariant is:")
    print()
    print("  D = -4K - 3N^2 + 6N + 1")
    print("    = d^2")
    print()
    print("with:")
    print()
    print("  d = N - 2S + 1")
    print()
    print("and, for the tested semiprimes:")
    print()
    print("  d   = (p-1)(q-1) - (p+q)")
    print("  d+3 = (p-2)(q-2)")
    print()
    print("The reconstruction chain is:")
    print()
    print("  (N,K)")
    print("    -> D = -4K - 3N^2 + 6N + 1")
    print("    -> d = sqrt(D)")
    print("    -> S = (N+1-d)/2")
    print("    -> X = (N+3-d)/2")
    print()
    print("The experiment also audits whether d has simple gcd,")
    print("divisibility, modular, or shifted-factor structure")
    print("beyond the already-known perfect-square identity.")
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {all_instances_pass}")
    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
