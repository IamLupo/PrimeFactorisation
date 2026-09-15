import sympy as sp


print("=" * 110)
print("EXPERIMENT 396 START")
print("=" * 110)
print()
print("EXACT DISCRIMINANT COMPRESSION / SHIFTED FACTORIZATION")
print()
print("Objective:")
print("  1. Introduce the shifted discriminant variable:")
print("       d = N - 2S + 1")
print("  2. Eliminate S exactly from K.")
print("  3. Verify the compressed identity:")
print("       4K + 3N^2 - 6N - 1 = -d^2")
print("  4. Verify:")
print("       d = (p-1)(q-1) - (p+q)")
print("  5. Test the shifted factor form:")
print("       d + 3 = (p-2)(q-2)")
print("  6. Recover d exactly from N and K.")
print("  7. Recover S exactly from N and d.")
print("  8. Recover X = S+1 exactly.")
print("  9. Verify the complete N,K -> d -> S -> X chain.")
print(" 10. Determine whether the new expression is simpler than")
print("     the previous quadratic/discriminant representation.")
print()
print("=" * 110)


# =============================================================================
# TEST INSTANCES
# =============================================================================

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


# =============================================================================
# HELPERS
# =============================================================================

def digits(v):
    return len(str(abs(int(v))))


def bits(v):
    return abs(int(v)).bit_length()


# =============================================================================
# MAIN
# =============================================================================

def main():

    all_k_compression = True
    all_d_factor = True
    all_shift_factor = True
    all_square_checks = True
    all_d_recovery = True
    all_s_recovery = True
    all_x_recovery = True
    all_chain = True
    all_parity = True

    rows = []

    for idx, (p_raw, q_raw) in enumerate(INSTANCES, start=1):

        p = sp.Integer(p_raw)
        q = sp.Integer(q_raw)

        N = sp.expand(p * q)
        S = sp.expand(p + q)
        X = sp.expand(S + 1)

        K = sp.expand(
            1
            - (p**2 - p + 1)
            * (q**2 - q + 1)
        )

        # -------------------------------------------------------------------------
        # Shifted discriminant variable
        # -------------------------------------------------------------------------

        d = sp.expand(
            N - 2*S + 1
        )

        d_factor = sp.expand(
            (p - 1) * (q - 1)
            - (p + q)
        )

        d_shift = sp.expand(
            d + 3
        )

        shifted_factor = sp.expand(
            (p - 2) * (q - 2)
        )

        # -------------------------------------------------------------------------
        # Compressed K identity
        #
        # K = -(3N^2 - 6N + d^2 - 1)/4
        #
        # equivalently:
        # 4K + 3N^2 - 6N - 1 = -d^2
        # -------------------------------------------------------------------------

        compressed_lhs = sp.expand(
            4*K
            + 3*N**2
            - 6*N
            - 1
        )

        compressed_rhs = sp.expand(
            -d**2
        )

        compression_check = (
            compressed_lhs == compressed_rhs
        )

        # -------------------------------------------------------------------------
        # Recover d^2 from N,K
        # -------------------------------------------------------------------------

        d_squared_from_NK = sp.expand(
            -4*K
            - 3*N**2
            + 6*N
            + 1
        )

        d_square_check = (
            d_squared_from_NK == d**2
        )

        # -------------------------------------------------------------------------
        # Exact d recovery
        # -------------------------------------------------------------------------

        d_recovered = sp.sqrt(d_squared_from_NK)

        # Everything is exact, so verify the positive square root explicitly.
        d_recovered_exact = (
            sp.expand(d_recovered**2 - d**2) == 0
            and d > 0
        )

        # -------------------------------------------------------------------------
        # Recover S from d
        #
        # d = N - 2S + 1
        #
        # S = (N+1-d)/2
        # -------------------------------------------------------------------------

        numerator_S = sp.expand(
            N + 1 - d
        )

        parity_check = (
            numerator_S % 2 == 0
        )

        S_recovered = sp.expand(
            numerator_S / 2
        )

        S_recovery_check = (
            S_recovered == S
        )

        # -------------------------------------------------------------------------
        # Recover X
        # -------------------------------------------------------------------------

        X_recovered = sp.expand(
            S_recovered + 1
        )

        X_recovery_check = (
            X_recovered == X
        )

        # -------------------------------------------------------------------------
        # Verify d-factor relations
        # -------------------------------------------------------------------------

        d_factor_check = (
            d == d_factor
        )

        shifted_factor_check = (
            d_shift == shifted_factor
        )

        # -------------------------------------------------------------------------
        # Reconstruct K directly from N,d
        # -------------------------------------------------------------------------

        K_reconstructed = sp.expand(
            -(
                3*N**2
                - 6*N
                + d**2
                - 1
            ) / 4
        )

        K_reconstruction_check = (
            K_reconstructed == K
        )

        # -------------------------------------------------------------------------
        # Verify equivalent X/S chain
        # -------------------------------------------------------------------------

        X_from_d = sp.expand(
            (N + 3 - d) / 2
        )

        X_from_d_check = (
            X_from_d == X
        )

        # -------------------------------------------------------------------------
        # Verify that d is exactly sqrt(Delta_X) = sqrt(Delta_S)
        # -------------------------------------------------------------------------

        Delta_X = sp.expand(
            (N + 3)**2
            - 4*(K + N**2 + 2)
        )

        Delta_S = sp.expand(
            (N + 1)**2
            - 4*(K + N**2 - N)
        )

        discriminant_match = (
            Delta_X == Delta_S
        )

        discriminant_square = (
            Delta_X == d**2
        )

        # -------------------------------------------------------------------------
        # Verify original X quadratic
        # -------------------------------------------------------------------------

        X_quadratic = sp.expand(
            X**2
            - (N + 3)*X
            + (K + N**2 + 2)
        )

        x_root_check = (
            X_quadratic == 0
        )

        # -------------------------------------------------------------------------
        # Verify corrected S quadratic
        # -------------------------------------------------------------------------

        S_quadratic = sp.expand(
            S**2
            - (N + 1)*S
            + (K + N**2 - N)
        )

        s_root_check = (
            S_quadratic == 0
        )

        # -------------------------------------------------------------------------
        # Print
        # -------------------------------------------------------------------------

        print()
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

        print("SHIFTED DISCRIMINANT VARIABLE")
        print("  d = N - 2S + 1")
        print(f"  d = {d}")
        print(f"  d digits = {digits(d)}")
        print(f"  d bit-length = {bits(d)}")
        print()

        print("K COMPRESSION")
        print("  4K + 3N^2 - 6N - 1 = -d^2")
        print(
            "  exact identity = "
            f"{compression_check}"
        )
        print(
            "  K reconstruction from N,d = "
            f"{K_reconstruction_check}"
        )
        print()

        print("FACTOR REPRESENTATION OF d")
        print(
            "  d = (p-1)(q-1) - (p+q)"
        )
        print(
            "  factor relation = "
            f"{d_factor_check}"
        )
        print()

        print("SHIFTED FACTOR REPRESENTATION")
        print("  d + 3 = (p-2)(q-2)")
        print(
            "  shifted factor relation = "
            f"{shifted_factor_check}"
        )
        print(
            f"  d+3 digits = {digits(d_shift)}"
        )
        print()

        print("N,K -> d^2")
        print("  d^2 = -4K - 3N^2 + 6N + 1")
        print(
            "  exact square identity = "
            f"{d_square_check}"
        )
        print(
            "  positive exact square root = "
            f"{d_recovered_exact}"
        )
        print()

        print("d -> S RECONSTRUCTION")
        print("  S = (N+1-d)/2")
        print(
            "  parity / integrality = "
            f"{parity_check}"
        )
        print(
            "  reconstructed S = "
            f"{S_recovered}"
        )
        print(
            "  S reconstruction = "
            f"{S_recovery_check}"
        )
        print()

        print("d -> X RECONSTRUCTION")
        print("  X = (N+3-d)/2")
        print(
            "  reconstructed X = "
            f"{X_from_d}"
        )
        print(
            "  X reconstruction = "
            f"{X_from_d_check}"
        )
        print()

        print("DISCRIMINANT TRANSFER")
        print(
            "  Delta_X == Delta_S = "
            f"{discriminant_match}"
        )
        print(
            "  Delta_X == d^2 = "
            f"{discriminant_square}"
        )
        print(
            "  d = sqrt(Delta_X) = "
            f"{d}"
        )
        print()

        print("ORIGINAL QUADRATIC CHECKS")
        print(
            "  S quadratic root = "
            f"{s_root_check}"
        )
        print(
            "  X quadratic root = "
            f"{x_root_check}"
        )
        print()

        print("COMPRESSION SIZE DATA")
        print(
            f"  |K| digits                 = {digits(K)}"
        )
        print(
            f"  |4K+3N^2-6N-1| digits      = "
            f"{digits(compressed_lhs)}"
        )
        print(
            f"  d^2 digits                 = "
            f"{digits(d**2)}"
        )
        print(
            f"  d digits                   = "
            f"{digits(d)}"
        )
        print()

        rows.append(
            (
                digits(K),
                digits(d_squared_from_NK),
                digits(d),
                digits(d_shift),
            )
        )

        # -------------------------------------------------------------------------
        # Global state
        # -------------------------------------------------------------------------

        all_k_compression &= compression_check
        all_d_factor &= d_factor_check
        all_shift_factor &= shifted_factor_check
        all_square_checks &= (
            d_square_check
            and d_recovered_exact
        )
        all_d_recovery &= d_square_check
        all_s_recovery &= S_recovery_check
        all_x_recovery &= (
            X_recovery_check
            and X_from_d_check
        )
        all_chain &= (
            compression_check
            and d_factor_check
            and shifted_factor_check
            and S_recovery_check
            and X_from_d_check
            and discriminant_match
            and discriminant_square
            and s_root_check
            and x_root_check
        )
        all_parity &= parity_check

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all compressed K identities       = "
        f"{all_k_compression}"
    )

    print(
        "  all d factor identities           = "
        f"{all_d_factor}"
    )

    print(
        "  all shifted factor identities     = "
        f"{all_shift_factor}"
    )

    print(
        "  all exact square-root checks      = "
        f"{all_square_checks}"
    )

    print(
        "  all d recoveries                  = "
        f"{all_d_recovery}"
    )

    print(
        "  all S recoveries                  = "
        f"{all_s_recovery}"
    )

    print(
        "  all X recoveries                  = "
        f"{all_x_recovery}"
    )

    print(
        "  all parity/integrality checks     = "
        f"{all_parity}"
    )

    print(
        "  complete N,K -> d -> S -> X chain = "
        f"{all_chain}"
    )

    print()

    print("SIZE RANGES")
    print(
        "  K minimum digits = "
        f"{min(r[0] for r in rows)}"
    )
    print(
        "  K maximum digits = "
        f"{max(r[0] for r in rows)}"
    )
    print(
        "  d^2 minimum digits = "
        f"{min(r[1] for r in rows)}"
    )
    print(
        "  d^2 maximum digits = "
        f"{max(r[1] for r in rows)}"
    )
    print(
        "  d minimum digits = "
        f"{min(r[2] for r in rows)}"
    )
    print(
        "  d maximum digits = "
        f"{max(r[2] for r in rows)}"
    )

    print()
    print("=" * 110)
    print("EXPERIMENT 396 FINAL STATUS")
    print("=" * 110)
    print()

    print("The corrected S-level structure admits a further exact")
    print("compression by defining:")
    print()
    print("  d = N - 2S + 1")
    print()
    print("Then:")
    print()
    print("  K = -(3N^2 - 6N + d^2 - 1)/4")
    print()
    print("equivalently:")
    print()
    print("  4K + 3N^2 - 6N - 1 = -d^2")
    print()
    print("and therefore:")
    print()
    print("  d = sqrt(-4K - 3N^2 + 6N + 1)")
    print()
    print("with the exact branch reconstruction:")
    print()
    print("  S = (N+1-d)/2")
    print("  X = (N+3-d)/2")
    print()
    print("For the tested semiprimes:")
    print()
    print("  d = (p-1)(q-1) - (p+q)")
    print("    = (p-2)(q-2) - 3")
    print()
    print("Thus the remaining structural question is whether")
    print("this compressed square relation reveals a simpler")
    print("N,K invariant or merely repackages the same quadratic.")
    print()
    print("=" * 110)
    print("EXPERIMENT 396 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
