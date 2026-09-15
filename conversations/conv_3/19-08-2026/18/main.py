import sympy as sp


print("=" * 110)
print("EXPERIMENT 395 START")
print("=" * 110)
print()
print("EXACT S-QUADRATIC / CORRECTED K-ELIMINATION")
print()
print("Purpose:")
print("  1. Correct the S = X-1 substitution from Experiment 394.")
print("  2. Verify:")
print("       K = -N^2 + N*S - S^2 + S + N")
print("  3. Verify the exact S quadratic:")
print("       S^2 - (N+1)S + (K+N^2-N) = 0")
print("  4. Recover the second S root.")
print("  5. Verify:")
print("       K + N^2 - N = S*(N-S+1)")
print("  6. Verify the factor form:")
print("       S*(N-S+1) = (p+q)(p-1)(q-1)")
print("  7. Verify the corrected S discriminant.")
print("  8. Verify that the S and X discriminants are identical.")
print("  9. Determine whether S=p+q is always the smaller branch.")
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
    v = int(v)
    return len(str(abs(v)))


def bit_length(v):
    return abs(int(v)).bit_length()


# =============================================================================
# MAIN
# =============================================================================

def main():

    all_k_identity = True
    all_s_quadratic = True
    all_other_root = True
    all_factor_identity = True
    all_discriminant = True
    all_disc_match_x = True
    all_smaller_branch = True
    all_factorized_form = True

    size_rows = []

    for idx, (p_raw, q_raw) in enumerate(INSTANCES, start=1):

        p = sp.Integer(p_raw)
        q = sp.Integer(q_raw)

        N = sp.expand(p * q)
        S = sp.expand(p + q)
        X_TRUE = sp.expand(S + 1)

        K = sp.expand(
            1
            - (p**2 - p + 1)
            * (q**2 - q + 1)
        )

        X_OTHER = sp.expand(
            N + 3 - X_TRUE
        )

        S_OTHER = sp.expand(
            N + 1 - S
        )

        print()
        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {N}")
        print(f"  S_TRUE  = {S}")
        print(f"  X_TRUE  = {X_TRUE}")
        print(f"  K_TRUE  = {K}")
        print()

        # =========================================================================
        # 1. CORRECTED K-IN-S IDENTITY
        # =========================================================================

        K_from_S = sp.expand(
            -N**2
            + N*S
            - S**2
            + S
            + N
        )

        k_identity = (K_from_S == K)

        print("CORRECTED K-IN-S IDENTITY")
        print("  K = -N^2 + N*S - S^2 + S + N")
        print(
            "  exact identity = "
            f"{k_identity}"
        )
        print()

        # =========================================================================
        # 2. NATURAL FACTORIZATION
        #
        # K + N^2 - N = S*(N-S+1)
        # =========================================================================

        constant = sp.expand(
            K + N**2 - N
        )

        natural_product = sp.expand(
            S * (N - S + 1)
        )

        factor_identity = (
            constant == natural_product
        )

        pq_factor = sp.expand(
            (p + q)
            * (p - 1)
            * (q - 1)
        )

        factorized_form = (
            constant == pq_factor
        )

        print("NATURAL CONSTANT FACTORIZATION")
        print("  K + N^2 - N = S*(N-S+1)")
        print(
            "  exact identity = "
            f"{factor_identity}"
        )
        print(
            "  factorized form = "
            "(p+q)(p-1)(q-1)"
        )
        print(
            "  exact factorization = "
            f"{factorized_form}"
        )
        print(
            "  constant digits = "
            f"{digits(constant)}"
        )
        print()

        # =========================================================================
        # 3. CORRECTED S QUADRATIC
        # =========================================================================

        Q_S = sp.expand(
            sp.Symbol("s")**2
            - (N + 1) * sp.Symbol("s")
            + constant
        )

        s = sp.Symbol("s")

        Q_S = sp.expand(
            s**2
            - (N + 1) * s
            + constant
        )

        true_root_check = (
            sp.expand(
                Q_S.subs(s, S)
            ) == 0
        )

        other_root_check = (
            sp.expand(
                Q_S.subs(s, S_OTHER)
            ) == 0
        )

        print("CORRECTED S-QUADRATIC")
        print(
            "  Q_S(S) = "
            "S^2 - (N+1)S + (K+N^2-N)"
        )
        print(
            "  S_TRUE root  = "
            f"{true_root_check}"
        )
        print(
            "  S_OTHER root = "
            f"{other_root_check}"
        )
        print()
        print("S ROOTS")
        print(f"  S_TRUE  = {S}")
        print(f"  S_OTHER = {S_OTHER}")
        print(
            "  S_TRUE < S_OTHER = "
            f"{S < S_OTHER}"
        )
        print()

        # =========================================================================
        # 4. S VIETA
        # =========================================================================

        sum_check = (
            sp.expand(S + S_OTHER)
            == N + 1
        )

        product_check = (
            sp.expand(S * S_OTHER)
            == constant
        )

        print("S VIETA")
        print(
            "  S1 + S2 = N+1 : "
            f"{sum_check}"
        )
        print(
            "  S1 * S2 = K+N^2-N : "
            f"{product_check}"
        )
        print()

        # =========================================================================
        # 5. S DISCRIMINANT
        # =========================================================================

        Delta_S = sp.expand(
            (N + 1)**2
            - 4 * constant
        )

        Delta_S_expected = sp.expand(
            (S_OTHER - S)**2
        )

        discriminant_check = (
            Delta_S == Delta_S_expected
        )

        sqrt_candidate = sp.expand(
            S_OTHER - S
        )

        square_check = (
            sp.expand(
                sqrt_candidate**2
                - Delta_S
            )
            == 0
        )

        print("S DISCRIMINANT")
        print(
            "  Delta_S = "
            "(N+1)^2 - 4(K+N^2-N)"
        )
        print(
            "  Delta_S digits = "
            f"{digits(Delta_S)}"
        )
        print(
            "  Delta_S = (S_OTHER-S_TRUE)^2 : "
            f"{discriminant_check}"
        )
        print(
            "  sqrt candidate digits = "
            f"{digits(sqrt_candidate)}"
        )
        print(
            "  Delta_S = candidate^2 : "
            f"{square_check}"
        )
        print()

        # =========================================================================
        # 6. EXACT FACTOR FORM OF THE S DISCRIMINANT
        #
        # sqrt(Delta_S)
        #   = N - 2S + 1
        #   = (p-1)(q-1) - (p+q)
        # =========================================================================

        factor_sqrt = sp.expand(
            (p - 1) * (q - 1)
            - (p + q)
        )

        simple_sqrt = sp.expand(
            N - 2*S + 1
        )

        factor_sqrt_check = (
            sqrt_candidate == factor_sqrt
        )

        simple_sqrt_check = (
            sqrt_candidate == simple_sqrt
        )

        print("DISCRIMINANT FACTOR FORM")
        print(
            "  sqrt(Delta_S) = "
            "(p-1)(q-1) - (p+q)"
        )
        print(
            "  factor form check = "
            f"{factor_sqrt_check}"
        )
        print(
            "  N - 2S + 1 check = "
            f"{simple_sqrt_check}"
        )
        print()

        # =========================================================================
        # 7. RELATION TO X DISCRIMINANT
        #
        # X = S+1
        #
        # Delta_X must equal Delta_S.
        # =========================================================================

        Delta_X = sp.expand(
            (N + 3)**2
            - 4 * (K + N**2 + 2)
        )

        x_disc_match = (
            Delta_X == Delta_S
        )

        x_root_check = (
            sp.expand(
                X_TRUE**2
                - (N + 3) * X_TRUE
                + (K + N**2 + 2)
            ) == 0
        )

        print("X / S DISCRIMINANT TRANSFER")
        print(
            "  Delta_X = "
            "(N+3)^2 - 4(K+N^2+2)"
        )
        print(
            "  Delta_X == Delta_S : "
            f"{x_disc_match}"
        )
        print(
            "  X_TRUE root check = "
            f"{x_root_check}"
        )
        print()

        # =========================================================================
        # 8. QUADRATIC FORMULA
        # =========================================================================

        S_MINUS = sp.expand(
            (N + 1 - sqrt_candidate) / 2
        )

        S_PLUS = sp.expand(
            (N + 1 + sqrt_candidate) / 2
        )

        formula_minus = (
            S_MINUS == S
        )

        formula_plus = (
            S_PLUS == S_OTHER
        )

        integer_branch = (
            sp.Integer(N + 1 - sqrt_candidate) % 2 == 0
        )

        print("EXACT S QUADRATIC FORMULA")
        print(
            "  S_MINUS = "
            "(N+1-sqrt(Delta_S))/2"
        )
        print(
            "  S_PLUS = "
            "(N+1+sqrt(Delta_S))/2"
        )
        print(
            "  S_MINUS == S_TRUE  : "
            f"{formula_minus}"
        )
        print(
            "  S_PLUS  == S_OTHER : "
            f"{formula_plus}"
        )
        print(
            "  parity / integrality check = "
            f"{integer_branch}"
        )
        print()

        # =========================================================================
        # 9. SMALLER-BRANCH TEST
        # =========================================================================

        smaller_branch = (
            S < S_OTHER
        )

        print("BRANCH TEST")
        print(
            "  S_TRUE < S_OTHER = "
            f"{smaller_branch}"
        )
        print(
            "  X_TRUE < X_OTHER = "
            f"{X_TRUE < X_OTHER}"
        )
        print()

        # =========================================================================
        # 10. DIRECT RELATION BETWEEN S AND X
        # =========================================================================

        x_other_expected = sp.expand(
            S_OTHER + 1
        )

        x_transfer_check = (
            X_OTHER == x_other_expected
        )

        print("S -> X BRANCH TRANSFER")
        print(
            "  X_TRUE = S_TRUE + 1 : "
            f"{X_TRUE == S + 1}"
        )
        print(
            "  X_OTHER = S_OTHER + 1 : "
            f"{x_transfer_check}"
        )
        print()

        # =========================================================================
        # 11. SYMMETRIC PRODUCT STRUCTURE
        # =========================================================================

        product_structure_1 = sp.expand(
            constant
        )

        product_structure_2 = sp.expand(
            S * (N - S + 1)
        )

        product_structure_3 = sp.expand(
            (p + q)
            * ((p - 1) * (q - 1))
        )

        print("THREE-WAY FACTORIZATION CHECK")
        print(
            "  K+N^2-N"
        )
        print(
            "    = S*(N-S+1)"
        )
        print(
            "    = (p+q)(p-1)(q-1)"
        )
        print(
            "  all three exact = "
            f"{product_structure_1 == product_structure_2 == product_structure_3}"
        )
        print()

        # =========================================================================
        # 12. SIZE DATA
        # =========================================================================

        size_rows.append(
            (
                digits(K),
                digits(constant),
                digits(Delta_S),
                digits(sqrt_candidate),
            )
        )

        # =========================================================================
        # 13. GLOBAL ACCUMULATION
        # =========================================================================

        all_k_identity &= k_identity
        all_s_quadratic &= (
            true_root_check
            and sum_check
            and product_check
        )
        all_other_root &= other_root_check
        all_factor_identity &= factor_identity
        all_factorized_form &= factorized_form
        all_discriminant &= (
            discriminant_check
            and square_check
            and factor_sqrt_check
            and simple_sqrt_check
        )
        all_disc_match_x &= (
            x_disc_match
            and x_root_check
        )
        all_smaller_branch &= smaller_branch

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all corrected K-in-S identities      = "
        f"{all_k_identity}"
    )

    print(
        "  all S-quadratic checks               = "
        f"{all_s_quadratic}"
    )

    print(
        "  all S-other-root checks              = "
        f"{all_other_root}"
    )

    print(
        "  all K+N^2-N factor checks            = "
        f"{all_factor_identity}"
    )

    print(
        "  all (p+q)(p-1)(q-1) checks           = "
        f"{all_factorized_form}"
    )

    print(
        "  all corrected S discriminant checks  = "
        f"{all_discriminant}"
    )

    print(
        "  all X/S discriminant transfers       = "
        f"{all_disc_match_x}"
    )

    print(
        "  smaller S branch always true         = "
        f"{all_smaller_branch}"
    )

    print()

    print("SIZE RANGES")

    print(
        "  K minimum digits = "
        f"{min(v[0] for v in size_rows)}"
    )

    print(
        "  K maximum digits = "
        f"{max(v[0] for v in size_rows)}"
    )

    print(
        "  K+N^2-N minimum digits = "
        f"{min(v[1] for v in size_rows)}"
    )

    print(
        "  K+N^2-N maximum digits = "
        f"{max(v[1] for v in size_rows)}"
    )

    print(
        "  Delta_S minimum digits = "
        f"{min(v[2] for v in size_rows)}"
    )

    print(
        "  Delta_S maximum digits = "
        f"{max(v[2] for v in size_rows)}"
    )

    print(
        "  sqrt(Delta_S) minimum digits = "
        f"{min(v[3] for v in size_rows)}"
    )

    print(
        "  sqrt(Delta_S) maximum digits = "
        f"{max(v[3] for v in size_rows)}"
    )

    print()
    print("=" * 110)
    print("EXPERIMENT 395 FINAL STATUS")
    print("=" * 110)
    print()

    print("The corrected substitution is:")
    print()
    print("  X = S + 1")
    print()
    print("which transforms")
    print()
    print("  K = -N^2 + N*X - X^2 + 3X - 2")
    print()
    print("into:")
    print()
    print("  K = -N^2 + N*S - S^2 + S + N")
    print()
    print("Therefore:")
    print()
    print("  S^2 - (N+1)S + (K+N^2-N) = 0")
    print()
    print("and the constant term has the exact factorization:")
    print()
    print("  K + N^2 - N")
    print("      = S(N-S+1)")
    print("      = (p+q)(p-1)(q-1)")
    print()
    print("The discriminant is:")
    print()
    print("  Delta_S")
    print("      = (N+1)^2 - 4(K+N^2-N)")
    print("      = (S_OTHER-S_TRUE)^2")
    print()
    print("This is the corrected S-level form of the same")
    print("two-branch structure found for X.")
    print()
    print("=" * 110)
    print("EXPERIMENT 395 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
