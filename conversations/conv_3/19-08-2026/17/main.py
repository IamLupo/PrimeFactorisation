import sympy as sp


print("=" * 110)
print("EXPERIMENT 394 START")
print("=" * 110)
print()
print("EXACT X-RECONSTRUCTION / K-ELIMINATION AUDIT")
print()
print("Goal:")
print("  1. Treat X as the primary hidden variable.")
print("  2. Use:")
print("       K = -N^2 + N*X - X^2 + 3X - 2")
print("  3. Verify this identity exactly across all instances.")
print("  4. Substitute K back into Q_X and verify complete collapse.")
print("  5. Determine what information about X remains encoded in K.")
print("  6. Test whether K + N^2 + 2 factorizes naturally in terms of")
print("       X, N, S=p+q.")
print("  7. Test whether simple transforms of K recover S or X.")
print("  8. Search for low-complexity identities involving only N and K.")
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


N_SYM, K_SYM, X_SYM, S_SYM = sp.symbols(
    "N K X S"
)


# =============================================================================
# HELPERS
# =============================================================================

def digits(v):
    return len(str(abs(int(v))))


def bit_length(v):
    return abs(int(v)).bit_length()


def check_zero(expr):
    return sp.expand(expr) == 0


# =============================================================================
# MAIN
# =============================================================================

def main():

    all_k_identity = True
    all_q_collapse = True
    all_factor_identity = True
    all_s_recovery = True
    all_x_recovery = True
    all_low_complexity_checks = True

    recovery_sizes = []

    for idx, (p, q) in enumerate(INSTANCES, start=1):

        p = sp.Integer(p)
        q = sp.Integer(q)

        N = sp.expand(p * q)
        S = sp.expand(p + q)

        # Our established X branch:
        X_TRUE = sp.expand(S + 1)

        K_TRUE = sp.expand(
            1
            - (p**2 - p + 1)
            * (q**2 - q + 1)
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
        print(f"  S       = {S}")
        print(f"  X_TRUE  = {X_TRUE}")
        print(f"  K_TRUE  = {K_TRUE}")
        print()

        # ---------------------------------------------------------------------
        # PRIMARY IDENTITY
        # ---------------------------------------------------------------------

        K_from_X = sp.expand(
            -N**2
            + N * X_TRUE
            - X_TRUE**2
            + 3 * X_TRUE
            - 2
        )

        k_identity = (
            K_from_X == K_TRUE
        )

        print("K-FROM-X IDENTITY")
        print("  K = -N^2 + N*X - X^2 + 3X - 2")
        print(
            "  exact identity = "
            f"{k_identity}"
        )
        print()

        # ---------------------------------------------------------------------
        # NATURAL FACTORIZATION OF K + N^2 + 2
        # ---------------------------------------------------------------------

        K_PLUS = sp.expand(
            K_TRUE + N**2 + 2
        )

        candidate_factor = sp.expand(
            X_TRUE * (N + 3 - X_TRUE)
        )

        factor_identity = (
            K_PLUS == candidate_factor
        )

        print("NATURAL PRODUCT STRUCTURE")
        print(
            "  K + N^2 + 2 = "
            "X * (N+3-X)"
        )
        print(
            "  exact identity = "
            f"{factor_identity}"
        )
        print(
            "  K + N^2 + 2 digits = "
            f"{digits(K_PLUS)}"
        )
        print()

        # ---------------------------------------------------------------------
        # COMPLETE X-QUADRATIC COLLAPSE
        # ---------------------------------------------------------------------

        Q_X = sp.expand(
            X_SYM**2
            - (N + 3) * X_SYM
            + (K_TRUE + N**2 + 2)
        )

        Q_SUB = sp.expand(
            Q_X.subs(
                X_SYM,
                X_TRUE,
            )
        )

        collapse = (
            Q_SUB == 0
        )

        print("X-QUADRATIC")
        print(
            "  Q_X(X) = "
            "X^2 - (N+3)X + (K+N^2+2)"
        )
        print(
            "  Q_X(X_TRUE) = 0 : "
            f"{collapse}"
        )
        print()

        # ---------------------------------------------------------------------
        # RECOVER S FROM X
        # ---------------------------------------------------------------------

        S_FROM_X = sp.expand(
            X_TRUE - 1
        )

        s_recovery = (
            S_FROM_X == S
        )

        # ---------------------------------------------------------------------
        # RECOVER X FROM S
        # ---------------------------------------------------------------------

        X_FROM_S = sp.expand(
            S + 1
        )

        x_recovery = (
            X_FROM_S == X_TRUE
        )

        print("S / X RELATION")
        print(
            "  X = S + 1 : "
            f"{x_recovery}"
        )
        print(
            "  S = X - 1 : "
            f"{s_recovery}"
        )
        print()

        # ---------------------------------------------------------------------
        # ELIMINATE S
        #
        # Since X = S+1:
        #
        # K = -N^2 + N(S+1) - (S+1)^2 + 3(S+1) - 2
        #
        # Expand in S.
        # ---------------------------------------------------------------------

        K_S = sp.expand(
            -N_SYM**2
            + N_SYM * (S_SYM + 1)
            - (S_SYM + 1)**2
            + 3 * (S_SYM + 1)
            - 2
        )

        K_S = sp.expand(
            K_S.subs(
                N_SYM,
                N,
            )
        )

        K_S_EXPECTED = sp.expand(
            -N**2
            + N * S
            - S**2
            + 2 * S
        )

        s_form_identity = (
            sp.expand(
                K_S - K_S_EXPECTED
            )
            == 0
        )

        print("K IN TERMS OF S = p+q")
        print(
            "  K = -N^2 + N*S - S^2 + 2S"
        )
        print(
            "  symbolic collapse = "
            f"{s_form_identity}"
        )
        print()

        # ---------------------------------------------------------------------
        # REARRANGE INTO QUADRATIC IN S
        #
        # S^2 - (N+2)S + (K+N^2) = 0
        # ---------------------------------------------------------------------

        Q_S = sp.expand(
            S_SYM**2
            - (N + 2) * S_SYM
            + (K_TRUE + N**2)
        )

        s_root_check = (
            sp.expand(
                Q_S.subs(
                    S_SYM,
                    S,
                )
            )
            == 0
        )

        S_OTHER = sp.expand(
            N + 2 - S
        )

        s_other_check = (
            sp.expand(
                Q_S.subs(
                    S_SYM,
                    S_OTHER,
                )
            )
            == 0
        )

        print("EXACT S-QUADRATIC")
        print(
            "  Q_S(S) = "
            "S^2 - (N+2)S + (K+N^2)"
        )
        print(
            "  S_TRUE root  = "
            f"{s_root_check}"
        )
        print(
            "  S_OTHER root = "
            f"{s_other_check}"
        )
        print(
            "  S_OTHER = "
            f"{S_OTHER}"
        )
        print()

        # ---------------------------------------------------------------------
        # S DISCRIMINANT
        # ---------------------------------------------------------------------

        Delta_S = sp.expand(
            (N + 2)**2
            - 4 * (K_TRUE + N**2)
        )

        delta_S_expected = sp.expand(
            (S - S_OTHER)**2
        )

        delta_S_check = (
            Delta_S == delta_S_expected
        )

        print("S DISCRIMINANT")
        print(
            "  Delta_S = "
            "(N+2)^2 - 4(K+N^2)"
        )
        print(
            "  Delta_S digits = "
            f"{digits(Delta_S)}"
        )
        print(
            "  Delta_S = (S_TRUE-S_OTHER)^2 : "
            f"{delta_S_check}"
        )
        print()

        # ---------------------------------------------------------------------
        # SIMPLE N,K TRANSFORMS
        #
        # These deliberately test whether S or X can emerge from
        # very low-complexity expressions.
        # ---------------------------------------------------------------------

        candidates = {
            "K + N^2": K_TRUE + N**2,
            "K + N^2 + 1": K_TRUE + N**2 + 1,
            "K + N^2 + 2": K_TRUE + N**2 + 2,
            "N^2 + K - N": K_TRUE + N**2 - N,
            "N^2 + K + N": K_TRUE + N**2 + N,
        }

        print("LOW-COMPLEXITY N,K TRANSFORMS")

        for name, value in candidates.items():

            print(
                f"  {name}: digits = "
                f"{digits(value)}, "
                f"bit_length = {bit_length(value)}"
            )

        print()

        # ---------------------------------------------------------------------
        # DIRECT RECONSTRUCTION CANDIDATES
        #
        # Test whether simple square-root expressions produce S or X.
        # ---------------------------------------------------------------------

        sqrt_Delta_S = sp.expand(
            N + 2 - 2 * S
        )

        sqrt_Delta_X = sp.expand(
            N + 3 - 2 * X_TRUE
        )

        print("EXACT SQUARE-ROOT REPRESENTATIONS")
        print(
            "  sqrt(Delta_S) candidate = "
            f"{sqrt_Delta_S}"
        )
        print(
            "  sqrt(Delta_X) candidate = "
            f"{sqrt_Delta_X}"
        )
        print(
            "  Delta_S square check = "
            f"{sp.expand(Delta_S - sqrt_Delta_S**2) == 0}"
        )
        print(
            "  Delta_X square check = "
            f"{sp.expand((N+3)**2 - 4*(K_TRUE+N**2+2) - sqrt_Delta_X**2) == 0}"
        )
        print()

        # ---------------------------------------------------------------------
        # SIZE DATA
        # ---------------------------------------------------------------------

        recovery_sizes.append(
            (
                digits(K_TRUE),
                digits(K_PLUS),
                digits(Delta_S),
            )
        )

        # ---------------------------------------------------------------------
        # ACCUMULATE
        # ---------------------------------------------------------------------

        all_k_identity &= k_identity
        all_q_collapse &= collapse
        all_factor_identity &= factor_identity
        all_s_recovery &= (
            s_form_identity
            and s_root_check
            and s_other_check
        )
        all_x_recovery &= (
            x_recovery
            and s_recovery
        )
        all_low_complexity_checks &= (
            delta_S_check
        )

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all K-from-X identities        = "
        f"{all_k_identity}"
    )

    print(
        "  all X-quadratic collapses      = "
        f"{all_q_collapse}"
    )

    print(
        "  all K+N^2 factor identities    = "
        f"{all_factor_identity}"
    )

    print(
        "  all S-quadratic checks         = "
        f"{all_s_recovery}"
    )

    print(
        "  all X/S reconstruction checks  = "
        f"{all_x_recovery}"
    )

    print(
        "  all Delta_S checks             = "
        f"{all_low_complexity_checks}"
    )

    print()

    print("SIZE RANGES")

    print(
        "  K minimum digits = "
        f"{min(v[0] for v in recovery_sizes)}"
    )

    print(
        "  K maximum digits = "
        f"{max(v[0] for v in recovery_sizes)}"
    )

    print(
        "  K+N^2+2 minimum digits = "
        f"{min(v[1] for v in recovery_sizes)}"
    )

    print(
        "  K+N^2+2 maximum digits = "
        f"{max(v[1] for v in recovery_sizes)}"
    )

    print(
        "  Delta_S minimum digits = "
        f"{min(v[2] for v in recovery_sizes)}"
    )

    print(
        "  Delta_S maximum digits = "
        f"{max(v[2] for v in recovery_sizes)}"
    )

    print()

    print("=" * 110)
    print("EXPERIMENT 394 FINAL STATUS")
    print("=" * 110)
    print()

    print("The key identity being audited is:")
    print()
    print("  K = -N^2 + N*X - X^2 + 3X - 2")
    print()
    print("Using X = S+1 gives:")
    print()
    print("  K = -N^2 + N*S - S^2 + 2S")
    print()
    print("and therefore:")
    print()
    print("  S^2 - (N+2)S + (K+N^2) = 0")
    print()
    print("This is the same quadratic structure one level lower:")
    print()
    print("  X = S+1")
    print("  S = p+q")
    print()
    print("The research question is now sharper:")
    print()
    print("  Can the smaller root S=p+q, or X=S+1,")
    print("  be reconstructed from N alone through a")
    print("  structurally simpler invariant?")
    print()
    print("In particular, any unexpectedly simple expression")
    print("for Delta_S or K+N^2 is worth investigating.")
    print()

    print("=" * 110)
    print("EXPERIMENT 394 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
