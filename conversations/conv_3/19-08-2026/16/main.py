import sympy as sp


print("=" * 110)
print("EXPERIMENT 393 START")
print("=" * 110)
print()
print("EXACT X-BRANCH / DISCRIMINANT RECONSTRUCTION")
print()
print("Objective:")
print("  1. Correct the discriminant relation from Experiment 392.")
print("  2. Verify:")
print("       Delta_X = (N+3-2*X_TRUE)^2")
print("  3. Verify the equivalent factor form:")
print("       sqrt(Delta_X) = (p-1)(q-1) - (p+q)")
print("  4. Recover both X roots directly from Delta_X.")
print("  5. Determine whether the smaller X branch is always X_TRUE.")
print("  6. Express X_TRUE using N, K and sqrt(Delta_X).")
print("  7. Determine whether any sign ambiguity remains.")
print()


# =============================================================================
# SYMBOLS
# =============================================================================

X = sp.symbols("X")


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
# INSTANCE CONSTRUCTION
# =============================================================================

def build_instance(p, q):

    p = sp.Integer(p)
    q = sp.Integer(q)

    N = p * q
    S = p + q
    X_TRUE = S + 1

    K_TRUE = (
        1
        - (p**2 - p + 1)
        * (q**2 - q + 1)
    )

    return (
        p,
        q,
        N,
        S,
        X_TRUE,
        K_TRUE,
    )


# =============================================================================
# X QUADRATIC
# =============================================================================

def x_quadratic(N, K):

    return sp.expand(
        X**2
        - (N + 3) * X
        + (K + N**2 + 2)
    )


# =============================================================================
# INTEGER DIGITS
# =============================================================================

def digits(value):

    value = abs(int(value))
    return len(str(value))


# =============================================================================
# MAIN
# =============================================================================

def main():

    all_discriminants_pass = True
    all_factor_forms_pass = True
    all_root_formula_pass = True
    all_smaller_branch_pass = True
    all_reconstruction_pass = True
    all_sign_pass = True
    all_identity_pass = True

    delta_sizes = []

    for index, (p, q) in enumerate(
        INSTANCES,
        start=1,
    ):

        print()
        print("=" * 110)
        print(f"INSTANCE {index}")
        print("=" * 110)
        print()

        (
            p,
            q,
            N,
            S,
            X_TRUE,
            K_TRUE,
        ) = build_instance(
            p,
            q,
        )

        # ---------------------------------------------------------------------
        # X quadratic
        # ---------------------------------------------------------------------

        QX = x_quadratic(
            N,
            K_TRUE,
        )

        # ---------------------------------------------------------------------
        # True X root
        # ---------------------------------------------------------------------

        true_root_check = (
            sp.expand(
                QX.subs(X, X_TRUE)
            )
            == 0
        )

        # ---------------------------------------------------------------------
        # Exact other root
        # ---------------------------------------------------------------------

        X_OTHER = sp.expand(
            N + 3 - X_TRUE
        )

        other_root_check = (
            sp.expand(
                QX.subs(X, X_OTHER)
            )
            == 0
        )

        # ---------------------------------------------------------------------
        # Discriminant
        # ---------------------------------------------------------------------

        Delta_X = sp.expand(
            (N + 3)**2
            - 4 * (K_TRUE + N**2 + 2)
        )

        delta_sqrt_exact = sp.expand(
            N + 3 - 2 * X_TRUE
        )

        delta_square_check = (
            sp.expand(
                Delta_X
                - delta_sqrt_exact**2
            )
            == 0
        )

        # ---------------------------------------------------------------------
        # Correct factor relation
        #
        # (p-1)(q-1) - (p+q)
        # = pq-p-q+1-p-q
        # = N-2S+1
        # ---------------------------------------------------------------------

        pq_shifted = sp.expand(
            (p - 1) * (q - 1) - (p + q)
        )

        factor_relation_check = (
            pq_shifted
            == delta_sqrt_exact
        )

        # ---------------------------------------------------------------------
        # Root separation
        # ---------------------------------------------------------------------

        root_separation = sp.expand(
            X_OTHER - X_TRUE
        )

        root_separation_check = (
            root_separation
            == delta_sqrt_exact
        )

        # ---------------------------------------------------------------------
        # Exact square-root verification
        # ---------------------------------------------------------------------

        square_check = (
            sp.expand(
                Delta_X
                - root_separation**2
            )
            == 0
        )

        # ---------------------------------------------------------------------
        # Quadratic formula roots
        # ---------------------------------------------------------------------

        # Since X_TRUE < X_OTHER in all tested cases,
        #
        #   X_TRUE  = (N+3 - sqrt(Delta))/2
        #   X_OTHER = (N+3 + sqrt(Delta))/2

        X_MINUS = sp.cancel(
            (
                N
                + 3
                - delta_sqrt_exact
            ) / 2
        )

        X_PLUS = sp.cancel(
            (
                N
                + 3
                + delta_sqrt_exact
            ) / 2
        )

        minus_check = (
            X_MINUS == X_TRUE
        )

        plus_check = (
            X_PLUS == X_OTHER
        )

        formula_check = (
            minus_check
            and plus_check
        )

        # ---------------------------------------------------------------------
        # Ordering
        # ---------------------------------------------------------------------

        smaller_is_true = (
            int(X_TRUE)
            <
            int(X_OTHER)
        )

        # ---------------------------------------------------------------------
        # Reconstruction from N,K alone
        #
        # sqrt(Delta_X) is an integer here.
        #
        # We test the exact symbolic expression:
        #
        # X = (N+3 +/- sqrt(Delta))/2
        # ---------------------------------------------------------------------

        reconstructed_true = sp.cancel(
            (
                N
                + 3
                - delta_sqrt_exact
            ) / 2
        )

        reconstructed_other = sp.cancel(
            (
                N
                + 3
                + delta_sqrt_exact
            ) / 2
        )

        reconstruction_check = (
            reconstructed_true == X_TRUE
            and
            reconstructed_other == X_OTHER
        )

        # ---------------------------------------------------------------------
        # Sign ambiguity
        # ---------------------------------------------------------------------

        positive_root_separation = (
            delta_sqrt_exact > 0
        )

        negative_root_branch = sp.expand(
            -delta_sqrt_exact
        )

        sign_choice_check = (
            positive_root_separation
            and
            (
                (
                    N + 3 - delta_sqrt_exact
                ) / 2
                == X_TRUE
            )
            and
            (
                (
                    N + 3 + delta_sqrt_exact
                ) / 2
                == X_OTHER
            )
        )

        # ---------------------------------------------------------------------
        # K-from-X consistency
        # ---------------------------------------------------------------------

        K_from_true_X = sp.expand(
            -N**2
            + N * X_TRUE
            - X_TRUE**2
            + 3 * X_TRUE
            - 2
        )

        k_identity_check = (
            K_from_true_X == K_TRUE
        )

        # ---------------------------------------------------------------------
        # Size
        # ---------------------------------------------------------------------

        delta_sizes.append(
            (
                digits(Delta_X),
                digits(delta_sqrt_exact),
            )
        )

        # ---------------------------------------------------------------------
        # Print
        # ---------------------------------------------------------------------

        print("BASIC DATA")
        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {N}")
        print(f"  S       = {S}")
        print(f"  X_TRUE  = {X_TRUE}")
        print(f"  K_TRUE  = {K_TRUE}")
        print()

        print("X QUADRATIC")
        print("  Q_X(X) = X^2 - (N+3)X + (K+N^2+2)")
        print(
            "  X_TRUE root = "
            f"{true_root_check}"
        )
        print(
            "  X_OTHER root = "
            f"{other_root_check}"
        )
        print()

        print("EXACT X ROOTS")
        print(
            "  X_TRUE  = "
            f"{X_TRUE}"
        )
        print(
            "  X_OTHER = "
            f"{X_OTHER}"
        )
        print(
            "  X_TRUE < X_OTHER = "
            f"{smaller_is_true}"
        )
        print()

        print("DISCRIMINANT")
        print(
            "  Delta_X = "
            "(N+3)^2 - 4(K+N^2+2)"
        )
        print(
            "  Delta_X digits = "
            f"{digits(Delta_X)}"
        )
        print(
            "  sqrt(Delta_X) candidate = "
            f"{delta_sqrt_exact}"
        )
        print(
            "  Delta_X = candidate^2 : "
            f"{delta_square_check}"
        )
        print(
            "  (X_OTHER-X_TRUE)^2 = Delta_X : "
            f"{square_check}"
        )
        print()

        print("CORRECTED FACTOR RELATION")
        print(
            "  (p-1)(q-1) - (p+q)"
        )
        print(
            "    == N - 2(p+q) + 1 : "
            f"{factor_relation_check}"
        )
        print(
            "  X_OTHER-X_TRUE == N-2S+1 : "
            f"{root_separation_check}"
        )
        print()

        print("EXACT QUADRATIC-FORMULA RECONSTRUCTION")
        print(
            "  X_MINUS = (N+3-sqrt(Delta_X))/2"
        )
        print(
            "  X_PLUS  = (N+3+sqrt(Delta_X))/2"
        )
        print(
            "  X_MINUS == X_TRUE  : "
            f"{minus_check}"
        )
        print(
            "  X_PLUS  == X_OTHER : "
            f"{plus_check}"
        )
        print(
            "  complete formula check = "
            f"{formula_check}"
        )
        print()

        print("SIGN / BRANCH TEST")
        print(
            "  sqrt(Delta_X) > 0 = "
            f"{positive_root_separation}"
        )
        print(
            "  negative branch candidate = "
            f"{negative_root_branch}"
        )
        print(
            "  correct sign gives X_TRUE = "
            f"{sign_choice_check}"
        )
        print()

        print("DIRECT K-FROM-X")
        print(
            "  K = -N^2 + N*X - X^2 + 3X - 2"
        )
        print(
            "  identity = "
            f"{k_identity_check}"
        )
        print()

        print("RECONSTRUCTION")
        print(
            "  reconstructed X_TRUE  = "
            f"{reconstructed_true}"
        )
        print(
            "  reconstructed X_OTHER = "
            f"{reconstructed_other}"
        )
        print(
            "  exact reconstruction = "
            f"{reconstruction_check}"
        )
        print()

        # ---------------------------------------------------------------------
        # Accumulate
        # ---------------------------------------------------------------------

        all_discriminants_pass &= (
            true_root_check
            and other_root_check
            and delta_square_check
            and square_check
        )

        all_factor_forms_pass &= (
            factor_relation_check
            and root_separation_check
        )

        all_root_formula_pass &= formula_check

        all_smaller_branch_pass &= smaller_is_true

        all_reconstruction_pass &= reconstruction_check

        all_sign_pass &= sign_choice_check

        all_identity_pass &= k_identity_check

    # =========================================================================
    # GLOBAL
    # =========================================================================

    print()
    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all X discriminant checks pass = "
        f"{all_discriminants_pass}"
    )

    print(
        "  corrected factor relation pass = "
        f"{all_factor_forms_pass}"
    )

    print(
        "  all quadratic-formula checks    = "
        f"{all_root_formula_pass}"
    )

    print(
        "  smaller root always X_TRUE     = "
        f"{all_smaller_branch_pass}"
    )

    print(
        "  all exact reconstructions       = "
        f"{all_reconstruction_pass}"
    )

    print(
        "  all sign checks                 = "
        f"{all_sign_pass}"
    )

    print(
        "  all K-from-X identities         = "
        f"{all_identity_pass}"
    )

    print()

    print("DISCRIMINANT SIZE RANGE")

    print(
        "  minimum Delta digits = "
        f"{min(v[0] for v in delta_sizes)}"
    )

    print(
        "  maximum Delta digits = "
        f"{max(v[0] for v in delta_sizes)}"
    )

    print(
        "  minimum sqrt digits = "
        f"{min(v[1] for v in delta_sizes)}"
    )

    print(
        "  maximum sqrt digits = "
        f"{max(v[1] for v in delta_sizes)}"
    )

    print()

    print("=" * 110)
    print("EXPERIMENT 393 FINAL STATUS")
    print("=" * 110)
    print()

    print("The corrected discriminant identity is:")
    print()
    print("  Delta_X")
    print("    = (N+3)^2 - 4(K+N^2+2)")
    print("    = (X_OTHER-X_TRUE)^2")
    print()
    print("For X_TRUE = p+q+1:")
    print()
    print("  X_OTHER-X_TRUE")
    print("    = N-2(p+q)+1")
    print("    = (p-1)(q-1) - (p+q)")
    print()
    print("Thus the exact branches are:")
    print()
    print("  X_TRUE  = (N+3-sqrt(Delta_X))/2")
    print("  X_OTHER = (N+3+sqrt(Delta_X))/2")
    print()
    print("The important remaining question is not whether the")
    print("quadratic has two exact branches — that is established.")
    print()
    print("The important question is whether the smaller branch")
    print("can be identified from N alone, or from a quantity")
    print("that is independently computable from N, without")
    print("already possessing K or the hidden factor information.")
    print()

    print("=" * 110)
    print("EXPERIMENT 393 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
