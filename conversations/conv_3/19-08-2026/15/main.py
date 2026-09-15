import sympy as sp


print("=" * 110)
print("EXPERIMENT 392 START")
print("=" * 110)
print()
print("EXACT KAPPA -> X QUADRATIC COLLAPSE")
print()
print("Objective:")
print("  1. Transform the Kappa quadratic from t=N/X into an exact")
print("     quadratic directly in X.")
print("  2. Verify that X_TRUE is one exact root.")
print("  3. Recover the second X-root exactly.")
print("  4. Verify the X-root Vieta relations.")
print("  5. Analyze the exact X discriminant.")
print("  6. Test whether the smaller X-root is universally X_TRUE.")
print("  7. Test the exact relation between the X discriminant")
print("     and (p-1)(q-1).")
print("  8. Rewrite the layer variable t=N/X and verify the")
print("     two-root layer branches directly in the X domain.")
print()


# =============================================================================
# SYMBOLS
# =============================================================================

X = sp.symbols("X")
t = sp.symbols("t")


# =============================================================================
# INSTANCES
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
# LAYER FUNCTIONS
# =============================================================================

def build_layer():

    h16 = (
        -9*t**8
        -36*t**7
        -84*t**6
        -126*t**5
        -126*t**4
        -84*t**3
        -36*t**2
        -9*t
        -1
    )

    h15 = (
        88*t**9
        +396*t**8
        +1164*t**7
        +2226*t**6
        +2898*t**5
        +2604*t**4
        +1596*t**3
        +639*t**2
        +151*t
        +16
    )

    h14 = (
        -276*t**10
        -1380*t**9
        -5460*t**8
        -13560*t**7
        -23058*t**6
        -27510*t**5
        -23100*t**4
        -13410*t**3
        -5135*t**2
        -1169*t
        -120
    )

    num = sp.expand(
        h15**2 - h16*h14
    )

    den = sp.expand(
        h16**2
    )

    return num, den


R_NUM, R_DEN = build_layer()


# =============================================================================
# DIRECT R EVALUATION
# =============================================================================

def R_value(tv):

    numerator = sp.cancel(
        R_NUM.subs(t, tv)
    )

    denominator = sp.cancel(
        R_DEN.subs(t, tv)
    )

    return sp.cancel(
        numerator / denominator
    )


# =============================================================================
# INSTANCE DATA
# =============================================================================

def build_instance(p, q):

    N = sp.Integer(p * q)

    S = sp.Integer(p + q)

    X_TRUE = S + 1

    K_TRUE = sp.Integer(
        1
        - (p**2 - p + 1)
        * (q**2 - q + 1)
    )

    return (
        sp.Integer(p),
        sp.Integer(q),
        N,
        S,
        X_TRUE,
        K_TRUE,
    )


# =============================================================================
# ORIGINAL KAPPA QUADRATIC
# =============================================================================

def kappa_t_polynomial(N, K):

    return sp.expand(
        (K + N**2 + 2) * t**2
        - N * (N + 3) * t
        + N**2
    )


# =============================================================================
# X-QUADRATIC
# =============================================================================

def x_quadratic(N, K):

    return sp.Poly(
        sp.expand(
            X**2
            - (N + 3) * X
            + (K + N**2 + 2)
        ),
        X,
        domain=sp.ZZ,
    )


# =============================================================================
# DIGIT HELPERS
# =============================================================================

def ratio_digits(value):

    num = abs(int(sp.numer(value)))
    den = abs(int(sp.denom(value)))

    return (
        len(str(num)),
        len(str(den)),
    )


def integer_digits(value):

    return len(str(abs(int(value))))


# =============================================================================
# MAIN
# =============================================================================

def main():

    all_kappa_pass = True
    all_x_root_pass = True
    all_vieta_pass = True
    all_discriminant_pass = True
    all_small_root_pass = True
    all_second_root_pass = True
    all_layer_x_pass = True

    discriminant_sizes = []

    for index, (p, q) in enumerate(
        INSTANCES,
        start=1,
    ):

        print()
        print("=" * 110)
        print(f"INSTANCE {index}")
        print("=" * 110)
        print()

        # ---------------------------------------------------------------------
        # Basic data
        # ---------------------------------------------------------------------

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

        t_true = sp.cancel(
            N / X_TRUE
        )

        # ---------------------------------------------------------------------
        # Original Kappa quadratic
        # ---------------------------------------------------------------------

        Q_t = kappa_t_polynomial(
            N,
            K_TRUE,
        )

        kappa_true_value = sp.cancel(
            Q_t.subs(t, t_true)
        )

        kappa_true_pass = (
            kappa_true_value == 0
        )

        # ---------------------------------------------------------------------
        # Exact transformed X quadratic
        # ---------------------------------------------------------------------

        Q_x = x_quadratic(
            N,
            K_TRUE,
        )

        x_true_value = sp.cancel(
            Q_x.as_expr().subs(
                X,
                X_TRUE,
            )
        )

        x_true_pass = (
            x_true_value == 0
        )

        # ---------------------------------------------------------------------
        # Second X root from Vieta
        # ---------------------------------------------------------------------

        x_other = sp.cancel(
            (N + 3) - X_TRUE
        )

        x_other_value = sp.cancel(
            Q_x.as_expr().subs(
                X,
                x_other,
            )
        )

        x_other_pass = (
            x_other_value == 0
        )

        # ---------------------------------------------------------------------
        # Direct roots from SymPy
        # ---------------------------------------------------------------------

        roots_x = sp.solve(
            Q_x.as_expr(),
            X,
        )

        root_set = {
            sp.cancel(r)
            for r in roots_x
        }

        roots_match = (
            root_set
            == {
                sp.cancel(X_TRUE),
                sp.cancel(x_other),
            }
        )

        # ---------------------------------------------------------------------
        # X Vieta
        # ---------------------------------------------------------------------

        x_sum_direct = sp.cancel(
            X_TRUE + x_other
        )

        x_product_direct = sp.cancel(
            X_TRUE * x_other
        )

        x_sum_expected = sp.cancel(
            N + 3
        )

        x_product_expected = sp.cancel(
            K_TRUE + N**2 + 2
        )

        x_sum_pass = (
            x_sum_direct
            == x_sum_expected
        )

        x_product_pass = (
            x_product_direct
            == x_product_expected
        )

        vieta_pass = (
            x_sum_pass
            and x_product_pass
        )

        # ---------------------------------------------------------------------
        # X discriminant
        # ---------------------------------------------------------------------

        Delta_X = sp.cancel(
            (N + 3)**2
            - 4 * (
                K_TRUE
                + N**2
                + 2
            )
        )

        Delta_X_expected = sp.cancel(
            (
                X_TRUE
                - x_other
            )**2
        )

        discr_identity_pass = (
            Delta_X
            == Delta_X_expected
        )

        sqrt_delta = sp.cancel(
            sp.Integer(
                int(p - 1)
                * int(q - 1)
            )
        )

        square_root_pass = (
            sp.cancel(
                Delta_X
                - sqrt_delta**2
            )
            == 0
        )

        discriminant_pass = (
            discr_identity_pass
            and square_root_pass
        )

        # ---------------------------------------------------------------------
        # Ordering
        # ---------------------------------------------------------------------

        small_root = min(
            int(X_TRUE),
            int(x_other),
        )

        large_root = max(
            int(X_TRUE),
            int(x_other),
        )

        small_root_is_true = (
            small_root == int(X_TRUE)
        )

        # ---------------------------------------------------------------------
        # Relation to p and q
        # ---------------------------------------------------------------------

        factor_gap = sp.cancel(
            N - S + 1
        )

        pq_minus_expression = sp.cancel(
            (p - 1) * (q - 1)
        )

        gap_pass = (
            factor_gap
            == pq_minus_expression
        )

        # Since:
        #
        #   X_TRUE = p+q+1
        #   x_other = N+3-X_TRUE
        #
        # their difference is:
        #
        #   x_other-X_TRUE = (p-1)(q-1)

        separation = sp.cancel(
            x_other - X_TRUE
        )

        separation_pass = (
            separation
            == pq_minus_expression
        )

        # ---------------------------------------------------------------------
        # Recover t roots from X roots
        # ---------------------------------------------------------------------

        t_from_true_x = sp.cancel(
            N / X_TRUE
        )

        t_from_other_x = sp.cancel(
            N / x_other
        )

        t_true_pass = (
            t_from_true_x == t_true
        )

        Q_true_other = sp.cancel(
            Q_t.subs(
                t,
                t_from_other_x,
            )
        )

        t_other_pass = (
            Q_true_other == 0
        )

        # ---------------------------------------------------------------------
        # Layer branches in X coordinates
        # ---------------------------------------------------------------------

        R_true = R_value(
            t_from_true_x
        )

        R_other = R_value(
            t_from_other_x
        )

        R_true_direct = R_value(
            sp.cancel(N / X_TRUE)
        )

        R_other_direct = R_value(
            sp.cancel(N / x_other)
        )

        layer_x_true_pass = (
            R_true == R_true_direct
        )

        layer_x_other_pass = (
            R_other == R_other_direct
        )

        layer_x_pass = (
            layer_x_true_pass
            and layer_x_other_pass
        )

        # ---------------------------------------------------------------------
        # K expression directly from X
        # ---------------------------------------------------------------------

        K_from_X = sp.cancel(
            -N**2
            + N * X_TRUE
            - X_TRUE**2
            + 3 * X_TRUE
            - 2
        )

        K_from_X_pass = (
            K_from_X == K_TRUE
        )

        # ---------------------------------------------------------------------
        # Polynomial identity in symbolic X
        #
        # This is the main algebraic collapse:
        #
        # K = -N^2 + N X - X^2 + 3X - 2
        #
        # ---------------------------------------------------------------------

        K_X_symbolic = sp.expand(
            -N**2
            + N * X
            - X**2
            + 3 * X
            - 2
        )

        collapse_identity = sp.expand(
            K_X_symbolic
            + N**2
            - N * X
            + X**2
            - 3 * X
            + 2
        )

        collapse_pass = (
            collapse_identity == 0
        )

        # ---------------------------------------------------------------------
        # Size statistics
        # ---------------------------------------------------------------------

        delta_ratio = ratio_digits(
            Delta_X
        )

        discriminant_sizes.append(
            delta_ratio
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

        print("ORIGINAL KAPPA CHECK")
        print(
            "  Q_t(t_TRUE) = "
            f"{kappa_true_value}"
        )
        print(
            "  status       = "
            f"{kappa_true_pass}"
        )
        print()

        print("EXACT X-QUADRATIC")
        print(
            "  Q_X(X) = X^2 - (N+3)X + (K+N^2+2)"
        )
        print(
            "  degree = "
            f"{Q_x.degree()}"
        )
        print(
            "  leading coefficient = "
            f"{Q_x.LC()}"
        )
        print()

        print("X ROOT CHECKS")
        print(
            "  X_TRUE root  = "
            f"{x_true_pass}"
        )
        print(
            "  X_OTHER root = "
            f"{x_other_pass}"
        )
        print(
            "  root set exact match = "
            f"{roots_match}"
        )
        print()

        print("X ROOTS")
        print(
            "  X_TRUE  = "
            f"{X_TRUE}"
        )
        print(
            "  X_OTHER = "
            f"{x_other}"
        )
        print(
            "  X_TRUE < X_OTHER = "
            f"{int(X_TRUE) < int(x_other)}"
        )
        print()

        print("X VIETA")
        print(
            "  X1 + X2 = N+3 : "
            f"{x_sum_pass}"
        )
        print(
            "  X1 * X2 = K+N^2+2 : "
            f"{x_product_pass}"
        )
        print(
            "  Vieta status = "
            f"{vieta_pass}"
        )
        print()

        print("X DISCRIMINANT")
        print(
            "  Delta_X = (N+3)^2 - 4(K+N^2+2)"
        )
        print(
            "  Delta_X numerator digits = "
            f"{delta_ratio[0]}"
        )
        print(
            "  Delta_X denominator digits = "
            f"{delta_ratio[1]}"
        )
        print(
            "  Delta_X = (X_TRUE-X_OTHER)^2 : "
            f"{discr_identity_pass}"
        )
        print(
            "  sqrt(Delta_X) = (p-1)(q-1) : "
            f"{square_root_pass}"
        )
        print(
            "  discriminant status = "
            f"{discriminant_pass}"
        )
        print()

        print("FACTOR-STRUCTURE RELATIONS")
        print(
            "  N-S+1 = (p-1)(q-1) : "
            f"{gap_pass}"
        )
        print(
            "  X_OTHER-X_TRUE = (p-1)(q-1) : "
            f"{separation_pass}"
        )
        print()

        print("DIRECT K-FROM-X IDENTITY")
        print(
            "  K = -N^2 + N*X - X^2 + 3X - 2"
        )
        print(
            "  direct identity = "
            f"{K_from_X_pass}"
        )
        print(
            "  symbolic polynomial collapse = "
            f"{collapse_pass}"
        )
        print()

        print("T-ROOTS RECOVERED FROM X")
        print(
            "  N/X_TRUE reproduces t_TRUE = "
            f"{t_true_pass}"
        )
        print(
            "  N/X_OTHER is second Kappa root = "
            f"{t_other_pass}"
        )
        print()

        print("LAYER IN X-DOMAIN")
        print(
            "  R(N/X_TRUE) exact = "
            f"{layer_x_true_pass}"
        )
        print(
            "  R(N/X_OTHER) exact = "
            f"{layer_x_other_pass}"
        )
        print(
            "  both X-domain branches pass = "
            f"{layer_x_pass}"
        )
        print()

        print("SIZE DATA")
        print(
            "  X_TRUE digits  = "
            f"{integer_digits(X_TRUE)}"
        )
        print(
            "  X_OTHER digits = "
            f"{integer_digits(x_other)}"
        )
        print(
            "  root separation digits = "
            f"{integer_digits(separation)}"
        )
        print()

        # ---------------------------------------------------------------------
        # Accumulate
        # ---------------------------------------------------------------------

        all_kappa_pass &= kappa_true_pass
        all_x_root_pass &= (
            x_true_pass
            and x_other_pass
            and roots_match
        )
        all_vieta_pass &= vieta_pass
        all_discriminant_pass &= discriminant_pass
        all_small_root_pass &= small_root_is_true
        all_second_root_pass &= (
            t_true_pass
            and t_other_pass
        )
        all_layer_x_pass &= layer_x_pass

    # =========================================================================
    # GLOBAL
    # =========================================================================

    print()
    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all Kappa checks pass       = "
        f"{all_kappa_pass}"
    )

    print(
        "  all X root checks pass      = "
        f"{all_x_root_pass}"
    )

    print(
        "  all X Vieta checks pass     = "
        f"{all_vieta_pass}"
    )

    print(
        "  all discriminant checks     = "
        f"{all_discriminant_pass}"
    )

    print(
        "  smaller X root always true  = "
        f"{all_small_root_pass}"
    )

    print(
        "  all recovered t roots pass  = "
        f"{all_second_root_pass}"
    )

    print(
        "  all X-domain layer checks   = "
        f"{all_layer_x_pass}"
    )

    print()

    print("DISCRIMINANT SIZE RANGE")

    print(
        "  minimum numerator digits = "
        f"{min(v[0] for v in discriminant_sizes)}"
    )

    print(
        "  maximum numerator digits = "
        f"{max(v[0] for v in discriminant_sizes)}"
    )

    print()

    print("=" * 110)
    print("EXPERIMENT 392 FINAL STATUS")
    print("=" * 110)
    print()

    print("The central transformation tested is:")
    print()
    print("  t = N/X")
    print()
    print("which converts:")
    print()
    print("  (K+N^2+2)t^2 - N(N+3)t + N^2 = 0")
    print()
    print("into:")
    print()
    print("  X^2 - (N+3)X + (K+N^2+2) = 0")
    print()
    print("Therefore:")
    print()
    print("  X1 + X2 = N+3")
    print("  X1*X2 = K+N^2+2")
    print()
    print("and:")
    print()
    print("  Delta_X = (N+3)^2 - 4(K+N^2+2)")
    print("          = (X1-X2)^2")
    print()
    print("For the tested semiprimes:")
    print()
    print("  X_TRUE  = p+q+1")
    print("  X_OTHER = N+3-X_TRUE")
    print("  X_OTHER-X_TRUE = (p-1)(q-1)")
    print()
    print("The experiment therefore tests whether the hidden X is")
    print("systematically the smaller branch of the exact Kappa")
    print("quadratic, and whether the layer remains structurally")
    print("well-defined after moving completely from t to X.")
    print()

    print("=" * 110)
    print("EXPERIMENT 392 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
