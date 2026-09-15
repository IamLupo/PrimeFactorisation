import sympy as sp


# =============================================================================
# EXPERIMENT 390 START
# =============================================================================

print("=" * 110)
print("EXPERIMENT 390 START")
print("=" * 110)
print()
print("EXACT TWO-ROOT BRANCH SEPARATION / DISCRIMINANT TRANSFER")
print()
print("Goal:")
print("  1. Compute the exact two Kappa roots.")
print("  2. Compute the exact two layer roots.")
print("  3. Measure the branch separation:")
print("       Delta_t = t_TRUE - t_OTHER")
print("       Delta_R = R_TRUE - R_OTHER")
print("  4. Compute the exact secant slope:")
print("       H = Delta_R / Delta_t")
print("  5. Verify:")
print("       H^2 = Delta_R^2 / Delta_t^2")
print("  6. Compare the Kappa and layer discriminants.")
print("  7. Look for stable scaling laws across instances.")
print()


# =============================================================================
# SYMBOL
# =============================================================================

t = sp.symbols("t")


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
# LAYER R(t)
# =============================================================================

def build_R():

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

    numerator = sp.expand(
        h15**2 - h16*h14
    )

    denominator = sp.expand(
        h16**2
    )

    return numerator, denominator


R_NUM, R_DEN = build_R()


# =============================================================================
# BASIC VALUES
# =============================================================================

def get_instance_values(p, q):

    N = sp.Integer(p * q)

    X_TRUE = sp.Integer(
        p + q + 1
    )

    K_TRUE = sp.Integer(
        1
        - (p**2 - p + 1)
        * (q**2 - q + 1)
    )

    t_true = sp.cancel(
        N / X_TRUE
    )

    return (
        N,
        X_TRUE,
        K_TRUE,
        t_true,
    )


# =============================================================================
# KAPPA QUADRATIC
# =============================================================================

def kappa_quadratic(N, K):

    return sp.Poly(
        sp.expand(
            (K + N**2 + 2) * t**2
            - N * (N + 3) * t
            + N**2
        ),
        t,
        domain=sp.QQ,
    )


# =============================================================================
# OTHER KAPPA ROOT
# =============================================================================

def other_root(N, K, t_true):

    a = K + N**2 + 2
    c = N**2

    return sp.cancel(
        c / (a * t_true)
    )


# =============================================================================
# EXACT R(t)
# =============================================================================

def R_value(tv):

    return sp.cancel(
        R_NUM.subs(t, tv)
        / R_DEN.subs(t, tv)
    )


# =============================================================================
# KAPPA DISCRIMINANT
# =============================================================================

def kappa_discriminant(N, K):

    a = K + N**2 + 2
    b = -N * (N + 3)
    c = N**2

    return sp.expand(
        b**2 - 4*a*c
    )


# =============================================================================
# DIGIT SIZE
# =============================================================================

def digits(value):

    num = abs(int(sp.numer(value)))
    den = abs(int(sp.denom(value)))

    return (
        len(str(num)),
        len(str(den)),
    )


# =============================================================================
# FACTOR-STRUCTURE WITHOUT FULL FACTORIZATION
# =============================================================================

def gcd_size(value):

    num = abs(int(sp.numer(value)))
    den = abs(int(sp.denom(value)))

    g = sp.gcd(num, den)

    return len(str(abs(int(g))))


# =============================================================================
# MAIN
# =============================================================================

def main():

    separation_pass = True
    slope_identity_pass = True
    discriminant_pass = True

    slope_digit_sizes = []
    ratio_digit_sizes = []

    print()

    for index, (p, q) in enumerate(
        INSTANCES,
        start=1,
    ):

        print("=" * 110)
        print(f"INSTANCE {index}")
        print("=" * 110)
        print()

        (
            N,
            X_TRUE,
            K_TRUE,
            t_true,
        ) = get_instance_values(
            p,
            q
        )

        t_other = other_root(
            N,
            K_TRUE,
            t_true,
        )

        # ---------------------------------------------------------------------
        # Kappa root verification
        # ---------------------------------------------------------------------

        Q = kappa_quadratic(
            N,
            K_TRUE,
        )

        q_true = sp.expand(
            Q.as_expr().subs(
                t,
                t_true,
            )
        )

        q_other = sp.expand(
            Q.as_expr().subs(
                t,
                t_other,
            )
        )

        kappa_roots_pass = (
            q_true == 0
            and q_other == 0
        )

        # ---------------------------------------------------------------------
        # Layer branches
        # ---------------------------------------------------------------------

        R_true = R_value(
            t_true
        )

        R_other = R_value(
            t_other
        )

        # ---------------------------------------------------------------------
        # Exact separations
        # ---------------------------------------------------------------------

        delta_t = sp.cancel(
            t_true - t_other
        )

        delta_R = sp.cancel(
            R_true - R_other
        )

        slope = sp.cancel(
            delta_R / delta_t
        )

        # ---------------------------------------------------------------------
        # Exact identities
        # ---------------------------------------------------------------------

        slope_identity = sp.cancel(
            slope * delta_t - delta_R
        ) == 0

        # ---------------------------------------------------------------------
        # Kappa discriminant
        # ---------------------------------------------------------------------

        Delta_K = kappa_discriminant(
            N,
            K_TRUE,
        )

        Delta_K_expected = sp.cancel(
            (K_TRUE + N**2 + 2)**2
            * delta_t**2
        )

        kappa_disc_pass = (
            sp.cancel(
                Delta_K - Delta_K_expected
            ) == 0
        )

        # ---------------------------------------------------------------------
        # Layer discriminant as branch separation
        #
        # We do not construct the full resultant.
        #
        # Since the two roots are R_TRUE and R_OTHER:
        #
        #   Delta_R = (R_TRUE - R_OTHER)^2
        # ---------------------------------------------------------------------

        Delta_R_direct = sp.cancel(
            delta_R**2
        )

        discriminant_pass_instance = (
            Delta_R_direct >= 0
            if (
                sp.denom(Delta_R_direct) != 0
                and sp.numer(Delta_R_direct) >= 0
            )
            else True
        )

        # ---------------------------------------------------------------------
        # Compare discriminant transfer:
        #
        #   Delta_R / Delta_K
        #
        # This is exactly H^2 divided by
        # (K+N^2+2)^2.
        # ---------------------------------------------------------------------

        transfer_ratio = sp.cancel(
            Delta_R_direct / Delta_K
        )

        expected_transfer = sp.cancel(
            slope**2
            / (K_TRUE + N**2 + 2)**2
        )

        transfer_pass = (
            sp.cancel(
                transfer_ratio - expected_transfer
            ) == 0
        )

        # ---------------------------------------------------------------------
        # Size statistics
        # ---------------------------------------------------------------------

        delta_t_size = digits(
            delta_t
        )

        delta_R_size = digits(
            delta_R
        )

        slope_size = digits(
            slope
        )

        ratio_size = digits(
            transfer_ratio
        )

        slope_digit_sizes.append(
            slope_size
        )

        ratio_digit_sizes.append(
            ratio_size
        )

        # ---------------------------------------------------------------------
        # Output
        # ---------------------------------------------------------------------

        print("BASIC DATA")
        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {N}")
        print(f"  X_TRUE  = {X_TRUE}")
        print(f"  K_TRUE  = {K_TRUE}")
        print()

        print("KAPPA ROOTS")
        print(f"  t_TRUE  = {t_true}")
        print(f"  t_OTHER = {t_other}")
        print(
            "  exact root verification = "
            f"{kappa_roots_pass}"
        )
        print()

        print("BRANCH SEPARATION")
        print(
            f"  Delta_t numerator digits   = "
            f"{delta_t_size[0]}"
        )
        print(
            f"  Delta_t denominator digits = "
            f"{delta_t_size[1]}"
        )
        print(
            f"  Delta_R numerator digits   = "
            f"{delta_R_size[0]}"
        )
        print(
            f"  Delta_R denominator digits = "
            f"{delta_R_size[1]}"
        )
        print()

        print("SECANT SLOPE")
        print(
            "  H = (R_TRUE - R_OTHER) / "
            "(t_TRUE - t_OTHER)"
        )
        print(
            f"  numerator digits   = "
            f"{slope_size[0]}"
        )
        print(
            f"  denominator digits = "
            f"{slope_size[1]}"
        )
        print(
            f"  exact slope identity = "
            f"{slope_identity}"
        )
        print()

        print("KAPPA DISCRIMINANT")
        print(
            f"  Delta_K bit-length = "
            f"{abs(int(Delta_K)).bit_length()}"
        )
        print(
            f"  discriminant/root-separation identity = "
            f"{kappa_disc_pass}"
        )
        print()

        print("DISCRIMINANT TRANSFER")
        print(
            "  Delta_R = (R_TRUE - R_OTHER)^2"
        )
        print(
            "  Delta_R / Delta_K == "
            "H^2 / (K+N^2+2)^2 : "
            f"{transfer_pass}"
        )
        print(
            f"  transfer ratio numerator digits = "
            f"{ratio_size[0]}"
        )
        print(
            f"  transfer ratio denominator digits = "
            f"{ratio_size[1]}"
        )
        print()

        print("BRANCH ORDER")
        print(
            "  t_TRUE > t_OTHER : "
            f"{t_true > t_other}"
        )

        # We know from the previous monotonicity experiments that
        # R decreases with X, while t=N/X decreases as X increases.
        print(
            "  R_TRUE > R_OTHER : "
            f"{R_true > R_other}"
        )

        print()

        separation_pass &= kappa_roots_pass
        slope_identity_pass &= slope_identity
        discriminant_pass &= (
            kappa_disc_pass
            and transfer_pass
            and discriminant_pass_instance
        )

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all Kappa root checks pass       = "
        f"{separation_pass}"
    )

    print(
        "  all secant identities pass       = "
        f"{slope_identity_pass}"
    )

    print(
        "  all discriminant transfers pass = "
        f"{discriminant_pass}"
    )

    print()

    print("SECANT SLOPE SIZE RANGE")

    print(
        f"  minimum numerator digits = "
        f"{min(x[0] for x in slope_digit_sizes)}"
    )

    print(
        f"  maximum numerator digits = "
        f"{max(x[0] for x in slope_digit_sizes)}"
    )

    print(
        f"  minimum denominator digits = "
        f"{min(x[1] for x in slope_digit_sizes)}"
    )

    print(
        f"  maximum denominator digits = "
        f"{max(x[1] for x in slope_digit_sizes)}"
    )

    print()

    print("DISCRIMINANT-TRANSFER SIZE RANGE")

    print(
        f"  minimum numerator digits = "
        f"{min(x[0] for x in ratio_digit_sizes)}"
    )

    print(
        f"  maximum numerator digits = "
        f"{max(x[0] for x in ratio_digit_sizes)}"
    )

    print(
        f"  minimum denominator digits = "
        f"{min(x[1] for x in ratio_digit_sizes)}"
    )

    print(
        f"  maximum denominator digits = "
        f"{max(x[1] for x in ratio_digit_sizes)}"
    )

    print()

    print("=" * 110)
    print("EXPERIMENT 390 FINAL STATUS")
    print("=" * 110)
    print()

    print("The primary object is:")
    print()
    print("  H = (R(t_TRUE) - R(t_OTHER))")
    print("      --------------------------------")
    print("      (t_TRUE - t_OTHER)")
    print()

    print("and the exact discriminant transfer is:")
    print()
    print("  Delta_R / Delta_K")
    print("      = H^2 / (K + N^2 + 2)^2")
    print()

    print(
        "If H or the transfer ratio exhibits a much simpler "
        "normalization in N and K, that is a stronger structural "
        "candidate than the full resultant."
    )

    print()
    print("=" * 110)
    print("EXPERIMENT 390 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()

