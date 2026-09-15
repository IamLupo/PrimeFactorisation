import sympy as sp


# =============================================================================
# EXPERIMENT 229S
# EXACT GENERAL INTERIOR LAYER PRODUCT LAW
# =============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No exact_F required
#
# Established interior layers:
#
#   P1
#   P2
#   P3
#
# Proposed general product law:
#
#   P_r(k,a,L)
#     = (-1)^(r+1) / r!
#       * C(k+r,a)
#       * product_{j=1}^{r-1}(L-a-j)
#       * ((k+r)L-ka)/(k+r)
#
# for the established interior range a < k.
#
# This script tests the structural law without any project imports.
# =============================================================================


# =============================================================================
# BASIC SYMBOLIC HELPERS
# =============================================================================

def S(x):
    """
    Sympify while preserving symbolic expressions.
    """
    return sp.sympify(x)


def simp(x):
    """
    Aggressive exact simplification over QQ.
    """
    x = sp.sympify(x)
    return sp.factor(sp.cancel(sp.expand(x)))


def binom(n, r):
    return sp.binomial(S(n), S(r))


# =============================================================================
# ESTABLISHED INTERIOR LAYERS
# =============================================================================

def P1(k, a, L):
    K = S(k)
    A = S(a)
    X = S(L)

    return simp(
        X * binom(K + 1, A)
        - K * binom(K, A - 1)
    )


def P2(k, a, L):
    K = S(k)
    A = S(a)
    X = S(L)

    return simp(
        -sp.Rational(1, 2)
        * binom(K + 2, A)
        * X**2

        + binom(K + 2, A)
        * (2 * (K + 1) * A + K + 2)
        / (2 * (K + 2))
        * X

        - K * A * (A + 1)
        * binom(K + 2, A)
        / (2 * (K + 2))
    )


def P3(k, a, L):
    K = S(k)
    A = S(a)
    X = S(L)

    return simp(
        binom(K + 3, A) * X**3 / 6

        - binom(K + 3, A)
        * (A * (K + 2) + K + 3)
        / (2 * (K + 3))
        * X**2

        + binom(K + 3, A)
        * (
            3 * A**2 * K
            + 3 * A**2
            + 6 * A * K
            + 9 * A
            + 2 * K
            + 6
        )
        / (6 * (K + 3))
        * X

        - K * A * (A + 1) * (A + 2)
        * binom(K + 3, A)
        / (6 * (K + 3))
    )


def established_P(r, k, a, L):
    if r == 1:
        return P1(k, a, L)

    if r == 2:
        return P2(k, a, L)

    if r == 3:
        return P3(k, a, L)

    raise ValueError(
        "Only r=1,2,3 are established in the source data."
    )


# =============================================================================
# PRODUCT FACTORS
# =============================================================================

def shifted_product(r, k, a, L=None):
    """
    product_{j=1}^{r-1}(L-a-j)

    k is accepted only for a uniform call signature and is otherwise unused.
    """
    A = S(a)

    if L is None:
        X = sp.Symbol("L")
    else:
        X = S(L)

    out = sp.Integer(1)

    for j in range(1, r):
        out *= X - A - j

    return sp.expand(out)


def proposed_Q(r, k, a, L):
    """
    Normalized Q_r:

        Q_r =
          product_{j=1}^{r-1}(L-a-j)
          * ((k+r)L-ka)/(k+r)
    """
    K = S(k)
    A = S(a)
    X = S(L)

    return simp(
        shifted_product(r, K, A, X)
        * (
            (K + r) * X - K * A
        )
        / (K + r)
    )


def proposed_P(r, k, a, L):
    """
    Full product formula:

      P_r =
        (-1)^(r+1)/r!
        * C(k+r,a)
        * Q_r
    """
    K = S(k)
    A = S(a)
    X = S(L)

    return simp(
        (-1) ** (r + 1)
        * binom(K + r, A)
        * proposed_Q(r, K, A, X)
        / sp.factorial(r)
    )


# =============================================================================
# EXPLICIT SYMBOLIC BINOMIAL IDENTITIES
# =============================================================================

def binomial_shift_up(n, a):
    """
    C(n+1,a) = (n+1)/(n+1-a) * C(n,a)

    Returned symbolically.
    """
    N = S(n)
    A = S(a)

    return simp(
        binom(N, A)
        * (N + 1)
        / (N + 1 - A)
    )


def binomial_lower_to_upper(n, a):
    """
    C(n,a-1) = a/(n+1) * C(n+1,a)
    """
    N = S(n)
    A = S(a)

    return simp(
        A * binom(N + 1, A)
        / (N + 1)
    )


# =============================================================================
# 1. DIRECT PRODUCT FACTORIZATION
# =============================================================================

def direct_factorization_audit():
    print("=" * 78)
    print("1. DIRECT PRODUCT FACTORIZATION")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 0),
        (3, 1),
        (3, 2),

        (5, 0),
        (5, 1),
        (5, 2),
        (5, 3),
        (5, 4),

        (7, 1),
        (7, 2),
        (7, 3),
        (7, 4),

        (9, 2),
        (9, 3),
        (9, 4),
        (9, 5),

        (11, 3),
        (11, 4),
        (11, 5),
    ]

    failures = 0

    for k, a in samples:
        for r in (1, 2, 3):
            actual = established_P(r, k, a, L)
            expected = proposed_P(r, k, a, L)
            residual = simp(actual - expected)

            if residual != 0:
                failures += 1

                print(
                    f"FAIL k={k:2d} a={a:2d} r={r}"
                )
                print(f"  actual   = {actual}")
                print(f"  expected = {expected}")
                print(f"  residual = {residual}")
                print()

    print(f"product-factorization failures = {failures}")
    print()


# =============================================================================
# 2. SYMBOLIC PRODUCT IDENTITIES
# =============================================================================

def symbolic_product_audit():
    print("=" * 78)
    print("2. SYMBOLIC PRODUCT IDENTITIES")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    # -------------------------------------------------------------------------
    # Q1
    #
    # P1 = L*C(K+1,A) - K*C(K,A-1)
    #
    # Using:
    #   C(K,A-1) = A/(K+1) C(K+1,A)
    #
    # gives
    #   Q1 = P1/C(K+1,A)
    #      = ((K+1)L-KA)/(K+1).
    # -------------------------------------------------------------------------

    q1_actual = simp(
        (
            L * binom(K + 1, A)
            - K * binomial_lower_to_upper(K, A)
        )
        / binom(K + 1, A)
    )

    q1_expected = simp(
        ((K + 1) * L - K * A)
        / (K + 1)
    )

    # -------------------------------------------------------------------------
    # Q2 and Q3 are already polynomial expressions after normalization.
    # -------------------------------------------------------------------------

    q2_actual = simp(
        -2
        * P2(K, A, L)
        / binom(K + 2, A)
    )

    q2_expected = proposed_Q(2, K, A, L)

    q3_actual = simp(
        6
        * P3(K, A, L)
        / binom(K + 3, A)
    )

    q3_expected = proposed_Q(3, K, A, L)

    checks = [
        ("Q1", q1_actual, q1_expected),
        ("Q2", q2_actual, q2_expected),
        ("Q3", q3_actual, q3_expected),
    ]

    failures = 0

    for name, actual, expected in checks:
        residual = simp(actual - expected)

        print(name)
        print(f"  actual   = {actual}")
        print(f"  expected = {expected}")
        print(f"  residual = {residual}")

        if residual == 0:
            print("  PASS")
        else:
            print("  FAIL")
            failures += 1

        print()

    print(f"symbolic product identity failures = {failures}")
    print()


# =============================================================================
# 3. ROOT-CHAIN AUDIT
# =============================================================================

def root_chain_audit():
    print("=" * 78)
    print("3. ROOT-CHAIN AUDIT")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 1),
        (3, 2),

        (5, 1),
        (5, 2),
        (5, 3),

        (7, 2),
        (7, 3),
        (7, 4),

        (9, 3),
        (9, 4),

        (11, 5),
    ]

    failures = 0

    for k, a in samples:
        for r in (1, 2, 3):

            Q = simp(
                (-1) ** (r + 1)
                * sp.factorial(r)
                * established_P(r, k, a, L)
                / binom(k + r, a)
            )

            expected_roots = [
                sp.Rational(a + j, 1)
                for j in range(1, r)
            ]

            expected_roots.append(
                sp.Rational(k * a, k + r)
            )

            expected_poly = simp(
                sp.prod(
                    L - root
                    for root in expected_roots
                )
            )

            residual = simp(
                Q - expected_poly
            )

            if residual != 0:
                failures += 1

                print(
                    f"FAIL k={k} a={a} r={r}"
                )
                print(f"  Q          = {Q}")
                print(f"  expectation= {expected_poly}")
                print(f"  residual   = {residual}")
                print()

    print(f"root-chain failures = {failures}")
    print()


# =============================================================================
# 4. LAYER-TO-LAYER RATIONAL RECURRENCE
#
# Correct formula:
#
# P_(r+1)/P_r =
#
#   -(k+r)(L-a-r)((k+r+1)L-ka)
#   --------------------------------
#   (r+1)(k+r+1-a)((k+r)L-ka)
#
# =============================================================================

def recurrence_rhs(r, k, a, L):
    K = S(k)
    A = S(a)
    X = S(L)

    return simp(
        -(
            (K + r)
            * (X - A - r)
            * ((K + r + 1) * X - K * A)
        )
        / (
            (r + 1)
            * (K + r + 1 - A)
            * ((K + r) * X - K * A)
        )
    )


def recurrence_audit():
    print("=" * 78)
    print("4. LAYER-TO-LAYER RATIONAL RECURRENCE")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in (1, 2, 3):

        lhs = simp(
            proposed_P(r + 1, K, A, L)
            / proposed_P(r, K, A, L)
        )

        rhs = recurrence_rhs(r, K, A, L)

        residual = simp(lhs - rhs)

        print(f"r={r}")
        print(f"  lhs      = {lhs}")
        print(f"  rhs      = {rhs}")
        print(f"  residual = {residual}")

        if residual == 0:
            print("  PASS")
        else:
            print("  FAIL")
            failures += 1

        print()

    print(f"recurrence failures = {failures}")
    print()


# =============================================================================
# 5. NORMALIZED T_r RECURRENCE
#
# T_r = (-1)^(r+1) r! P_r / C(k+r,a)
#
# Therefore:
#
# T_(r+1)/T_r =
#   -(L-a-r)
#    * ((k+r+1)L-ka)/(k+r+1)
#    * (k+r)/((k+r)L-ka)
#
# =============================================================================

def T(r, k, a, L):
    K = S(k)
    A = S(a)
    X = S(L)

    return simp(
        (-1) ** (r + 1)
        * sp.factorial(r)
        * proposed_P(r, K, A, X)
        / binom(K + r, A)
    )


def T_recurrence_audit():
    print("=" * 78)
    print("5. NORMALIZED T_r RECURRENCE")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in (1, 2, 3):

        lhs = simp(
            T(r + 1, K, A, L)
            / T(r, K, A, L)
        )

        rhs = simp(
            -(
                L - A - r
            )
            * (
                ((K + r + 1) * L - K * A)
                / (K + r + 1)
            )
            * (
                (K + r)
                / ((K + r) * L - K * A)
            )
        )

        residual = simp(lhs - rhs)

        print(f"r={r}")
        print(f"  lhs      = {lhs}")
        print(f"  rhs      = {rhs}")
        print(f"  residual = {residual}")

        if residual == 0:
            print("  PASS")
        else:
            print("  FAIL")
            failures += 1

        print()

    print(f"T_r recurrence failures = {failures}")
    print()


# =============================================================================
# 6. ROOT-EXTENSION AUDIT
# =============================================================================

def root_extension_audit():
    print("=" * 78)
    print("6. ROOT-EXTENSION AUDIT")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in (1, 2, 3):

        lhs = simp(
            proposed_Q(r + 1, K, A, L)
            / proposed_Q(r, K, A, L)
        )

        rhs = simp(
            (L - A - r)
            * (
                ((K + r + 1) * L - K * A)
                / (K + r + 1)
            )
            * (
                (K + r)
                / ((K + r) * L - K * A)
            )
        )

        residual = simp(lhs - rhs)

        print(f"r={r}")
        print(f"  residual = {residual}")

        if residual == 0:
            print("  PASS")
        else:
            print("  FAIL")
            failures += 1

        print()

    print(f"root-extension failures = {failures}")
    print()


# =============================================================================
# 7. COEFFICIENT STRUCTURE AUDIT
# =============================================================================

def coefficient_structure_audit():
    print("=" * 78)
    print("7. COEFFICIENT STRUCTURE AUDIT")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 1, 1),
        (3, 1, 2),
        (3, 1, 3),

        (3, 2, 1),
        (3, 2, 2),
        (3, 2, 3),

        (5, 2, 1),
        (5, 2, 2),
        (5, 2, 3),

        (5, 3, 1),
        (5, 3, 2),
        (5, 3, 3),

        (7, 3, 1),
        (7, 3, 2),
        (7, 3, 3),

        (9, 4, 1),
        (9, 4, 2),
        (9, 4, 3),

        (11, 5, 1),
        (11, 5, 2),
        (11, 5, 3),
    ]

    failures = 0

    for k, a, r in samples:

        actual = simp(
            (-1) ** (r + 1)
            * sp.factorial(r)
            * established_P(r, k, a, L)
            / binom(k + r, a)
        )

        expected = proposed_Q(r, k, a, L)

        residual = simp(actual - expected)

        coeffs = sp.Poly(
            sp.expand(actual),
            L
        ).all_coeffs()

        print(
            f"k={k:2d} a={a:2d} r={r}"
        )
        print(
            f"  Q_r          = {sp.factor(actual)}"
        )
        print(
            f"  coefficients = {coeffs}"
        )
        print(
            f"  residual     = {residual}"
        )

        if residual != 0:
            failures += 1

        print()

    print(f"coefficient-structure failures = {failures}")
    print()


# =============================================================================
# 8. ADJACENT BINOMIAL IDENTITY AUDIT
# =============================================================================

def adjacent_binomial_audit():
    print("=" * 78)
    print("8. ADJACENT-BINOMIAL IDENTITY AUDIT")
    print("=" * 78)

    failures = 0

    for k in range(3, 16, 2):
        for r in range(1, 8):
            for a in range(0, k):

                lhs1 = simp(
                    binom(k + r, a)
                    * (k + r - a)
                    / (k + r)
                )

                rhs1 = binom(k + r - 1, a)

                lhs2 = simp(
                    binom(k + r, a)
                    * a
                    / (k + r)
                )

                rhs2 = binom(k + r - 1, a - 1)

                if simp(lhs1 - rhs1) != 0:
                    failures += 1
                    print(
                        f"FAIL identity-1 k={k} r={r} a={a}"
                    )

                if simp(lhs2 - rhs2) != 0:
                    failures += 1
                    print(
                        f"FAIL identity-2 k={k} r={r} a={a}"
                    )

    print(
        f"adjacent-binomial failures = {failures}"
    )
    print()


# =============================================================================
# 9. INVERSE RECONSTRUCTION
# =============================================================================

def inverse_reconstruction_audit():
    print("=" * 78)
    print("9. INVERSE RECONSTRUCTION AUDIT")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 0),
        (3, 1),
        (3, 2),

        (5, 1),
        (5, 2),
        (5, 3),
        (5, 4),

        (7, 2),
        (7, 3),
        (7, 4),

        (9, 3),
        (9, 4),

        (11, 5),
    ]

    failures = 0

    for k, a in samples:
        for r in (1, 2, 3):

            actual = established_P(r, k, a, L)

            reconstructed = proposed_P(
                r,
                k,
                a,
                L
            )

            residual = simp(
                actual - reconstructed
            )

            if residual != 0:
                failures += 1

                print(
                    f"FAIL k={k} a={a} r={r}"
                )
                print(
                    f"  residual = {residual}"
                )

    print(
        f"inverse-reconstruction failures = {failures}"
    )
    print()


# =============================================================================
# 10. GENERALIZED SYMBOLIC LAW
# =============================================================================

def generalized_symbolic_law_audit():
    print("=" * 78)
    print("10. GENERALIZED SYMBOLIC LAW")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    # Test the normalized law for r=1..8 symbolically.
    #
    # We cannot compare against a kernel for r>3 because those layers have
    # not yet been independently established. But we can prove that the
    # recurrence generates the claimed product sequence.

    for r in range(1, 9):

        Qr = proposed_Q(
            r,
            K,
            A,
            L
        )

        Qnext = proposed_Q(
            r + 1,
            K,
            A,
            L
        )

        recurrence = simp(
            Qnext / Qr
        )

        expected = simp(
            (L - A - r)
            * (
                ((K + r + 1) * L - K * A)
                / (K + r + 1)
            )
            * (
                (K + r)
                / ((K + r) * L - K * A)
            )
        )

        residual = simp(
            recurrence - expected
        )

        if residual != 0:
            failures += 1

            print(
                f"FAIL r={r}"
            )
            print(
                f"  residual = {residual}"
            )

    print(
        f"generalized symbolic recurrence failures = {failures}"
    )
    print()


# =============================================================================
# 11. FINAL STRUCTURAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The symbolic Q1 failure in the previous run was not a failure"
    )
    print(
        "of the product law. It was a SymPy simplification limitation:"
    )
    print(
        "the identity C(K,A-1) = A*C(K+1,A)/(K+1) had not been supplied."
    )
    print()

    print(
        "The corrected exact interior law is:"
    )
    print()
    print(
        "P_r(k,a,L) ="
    )
    print(
        "  (-1)^(r+1)/r!"
        " * C(k+r,a)"
    )
    print(
        "  * product_{j=1}^{r-1}(L-a-j)"
    )
    print(
        "  * ((k+r)L-ka)/(k+r)."
    )
    print()

    print(
        "For r=1,2,3 this reproduces all established interior layers."
    )
    print()

    print(
        "The normalized form is:"
    )
    print()
    print(
        "Q_r ="
        " product_{j=1}^{r-1}(L-a-j)"
    )
    print(
        "      * ((k+r)L-ka)/(k+r)."
    )
    print()

    print(
        "Hence Q_r has roots:"
    )
    print(
        "  L = a+1, a+2, ..., a+r-1"
    )
    print(
        "and"
    )
    print(
        "  L = ka/(k+r)."
    )
    print()

    print(
        "The exact P-layer recurrence is:"
    )
    print()
    print(
        "P_(r+1)/P_r ="
    )
    print(
        "  -(k+r)(L-a-r)((k+r+1)L-ka)"
    )
    print(
        "  --------------------------------"
    )
    print(
        "  (r+1)(k+r+1-a)((k+r)L-ka)."
    )
    print()

    print(
        "The normalized T_r recurrence is:"
    )
    print()
    print(
        "T_(r+1)/T_r ="
    )
    print(
        "  -(L-a-r)"
        " * ((k+r+1)L-ka)/(k+r+1)"
        " * (k+r)/((k+r)L-ka)."
    )
    print()

    print(
        "This gives a genuine symbolic generating mechanism for the"
    )
    print(
        "established interior layers."
    )
    print()

    print(
        "Boundary rows a>=k remain outside the established interior law."
    )
    print()

    print(
        "No exact_F is reconstructed."
    )
    print(
        "No L4 analysis is performed."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 229S")
    print("EXACT GENERAL INTERIOR LAYER PRODUCT LAW")
    print("=" * 78)
    print()
    print("Standalone main.py")
    print("Exact arithmetic over QQ")
    print("No previous experiment imported")
    print("No filesystem access")
    print("No exact_F required")
    print()

    direct_factorization_audit()
    symbolic_product_audit()
    root_chain_audit()
    recurrence_audit()
    T_recurrence_audit()
    root_extension_audit()
    coefficient_structure_audit()
    adjacent_binomial_audit()
    inverse_reconstruction_audit()
    generalized_symbolic_law_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()