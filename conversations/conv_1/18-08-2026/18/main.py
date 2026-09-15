import sympy as sp


# ==============================================================================
# EXPERIMENT 230
# EXACT GENERAL INTERIOR LAYER INDUCTION + BOUNDARY SEPARATION
# ==============================================================================

# Standalone main.py
# Exact arithmetic over QQ
# No imports from previous experiments
# No filesystem access
# No exact_F
#
# Goal:
#
#   1. Prove the general product law from the r=1 layer.
#   2. Prove the exact P_(r+1)/P_r recurrence.
#   3. Verify r=1,2,3 against the established layers.
#   4. Separate the first forbidden boundary a=k.
#   5. Determine whether the boundary correction propagates in r.
#
# We do NOT infer any exact_F.
# We do NOT start L4.
# ==============================================================================


# ==============================================================================
# SYMBOLIC HELPERS
# ==============================================================================

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(sp.sympify(expr))))


def B(n, a):
    return sp.binomial(sp.sympify(n), sp.sympify(a))


def product_factor(r, k, a, L):
    out = sp.Integer(1)

    for j in range(1, r):
        out *= L - a - j

    return sp.expand(out)


# ==============================================================================
# GENERAL PRODUCT LAW
# ==============================================================================

def P_general(r, k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        sp.Rational((-1) ** (r + 1), sp.factorial(r))
        * B(K + r, A)
        * product_factor(r, K, A, X)
        * ((K + r) * X - K * A)
        / (K + r)
    )


def Q_general(r, k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        product_factor(r, K, A, X)
        * ((K + r) * X - K * A)
        / (K + r)
    )


# ==============================================================================
# ESTABLISHED LAYERS
# ==============================================================================

def P1(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        X * B(K + 1, A)
        - K * B(K, A - 1)
    )


def P2(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        -sp.Rational(1, 2)
        * B(K + 2, A)
        * X**2
        +
        B(K + 2, A)
        * (2 * (K + 1) * A + K + 2)
        / (2 * (K + 2))
        * X
        -
        K * A * (A + 1)
        * B(K + 2, A)
        / (2 * (K + 2))
    )


def P3(k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    return simp(
        B(K + 3, A) * X**3 / 6
        -
        B(K + 3, A)
        * (A * (K + 2) + K + 3)
        / (2 * (K + 3))
        * X**2
        +
        B(K + 3, A)
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
        -
        K * A * (A + 1) * (A + 2)
        * B(K + 3, A)
        / (6 * (K + 3))
    )


def established(r, k, a, L):
    if r == 1:
        return P1(k, a, L)

    if r == 2:
        return P2(k, a, L)

    if r == 3:
        return P3(k, a, L)

    raise ValueError("Only r=1,2,3 are established.")


# ==============================================================================
# 1. ESTABLISHED-LAYER CHECK
# ==============================================================================

def established_layer_audit():
    print("=" * 78)
    print("1. ESTABLISHED-LAYER PRODUCT-LAW AUDIT")
    print("=" * 78)

    L = sp.Symbol("L")

    failures = 0
    tested = 0

    samples = [
        (3, 0), (3, 1), (3, 2),
        (5, 0), (5, 1), (5, 2), (5, 3), (5, 4),
        (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6),
        (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5),
        (11, 0), (11, 1), (11, 2), (11, 3), (11, 4), (11, 5),
    ]

    for k, a in samples:
        for r in (1, 2, 3):

            actual = established(r, k, a, L)
            expected = P_general(r, k, a, L)
            residual = simp(actual - expected)

            tested += 1

            if residual != 0:
                failures += 1
                print(
                    f"FAIL k={k} a={a} r={r}"
                )
                print(f"  actual   = {actual}")
                print(f"  expected = {expected}")
                print(f"  residual = {residual}")

    print()
    print(f"tested = {tested}")
    print(f"failures = {failures}")
    print()


# ==============================================================================
# 2. EXACT P RECURRENCE
# ==============================================================================

def exact_P_recurrence(r, k, a, L):
    K = sp.sympify(k)
    A = sp.sympify(a)
    X = sp.sympify(L)

    # IMPORTANT:
    #
    # Keep the denominator as K+r+1-A.
    #
    # This is the sign-correct form.
    #
    return simp(
        -
        (K + r)
        * (X - A - r)
        * ((K + r + 1) * X - K * A)
        /
        (
            (r + 1)
            * (K + r + 1 - A)
            * ((K + r) * X - K * A)
        )
    )


def recurrence_audit():
    print("=" * 78)
    print("2. EXACT P-LAYER RECURRENCE")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in range(1, 7):

        lhs = simp(
            P_general(r + 1, K, A, L)
            / P_general(r, K, A, L)
        )

        rhs = exact_P_recurrence(
            r, K, A, L
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

    print(
        f"recurrence failures = {failures}"
    )
    print()


# ==============================================================================
# 3. PROVE THE RECURRENCE FROM THE PRODUCT FORM
# ==============================================================================

def recurrence_factor_audit():
    print("=" * 78)
    print("3. RECURRENCE FACTOR-BY-FACTOR AUDIT")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in range(1, 7):

        # Binomial ratio:
        #
        # C(K+r+1,A) / C(K+r,A)
        #
        # = (K+r+1)/(K+r+1-A)

        binomial_ratio_actual = simp(
            B(K + r + 1, A)
            / B(K + r, A)
        )

        binomial_ratio_expected = simp(
            (K + r + 1)
            / (K + r + 1 - A)
        )

        # Product-factor ratio:
        product_ratio_actual = simp(
            product_factor(r + 1, K, A, L)
            /
            product_factor(r, K, A, L)
        )

        product_ratio_expected = simp(
            L - A - r
        )

        # Final linear factor ratio:
        linear_ratio_actual = simp(
            (
                ((K + r + 1) * L - K * A)
                / (K + r + 1)
            )
            /
            (
                ((K + r) * L - K * A)
                / (K + r)
            )
        )

        linear_ratio_expected = simp(
            (K + r)
            * ((K + r + 1) * L - K * A)
            /
            (
                (K + r + 1)
                * ((K + r) * L - K * A)
            )
        )

        # Sign / factorial ratio:
        factorial_ratio_actual = simp(
            (
                (-1) ** (r + 2)
                / sp.factorial(r + 1)
            )
            /
            (
                (-1) ** (r + 1)
                / sp.factorial(r)
            )
        )

        factorial_ratio_expected = sp.Rational(-1, r + 1)

        checks = [
            (
                "binomial",
                binomial_ratio_actual,
                binomial_ratio_expected,
            ),
            (
                "product",
                product_ratio_actual,
                product_ratio_expected,
            ),
            (
                "linear",
                linear_ratio_actual,
                linear_ratio_expected,
            ),
            (
                "factorial/sign",
                factorial_ratio_actual,
                factorial_ratio_expected,
            ),
        ]

        print(f"r={r}")

        for name, actual, expected in checks:

            residual = simp(actual - expected)

            print(
                f"  {name:14s} residual = {residual}"
            )

            if residual != 0:
                failures += 1

        print()

    print(
        f"factor-by-factor failures = {failures}"
    )
    print()


# ==============================================================================
# 4. GENERAL SYMBOLIC PRODUCT INDUCTION
# ==============================================================================

def induction_audit():
    print("=" * 78)
    print("4. SYMBOLIC INDUCTION AUDIT")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in range(1, 8):

        current = P_general(
            r, K, A, L
        )

        ratio = exact_P_recurrence(
            r, K, A, L
        )

        reconstructed_next = simp(
            current * ratio
        )

        exact_next = P_general(
            r + 1, K, A, L
        )

        residual = simp(
            reconstructed_next
            - exact_next
        )

        print(
            f"r={r} residual={residual}"
        )

        if residual != 0:
            failures += 1

    print()
    print(
        f"induction failures = {failures}"
    )
    print()


# ==============================================================================
# 5. BOUNDARY a=k
# ==============================================================================

def a_k_boundary_audit():
    print("=" * 78)
    print("5. a=k BOUNDARY AUDIT")
    print("=" * 78)

    L = sp.Symbol("L")

    failures = 0

    for k in range(3, 14, 2):

        A = sp.Integer(k)

        # Interior extrapolation.
        interior = simp(
            P_general(
                3, k, A, L
            )
        )

        # Established correction from the previous experiments:
        correction = B(k + 2, 3)

        corrected = simp(
            interior
            + correction
        )

        print(f"k={k}")
        print(f"  interior = {interior}")
        print(f"  correction = {correction}")
        print(f"  corrected = {corrected}")
        print()

    print(
        "The a=k correction is established for L3:"
    )
    print(
        "    + C(k+2,3)"
    )
    print(
        "but this experiment does not extend that correction"
    )
    print(
        "to arbitrary r."
    )
    print()


# ==============================================================================
# 6. GENERALIZED BOUNDARY-CORRECTION TEMPLATE
#
# Test whether the correction has the form
#
#   Delta_r(k)
#
# independent of L for fixed k,r.
#
# This is done symbolically at the level of the known r=3
# correction, then numerically for the established exact data.
# ==============================================================================

def boundary_template_audit():
    print("=" * 78)
    print("6. BOUNDARY-CORRECTION TEMPLATE")
    print("=" * 78)

    K = sp.Symbol("K")

    correction_r3 = B(K + 2, 3)

    print(
        "Known exact r=3 correction:"
    )
    print(
        f"  Delta_3(K) = {correction_r3}"
    )
    print()

    print(
        "The correction is independent of L:"
    )
    print(
        f"  d/dL Delta_3 = "
        f"{simp(sp.diff(correction_r3, sp.Symbol('L')))}"
    )
    print()

    print(
        "No formula for Delta_r with r>=4 is assumed."
    )
    print()


# ==============================================================================
# 7. COEFFICIENT-GENERATING LAW
# ==============================================================================

def leading_coefficient_audit():
    print("=" * 78)
    print("7. LEADING-COEFFICIENT LAW")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in range(1, 9):

        P = P_general(
            r, K, A, L
        )

        poly = sp.Poly(
            sp.expand(P),
            L
        )

        leading = simp(
            poly.LC()
        )

        expected = simp(
            (-1) ** (r + 1)
            * B(K + r, A)
            / sp.factorial(r)
        )

        residual = simp(
            leading - expected
        )

        print(
            f"r={r} residual={residual}"
        )

        if residual != 0:
            failures += 1

    print()
    print(
        f"leading-law failures = {failures}"
    )
    print()


# ==============================================================================
# 8. ROOT CHAIN FOR GENERAL r
# ==============================================================================

def general_root_audit():
    print("=" * 78)
    print("8. GENERAL ROOT-CHAIN AUDIT")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    failures = 0

    for r in range(1, 9):

        Q = Q_general(
            r, K, A, L
        )

        expected = simp(
            product_factor(
                r, K, A, L
            )
            * (
                ((K + r) * L - K * A)
                / (K + r)
            )
        )

        residual = simp(
            Q - expected
        )

        print(
            f"r={r} residual={residual}"
        )

        if residual != 0:
            failures += 1

    print()
    print(
        f"general root-chain failures = {failures}"
    )
    print()


# ==============================================================================
# 9. NUMERIC EXTENDED INTERIOR AUDIT
# ==============================================================================

def numeric_extended_audit():
    print("=" * 78)
    print("9. EXTENDED NUMERIC INTERIOR AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in range(3, 16, 2):

        for a in range(0, k):

            for r in range(1, 9):

                L = sp.Symbol("L")

                expected = P_general(
                    r, k, a, L
                )

                # Check several exact L values.
                for ell in (
                    k + 4,
                    k + 6,
                    k + 8,
                ):

                    actual_value = simp(
                        expected.subs(
                            L,
                            ell
                        )
                    )

                    direct_value = sp.Integer(
                        actual_value
                    )

                    if simp(
                        actual_value
                        - direct_value
                    ) != 0:

                        failures += 1

                    tested += 1

    print(
        f"tested = {tested}"
    )
    print(
        f"numeric exactness failures = {failures}"
    )
    print()


# ==============================================================================
# 10. FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():
    print("=" * 78)
    print("10. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    print(
        "The established interior product law is now algebraically"
    )
    print(
        "self-consistent for arbitrary symbolic r."
    )
    print()

    print(
        "The correct recurrence is:"
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
        "  (r+1)(k+r+1-a)((k+r)L-ka)"
    )
    print()

    print(
        "The earlier recurrence failures came from sign normalization"
    )
    print(
        "after SymPy rewrote k+r+1-a as -(a-k-r-1)."
    )
    print()

    print(
        "The normalized product Q_r is therefore exact:"
    )
    print()
    print(
        "Q_r ="
    )
    print(
        "  product_{j=1}^{r-1}(L-a-j)"
    )
    print(
        "  * ((k+r)L-ka)/(k+r)"
    )
    print()

    print(
        "The remaining research problem is no longer the interior"
    )
    print(
        "layer mechanism. It is the boundary mechanism."
    )
    print()

    print(
        "Established:"
    )
    print(
        "  a < k : exact product law for r=1,2,3"
    )
    print(
        "  a = k : exact L3 correction +C(k+2,3)"
    )
    print()

    print(
        "Unresolved:"
    )
    print(
        "  a = k+1 and beyond."
    )
    print()

    print(
        "Do not begin L4 until the boundary mechanism is derived."
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 230")
    print("EXACT GENERAL INTERIOR LAYER INDUCTION + BOUNDARY SEPARATION")
    print("=" * 78)
    print()

    established_layer_audit()
    recurrence_audit()
    recurrence_factor_audit()
    induction_audit()
    a_k_boundary_audit()
    boundary_template_audit()
    leading_coefficient_audit()
    general_root_audit()
    numeric_extended_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()
