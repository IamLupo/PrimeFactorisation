import sympy as sp


# ============================================================================
# EXPERIMENT 227
# KERNEL RECONSTRUCTION FROM ESTABLISHED HOMOGENEOUS-LAYER DATA
# ============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No filesystem access
# No previous experiment imported
#
# Goal:
#   Recover the structural form of the missing exact_F(k, ell)
#   from identities that have already been established.
#
# We do NOT fit theorem-anchor values.
# We do NOT search arbitrary affine transformations.
# We do NOT assume a candidate exact_F.
#
# Instead we test the following consequences:
#
#   N = p*q
#   X = p+q+1
#   N+X = (p+1)(q+1)
#
#   G_top = -X^(ell-k) * ((N+X)^k - N^k)
#
#   L1 coefficient:
#       [N^a X^(ell-1-a)] L1
#         = ell*C(k+1,a) - k*C(k,a-1)
#
#   L2 interior coefficient:
#       A = -C(k+2,a)/2
#       B = C(k+2,a)*(2(k+1)a+k+2)/(2(k+2))
#       C = -k*a*(a+1)*C(k+2,a)/(2*(k+2))
#
#   L3 interior coefficient:
#       A = C(k+3,a)/6
#       B = -(a*k+2*a+k+3)*C(k+3,a)/(2*(k+3))
#       C = (3*a^2*k+3*a^2+6*a*k+9*a+2*k+6)
#           *C(k+3,a)/(6*(k+3))
#       D = -a*k*(a+1)*(a+2)*C(k+3,a)/(6*(k+3))
#
# The experiment then asks:
#
#   1. What universal coefficient-generating expression reproduces
#      the top layer?
#   2. What multiplicative / binomial correction generates L1?
#   3. Whether the L0,L1,L2,L3 coefficient arrays are consistent with
#      a single finite product / logarithmic derivative construction.
#   4. Whether simple candidates built only from
#         X, N, N+X, ell, k
#      can reproduce all four layers.
#
# The point is to recover the algebraic architecture of F rather than
# fabricate an exact_F from finitely many numerical values.
# ============================================================================


p, q = sp.symbols("p q")
N, X = sp.symbols("N X")
k, ell = sp.symbols("k ell", integer=True, positive=True)


# ----------------------------------------------------------------------------
# Basic helpers
# ----------------------------------------------------------------------------

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def C(n, r):
    return sp.binomial(n, r)


def coeff(expr, a, b):
    """
    Coefficient of N^a X^b in a polynomial.
    """
    poly = sp.Poly(sp.expand(expr), N, X)
    return sp.expand(poly.coeff_monomial(N**a * X**b))


def homogeneous(expr, total_degree):
    """
    Extract homogeneous component of total degree = total_degree
    in N,X.
    """
    poly = sp.Poly(sp.expand(expr), N, X)
    out = sp.Integer(0)

    for (a, b), c in poly.terms():
        if a + b == total_degree:
            out += c * N**a * X**b

    return sp.expand(out)


def substitute_x(expr, x_value):
    return sp.expand(expr.subs(X, x_value))


# ----------------------------------------------------------------------------
# Established layers
# ----------------------------------------------------------------------------

def L0_exact(kv, Lv):
    """
    Established highest homogeneous layer.
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    return sp.expand(
        -X**(L - K) * ((N + X)**K - N**K)
    )


def L1_coefficient(kv, Lv, a):
    K = sp.Integer(kv)
    L = sp.Integer(Lv)
    aa = sp.Integer(a)

    return sp.simplify(
        L * C(K + 1, aa)
        - K * C(K, aa - 1)
    )


def L1_exact(kv, Lv):
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    degree = L - 1
    out = sp.Integer(0)

    for a in range(0, K + 1):
        b = degree - a
        if b < 0:
            continue

        c = L1_coefficient(kv, Lv, a)
        if c != 0:
            out += c * N**a * X**b

    return sp.expand(out)


def L2_interior_coefficient(kv, Lv, a):
    K = sp.Integer(kv)
    L = sp.Integer(Lv)
    aa = sp.Integer(a)

    A = -C(K + 2, aa) / sp.Integer(2)

    B = (
        C(K + 2, aa)
        * (2 * (K + 1) * aa + K + 2)
        / (2 * (K + 2))
    )

    C0 = (
        -K * aa * (aa + 1) * C(K + 2, aa)
        / (2 * (K + 2))
    )

    return sp.expand(
        A * L**2 + B * L + C0
    )


def L2_boundary_k(kv, Lv):
    """
    a=k boundary:
      interior extrapolation + C(k+2,3)
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    return sp.expand(
        L2_interior_coefficient(kv, Lv, K)
        + C(K + 2, 3)
    )


def L2_boundary_k1(kv, Lv):
    """
    a=k+1.
    Established relation from the exact boundary data:
        D_(k+1) = 2 D_k / (k+1)
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    return sp.simplify(
        sp.Rational(2, K + 1) * L2_boundary_k(kv, Lv)
    )


def L2_exact(kv, Lv):
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    degree = L - 2
    out = sp.Integer(0)

    for a in range(0, K + 2):
        b = degree - a
        if b < 0:
            continue

        if a < K:
            c = L2_interior_coefficient(kv, Lv, a)
        elif a == K:
            c = L2_boundary_k(kv, Lv)
        else:
            c = L2_boundary_k1(kv, Lv)

        if c != 0:
            out += c * N**a * X**b

    return sp.expand(out)


def L3_interior_coefficient(kv, Lv, a):
    K = sp.Integer(kv)
    L = sp.Integer(Lv)
    aa = sp.Integer(a)

    A = C(K + 3, aa) / sp.Integer(6)

    B = (
        -(aa * K + 2 * aa + K + 3)
        * C(K + 3, aa)
        / (2 * (K + 3))
    )

    C0 = (
        (
            3 * aa**2 * K
            + 3 * aa**2
            + 6 * aa * K
            + 9 * aa
            + 2 * K
            + 6
        )
        * C(K + 3, aa)
        / (6 * (K + 3))
    )

    D = (
        -aa * K * (aa + 1) * (aa + 2)
        * C(K + 3, aa)
        / (6 * (K + 3))
    )

    return sp.expand(
        A * L**3
        + B * L**2
        + C0 * L
        + D
    )


def L3_boundary_k(kv, Lv):
    """
    Established a=k correction:
        interior + C(k+2,3)
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    return sp.expand(
        L3_interior_coefficient(kv, Lv, K)
        + C(K + 2, 3)
    )


# ----------------------------------------------------------------------------
# Structural candidate family
# ----------------------------------------------------------------------------

def candidate_top(kv, Lv):
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    return sp.expand(
        -X**(L - K) * ((N + X)**K - N**K)
    )


def candidate_first_order_operator(kv, Lv):
    """
    The established L1 coefficients strongly suggest differentiating
    the basic degree-k object while also differentiating X^(ell-k).

    This gives the structural expression:

      X^(ell-k-1) *
      [ ell (N+X)^k - k N (N+X)^(k-1) ]

    We test this identity directly.
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    return sp.expand(
        X**(L - K - 1)
        * (
            L * (N + X)**K
            - K * N * (N + X)**(K - 1)
        )
    )


def candidate_second_order_family(kv, Lv):
    """
    Build the most natural second-order candidate from repeated
    differentiation of X^(ell-k)*(N+X)^k.

    This is deliberately a structural test, not a fitted formula.
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    base = X**(L - K) * (N + X)**K

    # Euler-type operator acting on X.
    E = X * sp.diff(base, X)

    E2 = X * sp.diff(E, X)

    return sp.expand(
        E2
        - K * (X + N) * sp.diff(base, N)
    )


def candidate_shift_factor(kv, Lv, order):
    """
    Generic exact finite-difference-style structural candidate:

      X^(ell-k-order) * (N+X)^k

    with order = 0,1,2,3.

    This helps identify whether lower layers are generated by
    repeated shifts of the same basic binomial object.
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    power = L - K - order

    if power < 0:
        return sp.Integer(0)

    return sp.expand(
        X**power * (N + X)**K
    )


# ----------------------------------------------------------------------------
# Layer reconstruction from an unknown kernel
#
# We do not have exact_F. Instead, we encode ONLY the exact structural
# consequences already established, and ask what algebraic basis they imply.
# ----------------------------------------------------------------------------

def reconstructed_known_layers(kv, Lv):
    """
    Sum of established L0+L1+L2+known L3 interior/boundary pieces.
    This is NOT exact_F.
    It is a diagnostic reconstruction through degree ell-3.
    """
    K = sp.Integer(kv)
    L = sp.Integer(Lv)

    result = L0_exact(kv, Lv)
    result += L1_exact(kv, Lv)
    result += L2_exact(kv, Lv)

    degree = L - 3

    if degree >= 0:
        for a in range(0, K + 2):
            b = degree - a

            if b < 0:
                continue

            if a < K:
                c = L3_interior_coefficient(kv, Lv, a)
            elif a == K:
                c = L3_boundary_k(kv, Lv)
            else:
                # Observed exact support has no stable accepted formula here.
                continue

            if c != 0:
                result += c * N**a * X**b

    return sp.expand(result)


# ----------------------------------------------------------------------------
# Exact structural checks
# ----------------------------------------------------------------------------

def top_layer_audit():
    print("=" * 78)
    print("1. TOP-LAYER STRUCTURAL IDENTITY")
    print("=" * 78)

    pairs = [
        (1, 3),
        (1, 10),
        (3, 7),
        (3, 20),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for kv, Lv in pairs:
        top = candidate_top(kv, Lv)
        degree = Lv

        ok = True

        for a in range(0, kv):
            b = degree - a
            actual = coeff(top, a, b)
            expected = -C(kv, a)

            if sp.simplify(actual - expected) != 0:
                ok = False
                print(
                    f"FAIL k={kv} ell={Lv} a={a} "
                    f"actual={actual} expected={expected}"
                )

        print(
            f"k={kv:2d} ell={Lv:2d} "
            f"top-layer identity = {'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"top-layer failures = {failures}")
    print()


def first_layer_structural_audit():
    print("=" * 78)
    print("2. FIRST-LAYER STRUCTURAL DERIVATION")
    print("=" * 78)

    pairs = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for kv, Lv in pairs:
        expr = candidate_first_order_operator(kv, Lv)
        expected_expr = L1_exact(kv, Lv)

        delta = sp.expand(expr - expected_expr)

        ok = delta == 0

        print(
            f"k={kv:2d} ell={Lv:2d} "
            f"L1 operator identity = {'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            print("  residual =", sp.factor(delta))
            failures += 1

    print()
    print(f"L1 operator failures = {failures}")
    print()


def second_layer_operator_audit():
    print("=" * 78)
    print("3. SECOND-LAYER STRUCTURAL TEST")
    print("=" * 78)

    pairs = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    for kv, Lv in pairs:
        candidate = candidate_second_order_family(kv, Lv)

        degree = Lv - 2
        actual_layer = L2_exact(kv, Lv)

        candidate_layer = homogeneous(candidate, degree)

        print(f"k={kv:2d} ell={Lv:2d}")
        print("  candidate degree-(ell-2) layer:")
        print("   ", sp.factor(candidate_layer))
        print("  established L2:")
        print("   ", sp.factor(actual_layer))
        print("  residual:")
        print("   ", sp.factor(sp.expand(candidate_layer - actual_layer)))
        print()


def shift_basis_audit():
    print("=" * 78)
    print("4. SHIFTED BINOMIAL BASIS AUDIT")
    print("=" * 78)

    orders = [0, 1, 2, 3]

    pairs = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
    ]

    for kv, Lv in pairs:
        print(f"k={kv:2d} ell={Lv:2d}")

        for order in orders:
            expr = candidate_shift_factor(kv, Lv, order)
            print(
                f"  order={order} "
                f"degree={sp.Poly(expr, N, X).total_degree() if expr != 0 else '-'}"
            )

        print()


def endpoint_symbolic_audit():
    print("=" * 78)
    print("5. SYMBOLIC BOUNDARY AUDIT")
    print("=" * 78)

    K, L = sp.symbols("K L", integer=True, positive=True)

    # Interior L3 evaluated at a=k.
    aK = K

    A = C(K + 3, aK) / 6
    B = (
        -(aK * K + 2 * aK + K + 3)
        * C(K + 3, aK)
        / (2 * (K + 3))
    )
    C0 = (
        (
            3 * aK**2 * K
            + 3 * aK**2
            + 6 * aK * K
            + 9 * aK
            + 2 * K
            + 6
        )
        * C(K + 3, aK)
        / (6 * (K + 3))
    )
    D = (
        -aK * K * (aK + 1) * (aK + 2)
        * C(K + 3, aK)
        / (6 * (K + 3))
    )

    Dk_int = sp.factor(
        A * L**3 + B * L**2 + C0 * L + D
    )

    correction = C(K + 2, 3)
    Dk_exact = sp.factor(Dk_int + correction)

    print("Interior extrapolation at a=k:")
    print(" ", Dk_int)

    print()
    print("Boundary correction:")
    print(" ", correction)

    print()
    print("Exact symbolic boundary candidate:")
    print(" ", Dk_exact)

    print()
    print("Check correction expansion:")
    print(
        " ",
        sp.factor(
            sp.expand(Dk_exact - Dk_int - correction)
        ),
    )

    print()


def boundary_candidate_generation():
    print("=" * 78)
    print("6. a=k+1 / a=k+2 BOUNDARY CANDIDATE GENERATION")
    print("=" * 78)

    K, L = sp.symbols("K L", integer=True, positive=True)

    # ------------------------------------------------------------------
    # Interior extrapolation at a=k+1.
    # ------------------------------------------------------------------

    a = K + 1

    A = C(K + 3, a) / 6

    B = (
        -(a * K + 2 * a + K + 3)
        * C(K + 3, a)
        / (2 * (K + 3))
    )

    C0 = (
        (
            3 * a**2 * K
            + 3 * a**2
            + 6 * a * K
            + 9 * a
            + 2 * K
            + 6
        )
        * C(K + 3, a)
        / (6 * (K + 3))
    )

    D = (
        -a * K * (a + 1) * (a + 2)
        * C(K + 3, a)
        / (6 * (K + 3))
    )

    Dkp1_int = sp.factor(
        sp.expand(A * L**3 + B * L**2 + C0 * L + D)
    )

    print("Interior extrapolation at a=k+1:")
    print(" ", Dkp1_int)

    print()
    print("Expanded:")
    print(" ", sp.expand(Dkp1_int))

    print()
    print("Important:")
    print(
        "  This is an extrapolation only."
        "\n  It is NOT asserted to equal the exact kernel boundary row."
    )

    print()

    # ------------------------------------------------------------------
    # Interior extrapolation at a=k+2.
    # ------------------------------------------------------------------

    a = K + 2

    Dkp2_int = sp.factor(
        sp.expand(
            (
                C(K + 3, a) / 6 * L**3
                - (
                    a * K + 2 * a + K + 3
                )
                * C(K + 3, a)
                / (2 * (K + 3))
                * L**2
                + (
                    (
                        3 * a**2 * K
                        + 3 * a**2
                        + 6 * a * K
                        + 9 * a
                        + 2 * K
                        + 6
                    )
                    * C(K + 3, a)
                    / (6 * (K + 3))
                ) * L
                - (
                    a * K * (a + 1) * (a + 2)
                    * C(K + 3, a)
                    / (6 * (K + 3))
                )
            )
        )
    )

    print("Interior extrapolation at a=k+2:")
    print(" ", Dkp2_int)

    print()
    print("Required cancellation for zero tail:")
    print(" ", sp.factor(-Dkp2_int))

    print()


def coefficient_generating_audit():
    print("=" * 78)
    print("7. COEFFICIENT-GENERATING PATTERN AUDIT")
    print("=" * 78)

    # Compare the leading coefficients of L0,L1,L2,L3.
    #
    # They are:
    #
    # L0: -C(k,a)
    # L1:  ell*C(k+1,a) + ...
    # L2: -ell^2*C(k+2,a)/2 + ...
    # L3: +ell^3*C(k+3,a)/6 + ...
    #
    # This is the central structural signal.

    for kv in [3, 5, 7, 9]:
        print(f"k={kv}")

        for a in range(kv):
            c0 = -C(kv, a)
            c1 = C(kv + 1, a)
            c2 = -C(kv + 2, a) / 2
            c3 = C(kv + 3, a) / 6

            print(
                f"  a={a:2d} "
                f"L0={c0} "
                f"L1_lead={c1} "
                f"L2_lead={c2} "
                f"L3_lead={c3}"
            )

        print()


def reconstruction_from_generic_ansatz():
    print("=" * 78)
    print("8. GENERIC BINOMIAL-POWER ANSATZ AUDIT")
    print("=" * 78)

    print(
        "The highest-degree-in-ell coefficients are exactly"
    )
    print(
        "(-1)^r / r! * C(k+r,a)"
    )
    print(
        "for r=0,1,2,3."
    )

    print()

    failures = 0

    for kv in [3, 5, 7, 9, 11]:
        for a in range(kv):
            observed = [
                -C(kv, a),
                C(kv + 1, a),
                -C(kv + 2, a) / 2,
                C(kv + 3, a) / 6,
            ]

            expected = [
                (-1)**0 * C(kv + 0, a) / sp.factorial(0),
                (-1)**1 * C(kv + 1, a) / sp.factorial(1),
                (-1)**2 * C(kv + 2, a) / sp.factorial(2),
                (-1)**3 * C(kv + 3, a) / sp.factorial(3),
            ]

            if any(
                sp.simplify(observed[i] - expected[i]) != 0
                for i in range(4)
            ):
                failures += 1

    print(
        "highest-degree coefficient law =",
        "PASS" if failures == 0 else "FAIL",
    )
    print("failures =", failures)
    print()


def final_diagnostic():
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The strongest recovered structural law is:"
    )
    print()
    print(
        "  coefficient of ell^r in the degree-(ell-r) layer"
    )
    print(
        "      = (-1)^r / r! * C(k+r,a)"
    )
    print()
    print(
        "for the established interior range."
    )
    print()
    print(
        "This strongly suggests that the missing exact kernel"
    )
    print(
        "contains a finite binomial / falling-factorial mechanism"
    )
    print(
        "whose r-th homogeneous correction increments k -> k+r."
    )
    print()
    print(
        "The next derivation target should therefore be:"
    )
    print(
        "  identify the exact generating mechanism responsible"
    )
    print(
        "  for the sequence"
    )
    print(
        "      C(k,a), C(k+1,a), C(k+2,a), C(k+3,a), ..."
    )
    print()
    print(
        "rather than fitting another boundary polynomial."
    )
    print()
    print(
        "In particular, investigate whether F can be represented"
    )
    print(
        "as a finite difference / binomial expansion in an object"
    )
    print(
        "whose natural bases are"
    )
    print(
        "      X = p+q+1,"
    )
    print(
        "      N = pq,"
    )
    print(
        "      N+X = (p+1)(q+1)."
    )
    print()
    print(
        "No exact_F has been fabricated in this experiment."
    )
    print()


def main():
    print("=" * 78)
    print("EXPERIMENT 227")
    print("KERNEL RECONSTRUCTION FROM ESTABLISHED LAYER STRUCTURE")
    print("=" * 78)
    print()
    print("Standalone main.py")
    print("Exact arithmetic over QQ")
    print("No previous experiment imported")
    print("No filesystem access")
    print()
    print(
        "This experiment does NOT claim to reconstruct exact_F."
    )
    print(
        "It isolates the strongest algebraic constraints on exact_F."
    )
    print()

    top_layer_audit()
    first_layer_structural_audit()
    second_layer_operator_audit()
    shift_basis_audit()
    endpoint_symbolic_audit()
    boundary_candidate_generation()
    coefficient_generating_audit()
    reconstruction_from_generic_ansatz()
    final_diagnostic()


if __name__ == "__main__":
    main()