import sympy as sp

# ==============================================================================
# EXPERIMENT 251
# DIRECT EXACT-KERNEL r=5 UNIVERSAL BOUNDARY VALIDATION
# ==============================================================================

# Exact arithmetic over QQ.
# Standalone main.py.
# No previous experiment imported.
# No filesystem access.
# r=6 forbidden.
#
# Goal:
#   Test the universal boundary law for r=5 directly against
#   coefficients extracted from the exact pq kernel.
#
# Universal candidate:
#
#   Delta_(r,0)(k,d) = C(k+r-1, r)
#
#   Delta_(r,j)(k,d)
#     = (-1)^j / j!
#       * C(k+r-1, r-j)
#       * product_{m=1}^{j-1}(d-r-m)
#       * (d - (r-j)k/(k+j))
#
# where d = ell-k.
#
# This experiment does NOT fit polynomials.
# It does NOT infer the formula from the data.
# It computes the exact kernel coefficient first, subtracts
# the exact interior product-law extrapolation, then compares
# directly with the universal boundary candidate.
#
# IMPORTANT:
#   We use coeff extraction by Poly(...).coeff_monomial(N**a*X**b)
#   only when b >= 0.
#   This avoids the earlier 1/X PolynomialError.
# ==============================================================================


p, q = sp.symbols("p q")
N, X = sp.symbols("N X")


# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def pq_to_NX(expr):
    """
    Substitute:
        N = p + q
        X = p*q

    via exact symmetric reconstruction.

    The kernel is symmetric in p,q, so every homogeneous
    component can be rewritten in N,X.
    """
    expr = sp.expand(expr)

    poly = sp.Poly(expr, p, q)

    out = 0

    for (i, j), coeff in poly.terms():

        if i < j:
            continue

        # Symmetric combination:
        #
        # p^i q^j + p^j q^i
        #
        # For i=j:
        #   X^i
        #
        # Otherwise:
        #   X^j * (p^(i-j) + q^(i-j))
        #
        # and p^m + q^m is generated recursively from
        # N=p+q and X=pq.
        m = i - j

        if m == 0:
            sym_power = X**j
            out += coeff * sym_power
            continue

        # Power sum S_m = p^m + q^m.
        S0 = sp.Integer(2)
        S1 = N

        if m == 1:
            Sm = S1
        else:
            prev2 = S0
            prev1 = S1

            for n in range(2, m + 1):
                curr = sp.expand(N * prev1 - X * prev2)
                prev2, prev1 = prev1, curr

            Sm = prev1

        out += coeff * X**j * Sm

    return sp.expand(out)


def exact_G(k, ell):
    return sp.expand(
        pq_to_NX(
            exact_F(k, ell)
        )
    )


def coefficient(G, a, ell, r):
    """
    Extract coefficient of:
        N^a X^(ell-r-a)

    Return 0 when the X exponent is negative.
    """
    b = ell - r - a

    if b < 0:
        return sp.Integer(0)

    poly = sp.Poly(sp.expand(G), N, X)

    return sp.expand(
        poly.coeff_monomial(
            N**a * X**b
        )
    )


# ------------------------------------------------------------------------------
# Interior product law
# ------------------------------------------------------------------------------

def interior_P(r, k, a, ell):
    """
    Exact established interior product law:

      P_r(k,a,L)
        = (-1)^(r+1)/r!
          * C(k+r,a)
          * product_{j=1}^{r-1}(L-a-j)
          * ((k+r)L-ka)/(k+r)
    """
    k = sp.sympify(k)
    a = sp.sympify(a)
    L = sp.sympify(ell)

    prod = sp.Integer(1)

    for j in range(1, r):
        prod *= L - a - j

    return simp(
        sp.Rational((-1) ** (r + 1), sp.factorial(r))
        * sp.binomial(k + r, a)
        * prod
        * ((k + r) * L - k * a)
        / (k + r)
    )


# ------------------------------------------------------------------------------
# Universal boundary law
# ------------------------------------------------------------------------------

def universal_delta(r, j, k, d):
    k = sp.sympify(k)
    d = sp.sympify(d)

    if j == 0:
        return sp.binomial(k + r - 1, r)

    prod = sp.Integer(1)

    for m in range(1, j):
        prod *= d - r - m

    moving = (
        d
        - sp.Rational(r - j, 1) * k / (k + j)
    )

    return simp(
        sp.Rational((-1) ** j, sp.factorial(j))
        * sp.binomial(k + r - 1, r - j)
        * prod
        * moving
    )


# ------------------------------------------------------------------------------
# Direct exact-kernel boundary extraction
# ------------------------------------------------------------------------------

def exact_boundary_delta(r, k, ell, j):
    """
    Boundary index:
        a = k+j
    d = ell-k

    Delta = exact coefficient - interior extrapolation.
    """
    a = k + j
    G = exact_G(k, ell)

    actual = coefficient(
        G,
        a,
        ell,
        r
    )

    interior = interior_P(
        r,
        k,
        a,
        ell
    )

    return simp(actual - interior)


# ------------------------------------------------------------------------------
# 1. KERNEL SANITY
# ------------------------------------------------------------------------------

def kernel_sanity():
    print("=" * 78)
    print("1. EXACT KERNEL SANITY")
    print("=" * 78)

    failures = 0

    tests = [
        (3, 9),
        (5, 11),
        (7, 13),
        (9, 15),
        (11, 17),
        (13, 19),
    ]

    for k, ell in tests:
        F = exact_F(k, ell)

        symmetric = sp.expand(
            F - F.xreplace({p: q, q: p})
        ) == 0

        diagonal_zero = sp.expand(
            F.subs(q, p)
        ) == 0

        ok = symmetric and diagonal_zero

        print(
            f"k={k:2d} ell={ell:2d} "
            f"symmetric={symmetric} "
            f"diagonal={diagonal_zero} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"kernel sanity failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 2. r=5 DIRECT BOUNDARY TEST
# ------------------------------------------------------------------------------

def direct_r5_test():
    print("=" * 78)
    print("2. DIRECT EXACT-KERNEL r=5 BOUNDARY TEST")
    print("=" * 78)

    failures = 0
    tested = 0

    # These are intentionally small.
    #
    # d >= 6 is the stable domain.
    #
    # Multiple k values and multiple d values are used,
    # but the total number of exact expansions remains modest.
    test_points = [
        (3, 9),
        (3, 11),
        (3, 13),

        (5, 11),
        (5, 13),
        (5, 15),

        (7, 13),
        (7, 15),
        (7, 17),

        (9, 15),
        (9, 17),

        (11, 17),
        (11, 19),

        (13, 19),
        (13, 21),
    ]

    for k, ell in test_points:

        d = ell - k

        print(f"k={k:2d} ell={ell:2d} d={d:2d}")

        for j in range(0, 6):

            actual_delta = exact_boundary_delta(
                5,
                k,
                ell,
                j
            )

            expected_delta = universal_delta(
                5,
                j,
                k,
                d
            )

            residual = simp(
                actual_delta - expected_delta
            )

            ok = residual == 0

            print(
                f"  j={j} "
                f"actual={str(actual_delta):>12} "
                f"expected={str(expected_delta):>12} "
                f"residual={str(residual):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            tested += 1

            if not ok:
                failures += 1

        print()

    print(f"tested = {tested}")
    print(f"r=5 boundary failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 3. r=5 DEGREE LADDER FROM EXACT KERNEL
# ------------------------------------------------------------------------------

def r5_degree_audit():
    print("=" * 78)
    print("3. r=5 EXACT-KERNEL DEGREE LADDER")
    print("=" * 78)

    failures = 0

    # Use several k values with at least 7 stable d values.
    #
    # ell = k+d, d=6,8,...,18
    #
    # This provides enough points to identify degree <=5 without
    # asking the kernel to expand unnecessarily large cases.

    ks = [3, 5, 7, 9, 11, 13]
    ds = [6, 8, 10, 12, 14, 16, 18]

    for k in ks:
        print(f"k={k}")

        data = {}

        for j in range(0, 6):

            values = []

            for d in ds:
                ell = k + d

                delta = exact_boundary_delta(
                    5,
                    k,
                    ell,
                    j
                )

                values.append(
                    (d, delta)
                )

            data[j] = values

            # Interpolate only as an audit of degree.
            poly = sp.interpolate(
                values,
                D
            )

            poly = sp.Poly(
                sp.expand(poly),
                D
            )

            degree = poly.degree()

            ok = degree == j

            print(
                f"  j={j} "
                f"degree={degree} "
                f"expected={j} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

        print()

    print(
        f"degree-ladder failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 4. EXACT FACTORIZATION COMPARISON
# ------------------------------------------------------------------------------

def r5_factorization_audit():
    print("=" * 78)
    print("4. r=5 FACTORIZATION AUDIT")
    print("=" * 78)

    failures = 0

    ks = [3, 5, 7, 11, 13]
    ds = [6, 8, 10, 12, 14, 16]

    for k in ks:

        # Collect exact-kernel values and interpolate only after
        # direct Delta extraction.
        print(f"k={k}")

        for j in range(1, 6):

            values = []

            for d in ds:
                ell = k + d

                delta = exact_boundary_delta(
                    5,
                    k,
                    ell,
                    j
                )

                values.append(
                    (d, delta)
                )

            poly = sp.factor(
                sp.interpolate(values, D)
            )

            candidate = sp.factor(
                universal_delta(5, j, k, D)
            )

            residual = simp(
                poly - candidate
            )

            ok = residual == 0

            print(
                f"  j={j}"
            )
            print(
                f"    extracted = {poly}"
            )
            print(
                f"    candidate = {candidate}"
            )
            print(
                f"    residual  = {residual}"
            )
            print(
                f"    {'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

        print()

    print(
        f"factorization failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 5. MOVING ROOT AUDIT FROM DIRECT DATA
# ------------------------------------------------------------------------------

def moving_root_audit():
    print("=" * 78)
    print("5. r=5 MOVING-ROOT AUDIT")
    print("=" * 78)

    failures = 0

    # For each j, use the candidate moving root.
    # Then verify that the interpolated exact-kernel polynomial
    # has the same root.

    ks = [3, 5, 7, 11, 13]

    for k in ks:
        print(f"k={k}")

        for j in range(1, 6):

            ds = [
                sp.Integer(6),
                sp.Integer(8),
                sp.Integer(10),
                sp.Integer(12),
                sp.Integer(14),
                sp.Integer(16),
            ]

            values = []

            for d in ds:
                ell = k + int(d)

                delta = exact_boundary_delta(
                    5,
                    k,
                    ell,
                    j
                )

                values.append(
                    (d, delta)
                )

            poly = sp.interpolate(
                values,
                D
            )

            poly = sp.factor(poly)

            moving_root = simp(
                sp.Rational(5 - j, k + j)
                * k
            )

            residual = simp(
                poly.subs(D, moving_root)
            )

            ok = residual == 0

            print(
                f"  j={j} "
                f"moving_root={moving_root} "
                f"residual={residual} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

        print()

    print(
        f"moving-root failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 6. FIXED ROOT AUDIT
# ------------------------------------------------------------------------------

def fixed_root_audit():
    print("=" * 78)
    print("6. r=5 FIXED-ROOT AUDIT")
    print("=" * 78)

    failures = 0

    ks = [3, 5, 7, 11, 13]

    for k in ks:
        print(f"k={k}")

        for j in range(2, 6):

            ds = [
                sp.Integer(6),
                sp.Integer(8),
                sp.Integer(10),
                sp.Integer(12),
                sp.Integer(14),
                sp.Integer(16),
            ]

            values = []

            for d in ds:
                ell = k + int(d)

                delta = exact_boundary_delta(
                    5,
                    k,
                    ell,
                    j
                )

                values.append(
                    (d, delta)
                )

            poly = sp.factor(
                sp.interpolate(values, D)
            )

            print(f"  j={j}")

            for root in range(6, 6 + j - 1):

                residual = simp(
                    poly.subs(D, root)
                )

                ok = residual == 0

                print(
                    f"    D={root} "
                    f"residual={residual} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

                if not ok:
                    failures += 1

        print()

    print(
        f"fixed-root failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 7. FINAL SUMMARY
# ------------------------------------------------------------------------------

def final_summary():
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The experiment tests r=5 directly against the exact pq kernel."
    )
    print()

    print(
        "Universal candidate:"
    )
    print(
        "  Delta_(r,0)(k,d) = C(k+r-1,r)"
    )
    print()
    print(
        "  Delta_(r,j)(k,d) = "
        "(-1)^j/j! * C(k+r-1,r-j)"
    )
    print(
        "                      * product_{m=1}^{j-1}(d-r-m)"
    )
    print(
        "                      * (d-(r-j)k/(k+j))"
    )
    print()

    print(
        "For r=5:"
    )
    print()

    for j in range(6):
        print(
            f"  j={j}: "
            f"{sp.factor(universal_delta(5, j, K, D))}"
        )

    print()

    print(
        "No r=6 analysis is performed."
    )


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

def main():
    kernel_sanity()
    direct_r5_test()
    r5_degree_audit()
    r5_factorization_audit()
    moving_root_audit()
    fixed_root_audit()
    final_summary()


if __name__ == "__main__":
    main()
