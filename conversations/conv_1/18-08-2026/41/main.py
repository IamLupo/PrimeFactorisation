import sympy as sp

# ==============================================================================
# EXPERIMENT 251
# r=5 INTERIOR-LAW / BOUNDARY DELTA ISOLATION AUDIT
# ==============================================================================

# Exact arithmetic only.
# No previous experiment imported.
# No filesystem access.
# No r=6.

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")
K, A, L, D = sp.symbols("K A L D")


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


# ------------------------------------------------------------------------------
# Exact kernel
# ------------------------------------------------------------------------------

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ------------------------------------------------------------------------------
# Exact symmetric reduction
# ------------------------------------------------------------------------------

def power_sum(m):
    if m == 0:
        return sp.Integer(2)
    if m == 1:
        return N

    s0 = sp.Integer(2)
    s1 = N

    for _ in range(2, m + 1):
        s2 = sp.expand(N*s1 - X*s0)
        s0, s1 = s1, s2

    return s1


def pq_to_NX(expr):
    poly = sp.Poly(sp.expand(expr), p, q)

    out = sp.Integer(0)

    for (i, j), coeff in poly.terms():

        if i < j:
            continue

        if i == j:
            out += coeff * X**i
        else:
            out += coeff * X**j * power_sum(i - j)

    return sp.expand(out)


def coefficient(G, a, ell, r):
    b = ell - r - a

    if b < 0:
        return sp.Integer(0)

    return sp.Poly(
        sp.expand(G),
        N,
        X
    ).coeff_monomial(
        N**a * X**b
    )


# ------------------------------------------------------------------------------
# GENERAL interior product law
#
# P_r =
#   (-1)^(r+1) / r!
#   * C(k+r,a)
#   * product_{m=1}^{r-1}(L-a-m)
#   * ((k+r)L-ka)/(k+r)
# ------------------------------------------------------------------------------

def interior_general(r, k, a, ell):

    product_term = sp.Integer(1)

    for m in range(1, r):
        product_term *= ell - a - m

    return simp(
        sp.Rational((-1)**(r + 1), sp.factorial(r))
        * sp.binomial(k + r, a)
        * product_term
        * ((k + r)*ell - k*a)
        / (k + r)
    )


# ------------------------------------------------------------------------------
# Universal boundary law
# ------------------------------------------------------------------------------

def delta_universal(r, j, k, d):

    if j == 0:
        return sp.binomial(k + r - 1, r)

    product_term = sp.Integer(1)

    for m in range(1, j):
        product_term *= d - r - m

    moving_factor = (
        d - sp.Rational((r - j) * k, k + j)
    )

    return simp(
        sp.Rational((-1)**j, sp.factorial(j))
        * sp.binomial(k + r - 1, r - j)
        * product_term
        * moving_factor
    )


# ------------------------------------------------------------------------------
# 1. SYMBOLIC r=5 INTERIOR LAW
# ------------------------------------------------------------------------------

def symbolic_interior():

    print("=" * 78)
    print("1. SYMBOLIC r=5 INTERIOR LAW")
    print("=" * 78)

    actual = interior_general(5, K, A, L)

    expected = simp(
        sp.Rational(1, 120)
        * sp.binomial(K + 5, A)
        * (L - A - 1)
        * (L - A - 2)
        * (L - A - 3)
        * (L - A - 4)
        * ((K + 5)*L - K*A)
        / (K + 5)
    )

    residual = simp(actual - expected)

    print("actual   =", actual)
    print("expected =", expected)
    print("residual =", residual)
    print("status   =", "PASS" if residual == 0 else "FAIL")
    print()


# ------------------------------------------------------------------------------
# 2. INDIVIDUAL NUMERICAL INTERIOR CHECKS
# ------------------------------------------------------------------------------

def numerical_interior():

    print("=" * 78)
    print("2. NUMERICAL r=5 INTERIOR CHECK")
    print("=" * 78)

    tests = [
        (3, 9, 0),
        (3, 9, 1),
        (3, 9, 2),
        (3, 9, 3),
        (5, 11, 3),
        (5, 11, 4),
        (5, 11, 5),
    ]

    failures = 0

    for k, ell, a in tests:

        F = exact_F(k, ell)
        G = pq_to_NX(F)

        raw = coefficient(G, a, ell, 5)

        interior = interior_general(
            5, k, a, ell
        )

        delta = simp(
            raw - interior
        )

        print(
            f"k={k:2d} ell={ell:2d} a={a:2d}"
        )
        print(
            f"  raw      = {raw}"
        )
        print(
            f"  interior = {interior}"
        )
        print(
            f"  delta    = {delta}"
        )

        # Independent reconstruction of the displayed
        # interior value using the literal factors.
        literal = simp(
            sp.Rational(1, 120)
            * sp.binomial(k + 5, a)
            * (
                (ell - a - 1)
                * (ell - a - 2)
                * (ell - a - 3)
                * (ell - a - 4)
            )
            * (
                ((k + 5)*ell - k*a)
                / (k + 5)
            )
        )

        literal_residual = simp(
            interior - literal
        )

        print(
            f"  literal  = {literal}"
        )
        print(
            f"  interior-vs-literal = {literal_residual}"
        )

        ok = literal_residual == 0

        print(
            f"  {'PASS' if ok else 'FAIL'}"
        )
        print()

        if not ok:
            failures += 1

    print(
        f"interior arithmetic failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 3. FULL DELTA DECOMPOSITION
# ------------------------------------------------------------------------------

def delta_decomposition():

    print("=" * 78)
    print("3. RAW / INTERIOR / DELTA / PREDICTED DELTA")
    print("=" * 78)

    tests = [
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 11),
        (5, 13),
        (5, 15),
        (7, 13),
        (7, 15),
    ]

    failures = 0

    for k, ell in tests:

        d = ell - k

        F = exact_F(k, ell)
        G = pq_to_NX(F)

        print(
            f"k={k:2d} ell={ell:2d} d={d:2d}"
        )

        for j in range(6):

            a = k + j

            raw = coefficient(
                G, a, ell, 5
            )

            interior = interior_general(
                5, k, a, ell
            )

            delta = simp(
                raw - interior
            )

            predicted = delta_universal(
                5, j, k, d
            )

            residual = simp(
                delta - predicted
            )

            ok = residual == 0

            print(
                f"  j={j}"
                f" raw={raw}"
                f" interior={interior}"
                f" delta={delta}"
                f" predicted={predicted}"
                f" residual={residual}"
                f" {'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

        print()

    print(
        f"boundary-law failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 4. DIRECT LAYER RECONSTRUCTION
# ------------------------------------------------------------------------------

def reconstruction():

    print("=" * 78)
    print("4. DIRECT r=5 LAYER RECONSTRUCTION")
    print("=" * 78)

    tests = [
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 13),
    ]

    failures = 0

    for k, ell in tests:

        F = exact_F(k, ell)
        G = pq_to_NX(F)

        reconstructed = sp.Integer(0)

        for a in range(0, ell - 4):

            b = ell - 5 - a

            if b < 0:
                continue

            raw = coefficient(
                G, a, ell, 5
            )

            reconstructed += (
                raw
                * N**a
                * X**b
            )

        # Extract only the r=5 homogeneous layer
        # from the transformed kernel.
        direct_layer = sp.Integer(0)

        Gpoly = sp.Poly(
            sp.expand(G),
            N,
            X
        )

        for (a, b), coeff in Gpoly.terms():

            if a + 2*b == ell:
                pass

            # For this experiment the intended layer is
            # N^a X^(ell-5-a), so:
            if b == ell - 5 - a if b >= 0 else False:
                direct_layer += (
                    coeff
                    * N**a
                    * X**b
                )

        residual = sp.expand(
            reconstructed - direct_layer
        )

        ok = residual == 0

        print(
            f"k={k:2d} ell={ell:2d} "
            f"residual={residual} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"reconstruction failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 5. TERMINAL CHECK
# ------------------------------------------------------------------------------

def terminal_check():

    print("=" * 78)
    print("5. TERMINAL r=5 CHECK")
    print("=" * 78)

    tests = [
        (3, 13),
        (5, 15),
        (7, 17),
        (9, 19),
    ]

    failures = 0

    for k, ell in tests:

        d = ell - k

        F = exact_F(k, ell)
        G = pq_to_NX(F)

        j = 5
        a = k + j

        raw = coefficient(
            G, a, ell, 5
        )

        interior = interior_general(
            5, k, a, ell
        )

        delta = simp(
            raw - interior
        )

        predicted = delta_universal(
            5, j, k, d
        )

        residual = simp(
            delta - predicted
        )

        print(
            f"k={k:2d} ell={ell:2d} d={d:2d} "
            f"raw={raw} "
            f"interior={interior} "
            f"delta={delta} "
            f"predicted={predicted} "
            f"residual={residual}"
        )

        if residual != 0:
            failures += 1

    print()
    print(
        f"terminal failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    symbolic_interior()

    numerical_interior()

    delta_decomposition()

    reconstruction()

    terminal_check()

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "This run isolates the r=5 discrepancy."
    )
    print(
        "Do not perform a large r=5 sweep until:"
    )
    print(
        "  1. symbolic r=5 interior law passes;"
    )
    print(
        "  2. raw/interior/delta arithmetic is internally consistent;"
    )
    print(
        "  3. the direct r=5 layer reconstruction passes."
    )


if __name__ == "__main__":
    main()
