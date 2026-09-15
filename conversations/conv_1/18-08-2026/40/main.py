import sympy as sp

# ==============================================================================
# EXPERIMENT 250
# r=5 KERNEL EXTRACTION PIPELINE AUDIT
# ==============================================================================

# Exact arithmetic only.
# No previous experiment imported.
# No filesystem access.
# No r=6.

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")
D = sp.symbols("D")


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


# ------------------------------------------------------------------------------
# Exact pq kernel
# ------------------------------------------------------------------------------

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ------------------------------------------------------------------------------
# Method A: custom symmetric reduction
# ------------------------------------------------------------------------------

def power_sum(m):
    """
    p^m + q^m in terms of N=p+q, X=pq.
    """
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


def custom_pq_to_NX(expr):
    poly = sp.Poly(sp.expand(expr), p, q)

    out = sp.Integer(0)

    for (i, j), coeff in poly.terms():

        # Keep one member of each symmetric pair.
        if i < j:
            continue

        if i == j:
            out += coeff * X**i
            continue

        m = i - j
        out += coeff * X**j * power_sum(m)

    return sp.expand(out)


# ------------------------------------------------------------------------------
# Method B: SymPy symmetrize
# ------------------------------------------------------------------------------

def sympy_pq_to_NX(expr):
    """
    SymPy's independent elementary-symmetric reduction.

    For two variables:
        s1 = p+q = N
        s2 = pq  = X
    """
    symmetric_part, remainder, mapping = sp.symmetrize(
        sp.expand(expr),
        [p, q],
        formal=True
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Unexpected nonsymmetric remainder: {remainder}"
        )

    # mapping is normally:
    # [(s1, p + q), (s2, p*q)]
    s1, s2 = mapping[0][0], mapping[1][0]

    return sp.expand(
        symmetric_part.subs({
            s1: N,
            s2: X,
        })
    )


# ------------------------------------------------------------------------------
# Independent identity check
# ------------------------------------------------------------------------------

def direct_substitution_check(expr_nx, expr_pq):
    """
    Verify the N,X expression by substituting:

        N -> p+q
        X -> p*q
    """
    return sp.expand(
        expr_nx.subs({
            N: p + q,
            X: p*q,
        }) - expr_pq
    )


# ------------------------------------------------------------------------------
# Coefficient extraction
# ------------------------------------------------------------------------------

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
# r=5 interior product law
# ------------------------------------------------------------------------------

def interior_P5(k, a, ell):
    L = sp.sympify(ell)

    prod = sp.Integer(1)

    for j in range(1, 5):
        prod *= L - a - j

    return simp(
        sp.Rational(1, 120)
        * sp.binomial(k + 5, a)
        * prod
        * ((k + 5)*L - k*a)
        / (k + 5)
    )


# ------------------------------------------------------------------------------
# Universal r=5 boundary law
# ------------------------------------------------------------------------------

def delta5(j, k, d):
    if j == 0:
        return sp.binomial(k + 4, 5)

    prod = sp.Integer(1)

    for m in range(1, j):
        prod *= d - 5 - m

    moving = d - sp.Rational(5 - j, k + j)*k

    return simp(
        sp.Rational((-1)**j, sp.factorial(j))
        * sp.binomial(k + 4, 5 - j)
        * prod
        * moving
    )


# ------------------------------------------------------------------------------
# 1. Kernel symmetry audit
# ------------------------------------------------------------------------------

def kernel_audit():
    print("=" * 78)
    print("1. KERNEL SYMMETRY AUDIT")
    print("=" * 78)

    failures = 0

    tests = [
        (3, 9),
        (5, 11),
        (7, 13),
        (9, 15),
    ]

    for k, ell in tests:
        F = exact_F(k, ell)

        symmetric = sp.expand(
            F - F.xreplace({p: q, q: p})
        ) == 0

        # Correct diagonal statement:
        diagonal = sp.expand(F.subs(q, p))

        # Only symmetry is required.
        ok = symmetric

        print(
            f"k={k:2d} ell={ell:2d} "
            f"symmetric={symmetric} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        print(
            f"  F(p,p) = {sp.factor(diagonal)}"
        )

        if not ok:
            failures += 1

    print()
    print(f"kernel symmetry failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 2. CUSTOM vs SYMPY symmetric reduction
# ------------------------------------------------------------------------------

def conversion_audit():
    print("=" * 78)
    print("2. CUSTOM N,X CONVERSION AUDIT")
    print("=" * 78)

    failures = 0

    tests = [
        (1, 3),
        (2, 5),
        (3, 7),
        (3, 9),
        (5, 11),
        (7, 13),
    ]

    for k, ell in tests:
        F = exact_F(k, ell)

        A = custom_pq_to_NX(F)
        B = sympy_pq_to_NX(F)

        residual_custom = direct_substitution_check(A, F)
        residual_sympy = direct_substitution_check(B, F)
        residual_cross = sp.expand(A - B)

        ok = (
            residual_custom == 0
            and residual_sympy == 0
            and residual_cross == 0
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"custom={residual_custom == 0} "
            f"sympy={residual_sympy == 0} "
            f"cross={residual_cross == 0} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            print("  custom residual =", residual_custom)
            print("  sympy residual  =", residual_sympy)
            print("  cross residual  =", residual_cross)
            failures += 1

    print()
    print(f"conversion failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 3. r=2,r=3,r=4 coefficient regression
# ------------------------------------------------------------------------------

def known_layer_regression():
    print("=" * 78)
    print("3. KNOWN r=2/r=3/r=4 EXTRACTION REGRESSION")
    print("=" * 78)

    failures = 0

    # These are coefficients already established by the previous
    # experiments. The purpose is to make sure the extraction
    # machinery still reproduces them before attempting r=5.

    tests = [
        (2, 3, 9, 3),
        (2, 5, 11, 5),
        (3, 3, 9, 3),
        (3, 5, 11, 5),
        (4, 3, 9, 3),
        (4, 5, 11, 5),
        (4, 7, 13, 7),
    ]

    for r, k, ell, a in tests:

        F = exact_F(k, ell)
        G1 = custom_pq_to_NX(F)
        G2 = sympy_pq_to_NX(F)

        c1 = coefficient(G1, a, ell, r)
        c2 = coefficient(G2, a, ell, r)

        residual = simp(c1 - c2)

        ok = residual == 0

        print(
            f"r={r} k={k} ell={ell} a={a} "
            f"custom={c1} "
            f"sympy={c2} "
            f"residual={residual} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"known-layer extraction failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 4. r=5 LOCAL BOUNDARY CHECK
# ------------------------------------------------------------------------------

def r5_local_check():
    print("=" * 78)
    print("4. r=5 LOCAL BOUNDARY CHECK")
    print("=" * 78)

    failures_custom = 0
    failures_sympy = 0

    tests = [
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 13),
        (7, 13),
        (7, 15),
    ]

    for k, ell in tests:

        d = ell - k

        F = exact_F(k, ell)

        G_custom = custom_pq_to_NX(F)
        G_sympy = sympy_pq_to_NX(F)

        print(
            f"k={k:2d} ell={ell:2d} d={d:2d}"
        )

        for j in range(6):

            a = k + j

            actual_custom = coefficient(
                G_custom,
                a,
                ell,
                5
            )

            actual_sympy = coefficient(
                G_sympy,
                a,
                ell,
                5
            )

            interior = interior_P5(
                k,
                a,
                ell
            )

            delta_custom = simp(
                actual_custom - interior
            )

            delta_sympy = simp(
                actual_sympy - interior
            )

            expected = delta5(
                j,
                k,
                d
            )

            rc = simp(
                delta_custom - expected
            )

            rs = simp(
                delta_sympy - expected
            )

            okc = rc == 0
            oks = rs == 0

            print(
                f"  j={j} "
                f"custom={delta_custom} "
                f"sympy={delta_sympy} "
                f"expected={expected} "
                f"rc={rc} "
                f"rs={rs}"
            )

            print(
                f"       custom={'PASS' if okc else 'FAIL'} "
                f"sympy={'PASS' if oks else 'FAIL'}"
            )

            if not okc:
                failures_custom += 1

            if not oks:
                failures_sympy += 1

        print()

    print(
        f"custom r=5 failures = {failures_custom}"
    )
    print(
        f"SymPy r=5 failures  = {failures_sympy}"
    )
    print()


# ------------------------------------------------------------------------------
# 5. r=5 exact reconstruction comparison
# ------------------------------------------------------------------------------

def reconstruction_audit():
    print("=" * 78)
    print("5. r=5 FULL LOCAL RECONSTRUCTION")
    print("=" * 78)

    failures_custom = 0
    failures_sympy = 0

    tests = [
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 13),
        (7, 13),
    ]

    for k, ell in tests:

        F = exact_F(k, ell)

        Gc = custom_pq_to_NX(F)
        Gs = sympy_pq_to_NX(F)

        rc = sp.Integer(0)
        rs = sp.Integer(0)

        # Only reconstruct the r=5 layer.
        for a in range(0, 5 + ell):

            b = ell - 5 - a

            if b < 0:
                continue

            cc = coefficient(Gc, a, ell, 5)
            cs = coefficient(Gs, a, ell, 5)

            rc += (
                cc
                * N**a
                * X**b
            )

            rs += (
                cs
                * N**a
                * X**b
            )

        original_component_custom = sp.expand(rc)
        original_component_sympy = sp.expand(rs)

        residual = sp.expand(
            original_component_custom
            - original_component_sympy
        )

        ok = residual == 0

        print(
            f"k={k:2d} ell={ell:2d} "
            f"custom-vs-sympy residual={residual} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures_custom += 1

    print()
    print(
        f"reconstruction comparison failures = "
        f"{failures_custom}"
    )
    print()


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    kernel_audit()

    conversion_audit()

    known_layer_regression()

    r5_local_check()

    reconstruction_audit()

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "Do not run a large r=5 sweep until the custom and SymPy"
    )
    print(
        "symmetric reductions agree exactly and the known r=2/r=3/r=4"
    )
    print(
        "coefficient regressions pass."
    )


if __name__ == "__main__":
    main()

