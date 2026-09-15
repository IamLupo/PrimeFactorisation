import sympy as sp

# ==============================================================================
# EXPERIMENT 252
# EXACT r=5 ACTUAL BOUNDARY LADDER DERIVATION
# ==============================================================================

# Exact arithmetic only.
# Standalone.
# No filesystem access.
# No r=6.
# Stable domain: d = ell-k >= 6.

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")

K, A, L, D = sp.symbols("K A L D")


# ------------------------------------------------------------------------------
# Helpers
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


def power_sum(m):
    """
    p^m + q^m expressed in N=p+q, X=pq.
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


def pq_to_NX(expr):
    """
    Symmetric reduction in the basis N=p+q, X=pq.
    """
    poly = sp.Poly(sp.expand(expr), p, q)

    out = sp.Integer(0)

    for (i, j), coeff in poly.terms():

        # Use symmetry:
        # p^i q^j + p^j q^i
        # = X^min(i,j) * (p^(i-j)+q^(i-j)).
        #
        # Because F is symmetric, it is enough to retain i >= j.
        if i < j:
            continue

        if i == j:
            out += coeff * X**i
        else:
            out += coeff * X**j * power_sum(i - j)

    return sp.expand(out)


def coeff_r5(G, a, ell):
    """
    Coefficient of N^a X^(ell-5-a).
    """
    b = ell - 5 - a

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
# Interior r=5 law
# ------------------------------------------------------------------------------

def interior_r5(k, a, ell):

    return simp(
        sp.Rational(1, 120)
        * sp.binomial(k + 5, a)
        * (ell - a - 1)
        * (ell - a - 2)
        * (ell - a - 3)
        * (ell - a - 4)
        * ((k + 5)*ell - k*a)
        / (k + 5)
    )


# ------------------------------------------------------------------------------
# Actual Delta_5,j
# ------------------------------------------------------------------------------

def actual_delta5(k, ell, j):

    a = k + j

    F = exact_F(k, ell)
    G = pq_to_NX(F)

    raw = coeff_r5(G, a, ell)

    interior = interior_r5(
        k, a, ell
    )

    return simp(raw - interior)


# ------------------------------------------------------------------------------
# 1. RAW r=5 DATA
# ------------------------------------------------------------------------------

def raw_data():

    print("=" * 78)
    print("1. ACTUAL r=5 BOUNDARY DELTAS")
    print("=" * 78)

    ks = [3, 5, 7, 9, 11, 13]

    for k in ks:

        print(f"k={k}")

        for j in range(6):

            values = []

            for d in range(6, 18, 2):

                ell = k + d

                delta = actual_delta5(
                    k, ell, j
                )

                values.append(
                    (d, sp.factor(delta))
                )

            print(
                f"  j={j}: {values}"
            )

        print()


# ------------------------------------------------------------------------------
# 2. FINITE-DIFFERENCE DEGREE
# ------------------------------------------------------------------------------

def finite_differences(values):

    rows = [list(values)]

    while len(rows[-1]) > 1:
        prev = rows[-1]

        nxt = [
            sp.simplify(
                prev[i + 1] - prev[i]
            )
            for i in range(len(prev) - 1)
        ]

        rows.append(nxt)

    return rows


def detected_degree(values):

    rows = finite_differences(values)

    for order in range(len(rows)):

        if len(rows[order]) <= 1:
            return order

        if all(
            sp.simplify(x - rows[order][0]) == 0
            for x in rows[order]
        ):
            return order

    return len(rows) - 1


def degree_audit():

    print("=" * 78)
    print("2. ACTUAL FINITE-DIFFERENCE DEGREE")
    print("=" * 78)

    ks = [3, 5, 7, 9, 11, 13]

    for k in ks:

        print(f"k={k}")

        for j in range(6):

            ds = list(range(6, 18, 2))

            vals = [
                actual_delta5(
                    k,
                    k + d,
                    j
                )
                for d in ds
            ]

            degree = detected_degree(vals)

            print(
                f"  j={j} degree={degree}"
            )

        print()


# ------------------------------------------------------------------------------
# 3. INTERPOLATE ACTUAL STABLE POLYNOMIAL
# ------------------------------------------------------------------------------

def interpolate_actual():

    print("=" * 78)
    print("3. ACTUAL STABLE POLYNOMIALS")
    print("=" * 78)

    ds = list(range(6, 18, 2))
    Dsym = D

    all_polys = {}

    for k in [3, 5, 7, 9, 11, 13]:

        print(f"k={k}")

        for j in range(6):

            pts = []

            for d in ds:

                delta = actual_delta5(
                    k,
                    k + d,
                    j
                )

                pts.append(
                    (sp.Integer(d), delta)
                )

            poly = sp.interpolate(
                pts,
                Dsym
            )

            poly = sp.factor(
                sp.cancel(
                    sp.expand(poly)
                )
            )

            all_polys[(k, j)] = poly

            print(
                f"  j={j}: {poly}"
            )

        print()

    return all_polys


# ------------------------------------------------------------------------------
# 4. CROSS-k COEFFICIENT COMPARISON
# ------------------------------------------------------------------------------

def coefficient_table(all_polys):

    print("=" * 78)
    print("4. CROSS-k COEFFICIENT TABLE")
    print("=" * 78)

    for j in range(6):

        print(f"j={j}")

        for k in [3, 5, 7, 9, 11, 13]:

            poly = sp.Poly(
                all_polys[(k, j)],
                D
            )

            coeffs = [
                sp.factor(
                    poly.nth(n)
                )
                for n in range(poly.degree(), -1, -1)
            ]

            print(
                f"  k={k}: {coeffs}"
            )

        print()


# ------------------------------------------------------------------------------
# 5. ROOT AUDIT
# ------------------------------------------------------------------------------

def root_audit(all_polys):

    print("=" * 78)
    print("5. ACTUAL ROOT AUDIT")
    print("=" * 78)

    for j in range(6):

        print(f"j={j}")

        for k in [3, 5, 7, 9, 11, 13]:

            poly = all_polys[(k, j)]

            print(
                f"  k={k}: "
                f"{sp.factor(poly)}"
            )

            try:
                roots = sp.solve(
                    sp.Eq(poly, 0),
                    D
                )

                print(
                    f"       roots={roots}"
                )

            except Exception as exc:
                print(
                    f"       root extraction failed: {exc}"
                )

        print()


# ------------------------------------------------------------------------------
# 6. TEST THE OLD UNIVERSAL CANDIDATE
# ------------------------------------------------------------------------------

def proposed_universal(r, j, k, d):

    if j == 0:
        return sp.binomial(
            k + r - 1,
            r
        )

    product = sp.Integer(1)

    for m in range(1, j):
        product *= (
            d - r - m
        )

    return simp(
        sp.Rational(
            (-1)**j,
            sp.factorial(j)
        )
        * sp.binomial(
            k + r - 1,
            r - j
        )
        * product
        * (
            d
            - sp.Rational(
                (r - j)*k,
                k + j
            )
        )
    )


def compare_old_law():

    print("=" * 78)
    print("6. COMPARISON WITH PROPOSED UNIVERSAL LAW")
    print("=" * 78)

    failures = 0

    for k in [3, 5, 7, 9, 11, 13]:

        for j in range(6):

            for d in range(6, 18, 2):

                ell = k + d

                actual = actual_delta5(
                    k,
                    ell,
                    j
                )

                candidate = proposed_universal(
                    5,
                    j,
                    k,
                    d
                )

                residual = simp(
                    actual - candidate
                )

                if residual != 0:

                    print(
                        f"FAIL "
                        f"k={k} j={j} d={d} "
                        f"actual={actual} "
                        f"candidate={candidate} "
                        f"residual={residual}"
                    )

                    failures += 1

    print()
    print(
        f"old universal-law failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# 7. TEST WHETHER j=0,1,2 DIFFER FROM j=3,4,5
# ------------------------------------------------------------------------------

def structural_split():

    print("=" * 78)
    print("7. STRUCTURAL SPLIT")
    print("=" * 78)

    for j in range(6):

        print(f"j={j}")

        for k in [3, 5, 7, 9, 11, 13]:

            ds = list(range(6, 18, 2))

            vals = [
                actual_delta5(
                    k,
                    k + d,
                    j
                )
                for d in ds
            ]

            print(
                f"  k={k}: {vals}"
            )

        print()


# ------------------------------------------------------------------------------
# 8. DIRECT CHECK OF THE r=5 TERMINAL LAW
# ------------------------------------------------------------------------------

def terminal_check():

    print("=" * 78)
    print("8. TERMINAL j=5 CHECK")
    print("=" * 78)

    failures = 0

    for k in [3, 5, 7, 9, 11, 13]:

        for d in [6, 8, 10, 12, 14, 16]:

            actual = actual_delta5(
                k,
                k + d,
                5
            )

            expected = proposed_universal(
                5,
                5,
                k,
                d
            )

            residual = simp(
                actual - expected
            )

            print(
                f"k={k:2d} d={d:2d} "
                f"actual={actual} "
                f"expected={expected} "
                f"residual={residual} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
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

    raw_data()

    degree_audit()

    polys = interpolate_actual()

    coefficient_table(polys)

    root_audit(polys)

    compare_old_law()

    structural_split()

    terminal_check()

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The experiment derives the ACTUAL r=5 boundary ladder"
    )
    print(
        "directly from the exact pq kernel."
    )
    print()
    print(
        "No universal r,j formula is assumed."
    )
    print(
        "The old r<=4 candidate is tested only after"
    )
    print(
        "the actual r=5 polynomials have been extracted."
    )
    print()
    print(
        "Pay particular attention to:"
    )
    print(
        "  j=0, j=1, j=2"
    )
    print(
        "versus"
    )
    print(
        "  j=3, j=4, j=5."
    )


if __name__ == "__main__":
    main()

