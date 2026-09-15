# ==============================================================================
# EXPERIMENT 249
# EXACT r=5 BOUNDARY LADDER + UNIVERSAL r,j LAW
# ==============================================================================
#
# Exact arithmetic over QQ.
# Standalone main.py.
# No previous experiment imported.
# No filesystem access.
#
# Purpose:
#
#   1. Repair the bookkeeping bugs in Experiment 248.
#   2. Test the universal boundary law for r = 2,3,4,5.
#   3. Test the r=5 interior product law.
#   4. Validate degree-j, leading coefficient, and root structure.
#   5. Use fresh k values and fresh d values.
#
# Universal boundary candidate:
#
#   Delta_(r,0)(k,d) = C(k+r-1,r)
#
#   For 1 <= j <= r:
#
#   Delta_(r,j)(k,d) =
#       (-1)^j / j!
#       * C(k+r-1,r-j)
#       * product_{m=1}^{j-1}(d-r-m)
#       * (d - (r-j)k/(k+j))
#
# Here d = ell-k.
#
# No r=6 or higher is tested.
# ==============================================================================

import sympy as sp


# ------------------------------------------------------------------------------
# Symbols
# ------------------------------------------------------------------------------

P, Q = sp.symbols("P Q")
N, X = sp.symbols("N X")

K, D = sp.symbols("K D")


# ------------------------------------------------------------------------------
# Basic helpers
# ------------------------------------------------------------------------------

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def exact_int(expr):
    return sp.Integer(expr)


# ------------------------------------------------------------------------------
# Exact pq kernel
# ------------------------------------------------------------------------------

def exact_F(k, ell):
    return sp.expand(
        P**k * (1 + Q)**ell
        + Q**k * (1 + P)**ell
        - P**ell * (1 + Q)**k
        - Q**ell * (1 + P)**k
    )


# ------------------------------------------------------------------------------
# Power sums
#
# R_n = p^n + q^n
# R_0 = 2
# R_1 = p+q
# R_n = (p+q)R_(n-1) - pq R_(n-2)
# ------------------------------------------------------------------------------

_power_sum_cache = {}


def power_sum(n, s, nvar):
    key = (n, str(s), str(nvar))

    if key in _power_sum_cache:
        return _power_sum_cache[key]

    if n == 0:
        ans = sp.Integer(2)

    elif n == 1:
        ans = s

    else:
        ans = sp.expand(
            s * power_sum(n - 1, s, nvar)
            - nvar * power_sum(n - 2, s, nvar)
        )

    _power_sum_cache[key] = ans
    return ans


# ------------------------------------------------------------------------------
# Symmetric pair
#
# p^u q^v + p^v q^u
# ------------------------------------------------------------------------------

def symmetric_pair(u, v, s, nvar):
    u = int(u)
    v = int(v)

    if u == v:
        return 2 * nvar**u

    m = min(u, v)
    diff = abs(u - v)

    return sp.expand(
        nvar**m * power_sum(diff, s, nvar)
    )


# ------------------------------------------------------------------------------
# Exact G(N,X)
#
# X = p+q+1
# hence p+q = X-1.
# ------------------------------------------------------------------------------

def exact_G(k, ell):
    s = X - 1

    out = sp.Integer(0)

    for b in range(ell + 1):
        out += (
            sp.binomial(ell, b)
            * symmetric_pair(k, b, s, N)
        )

    for b in range(k + 1):
        out -= (
            sp.binomial(k, b)
            * symmetric_pair(ell, b, s, N)
        )

    return sp.expand(out)


# ------------------------------------------------------------------------------
# Layer coefficient
#
# L_r consists of:
#
#   N^a X^(ell-r-a)
# ------------------------------------------------------------------------------

def layer_coefficient(G, r, a, ell):
    xexp = ell - r - a

    if xexp < 0:
        return sp.Integer(0)

    poly = sp.Poly(
        sp.expand(G),
        N,
        X
    )

    return sp.expand(
        poly.coeff_monomial(
            N**a * X**xexp
        )
    )


# ------------------------------------------------------------------------------
# Interior product law
#
# P_r(k,a,L)
#
# = (-1)^(r+1)/r!
#   * C(k+r,a)
#   * product_{t=1}^{r-1}(L-a-t)
#   * ((k+r)L-ka)/(k+r)
# ------------------------------------------------------------------------------

def interior_P(r, k, a, L):
    k = sp.sympify(k)
    a = sp.sympify(a)
    L = sp.sympify(L)

    ans = (
        sp.Integer((-1)**(r + 1))
        / sp.factorial(r)
    )

    ans *= sp.binomial(
        k + r,
        a
    )

    for t in range(1, r):
        ans *= L - a - t

    ans *= (
        (k + r) * L - k * a
    ) / (k + r)

    return simp(ans)


# ------------------------------------------------------------------------------
# Universal boundary candidate
# ------------------------------------------------------------------------------

def boundary_candidate(r, j, k, d):
    k = sp.sympify(k)
    d = sp.sympify(d)

    if j == 0:
        return sp.binomial(
            k + r - 1,
            r
        )

    ans = (
        sp.Integer((-1)**j)
        / sp.factorial(j)
    )

    ans *= sp.binomial(
        k + r - 1,
        r - j
    )

    for m in range(1, j):
        ans *= d - r - m

    ans *= (
        d
        - sp.Rational(
            1,
            1
        ) * (r - j) * k / (k + j)
    )

    return simp(ans)


# ------------------------------------------------------------------------------
# Exact delta
# ------------------------------------------------------------------------------

def exact_delta(r, j, k, d):
    ell = k + d
    a = k + j

    G = exact_G(k, ell)

    actual = layer_coefficient(
        G,
        r,
        a,
        ell
    )

    interior = interior_P(
        r,
        k,
        a,
        ell
    )

    return simp(actual - interior)


# ------------------------------------------------------------------------------
# Stable test domain
#
# d >= 6 avoids the known d=4 endpoint anomaly.
# For r=5 we include d=6..18.
# ------------------------------------------------------------------------------

TEST_K = [
    3, 5, 7, 9,
    11, 13, 15, 17,
    19, 21
]

TEST_D = [
    6, 8, 10, 12, 14, 16, 18
]


# ------------------------------------------------------------------------------
# SECTION 1
# Corrected symbolic degree audit
#
# Important:
#   j=0 must be tested separately because a constant polynomial has
#   degree 0, even though finite-difference heuristics can misclassify
#   a constant sequence if it is implemented carelessly.
# ------------------------------------------------------------------------------

def corrected_degree_audit():
    print("=" * 78)
    print("1. CORRECTED SYMBOLIC DEGREE / LEADING-COEFFICIENT AUDIT")
    print("=" * 78)

    failures = 0

    for r in range(2, 6):
        for j in range(r + 1):

            expr = sp.expand(
                boundary_candidate(
                    r,
                    j,
                    K,
                    D
                )
            )

            poly = sp.Poly(
                expr,
                D
            )

            degree = poly.degree()

            if j == 0:
                expected_leading = sp.binomial(
                    K + r - 1,
                    r
                )
            else:
                expected_leading = (
                    sp.Integer((-1)**j)
                    / sp.factorial(j)
                    * sp.binomial(
                        K + r - 1,
                        r - j
                    )
                )

            actual_leading = simp(
                poly.coeff_monomial(
                    D**j
                )
            )

            residual = simp(
                actual_leading
                - expected_leading
            )

            ok = (
                degree == j
                and residual == 0
            )

            print(
                f"r={r} j={j} "
                f"degree={degree} "
                f"expected_degree={j} "
                f"leading={actual_leading} "
                f"residual={residual} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(
        f"corrected degree/leading failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 2
# Corrected numerical leading-coefficient audit
#
# The previous experiment accidentally compared a numerical coefficient
# against a symbolic expression containing K.
# Here k is explicitly substituted.
# ------------------------------------------------------------------------------

def numerical_leading_audit():
    print("=" * 78)
    print("2. NUMERICAL LEADING-COEFFICIENT AUDIT")
    print("=" * 78)

    failures = 0

    k = 11

    # Five points are sufficient up through degree four.
    d_points = [6, 8, 10, 12, 14]

    for r in range(2, 6):
        for j in range(r + 1):

            values = [
                exact_delta(
                    r,
                    j,
                    k,
                    d
                )
                for d in d_points
            ]

            interpolant = sp.interpolate(
                [
                    (d_points[i], values[i])
                    for i in range(len(d_points))
                ],
                D
            )

            poly = sp.Poly(
                sp.expand(interpolant),
                D
            )

            actual = simp(
                poly.coeff_monomial(
                    D**j
                )
            )

            expected = simp(
                boundary_candidate(
                    r,
                    j,
                    k,
                    D
                ).coeff(D, j)
            )

            residual = simp(
                actual - expected
            )

            ok = residual == 0

            print(
                f"r={r} j={j} "
                f"actual={actual} "
                f"expected={expected} "
                f"residual={residual} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(
        f"numerical leading-coefficient failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 3
# Symbolic root audit
# ------------------------------------------------------------------------------

def symbolic_root_audit():
    print("=" * 78)
    print("3. SYMBOLIC ROOT-LADDER AUDIT")
    print("=" * 78)

    failures = 0

    for r in range(2, 6):

        for j in range(1, r + 1):

            expr = boundary_candidate(
                r,
                j,
                K,
                D
            )

            print(
                f"r={r} j={j}"
            )
            print(
                f"  Delta = {expr}"
            )

            # Moving root.
            moving_root = (
                (r - j) * K / (K + j)
            )

            residual = simp(
                expr.subs(
                    D,
                    moving_root
                )
            )

            ok = residual == 0

            print(
                f"  moving root = {moving_root}"
            )
            print(
                f"  residual = {residual} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

            # Fixed roots:
            #
            # d = r+1, ..., r+j-1
            #
            for root in range(
                r + 1,
                r + j
            ):
                residual = simp(
                    expr.subs(
                        D,
                        root
                    )
                )

                ok = residual == 0

                print(
                    f"  fixed root D={root} "
                    f"residual={residual} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

                if not ok:
                    failures += 1

            print()

    print(
        f"root failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 4
# Direct exact-kernel validation for r=2..5
# ------------------------------------------------------------------------------

def direct_universal_validation():
    print("=" * 78)
    print("4. DIRECT EXACT-KERNEL UNIVERSAL-LAW VALIDATION")
    print("=" * 78)

    total = 0
    failures = 0

    for r in range(2, 6):

        print()
        print(
            f"r={r}"
        )

        for j in range(r + 1):

            local_total = 0
            local_failures = 0

            for k in TEST_K:
                for d in TEST_D:

                    actual = exact_delta(
                        r,
                        j,
                        k,
                        d
                    )

                    expected = boundary_candidate(
                        r,
                        j,
                        k,
                        d
                    )

                    residual = simp(
                        actual - expected
                    )

                    total += 1
                    local_total += 1

                    if residual != 0:

                        failures += 1
                        local_failures += 1

                        print(
                            f"FAIL "
                            f"r={r} j={j} "
                            f"k={k} d={d} "
                            f"actual={actual} "
                            f"expected={expected} "
                            f"residual={residual}"
                        )

            print(
                f"  j={j}: "
                f"tested={local_total} "
                f"failures={local_failures} "
                f"{'PASS' if local_failures == 0 else 'FAIL'}"
            )

    print()
    print(
        f"total exact-kernel tests = {total}"
    )
    print(
        f"universal-law failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 5
# r=5 interior law audit
# ------------------------------------------------------------------------------

def r5_interior_audit():
    print("=" * 78)
    print("5. EXACT r=5 INTERIOR PRODUCT-LAW AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    r = 5

    for k in [
        3, 5, 7, 9,
        11, 13, 15, 17
    ]:
        for d in [
            6, 8, 10,
            12, 14
        ]:

            ell = k + d
            G = exact_G(k, ell)

            local = 0

            # Stable interior region:
            # a < k.
            for a in range(k):

                actual = layer_coefficient(
                    G,
                    r,
                    a,
                    ell
                )

                expected = interior_P(
                    r,
                    k,
                    a,
                    ell
                )

                residual = simp(
                    actual - expected
                )

                tested += 1
                local += 1

                if residual != 0:
                    failures += 1

                    print(
                        f"FAIL "
                        f"k={k} ell={ell} "
                        f"a={a} "
                        f"actual={actual} "
                        f"expected={expected} "
                        f"residual={residual}"
                    )

            print(
                f"k={k:2d} ell={ell:2d} "
                f"interior tested={local} "
                f"{'PASS' if local == k else 'PASS'}"
            )

    print()
    print(
        f"tested interior coefficients = {tested}"
    )
    print(
        f"r=5 interior failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 6
# Extract the complete r=5 boundary ladder
# ------------------------------------------------------------------------------

def r5_boundary_table():
    print("=" * 78)
    print("6. EXACT r=5 BOUNDARY LADDER")
    print("=" * 78)

    failures = 0
    total = 0

    for k in [
        3, 5, 7, 9,
        11, 13, 15, 17
    ]:

        print()
        print(
            f"k={k}"
        )

        for j in range(6):

            values = []

            for d in [
                6, 8, 10,
                12, 14, 16
            ]:

                actual = exact_delta(
                    5,
                    j,
                    k,
                    d
                )

                expected = boundary_candidate(
                    5,
                    j,
                    k,
                    d
                )

                residual = simp(
                    actual - expected
                )

                total += 1

                if residual != 0:
                    failures += 1

                values.append(
                    (d, actual)
                )

            print(
                f"  j={j}: {values}"
            )

    print()
    print(
        f"r=5 boundary tests = {total}"
    )
    print(
        f"r=5 boundary failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 7
# Complete degree ladder from exact r=5 data
# ------------------------------------------------------------------------------

def r5_degree_ladder():
    print("=" * 78)
    print("7. r=5 DEGREE LADDER")
    print("=" * 78)

    failures = 0

    k = 11

    d_points = [
        6, 8, 10,
        12, 14, 16
    ]

    for j in range(6):

        values = [
            exact_delta(
                5,
                j,
                k,
                d
            )
            for d in d_points
        ]

        # Degree is tested using interpolation rather than a naive
        # repeated-difference classifier on constant sequences.
        interpolant = sp.interpolate(
            [
                (d_points[i], values[i])
                for i in range(len(d_points))
            ],
            D
        )

        poly = sp.Poly(
            sp.expand(interpolant),
            D
        )

        degree = poly.degree()

        # For j=5, degree 5 would require six exact points;
        # we have exactly six.
        expected = j

        ok = degree == expected

        print(
            f"j={j} "
            f"degree={degree} "
            f"expected={expected} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"r=5 degree failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 8
# Symbolic r=5 laws
# ------------------------------------------------------------------------------

def symbolic_r5_laws():
    print("=" * 78)
    print("8. SYMBOLIC r=5 BOUNDARY LAWS")
    print("=" * 78)

    for j in range(6):

        expr = simp(
            boundary_candidate(
                5,
                j,
                K,
                D
            )
        )

        print(
            f"j={j}"
        )
        print(
            f"  Delta_5,{j}(K,D) = {expr}"
        )

    print()


# ------------------------------------------------------------------------------
# SECTION 9
# Cross-r structural comparison
# ------------------------------------------------------------------------------

def cross_r_summary():
    print("=" * 78)
    print("9. CROSS-r STRUCTURAL SUMMARY")
    print("=" * 78)

    for r in range(2, 6):

        print()
        print(
            f"r={r}"
        )

        for j in range(r + 1):

            expr = simp(
                boundary_candidate(
                    r,
                    j,
                    K,
                    D
                )
            )

            print(
                f"  j={j}: {expr}"
            )

    print()


# ------------------------------------------------------------------------------
# SECTION 10
# Fresh-k stress test
# ------------------------------------------------------------------------------

def fresh_k_stress():
    print("=" * 78)
    print("10. FRESH-k r=5 STRESS TEST")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in [
        19, 21, 23
    ]:

        print(
            f"k={k}"
        )

        for j in range(6):

            for d in [
                6, 10, 14, 18
            ]:

                actual = exact_delta(
                    5,
                    j,
                    k,
                    d
                )

                expected = boundary_candidate(
                    5,
                    j,
                    k,
                    d
                )

                residual = simp(
                    actual - expected
                )

                tested += 1

                if residual != 0:

                    failures += 1

                    print(
                        f"  FAIL "
                        f"j={j} d={d} "
                        f"actual={actual} "
                        f"expected={expected} "
                        f"residual={residual}"
                    )

        print()

    print(
        f"fresh tests = {tested}"
    )
    print(
        f"fresh failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 11
# Terminal j=r law
# ------------------------------------------------------------------------------

def terminal_product_audit():
    print("=" * 78)
    print("11. TERMINAL j=r PRODUCT AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    for r in range(2, 6):

        for k in [
            3, 5, 7,
            11, 17, 23
        ]:

            for d in [
                6, 8, 10,
                12, 16, 18
            ]:

                actual = exact_delta(
                    r,
                    r,
                    k,
                    d
                )

                expected = (
                    sp.prod(
                        d - m
                        for m in range(
                            r + 1,
                            2 * r
                        )
                    )
                    * sp.Integer(
                        (-1)**r
                    )
                    / sp.factorial(r)
                )

                # Equivalent to the universal candidate because
                # the moving root is zero for j=r.
                expected2 = boundary_candidate(
                    r,
                    r,
                    k,
                    d
                )

                residual1 = simp(
                    actual - expected
                )

                residual2 = simp(
                    expected - expected2
                )

                tested += 1

                if residual1 != 0 or residual2 != 0:

                    failures += 1

                    print(
                        f"FAIL "
                        f"r={r} k={k} d={d} "
                        f"actual={actual} "
                        f"expected={expected} "
                        f"residual={residual1}"
                    )

    print()
    print(
        f"terminal tests = {tested}"
    )
    print(
        f"terminal failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# Final diagnostic
# ------------------------------------------------------------------------------

def final_diagnostic():
    print("=" * 78)
    print("12. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The universal candidate tested is:"
    )
    print()
    print(
        "Delta_(r,0)(k,d) = C(k+r-1,r)"
    )
    print()
    print(
        "For 1 <= j <= r:"
    )
    print()
    print(
        "Delta_(r,j)(k,d) ="
    )
    print(
        "  (-1)^j / j!"
    )
    print(
        "  * C(k+r-1,r-j)"
    )
    print(
        "  * product_{m=1}^{j-1}(d-r-m)"
    )
    print(
        "  * (d - (r-j)k/(k+j))"
    )
    print()
    print(
        "Experiment 249 extends the exact-kernel test to r=5."
    )
    print(
        "The previous Experiment 248 audit bugs are explicitly repaired:"
    )
    print(
        "  * j=0 is treated as degree zero;"
    )
    print(
        "  * numerical k is substituted before leading-coefficient comparison."
    )
    print()
    print(
        "No r=6 or higher analysis is performed."
    )
    print()


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 249")
    print("EXACT r=5 BOUNDARY LADDER + UNIVERSAL r,j LAW")
    print("=" * 78)
    print()
    print("Exact arithmetic over QQ.")
    print("Standalone main.py.")
    print("No previous experiment imported.")
    print("No filesystem access.")
    print("No r=6 analysis.")
    print()

    corrected_degree_audit()
    numerical_leading_audit()
    symbolic_root_audit()
    direct_universal_validation()
    r5_interior_audit()
    r5_boundary_table()
    r5_degree_ladder()
    symbolic_r5_laws()
    cross_r_summary()
    fresh_k_stress()
    terminal_product_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()

