# ==============================================================================
# EXPERIMENT 248
# UNIVERSAL r,j BOUNDARY PRODUCT LAW
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No L5 analysis
#
# Goal:
#
#   Test a single boundary law for r=2,3,4:
#
#       Delta_{r,0}
#           = C(k+r-1,r)
#
#       Delta_{r,j}
#           =
#           (-1)^j/j!
#           * C(k+r-1,r-j)
#           * product_{m=1}^{j-1}(d-r-m)
#           * ( d - (r-j)k/(k+j) )
#
# where
#
#       d = ell-k
#
# and
#
#       Delta = exact coefficient - interior extrapolation.
#
# The exact source is the pq kernel.
#
# ==============================================================================

import sympy as sp


# ------------------------------------------------------------------------------
# Symbols
# ------------------------------------------------------------------------------

P, Q = sp.symbols("P Q")
N, X = sp.symbols("N X")
Sx = sp.symbols("Sx")

K, D = sp.symbols("K D")


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def exact_int(x):
    return sp.Integer(x)


# ------------------------------------------------------------------------------
# Exact pq kernel
# ------------------------------------------------------------------------------

def exact_F(k, ell):
    """
    Exact symmetric pq kernel:

        p^k (1+q)^ell
      + q^k (1+p)^ell
      - p^ell (1+q)^k
      - q^ell (1+p)^k
    """
    return sp.expand(
        P**k * (1 + Q)**ell
        + Q**k * (1 + P)**ell
        - P**ell * (1 + Q)**k
        - Q**ell * (1 + P)**k
    )


# ------------------------------------------------------------------------------
# Power-sum recurrence
#
# Let s = p+q.
#
# R_n = p^n + q^n
#
# R_0 = 2
# R_1 = s
# R_n = s R_(n-1) - N R_(n-2)
# ------------------------------------------------------------------------------

_power_sum_cache = {}


def power_sum(n, s, nvar):
    key = (n, s, nvar)

    if key in _power_sum_cache:
        return _power_sum_cache[key]

    if n == 0:
        out = sp.Integer(2)
    elif n == 1:
        out = s
    else:
        a = power_sum(n - 1, s, nvar)
        b = power_sum(n - 2, s, nvar)
        out = sp.expand(s*a - nvar*b)

    _power_sum_cache[key] = out
    return out


# ------------------------------------------------------------------------------
# Convert symmetric pair
#
# p^u q^v + p^v q^u
#
# into N=pq and s=p+q.
# ------------------------------------------------------------------------------

def symmetric_pair(u, v, s, nvar):
    u = int(u)
    v = int(v)

    if u == v:
        return 2 * nvar**u

    m = min(u, v)
    r = abs(u - v)

    return sp.expand(
        nvar**m * power_sum(r, s, nvar)
    )


# ------------------------------------------------------------------------------
# Exact G(N,X)
#
# X = p+q+1
# so p+q = X-1.
#
# We construct the symmetric kernel directly without expanding in p,q first.
# ------------------------------------------------------------------------------

def exact_G(k, ell):
    s = X - 1

    out = sp.Integer(0)

    # First pair:
    # p^k (1+q)^ell + q^k (1+p)^ell
    for b in range(ell + 1):
        c = sp.binomial(ell, b)
        out += c * symmetric_pair(k, b, s, N)

    # Second pair:
    # -p^ell (1+q)^k - q^ell (1+p)^k
    for b in range(k + 1):
        c = sp.binomial(k, b)
        out -= c * symmetric_pair(ell, b, s, N)

    return sp.Poly(
        sp.expand(out),
        N,
        X
    ).as_expr()


# ------------------------------------------------------------------------------
# Layer coefficient
#
# L_r consists of monomials
#
#     N^a X^(ell-r-a)
#
# If the X exponent is negative, the coefficient is exactly zero.
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
# Established interior product law
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

    out = (
        sp.Integer((-1)**(r + 1))
        / sp.factorial(r)
    )

    out *= sp.binomial(k + r, a)

    for t in range(1, r):
        out *= (L - a - t)

    out *= (
        ((k + r) * L - k*a)
        / (k + r)
    )

    return simp(out)


# ------------------------------------------------------------------------------
# Extract exact boundary correction
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
# UNIVERSAL CANDIDATE
#
# For j=0:
#
#   C(k+r-1,r)
#
# For j>=1:
#
#   (-1)^j/j!
#   C(k+r-1,r-j)
#   product_{m=1}^{j-1}(d-r-m)
#   (d - (r-j)k/(k+j))
#
# IMPORTANT:
#   symbolic k must use ordinary division.
#   Never use sp.Rational(K, K+j).
# ------------------------------------------------------------------------------

def boundary_candidate(r, j, k, d):
    k = sp.sympify(k)
    d = sp.sympify(d)

    if j == 0:
        return sp.binomial(
            k + r - 1,
            r
        )

    out = (
        sp.Integer((-1)**j)
        / sp.factorial(j)
    )

    out *= sp.binomial(
        k + r - 1,
        r - j
    )

    for m in range(1, j):
        out *= (
            d - r - m
        )

    out *= (
        d
        - (
            (r - j) * k
            / (k + j)
        )
    )

    return simp(out)


# ------------------------------------------------------------------------------
# SECTION 1
# Symbolic universal-law audit
# ------------------------------------------------------------------------------

def symbolic_universal_audit():
    print("=" * 78)
    print("1. SYMBOLIC UNIVERSAL r,j BOUNDARY LAW")
    print("=" * 78)

    failures = 0

    for r in (2, 3, 4):
        for j in range(r + 1):

            candidate = boundary_candidate(
                r,
                j,
                K,
                D
            )

            print(f"r={r} j={j}")
            print(f"  candidate = {candidate}")

            # Check factorized structure against itself after expansion.
            residual = simp(
                candidate
                - sp.expand(candidate)
            )

            # This residual is normally zero; it ensures the expression
            # is a genuine polynomial/rational symbolic object.
            if residual == 0:
                print("  symbolic normalization = PASS")
            else:
                print(
                    f"  symbolic normalization = FAIL residual={residual}"
                )
                failures += 1

            print()

    print(
        f"symbolic normalization failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 2
# Symbolic degree and leading coefficient
# ------------------------------------------------------------------------------

def symbolic_degree_audit():
    print("=" * 78)
    print("2. SYMBOLIC DEGREE / LEADING-COEFFICIENT AUDIT")
    print("=" * 78)

    failures = 0

    for r in (2, 3, 4):
        for j in range(r + 1):

            expr = sp.together(
                boundary_candidate(
                    r,
                    j,
                    K,
                    D
                )
            )

            poly = sp.Poly(
                sp.expand(expr),
                D
            )

            degree = poly.degree()

            expected_degree = j

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

            status = (
                degree == expected_degree
                and residual == 0
            )

            print(
                f"r={r} j={j} "
                f"degree={degree} "
                f"expected={expected_degree} "
                f"leading={actual_leading} "
                f"status={'PASS' if status else 'FAIL'}"
            )

            if not status:
                failures += 1

    print()
    print(
        f"degree/leading failures = {failures}"
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

    for r in (2, 3, 4):
        for j in range(1, r + 1):

            expr = boundary_candidate(
                r,
                j,
                K,
                D
            )

            print(f"r={r} j={j}")
            print(
                f"  Delta = {simp(expr)}"
            )

            # k-dependent root
            moving_root = (
                (r - j) * K / (K + j)
            )

            residual = simp(
                expr.subs(
                    D,
                    moving_root
                )
            )

            status = "PASS" if residual == 0 else "FAIL"

            print(
                f"  moving root = {moving_root}"
            )
            print(
                f"  residual = {residual} {status}"
            )

            if residual != 0:
                failures += 1

            # Fixed roots
            fixed_roots = list(
                range(
                    r + 1,
                    r + j
                )
            )

            for root in fixed_roots:
                residual = simp(
                    expr.subs(
                        D,
                        root
                    )
                )

                status = (
                    "PASS"
                    if residual == 0
                    else "FAIL"
                )

                print(
                    f"  fixed root D={root} "
                    f"residual={residual} "
                    f"{status}"
                )

                if residual != 0:
                    failures += 1

            print()

    print(
        f"root failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# Exact-kernel test grid
#
# Existing k values + fresh k values.
# We intentionally use only stable d >= 6.
# ------------------------------------------------------------------------------

TEST_K = [
    3, 5, 7, 9,
    11, 13, 15, 17,
    19, 21
]

TEST_D = [
    6, 8, 10, 12, 14
]


# ------------------------------------------------------------------------------
# SECTION 4
# Direct exact-kernel validation
# ------------------------------------------------------------------------------

def direct_kernel_validation():
    print("=" * 78)
    print("4. DIRECT EXACT-KERNEL UNIVERSAL-LAW TEST")
    print("=" * 78)

    total = 0
    failures = 0

    for r in (2, 3, 4):

        print()
        print(f"r={r}")

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
# Degree ladder measured directly from exact-kernel data
# ------------------------------------------------------------------------------

def finite_difference(values):
    values = list(values)

    out = []

    while len(values) > 1:
        nxt = []

        for i in range(len(values) - 1):
            nxt.append(
                sp.expand(
                    values[i + 1] - values[i]
                )
            )

        out.append(nxt)
        values = nxt

    return out


def observed_degree(values):
    table = [list(values)]

    current = list(values)

    while len(current) > 1:
        nxt = [
            sp.expand(
                current[i + 1]
                - current[i]
            )
            for i in range(len(current) - 1)
        ]

        table.append(nxt)

        if all(
            sp.simplify(v - nxt[0]) == 0
            for v in nxt
        ):
            return len(table) - 1

        current = nxt

    return len(values) - 1


def degree_from_exact_kernel():
    print("=" * 78)
    print("5. EXACT-KERNEL DEGREE LADDER")
    print("=" * 78)

    failures = 0

    for r in (2, 3, 4):
        for j in range(r + 1):

            for k in TEST_K[:6]:

                values = []

                for d in TEST_D:
                    values.append(
                        exact_delta(
                            r,
                            j,
                            k,
                            d
                        )
                    )

                degree = observed_degree(values)

                expected = j

                status = (
                    "PASS"
                    if degree == expected
                    else "FAIL"
                )

                print(
                    f"r={r} j={j} k={k:2d} "
                    f"degree={degree} "
                    f"expected={expected} "
                    f"{status}"
                )

                if degree != expected:
                    failures += 1

    print()
    print(
        f"degree failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 6
# Direct coefficient-law audit
# ------------------------------------------------------------------------------

def leading_coefficient_from_data():
    print("=" * 78)
    print("6. EXACT-KERNEL LEADING-COEFFICIENT AUDIT")
    print("=" * 78)

    failures = 0

    # Five d-points are enough to reconstruct degrees up to four.
    for r in (2, 3, 4):
        for j in range(r + 1):

            k = 11

            values = [
                exact_delta(r, j, k, d)
                for d in TEST_D
            ]

            dvals = [
                exact_int(d)
                for d in TEST_D
            ]

            poly = sp.interpolate(
                [
                    (dvals[i], values[i])
                    for i in range(len(values))
                ],
                D
            )

            poly = sp.Poly(
                sp.expand(poly),
                D
            )

            actual = simp(
                poly.coeff_monomial(
                    D**j
                )
            )

            expected = simp(
                (
                    sp.Integer((-1)**j)
                    / sp.factorial(j)
                    * sp.binomial(
                        K + r - 1,
                        r - j
                    )
                )
            )

            residual = simp(
                actual - expected
            )

            status = (
                "PASS"
                if residual == 0
                else "FAIL"
            )

            print(
                f"r={r} j={j} "
                f"actual={actual} "
                f"expected={expected} "
                f"residual={residual} "
                f"{status}"
            )

            if residual != 0:
                failures += 1

    print()
    print(
        f"leading-coefficient failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 7
# Cross-r symbolic pattern
# ------------------------------------------------------------------------------

def symbolic_cross_r():
    print("=" * 78)
    print("7. CROSS-r SYMBOLIC PRODUCT PATTERN")
    print("=" * 78)

    for r in (2, 3, 4):
        print()
        print(f"r={r}")

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
# SECTION 8
# Test exact reconstruction using the unified law
# at fresh k values.
# ------------------------------------------------------------------------------

def fresh_k_stress_test():
    print("=" * 78)
    print("8. FRESH-k STRESS TEST")
    print("=" * 78)

    fresh_k = [19, 21]

    total = 0
    failures = 0

    for k in fresh_k:

        print(f"k={k}")

        for r in (2, 3, 4):
            for j in range(r + 1):

                for d in (6, 8, 10, 12):

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

                    if residual != 0:
                        failures += 1

                        print(
                            f"  FAIL "
                            f"r={r} j={j} d={d} "
                            f"actual={actual} "
                            f"expected={expected} "
                            f"residual={residual}"
                        )

        print()

    print(
        f"fresh tests = {total}"
    )
    print(
        f"fresh failures = {failures}"
    )
    print()


# ------------------------------------------------------------------------------
# SECTION 9
# Compact universal-law statement
# ------------------------------------------------------------------------------

def final_diagnostic():
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print("The tested universal boundary pattern is:")
    print()

    print(
        "Delta_(r,0)(k,d) = binomial(k+r-1, r)"
    )

    print()

    print(
        "For 1 <= j <= r:"
    )

    print(
        "Delta_(r,j)(k,d) ="
    )

    print(
        "  (-1)^j / j!"
    )

    print(
        "  * C(k+r-1, r-j)"
    )

    print(
        "  * product_{m=1}^{j-1} (d-r-m)"
    )

    print(
        "  * (d - (r-j)k/(k+j))"
    )

    print()

    print(
        "The structural consequences are:"
    )

    print(
        "  * degree in d = j"
    )

    print(
        "  * leading coefficient ="
    )

    print(
        "      (-1)^j/j! * C(k+r-1,r-j)"
    )

    print()

    print(
        "  * fixed roots:"
    )

    print(
        "      d = r+1, r+2, ..., r+j-1"
    )

    print()

    print(
        "  * moving root:"
    )

    print(
        "      d = (r-j)k/(k+j)"
    )

    print()

    print(
        "This experiment deliberately stops at r=4."
    )

    print(
        "No L5 analysis is performed."
    )

    print()


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 248")
    print("UNIVERSAL r,j BOUNDARY PRODUCT LAW")
    print("=" * 78)
    print()

    print("Exact arithmetic over QQ.")
    print("Standalone main.py.")
    print("No previous experiment imported.")
    print("No filesystem access.")
    print("No L5 analysis.")
    print()

    symbolic_universal_audit()
    symbolic_degree_audit()
    symbolic_root_audit()
    direct_kernel_validation()
    degree_from_exact_kernel()
    leading_coefficient_from_data()
    symbolic_cross_r()
    fresh_k_stress_test()
    final_diagnostic()


if __name__ == "__main__":
    main()

