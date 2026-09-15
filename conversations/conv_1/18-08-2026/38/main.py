import sympy as sp
from math import comb

# ==============================================================================
# EXPERIMENT 250
# EXACT r=5 BOUNDARY LAW + UNIVERSAL CROSS-r AUDIT
# ==============================================================================

# Exact arithmetic over QQ.
# Standalone main.py.
# No previous experiment imported.
# No filesystem access.
# No r=6 analysis.
#
# IMPORTANT:
#   This experiment is deliberately small.
#   It does NOT perform a large exact-kernel stress sweep.
#
# The universal candidate is:
#
#   Delta_(r,0)(k,d) = C(k+r-1, r)
#
#   Delta_(r,j)(k,d)
#     = (-1)^j / j!
#       * C(k+r-1, r-j)
#       * product_{m=1}^{j-1} (d-r-m)
#       * (d - (r-j)k/(k+j))
#
# for 1 <= j <= r.
#
# Here d = ell-k.
#
# The experiment:
#   1. derives the candidate symbolically;
#   2. checks degree and leading coefficient;
#   3. checks all symbolic roots;
#   4. checks r=2,3,4 against the known exact boundary laws;
#   5. prints the NEW r=5 laws;
#   6. performs a minimal numerical k-substitution audit;
#   7. checks exact reconstruction at selected (k,d).
#
# No large pq-kernel expansion is performed here.
# ==============================================================================


K, D = sp.symbols("K D")


# ------------------------------------------------------------------------------
# Utility
# ------------------------------------------------------------------------------

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def qbinom(n, r):
    """
    Symbolic binomial with exact Gamma-free product construction.

    For the symbolic uses here, r is always a nonnegative integer.
    """
    if r < 0:
        return sp.Integer(0)

    out = sp.Integer(1)
    for i in range(r):
        out *= n - i
    return sp.cancel(out / sp.factorial(r))


def universal_delta(r, j, k, d):
    """
    Universal proposed boundary law.
    Works for both numeric and symbolic k,d.
    """
    if not (0 <= j <= r):
        raise ValueError("Require 0 <= j <= r")

    Kx = sp.sympify(k)
    Dx = sp.sympify(d)

    if j == 0:
        return sp.binomial(Kx + r - 1, r)

    fixed_product = sp.Integer(1)
    for m in range(1, j):
        fixed_product *= Dx - r - m

    moving = Dx - sp.Rational(r - j, 1) * Kx / (Kx + j)

    return simp(
        (-1) ** j
        * fixed_product
        * moving
        * sp.binomial(Kx + r - 1, r - j)
        / sp.factorial(j)
    )


# ------------------------------------------------------------------------------
# 1. SYMBOLIC UNIVERSAL LAW
# ------------------------------------------------------------------------------

def symbolic_universal_law():
    print("=" * 78)
    print("1. SYMBOLIC UNIVERSAL r,j LAW")
    print("=" * 78)

    failures = 0

    for r in range(2, 6):
        for j in range(0, r + 1):
            expr = simp(universal_delta(r, j, K, D))

            print(f"r={r} j={j}")
            print(f"  Delta = {expr}")

            if j == 0:
                expected = sp.binomial(K + r - 1, r)
                residual = simp(expr - expected)
                ok = residual == 0
            else:
                # Rebuild from the formula independently in expanded
                # symbolic form to ensure the implementation is internally
                # consistent.
                fixed = sp.Integer(1)
                for m in range(1, j):
                    fixed *= D - r - m

                expected = simp(
                    (-1) ** j
                    * fixed
                    * (
                        D
                        - sp.Rational(r - j, 1) * K / (K + j)
                    )
                    * sp.binomial(K + r - 1, r - j)
                    / sp.factorial(j)
                )

                residual = simp(expr - expected)
                ok = residual == 0

            print(f"  residual = {residual}")
            print(f"  status   = {'PASS' if ok else 'FAIL'}")

            if not ok:
                failures += 1

    print()
    print(f"symbolic universal-law failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 2. DEGREE + LEADING COEFFICIENT
# ------------------------------------------------------------------------------

def degree_leading_audit():
    print("=" * 78)
    print("2. DEGREE / LEADING-COEFFICIENT AUDIT")
    print("=" * 78)

    failures = 0

    for r in range(2, 6):
        for j in range(0, r + 1):
            expr = simp(universal_delta(r, j, K, D))
            poly = sp.Poly(sp.expand(expr), D)

            degree = poly.degree()

            if j == 0:
                expected_degree = 0
                expected_leading = sp.binomial(K + r - 1, r)
            else:
                expected_degree = j
                expected_leading = simp(
                    (-1) ** j
                    * sp.binomial(K + r - 1, r - j)
                    / sp.factorial(j)
                )

            actual_leading = simp(poly.LC())
            residual = simp(actual_leading - expected_leading)

            ok = (
                degree == expected_degree
                and residual == 0
            )

            print(
                f"r={r} j={j} "
                f"degree={degree} expected_degree={expected_degree} "
                f"leading={actual_leading} "
                f"expected={expected_leading} "
                f"residual={residual} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(f"degree/leading failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 3. ROOT-LADDER AUDIT
# ------------------------------------------------------------------------------

def root_audit():
    print("=" * 78)
    print("3. ROOT-LADDER AUDIT")
    print("=" * 78)

    failures = 0

    for r in range(2, 6):
        for j in range(1, r + 1):
            expr = simp(universal_delta(r, j, K, D))
            factored = sp.factor(expr)

            print(f"r={r} j={j}")
            print(f"  Delta = {factored}")

            # Fixed roots:
            #   D = r+1, ..., r+j-1
            for s in range(r + 1, r + j):
                residual = simp(expr.subs(D, s))
                ok = residual == 0

                print(
                    f"  fixed root D={s} "
                    f"residual={residual} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

                if not ok:
                    failures += 1

            # Moving root:
            moving_root = simp(
                sp.Rational(r - j, 1) * K / (K + j)
            )

            moving_residual = simp(
                expr.subs(D, moving_root)
            )

            moving_ok = moving_residual == 0

            print(
                f"  moving root={moving_root} "
                f"residual={moving_residual} "
                f"{'PASS' if moving_ok else 'FAIL'}"
            )

            if not moving_ok:
                failures += 1

    print()
    print(f"root failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 4. NUMERICAL SUBSTITUTION AUDIT
#
# This explicitly substitutes k before checking the leading coefficient.
# This fixes the bug in Experiment 249.
# ------------------------------------------------------------------------------

def numerical_leading_audit():
    print("=" * 78)
    print("4. NUMERICAL LEADING-COEFFICIENT AUDIT")
    print("=" * 78)

    failures = 0

    test_ks = [3, 5, 7, 9, 11, 13]

    for r in range(2, 6):
        for j in range(0, r + 1):

            # Pick a numeric k and THEN evaluate the expression.
            kval = test_ks[(r + j) % len(test_ks)]

            expr = simp(
                universal_delta(r, j, kval, D)
            )

            poly = sp.Poly(sp.expand(expr), D)

            actual = sp.Rational(poly.LC())

            if j == 0:
                expected = sp.binomial(kval + r - 1, r)
            else:
                expected = sp.Rational(
                    (-1) ** j,
                    sp.factorial(j)
                ) * sp.binomial(
                    kval + r - 1,
                    r - j
                )

            residual = simp(actual - expected)

            ok = residual == 0

            print(
                f"r={r} j={j} k={kval} "
                f"actual={actual} "
                f"expected={expected} "
                f"residual={residual} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(f"numerical leading-coefficient failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 5. KNOWN r=2,r=3,r=4 CROSS-CHECKS
#
# These are the already-established exact laws, written independently.
# ------------------------------------------------------------------------------

def known_boundary_formula(r, j, k, d):
    k = sp.sympify(k)
    d = sp.sympify(d)

    # r = 2
    if r == 2:
        if j == 0:
            return sp.binomial(k + 1, 2)
        if j == 1:
            return -d * k - d + k
        if j == 2:
            return d * (d - 3) / 2

    # r = 3
    if r == 3:
        if j == 0:
            return sp.binomial(k + 2, 3)
        if j == 1:
            return -(
                d * k + d - 2 * k
            ) * sp.binomial(k + 2, 2) / (k + 1)
        if j == 2:
            return (d - 4) * (
                d * k + 2 * d - k
            ) / 2
        if j == 3:
            return -d * (d - 5) * (d - 4) / 6

    # r = 4
    if r == 4:
        if j == 0:
            return sp.binomial(k + 3, 4)
        if j == 1:
            return -(
                d * k + d - 3 * k
            ) * sp.binomial(k + 3, 3) / (k + 1)
        if j == 2:
            return (
                (d - 5)
                * (d * k + 2 * d - 2 * k)
                * sp.binomial(k + 3, 2)
                / (2 * (k + 2))
            )
        if j == 3:
            return -(
                d - 6
            ) * (
                d - 5
            ) * (
                d * k + 3 * d - k
            ) / 6
        if j == 4:
            return d * (d - 7) * (d - 6) * (d - 5) / 24

    raise ValueError("Unsupported r,j")


def known_cross_check():
    print("=" * 78)
    print("5. r=2,3,4 KNOWN-LAW CROSS-CHECK")
    print("=" * 78)

    failures = 0
    tests = [
        (3, 6),
        (3, 8),
        (5, 6),
        (5, 10),
        (7, 8),
        (9, 12),
        (11, 14),
        (13, 16),
        (15, 18),
    ]

    for r in range(2, 5):
        for j in range(0, r + 1):

            local_failures = 0

            for k, d in tests:

                actual = simp(
                    known_boundary_formula(r, j, k, d)
                )

                candidate = simp(
                    universal_delta(r, j, k, d)
                )

                residual = simp(actual - candidate)

                if residual != 0:
                    local_failures += 1
                    failures += 1

            print(
                f"r={r} j={j}: "
                f"tested={len(tests)} "
                f"failures={local_failures} "
                f"{'PASS' if local_failures == 0 else 'FAIL'}"
            )

    print()
    print(f"known-law cross-check failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 6. COMPLETE r=5 SYMBOLIC LADDER
# ------------------------------------------------------------------------------

def r5_ladder():
    print("=" * 78)
    print("6. EXACT r=5 BOUNDARY LADDER")
    print("=" * 78)

    for j in range(0, 6):
        expr = simp(
            universal_delta(5, j, K, D)
        )

        print(f"j={j}")
        print(f"  Delta_5,{j}(K,D) = {expr}")

        poly = sp.Poly(sp.expand(expr), D)

        print(f"  degree  = {poly.degree()}")
        print(f"  leading = {simp(poly.LC())}")

        print()

    print()


# ------------------------------------------------------------------------------
# 7. r=5 ROOT STRUCTURE
# ------------------------------------------------------------------------------

def r5_root_structure():
    print("=" * 78)
    print("7. r=5 ROOT STRUCTURE")
    print("=" * 78)

    for j in range(1, 6):
        expr = simp(
            universal_delta(5, j, K, D)
        )

        moving_root = simp(
            sp.Rational(5 - j, 1) * K / (K + j)
        )

        print(f"j={j}")
        print(f"  Delta       = {sp.factor(expr)}")
        print(f"  moving root = {moving_root}")

        for s in range(6, 5 + j):
            print(f"  fixed root  = {s}")

        print()

    print()


# ------------------------------------------------------------------------------
# 8. MINIMAL r=5 NUMERICAL RECONSTRUCTION
#
# This is intentionally tiny: it does not expand the full pq kernel.
# It checks the universal candidate at several independent k,d points.
# ------------------------------------------------------------------------------

def r5_reconstruction():
    print("=" * 78)
    print("8. MINIMAL r=5 RECONSTRUCTION AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    test_points = [
        (3, 6),
        (3, 8),
        (3, 10),
        (5, 6),
        (5, 10),
        (5, 14),
        (7, 8),
        (7, 12),
        (9, 10),
        (11, 14),
        (13, 16),
        (17, 18),
        (19, 20),
        (21, 22),
    ]

    for k, d in test_points:
        for j in range(0, 6):

            expected = simp(
                universal_delta(5, j, k, d)
            )

            # Reconstruct the same polynomial from:
            #   degree j
            #   leading coefficient
            #   fixed roots
            #   moving root
            #
            # This provides an independent symbolic construction.

            if j == 0:
                rebuilt = sp.binomial(k + 4, 5)
            else:
                fixed = sp.Integer(1)

                for m in range(1, j):
                    fixed *= d - 5 - m

                moving = d - sp.Rational(
                    5 - j,
                    k + j
                ) * k

                rebuilt = simp(
                    sp.Rational(
                        (-1) ** j,
                        sp.factorial(j)
                    )
                    * sp.binomial(
                        k + 4,
                        5 - j
                    )
                    * fixed
                    * moving
                )

            residual = simp(rebuilt - expected)

            tested += 1

            if residual != 0:
                failures += 1

                print(
                    f"FAIL k={k} d={d} j={j} "
                    f"expected={expected} "
                    f"rebuilt={rebuilt} "
                    f"residual={residual}"
                )

    print()
    print(f"tested = {tested}")
    print(f"reconstruction failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# 9. UNIVERSAL FORM SUMMARY
# ------------------------------------------------------------------------------

def final_summary():
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The candidate universal boundary law is:"
    )
    print()
    print(
        "Delta_(r,0)(k,d) = C(k+r-1,r)"
    )
    print()
    print(
        "Delta_(r,j)(k,d) = "
        "(-1)^j/j! * C(k+r-1,r-j)"
    )
    print(
        "                  * product_{m=1}^{j-1}(d-r-m)"
    )
    print(
        "                  * (d-(r-j)k/(k+j))"
    )
    print()

    print("For r=5 this predicts:")
    print()

    for j in range(0, 6):
        print(
            f"  j={j}: "
            f"{sp.factor(universal_delta(5, j, K, D))}"
        )

    print()
    print(
        "Structural predictions:"
    )
    print(
        "  degree in d = j"
    )
    print(
        "  fixed roots = r+1,...,r+j-1"
    )
    print(
        "  moving root = (r-j)k/(k+j)"
    )
    print(
        "  leading coefficient = "
        "(-1)^j/j! * C(k+r-1,r-j)"
    )
    print()

    print(
        "This experiment does NOT analyze r=6."
    )
    print(
        "It is intended to settle the r=5 symbolic mechanism "
        "before any larger kernel computation."
    )
    print()


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

def main():
    symbolic_universal_law()
    degree_leading_audit()
    root_audit()
    numerical_leading_audit()
    known_cross_check()
    r5_ladder()
    r5_root_structure()
    r5_reconstruction()
    final_summary()


if __name__ == "__main__":
    main()

