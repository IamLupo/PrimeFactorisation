# ==============================================================================
# EXPERIMENT 247
# EXACT r=4 BOUNDARY LADDER + CROSS-r STRUCTURE AUDIT
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No L5 analysis
#
# IMPORTANT:
#   Symbolic expressions use ordinary SymPy division.
#   Do NOT use sp.Rational(K, ...) for symbolic K.
#
# Stable region:
#   d >= 6
#
# The r=4 stable boundary laws recovered from Experiment 245/246 are tested
# exactly, then compared structurally with the already-established r=2 and
# r=3 boundary ladders.
#
# ==============================================================================

import sympy as sp


# ------------------------------------------------------------------------------
# Symbols
# ------------------------------------------------------------------------------

K, D = sp.symbols("K D")


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

def S(x):
    return sp.factor(sp.cancel(sp.expand(x)))


def exact_int(x):
    return sp.Integer(x)


# ------------------------------------------------------------------------------
# Stable r=4 data from the exact-kernel experiments
# ------------------------------------------------------------------------------
#
# Keyed by k, then j.
# Each item is (d, Delta_4,j).
#
# Only stable d >= 6 values are included here.
# d=4 is deliberately excluded because it is exceptional.
# ------------------------------------------------------------------------------

DATA_R4 = {
    3: {
        0: [(6,15),(8,15),(10,15),(12,15),(14,15),(16,15)],
        1: [(6,-75),(8,-115),(10,-155),(12,-195),(14,-235),(16,-275)],
        2: [(6,36),(8,153),(10,330),(12,567),(14,864),(16,1221)],
        3: [(6,0),(8,-45),(10,-190),(12,-483),(14,-972),(16,-1705)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },

    5: {
        0: [(6,70),(8,70),(10,70),(12,70),(14,70),(16,70)],
        1: [(6,-196),(8,-308),(10,-420),(12,-532),(14,-644),(16,-756)],
        2: [(6,64),(8,276),(10,600),(12,1036),(14,1584),(16,2244)],
        3: [(6,0),(8,-59),(10,-250),(12,-637),(14,-1284),(16,-2255)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },

    7: {
        0: [(6,210),(8,210),(10,210),(12,210),(14,210),(16,210)],
        1: [(6,-405),(8,-645),(10,-885),(12,-1125),(14,-1365),(16,-1605)],
        2: [(6,100),(8,435),(10,950),(12,1645),(14,2520),(16,3575)],
        3: [(6,0),(8,-73),(10,-310),(12,-791),(14,-1596),(16,-2805)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },

    9: {
        0: [(6,495),(8,495),(10,495),(12,495),(14,495),(16,495)],
        1: [(6,-726),(8,-1166),(10,-1606),(12,-2046),(14,-2486),(16,-2926)],
        2: [(6,144),(8,630),(10,1380),(12,2394),(14,3672),(16,5214)],
        3: [(6,0),(8,-87),(10,-370),(12,-945),(14,-1908),(16,-3355)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },

    11: {
        0: [(6,1001),(8,1001),(10,1001),(12,1001),(14,1001),(16,1001)],
        1: [(6,-1183),(8,-1911),(10,-2639),(12,-3367),(14,-4095),(16,-4823)],
        2: [(6,196),(8,861),(10,1890),(12,3283),(14,5040),(16,7161)],
        3: [(6,0),(8,-101),(10,-430),(12,-1099),(14,-2220),(16,-3905)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },

    13: {
        0: [(6,1820),(8,1820),(10,1820),(12,1820),(14,1820),(16,1820)],
        1: [(6,-1800),(8,-2920),(10,-4040),(12,-5160),(14,-6280),(16,-7400)],
        2: [(6,256),(8,1128),(10,2480),(12,4312),(14,6624),(16,9416)],
        3: [(6,0),(8,-115),(10,-490),(12,-1253),(14,-2532),(16,-4455)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },

    15: {
        0: [(6,3060),(8,3060),(10,3060),(12,3060),(14,3060),(16,3060)],
        1: [(6,-2601),(8,-4233),(10,-5865),(12,-7497),(14,-9129),(16,-10761)],
        2: [(6,324),(8,1431),(10,3150),(12,5481),(14,8424),(16,11979)],
        3: [(6,0),(8,-129),(10,-550),(12,-1407),(14,-2844),(16,-5005)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },

    17: {
        0: [(6,4845),(8,4845),(10,4845),(12,4845),(14,4845),(16,4845)],
        1: [(6,-3610),(8,-5890),(10,-8170),(12,-10450),(14,-12730),(16,-15010)],
        2: [(6,400),(8,1770),(10,3900),(12,6790),(14,10440),(16,14850)],
        3: [(6,0),(8,-143),(10,-610),(12,-1561),(14,-3156),(16,-5555)],
        4: [(6,0),(8,2),(10,25),(12,105),(14,294),(16,660)],
    },
}


# ------------------------------------------------------------------------------
# EXACT r=4 STABLE BOUNDARY LAWS
#
# These are now written in a symbolic-safe way.
# ------------------------------------------------------------------------------

def delta4(j, k, d):
    k = sp.sympify(k)
    d = sp.sympify(d)

    if j == 0:
        return sp.binomial(k + 3, 4)

    if j == 1:
        return (
            -sp.binomial(k + 3, 3)
            * (d - 3 * k / (k + 1))
        )

    if j == 2:
        return (
            (k + 3) / 4
            * (d - 5)
            * ((k + 2) * d - 2 * k)
        )

    if j == 3:
        return (
            -(k + 3) / 6
            * (d - 5)
            * (d - 6)
            * (d - k / (k + 3))
        )

    if j == 4:
        return (
            d
            * (d - 5)
            * (d - 6)
            * (d - 7)
            / 24
        )

    raise ValueError("j must be 0,1,2,3,4")


# ------------------------------------------------------------------------------
# Section 1
# Direct validation
# ------------------------------------------------------------------------------

def direct_validation():
    print("=" * 78)
    print("1. DIRECT r=4 STABLE-LAW VALIDATION")
    print("=" * 78)

    tested = 0
    failures = 0

    for k in sorted(DATA_R4):
        for j in range(5):
            for d0, actual0 in DATA_R4[k][j]:
                actual = exact_int(actual0)
                expected = S(delta4(j, exact_int(k), exact_int(d0)))
                residual = S(actual - expected)

                tested += 1

                if residual == 0:
                    status = "PASS"
                else:
                    status = "FAIL"
                    failures += 1

                print(
                    f"k={k:2d} j={j} d={d0:2d} "
                    f"actual={str(actual):>7} "
                    f"expected={str(expected):>7} "
                    f"residual={str(residual):>7} {status}"
                )

    print()
    print(f"tested = {tested}")
    print(f"stable-law failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 2
# Finite-difference degree
# ------------------------------------------------------------------------------

def finite_difference_degree(values):
    current = list(values)

    for order in range(len(current)):
        if len(current) == 1:
            return order

        if all(sp.simplify(current[i] - current[0]) == 0
               for i in range(len(current))):
            return order

        current = [
            sp.simplify(current[i + 1] - current[i])
            for i in range(len(current) - 1)
        ]

    return len(values) - 1


def degree_audit():
    print("=" * 78)
    print("2. STABLE DEGREE LADDER")
    print("=" * 78)

    failures = 0

    for k in sorted(DATA_R4):
        for j in range(5):
            values = [
                exact_int(v)
                for _, v in DATA_R4[k][j]
            ]

            degree = finite_difference_degree(values)

            status = "PASS" if degree == j else "FAIL"

            print(
                f"k={k:2d} j={j} "
                f"degree={degree} expected={j} {status}"
            )

            if degree != j:
                failures += 1

    print()
    print(f"degree failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 3
# Symbolic formulas
# ------------------------------------------------------------------------------

def symbolic_laws():
    print("=" * 78)
    print("3. SYMBOLIC r=4 BOUNDARY LAWS")
    print("=" * 78)

    for j in range(5):
        expression = S(delta4(j, K, D))

        print(f"j={j}")
        print(f"  Delta_4,{j}(K,D) = {expression}")
        print()

    print()


# ------------------------------------------------------------------------------
# Section 4
# Symbolic degree and coefficient audit
# ------------------------------------------------------------------------------

def coefficient_audit():
    print("=" * 78)
    print("4. SYMBOLIC DEGREE / COEFFICIENT AUDIT")
    print("=" * 78)

    failures = 0

    for j in range(5):
        expr = sp.Poly(
            sp.expand(delta4(j, K, D)),
            D
        )

        degree = expr.degree()

        print(f"j={j}")
        print(f"  degree = {degree}")

        if degree != j:
            failures += 1

        for p in range(j, -1, -1):
            coeff = S(expr.coeff_monomial(D**p))
            print(f"  [D^{p}] = {coeff}")

        print()

    print(f"symbolic degree failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 5
# Symbolic roots
# ------------------------------------------------------------------------------

def root_audit():
    print("=" * 78)
    print("5. SYMBOLIC ROOT AUDIT")
    print("=" * 78)

    failures = 0

    expected_roots = {
        0: [],
        1: [3 * K / (K + 1)],
        2: [5, 2 * K / (K + 2)],
        3: [5, 6, K / (K + 3)],
        4: [0, 5, 6, 7],
    }

    for j in range(5):
        expr = S(delta4(j, K, D))

        print(f"j={j}")
        print(f"  Delta = {expr}")

        for root in expected_roots[j]:
            residual = S(expr.subs(D, root))

            print(
                f"  root={root} "
                f"residual={residual} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

            if residual != 0:
                failures += 1

        print()

    print(f"root failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 6
# d = 6 normalization audit
#
# This is useful because d=6 is the first point in the stable region.
# ------------------------------------------------------------------------------

def d6_audit():
    print("=" * 78)
    print("6. d=6 INITIAL-VALUE AUDIT")
    print("=" * 78)

    for k in sorted(DATA_R4):
        print(f"k={k}")

        for j in range(5):
            value = S(delta4(j, exact_int(k), exact_int(6)))
            print(f"  j={j}: {value}")

        print()


# ------------------------------------------------------------------------------
# Section 7
# Cross-k normalization
#
# Test whether the k-dependence itself has a simple binomial structure.
# ------------------------------------------------------------------------------

def cross_k_audit():
    print("=" * 78)
    print("7. CROSS-k NORMALIZATION AUDIT")
    print("=" * 78)

    for j in range(5):
        print(f"j={j}")

        for k in sorted(DATA_R4):
            expr = S(delta4(j, exact_int(k), D))

            degree = sp.Poly(sp.expand(expr), D).degree()

            leading = S(
                sp.Poly(sp.expand(expr), D)
                .coeff_monomial(D**degree)
            )

            print(
                f"  k={k:2d} "
                f"degree={degree} "
                f"leading={leading}"
            )

        print()


# ------------------------------------------------------------------------------
# Section 8
# Candidate general pattern detector
#
# We compare the observed root pattern against:
#
#   roots = 5,6,...,j+4
# plus one k-dependent root for j<4.
#
# The last row j=4 has exactly:
#   D(D-5)(D-6)(D-7)/24.
#
# ------------------------------------------------------------------------------

def pattern_audit():
    print("=" * 78)
    print("8. ROOT-LADDER PATTERN AUDIT")
    print("=" * 78)

    print("Observed stable root pattern:")
    print()
    print("j=0: no D-root")
    print("j=1: D = 3k/(k+1)")
    print("j=2: D = 5, 2k/(k+2)")
    print("j=3: D = 5, 6, k/(k+3)")
    print("j=4: D = 0, 5, 6, 7")
    print()

    print(
        "The fixed roots form the sequence 5,6,7,... while the "
        "remaining root moves with k."
    )
    print()


# ------------------------------------------------------------------------------
# Section 9
# Compare r=4 with known r=3 boundary structure.
#
# r=3 stable laws, written in d = ell-k.
# These are extracted from the established data.
# ------------------------------------------------------------------------------

def delta3(j, k, d):
    k = sp.sympify(k)
    d = sp.sympify(d)

    if j == 0:
        return sp.binomial(k + 2, 3)

    if j == 1:
        # Delta_3,1 = -C(k+2,2)*L + k(k+2)(k+3)/2
        # L = k+d
        L = k + d
        return (
            -sp.binomial(k + 2, 2) * L
            + k * (k + 2) * (k + 3) / 2
        )

    if j == 2:
        # Stable r=3 law recovered in d:
        # (d-4) * ((2k-1)d - 2k) / 2
        #
        # This is the corrected stable form visible from the exact data.
        return (
            (d - 4)
            * ((2 * k - 1) * d - 2 * k)
            / 2
        )

    if j == 3:
        return (
            -(d - 5)
            * (d - 4)
            * (d - k / (k + 3))
            / 6
            * (k + 3)
        )

    raise ValueError("j must be 0..3")


def cross_r_comparison():
    print("=" * 78)
    print("9. r=3 / r=4 BOUNDARY COMPARISON")
    print("=" * 78)

    for j in range(4):
        print(f"j={j}")
        print(f"  r=3: {S(delta3(j, K, D))}")
        print(f"  r=4: {S(delta4(j, K, D))}")
        print()

    print()


# ------------------------------------------------------------------------------
# Section 10
# Compact factor comparison
# ------------------------------------------------------------------------------

def compact_factor_audit():
    print("=" * 78)
    print("10. COMPACT FACTOR AUDIT")
    print("=" * 78)

    formulas = {
        0: sp.binomial(K + 3, 4),

        1: -sp.binomial(K + 3, 3)
           * (D - 3 * K / (K + 1)),

        2: (K + 3) / 4
           * (D - 5)
           * ((K + 2) * D - 2 * K),

        3: -(K + 3) / 6
           * (D - 5)
           * (D - 6)
           * (D - K / (K + 3)),

        4: D * (D - 5) * (D - 6) * (D - 7) / 24,
    }

    for j in range(5):
        residual = S(
            formulas[j] - delta4(j, K, D)
        )

        print(
            f"j={j} residual={residual} "
            f"{'PASS' if residual == 0 else 'FAIL'}"
        )

    print()


# ------------------------------------------------------------------------------
# Section 11
# Stable reconstruction at all available data points.
# ------------------------------------------------------------------------------

def reconstruction_audit():
    print("=" * 78)
    print("11. COMPLETE STABLE RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = 0

    for k in sorted(DATA_R4):
        for j in range(5):
            for d0, actual0 in DATA_R4[k][j]:
                actual = exact_int(actual0)
                predicted = S(
                    delta4(j, exact_int(k), exact_int(d0))
                )
                residual = S(actual - predicted)

                tested += 1

                if residual != 0:
                    failures += 1
                    print(
                        f"FAIL k={k} j={j} d={d0} "
                        f"actual={actual} predicted={predicted} "
                        f"residual={residual}"
                    )

    print()
    print(f"tested = {tested}")
    print(f"reconstruction failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 12
# Final diagnostic
# ------------------------------------------------------------------------------

def final_diagnostic():
    print("=" * 78)
    print("12. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print("The r=4 stable boundary ladder is exactly:")
    print()
    print("Delta_4,0 = C(k+3,4)")
    print()
    print(
        "Delta_4,1 = -C(k+3,3)"
        " * (d - 3k/(k+1))"
    )
    print()
    print(
        "Delta_4,2 = (k+3)/4"
        " * (d-5)"
        " * ((k+2)d - 2k)"
    )
    print()
    print(
        "Delta_4,3 = -(k+3)/6"
        " * (d-5)(d-6)"
        " * (d-k/(k+3))"
    )
    print()
    print(
        "Delta_4,4 = d(d-5)(d-6)(d-7)/24"
    )
    print()

    print("Verified:")
    print("  * 240 stable exact data points.")
    print("  * zero stable-law failures.")
    print("  * degree-j ladder for j=0,...,4.")
    print("  * symbolic root identities.")
    print("  * complete stable reconstruction.")
    print()

    print("Important boundary fact:")
    print("  d=4 remains an exceptional endpoint.")
    print("  It is not part of the stable polynomial ladder.")
    print()

    print("Next mathematical target:")
    print(
        "  derive the GENERAL r,j boundary formula from the"
    )
    print(
        "  exact pq kernel, using the r=2, r=3 and r=4 ladders"
    )
    print(
        "  as constraints."
    )
    print()

    print(
        "Do NOT start L5 merely by interpolation."
    )
    print(
        "First attempt to identify the common combinatorial"
    )
    print(
        "boundary mechanism across r."
    )
    print()


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("EXPERIMENT 247")
    print("EXACT r=4 BOUNDARY LADDER + CROSS-r STRUCTURE AUDIT")
    print("=" * 78)
    print()
    print("Exact arithmetic over QQ.")
    print("Standalone main.py.")
    print("No previous experiment imported.")
    print("No filesystem access.")
    print("No L5 analysis.")
    print()

    direct_validation()
    degree_audit()
    symbolic_laws()
    coefficient_audit()
    root_audit()
    d6_audit()
    cross_k_audit()
    pattern_audit()
    cross_r_comparison()
    compact_factor_audit()
    reconstruction_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()

