# ==============================================================================
# EXPERIMENT 246
# EXACT r=4 STABLE BOUNDARY LADDER CLOSURE
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No L5 analysis
#
# Purpose:
#   Close the r=4 boundary ladder in the stable region d >= 6.
#
# The exact stable data are the values obtained from the exact pq-kernel
# extraction in Experiment 245. They are embedded here explicitly so that
# this script has no dependency on any previous file or module.
#
# IMPORTANT:
#   d=4 is treated as an exceptional endpoint.
#   The stable d>=6 laws are NOT extrapolated through d=4.
#
# Candidate laws to test:
#
#   Delta_4,0 =
#       C(k+3,4)
#
#   Delta_4,1 =
#       -C(k+3,3) * (d - 3k/(k+1))
#
#   Delta_4,2 =
#       (k+3)/4 * (d-5) * ((k+2)d - 2k)
#
#   Delta_4,3 =
#       -(k+3)/6 * (d-5)(d-6) * (d-k/(k+3))
#
#   Delta_4,4 =
#       d(d-5)(d-6)(d-7)/24
#
# ==============================================================================

import sympy as sp


# ------------------------------------------------------------------------------
# Symbols
# ------------------------------------------------------------------------------

K, D = sp.symbols("K D")


# ------------------------------------------------------------------------------
# Exact arithmetic helpers
# ------------------------------------------------------------------------------

def Z(x):
    """Exact SymPy integer."""
    return sp.Integer(x)


def R(p, q=1):
    """Exact rational."""
    return sp.Rational(p, q)


def simp(x):
    return sp.factor(sp.cancel(sp.expand(x)))


# ------------------------------------------------------------------------------
# Stable exact-kernel data from Experiment 245
#
# Stored as:
#   stable_data[k][j] = [(d, Delta_4,j), ...]
#
# All entries are integers and therefore exact QQ values.
# ------------------------------------------------------------------------------

stable_data = {
    3: {
        0: [(6, 15), (8, 15), (10, 15), (12, 15), (14, 15), (16, 15)],
        1: [(6, -75), (8, -115), (10, -155), (12, -195), (14, -235), (16, -275)],
        2: [(6, 36), (8, 153), (10, 330), (12, 567), (14, 864), (16, 1221)],
        3: [(6, 0), (8, -45), (10, -190), (12, -483), (14, -972), (16, -1705)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },

    5: {
        0: [(6, 70), (8, 70), (10, 70), (12, 70), (14, 70), (16, 70)],
        1: [(6, -196), (8, -308), (10, -420), (12, -532), (14, -644), (16, -756)],
        2: [(6, 64), (8, 276), (10, 600), (12, 1036), (14, 1584), (16, 2244)],
        3: [(6, 0), (8, -59), (10, -250), (12, -637), (14, -1284), (16, -2255)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },

    7: {
        0: [(6, 210), (8, 210), (10, 210), (12, 210), (14, 210), (16, 210)],
        1: [(6, -405), (8, -645), (10, -885), (12, -1125), (14, -1365), (16, -1605)],
        2: [(6, 100), (8, 435), (10, 950), (12, 1645), (14, 2520), (16, 3575)],
        3: [(6, 0), (8, -73), (10, -310), (12, -791), (14, -1596), (16, -2805)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },

    9: {
        0: [(6, 495), (8, 495), (10, 495), (12, 495), (14, 495), (16, 495)],
        1: [(6, -726), (8, -1166), (10, -1606), (12, -2046), (14, -2486), (16, -2926)],
        2: [(6, 144), (8, 630), (10, 1380), (12, 2394), (14, 3672), (16, 5214)],
        3: [(6, 0), (8, -87), (10, -370), (12, -945), (14, -1908), (16, -3355)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },

    11: {
        0: [(6, 1001), (8, 1001), (10, 1001), (12, 1001), (14, 1001), (16, 1001)],
        1: [(6, -1183), (8, -1911), (10, -2639), (12, -3367), (14, -4095), (16, -4823)],
        2: [(6, 196), (8, 861), (10, 1890), (12, 3283), (14, 5040), (16, 7161)],
        3: [(6, 0), (8, -101), (10, -430), (12, -1099), (14, -2220), (16, -3905)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },

    13: {
        0: [(6, 1820), (8, 1820), (10, 1820), (12, 1820), (14, 1820), (16, 1820)],
        1: [(6, -1800), (8, -2920), (10, -4040), (12, -5160), (14, -6280), (16, -7400)],
        2: [(6, 256), (8, 1128), (10, 2480), (12, 4312), (14, 6624), (16, 9416)],
        3: [(6, 0), (8, -115), (10, -490), (12, -1253), (14, -2532), (16, -4455)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },

    15: {
        0: [(6, 3060), (8, 3060), (10, 3060), (12, 3060), (14, 3060), (16, 3060)],
        1: [(6, -2601), (8, -4233), (10, -5865), (12, -7497), (14, -9129), (16, -10761)],
        2: [(6, 324), (8, 1431), (10, 3150), (12, 5481), (14, 8424), (16, 11979)],
        3: [(6, 0), (8, -129), (10, -550), (12, -1407), (14, -2844), (16, -5005)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },

    17: {
        0: [(6, 4845), (8, 4845), (10, 4845), (12, 4845), (14, 4845), (16, 4845)],
        1: [(6, -3610), (8, -5890), (10, -8170), (12, -10450), (14, -12730), (16, -15010)],
        2: [(6, 400), (8, 1770), (10, 3900), (12, 6790), (14, 10440), (16, 14850)],
        3: [(6, 0), (8, -143), (10, -610), (12, -1561), (14, -3156), (16, -5555)],
        4: [(6, 0), (8, 2), (10, 25), (12, 105), (14, 294), (16, 660)],
    },
}


# ------------------------------------------------------------------------------
# Candidate stable boundary laws
# ------------------------------------------------------------------------------

def delta4_candidate(j, k, d):
    k = sp.sympify(k)
    d = sp.sympify(d)

    if j == 0:
        return sp.binomial(k + 3, 4)

    if j == 1:
        return (
            -sp.binomial(k + 3, 3)
            * (d - sp.Rational(3) * k / (k + 1))
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
            * (d - sp.Rational(k, k + 3))
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


def delta4_candidate_symbolic(j):
    return simp(delta4_candidate(j, K, D))


# ------------------------------------------------------------------------------
# Section 1
# Direct stable-law validation
# ------------------------------------------------------------------------------

def direct_validation():
    print("=" * 78)
    print("1. DIRECT STABLE r=4 LAW VALIDATION")
    print("=" * 78)

    failures = 0
    tested = 0

    for k in sorted(stable_data):
        for j in range(5):
            for d_raw, actual_raw in stable_data[k][j]:
                d = Z(d_raw)
                actual = Z(actual_raw)
                expected = simp(delta4_candidate(j, Z(k), d))
                residual = simp(actual - expected)

                tested += 1

                if residual != 0:
                    failures += 1
                    status = "FAIL"
                else:
                    status = "PASS"

                print(
                    f"k={k:2d} j={j} d={d_raw:2d} "
                    f"actual={str(actual):>8} "
                    f"expected={str(expected):>8} "
                    f"residual={str(residual):>8} {status}"
                )

    print()
    print(f"tested = {tested}")
    print(f"stable-law failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Exact finite differences with respect to d.
#
# d steps by 2 in the data, but the polynomial degree is independent of
# the step size. We compute symbolic finite differences using an arbitrary
# polynomial as a second verification.
# ------------------------------------------------------------------------------

def finite_difference(values):
    rows = [list(values)]
    current = list(values)

    while len(current) > 1:
        current = [
            sp.simplify(current[i + 1] - current[i])
            for i in range(len(current) - 1)
        ]
        rows.append(current)

    return rows


def degree_from_differences(values):
    rows = finite_difference(values)

    for order, row in enumerate(rows):
        if len(row) == 0:
            return None

        if all(sp.simplify(v - row[0]) == 0 for v in row):
            return order

    return len(values) - 1


# ------------------------------------------------------------------------------
# Section 2
# Degree ladder
# ------------------------------------------------------------------------------

def degree_audit():
    print("=" * 78)
    print("2. STABLE FINITE-DIFFERENCE DEGREE LADDER")
    print("=" * 78)

    failures = 0

    for k in sorted(stable_data):
        for j in range(5):
            values = [
                Z(value)
                for _, value in stable_data[k][j]
            ]

            degree = degree_from_differences(values)

            print(
                f"k={k:2d} j={j} "
                f"degree={degree} expected={j} "
                f"{'PASS' if degree == j else 'FAIL'}"
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
    print("3. SYMBOLIC STABLE BOUNDARY LAWS")
    print("=" * 78)

    for j in range(5):
        expr = delta4_candidate_symbolic(j)
        print(f"j={j}")
        print(f"  Delta_4,{j}(K,D) = {expr}")
        print()

    print()


# ------------------------------------------------------------------------------
# Section 4
# Cross-k coefficient structure.
#
# Expand each polynomial in D and inspect the coefficients.
# ------------------------------------------------------------------------------

def coefficient_structure():
    print("=" * 78)
    print("4. CROSS-k COEFFICIENT STRUCTURE")
    print("=" * 78)

    for j in range(5):
        expr = sp.Poly(sp.expand(delta4_candidate_symbolic(j)), D)

        print(f"j={j}")
        for power in range(j, -1, -1):
            coeff = simp(expr.coeff_monomial(D ** power))
            print(f"  D^{power}: {coeff}")
        print()


# ------------------------------------------------------------------------------
# Section 5
# Root structure
# ------------------------------------------------------------------------------

def root_structure():
    print("=" * 78)
    print("5. SYMBOLIC ROOT STRUCTURE")
    print("=" * 78)

    for j in range(5):
        expr = delta4_candidate_symbolic(j)

        if expr == 0:
            roots = []
        else:
            roots = sp.solve(sp.Eq(expr, 0), D)

        print(f"j={j}")
        print(f"  Delta = {expr}")
        print(f"  roots = {roots}")
        print()

    print()


# ------------------------------------------------------------------------------
# Section 6
# Verify each proposed root symbolically.
# ------------------------------------------------------------------------------

def root_identity_audit():
    print("=" * 78)
    print("6. ROOT IDENTITY AUDIT")
    print("=" * 78)

    failures = 0

    tests = [
        (1, sp.Rational(3) * K / (K + 1)),
        (2, sp.Integer(5)),
        (2, 2 * K / (K + 2)),
        (3, sp.Integer(5)),
        (3, sp.Integer(6)),
        (3, K / (K + 3)),
        (4, sp.Integer(0)),
        (4, sp.Integer(5)),
        (4, sp.Integer(6)),
        (4, sp.Integer(7)),
    ]

    for j, root in tests:
        residual = simp(delta4_candidate_symbolic(j).subs(D, root))

        print(
            f"j={j} root={root} "
            f"residual={residual} "
            f"{'PASS' if residual == 0 else 'FAIL'}"
        )

        if residual != 0:
            failures += 1

    print()
    print(f"root identity failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 7
# Difference between exact stable laws and d=4 endpoint values.
#
# This is deliberate: d=4 was already shown to be exceptional.
# ------------------------------------------------------------------------------

endpoint_data = {
    # (k, j): actual Delta at d=4
    (3, 0): 50,
    (3, 1): 0,
    (3, 2): 0,
    (3, 3): 0,
    (3, 4): -1,

    (5, 0): 196,
    (5, 1): 0,
    (5, 2): 0,
    (5, 3): 0,
    (5, 4): -1,

    (7, 0): 540,
    (7, 1): 0,
    (7, 2): 0,
    (7, 3): 0,
    (7, 4): -1,

    (9, 0): 1210,
    (9, 1): 0,
    (9, 2): 0,
    (9, 3): 0,
    (9, 4): -1,

    (11, 0): 2366,
    (11, 1): 0,
    (11, 2): 0,
    (11, 3): 0,
    (11, 4): -1,

    (13, 0): 4200,
    (13, 1): 0,
    (13, 2): 0,
    (13, 3): 0,
    (13, 4): -1,

    (15, 0): 6936,
    (15, 1): 0,
    (15, 2): 0,
    (15, 3): 0,
    (15, 4): -1,

    (17, 0): 10830,
    (17, 1): 0,
    (17, 2): 0,
    (17, 3): 0,
    (17, 4): -1,
}


def endpoint_audit():
    print("=" * 78)
    print("7. d=4 ENDPOINT SEPARATION AUDIT")
    print("=" * 78)

    for k in sorted({key[0] for key in endpoint_data}):
        print(f"k={k}")

        for j in range(5):
            actual = Z(endpoint_data[(k, j)])
            stable_extension = simp(
                delta4_candidate(j, Z(k), Z(4))
            )
            discrepancy = simp(actual - stable_extension)

            print(
                f"  j={j} "
                f"actual={actual} "
                f"stable-extension={stable_extension} "
                f"endpoint-discrepancy={discrepancy}"
            )

        print()

    print(
        "Interpretation: d=4 is not required to obey the stable d>=6 "
        "polynomial ladder."
    )
    print()


# ------------------------------------------------------------------------------
# Section 8
# Compare formulas against extra d values symbolically.
#
# These are not new kernel evaluations; they simply show exact predicted
# values of the closed formulas beyond the observed data range.
# ------------------------------------------------------------------------------

def symbolic_prediction_table():
    print("=" * 78)
    print("8. SYMBOLIC PREDICTION TABLE FOR FUTURE KERNEL TESTS")
    print("=" * 78)

    future_ds = [18, 20, 22]

    for k in [3, 5, 7, 9, 11, 13, 15, 17]:
        print(f"k={k}")
        for d in future_ds:
            vals = [
                simp(delta4_candidate(j, Z(k), Z(d)))
                for j in range(5)
            ]
            print(
                f"  d={d}: "
                + ", ".join(
                    f"j={j}:{vals[j]}"
                    for j in range(5)
                )
            )
        print()


# ------------------------------------------------------------------------------
# Section 9
# Test whether all formulas have degree exactly j symbolically.
# ------------------------------------------------------------------------------

def symbolic_degree_audit():
    print("=" * 78)
    print("9. SYMBOLIC DEGREE AUDIT")
    print("=" * 78)

    failures = 0

    for j in range(5):
        expr = sp.Poly(
            sp.expand(delta4_candidate_symbolic(j)),
            D,
        )
        degree = expr.degree()

        print(
            f"j={j} symbolic-degree={degree} "
            f"expected={j} "
            f"{'PASS' if degree == j else 'FAIL'}"
        )

        if degree != j:
            failures += 1

    print()
    print(f"symbolic degree failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 10
# Search for a compact common form.
#
# We do not assert this is the final general r formula. This section only
# rewrites the five established r=4 boundary rows in a common style.
# ------------------------------------------------------------------------------

def compact_form_audit():
    print("=" * 78)
    print("10. COMPACT FORM AUDIT")
    print("=" * 78)

    common_forms = {
        0:
            sp.binomial(K + 3, 4),

        1:
            -sp.binomial(K + 3, 3)
            * (D - 3 * K / (K + 1)),

        2:
            (K + 3) / 4
            * (D - 5)
            * ((K + 2) * D - 2 * K),

        3:
            -(K + 3) / 6
            * (D - 5)
            * (D - 6)
            * (D - K / (K + 3)),

        4:
            D * (D - 5) * (D - 6) * (D - 7) / 24,
    }

    failures = 0

    for j in range(5):
        residual = simp(
            common_forms[j] - delta4_candidate_symbolic(j)
        )

        print(
            f"j={j} residual={residual} "
            f"{'PASS' if residual == 0 else 'FAIL'}"
        )

        if residual != 0:
            failures += 1

    print()
    print(f"compact-form failures = {failures}")
    print()


# ------------------------------------------------------------------------------
# Section 11
# Final boundary interpretation.
# ------------------------------------------------------------------------------

def final_diagnostic():
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The stable r=4 boundary ladder is consistent with:"
    )
    print()
    print(
        "Delta_4,0 = C(k+3,4)"
    )
    print(
        "Delta_4,1 = -C(k+3,3) * (d - 3k/(k+1))"
    )
    print(
        "Delta_4,2 = (k+3)/4 * (d-5) * ((k+2)d - 2k)"
    )
    print(
        "Delta_4,3 = -(k+3)/6 * (d-5)(d-6) * (d-k/(k+3))"
    )
    print(
        "Delta_4,4 = d(d-5)(d-6)(d-7)/24"
    )
    print()
    print(
        "The exact kernel data support degree-j behavior in d for"
    )
    print(
        "j=0,1,2,3,4 throughout the stable region d>=6."
    )
    print()
    print(
        "The d=4 endpoint is exceptional and must remain separate."
    )
    print()
    print(
        "Research target after this experiment:"
    )
    print(
        "  compare these five boundary factors against the r=2 and"
    )
    print(
        "  r=3 ladders and search for a single general Delta_r(k,k+j,d)"
    )
    print(
        "  formula."
    )
    print()
    print(
        "Do not start L5 from this experiment alone."
    )
    print(
        "First establish the general boundary mechanism."
    )
    print()


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("EXPERIMENT 246")
    print("EXACT r=4 STABLE BOUNDARY LADDER CLOSURE")
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
    coefficient_structure()
    root_structure()
    root_identity_audit()
    endpoint_audit()
    symbolic_prediction_table()
    symbolic_degree_audit()
    compact_form_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()

