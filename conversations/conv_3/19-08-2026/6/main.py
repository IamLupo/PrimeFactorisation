import sympy as sp


# =============================================================================
# EXPERIMENT 381 START
# EXACT KAPPA / LAYER ELIMINATION — CORRECTED KAPPA POLYNOMIAL
# =============================================================================

print("=" * 110)
print("EXPERIMENT 381 START")
print("=" * 110)
print()
print("EXACT KAPPA / LAYER ELIMINATION")
print("Corrected identity:")
print()
print("  (K + N^2 + 2) t^2 - N(N+3)t + N^2 = 0")
print()
print("The experiment verifies:")
print("  1. exact Kappa root")
print("  2. exact layer root")
print("  3. exact resultant compatibility")
print("  4. resultant degree structure")
print("  5. whether a common structural factor occurs across instances")
print()


# =============================================================================
# SYMBOLIC VARIABLE
# =============================================================================

t = sp.symbols("t")
K = sp.symbols("K")
R = sp.symbols("R")


# =============================================================================
# EXACT LAYER POLYNOMIALS
# =============================================================================

def h16(t):
    return (
        -9*t**8
        -36*t**7
        -84*t**6
        -126*t**5
        -126*t**4
        -84*t**3
        -36*t**2
        -9*t
        -1
    )


def h15(t):
    return (
        88*t**9
        +396*t**8
        +1164*t**7
        +2226*t**6
        +2898*t**5
        +2604*t**4
        +1596*t**3
        +639*t**2
        +151*t
        +16
    )


def h14(t):
    return (
        -276*t**10
        -1380*t**9
        -5460*t**8
        -13560*t**7
        -23058*t**6
        -27510*t**5
        -23100*t**4
        -13410*t**3
        -5135*t**2
        -1169*t
        -120
    )


H16 = sp.expand(h16(t))
H15 = sp.expand(h15(t))
H14 = sp.expand(h14(t))


NUM_R = sp.expand(H15**2)
DEN_R = sp.expand(H16 * H14)

LAYER_POLY_GENERIC = sp.expand(
    R * DEN_R - NUM_R
)


# =============================================================================
# TEST DATA
# =============================================================================

TEST_PAIRS = (
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (10007, 10009),
    (50021, 50047),
    (100003, 100019),
    (200003, 200009),
    (300007, 900001),
    (500009, 700001),
    (1000003, 1000033),
    (2000003, 3000017),
)


# =============================================================================
# EXACT KAPPA
# =============================================================================

def kappa_from_x(N, X):
    S = X - 1

    return (
        -S**2
        + (N + 1)*S
        - N**2
        + N
    )


# =============================================================================
# EXACT LAYER INVARIANT
# =============================================================================

def layer_from_x(N, X):
    t_value = sp.Rational(N, X)

    return sp.cancel(
        NUM_R.subs(t, t_value)
        /
        DEN_R.subs(t, t_value)
    )


# =============================================================================
# STORAGE FOR CROSS-INSTANCE ANALYSIS
# =============================================================================

all_pass = True

resultant_degrees = []
resultant_content_values = []
resultant_specialized_values = []

# =============================================================================
# INSTANCE LOOP
# =============================================================================

for instance_id, (p, q) in enumerate(TEST_PAIRS, start=1):

    print()
    print("=" * 110)
    print(f"INSTANCE {instance_id}")
    print("=" * 110)

    N = p * q
    X_true = p + q + 1
    t_true = sp.Rational(N, X_true)

    K_true = sp.Integer(
        kappa_from_x(N, X_true)
    )

    R_true = layer_from_x(N, X_true)

    print()
    print("BASIC DATA")
    print("  p       =", p)
    print("  q       =", q)
    print("  N       =", N)
    print("  X_TRUE  =", X_true)
    print("  t_TRUE  =", t_true)

    print()
    print("TRUE INVARIANTS")
    print("  K_TRUE  =", K_true)
    print("  R_TRUE  =", R_true)

    # =========================================================================
    # CORRECTED KAPPA POLYNOMIAL
    # =========================================================================

    kappa_poly = sp.expand(
        (K + N**2 + 2)*t**2
        - N*(N + 3)*t
        + N**2
    )

    print()
    print("CORRECTED KAPPA POLYNOMIAL")
    print("  degree =", sp.degree(kappa_poly, t))

    # =========================================================================
    # DIRECT KAPPA ROOT CHECK
    # =========================================================================

    kappa_root_check = sp.cancel(
        kappa_poly.subs({
            K: K_true,
            t: t_true
        })
    )

    print()
    print("KAPPA ROOT CHECK")
    print("  polynomial(t_TRUE) =", kappa_root_check)

    if kappa_root_check != 0:
        print("  STATUS = FAIL")
        all_pass = False
    else:
        print("  STATUS = PASS")

    # =========================================================================
    # DIRECT LAYER ROOT CHECK
    # =========================================================================

    layer_poly = sp.expand(
        LAYER_POLY_GENERIC
    )

    layer_root_check = sp.cancel(
        layer_poly.subs({
            R: R_true,
            t: t_true
        })
    )

    print()
    print("LAYER ROOT CHECK")
    print("  polynomial(t_TRUE) =", layer_root_check)

    if layer_root_check != 0:
        print("  STATUS = FAIL")
        all_pass = False
    else:
        print("  STATUS = PASS")

    # =========================================================================
    # RESULTANT
    # =========================================================================

    print()
    print("COMPUTING EXACT RESULTANT")

    resultant = sp.resultant(
        kappa_poly,
        layer_poly,
        t
    )

    resultant = sp.factor(resultant)

    resultant_poly = sp.Poly(
        sp.expand(resultant),
        K,
        R
    )

    total_degree = resultant_poly.total_degree()
    k_degree = resultant_poly.degree(K)
    r_degree = resultant_poly.degree(R)

    resultant_degrees.append(
        (total_degree, k_degree, r_degree)
    )

    print()
    print("RESULTANT STRUCTURE")
    print("  total degree =", total_degree)
    print("  K degree     =", k_degree)
    print("  R degree     =", r_degree)

    # =========================================================================
    # CONTENT
    # =========================================================================

    content = sp.Poly(
        resultant,
        K,
        R
    ).content()

    resultant_content_values.append(content)

    print()
    print("RESULTANT CONTENT")
    print("  content =", content)

    # =========================================================================
    # RESULTANT AT TRUE POINT
    # =========================================================================

    specialized = sp.cancel(
        resultant.subs({
            K: K_true,
            R: R_true
        })
    )

    resultant_specialized_values.append(specialized)

    print()
    print("RESULTANT TRUE-POINT CHECK")

    if specialized == 0:
        print("  resultant(K_TRUE, R_TRUE) = 0")
        print("  STATUS = PASS")
    else:
        print("  resultant(K_TRUE, R_TRUE) != 0")
        print("  VALUE =")
        print("  ", specialized)
        print("  STATUS = FAIL")
        all_pass = False

    # =========================================================================
    # FACTOR STRUCTURE WITHOUT PRINTING HUGE FULL RESULTANT
    # =========================================================================

    print()
    print("RESULTANT FACTOR STRUCTURE")

    factor_list = sp.factor_list(resultant)

    constant_factor = factor_list[0]
    factors = factor_list[1]

    print("  number of nonconstant factors =", len(factors))
    print("  constant factor bit-length    =",
          abs(int(constant_factor)).bit_length()
          if constant_factor != 0 and constant_factor.is_Integer
          else "non-integer/large")

    for factor_id, (factor, multiplicity) in enumerate(
        factors,
        start=1
    ):

        poly_factor = sp.Poly(
            factor,
            K,
            R
        )

        print()
        print(f"  FACTOR {factor_id}")
        print("    multiplicity =", multiplicity)
        print("    total degree =", poly_factor.total_degree())
        print("    K degree     =", poly_factor.degree(K))
        print("    R degree     =", poly_factor.degree(R))

        # Only print factors that are structurally small.
        term_count = len(poly_factor.terms())

        print("    term count   =", term_count)

        if term_count <= 12:
            print("    expression   =")
            print("      ", factor)
        else:
            print("    expression   = ** TOO LARGE **")

    # =========================================================================
    # CHECK WHETHER RESULTANT IS QUADRATIC IN R
    # =========================================================================

    print()
    print("LOW-DEGREE CHECKS")

    if r_degree == 2:
        print("  Resultant is quadratic in R: PASS")
    else:
        print("  Resultant is quadratic in R: FAIL")
        all_pass = False

    if k_degree == 18:
        print("  Resultant K-degree = 18: PASS")
    else:
        print("  Resultant K-degree = 18: NOT OBSERVED")

    if total_degree == 20:
        print("  Resultant total degree = 20: PASS")
    else:
        print("  Resultant total degree = 20: NOT OBSERVED")


# =============================================================================
# CROSS-INSTANCE SUMMARY
# =============================================================================

print()
print("=" * 110)
print("CROSS-INSTANCE SUMMARY")
print("=" * 110)

print()

for index, degrees in enumerate(
    resultant_degrees,
    start=1
):
    total_degree, k_degree, r_degree = degrees

    print(
        f"  INSTANCE {index:2d}: "
        f"total={total_degree:<3d} "
        f"K={k_degree:<3d} "
        f"R={r_degree:<3d}"
    )

print()
print("All degree triples identical:",
      len(set(resultant_degrees)) == 1)

print()
print("All resultant true-point checks passed:",
      all(value == 0 for value in resultant_specialized_values))

print()
print("=" * 110)
print("EXPERIMENT 381 FINAL STATUS")
print("=" * 110)

print()
print("  ALL EXACT CHECKS PASS:", all_pass)

print()
print("Interpretation:")
print()
print("  If PASS:")
print("    Kappa and the layer invariant are exactly compatible")
print("    through the common hidden variable t = N/X.")
print()
print("    The next question is whether the resultant itself")
print("    has a useful normalization or universal closed form")
print("    in N, K, and R.")
print()
print("  If FAIL:")
print("    inspect the direct root checks first; this indicates")
print("    an algebraic transcription problem before interpreting")
print("    the resultant.")

print()
print("=" * 110)
print("EXPERIMENT 381 FINISHED")
print("=" * 110)
