import sympy as sp


# ============================================================
# EXPERIMENT
#
# EXACT KAPPA / LAYER COMPATIBILITY RELATION
#
# Hidden parameter:
#
#       t = N / X
#
# Kappa:
#
#       K = -X^2 + (N+3)X - N^2 - 1
#
# substituting X=N/t gives:
#
#       (K + N^2 + 1)t^2
#       - N(N+3)t
#       + N^2 = 0
#
# Layer invariant:
#
#       R = R(t)
#
# Goal:
#
#   Eliminate t exactly and determine whether K and R
#   satisfy a surprisingly low-complexity algebraic relation.
#
# ============================================================


t = sp.symbols("t")
K = sp.symbols("K")
R = sp.symbols("R")


# ============================================================
# EXACT LAYER POLYNOMIALS
# ============================================================

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


# ============================================================
# R(t) = NUM_R / DEN_R
# ============================================================

NUM_R = sp.expand(H15**2)
DEN_R = sp.expand(H16 * H14)

print("=" * 110)
print("EXPERIMENT: EXACT KAPPA/LAYER ELIMINATION")
print("=" * 110)

print()
print("1. BASIC SYMBOLIC OBJECTS")
print("=" * 110)

print()
print("NUM_R degree =", sp.degree(NUM_R, t))
print("DEN_R degree =", sp.degree(DEN_R, t))


# ============================================================
# TEST INSTANCES
# ============================================================

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


# ============================================================
# EXACT KAPPA
# ============================================================

def kappa_value(N, X):
    S = X - 1

    return (
        -S*S
        + (N+1)*S
        - N*N
        + N
    )


# ============================================================
# EXACT R
# ============================================================

def layer_value(N, X):
    tx = sp.Rational(N, X)

    return sp.cancel(
        NUM_R.subs(t, tx)
        /
        DEN_R.subs(t, tx)
    )


# ============================================================
# MAIN INSTANCE LOOP
# ============================================================

all_pass = True

for idx, (p, q) in enumerate(TEST_PAIRS, start=1):

    N = p*q
    X_true = p + q + 1

    K_true = sp.Integer(
        kappa_value(N, X_true)
    )

    R_true = sp.cancel(
        layer_value(N, X_true)
    )

    print()
    print("=" * 110)
    print(f"INSTANCE {idx}")
    print("=" * 110)

    print()
    print("N      =", N)
    print("X_TRUE =", X_true)
    print("K_TRUE =", K_true)
    print("R_TRUE =", R_true)

    # --------------------------------------------------------
    # Kappa equation in t
    #
    # (K + N^2 + 1)t^2 - N(N+3)t + N^2 = 0
    # --------------------------------------------------------

    kappa_poly = sp.expand(
        (K + N*N + 1)*t**2
        - N*(N+3)*t
        + N*N
    )

    print()
    print("KAPPA POLYNOMIAL IN t:")
    print(kappa_poly)

    # --------------------------------------------------------
    # Layer equation
    #
    # R * DEN_R - NUM_R = 0
    # --------------------------------------------------------

    layer_poly = sp.expand(
        R * DEN_R
        - NUM_R
    )

    print()
    print("LAYER POLYNOMIAL DEGREE =", sp.degree(layer_poly, t))

    # --------------------------------------------------------
    # Eliminate t.
    #
    # Because kappa_poly is degree 2 and layer_poly degree 18,
    # this resultant is much cheaper than eliminating N,K,R
    # symbolically all at once.
    # --------------------------------------------------------

    print()
    print("COMPUTING EXACT RESULTANT...")

    resultant = sp.resultant(
        kappa_poly,
        layer_poly,
        t
    )

    resultant = sp.factor(resultant)

    print()
    print("RESULTANT FACTORIZATION:")
    print("** TO BIG RESULT **")

    # --------------------------------------------------------
    # Remove obvious powers/constants.
    # --------------------------------------------------------

    resultant_poly = sp.Poly(
        sp.expand(resultant),
        K,
        R
    )

    print()
    print("RESULTANT TOTAL DEGREE =",
          resultant_poly.total_degree())

    print("RESULTANT K-DEGREE =",
          resultant_poly.degree(K))

    print("RESULTANT R-DEGREE =",
          resultant_poly.degree(R))

    # --------------------------------------------------------
    # Verify TRUE (K,R) lies exactly on relation.
    # --------------------------------------------------------

    compatibility = sp.cancel(
        resultant.subs({
            K: K_true,
            R: R_true
        })
    )

    print()
    print("TRUE POINT COMPATIBILITY:")
    print("  resultant(K_TRUE, R_TRUE) =", compatibility)

    if compatibility != 0:
        all_pass = False

    # --------------------------------------------------------
    # Check whether relation is linear in R or K.
    # --------------------------------------------------------

    degree_R = resultant_poly.degree(R)
    degree_K = resultant_poly.degree(K)

    print()
    print("LOW-COMPLEXITY STRUCTURE CHECK")

    if degree_R == 1:
        print("  RESULTANT IS LINEAR IN R")
    else:
        print("  Resultant is NOT linear in R")

    if degree_K == 1:
        print("  RESULTANT IS LINEAR IN K")
    else:
        print("  Resultant is NOT linear in K")

    if degree_R <= 2:
        print("  R-degree <= 2")
    else:
        print("  R-degree > 2")

    if degree_K <= 2:
        print("  K-degree <= 2")
    else:
        print("  K-degree > 2")


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("=" * 110)
print("GLOBAL RESULT")
print("=" * 110)

print()
print("All exact compatibility checks passed:", all_pass)

print()
print("INTERPRETATION")
print("=" * 110)

print()
print("The experiment asks whether Kappa and the layer invariant")
print("are algebraically coupled through the same hidden variable")
print("t = N/X in a low-complexity way.")

print()
print("If the resultant is surprisingly low-degree or factors into")
print("simple pieces, that may provide a new route to eliminating X.")

print()
print("If the resultant is extremely complicated, that is also")
print("informative: it indicates that Kappa and R are independent")
print("coordinates of X rather than a simple closed-form transform.")

print()
print("Most important:")
print("  K and R are NOT being treated as independent measurements.")
print("  Both are generated from the same hidden t.")
print()

print("=" * 110)
print("END")
print("=" * 110)
