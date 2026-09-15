#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 532 START")
print("=" * 78)
print("ABSTRACT OPERATOR CONJUGACY / INVARIANT-IDEAL TEST")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

U, V = sp.symbols("U V")
c = sp.symbols("c", nonzero=True)
h = sp.symbols("h")
z = sp.symbols("z")

failures = 0


# ============================================================================
# CERTIFICATE
# ============================================================================

def cert(label, expr, expected=0):
    global failures

    difference = sp.factor(
        sp.expand(expr - expected)
    )

    ok = difference == 0

    print(f"  {label}")
    print(f"    difference = {difference}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# OPERATOR
# ============================================================================

def T(expr, step=1):
    """
    Distinguished translation operator:

        U -> U + step*c
        V -> V
    """
    return sp.expand(
        expr.subs(U, U + step*c)
    )


# ============================================================================
# [1] CANONICAL OPERATOR
# ============================================================================

print()
print("[1] CANONICAL OPERATOR")
print("-" * 78)

print("""
  T_h:

      U -> U + h*c
      V -> V
""")

cert(
    "T_h(V)-V",
    T(V, h) - V,
    0,
)

cert(
    "T_h(U)-U",
    T(U, h) - U,
    h*c,
)


# ============================================================================
# [2] POLYNOMIAL INVARIANT ALGEBRA
# ============================================================================

print()
print("[2] POLYNOMIAL INVARIANT ALGEBRA")
print("-" * 78)

coefficients = {}
terms = []

for i in range(4):
    for j in range(4 - i):
        coeff = sp.symbols(f"c_{i}_{j}")
        coefficients[(i, j)] = coeff
        terms.append(
            coeff * U**i * V**j
        )

generic_poly = sp.expand(
    sum(terms)
)

invariant_difference = sp.Poly(
    sp.expand(T(generic_poly) - generic_poly),
    U,
    V,
)

equations = [
    sp.Eq(coef, 0)
    for _, coef in invariant_difference.terms()
]

solutions = sp.solve(
    equations,
    list(coefficients.values()),
    dict=True,
)

print("  generic polynomial:")
print("   ", generic_poly)

print()
print("  invariant solutions:")

for solution in solutions:
    print("   ", solution)

print()
print("""
  Expected invariant algebra:

      Q[V]

  i.e. every polynomial invariant is independent of U.
""")


# ============================================================================
# [3] GENERAL AFFINE OBSERVABLE
# ============================================================================

print()
print("[3] GENERAL AFFINE OBSERVABLE")
print("-" * 78)

aU, aV, a0 = sp.symbols(
    "aU aV a0"
)

H = sp.expand(
    aU*U + aV*V + a0
)

dH = sp.factor(
    T(H, h) - H
)

print("  H =", H)
print()
print("  T_h(H)-H =", dH)

cert(
    "constant translation response",
    dH,
    aU*h*c,
)

cert(
    "no V dependence in response",
    sp.diff(dH, V),
    0,
)


# ============================================================================
# [4] GENERAL QUADRATIC OBSERVABLE
# ============================================================================

print()
print("[4] GENERAL QUADRATIC OBSERVABLE")
print("-" * 78)

k20, k11, k02, k10, k01, k00 = sp.symbols(
    "k20 k11 k02 k10 k01 k00"
)

K = sp.expand(
    k20*U**2
    + k11*U*V
    + k02*V**2
    + k10*U
    + k01*V
    + k00
)

dK = sp.factor(
    T(K, h) - K
)

print("  K =")
print("   ", K)

print()
print("  T_h(K)-K =")
print("   ", dK)


# ============================================================================
# [5] PURE MOVING CHANNEL
# ============================================================================

print()
print("[5] PURE MOVING CHANNEL")
print("-" * 78)

print("""
  A genuine KAPPA quadratic must have no independent V term
  in its translation response.

  Therefore:

      k11 = 0.
""")

cert(
    "V-dependent quadratic cross term condition",
    sp.diff(dK, V).subs({
        k11: 0
    }),
    0,
)


# ============================================================================
# [6] CANONICAL QUADRATIC KERNEL
# ============================================================================

print()
print("[6] CANONICAL QUADRATIC KERNEL")
print("-" * 78)

Kcanon = sp.expand(
    U**2 - V**2
)

dKcanon = sp.factor(
    T(Kcanon, h) - Kcanon
)

print("  K =", Kcanon)
print()
print("  T_h(K)-K =", dKcanon)

cert(
    "canonical translation response",
    dKcanon,
    2*c*h*U + c**2*h**2,
)


# ============================================================================
# [7] MOVING COORDINATE EXTRACTION
# ============================================================================

print()
print("[7] MOVING COORDINATE EXTRACTION")
print("-" * 78)

H_recovered = sp.factor(
    (
        dKcanon
        - c**2*h**2
    ) /
    (2*c*h)
)

print("  H_recovered =", H_recovered)

cert(
    "recover U",
    H_recovered,
    U,
)


# ============================================================================
# [8] GAP-INVARIANT EXTRACTION
# ============================================================================

print()
print("[8] GAP-INVARIANT EXTRACTION")
print("-" * 78)

V2_recovered = sp.expand(
    U**2 - Kcanon
)

print("  V^2_recovered =", V2_recovered)

cert(
    "recover V^2",
    V2_recovered,
    V**2,
)


# ============================================================================
# [9] DISCRETE DIFFERENCE STRUCTURE
# ============================================================================

print()
print("[9] DISCRETE DIFFERENCE STRUCTURE")
print("-" * 78)

D1 = sp.expand(
    T(Kcanon) - Kcanon
)

D2 = sp.expand(
    T(Kcanon, 2)
    - 2*T(Kcanon)
    + Kcanon
)

D3 = sp.expand(
    T(Kcanon, 3)
    - 3*T(Kcanon, 2)
    + 3*T(Kcanon)
    - Kcanon
)

cert(
    "D(K)",
    D1,
    2*c*U + c**2,
)

cert(
    "D^2(K)",
    D2,
    2*c**2,
)

cert(
    "D^3(K)",
    D3,
    0,
)


# ============================================================================
# [10] AFFINE MOVING COORDINATE
# ============================================================================

print()
print("[10] AFFINE MOVING COORDINATE")
print("-" * 78)

B, I = sp.symbols("B I")

H_general = U + B*V + I

dH_general = sp.expand(
    T(H_general) - H_general
)

print("  H =", H_general)
print()
print("  D(H) =", dH_general)

cert(
    "affine gauge response",
    dH_general,
    c,
)


# ============================================================================
# [11] INVARIANT GAUGE FOR K
# ============================================================================

print()
print("[11] INVARIANT GAUGE FOR K")
print("-" * 78)

J0, J1, J2 = sp.symbols(
    "J0 J1 J2"
)

# Arbitrary invariant polynomial in V.
J = sp.expand(
    J2*V**2
    + J1*V
    + J0
)

K_gauge = sp.expand(
    Kcanon + J
)

cert(
    "K invariant-gauge response",
    T(K_gauge) - K_gauge,
    2*c*U + c**2,
)


# ============================================================================
# [12] NORMALIZED KAPPA COORDINATES
# ============================================================================

print()
print("[12] NORMALIZED KAPPA COORDINATES")
print("-" * 78)

N = sp.expand(
    Kcanon / c**2
)

S = sp.expand(
    2*U / c
)

Delta = sp.expand(
    4*V**2 / c**2
)

print("  N =", N)
print("  S =", S)
print("  Delta =", Delta)

cert(
    "N",
    N,
    (U**2 - V**2)/c**2,
)

cert(
    "S",
    S,
    2*U/c,
)

cert(
    "Delta",
    Delta,
    4*V**2/c**2,
)


# ============================================================================
# [13] SPECTRAL POLYNOMIAL
# ============================================================================

print()
print("[13] SPECTRAL POLYNOMIAL")
print("-" * 78)

Q = sp.expand(
    z**2 - S*z + N
)

print("  Q(z) =")
print("   ", sp.factor(Q))

Q_discriminant = sp.factor(
    sp.discriminant(Q, z)
)

print()
print("  discriminant(Q) =", Q_discriminant)

cert(
    "Q discriminant",
    Q_discriminant,
    Delta,
)


# ============================================================================
# [14] ROOT FACTORIZATION
# ============================================================================

print()
print("[14] ROOT FACTORIZATION")
print("-" * 78)

Q_expected = sp.expand(
    (-U + c*z - V)
    * (-U + c*z + V)
    / c**2
)

print("  expected factor form:")
print("   ", sp.factor(Q_expected))

cert(
    "Q factorization",
    Q,
    Q_expected,
)


# ============================================================================
# [15] +1 TRANSLATED FACTOR
# ============================================================================

print()
print("[15] +1 TRANSLATED FACTOR")
print("-" * 78)

Q_shift = sp.expand(
    Q.subs(z, z - 1)
)

Q_shift_expected = sp.expand(
    (-U + c*(z - 1) - V)
    * (-U + c*(z - 1) + V)
    / c**2
)

print("  Q(z-1) =")
print("   ", sp.factor(Q_shift))

cert(
    "Q(z-1)",
    Q_shift,
    Q_shift_expected,
)


# ============================================================================
# [16] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ============================================================================

print()
print("[16] FOUR-BASE CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi = sp.expand(
    Q * Q_shift
)

chi_expected = sp.expand(
    (
        (-U + c*z - V)
        * (-U + c*z + V)
        * (-U + c*(z - 1) - V)
        * (-U + c*(z - 1) + V)
    )
    / c**4
)

print("  chi(z) =")
print("   ", sp.factor(chi))

cert(
    "four-base characteristic polynomial",
    chi,
    chi_expected,
)


# ============================================================================
# [17] ORBIT REPRESENTATION
# ============================================================================

print()
print("[17] TRANSLATION ORBIT")
print("-" * 78)

m = sp.symbols("m")

L_m = sp.expand(
    (U - V + m*c)
    * (U + V + m*c)
    / c**2
)

L_target = sp.expand(
    m**2 + S*m + N
)

print("  L(m) =")
print("   ", sp.factor(L_m))

cert(
    "L(m)=m^2+S*m+N",
    L_m,
    L_target,
)


# ============================================================================
# [18] ORBIT DIFFERENCES
# ============================================================================

print()
print("[18] ORBIT DIFFERENCES")
print("-" * 78)

D_L = sp.expand(
    L_m.subs(m, m + 1)
    - L_m
)

D2_L = sp.expand(
    L_m.subs(m, m + 2)
    - 2*L_m.subs(m, m + 1)
    + L_m
)

print("  ΔL =", sp.factor(D_L))
print()
print("  Δ²L =", sp.factor(D2_L))

cert(
    "first difference",
    D_L,
    S + 2*m + 1,
)

cert(
    "second difference",
    D2_L,
    2,
)


# ============================================================================
# [19] RECOVERY FROM ORBIT
# ============================================================================

print()
print("[19] ORBIT -> KAPPA COORDINATES")
print("-" * 78)

S_from_orbit = sp.factor(
    D_L - (2*m + 1)
)

N_from_orbit = sp.expand(
    L_m - m*S_from_orbit - m**2
)

print("  S_from_orbit =", S_from_orbit)
print("  N_from_orbit =", N_from_orbit)

cert(
    "S from orbit",
    S_from_orbit,
    S,
)

cert(
    "N from orbit",
    N_from_orbit,
    N,
)


# ============================================================================
# [20] OPERATOR CONJUGACY RECOGNITION
# ============================================================================

print()
print("[20] OPERATOR CONJUGACY RECOGNITION")
print("-" * 78)

print("""
  An unknown upstream operator T is a direct candidate when
  there exist observables J, H, K and a nonzero constant c0
  satisfying:

      T(J)-J = 0

      T(H)-H = c0

      T(K)-K = 2*c0*H + c0^2

      T^2(K)-2*T(K)+K = 2*c0^2.

  Gauge freedom is allowed:

      H -> H + I(V)

      K -> K + J(V)

  for T-invariant I,J.

  Canonical realization:

      J = V
      H = U
      K = U^2-V^2.
""")


# ============================================================================
# [21] PURE OPERATOR INVARIANT
# ============================================================================

print()
print("[21] INVARIANT / MOVING DECOMPOSITION")
print("-" * 78)

print("""
  In the canonical representation:

      invariant coordinate:
          V

      moving coordinate:
          U

      quadratic kernel:
          U^2-V^2.

  Consequently:

      V  -> invariant
      U  -> affine
      K  -> quadratic orbit.
""")


# ============================================================================
# [22] FINAL EXACT AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 532 FINISHED")
print("=" * 78)

print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
This experiment only establishes the recognition normal form.

The next step is to substitute an ACTUAL transformation from
the original homogeneous-layer construction.

For that operator T, search for observables satisfying:

    T(J)-J = 0

    T(H)-H = c

    T(K)-K-2*c*H-c^2 = invariant

    T^2(K)-2*T(K)+K = constant.

Do not fit numerical data.

Do not enumerate factors.

Do not introduce p or q.

The desired result is an exact operator-level conjugacy:

    original homogeneous operator
              |
              v
        invariant algebra
              +
        moving coordinate
              +
        quadratic kernel
              |
              v
        KAPPA translation system.
""")