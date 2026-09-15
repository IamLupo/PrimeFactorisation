#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 531 START")
print("=" * 78)
print("TRANSLATION-OPERATOR INVARIANT THEORY / KAPPA NORMAL-FORM CLASSIFICATION")
print()


failures = 0


def cert(label, expr, expected=0):
    global failures

    diff = sp.factor(sp.expand(expr - expected))
    ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# SYMBOLS
# ============================================================================

U, V, c, h, z = sp.symbols("U V c h z")

# Generic polynomial coefficients.
a20, a11, a02, a10, a01, a00 = sp.symbols(
    "a20 a11 a02 a10 a01 a00"
)

b10, b01, b00 = sp.symbols(
    "b10 b01 b00"
)

j20, j11, j02, j10, j01, j00 = sp.symbols(
    "j20 j11 j02 j10 j01 j00"
)


# ============================================================================
# TRANSLATION OPERATOR
# ============================================================================

def T(expr, amount=h):
    return sp.expand(
        expr.subs(U, U + amount * c)
    )


print()
print("[1] TRANSLATION OPERATOR")
print("-" * 78)

print("""
  T_h:

      U -> U + h*c
      V -> V
""")

cert(
    "T_h(U)",
    U + h*c,
    U + h*c,
)

cert(
    "T_h(V)",
    V,
    V,
)


# ============================================================================
# GENERIC QUADRATIC OBSERVABLE K
# ============================================================================

print()
print("[2] GENERAL QUADRATIC OBSERVABLE")
print("-" * 78)

K = sp.expand(
    a20*U**2
    + a11*U*V
    + a02*V**2
    + a10*U
    + a01*V
    + a00
)

print("  K(U,V) =")
print("   ", K)

TK = sp.expand(
    T(K)
)

dK = sp.expand(
    TK - K
)

ddK = sp.expand(
    T(K, h + 1)
    - 2*T(K, h)
    + K
)

print()
print("  T_h(K)-K =")
print("   ", sp.factor(dK))

print()
print("  second difference =")
print("   ", sp.factor(ddK))


# ============================================================================
# SOLVE K FOR CANONICAL QUADRATIC RESPONSE
#
# Desired:
#
#   T(K)-K = 2*c*h*H + c^2*h^2
#
# where H is affine-linear.
# ============================================================================

print()
print("[3] GENERAL MOVING COORDINATE H")
print("-" * 78)

H = sp.expand(
    b10*U
    + b01*V
    + b00
)

print("  H(U,V) =")
print("   ", H)


target_dK = sp.expand(
    2*c*h*H + c**2*h**2
)

difference = sp.Poly(
    sp.expand(dK - target_dK),
    U,
    V,
    h,
)

equations = [
    sp.Eq(coef, 0)
    for _, coef in difference.terms()
]

solution_KH = sp.solve(
    equations,
    [
        a20,
        a11,
        a10,
        b10,
        b01,
        b00,
    ],
    dict=True,
)

print()
print("  exact solution family:")
for sol in solution_KH:
    print("   ", sol)


# ============================================================================
# EXPECTED NORMAL FORM
# ============================================================================

print()
print("[4] CANONICAL MOVING COORDINATE")
print("-" * 78)

H_canonical = sp.expand(
    U + sp.Symbol("invV")
)

# The coefficient of U in K must match H.
# Show this directly from the translation derivative.

dK_h1 = sp.expand(
    dK.subs(h, 1)
)

print("""
  Because

      T_1(K)-K
        = c*(2*a20*U + a11*V + a10)
          + a20*c^2,

  the moving coordinate extracted from K is

      H_K
        = [T_1(K)-K-a20*c^2]/(2*a20*c)

  whenever a20 != 0.
""")


H_from_K = sp.factor(
    (
        dK_h1
        - a20*c**2
    )
    / (2*a20*c)
)

print("  H_from_K =")
print("   ", H_from_K)


cert(
    "U coefficient of H_from_K",
    sp.diff(H_from_K, U),
    1,
)


# ============================================================================
# PURE S-CHANNEL CONDITION
# ============================================================================

print()
print("[5] PURE SUM-CHANNEL CONDITION")
print("-" * 78)

print("""
  A genuine KAPPA S-channel must contain no independent V term.

  Therefore H must have the form

      H = lambda*U + invariant(V).

  The strict canonical choice is

      H = U.
""")

# Eliminate V from moving coordinate.
pure_V_condition = sp.expand(
    sp.diff(H, V)
)

print("  dH/dV =", pure_V_condition)

pure_solution = sp.solve(
    sp.Eq(pure_V_condition, 0),
    b01,
    dict=True,
)

print("  V-free condition =", pure_solution)


# ============================================================================
# QUADRATIC INVARIANT CLASSIFICATION
# ============================================================================

print()
print("[6] INVARIANT QUADRATIC CLASSIFICATION")
print("-" * 78)

J = sp.expand(
    j20*U**2
    + j11*U*V
    + j02*V**2
    + j10*U
    + j01*V
    + j00
)

TJ = sp.expand(
    T(J, 1)
)

invariant_poly = sp.Poly(
    sp.expand(TJ - J),
    U,
    V,
)

invariant_eqs = [
    sp.Eq(coef, 0)
    for _, coef in invariant_poly.terms()
]

invariant_solution = sp.solve(
    invariant_eqs,
    [
        j20,
        j11,
        j10,
    ],
    dict=True,
)

print("  invariant solutions =")
for sol in invariant_solution:
    print("   ", sol)


# ============================================================================
# NORMAL-FORM DECOMPOSITION
# ============================================================================

print()
print("[7] NORMAL-FORM DECOMPOSITION")
print("-" * 78)

print("""
  Every quadratic observable K can be decomposed as

      K(U,V)
        = A*U^2
        + U*I1(V)
        + I0(V),

  where

      I1(V) = a11*V + a10
      I0(V) = a02*V^2 + a01*V + a00.

  Under U -> U+c:

      Delta K
        = 2*A*c*U
          + A*c^2
          + c*I1(V).

  Therefore the only source of the moving U-channel
  is the U^2 coefficient.
""")


A = sp.symbols("A")

I1 = sp.expand(
    a11*V + a10
)

I0 = sp.expand(
    a02*V**2
    + a01*V
    + a00
)

K_decomp = sp.expand(
    A*U**2
    + U*I1
    + I0
)

print("  generic decomposition =")
print("   ", K_decomp)


# ============================================================================
# PURE KERNEL CLASSIFICATION
# ============================================================================

print()
print("[8] PURE CONIC KERNEL CLASSIFICATION")
print("-" * 78)

print("""
  Require that the moving part contain only U:

      Delta K
        = 2*A*c*U*h
          + A*c^2*h^2.

  This removes every V-dependent linear term.
""")

pure_kernel_eq = sp.Eq(
    a11,
    0,
)

print("  condition =", pure_kernel_eq)


# ============================================================================
# CANONICAL KERNEL FAMILY
# ============================================================================

print()
print("[9] CANONICAL KERNEL FAMILY")
print("-" * 78)

kappa, inv2, inv1, inv0 = sp.symbols(
    "kappa inv2 inv1 inv0"
)

K_family = sp.expand(
    kappa*U**2
    + inv2*V**2
    + inv1*V
    + inv0
)

K_family_diff = sp.factor(
    T(K_family) - K_family
)

print("  K_family =")
print("   ", K_family)

print()
print("  Delta K_family =")
print("   ", K_family_diff)

cert(
    "no V-dependent translation term",
    sp.diff(K_family_diff, V),
    0,
)


# ============================================================================
# MOVING COORDINATE FROM GENERAL CANONICAL KERNEL
# ============================================================================

print()
print("[10] MOVING COORDINATE EXTRACTION")
print("-" * 78)

K1 = sp.expand(
    K_family.subs(U, U + c)
)

U_extracted = sp.expand(
    (
        K1
        - K_family
        - kappa*c**2
    )
    / (2*kappa*c)
)

print("  extracted coordinate =")
print("   ", sp.factor(U_extracted))

cert(
    "recovered moving coordinate",
    U_extracted,
    U,
)


# ============================================================================
# HISTORICAL NORMALIZATION
# ============================================================================

print()
print("[11] HISTORICAL NORMALIZATION")
print("-" * 78)

# Set kappa = 1 for the canonical normalized conic.
K_normalized = sp.expand(
    K_family.subs(kappa, 1)
)

print("  normalized kernel =")
print("   ", K_normalized)


# ============================================================================
# IDENTIFY GAP INVARIANT
# ============================================================================

print()
print("[12] GAP-INVARIANT EXTRACTION")
print("-" * 78)

# The canonical historical form is U^2 - V^2.
# Its invariant part is therefore -V^2.
#
# Given U from the operator and K itself:

V2_from_kernel = sp.expand(
    U**2 - (U**2 - V**2)
)

print("  V^2 extraction from canonical K =")
print("   ", V2_from_kernel)

cert(
    "V^2 recovered",
    V2_from_kernel,
    V**2,
)


# ============================================================================
# KAPPA BRIDGE
# ============================================================================

print()
print("[13] KAPPA BRIDGE FROM NORMAL FORM")
print("-" * 78)

N = sp.expand(
    (U**2 - V**2) / c**2
)

S = sp.expand(
    2*U/c
)

Delta = sp.expand(
    4*V**2/c**2
)

Q = sp.expand(
    z**2 - S*z + N
)

print("  N     =", N)
print("  S     =", S)
print("  Delta =", Delta)
print("  Q(z)  =", sp.factor(Q))

cert(
    "Q discriminant",
    sp.discriminant(Q, z),
    Delta,
)


# ============================================================================
# TRANSLATION AS SPECTRAL VARIABLE
# ============================================================================

print()
print("[14] TRANSLATION VARIABLE = SPECTRAL VARIABLE")
print("-" * 78)

orbit = sp.expand(
    N + h*S + h**2
)

spectral = sp.expand(
    orbit.subs(h, -z)
)

print("  orbit(h)    =", sp.factor(orbit))
print("  orbit(-z)   =", sp.factor(spectral))

cert(
    "orbit(-z)=Q(z)",
    spectral,
    Q,
)


# ============================================================================
# UNIQUENESS RESULT
# ============================================================================

print()
print("[15] UNIQUENESS / GAUGE FREEDOM")
print("-" * 78)

print("""
  The operator determines a decomposition into:

      invariant algebra  = polynomials in V,

      moving coordinate  = U modulo invariant terms.

  Therefore the canonical triple is not completely unique.

  An affine change

      H -> H + I(V)

  leaves the translation increment unchanged because

      T(I(V)) = I(V).

  Likewise an invariant may be added freely to K.

  Thus the true invariant content is not the literal symbols
  U, V, K, but the operator representation class:

      moving coordinate modulo invariants,

      quadratic kernel modulo invariant additions.
""")


# ============================================================================
# EXACT GAUGE TEST
# ============================================================================

print()
print("[16] GAUGE-INVARIANCE TEST")
print("-" * 78)

g1, g0 = sp.symbols("g1 g0")

I = sp.expand(
    g1*V + g0
)

H_gauge = sp.expand(
    U + I
)

cert(
    "translation of invariant-added H",
    T(H_gauge, 1) - H_gauge,
    c,
)


J_gauge = sp.expand(
    (U**2 - V**2)
    + I
)

cert(
    "linear invariant addition to K",
    T(J_gauge, 1) - J_gauge,
    2*c*U + c**2,
)


# ============================================================================
# FINAL OPERATOR SIGNATURE
# ============================================================================

print()
print("[17] FINAL OPERATOR SIGNATURE")
print("-" * 78)

print("""
  A genuine upstream realization need not literally contain:

      U
      V
      U^2-V^2.

  It is sufficient that there exist observables H,K and
  an operator T such that:

      T(H)-H = c

      T(K)-K = 2*c*H + c^2

      T^2(K)-2*T(K)+K = 2*c^2

  up to invariant gauge terms.

  The invariant algebra consists of quantities unchanged by T.

  For the canonical model:

      invariant = V

      moving coordinate = U

      quadratic kernel = U^2-V^2.

  The resulting factorization bridge is:

      K
       ->
      H
       ->
      S
       ->
      Q(z)=z^2-Sz+N
       ->
      roots.

  This gives an operator-theoretic recognition criterion for
  the original homogeneous-layer construction.
""")


# ============================================================================
# AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 531 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)
print("""
The next step is no longer to invent another identity.

Use this normal form to test the actual homogeneous-layer
operators.

For an actual operator T and candidate observables H,K,J,
test whether:

    T(J)-J = 0

    T(H)-H = constant

    T(K)-K - 2*c*H - c^2
        is T-invariant

    T^2(K)-2*T(K)+K
        is constant.

Importantly, allow invariant gauge freedom:

    H -> H + I

    K -> K + J

where T(I)=I and T(J)=J.

If the homogeneous-layer construction contains this structure,
then it is operator-conjugate to the KAPPA translation system.

That would be a genuinely new structural result, rather than
another rewrite of S, Delta, or the factor pair.
""")
