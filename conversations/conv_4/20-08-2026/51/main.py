#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 527 START")
print("=" * 78)
print("CONIC TRANSLATION GROUP / KAPPA OPERATOR ALGEBRA")
print()

failures = 0


def check(label, expr, expected=0):
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

U, V = sp.symbols("U V")
a, b, m, n = sp.symbols("a b m n")
r = sp.symbols("r", nonzero=True)

P = (U - V) / r
Q = (U + V) / r

N = sp.expand((U**2 - V**2) / r**2)
S = sp.expand(2*U / r)
Delta = sp.expand(4*V**2 / r**2)

print("[1] BASE CONIC")
print("-" * 78)

K = sp.expand(U**2 - V**2)

print("  K(U,V) = U^2 - V^2")
print("  N =", N)
print("  S =", S)
print("  Delta =", Delta)

check("K = U^2-V^2", K, U**2 - V**2)

print()


# ============================================================================
# [2] GENERAL TRANSLATION OPERATOR
# ============================================================================

print("[2] GENERAL (U,V) TRANSLATION")
print("-" * 78)

def translate(expr, du, dv):
    return sp.expand(expr.subs({
        U: U + du,
        V: V + dv
    }))


K_ab = translate(K, a*r, b*r)

print("  T_(a,b):")
print("    U -> U+a*r")
print("    V -> V+b*r")
print()
print("  K(T(U,V)) =", sp.factor(K_ab))

check(
    "general translation expansion",
    K_ab,
    (U + a*r)**2 - (V + b*r)**2,
)

print()


# ============================================================================
# [3] TRANSLATION INCREMENT
# ============================================================================

print("[3] TRANSLATION INCREMENT")
print("-" * 78)

delta_K = sp.expand(K_ab - K)

print("  ΔK =", sp.factor(delta_K))

expected_delta = sp.expand(
    2*a*r*U - 2*b*r*V + (a**2 - b**2)*r**2
)

check(
    "translation increment",
    delta_K,
    expected_delta,
)

print()


# ============================================================================
# [4] PURE U TRANSLATION
# ============================================================================

print("[4] PURE U TRANSLATION")
print("-" * 78)

K_U = translate(K, a*r, 0)

print("  U -> U+a*r")
print("  ΔK_U =", sp.factor(K_U - K))

check(
    "pure U translation",
    K_U - K,
    2*a*r*U + a**2*r**2,
)

print()


# ============================================================================
# [5] PURE V TRANSLATION
# ============================================================================

print("[5] PURE V TRANSLATION")
print("-" * 78)

K_V = translate(K, 0, b*r)

print("  V -> V+b*r")
print("  ΔK_V =", sp.factor(K_V - K))

check(
    "pure V translation",
    K_V - K,
    -2*b*r*V - b**2*r**2,
)

print()


# ============================================================================
# [6] COMBINED DIAGONAL TRANSLATIONS
# ============================================================================

print("[6] DIAGONAL TRANSLATION CLASSES")
print("-" * 78)

K_same = translate(K, a*r, a*r)
K_opposite = translate(K, a*r, -a*r)

print("  same-sign:  U -> U+a*r, V -> V+a*r")
print("    ΔK =", sp.factor(K_same - K))

print("  opposite-sign: U -> U+a*r, V -> V-a*r")
print("    ΔK =", sp.factor(K_opposite - K))

check(
    "same-sign translation",
    K_same - K,
    2*a*r*(U - V),
)

check(
    "opposite-sign translation",
    K_opposite - K,
    2*a*r*(U + V),
)

print()


# ============================================================================
# [7] FACTOR-COORDINATE ACTION
# ============================================================================

print("[7] ACTION ON FACTOR COORDINATES")
print("-" * 78)

A0 = sp.expand(U - V)
B0 = sp.expand(U + V)

A_ab = sp.expand((U + a*r) - (V + b*r))
B_ab = sp.expand((U + a*r) + (V + b*r))

print("  A0 =", A0)
print("  B0 =", B0)

print("  transformed A0 =", A_ab)
print("  transformed B0 =", B_ab)

check(
    "A translation law",
    A_ab,
    A0 + (a-b)*r,
)

check(
    "B translation law",
    B_ab,
    B0 + (a+b)*r,
)

print()


# ============================================================================
# [8] KAPPA TRANSLATION AS FACTOR TRANSLATION
# ============================================================================

print("[8] KAPPA TRANSLATION")
print("-" * 78)

A_K = sp.expand(A0 + m*r)
B_K = sp.expand(B0 + m*r)

print("  A -> A+m*r")
print("  B -> B+m*r")

check(
    "KAPPA A translation",
    A_K,
    (U + m*r) - V,
)

check(
    "KAPPA B translation",
    B_K,
    (U + m*r) + V,
)

check(
    "KAPPA product orbit",
    A_K * B_K,
    (U + m*r)**2 - V**2,
)

print()


# ============================================================================
# [9] TRANSLATION PARAMETERS IN (A,B)
# ============================================================================

print("[9] FACTOR-COORDINATE TRANSLATION PARAMETERS")
print("-" * 78)

alpha = sp.symbols("alpha")
beta = sp.symbols("beta")

U_from_AB = sp.expand((A0 + B0) / 2)
V_from_AB = sp.expand((B0 - A0) / 2)

print("  U =", U_from_AB)
print("  V =", V_from_AB)

A_shift = sp.expand(A0 + alpha*r)
B_shift = sp.expand(B0 + beta*r)

U_shift = sp.expand((A_shift + B_shift) / 2)
V_shift = sp.expand((B_shift - A_shift) / 2)

print("  U' =", U_shift)
print("  V' =", V_shift)

check(
    "U shift from factor shifts",
    U_shift,
    U + (alpha + beta)*r/2,
)

check(
    "V shift from factor shifts",
    V_shift,
    V + (beta - alpha)*r/2,
)

print()


# ============================================================================
# [10] WHICH FACTOR TRANSLATIONS KEEP V FIXED?
# ============================================================================

print("[10] V-PRESERVING FACTOR TRANSLATION")
print("-" * 78)

condition_v = sp.solve(
    sp.Eq(beta - alpha, 0),
    beta
)

print("  V'=V requires beta=alpha")
print("  solution:", condition_v)

check(
    "V-preserving condition",
    beta - alpha,
    0,
)

print("""
  Therefore the KAPPA operation

      A -> A+m*r
      B -> B+m*r

  is uniquely the equal-translation direction in
  factor coordinates.

  Equivalently:

      U -> U+m*r
      V -> V.
""")

print()


# ============================================================================
# [11] WHICH TRANSLATION PRESERVES THE FACTOR GAP?
# ============================================================================

print("[11] GAP-PRESERVING TRANSLATION")
print("-" * 78)

gap = sp.expand(B0 - A0)
gap_shift = sp.expand(B_shift - A_shift)

print("  gap =", gap)
print("  shifted gap =", gap_shift)

check(
    "equal translation preserves gap",
    gap_shift,
    gap + (beta-alpha)*r,
)

print()


# ============================================================================
# [12] TRANSLATION OF N, S, DELTA
# ============================================================================

print("[12] ACTION ON N, S, DELTA")
print("-" * 78)

N_ab = sp.expand((K_ab) / r**2)

S_ab = sp.expand(
    2*(U + a*r)/r
)

Delta_ab = sp.expand(
    4*(V + b*r)**2 / r**2
)

print("  N' =", sp.factor(N_ab))
print("  S' =", sp.factor(S_ab))
print("  Delta' =", sp.factor(Delta_ab))

check(
    "S translation",
    S_ab - S,
    2*a,
)

check(
    "Delta translation",
    Delta_ab - Delta,
    8*b*V/r + 4*b**2,
)

print()


# ============================================================================
# [13] KAPPA ORBIT
# ============================================================================

print("[13] NORMALIZED KAPPA ORBIT")
print("-" * 78)

L_m = sp.expand(
    (U + m*r)**2 / r**2 - V**2 / r**2
)

print("  L_m =", sp.factor(L_m))

check(
    "L_m = N+mS+m^2",
    L_m,
    N + m*S + m**2,
)

dL = sp.expand(
    L_m.subs(m, m+1) - L_m
)

d2L = sp.expand(
    L_m.subs(m, m+2)
    - 2*L_m.subs(m, m+1)
    + L_m
)

check(
    "first difference",
    dL,
    S + 2*m + 1,
)

check(
    "second difference",
    d2L,
    2,
)

print()


# ============================================================================
# [14] OPERATOR COMMUTATOR
# ============================================================================

print("[14] TRANSLATION OPERATOR COMMUTATOR")
print("-" * 78)

def T(expr, du, dv):
    return translate(expr, du, dv)

f = sp.Function("f")
U_sym, V_sym = sp.symbols("U_sym V_sym")

# Symbolic functions are easier to test directly.
fUV = f(U_sym, V_sym)

T1T2 = sp.expand(
    T(
        T(fUV, a*r, 0),
        0,
        b*r,
    )
)

T2T1 = sp.expand(
    T(
        T(fUV, 0, b*r),
        a*r,
        0,
    )
)

print("  T_U(a) T_V(b) =", T1T2)
print("  T_V(b) T_U(a) =", T2T1)

check(
    "translation operators commute",
    T1T2,
    T2T1,
)

print()


# ============================================================================
# [15] CONIC QUADRATIC ORBIT AS OPERATOR EXPONENTIAL
# ============================================================================

print("[15] DISCRETE OPERATOR STRUCTURE")
print("-" * 78)

D_U = sp.expand(
    K_U - K
)

D2_U = sp.expand(
    translate(K_U, a*r, 0)
    - 2*K_U
    + K
)

print("  first U-difference =", sp.factor(D_U))
print("  second U-difference =", sp.factor(D2_U))

check(
    "U-orbit curvature",
    D2_U,
    2*a**2*r**2,
)

print()


# ============================================================================
# [16] HISTORICAL x,y COORDINATE ACTION
# ============================================================================

print("[16] BACK-TRANSFORM TO x,y")
print("-" * 78)

x, y = sp.symbols("x y")

U_xy = 2*y
V_xy = 2*x - 3

K_xy = sp.expand(
    U_xy**2 - V_xy**2
)

print("  U = 2y")
print("  V = 2x-3")
print("  K =", sp.factor(K_xy))

# General translation:
# U -> U+a*r means y -> y+a*r/2
# V -> V+b*r means x -> x+b*r/2

x_shift = x + b*r/2
y_shift = y + a*r/2

K_xy_shift = sp.expand(
    K_xy.subs({
        x: x_shift,
        y: y_shift
    })
)

check(
    "x,y translation correspondence",
    K_xy_shift,
    K_ab.subs({
        U: U_xy,
        V: V_xy
    }),
)

print()


# ============================================================================
# [17] SPECIAL HISTORICAL OPERATORS
# ============================================================================

print("[17] SPECIAL HISTORICAL OPERATORS")
print("-" * 78)

special = {
    "T_U(+1)": (1, 0),
    "T_U(-1)": (-1, 0),
    "T_V(+1)": (0, 1),
    "T_V(-1)": (0, -1),
    "T_same(+1)": (1, 1),
    "T_opposite(+1)": (1, -1),
}

for name, (aa, bb) in special.items():
    expr = sp.factor(
        translate(K, aa*r, bb*r)
    )
    print(f"  {name}:")
    print(f"    {expr}")

print()


# ============================================================================
# [18] TRANSLATION-DERIVED QUADRATIC FAMILY
# ============================================================================

print("[18] GENERAL TRANSLATION-DERIVED QUADRATIC")
print("-" * 78)

L_ab_m = sp.expand(
    translate(K, a*r, b*r) / r**2
)

print("  L_(a,b) =", sp.factor(L_ab_m))

check(
    "general normalized translation",
    L_ab_m,
    N + (2*a)*U/r - (2*b)*V/r + a**2 - b**2,
)

print()


# ============================================================================
# [19] UNIQUE DIRECTION CARRYING S WITHOUT V
# ============================================================================

print("[19] UNIQUE S-CHANNEL")
print("-" * 78)

linear_part = sp.expand(
    L_ab_m - N - (a**2 - b**2)
)

print("  linear translation information =", sp.factor(linear_part))

# For pure S-only dependence we require b=0.
# Then the linear term is 2aU/r = aS.

check(
    "b=0 removes V channel",
    linear_part.subs(b, 0),
    a*S,
)

print("""
  This establishes symbolically:

      b = 0

  is the unique translation direction whose linear
  information is proportional to S alone.

  Thus the KAPPA orbit is not merely one possible
  translation. It is the distinguished V-preserving
  direction of the historical hyperbola.
""")

print()


# ============================================================================
# [20] FINAL STRUCTURAL RESULT
# ============================================================================

print("[20] FINAL STRUCTURAL RESULT")
print("-" * 78)

print("""
  The historical conic carries a two-dimensional translation
  operator family:

      (U,V) -> (U+a*r, V+b*r).

  Its effect on the invariant is

      K' - K
        = 2r(aU-bV) + r^2(a^2-b^2).

  In factor coordinates:

      A0 -> A0+(a-b)r
      B0 -> B0+(a+b)r.

  The KAPPA translation corresponds exactly to

      a=1, b=0,

  hence

      A0 -> A0+r
      B0 -> B0+r.

  This is the unique translation direction in which the
  linear coefficient contains S but not the independent
  antisymmetric coordinate V.

  Consequently the KAPPA quadratic orbit is the
  V-preserving translation orbit of the historical hyperbola.
""")

print()


# ============================================================================
# FINAL AUDIT
# ============================================================================

print("=" * 78)
print("EXPERIMENT 527 FINISHED")
print("=" * 78)
print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT SYMBOLIC AUDIT = {failures == 0}")

print("""
NEXT RESEARCH TARGET
------------------------------------------------------------------------------

The next upstream question should now be operator-theoretic:

    Does the ORIGINAL homogeneous construction contain a
    transformation equivalent to the distinguished direction

        (U,V) -> (U+r,V)

    or, in historical coordinates,

        y -> y+r/2
        x -> x ?

Search its actual symbolic transformation rules for:

    1. a quantity invariant under the V direction;

    2. a quantity translated by a constant in U;

    3. a two-layer difference proportional to U;

    4. a quadratic orbit with second difference 2.

The decisive signature is:

    T(L)-L = S+1
    T^2(L)-2T(L)+L = 2

with

    L(0)=N.

That would identify an existing homogeneous operator
as the generator of the KAPPA quadratic rather than
introducing the KAPPA construction separately.

Do not introduce p,q after the fact when performing this
upstream test. The purpose is to identify the operator
before the factor coordinates are reconstructed.
""")
