#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 536 START")
print("=" * 78)
print("AFFINE TRANSLATION + GENERAL QUADRATIC FORM / MINKOWSKI-KAPPA CLASSIFICATION")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

U, V = sp.symbols("U V")
Up, Vp = sp.symbols("Up Vp")

a, b = sp.symbols("a b")
alpha, beta = sp.symbols("alpha beta")

A, B, C, D, E, F = sp.symbols(
    "A B C D E F"
)

h = sp.symbols("h")
z = sp.symbols("z")

failures = 0


# ============================================================================
# CERTIFICATE
# ============================================================================

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
# GENERAL TRANSLATION
# ============================================================================

def T(expr, step=1):
    return sp.expand(
        expr.subs({
            U: U + step * alpha,
            V: V + step * beta,
        })
    )


print()
print("[1] GENERAL TRANSLATION")
print("-" * 78)

print("""
  T_h:

      U -> U + h*alpha
      V -> V + h*beta

  Translation direction:

      d = (alpha,beta)
""")


# ============================================================================
# GENERAL QUADRATIC FORM
# ============================================================================

print()
print("[2] GENERAL QUADRATIC FORM")
print("-" * 78)

K = sp.expand(
    A*U**2
    + B*U*V
    + C*V**2
    + D*U
    + E*V
    + F
)

print("  K(U,V) =")
print("   ", K)


# ============================================================================
# TRANSLATION RESPONSE
# ============================================================================

print()
print("[3] TRANSLATION RESPONSE")
print("-" * 78)

TK = sp.expand(T(K, 1))
dK = sp.expand(TK - K)

print("  T(K)-K =")
print("   ", sp.factor(dK))


# ============================================================================
# SECOND DIFFERENCE
# ============================================================================

print()
print("[4] SECOND DIFFERENCE")
print("-" * 78)

T2K = sp.expand(T(K, 2))

D2K = sp.expand(
    T2K - 2*TK + K
)

print("  T^2(K)-2T(K)+K =")
print("   ", sp.factor(D2K))

# For a quadratic form, the second difference must be constant.
second_expected = sp.expand(
    2*(A*alpha**2 + B*alpha*beta + C*beta**2)
)

cert(
    "quadratic second-difference law",
    D2K,
    second_expected,
)


# ============================================================================
# PURE MOVING-CHANNEL CONDITION
# ============================================================================

print()
print("[5] PURE MOVING CHANNEL")
print("-" * 78)

print("""
  We want:

      T(K)-K

  to depend on one linear moving coordinate only.

  Therefore its linear U,V part must define a single direction.
""")


linear_U = sp.expand(
    sp.diff(dK, U)
)

linear_V = sp.expand(
    sp.diff(dK, V)
)

print("  linear U coefficient =", sp.factor(linear_U))
print("  linear V coefficient =", sp.factor(linear_V))


# ============================================================================
# TRANSLATION-DIRECTION QUADRATIC
# ============================================================================

print()
print("[6] QUADRATIC FORM ALONG THE TRANSLATION")
print("-" * 78)

direction_value = sp.expand(
    A*alpha**2
    + B*alpha*beta
    + C*beta**2
)

print("  q(d,d) =", sp.factor(direction_value))

print("""
  The second difference is controlled only by:

      q(d,d).

  Thus the translation orbit curvature is a property of the
  quadratic form evaluated on the translation direction.
""")


# ============================================================================
# HISTORICAL MINKOWSKI FORM
# ============================================================================

print()
print("[7] HISTORICAL MINKOWSKI FORM")
print("-" * 78)

K_minkowski = sp.expand(
    U**2 - V**2
)

print("  K_M = U^2 - V^2")

TKM = sp.expand(
    T(K_minkowski, 1)
)

dKM = sp.expand(
    TKM - K_minkowski
)

D2KM = sp.expand(
    T(K_minkowski, 2)
    - 2*TKM
    + K_minkowski
)

print("  T(K_M)-K_M =")
print("   ", sp.factor(dKM))

print("  T^2(K_M)-2T(K_M)+K_M =")
print("   ", sp.factor(D2KM))


# ============================================================================
# PURE U TRANSLATION
# ============================================================================

print()
print("[8] HISTORICAL PURE-U TRANSLATION")
print("-" * 78)

historical_subs = {
    alpha: 1,
    beta: 0,
}

dKM_U = sp.factor(
    dKM.subs(historical_subs)
)

D2KM_U = sp.factor(
    D2KM.subs(historical_subs)
)

print("  U -> U+1, V -> V")

cert(
    "Minkowski first response",
    dKM_U,
    2*U + 1,
)

cert(
    "Minkowski second difference",
    D2KM_U,
    2,
)


# ============================================================================
# GENERAL STEP c
# ============================================================================

print()
print("[9] GENERAL MINKOWSKI STEP")
print("-" * 78)

c = sp.symbols("c", nonzero=True)

KM_c = sp.expand(
    (U + c*h)**2 - V**2
)

dKM_c = sp.expand(
    KM_c.subs(h, h + 1)
    - KM_c
)

D2KM_c = sp.expand(
    KM_c.subs(h, h + 2)
    - 2*KM_c.subs(h, h + 1)
    + KM_c
)

print("  L_h = (U+c*h)^2 - V^2")

cert(
    "first difference",
    dKM_c,
    2*c*U + 2*c**2*h + c**2,
)

cert(
    "second difference",
    D2KM_c,
    2*c**2,
)


# ============================================================================
# NORMALIZED MINKOWSKI ORBIT
# ============================================================================

print()
print("[10] NORMALIZED MINKOWSKI ORBIT")
print("-" * 78)

Lh = sp.expand(
    KM_c / c**2
)

L_expected = sp.expand(
    h**2
    + (2*U/c)*h
    + (U**2 - V**2)/c**2
)

cert(
    "L_h = h^2 + S*h + N",
    Lh,
    L_expected,
)


# ============================================================================
# FIRST-DIFFERENCE EXTRACTION
# ============================================================================

print()
print("[11] FIRST-DIFFERENCE S EXTRACTION")
print("-" * 78)

DLh = sp.expand(
    Lh.subs(h, h + 1)
    - Lh
)

S_expected = sp.expand(
    2*U/c
)

print("  Delta L_h =")
print("   ", sp.factor(DLh))

cert(
    "Delta L_h = S + 2h + 1",
    DLh,
    S_expected + 2*h + 1,
)

S_from_difference = sp.expand(
    DLh - (2*h + 1)
)

cert(
    "S extracted from first difference",
    S_from_difference,
    S_expected,
)


# ============================================================================
# KAPPA QUADRATIC
# ============================================================================

print()
print("[12] KAPPA GENERATING QUADRATIC")
print("-" * 78)

Q = sp.expand(
    Lh.subs(h, -z)
)

Q_expected = sp.expand(
    z**2
    - S_expected*z
    + (U**2 - V**2)/c**2
)

print("  Q(z) =")
print("   ", sp.factor(Q))

cert(
    "Q(z)=z^2-Sz+N",
    Q,
    Q_expected,
)


# ============================================================================
# DISCRIMINANT
# ============================================================================

print()
print("[13] DISCRIMINANT CHANNEL")
print("-" * 78)

discQ = sp.factor(
    sp.discriminant(Q, z)
)

Delta_expected = sp.expand(
    4*V**2/c**2
)

print("  discriminant(Q) =", discQ)

cert(
    "discriminant = 4V^2/c^2",
    discQ,
    Delta_expected,
)


# ============================================================================
# CENTER / MOVING / INVARIANT DECOMPOSITION
# ============================================================================

print()
print("[14] OPERATOR COORDINATE DECOMPOSITION")
print("-" * 78)

H = sp.expand(U/c)
J2 = sp.expand(V**2/c**2)

print("  moving coordinate H =", H)
print("  invariant quadratic J^2 =", J2)

cert(
    "T(H)-H = 1",
    sp.expand(
        H.subs(U, U + c) - H
    ),
    1,
)

cert(
    "T(J^2)-J^2 = 0",
    sp.expand(
        J2.subs(V, V)
        - J2
    ),
    0,
)


# ============================================================================
# HISTORICAL x,y PULLBACK
# ============================================================================

print()
print("[15] HISTORICAL x,y PULLBACK")
print("-" * 78)

x, y, r = sp.symbols("x y r", nonzero=True)

U_xy = sp.expand(2*y)
V_xy = sp.expand(2*x - 3)

K_xy = sp.expand(
    U_xy**2 - V_xy**2
)

print("  U = 2y")
print("  V = 2x-3")

print("  K =")
print("   ", sp.factor(K_xy))

v_xy = sp.expand(
    -x**2 + 3*x + y**2 - 2
)

cert(
    "K = 4v-1",
    K_xy,
    4*v_xy - 1,
)


# ============================================================================
# HISTORICAL y TRANSLATION
# ============================================================================

print()
print("[16] HISTORICAL y TRANSLATION")
print("-" * 78)

yshift = sp.expand(
    y + r*h/2
)

v_shift = sp.expand(
    v_xy.subs(y, yshift)
)

K_shift = sp.expand(
    4*v_shift - 1
)

K_historical = sp.expand(
    (U_xy + r*h)**2 - V_xy**2
)

cert(
    "y -> y+r*h/2 reproduces U -> U+r*h",
    K_shift,
    K_historical,
)


# ============================================================================
# HISTORICAL NORMALIZED ORBIT
# ============================================================================

print()
print("[17] HISTORICAL NORMALIZED ORBIT")
print("-" * 78)

L_hist = sp.expand(
    K_historical / r**2
)

N_hist = sp.expand(
    (U_xy**2 - V_xy**2) / r**2
)

S_hist = sp.expand(
    2*U_xy/r
)

cert(
    "historical orbit",
    L_hist,
    h**2 + S_hist*h + N_hist,
)


# ============================================================================
# HISTORICAL FIRST DIFFERENCE
# ============================================================================

print()
print("[18] HISTORICAL FIRST-DIFFERENCE CHANNEL")
print("-" * 78)

DL_hist = sp.expand(
    L_hist.subs(h, h + 1)
    - L_hist
)

cert(
    "Delta L - (2h+1) = S",
    DL_hist - (2*h + 1),
    S_hist,
)


# ============================================================================
# CRITICAL COMPARISON
# ============================================================================

print()
print("[19] EUCLIDEAN VS MINKOWSKI STRUCTURE")
print("-" * 78)

print("""
  The earlier generic experiment used a translation-adapted
  Euclidean construction.

  The historical KAPPA system instead uses the indefinite form:

      K = U^2 - V^2.

  This matters because the invariant channel is not produced
  by an ordinary Euclidean norm.

  The exact historical decomposition is:

      moving direction:
          U

      invariant direction:
          V

      quadratic form:
          U^2 - V^2.

  Therefore the correct upstream search must preserve the
  SIGNATURE of the quadratic form, not merely the existence
  of a translation orbit.
""")


# ============================================================================
# OPERATOR SIGNATURE
# ============================================================================

print()
print("[20] FINAL OPERATOR SIGNATURE")
print("-" * 78)

print("""
  A genuine upstream realization should admit observables
  (H,J,K) and a step c satisfying:

      T(H) = H + c

      T(J) = J

      K = H^2 - J^2
          up to invariant gauge terms

      T(K)-K = 2cH + c^2

      T^2(K)-2T(K)+K = 2c^2.

  After normalization by c:

      L_h = h^2 + S*h + N

      S = 2H/c

      N = K/c^2

      Delta = 4J^2/c^2.

  The important structural requirement is now stronger than
  "there is a translation":

      there must be an INDEFINITE quadratic form whose
      translation direction is one coordinate and whose
      transverse coordinate is invariant.
""")


# ============================================================================
# FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 536 FINISHED")
print("=" * 78)

print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The generic Euclidean affine construction is now separated from
the historical Minkowski construction.

The next upstream experiment should therefore search the ORIGINAL
homogeneous-layer algebra for a quadratic form Q and transformation
T satisfying:

    rank(Q) = 2

    discriminant(Q) != 0

    T has one invariant linear direction

    T has one translated linear direction

    Q is equivalent to:

        H^2 - J^2

    and:

        T(Q)-Q = 2cH+c^2.

The decisive issue is the quadratic-form SIGNATURE.

Do not merely search for:

    second difference = 2.

Search for the full triple:

    invariant direction
    moving direction
    indefinite quadratic form.

If the original homogeneous construction contains that triple,
then it can be conjugated into the historical/KAPPA normal form.

That is the next structural bridge.
""")
