#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 535 START")
print("=" * 78)
print("AFFINE TRANSLATION CONJUGACY / SCALE-CORRECTED KAPPA NORMAL FORM")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

x, y = sp.symbols("x y")
xp, yp = sp.symbols("xp yp")

alpha, beta = sp.symbols("alpha beta", nonzero=True)
n = sp.symbols("n", integer=True)

c = sp.symbols("c", nonzero=True)

h, z = sp.symbols("h z")

A, B, C, D, E, F = sp.symbols("A B C D E F")

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
# GENERAL PURE TRANSLATION
# ============================================================================

def T(expr, k=1):
    return sp.expand(
        expr.subs({
            x: x + k * alpha,
            y: y + k * beta,
        })
    )


print()
print("[1] GENERAL AFFINE TRANSLATION")
print("-" * 78)

print("""
  T(x) = x + alpha
  T(y) = y + beta

  translation vector:

      t = (alpha,beta)

  define the Euclidean translation norm:

      rho^2 = alpha^2 + beta^2.
""")

rho2 = sp.expand(alpha**2 + beta**2)

print("  rho^2 =", rho2)


# ============================================================================
# CANONICAL MOVING / INVARIANT LINEAR FORMS
# ============================================================================

print()
print("[2] CANONICAL MOVING / INVARIANT BASIS")
print("-" * 78)

H = sp.expand(
    alpha*x + beta*y
)

J = sp.expand(
    -beta*x + alpha*y
)

print("  H =", H)
print("  J =", J)

cert(
    "T(H)-H = rho^2",
    T(H) - H,
    rho2,
)

cert(
    "T(J)-J = 0",
    T(J) - J,
    0,
)

print("""
  Thus:

      H  = moving coordinate
      J  = invariant coordinate

  with

      T(H) = H + rho^2
      T(J) = J.
""")


# ============================================================================
# NORMALIZED MOVING STEP
# ============================================================================

print()
print("[3] UNIT MOVING COORDINATE")
print("-" * 78)

Hn = sp.expand(H / rho2)

print("  H_normalized =", Hn)

cert(
    "T(H_normalized)-H_normalized = 1",
    T(Hn) - Hn,
    1,
)


# ============================================================================
# QUADRATIC CONIC ADAPTED TO TRANSLATION
# ============================================================================

print()
print("[4] CANONICAL QUADRATIC KERNEL")
print("-" * 78)

# The natural quadratic with moving coordinate H and invariant J is:
#
#     K = H^2 - rho^2*J^2
#
# This normalization is chosen so that after division by rho^4
# the orbit has unit quadratic curvature.

K = sp.expand(
    H**2 - rho2 * J**2
)

print("  K =")
print("   ", sp.factor(K))


# ============================================================================
# TRANSLATION RESPONSE
# ============================================================================

print()
print("[5] QUADRATIC TRANSLATION RESPONSE")
print("-" * 78)

dK = sp.expand(
    T(K) - K
)

print("  T(K)-K =")
print("   ", sp.factor(dK))

cert(
    "quadratic response",
    dK,
    2 * rho2 * H + rho2**2,
)


# ============================================================================
# SECOND DIFFERENCE
# ============================================================================

print()
print("[6] SECOND DIFFERENCE")
print("-" * 78)

K0 = K
K1 = T(K, 1)
K2 = T(K, 2)

D2K = sp.expand(
    K2 - 2*K1 + K0
)

print("  T^2(K)-2T(K)+K =")
print("   ", sp.factor(D2K))

cert(
    "second difference",
    D2K,
    2 * rho2**2,
)


# ============================================================================
# CORRECT NORMALIZATION
# ============================================================================

print()
print("[7] CORRECT NORMALIZED ORBIT")
print("-" * 78)

# Divide by rho^4.
#
# The normalized moving coordinate is H/rho^2.
#
# L_n = T^n(K)/rho^4
#
# Then:
#
#     L_n = n^2 + (2H/rho^2)n + K/rho^4.

L_n = sp.expand(
    T(K, n) / rho2**2
)

S_eff = sp.expand(
    2 * H / rho2
)

N_eff = sp.expand(
    K / rho2**2
)

L_expected = sp.expand(
    n**2 + S_eff*n + N_eff
)

cert(
    "L_n quadratic orbit",
    L_n,
    L_expected,
)


# ============================================================================
# FIRST DIFFERENCE
# ============================================================================

print()
print("[8] FIRST-DIFFERENCE CHANNEL")
print("-" * 78)

L0 = sp.expand(T(K, n) / rho2**2)
L1 = sp.expand(T(K, n + 1) / rho2**2)

dL = sp.expand(L1 - L0)

print("  Delta L_n =")
print("   ", sp.factor(dL))

cert(
    "Delta L_n",
    dL,
    S_eff + 2*n + 1,
)


# ============================================================================
# SECOND DIFFERENCE OF NORMALIZED ORBIT
# ============================================================================

print()
print("[9] NORMALIZED CURVATURE")
print("-" * 78)

L2 = sp.expand(T(K, n + 2) / rho2**2)

D2L = sp.expand(
    L2 - 2*L1 + L0
)

print("  Delta^2 L_n =", sp.factor(D2L))

cert(
    "normalized curvature = 2",
    D2L,
    2,
)


# ============================================================================
# GENERATING POLYNOMIAL
# ============================================================================

print()
print("[10] GENERATING QUADRATIC")
print("-" * 78)

Q = sp.expand(
    L_expected.subs(n, -z)
)

print("  Q(z) = L(-z)")
print("   ", sp.factor(Q))

Q_expected = sp.expand(
    z**2 - S_eff*z + N_eff
)

cert(
    "Q(z)=z^2-Sz+N",
    Q,
    Q_expected,
)


# ============================================================================
# DISCRIMINANT
# ============================================================================

print()
print("[11] DISCRIMINANT")
print("-" * 78)

discQ = sp.factor(
    sp.discriminant(Q, z)
)

print("  discriminant(Q) =")
print("   ", discQ)

# Since Q = z^2 - S_eff z + N_eff,
# discriminant = 4 J^2 / rho^2.

Delta_eff = sp.expand(
    4 * J**2 / rho2
)

cert(
    "discriminant = 4 J^2/rho^2",
    discQ,
    Delta_eff,
)


# ============================================================================
# RECONSTRUCT J^2
# ============================================================================

print()
print("[12] INVARIANT GAP CHANNEL")
print("-" * 78)

J2_recovered = sp.expand(
    discQ * rho2 / 4
)

cert(
    "J^2 recovery",
    J2_recovered,
    J**2,
)


# ============================================================================
# RECONSTRUCT MOVING COORDINATE
# ============================================================================

print()
print("[13] MOVING COORDINATE RECOVERY")
print("-" * 78)

H_recovered = sp.expand(
    (T(K) - K - rho2**2) / (2*rho2)
)

print("  H_recovered =")
print("   ", sp.factor(H_recovered))

cert(
    "recover H from one translated K",
    H_recovered,
    H,
)


# ============================================================================
# RECONSTRUCT EFFECTIVE S
# ============================================================================

print()
print("[14] EFFECTIVE S RECOVERY")
print("-" * 78)

S_recovered = sp.expand(
    2 * H_recovered / rho2
)

cert(
    "S_eff recovery",
    S_recovered,
    S_eff,
)


# ============================================================================
# RECONSTRUCT EFFECTIVE N
# ============================================================================

print()
print("[15] EFFECTIVE N RECOVERY")
print("-" * 78)

N_recovered = sp.expand(
    K / rho2**2
)

cert(
    "N_eff recovery",
    N_recovered,
    N_eff,
)


# ============================================================================
# DIRECT OPERATOR FORM
# ============================================================================

print()
print("[16] OPERATOR NORMAL FORM")
print("-" * 78)

print("""
  After the scale correction:

      H' = H/rho^2
      J' = J/rho

  satisfy:

      T(H') = H' + 1
      T(J') = J'

  and the normalized quadratic is:

      L = H'^2 - J'^2.

  Therefore:

      T(L)-L = 2H' + 1

      T^2(L)-2T(L)+L = 2.
""")


cert(
    "T(H')-H' = 1",
    T(Hn) - Hn,
    1,
)


J_normalized = sp.expand(
    J / sp.sqrt(rho2)
)

# We do not need a square-root normalization in the symbolic
# structural test. The invariant contribution is already captured
# by J^2/rho^2.
print("  invariant channel J^2/rho^2 =", sp.factor(J**2 / rho2))


# ============================================================================
# SPECIALIZATION TO HISTORICAL MINKOWSKI COORDINATES
# ============================================================================

print()
print("[17] HISTORICAL MINKOWSKI SPECIALIZATION")
print("-" * 78)

U, V, r = sp.symbols("U V r", nonzero=True)

# Historical translation:
#
#     U -> U+r
#     V -> V
#
# Corresponds to ordinary coordinate translation alpha=r, beta=0.

historical_subs = {
    alpha: r,
    beta: 0,
    x: U,
    y: V,
}

rho2_hist = sp.expand(
    rho2.subs(historical_subs)
)

H_hist = sp.expand(
    H.subs(historical_subs)
)

J_hist = sp.expand(
    J.subs(historical_subs)
)

K_hist = sp.expand(
    K.subs(historical_subs)
)

print("  rho^2 =", rho2_hist)
print("  H     =", H_hist)
print("  J     =", J_hist)
print("  K     =", sp.factor(K_hist))

cert(
    "historical rho^2 = r^2",
    rho2_hist,
    r**2,
)

cert(
    "historical H = r*U",
    H_hist,
    r*U,
)

cert(
    "historical J = -r*V",
    J_hist,
    -r*V,
)

cert(
    "historical K = r^4*(U^2-V^2)?",
    K_hist,
    r**2 * (r**2) * (U**2 - V**2),
)


# ============================================================================
# DIRECT HISTORICAL NORMALIZATION
# ============================================================================

print()
print("[18] HISTORICAL NORMALIZED RECOVERY")
print("-" * 78)

# K_hist / r^4 = U^2/r^2 - V^2/r^2.
#
# H_hist/r^2 = U/r.
# Therefore effective S = 2U/r.

N_hist = sp.expand(
    K_hist / r**4
)

S_hist = sp.expand(
    2 * H_hist / r**3
)

Delta_hist = sp.expand(
    4 * J_hist**2 / r**4
)

print("  N_hist =", sp.factor(N_hist))
print("  S_hist =", sp.factor(S_hist))
print("  Delta_hist =", sp.factor(Delta_hist))

cert(
    "N_hist = (U^2-V^2)/r^2",
    N_hist,
    (U**2 - V**2) / r**2,
)

cert(
    "S_hist = 2U/r",
    S_hist,
    2*U/r,
)

cert(
    "Delta_hist = 4V^2/r^2",
    Delta_hist,
    4*V**2/r**2,
)


# ============================================================================
# AFFINE-CONJUGACY CRITERION
# ============================================================================

print()
print("[19] AFFINE-CONJUGACY CRITERION")
print("-" * 78)

print("""
  The scale-corrected canonical structure is:

      T(H') = H' + 1
      T(J') = J'

      L = H'^2 - J'^2

  with:

      T(L)-L = 2H' + 1
      T^2(L)-2T(L)+L = 2.

  Therefore an upstream operator is a candidate whenever,
  after an invertible affine change of observables, it admits:

      one invariant linear coordinate,
      one translated linear coordinate,
      one quadratic whose pure moving part is square-minus-invariant.
""")


# ============================================================================
# GENERAL QUADRATIC GAUGE
# ============================================================================

print()
print("[20] QUADRATIC GAUGE FREEDOM")
print("-" * 78)

g2, g1, g0 = sp.symbols("g2 g1 g0")

I_J = sp.expand(
    g2*J**2 + g1*J + g0
)

K_gauged = sp.expand(
    K + I_J
)

cert(
    "invariant quadratic gauge preserves response",
    T(K_gauged) - K_gauged,
    T(K) - K,
)


# ============================================================================
# LINEAR GAUGE
# ============================================================================

print()
print("[21] MOVING-COORDINATE GAUGE")
print("-" * 78)

q1, q0 = sp.symbols("q1 q0")

H_gauged = sp.expand(
    H + q1*J + q0
)

cert(
    "H gauge preserves translation increment",
    T(H_gauged) - H_gauged,
    T(H) - H,
)


# ============================================================================
# WHAT IS INTRINSIC?
# ============================================================================

print()
print("[22] INTRINSIC OPERATOR DATA")
print("-" * 78)

print("""
  The following quantities survive the coordinate-gauge freedom:

      rho^2
          translation scale

      invariant algebra generated by J
          unchanged observables

      moving quotient H mod invariant functions
          translated coordinate

      quadratic orbit curvature
          2

      discriminant channel
          proportional to J^2.

  Thus the real upstream search target is not a literal U/V
  representation.

  It is an affine operator possessing:

      invariant subspace
          +
      translated quotient
          +
      quadratic lift.
""")


# ============================================================================
# FINAL OPERATOR CERTIFICATE
# ============================================================================

print()
print("[23] FINAL OPERATOR CERTIFICATE")
print("-" * 78)

print("""
  SCALE-CORRECTED KAPPA NORMAL FORM:

      T(H') = H' + 1
      T(J') = J'

      L = H'^2 - J'^2

      T(L)-L = 2H' + 1

      T^2(L)-2T(L)+L = 2

      Q(z)=L(-z)

      Q(z)=z^2-Sz+N

      Delta = discriminant(Q).

  The previous apparent failure for arbitrary (alpha,beta)
  was a normalization error rather than a structural failure.

  The operator itself carries an intrinsic step-size rho^2.
  Once that step is normalized, the curvature is universally 2.
""")


# ============================================================================
# FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 535 FINISHED")
print("=" * 78)
print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The next experiment should stop using an invented translation
operator.

Instead, represent the ORIGINAL homogeneous-layer transformation
as an actual map on its lowest-degree observable space.

Then compute:

    T - I

and classify:

    ker(T-I)

    image(T-I)

    generalized eigenspaces.

The desired signature is:

    dim ker(T-I) >= 1

    plus a generalized moving observable H with

        (T-I)H = constant

    plus a quadratic K satisfying

        (T-I)K = 2H+1 + invariant.

The decisive question is now:

    Does the ORIGINAL homogeneous-layer operator have a
    one-dimensional affine translation quotient with a
    quadratic lift?

If yes, construct its conjugacy explicitly.

If no, that is equally useful: it would prove that the
KAPPA translation structure is not generated by that operator.

Do not add p or q to the upstream search.
Do not fit numerical values.
Do not search continued fractions.
The next experiment is an operator/Jordan classification.
""")
