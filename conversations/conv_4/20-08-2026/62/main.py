#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 538 START")
print("=" * 78)
print("COORDINATE-FREE QUADRATIC ORBIT -> KAPPA NORMAL FORM RECONSTRUCTION")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

x1, x2 = sp.symbols("x1 x2")
t1, t2 = sp.symbols("t1 t2", nonzero=True)

A, B, C = sp.symbols("A B C")
D, E, F = sp.symbols("D E F")

h, z = sp.symbols("h z")
c = sp.symbols("c", nonzero=True)

rho = sp.symbols("rho", nonzero=True)


# ============================================================================
# HELPERS
# ============================================================================

failures = 0


def simp(expr):
    return sp.factor(sp.expand(sp.simplify(expr)))


def cert(label, expr, expected=0):
    global failures

    diff = simp(expr - expected)

    # Handle matrices explicitly.
    if isinstance(diff, sp.MatrixBase):
        ok = all(simp(v) == 0 for v in diff)
    else:
        ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# 1. GENERAL QUADRATIC FORM
# ============================================================================

print("[1] GENERAL QUADRATIC FORM")
print("-" * 78)

Qmat = sp.Matrix([
    [A, B / 2],
    [B / 2, C],
])

x = sp.Matrix([x1, x2])
t = sp.Matrix([t1, t2])

linear = sp.Matrix([D, E])

K = sp.expand(
    (x.T * Qmat * x)[0]
    + (linear.T * x)[0]
    + F
)

print("  Q =")
sp.pprint(Qmat)

print()
print("  K(x) =")
sp.pprint(K)


# ============================================================================
# 2. TRANSLATION OPERATOR
# ============================================================================

print()
print("[2] TRANSLATION OPERATOR")
print("-" * 78)

def T(expr, n=1):
    return sp.expand(
        expr.subs(
            {
                x1: x1 + n * t1,
                x2: x2 + n * t2,
            },
            simultaneous=True,
        )
    )


K1 = T(K)
K2 = T(K, 2)

dK = sp.expand(K1 - K)
d2K = sp.expand(K2 - 2 * K1 + K)

print("  T(K)-K =")
sp.pprint(sp.factor(dK))

print()
print("  T²(K)-2T(K)+K =")
sp.pprint(sp.factor(d2K))


# ============================================================================
# 3. TRANSLATION QUADRATIC NORM
# ============================================================================

print()
print("[3] TRANSLATION QUADRATIC NORM")
print("-" * 78)

qtt = sp.expand(
    (t.T * Qmat * t)[0]
)

print("  q(t,t) =")
sp.pprint(sp.factor(qtt))

cert(
    "second difference = 2 q(t,t)",
    d2K,
    2 * qtt,
)


# ============================================================================
# 4. FIRST DIFFERENCE AS A LINEAR FUNCTION
# ============================================================================

print()
print("[4] FIRST-DIFFERENCE LINEAR CHANNEL")
print("-" * 78)

response = sp.expand(
    dK - qtt - (linear.T * t)[0]
)

expected_linear_response = sp.expand(
    2 * (Qmat * t).dot(x)
)

cert(
    "linear part = 2 q(t,x)",
    response,
    expected_linear_response,
)


# ============================================================================
# 5. MOVING COVECTOR
# ============================================================================

print()
print("[5] MOVING COVECTOR")
print("-" * 78)

moving_covector = sp.simplify(
    2 * Qmat * t
)

print("  m = 2 Q t")
sp.pprint(moving_covector)

moving_value = sp.expand(
    (moving_covector.T * x)[0]
)

print()
print("  moving observable M(x) =")
sp.pprint(moving_value)


# ============================================================================
# 6. NORMALIZED MOVING COORDINATE
# ============================================================================

print()
print("[6] NORMALIZED MOVING COORDINATE")
print("-" * 78)

# Since:
#
#   T(M)-M = 2 q(t,t)
#
# define H = M/(2 q(t,t)).
#
# Then:
#
#   T(H)-H = 1.

H = sp.cancel(
    moving_value / (2 * qtt)
)

print("  H =")
sp.pprint(H)

cert(
    "T(H)-H = 1",
    T(H) - H,
    1,
)


# ============================================================================
# 7. INVARIANT COVECTOR
# ============================================================================

print()
print("[7] TRANSVERSE INVARIANT COVECTOR")
print("-" * 78)

# In two dimensions, a covector annihilating t is:
#
#   j = (-t2, t1)

j = sp.Matrix([
    -t2,
    t1,
])

J = sp.expand(
    (j.T * x)[0]
)

print("  J = -t2*x1 + t1*x2")

cert(
    "J(t) = 0",
    (j.T * t)[0],
    0,
)

cert(
    "T(J)-J = 0",
    T(J) - J,
    0,
)


# ============================================================================
# 8. COORDINATE TRANSFORMATION
# ============================================================================

print()
print("[8] MOVING / INVARIANT COORDINATE MAP")
print("-" * 78)

Hraw = sp.expand(
    moving_value
)

Jraw = sp.expand(J)

print("  H_raw =")
sp.pprint(Hraw)

print()
print("  J_raw =")
sp.pprint(Jraw)

# Jacobian determinant of (Hraw,Jraw) with respect to (x1,x2).
jac = sp.Matrix([
    [sp.diff(Hraw, x1), sp.diff(Hraw, x2)],
    [sp.diff(Jraw, x1), sp.diff(Jraw, x2)],
])

jac_det = sp.factor(
    sp.det(jac)
)

print()
print("  Jacobian determinant =")
sp.pprint(jac_det)


# ============================================================================
# 9. DISCRIMINANT OF THE QUADRATIC ORBIT
# ============================================================================

print()
print("[9] GENERATING QUADRATIC")
print("-" * 78)

L = sp.expand(
    T(K, 1)  # keep the construction visibly tied to T
)

# Instead construct the normalized orbit from K itself.
#
#   L_h = [K(x+h t) - invariant linear/constant gauges] / q(t,t)
#
# The coefficients relevant to the orbit are obtained directly from
# K(x+h t).

Kh = sp.expand(
    K.subs(
        {
            x1: x1 + h * t1,
            x2: x2 + h * t2,
        },
        simultaneous=True,
    )
)

Lh = sp.cancel(Kh / qtt)

Lpoly = sp.Poly(
    sp.expand(Lh),
    h,
)

c2 = sp.expand(
    Lpoly.coeff_monomial(h**2)
)

c1 = sp.expand(
    Lpoly.coeff_monomial(h)
)

c0 = sp.expand(
    Lpoly.coeff_monomial(1)
)

print("  normalized orbit coefficient h² =")
sp.pprint(c2)

print()
print("  normalized orbit coefficient h =")
sp.pprint(c1)

print()
print("  normalized orbit constant =")
sp.pprint(c0)

cert(
    "orbit h² coefficient = 1",
    c2,
    1,
)


# ============================================================================
# 10. GENERATING POLYNOMIAL
# ============================================================================

print()
print("[10] SPECTRAL GENERATING POLYNOMIAL")
print("-" * 78)

Qorbit = sp.expand(
    Lh.subs(h, -z)
)

Qpoly = sp.Poly(
    Qorbit,
    z,
)

q_z2 = sp.expand(
    Qpoly.coeff_monomial(z**2)
)

q_z1 = sp.expand(
    Qpoly.coeff_monomial(z)
)

q_z0 = sp.expand(
    Qpoly.coeff_monomial(1)
)

print("  Q(z)=L(-z)")
sp.pprint(Qorbit)

print()
print("  z² coefficient =")
sp.pprint(q_z2)

print()
print("  z coefficient =")
sp.pprint(q_z1)

print()
print("  constant =")
sp.pprint(q_z0)

cert(
    "Q is monic",
    q_z2,
    1,
)


# ============================================================================
# 11. DISCRIMINANT
# ============================================================================

print()
print("[11] DISCRIMINANT")
print("-" * 78)

disc_Q = sp.factor(
    sp.discriminant(Qorbit, z)
)

print("  discriminant(Q) =")
sp.pprint(disc_Q)


# ============================================================================
# 12. GENERAL DISCRIMINANT INVARIANT
# ============================================================================

print()
print("[12] TRANSVERSE INVARIANT CHANNEL")
print("-" * 78)

# The key coordinate-free object is J².
#
# In the canonical Minkowski realization this becomes V²
# and the discriminant becomes 4J² after normalization.

J2 = sp.expand(
    J**2
)

print("  J² =")
sp.pprint(J2)

# Verify that J² is invariant.
cert(
    "T(J²)-J² = 0",
    T(J2) - J2,
    0,
)


# ============================================================================
# 13. PURE MINKOWSKI CONDITION
# ============================================================================

print()
print("[13] MINKOWSKI / INDEFINITE SPECIALIZATION")
print("-" * 78)

U, V = sp.symbols("U V")

mink_subs = {
    A: 1,
    B: 0,
    C: -1,
    D: 0,
    E: 0,
    F: 0,
    x1: U,
    x2: V,
    t1: c,
    t2: 0,
}

Km = sp.expand(
    K.subs(mink_subs)
)

Hm = sp.expand(
    H.subs(mink_subs)
)

Jm = sp.expand(
    J.subs(mink_subs)
)

qmm = sp.expand(
    qtt.subs(mink_subs)
)

print("  K =")
sp.pprint(Km)

print()
print("  H =")
sp.pprint(Hm)

print()
print("  J =")
sp.pprint(Jm)

print()
print("  q(t,t) =")
sp.pprint(qmm)

cert(
    "K = U²-V²",
    Km,
    U**2 - V**2,
)

cert(
    "H = U/c",
    Hm,
    U / c,
)

cert(
    "J = c V",
    Jm,
    c * V,
)

cert(
    "q(t,t)=c²",
    qmm,
    c**2,
)


# ============================================================================
# 14. MINKOWSKI NORMALIZED ORBIT
# ============================================================================

print()
print("[14] MINKOWSKI NORMAL FORM")
print("-" * 78)

Lm = sp.expand(
    ((U + h * c)**2 - V**2) / c**2
)

expected_Lm = sp.expand(
    h**2
    + (2 * U / c) * h
    + (U**2 - V**2) / c**2
)

cert(
    "L(h)=h²+S*h+N",
    Lm,
    expected_Lm,
)


# ============================================================================
# 15. MINKOWSKI KAPPA POLYNOMIAL
# ============================================================================

print()
print("[15] MINKOWSKI KAPPA POLYNOMIAL")
print("-" * 78)

Qm = sp.expand(
    Lm.subs(h, -z)
)

expected_Qm = sp.expand(
    z**2
    - (2 * U / c) * z
    + (U**2 - V**2) / c**2
)

cert(
    "Q(z)=z²-Sz+N",
    Qm,
    expected_Qm,
)

disc_m = sp.factor(
    sp.discriminant(Qm, z)
)

print()
print("  discriminant =")
sp.pprint(disc_m)

cert(
    "Delta = 4V²/c²",
    disc_m,
    4 * V**2 / c**2,
)


# ============================================================================
# 16. HISTORICAL CONIC PULLBACK
# ============================================================================

print()
print("[16] HISTORICAL CONIC PULLBACK")
print("-" * 78)

xh, yh = sp.symbols("x y")

Uh = 2 * yh
Vh = 2 * xh - 3

v = sp.expand(
    yh**2 - xh**2 + 3*xh - 2
)

Kh = sp.expand(
    Uh**2 - Vh**2
)

cert(
    "U²-V² = 4v-1",
    Kh,
    4 * v - 1,
)


# ============================================================================
# 17. HISTORICAL TRANSLATION
# ============================================================================

print()
print("[17] HISTORICAL TRANSLATION OPERATOR")
print("-" * 78)

alpha = sp.symbols("alpha", nonzero=True)

v_shift = sp.expand(
    v.subs(
        yh,
        yh + alpha * c / 2,
    )
)

K_shift = sp.expand(
    4 * v_shift - 1
)

expected_shift = sp.expand(
    (Uh + alpha*c)**2 - Vh**2
)

cert(
    "translated historical kernel",
    K_shift,
    expected_shift,
)


# ============================================================================
# 18. DISTINCTIVE KAPPA SIGNATURE
# ============================================================================

print()
print("[18] DISTINCTIVE KAPPA SIGNATURE")
print("-" * 78)

print("""
  A generic quadratic-plus-translation system gives a quadratic
  orbit automatically.

  The KAPPA-specific structure requires more:

      1. a non-null translation direction;
      2. a moving linear observable;
      3. an independent invariant linear observable;
      4. a quadratic form with nonzero discriminant;
      5. after normalization, the quadratic must reduce to

             H^2 - J^2;

      6. the orbit must therefore become

             L_m = m^2 + S*m + N;

      7. L(-z) must become

             z^2 - S*z + N;

      8. the discriminant must be carried entirely by the
         invariant transverse coordinate.

  Thus the next upstream search should not stop at

      Delta^2(K) = constant.

  It must identify the full triple:

      moving observable H
      invariant observable J
      indefinite quadratic K.
""")


# ============================================================================
# 19. RECOGNITION CONDITIONS
# ============================================================================

print()
print("[19] COORDINATE-FREE RECOGNITION CONDITIONS")
print("-" * 78)

print("""
  Given an actual homogeneous-layer operator T:

      A. Find a quadratic observable K.

      B. Compute

             D1(K) = T(K)-K

         and verify it is affine-linear.

      C. Extract its moving linear component H.

      D. Verify

             T(H)-H = constant != 0.

      E. Find an independent invariant observable J:

             T(J)-J = 0.

      F. Check whether K is gauge-equivalent to

             H^2 - J^2.

  The key negative test is also important:

      if the transverse invariant enters with the wrong signature,
      or if K cannot be reduced to H^2-J^2,

  then the operator is not in the historical/KAPPA conjugacy class.
""")


# ============================================================================
# 20. FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 538 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The generic calculation has now been reduced to an intrinsic
recognition problem.

The next experiment should use the ACTUAL homogeneous-layer
transformation, not another synthetic translation.

For that transformation T, search its low-degree observable space
for:

    (T-I)H = c

    (T-I)J = 0

    (T-I)K = 2*c*H + c² + invariant

and then test whether the quadratic part of K has indefinite
signature relative to the moving/invariant decomposition.

The decisive object is:

    invariant direction
          +
    moving direction
          +
    indefinite quadratic lift.

Only that full structure can establish a genuine conjugacy to the
historical/KAPPA system.

Do not introduce p or q.
Do not enumerate factors.
Do not use continued fractions.
Do not fit numerical data.
""")
