#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 537 START")
print("=" * 78)
print("COORDINATE-FREE QUADRATIC FORM / TRANSLATION CONJUGACY TEST")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

x1, x2 = sp.symbols("x1 x2")
t1, t2 = sp.symbols("t1 t2", nonzero=True)

A, B, C = sp.symbols("A B C")
D, E, F = sp.symbols("D E F")

z, h = sp.symbols("z h")
c = sp.symbols("c", nonzero=True)


# ============================================================================
# CERTIFICATE
# ============================================================================

failures = 0


def simplify_expr(expr):
    return sp.factor(sp.expand(sp.simplify(expr)))


def cert(label, expr, expected=0):
    global failures

    diff = simplify_expr(expr - expected)
    ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# VARIABLES / QUADRATIC FORM
# ============================================================================

x = sp.Matrix([x1, x2])
t = sp.Matrix([t1, t2])

Q = sp.Matrix([
    [A, B / 2],
    [B / 2, C],
])

linear = sp.Matrix([D, E])

K = sp.expand(
    (x.T * Q * x)[0]
    + (linear.T * x)[0]
    + F
)


def translate(expr, n=1):
    return sp.expand(
        expr.subs(
            {
                x1: x1 + n * t1,
                x2: x2 + n * t2,
            },
            simultaneous=True,
        )
    )


# ============================================================================
# 1. GENERAL QUADRATIC FORM
# ============================================================================

print("[1] GENERAL QUADRATIC FORM")
print("-" * 78)

print("  Q-matrix =")
sp.pprint(Q)

print()
print("  K(x) =")
sp.pprint(K)


# ============================================================================
# 2. GENERAL TRANSLATION
# ============================================================================

print()
print("[2] GENERAL TRANSLATION RESPONSE")
print("-" * 78)

TK = translate(K, 1)
dK = sp.expand(TK - K)

print("  T(K)-K =")
sp.pprint(sp.factor(dK))


# ============================================================================
# 3. LINEAR RESPONSE
# ============================================================================

print()
print("[3] TRANSLATION GRADIENT CHANNEL")
print("-" * 78)

response_vector = sp.Matrix([
    sp.diff(dK, x1),
    sp.diff(dK, x2),
])

expected_response_vector = 2 * Q * t

print("  coefficient vector of x:")
sp.pprint(response_vector)

print()
print("  expected 2*Q*t:")
sp.pprint(expected_response_vector)

cert(
    "linear response = 2*Q*t",
    response_vector,
    expected_response_vector,
)


# ============================================================================
# 4. QUADRATIC FORM ALONG TRANSLATION DIRECTION
# ============================================================================

print()
print("[4] QUADRATIC FORM ON TRANSLATION DIRECTION")
print("-" * 78)

q_tt = sp.expand(
    (t.T * Q * t)[0]
)

print("  q(t,t) =")
sp.pprint(sp.factor(q_tt))


# ============================================================================
# 5. CONSTANT PART OF FIRST DIFFERENCE
# ============================================================================

print()
print("[5] CONSTANT TRANSLATION RESPONSE")
print("-" * 78)

linear_part = sp.expand(
    (response_vector.T * x)[0]
)

constant_part = sp.expand(
    dK - linear_part
)

expected_constant = sp.expand(
    q_tt + (linear.T * t)[0]
)

cert(
    "constant response",
    constant_part,
    expected_constant,
)


# ============================================================================
# 6. SECOND DIFFERENCE
# ============================================================================

print()
print("[6] SECOND DIFFERENCE")
print("-" * 78)

D2K = sp.expand(
    translate(K, 2)
    - 2 * translate(K, 1)
    + K
)

cert(
    "second difference = 2*q(t,t)",
    D2K,
    2 * q_tt,
)


# ============================================================================
# 7. MOVING COORDINATE
# ============================================================================

print()
print("[7] MOVING COORDINATE EXTRACTION")
print("-" * 78)

print("""
  When q(t,t) != 0, define the normalized moving coordinate

      H = (linear response)/(2*q(t,t)).
""")

H = sp.expand(
    linear_part / (2 * q_tt)
)

print("  H =")
sp.pprint(H)

cert(
    "T(H)-H = 1",
    translate(H, 1) - H,
    1,
)


# ============================================================================
# 8. INVARIANT LINEAR COORDINATE
# ============================================================================

print()
print("[8] INVARIANT LINEAR COORDINATE")
print("-" * 78)

j1, j2 = sp.symbols("j1 j2")

J = j1 * x1 + j2 * x2

print("  J = j1*x1 + j2*x2")

cert(
    "T(J)-J = j.t",
    translate(J, 1) - J,
    j1 * t1 + j2 * t2,
)

print("""
  Therefore an invariant linear form satisfies

      j1*t1 + j2*t2 = 0.
""")


# ============================================================================
# 9. CANONICAL INVARIANT
# ============================================================================

print()
print("[9] CANONICAL TRANSVERSE INVARIANT")
print("-" * 78)

Jcanonical = sp.expand(
    -t2 * x1 + t1 * x2
)

print("  J = -t2*x1 + t1*x2")

cert(
    "T(J)-J = 0",
    translate(Jcanonical, 1) - Jcanonical,
    0,
)


# ============================================================================
# 10. NORMALIZED QUADRATIC ORBIT
# ============================================================================

print()
print("[10] NORMALIZED QUADRATIC ORBIT")
print("-" * 78)

Kn = sp.expand(K / q_tt)
Lh = sp.expand(translate(Kn, h))

print("  L(h) =")
sp.pprint(sp.factor(Lh))


# ============================================================================
# 11. ORBIT COEFFICIENTS
# ============================================================================

print()
print("[11] ORBIT COEFFICIENT EXTRACTION")
print("-" * 78)

Lpoly = sp.Poly(sp.expand(Lh), h)

h2 = sp.expand(Lpoly.coeff_monomial(h**2))
h1 = sp.expand(Lpoly.coeff_monomial(h))
h0 = sp.expand(Lpoly.coeff_monomial(1))

print("  coefficient h^2 =")
sp.pprint(sp.factor(h2))

print("  coefficient h   =")
sp.pprint(sp.factor(h1))

print("  coefficient 1   =")
sp.pprint(sp.factor(h0))

cert(
    "quadratic orbit coefficient",
    h2,
    1,
)


# ============================================================================
# 12. ORBIT FINITE DIFFERENCES
# ============================================================================

print()
print("[12] ORBIT FINITE DIFFERENCES")
print("-" * 78)

D_L = sp.expand(
    Lh.subs(h, h + 1) - Lh
)

D2_L = sp.expand(
    Lh.subs(h, h + 2)
    - 2 * Lh.subs(h, h + 1)
    + Lh
)

print("  ΔL(h) =")
sp.pprint(sp.factor(D_L))

print()
print("  Δ²L(h) =")
sp.pprint(sp.factor(D2_L))

cert(
    "Δ²L = 2",
    D2_L,
    2,
)


# ============================================================================
# 13. GENERATING POLYNOMIAL
# ============================================================================

print()
print("[13] GENERATING POLYNOMIAL")
print("-" * 78)

Qorbit = sp.expand(
    Lh.subs(h, -z)
)

print("  Q(z)=L(-z):")
sp.pprint(sp.factor(Qorbit))

Qpoly = sp.Poly(Qorbit, z)

q2 = sp.expand(Qpoly.coeff_monomial(z**2))
q1 = sp.expand(Qpoly.coeff_monomial(z))
q0 = sp.expand(Qpoly.coeff_monomial(1))

print()
print("  z^2 coefficient =")
sp.pprint(sp.factor(q2))

print("  z coefficient =")
sp.pprint(sp.factor(q1))

print("  constant =")
sp.pprint(sp.factor(q0))

cert(
    "Q quadratic coefficient",
    q2,
    1,
)


# ============================================================================
# 14. DISCRIMINANT
# ============================================================================

print()
print("[14] SPECTRAL DISCRIMINANT")
print("-" * 78)

disc_Q = sp.factor(
    sp.discriminant(Qorbit, z)
)

print("  discriminant(Q) =")
sp.pprint(disc_Q)


# ============================================================================
# 15. TRANSLATION COMPOSITION
# ============================================================================

print()
print("[15] TRANSLATION COMPOSITION")
print("-" * 78)

composition_difference = sp.expand(
    translate(translate(K, 1), 1)
    - translate(K, 2)
)

cert(
    "T_1(T_1(K)) = T_2(K)",
    composition_difference,
    0,
)


# ============================================================================
# 16. MINKOWSKI SPECIALIZATION
# ============================================================================

print()
print("[16] MINKOWSKI SPECIALIZATION")
print("-" * 78)

U, V = sp.symbols("U V")

minkowski_subs = {
    A: 1,
    B: 0,
    C: -1,
    D: 0,
    E: 0,
    F: 0,
    t1: c,
    t2: 0,
    x1: U,
    x2: V,
}

K_m = sp.expand(
    K.subs(minkowski_subs)
)

qtt_m = sp.expand(
    q_tt.subs(minkowski_subs)
)

H_m = sp.expand(
    H.subs(minkowski_subs)
)

J_m = sp.expand(
    Jcanonical.subs(minkowski_subs)
)

print("  K =")
sp.pprint(K_m)

print()
print("  q(t,t) =")
sp.pprint(qtt_m)

print()
print("  H =")
sp.pprint(H_m)

print()
print("  J =")
sp.pprint(J_m)

cert(
    "Minkowski kernel",
    K_m,
    U**2 - V**2,
)

cert(
    "Minkowski q(t,t)",
    qtt_m,
    c**2,
)

cert(
    "Minkowski moving coordinate",
    H_m,
    U / c,
)

cert(
    "Minkowski invariant coordinate",
    J_m,
    c * V,
)


# ============================================================================
# 17. MINKOWSKI NORMALIZED ORBIT
# ============================================================================

print()
print("[17] MINKOWSKI NORMALIZED ORBIT")
print("-" * 78)

L_m = sp.expand(
    (U + c*h)**2 / c**2
    - V**2 / c**2
)

print("  L(h) =")
sp.pprint(sp.expand(L_m))

expected_L = sp.expand(
    h**2 + (2*U/c)*h + (U**2 - V**2)/c**2
)

cert(
    "L(h)=h^2+S*h+N",
    L_m,
    expected_L,
)


# ============================================================================
# 18. MINKOWSKI GENERATING QUADRATIC
# ============================================================================

print()
print("[18] MINKOWSKI KAPPA QUADRATIC")
print("-" * 78)

Q_m = sp.expand(
    L_m.subs(h, -z)
)

expected_Q = sp.expand(
    z**2
    - (2*U/c)*z
    + (U**2 - V**2)/c**2
)

print("  Q(z) =")
sp.pprint(Q_m)

cert(
    "Q(z)=z^2-Sz+N",
    Q_m,
    expected_Q,
)


# ============================================================================
# 19. MINKOWSKI DISCRIMINANT
# ============================================================================

print()
print("[19] MINKOWSKI DISCRIMINANT CHANNEL")
print("-" * 78)

disc_m = sp.factor(
    sp.discriminant(Q_m, z)
)

print("  discriminant =")
sp.pprint(disc_m)

cert(
    "discriminant = 4V^2/c^2",
    disc_m,
    4 * V**2 / c**2,
)


# ============================================================================
# 20. HISTORICAL x,y SPECIALIZATION
# ============================================================================

print()
print("[20] HISTORICAL x,y PULLBACK")
print("-" * 78)

x_hist, y_hist = sp.symbols("x_hist y_hist")

U_hist = 2 * y_hist
V_hist = 2 * x_hist - 3

K_hist = sp.expand(
    U_hist**2 - V_hist**2
)

v_hist = sp.expand(
    y_hist**2
    - x_hist**2
    + 3*x_hist
    - 2
)

print("  U = 2y")
print("  V = 2x-3")

print()
print("  K =")
sp.pprint(sp.factor(K_hist))

print()
print("  4v-1 =")
sp.pprint(sp.factor(4*v_hist - 1))

cert(
    "K = 4v-1",
    K_hist,
    4*v_hist - 1,
)


# ============================================================================
# 21. HISTORICAL TRANSLATION
# ============================================================================

print()
print("[21] HISTORICAL TRANSLATION")
print("-" * 78)

historical_shift = sp.expand(
    v_hist.subs(
        y_hist,
        y_hist + c*h/2,
    )
)

historical_K_shift = sp.expand(
    4*historical_shift - 1
)

expected_hist_K = sp.expand(
    (U_hist + c*h)**2 - V_hist**2
)

cert(
    "4*T(v)-1 = translated Minkowski kernel",
    historical_K_shift,
    expected_hist_K,
)


# ============================================================================
# 22. HISTORICAL NORMALIZED ORBIT
# ============================================================================

print()
print("[22] HISTORICAL NORMALIZED ORBIT")
print("-" * 78)

hist_L = sp.expand(
    historical_K_shift / c**2
)

expected_hist_L = sp.expand(
    h**2
    + (4*y_hist/c)*h
    + (4*v_hist - 1)/c**2
)

cert(
    "historical orbit is quadratic",
    hist_L,
    expected_hist_L,
)


# ============================================================================
# 23. COORDINATE-FREE SIGNATURE
# ============================================================================

print()
print("[23] COORDINATE-FREE KAPPA SIGNATURE")
print("-" * 78)

print("""
  A translation-plus-quadratic system is a candidate for the
  KAPPA mechanism when it contains:

      translation vector t

      quadratic form q

      q(t,t) != 0

  and therefore has:

      T^2(K)-2T(K)+K = constant,

  while the first difference is linear.

  In the distinguished Minkowski realization:

      moving coordinate = U/c

      invariant coordinate = V

      quadratic kernel = (U^2-V^2)/c^2

      orbit = h^2 + S*h + N

      Q(z) = orbit(-z)

      discriminant(Q) = Delta.
""")


# ============================================================================
# 24. UPSTREAM RECOGNITION TEST
# ============================================================================

print()
print("[24] UPSTREAM RECOGNITION CRITERION")
print("-" * 78)

print("""
  For the ORIGINAL homogeneous-layer construction, the next
  test should look for an actual transformation T and observable K
  satisfying:

      T^2(K)-2T(K)+K = nonzero constant

  and:

      T(K)-K = linear_observable + constant.

  Then determine whether the linear response itself transforms
  by a nonzero constant.

  Finally identify an independent invariant observable J.

  The desired hierarchy is:

      K
      |
      +--> first difference
      |
      +--> moving linear observable H
      |
      +--> invariant transverse observable J
      |
      v
      normalized orbit
      |
      v
      Q(z)
      |
      v
      discriminant.


  This is the operator-level signature to search for upstream.
""")


# ============================================================================
# FINAL
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 537 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT STRUCTURAL TARGET")
print("=" * 78)

print("""
The remaining problem is not another downstream identity.

The coordinate-free construction shows that a quadratic form plus
a non-null translation direction automatically produces:

    a quadratic translation orbit,
    a linear first-difference channel,
    an invariant transverse direction,
    and a generating quadratic.

The next experiment should therefore take the ACTUAL homogeneous-layer
transformation rules from the earlier construction and compute their
action on its lowest-degree observables.

The decisive test is:

    T^2(K)-2T(K)+K = constant

and whether:

    T(K)-K

produces a linear observable whose own translation response is
constant.

If such a chain exists, identify its invariant transverse channel
and construct the corresponding Q(z).

Do not insert p or q.
Do not fit numerical data.
Do not enumerate divisors.
Do not use continued fractions.

The objective is an exact operator-level identification.
""")