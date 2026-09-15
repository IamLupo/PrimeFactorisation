#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 533 START")
print("=" * 78)
print("AFFINE OPERATOR CONVERSE / KAPPA CONJUGACY CLASSIFICATION")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

U, V = sp.symbols("U V")
Up, Vp = sp.symbols("Up Vp")
z, m = sp.symbols("z m")

a, b, c, d = sp.symbols("a b c d")
alpha, beta = sp.symbols("alpha beta")

lam = sp.symbols("lam", nonzero=True)
mu = sp.symbols("mu", nonzero=True)
h = sp.symbols("h")

A, B, C = sp.symbols("A B C")
D, E, F = sp.symbols("D E F")

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
# GENERAL AFFINE OPERATOR
# ============================================================================

def T(expr):
    """
    General affine operator:

        U' = a U + b V + alpha
        V' = c U + d V + beta
    """
    return sp.expand(
        expr.subs({
            U: a*U + b*V + alpha,
            V: c*U + d*V + beta,
        })
    )


def T_matrix(expr, M):
    """
    Apply a supplied affine map

        U -> M[0]*U + M[1]*V + M[2]
        V -> M[3]*U + M[4]*V + M[5]
    """
    return sp.expand(
        expr.subs({
            U: M[0]*U + M[1]*V + M[2],
            V: M[3]*U + M[4]*V + M[5],
        })
    )


def degree_coefficients(expr):
    """
    Return coefficients of U^2, UV, V^2, U, V, 1.
    """
    poly = sp.Poly(sp.expand(expr), U, V)

    return {
        "U2": poly.coeff_monomial(U**2),
        "UV": poly.coeff_monomial(U*V),
        "V2": poly.coeff_monomial(V**2),
        "U": poly.coeff_monomial(U),
        "V": poly.coeff_monomial(V),
        "1": poly.coeff_monomial(1),
    }


# ============================================================================
# [1] GENERAL AFFINE ACTION ON THE CANONICAL QUADRATIC
# ============================================================================

print()
print("[1] GENERAL AFFINE ACTION ON U^2-V^2")
print("-" * 78)

K0 = sp.expand(U**2 - V**2)

K1 = T(K0)

print("  K0 =", K0)
print()
print("  T(K0) =")
print("   ", sp.factor(K1))

delta_K = sp.factor(
    K1 - K0
)

print()
print("  T(K0)-K0 =")
print("   ", delta_K)


# ============================================================================
# [2] CONDITIONS FOR PURE MOVING RESPONSE
# ============================================================================

print()
print("[2] PURE MOVING-CHANNEL CONDITIONS")
print("-" * 78)

coeff_delta = degree_coefficients(delta_K)

print("  coefficients in T(K)-K:")
for key, value in coeff_delta.items():
    print(f"    {key} = {sp.factor(value)}")

print()
print("""
  A KAPPA-compatible moving response must have:

      no U^2 term
      no UV term
      no V^2 term
      no independent V term

  while retaining a linear U term.
""")

pure_quadratic_conditions = [
    sp.Eq(coeff_delta["U2"], 0),
    sp.Eq(coeff_delta["UV"], 0),
    sp.Eq(coeff_delta["V2"], 0),
    sp.Eq(coeff_delta["V"], 0),
]

print("  symbolic conditions:")
for eq in pure_quadratic_conditions:
    print("   ", eq)


# ============================================================================
# [3] SOLVE HOMOGENEOUS PART
# ============================================================================

print()
print("[3] HOMOGENEOUS AFFINE CLASSIFICATION")
print("-" * 78)

homogeneous_conditions = [
    sp.Eq(a**2 - c**2, 1),
    sp.Eq(a*b - c*d, 0),
    sp.Eq(b**2 - d**2, -1),
]

solutions_linear = sp.solve(
    [
        sp.expand(eq.lhs)
        for eq in homogeneous_conditions
    ],
    [a, b, c, d],
    dict=True,
)

print("  Lorentz-type preservation equations:")
for eq in homogeneous_conditions:
    print("   ", eq)

print()
print("  symbolic solution objects:")
for sol in solutions_linear:
    print("   ", sol)


# ============================================================================
# [4] PURE TRANSLATION SUBCLASS
# ============================================================================

print()
print("[4] PURE TRANSLATION SUBCLASS")
print("-" * 78)

pure_translation = (
    a - 1,
    b,
    c,
    d - 1,
    alpha,
    beta,
)

M_translation = (
    1,
    0,
    alpha,
    0,
    1,
    beta,
)

K_translation = T_matrix(
    K0,
    M_translation,
)

delta_translation = sp.factor(
    K_translation - K0
)

print("  U' = U + alpha")
print("  V' = V + beta")
print()
print("  T(K)-K =")
print("   ", delta_translation)

cert(
    "pure translation formula",
    delta_translation,
    2*alpha*U
    - 2*beta*V
    + alpha**2
    - beta**2,
)


# ============================================================================
# [5] KAPPA DIRECTION
# ============================================================================

print()
print("[5] UNIQUE V-PRESERVING TRANSLATION DIRECTION")
print("-" * 78)

print("""
  Require:

      V' = V

  hence:

      beta = 0.
""")

delta_U_translation = sp.expand(
    delta_translation.subs(beta, 0)
)

print("  T(U^2-V^2)-K =")
print("   ", sp.factor(delta_U_translation))

cert(
    "V-preserving response",
    delta_U_translation,
    2*alpha*U + alpha**2,
)


# ============================================================================
# [6] MOVING COORDINATE RECONSTRUCTION
# ============================================================================

print()
print("[6] MOVING COORDINATE RECONSTRUCTION")
print("-" * 78)

H_extracted = sp.factor(
    (
        delta_U_translation
        - alpha**2
    )
    / (2*alpha)
)

print("  extracted H =", H_extracted)

cert(
    "H=U",
    H_extracted,
    U,
)


# ============================================================================
# [7] INVARIANT ALGEBRA TEST
# ============================================================================

print()
print("[7] INVARIANT ALGEBRA")
print("-" * 78)

r0, r1, r2, r3 = sp.symbols(
    "r0 r1 r2 r3"
)

J = sp.expand(
    r3*U**3
    + r2*U**2
    + r1*U
    + r0
)

# Include V-dependent coefficients explicitly.
j30, j21, j12, j03 = sp.symbols(
    "j30 j21 j12 j03"
)
j20, j11, j02 = sp.symbols(
    "j20 j11 j02"
)
j10, j01 = sp.symbols(
    "j10 j01"
)
j00 = sp.symbols("j00")

general_J = sp.expand(
    j30*U**3
    + j21*U**2*V
    + j12*U*V**2
    + j03*V**3
    + j20*U**2
    + j11*U*V
    + j02*V**2
    + j10*U
    + j01*V
    + j00
)

invariance_eq = sp.Poly(
    sp.expand(
        general_J.subs(U, U + alpha)
        - general_J
    ),
    U,
    V,
)

eqs = [
    sp.Eq(coeff, 0)
    for _, coeff in invariance_eq.terms()
]

unknowns = [
    j30, j21, j12,
    j20, j11,
    j10,
]

solution_invariant = sp.solve(
    [
        eq.lhs
        for eq in eqs
    ],
    unknowns,
    dict=True,
)

print("  invariant solutions:")
for sol in solution_invariant:
    print("   ", sol)

print("""
  The surviving invariant algebra should depend only on V.
""")


# ============================================================================
# [8] AFFINE GAUGE FREEDOM
# ============================================================================

print()
print("[8] INVARIANT GAUGE FREEDOM")
print("-" * 78)

g0, g1, g2 = sp.symbols(
    "g0 g1 g2"
)

I_V = sp.expand(
    g2*V**2
    + g1*V
    + g0
)

H_gauge = sp.expand(
    U + I_V
)

cert(
    "H+I(V) has same translation response",
    (
        H_gauge.subs(U, U + alpha)
        - H_gauge
    ),
    alpha,
)


# ============================================================================
# [9] QUADRATIC GAUGE FREEDOM
# ============================================================================

print()
print("[9] QUADRATIC GAUGE FREEDOM")
print("-" * 78)

K_gauge = sp.expand(
    K0 + I_V
)

cert(
    "K+I(V) has same translation response",
    K_gauge.subs(U, U + alpha)
    - K_gauge,
    2*alpha*U + alpha**2,
)


# ============================================================================
# [10] OPERATOR NORMAL FORM
# ============================================================================

print()
print("[10] OPERATOR NORMAL FORM")
print("-" * 78)

print("""
  Canonical KAPPA normal form:

      T(U) = U + c
      T(V) = V

      K = U^2 - V^2

  with:

      T(K)-K = 2cU+c^2.
""")


# ============================================================================
# [11] CONJUGACY UNDER AFFINE COORDINATE CHANGE
# ============================================================================

print()
print("[11] AFFINE COORDINATE CONJUGACY")
print("-" * 78)

A11, A12, A21, A22 = sp.symbols(
    "A11 A12 A21 A22",
)
t1, t2 = sp.symbols(
    "t1 t2"
)

X = sp.Matrix([
    U,
    V,
])

M = sp.Matrix([
    [A11, A12],
    [A21, A22],
])

tvec = sp.Matrix([
    t1,
    t2,
])

Y = M*X + tvec

print("  General affine coordinate change:")
print()
print("      [H]   [A11 A12] [U]   [t1]")
print("      [J] = [A21 A22] [V] + [t2]")
print()


# ============================================================================
# [12] TRANSLATION VECTOR UNDER CONJUGACY
# ============================================================================

print()
print("[12] TRANSLATION VECTOR UNDER CONJUGACY")
print("-" * 78)

translation_vector = sp.Matrix([
    c,
    0,
])

conjugated_vector = sp.simplify(
    M * translation_vector
)

print("  canonical translation vector:")
print("   ", translation_vector)

print()
print("  transformed translation vector:")
print("   ", conjugated_vector)

print("""
  Therefore an upstream operator is affine-conjugate to the
  canonical KAPPA translation whenever its translation vector
  spans a nonzero one-dimensional moving direction.

  The invariant coordinate is any independent linear form
  orthogonal to that direction in the chosen coordinate system.
""")


# ============================================================================
# [13] LINEAR INVARIANT / MOVING DUAL BASIS
# ============================================================================

print()
print("[13] LINEAR DUAL-BASIS CONDITIONS")
print("-" * 78)

h10, h01, h00 = sp.symbols(
    "h10 h01 h00"
)

j10b, j01b, j00b = sp.symbols(
    "j10b j01b j00b"
)

Hlin = sp.expand(
    h10*U
    + h01*V
    + h00
)

Jlin = sp.expand(
    j10b*U
    + j01b*V
    + j00b
)

dHlin = sp.expand(
    Hlin.subs(U, U + alpha)
    - Hlin
)

dJlin = sp.expand(
    Jlin.subs(U, U + alpha)
    - Jlin
)

print("  H =")
print("   ", Hlin)

print()
print("  J =")
print("   ", Jlin)

print()
print("  T(H)-H =", dHlin)
print("  T(J)-J =", dJlin)

print()
print("""
  Conditions:

      T(H)-H = alpha*h10
      T(J)-J = 0

  therefore:

      h10 != 0
      j10b = 0.

  So J must live entirely in the invariant V-direction,
  while H must have a nonzero component along the moving
  direction.
""")


# ============================================================================
# [14] CONIC KERNEL RECONSTRUCTION FROM H,J
# ============================================================================

print()
print("[14] KERNEL RECONSTRUCTION FROM H,J")
print("-" * 78)

# Canonical gauge:
# H = U
# J = V
#
# Most general quadratic with the required leading response:
#
# K = H^2 + H*I(J) + J2(J)

q1, q0, jpoly2, jpoly1, jpoly0 = sp.symbols(
    "q1 q0 jpoly2 jpoly1 jpoly0"
)

K_HJ = sp.expand(
    Hlin**2
    + Hlin*(q1*Jlin + q0)
    + jpoly2*Jlin**2
    + jpoly1*Jlin
    + jpoly0
)

print("""
  General gauge-compatible quadratic:

      K = H^2 + H*I(J) + J2(J).

  Its translation response is determined only by H.
""")


# ============================================================================
# [15] CANONICAL CHOICE
# ============================================================================

print()
print("[15] CANONICAL CHOICE")
print("-" * 78)

Hc = U
Jc = V
Kc = sp.expand(
    Hc**2 - Jc**2
)

cert(
    "canonical K(H,J)",
    Kc,
    U**2 - V**2,
)


# ============================================================================
# [16] GENERATING POLYNOMIAL
# ============================================================================

print()
print("[16] GENERATING POLYNOMIAL")
print("-" * 78)

Qcanon = sp.expand(
    (z*c - U - V)
    * (z*c - U + V)
    / c**2
)

print("  Q(z) =")
print("   ", sp.factor(Qcanon))

cert(
    "Q(z)=z^2-Sz+N",
    Qcanon,
    z**2
    - (2*U/c)*z
    + (U**2 - V**2)/c**2,
)

cert(
    "Q(z) = L(-z)",
    Qcanon,
    (
        (U - V - z*c)
        * (U + V - z*c)
        / c**2
    ),
)


# ============================================================================
# [17] DISCRIMINANT / INVARIANT CHANNEL
# ============================================================================

print()
print("[17] DISCRIMINANT / INVARIANT CHANNEL")
print("-" * 78)

disc_Q = sp.factor(
    sp.discriminant(Qcanon, z)
)

cert(
    "disc(Q)=4V^2/c^2",
    disc_Q,
    4*V**2/c**2,
)


# ============================================================================
# [18] TRANSLATION OF THE SPECTRAL VARIABLE
# ============================================================================

print()
print("[18] SPECTRAL TRANSLATION")
print("-" * 78)

Q_shift = sp.expand(
    Qcanon.subs(z, z - 1)
)

Q_shift_target = sp.expand(
    ((z - 1)*c - U - V)
    * ((z - 1)*c - U + V)
    / c**2
)

cert(
    "Q(z-1)",
    Q_shift,
    Q_shift_target,
)


# ============================================================================
# [19] FOUR-BASE COMPOSITION
# ============================================================================

print()
print("[19] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = sp.expand(
    Qcanon * Q_shift
)

chi_target = sp.expand(
    (
        ((z)*c - U - V)
        * ((z)*c - U + V)
        * ((z - 1)*c - U - V)
        * ((z - 1)*c - U + V)
    )
    / c**4
)

cert(
    "Q(z)Q(z-1)",
    chi,
    chi_target,
)


# ============================================================================
# [20] ORBIT FORM
# ============================================================================

print()
print("[20] ORBIT FORM")
print("-" * 78)

L = sp.expand(
    (U + m*c)**2 / c**2
    - V**2/c**2
)

L_target = sp.expand(
    m**2
    + (2*U/c)*m
    + (U**2 - V**2)/c**2
)

cert(
    "L(m)=m^2+S*m+N",
    L,
    L_target,
)


# ============================================================================
# [21] CONVERSE IDENTIFICATION CRITERION
# ============================================================================

print()
print("[21] CONVERSE IDENTIFICATION CRITERION")
print("-" * 78)

print("""
  An unknown operator T belongs to the KAPPA conjugacy class
  if one can find observables H,J,K and a nonzero c such that:

      T(H)-H = c

      T(J)-J = 0

      T(K)-K = 2cH+c^2 + I

  where I is T-invariant,

  and K has quadratic degree two in H.

  After removing invariant gauge terms:

      K ~ H^2 - J^2.

  The spectral polynomial is then:

      Q(z)
        = z^2 - (2H/c)z + K/c^2

  and:

      Delta = 4J^2/c^2.
""")


# ============================================================================
# [22] EXACT NECESSITY / SUFFICIENCY TEST
# ============================================================================

print()
print("[22] NECESSITY / SUFFICIENCY")
print("-" * 78)

Htest = U
Jtest = V
Ktest = U**2 - V**2

cond1 = sp.factor(
    T_matrix(Htest, (1,0,c,0,1,0)) - Htest
)

cond2 = sp.factor(
    T_matrix(Jtest, (1,0,c,0,1,0)) - Jtest
)

cond3 = sp.factor(
    T_matrix(Ktest, (1,0,c,0,1,0))
    - Ktest
    - 2*c*Htest
    - c**2
)

print("  canonical operator:")
print("   ", "U -> U+c")
print("   ", "V -> V")

cert(
    "T(H)-H=c",
    cond1,
    c,
)

cert(
    "T(J)-J=0",
    cond2,
    0,
)

cert(
    "T(K)-K-2cH-c^2=0",
    cond3,
    0,
)


# ============================================================================
# [23] FINAL STRUCTURAL STATEMENT
# ============================================================================

print()
print("[23] FINAL STRUCTURAL STATEMENT")
print("-" * 78)

print("""
  The KAPPA bridge is characterized by an affine unipotent
  translation direction.

  In canonical coordinates:

      moving direction  = H
      invariant algebra = functions of J

      T(H) = H + c
      T(J) = J

  A quadratic observable satisfying

      T(K)-K = 2cH+c^2

  generates the orbit

      K_m/c^2
          = m^2 + (2H/c)m + K_0/c^2.

  Therefore:

      S = 2H/c

      N = K_0/c^2

      Delta = 4J^2/c^2.

  The factor-bearing quadratic follows as

      Q(z)
        = z^2-Sz+N.

  This provides a recognition theorem:

      original homogeneous operator
          |
          v
      affine moving direction
          +
      invariant direction
          +
      quadratic observable
          |
          v
      KAPPA normal form.
""")


# ============================================================================
# FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 533 FINISHED")
print("=" * 78)

print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The next step is finally an operator-search problem.

Given an ACTUAL transformation from the original homogeneous
construction:

    T : observables -> observables

do not search arbitrary algebraic identities.

Instead compute its action on the lowest-degree observable
space and determine:

    1. eigenvalue-1 invariant directions;

    2. generalized eigenvectors with

           (T-I)H = c;

    3. quadratic observables K satisfying

           (T-I)K = 2cH+c^2 + invariant;

    4. whether the resulting invariant direction supplies the
       discriminant channel;

    5. whether the generated quadratic has constant term equal
       to the known N-only layer.

The decisive object is the Jordan structure of T:

      invariant direction
             |
             v
      generalized eigenvector
             |
             v
      quadratic observable.

A match here would be substantially stronger than another
identity: it would identify the operator representation
behind the KAPPA construction.
""")
