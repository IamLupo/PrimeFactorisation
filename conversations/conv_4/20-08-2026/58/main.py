#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 534 START")
print("=" * 78)
print("AFFINE JORDAN OPERATOR -> INVARIANT / MOVING BASIS -> KAPPA NORMAL FORM")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

x, y = sp.symbols("x y")
X, Y = sp.symbols("X Y")

a, b, c, d = sp.symbols("a b c d")
alpha, beta = sp.symbols("alpha beta")

h, z, m = sp.symbols("h z m")
step = sp.symbols("step", nonzero=True)

# General linear observables
l1, l2, l0 = sp.symbols("l1 l2 l0")
j1, j2, j0 = sp.symbols("j1 j2 j0")

# Quadratic coefficients
k20, k11, k02 = sp.symbols("k20 k11 k02")
k10, k01, k00 = sp.symbols("k10 k01 k00")

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
# GENERAL AFFINE OPERATOR
# ============================================================================

def T_expr(expr):
    """
    General affine operator:

        x' = a*x + b*y + alpha
        y' = c*x + d*y + beta
    """
    return sp.expand(
        expr.subs({
            x: a*x + b*y + alpha,
            y: c*x + d*y + beta,
        })
    )


print()
print("[1] GENERAL AFFINE OPERATOR")
print("-" * 78)

print("""
  x' = a*x + b*y + alpha
  y' = c*x + d*y + beta

  Linear part:

      M = [[a,b],
           [c,d]]

  Translation vector:

      t = [alpha,beta]^T
""")


# ============================================================================
# LINEAR OPERATOR MATRIX
# ============================================================================

print()
print("[2] LINEAR JORDAN DATA")
print("-" * 78)

M = sp.Matrix([
    [a, b],
    [c, d],
])

I2 = sp.eye(2)

trace_M = sp.expand(sp.trace(M))
det_M = sp.expand(M.det())

char_M = sp.expand(
    sp.det(sp.Symbol("lambda") * I2 - M)
)

print("  trace(M) =", trace_M)
print("  det(M)   =", det_M)
print("  characteristic polynomial:")
print("   ", sp.factor(char_M))

print()
print("""
  The KAPPA translation operator has linear matrix

      I.

  Its nontrivial behavior is therefore entirely in the
  affine translation vector.

  A genuinely affine-conjugate translation requires:

      M = I

  after choosing the appropriate observable coordinates.
""")


# ============================================================================
# UNIPOTENT / IDENTITY CONDITIONS
# ============================================================================

print()
print("[3] IDENTITY-LINEAR-PART TEST")
print("-" * 78)

identity_conditions = [
    sp.Eq(a, 1),
    sp.Eq(b, 0),
    sp.Eq(c, 0),
    sp.Eq(d, 1),
]

for eq in identity_conditions:
    print("   ", eq)

print()
print("""
  More generally, an upstream realization may first require
  conjugating M into the identity action on the relevant
  observable quotient.

  This experiment therefore separates:

      linear-action obstruction
      from
      affine-translation structure.
""")


# ============================================================================
# GENERAL LINEAR OBSERVABLE
# ============================================================================

print()
print("[4] ACTION ON A GENERAL LINEAR OBSERVABLE")
print("-" * 78)

H = sp.expand(
    l1*x + l2*y + l0
)

J = sp.expand(
    j1*x + j2*y + j0
)

TH_minus_H = sp.expand(
    T_expr(H) - H
)

TJ_minus_J = sp.expand(
    T_expr(J) - J
)

print("  H =", H)
print()
print("  T(H)-H =")
print("   ", sp.factor(TH_minus_H))

print()
print("  J =", J)
print()
print("  T(J)-J =")
print("   ", sp.factor(TJ_minus_J))


# ============================================================================
# LINEAR OBSERVABLE COEFFICIENT SYSTEM
# ============================================================================

print()
print("[5] INVARIANT LINEAR OBSERVABLE EQUATIONS")
print("-" * 78)

poly_H = sp.Poly(
    sp.expand(TH_minus_H),
    x,
    y,
)

poly_J = sp.Poly(
    sp.expand(TJ_minus_J),
    x,
    y,
)

eq_H_zero = [
    sp.Eq(coeff, 0)
    for _, coeff in poly_H.terms()
]

eq_J_zero = [
    sp.Eq(coeff, 0)
    for _, coeff in poly_J.terms()
]

print("  J-invariance requires:")
for eq in eq_J_zero:
    print("   ", eq)

print()
print("""
  A nonconstant invariant linear observable must satisfy
  the left-eigenvector condition

      [j1,j2] M = [j1,j2].

  A moving coordinate with constant translation response
  must satisfy

      [l1,l2] M = [l1,l2],

  while its affine part is shifted by the translation vector.
""")


# ============================================================================
# IDENTITY-LINEAR-PART SPECIALIZATION
# ============================================================================

print()
print("[6] IDENTITY-LINEAR-PART SPECIALIZATION")
print("-" * 78)

T_identity = lambda expr: sp.expand(
    expr.subs({
        x: x + alpha,
        y: y + beta,
    })
)

H_id = sp.expand(
    l1*x + l2*y + l0
)

J_id = sp.expand(
    j1*x + j2*y + j0
)

dH_id = sp.expand(
    T_identity(H_id) - H_id
)

dJ_id = sp.expand(
    T_identity(J_id) - J_id
)

print("  T(H)-H =", dH_id)
print("  T(J)-J =", dJ_id)


# ============================================================================
# TRANSLATION DIRECTION
# ============================================================================

print()
print("[7] MOVING / INVARIANT DUAL BASIS")
print("-" * 78)

print("""
  Translation vector:

      t = (alpha,beta).

  A linear invariant J must obey:

      j1*alpha + j2*beta = 0.

  A moving coordinate H obeys:

      h1*alpha + h2*beta != 0.

  Thus the observable plane splits into:

      invariant direction  perpendicular to t
      moving direction     parallel to t in the dual space.
""")


# ============================================================================
# EXPLICIT DUAL BASIS
# ============================================================================

print()
print("[8] CANONICAL DUAL BASIS FOR NONZERO TRANSLATION")
print("-" * 78)

# Avoid dividing by alpha/beta by using polynomial expressions.
H_canonical = sp.expand(
    alpha*x + beta*y
)

J_canonical = sp.expand(
    -beta*x + alpha*y
)

THc = sp.expand(
    T_identity(H_canonical)
    - H_canonical
)

TJc = sp.expand(
    T_identity(J_canonical)
    - J_canonical
)

print("  H =", H_canonical)
print("  J =", J_canonical)

print()
print("  T(H)-H =", sp.factor(THc))
print("  T(J)-J =", sp.factor(TJc))

cert(
    "canonical moving coordinate response",
    THc,
    alpha**2 + beta**2,
)

cert(
    "canonical invariant coordinate response",
    TJc,
    0,
)


# ============================================================================
# QUADRATIC OBSERVABLE
# ============================================================================

print()
print("[9] GENERAL QUADRATIC OBSERVABLE")
print("-" * 78)

K = sp.expand(
    k20*x**2
    + k11*x*y
    + k02*y**2
    + k10*x
    + k01*y
    + k00
)

dK = sp.expand(
    T_identity(K) - K
)

print("  K =")
print("   ", K)

print()
print("  T(K)-K =")
print("   ", sp.factor(dK))


# ============================================================================
# KAPPA-TYPE QUADRATIC CLASS
# ============================================================================

print()
print("[10] KAPPA-TYPE QUADRATIC CLASS")
print("-" * 78)

# Build an invariant quadratic using J_canonical.
#
# The canonical moving coordinate is H_canonical.
# The most general quadratic with pure moving response is
#
#     K = A H^2 + H*I(J) + J2(J)
#
# Here we test the basic canonical representative.

K_basic = sp.expand(
    H_canonical**2
)

dK_basic = sp.expand(
    T_identity(K_basic) - K_basic
)

print("  K_basic = H^2")
print("  T(K_basic)-K_basic =")
print("   ", sp.factor(dK_basic))

cert(
    "quadratic translation response",
    dK_basic,
    2*(alpha**2 + beta**2)*H_canonical
    + (alpha**2 + beta**2)**2,
)


# ============================================================================
# CONIC-LIKE NORMALIZATION
# ============================================================================

print()
print("[11] INVARIANT SUBTRACTION")
print("-" * 78)

# The natural quadratic invariant of the translation direction is
# J_canonical^2.

K_conic = sp.expand(
    H_canonical**2
    - J_canonical**2
)

dK_conic = sp.expand(
    T_identity(K_conic) - K_conic
)

print("  K_conic = H^2-J^2")
print()
print("  T(K_conic)-K_conic =")
print("   ", sp.factor(dK_conic))

cert(
    "pure moving response after invariant subtraction",
    dK_conic,
    2*(alpha**2 + beta**2)*H_canonical
    + (alpha**2 + beta**2)**2,
)


# ============================================================================
# SECOND DIFFERENCE
# ============================================================================

print()
print("[12] SECOND OPERATOR DIFFERENCE")
print("-" * 78)

def T_integer(expr, n):
    return sp.expand(
        expr.subs({
            x: x + n*alpha,
            y: y + n*beta,
        })
    )


K0 = K_conic
K1 = T_integer(K_conic, 1)
K2 = T_integer(K_conic, 2)

second_difference = sp.expand(
    K2 - 2*K1 + K0
)

print("  T^2(K)-2T(K)+K =")
print("   ", sp.factor(second_difference))

cert(
    "second difference constant",
    second_difference,
    2*(alpha**2 + beta**2),
)


# ============================================================================
# NORMALIZED ORBIT
# ============================================================================

print()
print("[13] NORMALIZED ORBIT")
print("-" * 78)

rho2 = sp.expand(
    alpha**2 + beta**2
)

L_n = sp.expand(
    T_integer(K_conic, m) / rho2
)

print("  normalized orbit =")
print("   ", sp.factor(L_n))

# Rewrite expected form using H_canonical/rho2.
S_abstract = sp.expand(
    2*H_canonical/rho2
)

N_abstract = sp.expand(
    K_conic/rho2
)

L_expected = sp.expand(
    m**2
    + S_abstract*m
    + N_abstract
)

cert(
    "quadratic orbit structure",
    L_n,
    L_expected,
)


# ============================================================================
# GENERATING QUADRATIC
# ============================================================================

print()
print("[14] GENERATING QUADRATIC")
print("-" * 78)

Q_operator = sp.expand(
    L_expected.subs(m, -z)
)

print("  Q(z) = L(-z)")
print("   ", sp.factor(Q_operator))

Q_expected = sp.expand(
    z**2
    - S_abstract*z
    + N_abstract
)

cert(
    "L(-z)=z^2-Sz+N",
    Q_operator,
    Q_expected,
)


# ============================================================================
# DISCRIMINANT
# ============================================================================

print()
print("[15] DISCRIMINANT CHANNEL")
print("-" * 78)

disc_operator = sp.factor(
    sp.discriminant(Q_operator, z)
)

print("  discriminant(Q) =")
print("   ", disc_operator)

print()
print("""
  In the canonical historical representation:

      H = U
      J = V

  and:

      Delta = 4 V^2/c^2.

  Here the same structure appears abstractly as the invariant
  quadratic coordinate J^2.
""")


# ============================================================================
# GAUGE TEST
# ============================================================================

print()
print("[16] INVARIANT GAUGE TEST")
print("-" * 78)

g1, g0 = sp.symbols("g1 g0")

I_linear = sp.expand(
    g1*J_canonical + g0
)

H_gauged = sp.expand(
    H_canonical + I_linear
)

dH_gauged = sp.expand(
    T_identity(H_gauged) - H_gauged
)

cert(
    "H -> H + invariant has same response",
    dH_gauged,
    alpha**2 + beta**2,
)


# ============================================================================
# OPERATOR-LEVEL RECOGNITION CONDITIONS
# ============================================================================

print()
print("[17] OPERATOR RECOGNITION CONDITIONS")
print("-" * 78)

print("""
  For a genuine upstream homogeneous operator, seek:

      T(J) = J

      T(H) = H + c

      T(K) = K + 2*c*H + c^2 + I(J)

  where I(J) is invariant.

  Then:

      T^2(K)-2T(K)+K

  must be constant.

  The corresponding normalized orbit has the form

      L_m = m^2 + S*m + N

  after a suitable scaling and invariant gauge.

  The important distinction is:

      the operator supplies the Jordan direction;

      the quadratic supplies the factor-bearing invariant.
""")


# ============================================================================
# CONVERSE THEOREM IN SYMBOLIC FORM
# ============================================================================

print()
print("[18] CONVERSE NORMAL-FORM CERTIFICATE")
print("-" * 78)

# Use the abstract H,J,c directly.
c0 = sp.symbols("c0", nonzero=True)
H0, J0 = sp.symbols("H0 J0")

# Canonical abstract orbit.
K_abstract = sp.expand(
    H0**2 - J0**2
)

K_abstract_shift = sp.expand(
    (H0 + c0)**2 - J0**2
)

abstract_increment = sp.expand(
    K_abstract_shift
    - K_abstract
)

print("  K = H^2-J^2")
print("  T(K)-K =")
print("   ", abstract_increment)

cert(
    "abstract KAPPA increment",
    abstract_increment,
    2*c0*H0 + c0**2,
)


# ============================================================================
# JORDAN STRUCTURE SUMMARY
# ============================================================================

print()
print("[19] JORDAN-STRUCTURE SUMMARY")
print("-" * 78)

print("""
  The affine operator has a useful KAPPA representation precisely
  when its observable space contains:

      1. an invariant linear form J;

      2. a generalized translation coordinate H;

      3. a quadratic K with the same moving direction.

  In operator language:

      (T-I)J = 0

      (T-I)H = c

      (T-I)K = 2cH + c^2 + invariant.

  The quadratic orbit then integrates the linear Jordan channel:

      H
      |
      v
      K
      |
      v
      L_m = m^2 + S*m + N.

  This is the operator-level structure that should be searched
  for in the actual homogeneous-layer construction.
""")


# ============================================================================
# FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 534 FINISHED")
print("=" * 78)

print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The canonical KAPPA form is now characterized without assuming
U,V as primitive coordinates.

The next genuine upstream test is:

    take each ACTUAL homogeneous-layer operator T

    restrict T to the lowest-degree observable space

    compute:

        (T-I) on linear observables

        (T-I) on quadratic observables

    and determine its Jordan structure.

Specifically search for:

    ker(T-I)                 -> invariant observables

    generalized kernel       -> moving observable H

    quadratic lift           -> K

such that:

    (T-I)H = c

    (T-I)K = 2cH+c^2+I

and:

    (T-I)^2 K = 2c^2.

The decisive result would be an EXISTING homogeneous-layer
operator whose Jordan chain has exactly this structure.

That would establish an operator conjugacy to the KAPPA
translation system before introducing p, q, S, Delta, or
continued fractions.
""")
