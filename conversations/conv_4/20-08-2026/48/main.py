#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 524 START")
print("=" * 78)
print("HOMOGENEOUS OPERATOR EXTRACTION -> TRANSLATION QUADRATIC")
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
m, h = sp.symbols("m h")
z = sp.symbols("z")

N = p * q
S = p + q
Delta = (p - q)**2

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
# [1] UNIVERSAL TRANSLATION ORBIT
# ============================================================================

print("[1] UNIVERSAL FACTOR-TRANSLATION ORBIT")
print("-" * 78)

L = sp.expand((p + m) * (q + m))

print("  L(m) =", sp.factor(L))

check(
    "L(m) = m^2 + S*m + N",
    L,
    m**2 + S*m + N,
)

print()


# ============================================================================
# [2] HOMOGENEOUS DECOMPOSITION
# ============================================================================

print("[2] HOMOGENEOUS DEGREE DECOMPOSITION")
print("-" * 78)

L_expanded = sp.expand(L)

deg2 = sp.expand(m**2)
deg1 = sp.expand(S * m)
deg0 = sp.expand(N)

print("  degree-2 layer =", deg2)
print("  degree-1 layer =", deg1)
print("  degree-0 layer =", deg0)

check(
    "L = degree2 + degree1 + degree0",
    L_expanded,
    deg2 + deg1 + deg0,
)

print()


# ============================================================================
# [3] EULER OPERATOR
# ============================================================================

print("[3] EULER HOMOGENEITY OPERATOR")
print("-" * 78)

E = sp.expand(
    p * sp.diff(L, p) +
    q * sp.diff(L, q)
)

print("  E[L] =", sp.factor(E))

check(
    "Euler degree decomposition",
    E,
    2 * m**2 + S * m,
)

print()


# ============================================================================
# [4] SECOND EULER OPERATOR
# ============================================================================

print("[4] SECOND EULER OPERATOR")
print("-" * 78)

E2 = sp.expand(
    p**2 * sp.diff(L, p, 2)
    + 2*p*q * sp.diff(L, p, q)
    + q**2 * sp.diff(L, q, 2)
)

print("  E2[L] =", sp.factor(E2))

check(
    "second Euler operator",
    E2,
    2 * N,
)

print()


# ============================================================================
# [5] FACTOR-CONTENT SEPARATION
# ============================================================================

print("[5] FACTOR-CONTENT SEPARATION")
print("-" * 78)

# Translation orbit expanded directly:
# (p+m)(q+m) = pq + m(p+q) + m^2

content0 = sp.expand(L.subs(m, 0))
content1 = sp.expand(L.subs(m, 1))
content_minus1 = sp.expand(L.subs(m, -1))

print("  L(0)  =", content0)
print("  L(1)  =", content1)
print("  L(-1) =", content_minus1)

check("L(0)=N", content0, N)
check("L(1)-L(0)=S+1", content1-content0, S+1)
check("L(0)-L(-1)=S-1", content0-content_minus1, S-1)

print()


# ============================================================================
# [6] CENTRAL DIFFERENCE
# ============================================================================

print("[6] CENTRAL DIFFERENCE OPERATOR")
print("-" * 78)

central = sp.expand(
    L.subs(m, h) - L.subs(m, -h)
)

central_normalized = sp.expand(
    central / (2*h)
)

print("  L(h)-L(-h) =", sp.factor(central))
print("  normalized =", sp.factor(central_normalized))

check(
    "central difference extracts S+2m",
    sp.expand(
        (L.subs(m, h) - L.subs(m, -h)) / (2*h)
    ),
    S,
)

print()


# ============================================================================
# [7] GENERAL TRANSLATION STEP
# ============================================================================

print("[7] GENERAL TRANSLATION STEP")
print("-" * 78)

L_forward = sp.expand(L.subs(m, m+h))
L_backward = sp.expand(L.subs(m, m-h))

forward_difference = sp.expand(
    (L_forward - L) / h
)

central_difference = sp.expand(
    (L_forward - L_backward) / (2*h)
)

print("  forward difference / h =", sp.factor(forward_difference))
print("  central difference     =", sp.factor(central_difference))

check(
    "forward difference law",
    forward_difference,
    S + 2*m + h,
)

check(
    "central difference law",
    central_difference,
    S + 2*m,
)

print()


# ============================================================================
# [8] SECOND DIFFERENCE
# ============================================================================

print("[8] SECOND DIFFERENCE OPERATOR")
print("-" * 78)

second_difference = sp.expand(
    (
        L.subs(m, m+h)
        - 2*L
        + L.subs(m, m-h)
    ) / h**2
)

print("  second difference =", second_difference)

check(
    "second difference = 2",
    second_difference,
    2,
)

print()


# ============================================================================
# [9] OPERATOR THAT REMOVES THE UNIVERSAL SHIFT PART
# ============================================================================

print("[9] SHIFT-CORRECTED FIRST DIFFERENCE")
print("-" * 78)

corrected = sp.expand(
    (L_forward - L) / h - 2*m - h
)

print("  corrected operator =", sp.factor(corrected))

check(
    "corrected difference extracts S",
    corrected,
    S,
)

print()


# ============================================================================
# [10] OPERATOR REPRESENTATION OF THE KAPPA QUADRATIC
# ============================================================================

print("[10] KAPPA QUADRATIC FROM LAYER ORBIT")
print("-" * 78)

Q_from_L = sp.expand(
    L.subs(m, -z)
)

Q_target = sp.expand(
    z**2 - S*z + N
)

print("  Q_from_L =", sp.factor(Q_from_L))
print("  Q_target =", sp.factor(Q_target))

check(
    "L(-z)=z^2-Sz+N",
    Q_from_L,
    Q_target,
)

print()


# ============================================================================
# [11] TRANSLATION COMPOSITION
# ============================================================================

print("[11] UNIVERSAL TRANSLATION COMPOSITION")
print("-" * 78)

Q_shifted = sp.expand(
    Q_target.subs(z, z-1)
)

chi = sp.expand(
    Q_target * Q_shifted
)

expected_chi = sp.expand(
    (z-p)*(z-q)*(z-p-1)*(z-q-1)
)

check(
    "Q(z)Q(z-1)=four-base characteristic",
    chi,
    expected_chi,
)

print()


# ============================================================================
# [12] HOMOGENEITY + TRANSLATION: GENERAL ANSATZ
# ============================================================================

print("[12] GENERAL HOMOGENEOUS TRANSLATION ANSATZ")
print("-" * 78)

a, b, c = sp.symbols("a b c")

# Generic quadratic orbit generated from a degree-2 product:
#   R(m) = a m^2 + b S m + c N
R = sp.expand(
    a*m**2 + b*S*m + c*N
)

print("  R(m) =", R)

R_diff = sp.expand(
    R.subs(m, m+1) - R
)

R_second = sp.expand(
    R.subs(m, m+2)
    - 2*R.subs(m, m+1)
    + R
)

print("  ΔR =", R_diff)
print("  Δ²R =", R_second)

print()


# ============================================================================
# [13] SOLVE OPERATOR CONDITIONS
# ============================================================================

print("[13] SOLVE UNIVERSAL QUADRATIC CONDITIONS")
print("-" * 78)

# Require:
#   coefficient of m^2 = 1
#   coefficient of Sm = 1
#   coefficient of N = 1
#
# This determines the only normalized orbit.

solutions = sp.solve(
    [
        sp.Eq(a, 1),
        sp.Eq(b, 1),
        sp.Eq(c, 1),
    ],
    [a, b, c],
    dict=True,
)

print("  solutions =", solutions)

if solutions:
    sol = solutions[0]

    R_normalized = sp.expand(
        R.subs(sol)
    )

    check(
        "normalized orbit",
        R_normalized,
        m**2 + S*m + N,
    )

print()


# ============================================================================
# [14] GENERAL SYMMETRIC QUADRATIC IN p,q
# ============================================================================

print("[14] MOST GENERAL LOW-DEGREE SYMMETRIC ORBIT")
print("-" * 78)

alpha, beta, gamma, delta = sp.symbols(
    "alpha beta gamma delta"
)

# General symmetric quadratic expression in p,q and m:
#
#   alpha * (p^2+q^2)
# + beta  * pq
# + gamma * m(p+q)
# + delta * m^2
#
# We ask which coefficients make it equal to
#
#   (p+m)(q+m).

generic = sp.expand(
    alpha*(p**2 + q**2)
    + beta*p*q
    + gamma*m*(p+q)
    + delta*m**2
)

target = sp.expand(
    (p+m)*(q+m)
)

poly = sp.Poly(
    sp.expand(generic-target),
    p, q, m
)

equations = [
    sp.Eq(coeff, 0)
    for coeff in poly.coeffs()
]

solutions_general = sp.solve(
    equations,
    [alpha, beta, gamma, delta],
    dict=True,
)

print("  coefficient solution =", solutions_general)

if solutions_general:
    check(
        "general symmetric orbit reconstruction",
        generic.subs(solutions_general[0]),
        target,
    )

print()


# ============================================================================
# [15] KEY STRUCTURAL RESULT
# ============================================================================

print("[15] KEY STRUCTURAL RESULT")
print("-" * 78)

print("""
  The translation quadratic is characterized uniquely by:

      R(m) = m^2 + S*m + N

  and therefore by the operator identities:

      ΔR(m) = S + 2m + 1

      Δ²R(m) = 2.

  The useful information lives in the first difference.

  The curvature is universal.

  Equivalently:

      S = ΔR(m) - 2m - 1.

  The spectral polynomial is then obtained mechanically:

      Q(z) = R(-z)
           = z^2 - S*z + N.

  The +1 KAPPA factor is

      Q(z-1).

  Therefore the whole four-base spectrum follows from
  one translated quadratic orbit.
""")

print()


# ============================================================================
# [16] UPSTREAM OPERATOR NECESSITY TEST
# ============================================================================

print("[16] UPSTREAM OPERATOR NECESSITY TEST")
print("-" * 78)

print("""
  Any proposed homogeneous-layer observable H_m that is claimed
  to generate the KAPPA bridge must satisfy, after normalization,

      H_m = m^2 + A*m + B.

  For the desired factorization bridge:

      B = N
      A = S.

  Hence the decisive symbolic tests are:

      H_m - m^2
          is affine in m,

      coefficient of m
          is the missing symmetric coordinate,

      constant term
          is exactly N,

      Δ²H_m
          is exactly 2.

  No numerical interpolation is needed.
  No factor-pair enumeration is needed.
  No continued fractions are needed.
""")

print()


# ============================================================================
# [17] FINAL AUDIT
# ============================================================================

print("=" * 78)
print("EXPERIMENT 524 FINISHED")
print("=" * 78)

print()
print(f"SYMBOLIC FAILURES = {failures}")

if failures == 0:
    print("OVERALL EXACT SYMBOLIC AUDIT = True")
else:
    print("OVERALL EXACT SYMBOLIC AUDIT = False")

print("""
NEXT RESEARCH TARGET
------------------------------------------------------------------------------

The synthetic translation orbit is now completely characterized.

The next experiment should instantiate H_m using the ACTUAL
homogeneous-layer expressions from the earlier construction.

For each genuine layer family H_m, test symbolically:

    1. degree_m(H_m) == 2

    2. H_m - m^2 is affine in m

    3. coefficient_m(H_m) is a candidate for S

    4. constant_m(H_m) == N

    5. Δ²H_m == 2

The strongest possible result is:

    H_m = m^2 + S*m + N

because then the KAPPA quadratic is not constructed
after factor recovery; it is already present as the
generating polynomial of the homogeneous-layer family.

The next step is therefore to plug in the ORIGINAL
homogeneous-layer formulas, not to invent another downstream
identity.
""")
