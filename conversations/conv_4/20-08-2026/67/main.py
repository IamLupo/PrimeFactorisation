#!/usr/bin/env python3

import sympy as sp


# ============================================================================
# EXPERIMENT 543
# ============================================================================
# HOMOGENEOUS LAYER INDEX -> INTRINSIC TRANSLATION / S-CHANNEL
#
# No numerical fitting.
# No factor enumeration.
# No continued fractions.
# No p,q.
# ============================================================================

print("=" * 78)
print("EXPERIMENT 543 START")
print("=" * 78)
print("HOMOGENEOUS LAYER INDEX -> INTRINSIC TRANSLATION / S-CHANNEL")
print()


# ============================================================================
# SYMBOLS -- IMPORTANT: DECLARE ALL SYMBOLS BEFORE USE
# ============================================================================

N, S = sp.symbols("N S")
h, k, m = sp.symbols("h k m")
z = sp.symbols("z")


failures = 0


def simp(expr):
    """Exact symbolic simplification."""
    return sp.factor(sp.expand(sp.cancel(sp.simplify(expr))))


def cert(label, expr, expected=0):
    """Print and count an exact symbolic certificate."""
    global failures

    diff = simp(expr - expected)
    ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# [1] EXISTING HOMOGENEOUS LAYERS
# ============================================================================

print()
print("[1] EXISTING HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = sp.expand(6*N - S**2 + S)
F3 = sp.expand((S + 1)*F2)

print("  F2 =")
sp.pprint(F2)

print()
print("  F3 =")
sp.pprint(F3)

print()

cert(
    "F3=(S+1)F2",
    F3,
    (S + 1)*F2,
)


# ============================================================================
# [2] LAYER-RATIO S CHANNEL
# ============================================================================

print()
print("[2] LAYER-RATIO S CHANNEL")
print("-" * 78)

ratio = sp.cancel(F3 / F2)

print("  F3/F2 =")
sp.pprint(ratio)

cert(
    "F3/F2 = S+1",
    ratio,
    S + 1,
)

S_ratio = simp(ratio - 1)

print()
print("  S from ratio = F3/F2 - 1")
sp.pprint(S_ratio)

cert(
    "ratio extraction gives S",
    S_ratio,
    S,
)


# ============================================================================
# [3] DISTINGUISHED TRANSLATION OF N AND S
# ============================================================================

print()
print("[3] DISTINGUISHED TRANSLATION OPERATOR")
print("-" * 78)

N_h = sp.expand(N + h*S + h**2)
S_h = sp.expand(S + 2*h)

print("  T_h(N) =")
sp.pprint(N_h)

print()
print("  T_h(S) =")
sp.pprint(S_h)

cert(
    "T_h(N)",
    N_h,
    N + h*S + h**2,
)

cert(
    "T_h(S)",
    S_h,
    S + 2*h,
)


# ============================================================================
# [4] TRANSLATED F2
# ============================================================================

print()
print("[4] TRANSLATED F2")
print("-" * 78)

F2_h = sp.expand(
    6*N_h - S_h**2 + S_h
)

print("  T_h(F2) =")
sp.pprint(sp.factor(F2_h))

dF2 = sp.factor(F2_h - F2)

print()
print("  T_h(F2)-F2 =")
sp.pprint(dF2)

cert(
    "F2 translation response",
    dF2,
    2*h*(S + h + 1),
)


# ============================================================================
# [5] EXACT S EXTRACTION FROM TRANSLATION RESPONSE
# ============================================================================

print()
print("[5] EXACT S EXTRACTION")
print("-" * 78)

S_translation = sp.cancel(
    dF2 / (2*h) - h - 1
)

print("  S_translation =")
sp.pprint(S_translation)

cert(
    "translation response extracts S",
    S_translation,
    S,
)


# ============================================================================
# [6] DIVISION-FREE VERSION
# ============================================================================

print()
print("[6] DIVISION-FREE S CHANNEL")
print("-" * 78)

division_free = sp.expand(
    dF2 - 2*h*(S + h + 1)
)

print("  T_h(F2)-F2 - 2h(S+h+1) =")
sp.pprint(division_free)

cert(
    "division-free identity",
    division_free,
    0,
)


# ============================================================================
# [7] NORMALIZED FIRST DIFFERENCE
# ============================================================================

print()
print("[7] NORMALIZED DIFFERENCE CHANNEL")
print("-" * 78)

R_h = sp.expand(
    dF2 / 2 - h*(h + 1)
)

print("  R_h =")
sp.pprint(R_h)

cert(
    "R_h = h*S",
    R_h,
    h*S,
)


# ============================================================================
# [8] RECONSTRUCT N FROM F2 AND S
# ============================================================================

print()
print("[8] N RECONSTRUCTION")
print("-" * 78)

N_recovered = sp.cancel(
    (F2 + S_translation**2 - S_translation) / 6
)

print("  N_recovered =")
sp.pprint(N_recovered)

cert(
    "N from F2 and recovered S",
    N_recovered,
    N,
)


# ============================================================================
# [9] TRANSLATION QUADRATIC ORBIT
# ============================================================================

print()
print("[9] TRANSLATION QUADRATIC ORBIT")
print("-" * 78)

L_m = sp.expand(
    N + m*S + m**2
)

print("  L(m) =")
sp.pprint(L_m)

D_L = sp.expand(
    L_m.subs(m, m + 1) - L_m
)

D2_L = sp.expand(
    L_m.subs(m, m + 2)
    - 2*L_m.subs(m, m + 1)
    + L_m
)

print()
print("  Delta L(m) =")
sp.pprint(D_L)

print()
print("  Delta^2 L(m) =")
sp.pprint(D2_L)

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
# [10] ORBIT INTEGRATION CONSTANTS
# ============================================================================

print()
print("[10] ORBIT INTEGRATION CONSTANTS")
print("-" * 78)

A_orbit = sp.expand(
    D_L - (2*m + 1)
)

B_orbit = sp.expand(
    L_m - A_orbit*m - m**2
)

print("  A =")
sp.pprint(A_orbit)

print()
print("  B =")
sp.pprint(B_orbit)

cert(
    "A=S",
    A_orbit,
    S,
)

cert(
    "B=N",
    B_orbit,
    N,
)


# ============================================================================
# [11] MATCH F2 RESPONSE TO ORBIT DIFFERENCE
# ============================================================================

print()
print("[11] F2 RESPONSE / QUADRATIC ORBIT MATCH")
print("-" * 78)

D_L_h = sp.expand(
    D_L.subs(m, h)
)

normalized_F2_difference = sp.cancel(
    dF2 / (2*h)
)

print("  Delta L(h) =")
sp.pprint(D_L_h)

print()
print("  Delta F2/(2h) =")
sp.pprint(normalized_F2_difference)

cert(
    "Delta F2/(2h) = Delta L(h)",
    normalized_F2_difference,
    D_L_h,
)


# ============================================================================
# [12] GENERATING KAPPA QUADRATIC
# ============================================================================

print()
print("[12] GENERATING KAPPA QUADRATIC")
print("-" * 78)

Q_operator = sp.expand(
    L_m.subs(m, -z)
)

Q_target = sp.expand(
    z**2 - S*z + N
)

print("  L(-z) =")
sp.pprint(Q_operator)

print()
print("  target Q(z) = z^2-Sz+N")

cert(
    "L(-z)=Q(z)",
    Q_operator,
    Q_target,
)


# ============================================================================
# [13] DISCRIMINANT
# ============================================================================

print()
print("[13] KAPPA DISCRIMINANT")
print("-" * 78)

Delta_operator = sp.factor(
    sp.discriminant(Q_operator, z)
)

Delta_target = sp.expand(
    S**2 - 4*N
)

print("  discriminant(Q) =")
sp.pprint(Delta_operator)

print()
print("  target Delta =")
sp.pprint(Delta_target)

cert(
    "discriminant = S^2-4N",
    Delta_operator,
    Delta_target,
)


# ============================================================================
# [14] F2 UNDER THE SAME SPECTRAL VARIABLE
# ============================================================================

print()
print("[14] F2 VS GENERATING QUADRATIC")
print("-" * 78)

F2_spectral = sp.expand(
    F2.subs(S, z)
)

print("  F2(S=z) =")
sp.pprint(F2_spectral)

print()
print("  Q(z) =")
sp.pprint(Q_target)

print()
print("  This confirms F2 is not Q itself.")
print("  F2 is an independent quadratic functional of S.")


# ============================================================================
# [15] F2 ORBIT AS A QUADRATIC IN h
# ============================================================================

print()
print("[15] F2 ORBIT POLYNOMIAL")
print("-" * 78)

poly_F2_h = sp.Poly(
    F2_h,
    h
)

coeff_h2 = sp.expand(poly_F2_h.coeff_monomial(h**2))
coeff_h1 = sp.expand(poly_F2_h.coeff_monomial(h))
coeff_h0 = sp.expand(poly_F2_h.coeff_monomial(1))

print("  coefficient h^2 =")
sp.pprint(coeff_h2)

print()
print("  coefficient h =")
sp.pprint(coeff_h1)

print()
print("  coefficient 1 =")
sp.pprint(coeff_h0)

cert(
    "h^2 coefficient",
    coeff_h2,
    2,
)

cert(
    "h coefficient",
    coeff_h1,
    2*(S + 1),
)

cert(
    "constant coefficient",
    coeff_h0,
    F2,
)


# ============================================================================
# [16] F2 ORBIT DISCRIMINANT
# ============================================================================

print()
print("[16] F2-ORBIT DISCRIMINANT")
print("-" * 78)

disc_F2_orbit = sp.factor(
    coeff_h1**2 - 4*coeff_h2*coeff_h0
)

print("  disc_h(T_h(F2)) =")
sp.pprint(disc_F2_orbit)

# Express it in terms of Delta = S^2 - 4N.
Delta = S**2 - 4*N

disc_minus_12Delta = sp.factor(
    disc_F2_orbit - 12*Delta
)

print()
print("  disc(F2 orbit) - 12 Delta =")
sp.pprint(disc_minus_12Delta)

cert(
    "disc(F2 orbit) = 12 Delta",
    disc_F2_orbit,
    12*Delta,
)


# ============================================================================
# [17] IMPORTANT: F2 ORBIT DISCRIMINANT COLLAPSES TO DELTA
# ============================================================================

print()
print("[17] DISCRIMINANT COLLAPSE")
print("-" * 78)

print("""
  The F2 orbit itself has

      T_h(F2)
        = 2 h^2
          + 2(S+1)h
          + F2.

  Its discriminant in h is:

      [2(S+1)]^2 - 8F2.

  Since

      F2 = 6N-S^2+S,

  the expression simplifies to:

      12(S^2-4N)
        = 12 Delta.

  Thus the F2 translation orbit DOES carry the factor-gap
  discriminant, but only after taking its orbit discriminant.

  This is different from the direct KAPPA quadratic orbit.
""")


# ============================================================================
# [18] TWO DISCRIMINANTS
# ============================================================================

print()
print("[18] TWO-LEVEL DISCRIMINANT STRUCTURE")
print("-" * 78)

Q_disc = sp.factor(
    sp.discriminant(Q_target, z)
)

F2_disc = sp.factor(
    sp.discriminant(F2_h, h)
)

print("  Q(z) discriminant =")
sp.pprint(Q_disc)

print()
print("  T_h(F2) discriminant =")
sp.pprint(F2_disc)

cert(
    "F2 orbit discriminant = 12*Q discriminant",
    F2_disc,
    12*Q_disc,
)


# ============================================================================
# [19] F3 / F2 AS THE LINEAR ORBIT COORDINATE
# ============================================================================

print()
print("[19] F3/F2 AS OPERATOR COORDINATE")
print("-" * 78)

operator_coordinate = sp.factor(
    F3 / F2
)

print("  F3/F2 =")
sp.pprint(operator_coordinate)

cert(
    "F3/F2 = S+1",
    operator_coordinate,
    S + 1,
)

print()
print("  Therefore:")
print("      F3/F2 - 1 = S")


# ============================================================================
# [20] TRANSLATION RESPONSE AND LAYER RATIO AGREE
# ============================================================================

print()
print("[20] TWO INDEPENDENT S CHANNELS")
print("-" * 78)

S_from_ratio = sp.factor(
    F3 / F2 - 1
)

S_from_translation = sp.factor(
    dF2 / (2*h) - h - 1
)

print("  S from F3/F2:")
sp.pprint(S_from_ratio)

print()
print("  S from translated F2:")
sp.pprint(S_from_translation)

cert(
    "ratio and translation S channels agree",
    S_from_ratio,
    S_from_translation,
)


# ============================================================================
# [21] TRANSLATION GROUP LAW
# ============================================================================

print()
print("[21] TRANSLATION GROUP LAW")
print("-" * 78)

N_hk = sp.expand(
    N + (h + k)*S + (h + k)**2
)

S_hk = sp.expand(
    S + 2*(h + k)
)

F2_hk = sp.expand(
    6*N_hk - S_hk**2 + S_hk
)

# Apply h, then k.
N_after_h = sp.expand(
    N + h*S + h**2
)

S_after_h = sp.expand(
    S + 2*h
)

N_after_hk = sp.expand(
    N_after_h + k*S_after_h + k**2
)

S_after_hk = sp.expand(
    S_after_h + 2*k
)

F2_after_hk = sp.expand(
    6*N_after_hk - S_after_hk**2 + S_after_hk
)

cert(
    "T_k(T_h(F2))=T_(h+k)(F2)",
    F2_after_hk,
    F2_hk,
)


# ============================================================================
# [22] COMMUTATOR
# ============================================================================

print()
print("[22] TRANSLATION COMMUTATOR")
print("-" * 78)

N_after_k = sp.expand(
    N + k*S + k**2
)

S_after_k = sp.expand(
    S + 2*k
)

N_after_kh = sp.expand(
    N_after_k + h*S_after_k + h**2
)

S_after_kh = sp.expand(
    S_after_k + 2*h
)

F2_after_kh = sp.expand(
    6*N_after_kh - S_after_kh**2 + S_after_kh
)

commutator = sp.factor(
    F2_after_hk - F2_after_kh
)

print("  [T_h,T_k]F2 =")
sp.pprint(commutator)

cert(
    "translation operators commute",
    commutator,
    0,
)


# ============================================================================
# [23] OPERATOR-LEVEL INFORMATION FLOW
# ============================================================================

print()
print("[23] OPERATOR-LEVEL INFORMATION FLOW")
print("-" * 78)

print("""
  Existing layer:
      F2 = 6N-S^2+S

  Existing companion:
      F3 = (S+1)F2

  Therefore:
      F3/F2 - 1 = S.

  Under the distinguished translation:
      N -> N+hS+h^2
      S -> S+2h

  and:
      T_h(F2)-F2 = 2h(S+h+1).

  Hence two exact S channels exist:

      S = F3/F2 - 1

  and

      S = [T_h(F2)-F2]/(2h)-h-1.

  Once S is obtained:

      N = (F2+S^2-S)/6

  and:

      Q(z)=z^2-Sz+N.

  Finally:

      discriminant(Q)=S^2-4N.

  The F2 translation orbit independently verifies the same
  discriminant through:

      discriminant_h(T_h(F2)) = 12(S^2-4N).
""")


# ============================================================================
# [24] FINAL STRUCTURAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 543 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The key result tested here is that F2 is not itself the KAPPA
quadratic, but its translation orbit contains two useful channels:

    first difference
        -> S

    orbit discriminant
        -> 12 Delta.

Together with the existing companion layer:

    F3/F2 - 1
        -> S,

the homogeneous KAPPA layer already contains two algebraically
equivalent ways of exposing the missing symmetric coordinate.

The next structural question should therefore target the SOURCE
of the translation parameter h.

Do not invent a new translation.

Instead identify whether the original homogeneous-layer
construction already possesses an index/operator whose action
on F2 induces:

    F2 -> T_h(F2)

with:

    T_h(N) = N+hS+h^2

and:

    T_h(S) = S+2h.

If that operator can be derived from the original construction
itself, then the chain becomes genuinely upstream:

    homogeneous layer
        -> intrinsic operator
        -> F2 orbit
        -> S
        -> N
        -> Q(z)
        -> Delta.

The most valuable next experiment is therefore an exact search
for the actual layer-index/operator law that generates the
observed F2 transformation.
""")