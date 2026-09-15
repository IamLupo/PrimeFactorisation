#!/usr/bin/env python3

import sympy as sp


# ============================================================================
# EXPERIMENT 544
# ============================================================================
# CORRECTED F2 ORBIT DISCRIMINANT -> DELTA -> N CLOSURE
#
# No numerical data.
# No factor enumeration.
# No continued fractions.
# No p,q.
# ============================================================================


print("=" * 78)
print("EXPERIMENT 544 START")
print("=" * 78)
print("CORRECTED F2 ORBIT DISCRIMINANT -> DELTA -> N CLOSURE")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

N, S, h, k, z = sp.symbols("N S h k z")


failures = 0


def simp(expr):
    return sp.factor(
        sp.expand(
            sp.cancel(
                sp.simplify(expr)
            )
        )
    )


def cert(label, expr, expected=0):
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
# [1] BASIC KAPPA INVARIANTS
# ============================================================================

print()
print("[1] BASIC KAPPA INVARIANTS")
print("-" * 78)

Delta = sp.expand(S**2 - 4*N)

F2 = sp.expand(
    6*N - S**2 + S
)

F3 = sp.expand(
    (S + 1)*F2
)

print("  F2 =")
sp.pprint(F2)

print()
print("  F3 =")
sp.pprint(F3)

print()
print("  Delta =")
sp.pprint(Delta)

cert(
    "F3=(S+1)F2",
    F3,
    (S + 1)*F2,
)


# ============================================================================
# [2] FIRST IMPORTANT REWRITE OF F2
# ============================================================================

print()
print("[2] F2 IN TERMS OF N, S, DELTA")
print("-" * 78)

F2_rewrite = sp.expand(
    2*N + S - Delta
)

print("  proposed:")
sp.pprint(F2_rewrite)

cert(
    "F2 = 2N + S - Delta",
    F2,
    F2_rewrite,
)


# ============================================================================
# [3] TRANSLATION OPERATOR
# ============================================================================

print()
print("[3] DISTINGUISHED TRANSLATION")
print("-" * 78)

N_h = sp.expand(
    N + h*S + h**2
)

S_h = sp.expand(
    S + 2*h
)

Delta_h = sp.expand(
    S_h**2 - 4*N_h
)

print("  T_h(N) =")
sp.pprint(N_h)

print()
print("  T_h(S) =")
sp.pprint(S_h)

print()
print("  T_h(Delta) =")
sp.pprint(Delta_h)

cert(
    "N translation",
    N_h,
    N + h*S + h**2,
)

cert(
    "S translation",
    S_h,
    S + 2*h,
)

cert(
    "Delta is translation invariant",
    Delta_h,
    Delta,
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

dF2 = sp.expand(
    F2_h - F2
)

print("  T_h(F2) =")
sp.pprint(sp.factor(F2_h))

print()
print("  T_h(F2)-F2 =")
sp.pprint(sp.factor(dF2))

cert(
    "F2 translation law",
    dF2,
    2*h*(S + h + 1),
)


# ============================================================================
# [5] CORRECT FIRST-DIFFERENCE CHANNEL
# ============================================================================

print()
print("[5] CORRECT FIRST-DIFFERENCE CHANNEL")
print("-" * 78)

first_channel = sp.cancel(
    dF2 / (2*h)
)

print("  Delta(F2)/(2h) =")
sp.pprint(first_channel)

cert(
    "Delta(F2)/(2h) = S+h+1",
    first_channel,
    S + h + 1,
)

S_from_first_difference = sp.cancel(
    first_channel - h - 1
)

print()
print("  recovered S =")
sp.pprint(S_from_first_difference)

cert(
    "S recovery",
    S_from_first_difference,
    S,
)


# ============================================================================
# [6] COMPARE F2 DIFFERENCE WITH KAPPA ORBIT DIFFERENCE
# ============================================================================

print()
print("[6] F2 DIFFERENCE VS KAPPA ORBIT")
print("-" * 78)

L_h = sp.expand(
    N + h*S + h**2
)

Delta_L_h = sp.expand(
    L_h.subs(h, h + 1) - L_h
)

print("  L(h) =")
sp.pprint(L_h)

print()
print("  one-step Delta L(h) =")
sp.pprint(Delta_L_h)

print()
print("  Delta(F2)/(2h) =")
sp.pprint(first_channel)

cert(
    "Delta(F2)/(2h) - (Delta L(h)-h)",
    first_channel,
    Delta_L_h - h,
)


# ============================================================================
# [7] F2 ORBIT POLYNOMIAL
# ============================================================================

print()
print("[7] F2 ORBIT POLYNOMIAL")
print("-" * 78)

poly_F2 = sp.Poly(
    F2_h,
    h
)

a = sp.expand(
    poly_F2.coeff_monomial(h**2)
)

b = sp.expand(
    poly_F2.coeff_monomial(h)
)

c = sp.expand(
    poly_F2.coeff_monomial(1)
)

print("  a = coefficient(h^2) =")
sp.pprint(a)

print()
print("  b = coefficient(h) =")
sp.pprint(b)

print()
print("  c = coefficient(1) =")
sp.pprint(c)

cert(
    "quadratic coefficient",
    a,
    2,
)

cert(
    "linear coefficient",
    b,
    2*(S + 1),
)

cert(
    "constant coefficient",
    c,
    F2,
)


# ============================================================================
# [8] CORRECTED ORBIT DISCRIMINANT
# ============================================================================

print()
print("[8] CORRECTED F2-ORBIT DISCRIMINANT")
print("-" * 78)

D_F2 = sp.factor(
    b**2 - 4*a*c
)

print("  discriminant_h(T_h(F2)) =")
sp.pprint(D_F2)

cert(
    "D_F2 = 12*Delta + 4",
    D_F2,
    12*Delta + 4,
)


# ============================================================================
# [9] DIRECT DELTA RECOVERY FROM F2 ORBIT
# ============================================================================

print()
print("[9] DELTA RECOVERY FROM F2 ORBIT")
print("-" * 78)

Delta_from_orbit = sp.factor(
    (D_F2 - 4) / 12
)

print("  Delta_candidate =")
sp.pprint(Delta_from_orbit)

cert(
    "Delta=(D_F2-4)/12",
    Delta_from_orbit,
    Delta,
)


# ============================================================================
# [10] N RECOVERY USING F2, S, DELTA
# ============================================================================

print()
print("[10] N RECOVERY FROM F2 + S + DELTA")
print("-" * 78)

N_from_F2_S_Delta = sp.factor(
    (F2 + S + Delta*0 + Delta - S) / 2
)

print("  N_candidate =")
sp.pprint(N_from_F2_S_Delta)

cert(
    "N=(F2+Delta-S)/2",
    N_from_F2_S_Delta,
    N,
)


# ============================================================================
# [11] ALTERNATIVE N RECOVERY
# ============================================================================

print()
print("[11] ALTERNATIVE N RECOVERY")
print("-" * 78)

N_from_F2_S = sp.factor(
    (F2 + S**2 - S) / 6
)

N_from_F2_Delta_S = sp.factor(
    (F2 + Delta - S) / 2
)

print("  Route A:")
print("    N = (F2+S^2-S)/6")
sp.pprint(N_from_F2_S)

print()
print("  Route B:")
print("    N = (F2+Delta-S)/2")
sp.pprint(N_from_F2_Delta_S)

cert(
    "Route A = N",
    N_from_F2_S,
    N,
)

cert(
    "Route B = N",
    N_from_F2_Delta_S,
    N,
)

cert(
    "Route A = Route B",
    N_from_F2_S,
    N_from_F2_Delta_S,
)


# ============================================================================
# [12] F3 COMPANION CHANNEL
# ============================================================================

print()
print("[12] F3/F2 S CHANNEL")
print("-" * 78)

ratio = sp.cancel(
    F3 / F2
)

print("  F3/F2 =")
sp.pprint(ratio)

cert(
    "F3/F2 = S+1",
    ratio,
    S + 1,
)

S_ratio = sp.factor(
    ratio - 1
)

cert(
    "F3/F2-1 = S",
    S_ratio,
    S,
)


# ============================================================================
# [13] TWO COMPLETE RECOVERY ROUTES
# ============================================================================

print()
print("[13] TWO COMPLETE SYMBOLIC RECOVERY ROUTES")
print("-" * 78)

print("""
  ROUTE A:

      F2, F3
        |
        v
      S = F3/F2 - 1
        |
        v
      Delta = S^2 - 4N
        |
        v
      Q(z)

  ROUTE B:

      F2 + translation operator
        |
        v
      orbit discriminant D_F2
        |
        v
      Delta = (D_F2-4)/12
        |
        +---- first difference -> S
        |
        v
      N = (F2+Delta-S)/2
        |
        v
      Q(z)
""")


# ============================================================================
# [14] OPERATOR-INVARIANT DISCRIMINANT CHANNEL
# ============================================================================

print()
print("[14] DISCRIMINANT AS TRANSLATION INVARIANT")
print("-" * 78)

D_F2_h_symbolic = sp.factor(
    sp.discriminant(F2_h, h)
)

D_F2_under_shift = sp.factor(
    D_F2_h_symbolic.subs(
        {
            N: N_h,
            S: S_h,
        }
    )
)

cert(
    "orbit discriminant invariant",
    D_F2_under_shift,
    D_F2_h_symbolic,
)


# ============================================================================
# [15] DELTA INVARIANCE
# ============================================================================

print()
print("[15] DELTA INVARIANCE THROUGH F2 ORBIT")
print("-" * 78)

Delta_recovered_shift = sp.factor(
    (D_F2_h_symbolic - 4) / 12
)

cert(
    "Delta from orbit discriminant",
    Delta_recovered_shift,
    Delta,
)

cert(
    "Delta remains invariant",
    Delta_h,
    Delta,
)


# ============================================================================
# [16] F2 AS AN AFFINE FUNCTION OF DELTA
# ============================================================================

print()
print("[16] F2 / DELTA LINEAR DECOMPOSITION")
print("-" * 78)

F2_delta_form = sp.expand(
    2*N + S - Delta
)

print("  F2 = 2N + S - Delta")

cert(
    "F2 delta decomposition",
    F2,
    F2_delta_form,
)


# ============================================================================
# [17] ORBIT DISCRIMINANT AS A PURE INVARIANT
# ============================================================================

print()
print("[17] PURE INVARIANT FORM")
print("-" * 78)

pure_invariant = sp.expand(
    12*(S**2 - 4*N) + 4
)

print("  D_F2 =")
sp.pprint(pure_invariant)

cert(
    "D_F2 pure invariant",
    D_F2,
    pure_invariant,
)


# ============================================================================
# [18] S-DISCRIMINANT COMBINATION
# ============================================================================

print()
print("[18] S CHANNEL + ORBIT DISCRIMINANT")
print("-" * 78)

S_from_orbit_difference = sp.factor(
    b/2 - h*0 - 1
)

print("  S from orbit coefficient b:")
sp.pprint(S_from_orbit_difference)

cert(
    "b/2-1 = S",
    S_from_orbit_difference,
    S,
)


# ============================================================================
# [19] QUADRATIC RECOVERY WITHOUT EXPLICIT FACTORS
# ============================================================================

print()
print("[19] SPECTRAL QUADRATIC RECOVERY")
print("-" * 78)

S_eff = sp.factor(
    b/2 - 1
)

Delta_eff = sp.factor(
    (D_F2 - 4) / 12
)

N_eff = sp.factor(
    (F2 + Delta_eff - S_eff) / 2
)

Q_eff = sp.expand(
    z**2 - S_eff*z + N_eff
)

Q_target = sp.expand(
    z**2 - S*z + N
)

print("  S_eff =")
sp.pprint(S_eff)

print()
print("  Delta_eff =")
sp.pprint(Delta_eff)

print()
print("  N_eff =")
sp.pprint(N_eff)

print()
print("  Q_eff =")
sp.pprint(Q_eff)

cert(
    "S_eff=S",
    S_eff,
    S,
)

cert(
    "Delta_eff=Delta",
    Delta_eff,
    Delta,
)

cert(
    "N_eff=N",
    N_eff,
    N,
)

cert(
    "Q_eff=Q",
    Q_eff,
    Q_target,
)


# ============================================================================
# [20] CHARACTERISTIC IDENTITY
# ============================================================================

print()
print("[20] SPECTRAL DISCRIMINANT")
print("-" * 78)

Q_disc = sp.factor(
    sp.discriminant(Q_eff, z)
)

print("  discriminant(Q_eff) =")
sp.pprint(Q_disc)

cert(
    "discriminant(Q_eff)=Delta",
    Q_disc,
    Delta,
)


# ============================================================================
# [21] TRANSLATED KAPPA FACTOR
# ============================================================================

print()
print("[21] +1 KAPPA TRANSLATION")
print("-" * 78)

Q_shift = sp.expand(
    Q_eff.subs(z, z - 1)
)

Q_shift_target = sp.expand(
    (z - (S + 1)) * (z - (0 + 1))
    - (z - (S + 1))
)

#
# Better exact target:
#
# Q(z-1) = (z-1)^2 - S(z-1) + N
#

Q_shift_target = sp.expand(
    (z - 1)**2 - S*(z - 1) + N
)

cert(
    "Q_eff(z-1)=Q(z-1)",
    Q_shift,
    Q_shift_target,
)


# ============================================================================
# [22] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ============================================================================

print()
print("[22] FOUR-BASE COMPOSITION")
print("-" * 78)

chi_eff = sp.expand(
    Q_eff * Q_shift
)

chi_target = sp.expand(
    Q_target * Q_shift_target
)

cert(
    "Q_eff*Q_eff(z-1)",
    chi_eff,
    chi_target,
)


# ============================================================================
# [23] KEY CORRECTION CERTIFICATE
# ============================================================================

print()
print("[23] CORRECTED IDENTITIES")
print("-" * 78)

print("""
  Correct:

      Delta(F2)/(2h) = S+h+1

  not:

      Delta(F2)/(2h) = Delta L(h).

  Also:

      disc_h(T_h(F2))
        = 12 Delta + 4

  not:

      12 Delta.

  Therefore:

      Delta
        = [disc_h(T_h(F2))-4]/12.

  And since:

      F2 = 2N + S - Delta,

  we obtain:

      N = [F2 + Delta - S]/2.

  This gives an exact closure from the translated F2 orbit.
""")


# ============================================================================
# [24] OPERATOR-LEVEL INTERPRETATION
# ============================================================================

print()
print("[24] OPERATOR-LEVEL INTERPRETATION")
print("-" * 78)

print("""
  The corrected structure is:

      F2
       |
       +---- coefficient of h
       |        |
       |        v
       |      S+1
       |
       +---- discriminant in h
                |
                v
            12 Delta + 4

  Therefore the orbit carries two separate channels:

      linear coefficient
          -> S

      orbit discriminant
          -> Delta.

  The constant layer value then supplies N:

      N = (F2 + Delta - S)/2.

  Consequently the complete quadratic

      Q(z)=z^2-Sz+N

  can be reconstructed from the operator orbit coefficients.

  The unresolved upstream problem is now very narrow:

      where does the actual translation parameter/operator
      acting on F2 originate inside the original homogeneous
      construction?

  This experiment does NOT assume such an operator exists
  upstream. It only identifies exactly what its observable
  output would have to be.
""")


# ============================================================================
# FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 544 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The corrected downstream invariant package is now:

    F2 = 6N-S^2+S

    coefficient_h(T_h(F2)) = 2(S+1)

    discriminant_h(T_h(F2)) = 12 Delta + 4

    Delta = (D_F2-4)/12

    N = (F2+Delta-S)/2

    Q(z)=z^2-Sz+N.

The next experiment should therefore attack only the missing
SOURCE of the translation parameter.

Do not invent another quadratic.

Search the existing homogeneous-layer index itself for an
operator T_m such that its action on F2 has coefficients:

    h^2 coefficient = 2

    h coefficient = 2(S+1)

    constant = F2.

Equivalently, look for an intrinsic layer transformation for
which the discriminant of the resulting quadratic orbit is
translation-invariant.

The target is now the operator origin, not another downstream
identity.
""")
