#!/usr/bin/env python3

import sympy as sp


# ============================================================================
# EXPERIMENT 545
# ============================================================================
# F2/F3-ONLY TRANSLATION OPERATOR
#
# Goal:
#   Remove N and S completely from the observable interface.
#
# Starting observables:
#   F2 = 6N - S^2 + S
#   F3 = (S+1)F2
#
# Since:
#   S = F3/F2 - 1
#
# derive the distinguished translation operator directly on (F2,F3).
#
# No p,q.
# No numerical data.
# No factor enumeration.
# No continued fractions.
# ============================================================================


print("=" * 78)
print("EXPERIMENT 545 START")
print("=" * 78)
print("F2/F3-ONLY TRANSLATION OPERATOR")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

N, S = sp.symbols("N S")
h, k = sp.symbols("h k")
z = sp.symbols("z")

A, B = sp.symbols("A B")


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
# [1] ORIGINAL HOMOGENEOUS LAYERS
# ============================================================================

print()
print("[1] ORIGINAL HOMOGENEOUS LAYERS")
print("-" * 78)

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

cert(
    "F3=(S+1)F2",
    F3,
    (S + 1)*F2,
)


# ============================================================================
# [2] OBSERVABLE S CHANNEL
# ============================================================================

print()
print("[2] OBSERVABLE S CHANNEL")
print("-" * 78)

S_obs = sp.cancel(
    F3/F2 - 1
)

print("  S_obs = F3/F2 - 1")
sp.pprint(S_obs)

cert(
    "S_obs=S",
    S_obs,
    S,
)


# ============================================================================
# [3] OBSERVABLE N CHANNEL
# ============================================================================

print()
print("[3] OBSERVABLE N CHANNEL")
print("-" * 78)

N_obs = sp.cancel(
    (F2 + S_obs**2 - S_obs)/6
)

print("  N_obs =")
sp.pprint(N_obs)

cert(
    "N_obs=N",
    N_obs,
    N,
)


# ============================================================================
# [4] DISTINGUISHED TRANSLATION IN N,S
# ============================================================================

print()
print("[4] REFERENCE TRANSLATION")
print("-" * 78)

N_h = sp.expand(
    N + h*S + h**2
)

S_h = sp.expand(
    S + 2*h
)

F2_h_reference = sp.expand(
    6*N_h - S_h**2 + S_h
)

F3_h_reference = sp.expand(
    (S_h + 1)*F2_h_reference
)

print("  T_h(N) =")
sp.pprint(N_h)

print()
print("  T_h(S) =")
sp.pprint(S_h)

print()
print("  reference T_h(F2) =")
sp.pprint(sp.factor(F2_h_reference))

print()
print("  reference T_h(F3) =")
sp.pprint(sp.factor(F3_h_reference))


# ============================================================================
# [5] ELIMINATE N,S ENTIRELY
# ============================================================================

print()
print("[5] ELIMINATE N AND S")
print("-" * 78)

# A represents F2
# B represents F3
#
# S = B/A - 1
#
# Therefore:
#
# T_h(F2)
#   = A + 2h(S+h+1)
#   = A + 2h(B/A+h).

T_A = sp.factor(
    A + 2*h*(B/A + h)
)

print("  A = F2")
print("  B = F3")

print()
print("  T_h(A) =")
sp.pprint(T_A)

cert(
    "observable T_h(F2) formula",
    T_A.subs({A: F2, B: F3}) - F2_h_reference,
    0,
)


# ============================================================================
# [6] OBSERVABLE TRANSLATION FORM
# ============================================================================

print()
print("[6] F2-ONLY/F3-ONLY TRANSLATION RESPONSE")
print("-" * 78)

Delta_A = sp.factor(
    T_A - A
)

print("  T_h(A)-A =")
sp.pprint(Delta_A)

cert(
    "T_h(A)-A = 2h(B/A+h)",
    Delta_A,
    2*h*(B/A + h),
)


# ============================================================================
# [7] RECONSTRUCT F3 TRANSFORMATION
# ============================================================================

print()
print("[7] F3 TRANSFORMATION")
print("-" * 78)

#
# Since:
#
#   F3 = (S+1) F2
#
# and:
#
#   S+1 = B/A,
#
# after translation:
#
#   S'+1 = B/A + 2h.
#
# Hence:
#
#   B' = (B/A + 2h) A'.
#

T_B = sp.factor(
    (B/A + 2*h)*T_A
)

print("  T_h(B) =")
sp.pprint(T_B)

cert(
    "observable T_h(F3) formula",
    T_B.subs({A: F2, B: F3}) - F3_h_reference,
    0,
)


# ============================================================================
# [8] FIRST DIFFERENCE OF F3
# ============================================================================

print()
print("[8] F3 TRANSLATION RESPONSE")
print("-" * 78)

Delta_B = sp.factor(
    T_B - B
)

print("  T_h(B)-B =")
sp.pprint(Delta_B)

reference_Delta_B = sp.factor(
    F3_h_reference - F3
)

cert(
    "F3 response matches reference",
    Delta_B.subs({A: F2, B: F3}) - reference_Delta_B,
    0,
)


# ============================================================================
# [9] QUOTIENT EVOLUTION
# ============================================================================

print()
print("[9] QUOTIENT EVOLUTION")
print("-" * 78)

ratio = sp.cancel(B/A)

ratio_translated = sp.cancel(
    T_B/T_A
)

print("  B/A =")
sp.pprint(ratio)

print()
print("  T_h(B)/T_h(A) =")
sp.pprint(ratio_translated)

cert(
    "ratio transforms as S+1 -> S+1+2h",
    ratio_translated,
    B/A + 2*h,
)


# ============================================================================
# [10] DIRECT S RECOVERY AFTER TRANSLATION
# ============================================================================

print()
print("[10] DIRECT S RECOVERY FROM TRANSLATED LAYERS")
print("-" * 78)

S_after_observables = sp.cancel(
    T_B/T_A - 1
)

print("  S' = T_h(F3)/T_h(F2)-1")
sp.pprint(S_after_observables)

cert(
    "translated S channel",
    S_after_observables,
    S_obs + 2*h,
)


# ============================================================================
# [11] DIRECT N RECOVERY AFTER TRANSLATION
# ============================================================================

print()
print("[11] DIRECT N RECOVERY AFTER TRANSLATION")
print("-" * 78)

N_after_observables = sp.cancel(
    (
        T_A
        + S_after_observables**2
        - S_after_observables
    ) / 6
)

print("  N' recovered from translated F2/F3 =")
sp.pprint(sp.factor(N_after_observables))

cert(
    "translated N channel",
    N_after_observables,
    N_obs + h*S_obs + h**2,
)


# ============================================================================
# [12] PURE OBSERVABLE TRANSFORMATION LAW
# ============================================================================

print()
print("[12] CLOSED TWO-OBSERVABLE OPERATOR")
print("-" * 78)

print("""
  With:

      A = F2
      B = F3

  the translation closes entirely on (A,B):

      A' = A + 2h(B/A + h)

      B' = (B/A + 2h) A'.

  Equivalently:

      A' = A + 2h B/A + 2h^2

      B' = (B+2hA)A'/A.
""")

# Verify expanded closed forms.

T_A_closed = sp.factor(
    A + 2*h*B/A + 2*h**2
)

T_B_closed = sp.factor(
    (B + 2*h*A)*T_A_closed/A
)

cert(
    "closed A transformation",
    T_A_closed,
    T_A,
)

cert(
    "closed B transformation",
    T_B_closed,
    T_B,
)


# ============================================================================
# [13] GROUP COMPOSITION
# ============================================================================

print()
print("[13] GROUP COMPOSITION IN F2/F3 SPACE")
print("-" * 78)

def T_A_of(A0, B0, hh):
    return sp.factor(
        A0 + 2*hh*B0/A0 + 2*hh**2
    )


def T_B_of(A0, B0, hh):
    A1 = T_A_of(A0, B0, hh)
    return sp.factor(
        (B0 + 2*hh*A0)*A1/A0
    )


A_h = T_A_of(A, B, h)
B_h = T_B_of(A, B, h)

A_hk = T_A_of(A_h, B_h, k)
B_hk = T_B_of(A_h, B_h, k)

A_sum = T_A_of(A, B, h + k)
B_sum = T_B_of(A, B, h + k)

cert(
    "T_k(T_h(A))=T_(h+k)(A)",
    A_hk,
    A_sum,
)

cert(
    "T_k(T_h(B))=T_(h+k)(B)",
    B_hk,
    B_sum,
)


# ============================================================================
# [14] COMMUTATOR
# ============================================================================

print()
print("[14] OBSERVABLE OPERATOR COMMUTATOR")
print("-" * 78)

A_k = T_A_of(A, B, k)
B_k = T_B_of(A, B, k)

A_kh = T_A_of(A_k, B_k, h)
B_kh = T_B_of(A_k, B_k, h)

cert(
    "[T_h,T_k]A=0",
    A_hk,
    A_kh,
)

cert(
    "[T_h,T_k]B=0",
    B_hk,
    B_kh,
)


# ============================================================================
# [15] INVARIANT DISCRIMINANT CHANNEL
# ============================================================================

print()
print("[15] INVARIANT DISCRIMINANT CHANNEL")
print("-" * 78)

#
# The observable expression for Delta is:
#
#   S = B/A - 1
#   N = (A+S^2-S)/6
#
#   Delta = S^2 - 4N
#

Delta_obs = sp.factor(
    S_obs**2 - 4*N_obs
)

Delta_translated = sp.factor(
    S_after_observables**2 - 4*N_after_observables
)

print("  Delta(F2,F3) =")
sp.pprint(Delta_obs)

print()
print("  Delta(T_h(F2),T_h(F3)) =")
sp.pprint(Delta_translated)

cert(
    "observable Delta invariant",
    Delta_translated,
    Delta_obs,
)


# ============================================================================
# [16] OBSERVABLE DISCRIMINANT WITHOUT N OR S
# ============================================================================

print()
print("[16] DIRECT DELTA FORMULA IN F2,F3")
print("-" * 78)

Delta_AB = sp.factor(
    Delta_obs.subs(
        {
            F2: A,
            F3: B,
        }
    )
)

print("  Delta(A,B) =")
sp.pprint(Delta_AB)

# Independently derive it without symbolic substitution from N,S.

S_AB = sp.cancel(
    B/A - 1
)

N_AB = sp.cancel(
    (A + S_AB**2 - S_AB)/6
)

Delta_AB_direct = sp.factor(
    S_AB**2 - 4*N_AB
)

cert(
    "Delta_AB direct construction",
    Delta_AB_direct,
    Delta_AB,
)


# ============================================================================
# [17] GENERATING QUADRATIC DIRECTLY FROM F2/F3
# ============================================================================

print()
print("[17] KAPPA QUADRATIC FROM F2/F3 ONLY")
print("-" * 78)

Q_AB = sp.factor(
    z**2
    - S_AB*z
    + N_AB
)

print("  Q_AB(z) =")
sp.pprint(Q_AB)

# Substitute back the actual F2/F3 definitions.

Q_NS = sp.expand(
    z**2 - S*z + N
)

cert(
    "Q(F2,F3)=z^2-Sz+N",
    Q_AB.subs({A: F2, B: F3}),
    Q_NS,
)


# ============================================================================
# [18] F2 ORBIT DISCRIMINANT CORRECTION
# ============================================================================

print()
print("[18] F2 ORBIT DISCRIMINANT")
print("-" * 78)

F2_h_poly = sp.Poly(
    F2_h_reference,
    h
)

f2_a = sp.expand(
    F2_h_poly.coeff_monomial(h**2)
)

f2_b = sp.expand(
    F2_h_poly.coeff_monomial(h)
)

f2_c = sp.expand(
    F2_h_poly.coeff_monomial(1)
)

f2_disc = sp.factor(
    f2_b**2 - 4*f2_a*f2_c
)

print("  disc_h(T_h(F2)) =")
sp.pprint(f2_disc)

cert(
    "F2 orbit discriminant",
    f2_disc,
    12*(S**2 - 4*N) + 4,
)

Delta_from_F2_orbit = sp.factor(
    (f2_disc - 4)/12
)

cert(
    "Delta=(disc-4)/12",
    Delta_from_F2_orbit,
    S**2 - 4*N,
)


# ============================================================================
# [19] F2/F3 OBSERVABLE RECOVERY OF ORBIT DISCRIMINANT
# ============================================================================

print()
print("[19] ORBIT DISCRIMINANT AS OBSERVABLE")
print("-" * 78)

S_from_AB = S_AB
N_from_AB = N_AB

orbit_disc_AB = sp.factor(
    12*(S_from_AB**2 - 4*N_from_AB) + 4
)

print("  D_orbit(F2,F3) =")
sp.pprint(orbit_disc_AB)

cert(
    "observable orbit discriminant",
    orbit_disc_AB,
    f2_disc.subs(
        {
            N: N_from_AB,
            S: S_from_AB
        }
    ),
)


# ============================================================================
# [20] MINIMAL OBSERVABLE CLOSURE
# ============================================================================

print()
print("[20] MINIMAL OBSERVABLE CLOSURE")
print("-" * 78)

print("""
  The pair:

      (F2,F3)

  is closed under the distinguished translation.

  Define:

      A = F2
      B = F3.

  Then:

      S = B/A - 1

      N = [A + S^2 - S]/6

      Delta = S^2 - 4N

  and the translation becomes:

      A' = A + 2h(B/A + h)

      B' = (B/A + 2h) A'.

  Therefore the entire (N,S,Delta) state is encoded by
  the two existing homogeneous observables F2 and F3.

  No explicit p or q is required at this stage.
""")


# ============================================================================
# [21] OPERATOR INVARIANTS / MOVING COORDINATE
# ============================================================================

print()
print("[21] OPERATOR COORDINATES IN OBSERVABLE SPACE")
print("-" * 78)

print("""
  Observable quotient coordinate:

      R = F3/F2 - 1

  Under translation:

      R -> R + 2h.

  Therefore:

      R/2

  is the moving coordinate.

  The invariant is:

      Delta(F2,F3).

  Thus the observable space already has the same
  moving/invariant decomposition found in the historical
  conic system:

      moving:
          S/2

      invariant:
          Delta.
""")


# ============================================================================
# [22] FINAL STRUCTURAL CERTIFICATE
# ============================================================================

print()
print("[22] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print("""
  EXISTING LAYERS:

      F2 = 6N-S^2+S

      F3 = (S+1)F2

  OBSERVABLE COORDINATES:

      R = F3/F2 - 1 = S

      N = (F2+R^2-R)/6

  KAPPA QUADRATIC:

      Q(z)=z^2-Rz+N

  DISCRIMINANT:

      Delta = R^2-4N.

  TRANSLATION:

      R -> R+2h

      N -> N+hR+h^2

      Delta -> Delta.

  F2 response:

      F2' - F2 = 2h(R+h+1).

  F2 orbit discriminant:

      disc_h(F2') = 12Delta+4.

  Therefore the observable pair (F2,F3) carries:

      moving symmetric coordinate R,

      invariant discriminant Delta,

      product coordinate N,

      generating polynomial Q.

  The remaining question is now extremely specific:

      Is there an actual INDEX SHIFT in the original
      homogeneous construction whose action on its existing
      layer pair is exactly the rational map above?

  This is stronger and more concrete than searching
  arbitrary polynomial identities.
""")


# ============================================================================
# FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 545 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The next experiment should no longer manufacture a new
quadratic or conic.

The exact rational transformation to search for in the
ORIGINAL homogeneous-layer index is:

    A = F2
    B = F3

    A' = A + 2h(B/A + h)

    B' = (B/A + 2h) A'.

Equivalently, the ratio obeys the exceptionally simple law:

    B'/A' = B/A + 2h.

This gives:

    S' = S + 2h.

The invariant channel is:

    Delta(F2,F3) = invariant.

The decisive upstream question is therefore:

    Does the genuine homogeneous-layer recurrence/index
    already implement a shift for which

        F3/F2 -> F3/F2 + 2h ?

If yes, then the KAPPA symmetric coordinate is an intrinsic
coordinate of the layer-index action:

    layer index
        -> F3/F2
        -> S
        -> N
        -> Q(z).

The strongest next test is to derive the actual index recurrence
for F2,F3 and compare it directly with this rational operator.
""")
