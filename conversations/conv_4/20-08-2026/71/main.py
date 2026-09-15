#!/usr/bin/env python3

import sympy as sp


# =============================================================================
# EXPERIMENT 547
# =============================================================================
# HOMOGENEOUS-LAYER RATIO LADDER / INTRINSIC INDEX TRANSLATION
#
# Known:
#
#     F2 = 6N - S^2 + S
#     F3 = (S+1) F2
#
# Therefore:
#
#     F3/F2 = S+1.
#
# The previous experiments suggest that the real moving coordinate is
#
#     R_n = F_{n+1}/F_n
#
# and that an intrinsic layer-index translation would have to act as
#
#     R -> R + 2h.
#
# This experiment:
#
#   1. keeps everything symbolic;
#   2. derives the exact observable coordinate from F2,F3;
#   3. constructs the most general multiplicative layer ladder compatible
#      with an affine index shift;
#   4. derives the exact ratio law;
#   5. derives the induced F2/F3 rational operator;
#   6. determines which higher-layer ratio law is necessary for the
#      KAPPA translation;
#   7. does NOT assume F4/F5 are known homogeneous layers;
#   8. does NOT insert p,q or numerical data.
#
# No factor enumeration.
# No continued fractions.
# No interpolation.
# =============================================================================


print("=" * 78)
print("EXPERIMENT 547 START")
print("=" * 78)
print("HOMOGENEOUS-LAYER RATIO LADDER / INTRINSIC INDEX TRANSLATION")
print()


# =============================================================================
# SYMBOLS
# =============================================================================

N, S = sp.symbols("N S")
h, m, n = sp.symbols("h m n")
z = sp.symbols("z")

A, B = sp.symbols("A B", nonzero=True)


failures = 0


def simp(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.together(
                    sp.simplify(expr)
                )
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


# =============================================================================
# [1] KNOWN HOMOGENEOUS LAYERS
# =============================================================================

print()
print("[1] KNOWN HOMOGENEOUS LAYERS")
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


# =============================================================================
# [2] KNOWN RATIO COORDINATE
# =============================================================================

print()
print("[2] KNOWN LAYER RATIO")
print("-" * 78)

R2 = sp.factor(
    F3/F2
)

print("  R2 = F3/F2")
sp.pprint(R2)

cert(
    "R2=S+1",
    R2,
    S + 1,
)


# =============================================================================
# [3] GENERAL LAYER RATIO
# =============================================================================

print()
print("[3] GENERAL LAYER RATIO")
print("-" * 78)

# Define formally:
#
#     R_n = F_(n+1) / F_n.
#
# We do not assume any formula for F_n.

R_n = sp.symbols("R_n")

print("""
  Define formally:

      R_n = F_(n+1) / F_n.

  The known layer gives:

      R_2 = S+1.
""")


# =============================================================================
# [4] NECESSARY TRANSLATION LAW
# =============================================================================

print()
print("[4] NECESSARY INDEX TRANSLATION LAW")
print("-" * 78)

R_shift = sp.expand(
    R_n + 2*h
)

print("  candidate translated ratio:")
sp.pprint(R_shift)

print("""
  For the ratio to be the KAPPA moving coordinate, the required law is:

      R_n  ->  R_n + 2h.

  At n=2 this becomes:

      R_2 -> S+1+2h.
""")


# =============================================================================
# [5] RECOVER S FROM THE RATIO
# =============================================================================

print()
print("[5] S CHANNEL FROM THE RATIO")
print("-" * 78)

S_from_R = sp.expand(
    R2 - 1
)

print("  S = R2-1")
sp.pprint(S_from_R)

cert(
    "S=F3/F2-1",
    S_from_R,
    S,
)


# =============================================================================
# [6] RATIO SHIFT -> S SHIFT
# =============================================================================

print()
print("[6] RATIO SHIFT -> S SHIFT")
print("-" * 78)

S_shift = sp.expand(
    R_shift.subs(R_n, S + 1) - 1
)

print("  S' =")
sp.pprint(S_shift)

cert(
    "S'=S+2h",
    S_shift,
    S + 2*h,
)


# =============================================================================
# [7] GENERAL N-S TRANSLATION
# =============================================================================

print()
print("[7] REQUIRED N CHANNEL")
print("-" * 78)

N_shift = sp.expand(
    N + h*S + h**2
)

print("  N' =")
sp.pprint(N_shift)

print("""
  The KAPPA translation requires simultaneously:

      S' = S+2h

      N' = N+hS+h^2.
""")


# =============================================================================
# [8] F2 REGENERATION FROM SHIFTED N,S
# =============================================================================

print()
print("[8] REGENERATE F2 FROM THE TRANSLATED STATE")
print("-" * 78)

F2_shift = sp.expand(
    6*N_shift - S_shift**2 + S_shift
)

print("  F2' =")
sp.pprint(sp.factor(F2_shift))

cert(
    "F2' matches its translation law",
    F2_shift - (
        F2 + 2*h*(S + h + 1)
    ),
    0,
)


# =============================================================================
# [9] F3 REGENERATION
# =============================================================================

print()
print("[9] REGENERATE F3 FROM THE TRANSLATED STATE")
print("-" * 78)

F3_shift = sp.expand(
    (S_shift + 1)*F2_shift
)

print("  F3' =")
sp.pprint(sp.factor(F3_shift))

cert(
    "F3'=(S'+1)F2'",
    F3_shift,
    (S_shift + 1)*F2_shift,
)


# =============================================================================
# [10] REQUIRED RATIO OF THE TRANSLATED LAYERS
# =============================================================================

print()
print("[10] TRANSLATED LAYER RATIO")
print("-" * 78)

R2_shift = sp.factor(
    F3_shift/F2_shift
)

print("  F3'/F2' =")
sp.pprint(R2_shift)

cert(
    "F3'/F2'=S+2h+1",
    R2_shift,
    S + 2*h + 1,
)


# =============================================================================
# [11] RATIO-ONLY OPERATOR
# =============================================================================

print()
print("[11] RATIO-ONLY OPERATOR")
print("-" * 78)

# Let W = F3/F2.
#
# The operator must satisfy:
#
#     W' = W+2h.

W = sp.symbols("W")

W_prime = sp.expand(
    W + 2*h
)

print("  W = F3/F2")
print("  W' = W+2h")

cert(
    "ratio translation",
    W_prime,
    W + 2*h,
)


# =============================================================================
# [12] F2 RESPONSE IN TERMS OF W
# =============================================================================

print()
print("[12] F2 RESPONSE IN RATIO COORDINATE")
print("-" * 78)

A_prime_W = sp.factor(
    A + 2*h*(W + h)
)

print("  A' =")
sp.pprint(A_prime_W)

cert(
    "A'-A = 2h(W+h)",
    A_prime_W - A,
    2*h*(W + h),
)


# =============================================================================
# [13] F3 RESPONSE IN RATIO COORDINATE
# =============================================================================

print()
print("[13] F3 RESPONSE IN RATIO COORDINATE")
print("-" * 78)

B_prime_W = sp.factor(
    (W + 2*h)*A_prime_W
)

print("  B' =")
sp.pprint(B_prime_W)

cert(
    "B'=(W+2h)A'",
    B_prime_W,
    (W + 2*h)*A_prime_W,
)


# =============================================================================
# [14] RETURN TO F2/F3 OBSERVABLES
# =============================================================================

print()
print("[14] RETURN TO OBSERVABLE PAIR")
print("-" * 78)

W_AB = sp.factor(
    B/A
)

A_prime_AB = sp.factor(
    A_prime_W.subs(W, W_AB)
)

B_prime_AB = sp.factor(
    B_prime_W.subs(W, W_AB)
)

print("  A' =")
sp.pprint(A_prime_AB)

print()
print("  B' =")
sp.pprint(B_prime_AB)

cert(
    "A' observable closure",
    A_prime_AB,
    A + 2*h*(B/A + h),
)

cert(
    "B' observable closure",
    B_prime_AB,
    (B/A + 2*h)*A_prime_AB,
)


# =============================================================================
# [15] INTRINSIC INDEX CONDITION
# =============================================================================

print()
print("[15] INTRINSIC INDEX CONDITION")
print("-" * 78)

print("""
  Suppose the genuine homogeneous sequence has:

      R_n = F_(n+1)/F_n.

  If one index step corresponds to the unit KAPPA translation,
  then the exact required condition is:

      R_(n+1) - R_n = 2.

  More generally, a step of size h requires:

      R_(n+h) - R_n = 2h.

  This is now an exact condition on the ORIGINAL layer sequence,
  independent of N, S, p, q, or any factorization.
""")


# =============================================================================
# [16] GENERAL AFFINE RATIO LADDER
# =============================================================================

print()
print("[16] GENERAL AFFINE RATIO LADDER")
print("-" * 78)

alpha, beta = sp.symbols("alpha beta")

# General candidate:
#
#     R_n = alpha*n + beta + S.
#
# Requiring one index step to shift the ratio by 2 determines alpha.

R_ladder = sp.expand(
    S + alpha*n + beta
)

R_ladder_next = sp.expand(
    R_ladder.subs(n, n + 1)
)

difference_R = sp.expand(
    R_ladder_next - R_ladder
)

print("  R_n =")
sp.pprint(R_ladder)

print()
print("  R_(n+1)-R_n =")
sp.pprint(difference_R)

sol_alpha = sp.solve(
    sp.Eq(difference_R, 2),
    alpha,
    dict=True,
)

print()
print("  solutions for unit KAPPA step:")
print(f"    {sol_alpha}")

cert(
    "unit step requires alpha=2",
    difference_R.subs(alpha, 2),
    2,
)


# =============================================================================
# [17] CALIBRATE AGAINST THE KNOWN F2/F3 RATIO
# =============================================================================

print()
print("[17] INDEX CALIBRATION")
print("-" * 78)

# Known:
#
#     R_2 = S+1.
#
# For:
#
#     R_n = S+2n+beta
#
# require:
#
#     S+4+beta = S+1
#
# hence beta=-3.

beta_solution = sp.solve(
    sp.Eq(
        (S + 2*2 + beta),
        S + 1
    ),
    beta,
    dict=True,
)

print("  calibration condition:")
sp.pprint(sp.Eq(S + 4 + beta, S + 1))

print()
print("  beta solution:")
print(f"    {beta_solution}")


R_canonical = sp.expand(
    S + 2*n - 3
)

print()
print("  calibrated candidate:")
sp.pprint(R_canonical)

cert(
    "R_2=S+1",
    R_canonical.subs(n, 2),
    S + 1,
)

cert(
    "R_(n+1)-R_n=2",
    R_canonical.subs(n, n + 1) - R_canonical,
    2,
)


# =============================================================================
# [18] HYPOTHETICAL HOMOGENEOUS RATIO LAW
# =============================================================================

print()
print("[18] HYPOTHETICAL INDEX LAW")
print("-" * 78)

print("""
  The UNIQUE affine ratio law calibrated to the known layer is:

      R_n = S + 2n - 3.

  Therefore it predicts:

      F_(n+1) / F_n = S + 2n - 3.

  The first known case is:

      F3/F2 = S+1.

  The next would necessarily be:

      F4/F3 = S+3.

  Then:

      F5/F4 = S+5,

      F6/F5 = S+7,

  etc.

  This is NOT assumed to be true.

  It is the exact prediction that the original homogeneous-layer
  construction must satisfy if its layer index is the same
  translation coordinate discovered above.
""")


# =============================================================================
# [19] IMPLIED HOMOGENEOUS LAYER RECURRENCE
# =============================================================================

print()
print("[19] IMPLIED LAYER RECURRENCE")
print("-" * 78)

# Symbolically define the normalized recurrence:
#
#     F_(n+1) = (S+2n-3) F_n.

Fn = sp.symbols("F_n")

Fn1_pred = sp.expand(
    R_canonical * Fn
)

print("  F_(n+1) =")
sp.pprint(Fn1_pred)


# =============================================================================
# [20] EXPLICIT FIRST TERMS
# =============================================================================

print()
print("[20] FIRST TERMS IMPLIED BY THE INDEX LAW")
print("-" * 78)

F2_symbol = sp.symbols("F2_symbol")

F = {
    2: F2_symbol
}

for idx in range(2, 7):
    ratio_idx = sp.expand(
        R_canonical.subs(n, idx)
    )
    F[idx + 1] = sp.factor(
        ratio_idx * F[idx]
    )

for idx in range(2, 8):
    print(f"  F{idx} =")
    sp.pprint(F[idx])
    print()


# =============================================================================
# [21] VERIFY KNOWN F3
# =============================================================================

print()
print("[21] KNOWN F3 CONSISTENCY")
print("-" * 78)

cert(
    "predicted F3=(S+1)F2",
    F[3],
    (S + 1)*F2_symbol,
)


# =============================================================================
# [22] PREDICTED F4/F3 RATIO
# =============================================================================

print()
print("[22] FIRST UNKNOWN LAYER PREDICTION")
print("-" * 78)

R3_pred = sp.factor(
    F[4] / F[3]
)

print("  predicted F4/F3 =")
sp.pprint(R3_pred)

cert(
    "F4/F3=S+3",
    R3_pred,
    S + 3,
)


# =============================================================================
# [23] TRANSLATION ACTION ON THE RATIO LADDER
# =============================================================================

print()
print("[23] TRANSLATION ACTION ON RATIO LADDER")
print("-" * 78)

R_n_shifted = sp.expand(
    R_canonical + 2*h
)

print("  R_n' =")
sp.pprint(R_n_shifted)

cert(
    "ratio ladder shift",
    R_n_shifted - R_canonical,
    2*h,
)


# =============================================================================
# [24] INDEX SHIFT EQUIVALENCE
# =============================================================================

print()
print("[24] INDEX SHIFT EQUIVALENCE")
print("-" * 78)

# Since:
#
#     R_(n+h) = S + 2(n+h)-3
#
# this is exactly:
#
#     R_n + 2h.

cert(
    "R_(n+h)=R_n+2h",
    R_canonical.subs(n, n + h),
    R_canonical + 2*h,
)


# =============================================================================
# [25] TRANSLATION AS INDEX SHIFT
# =============================================================================

print()
print("[25] TRANSLATION AS INDEX SHIFT")
print("-" * 78)

print("""
  Under the calibrated ratio law:

      R_n = S+2n-3,

  the discovered translation

      R -> R+2h

  is exactly the ordinary layer-index shift

      n -> n+h.

  Therefore:

      KAPPA translation parameter h
          =
      homogeneous-layer index displacement h.

  This is the precise upstream identification being tested.
""")


# =============================================================================
# [26] IMPORTANT NEGATIVE CONTROL
# =============================================================================

print()
print("[26] INDEX-LAW OBSTRUCTION TEST")
print("-" * 78)

#
# We deliberately do NOT claim F4/F3=S+3.
#
# Instead we state the condition that must be verified from the
# actual homogeneous-layer construction.
#

print("""
  The following has NOT been proven by symbolic reconstruction:

      F4/F3 = S+3.

  It is only a necessary prediction of the translation-conjugacy
  hypothesis.

  Therefore the real next check against the original construction is:

      actual(F4/F3) - (S+3) == 0.

  Likewise:

      actual(F5/F4) - (S+5) == 0.

  and generally:

      actual(F_(n+1)/F_n) - (S+2n-3) == 0.
""")


# =============================================================================
# [27] IF THE INDEX LAW HOLDS
# =============================================================================

print()
print("[27] CONSEQUENCE IF THE INDEX LAW HOLDS")
print("-" * 78)

print("""
  If the genuine homogeneous construction satisfies

      F_(n+1)/F_n = S+2n-3,

  then:

      F3/F2-1 = S

  is not merely an algebraic accident.

  It is the n=2 member of a translated ratio ladder.

  The layer index itself becomes the KAPPA translation coordinate.

  Then:

      n -> n+h

  induces:

      S -> S+2h

  and the associated product coordinate must evolve as:

      N -> N+hS+h^2.

  At that point the homogeneous-layer index has been identified
  with the historical/KAPPA translation operator.
""")


# =============================================================================
# [28] FINAL STRUCTURAL CERTIFICATE
# =============================================================================

print()
print("[28] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print("""
  Known exact fact:

      F3/F2 = S+1.

  Exact observable coordinate:

      R_n = F_(n+1)/F_n.

  Exact required unit-step law:

      R_(n+1)-R_n = 2.

  Unique affine calibrated law:

      R_n = S+2n-3.

  Therefore the decisive unresolved identity is:

      F_(n+1)/F_n = S+2n-3.

  The first unknown test case is:

      F4/F3 = S+3.

  This experiment does not assume that identity.

  It isolates it as the exact statement that the original
  homogeneous-layer construction must satisfy for the layer
  index itself to be the KAPPA translation coordinate.
""")


# =============================================================================
# FINAL AUDIT
# =============================================================================

print()
print("=" * 78)
print("EXPERIMENT 547 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
This experiment has reduced the upstream problem to one concrete
identity involving the ACTUAL homogeneous layers:

    R_n = F_(n+1)/F_n.

The KAPPA translation hypothesis requires:

    R_(n+1)-R_n = 2.

Using the already verified n=2 relation:

    F3/F2 = S+1,

the unique affine law is:

    F_(n+1)/F_n = S+2n-3.

Therefore the next experiment should use the ORIGINAL formulas
for F4 and F5, if they exist in the earlier homogeneous-layer
construction, and test exactly:

    F4/F3 == S+3

    F5/F4 == S+5

without introducing p,q or numerical fitting.

If those hold, the layer index itself is the translation operator.
If they fail, the failure identifies exactly where the homogeneous
layer hierarchy departs from the KAPPA translation orbit.
""")
