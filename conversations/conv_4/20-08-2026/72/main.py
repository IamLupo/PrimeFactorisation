#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 548
# ==============================================================================
# FIRST UNKNOWN HOMOGENEOUS-LAYER RATIO AUDIT
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 548 START")
print("=" * 78)
print("FIRST UNKNOWN HOMOGENEOUS-LAYER RATIO AUDIT")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

N, S = sp.symbols("N S")
n, h = sp.symbols("n h")
z = sp.symbols("z")

# These are intentionally placeholders.
# Replace ONLY these two definitions with the genuine historical formulas
# for F4 and F5 when those formulas are available.
F4_actual = sp.Symbol("F4_actual", nonzero=True)
F5_actual = sp.Symbol("F5_actual", nonzero=True)

A, B = sp.symbols("A B", nonzero=True)

failures = 0


# ==============================================================================
# HELPERS
# ==============================================================================

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


# ==============================================================================
# [1] ESTABLISHED HOMOGENEOUS LAYERS
# ==============================================================================

print("[1] ESTABLISHED HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = sp.expand(
    6 * N - S**2 + S
)

F3 = sp.expand(
    (S + 1) * F2
)

print("  F2 =")
sp.pprint(F2)

print()
print("  F3 =")
sp.pprint(F3)

cert(
    "F3=(S+1)F2",
    F3,
    (S + 1) * F2,
)


# ==============================================================================
# [2] ESTABLISHED RATIO
# ==============================================================================

print()
print("[2] ESTABLISHED RATIO")
print("-" * 78)

R2 = sp.factor(
    F3 / F2
)

print("  R2 = F3/F2")
sp.pprint(R2)

cert(
    "R2=S+1",
    R2,
    S + 1,
)


# ==============================================================================
# [3] CALIBRATED INDEX LADDER
# ==============================================================================

print()
print("[3] CALIBRATED INDEX LADDER")
print("-" * 78)

R_law = sp.expand(
    S + 2*n - 3
)

print("  proposed R_n =")
sp.pprint(R_law)

cert(
    "R2=S+1",
    R_law.subs(n, 2),
    S + 1,
)

cert(
    "R_(n+1)-R_n=2",
    R_law.subs(n, n + 1) - R_law,
    2,
)


# ==============================================================================
# [4] PREDICTED FIRST UNKNOWN RATIO
# ==============================================================================

print()
print("[4] FIRST UNKNOWN RATIO")
print("-" * 78)

R3_pred = sp.expand(
    R_law.subs(n, 3)
)

print("  predicted F4/F3 =")
sp.pprint(R3_pred)

cert(
    "R3=S+3",
    R3_pred,
    S + 3,
)


# ==============================================================================
# [5] PREDICTED SECOND UNKNOWN RATIO
# ==============================================================================

print()
print("[5] SECOND UNKNOWN RATIO")
print("-" * 78)

R4_pred = sp.expand(
    R_law.subs(n, 4)
)

print("  predicted F5/F4 =")
sp.pprint(R4_pred)

cert(
    "R4=S+5",
    R4_pred,
    S + 5,
)


# ==============================================================================
# [6] PREDICTED F4
# ==============================================================================

print()
print("[6] PREDICTED F4")
print("-" * 78)

F4_pred = sp.factor(
    F3 * R3_pred
)

print("  F4_pred =")
sp.pprint(F4_pred)

cert(
    "F4 prediction",
    F4_pred,
    (S + 3) * F3,
)


# ==============================================================================
# [7] PREDICTED F5
# ==============================================================================

print()
print("[7] PREDICTED F5")
print("-" * 78)

F5_pred = sp.factor(
    F4_pred * R4_pred
)

print("  F5_pred =")
sp.pprint(F5_pred)

cert(
    "F5 prediction",
    F5_pred,
    (S + 5) * F4_pred,
)


# ==============================================================================
# [8] ACTUAL F4 PLACEHOLDER AUDIT
# ==============================================================================

print()
print("[8] ACTUAL F4 AUDIT")
print("-" * 78)

R3_actual = sp.factor(
    F4_actual / F3
)

print("  R3_actual = F4/F3")
sp.pprint(R3_actual)

E4 = sp.factor(
    R3_actual - (S + 3)
)

print()
print("  residual E4 =")
sp.pprint(E4)

print("""
  Required exact identity:

      E4 = 0.

  This remains untested until F4_actual is replaced by the
  genuine homogeneous-layer formula.
""")


# ==============================================================================
# [9] ACTUAL F5 PLACEHOLDER AUDIT
# ==============================================================================

print()
print("[9] ACTUAL F5 AUDIT")
print("-" * 78)

R4_actual = sp.factor(
    F5_actual / F4_actual
)

print("  R4_actual = F5/F4")
sp.pprint(R4_actual)

E5 = sp.factor(
    R4_actual - (S + 5)
)

print()
print("  residual E5 =")
sp.pprint(E5)

print("""
  Required exact identity:

      E5 = 0.

  This remains untested until F5_actual is replaced by the
  genuine homogeneous-layer formula.
""")


# ==============================================================================
# [10] RATIO DIFFERENCE FORM
# ==============================================================================

print()
print("[10] RATIO DIFFERENCE FORM")
print("-" * 78)

D3 = sp.factor(
    R3_actual - R2
)

D4 = sp.factor(
    R4_actual - R3_actual
)

print("  R3-R2 =")
sp.pprint(D3)

print()
print("  R4-R3 =")
sp.pprint(D4)

print("""
  The translation hypothesis is exactly:

      R3-R2 = 2

      R4-R3 = 2.
""")


# ==============================================================================
# [11] SECOND DIFFERENCE OF THE RATIO LADDER
# ==============================================================================

print()
print("[11] RATIO-LADDER CURVATURE")
print("-" * 78)

ratio_second_difference = sp.factor(
    R4_actual - 2*R3_actual + R2
)

print("  R4 - 2 R3 + R2 =")
sp.pprint(ratio_second_difference)

print("""
  Under the proposed affine ratio ladder:

      R_(n+2) - 2R_(n+1) + R_n = 0.
""")


# ==============================================================================
# [12] OBSERVABLE-ONLY PREDICTIONS
# ==============================================================================

print()
print("[12] OBSERVABLE-ONLY PREDICTIONS")
print("-" * 78)

W = sp.factor(
    B / A
)

F4_obs_pred = sp.factor(
    B * (W + 2)
)

F5_obs_pred = sp.factor(
    F4_obs_pred * (W + 4)
)

print("  W = F3/F2")
sp.pprint(W)

print()
print("  F4 predicted from F2,F3:")
sp.pprint(F4_obs_pred)

print()
print("  F5 predicted from F2,F3:")
sp.pprint(F5_obs_pred)

cert(
    "observable F4 formula",
    F4_obs_pred,
    B * (B/A + 2),
)

cert(
    "observable F5 formula",
    F5_obs_pred,
    F4_obs_pred * (B/A + 4),
)


# ==============================================================================
# [13] INDEX-SHIFT OPERATOR
# ==============================================================================

print()
print("[13] INDEX-SHIFT OPERATOR")
print("-" * 78)

R_n_shift = sp.expand(
    R_law.subs(n, n + h)
)

print("  R_(n+h) =")
sp.pprint(R_n_shift)

cert(
    "R_(n+h)=R_n+2h",
    R_n_shift,
    R_law + 2*h,
)


# ==============================================================================
# [14] PREDICTED LAYER RECURRENCE
# ==============================================================================

print()
print("[14] PREDICTED LAYER RECURRENCE")
print("-" * 78)

F_symbol = sp.Symbol("F_n")

F_next = sp.factor(
    R_law * F_symbol
)

print("  F_(n+1) =")
sp.pprint(F_next)

print("""
  The hypothesis predicts:

      F_(n+1) = (S+2n-3) F_n.
""")


# ==============================================================================
# [15] FIRST TERMS OF THE PREDICTED LADDER
# ==============================================================================

print()
print("[15] PREDICTED LAYER LADDER")
print("-" * 78)

F = {
    2: sp.Symbol("F2")
}

for k in range(2, 7):
    ratio = sp.expand(
        R_law.subs(n, k)
    )
    F[k + 1] = sp.factor(
        ratio * F[k]
    )

for k in range(2, 8):
    print(f"  F{k} =")
    sp.pprint(F[k])
    print()


# ==============================================================================
# [16] TRANSLATION OF RATIO COORDINATE
# ==============================================================================

print()
print("[16] RATIO TRANSLATION")
print("-" * 78)

R_shift = sp.expand(
    R_law + 2*h
)

print("  R_n' =")
sp.pprint(R_shift)

cert(
    "R_n' - R_n = 2h",
    R_shift - R_law,
    2*h,
)


# ==============================================================================
# [17] MOVING COORDINATE
# ==============================================================================

print()
print("[17] MOVING COORDINATE")
print("-" * 78)

W_n = sp.expand(
    R_law
)

print("  W_n =")
sp.pprint(W_n)

print("""
  Since:

      W_(n+h) - W_n = 2h,

  the layer index is an affine coordinate in ratio space.
""")


# ==============================================================================
# [18] S RECOVERY AT n=2
# ==============================================================================

print()
print("[18] S RECOVERY")
print("-" * 78)

S_from_ratio = sp.expand(
    R2 - 1
)

cert(
    "S=R2-1",
    S_from_ratio,
    S,
)


# ==============================================================================
# [19] N CHANNEL
# ==============================================================================

print()
print("[19] N CHANNEL")
print("-" * 78)

N_from_F2S = sp.factor(
    (F2 + S**2 - S) / 6
)

print("  N =")
sp.pprint(N_from_F2S)

cert(
    "N recovery from F2,S",
    N_from_F2S,
    N,
)


# ==============================================================================
# [20] GENERATING QUADRATIC
# ==============================================================================

print()
print("[20] GENERATING KAPPA QUADRATIC")
print("-" * 78)

Q = sp.expand(
    z**2 - S*z + N
)

print("  Q(z) =")
sp.pprint(Q)

disc = sp.factor(
    sp.discriminant(Q, z)
)

print()
print("  discriminant(Q) =")
sp.pprint(disc)

cert(
    "discriminant = S^2-4N",
    disc,
    S**2 - 4*N,
)


# ==============================================================================
# [21] HISTORICAL TRANSLATION CONNECTION
# ==============================================================================

print()
print("[21] HISTORICAL TRANSLATION CONNECTION")
print("-" * 78)

print("""
  If the ratio ladder is genuine:

      R_n = S+2n-3

  then an index displacement h gives:

      R_(n+h)=R_n+2h.

  This is exactly the same affine motion as the historical/KAPPA
  translation coordinate.

  Therefore:

      homogeneous layer index
          <-> KAPPA translation parameter.
""")


# ==============================================================================
# [22] IMPORTANT VALIDATION BOUNDARY
# ==============================================================================

print()
print("[22] VALIDATION BOUNDARY")
print("-" * 78)

print("""
  This experiment deliberately does NOT claim that:

      F4/F3 = S+3

  or:

      F5/F4 = S+5.

  Those are still hypotheses.

  To validate them, the placeholders:

      F4_actual
      F5_actual

  must be replaced with the genuine formulas generated by the
  original homogeneous-layer construction.

  Once replaced, the decisive tests are:

      F4_actual/F3 - F3/F2 = 2

      F5_actual/F4_actual - F4_actual/F3 = 2

  These tests do not require p, q, or numerical fitting.
""")


# ==============================================================================
# [23] FINAL STRUCTURAL CERTIFICATE
# ==============================================================================

print()
print("[23] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print("""
  Established:

      F3/F2 = S+1.

  Proposed intrinsic index law:

      R_n = F_(n+1)/F_n

      R_(n+1)-R_n = 2.

  Therefore necessarily:

      R_n = S+2n-3.

  First unknown consequences:

      F4/F3 = S+3

      F5/F4 = S+5.

  Observable-only form:

      F4 = F3(F3/F2+2)

      F5 = F4(F3/F2+4).

  The strongest upstream test is therefore:

      F4/F3 - F3/F2 = 2

      F5/F4 - F4/F3 = 2.

  No further downstream identities are required until the actual
  F4 and F5 formulas are tested.
""")


# ==============================================================================
# FINAL
# ==============================================================================

print()
print("=" * 78)
print("EXPERIMENT 548 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The next experiment should NOT derive another synthetic identity.

It should use the genuine original homogeneous-layer formulas for
F4 and F5.

Then test only:

    F4/F3 - F3/F2

and

    F5/F4 - F4/F3.

The predicted result is:

    2

and:

    2.

If both equal 2 exactly, the homogeneous-layer index itself has
the affine translation law:

    R_(n+1)=R_n+2.

That would identify the layer index as the same translation
coordinate underlying the historical conic and KAPPA quadratic.

If either fails, the residual tells us exactly how the actual
homogeneous hierarchy differs from the proposed translation
representation.
""")