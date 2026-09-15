#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 549
# ==============================================================================
# DIVISION-FREE HOMOGENEOUS RATIO LADDER / INTRINSIC INDEX OPERATOR
#
# Known:
#
#     F2 = 6N - S^2 + S
#     F3 = (S+1)F2
#
# Hypothesis:
#
#     R_n = F_(n+1)/F_n
#
#     R_(n+1)-R_n = 2
#
# Instead of using ratios, this experiment derives the completely
# division-free polynomial identities equivalent to that statement.
#
# This is important because the genuine homogeneous construction may
# have polynomial formulas even when F_(n+1)/F_n is rational.
#
# Main identities:
#
#   1. F_(n+2) F_n - F_(n+1)^2 - 2 F_(n+1) F_n = 0
#
#   2. H_n = F_(n+1) - (2n-3) F_n
#
#      H_n/F_n = S
#
#   3. H_(n+1) F_n - H_n F_(n+1) = 0
#
#      This expresses constancy of the hidden S-channel without
#      explicitly inserting S.
#
# The experiment also derives the first constraints on F4 and F5.
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 549 START")
print("=" * 78)
print("DIVISION-FREE HOMOGENEOUS RATIO LADDER / INTRINSIC INDEX OPERATOR")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

N, S = sp.symbols("N S")
F2, F3 = sp.symbols("F2 F3", nonzero=True)
F4, F5, F6 = sp.symbols("F4 F5 F6", nonzero=True)

n = sp.symbols("n", integer=True)
k = sp.symbols("k", integer=True)

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


# ==============================================================================
# [1] KNOWN HOMOGENEOUS LAYERS
# ==============================================================================

print("[1] KNOWN HOMOGENEOUS LAYERS")
print("-" * 78)

F2_exact = 6*N - S**2 + S
F3_exact = (S + 1)*F2_exact

print("  F2 =")
sp.pprint(F2_exact)

print()
print("  F3 =")
sp.pprint(sp.expand(F3_exact))

cert(
    "F3=(S+1)F2",
    F3_exact,
    (S + 1)*F2_exact,
)


# ==============================================================================
# [2] KNOWN RATIO
# ==============================================================================

print()
print("[2] KNOWN RATIO")
print("-" * 78)

R2 = sp.factor(F3_exact / F2_exact)

print("  R2 = F3/F2")
sp.pprint(R2)

cert(
    "R2=S+1",
    R2,
    S + 1,
)


# ==============================================================================
# [3] DIVISION-FREE FORM OF THE UNIT RATIO STEP
# ==============================================================================

print()
print("[3] DIVISION-FREE UNIT RATIO STEP")
print("-" * 78)

#
# R_(n+1)-R_n = 2
#
# becomes:
#
# F_(n+2)/F_(n+1) - F_(n+1)/F_n = 2
#
# multiply by F_(n+1)F_n:
#
# F_(n+2)F_n - F_(n+1)^2 = 2F_(n+1)F_n.
#
# Therefore:
#
# F_(n+2)F_n - F_(n+1)^2 - 2F_(n+1)F_n = 0.
#

F_n, F_np1, F_np2 = sp.symbols(
    "F_n F_np1 F_np2",
    nonzero=True,
)

ladder_identity = sp.expand(
    F_np2*F_n
    - F_np1**2
    - 2*F_np1*F_n
)

print("  universal division-free identity:")
print()
print("      F_(n+2) F_n - F_(n+1)^2 - 2F_(n+1)F_n = 0")
print()
sp.pprint(ladder_identity)


# ==============================================================================
# [4] FIRST UNKNOWN DIVISION-FREE TEST
# ==============================================================================

print()
print("[4] FIRST UNKNOWN DIVISION-FREE TEST")
print("-" * 78)

test_F4 = sp.expand(
    F4*F2 - F3**2 - 2*F3*F2
)

print("  F4*F2 - F3^2 - 2F3F2 =")
sp.pprint(test_F4)

print()
print("""
  Genuine F4 must satisfy exactly:

      F4*F2 - F3^2 - 2F3F2 = 0.
""")


# ==============================================================================
# [5] SECOND UNKNOWN DIVISION-FREE TEST
# ==============================================================================

print()
print("[5] SECOND UNKNOWN DIVISION-FREE TEST")
print("-" * 78)

test_F5 = sp.expand(
    F5*F3 - F4**2 - 2*F4*F3
)

print("  F5*F3 - F4^2 - 2F4F3 =")
sp.pprint(test_F5)

print()
print("""
  Genuine F5 must satisfy exactly:

      F5*F3 - F4^2 - 2F4F3 = 0.
""")


# ==============================================================================
# [6] THIRD UNKNOWN DIVISION-FREE TEST
# ==============================================================================

print()
print("[6] THIRD UNKNOWN DIVISION-FREE TEST")
print("-" * 78)

test_F6 = sp.expand(
    F6*F4 - F5**2 - 2*F5*F4
)

print("  F6*F4 - F5^2 - 2F5F4 =")
sp.pprint(test_F6)

print()
print("""
  Once F6 exists, the same identity continues:

      F6*F4 - F5^2 - 2F5F4 = 0.
""")


# ==============================================================================
# [7] S-CHANNEL WITHOUT S
# ==============================================================================

print()
print("[7] DIVISION-FREE S-CHANNEL")
print("-" * 78)

#
# Under:
#
#   R_n = S + 2n - 3
#
# we have:
#
#   F_(n+1) = (S+2n-3)F_n.
#
# Rearrange:
#
#   F_(n+1) - (2n-3)F_n = S F_n.
#
# Define:
#
#   H_n = F_(n+1) - (2n-3)F_n.
#
# Then:
#
#   H_n/F_n = S.
#

H_n = sp.expand(
    F_np1 - (2*n - 3)*F_n
)

print("  H_n = F_(n+1) - (2n-3)F_n")
sp.pprint(H_n)

print()
print("  expected:")
print("      H_n = S*F_n")


# ==============================================================================
# [8] KNOWN n=2 CHANNEL
# ==============================================================================

print()
print("[8] KNOWN n=2 CHANNEL")
print("-" * 78)

H2_exact = sp.expand(
    F3_exact - (2*2 - 3)*F2_exact
)

print("  H2 = F3-F2")
sp.pprint(H2_exact)

cert(
    "H2=S*F2",
    H2_exact,
    S*F2_exact,
)


# ==============================================================================
# [9] CHANNEL CONSTANCY WITHOUT S
# ==============================================================================

print()
print("[9] CHANNEL CONSTANCY WITHOUT S")
print("-" * 78)

#
# H_n/F_n = H_(n+1)/F_(n+1)
#
# is equivalent to:
#
# H_(n+1)F_n - H_nF_(n+1) = 0.
#

H_np1 = sp.expand(
    F_np2 - (2*(n+1) - 3)*F_np1
)

channel_constancy = sp.expand(
    H_np1*F_n - H_n*F_np1
)

print("  H_(n+1)F_n - H_nF_(n+1) =")
sp.pprint(channel_constancy)

print()
print("""
  This is a completely division-free expression of:

      H_n/F_n = constant.

  Under the KAPPA translation hypothesis that constant is S.
""")


# ==============================================================================
# [10] EQUIVALENCE OF THE TWO CONDITIONS
# ==============================================================================

print()
print("[10] EQUIVALENCE OF RATIO LAW AND CHANNEL LAW")
print("-" * 78)

R_general = sp.symbols("R_general")
R_next = sp.symbols("R_next")

ratio_condition = sp.expand(
    R_next - R_general - 2
)

print("  ratio condition:")
sp.pprint(ratio_condition)

H_ratio = sp.expand(
    (R_general - (2*n - 3))
)

H_next_ratio = sp.expand(
    (R_next - (2*(n+1) - 3))
)

channel_condition = sp.expand(
    H_next_ratio - H_ratio
)

print()
print("  shifted channel condition:")
sp.pprint(channel_condition)

cert(
    "equivalent increment condition",
    channel_condition,
    ratio_condition,
)


# ==============================================================================
# [11] PREDICTED F4 FROM DIVISION-FREE CONDITION
# ==============================================================================

print()
print("[11] F4 CONSTRAINT")
print("-" * 78)

F4_solution = sp.factor(
    (F3**2 + 2*F3*F2) / F2
)

print("  Solving the first identity for F4:")
print()
print("      F4 =")
sp.pprint(F4_solution)

print()
print("  This is equivalent to:")
print()
print("      F4 = F3(F3/F2 + 2)")


# ==============================================================================
# [12] SUBSTITUTE THE KNOWN F2,F3
# ==============================================================================

print()
print("[12] F4 PREDICTION FROM KNOWN LAYERS")
print("-" * 78)

F4_pred = sp.factor(
    F4_solution.subs({
        F2: F2_exact,
        F3: F3_exact,
    })
)

print("  F4_pred =")
sp.pprint(sp.expand(F4_pred))

F4_expected = sp.factor(
    (S + 3)*F3_exact
)

cert(
    "F4=(S+3)F3",
    F4_pred,
    F4_expected,
)


# ==============================================================================
# [13] F5 CONSTRAINT
# ==============================================================================

print()
print("[13] F5 CONSTRAINT")
print("-" * 78)

F5_solution = sp.factor(
    (F4**2 + 2*F4*F3) / F3
)

print("  Solving the second identity for F5:")
print()
print("      F5 =")
sp.pprint(F5_solution)


# ==============================================================================
# [14] SUBSTITUTE F4 PREDICTION
# ==============================================================================

print()
print("[14] F5 PREDICTION FROM KNOWN LAYERS")
print("-" * 78)

F5_pred = sp.factor(
    F5_solution.subs({
        F3: F3_exact,
        F4: F4_pred,
    })
)

print("  F5_pred =")
sp.pprint(sp.expand(F5_pred))

F5_expected = sp.factor(
    (S + 5)*F4_expected
)

cert(
    "F5=(S+5)F4",
    F5_pred,
    F5_expected,
)


# ==============================================================================
# [15] CLOSED OBSERVABLE FORM
# ==============================================================================

print()
print("[15] CLOSED OBSERVABLE FORM")
print("-" * 78)

A, B = sp.symbols("A B", nonzero=True)

F4_AB = sp.factor(
    (B**2 + 2*A*B) / A
)

F5_AB = sp.factor(
    (F4_AB**2 + 2*B*F4_AB) / B
)

print("  F4(F2,F3) =")
sp.pprint(F4_AB)

print()
print("  F5(F2,F3) =")
sp.pprint(F5_AB)

print()
print("""
  The division-free recurrence is therefore equivalent to the
  rational observable evolution:

      F4 = F3(F3/F2 + 2)

      F5 = F4(F4/F3 + 2).
""")


# ==============================================================================
# [16] CASORATIAN-LIKE DEFECT
# ==============================================================================

print()
print("[16] LADDER DEFECT")
print("-" * 78)

#
# Define:
#
#   E_n = F_(n+2)F_n - F_(n+1)^2 - 2F_(n+1)F_n.
#
# Exact translation law:
#
#   E_n = 0.
#
# This is a natural defect observable for the real sequence.
#

E_n = sp.Function("E")(n)

print("""
  Define the intrinsic defect:

      E_n =
        F_(n+2)F_n
        - F_(n+1)^2
        - 2F_(n+1)F_n.

  The KAPPA translation hypothesis is exactly:

      E_n = 0.
""")


# ==============================================================================
# [17] CHANNEL DEFECT
# ==============================================================================

print()
print("[17] S-CHANNEL DEFECT")
print("-" * 78)

#
# Define:
#
#   C_n = H_(n+1)F_n - H_nF_(n+1).
#
# Then:
#
#   C_n = 0
#
# is equivalent to the constancy of H_n/F_n.
#

print("""
  Define:

      H_n = F_(n+1) - (2n-3)F_n

      C_n =
        H_(n+1)F_n
        - H_nF_(n+1).

  Exact translation structure requires:

      C_n = 0.
""")


# ==============================================================================
# [18] RELATION BETWEEN THE TWO DEFECTS
# ==============================================================================

print()
print("[18] DEFECT EQUIVALENCE")
print("-" * 78)

#
# Explicitly construct C_n in terms of F_n,F_(n+1),F_(n+2).
#

Fn, Fnp1, Fnp2 = sp.symbols(
    "Fn Fnp1 Fnp2",
    nonzero=True,
)

Hn_generic = sp.expand(
    Fnp1 - (2*n - 3)*Fn
)

Hnp1_generic = sp.expand(
    Fnp2 - (2*n - 1)*Fnp1
)

C_generic = sp.factor(
    Hnp1_generic*Fn
    - Hn_generic*Fnp1
)

E_generic = sp.factor(
    Fnp2*Fn
    - Fnp1**2
    - 2*Fnp1*Fn
)

print("  C_n =")
sp.pprint(C_generic)

print()
print("  E_n =")
sp.pprint(E_generic)

cert(
    "C_n-E_n=0",
    C_generic,
    E_generic,
)


# ==============================================================================
# [19] THIS IS THE KEY RESULT
# ==============================================================================

print()
print("[19] KEY OPERATOR IDENTITY")
print("-" * 78)

print("""
  The two apparently different tests are EXACTLY THE SAME:

      ratio increment:

          R_(n+1)-R_n = 2

      channel constancy:

          H_(n+1)/F_(n+1) = H_n/F_n

      division-free defect:

          E_n = 0.

  Therefore only ONE genuine upstream identity has to be
  established from the original homogeneous-layer construction:

      F_(n+2)F_n
        - F_(n+1)^2
        - 2F_(n+1)F_n
        = 0.

  This is stronger experimentally because it contains no
  division and no explicit S.
""")


# ==============================================================================
# [20] INDEX TRANSLATION CONSEQUENCE
# ==============================================================================

print()
print("[20] INDEX TRANSLATION CONSEQUENCE")
print("-" * 78)

print("""
  If E_n = 0 for the genuine homogeneous sequence, then:

      R_(n+1)=R_n+2.

  Hence:

      R_n-2n

  is invariant under the layer-index shift.

  Since:

      R_2=S+1,

  the invariant is:

      R_n-2n = S-3.

  Therefore:

      R_n = S+2n-3.
""")


# ==============================================================================
# [21] SECOND-ORDER NONLINEAR RECURRENCE
# ==============================================================================

print()
print("[21] SECOND-ORDER NONLINEAR RECURRENCE")
print("-" * 78)

print("""
  The upstream identity can be written as:

      F_(n+2)F_n
        = F_(n+1)(F_(n+1)+2F_n).

  Thus:

      F_(n+2)
        = F_(n+1)^2/F_n
          + 2F_(n+1).

  This is a nonlinear second-order recurrence.

  The important point is that it requires no explicit S.

  If the genuine homogeneous hierarchy obeys this recurrence,
  its ratio ladder automatically has the KAPPA translation law.
""")


# ==============================================================================
# [22] INTEGRATED FORM
# ==============================================================================

print()
print("[22] INTEGRATED RATIO FORM")
print("-" * 78)

print("""
  Let:

      R_n = F_(n+1)/F_n.

  Then the nonlinear recurrence becomes:

      R_(n+1)=R_n+2.

  Therefore:

      R_n = 2n + C

  for a layer-independent constant C.

  Calibration at n=2 gives:

      C = S-3.

  Hence:

      R_n=S+2n-3.
""")


# ==============================================================================
# [23] WHAT MUST NOW BE CHECKED AGAINST THE REAL CONSTRUCTION
# ==============================================================================

print()
print("[23] ACTUAL-CONSTRUCTION TEST")
print("-" * 78)

print("""
  Once the genuine formulas are available, DO NOT first compute
  S or p,q.

  Compute directly:

      E_2 =
        F4*F2
        - F3^2
        - 2F3F2.

  Then:

      E_3 =
        F5*F3
        - F4^2
        - 2F4F3.

  Then:

      E_4 =
        F6*F4
        - F5^2
        - 2F5F4.

  The ideal result is:

      E_2 = E_3 = E_4 = ... = 0.

  This would establish the layer-index translation law without
  ever introducing S as an input.
""")


# ==============================================================================
# [24] VALIDATION BOUNDARY
# ==============================================================================

print()
print("[24] VALIDATION BOUNDARY")
print("-" * 78)

print("""
  This experiment proves the equivalence of the proposed
  translation law and its division-free recurrence form.

  It does NOT yet prove that the original homogeneous-layer
  construction satisfies E_n=0.

  That requires the genuine F4/F5/F6 formulas.

  Therefore this experiment deliberately separates:

      algebraic consequence
          from
      empirical/original-construction validation.
""")


# ==============================================================================
# FINAL AUDIT
# ==============================================================================

print()
print("=" * 78)
print("EXPERIMENT 549 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The next experiment should finally stop predicting F4/F5 individually.

The strongest upstream identity is now:

    F_(n+2)F_n
      - F_(n+1)^2
      - 2F_(n+1)F_n
      = 0.

This is equivalent to:

    F_(n+2)/F_(n+1)
      - F_(n+1)/F_n
      = 2,

but is completely division-free.

Therefore the next genuine homogeneous-layer formulas should be
substituted directly into:

    E_2 = F4F2-F3^2-2F3F2

    E_3 = F5F3-F4^2-2F4F3

    E_4 = F6F4-F5^2-2F5F4.

If these vanish identically, then the homogeneous-layer index
itself has the affine KAPPA translation law.

That would be the first real upstream test of the hypothesis.
""")
