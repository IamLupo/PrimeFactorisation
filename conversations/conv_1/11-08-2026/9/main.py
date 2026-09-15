import sympy as sp
from itertools import combinations
from sympy import primerange


# =============================================================================
# KAPPA 3-BODY / 4-BODY — PAPER FORMULA DIAGNOSTIC
# =============================================================================

print("=" * 80)
print("KAPPA STRUCTURE — EXACT PAPER-FORMULA FAILURE ANALYSIS")
print("=" * 80)


# -----------------------------------------------------------------------------
# SYMBOLS
# -----------------------------------------------------------------------------

e1, e2, e3, e4, n = sp.symbols("e1 e2 e3 e4 n")
p, q, r, s = sp.symbols("p q r s")


# =============================================================================
# TEST 1: VERIFIED 3-BODY KAPPA
# =============================================================================

print("\n" + "=" * 80)
print("TEST 1: 3-BODY KAPPA")
print("=" * 80)

kappa3 = (
    -e1**2
    + e1*e2
    + e1*e3
    + e1
    - e2**2
    + e2*e3
    + e2
    - 3*e3
)

print("\nkappa3 =")
print(kappa3)


# =============================================================================
# TEST 2: VERIFIED 4-BODY KAPPA
# =============================================================================

print("\n" + "=" * 80)
print("TEST 2: 4-BODY KAPPA")
print("=" * 80)

kappa4 = (
    -e1**2
    + e1*e2
    + e1*e3
    - 2*e1*e4
    + e1
    - e2**2
    + e2*e3
    + e2*e4
    + e2
    - e3**2
    + e3*e4
    - 2*e3
)

print("\nkappa4 =")
print(kappa4)

kappa4_n = sp.expand(kappa4.subs(e4, n))

print("\nkappa4(e4=n) =")
print(kappa4_n)


# =============================================================================
# TEST 3: EXACT 3 -> 4 CORRECTION
# =============================================================================

print("\n" + "=" * 80)
print("TEST 3: EXACT 3-BODY -> 4-BODY CORRECTION")
print("=" * 80)

correction = sp.expand(kappa4_n - kappa3)

print("\nkappa4 - kappa3 =")
print(correction)

correction_expected = (
    n * (-2*e1 + e2 + e3)
    - e3 * (e3 - 1)
)

print("\nExpected compact correction =")
print(correction_expected)

check = sp.simplify(correction - correction_expected)

print("\nCheck:")
print(check)

assert check == 0

print("\nPASS: exact correction confirmed.")


# =============================================================================
# TEST 4: SEPARATE n-DEPENDENT AND n-INDEPENDENT PARTS
# =============================================================================

print("\n" + "=" * 80)
print("TEST 4: CORRECTION DECOMPOSITION")
print("=" * 80)

correction_poly = sp.Poly(correction, n)

c0 = correction_poly.coeff_monomial(n**0)
c1 = correction_poly.coeff_monomial(n**1)

print("\nn^0 coefficient:")
print(sp.factor(c0))

print("\nn^1 coefficient:")
print(sp.factor(c1))

assert sp.expand(c0 + c1*n - correction) == 0

print("\nTherefore:")
print("  kappa4 - kappa3 =")
print("      n*(-2*e1 + e2 + e3)")
print("      - e3*(e3 - 1)")


# =============================================================================
# TEST 5: NORMALIZED CUBIC PRODUCT — 3 BODY
# =============================================================================

print("\n" + "=" * 80)
print("TEST 5: NORMALIZED CUBIC PRODUCT — 3 BODY")
print("=" * 80)

X3 = (p**3 + 1) * (q**3 + 1) * (r**3 + 1)
M13 = (p + 1) * (q + 1) * (r + 1)

X3_norm_raw = sp.cancel(X3 / M13)

e1_3 = p + q + r
e2_3 = p*q + p*r + q*r
e3_3 = p*q*r

X3_norm_expected = (
    e1_3**2
    - e1_3*e2_3
    - e1_3*e3_3
    - e1_3
    + e2_3**2
    - e2_3*e3_3
    - e2_3
    + e3_3**2
    + 2*e3_3
    + 1
)

print("\nX3/M1 =")
print(sp.expand(X3_norm_raw))

print("\nCheck:")
print(sp.simplify(X3_norm_raw - X3_norm_expected))

assert sp.simplify(X3_norm_raw - X3_norm_expected) == 0

print("\nPASS: 3-body normalized cubic product confirmed.")


# =============================================================================
# TEST 6: NORMALIZED CUBIC PRODUCT — 4 BODY
# =============================================================================

print("\n" + "=" * 80)
print("TEST 6: NORMALIZED CUBIC PRODUCT — 4 BODY")
print("=" * 80)

X4 = (
    (p**3 + 1)
    * (q**3 + 1)
    * (r**3 + 1)
    * (s**3 + 1)
)

M14 = (p + 1) * (q + 1) * (r + 1) * (s + 1)

X4_norm_raw = sp.cancel(X4 / M14)

e1_4 = p + q + r + s
e2_4 = (
    p*q + p*r + p*s
    + q*r + q*s + r*s
)
e3_4 = (
    p*q*r + p*q*s
    + p*r*s + q*r*s
)
e4_4 = p*q*r*s

X4_norm_expected = (
    e1_4**2
    - e1_4*e2_4
    - e1_4*e3_4
    + 2*e1_4*e4_4
    - e1_4
    + e2_4**2
    - e2_4*e3_4
    - e2_4*e4_4
    - e2_4
    + e3_4**2
    - e3_4*e4_4
    + 2*e3_4
    + e4_4**2
    - e4_4
    + 1
)

print("\nCheck:")
print(sp.expand(X4_norm_raw - X4_norm_expected))

assert sp.expand(X4_norm_raw - X4_norm_expected) == 0

print("\nPASS: 4-body normalized cubic product confirmed.")


# =============================================================================
# TEST 7: THE NAIVE PAPER EXTENSION
#
# This is the important diagnostic.
#
# If the 3-body formula is mechanically extended to four bodies while
# ignoring the new elementary-symmetric variable e4, the normalized product
# would incorrectly look like:
#
#   e1^2 - e1 e2 - e1 e4 - e1
#   + e2^2 - e2 e4 - e2
#   + e4^2 + 2 e4 + 1
#
# This is exactly the kind of "body-3 formula transplanted to body-4"
# structure we want to test.
# =============================================================================

print("\n" + "=" * 80)
print("TEST 7: NAIVE PAPER BODY-3 -> BODY-4 EXTENSION")
print("=" * 80)

paper_normalized = (
    e1**2
    - e1*e2
    - e1*e4
    - e1
    + e2**2
    - e2*e4
    - e2
    + e4**2
    + 2*e4
    + 1
)

actual_normalized = (
    e1**2
    - e1*e2
    - e1*e3
    + 2*e1*e4
    - e1
    + e2**2
    - e2*e3
    - e2*e4
    - e2
    + e3**2
    - e3*e4
    + 2*e3
    + e4**2
    - e4
    + 1
)

paper_residual = sp.expand(actual_normalized - paper_normalized)

print("\nPaper normalized expression:")
print(paper_normalized)

print("\nActual normalized expression:")
print(actual_normalized)

print("\nExact residual = actual - paper:")
print(sp.factor(paper_residual))

print("\nExpanded residual:")
print(sp.expand(paper_residual))


# =============================================================================
# TEST 8: PAPER KAPPA VS ACTUAL KAPPA
# =============================================================================

print("\n" + "=" * 80)
print("TEST 8: PAPER KAPPA VS ACTUAL KAPPA")
print("=" * 80)

# Since
#
#   X/M1 = e4^2 - e4 + 1 - kappa4
#
# the paper's normalized product implies a corresponding paper kappa.

paper_kappa4 = sp.expand(
    e4**2 - e4 + 1 - paper_normalized
)

print("\nPaper-implied kappa4:")
print(sp.factor(paper_kappa4))

print("\nActual kappa4:")
print(kappa4)

kappa_residual = sp.expand(kappa4 - paper_kappa4)

print("\nExact kappa residual:")
print(sp.factor(kappa_residual))

print("\nExpanded:")
print(kappa_residual)


# =============================================================================
# TEST 9: FACTOR THE PAPER FAILURE
# =============================================================================

print("\n" + "=" * 80)
print("TEST 9: FACTOR THE PAPER FAILURE")
print("=" * 80)

print("\nkappa4 - kappa4_paper =")
print(sp.factor(kappa_residual))

print("\nThis identifies the missing structural terms.")


# =============================================================================
# TEST 10: SUBSTITUTE e4 = n
# =============================================================================

print("\n" + "=" * 80)
print("TEST 10: PAPER FAILURE AFTER e4 = n")
print("=" * 80)

actual_n = sp.expand(kappa4.subs(e4, n))
paper_n = sp.expand(paper_kappa4.subs(e4, n))

residual_n = sp.factor(actual_n - paper_n)

print("\nActual kappa4:")
print(actual_n)

print("\nPaper kappa4:")
print(paper_n)

print("\nResidual:")
print(residual_n)


# =============================================================================
# TEST 11: NUMERICAL TEST
# =============================================================================

print("\n" + "=" * 80)
print("TEST 11: NUMERICAL PAPER DIAGNOSTIC")
print("=" * 80)

primes = list(primerange(2, 100))

cases = 0
failures = []

for a, b, c, d in combinations(primes, 4):

    values = {
        e1: a + b + c + d,
        e2: a*b + a*c + a*d + b*c + b*d + c*d,
        e3: a*b*c + a*b*d + a*c*d + b*c*d,
        e4: a*b*c*d,
    }

    actual = int(kappa4.subs(values))
    paper = int(paper_kappa4.subs(values))

    cases += 1

    if actual == paper:
        failures.append(
            ((a, b, c, d), actual, paper)
        )

print("TEST CASES =", cases)

if failures:
    print("WARNING: paper formula accidentally agrees in:")
    for item in failures[:10]:
        print(item)
else:
    print("PASS: paper formula fails on every tested 4-body case.")


# =============================================================================
# TEST 12: NUMERICAL CORRECTION CHECK
# =============================================================================

print("\n" + "=" * 80)
print("TEST 12: NUMERICAL CHECK OF COMPACT CORRECTION")
print("=" * 80)

correction_failures = []

for a, b, c, d in combinations(primes[:20], 4):

    E1 = a + b + c + d
    E2 = a*b + a*c + a*d + b*c + b*d + c*d
    E3 = a*b*c + a*b*d + a*c*d + b*c*d
    N = a*b*c*d

    k3 = int(
        kappa3.subs({
            e1: E1,
            e2: E2,
            e3: E3
        })
    )

    k4 = int(
        kappa4.subs({
            e1: E1,
            e2: E2,
            e3: E3,
            e4: N
        })
    )

    compact = (
        N * (-2*E1 + E2 + E3)
        - E3 * (E3 - 1)
    )

    if k4 - k3 != compact:
        correction_failures.append(
            ((a, b, c, d), k4, k3, compact)
        )

if correction_failures:
    print("FAILURES:")
    for item in correction_failures[:10]:
        print(item)
else:
    print("PASS: compact correction holds numerically.")


# =============================================================================
# TEST 13: WHY THE BODY-4 FORMULA CANNOT BE THE BODY-3 FORMULA
# =============================================================================

print("\n" + "=" * 80)
print("TEST 13: STRUCTURAL REASON")
print("=" * 80)

print("""
The 3-body normalized cubic product contains:

    e3^2 + 2*e3 + 1

and therefore:

    X3/M1 = e3^2 - e3 + 1 - kappa3.

For four bodies the new elementary symmetric variable e4 appears
and the actual normalized product contains:

    e3^2
    - e3*e4
    + 2*e3
    + e4^2
    - e4
    + 1

as well as the cross terms involving e1,e2.

Therefore the 4-body transition is not obtained by simply replacing
e3 with e4.

The missing information is the interaction between the old e3 layer
and the new e4 layer.

The exact kappa correction is:

    kappa4 - kappa3
      = n*(-2*e1 + e2 + e3)
        - e3*(e3 - 1)

with n = e4.
""")


# =============================================================================
# TEST 14: FINAL COMPACT RESULT
# =============================================================================

print("\n" + "=" * 80)
print("FINAL STRUCTURAL RESULT")
print("=" * 80)

print("""
3-BODY:

    R3 = M1_3 * kappa3

    kappa3 =
        -e1^2
        + e1*e2
        + e1*e3
        + e1
        - e2^2
        + e2*e3
        + e2
        - 3*e3


4-BODY:

    R4 = M1_4 * kappa4

    kappa4 =
        -e1^2
        + e1*e2
        + e1*e3
        - 2*e1*e4
        + e1
        - e2^2
        + e2*e3
        + e2*e4
        + e2
        - e3^2
        + e3*e4
        - 2*e3


With e4 = n:

    kappa4 =
        -e1^2
        + e1*e2
        + e1*e3
        - 2*e1*n
        + e1
        - e2^2
        + e2*e3
        + e2*n
        + e2
        - e3^2
        + e3*n
        - 2*e3


Exact body transition:

    kappa4 - kappa3
      = n*(-2*e1 + e2 + e3)
        - e3*(e3 - 1)


This is the quantity to compare directly against the paper formula.
""")

print("=" * 80)
print("DONE")
print("=" * 80)

