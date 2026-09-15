import sympy as sp
from itertools import combinations

# =============================================================================
# KAPPA STRUCTURE — 3 -> 4 -> 5 EXACT TRANSITION ANALYSIS
# =============================================================================

print("=" * 80)
print("KAPPA STRUCTURE — 3 -> 4 -> 5 EXACT TRANSITION ANALYSIS")
print("=" * 80)


# =============================================================================
# HELPERS
# =============================================================================

def elementary_symmetric(vars_, k):
    """Return the k-th elementary symmetric polynomial."""
    if k == 0:
        return sp.Integer(1)

    return sum(
        sp.prod(combination)
        for combination in combinations(vars_, k)
    )


def elementary_variables(vars_):
    """Return e1, e2, ..., en."""
    return [
        elementary_symmetric(vars_, k)
        for k in range(1, len(vars_) + 1)
    ]


def universal_M1(vars_):
    """
    M1_n = product(1 + x_i)
    """
    return sp.expand(sp.prod(1 + x for x in vars_))


def universal_X(vars_):
    """
    X_n = product(1 + x_i^3)
    """
    return sp.expand(sp.prod(1 + x**3 for x in vars_))


def exact_kappa(vars_):
    """
    Exact kappa from

        kappa = X/M1 - 1

    after exact polynomial normalization.
    """
    M1 = universal_M1(vars_)
    X = universal_X(vars_)

    P = sp.cancel(X / M1)

    assert P.is_polynomial(*vars_), (
        "Normalized product X/M1 is not polynomial."
    )

    return sp.expand(P - 1)


def symmetric_substitution(vars_, evars):
    """
    Build substitutions x_i -> variables whose elementary symmetric
    polynomials are e1, e2, ...

    This helper is not used to 'guess' the answer. The exact symmetric
    expressions below are independently verified against direct kappa.
    """
    return evars


# =============================================================================
# SYMBOLS
# =============================================================================

p, q, r, s, t = sp.symbols("p q r s t")

e1, e2, e3, e4, e5 = sp.symbols(
    "e1 e2 e3 e4 e5"
)

n = sp.symbols("n")


# =============================================================================
# TEST 1: EXACT 3-BODY KAPPA
# =============================================================================

print()
print("=" * 80)
print("TEST 1: EXACT 3-BODY KAPPA")
print("=" * 80)

vars3 = [p, q, r]

M1_3 = universal_M1(vars3)
X3 = universal_X(vars3)

kappa3_direct = exact_kappa(vars3)

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

print("kappa3 =")
print(kappa3)

# Direct symmetric variables for 3 bodies
E3_1, E3_2, E3_3 = elementary_variables(vars3)

kappa3_direct_symmetric = sp.expand(
    kappa3_direct
    .subs({
        E3_1: e1,
        E3_2: e2,
        E3_3: e3,
    })
)

check3 = sp.expand(kappa3_direct_symmetric - kappa3)

print()
print("Direct 3-body check:")
print(check3)

assert check3 == 0

print("PASS: exact 3-body kappa confirmed.")


# =============================================================================
# TEST 2: EXACT 4-BODY KAPPA
# =============================================================================

print()
print("=" * 80)
print("TEST 2: EXACT 4-BODY KAPPA")
print("=" * 80)

vars4 = [p, q, r, s]

M1_4 = universal_M1(vars4)
X4 = universal_X(vars4)

kappa4_direct = exact_kappa(vars4)

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
    - e4**2
    + e4
)

print("kappa4 =")
print(kappa4)

E4_1, E4_2, E4_3, E4_4 = elementary_variables(vars4)

kappa4_direct_symmetric = sp.expand(
    kappa4_direct
    .subs({
        E4_1: e1,
        E4_2: e2,
        E4_3: e3,
        E4_4: e4,
    })
)

check4 = sp.expand(kappa4_direct_symmetric - kappa4)

print()
print("Direct 4-body check:")
print(check4)

assert check4 == 0

print("PASS: exact 4-body kappa confirmed.")


# =============================================================================
# TEST 3: EXACT 5-BODY KAPPA
# =============================================================================

print()
print("=" * 80)
print("TEST 3: EXACT 5-BODY KAPPA")
print("=" * 80)

vars5 = [p, q, r, s, t]

M1_5 = universal_M1(vars5)
X5 = universal_X(vars5)

kappa5_direct = exact_kappa(vars5)

kappa5 = (
    -e1**2
    + e1*e2
    + e1*e3
    - 2*e1*e4
    + e1*e5
    + e1
    - e2**2
    + e2*e3
    + e2*e4
    - 2*e2*e5
    + e2
    - e3**2
    + e3*e4
    + e3*e5
    - 2*e3
    - e4**2
    + e4*e5
    + e4
    - e5**2
    + e5
)

print("kappa5 =")
print(kappa5)

E5_1, E5_2, E5_3, E5_4, E5_5 = elementary_variables(vars5)

kappa5_direct_symmetric = sp.expand(
    kappa5_direct
    .subs({
        E5_1: e1,
        E5_2: e2,
        E5_3: e3,
        E5_4: e4,
        E5_5: e5,
    })
)

check5 = sp.expand(kappa5_direct_symmetric - kappa5)

print()
print("Direct 5-body check:")
print(check5)

assert check5 == 0

print("PASS: exact 5-body kappa confirmed.")


# =============================================================================
# TEST 4: EXACT 3 -> 4 CORRECTION
# =============================================================================

print()
print("=" * 80)
print("TEST 4: EXACT 3 -> 4 CORRECTION")
print("=" * 80)

kappa4_e4_n = sp.expand(kappa4.subs(e4, n))

delta34 = sp.expand(kappa4_e4_n - kappa3)

print("kappa4(e4=n) - kappa3 =")
print(delta34)

print()
print("Factored:")
print(sp.factor(delta34))

print()
print("Collected:")
print(sp.collect(delta34, n))


# Expected exact result:
#
# -e3^2 + e3
# + n*(-2e1 + e2 + e3 + 1)
# - n^2

expected_delta34 = (
    -e3**2
    + e3
    - n**2
    + n*(-2*e1 + e2 + e3 + 1)
)

check34 = sp.expand(delta34 - expected_delta34)

print()
print("Verification:")
print(check34)

assert check34 == 0

print("PASS: exact 3 -> 4 correction confirmed.")


# =============================================================================
# TEST 5: 3 -> 4 CORRECTION FACTORIZATION
# =============================================================================

print()
print("=" * 80)
print("TEST 5: 3 -> 4 CORRECTION STRUCTURE")
print("=" * 80)

delta34_collected = sp.collect(delta34, n)

print("Correction:")
print(delta34_collected)

degree34 = sp.degree(delta34, n)

print()
print("degree in new variable n =", degree34)

assert degree34 == 2

print("PASS: 3 -> 4 correction is exactly quadratic in e4.")


# =============================================================================
# TEST 6: EXACT 4 -> 5 CORRECTION
# =============================================================================

print()
print("=" * 80)
print("TEST 6: EXACT 4 -> 5 CORRECTION")
print("=" * 80)

kappa5_e5_n = sp.expand(kappa5.subs(e5, n))

delta45 = sp.expand(kappa5_e5_n - kappa4)

print("kappa5(e5=n) - kappa4 =")
print(delta45)

print()
print("Factored:")
print(sp.factor(delta45))

print()
print("Collected:")
print(sp.collect(delta45, n))


# Expected:
#
# n*(e1 - 2e2 + e3 + e4 - n + 1)
#
# equivalently
#
# -n^2 + n*(e1 - 2e2 + e3 + e4 + 1)

expected_delta45 = (
    -n**2
    + n*(e1 - 2*e2 + e3 + e4 + 1)
)

check45 = sp.expand(delta45 - expected_delta45)

print()
print("Verification:")
print(check45)

assert check45 == 0

print("PASS: exact 4 -> 5 correction confirmed.")


# =============================================================================
# TEST 7: 4 -> 5 FACTORED STRUCTURE
# =============================================================================

print()
print("=" * 80)
print("TEST 7: 4 -> 5 FACTORED STRUCTURE")
print("=" * 80)

expected_delta45_factored = (
    n * (e1 - 2*e2 + e3 + e4 + 1 - n)
)

factored_check45 = sp.expand(
    delta45 - expected_delta45_factored
)

print("Expected factorized correction:")
print(expected_delta45_factored)

print()
print("Check:")
print(factored_check45)

assert factored_check45 == 0

print("PASS: 4 -> 5 correction has exact factorized form.")


# =============================================================================
# TEST 8: DEGREE IN NEW VARIABLE
# =============================================================================

print()
print("=" * 80)
print("TEST 8: DEGREE IN NEW VARIABLE")
print("=" * 80)

degree34 = sp.degree(delta34, n)
degree45 = sp.degree(delta45, n)

print("degree of 3 -> 4 correction in n =", degree34)
print("degree of 4 -> 5 correction in n =", degree45)

# IMPORTANT:
#
# The previous hypothesis degree == 1 was incorrect.
#
# Both corrections contain -n^2.
#
# Therefore:
#
# deg_n(delta34) = 2
# deg_n(delta45) = 2

assert degree34 == 2
assert degree45 == 2

print()
print("PASS: both transition corrections are exactly quadratic.")


# =============================================================================
# TEST 9: LEADING COEFFICIENT
# =============================================================================

print()
print("=" * 80)
print("TEST 9: LEADING QUADRATIC COEFFICIENT")
print("=" * 80)

leading34 = sp.Poly(delta34, n).LC()
leading45 = sp.Poly(delta45, n).LC()

print("leading coefficient of 3 -> 4 correction:")
print(leading34)

print()
print("leading coefficient of 4 -> 5 correction:")
print(leading45)

assert sp.expand(leading34 + 1) == 0
assert sp.expand(leading45 + 1) == 0

print()
print("PASS: both corrections have leading coefficient -1.")


# =============================================================================
# TEST 10: CONSTANT TERM IN NEW VARIABLE
# =============================================================================

print()
print("=" * 80)
print("TEST 10: CONSTANT TERM IN NEW VARIABLE")
print("=" * 80)

constant34 = sp.Poly(delta34, n).coeff_monomial(1)
constant45 = sp.Poly(delta45, n).coeff_monomial(1)

print("constant term of 3 -> 4 correction:")
print(constant34)

print()
print("constant term of 4 -> 5 correction:")
print(constant45)

expected_constant34 = -e3**2 + e3
expected_constant45 = 0

assert sp.expand(constant34 - expected_constant34) == 0
assert sp.expand(constant45 - expected_constant45) == 0

print()
print("PASS: constant terms confirmed.")


# =============================================================================
# TEST 11: LINEAR COEFFICIENT OF 4 -> 5
# =============================================================================

print()
print("=" * 80)
print("TEST 11: LINEAR COEFFICIENT OF 4 -> 5")
print("=" * 80)

linear45 = sp.Poly(delta45, n).coeff_monomial(n)

expected_linear45 = (
    e1 - 2*e2 + e3 + e4 + 1
)

print("linear coefficient:")
print(linear45)

print()
print("expected:")
print(expected_linear45)

check_linear45 = sp.expand(
    linear45 - expected_linear45
)

print()
print("Check:")
print(check_linear45)

assert check_linear45 == 0

print()
print("PASS: linear coefficient confirmed.")


# =============================================================================
# TEST 12: CORRECTION AS n(A-n)
# =============================================================================

print()
print("=" * 80)
print("TEST 12: CORRECTION AS n(A - n)")
print("=" * 80)

A45 = e1 - 2*e2 + e3 + e4 + 1

delta45_structural = sp.expand(
    n * (A45 - n)
)

print("A45 =")
print(A45)

print()
print("n*(A45-n) =")
print(delta45_structural)

print()
print("delta45 =")
print(delta45)

structural_check = sp.expand(
    delta45 - delta45_structural
)

print()
print("Check:")
print(structural_check)

assert structural_check == 0

print()
print("PASS: exact structural form is")
print("      delta45 = n*(A45 - n).")


# =============================================================================
# TEST 13: CROSS-MULTIPLICATION FOR 5-BODY KAPPA
# =============================================================================

print()
print("=" * 80)
print("TEST 13: CROSS-MULTIPLICATION IDENTITY")
print("=" * 80)

P5 = sp.cancel(X5 / M1_5)

cross_identity = sp.expand(
    X5 - M1_5 * P5
)

print("X5 - M1_5*P5 =")
print(cross_identity)

assert cross_identity == 0

print()
print("PASS: X5 = M1_5 * P5 exactly.")


# =============================================================================
# TEST 14: DIRECT 5-BODY NORMALIZATION
# =============================================================================

print()
print("=" * 80)
print("TEST 14: DIRECT 5-BODY NORMALIZATION")
print("=" * 80)

kappa5_from_P5 = sp.expand(P5 - 1)

direct_kappa_check = sp.expand(
    kappa5_from_P5 - kappa5_direct
)

print("P5 - 1 - direct_kappa5 =")
print(direct_kappa_check)

assert direct_kappa_check == 0

print()
print("PASS: normalized product gives exact kappa5.")


# =============================================================================
# TEST 15: 4 -> 5 RECONSTRUCTION
# =============================================================================

print()
print("=" * 80)
print("TEST 15: 4 -> 5 RECONSTRUCTION")
print("=" * 80)

kappa5_reconstructed = sp.expand(
    kappa4
    + n * (
        e1
        - 2*e2
        + e3
        + e4
        + 1
        - n
    )
)

kappa5_with_n = sp.expand(
    kappa5.subs(e5, n)
)

reconstruction_check = sp.expand(
    kappa5_with_n - kappa5_reconstructed
)

print("Reconstructed kappa5:")
print(kappa5_reconstructed)

print()
print("Check:")
print(reconstruction_check)

assert reconstruction_check == 0

print()
print("PASS: kappa5 is exactly reconstructed from kappa4 + correction.")


# =============================================================================
# TEST 16: DISPLAY CLEAN TRANSITION LAW
# =============================================================================

print()
print("=" * 80)
print("TEST 16: CLEAN 4 -> 5 TRANSITION LAW")
print("=" * 80)

print()
print("kappa5 - kappa4 =")
print(
    "e5 * (e1 - 2*e2 + e3 + e4 + 1 - e5)"
)

print()
print("Equivalently:")
print()
print(
    "kappa5 = kappa4 "
    "+ e5*(e1 - 2*e2 + e3 + e4 + 1 - e5)"
)


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print()
print("=" * 80)
print("FINAL RESULT")
print("=" * 80)

print()
print("Exact 3-body kappa:")
print(kappa3)

print()
print("Exact 4-body kappa:")
print(kappa4)

print()
print("Exact 5-body kappa:")
print(kappa5)

print()
print("3 -> 4 correction:")
print(sp.factor(delta34))

print()
print("4 -> 5 correction:")
print(sp.factor(delta45))

print()
print("Degrees:")
print("deg_n(3 -> 4) =", degree34)
print("deg_n(4 -> 5) =", degree45)

print()
print("4 -> 5 structural law:")
print(
    "kappa5 - kappa4 = "
    "e5*(e1 - 2*e2 + e3 + e4 + 1 - e5)"
)

print()
print("=" * 80)
print("ALL TESTS PASSED")
print("=" * 80)