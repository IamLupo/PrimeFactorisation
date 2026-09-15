import sympy as sp
import itertools


# =============================================================================
# KAPPA STRUCTURE — 3 -> 4 -> 5 EXACT TRANSITION ANALYSIS
#
# Definitions:
#
#   M1_n = sum of all square-free monomials in the n variables
#           = product(1 + x_i)
#
#   X_n  = sum of cubes of all square-free monomials
#           = product(1 + x_i^3)
#
#   P_n  = X_n / M1_n
#
#   kappa_n = 1 - P_n
#
# Everything below is derived directly from these definitions.
# No guessed kappa formulas are used.
# =============================================================================


sp.init_printing()


# =============================================================================
# HELPERS
# =============================================================================

def section(title):
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)


def elementary_symmetric_variables(n):
    """
    Return formal elementary symmetric variables e1,...,en.
    """
    return sp.symbols("e1:" + str(n + 1))


def all_square_free_monomials(vars_):
    """
    Generate every square-free monomial, including 1.
    """
    terms = []

    for k in range(len(vars_) + 1):
        for subset in itertools.combinations(vars_, k):
            terms.append(sp.prod(subset))

    return terms


def build_M1(vars_):
    """
    M1_n = sum of all square-free monomials.
    """
    return sp.expand(sum(all_square_free_monomials(vars_)))


def build_X(vars_):
    """
    X_n = sum of cubes of all square-free monomials.
    """
    return sp.expand(
        sum(term**3 for term in all_square_free_monomials(vars_))
    )


def build_exact_objects(vars_):
    """
    Construct M1, X, exact normalized product P, and exact kappa.
    """
    M1 = build_M1(vars_)
    X = build_X(vars_)

    P = sp.cancel(X / M1)

    # kappa = 1 - normalized product
    kappa = sp.cancel(1 - P)

    return M1, X, P, kappa


def symmetric_reduce(expr, vars_):
    """
    Convert a symmetric polynomial in vars_ into elementary symmetric
    variables using SymPy's exact symmetrization.
    """
    result = sp.symmetrize(
        sp.expand(expr),
        vars_,
        formal=True
    )

    symmetric_expr = sp.expand(result[0])
    remainder = sp.expand(result[1])
    substitutions = result[2]

    return symmetric_expr, remainder, substitutions


def make_elementary_substitution(substitutions):
    """
    Convert SymPy's formal symbols s1,s2,... into e1,e2,...
    """
    formal_symbols = [item[0] for item in substitutions]

    e_symbols = sp.symbols(
        "e1:" + str(len(formal_symbols) + 1)
    )

    return {
        formal_symbols[i]: e_symbols[i]
        for i in range(len(formal_symbols))
    }


# =============================================================================
# VARIABLES
# =============================================================================

p, q, r, s, t = sp.symbols("p q r s t")

vars3 = (p, q, r)
vars4 = (p, q, r, s)
vars5 = (p, q, r, s, t)


# =============================================================================
# TEST 1 — EXACT 3-BODY OBJECTS
# =============================================================================

section("TEST 1: EXACT 3-BODY DEFINITIONS")

M1_3, X3, P3, kappa3_direct = build_exact_objects(vars3)

print("M1_3 =")
print(M1_3)

print()
print("X3 =")
print(X3)

print()
print("P3 = cancel(X3/M1_3)")
print(P3)

print()
print("kappa3 = 1 - P3")
print(kappa3_direct)

assert sp.expand(X3 - M1_3 * P3) == 0
assert sp.expand(kappa3_direct - (1 - P3)) == 0

print()
print("PASS: exact 3-body normalization confirmed.")


# =============================================================================
# TEST 2 — EXACT 4-BODY OBJECTS
# =============================================================================

section("TEST 2: EXACT 4-BODY DEFINITIONS")

M1_4, X4, P4, kappa4_direct = build_exact_objects(vars4)

print("M1_4 =")
print(M1_4)

print()
print("X4 =")
print(X4)

print()
print("P4 = cancel(X4/M1_4)")
print(P4)

print()
print("kappa4 = 1 - P4")
print(kappa4_direct)

assert sp.expand(X4 - M1_4 * P4) == 0
assert sp.expand(kappa4_direct - (1 - P4)) == 0

print()
print("PASS: exact 4-body normalization confirmed.")


# =============================================================================
# TEST 3 — EXACT 5-BODY OBJECTS
# =============================================================================

section("TEST 3: EXACT 5-BODY DEFINITIONS")

M1_5, X5, P5, kappa5_direct = build_exact_objects(vars5)

print("M1_5 =")
print(M1_5)

print()
print("X5 =")
print(X5)

print()
print("P5 = cancel(X5/M1_5)")
print(P5)

print()
print("kappa5 = 1 - P5")
print(kappa5_direct)

assert sp.expand(X5 - M1_5 * P5) == 0
assert sp.expand(kappa5_direct - (1 - P5)) == 0

print()
print("PASS: exact 5-body normalization confirmed.")


# =============================================================================
# TEST 4 — ELEMENTARY SYMMETRIC REDUCTION
# =============================================================================

section("TEST 4: EXACT ELEMENTARY-SYMMETRIC REDUCTION")


# -----------------------------------------------------------------------------
# 3-body
# -----------------------------------------------------------------------------

kappa3_sym_raw, rem3, subs3 = symmetric_reduce(
    kappa3_direct,
    vars3
)

sub_map3 = make_elementary_substitution(subs3)

kappa3 = sp.expand(
    kappa3_sym_raw.subs(sub_map3)
)

print("kappa3 symmetric =")
print(kappa3)

print()
print("3-body symmetrization remainder =")
print(rem3)

assert rem3 == 0


# -----------------------------------------------------------------------------
# 4-body
# -----------------------------------------------------------------------------

kappa4_sym_raw, rem4, subs4 = symmetric_reduce(
    kappa4_direct,
    vars4
)

sub_map4 = make_elementary_substitution(subs4)

kappa4 = sp.expand(
    kappa4_sym_raw.subs(sub_map4)
)

print()
print("kappa4 symmetric =")
print(kappa4)

print()
print("4-body symmetrization remainder =")
print(rem4)

assert rem4 == 0


# -----------------------------------------------------------------------------
# 5-body
# -----------------------------------------------------------------------------

kappa5_sym_raw, rem5, subs5 = symmetric_reduce(
    kappa5_direct,
    vars5
)

sub_map5 = make_elementary_substitution(subs5)

kappa5 = sp.expand(
    kappa5_sym_raw.subs(sub_map5)
)

print()
print("kappa5 symmetric =")
print(kappa5)

print()
print("5-body symmetrization remainder =")
print(rem5)

assert rem5 == 0


# =============================================================================
# TEST 5 — DIRECT VS SYMMETRIC KAPPA
# =============================================================================

section("TEST 5: DIRECT KAPPA VS SYMMETRIC KAPPA")


# 3-body substitution back
e1, e2, e3 = sp.symbols("e1 e2 e3")

elementary3 = {
    e1: p + q + r,
    e2: p*q + p*r + q*r,
    e3: p*q*r,
}

check3 = sp.expand(
    kappa3.subs(elementary3) - kappa3_direct
)

print("3-body check:")
print(check3)

assert check3 == 0

# 4-body substitution back
e1, e2, e3, e4 = sp.symbols("e1 e2 e3 e4")

elementary4 = {
    e1: p + q + r + s,
    e2: (
        p*q + p*r + p*s +
        q*r + q*s +
        r*s
    ),
    e3: (
        p*q*r +
        p*q*s +
        p*r*s +
        q*r*s
    ),
    e4: p*q*r*s,
}

check4 = sp.expand(
    kappa4.subs(elementary4) - kappa4_direct
)

print()
print("4-body check:")
print(check4)

assert check4 == 0

# 5-body substitution back
e1, e2, e3, e4, e5 = sp.symbols(
    "e1 e2 e3 e4 e5"
)

elementary5 = {
    e1: p + q + r + s + t,

    e2: (
        p*q + p*r + p*s + p*t +
        q*r + q*s + q*t +
        r*s + r*t +
        s*t
    ),

    e3: (
        p*q*r +
        p*q*s +
        p*q*t +
        p*r*s +
        p*r*t +
        p*s*t +
        q*r*s +
        q*r*t +
        q*s*t +
        r*s*t
    ),

    e4: (
        p*q*r*s +
        p*q*r*t +
        p*q*s*t +
        p*r*s*t +
        q*r*s*t
    ),

    e5: p*q*r*s*t,
}

check5 = sp.expand(
    kappa5.subs(elementary5) - kappa5_direct
)

print()
print("5-body check:")
print(check5)

assert check5 == 0

print()
print("PASS: all direct kappa formulas exactly match")
print("their elementary-symmetric representations.")


# =============================================================================
# TEST 6 — CLEAN EXACT KAPPA FORMULAS
# =============================================================================

section("TEST 6: CLEAN EXACT KAPPA FORMULAS")

print("kappa3(e1,e2,e3) =")
print(kappa3)

print()
print("kappa4(e1,e2,e3,e4) =")
print(kappa4)

print()
print("kappa5(e1,e2,e3,e4,e5) =")
print(kappa5)


# =============================================================================
# TEST 7 — EXACT 3 -> 4 CORRECTION
# =============================================================================

section("TEST 7: EXACT 3 -> 4 CORRECTION")

n = sp.symbols("n")

kappa4_e4_n = sp.expand(
    kappa4.subs(e4, n)
)

correction34 = sp.expand(
    kappa4_e4_n - kappa3
)

print("kappa4(e4=n) - kappa3 =")
print(correction34)

print()
print("Factored:")
print(sp.factor(correction34))

print()
print("Collected:")
print(sp.collect(correction34, n))

# Expected exact correction:
expected34 = sp.expand(
    -e3**2
    + e3
    + n * (-2*e1 + e2 + e3 + 1)
    - n**2
)

print()
print("Expected exact 3 -> 4 correction =")
print(expected34)

check34 = sp.expand(
    correction34 - expected34
)

print()
print("Check:")
print(check34)

assert check34 == 0

print()
print("PASS: exact 3 -> 4 correction confirmed.")


# =============================================================================
# TEST 8 — EXACT 4 -> 5 CORRECTION
# =============================================================================

section("TEST 8: EXACT 4 -> 5 CORRECTION")

kappa5_e5_n = sp.expand(
    kappa5.subs(e5, n)
)

correction45 = sp.expand(
    kappa5_e5_n - kappa4
)

print("kappa5(e5=n) - kappa4 =")
print(correction45)

print()
print("Factored:")
print(sp.factor(correction45))

print()
print("Collected:")
print(sp.collect(correction45, n))

# Expected exact correction:
expected45 = sp.expand(
    n * (e1 - 2*e2 + e3 + e4 + 1)
    - n**2
)

print()
print("Expected exact 4 -> 5 correction =")
print(expected45)

check45 = sp.expand(
    correction45 - expected45
)

print()
print("Check:")
print(check45)

assert check45 == 0

print()
print("PASS: exact 4 -> 5 correction confirmed.")


# =============================================================================
# TEST 9 — DEGREE OF 3 -> 4 CORRECTION
# =============================================================================

section("TEST 9: DEGREE IN NEW VARIABLE — 3 -> 4")

degree34 = sp.Poly(
    correction34,
    n
).degree()

print("degree of 3 -> 4 correction in n =", degree34)

# Important:
# The correction is genuinely quadratic.
#
# Therefore the correct assertion is degree == 2,
# not degree == 1.

assert degree34 == 2

print()
print("PASS: 3 -> 4 correction is exactly quadratic in n.")


# =============================================================================
# TEST 10 — DEGREE OF 4 -> 5 CORRECTION
# =============================================================================

section("TEST 10: DEGREE IN NEW VARIABLE — 4 -> 5")

degree45 = sp.Poly(
    correction45,
    n
).degree()

print("degree of 4 -> 5 correction in n =", degree45)

assert degree45 == 2

print()
print("PASS: 4 -> 5 correction is exactly quadratic in n.")


# =============================================================================
# TEST 11 — COEFFICIENT STRUCTURE OF 3 -> 4
# =============================================================================

section("TEST 11: COEFFICIENT STRUCTURE — 3 -> 4")

coeff34_n2 = sp.expand(
    sp.Poly(correction34, n).coeff_monomial(n**2)
)

coeff34_n1 = sp.expand(
    sp.Poly(correction34, n).coeff_monomial(n)
)

coeff34_n0 = sp.expand(
    sp.Poly(correction34, n).coeff_monomial(1)
)

print("n^2 coefficient:")
print(coeff34_n2)

print()
print("n^1 coefficient:")
print(coeff34_n1)

print()
print("n^0 coefficient:")
print(coeff34_n0)

assert coeff34_n2 == -1
assert coeff34_n1 == -2*e1 + e2 + e3 + 1
assert coeff34_n0 == -e3**2 + e3

print()
print("PASS: exact coefficient structure confirmed.")


# =============================================================================
# TEST 12 — COEFFICIENT STRUCTURE OF 4 -> 5
# =============================================================================

section("TEST 12: COEFFICIENT STRUCTURE — 4 -> 5")

coeff45_n2 = sp.expand(
    sp.Poly(correction45, n).coeff_monomial(n**2)
)

coeff45_n1 = sp.expand(
    sp.Poly(correction45, n).coeff_monomial(n)
)

coeff45_n0 = sp.expand(
    sp.Poly(correction45, n).coeff_monomial(1)
)

print("n^2 coefficient:")
print(coeff45_n2)

print()
print("n^1 coefficient:")
print(coeff45_n1)

print()
print("n^0 coefficient:")
print(coeff45_n0)

assert coeff45_n2 == -1
assert coeff45_n1 == e1 - 2*e2 + e3 + e4 + 1
assert coeff45_n0 == 0

print()
print("PASS: exact coefficient structure confirmed.")


# =============================================================================
# TEST 13 — FACTORED 4 -> 5 FORM
# =============================================================================

section("TEST 13: CLEAN FACTORED 4 -> 5 FORM")

factored45 = sp.factor(correction45)

expected_factored45 = (
    n * (
        e1
        - 2*e2
        + e3
        + e4
        - n
        + 1
    )
)

print("Exact factored correction:")
print(factored45)

print()
print("Expected:")
print(expected_factored45)

assert sp.expand(
    factored45 - expected_factored45
) == 0

print()
print("PASS: 4 -> 5 correction has exact factorization.")


# =============================================================================
# TEST 14 — FACTORED 3 -> 4 FORM
# =============================================================================

section("TEST 14: CLEAN FACTORED 3 -> 4 FORM")

factored34 = sp.factor(correction34)

expected_factored34 = (
    -e3**2
    + e3
    + n * (
        -2*e1
        + e2
        + e3
        + 1
        - n
    )
)

print("Exact factored correction:")
print(factored34)

print()
print("Expected:")
print(expected_factored34)

assert sp.expand(
    factored34 - expected_factored34
) == 0

print()
print("PASS: 3 -> 4 correction has exact quadratic structure.")


# =============================================================================
# TEST 15 — CROSS-MULTIPLICATION IDENTITIES
# =============================================================================

section("TEST 15: CROSS-MULTIPLICATION IDENTITIES")

cross3 = sp.expand(
    X3 - M1_3 * P3
)

cross4 = sp.expand(
    X4 - M1_4 * P4
)

cross5 = sp.expand(
    X5 - M1_5 * P5
)

print("X3 - M1_3*P3 =")
print(cross3)

print()
print("X4 - M1_4*P4 =")
print(cross4)

print()
print("X5 - M1_5*P5 =")
print(cross5)

assert cross3 == 0
assert cross4 == 0
assert cross5 == 0

print()
print("PASS: all normalized products satisfy exact")
print("cross-multiplication identities.")


# =============================================================================
# TEST 16 — POLYNOMIAL NORMALIZATION
# =============================================================================

section("TEST 16: EXACT NORMALIZED PRODUCTS ARE POLYNOMIALS")

for label, P, vars_ in [
    ("P3", P3, vars3),
    ("P4", P4, vars4),
    ("P5", P5, vars5),
]:

    numerator, denominator = sp.fraction(
        sp.cancel(P)
    )

    denominator = sp.expand(denominator)

    print(label, "denominator after cancel =", denominator)

    assert denominator == 1

    print("PASS:", label, "is polynomial.")

    print()


# =============================================================================
# TEST 17 — VERIFY THE 3-BODY FORMULA EXPLICITLY
# =============================================================================

section("TEST 17: CLEAN EXACT 3-BODY KAPPA")

expected_kappa3 = sp.expand(
    -e1**2
    + e1*e2
    + e1*e3
    + e1
    - e2**2
    + e2*e3
    + e2
    - e3**2
    - 2*e3
)

print("Derived kappa3:")
print(kappa3)

print()
print("Expected exact kappa3:")
print(expected_kappa3)

check_kappa3 = sp.expand(
    kappa3 - expected_kappa3
)

print()
print("Check:")
print(check_kappa3)

assert check_kappa3 == 0

print()
print("PASS: exact 3-body kappa formula confirmed.")


# =============================================================================
# TEST 18 — VERIFY THE 4-BODY FORMULA EXPLICITLY
# =============================================================================

section("TEST 18: CLEAN EXACT 4-BODY KAPPA")

expected_kappa4 = sp.expand(
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

print("Derived kappa4:")
print(kappa4)

print()
print("Expected exact kappa4:")
print(expected_kappa4)

check_kappa4 = sp.expand(
    kappa4 - expected_kappa4
)

print()
print("Check:")
print(check_kappa4)

assert check_kappa4 == 0

print()
print("PASS: exact 4-body kappa formula confirmed.")


# =============================================================================
# TEST 19 — VERIFY THE 5-BODY FORMULA EXPLICITLY
# =============================================================================

section("TEST 19: CLEAN EXACT 5-BODY KAPPA")

expected_kappa5 = sp.expand(
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

print("Derived kappa5:")
print(kappa5)

print()
print("Expected exact kappa5:")
print(expected_kappa5)

check_kappa5 = sp.expand(
    kappa5 - expected_kappa5
)

print()
print("Check:")
print(check_kappa5)

assert check_kappa5 == 0

print()
print("PASS: exact 5-body kappa formula confirmed.")


# =============================================================================
# FINAL STRUCTURAL SUMMARY
# =============================================================================

section("FINAL EXACT STRUCTURAL RESULT")

print("EXACT kappa3:")
print(kappa3)

print()
print("EXACT kappa4:")
print(kappa4)

print()
print("EXACT kappa5:")
print(kappa5)

print()
print("-" * 80)

print("EXACT 3 -> 4 CORRECTION:")
print(
    sp.factor(correction34)
)

print()
print("EXACT 4 -> 5 CORRECTION:")
print(
    sp.factor(correction45)
)

print()
print("-" * 80)

print("DEGREE(3 -> 4) =", degree34)
print("DEGREE(4 -> 5) =", degree45)

print()
print("The exact 3 -> 4 correction is:")
print(
    "-e3^2 + e3 + n*(-2*e1 + e2 + e3 + 1) - n^2"
)

print()
print("The exact 4 -> 5 correction is:")
print(
    "n*(e1 - 2*e2 + e3 + e4 + 1 - n)"
)

print()
print("=" * 80)
print("ALL EXACT TESTS PASSED")
print("=" * 80)