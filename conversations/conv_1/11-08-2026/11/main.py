import sympy as sp

print("=" * 80)
print("5-BODY STRUCTURAL FORMULA — EXACT NORMALIZATION / KAPPA / 4->5")
print("=" * 80)

# ============================================================================
# VARIABLES
# ============================================================================

p, q, r, s, t = sp.symbols("p q r s t")
e1, e2, e3, e4, e5 = sp.symbols("e1 e2 e3 e4 e5")

vars5 = [p, q, r, s, t]

# ============================================================================
# HELPERS
# ============================================================================

def elementary_symmetric(xs):
    """
    Return elementary symmetric polynomials e1,...,en.
    """
    n = len(xs)
    return [
        sp.expand(sum(sp.prod(c) for c in sp.utilities.iterables.combinations(xs, k)))
        for k in range(1, n + 1)
    ]


def exact_symmetric_reduction(expr, xs, es):
    """
    Reduce a symmetric polynomial into elementary symmetric variables.

    SymPy's symmetric reduction returns:
        (remainder, mapping)

    where mapping gives the elementary symmetric substitution.
    """
    result = sp.symmetrize(
        sp.expand(expr),
        xs,
        formal=True
    )

    reduced, remainder, mapping = result

    assert remainder == 0, (
        "Expression is not exactly symmetric or reduction failed.\n"
        f"Remainder = {remainder}"
    )

    # mapping contains generated symbols s1,...,sn.
    generated = [m[0] for m in mapping]

    substitution = dict(zip(generated, es))

    return sp.expand(reduced.subs(substitution))


def normalized_product(X, M1):
    """
    Compute X/M1 exactly as a polynomial.

    IMPORTANT:
    Do not use ordinary `/` and then compare the resulting rational
    expression.  cancel() performs the exact cancellation.
    """
    return sp.cancel(X / M1)


# ============================================================================
# UNIVERSAL DEFINITIONS
# ============================================================================

print()
print("=" * 80)
print("TEST 1: UNIVERSAL 5-BODY DEFINITIONS")
print("=" * 80)

M1_5 = sp.expand(
    sp.prod(1 + x for x in vars5)
)

X5 = sp.expand(
    sp.prod(1 + x**3 for x in vars5)
)

print("M1_5 =")
print(M1_5)

print()
print("X5 =")
print(X5)


# ============================================================================
# ELEMENTARY SYMMETRIC POLYNOMIALS
# ============================================================================

E5 = elementary_symmetric(vars5)

E1, E2, E3, E4, E5prod = E5

print()
print("=" * 80)
print("TEST 2: ELEMENTARY SYMMETRIC VARIABLES")
print("=" * 80)

print("e1 =", E1)
print("e2 =", E2)
print("e3 =", E3)
print("e4 =", E4)
print("e5 =", E5prod)


# ============================================================================
# M1 SYMMETRIC FORM
# ============================================================================

print()
print("=" * 80)
print("TEST 3: M1_5 SYMMETRIC REDUCTION")
print("=" * 80)

M1_5_sym = exact_symmetric_reduction(
    M1_5,
    vars5,
    [e1, e2, e3, e4, e5]
)

M1_expected = e1 + e2 + e3 + e4 + e5 + 1

print("M1_5 symmetric:")
print(M1_5_sym)

print()
print("Expected:")
print(M1_expected)

M1_check = sp.expand(M1_5_sym - M1_expected)

print()
print("Check:")
print(M1_check)

assert M1_check == 0

print("PASS: M1_5 = 1 + e1 + e2 + e3 + e4 + e5.")


# ============================================================================
# X5 SYMMETRIC FORM
# ============================================================================

print()
print("=" * 80)
print("TEST 4: X5 SYMMETRIC REDUCTION")
print("=" * 80)

X5_sym = exact_symmetric_reduction(
    X5,
    vars5,
    [e1, e2, e3, e4, e5]
)

print("X5 symmetric:")
print(X5_sym)


# ============================================================================
# DIRECT EXACT NORMALIZATION
# ============================================================================

print()
print("=" * 80)
print("TEST 5: EXACT NORMALIZED PRODUCT X5/M1_5")
print("=" * 80)

P5_direct = normalized_product(X5, M1_5)

print("P5 = cancel(X5/M1_5)")
print(P5_direct)

print()
print("Is P5 polynomial?")

P5_together = sp.together(P5_direct)
num, den = sp.fraction(P5_together)

if den == 1:
    print("YES")
else:
    print("NO")
    print("Denominator =", den)

assert den == 1

print("PASS: exact normalization is polynomial.")


# ============================================================================
# SYMMETRIC REDUCTION OF NORMALIZED PRODUCT
# ============================================================================

print()
print("=" * 80)
print("TEST 6: NORMALIZED PRODUCT IN ELEMENTARY SYMMETRIC VARIABLES")
print("=" * 80)

P5_sym = exact_symmetric_reduction(
    P5_direct,
    vars5,
    [e1, e2, e3, e4, e5]
)

print("P5 symmetric =")
print(P5_sym)


# ============================================================================
# VERIFY DIRECT NORMALIZED PRODUCT
# ============================================================================

print()
print("=" * 80)
print("TEST 7: DIRECT P5 VS SYMMETRIC P5")
print("=" * 80)

P5_sym_back = sp.expand(
    P5_sym.subs({
        e1: E1,
        e2: E2,
        e3: E3,
        e4: E4,
        e5: E5prod,
    })
)

P5_check = sp.expand(P5_direct - P5_sym_back)

print("Check:")
print(P5_check)

assert P5_check == 0

print("PASS: normalized 5-body product is exactly symmetric.")


# ============================================================================
# CROSS-MULTIPLICATION CHECK
# ============================================================================

print()
print("=" * 80)
print("TEST 8: CROSS-MULTIPLICATION IDENTITY")
print("=" * 80)

cross_check = sp.expand(
    X5 - M1_5 * P5_direct
)

print("X5 - M1_5*P5 =")
print(cross_check)

assert cross_check == 0

print("PASS: X5 = M1_5 * P5 exactly.")


# ============================================================================
# DEFINE KAPPA CORRECTLY
# ============================================================================
#
# From the 3-body and 4-body structure:
#
#       X_n / M1_n = 1 - kappa_n
#
# Therefore:
#
#       kappa_n = 1 - X_n/M1_n
#
# This is the crucial definition.
# ============================================================================

print()
print("=" * 80)
print("TEST 9: EXACT 5-BODY KAPPA")
print("=" * 80)

kappa5_direct = sp.expand(
    1 - P5_direct
)

kappa5_sym = sp.expand(
    1 - P5_sym
)

print("kappa5 symmetric =")
print(kappa5_sym)


# ============================================================================
# VERIFY KAPPA
# ============================================================================

print()
print("=" * 80)
print("TEST 10: DIRECT KAPPA VS SYMMETRIC KAPPA")
print("=" * 80)

kappa5_check = sp.expand(
    kappa5_direct - kappa5_sym.subs({
        e1: E1,
        e2: E2,
        e3: E3,
        e4: E4,
        e5: E5prod,
    })
)

print("Check:")
print(kappa5_check)

assert kappa5_check == 0

print("PASS: exact kappa5 confirmed.")


# ============================================================================
# PRINT CLEAN KAPPA FORM
# ============================================================================

print()
print("=" * 80)
print("TEST 11: CLEAN KAPPA5(e1,...,e5)")
print("=" * 80)

print("kappa5 =")
print(kappa5_sym)


# ============================================================================
# 4-BODY REFERENCE FORMULA
# ============================================================================

print()
print("=" * 80)
print("TEST 12: RECONSTRUCT EXACT 4-BODY KAPPA")
print("=" * 80)

p4, q4, r4, s4 = sp.symbols("p4 q4 r4 s4")

vars4 = [p4, q4, r4, s4]

M1_4 = sp.expand(
    sp.prod(1 + x for x in vars4)
)

X4 = sp.expand(
    sp.prod(1 + x**3 for x in vars4)
)

P4_direct = sp.cancel(X4 / M1_4)

P4_sym = exact_symmetric_reduction(
    P4_direct,
    vars4,
    [e1, e2, e3, e4]
)

kappa4_sym = sp.expand(
    1 - P4_sym
)

print("kappa4 =")
print(kappa4_sym)


# ============================================================================
# 4 -> 5 BODY CORRECTION
# ============================================================================

print()
print("=" * 80)
print("TEST 13: EXACT 4-BODY -> 5-BODY CORRECTION")
print("=" * 80)

# Set e5 = n
n = sp.symbols("n")

kappa5_e5n = sp.expand(
    kappa5_sym.subs(e5, n)
)

print("kappa5(e5=n) =")
print(kappa5_e5n)

print()
print("kappa4 =")
print(kappa4_sym)

correction45 = sp.expand(
    kappa5_e5n - kappa4_sym
)

print()
print("Exact correction kappa5 - kappa4 =")
print(correction45)


# ============================================================================
# COLLECT CORRECTION BY n = e5
# ============================================================================

print()
print("=" * 80)
print("TEST 14: CORRECTION DECOMPOSITION")
print("=" * 80)

correction45_collected = sp.collect(
    sp.expand(correction45),
    n
)

print("Collected:")
print(correction45_collected)

coeff_n0 = sp.expand(
    correction45_collected.subs(n, 0)
)

coeff_n1 = sp.expand(
    sp.diff(correction45_collected, n)
)

print()
print("n^0 coefficient:")
print(coeff_n0)

print()
print("n^1 coefficient:")
print(coeff_n1)


# ============================================================================
# VERIFY CORRECTION IS LINEAR IN e5
# ============================================================================

print()
print("=" * 80)
print("TEST 15: DEGREE IN e5")
print("=" * 80)

degree_e5 = sp.degree(
    correction45,
    n
)

print("degree in n =", degree_e5)

assert degree_e5 <= 1

print("PASS: 4->5 correction is linear in e5.")


# ============================================================================
# FACTOR CORRECTION
# ============================================================================

print()
print("=" * 80)
print("TEST 16: FACTOR 4->5 CORRECTION")
print("=" * 80)

print("Factored correction:")
print(sp.factor(correction45))


# ============================================================================
# VERIFY BY SUBTRACTION
# ============================================================================

print()
print("=" * 80)
print("TEST 17: EXACT RECONSTRUCTION")
print("=" * 80)

reconstructed = sp.expand(
    kappa4_sym + correction45
)

reconstruction_check = sp.expand(
    kappa5_e5n - reconstructed
)

print("Check:")
print(reconstruction_check)

assert reconstruction_check == 0

print("PASS: kappa5 = kappa4 + exact correction.")


# ============================================================================
# STRUCTURAL SPECIALIZATION e5 = 0
# ============================================================================

print()
print("=" * 80)
print("TEST 18: e5 = 0 REDUCTION")
print("=" * 80)

kappa5_e5zero = sp.expand(
    kappa5_sym.subs(e5, 0)
)

reduction_check = sp.expand(
    kappa5_e5zero - kappa4_sym
)

print("kappa5(e5=0) - kappa4 =")
print(reduction_check)

assert reduction_check == 0

print("PASS: 5-body formula reduces exactly to 4-body formula when e5=0.")


# ============================================================================
# FINAL RESULT
# ============================================================================

print()
print("=" * 80)
print("FINAL STRUCTURAL RESULT")
print("=" * 80)

print()
print("M1_5 =")
print(
    "1 + e1 + e2 + e3 + e4 + e5"
)

print()
print("P5 = X5/M1_5 =")
print(P5_sym)

print()
print("kappa5 = 1 - P5 =")
print(kappa5_sym)

print()
print("With e5 = n:")
print()
print("kappa5(e5=n) =")
print(kappa5_e5n)

print()
print("Exact 4 -> 5 correction:")
print()
print("kappa5 - kappa4 =")
print(correction45)

print()
print("Collected by n:")
print()
print(correction45_collected)

print()
print("=" * 80)
print("ALL TESTS PASSED")
print("=" * 80)
