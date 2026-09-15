import sympy as sp
from itertools import combinations
from sympy import primerange

from m import calc_M1, calc_M2, calc_x, get_R


# =============================================================================
# KAPPA STRUCTURAL INVESTIGATION — NEXT STAGE
# =============================================================================

print("=" * 80)
print("KAPPA STRUCTURE — PRODUCT / CORRECTION / PAPER DIAGNOSTIC")
print("=" * 80)


# -----------------------------------------------------------------------------
# SYMBOLS
# -----------------------------------------------------------------------------

p, q, r, s = sp.symbols("p q r s")

e1, e2, e3, e4, n = sp.symbols(
    "e1 e2 e3 e4 n"
)


# =============================================================================
# VERIFIED KAPPAS
# =============================================================================

kappa3 = sp.expand(
    -e1**2
    + e1*e2
    + e1*e3
    + e1
    - e2**2
    + e2*e3
    + e2
    - 3*e3
)

kappa4 = sp.expand(
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

kappa4_n = sp.expand(kappa4.subs(e4, n))


# =============================================================================
# TEST 1: CLEAN KAPPA FORM
# =============================================================================

print("\n" + "=" * 80)
print("TEST 1: CLEAN KAPPA FORM")
print("=" * 80)

print("\nkappa3 =")
print(kappa3)

print("\nkappa4 =")
print(kappa4)

print("\nkappa4(e4=n) =")
print(kappa4_n)


# =============================================================================
# TEST 2: KAPPA CORRECTION
# =============================================================================

print("\n" + "=" * 80)
print("TEST 2: KAPPA4 - KAPPA3")
print("=" * 80)

correction = sp.factor(kappa4_n - kappa3)

print("\nkappa4 - kappa3 =")
print(correction)

print("\nCollected by n:")
print(sp.collect(sp.expand(correction), n))

expected_correction = (
    n * (-2*e1 + e2 + e3)
    - e3 * (e3 - 1)
)

check = sp.simplify(
    correction - expected_correction
)

print("\nCheck against compact correction:")
print(check)

if check == 0:
    print("\nPASS: compact kappa correction confirmed.")


# =============================================================================
# TEST 3: EXPRESS THE CORRECTION USING e1,e2,e3,n
# =============================================================================

print("\n" + "=" * 80)
print("TEST 3: CORRECTION STRUCTURE")
print("=" * 80)

print("""
kappa4 = kappa3
         + n*(-2*e1 + e2 + e3)
         - e3*(e3 - 1)
""")

print("\nCoefficient of n:")
print(sp.factor(sp.diff(correction, n)))

print("\nConstant correction:")
print(sp.factor(correction.subs(n, 0)))


# =============================================================================
# TEST 4: CUBIC PRODUCT — SYMMETRIC REDUCTION
# =============================================================================

print("\n" + "=" * 80)
print("TEST 4: CUBIC PRODUCT SYMMETRIC REDUCTION")
print("=" * 80)

X = sp.expand(
    (p**3 + 1)
    * (q**3 + 1)
    * (r**3 + 1)
    * (s**3 + 1)
)

sym_result = sp.symmetrize(
    X,
    [p, q, r, s],
    formal=True
)

# SymPy returns:
# (symmetric_expression, remainder, mapping)

X_sym = sym_result[0]
X_remainder = sym_result[1]
X_mapping = sym_result[2]

print("\nX symmetric:")
print(sp.factor(X_sym))

print("\nRemainder:")
print(X_remainder)

print("\nMapping:")
print(X_mapping)

if X_remainder == 0:
    print("\nPASS: X is completely symmetric.")


# =============================================================================
# TEST 5: REPLACE FORMAL s1..s4 WITH e1..e4
# =============================================================================

print("\n" + "=" * 80)
print("TEST 5: X(e1,e2,e3,e4)")
print("=" * 80)

formal_symbols = [item[0] for item in X_mapping]

formal_to_e = {
    formal_symbols[0]: e1,
    formal_symbols[1]: e2,
    formal_symbols[2]: e3,
    formal_symbols[3]: e4,
}

X_e = sp.expand(
    X_sym.subs(formal_to_e)
)

print("\nX(e1,e2,e3,e4) =")
print(sp.factor(X_e))


# =============================================================================
# TEST 6: STRUCTURAL R RELATION
# =============================================================================

print("\n" + "=" * 80)
print("TEST 6: X = M1*(n^2-n+1) - R")
print("=" * 80)

M1_4 = e1 + e2 + e3 + e4 + 1

R4 = sp.expand(
    M1_4 * kappa4
)

X_from_R = sp.expand(
    M1_4 * (e4**2 - e4 + 1)
    - R4
)

difference_X = sp.factor(
    X_e - X_from_R
)

print("\nX_direct - X_reconstructed =")
print(difference_X)

if difference_X == 0:
    print("\nPASS: cubic product reconstructed exactly from R.")


# =============================================================================
# TEST 7: X/M1
# =============================================================================

print("\n" + "=" * 80)
print("TEST 7: NORMALIZED CUBIC PRODUCT")
print("=" * 80)

# Since X itself contains M1 as a factor, determine the remaining factor.

quotient_X = sp.factor(
    sp.cancel(X_e / M1_4)
)

print("\nX / M1 =")
print(quotient_X)

print("\nCheck:")
print(
    sp.simplify(
        X_e - M1_4 * quotient_X
    )
)


# =============================================================================
# TEST 8: DIRECT RELATION BETWEEN X/M1 AND KAPPA4
# =============================================================================

print("\n" + "=" * 80)
print("TEST 8: X/M1 VS KAPPA4")
print("=" * 80)

normalized_R = sp.expand(
    e4**2 - e4 + 1 - kappa4
)

print("\nX/M1 predicted by R:")
print(sp.factor(normalized_R))

check_normalized = sp.factor(
    quotient_X - normalized_R
)

print("\nDifference:")
print(check_normalized)

if check_normalized == 0:
    print("\nPASS: X/M1 = e4^2 - e4 + 1 - kappa4.")


# =============================================================================
# TEST 9: ELIMINATE e4 = n
# =============================================================================

print("\n" + "=" * 80)
print("TEST 9: NORMALIZED PRODUCT AFTER e4=n")
print("=" * 80)

normalized_n = sp.expand(
    quotient_X.subs(e4, n)
)

print("\nX/M1 =")
print(sp.factor(normalized_n))

print("\nExpanded:")
print(normalized_n)


# =============================================================================
# TEST 10: COMPARE 3-BODY NORMALIZED PRODUCT
# =============================================================================

print("\n" + "=" * 80)
print("TEST 10: 3-BODY NORMALIZED PRODUCT")
print("=" * 80)

X3 = sp.expand(
    (p**3 + 1)
    * (q**3 + 1)
    * (r**3 + 1)
)

sym3_result = sp.symmetrize(
    X3,
    [p, q, r],
    formal=True
)

X3_sym = sym3_result[0]
X3_rem = sym3_result[1]
X3_mapping = sym3_result[2]

print("\nX3 symmetric:")
print(sp.factor(X3_sym))

print("\nRemainder:")
print(X3_rem)

formal3 = [item[0] for item in X3_mapping]

X3_e = sp.expand(
    X3_sym.subs({
        formal3[0]: e1,
        formal3[1]: e2,
        formal3[2]: e3,
    })
)

M1_3 = e1 + e2 + e3 + 1

normalized_X3 = sp.factor(
    sp.cancel(X3_e / M1_3)
)

print("\nX3/M1_3 =")
print(normalized_X3)


# =============================================================================
# TEST 11: 3-BODY KAPPA FROM PRODUCT
# =============================================================================

print("\n" + "=" * 80)
print("TEST 11: RECONSTRUCT KAPPA3 FROM PRODUCT")
print("=" * 80)

# For 3 body:
#
# R3 = M1*(n^2-n+1) - X3
#
# n=e3.
#
# Therefore:
#
# kappa3 = e3^2-e3+1-X3/M1.

kappa3_from_product = sp.expand(
    e3**2
    - e3
    + 1
    - normalized_X3
)

print("\nkappa3 reconstructed:")
print(sp.factor(kappa3_from_product))

print("\nDifference from verified kappa3:")
print(
    sp.factor(
        kappa3_from_product - kappa3
    )
)


# =============================================================================
# TEST 12: THE KEY COMPARISON
# =============================================================================

print("\n" + "=" * 80)
print("TEST 12: KEY STRUCTURAL COMPARISON")
print("=" * 80)

print("""
3-BODY:

    X3/M1_3 = e3^2 - e3 + 1 - kappa3

4-BODY:

    X4/M1_4 = e4^2 - e4 + 1 - kappa4

Therefore the difference between bodies is entirely encoded
in the change of the normalized cubic product.
""")

print("\n4-body normalized product:")
print(
    sp.factor(
        normalized_X3.subs(e3, e4)
    )
)

print("\nActual 4-body normalized product:")
print(
    sp.factor(
        quotient_X
    )
)

print("\nDifference:")
print(
    sp.factor(
        quotient_X
        - normalized_X3.subs(e3, e4)
    )
)


# =============================================================================
# TEST 13: NUMERICAL KAPPA TEST
# =============================================================================

print("\n" + "=" * 80)
print("TEST 13: NUMERICAL KAPPA VERIFICATION")
print("=" * 80)

MAX_PRIME = 50
MAX_N = 20000

primes = list(primerange(2, MAX_PRIME))

tests = [
    combo
    for combo in combinations(primes, 4)
    if combo[0] * combo[1] * combo[2] * combo[3] <= MAX_N
]

print(f"\nMAX_PRIME = {MAX_PRIME}")
print(f"MAX_N     = {MAX_N}")
print(f"TEST CASES = {len(tests)}")

failures = []

for pp, qq, rr, ss in tests:

    nn = pp * qq * rr * ss

    factors = {
        pp: 1,
        qq: 1,
        rr: 1,
        ss: 1,
    }

    M1 = calc_M1(factors)
    M2 = calc_M2(nn, factors, M1)

    R = get_R(nn, M1, M2)

    E1 = pp + qq + rr + ss

    E2 = (
        pp*qq
        + pp*rr
        + pp*ss
        + qq*rr
        + qq*ss
        + rr*ss
    )

    E3 = (
        pp*qq*rr
        + pp*qq*ss
        + pp*rr*ss
        + qq*rr*ss
    )

    E4 = nn

    kap = int(
        kappa4_n.subs({
            e1: E1,
            e2: E2,
            e3: E3,
            n: E4,
        })
    )

    predicted = M1 * kap

    if R != predicted:
        failures.append(
            (pp, qq, rr, ss, R, M1, kap, predicted)
        )

if not failures:
    print("\nPASS: every numerical test satisfies R=M1*kappa4.")
else:
    print(f"\nFAILURES: {len(failures)}")
    for item in failures[:20]:
        print(item)


# =============================================================================
# TEST 14: EXPLICIT VALUES
# =============================================================================

print("\n" + "=" * 80)
print("TEST 14: EXPLICIT KAPPA VALUES")
print("=" * 80)

examples = [
    (2, 3, 5, 7),
    (2, 3, 5, 11),
    (2, 3, 7, 11),
    (2, 5, 7, 11),
    (3, 5, 7, 11),
    (5, 7, 11, 13),
]

for pp, qq, rr, ss in examples:

    nn = pp * qq * rr * ss

    E1 = pp + qq + rr + ss

    E2 = (
        pp*qq
        + pp*rr
        + pp*ss
        + qq*rr
        + qq*ss
        + rr*ss
    )

    E3 = (
        pp*qq*rr
        + pp*qq*ss
        + pp*rr*ss
        + qq*rr*ss
    )

    kap = int(
        kappa4_n.subs({
            e1: E1,
            e2: E2,
            e3: E3,
            n: nn,
        })
    )

    factors = {
        pp: 1,
        qq: 1,
        rr: 1,
        ss: 1,
    }

    M1 = calc_M1(factors)
    M2 = calc_M2(nn, factors, M1)
    R = get_R(nn, M1, M2)

    print(f"\nprimes = {(pp, qq, rr, ss)}")
    print(f"n      = {nn}")
    print(f"e1     = {E1}")
    print(f"e2     = {E2}")
    print(f"e3     = {E3}")
    print(f"e4     = {nn}")
    print(f"M1     = {M1}")
    print(f"kappa4 = {kap}")
    print(f"R      = {R}")
    print(f"M1*k   = {M1 * kap}")


# =============================================================================
# FINAL STRUCTURAL STATEMENT
# =============================================================================

print("\n" + "=" * 80)
print("FINAL STRUCTURAL RESULT")
print("=" * 80)

print("""
VERIFIED:

    R4 = M1 * kappa4

where

    M1 = e1 + e2 + e3 + e4 + 1

and

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

After e4=n:

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

And:

    kappa4 - kappa3
      = n*(-2*e1 + e2 + e3)
        - e3*(e3 - 1)

This is the structural correction that must be explained
when comparing the 3-body and 4-body formulas.
""")

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)