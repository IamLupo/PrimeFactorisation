import sympy as sp
from itertools import combinations
from sympy import primerange


# =============================================================================
# 4-BODY STRUCTURAL FORMULA — CLEAN SYMBOLIC DERIVATION
# =============================================================================

print("=" * 80)
print("4-BODY STRUCTURAL FORMULA — CLEAN DERIVATION")
print("=" * 80)


# =============================================================================
# SYMBOLS
# =============================================================================

p, q, r, s = sp.symbols("p q r s")

e1, e2, e3, e4 = sp.symbols(
    "e1 e2 e3 e4"
)

n = sp.symbols("n")


# =============================================================================
# ELEMENTARY SYMMETRIC POLYNOMIALS
# =============================================================================

E1 = p + q + r + s

E2 = (
    p*q + p*r + p*s
    + q*r + q*s + r*s
)

E3 = (
    p*q*r
    + p*q*s
    + p*r*s
    + q*r*s
)

E4 = p*q*r*s


# =============================================================================
# STRUCTURAL INGREDIENTS
# =============================================================================

M1 = (
    (p + 1)
    * (q + 1)
    * (r + 1)
    * (s + 1)
)

X = (
    (p**3 + 1)
    * (q**3 + 1)
    * (r**3 + 1)
    * (s**3 + 1)
)

N = E4


# =============================================================================
# STRUCTURAL R
#
# R = M1 (n^2 - n + 1) - X
# =============================================================================

R = sp.expand(
    M1 * (N**2 - N + 1) - X
)


# =============================================================================
# TEST 1
# DIRECT SYMMETRIC REDUCTION
# =============================================================================

print()
print("=" * 80)
print("TEST 1: SYMMETRIC REDUCTION")
print("=" * 80)

sym_result = sp.symmetrize(
    R,
    [p, q, r, s],
    formal=True
)

R_formal = sp.expand(sym_result[0])
remainder = sp.expand(sym_result[1])

formal_symbols = sym_result[2]

S1, S2, S3, S4 = formal_symbols


print()
print("Formal symmetric variables:")
print(S1, S2, S3, S4)

print()
print("R in formal elementary symmetric variables:")
print(sp.factor(R_formal))

print()
print("Remainder:")
print(remainder)


# =============================================================================
# MAP SymPy'S FORMAL VARIABLES TO OUR VARIABLES
# =============================================================================

R_e = sp.expand(
    R_formal.subs({
        S1: e1,
        S2: e2,
        S3: e3,
        S4: e4
    })
)

print()
print("=" * 80)
print("TEST 2: R(e1,e2,e3,e4)")
print("=" * 80)

print()
print(sp.factor(R_e))


# =============================================================================
# M1 IN ELEMENTARY SYMMETRIC VARIABLES
# =============================================================================

M1_e = (
    1 + e1 + e2 + e3 + e4
)


# =============================================================================
# TEST 3
# FACTOR BY M1
# =============================================================================

print()
print("=" * 80)
print("TEST 3: FACTOR BY M1")
print("=" * 80)

Q4 = sp.factor(
    sp.cancel(
        R_e / M1_e
    )
)

print()
print("M1 =")
print(M1_e)

print()
print("Q4 = R/M1 =")
print(Q4)

print()
print("Check R - M1*Q4:")
print(
    sp.simplify(
        R_e - M1_e * Q4
    )
)


# =============================================================================
# TEST 4
# SUBSTITUTE e4 = n
# =============================================================================

print()
print("=" * 80)
print("TEST 4: SUBSTITUTE e4 = n")
print("=" * 80)

R_n = sp.factor(
    R_e.subs(e4, n)
)

Q4_n = sp.factor(
    Q4.subs(e4, n)
)

M1_n = (
    1 + e1 + e2 + e3 + n
)

print()
print("M1(n) =")
print(M1_n)

print()
print("R(e1,e2,e3,n) =")
print(R_n)

print()
print("Q4(e1,e2,e3,n) =")
print(Q4_n)


# =============================================================================
# TEST 5
# EXPANDED Q4
# =============================================================================

print()
print("=" * 80)
print("TEST 5: EXPANDED Q4")
print("=" * 80)

Q4_expanded = sp.expand(Q4_n)

print()
print(Q4_expanded)


# =============================================================================
# TEST 6
# COLLECT BY n
# =============================================================================

print()
print("=" * 80)
print("TEST 6: COLLECT Q4 BY n")
print("=" * 80)

print()
print(
    sp.collect(
        Q4_expanded,
        n
    )
)


# =============================================================================
# TEST 7
# 3-BODY FORMULA
# =============================================================================

print()
print("=" * 80)
print("TEST 7: 3-BODY FORMULA")
print("=" * 80)

a, b, c = sp.symbols(
    "a b c"
)

f1 = a + b + c
f2 = a*b + a*c + b*c
f3 = a*b*c

Q3 = sp.expand(
    -f1**2
    + f1*f2
    + f1*f3
    + f1
    - f2**2
    + f2*f3
    + f2
    - 3*f3
)

R3 = sp.factor(
    (1 + f1 + f2 + f3) * Q3
)

print()
print("Q3 =")
print(Q3)

print()
print("R3 =")
print(R3)


# =============================================================================
# TEST 8
# LOOK FOR A GENERAL PATTERN
# =============================================================================

print()
print("=" * 80)
print("TEST 8: 4-BODY CORRECTION RELATIVE TO 3-BODY")
print("=" * 80)

# We cannot literally subtract Q3 from Q4 because the elementary
# symmetric variables have different meanings.
#
# Instead, embed a 3-body system (p,q,r) into the 4-body system
# by setting s = 1.
#
# Then:
#
# e1(4) = f1 + 1
# e2(4) = f2 + f1
# e3(4) = f3 + f2
# e4(4) = f3
#
# This is a meaningful structural specialization.

Q4_s1 = sp.expand(
    Q4.subs({
        e1: f1 + 1,
        e2: f2 + f1,
        e3: f3 + f2,
        e4: f3
    })
)

print()
print("Q4 under s = 1:")
print(
    sp.factor(Q4_s1)
)

print()
print("Difference from Q3:")
print(
    sp.factor(
        sp.expand(Q4_s1 - Q3)
    )
)


# =============================================================================
# TEST 9
# DERIVE PRODUCT(p^3+1) DIRECTLY
# =============================================================================

print()
print("=" * 80)
print("TEST 9: CUBIC PRODUCT STRUCTURE")
print("=" * 80)

X_sym = sp.symmetrize(
    X,
    [p, q, r, s],
    formal=True
)

X_formal = sp.expand(
    X_sym[0]
)

X_S1, X_S2, X_S3, X_S4 = X_sym[2]

X_e = sp.expand(
    X_formal.subs({
        X_S1: e1,
        X_S2: e2,
        X_S3: e3,
        X_S4: e4
    })
)

print()
print("Product (p^3+1)(q^3+1)(r^3+1)(s^3+1):")
print()
print(
    sp.factor(X_e)
)


# =============================================================================
# TEST 10
# COMPARE X WITH M1(n^2-n+1)-R
# =============================================================================

print()
print("=" * 80)
print("TEST 10: RECONSTRUCT X")
print("=" * 80)

X_reconstructed = sp.expand(
    M1_e * (e4**2 - e4 + 1)
    - R_e
)

print()
print("X reconstructed from R:")
print(
    sp.factor(X_reconstructed)
)

print()
print("Difference:")
print(
    sp.simplify(
        X_e - X_reconstructed
    )
)


# =============================================================================
# TEST 11
# NUMERICAL VERIFICATION
# =============================================================================

print()
print("=" * 80)
print("TEST 11: NUMERICAL VERIFICATION")
print("=" * 80)


def elementary(values):

    p0, q0, r0, s0 = values

    E1v = p0 + q0 + r0 + s0

    E2v = (
        p0*q0
        + p0*r0
        + p0*s0
        + q0*r0
        + q0*s0
        + r0*s0
    )

    E3v = (
        p0*q0*r0
        + p0*q0*s0
        + p0*r0*s0
        + q0*r0*s0
    )

    E4v = p0*q0*r0*s0

    return E1v, E2v, E3v, E4v


def R_direct(values):

    p0, q0, r0, s0 = values

    nv = p0*q0*r0*s0

    m1 = (
        (p0 + 1)
        * (q0 + 1)
        * (r0 + 1)
        * (s0 + 1)
    )

    x = (
        (p0**3 + 1)
        * (q0**3 + 1)
        * (r0**3 + 1)
        * (s0**3 + 1)
    )

    return (
        m1 * (nv**2 - nv + 1)
        - x
    )


def R_formula(values):

    E1v, E2v, E3v, E4v = elementary(values)

    return int(
        R_e.subs({
            e1: E1v,
            e2: E2v,
            e3: E3v,
            e4: E4v
        })
    )


primes = list(
    primerange(2, 100)
)

tested = 0
failures = []

for values in combinations(
    primes,
    4
):

    nv = 1

    for x in values:
        nv *= x

    if nv > 1000000:
        continue

    tested += 1

    direct = R_direct(values)
    formula = R_formula(values)

    if direct != formula:

        failures.append(
            (
                values,
                direct,
                formula,
                direct - formula
            )
        )


print()
print("Tested:", tested)
print("Failures:", len(failures))

if failures:

    print()

    for failure in failures[:20]:
        print(failure)

else:

    print()
    print(
        "PASS: exact 4-body polynomial verified."
    )


# =============================================================================
# TEST 12
# FINAL STRUCTURAL FORM
# =============================================================================

print()
print("=" * 80)
print("TEST 12: FINAL STRUCTURAL FORM")
print("=" * 80)

print()
print("M1 =")
print(
    sp.factor(M1_e)
)

print()
print("Q4 =")
print(
    sp.factor(Q4)
)

print()
print("R4 =")
print(
    sp.factor(
        M1_e * Q4
    )
)

print()
print("=" * 80)
print("DONE")
print("=" * 80)

