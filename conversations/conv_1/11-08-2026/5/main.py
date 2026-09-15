import sympy as sp
from sympy import primerange
from itertools import combinations


# ============================================================================
# 4-BODY STRUCTURAL FORMULA INVESTIGATION
# ============================================================================

print("=" * 80)
print("4-BODY STRUCTURAL FORMULA INVESTIGATION")
print("=" * 80)


# ----------------------------------------------------------------------------
# SYMBOLIC VARIABLES
# ----------------------------------------------------------------------------

p, q, r, s = sp.symbols("p q r s")
e1, e2, e3, e4, n = sp.symbols("e1 e2 e3 e4 n")


# ----------------------------------------------------------------------------
# ELEMENTARY SYMMETRIC POLYNOMIALS
# ----------------------------------------------------------------------------

E1 = p + q + r + s

E2 = (
    p*q + p*r + p*s
    + q*r + q*s
    + r*s
)

E3 = (
    p*q*r + p*q*s
    + p*r*s + q*r*s
)

E4 = p*q*r*s


# ----------------------------------------------------------------------------
# STRUCTURAL R4
#
# From:
#
# R = M1 * (n^2 - 3n + 2) - 8 M2
#
# and
#
# M2 = (X - 2 n M1 + M1) / 8
#
# where
#
# X = product(p_i^3 + 1)
#
# therefore
#
# R = M1 * (n^2 - n + 1) - X
# ----------------------------------------------------------------------------

M1 = (p + 1) * (q + 1) * (r + 1) * (s + 1)

X = (
    (p**3 + 1)
    * (q**3 + 1)
    * (r**3 + 1)
    * (s**3 + 1)
)

R = sp.expand(
    M1 * (E4**2 - E4 + 1) - X
)


# ============================================================================
# TEST 1: SYMMETRIC REDUCTION
# ============================================================================

print()
print("=" * 80)
print("TEST 1: SYMMETRIC REDUCTION")
print("=" * 80)

# Substitute the known elementary symmetric variables.
#
# We ask SymPy to express R through symmetric reduction.

R_sym = sp.symmetrize(
    R,
    [p, q, r, s],
    formal=True
)

R_symmetric = sp.expand(R_sym[0])

print()
print("R in elementary symmetric variables:")
print(R_symmetric)

print()
print("Remainder:")
print(R_sym[1])


# ============================================================================
# TEST 2: SUBSTITUTE e4 = n
# ============================================================================

print()
print("=" * 80)
print("TEST 2: ELIMINATE e4 USING e4 = n")
print("=" * 80)

# Extract the symmetric expression and replace s1,s2,s3,s4
s1, s2, s3, s4 = R_sym[2]

R_e = R_symmetric.subs({
    s1: e1,
    s2: e2,
    s3: e3,
    s4: e4
})

R_n = sp.factor(
    sp.expand(R_e.subs(e4, n))
)

print()
print("R(e1,e2,e3,n) =")
print(R_n)


# ============================================================================
# TEST 3: FACTOR OUT M1
# ============================================================================

print()
print("=" * 80)
print("TEST 3: FACTOR OUT M1")
print("=" * 80)

M1_e = e1 + e2 + e3 + e4 + 1

Q4 = sp.factor(
    sp.cancel(R_e / M1_e)
)

print()
print("Q4 = R4 / M1:")
print(Q4)

print()
print("After e4 = n:")
print(
    sp.factor(
        Q4.subs(e4, n)
    )
)


# ============================================================================
# TEST 4: COMPLETE EXPANSION AFTER e4 = n
# ============================================================================

print()
print("=" * 80)
print("TEST 4: EXPANDED Q4(e1,e2,e3,n)")
print("=" * 80)

Q4_n = sp.expand(
    Q4.subs(e4, n)
)

print()
print(Q4_n)


# ============================================================================
# TEST 5: COLLECT BY n
# ============================================================================

print()
print("=" * 80)
print("TEST 5: COLLECT Q4 BY n")
print("=" * 80)

Q4_collected = sp.collect(
    Q4_n,
    n
)

print()
print(Q4_collected)


# ============================================================================
# TEST 6: COMPARE WITH 3-BODY FORMULA
# ============================================================================

print()
print("=" * 80)
print("TEST 6: 3-BODY VS 4-BODY STRUCTURE")
print("=" * 80)

# 3-body formula already established:
#
# R3 =
# -(e1 + e2 + e3 + 1)
# *
# (
#     e1^2 - e1 e2 - e1 e3 - e1
#     + e2^2 - e2 e3 - e2
#     + 3 e3
# )
#
# Therefore
#
# Q3 =
# -e1^2 + e1e2 + e1e3 + e1
# -e2^2 + e2e3 + e2
# -3e3

Q3 = (
    -e1**2
    + e1*e2
    + e1*e3
    + e1
    - e2**2
    + e2*e3
    + e2
    - 3*e3
)

print()
print("Q3 =")
print(sp.expand(Q3))

print()
print("Q4 =")
print(Q4_n)


# ============================================================================
# TEST 7: TRY TO IDENTIFY THE NEW 4-BODY TERMS
# ============================================================================

print()
print("=" * 80)
print("TEST 7: Q4 - Q3")
print("=" * 80)

# Q3 does not have the same variables as Q4, but formally comparing
# the first three elementary symmetric variables exposes the structural
# correction caused by the fourth body.

difference = sp.expand(Q4_n - Q3)

print()
print("Q4 - Q3 =")
print(sp.factor(difference))


# ============================================================================
# TEST 8: HOMOGENEOUS COMPONENTS
# ============================================================================

print()
print("=" * 80)
print("TEST 8: HOMOGENEOUS COMPONENTS OF Q4")
print("=" * 80)

poly_Q4 = sp.Poly(Q4_n, e1, e2, e3, n)

components = {}

for monomial, coefficient in poly_Q4.terms():

    degree = sum(monomial)

    components.setdefault(degree, 0)
    components[degree] += coefficient * (
        e1 ** monomial[0]
        * e2 ** monomial[1]
        * e3 ** monomial[2]
        * n  ** monomial[3]
    )

for degree in sorted(components, reverse=True):
    print()
    print(f"Degree {degree}:")
    print(sp.factor(components[degree]))


# ============================================================================
# TEST 9: NUMERICAL VERIFICATION AT LARGER RANGE
# ============================================================================

print()
print("=" * 80)
print("TEST 9: NUMERICAL VERIFICATION")
print("=" * 80)


def elementary_symmetric_4(values):
    p, q, r, s = values

    e1 = p + q + r + s

    e2 = (
        p*q + p*r + p*s
        + q*r + q*s + r*s
    )

    e3 = (
        p*q*r + p*q*s
        + p*r*s + q*r*s
    )

    e4 = p*q*r*s

    return e1, e2, e3, e4


def R_structural(values):
    p, q, r, s = values

    n = p*q*r*s

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

    return M1 * (n*n - n + 1) - X


def R_polynomial(values):
    e1v, e2v, e3v, e4v = elementary_symmetric_4(values)

    nv = e4v

    return int(
        R_n.subs({
            e1: e1v,
            e2: e2v,
            e3: e3v,
            n: nv
        })
    )


primes = list(primerange(2, 100))

tested = 0
failures = []

for values in combinations(primes, 4):

    p0, q0, r0, s0 = values
    nv = p0*q0*r0*s0

    if nv > 100000:
        continue

    tested += 1

    a = R_structural(values)
    b = R_polynomial(values)

    if a != b:
        failures.append(
            (values, a, b, a - b)
        )


print()
print(f"Tested: {tested}")
print(f"Failures: {len(failures)}")

if failures:
    for item in failures[:20]:
        print(item)
else:
    print("PASS: polynomial agrees with structural formula.")


# ============================================================================
# TEST 10: SPECIALIZATION e4 = n
# ============================================================================

print()
print("=" * 80)
print("TEST 10: SPECIALIZED FACTORIZATION")
print("=" * 80)

specialized = sp.factor(
    R_n
)

print()
print("R4 =")
print(specialized)

print()
print("R4 / (e1+e2+e3+n+1) =")
print(
    sp.factor(
        sp.cancel(
            R_n / (e1 + e2 + e3 + n + 1)
        )
    )
)


# ============================================================================
# DONE
# ============================================================================

print()
print("=" * 80)
print("DONE")
print("=" * 80)

