import sympy as sp
from itertools import combinations
from math import prod
from sympy import primerange

from m import (
    macmahon_M1,
    macmahon_M2,
    calc_M1,
    calc_M2,
    calc_x,
    get_R,
)


# =============================================================================
# 4-BODY STRUCTURAL FORMULA — DIRECT POLYNOMIAL TEST
# =============================================================================

print("=" * 80)
print("4-BODY STRUCTURAL FORMULA — DIRECT POLYNOMIAL TEST")
print("=" * 80)


# =============================================================================
# ELEMENTARY SYMMETRIC VARIABLES
# =============================================================================

e1, e2, e3, e4 = sp.symbols("e1 e2 e3 e4")
p, q, r, s = sp.symbols("p q r s")


def elementary_4(values):
    p, q, r, s = values

    e1 = p + q + r + s

    e2 = (
        p*q + p*r + p*s
        + q*r + q*s + r*s
    )

    e3 = (
        p*q*r
        + p*q*s
        + p*r*s
        + q*r*s
    )

    e4 = p*q*r*s

    return e1, e2, e3, e4


# =============================================================================
# CORRECT 4-BODY Q POLYNOMIAL
# =============================================================================

Q4 = (
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
)

M1_e = e1 + e2 + e3 + e4 + 1

R4_formula = sp.expand(-M1_e * Q4)


# =============================================================================
# TEST 1: SYMBOLIC FACTORIZATION
# =============================================================================

print("\n" + "=" * 80)
print("TEST 1: SYMBOLIC FACTORIZATION")
print("=" * 80)

print("\nQ4 =")
print(sp.factor(Q4))

print("\nM1 =")
print(M1_e)

print("\nR4 =")
print(sp.factor(R4_formula))

check_factor = sp.simplify(
    R4_formula + M1_e * Q4
)

print("\nCheck R4 + M1*Q4:")
print(check_factor)

assert check_factor == 0

print("\nPASS: R4 = -M1*Q4 exactly.")


# =============================================================================
# TEST 2: VERIFY Q4 DIRECTLY FROM p,q,r,s
# =============================================================================

print("\n" + "=" * 80)
print("TEST 2: SYMMETRIC POLYNOMIAL VS RAW VARIABLES")
print("=" * 80)

e1_raw, e2_raw, e3_raw, e4_raw = elementary_4((p, q, r, s))

Q4_raw = sp.expand(Q4.subs({
    e1: e1_raw,
    e2: e2_raw,
    e3: e3_raw,
    e4: e4_raw,
}))

print("\nQ4(p,q,r,s) =")
print(sp.factor(Q4_raw))

print("\nExpanded:")
print(Q4_raw)


# =============================================================================
# TEST 3: NUMERICAL STRUCTURAL VERIFICATION
# =============================================================================

print("\n" + "=" * 80)
print("TEST 3: NUMERICAL R VS POLYNOMIAL")
print("=" * 80)

primes = list(primerange(2, 100))

cases = []

for combo in combinations(primes, 4):

    n = prod(combo)

    if n <= 5000:
        cases.append(combo)

print(f"TEST CASES = {len(cases)}")

failures = []

for values in cases:

    p0, q0, r0, s0 = values

    n = p0 * q0 * r0 * s0

    # Factorized M1
    factors = {
        p0: 1,
        q0: 1,
        r0: 1,
        s0: 1,
    }

    M1 = calc_M1(factors)

    M2 = calc_M2(n, factors, M1)

    R_direct = get_R(n, M1, M2)

    a, b, c, d = elementary_4(values)

    R_poly = int(
        R4_formula.subs({
            e1: a,
            e2: b,
            e3: c,
            e4: d,
        })
    )

    if R_direct != R_poly:
        failures.append(
            (values, R_direct, R_poly, R_direct - R_poly)
        )


if failures:

    print("\nFAILURES:", len(failures))

    for failure in failures[:20]:
        print(failure)

else:

    print("\nPASS: every numerical test agrees.")


# =============================================================================
# TEST 4: VERIFY M1 FACTOR
# =============================================================================

print("\n" + "=" * 80)
print("TEST 4: M1 FACTOR")
print("=" * 80)

M1_expected = e1 + e2 + e3 + e4 + 1

factor_check = sp.factor(
    R4_formula / (-Q4)
)

print("\nR4 / (-Q4) =")
print(factor_check)

assert sp.simplify(
    factor_check - M1_expected
) == 0

print("\nPASS: exact M1 factor confirmed.")


# =============================================================================
# TEST 5: SUBSTITUTE e4 = n
# =============================================================================

print("\n" + "=" * 80)
print("TEST 5: SUBSTITUTE e4 = n")
print("=" * 80)

n = sp.symbols("n")

Q4_n = sp.expand(
    Q4.subs(e4, n)
)

R4_n = sp.expand(
    R4_formula.subs(e4, n)
)

M1_n = sp.expand(
    M1_e.subs(e4, n)
)

print("\nM1(e1,e2,e3,n) =")
print(M1_n)

print("\nQ4(e1,e2,e3,n) =")
print(sp.factor(Q4_n))

print("\nExpanded Q4 =")
print(Q4_n)

print("\nR4(e1,e2,e3,n) =")
print(sp.factor(R4_n))


# =============================================================================
# TEST 6: CHECK WHETHER Q4 IS POLYNOMIAL
# =============================================================================

print("\n" + "=" * 80)
print("TEST 6: Q4 POLYNOMIAL CHECK")
print("=" * 80)

poly_Q4 = sp.Poly(
    Q4_n,
    e1,
    e2,
    e3,
    n,
)

print("\nTotal degree:")
print(poly_Q4.total_degree())

print("\nNumber of terms:")
print(len(poly_Q4.terms()))

print("\nPASS: Q4 is a genuine polynomial.")


# =============================================================================
# TEST 7: HOMOGENEOUS COMPONENTS
# =============================================================================

print("\n" + "=" * 80)
print("TEST 7: HOMOGENEOUS COMPONENTS OF Q4")
print("=" * 80)

terms_by_degree = {}

for term, coefficient in poly_Q4.terms():

    degree = sum(term)

    terms_by_degree.setdefault(degree, 0)

    terms_by_degree[degree] += 1


for degree in sorted(terms_by_degree):

    component = 0

    for term, coefficient in poly_Q4.terms():

        if sum(term) == degree:

            monomial = (
                e1 ** term[0]
                * e2 ** term[1]
                * e3 ** term[2]
                * n ** term[3]
            )

            component += coefficient * monomial

    print(f"\nDegree {degree}:")
    print(sp.factor(component))


# =============================================================================
# TEST 8: CUBIC PRODUCT
# =============================================================================

print("\n" + "=" * 80)
print("TEST 8: CUBIC PRODUCT IDENTITY")
print("=" * 80)

X_raw = (
    (p**3 + 1)
    * (q**3 + 1)
    * (r**3 + 1)
    * (s**3 + 1)
)

X_sym = sp.expand(
    X_raw.subs({
        p: p,
        q: q,
        r: r,
        s: s,
    })
)

X_e = sp.expand(
    (e1 + e2 + e3 + e4 + 1)
    * (
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
)

X_e_raw = sp.expand(
    X_e.subs({
        e1: e1_raw,
        e2: e2_raw,
        e3: e3_raw,
        e4: e4_raw,
    })
)

print("\nCheck:")
print(
    sp.simplify(X_raw - X_e_raw)
)

assert sp.simplify(
    X_raw - X_e_raw
) == 0

print("\nPASS: cubic product has exact elementary-symmetric form.")


# =============================================================================
# TEST 9: LOOK FOR A SIMPLE RELATION BETWEEN X AND R
# =============================================================================

print("\n" + "=" * 80)
print("TEST 9: X - R RELATION")
print("=" * 80)

X_e = sp.expand(X_e)

difference = sp.factor(
    X_e - R4_formula
)

print("\nX - R =")
print(difference)

print("\nCollected by e4:")
print(
    sp.collect(
        sp.expand(difference),
        e4
    )
)


# =============================================================================
# TEST 10: EXPLICIT NUMERICAL EXAMPLES
# =============================================================================

print("\n" + "=" * 80)
print("TEST 10: EXPLICIT EXAMPLES")
print("=" * 80)

examples = [
    (2, 3, 5, 7),
    (2, 3, 5, 11),
    (2, 3, 7, 11),
    (2, 5, 7, 11),
    (3, 5, 7, 11),
    (5, 7, 11, 13),
]

for values in examples:

    a, b, c, d = values

    n0 = a * b * c * d

    factors = {
        a: 1,
        b: 1,
        c: 1,
        d: 1,
    }

    M1_0 = calc_M1(factors)
    M2_0 = calc_M2(n0, factors, M1_0)

    R0 = get_R(n0, M1_0, M2_0)

    vals = elementary_4(values)

    R_formula_0 = int(
        R4_formula.subs({
            e1: vals[0],
            e2: vals[1],
            e3: vals[2],
            e4: vals[3],
        })
    )

    Q4_0 = int(
        Q4.subs({
            e1: vals[0],
            e2: vals[1],
            e3: vals[2],
            e4: vals[3],
        })
    )

    print("\nprimes =", values)
    print("n      =", n0)
    print("e1     =", vals[0])
    print("e2     =", vals[1])
    print("e3     =", vals[2])
    print("e4     =", vals[3])
    print("M1     =", M1_0)
    print("M2     =", M2_0)
    print("Q4     =", Q4_0)
    print("R      =", R0)
    print("-M1*Q4 =", -M1_0 * Q4_0)

    assert R0 == -M1_0 * Q4_0
    assert R0 == R_formula_0


# =============================================================================
# DONE
# =============================================================================

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)

