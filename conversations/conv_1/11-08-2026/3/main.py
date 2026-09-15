import sympy as sp
from itertools import combinations
from sympy import primerange


# ============================================================
# BASIC FUNCTIONS
# ============================================================

def brute_M1(n):
    """Brute-force M1."""
    total = 0

    for s in range(1, n + 1):
        if n % s == 0:
            total += n // s

    return total


def brute_M2(n):
    """Brute-force M2."""
    total = 0

    for s1 in range(1, n):
        for s2 in range(s1 + 1, n + 1):
            for m1 in range(1, n // s1 + 1):
                remainder = n - m1 * s1

                if remainder > 0 and remainder % s2 == 0:
                    m2 = remainder // s2
                    total += m1 * m2

    return total


def closed_M1(primes):
    result = 1

    for p in primes:
        result *= p + 1

    return result


def closed_M2(primes):
    """
    Your closed M2 formula:

        8 M2 =
            prod(p^3 + 1)
            - 2 n M1
            + M1
    """

    n = 1
    x = 1

    for p in primes:
        n *= p
        x *= p**3 + 1

    M1 = closed_M1(primes)

    return (x - 2 * n * M1 + M1) // 8


def R_direct(n, M1, M2):
    return M1 * (n**2 - 3*n + 2) - 8*M2


# ============================================================
# ELEMENTARY SYMMETRIC POLYNOMIALS
# ============================================================

def elementary_symmetric(primes):
    """
    Return e1, e2, ..., ek.
    """

    k = len(primes)
    e = []

    for r in range(1, k + 1):
        total = 0

        for subset in combinations(primes, r):
            product = 1

            for x in subset:
                product *= x

            total += product

        e.append(total)

    return e


# ============================================================
# CORRECT 3-BODY STRUCTURAL FORMULA
# ============================================================

def R_three_body_structural(primes):
    """
    Correct formula for three distinct primes.

    e1 = p + q + r
    e2 = pq + pr + qr
    e3 = pqr = n
    """

    if len(primes) != 3:
        raise ValueError("Expected exactly 3 primes")

    e1, e2, e3 = elementary_symmetric(primes)

    return (
        -e1**3
        -e2**3
        +3*e1*e2*e3
        +3*e1*e2
        +e1*e3**2
        -e1*e3
        +e1
        +e2*e3**2
        -e2*e3
        +e2
        -3*e3**2
        -3*e3
    )


# ============================================================
# DIRECT STRUCTURAL IDENTITY
# ============================================================

def R_structural_direct(primes):
    """
    Most useful form:

        R =
          M1 * (n^2 - n + 1)
          - prod(p^3 + 1)
    """

    n = 1
    x = 1

    for p in primes:
        n *= p
        x *= p**3 + 1

    M1 = closed_M1(primes)

    return M1 * (n**2 - n + 1) - x


# ============================================================
# WRONG / SUSPECTED FORMULA
# ============================================================

def R_old_formula(e1, e2, e3):
    """
    The formula produced by the previous test.
    Kept here deliberately so we can demonstrate its failure.
    """

    return (
        -e1**3
        -e1**2*e2
        -e1**2*e3
        +2*e1*e2
        +e1*e3**2
        +2*e1
        +e2**2
        +e2*e3**2
        +2*e2
        +e3**3
        +1
    )


# ============================================================
# TEST 1
# BRUTE FORCE M1/M2
# ============================================================

print("=" * 80)
print("TEST 1: BRUTE-FORCE M1/M2 VS CLOSED FORM")
print("=" * 80)

primes = list(primerange(2, 100))

failures = []

for p, q, r in combinations(primes, 3):

    n = p * q * r

    if n >= 1000:
        continue

    M1_b = brute_M1(n)
    M1_c = closed_M1((p, q, r))

    M2_b = brute_M2(n)
    M2_c = closed_M2((p, q, r))

    if M1_b != M1_c or M2_b != M2_c:
        failures.append(
            ((p, q, r), M1_b, M1_c, M2_b, M2_c)
        )

if not failures:
    print("PASS: all M1/M2 tests agree.")
else:
    print("FAILURES:", len(failures))

    for x in failures[:10]:
        print(x)


# ============================================================
# TEST 2
# R DIRECT VS STRUCTURAL PRODUCT FORM
# ============================================================

print()
print("=" * 80)
print("TEST 2: R DIRECT VS STRUCTURAL PRODUCT FORM")
print("=" * 80)

failures = []

for p, q, r in combinations(primes, 3):

    n = p * q * r

    if n >= 1000:
        continue

    M1 = closed_M1((p, q, r))
    M2 = closed_M2((p, q, r))

    R1 = R_direct(n, M1, M2)
    R2 = R_structural_direct((p, q, r))

    if R1 != R2:
        failures.append(
            ((p, q, r), R1, R2, R1 - R2)
        )

if not failures:
    print("PASS: R structural identity holds for every test.")
else:
    print("FAILURES:", len(failures))

    for x in failures[:20]:
        print(x)


# ============================================================
# TEST 3
# R DIRECT VS ELEMENTARY-SYMMETRIC POLYNOMIAL
# ============================================================

print()
print("=" * 80)
print("TEST 3: R VS CORRECT e1/e2/e3 POLYNOMIAL")
print("=" * 80)

failures = []

for p, q, r in combinations(primes, 3):

    n = p * q * r

    if n >= 1000:
        continue

    M1 = closed_M1((p, q, r))
    M2 = closed_M2((p, q, r))

    R1 = R_direct(n, M1, M2)
    R2 = R_three_body_structural((p, q, r))

    if R1 != R2:
        e1, e2, e3 = elementary_symmetric((p, q, r))

        failures.append(
            ((p, q, r), e1, e2, e3, R1, R2, R1 - R2)
        )

if not failures:
    print("PASS: correct 3-body polynomial holds.")
else:
    print("FAILURES:", len(failures))

    for x in failures[:20]:
        print(x)


# ============================================================
# TEST 4
# OLD FORMULA VS CORRECT FORMULA
# ============================================================

print()
print("=" * 80)
print("TEST 4: OLD PAPER/SYMMETRIC FORMULA DIAGNOSTIC")
print("=" * 80)

failures = []

for p, q, r in combinations(primes, 3):

    n = p * q * r

    if n >= 1000:
        continue

    M1 = closed_M1((p, q, r))
    M2 = closed_M2((p, q, r))

    R_true = R_direct(n, M1, M2)

    e1, e2, e3 = elementary_symmetric((p, q, r))

    R_old = R_old_formula(e1, e2, e3)

    if R_true != R_old:
        failures.append(
            ((p, q, r), R_true, R_old, R_true - R_old)
        )

print("Old formula failures:", len(failures))

for x in failures[:10]:
    print(x)


# ============================================================
# TEST 5
# SYMBOLIC DERIVATION
# ============================================================

print()
print("=" * 80)
print("TEST 5: SYMBOLIC DERIVATION")
print("=" * 80)

e1, e2, e3 = sp.symbols("e1 e2 e3")

M1 = 1 + e1 + e2 + e3

sum_p3 = e1**3 - 3*e1*e2 + 3*e3

sum_pair_cubes = (
    e2**3
    - 3*e1*e2*e3
    + 3*e3**2
)

prod_p3_plus_1 = (
    1
    + sum_p3
    + sum_pair_cubes
    + e3**3
)

R = sp.expand(
    M1 * (e3**2 - e3 + 1)
    - prod_p3_plus_1
)

print()
print("Correct R(e1,e2,e3):")
print()
print(sp.factor(R))
print()
print("Expanded:")
print()
print(R)


# ============================================================
# TEST 6
# SHOW FIRST EXAMPLES
# ============================================================

print()
print("=" * 80)
print("TEST 6: EXPLICIT EXAMPLES")
print("=" * 80)

examples = [
    (2, 3, 5),
    (2, 3, 7),
    (2, 3, 11),
    (2, 5, 7),
    (3, 5, 7),
]

for primes3 in examples:

    p, q, r = primes3

    n = p * q * r
    e1, e2, e3 = elementary_symmetric(primes3)

    M1 = closed_M1(primes3)
    M2 = closed_M2(primes3)

    R = R_direct(n, M1, M2)

    print()
    print(f"(p,q,r) = {primes3}")
    print(f"n       = {n}")
    print(f"e1      = {e1}")
    print(f"e2      = {e2}")
    print(f"e3      = {e3}")
    print(f"M1      = {M1}")
    print(f"M2      = {M2}")
    print(f"R       = {R}")

    print()
    print("Structural:")
    print(
        "R = M1*(n^2-n+1) - product(p^3+1)"
    )

    print(
        "R =", R_structural_direct(primes3)
    )

    print()
    print("3-body polynomial:")
    print(
        "R =",
        R_three_body_structural(primes3)
    )

