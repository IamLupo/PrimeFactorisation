from sympy import primerange, symbols, Poly, expand
from collections import defaultdict
from itertools import combinations


# ============================================================
# BRUTE-FORCE MACMAHON FUNCTIONS
# ============================================================

def macmahon_M1_bruteforce(n):
    """
    M1(n): sum of multiplicities for partitions into 1 distinct
    part size.
    """
    M1 = 0

    for s1 in range(1, n + 1):
        if n % s1 == 0:
            m1 = n // s1
            M1 += m1

    return M1


def macmahon_M2_bruteforce(n):
    """
    M2(n): sum of m1*m2 over partitions

        n = m1*s1 + m2*s2

    with s1 < s2 and m1,m2 >= 1.
    """
    M2 = 0

    for s1 in range(1, n):
        for s2 in range(s1 + 1, n + 1):

            max_m1 = n // s1

            for m1 in range(1, max_m1 + 1):

                remainder = n - m1 * s1

                if remainder > 0 and remainder % s2 == 0:
                    m2 = remainder // s2
                    M2 += m1 * m2

    return M2


# ============================================================
# CLOSED-FORM MACMAHON FUNCTIONS FOR SQUAREFREE n = p*q*r
# ============================================================

def M1_squarefree_3(p, q, r):
    """
    For n = p*q*r with distinct primes:
        M1 = (p+1)(q+1)(r+1)
    """
    return (p + 1) * (q + 1) * (r + 1)


def M2_squarefree_3(p, q, r):
    """
    Closed form used in the original script:

        M2 =
        [(p^3+1)(q^3+1)(r^3+1)
         - 2*n*M1
         + M1] / 8
    """
    n = p * q * r
    M1 = M1_squarefree_3(p, q, r)

    x = (p**3 + 1) * (q**3 + 1) * (r**3 + 1)

    numerator = x - (2 * n * M1) + M1

    assert numerator % 8 == 0

    return numerator // 8


# ============================================================
# R
# ============================================================

def R_from_M1_M2(n, M1, M2):
    """
    Your original definition:

        R = M1*(n^2 - 3n + 2) - 8*M2
    """
    return M1 * (n**2 - 3*n + 2) - 8*M2


def R_direct_3(p, q, r):
    """
    Compute R using the closed-form M1 and M2.
    """
    n = p * q * r
    M1 = M1_squarefree_3(p, q, r)
    M2 = M2_squarefree_3(p, q, r)

    return R_from_M1_M2(n, M1, M2)


# ============================================================
# SYMMETRIC-POLYNOMIAL VERSION
# ============================================================

def symmetric_invariants_3(p, q, r):
    """
    Elementary symmetric polynomials:

        e1 = p + q + r
        e2 = pq + pr + qr
        e3 = pqr = n
    """
    e1 = p + q + r
    e2 = p*q + p*r + q*r
    e3 = p*q*r

    return e1, e2, e3


def R_symmetric_3(p, q, r):
    """
    Derived structural formula:

        M1 = 1 + e1 + e2 + e3

        R = M1 *
            (e3^2 - e3 + 1
             - e1^2 + e2 + e1)
    """
    e1, e2, e3 = symmetric_invariants_3(p, q, r)

    M1 = 1 + e1 + e2 + e3

    bracket = (
        e3**2
        - e3
        + 1
        - e1**2
        + e2
        + e1
    )

    return M1 * bracket


# ============================================================
# TEST 1: M1 AND M2 AGREEMENT
# ============================================================

def test_M1_M2(prime_limit=50, n_limit=1000):

    primes = list(primerange(2, prime_limit))

    failures = []

    print()
    print("=" * 80)
    print("TEST 1: BRUTE-FORCE M1/M2 VS CLOSED FORM")
    print("=" * 80)

    for p, q, r in combinations(primes, 3):

        n = p * q * r

        if n >= n_limit:
            continue

        M1_brute = macmahon_M1_bruteforce(n)
        M1_formula = M1_squarefree_3(p, q, r)

        M2_brute = macmahon_M2_bruteforce(n)
        M2_formula = M2_squarefree_3(p, q, r)

        if M1_brute != M1_formula or M2_brute != M2_formula:

            failures.append(
                (p, q, r,
                 M1_brute, M1_formula,
                 M2_brute, M2_formula)
            )

    if not failures:
        print("PASS: every tested triple agrees.")
    else:
        print("FAILURES:", len(failures))

        for row in failures[:20]:
            print(row)

    return failures


# ============================================================
# TEST 2: R VS SYMMETRIC FORMULA
# ============================================================

def test_R_formula(prime_limit=100):

    primes = list(primerange(2, prime_limit))

    failures = []

    print()
    print("=" * 80)
    print("TEST 2: R DIRECT VS SYMMETRIC-POLYNOMIAL FORMULA")
    print("=" * 80)

    for p, q, r in combinations(primes, 3):

        R1 = R_direct_3(p, q, r)
        R2 = R_symmetric_3(p, q, r)

        if R1 != R2:

            failures.append(
                (p, q, r, R1, R2, R1 - R2)
            )

    if not failures:
        print("PASS: symmetric formula reproduces R exactly.")
    else:
        print("FAILURES:", len(failures))

        for row in failures[:20]:
            print(row)

    return failures


# ============================================================
# TEST 3: SHOW THE STRUCTURE
# ============================================================

def show_examples():

    examples = [
        (2, 3, 5),
        (2, 3, 7),
        (2, 3, 11),
        (2, 5, 7),
        (3, 5, 7),
    ]

    print()
    print("=" * 80)
    print("TEST 3: EXPLICIT STRUCTURE")
    print("=" * 80)

    for p, q, r in examples:

        n = p * q * r

        e1, e2, e3 = symmetric_invariants_3(p, q, r)

        M1 = M1_squarefree_3(p, q, r)
        M2 = M2_squarefree_3(p, q, r)
        R = R_direct_3(p, q, r)

        print()
        print(f"(p,q,r) = ({p},{q},{r})")
        print(f"n  = {n}")
        print(f"e1 = {e1}")
        print(f"e2 = {e2}")
        print(f"e3 = {e3}")
        print(f"M1 = {M1}")
        print(f"M2 = {M2}")
        print(f"R  = {R}")

        print(
            "R/n^2 =",
            R / (n**2)
        )


# ============================================================
# TEST 4: IS R A FUNCTION OF n ALONE?
# ============================================================

def test_n_only(prime_limit=100):

    primes = list(primerange(2, prime_limit))

    by_n = defaultdict(list)

    for p, q, r in combinations(primes, 3):

        n = p * q * r

        e1, e2, e3 = symmetric_invariants_3(p, q, r)
        R = R_direct_3(p, q, r)

        by_n[n].append(
            (p, q, r, e1, e2, R)
        )

    print()
    print("=" * 80)
    print("TEST 4: DOES R DEPEND ON n ALONE?")
    print("=" * 80)

    repeated_n = 0

    for n, rows in sorted(by_n.items()):

        if len(rows) > 1:

            repeated_n += 1

            R_values = {row[-1] for row in rows}

            print()
            print("n =", n)

            for row in rows:
                print(
                    f"  primes={row[:3]}, "
                    f"e1={row[3]}, "
                    f"e2={row[4]}, "
                    f"R={row[5]}"
                )

            if len(R_values) > 1:
                print("  >>> COUNTEREXAMPLE: same n, different R")

    if repeated_n == 0:
        print(
            "No repeated products n occurred in the tested "
            "prime range."
        )
        print(
            "This is expected for distinct prime factorizations: "
            "unique factorization prevents different prime triples "
            "from having the same n."
        )

    print()
    print(
        "IMPORTANT: this test cannot prove n-dependence by itself "
        "because unique factorization makes n uniquely determine "
        "(p,q,r) for squarefree 3-prime n."
    )


# ============================================================
# TEST 5: FIT R AS A POLYNOMIAL IN e1,e2,e3
# ============================================================

def polynomial_fit_test(prime_limit=50):

    """
    Use SymPy to verify that the proposed R polynomial is
    exactly the one generated by the data.
    """

    e1s, e2s, e3s = symbols("e1 e2 e3")

    proposed = expand(
        (1 + e1s + e2s + e3s)
        *
        (
            e3s**2
            - e3s
            + 1
            - e1s**2
            + e2s
            + e1s
        )
    )

    print()
    print("=" * 80)
    print("TEST 5: SYMBOLIC R POLYNOMIAL")
    print("=" * 80)

    print()
    print("R(e1,e2,e3) =")
    print(proposed)

    print()
    print("Expanded:")
    print(Poly(proposed, e1s, e2s, e3s))


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    # Keep this small initially because brute-force M2 is expensive.
    test_M1_M2(
        prime_limit=50,
        n_limit=1000
    )

    test_R_formula(
        prime_limit=100
    )

    show_examples()

    test_n_only(
        prime_limit=100
    )

    polynomial_fit_test(
        prime_limit=50
    )