from math import isqrt


TEST_CASES = [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
    (37, 71),
    (41, 73),
    (43, 79),
    (47, 83),
    (53, 89),
]


# ============================================================
# NUMBER THEORY
# ============================================================

def sigma_and_sigma3_from_factorization(factors):
    """
    factors = [(p, exponent), ...]

    Returns:

        sigma(n)
        sigma_3(n)

    using

        sigma_k(p^a)
        =
        1 + p^k + p^(2k) + ... + p^(ak)
    """

    sigma = 1
    sigma3 = 1

    for p, a in factors:

        local_sigma = sum(
            p ** j
            for j in range(a + 1)
        )

        local_sigma3 = sum(
            p ** (3 * j)
            for j in range(a + 1)
        )

        sigma *= local_sigma
        sigma3 *= local_sigma3

    return sigma, sigma3


def divisors(n):
    result = []

    r = isqrt(n)

    for d in range(1, r + 1):

        if n % d != 0:
            continue

        result.append(d)

        other = n // d

        if other != d:
            result.append(other)

    return result


def sigma_direct(n):
    return sum(divisors(n))


def sigma3_direct(n):
    return sum(
        d ** 3
        for d in divisors(n)
    )


# ============================================================
# BRUTE FORCE MACMAHON M2
# ============================================================

def macmahon_M2_bruteforce(n):

    total = 0

    for s1 in range(1, n):

        for s2 in range(s1 + 1, n + 1):

            for m1 in range(1, n // s1 + 1):

                remainder = n - m1 * s1

                if remainder <= 0:
                    continue

                if remainder % s2 != 0:
                    continue

                m2 = remainder // s2

                if m2 <= 0:
                    continue

                total += m1 * m2

    return total


# ============================================================
# CLOSED M2 FORMULA
# ============================================================

def macmahon_M2_closed(n, sigma, sigma3):

    numerator = (
        sigma3
        + (1 - 2 * n) * sigma
    )

    assert numerator % 8 == 0

    return numerator // 8


# ============================================================
# ORIGINAL SEMIPRIME COORDINATE FORM
# ============================================================

def macmahon_M2_semiprime(p, q):

    n = p * q

    x = (
        (p ** 3 + 1)
        * (q ** 3 + 1)
    )

    M1 = (
        (p + 1)
        * (q + 1)
    )

    return (
        x
        - 2 * n * M1
        + M1
    ) // 8


# ============================================================
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 10")
print()

print("PART 1: GENERAL CLOSED FORM")
print("-" * 100)

print(
    "n\t"
    "sigma\t"
    "sigma3\t"
    "M2_closed\t"
    "M2_bruteforce\t"
    "OK"
)

GENERAL_LIMIT = 100

for n in range(2, GENERAL_LIMIT + 1):

    sigma = sigma_direct(n)
    sigma3 = sigma3_direct(n)

    M2_closed = macmahon_M2_closed(
        n,
        sigma,
        sigma3
    )

    M2_brute = macmahon_M2_bruteforce(n)

    ok = (
        M2_closed == M2_brute
    )

    print(
        f"{n}\t"
        f"{sigma}\t"
        f"{sigma3}\t"
        f"{M2_closed}\t"
        f"{M2_brute}\t"
        f"{ok}"
    )


print()
print("PART 2: SEMIPRIME SPECIALIZATION")
print("-" * 120)

print(
    "n\tp\tq\t"
    "sigma_factorized\t"
    "sigma3_factorized\t"
    "M2_closed\t"
    "M2_semiprime_formula\t"
    "M2_match"
)

for p, q in TEST_CASES:

    n = p * q

    #
    # Factorized expressions, used only to verify
    # the squarefree semiprime specialization.
    #
    sigma, sigma3 = (
        sigma_and_sigma3_from_factorization(
            [(p, 1), (q, 1)]
        )
    )

    M2_closed = macmahon_M2_closed(
        n,
        sigma,
        sigma3
    )

    M2_semiprime = (
        macmahon_M2_semiprime(p, q)
    )

    match = (
        M2_closed == M2_semiprime
    )

    print(
        f"{n}\t"
        f"{p}\t"
        f"{q}\t"
        f"{sigma}\t"
        f"{sigma3}\t"
        f"{M2_closed}\t"
        f"{M2_semiprime}\t"
        f"{match}"
    )


print()
print("PART 3: DIRECT DIVISOR COMPUTATION")
print("-" * 100)

print(
    "n\t"
    "sigma_direct\t"
    "sigma3_direct\t"
    "M2_closed\t"
    "M2_bruteforce\t"
    "OK"
)

for p, q in TEST_CASES:

    n = p * q

    #
    # IMPORTANT:
    #
    # These sigma values are obtained from n itself
    # by enumerating divisors, not from p and q.
    #
    sigma = sigma_direct(n)
    sigma3 = sigma3_direct(n)

    M2_closed = macmahon_M2_closed(
        n,
        sigma,
        sigma3
    )

    #
    # Brute force is only practical for these
    # relatively small test cases.
    #
    M2_brute = macmahon_M2_bruteforce(n)

    ok = (
        M2_closed == M2_brute
    )

    print(
        f"{n}\t"
        f"{sigma}\t"
        f"{sigma3}\t"
        f"{M2_closed}\t"
        f"{M2_brute}\t"
        f"{ok}"
    )


print()
print("FINISHED EXPERIMENT 10")
