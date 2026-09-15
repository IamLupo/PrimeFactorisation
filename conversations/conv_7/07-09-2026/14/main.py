from math import isqrt


TEST_CASES = [
    # Semiprimes
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

    # Repeated-prime composites
    (2, 2),      # 4
    (3, 3),      # 9
    (5, 5),      # 25
    (7, 7),      # 49
    (11, 11),    # 121

    # Three distinct prime factors
    (2, 3),      # 6 -- interpreted below as a semiprime
]


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


def sigma(n):
    """
    sigma(n) = sum of positive divisors.
    """
    return sum(divisors(n))


def recover_from_sigma(n, sigma_n):
    """
    Assume n=pq with p,q prime.

    Then:

        S = p+q = sigma(n)-n-1

    and:

        p,q = (S +/- sqrt(S^2-4n))/2
    """

    S = sigma_n - n - 1

    discriminant = S * S - 4 * n

    if discriminant < 0:
        return S, discriminant, None, None, False

    root = isqrt(discriminant)

    if root * root != discriminant:
        return S, discriminant, None, None, False

    if (S + root) % 2 != 0:
        return S, discriminant, None, None, False

    p = (S + root) // 2
    q = (S - root) // 2

    valid = (
        p > 1
        and q > 1
        and p * q == n
    )

    return S, discriminant, p, q, valid


print("START EXPERIMENT 11")
print()

print(
    "n\t"
    "sigma\t"
    "S=sigma-n-1\t"
    "discriminant\t"
    "recovered_p\t"
    "recovered_q\t"
    "valid_factorization"
)

for p_true, q_true in TEST_CASES:

    n = p_true * q_true

    sigma_n = sigma(n)

    S, discriminant, p, q, valid = (
        recover_from_sigma(n, sigma_n)
    )

    print(
        f"{n}\t"
        f"{sigma_n}\t"
        f"{S}\t"
        f"{discriminant}\t"
        f"{p}\t"
        f"{q}\t"
        f"{valid}"
    )


print()
print("CONTROL TESTS")
print("-" * 80)

# Numbers that are not semiprimes.
CONTROL_NUMBERS = [
    8,      # 2^3
    12,     # 2^2*3
    18,     # 2*3^2
    24,     # 2^3*3
    30,     # 2*3*5
    36,     # 2^2*3^2
    60,     # 2^2*3*5
    72,     # 2^3*3^2
    180,    # 2^2*3^2*5
]

print(
    "n\t"
    "sigma\t"
    "S_candidate\t"
    "discriminant\t"
    "candidate_p\t"
    "candidate_q\t"
    "valid"
)

for n in CONTROL_NUMBERS:

    sigma_n = sigma(n)

    S, discriminant, p, q, valid = (
        recover_from_sigma(n, sigma_n)
    )

    print(
        f"{n}\t"
        f"{sigma_n}\t"
        f"{S}\t"
        f"{discriminant}\t"
        f"{p}\t"
        f"{q}\t"
        f"{valid}"
    )


print()
print("IDENTITY CHECK FOR SEMIPRIMES")
print("-" * 80)

for p, q in [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
]:

    n = p * q
    sigma_n = sigma(n)

    expected_sigma = (
        (p + 1) * (q + 1)
    )

    expected_S = p + q

    recovered_S = (
        sigma_n - n - 1
    )

    print(
        f"n={n}: "
        f"sigma={sigma_n}, "
        f"expected_sigma={expected_sigma}, "
        f"S={recovered_S}, "
        f"expected_S={expected_S}, "
        f"OK={sigma_n == expected_sigma and recovered_S == expected_S}"
    )


print()
print("FINISHED EXPERIMENT 11")
