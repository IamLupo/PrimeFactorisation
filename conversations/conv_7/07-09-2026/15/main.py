from math import isqrt


# ============================================================
# NUMBER THEORY
# ============================================================

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
    return sum(divisors(n))


def sigma3(n):
    return sum(
        d ** 3
        for d in divisors(n)
    )


def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    r = isqrt(n)

    for d in range(3, r + 1, 2):
        if n % d == 0:
            return False

    return True


def factor_count_distinct(n):
    """
    Return the number of distinct prime factors.

    Used only for classification of the test data.
    """
    count = 0
    x = n

    if x % 2 == 0:
        count += 1
        while x % 2 == 0:
            x //= 2

    p = 3

    while p * p <= x:

        if x % p == 0:
            count += 1

            while x % p == 0:
                x //= p

        p += 2

    if x > 1:
        count += 1

    return count


def classify(n):
    """
    Classification for analysis only.
    """
    if is_prime(n):
        return "prime"

    count = factor_count_distinct(n)

    # Detect square semiprimes such as p^2.
    r = isqrt(n)

    if r * r == n and is_prime(r):
        return "prime_square"

    if count == 2:
        # Check whether both prime factors have exponent 1.
        # Equivalently: n is squarefree when the two
        # distinct prime factors occur once each.
        for p in range(2, r + 1):
            if n % p == 0:
                q = n // p

                if is_prime(p) and is_prime(q) and p != q:
                    return "semiprime"

                break

    if count >= 3:
        return "3plus_prime_factors"

    return "prime_power_or_other"


# ============================================================
# SEMIPRIME RECOVERY
# ============================================================

def recover_from_sigma(n, sigma_n):
    """
    Assume n=pq with p,q prime.

        S = p+q = sigma(n)-n-1

        Delta = S^2 - 4n

    """
    S = sigma_n - n - 1

    delta = S * S - 4 * n

    if delta < 0:
        return S, delta, None, None, False

    root = isqrt(delta)

    if root * root != delta:
        return S, delta, None, None, False

    if (S + root) % 2 != 0:
        return S, delta, None, None, False

    p = (S + root) // 2
    q = (S - root) // 2

    valid = (
        p > 1
        and q > 1
        and p * q == n
        and is_prime(p)
        and is_prime(q)
    )

    return S, delta, p, q, valid


# ============================================================
# EXPERIMENT 12
# ============================================================

print("START EXPERIMENT 12")
print()

MAX_N = 1000

print("PART 1: SEMIPRIME IDENTITY CHECK")
print("-" * 120)

semiprime_total = 0
semiprime_failures = 0

print(
    "n\t"
    "sigma\t"
    "sigma3\t"
    "S\t"
    "predicted_sigma3\t"
    "D\t"
    "recovered_p\t"
    "recovered_q\t"
    "OK"
)

for n in range(4, MAX_N + 1):

    if classify(n) != "semiprime":
        continue

    semiprime_total += 1

    sigma_n = sigma(n)
    sigma3_n = sigma3(n)

    S = sigma_n - n - 1

    predicted_sigma3 = (
        n**3
        + 1
        + S**3
        - 3 * n * S
    )

    D = sigma3_n - predicted_sigma3

    recovered_S, delta, p, q, valid = (
        recover_from_sigma(n, sigma_n)
    )

    ok = (
        D == 0
        and recovered_S == S
        and valid
    )

    if not ok:
        semiprime_failures += 1

    print(
        f"{n}\t"
        f"{sigma_n}\t"
        f"{sigma3_n}\t"
        f"{S}\t"
        f"{predicted_sigma3}\t"
        f"{D}\t"
        f"{p}\t"
        f"{q}\t"
        f"{ok}"
    )


print()
print(
    f"Semiprimes tested: {semiprime_total}"
)
print(
    f"Semiprime failures: {semiprime_failures}"
)

print()
print("PART 2: SEARCH FOR NON-SEMIPRIME FALSE POSITIVES")
print("-" * 120)

false_positives = []

print(
    "n\t"
    "classification\t"
    "sigma\t"
    "sigma3\t"
    "S_candidate\t"
    "D\t"
    "quadratic_p\t"
    "quadratic_q\t"
    "valid_prime_factorization"
)

for n in range(4, MAX_N + 1):

    classification = classify(n)

    if classification == "semiprime":
        continue

    sigma_n = sigma(n)
    sigma3_n = sigma3(n)

    S = sigma_n - n - 1

    predicted_sigma3 = (
        n**3
        + 1
        + S**3
        - 3 * n * S
    )

    D = sigma3_n - predicted_sigma3

    recovered_S, delta, p, q, valid = (
        recover_from_sigma(n, sigma_n)
    )

    #
    # A false positive would be a non-semiprime
    # satisfying the same sigma/sigma3 identity
    # AND producing two prime factors.
    #
    if D == 0 and valid:
        false_positives.append(n)

        print(
            f"{n}\t"
            f"{classification}\t"
            f"{sigma_n}\t"
            f"{sigma3_n}\t"
            f"{S}\t"
            f"{D}\t"
            f"{p}\t"
            f"{q}\t"
            f"{valid}"
        )


print()

print(
    f"Non-semiprime false positives found: "
    f"{len(false_positives)}"
)

if false_positives:
    print(
        "FALSE POSITIVES:",
        false_positives
    )
else:
    print(
        "No non-semiprime false positives found."
    )


print()
print("PART 3: SAMPLE RESIDUALS")
print("-" * 100)

sample_numbers = [
    6, 8, 9, 10, 12,
    15, 16, 18, 20, 24,
    30, 36, 42, 60,
    72, 120, 180,
]

print(
    "n\t"
    "classification\t"
    "sigma\t"
    "sigma3\t"
    "S_candidate\t"
    "D"
)

for n in sample_numbers:

    sigma_n = sigma(n)
    sigma3_n = sigma3(n)

    S = sigma_n - n - 1

    predicted_sigma3 = (
        n**3
        + 1
        + S**3
        - 3 * n * S
    )

    D = sigma3_n - predicted_sigma3

    print(
        f"{n}\t"
        f"{classify(n)}\t"
        f"{sigma_n}\t"
        f"{sigma3_n}\t"
        f"{S}\t"
        f"{D}"
    )


print()
print("FINISHED EXPERIMENT 12")

