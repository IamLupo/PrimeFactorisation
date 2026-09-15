import math
import random

from sympy import isprime


EXPERIMENT = 85


def direct_b(n, r):
    total = 0

    for j in range(r):
        coefficient = (
            math.comb(r - 1, j)
            * (2 ** (r - 1 - j))
        )

        index = r + j

        if index % 2 == 1:
            value = 1
        else:
            value = n - 1

        total += coefficient * value

    return total


def closed_b(n, r):
    sign = (-1) ** r

    numerator = (
        n * (3 ** (r - 1))
        + (n - 2) * sign
    )

    return numerator // 2


def recurrence_b(n, count):
    values = []

    if count >= 1:
        values.append(1)

    if count >= 2:
        values.append(2 * n - 1)

    while len(values) < count:
        values.append(
            2 * values[-1]
            + 3 * values[-2]
        )

    return values


def verify_direct_formula(max_n=100, max_r=20):
    failures = 0

    for n in range(2, max_n + 1):
        for r in range(1, max_r + 1):
            direct = direct_b(n, r)
            closed = closed_b(n, r)

            if direct != closed:
                failures += 1

    return failures


def verify_recurrence(max_n=100, max_r=20):
    failures = 0

    for n in range(2, max_n + 1):
        values = recurrence_b(n, max_r)

        for r in range(1, max_r - 1):
            expected = (
                2 * values[r]
                + 3 * values[r - 1]
            )

            if values[r + 1] != expected:
                failures += 1

    return failures


def verify_congruence(max_n=100, max_r=50):
    failures = 0

    for n in range(2, max_n + 1):
        for r in range(1, max_r + 1):
            value = closed_b(n, r)
            expected = (-1) ** (r + 1)

            if (value - expected) % n != 0:
                failures += 1

    return failures


def verify_quotient_identity(max_n=100, max_r=50):
    failures = 0

    for n in range(2, max_n + 1):
        for r in range(1, max_r + 1):
            value = closed_b(n, r)

            left = (
                value
                - ((-1) ** (r + 1))
            )

            right = (
                n
                * (
                    (3 ** (r - 1))
                    - ((-1) ** r)
                )
                // 2
            )

            if left != right:
                failures += 1

    return failures


def print_example(n, count=8):
    values = recurrence_b(n, count)

    print(f"N = {n}")
    print(f"VALUES = {values}")
    print()

    for r, value in enumerate(values, start=1):
        print(
            f"r={r:2d} "
            f"b_r={value} "
            f"prime={isprime(value)} "
            f"mod_N={value % n}"
        )

    print()


def prime_pattern_scan(max_n=100):
    prime_n_prime_bn = 0
    prime_n_composite_bn = 0

    composite_n_prime_bn = 0
    composite_n_composite_bn = 0

    for n in range(2, max_n + 1):
        value = closed_b(n, n)

        n_prime = isprime(n)
        b_prime = isprime(value)

        if n_prime and b_prime:
            prime_n_prime_bn += 1

        elif n_prime and not b_prime:
            prime_n_composite_bn += 1

        elif not n_prime and b_prime:
            composite_n_prime_bn += 1

        else:
            composite_n_composite_bn += 1

    return (
        prime_n_prime_bn,
        prime_n_composite_bn,
        composite_n_prime_bn,
        composite_n_composite_bn,
    )


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print("1. DIRECT TRANSFORM VS CLOSED FORM")
    print("-" * 60)

    failures = verify_direct_formula(
        max_n=100,
        max_r=20
    )

    print(
        f"IDENTITY FAILURES: {failures}"
    )
    print()

    print("2. RECURRENCE")
    print("-" * 60)

    failures = verify_recurrence(
        max_n=100,
        max_r=30
    )

    print(
        f"RECURRENCE FAILURES: {failures}"
    )
    print()

    print("3. MODULO n")
    print("-" * 60)

    failures = verify_congruence(
        max_n=100,
        max_r=50
    )

    print(
        f"CONGRUENCE FAILURES: {failures}"
    )

    print(
        "EXPECTED:"
    )
    print(
        "b_r == (-1)^(r+1) (mod n)"
    )
    print()

    print("4. QUOTIENT IDENTITY")
    print("-" * 60)

    failures = verify_quotient_identity(
        max_n=100,
        max_r=50
    )

    print(
        f"QUOTIENT IDENTITY FAILURES: {failures}"
    )
    print()

    print("5. EXAMPLES")
    print("-" * 60)

    print_example(
        n=7,
        count=8
    )

    print_example(
        n=11,
        count=8
    )

    print_example(
        n=15,
        count=8
    )

    print("6. PRIME PATTERN FOR b_n")
    print("-" * 60)

    (
        prime_n_prime_bn,
        prime_n_composite_bn,
        composite_n_prime_bn,
        composite_n_composite_bn,
    ) = prime_pattern_scan(
        max_n=100
    )

    print(
        "n prime, b_n prime: "
        f"{prime_n_prime_bn}"
    )

    print(
        "n prime, b_n composite: "
        f"{prime_n_composite_bn}"
    )

    print(
        "n composite, b_n prime: "
        f"{composite_n_prime_bn}"
    )

    print(
        "n composite, b_n composite: "
        f"{composite_n_composite_bn}"
    )

    print()

    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    random.seed(850011)
    run_experiment()


if __name__ == "__main__":
    main()
