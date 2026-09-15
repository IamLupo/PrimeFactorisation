import math
from sympy import isprime


EXPERIMENT = 86


def binomial_transform_value(sequence, r):
    total = 0

    for j in range(r):
        coefficient = (
            math.comb(r - 1, j)
            * (2 ** (r - 1 - j))
        )

        total += (
            coefficient
            * sequence(r + j)
        )

    return total


def make_periodic_sequence(values):
    period = len(values)

    def sequence(index):
        return values[(index - 1) % period]

    return sequence


def transformed_sequence(sequence, count):
    result = []

    for r in range(1, count + 1):
        result.append(
            binomial_transform_value(
                sequence,
                r
            )
        )

    return result


def characteristic_polynomial(period):
    """
    Characteristic polynomial:

        (lambda - 2)^period - 1

    Coefficients are returned low-degree first.
    """

    polynomial = [1]

    base = [-2, 1]

    for _ in range(period):
        new = [
            0
            for _ in range(
                len(polynomial) + len(base) - 1
            )
        ]

        for i, a in enumerate(polynomial):
            for j, b in enumerate(base):
                new[i + j] += a * b

        polynomial = new

    polynomial[0] -= 1

    return polynomial


def recurrence_failures(values, polynomial):
    """
    If

        c_0 + c_1 lambda + ... + c_k lambda^k = 0

    then

        b_n = -(c_0 b_{n-1} + ... + c_{k-1} b_{n-k})
    """

    degree = len(polynomial) - 1
    failures = 0

    for i in range(degree, len(values)):
        total = 0

        for j in range(degree):
            total += (
                polynomial[j]
                * values[i - j - 1]
            )

        expected = -total

        if values[i] != expected:
            failures += 1

    return failures


def alternating_values(n):
    return [
        1,
        n - 1
    ]


def period_three_values(n):
    return [
        1,
        n - 1,
        1
    ]


def period_four_values(n):
    return [
        1,
        n - 1,
        1,
        n - 1
    ]


def period_five_values(n):
    return [
        1,
        n - 1,
        1,
        n - 1,
        1
    ]


def verify_transform_and_recurrence(
    value_function,
    n,
    count=30
):
    period_values = value_function(n)

    sequence = make_periodic_sequence(
        period_values
    )

    transformed = transformed_sequence(
        sequence,
        count
    )

    polynomial = characteristic_polynomial(
        len(period_values)
    )

    failures = recurrence_failures(
        transformed,
        polynomial
    )

    return (
        period_values,
        transformed,
        polynomial,
        failures
    )


def prime_statistics(
    value_function,
    min_n=2,
    max_n=100
):
    prime_prime = 0
    prime_composite = 0
    composite_prime = 0
    composite_composite = 0

    for n in range(min_n, max_n + 1):
        period_values = value_function(n)

        sequence = make_periodic_sequence(
            period_values
        )

        transformed = transformed_sequence(
            sequence,
            n
        )

        b_n = transformed[n - 1]

        n_prime = isprime(n)
        b_prime = isprime(b_n)

        if n_prime and b_prime:
            prime_prime += 1

        elif n_prime and not b_prime:
            prime_composite += 1

        elif not n_prime and b_prime:
            composite_prime += 1

        else:
            composite_composite += 1

    return (
        prime_prime,
        prime_composite,
        composite_prime,
        composite_composite
    )


def print_period_example(
    name,
    value_function,
    n
):
    (
        period_values,
        transformed,
        polynomial,
        failures
    ) = verify_transform_and_recurrence(
        value_function,
        n,
        count=12
    )

    residues = [
        value % n
        for value in transformed
    ]

    print(name)
    print("-" * 60)
    print(f"N = {n}")
    print(f"PERIOD VALUES = {period_values}")
    print(f"B VALUES = {transformed}")
    print(f"B MOD N = {residues}")
    print(f"CHARACTERISTIC POLYNOMIAL = {polynomial}")
    print(f"RECURRENCE FAILURES = {failures}")
    print()


def print_prime_statistics(
    name,
    value_function
):
    (
        prime_prime,
        prime_composite,
        composite_prime,
        composite_composite
    ) = prime_statistics(
        value_function,
        2,
        100
    )

    print(name)
    print(
        f"  PRIME n, PRIME b_n: "
        f"{prime_prime}"
    )
    print(
        f"  PRIME n, COMPOSITE b_n: "
        f"{prime_composite}"
    )
    print(
        f"  COMPOSITE n, PRIME b_n: "
        f"{composite_prime}"
    )
    print(
        f"  COMPOSITE n, COMPOSITE b_n: "
        f"{composite_composite}"
    )
    print()


def verify_periodic_recurrences():
    tests = [
        (
            "PERIOD 2",
            alternating_values
        ),
        (
            "PERIOD 3",
            period_three_values
        ),
        (
            "PERIOD 4",
            period_four_values
        ),
        (
            "PERIOD 5",
            period_five_values
        )
    ]

    failures = 0

    for name, function in tests:
        (
            _,
            _,
            _,
            recurrence_failure_count
        ) = verify_transform_and_recurrence(
            function,
            17,
            count=40
        )

        print(
            f"{name}: "
            f"{recurrence_failure_count} failures"
        )

        failures += recurrence_failure_count

    return failures


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print("1. PERIODIC RECURRENCE VERIFICATION")
    print("-" * 60)

    total_failures = verify_periodic_recurrences()

    print()
    print(
        f"TOTAL RECURRENCE FAILURES: "
        f"{total_failures}"
    )
    print()

    print("2. PERIODIC EXAMPLES")
    print("-" * 60)

    print_period_example(
        "PERIOD 2",
        alternating_values,
        17
    )

    print_period_example(
        "PERIOD 3",
        period_three_values,
        17
    )

    print_period_example(
        "PERIOD 4",
        period_four_values,
        17
    )

    print_period_example(
        "PERIOD 5",
        period_five_values,
        17
    )

    print("3. PRIME STATISTICS")
    print("-" * 60)

    print_prime_statistics(
        "PERIOD 2",
        alternating_values
    )

    print_prime_statistics(
        "PERIOD 3",
        period_three_values
    )

    print_prime_statistics(
        "PERIOD 4",
        period_four_values
    )

    print_prime_statistics(
        "PERIOD 5",
        period_five_values
    )

    print("4. CHARACTERISTIC POLYNOMIALS")
    print("-" * 60)

    for period in range(2, 9):
        polynomial = characteristic_polynomial(
            period
        )

        print(
            f"PERIOD {period}: "
            f"{polynomial}"
        )

    print()

    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()