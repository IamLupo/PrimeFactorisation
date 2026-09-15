import math
from sympy import symbols, expand, resultant, isprime


EXPERIMENT = 88


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


def full_spectral_polynomial(period):
    omega = symbols("omega")
    lam = symbols("lam")

    relation_root = (
        lam
        - omega * (2 + omega)
    )

    root_of_unity = (
        omega ** period - 1
    )

    result = expand(
        resultant(
            root_of_unity,
            relation_root,
            omega
        )
    )

    polynomial = result.as_poly(lam)

    degree = polynomial.degree()

    coefficients = []

    for i in range(degree + 1):
        coefficients.append(
            int(
                polynomial.coeff_monomial(
                    lam ** i
                )
            )
        )

    return coefficients


def recurrence_failures(values, coefficients):
    degree = len(coefficients) - 1

    failures = 0

    for i in range(degree, len(values)):
        total = 0

        for j in range(degree):
            total += (
                coefficients[j]
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


def generic_period_values(n, period):
    values = []

    for i in range(period):
        if i % 2 == 0:
            values.append(1)
        else:
            values.append(n - 1)

    return values


def verify_period(
    n,
    period,
    count
):
    values = generic_period_values(
        n,
        period
    )

    sequence = make_periodic_sequence(
        values
    )

    transformed = transformed_sequence(
        sequence,
        count
    )

    polynomial = full_spectral_polynomial(
        period
    )

    failures = recurrence_failures(
        transformed,
        polynomial
    )

    return (
        values,
        transformed,
        polynomial,
        failures
    )


def print_period_result(
    n,
    period,
    count=50
):
    (
        values,
        transformed,
        polynomial,
        failures
    ) = verify_period(
        n,
        period,
        count
    )

    print(
        f"PERIOD {period}"
    )
    print("-" * 60)

    print(
        f"PERIOD VALUES: {values}"
    )

    print(
        f"B VALUES: {transformed[:12]}"
    )

    print(
        f"B MOD {n}: "
        f"{[x % n for x in transformed[:20]]}"
    )

    print(
        f"SPECTRAL POLYNOMIAL: "
        f"{polynomial}"
    )

    print(
        f"RECURRENCE FAILURES: "
        f"{failures}"
    )

    print()


def prime_statistics(
    period,
    min_n=2,
    max_n=100
):
    prime_prime = 0
    prime_composite = 0
    composite_prime = 0
    composite_composite = 0

    for n in range(
        min_n,
        max_n + 1
    ):
        values = generic_period_values(
            n,
            period
        )

        sequence = make_periodic_sequence(
            values
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


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print("1. FULL ROOT-OF-UNITY SPECTRUM")
    print("-" * 60)

    for period in range(2, 9):
        polynomial = full_spectral_polynomial(
            period
        )

        print(
            f"PERIOD {period}: "
            f"{polynomial}"
        )

    print()

    print("2. EXACT RECURRENCE VERIFICATION")
    print("-" * 60)

    total_failures = 0

    for period in range(2, 9):
        (
            _,
            _,
            _,
            failures
        ) = verify_period(
            17,
            period,
            60
        )

        print(
            f"PERIOD {period}: "
            f"{failures} failures"
        )

        total_failures += failures

    print()

    print(
        f"TOTAL RECURRENCE FAILURES: "
        f"{total_failures}"
    )
    print()

    print("3. EXPLICIT PERIODS")
    print("-" * 60)

    print_period_result(
        17,
        2
    )

    print_period_result(
        17,
        3
    )

    print_period_result(
        17,
        4
    )

    print_period_result(
        17,
        5
    )

    print("4. PRIME STATISTICS")
    print("-" * 60)

    for period in range(2, 9):
        (
            pp,
            pc,
            cp,
            cc
        ) = prime_statistics(
            period
        )

        print(
            f"PERIOD {period}"
        )
        print(
            f"  PRIME n, PRIME b_n: {pp}"
        )
        print(
            f"  PRIME n, COMPOSITE b_n: {pc}"
        )
        print(
            f"  COMPOSITE n, PRIME b_n: {cp}"
        )
        print(
            f"  COMPOSITE n, COMPOSITE b_n: {cc}"
        )
        print()

    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
