import math
from sympy import (
    symbols,
    expand,
    cyclotomic_poly,
    resultant,
    isprime,
)


EXPERIMENT = 87


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


def polynomial_coefficients(expr, variable):
    poly = expand(expr)

    degree = poly.as_poly(variable).degree()

    coefficients = []

    for i in range(degree + 1):
        coefficients.append(
            int(
                poly.as_poly(variable).coeff_monomial(
                    variable ** i
                )
            )
        )

    return coefficients


def spectral_polynomial(period):
    """
    The eigenvalues are

        lambda = omega * (2 + omega)

    where

        omega^period = 1.

    We obtain the polynomial by eliminating omega from

        Phi_k(omega) = 0
        lambda - omega*(2 + omega) = 0.
    """

    omega = symbols("omega")
    lam = symbols("lam")

    cyclotomic = cyclotomic_poly(
        period,
        omega
    )

    relation = (
        lam
        - omega * (2 + omega)
    )

    resultant_poly = expand(
        resultant(
            cyclotomic,
            relation,
            omega
        )
    )

    result = resultant_poly.as_poly(
        lam
    )

    return [
        int(
            result.coeff_monomial(
                lam ** i
            )
        )
        for i in range(
            result.degree() + 1
        )
    ]


def recurrence_failures(
    values,
    coefficients
):
    degree = len(coefficients) - 1

    failures = 0

    for i in range(
        degree,
        len(values)
    ):
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


def run_period(
    name,
    value_function,
    n,
    count=30
):
    period_values = value_function(n)

    sequence = make_periodic_sequence(
        period_values
    )

    values = transformed_sequence(
        sequence,
        count
    )

    coefficients = spectral_polynomial(
        len(period_values)
    )

    failures = recurrence_failures(
        values,
        coefficients
    )

    residues = [
        value % n
        for value in values[:20]
    ]

    print(name)
    print("-" * 60)
    print(
        f"PERIOD VALUES: "
        f"{period_values}"
    )
    print(
        f"VALUES: "
        f"{values[:12]}"
    )
    print(
        f"MOD {n}: "
        f"{residues}"
    )
    print(
        f"SPECTRAL POLYNOMIAL: "
        f"{coefficients}"
    )
    print(
        f"RECURRENCE FAILURES: "
        f"{failures}"
    )
    print()


def prime_statistics(
    value_function,
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
        values = value_function(n)

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


def print_prime_statistics(
    name,
    value_function
):
    (
        pp,
        pc,
        cp,
        cc
    ) = prime_statistics(
        value_function
    )

    print(name)
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


def verify_modular_pattern(
    value_function,
    n,
    count=30
):
    values = transformed_sequence(
        make_periodic_sequence(
            value_function(n)
        ),
        count
    )

    residues = [
        value % n
        for value in values
    ]

    return residues


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print("1. EXACT SPECTRAL RECURRENCES")
    print("-" * 60)

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

    for name, function in tests:
        run_period(
            name,
            function,
            17,
            count=40
        )

    print(
        "2. SPECTRAL POLYNOMIALS"
    )
    print("-" * 60)

    for period in range(2, 9):
        coefficients = spectral_polynomial(
            period
        )

        print(
            f"PERIOD {period}: "
            f"{coefficients}"
        )

    print()

    print("3. MODULAR PATTERNS")
    print("-" * 60)

    for name, function in tests:
        residues = verify_modular_pattern(
            function,
            17,
            count=20
        )

        print(
            f"{name}: {residues}"
        )

    print()

    print("4. PRIME STATISTICS")
    print("-" * 60)

    for name, function in tests:
        print_prime_statistics(
            name,
            function
        )

    print(
        "5. PERIOD-2 DIRECT CHECK"
    )
    print("-" * 60)

    values = transformed_sequence(
        make_periodic_sequence(
            alternating_values(17)
        ),
        20
    )

    failures = 0

    for i in range(
        len(values) - 2
    ):
        expected = (
            2 * values[i + 1]
            + 3 * values[i]
        )

        if values[i + 2] != expected:
            failures += 1

    print(
        f"PERIOD-2 RECURRENCE:"
    )
    print(
        f"b_(r+2) = 2*b_(r+1) + 3*b_r"
    )
    print(
        f"FAILURES: {failures}"
    )
    print()

    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
