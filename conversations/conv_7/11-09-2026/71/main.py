import math
import random

from sympy import factorint, isprime


EXPERIMENT = 89


def binomial_transform_value(n, r):
    total = 0

    for j in range(r):
        coefficient = (
            math.comb(r - 1, j)
            * (2 ** (r - 1 - j))
        )

        if (r + j) % 2 == 1:
            value = 1
        else:
            value = n - 1

        total += coefficient * value

    return total


def closed_form(n, r):
    return (
        n * (3 ** (r - 1))
        + (n - 2) * ((-1) ** r)
    ) // 2


def verify_closed_form(max_n=100, max_r=20):
    failures = 0

    for n in range(2, max_n + 1):
        for r in range(1, max_r + 1):
            direct = binomial_transform_value(
                n,
                r
            )

            formula = closed_form(
                n,
                r
            )

            if direct != formula:
                failures += 1

    return failures


def factorization_test(n, r):
    b = closed_form(n, r)

    g_minus = math.gcd(
        b - 1,
        n
    )

    g_plus = math.gcd(
        b + 1,
        n
    )

    sign = (-1) ** (r + 1)

    numerator = (
        2 * (b - sign)
    )

    denominator = (
        3 ** (r - 1)
        - ((-1) ** r)
    )

    if denominator == 0:
        quotient_gcd = None
    else:
        if numerator % denominator != 0:
            quotient_gcd = None
        else:
            quotient = numerator // denominator

            quotient_gcd = math.gcd(
                quotient,
                n
            )

    g_three_minus = math.gcd(
        3 ** r - 1,
        n
    )

    g_three_plus = math.gcd(
        3 ** r + 1,
        n
    )

    return {
        "b": b,
        "g_minus": g_minus,
        "g_plus": g_plus,
        "quotient_gcd": quotient_gcd,
        "g_three_minus": g_three_minus,
        "g_three_plus": g_three_plus,
    }


def orders_of_three(n):
    if math.gcd(3, n) != 1:
        return None

    factors = factorint(n)

    orders = {}

    for p in factors:
        if p == 3:
            orders[p] = None
            continue

        value = 1

        for r in range(1, p):
            value = (
                value * 3
            ) % p

            if value == 1:
                orders[p] = r
                break

    return orders


def run_example(n, max_r=12):
    print(f"N = {n}")
    print(
        f"PRIME = {isprime(n)}"
    )

    orders = orders_of_three(n)

    print(
        f"ORDERS OF 3: {orders}"
    )

    print()

    for r in range(1, max_r + 1):
        result = factorization_test(
            n,
            r
        )

        print(
            f"r={r:2d} "
            f"b={result['b']} "
            f"gcd(b-1,n)={result['g_minus']} "
            f"gcd(b+1,n)={result['g_plus']} "
            f"quotient_gcd={result['quotient_gcd']} "
            f"gcd(3^r-1,n)={result['g_three_minus']} "
            f"gcd(3^r+1,n)={result['g_three_plus']}"
        )

    print()


def random_semiprime():
    while True:
        p = random.randint(
            100,
            10000
        )

        q = random.randint(
            100,
            10000
        )

        if p == q:
            continue

        if not isprime(p):
            continue

        if not isprime(q):
            continue

        if p > q:
            p, q = q, p

        return p * q, p, q


def test_semiprimes(
    count=500,
    max_r=32
):
    nontrivial_b_minus = 0
    nontrivial_b_plus = 0
    nontrivial_quotient = 0
    nontrivial_three_minus = 0
    nontrivial_three_plus = 0

    for _ in range(count):
        n, p, q = random_semiprime()

        hit_minus = False
        hit_plus = False
        hit_quotient = False
        hit_three_minus = False
        hit_three_plus = False

        for r in range(
            1,
            max_r + 1
        ):
            result = factorization_test(
                n,
                r
            )

            for value in [
                result["g_minus"],
                result["g_plus"],
                result["quotient_gcd"],
                result["g_three_minus"],
                result["g_three_plus"],
            ]:
                if value in (None, 1, n):
                    continue

                if value == p or value == q:
                    if value == result["g_minus"]:
                        hit_minus = True

                    if value == result["g_plus"]:
                        hit_plus = True

                    if value == result["quotient_gcd"]:
                        hit_quotient = True

                    if value == result["g_three_minus"]:
                        hit_three_minus = True

                    if value == result["g_three_plus"]:
                        hit_three_plus = True

        if hit_minus:
            nontrivial_b_minus += 1

        if hit_plus:
            nontrivial_b_plus += 1

        if hit_quotient:
            nontrivial_quotient += 1

        if hit_three_minus:
            nontrivial_three_minus += 1

        if hit_three_plus:
            nontrivial_three_plus += 1

    return (
        nontrivial_b_minus,
        nontrivial_b_plus,
        nontrivial_quotient,
        nontrivial_three_minus,
        nontrivial_three_plus,
    )


def run_experiment():
    print("=" * 60)
    print(
        f"START EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)
    print()

    print("1. CLOSED FORM VERIFICATION")
    print("-" * 60)

    failures = verify_closed_form()

    print(
        f"FAILURES: {failures}"
    )
    print()

    print("2. PRIME EXAMPLE")
    print("-" * 60)

    run_example(
        17,
        max_r=12
    )

    print("3. COMPOSITE EXAMPLES")
    print("-" * 60)

    run_example(
        15,
        max_r=12
    )

    run_example(
        77,
        max_r=12
    )

    run_example(
        143,
        max_r=12
    )

    print("4. SEMIPRIME GCD TEST")
    print("-" * 60)

    (
        b_minus,
        b_plus,
        quotient,
        three_minus,
        three_plus
    ) = test_semiprimes(
        count=500,
        max_r=32
    )

    print(
        "CASES WITH NONTRIVIAL gcd(b_r-1,N): "
        f"{b_minus}/500"
    )

    print(
        "CASES WITH NONTRIVIAL gcd(b_r+1,N): "
        f"{b_plus}/500"
    )

    print(
        "CASES WITH NONTRIVIAL QUOTIENT GCD: "
        f"{quotient}/500"
    )

    print(
        "CASES WITH NONTRIVIAL gcd(3^r-1,N): "
        f"{three_minus}/500"
    )

    print(
        "CASES WITH NONTRIVIAL gcd(3^r+1,N): "
        f"{three_plus}/500"
    )

    print()

    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    random.seed(890011)

    run_experiment()


if __name__ == "__main__":
    main()
