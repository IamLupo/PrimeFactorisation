import math
import random


EXPERIMENT = 53


def is_prime(n):
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    ]

    for p in small_primes:
        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    r = 0

    while d % 2 == 0:
        d //= 2
        r += 1

    for a in [2, 3, 5, 7, 11, 13, 17]:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(r - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(bits):
    while True:
        n = random.getrandbits(bits)

        n |= 1
        n |= 1 << (bits - 1)

        if is_prime(n):
            return n


def generate_case(bits):
    while True:
        p = random_prime(bits)
        q = random_prime(bits)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q
        s = math.isqrt(n)

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": n - s * s,
        }


def diagonal_value(case, r, x):
    n = case["N"]
    s = case["s"]

    return (
        n
        - r * (s + 1 - x)
    )


def B_direct(case, r, n):
    total = 0

    for x in range(n + 1):
        total += (
            diagonal_value(
                case,
                r,
                x,
            )
            * math.comb(n, x)
            * (2 ** (n - x))
        )

    return total


def C_closed(case, r, n):
    """
    C_n(r) =
        [B_n(r) - 3^n N] / 3^(n-1)

    =
        r [n - 3(s+1)].
    """

    s = case["s"]

    return (
        r
        * (
            n
            - 3 * (s + 1)
        )
    )


def transform_difference(case, r, n):
    """
    C_(n+1)(r) - C_n(r) = r.
    """

    return (
        C_closed(
            case,
            r,
            n + 1,
        )
        - C_closed(
            case,
            r,
            n,
        )
    )


def test_difference_identity(cases):
    print("============================================================")
    print("TEST 1: CONSECUTIVE TRANSFORM DIFFERENCE")
    print("============================================================")

    failures = 0

    for case in cases:
        # Arbitrary r values, not just divisors.
        values = [
            1,
            2,
            case["s"],
            case["p"],
            case["q"],
            case["N"],
        ]

        for r in values:
            for n in range(1, 10):
                difference = transform_difference(
                    case,
                    r,
                    n,
                )

                if difference != r:
                    failures += 1

    print(
        f"identity failures = {failures}"
    )

    print()


def test_factor_signal(cases):
    print("============================================================")
    print("TEST 2: gcd OF TRANSFORM DIFFERENCE")
    print("============================================================")

    p_hits = 0
    q_hits = 0
    total_tests = 0

    for case in cases:
        for r in range(
            max(1, case["s"] - 100),
            case["s"] + 101,
        ):
            g = math.gcd(
                transform_difference(
                    case,
                    r,
                    5,
                ),
                case["N"],
            )

            expected = math.gcd(
                r,
                case["N"],
            )

            if g != expected:
                print(
                    "ERROR:",
                    case["N"],
                    r,
                    g,
                    expected,
                )
                return

            if g == case["p"]:
                p_hits += 1

            if g == case["q"]:
                q_hits += 1

            total_tests += 1

    print(
        f"tests = {total_tests}"
    )

    print(
        f"p recovered = {p_hits}"
    )

    print(
        f"q recovered = {q_hits}"
    )

    print()


def scan_z_direct(case, radius):
    """
    Directly test z values around 0.

    Since r = z+s, this is equivalent to testing
    candidate diagonal slopes near s.
    """

    hits = []

    for z in range(
        -radius,
        radius + 1,
    ):
        r = z + case["s"]

        if r <= 1:
            continue

        g = math.gcd(
            transform_difference(
                case,
                r,
                5,
            ),
            case["N"],
        )

        if 1 < g < case["N"]:
            hits.append(
                (
                    z,
                    r,
                    g,
                )
            )

    return hits


def test_search_equivalence(cases):
    print("============================================================")
    print("TEST 3: DIAGONAL SEARCH VS TRIAL-DIVISOR SEARCH")
    print("============================================================")

    for radius in [10, 100, 1000, 5000]:
        transform_hits = 0
        direct_hits = 0

        for case in cases:
            hits = scan_z_direct(
                case,
                radius,
            )

            if hits:
                transform_hits += 1

            # Directly test the same candidate slopes.
            found = False

            for z in range(
                -radius,
                radius + 1,
            ):
                r = z + case["s"]

                if r <= 1:
                    continue

                if (
                    1 < math.gcd(
                        r,
                        case["N"],
                    ) < case["N"]
                ):
                    found = True
                    break

            if found:
                direct_hits += 1

        print(
            f"radius={radius:5d}  "
            f"transform={transform_hits:4d}  "
            f"direct={direct_hits:4d}"
        )

    print()


def aggregate_product(case, n, radius):
    """
    Product of all transform differences over a z range:

        ∏ gcd-sensitive quantities

    We reduce modulo N to keep integers manageable.
    """

    N = case["N"]
    s = case["s"]

    value = 1

    for z in range(
        -radius,
        radius + 1,
    ):
        r = z + s

        value = (
            value
            * (
                transform_difference(
                    case,
                    r,
                    n,
                )
                % N
            )
        ) % N

    return value


def aggregate_gcd(case, n, radius):
    value = aggregate_product(
        case,
        n,
        radius,
    )

    return math.gcd(
        value,
        case["N"],
    )


def test_aggregate_products(cases):
    print("============================================================")
    print("TEST 4: AGGREGATED PRODUCT")
    print("============================================================")

    for radius in [1, 2, 4, 8, 16, 32]:
        factor_hits = 0
        nontrivial = 0

        for case in cases:
            g = aggregate_gcd(
                case,
                5,
                radius,
            )

            if 1 < g < case["N"]:
                nontrivial += 1

            if g in [
                case["p"],
                case["q"],
            ]:
                factor_hits += 1

        print(
            f"radius={radius:2d}  "
            f"nontrivial={nontrivial:4d}  "
            f"exact factor={factor_hits:4d}"
        )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 5: EXAMPLES")
    print("============================================================")

    for case in cases[:5]:
        print()

        print(
            f"N={case['N']} "
            f"p={case['p']} "
            f"q={case['q']} "
            f"s={case['s']}"
        )

        for r in [
            case["p"],
            case["q"],
        ]:
            print()
            print(f"r={r}")

            for n in [1, 2, 5, 10]:
                c = transform_difference(
                    case,
                    r,
                    n,
                )

                g = math.gcd(
                    c,
                    case["N"],
                )

                print(
                    f"n={n:2d}  "
                    f"C_(n+1)-C_n={c}  "
                    f"gcd={g}"
                )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(535353)

    CASE_COUNT = 200
    PRIME_BITS = 16

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = [
        generate_case(PRIME_BITS)
        for _ in range(CASE_COUNT)
    ]

    test_difference_identity(cases)
    test_factor_signal(cases)
    test_search_equivalence(cases)
    test_aggregate_products(cases)
    print_examples(cases)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
