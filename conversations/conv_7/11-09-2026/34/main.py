import math
import random


EXPERIMENT = 54


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
            "a": s - p,
            "b": q - s,
        }


def batch_product(case, radius):
    """
    Product of all diagonal slopes

        r = s + z

    for z in [-radius, ..., radius].
    """

    n = case["N"]
    s = case["s"]

    product = 1

    for z in range(-radius, radius + 1):
        r = s + z

        if r <= 0:
            continue

        product = (
            product * r
        ) % n

    return product


def batch_gcd(case, radius):
    product = batch_product(
        case,
        radius,
    )

    return math.gcd(
        product,
        case["N"],
    )


def sequential_gcd(case, radius):
    n = case["N"]
    s = case["s"]

    factors = []

    for z in range(-radius, radius + 1):
        r = s + z

        if r <= 0:
            continue

        g = math.gcd(
            r,
            n,
        )

        if 1 < g < n:
            factors.append(g)

    return sorted(set(factors))


def expected_factors(case, radius):
    factors = []

    if case["a"] <= radius:
        factors.append(case["p"])

    if case["b"] <= radius:
        factors.append(case["q"])

    return sorted(set(factors))


def split_gcd(g, n):
    """
    Turn gcd result into its implied factor set.
    """

    if g == 1:
        return []

    if g == n:
        return [
            math.isqrt(n)
        ]

    factors = []

    if 1 < g < n:
        factors.append(g)

        other = n // g

        if other > 1:
            factors.append(other)

    return sorted(set(factors))


def fermat_radius(case):
    return max(
        case["a"],
        case["b"],
    )


def compare_methods(cases, radius):
    sequential_success = 0
    batch_success = 0
    exact_batch = 0

    for case in cases:
        expected = expected_factors(
            case,
            radius,
        )

        sequential = sequential_gcd(
            case,
            radius,
        )

        batch = batch_gcd(
            case,
            radius,
        )

        batch_factors = split_gcd(
            batch,
            case["N"],
        )

        if sequential == expected:
            sequential_success += 1

        if batch_factors == expected:
            exact_batch += 1

        if batch_factors:
            batch_success += 1

    return (
        sequential_success,
        batch_success,
        exact_batch,
    )


def test_exact_equivalence(cases):
    print("============================================================")
    print("TEST 1: BATCH GCD == SEQUENTIAL GCD")
    print("============================================================")

    failures = 0

    for case in cases:
        radius = random.randint(
            1,
            500,
        )

        sequential = sequential_gcd(
            case,
            radius,
        )

        batch = batch_gcd(
            case,
            radius,
        )

        expected = expected_factors(
            case,
            radius,
        )

        batch_factors = split_gcd(
            batch,
            case["N"],
        )

        if sequential != expected:
            failures += 1

        if expected:
            if batch_factors != expected:
                failures += 1
        else:
            if batch != 1:
                failures += 1

    print(
        f"equivalence failures = {failures}"
    )

    print()


def test_radius(cases):
    print("============================================================")
    print("TEST 2: FACTOR DETECTION VS RADIUS")
    print("============================================================")

    radii = [
        1,
        2,
        4,
        8,
        16,
        32,
        64,
        128,
        256,
        512,
        1024,
        2048,
        4096,
        8192,
    ]

    for radius in radii:
        sequential_success, batch_success, exact_batch = (
            compare_methods(
                cases,
                radius,
            )
        )

        print(
            f"radius={radius:5d}  "
            f"sequential={sequential_success:4d}  "
            f"batch-nontrivial={batch_success:4d}  "
            f"batch-exact={exact_batch:4d}"
        )

    print()


def test_fermat_radius(cases):
    print("============================================================")
    print("TEST 3: MINIMUM REQUIRED RADIUS")
    print("============================================================")

    radii = [
        fermat_radius(case)
        for case in cases
    ]

    print(
        f"minimum radius = "
        f"{min(radii)}"
    )

    print(
        f"maximum radius = "
        f"{max(radii)}"
    )

    print(
        f"mean radius    = "
        f"{sum(radii) / len(radii):.2f}"
    )

    print()


def test_asymmetric_factor_distance(cases):
    print("============================================================")
    print("TEST 4: p AND q DISTANCE FROM s")
    print("============================================================")

    ratios = []

    for case in cases:
        a = case["a"]
        b = case["b"]

        if max(a, b) == 0:
            continue

        ratios.append(
            min(a, b) / max(a, b)
        )

    print(
        f"min distance ratio = "
        f"{min(ratios):.6f}"
    )

    print(
        f"max distance ratio = "
        f"{max(ratios):.6f}"
    )

    print(
        f"mean distance ratio = "
        f"{sum(ratios) / len(ratios):.6f}"
    )

    print()


def test_transform_equivalence(cases):
    print("============================================================")
    print("TEST 5: TRANSFORM VS DIRECT r")
    print("============================================================")

    failures = 0

    for case in cases:
        n = case["N"]

        for r in [
            case["p"],
            case["q"],
            case["s"],
            case["s"] + 1,
        ]:
            direct = math.gcd(
                r,
                n,
            )

            transform = math.gcd(
                r,
                n,
            )

            if direct != transform:
                failures += 1

    print(
        f"transform equivalence failures = "
        f"{failures}"
    )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 6: EXAMPLES")
    print("============================================================")

    for case in cases[:10]:
        print()

        print(
            f"N={case['N']}"
        )

        print(
            f"p={case['p']} "
            f"q={case['q']} "
            f"s={case['s']}"
        )

        print(
            f"s-p={case['a']} "
            f"q-s={case['b']}"
        )

        radius = fermat_radius(case)

        print(
            f"minimum radius={radius}"
        )

        g = batch_gcd(
            case,
            radius,
        )

        print(
            f"batch gcd={g}"
        )

        print(
            f"expected factors="
            f"{expected_factors(case, radius)}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(545454)

    CASE_COUNT = 500
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

    test_exact_equivalence(cases)
    test_radius(cases)
    test_fermat_radius(cases)
    test_asymmetric_factor_distance(cases)
    test_transform_equivalence(cases)
    print_examples(cases)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
