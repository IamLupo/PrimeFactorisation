import math
import random


EXPERIMENT = 49


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

        y = p + q - 2 * s - 1
        gap = q - p

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": n - s * s,
            "y": y,
            "gap": gap,
        }


def discriminant(case, y):
    s = case["s"]
    n = case["N"]

    z = 2 * s + 1 + y

    return z * z - 4 * n


def is_square(n):
    if n < 0:
        return False

    r = math.isqrt(n)

    return r * r == n


def fermat_start(n):
    return math.isqrt(n - 1) + 1


def fermat_iterations(n, p, q):
    """
    Standard Fermat search:

        A^2 - N = B^2

    with A starting at ceil(sqrt(N)).
    """

    a = fermat_start(n)

    iterations = 0

    while True:
        value = a * a - n

        if is_square(value):
            b = math.isqrt(value)

            if (
                a - b == p
                and a + b == q
            ):
                return iterations, a, b

        a += 1
        iterations += 1


def y_square_search(case):
    """
    Search y = 0,1,2,... until

        discriminant(y)

    becomes a square.
    """

    y = 0
    tested = 0

    while True:
        delta = discriminant(
            case,
            y,
        )

        if is_square(delta):
            return tested, y, math.isqrt(delta)

        y += 1
        tested += 1


def theoretical_fermat_A(case):
    return (
        case["p"] + case["q"]
    ) // 2


def theoretical_fermat_B(case):
    return (
        case["q"] - case["p"]
    ) // 2


def verify_exact_relation(case):
    y = case["y"]

    delta = discriminant(
        case,
        y,
    )

    expected = case["gap"] ** 2

    return delta == expected


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(494949)

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

    print("============================================================")
    print("TEST 1: DISCRIMINANT IDENTITY")
    print("============================================================")

    failures = 0

    for case in cases:
        if not verify_exact_relation(case):
            failures += 1

    print(
        f"failures = {failures}/{CASE_COUNT}"
    )

    print()

    print("============================================================")
    print("TEST 2: y-SEARCH VS FERMAT")
    print("============================================================")

    y_mismatches = 0
    iteration_mismatches = 0

    for case in cases:
        y_tests, found_y, found_gap = y_square_search(
            case
        )

        fermat_tests, a, b = fermat_iterations(
            case["N"],
            case["p"],
            case["q"],
        )

        expected_a = theoretical_fermat_A(
            case
        )

        expected_b = theoretical_fermat_B(
            case
        )

        if found_y != case["y"]:
            y_mismatches += 1

        if (
            a != expected_a
            or b != expected_b
        ):
            iteration_mismatches += 1

        if y_tests != case["y"]:
            iteration_mismatches += 1

    print(
        f"y-search mismatches      = "
        f"{y_mismatches}/{CASE_COUNT}"
    )

    print(
        f"Fermat/result mismatches = "
        f"{iteration_mismatches}/{CASE_COUNT}"
    )

    print()

    print("============================================================")
    print("TEST 3: RELATION BETWEEN y AND FERMAT A")
    print("============================================================")

    failures = 0

    for case in cases:
        a = theoretical_fermat_A(case)

        expected_a = case["s"] + (
            case["y"] + 1
        ) // 2

        if a != expected_a:
            failures += 1

    print(
        f"A = s + ceil(y/2) failures = "
        f"{failures}/{CASE_COUNT}"
    )

    print()

    print("============================================================")
    print("TEST 4: ITERATION COUNT VS y")
    print("============================================================")

    differences = []

    for case in cases:
        fermat_tests, _, _ = fermat_iterations(
            case["N"],
            case["p"],
            case["q"],
        )

        y_tests, _, _ = y_square_search(
            case
        )

        differences.append(
            fermat_tests - y_tests
        )

    print(
        f"minimum difference = "
        f"{min(differences)}"
    )

    print(
        f"maximum difference = "
        f"{max(differences)}"
    )

    print(
        f"mean difference    = "
        f"{sum(differences) / len(differences):.4f}"
    )

    print()

    print("============================================================")
    print("TEST 5: SQUARE-GAP FORM")
    print("============================================================")

    for case in cases[:10]:
        s = case["s"]
        y = case["y"]
        gap = case["gap"]

        z = 2 * s + 1 + y

        print()
        print(f"N       = {case['N']}")
        print(f"s       = {s}")
        print(f"y       = {y}")
        print(f"z       = {z}")
        print(f"q-p     = {gap}")
        print(f"z^2-4N  = {z*z - 4*case['N']}")
        print(f"(q-p)^2 = {gap*gap}")
        print(
            f"match   = "
            f"{z*z - 4*case['N'] == gap*gap}"
        )

    print()

    print("============================================================")
    print("TEST 6: FACTORIZATION FROM y")
    print("============================================================")

    failures = 0

    for case in cases:
        y = case["y"]
        delta = discriminant(
            case,
            y,
        )

        gap = math.isqrt(delta)

        z = 2 * case["s"] + 1 + y

        p_recovered = (z - gap) // 2
        q_recovered = (z + gap) // 2

        if (
            p_recovered != case["p"]
            or q_recovered != case["q"]
        ):
            failures += 1

    print(
        f"factor recovery failures = "
        f"{failures}/{CASE_COUNT}"
    )

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
