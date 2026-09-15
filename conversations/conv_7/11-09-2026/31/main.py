import math
import random


EXPERIMENT = 51


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
        d = n - s * s

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
        }


def P(n, s, d, x, y):
    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def diagonal_P(n, s, d, x, z):
    """
    P(x, z-x).
    """

    y = z - x

    return P(
        n,
        s,
        d,
        x,
        y,
    )


def diagonal_closed_form(n, s, d, x, z):
    """
    Exact identity:

        P(x,z-x)
        =
        (z+s)x + D-s-(s+1)z.
    """

    return (
        (z + s) * x
        + d
        - s
        - (s + 1) * z
    )


def diagonal_root(n, s, d, z):
    """
    Solve P(x,z-x)=0 for x.

        x = s+1 - N/(z+s)

    Returns None if z+s does not divide N.
    """

    r = z + s

    if r == 0:
        return None

    if n % r != 0:
        return None

    x = s + 1 - n // r

    return x


def divisor_parameterization(n, s, r):
    """
    Given a divisor r of N:

        z = r-s
        x = s+1-N/r
        y = z-x

    """

    if r == 0 or n % r != 0:
        return None

    z = r - s
    x = s + 1 - n // r
    y = z - x

    return x, y, z


def all_divisors(n):
    divisors = []

    limit = math.isqrt(n)

    for d in range(1, limit + 1):
        if n % d != 0:
            continue

        divisors.append(d)

        other = n // d

        if other != d:
            divisors.append(other)

    return sorted(divisors)


def test_diagonal_identity(cases):
    print("============================================================")
    print("TEST 1: DIAGONAL LINEARIZATION")
    print("============================================================")

    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]
        d = case["D"]

        for _ in range(50):
            z = random.randint(
                -2 * s,
                2 * s,
            )

            x = random.randint(
                -2 * s,
                2 * s,
            )

            direct = diagonal_P(
                n,
                s,
                d,
                x,
                z,
            )

            closed = diagonal_closed_form(
                n,
                s,
                d,
                x,
                z,
            )

            if direct != closed:
                failures += 1

    print(
        f"linearization failures = "
        f"{failures}"
    )

    print()


def test_divisor_equivalence(cases):
    print("============================================================")
    print("TEST 2: P(x,z-x)=0 <=> z+s DIVIDES N")
    print("============================================================")

    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]
        d = case["D"]

        for z in range(
            -s + 1,
            s + 1,
        ):
            x = diagonal_root(
                n,
                s,
                d,
                z,
            )

            has_divisor = (
                (z + s) != 0
                and n % (z + s) == 0
            )

            if has_divisor:
                if x is None:
                    failures += 1
                    continue

                y = z - x

                if P(
                    n,
                    s,
                    d,
                    x,
                    y,
                ) != 0:
                    failures += 1

            else:
                if x is not None:
                    failures += 1

    print(
        f"equivalence failures = "
        f"{failures}"
    )

    print()


def test_factor_points(cases):
    print("============================================================")
    print("TEST 3: FACTOR DIVISORS")
    print("============================================================")

    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]
        p = case["p"]
        q = case["q"]

        for r, expected_other in [
            (p, q),
            (q, p),
        ]:
            result = divisor_parameterization(
                n,
                s,
                r,
            )

            if result is None:
                failures += 1
                continue

            x, y, z = result

            expected_x = (
                s + 1 - expected_other
            )

            expected_z = r - s

            if (
                x != expected_x
                or z != expected_z
                or x + y != z
                or P(
                    n,
                    s,
                    n - s * s,
                    x,
                    y,
                ) != 0
            ):
                failures += 1

    print(
        f"factor-divisor failures = "
        f"{failures}/{2 * len(cases)}"
    )

    print()


def test_all_divisors(cases):
    print("============================================================")
    print("TEST 4: ALL DIVISORS <-> INTEGER ROOTS")
    print("============================================================")

    total_divisors = 0
    total_roots = 0
    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]
        d = case["D"]

        divisors = all_divisors(n)

        total_divisors += len(divisors)

        roots = []

        for r in divisors:
            z = r - s

            x = diagonal_root(
                n,
                s,
                d,
                z,
            )

            if x is None:
                failures += 1
                continue

            y = z - x

            if P(
                n,
                s,
                d,
                x,
                y,
            ) != 0:
                failures += 1
                continue

            roots.append(
                (x, y, z, r)
            )

        total_roots += len(roots)

        if len(roots) != len(divisors):
            failures += 1

    print(
        f"total divisors = {total_divisors}"
    )

    print(
        f"total roots    = {total_roots}"
    )

    print(
        f"failures       = {failures}"
    )

    print()


def test_diagonal_slope(cases):
    print("============================================================")
    print("TEST 5: DIAGONAL SLOPE")
    print("============================================================")

    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]
        d = case["D"]

        divisors = all_divisors(n)

        for r in divisors:
            z = r - s

            x1 = 0
            x2 = 1

            f1 = diagonal_closed_form(
                n,
                s,
                d,
                x1,
                z,
            )

            f2 = diagonal_closed_form(
                n,
                s,
                d,
                x2,
                z,
            )

            slope = f2 - f1

            if slope != r:
                failures += 1

    print(
        f"slope failures = "
        f"{failures}"
    )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 6: FACTOR EXAMPLES")
    print("============================================================")

    for case in cases[:10]:
        n = case["N"]
        s = case["s"]
        p = case["p"]
        q = case["q"]

        print()
        print(f"N = {n}")
        print(f"p = {p}")
        print(f"q = {q}")
        print(f"s = {s}")

        for r in [p, q]:
            x, y, z = divisor_parameterization(
                n,
                s,
                r,
            )

            print()
            print(f"r = {r}")
            print(f"x = {x}")
            print(f"y = {y}")
            print(f"z = {z}")
            print(f"x+y = {x+y}")
            print(f"z+s = {z+s}")
            print(f"P = {P(n, s, case['D'], x, y)}")

            direct_slope = (
                diagonal_closed_form(
                    n,
                    s,
                    case["D"],
                    1,
                    z,
                )
                -
                diagonal_closed_form(
                    n,
                    s,
                    case["D"],
                    0,
                    z,
                )
            )

            print(
                f"diagonal slope = "
                f"{direct_slope}"
            )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(515151)

    CASE_COUNT = 100
    PRIME_BITS = 14

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = [
        generate_case(PRIME_BITS)
        for _ in range(CASE_COUNT)
    ]

    test_diagonal_identity(cases)
    test_divisor_equivalence(cases)
    test_factor_points(cases)
    test_all_divisors(cases)
    test_diagonal_slope(cases)
    print_examples(cases)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
