import math
import random


EXPERIMENT = 52


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

        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "y": y,
        }


def diagonal_value(case, r, x):
    """
    Since z = r-s and y = z-x,

        P(x,z-x)
        = N - r(s+1-x).

    This is the linear diagonal form.
    """

    n = case["N"]
    s = case["s"]

    return (
        n
        - r * (s + 1 - x)
    )


def direct_diagonal_value(case, r, x):
    n = case["N"]
    s = case["s"]
    d = case["D"]

    z = r - s
    y = z - x

    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def binomial_transform_direct(case, r, n):
    """
    Directly evaluate

        sum_{x=0}^n
            P(x,z-x) C(n,x) 2^(n-x).
    """

    total = 0

    for x in range(n + 1):
        value = diagonal_value(
            case,
            r,
            x,
        )

        total += (
            value
            * math.comb(n, x)
            * (2 ** (n - x))
        )

    return total


def binomial_transform_closed(case, r, n):
    """
    Exact closed form:

        B_n(r)
        =
        3^(n-1)
        [
            3N + r(n - 3(s+1))
        ]
    """

    N = case["N"]
    s = case["s"]

    if n == 0:
        return (
            N
            - r * (s + 1)
        )

    return (
        3 ** (n - 1)
        * (
            3 * N
            + r * (
                n
                - 3 * (s + 1)
            )
        )
    )


def transformed_difference(case, r, n):
    """
    B_n(r) - 3^n N.
    """

    B = binomial_transform_closed(
        case,
        r,
        n,
    )

    return (
        B
        - (3 ** n) * case["N"]
    )


def predicted_gcd(case, r, n):
    N = case["N"]

    value = transformed_difference(
        case,
        r,
        n,
    )

    return math.gcd(
        value,
        N,
    )


def coefficient(case, n):
    """
    c_n = n - 3(s+1)
    """

    return (
        n
        - 3 * (
            case["s"] + 1
        )
    )


def test_transform_identity(cases):
    print("============================================================")
    print("TEST 1: DIRECT VS CLOSED BINOMIAL TRANSFORM")
    print("============================================================")

    failures = 0

    for case in cases:
        for r in [case["p"], case["q"]]:
            for n in range(1, 12):
                direct = (
                    binomial_transform_direct(
                        case,
                        r,
                        n,
                    )
                )

                closed = (
                    binomial_transform_closed(
                        case,
                        r,
                        n,
                    )
                )

                if direct != closed:
                    failures += 1

    print(
        f"closed-form failures = {failures}"
    )

    print()


def test_factor_gcd(cases):
    print("============================================================")
    print("TEST 2: FACTOR GCD SIGNAL")
    print("============================================================")

    for n in range(1, 17):
        p_hits = 0
        q_hits = 0
        gcd_n = 0

        for case in cases:
            p = case["p"]
            q = case["q"]
            N = case["N"]

            gp = predicted_gcd(
                case,
                p,
                n,
            )

            gq = predicted_gcd(
                case,
                q,
                n,
            )

            if gp == p:
                p_hits += 1

            if gq == q:
                q_hits += 1

            if gp == N:
                gcd_n += 1

            if gq == N:
                gcd_n += 1

        print(
            f"n={n:2d}  "
            f"p recovered={p_hits:4d}/"
            f"{len(cases)}  "
            f"q recovered={q_hits:4d}/"
            f"{len(cases)}  "
            f"gcd=N={gcd_n:4d}"
        )

    print()


def test_coprime_condition(cases):
    print("============================================================")
    print("TEST 3: COPRIME COEFFICIENT CONDITION")
    print("============================================================")

    for n in range(1, 17):
        coprime = 0
        p_exact = 0
        q_exact = 0

        for case in cases:
            N = case["N"]
            p = case["p"]
            q = case["q"]

            c = coefficient(
                case,
                n,
            )

            if math.gcd(c, N) == 1:
                coprime += 1

                if predicted_gcd(
                    case,
                    p,
                    n,
                ) == p:
                    p_exact += 1

                if predicted_gcd(
                    case,
                    q,
                    n,
                ) == q:
                    q_exact += 1

        print(
            f"n={n:2d}  "
            f"gcd(c,N)=1: {coprime:4d}/"
            f"{len(cases)}  "
            f"p exact={p_exact:4d}  "
            f"q exact={q_exact:4d}"
        )

    print()


def test_unknown_diagonal(cases):
    print("============================================================")
    print("TEST 4: ARBITRARY DIAGONALS")
    print("============================================================")

    examples = []

    for case in cases:
        N = case["N"]

        # Pick several arbitrary r values.
        values = [
            1,
            case["s"],
            case["p"],
            case["q"],
            N,
        ]

        for r in values:
            n = 5

            g = predicted_gcd(
                case,
                r,
                n,
            )

            expected = math.gcd(
                r * coefficient(
                    case,
                    n,
                ),
                N,
            )

            if g != expected:
                examples.append(
                    (
                        case,
                        r,
                        g,
                        expected,
                    )
                )

    print(
        f"arbitrary diagonal failures = "
        f"{len(examples)}"
    )

    print()


def test_minimal_n(cases):
    print("============================================================")
    print("TEST 5: FIRST n WITH gcd(c_n,N)=1")
    print("============================================================")

    distribution = {}

    for case in cases:
        found = None

        for n in range(1, 100):
            c = coefficient(
                case,
                n,
            )

            if math.gcd(
                c,
                case["N"],
            ) == 1:
                found = n
                break

        distribution[found] = (
            distribution.get(found, 0)
            + 1
        )

    for n in sorted(distribution):
        print(
            f"first n={n}: "
            f"{distribution[n]}"
        )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 6: TRUE FACTOR DIAGONALS")
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
            z = r - case["s"]

            print()
            print(f"r={r}")
            print(f"z=r-s={z}")

            for n in [1, 2, 3, 4, 5, 8]:
                c = coefficient(
                    case,
                    n,
                )

                g = predicted_gcd(
                    case,
                    r,
                    n,
                )

                print(
                    f"n={n:2d}  "
                    f"c={c:8d}  "
                    f"gcd={g}"
                )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(525252)

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

    test_transform_identity(cases)
    test_factor_gcd(cases)
    test_coprime_condition(cases)
    test_unknown_diagonal(cases)
    test_minimal_n(cases)
    print_examples(cases)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
