import math
import random

EXPERIMENT = 82


def is_prime(n):
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def random_prime(bits):
    low = 1 << (bits - 1)
    high = 1 << bits

    while True:
        n = random.randrange(low, high)
        n |= 1

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


def J(case, u):
    s = case["s"]
    d = case["D"]

    return (
        u * u
        + 4 * s * u
        - 4 * d
    )


def discriminant(case, y):
    return (
        2 * case["s"]
        + 1
        + y
    ) ** 2 - 4 * case["N"]


def main():
    random.seed(82001)

    bits = 18
    cases_count = 500

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"cases      = {cases_count}")
    print(f"prime bits = {bits}")
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: OFFSET-SUM IDENTITY")
    print("=" * 60)

    failures = 0

    for case in cases:
        u = case["y"] + 1

        expected = (
            case["q"] - case["p"]
        ) ** 2

        actual = J(case, u)

        if actual != expected:
            failures += 1

    print(
        f"identity failures = "
        f"{failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: EQUIVALENCE WITH DISCRIMINANT")
    print("=" * 60)

    failures = 0

    for case in cases:
        u = case["y"] + 1

        lhs = J(case, u)

        rhs = discriminant(
            case,
            case["y"],
        )

        if lhs != rhs:
            failures += 1

    print(
        f"equivalence failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: FIRST DIFFERENCE")
    print("=" * 60)

    failures = 0

    for case in cases:
        for u in range(
            0,
            min(case["y"] + 100, 1000),
        ):
            lhs = (
                J(case, u + 1)
                - J(case, u)
            )

            rhs = (
                2 * u
                + 1
                + 4 * case["s"]
            )

            if lhs != rhs:
                failures += 1

    print(
        f"difference failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: TRUE-VALUE gcds")
    print("=" * 60)

    expressions = {
        "J": [],
        "J-1": [],
        "J+1": [],
        "J-(q-p)^2": [],
        "first difference": [],
    }

    counts = {
        name: {
            "p": 0,
            "q": 0,
            "1": 0,
            "N": 0,
            "other": 0,
        }
        for name in expressions
    }

    for case in cases:
        n = case["N"]
        u = case["y"] + 1
        dgap = case["q"] - case["p"]

        values = {
            "J": J(case, u),
            "J-1": J(case, u) - 1,
            "J+1": J(case, u) + 1,
            "J-(q-p)^2": J(case, u) - dgap * dgap,
            "first difference": (
                J(case, u + 1)
                - J(case, u)
            ),
        }

        for name, value in values.items():
            g = math.gcd(
                value,
                n,
            )

            if g == case["p"]:
                counts[name]["p"] += 1
            elif g == case["q"]:
                counts[name]["q"] += 1
            elif g == 1:
                counts[name]["1"] += 1
            elif g == n:
                counts[name]["N"] += 1
            else:
                counts[name]["other"] += 1

    for name in counts:
        print()
        print(name)

        for key in [
            "p",
            "q",
            "1",
            "N",
            "other",
        ]:
            print(
                f"  gcd={key:5s} : "
                f"{counts[name][key]}"
            )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: RANDOM-u CONTROL")
    print("=" * 60)

    total = 0
    factor_hits = 0

    for case in cases:
        n = case["N"]

        for _ in range(100):
            u = random.randint(
                0,
                case["y"] + 100,
            )

            g = math.gcd(
                J(case, u),
                n,
            )

            total += 1

            if (
                g == case["p"]
                or g == case["q"]
            ):
                factor_hits += 1

    print(
        f"random u tested = {total}"
    )

    print(
        f"factor gcd hits = "
        f"{factor_hits}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: NEIGHBOURING J VALUES")
    print("=" * 60)

    k_values = [
        1,
        2,
        3,
        4,
        8,
        16,
        32,
    ]

    for k in k_values:
        hits = 0

        for case in cases:
            u = case["y"] + 1

            g = math.gcd(
                J(case, u),
                J(case, u + k),
            )

            g = math.gcd(
                g,
                case["N"],
            )

            if (
                g == case["p"]
                or g == case["q"]
            ):
                hits += 1

        print(
            f"k={k:2d} factor hits = "
            f"{hits}/{cases_count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: EXPLICIT CASES")
    print("=" * 60)

    for number, case in enumerate(
        cases[:10],
        start=1,
    ):
        u = case["y"] + 1

        print()
        print(f"CASE {number}")
        print(f"N = {case['N']}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {case['s']}")
        print(f"D = {case['D']}")
        print(f"y = {case['y']}")
        print(f"u = {u}")
        print()

        print(
            "J(u) =",
            J(case, u),
        )

        print(
            "(q-p)^2 =",
            (
                case["q"]
                - case["p"]
            ) ** 2,
        )

        print(
            "gcd(J,N) =",
            math.gcd(
                J(case, u),
                case["N"],
            ),
        )

        print(
            "J(u+1)-J(u) =",
            J(case, u + 1)
            - J(case, u),
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
