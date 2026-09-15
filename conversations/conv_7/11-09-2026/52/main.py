import math
import random

EXPERIMENT = 70


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
    high = (1 << bits) - 1

    while True:
        n = random.randint(low, high)
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

        xp = s + 1 - p
        xq = s + 1 - q
        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "xp": xp,
            "xq": xq,
            "y": y,
        }


def P(case, x, y):
    return (
        x * x
        + (y - 1) * x
        + (case["D"] - case["s"])
        - (case["s"] + 1) * y
    )


def slope(case, x, y):
    return P(case, x + 1, y) - P(case, x, y)


def second_difference(case, x, y):
    return (
        P(case, x + 1, y)
        - 2 * P(case, x, y)
        + P(case, x - 1, y)
    )


def main():
    random.seed(70001)

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
    print("TEST 1: FIRST DIFFERENCE IDENTITY")
    print("=" * 60)

    failures = 0

    for case in cases:
        for _ in range(100):
            x = random.randint(
                -100,
                case["s"] + 100,
            )

            y = random.randint(
                0,
                case["y"] + 100,
            )

            lhs = slope(case, x, y)

            rhs = 2 * x + y

            if lhs != rhs:
                failures += 1

    print(
        f"identity failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: TRUE x_p SLOPE")
    print("=" * 60)

    failures = 0

    for case in cases:
        xp = case["xp"]
        y = case["y"]

        value = slope(
            case,
            xp,
            y,
        )

        expected = (
            case["q"]
            - case["p"]
            + 1
        )

        if value != expected:
            failures += 1

    print(
        f"q-p+1 failures = "
        f"{failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE x_q SLOPE")
    print("=" * 60)

    failures = 0

    for case in cases:
        xq = case["xq"]
        y = case["y"]

        value = slope(
            case,
            xq,
            y,
        )

        expected = (
            case["p"]
            - case["q"]
            + 1
        )

        if value != expected:
            failures += 1

    print(
        f"p-q+1 failures = "
        f"{failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: SECOND DIFFERENCE")
    print("=" * 60)

    failures = 0

    for case in cases:
        for _ in range(100):
            x = random.randint(
                -100,
                case["s"] + 100,
            )

            y = random.randint(
                0,
                case["y"] + 100,
            )

            value = second_difference(
                case,
                x,
                y,
            )

            if value != 2:
                failures += 1

    print(
        f"second-difference failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 5
    #
    # At the true root:
    #
    # P(xp,y) = 0
    # P(xp+1,y) = q-p+1
    # P(xp-1,y) = -(q-p-1)
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: LOCAL VALUES AT x_p")
    print("=" * 60)

    failures = 0

    for case in cases:
        xp = case["xp"]
        y = case["y"]

        left = P(
            case,
            xp - 1,
            y,
        )

        center = P(
            case,
            xp,
            y,
        )

        right = P(
            case,
            xp + 1,
            y,
        )

        expected_left = (
            -(case["q"] - case["p"] - 1)
        )

        expected_center = 0

        expected_right = (
            case["q"] - case["p"] + 1
        )

        if left != expected_left:
            failures += 1

        if center != expected_center:
            failures += 1

        if right != expected_right:
            failures += 1

    print(
        f"local-value failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 6
    #
    # Test gcd combinations of the local slope.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: gcd OF LOCAL GAP EXPRESSIONS")
    print("=" * 60)

    combinations = {
        "gcd(d, N)": 0,
        "gcd(d+1, N)": 0,
        "gcd(d-1, N)": 0,
        "gcd(d^2-1, N)": 0,
        "gcd(d^2, N)": 0,
    }

    for case in cases:
        p = case["p"]
        q = case["q"]
        n = case["N"]

        d = q - p

        values = {
            "gcd(d, N)": d,
            "gcd(d+1, N)": d + 1,
            "gcd(d-1, N)": d - 1,
            "gcd(d^2-1, N)": d * d - 1,
            "gcd(d^2, N)": d * d,
        }

        for name, value in values.items():
            g = math.gcd(
                value,
                n,
            )

            if g != 1 and g != n:
                combinations[name] += 1

    for name, count in combinations.items():
        print(
            f"{name:20s} "
            f"nontrivial = {count}/{cases_count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 7
    #
    # The most interesting combination:
    #
    # (q-p+1)(q-p-1) = (q-p)^2 - 1
    #
    # Compare it with N.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: (q-p)^2 - 1")
    print("=" * 60)

    p_hits = 0
    q_hits = 0
    trivial = 0
    other = 0

    for case in cases:
        d = case["q"] - case["p"]

        value = d * d - 1

        g = math.gcd(
            value,
            case["N"],
        )

        if g == case["p"]:
            p_hits += 1
        elif g == case["q"]:
            q_hits += 1
        elif g == 1 or g == case["N"]:
            trivial += 1
        else:
            other += 1

    print(
        f"gcd = p : {p_hits}"
    )

    print(
        f"gcd = q : {q_hits}"
    )

    print(
        f"trivial  : {trivial}"
    )

    print(
        f"other    : {other}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 8
    #
    # Random x control.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 8: RANDOM SLOPE CONTROL")
    print("=" * 60)

    total = 0
    factor_hits = 0

    for case in cases:
        n = case["N"]

        for _ in range(100):
            x = random.randint(
                0,
                case["s"] + 1,
            )

            y = random.randint(
                0,
                case["y"] + 100,
            )

            d = slope(
                case,
                x,
                y,
            )

            g = math.gcd(
                d,
                n,
            )

            total += 1

            if g == case["p"] or g == case["q"]:
                factor_hits += 1

    print(
        f"random slopes = {total}"
    )

    print(
        f"factor gcd hits = {factor_hits}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 9
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 9: EXPLICIT CASES")
    print("=" * 60)

    for number, case in enumerate(
        cases[:10],
        start=1,
    ):
        xp = case["xp"]
        xq = case["xq"]
        y = case["y"]

        dp = slope(
            case,
            xp,
            y,
        )

        dq = slope(
            case,
            xq,
            y,
        )

        print()
        print(f"CASE {number}")
        print(f"N = {case['N']}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {case['s']}")
        print(f"x_p = {xp}")
        print(f"x_q = {xq}")
        print(f"y = {y}")
        print()
        print(
            "slope(x_p) =",
            dp,
        )
        print(
            "expected   =",
            case["q"] - case["p"] + 1,
        )
        print(
            "slope(x_q) =",
            dq,
        )
        print(
            "expected   =",
            case["p"] - case["q"] + 1,
        )
        print(
            "gcd(slope(x_p),N) =",
            math.gcd(
                dp,
                case["N"],
            ),
        )
        print(
            "gcd((q-p)^2-1,N) =",
            math.gcd(
                (case["q"] - case["p"]) ** 2 - 1,
                case["N"],
            ),
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
