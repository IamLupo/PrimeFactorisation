import math
import random


EXPERIMENT = 73


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

        xp = s + 1 - p
        xq = s + 1 - q

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "xp": xp,
            "xq": xq,
        }


def crt_two(a, p, b, q):
    t = (
        (b - a)
        * pow(p, -1, q)
    ) % q

    return (a + p * t) % (p * q)


def get_nontrivial_gap(case):
    p = case["p"]
    q = case["q"]
    n = case["N"]
    s = case["s"]

    xp = case["xp"] % p
    xq = case["xq"] % q

    roots_p = [
        xp,
        (1 - xp) % p,
    ]

    roots_q = [
        xq,
        (1 - xq) % q,
    ]

    roots = []

    for a in roots_p:
        for b in roots_q:
            roots.append(
                crt_two(
                    a,
                    p,
                    b,
                    q,
                )
            )

    roots = sorted(set(roots))

    known1 = s + 1
    known2 = n - s

    unknown = [
        r
        for r in roots
        if r != known1 and r != known2
    ]

    unknown.sort()

    assert len(unknown) == 2

    A = unknown[0] - known1

    return A


def main():
    random.seed(73001)

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
    print("TEST 1: COMPLEMENTARY FACTOR SPLIT")
    print("=" * 60)

    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]

        A = get_nontrivial_gap(case)

        g1 = math.gcd(
            A,
            n,
        )

        g2 = math.gcd(
            A + 2 * s + 1,
            n,
        )

        if g1 * g2 != n:
            failures += 1

    print(
        f"g1*g2=N failures = "
        f"{failures}/{cases_count}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: DIRECT FACTOR RECOVERY")
    print("=" * 60)

    correct = 0
    swapped = 0
    other = 0

    for case in cases:
        n = case["N"]
        p = case["p"]
        q = case["q"]
        s = case["s"]

        A = get_nontrivial_gap(case)

        g1 = math.gcd(
            A,
            n,
        )

        g2 = math.gcd(
            A + 2 * s + 1,
            n,
        )

        if (
            {g1, g2}
            == {p, q}
        ):
            correct += 1

            if g1 == q and g2 == p:
                swapped += 1

        else:
            other += 1

    print(
        f"factor pair recovered = "
        f"{correct}/{cases_count}"
    )

    print(
        f"orientation q,p = "
        f"{swapped}/{cases_count}"
    )

    print(
        f"other = {other}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: EXACT PRODUCT QUOTIENT")
    print("=" * 60)

    quotients = []

    for case in cases:
        n = case["N"]
        s = case["s"]

        A = get_nontrivial_gap(case)

        value = (
            A
            * (A + 2 * s + 1)
        )

        assert value % n == 0

        quotients.append(
            value // n
        )

    print(
        f"minimum quotient = "
        f"{min(quotients)}"
    )

    print(
        f"maximum quotient = "
        f"{max(quotients)}"
    )

    print(
        f"mean quotient = "
        f"{sum(quotients) / len(quotients):.6f}"
    )

    print(
        f"median quotient = "
        f"{sorted(quotients)[len(quotients)//2]}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: QUOTIENT AGAINST y")
    print("=" * 60)

    matches = {
        "q = y": 0,
        "q = y+1": 0,
        "q = 2y": 0,
        "q = y(y+1)": 0,
        "q = p": 0,
        "q = q": 0,
    }

    for case in cases:
        n = case["N"]
        s = case["s"]
        y = case["xq"]  # overwritten below
        p = case["p"]
        q = case["q"]

        # Recover y from p,q only for classification.
        y = (
            p + q
            - 2 * s
            - 1
        )

        A = get_nontrivial_gap(case)

        quotient = (
            A
            * (A + 2 * s + 1)
        ) // n

        tests = {
            "q = y": y,
            "q = y+1": y + 1,
            "q = 2y": 2 * y,
            "q = y(y+1)": y * (y + 1),
            "q = p": p,
            "q = q": q,
        }

        for name, value in tests.items():
            if quotient == value:
                matches[name] += 1

    for name, count in matches.items():
        print(
            f"{name:16s} = {count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: SECOND FACTOR EXPRESSION")
    print("=" * 60)

    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]
        p = case["p"]
        q = case["q"]

        A = get_nontrivial_gap(case)

        g1 = math.gcd(
            A,
            n,
        )

        g2 = math.gcd(
            A + 2 * s + 1,
            n,
        )

        if g1 not in (p, q):
            failures += 1

        if g2 not in (p, q):
            failures += 1

        if g1 == g2:
            failures += 1

    print(
        f"split failures = "
        f"{failures}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: CAN FACTOR SPLIT BE EXPRESSED USING A")
    print("=" * 60)

    failures = 0

    for case in cases:
        n = case["N"]
        s = case["s"]

        A = get_nontrivial_gap(case)

        # From:
        #
        # A(A + 2s+1) = QN
        #
        # check the equivalent congruence.
        #
        lhs = (
            A
            * (A + 2 * s + 1)
        )

        if lhs % n != 0:
            failures += 1

    print(
        f"congruence failures = "
        f"{failures}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: RANDOM A CONTROL")
    print("=" * 60)

    total = 0
    split_hits = 0

    for case in cases:
        n = case["N"]

        for _ in range(100):
            A = random.randint(
                1,
                n - 1,
            )

            g1 = math.gcd(
                A,
                n,
            )

            g2 = math.gcd(
                A + 2 * case["s"] + 1,
                n,
            )

            total += 1

            if (
                g1 != 1
                and g1 != n
                and g2 != 1
                and g2 != n
                and g1 * g2 == n
            ):
                split_hits += 1

    print(
        f"random A tested = {total}"
    )

    print(
        f"factor split hits = "
        f"{split_hits}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 8
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 8: EXPLICIT CASES")
    print("=" * 60)

    for number, case in enumerate(
        cases[:10],
        start=1,
    ):
        n = case["N"]
        s = case["s"]

        A = get_nontrivial_gap(case)

        g1 = math.gcd(
            A,
            n,
        )

        g2 = math.gcd(
            A + 2 * s + 1,
            n,
        )

        quotient = (
            A
            * (A + 2 * s + 1)
        ) // n

        print()
        print(f"CASE {number}")
        print(f"N = {n}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {s}")
        print(f"A = {A}")
        print(
            "A+2s+1 =",
            A + 2 * s + 1,
        )
        print(
            "gcd(A,N) =",
            g1,
        )
        print(
            "gcd(A+2s+1,N) =",
            g2,
        )
        print(
            "product quotient =",
            quotient,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
