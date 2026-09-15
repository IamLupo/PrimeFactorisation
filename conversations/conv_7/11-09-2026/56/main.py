import math
import random

EXPERIMENT = 74


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
        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "y": y,
        }


def crt_two(a, p, b, q):
    t = (
        (b - a)
        * pow(p, -1, q)
    ) % q

    return (a + p * t) % (p * q)


def get_gap(case):
    p = case["p"]
    q = case["q"]
    n = case["N"]
    s = case["s"]

    xp = (s + 1) % p
    xq = (s + 1) % q

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
                crt_two(a, p, b, q)
            )

    roots = sorted(set(roots))

    known = s + 1

    unknown = [
        r for r in roots
        if r != known and r != n - s
    ]

    unknown.sort()

    assert len(unknown) == 2

    return unknown[0] - known


def main():
    random.seed(74001)

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
    print("TEST 1: CRT MULTIPLIER FACTORIZATION")
    print("=" * 60)

    failures = 0

    for case in cases:
        s = case["s"]
        n = case["N"]

        A = get_gap(case)

        g1 = math.gcd(A, n)
        g2 = math.gcd(A + 2 * s + 1, n)

        r = A // g1
        t = (A + 2 * s + 1) // g2

        quotient = (
            A * (A + 2 * s + 1)
        ) // n

        if r * t != quotient:
            failures += 1

    print(
        f"r*t = quotient failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: CRT CONGRUENCE MULTIPLIERS")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]
        n = case["N"]

        A = get_gap(case)

        g1 = math.gcd(A, n)
        g2 = math.gcd(A + 2 * s + 1, n)

        r = A // g1
        t = (A + 2 * s + 1) // g2

        # A is divisible by g1 and A+2s+1 by g2.
        if A % g1 != 0:
            failures += 1

        if (A + 2 * s + 1) % g2 != 0:
            failures += 1

        # The complementary factors are different.
        if g1 == g2:
            failures += 1

        # Both CRT multipliers are positive.
        if r <= 0 or t <= 0:
            failures += 1

        # Product of the complementary prime divisors.
        if g1 * g2 != n:
            failures += 1

    print(
        f"structural failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: MULTIPLIER VALUES")
    print("=" * 60)

    r_values = []
    t_values = []
    q_values = []

    for case in cases:
        s = case["s"]
        n = case["N"]

        A = get_gap(case)

        g1 = math.gcd(A, n)
        g2 = math.gcd(A + 2 * s + 1, n)

        r = A // g1
        t = (A + 2 * s + 1) // g2

        r_values.append(r)
        t_values.append(t)
        q_values.append(r * t)

    print(
        f"r minimum    = {min(r_values)}"
    )
    print(
        f"r maximum    = {max(r_values)}"
    )
    print(
        f"r mean       = "
        f"{sum(r_values) / len(r_values):.6f}"
    )
    print(
        f"r median     = "
        f"{sorted(r_values)[len(r_values)//2]}"
    )

    print()

    print(
        f"t minimum    = {min(t_values)}"
    )
    print(
        f"t maximum    = {max(t_values)}"
    )
    print(
        f"t mean       = "
        f"{sum(t_values) / len(t_values):.6f}"
    )
    print(
        f"t median     = "
        f"{sorted(t_values)[len(t_values)//2]}"
    )

    print()

    print(
        f"Q minimum    = {min(q_values)}"
    )
    print(
        f"Q maximum    = {max(q_values)}"
    )
    print(
        f"Q mean       = "
        f"{sum(q_values) / len(q_values):.6f}"
    )
    print(
        f"Q median     = "
        f"{sorted(q_values)[len(q_values)//2]}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: MULTIPLIER AGAINST KNOWN QUANTITIES")
    print("=" * 60)

    matches = {
        "r=y": 0,
        "t=y": 0,
        "r=y+1": 0,
        "t=y+1": 0,
        "r=D": 0,
        "t=D": 0,
        "r=q-p": 0,
        "t=q-p": 0,
        "r=s": 0,
        "t=s": 0,
        "r=s+1": 0,
        "t=s+1": 0,
        "r=p": 0,
        "t=p": 0,
        "r=q": 0,
        "t=q": 0,
    }

    for case in cases:
        s = case["s"]
        p = case["p"]
        q = case["q"]
        y = case["y"]

        D = case["N"] - s * s

        A = get_gap(case)

        g1 = math.gcd(A, case["N"])
        g2 = math.gcd(
            A + 2 * s + 1,
            case["N"],
        )

        r = A // g1
        t = (
            A + 2 * s + 1
        ) // g2

        tests = {
            "r=y": r == y,
            "t=y": t == y,
            "r=y+1": r == y + 1,
            "t=y+1": t == y + 1,
            "r=D": r == D,
            "t=D": t == D,
            "r=q-p": r == q - p,
            "t=q-p": t == q - p,
            "r=s": r == s,
            "t=s": t == s,
            "r=s+1": r == s + 1,
            "t=s+1": t == s + 1,
            "r=p": r == p,
            "t=p": t == p,
            "r=q": r == q,
            "t=q": t == q,
        }

        for name, value in tests.items():
            if value:
                matches[name] += 1

    for name, count in matches.items():
        print(
            f"{name:12s} = {count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: LINEAR RELATIONS")
    print("=" * 60)

    relations = {
        "r+t = y": 0,
        "r+t = y+1": 0,
        "r+t = 2s+1": 0,
        "r+t = p+q": 0,
        "r+t = q-p": 0,
        "r-t = y": 0,
        "r-t = q-p": 0,
        "r*t = N": 0,
        "r*t = D": 0,
    }

    for case in cases:
        s = case["s"]
        p = case["p"]
        q = case["q"]
        y = case["y"]
        n = case["N"]
        D = n - s * s

        A = get_gap(case)

        g1 = math.gcd(A, n)
        g2 = math.gcd(
            A + 2 * s + 1,
            n,
        )

        r = A // g1
        t = (A + 2 * s + 1) // g2

        checks = {
            "r+t = y": r + t == y,
            "r+t = y+1": r + t == y + 1,
            "r+t = 2s+1": r + t == 2 * s + 1,
            "r+t = p+q": r + t == p + q,
            "r+t = q-p": r + t == q - p,
            "r-t = y": abs(r - t) == y,
            "r-t = q-p": abs(r - t) == q - p,
            "r*t = N": r * t == n,
            "r*t = D": r * t == D,
        }

        for name, value in checks.items():
            if value:
                relations[name] += 1

    for name, count in relations.items():
        print(
            f"{name:18s} = {count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: MODULAR-INVERSE FORM")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]
        n = case["N"]

        A = get_gap(case)

        g = math.gcd(A, n)

        other = n // g

        multiplier = (
            A // g
        )

        # A is 0 mod g and
        # A = -(2s+1) mod other.
        predicted = (
            -(
                2 * s + 1
            )
            * pow(
                g,
                -1,
                other,
            )
        ) % other

        if multiplier % other != predicted:
            failures += 1

    print(
        f"modular-inverse failures = "
        f"{failures}"
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
        s = case["s"]
        n = case["N"]

        A = get_gap(case)

        g1 = math.gcd(A, n)
        g2 = math.gcd(
            A + 2 * s + 1,
            n,
        )

        r = A // g1
        t = (
            A + 2 * s + 1
        ) // g2

        print()
        print(f"CASE {number}")
        print(f"N = {n}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {s}")
        print(f"y = {case['y']}")
        print()

        print(
            "A =",
            A,
        )

        print(
            "g1 =",
            g1,
        )

        print(
            "g2 =",
            g2,
        )

        print(
            "r = A/g1 =",
            r,
        )

        print(
            "t = (A+2s+1)/g2 =",
            t,
        )

        print(
            "r*t =",
            r * t,
        )

        print(
            "A(A+2s+1)/N =",
            (
                A
                * (A + 2 * s + 1)
            ) // n,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
