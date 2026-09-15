import math
import random

EXPERIMENT = 64


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
        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "xp": xp,
            "y": y,
        }


def P0(case, x):
    return (
        x * x
        - x
        + case["D"]
        - case["s"]
    )


def H(case, a, x):
    s = case["s"]
    d = case["D"]

    return (
        (a + 2) * x * x
        - (a + s + 2) * x
        + (a + 1) * (d - s)
    )


def direct_transform(case, a, x):
    total = 0

    for y in range(x + 1):
        weight = (
            math.comb(x, y)
            * (a ** (x - y))
        )

        p_value = (
            x * x
            + (y - 1) * x
            + (case["D"] - case["s"])
            - (case["s"] + 1) * y
        )

        total += p_value * weight

    return total


def closed_transform(case, a, x):
    if x == 0:
        return (
            (a + 1) ** 0
            * H(case, a, x)
        )

    return (
        (a + 1) ** (x - 1)
        * H(case, a, x)
    )


def main():
    random.seed(64001)

    bits = 18
    cases_count = 500

    bases = [
        1,
        2,
        3,
        4,
        5,
        7,
        11,
    ]

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"cases      = {cases_count}")
    print(f"prime bits = {bits}")
    print(f"bases      = {bases}")
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: CORRECTED GENERALIZED TRANSFORM")
    print("=" * 60)

    failures = 0

    for case in cases[:100]:
        for x in range(1, 41):
            for a in bases:
                direct = direct_transform(
                    case,
                    a,
                    x,
                )

                closed = closed_transform(
                    case,
                    a,
                    x,
                )

                if direct != closed:
                    failures += 1

    print(
        f"identity failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: BASE DIFFERENCE IDENTITY")
    print("=" * 60)

    failures = 0

    for case in cases:
        for x in range(0, min(case["s"] + 2, 200)):
            b_value = P0(case, x)

            for i in range(len(bases)):
                for j in range(i + 1, len(bases)):
                    a = bases[i]
                    b = bases[j]

                    lhs = (
                        H(case, b, x)
                        - H(case, a, x)
                    )

                    rhs = (
                        (b - a)
                        * b_value
                    )

                    if lhs != rhs:
                        failures += 1

    print(
        f"difference identity failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE x_p COLLAPSE")
    print("=" * 60)

    failures = 0

    for case in cases:
        xp = case["xp"]
        p = case["p"]
        y = case["y"]

        expected = p * y
        actual = P0(case, xp)

        if actual != expected:
            failures += 1

    print(
        f"P0(x_p) = p*y failures = "
        f"{failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: BASE DIFFERENCES AT TRUE x_p")
    print("=" * 60)

    failures = 0

    for case in cases:
        xp = case["xp"]
        p = case["p"]
        y = case["y"]

        for i in range(len(bases)):
            for j in range(i + 1, len(bases)):
                a = bases[i]
                b = bases[j]

                lhs = (
                    H(case, b, xp)
                    - H(case, a, xp)
                )

                rhs = (
                    (b - a)
                    * p
                    * y
                )

                if lhs != rhs:
                    failures += 1

    print(
        f"true-x difference failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: gcd INFORMATION")
    print("=" * 60)

    recovered_p = 0
    recovered_q = 0
    trivial = 0
    other = 0

    for case in cases:
        xp = case["xp"]
        n = case["N"]

        a = 1
        b = 2

        difference = (
            H(case, b, xp)
            - H(case, a, xp)
        )

        g = math.gcd(
            difference,
            n,
        )

        if g == case["p"]:
            recovered_p += 1

        elif g == case["q"]:
            recovered_q += 1

        elif g == 1 or g == n:
            trivial += 1

        else:
            other += 1

    print(
        f"recovered p = {recovered_p}/{cases_count}"
    )

    print(
        f"recovered q = {recovered_q}/{cases_count}"
    )

    print(
        f"trivial     = {trivial}/{cases_count}"
    )

    print(
        f"other       = {other}/{cases_count}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: MULTI-BASE gcd VS SINGLE P0")
    print("=" * 60)

    failures = 0

    for case in cases:
        xp = case["xp"]
        n = case["N"]

        difference = (
            H(case, 2, xp)
            - H(case, 1, xp)
        )

        direct = P0(case, xp)

        g1 = math.gcd(
            difference,
            n,
        )

        g2 = math.gcd(
            direct,
            n,
        )

        if g1 != g2:
            failures += 1

    print(
        f"gcd equivalence failures = "
        f"{failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: RANDOM x FALSE-POSITIVE TEST")
    print("=" * 60)

    random.seed(64002)

    total = 0
    factor_hits = 0
    trivial_hits = 0

    for case in cases:
        xp = case["xp"]
        s = case["s"]
        n = case["N"]

        # Test random x values other than x_p.
        for _ in range(100):
            x = random.randint(
                0,
                s + 1,
            )

            if x == xp:
                continue

            total += 1

            g = math.gcd(
                P0(case, x),
                n,
            )

            if g == case["p"] or g == case["q"]:
                factor_hits += 1

            if g == 1 or g == n:
                trivial_hits += 1

    print(
        f"random x tested = {total}"
    )

    print(
        f"factor gcd hits = {factor_hits}"
    )

    print(
        f"trivial gcd     = {trivial_hits}"
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
        xp = case["xp"]

        print()
        print(f"CASE {number}")
        print(f"N       = {case['N']}")
        print(f"p       = {case['p']}")
        print(f"q       = {case['q']}")
        print(f"s       = {case['s']}")
        print(f"D       = {case['D']}")
        print(f"x_p     = {xp}")
        print(f"y       = {case['y']}")
        print()
        print(
            "P0(x_p) =",
            P0(case, xp),
        )
        print(
            "p*y     =",
            case["p"] * case["y"],
        )

        for a in bases:
            print(
                f"H_{a}(x_p) =",
                H(case, a, xp),
            )

        print(
            "H_2 - H_1 =",
            H(case, 2, xp)
            - H(case, 1, xp),
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
