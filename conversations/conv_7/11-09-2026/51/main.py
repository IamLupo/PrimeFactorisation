import math
import random
import sympy as sp

EXPERIMENT = 69


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
    s = case["s"]
    d = case["D"]

    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def resultant_formula(case, k):
    return k * k * case["N"]


def direct_resultant_y(case, y, k):
    x = sp.symbols("x")

    P1 = (
        x * x
        + (y - 1) * x
        + (case["D"] - case["s"])
        - (case["s"] + 1) * y
    )

    y2 = y + k

    P2 = (
        x * x
        + (y2 - 1) * x
        + (case["D"] - case["s"])
        - (case["s"] + 1) * y2
    )

    return int(
        sp.resultant(P1, P2, x)
    )


def main():
    random.seed(69001)

    bits = 18
    cases_count = 500

    k_values = [
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
    print(f"k values   = {k_values}")
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: RESULTANT k=1")
    print("=" * 60)

    failures = 0

    for case in cases[:50]:
        y = random.randint(
            0,
            case["y"] + 100,
        )

        actual = direct_resultant_y(
            case,
            y,
            1,
        )

        expected = case["N"]

        if actual != expected:
            failures += 1

    print(
        f"resultant failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: GENERAL RESULTANT k^2 N")
    print("=" * 60)

    failures = 0

    for case in cases[:50]:
        y = random.randint(
            0,
            case["y"] + 100,
        )

        for k in k_values:
            actual = direct_resultant_y(
                case,
                y,
                k,
            )

            expected = resultant_formula(
                case,
                k,
            )

            if actual != expected:
                failures += 1

    print(
        f"resultant failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE x_p LOCAL FACTOR")
    print("=" * 60)

    failures = 0

    for case in cases:
        xp = case["xp"]
        y = case["y"]

        if P(case, xp, y) != 0:
            failures += 1
            continue

        value = P(case, xp, y + 1)

        if value != -case["p"]:
            failures += 1

    print(
        f"P(x_p,y+1) = -p failures = "
        f"{failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: TRUE x_q LOCAL FACTOR")
    print("=" * 60)

    failures = 0

    for case in cases:
        xq = case["xq"]
        y = case["y"]

        if P(case, xq, y) != 0:
            failures += 1
            continue

        value = P(case, xq, y + 1)

        if value != -case["q"]:
            failures += 1

    print(
        f"P(x_q,y+1) = -q failures = "
        f"{failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: DOUBLE gcd AT TRUE ROOT")
    print("=" * 60)

    p_hits = 0
    q_hits = 0
    trivial = 0

    for case in cases:
        xp = case["xp"]
        xq = case["xq"]
        y = case["y"]
        n = case["N"]

        gp = math.gcd(
            P(case, xp, y + 1),
            n,
        )

        gq = math.gcd(
            P(case, xq, y + 1),
            n,
        )

        if gp == case["p"]:
            p_hits += 1

        elif gp != 1:
            trivial += 1

        if gq == case["q"]:
            q_hits += 1

        elif gq != 1:
            trivial += 1

    print(
        f"p recovery = {p_hits}/{cases_count}"
    )

    print(
        f"q recovery = {q_hits}/{cases_count}"
    )

    print(
        f"unexpected = {trivial}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: RANDOM x,y DOUBLE gcd")
    print("=" * 60)

    factor_hits = 0
    trivial_hits = 0
    total = 0

    for case in cases:
        n = case["N"]

        for _ in range(100):
            x = random.randint(
                -case["s"],
                case["s"] + 1,
            )

            y = random.randint(
                0,
                case["y"] + 100,
            )

            k = random.choice(
                k_values
            )

            g = math.gcd(
                math.gcd(
                    P(case, x, y),
                    P(case, x, y + k),
                ),
                n,
            )

            total += 1

            if g == case["p"] or g == case["q"]:
                factor_hits += 1
            elif g == 1 or g == n:
                trivial_hits += 1

    print(
        f"total         = {total}"
    )

    print(
        f"factor hits   = {factor_hits}"
    )

    print(
        f"trivial hits  = {trivial_hits}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: ELIMINATING y")
    print("=" * 60)

    failures = 0

    for case in cases:
        for _ in range(100):
            x = random.randint(
                -case["s"],
                case["s"] + 1,
            )

            y = random.randint(
                0,
                case["y"] + 100,
            )

            difference = (
                P(case, x, y + 1)
                - P(case, x, y)
            )

            expected = (
                x - case["s"] - 1
            )

            if difference != expected:
                failures += 1

    print(
        f"difference failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 8
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 8: RESULTANT AS A FACTORIZATION CERTIFICATE")
    print("=" * 60)

    failures = 0

    for case in cases:
        n = case["N"]

        for k in k_values:
            r = resultant_formula(
                case,
                k,
            )

            if r % n != 0:
                failures += 1

            if r != k * k * n:
                failures += 1

    print(
        f"certificate failures = {failures}"
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

        print()
        print(f"CASE {number}")
        print(f"N       = {case['N']}")
        print(f"p       = {case['p']}")
        print(f"q       = {case['q']}")
        print(f"s       = {case['s']}")
        print(f"x_p     = {xp}")
        print(f"x_q     = {xq}")
        print(f"y       = {y}")
        print()

        print(
            "P(x_p,y)     =",
            P(case, xp, y),
        )

        print(
            "P(x_p,y+1)   =",
            P(case, xp, y + 1),
        )

        print(
            "P(x_q,y)     =",
            P(case, xq, y),
        )

        print(
            "P(x_q,y+1)   =",
            P(case, xq, y + 1),
        )

        print(
            "resultant k=1 =",
            resultant_formula(case, 1),
        )

        print(
            "resultant k=2 =",
            resultant_formula(case, 2),
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
