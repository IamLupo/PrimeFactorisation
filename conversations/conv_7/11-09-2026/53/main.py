import math
import random

EXPERIMENT = 71


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

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
        }


def P0(case, x):
    return (
        x * x
        - x
        + case["D"]
        - case["s"]
    )


def shifted_P0(case, k):
    s = case["s"]

    return P0(
        case,
        s + 1 + k,
    )


def shifted_closed(case, k):
    n = case["N"]
    s = case["s"]

    return (
        n
        + k * (2 * s + 1 + k)
    )


def factor_status(g, case):
    if g == case["p"]:
        return "p"

    if g == case["q"]:
        return "q"

    if g == case["N"]:
        return "N"

    if g == 1:
        return "1"

    return "other"


def main():
    random.seed(71001)

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
    print("TEST 1: SHIFTED P0 IDENTITY")
    print("=" * 60)

    failures = 0

    for case in cases:
        for k in range(-100, 101):
            direct = shifted_P0(
                case,
                k,
            )

            closed = shifted_closed(
                case,
                k,
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
    print("TEST 2: SMALL POSITIVE SHIFTS")
    print("=" * 60)

    hit_counts = {
        "p": 0,
        "q": 0,
        "N": 0,
        "other": 0,
        "1": 0,
    }

    smallest_hits = []

    for case in cases:
        first_hit = None

        for k in range(1, 65):
            value = (
                k
                * (
                    2 * case["s"]
                    + 1
                    + k
                )
            )

            g = math.gcd(
                value,
                case["N"],
            )

            status = factor_status(
                g,
                case,
            )

            hit_counts[status] += 1

            if status in ("p", "q"):
                if first_hit is None:
                    first_hit = k

        if first_hit is not None:
            smallest_hits.append(first_hit)

    for name in ["p", "q", "N", "other", "1"]:
        print(
            f"{name:5s} = {hit_counts[name]}"
        )

    print()

    print(
        f"cases with factor hit = "
        f"{len(smallest_hits)}/{cases_count}"
    )

    if smallest_hits:
        print(
            f"minimum k = "
            f"{min(smallest_hits)}"
        )

        print(
            f"maximum k = "
            f"{max(smallest_hits)}"
        )

        print(
            f"mean k    = "
            f"{sum(smallest_hits) / len(smallest_hits):.6f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: POWERS OF TWO")
    print("=" * 60)

    powers = [
        1 << j
        for j in range(0, 18)
    ]

    for k in powers:
        p_hits = 0
        q_hits = 0
        trivial = 0

        for case in cases:
            value = (
                k
                * (
                    2 * case["s"]
                    + 1
                    + k
                )
            )

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

        print(
            f"k={k:6d}  "
            f"p={p_hits:3d}  "
            f"q={q_hits:3d}  "
            f"trivial={trivial:3d}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: s MOD m")
    print("=" * 60)

    moduli = [
        2,
        3,
        4,
        5,
        7,
        8,
        11,
        13,
        16,
        17,
        19,
        23,
        31,
        32,
        64,
    ]

    for modulus in moduli:
        p_hits = 0
        q_hits = 0

        for case in cases:
            k = case["s"] % modulus

            if k == 0:
                k = modulus

            value = (
                k
                * (
                    2 * case["s"]
                    + 1
                    + k
                )
            )

            g = math.gcd(
                value,
                case["N"],
            )

            if g == case["p"]:
                p_hits += 1

            elif g == case["q"]:
                q_hits += 1

        print(
            f"m={modulus:2d}  "
            f"p={p_hits:3d}  "
            f"q={q_hits:3d}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: FLOOR(s/m)")
    print("=" * 60)

    divisors = [
        2,
        3,
        4,
        5,
        7,
        8,
        11,
        13,
        16,
        17,
        19,
        23,
        31,
        32,
    ]

    for m in divisors:
        p_hits = 0
        q_hits = 0

        for case in cases:
            k = case["s"] // m

            if k == 0:
                continue

            value = (
                k
                * (
                    2 * case["s"]
                    + 1
                    + k
                )
            )

            g = math.gcd(
                value,
                case["N"],
            )

            if g == case["p"]:
                p_hits += 1

            elif g == case["q"]:
                q_hits += 1

        print(
            f"m={m:2d}  "
            f"p={p_hits:3d}  "
            f"q={q_hits:3d}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: DIRECT gcd VS SHIFTED FORM")
    print("=" * 60)

    failures = 0

    for case in cases:
        n = case["N"]

        for k in range(1, 100):
            direct = math.gcd(
                shifted_P0(case, k),
                n,
            )

            closed = math.gcd(
                shifted_closed(case, k),
                n,
            )

            if direct != closed:
                failures += 1

    print(
        f"gcd equivalence failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: CONDITION FOR A p-HIT")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]

        # p divides shifted value iff
        #
        # k == 0 (mod p)
        #
        # or
        #
        # 2s+1+k == 0 (mod p).
        #
        # Test direct equivalence for small k.

        for k in range(1, 1000):
            value_hit = (
                k
                * (
                    2 * s
                    + 1
                    + k
                )
            ) % p == 0

            predicted = (
                k % p == 0
                or
                (2 * s + 1 + k) % p == 0
            )

            if value_hit != predicted:
                failures += 1

    print(
        f"p-divisibility equivalence failures = "
        f"{failures}"
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
        print()
        print(f"CASE {number}")
        print(f"N = {case['N']}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {case['s']}")
        print()

        for k in [
            1,
            2,
            3,
            4,
            8,
            16,
            32,
            64,
        ]:
            value = shifted_closed(
                case,
                k,
            )

            g = math.gcd(
                value,
                case["N"],
            )

            print(
                f"k={k:2d} "
                f"gcd={g} "
                f"status={factor_status(g, case)}"
            )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
