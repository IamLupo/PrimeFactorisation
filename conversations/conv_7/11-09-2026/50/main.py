import math
import random

EXPERIMENT = 68


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


def P0(case, x):
    return (
        x * x
        - x
        + case["D"]
        - case["s"]
    )


def P0_factored(case, x):
    n = case["N"]
    s = case["s"]

    return (
        n
        + (x - s - 1)
        * (x + s)
    )


def theoretical_hit_modulus(case):
    """
    Expected probability that two random evaluations
    share at least one prime factor.

    For prime r, P0(x) has two roots modulo r,
    so two independent x values both being roots
    has probability approximately 4/r^2.

    We use this only as a rough reference.
    """

    p = case["p"]
    q = case["q"]

    return (
        4.0 / (p * p)
        + 4.0 / (q * q)
    )


def main():
    random.seed(68001)

    bits = 18
    cases_count = 500

    x_samples_per_case = 200
    k_values = list(range(1, 65))

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"cases             = {cases_count}")
    print(f"prime bits        = {bits}")
    print(
        f"x samples/case   = "
        f"{x_samples_per_case}"
    )
    print(
        f"k values          = "
        f"{len(k_values)}"
    )
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: P0 FACTORIZATION IDENTITY")
    print("=" * 60)

    failures = 0

    for case in cases:
        for _ in range(100):
            x = random.randint(
                0,
                case["s"] + 1,
            )

            if (
                P0(case, x)
                != P0_factored(case, x)
            ):
                failures += 1

    print(
        f"identity failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: SINGLE-VALUE gcd")
    print("=" * 60)

    factor_hits = 0
    trivial_hits = 0

    for case in cases:
        n = case["N"]

        for _ in range(x_samples_per_case):
            x = random.randint(
                0,
                case["s"] + 1,
            )

            g = math.gcd(
                P0(case, x),
                n,
            )

            if g == case["p"] or g == case["q"]:
                factor_hits += 1

            elif g == 1 or g == n:
                trivial_hits += 1

    print(
        f"factor hits = {factor_hits}"
    )

    print(
        f"trivial hits = {trivial_hits}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: PAIRWISE gcd")
    print("=" * 60)

    pair_factor_hits = 0
    pair_trivial_hits = 0

    hit_by_k = {
        k: 0
        for k in k_values
    }

    hit_p_by_k = {
        k: 0
        for k in k_values
    }

    hit_q_by_k = {
        k: 0
        for k in k_values
    }

    total_pairs = 0

    for case in cases:
        n = case["N"]

        for _ in range(x_samples_per_case):
            x = random.randint(
                0,
                case["s"] + 1,
            )

            value_x = P0(
                case,
                x,
            )

            for k in k_values:
                x2 = x + k

                value_x2 = P0(
                    case,
                    x2,
                )

                g = math.gcd(
                    math.gcd(
                        value_x,
                        value_x2,
                    ),
                    n,
                )

                total_pairs += 1

                if (
                    g == case["p"]
                    or g == case["q"]
                ):
                    pair_factor_hits += 1
                    hit_by_k[k] += 1

                    if g == case["p"]:
                        hit_p_by_k[k] += 1

                    if g == case["q"]:
                        hit_q_by_k[k] += 1

                elif g == 1 or g == n:
                    pair_trivial_hits += 1

    print(
        f"total pairs = {total_pairs}"
    )

    print(
        f"factor hits = {pair_factor_hits}"
    )

    print(
        f"trivial hits = {pair_trivial_hits}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: HIT FREQUENCY BY k")
    print("=" * 60)

    for k in k_values:
        print(
            f"k={k:2d}  "
            f"hits={hit_by_k[k]:5d}  "
            f"p={hit_p_by_k[k]:5d}  "
            f"q={hit_q_by_k[k]:5d}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: DIFFERENCE gcd")
    print("=" * 60)

    difference_hits = 0
    difference_p = 0
    difference_q = 0
    difference_trivial = 0

    for case in cases:
        n = case["N"]

        for _ in range(x_samples_per_case):
            x = random.randint(
                0,
                case["s"] + 1,
            )

            for k in k_values:
                x2 = x + k

                d = (
                    P0(case, x2)
                    - P0(case, x)
                )

                g = math.gcd(
                    d,
                    n,
                )

                if g == case["p"]:
                    difference_hits += 1
                    difference_p += 1

                elif g == case["q"]:
                    difference_hits += 1
                    difference_q += 1

                elif g == 1 or g == n:
                    difference_trivial += 1

    print(
        f"factor hits = {difference_hits}"
    )

    print(
        f"p hits      = {difference_p}"
    )

    print(
        f"q hits      = {difference_q}"
    )

    print(
        f"trivial     = {difference_trivial}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: EXACT DIFFERENCE FORMULA")
    print("=" * 60)

    failures = 0

    for case in cases:
        n = case["N"]

        for _ in range(1000):
            x = random.randint(
                0,
                case["s"] + 1,
            )

            k = random.randint(
                1,
                1000,
            )

            lhs = (
                P0(case, x + k)
                - P0(case, x)
            )

            rhs = (
                k
                * (
                    2 * x
                    + k
                    - 1
                )
            )

            if lhs != rhs:
                failures += 1

    print(
        f"difference identity failures = "
        f"{failures}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: CAN SMALL k FORCE A FACTOR?")
    print("=" * 60)

    first_hit_k = []

    for case in cases:
        n = case["N"]

        case_first = None

        for k in k_values:
            found = False

            for _ in range(100):
                x = random.randint(
                    0,
                    case["s"] + 1,
                )

                g = math.gcd(
                    math.gcd(
                        P0(case, x),
                        P0(case, x + k),
                    ),
                    n,
                )

                if (
                    g == case["p"]
                    or g == case["q"]
                ):
                    found = True
                    break

            if found:
                case_first = k
                break

        if case_first is not None:
            first_hit_k.append(case_first)

    print(
        f"cases with a hit = "
        f"{len(first_hit_k)}/{cases_count}"
    )

    if first_hit_k:
        print(
            f"minimum first k = "
            f"{min(first_hit_k)}"
        )

        print(
            f"maximum first k = "
            f"{max(first_hit_k)}"
        )

        print(
            "mean first k = "
            f"{sum(first_hit_k) / len(first_hit_k):.6f}"
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
        print(f"N       = {case['N']}")
        print(f"p       = {case['p']}")
        print(f"q       = {case['q']}")
        print(f"s       = {case['s']}")
        print()

        for k in [1, 2, 3, 4, 8, 16, 32, 64]:
            best = 1

            best_x = None

            for x in range(
                0,
                min(case["s"] + 1, 1000),
            ):
                g = math.gcd(
                    math.gcd(
                        P0(case, x),
                        P0(case, x + k),
                    ),
                    case["N"],
                )

                if (
                    g != 1
                    and g != case["N"]
                ):
                    best = g
                    best_x = x
                    break

            print(
                f"k={k:2d} "
                f"first_factor_x={best_x} "
                f"gcd={best}"
            )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
