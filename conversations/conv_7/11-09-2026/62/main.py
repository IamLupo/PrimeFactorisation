import math
import random

EXPERIMENT = 80


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


def next_prime(n):
    n = max(3, n | 1)

    while not is_prime(n):
        n += 2

    return n


def generate_controlled_case(bits, gap_ratio):
    """
    Generate p and q with approximately

        q = p * (1 + gap_ratio).

    The actual gap is determined by the nearby prime q.
    """

    low = 1 << (bits - 1)
    high = 1 << bits

    p = random_prime_range(
        low,
        high,
    )

    target_q = int(
        p * (1.0 + gap_ratio)
    )

    q = next_prime(target_q)

    if q <= p:
        q = next_prime(p + 2)

    n = p * q
    s = math.isqrt(n)

    y = p + q - 2 * s - 1
    delta_gap = q - p

    return {
        "p": p,
        "q": q,
        "N": n,
        "s": s,
        "y": y,
        "gap": delta_gap,
        "gap_ratio": delta_gap / p,
    }


def random_prime_range(low, high):
    while True:
        n = random.randrange(
            low,
            high,
        )

        n |= 1

        if is_prime(n):
            return n


def generate_case(bits, gap_ratio):
    return generate_controlled_case(
        bits,
        gap_ratio,
    )


def delta(case, y):
    s = case["s"]
    n = case["N"]

    z = 2 * s + 1 + y

    return z * z - 4 * n


def build_square_tables(moduli):
    tables = {}

    for modulus in moduli:
        table = bytearray(modulus)

        for x in range(modulus):
            table[
                (x * x) % modulus
            ] = 1

        tables[modulus] = table

    return tables


def survives_sieve(case, y, tables):
    d = delta(case, y)

    for modulus, table in tables.items():
        if table[d % modulus] == 0:
            return False

    return True


def count_search(case, tables):
    target = case["y"]

    raw = target + 1
    survivors = 0
    square_tests = 0
    found = None

    for y in range(target + 1):
        if not survives_sieve(
            case,
            y,
            tables,
        ):
            continue

        survivors += 1
        square_tests += 1

        d = delta(case, y)

        root = math.isqrt(d)

        if root * root == d:
            found = y
            break

    if found != target:
        raise RuntimeError(
            "True y was not the first exact square"
        )

    return {
        "raw": raw,
        "survivors": survivors,
        "square_tests": square_tests,
    }


def theoretical_density(moduli):
    density = 1.0

    for modulus in moduli:
        residues = set()

        for x in range(modulus):
            residues.add(
                (x * x) % modulus
            )

        density *= (
            len(residues) / modulus
        )

    return density


def mean(values):
    if not values:
        return 0.0

    return sum(values) / len(values)


def median(values):
    if not values:
        return 0.0

    values = sorted(values)
    n = len(values)

    if n % 2:
        return float(values[n // 2])

    return (
        values[n // 2 - 1]
        + values[n // 2]
    ) / 2.0


def main():
    random.seed(80001)

    bits = 18
    cases_per_group = 50

    gap_ratios = [
        0.0005,
        0.001,
        0.002,
        0.005,
        0.01,
        0.02,
        0.05,
        0.10,
        0.20,
    ]

    moduli = [
        16,
        3,
        5,
        7,
        11,
        13,
        17,
        19,
        23,
    ]

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"prime bits        = {bits}")
    print(f"cases/group       = {cases_per_group}")
    print(f"moduli            = {moduli}")
    print()

    tables = build_square_tables(
        moduli
    )

    rho = theoretical_density(
        moduli
    )

    print(
        f"theoretical residue density = "
        f"{rho:.12f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: TRUE y SURVIVAL")
    print("=" * 60)

    total_failures = 0

    for gap_ratio in gap_ratios:
        failures = 0

        for _ in range(cases_per_group):
            case = generate_case(
                bits,
                gap_ratio,
            )

            if not survives_sieve(
                case,
                case["y"],
                tables,
            ):
                failures += 1

        total_failures += failures

        print(
            f"gap_ratio={gap_ratio:8.4f} "
            f"failures={failures}/{cases_per_group}"
        )

    print(
        f"total failures = {total_failures}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: SEARCH SIZE VS FACTOR GAP")
    print("=" * 60)

    print(
        "columns:"
        " ratio, mean gap, mean y,"
        " mean survivors, mean raw,"
        " reduction"
    )

    for gap_ratio in gap_ratios:
        gaps = []
        ys = []
        survivors = []
        raws = []

        for _ in range(
            cases_per_group
        ):
            case = generate_case(
                bits,
                gap_ratio,
            )

            result = count_search(
                case,
                tables,
            )

            gaps.append(
                case["gap"]
            )

            ys.append(
                case["y"]
            )

            survivors.append(
                result["survivors"]
            )

            raws.append(
                result["raw"]
            )

        mean_gap = mean(gaps)
        mean_y = mean(ys)
        mean_survivors = mean(
            survivors
        )
        mean_raw = mean(raws)

        reduction = (
            1.0
            - mean_survivors
            / mean_raw
        )

        print(
            f"{gap_ratio:8.4f} "
            f"{mean_gap:12.2f} "
            f"{mean_y:12.2f} "
            f"{mean_survivors:12.2f} "
            f"{mean_raw:12.2f} "
            f"{reduction:12.8f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: SIEVE SURVIVORS / RAW")
    print("=" * 60)

    ratios = []

    for gap_ratio in gap_ratios:
        for _ in range(
            cases_per_group
        ):
            case = generate_case(
                bits,
                gap_ratio,
            )

            result = count_search(
                case,
                tables,
            )

            ratios.append(
                result["survivors"]
                / result["raw"]
            )

    empirical_density = mean(
        ratios
    )

    print(
        f"empirical density = "
        f"{empirical_density:.12f}"
    )

    print(
        f"theoretical density = "
        f"{rho:.12f}"
    )

    print(
        f"ratio empirical/theory = "
        f"{empirical_density / rho:.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: NORMALIZED SIEVE WORK")
    print("=" * 60)

    print(
        "For a purely density-based sieve,"
        " survivors/y should stay approximately"
        " constant."
    )

    normalized = []

    for gap_ratio in gap_ratios:
        for _ in range(
            cases_per_group
        ):
            case = generate_case(
                bits,
                gap_ratio,
            )

            result = count_search(
                case,
                tables,
            )

            normalized.append(
                result["survivors"]
                / (
                    case["gap"] ** 2
                    / (
                        4
                        * math.isqrt(
                            case["N"]
                        )
                    )
                    + 1
                )
            )

    print(
        f"minimum normalized work = "
        f"{min(normalized):.12e}"
    )

    print(
        f"maximum normalized work = "
        f"{max(normalized):.12e}"
    )

    print(
        f"mean normalized work = "
        f"{mean(normalized):.12e}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: SQUARE-TEST REDUCTION")
    print("=" * 60)

    for gap_ratio in gap_ratios:
        square_counts = []
        raw_counts = []

        for _ in range(
            cases_per_group
        ):
            case = generate_case(
                bits,
                gap_ratio,
            )

            result = count_search(
                case,
                tables,
            )

            square_counts.append(
                result["square_tests"]
            )

            raw_counts.append(
                result["raw"]
            )

        reduction = (
            1.0
            - mean(square_counts)
            / mean(raw_counts)
        )

        print(
            f"gap_ratio={gap_ratio:8.4f} "
            f"square-test reduction="
            f"{reduction:.8f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: EXTREME GAP EXAMPLES")
    print("=" * 60)

    for gap_ratio in [
        gap_ratios[0],
        gap_ratios[3],
        gap_ratios[6],
        gap_ratios[-1],
    ]:
        case = generate_case(
            bits,
            gap_ratio,
        )

        result = count_search(
            case,
            tables,
        )

        print()
        print(
            f"gap_ratio target = "
            f"{gap_ratio}"
        )

        print(
            f"p = {case['p']}"
        )

        print(
            f"q = {case['q']}"
        )

        print(
            f"actual gap = {case['gap']}"
        )

        print(
            f"actual ratio = "
            f"{case['gap_ratio']:.8f}"
        )

        print(
            f"y = {case['y']}"
        )

        print(
            f"raw candidates = "
            f"{result['raw']}"
        )

        print(
            f"sieve survivors = "
            f"{result['survivors']}"
        )

        print(
            f"reduction = "
            f"{1.0 - result['survivors'] / result['raw']:.8f}"
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
