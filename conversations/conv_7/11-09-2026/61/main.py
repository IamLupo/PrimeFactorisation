import math
import random


EXPERIMENT = 79


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


def square_gap(value):
    root = math.isqrt(value)

    if root * root == value:
        return 0, root

    upper = root + 1

    return (
        upper * upper - value,
        upper,
    )


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


def percentile(values, value):
    if not values:
        return 0.0

    count = 0

    for x in values:
        if x <= value:
            count += 1

    return count / len(values)


def main():
    random.seed(79001)

    bits = 18
    cases_count = 500

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
    print(f"cases      = {cases_count}")
    print(f"prime bits = {bits}")
    print(f"moduli     = {moduli}")
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    tables = build_square_tables(
        moduli
    )

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: TRUE y SURVIVAL")
    print("=" * 60)

    failures = 0

    for case in cases:
        if not survives_sieve(
            case,
            case["y"],
            tables,
        ):
            failures += 1

    print(
        f"failures = {failures}/{cases_count}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 2
    #
    # Collect square gaps of all surviving candidates.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: GLOBAL SQUARE-GAP DISTRIBUTION")
    print("=" * 60)

    all_gaps = []
    all_square_hits = 0
    total_survivors = 0

    for case in cases:
        target = case["y"]

        for y in range(target + 1):
            if not survives_sieve(
                case,
                y,
                tables,
            ):
                continue

            total_survivors += 1

            d = delta(
                case,
                y,
            )

            gap, _ = square_gap(d)

            all_gaps.append(gap)

            if gap == 0:
                all_square_hits += 1

    print(
        f"total survivors = "
        f"{total_survivors}"
    )

    print(
        f"exact square hits = "
        f"{all_square_hits}"
    )

    print(
        f"minimum gap = "
        f"{min(all_gaps)}"
    )

    print(
        f"maximum gap = "
        f"{max(all_gaps)}"
    )

    print(
        f"mean gap = "
        f"{mean(all_gaps):.6f}"
    )

    print(
        f"median gap = "
        f"{median(all_gaps):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    #
    # Rank the true y's square gap among all survivors
    # before it. The true gap is always zero, but this tells
    # us how often an earlier square survivor exists.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: EARLIER EXACT-SQUARE SURVIVORS")
    print("=" * 60)

    earlier_square_cases = 0
    earlier_square_total = 0

    true_ranks = []

    for case in cases:
        target = case["y"]

        squares_before = 0
        survivors_before = 0

        for y in range(target):
            if not survives_sieve(
                case,
                y,
                tables,
            ):
                continue

            survivors_before += 1

            d = delta(
                case,
                y,
            )

            gap, _ = square_gap(d)

            if gap == 0:
                squares_before += 1

        if squares_before > 0:
            earlier_square_cases += 1
            earlier_square_total += (
                squares_before
            )

        true_ranks.append(
            survivors_before + 1
        )

    print(
        f"cases with earlier square = "
        f"{earlier_square_cases}/{cases_count}"
    )

    print(
        f"total earlier squares = "
        f"{earlier_square_total}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    #
    # Instead of asking whether the true y is square
    # (which is tautological), examine candidates in a
    # fixed neighborhood before y.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: PRE-TRUE SQUARE-GAP TREND")
    print("=" * 60)

    windows = [
        1,
        2,
        3,
        5,
        10,
        20,
        50,
    ]

    for window in windows:
        samples = []
        case_count = 0

        for case in cases:
            target = case["y"]

            start = max(
                0,
                target - window,
            )

            local = []

            for y in range(
                start,
                target,
            ):
                if not survives_sieve(
                    case,
                    y,
                    tables,
                ):
                    continue

                d = delta(
                    case,
                    y,
                )

                gap, _ = square_gap(d)

                local.append(gap)

            if local:
                case_count += 1
                samples.extend(local)

        if samples:
            print(
                f"window={window:2d}  "
                f"cases={case_count:3d}  "
                f"mean={mean(samples):.6f}  "
                f"median={median(samples):.6f}  "
                f"min={min(samples)}"
            )

    print()

    # ---------------------------------------------------------
    # TEST 5
    #
    # Compare square gaps of the true case's survivors
    # against survivors from an unrelated random y-position
    # in the same range.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: RANDOM-POSITION CONTROL")
    print("=" * 60)

    random_control_gaps = []

    for case in cases:
        target = case["y"]

        if target < 100:
            continue

        random_limit = target

        for _ in range(20):
            y = random.randint(
                0,
                random_limit,
            )

            if not survives_sieve(
                case,
                y,
                tables,
            ):
                continue

            d = delta(
                case,
                y,
            )

            gap, _ = square_gap(d)

            random_control_gaps.append(
                gap
            )

    print(
        f"random control samples = "
        f"{len(random_control_gaps)}"
    )

    print(
        f"mean = "
        f"{mean(random_control_gaps):.6f}"
    )

    print(
        f"median = "
        f"{median(random_control_gaps):.6f}"
    )

    print(
        f"minimum = "
        f"{min(random_control_gaps)}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    #
    # Examine how small the gap becomes before the true y.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: MINIMUM PRE-TRUE SQUARE GAP")
    print("=" * 60)

    minima = []

    for case in cases:
        target = case["y"]

        best = None

        for y in range(target):
            if not survives_sieve(
                case,
                y,
                tables,
            ):
                continue

            d = delta(
                case,
                y,
            )

            gap, _ = square_gap(d)

            if best is None or gap < best:
                best = gap

        if best is not None:
            minima.append(best)

    print(
        f"cases = {len(minima)}"
    )

    print(
        f"minimum = "
        f"{min(minima)}"
    )

    print(
        f"maximum = "
        f"{max(minima)}"
    )

    print(
        f"mean = "
        f"{mean(minima):.6f}"
    )

    print(
        f"median = "
        f"{median(minima):.6f}"
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
        target = case["y"]

        survivors = []

        for y in range(target + 1):
            if not survives_sieve(
                case,
                y,
                tables,
            ):
                continue

            d = delta(
                case,
                y,
            )

            gap, root = square_gap(d)

            survivors.append(
                (y, gap, root)
            )

        print()
        print(f"CASE {number}")
        print(f"N = {case['N']}")
        print(f"s = {case['s']}")
        print(f"true y = {target}")
        print(
            "survivor count =",
            len(survivors),
        )

        print(
            "last survivors =",
            survivors[-10:],
        )

        earlier = [
            item
            for item in survivors
            if item[0] < target
            and item[1] == 0
        ]

        print(
            "earlier square survivors =",
            earlier,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
