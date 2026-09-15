import math
import random

EXPERIMENT = 61


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

        N = p * q
        s = math.isqrt(N)
        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": N,
            "s": s,
            "y": y,
        }


def discriminant(case, y):
    s = case["s"]
    N = case["N"]

    z = 2 * s + 1 + y

    return z * z - 4 * N


def build_square_tables(moduli):
    tables = {}

    for modulus in moduli:
        table = bytearray(modulus)

        for x in range(modulus):
            table[(x * x) % modulus] = 1

        tables[modulus] = table

    return tables


def survives_sieve(case, y, tables):
    d = discriminant(case, y)

    for modulus, table in tables.items():
        if table[d % modulus] == 0:
            return False

    return True


def collect_survivors(case, tables, maximum_y):
    survivors = []

    for y in range(maximum_y + 1):
        if survives_sieve(case, y, tables):
            survivors.append(y)

    return survivors


def mean(values):
    if not values:
        return 0.0

    return sum(values) / len(values)


def median(values):
    if not values:
        return 0.0

    values = sorted(values)
    n = len(values)

    if n % 2 == 1:
        return float(values[n // 2])

    return 0.5 * (
        values[n // 2 - 1] +
        values[n // 2]
    )


def percentile(values, value):
    if not values:
        return 0.0

    count = 0

    for x in values:
        if x <= value:
            count += 1

    return count / len(values)


def main():
    random.seed(61001)

    bits = 18
    cases_count = 500
    random_samples_per_case = 1000

    forward_window = 5000

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
    print(f"cases                = {cases_count}")
    print(f"prime bits           = {bits}")
    print(
        f"random samples/case = "
        f"{random_samples_per_case}"
    )
    print(f"forward window       = {forward_window}")
    print()

    cases = []

    for _ in range(cases_count):
        cases.append(generate_case(bits))

    tables = build_square_tables(moduli)

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

    print(f"failures = {failures}/{cases_count}")
    print()

    # ---------------------------------------------------------
    # TEST 2
    #
    # Collect survivors around the true y.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: COLLECT SURVIVOR STRUCTURE")
    print("=" * 60)

    case_data = []

    for case in cases:
        target = case["y"]

        survivors = collect_survivors(
            case,
            tables,
            target + forward_window,
        )

        true_index = survivors.index(target)

        case_data.append({
            "case": case,
            "survivors": survivors,
            "true_index": true_index,
        })

    print("cases collected =", len(case_data))
    print()

    # ---------------------------------------------------------
    # TEST 3
    #
    # Direct true-y gap measurements.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE y GAPS")
    print("=" * 60)

    true_previous = []
    true_next = []

    true_ratios = []

    for data in case_data:
        survivors = data["survivors"]
        i = data["true_index"]

        if i > 0:
            gap_before = survivors[i] - survivors[i - 1]
            true_previous.append(gap_before)

        if i + 1 < len(survivors):
            gap_after = survivors[i + 1] - survivors[i]
            true_next.append(gap_after)

        if i > 0 and i + 1 < len(survivors):
            gap_before = survivors[i] - survivors[i - 1]
            gap_after = survivors[i + 1] - survivors[i]

            true_ratios.append(
                gap_before / gap_after
            )

    print(
        f"previous-gap cases = "
        f"{len(true_previous)}"
    )

    print(
        f"previous-gap mean = "
        f"{mean(true_previous):.6f}"
    )

    print(
        f"previous-gap median = "
        f"{median(true_previous):.6f}"
    )

    print()

    print(
        f"next-gap cases = "
        f"{len(true_next)}"
    )

    print(
        f"next-gap mean = "
        f"{mean(true_next):.6f}"
    )

    print(
        f"next-gap median = "
        f"{median(true_next):.6f}"
    )

    print()

    print(
        f"gap-ratio mean = "
        f"{mean(true_ratios):.6f}"
    )

    print(
        f"gap-ratio median = "
        f"{median(true_ratios):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    #
    # Random survivor controls.
    #
    # We choose a random interior survivor so that both
    # neighboring gaps exist.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: RANDOM SURVIVOR CONTROL")
    print("=" * 60)

    random_previous = []
    random_next = []
    random_ratios = []

    true_vs_random_previous = []
    true_vs_random_next = []
    true_vs_random_ratio = []

    random_percentiles_previous = []
    true_percentiles_previous = []

    random_percentiles_next = []
    true_percentiles_next = []

    for data in case_data:
        survivors = data["survivors"]
        true_index = data["true_index"]

        # Interior survivors only.
        if len(survivors) < 3:
            continue

        if (
            true_index == 0
            or true_index + 1 >= len(survivors)
        ):
            continue

        true_left = (
            survivors[true_index]
            - survivors[true_index - 1]
        )

        true_right = (
            survivors[true_index + 1]
            - survivors[true_index]
        )

        # Random samples.
        local_random_left = []
        local_random_right = []

        for _ in range(random_samples_per_case):
            while True:
                random_index = random.randint(
                    1,
                    len(survivors) - 2,
                )

                # Avoid selecting the true y itself.
                if random_index != true_index:
                    break

            left = (
                survivors[random_index]
                - survivors[random_index - 1]
            )

            right = (
                survivors[random_index + 1]
                - survivors[random_index]
            )

            local_random_left.append(left)
            local_random_right.append(right)

            random_previous.append(left)
            random_next.append(right)

            random_ratios.append(
                left / right
            )

        # True percentile inside this case's
        # random-control distribution.
        true_percentiles_previous.append(
            percentile(
                local_random_left,
                true_left,
            )
        )

        true_percentiles_next.append(
            percentile(
                local_random_right,
                true_right,
            )
        )

        # Median random controls.
        random_left_median = median(
            local_random_left
        )

        random_right_median = median(
            local_random_right
        )

        true_vs_random_previous.append(
            true_left / random_left_median
        )

        true_vs_random_next.append(
            true_right / random_right_median
        )

        random_ratio_median = median(
            [
                x / y
                for x, y in zip(
                    local_random_left,
                    local_random_right,
                )
                if y != 0
            ]
        )

        true_vs_random_ratio.append(
            (true_left / true_right)
            / random_ratio_median
        )

    print(
        f"random previous-gap mean = "
        f"{mean(random_previous):.6f}"
    )

    print(
        f"random previous-gap median = "
        f"{median(random_previous):.6f}"
    )

    print()

    print(
        f"random next-gap mean = "
        f"{mean(random_next):.6f}"
    )

    print(
        f"random next-gap median = "
        f"{median(random_next):.6f}"
    )

    print()

    print(
        f"random gap-ratio mean = "
        f"{mean(random_ratios):.6f}"
    )

    print(
        f"random gap-ratio median = "
        f"{median(random_ratios):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: TRUE / RANDOM GAP RATIOS")
    print("=" * 60)

    print(
        "previous-gap ratio:"
    )

    print(
        f"  mean   = "
        f"{mean(true_vs_random_previous):.6f}"
    )

    print(
        f"  median = "
        f"{median(true_vs_random_previous):.6f}"
    )

    print()

    print(
        "next-gap ratio:"
    )

    print(
        f"  mean   = "
        f"{mean(true_vs_random_next):.6f}"
    )

    print(
        f"  median = "
        f"{median(true_vs_random_next):.6f}"
    )

    print()

    print(
        "gap-asymmetry ratio:"
    )

    print(
        f"  mean   = "
        f"{mean(true_vs_random_ratio):.6f}"
    )

    print(
        f"  median = "
        f"{median(true_vs_random_ratio):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: TRUE y PERCENTILE AGAINST RANDOM SURVIVORS")
    print("=" * 60)

    print(
        "Previous-gap percentile:"
    )

    print(
        f"  minimum = "
        f"{min(true_percentiles_previous):.6f}"
    )

    print(
        f"  maximum = "
        f"{max(true_percentiles_previous):.6f}"
    )

    print(
        f"  mean    = "
        f"{mean(true_percentiles_previous):.6f}"
    )

    print(
        f"  median  = "
        f"{median(true_percentiles_previous):.6f}"
    )

    print()

    print(
        "Next-gap percentile:"
    )

    print(
        f"  minimum = "
        f"{min(true_percentiles_next):.6f}"
    )

    print(
        f"  maximum = "
        f"{max(true_percentiles_next):.6f}"
    )

    print(
        f"  mean    = "
        f"{mean(true_percentiles_next):.6f}"
    )

    print(
        f"  median  = "
        f"{median(true_percentiles_next):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: EXPLICIT CASES")
    print("=" * 60)

    for number, data in enumerate(
        case_data[:10],
        start=1,
    ):
        case = data["case"]
        survivors = data["survivors"]
        i = data["true_index"]

        print()
        print(f"CASE {number}")
        print(f"N       = {case['N']}")
        print(f"s       = {case['s']}")
        print(f"true y  = {case['y']}")
        print(f"rank    = {i + 1}")

        if i > 0:
            print(
                "previous =",
                survivors[i - 1],
            )

            print(
                "gap before =",
                survivors[i] - survivors[i - 1],
            )
        else:
            print(
                "previous = NONE"
            )

        if i + 1 < len(survivors):
            print(
                "next     =",
                survivors[i + 1],
            )

            print(
                "gap after  =",
                survivors[i + 1] - survivors[i],
            )
        else:
            print(
                "next     = NONE"
            )

        lo = max(0, i - 3)
        hi = min(len(survivors), i + 4)

        print(
            "local survivors =",
            survivors[lo:hi],
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
