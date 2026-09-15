import math
import random

EXPERIMENT = 60


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
            table[(x * x) % modulus] = 1

        tables[modulus] = table

    return tables


def survives_sieve(case, y, tables):
    d = delta(case, y)

    for modulus, table in tables.items():
        if table[d % modulus] == 0:
            return False

    return True


def collect_survivors(case, tables, stop_y):
    survivors = []

    for y in range(stop_y + 1):
        if survives_sieve(case, y, tables):
            survivors.append(y)

    return survivors


def percentile_rank(value, values):
    if not values:
        return 0.0

    values = sorted(values)

    count = 0

    for x in values:
        if x <= value:
            count += 1

    return count / len(values)


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

    return 0.5 * (values[n // 2 - 1] + values[n // 2])


def standard_deviation(values):
    if len(values) < 2:
        return 0.0

    m = mean(values)

    return math.sqrt(
        sum((x - m) ** 2 for x in values)
        / (len(values) - 1)
    )


def main():
    random.seed(60001)

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

    # We need a survivor after the true y.
    # 5000 is deliberately much larger than the usual sieve gap.
    forward_window = 5000

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"cases           = {cases_count}")
    print(f"prime bits      = {bits}")
    print(f"forward window  = {forward_window}")
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
        if not survives_sieve(case, case["y"], tables):
            failures += 1

    print(f"failures = {failures}/{cases_count}")
    print()

    # ---------------------------------------------------------
    # MAIN COLLECTION
    # ---------------------------------------------------------

    all_previous_gaps = []
    all_next_gaps = []
    all_ratio_values = []

    true_previous_gaps = []
    true_next_gaps = []

    previous_percentiles = []
    next_percentiles = []

    cases_without_next = 0

    for case in cases:
        target = case["y"]

        stop_y = target + forward_window

        survivors = collect_survivors(
            case,
            tables,
            stop_y,
        )

        target_index = survivors.index(target)

        # A predecessor is usually available because y is
        # normally not the first survivor.
        if target_index > 0:
            prev_y = survivors[target_index - 1]

            true_gap_before = target - prev_y
            true_previous_gaps.append(true_gap_before)

            # Ordinary gaps before the true y.
            ordinary_gaps = []

            for i in range(1, target_index + 1):
                ordinary_gaps.append(
                    survivors[i] - survivors[i - 1]
                )

            all_previous_gaps.extend(ordinary_gaps)

            previous_percentiles.append(
                percentile_rank(
                    true_gap_before,
                    ordinary_gaps,
                )
            )

        # Find next survivor after true y.
        if target_index + 1 < len(survivors):
            next_y = survivors[target_index + 1]

            true_gap_after = next_y - target
            true_next_gaps.append(true_gap_after)

            ordinary_next_gaps = []

            for i in range(target_index + 1, len(survivors) - 1):
                ordinary_next_gaps.append(
                    survivors[i + 1] - survivors[i]
                )

            all_next_gaps.extend(ordinary_next_gaps)

        else:
            cases_without_next += 1

        # A symmetric local quantity.
        if (
            target_index > 0
            and target_index + 1 < len(survivors)
        ):
            prev_y = survivors[target_index - 1]
            next_y = survivors[target_index + 1]

            gap_left = target - prev_y
            gap_right = next_y - target

            ratio = (
                gap_left / gap_right
                if gap_right != 0
                else float("inf")
            )

            all_ratio_values.append(ratio)

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: PREVIOUS-GAP DISTRIBUTION")
    print("=" * 60)

    print(
        f"true gap count = {len(true_previous_gaps)}"
    )

    print(
        f"ordinary gap count = {len(all_previous_gaps)}"
    )

    print(
        f"true minimum = {min(true_previous_gaps)}"
    )

    print(
        f"true maximum = {max(true_previous_gaps)}"
    )

    print(
        f"true mean    = {mean(true_previous_gaps):.6f}"
    )

    print(
        f"true median  = {median(true_previous_gaps):.6f}"
    )

    print()

    print(
        f"ordinary minimum = {min(all_previous_gaps)}"
    )

    print(
        f"ordinary maximum = {max(all_previous_gaps)}"
    )

    print(
        f"ordinary mean    = {mean(all_previous_gaps):.6f}"
    )

    print(
        f"ordinary median  = {median(all_previous_gaps):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE PREVIOUS-GAP PERCENTILE")
    print("=" * 60)

    print(
        "Each percentile is the fraction of prior "
        "sieve gaps <= the true-y gap."
    )

    print(
        f"minimum percentile = "
        f"{min(previous_percentiles):.6f}"
    )

    print(
        f"maximum percentile = "
        f"{max(previous_percentiles):.6f}"
    )

    print(
        f"mean percentile    = "
        f"{mean(previous_percentiles):.6f}"
    )

    print(
        f"median percentile  = "
        f"{median(previous_percentiles):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: NEXT-GAP DISTRIBUTION")
    print("=" * 60)

    print(
        f"cases without next survivor = "
        f"{cases_without_next}"
    )

    if true_next_gaps:
        print(
            f"true minimum = {min(true_next_gaps)}"
        )

        print(
            f"true maximum = {max(true_next_gaps)}"
        )

        print(
            f"true mean    = {mean(true_next_gaps):.6f}"
        )

        print(
            f"true median  = {median(true_next_gaps):.6f}"
        )

        print()

        print(
            f"ordinary minimum = {min(all_next_gaps)}"
        )

        print(
            f"ordinary maximum = {max(all_next_gaps)}"
        )

        print(
            f"ordinary mean    = {mean(all_next_gaps):.6f}"
        )

        print(
            f"ordinary median  = {median(all_next_gaps):.6f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: GAP RATIO")
    print("=" * 60)

    if all_ratio_values:
        finite_ratios = [
            x for x in all_ratio_values
            if math.isfinite(x)
        ]

        print(
            f"count  = {len(finite_ratios)}"
        )

        print(
            f"minimum = {min(finite_ratios):.6f}"
        )

        print(
            f"maximum = {max(finite_ratios):.6f}"
        )

        print(
            f"mean    = {mean(finite_ratios):.6f}"
        )

        print(
            f"median  = {median(finite_ratios):.6f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: EXPLICIT CASES")
    print("=" * 60)

    for index, case in enumerate(cases[:10], start=1):
        target = case["y"]

        stop_y = target + forward_window

        survivors = collect_survivors(
            case,
            tables,
            stop_y,
        )

        target_index = survivors.index(target)

        print()
        print(f"CASE {index}")
        print(f"N       = {case['N']}")
        print(f"s       = {case['s']}")
        print(f"true y  = {target}")

        if target_index > 0:
            previous = survivors[target_index - 1]
            gap_before = target - previous

            print(f"previous survivor = {previous}")
            print(f"gap before        = {gap_before}")
        else:
            print("previous survivor = NONE")
            print("gap before        = NONE")

        if target_index + 1 < len(survivors):
            following = survivors[target_index + 1]
            gap_after = following - target

            print(f"next survivor     = {following}")
            print(f"gap after         = {gap_after}")
        else:
            print("next survivor     = NONE")
            print("gap after         = NONE")

        start = max(0, target_index - 3)
        end = min(len(survivors), target_index + 4)

        print(
            "local survivors   =",
            survivors[start:end],
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
