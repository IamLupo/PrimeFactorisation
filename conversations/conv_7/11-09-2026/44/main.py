import math
import random

EXPERIMENT = 62


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


def percentile_rank(values, value):
    if not values:
        return 0.0

    count = 0

    for x in values:
        if x <= value:
            count += 1

    return count / len(values)


def main():
    random.seed(62001)

    bits = 18
    cases_count = 500

    # This is enough to obtain many survivors after y.
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
    print(f"cases          = {cases_count}")
    print(f"prime bits     = {bits}")
    print(f"forward window = {forward_window}")
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

        case_data.append(
            (case, survivors, true_index)
        )

    print(f"cases collected = {len(case_data)}")
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE y GAPS")
    print("=" * 60)

    true_left = []
    true_right = []
    true_ratio = []

    for case, survivors, i in case_data:

        if i > 0:
            left = survivors[i] - survivors[i - 1]
            true_left.append(left)

        if i + 1 < len(survivors):
            right = survivors[i + 1] - survivors[i]
            true_right.append(right)

        if i > 0 and i + 1 < len(survivors):
            left = survivors[i] - survivors[i - 1]
            right = survivors[i + 1] - survivors[i]

            true_ratio.append(left / right)

    print(f"previous-gap cases = {len(true_left)}")
    print(f"previous-gap mean  = {mean(true_left):.6f}")
    print(f"previous-gap median= {median(true_left):.6f}")
    print()

    print(f"next-gap cases = {len(true_right)}")
    print(f"next-gap mean  = {mean(true_right):.6f}")
    print(f"next-gap median= {median(true_right):.6f}")
    print()

    print(f"gap-ratio mean   = {mean(true_ratio):.6f}")
    print(f"gap-ratio median = {median(true_ratio):.6f}")
    print()

    # ---------------------------------------------------------
    # TEST 4
    #
    # ALL OTHER SURVIVORS
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: ALL-OTHER-SURVIVOR CONTROL")
    print("=" * 60)

    control_left = []
    control_right = []
    control_ratio = []

    true_left_percentiles = []
    true_right_percentiles = []
    true_ratio_percentiles = []

    cases_with_left_control = 0
    cases_with_right_control = 0
    cases_with_ratio_control = 0

    for case, survivors, true_index in case_data:

        if true_index == 0:
            continue

        if true_index + 1 >= len(survivors):
            continue

        # True gaps.
        left_true = (
            survivors[true_index]
            - survivors[true_index - 1]
        )

        right_true = (
            survivors[true_index + 1]
            - survivors[true_index]
        )

        ratio_true = left_true / right_true

        # Every other interior survivor is a control.
        local_left = []
        local_right = []
        local_ratio = []

        for i in range(1, len(survivors) - 1):

            if i == true_index:
                continue

            left = survivors[i] - survivors[i - 1]
            right = survivors[i + 1] - survivors[i]

            ratio = left / right

            local_left.append(left)
            local_right.append(right)
            local_ratio.append(ratio)

        if local_left:
            cases_with_left_control += 1

            control_left.extend(local_left)

            true_left_percentiles.append(
                percentile_rank(
                    local_left,
                    left_true,
                )
            )

        if local_right:
            cases_with_right_control += 1

            control_right.extend(local_right)

            true_right_percentiles.append(
                percentile_rank(
                    local_right,
                    right_true,
                )
            )

        if local_ratio:
            cases_with_ratio_control += 1

            control_ratio.extend(local_ratio)

            true_ratio_percentiles.append(
                percentile_rank(
                    local_ratio,
                    ratio_true,
                )
            )

    print(
        f"left-control cases = "
        f"{cases_with_left_control}"
    )

    print(
        f"control left mean = "
        f"{mean(control_left):.6f}"
    )

    print(
        f"control left median = "
        f"{median(control_left):.6f}"
    )

    print()

    print(
        f"right-control cases = "
        f"{cases_with_right_control}"
    )

    print(
        f"control right mean = "
        f"{mean(control_right):.6f}"
    )

    print(
        f"control right median = "
        f"{median(control_right):.6f}"
    )

    print()

    print(
        f"ratio-control cases = "
        f"{cases_with_ratio_control}"
    )

    print(
        f"control ratio mean = "
        f"{mean(control_ratio):.6f}"
    )

    print(
        f"control ratio median = "
        f"{median(control_ratio):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: TRUE y PERCENTILES")
    print("=" * 60)

    print("Previous gap:")
    print(
        f"  min    = "
        f"{min(true_left_percentiles):.6f}"
    )
    print(
        f"  max    = "
        f"{max(true_left_percentiles):.6f}"
    )
    print(
        f"  mean   = "
        f"{mean(true_left_percentiles):.6f}"
    )
    print(
        f"  median = "
        f"{median(true_left_percentiles):.6f}"
    )

    print()

    print("Next gap:")
    print(
        f"  min    = "
        f"{min(true_right_percentiles):.6f}"
    )
    print(
        f"  max    = "
        f"{max(true_right_percentiles):.6f}"
    )
    print(
        f"  mean   = "
        f"{mean(true_right_percentiles):.6f}"
    )
    print(
        f"  median = "
        f"{median(true_right_percentiles):.6f}"
    )

    print()

    print("Gap ratio:")
    print(
        f"  min    = "
        f"{min(true_ratio_percentiles):.6f}"
    )
    print(
        f"  max    = "
        f"{max(true_ratio_percentiles):.6f}"
    )
    print(
        f"  mean   = "
        f"{mean(true_ratio_percentiles):.6f}"
    )
    print(
        f"  median = "
        f"{median(true_ratio_percentiles):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: TRUE y VS CONTROL")
    print("=" * 60)

    print(
        "A value near 0.5 means the true y looks "
        "like an ordinary survivor."
    )

    print()

    print(
        "Expected neutral percentile:"
    )

    print("  approximately 0.50")
    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: EXPLICIT CASES")
    print("=" * 60)

    for number, (case, survivors, i) in enumerate(
        case_data[:10],
        start=1,
    ):
        print()
        print(f"CASE {number}")
        print(f"N       = {case['N']}")
        print(f"s       = {case['s']}")
        print(f"true y  = {case['y']}")
        print(f"rank    = {i + 1}")

        if i > 0:
            print(
                "left gap  =",
                survivors[i] - survivors[i - 1],
            )
        else:
            print("left gap  = NONE")

        if i + 1 < len(survivors):
            print(
                "right gap =",
                survivors[i + 1] - survivors[i],
            )
        else:
            print("right gap = NONE")

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
