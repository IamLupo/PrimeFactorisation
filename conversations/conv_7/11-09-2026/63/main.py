import math
import random

EXPERIMENT = 81


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


def random_prime_range(low, high):
    while True:
        n = random.randrange(low, high)
        n |= 1

        if is_prime(n):
            return n


def next_prime(n):
    n = max(3, n | 1)

    while not is_prime(n):
        n += 2

    return n


def generate_case(bits, gap_ratio):
    p = random_prime_range(
        1 << (bits - 1),
        1 << bits,
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

    for m in moduli:
        table = bytearray(m)

        for x in range(m):
            table[(x * x) % m] = 1

        tables[m] = table

    return tables


def conditional_density(case, modulus, square_table):
    count = 0

    for y in range(modulus):
        d = delta(case, y)

        if square_table[d % modulus]:
            count += 1

    return count / modulus


def exact_conditional_densities(case, moduli, tables):
    densities = {}

    for m in moduli:
        densities[m] = conditional_density(
            case,
            m,
            tables[m],
        )

    return densities


def combined_density(densities):
    result = 1.0

    for value in densities.values():
        result *= value

    return result


def survives_sieve(case, y, tables):
    d = delta(case, y)

    for m, table in tables.items():
        if table[d % m] == 0:
            return False

    return True


def count_survivors(case, tables):
    target = case["y"]
    count = 0

    for y in range(target + 1):
        if survives_sieve(
            case,
            y,
            tables,
        ):
            count += 1

    return count


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
    random.seed(81001)

    bits = 20
    cases_per_group = 50

    gap_ratios = [
        0.005,
        0.01,
        0.02,
        0.05,
        0.10,
        0.20,
        0.30,
        0.40,
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
    print(f"prime bits   = {bits}")
    print(f"cases/group  = {cases_per_group}")
    print(f"moduli       = {moduli}")
    print()

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: CONDITIONAL MODULAR DENSITIES")
    print("=" * 60)

    reference_case = generate_case(
        bits,
        0.10,
    )

    tables = build_square_tables(
        moduli
    )

    densities = exact_conditional_densities(
        reference_case,
        moduli,
        tables,
    )

    for m in moduli:
        print(
            f"m={m:2d}  "
            f"conditional density="
            f"{densities[m]:.12f}"
        )

    theoretical = combined_density(
        densities
    )

    print()
    print(
        f"combined density = "
        f"{theoretical:.12f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 2
    #
    # Verify whether conditional densities are stable
    # across different semiprimes.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: DENSITY STABILITY ACROSS CASES")
    print("=" * 60)

    density_variations = {
        m: []
        for m in moduli
    }

    combined_variations = []

    test_cases = []

    for _ in range(100):
        case = generate_case(
            bits,
            random.choice(
                gap_ratios
            ),
        )

        test_cases.append(case)

        local = exact_conditional_densities(
            case,
            moduli,
            tables,
        )

        for m in moduli:
            density_variations[m].append(
                local[m]
            )

        combined_variations.append(
            combined_density(local)
        )

    for m in moduli:
        values = density_variations[m]

        print(
            f"m={m:2d}  "
            f"min={min(values):.12f}  "
            f"max={max(values):.12f}  "
            f"mean={mean(values):.12f}"
        )

    print()

    print(
        f"combined min  = "
        f"{min(combined_variations):.12f}"
    )

    print(
        f"combined max  = "
        f"{max(combined_variations):.12f}"
    )

    print(
        f"combined mean = "
        f"{mean(combined_variations):.12f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE y SURVIVAL")
    print("=" * 60)

    failures = 0

    for case in test_cases:
        if not survives_sieve(
            case,
            case["y"],
            tables,
        ):
            failures += 1

    print(
        f"failures = "
        f"{failures}/{len(test_cases)}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    #
    # Correct aggregate empirical density.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: AGGREGATE EMPIRICAL DENSITY")
    print("=" * 60)

    total_raw = 0
    total_survivors = 0
    total_pre_true_survivors = 0

    per_case_ratios = []
    normalized_pre_true = []

    for case in test_cases:
        target = case["y"]

        survivors = count_survivors(
            case,
            tables,
        )

        total_raw += target + 1
        total_survivors += survivors
        total_pre_true_survivors += (
            survivors - 1
        )

        per_case_ratios.append(
            survivors / (target + 1)
        )

        if target > 0:
            normalized_pre_true.append(
                (survivors - 1) / target
            )

    aggregate_density = (
        total_survivors
        / total_raw
    )

    aggregate_pre_true_density = (
        total_pre_true_survivors
        / (
            total_raw
            - len(test_cases)
        )
    )

    print(
        f"total raw = "
        f"{total_raw}"
    )

    print(
        f"total survivors = "
        f"{total_survivors}"
    )

    print(
        f"aggregate survivor density = "
        f"{aggregate_density:.12f}"
    )

    print(
        f"aggregate PRE-TRUE density = "
        f"{aggregate_pre_true_density:.12f}"
    )

    print(
        f"mean per-case ratio = "
        f"{mean(per_case_ratios):.12f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: PREDICTED PRE-TRUE SURVIVORS")
    print("=" * 60)

    prediction_errors = []
    observed = []
    predicted = []

    for case in test_cases:
        target = case["y"]

        local = exact_conditional_densities(
            case,
            moduli,
            tables,
        )

        rho = combined_density(
            local
        )

        survivors = count_survivors(
            case,
            tables,
        )

        actual_pre = survivors - 1

        expected_pre = rho * target

        observed.append(
            actual_pre
        )

        predicted.append(
            expected_pre
        )

        prediction_errors.append(
            actual_pre
            - expected_pre
        )

    print(
        f"mean observed = "
        f"{mean(observed):.6f}"
    )

    print(
        f"mean predicted = "
        f"{mean(predicted):.6f}"
    )

    print(
        f"mean error = "
        f"{mean(prediction_errors):.6f}"
    )

    print(
        f"mean absolute error = "
        f"{mean([abs(x) for x in prediction_errors]):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: PRE-TRUE DENSITY RATIO")
    print("=" * 60)

    ratios = []

    for case in test_cases:
        target = case["y"]

        if target == 0:
            continue

        local = exact_conditional_densities(
            case,
            moduli,
            tables,
        )

        rho = combined_density(
            local
        )

        survivors = count_survivors(
            case,
            tables,
        )

        actual_pre = survivors - 1

        if rho * target == 0:
            continue

        ratios.append(
            actual_pre
            / (rho * target)
        )

    print(
        f"count  = {len(ratios)}"
    )

    print(
        f"minimum = "
        f"{min(ratios):.6f}"
    )

    print(
        f"maximum = "
        f"{max(ratios):.6f}"
    )

    print(
        f"mean    = "
        f"{mean(ratios):.6f}"
    )

    print(
        f"median  = "
        f"{median(ratios):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: DENSITY VS GAP SIZE")
    print("=" * 60)

    for ratio_target in gap_ratios:
        cases = [
            generate_case(
                bits,
                ratio_target,
            )
            for _ in range(cases_per_group)
        ]

        raw_total = 0
        pre_total = 0

        for case in cases:
            survivors = count_survivors(
                case,
                tables,
            )

            raw_total += case["y"]
            pre_total += survivors - 1

        empirical = (
            pre_total / raw_total
            if raw_total
            else 0
        )

        print(
            f"gap_ratio={ratio_target:7.3f}  "
            f"raw={raw_total:10d}  "
            f"pre_true={pre_total:8d}  "
            f"density={empirical:.10f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 8
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 8: EXPLICIT CASES")
    print("=" * 60)

    for number, case in enumerate(
        test_cases[:10],
        start=1,
    ):
        target = case["y"]

        local = exact_conditional_densities(
            case,
            moduli,
            tables,
        )

        rho = combined_density(
            local
        )

        survivors = count_survivors(
            case,
            tables,
        )

        print()
        print(f"CASE {number}")
        print(f"N = {case['N']}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {case['s']}")
        print(f"y = {target}")
        print()
        print(
            f"conditional density = "
            f"{rho:.12f}"
        )

        print(
            f"pre-true survivors = "
            f"{survivors - 1}"
        )

        print(
            f"expected pre-true = "
            f"{rho * target:.6f}"
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
