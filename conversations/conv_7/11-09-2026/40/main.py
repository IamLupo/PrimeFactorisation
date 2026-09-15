import math
import random

EXPERIMENT = 59


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


def is_quadratic_residue_mod(value, modulus, square_table):
    return square_table[value % modulus] != 0


def build_square_table(modulus):
    table = bytearray(modulus)

    for x in range(modulus):
        table[(x * x) % modulus] = 1

    return table


def build_square_tables(moduli):
    tables = {}

    for modulus in moduli:
        tables[modulus] = build_square_table(modulus)

    return tables


def survives_sieve(case, y, tables):
    d = delta(case, y)

    for modulus, table in tables.items():
        if table[d % modulus] == 0:
            return False

    return True


def collect_survivors(case, tables):
    survivors = []
    target = case["y"]

    for y in range(target + 1):
        if survives_sieve(case, y, tables):
            survivors.append(y)

    return survivors


def empirical_density(case, tables):
    target = case["y"]
    count = 0

    for y in range(target + 1):
        if survives_sieve(case, y, tables):
            count += 1

    return count / (target + 1)


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

    return 0.5 * (values[n // 2 - 1] + values[n // 2])


def standard_deviation(values):
    if len(values) < 2:
        return 0.0

    m = mean(values)

    return math.sqrt(
        sum((x - m) ** 2 for x in values) / (len(values) - 1)
    )


def percentile(values, fraction):
    if not values:
        return 0.0

    values = sorted(values)

    index = fraction * (len(values) - 1)

    lo = math.floor(index)
    hi = math.ceil(index)

    if lo == hi:
        return float(values[lo])

    weight = index - lo

    return (
        values[lo] * (1.0 - weight)
        + values[hi] * weight
    )


def expected_binomial_mean(y, density):
    return (y + 1) * density


def expected_binomial_variance(y, density):
    n = y + 1
    return n * density * (1.0 - density)


def main():
    random.seed(59001)

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
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: GLOBAL SIEVE DENSITY")
    print("=" * 60)

    all_survivors = 0
    all_candidates = 0

    per_case_densities = []

    for case in cases:
        target = case["y"]

        survivors = 0

        for y in range(target + 1):
            if survives_sieve(case, y, tables):
                survivors += 1

        all_survivors += survivors
        all_candidates += target + 1

        per_case_densities.append(
            survivors / (target + 1)
        )

    global_density = all_survivors / all_candidates

    print(f"raw candidates    = {all_candidates}")
    print(f"sieved candidates = {all_survivors}")
    print(f"global density    = {global_density:.10f}")
    print(f"minimum density   = {min(per_case_densities):.10f}")
    print(f"maximum density   = {max(per_case_densities):.10f}")
    print(f"mean density      = {mean(per_case_densities):.10f}")
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE y RANK")
    print("=" * 60)

    ranks = []
    normalized_ranks = []
    local_density_ratios = []

    for case in cases:
        target = case["y"]

        rank = 0

        for y in range(target + 1):
            if survives_sieve(case, y, tables):
                rank += 1

        ranks.append(rank)

        expected = (target + 1) * global_density

        if expected > 0:
            normalized_ranks.append(rank / expected)
        else:
            normalized_ranks.append(0.0)

        local_density = per_case_densities[len(ranks) - 1]

        if local_density > 0:
            local_density_ratios.append(
                rank / ((target + 1) * local_density)
            )

    print(f"minimum rank = {min(ranks)}")
    print(f"maximum rank = {max(ranks)}")
    print(f"mean rank    = {mean(ranks):.6f}")
    print(f"median rank  = {median(ranks):.6f}")
    print()
    print("Normalized by GLOBAL density:")
    print(f"minimum Q    = {min(normalized_ranks):.6f}")
    print(f"maximum Q    = {max(normalized_ranks):.6f}")
    print(f"mean Q       = {mean(normalized_ranks):.6f}")
    print(f"median Q     = {median(normalized_ranks):.6f}")
    print()

    print("Normalized by LOCAL density:")
    print(f"minimum Q    = {min(local_density_ratios):.6f}")
    print(f"maximum Q    = {max(local_density_ratios):.6f}")
    print(f"mean Q       = {mean(local_density_ratios):.6f}")
    print(f"median Q     = {median(local_density_ratios):.6f}")
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: RANK PERCENTILES")
    print("=" * 60)

    print(
        "Q percentiles:"
    )
    print(
        f"  1%  = {percentile(normalized_ranks, 0.01):.6f}"
    )
    print(
        f"  5%  = {percentile(normalized_ranks, 0.05):.6f}"
    )
    print(
        f" 25%  = {percentile(normalized_ranks, 0.25):.6f}"
    )
    print(
        f" 50%  = {percentile(normalized_ranks, 0.50):.6f}"
    )
    print(
        f" 75%  = {percentile(normalized_ranks, 0.75):.6f}"
    )
    print(
        f" 95%  = {percentile(normalized_ranks, 0.95):.6f}"
    )
    print(
        f" 99%  = {percentile(normalized_ranks, 0.99):.6f}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: RANK VS EXPECTED BINOMIAL")
    print("=" * 60)

    mean_error = []
    variance_error = []

    for i, case in enumerate(cases):
        target = case["y"]
        observed = ranks[i]

        expected_mean = expected_binomial_mean(
            target,
            global_density
        )

        expected_variance = expected_binomial_variance(
            target,
            global_density
        )

        mean_error.append(
            observed - expected_mean
        )

        # Compare squared deviation with expected variance.
        variance_error.append(
            (observed - expected_mean) ** 2
            - expected_variance
        )

    print(
        "mean(observed - expected) = "
        f"{mean(mean_error):.6f}"
    )

    print(
        "mean((observed-expected)^2 - variance) = "
        f"{mean(variance_error):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: TRUE y ORDER BIAS")
    print("=" * 60)

    bins = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    for q in normalized_ranks:
        if q <= 0:
            index = 0
        else:
            index = min(9, int(q))

        bins[index] += 1

    print("Q intervals:")
    print("0.0 <= Q < 1.0 :", bins[0])
    print("1.0 <= Q < 2.0 :", bins[1])
    print("2.0 <= Q < 3.0 :", bins[2])
    print("3.0 <= Q < 4.0 :", bins[3])
    print("4.0 <= Q < 5.0 :", bins[4])
    print("5.0 <= Q < 6.0 :", bins[5])
    print("6.0 <= Q < 7.0 :", bins[6])
    print("7.0 <= Q < 8.0 :", bins[7])
    print("8.0 <= Q < 9.0 :", bins[8])
    print("Q >= 9.0       :", bins[9])
    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: EXPLICIT CASES")
    print("=" * 60)

    for i, case in enumerate(cases[:10], start=1):
        target = case["y"]

        survivors = []

        for y in range(target + 1):
            if survives_sieve(case, y, tables):
                survivors.append(y)

        rank = survivors.index(target) + 1

        expected = (target + 1) * global_density
        q = rank / expected if expected > 0 else 0.0

        print()
        print(f"CASE {i}")
        print(f"N       = {case['N']}")
        print(f"s       = {case['s']}")
        print(f"true y  = {target}")
        print(f"rank    = {rank}")
        print(f"expected= {expected:.4f}")
        print(f"Q       = {q:.6f}")
        print(f"survivors before/at true y = {len(survivors)}")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
