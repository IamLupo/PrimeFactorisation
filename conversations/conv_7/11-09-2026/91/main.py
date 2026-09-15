import math
import random
from sympy import isprime


EXPERIMENT = 112

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384

PROBE_BUDGETS = (
    1,
    2,
    4,
    8,
    16,
    32,
    64,
)


def gcd_state(g, n, p, q):
    if g == 1:
        return "1"

    if g == n:
        return "n"

    if g == p:
        return "p"

    if g == q:
        return "q"

    return "other"


def value_gcd(n, k, z):
    """
    For d = 2z - 1:

        v_k = -C(k-2,d)

    whenever d < k-1.

    Otherwise v_k = 0 and gcd(v_k,n) = n.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d)
    value %= n
    value = (-value) % n

    return math.gcd(value, n)


def generate_semiprime():
    while True:
        p = random.randrange(
            P_MIN,
            P_MAX + 1,
        )

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, P_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        if p == q:
            continue

        return p, q


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def random_sequence(max_z, count, seed):
    """
    Uniform random unique positions.
    """

    count = min(
        count,
        max_z,
    )

    rng = random.Random(seed)

    return rng.sample(
        range(1, max_z + 1),
        count,
    )


def evenly_spaced_sequence(max_z, count):
    """
    Deterministic approximately-uniform spacing.
    """

    count = min(
        count,
        max_z,
    )

    if count == 1:
        return [1]

    result = []

    for i in range(count):
        numerator = i * (max_z - 1)
        denominator = count - 1

        z = 1 + (
            numerator // denominator
        )

        result.append(z)

    # Remove duplicates caused by very small ranges.
    result = list(dict.fromkeys(result))

    return result


def bit_reverse(value, bits):
    result = 0

    for _ in range(bits):
        result <<= 1
        result |= value & 1
        value >>= 1

    return result


def van_der_corput_sequence(max_z, count):
    """
    Generate positions using bit-reversal ordering.

    The sequence spreads early samples throughout the
    complete interval rather than clustering them.
    """

    count = min(
        count,
        max_z,
    )

    bits = max(
        1,
        (count - 1).bit_length(),
    )

    result = []
    used = set()

    index = 0

    while len(result) < count:
        reversed_value = bit_reverse(
            index,
            bits,
        )

        # Map [0, 2^bits-1] into [1, max_z].
        z = (
            reversed_value * max_z
            // (1 << bits)
        ) + 1

        if z < 1:
            z = 1

        if z > max_z:
            z = max_z

        if z not in used:
            used.add(z)
            result.append(z)

        index += 1

    return result


def evaluate_sequence(
    n,
    p,
    q,
    sequence,
    budgets,
):
    """
    Evaluate a sequence cumulatively.

    Returns hit information for every requested budget.
    """

    results = {}

    found = False
    found_z = None
    found_gcd = None
    found_state = None

    budget_set = set(budgets)

    for index, z in enumerate(
        sequence,
        start=1,
    ):
        g = value_gcd(
            n,
            len(sequence) * 0 + CURRENT_K,
            z,
        )

        state = gcd_state(
            g,
            n,
            p,
            q,
        )

        if (
            not found
            and (
                state == "p"
                or state == "q"
            )
        ):
            found = True
            found_z = z
            found_gcd = g
            found_state = state

        if index in budget_set:
            results[index] = {
                "found": found,
                "z": found_z,
                "gcd": found_gcd,
                "state": found_state,
            }

    return results


def evaluate_strategy(
    n,
    p,
    q,
    k,
    sequence,
):
    """
    Evaluate one candidate sequence for one (N,k) case.
    """

    results = {}

    found = False
    found_z = None
    found_gcd = None
    found_state = None

    budget_set = set(PROBE_BUDGETS)

    for index, z in enumerate(
        sequence,
        start=1,
    ):
        g = value_gcd(
            n,
            k,
            z,
        )

        state = gcd_state(
            g,
            n,
            p,
            q,
        )

        if (
            not found
            and (
                state == "p"
                or state == "q"
            )
        ):
            found = True
            found_z = z
            found_gcd = g
            found_state = state

        if index in budget_set:
            results[index] = {
                "found": found,
                "z": found_z,
                "gcd": found_gcd,
                "state": found_state,
            }

    return results


def run_experiment():
    global CURRENT_K

    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q MAX: {Q_MAX}")
    print(f"MAX Z: {MAX_Z}")

    print()
    print("PROBE BUDGETS:")
    print(
        " ".join(
            str(x)
            for x in PROBE_BUDGETS
        )
    )

    print()

    random.seed(112)

    strategies = (
        "random",
        "even",
        "van_der_corput",
    )

    hit_counts = {
        strategy: {
            budget: 0
            for budget in PROBE_BUDGETS
        }
        for strategy in strategies
    }

    total_k_cases = 0

    example_misses = {
        strategy: []
        for strategy in strategies
    }

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            max_z = min(
                k // 2,
                MAX_Z,
            )

            max_budget = min(
                max(PROBE_BUDGETS),
                max_z,
            )

            # -------------------------------------------------
            # RANDOM
            # -------------------------------------------------

            random_seq = random_sequence(
                max_z,
                max_budget,
                seed=(
                    case_index * 1000003
                    + n * 1009
                    + k * 9176
                ),
            )

            random_results = (
                evaluate_strategy(
                    n,
                    p,
                    q,
                    k,
                    random_seq,
                )
            )

            # -------------------------------------------------
            # EVENLY SPACED
            # -------------------------------------------------

            even_seq = evenly_spaced_sequence(
                max_z,
                max_budget,
            )

            even_results = (
                evaluate_strategy(
                    n,
                    p,
                    q,
                    k,
                    even_seq,
                )
            )

            # -------------------------------------------------
            # VAN DER CORPUT
            # -------------------------------------------------

            vdc_seq = van_der_corput_sequence(
                max_z,
                max_budget,
            )

            vdc_results = (
                evaluate_strategy(
                    n,
                    p,
                    q,
                    k,
                    vdc_seq,
                )
            )

            result_sets = {
                "random": random_results,
                "even": even_results,
                "van_der_corput": vdc_results,
            }

            for strategy in strategies:
                results = result_sets[
                    strategy
                ]

                for budget in PROBE_BUDGETS:
                    effective_budget = min(
                        budget,
                        max_z,
                    )

                    if (
                        effective_budget
                        not in results
                    ):
                        continue

                    if results[
                        effective_budget
                    ]["found"]:
                        hit_counts[
                            strategy
                        ][budget] += 1

                if not results[
                    min(
                        max(PROBE_BUDGETS),
                        max_z,
                    )
                ]["found"]:
                    if (
                        len(
                            example_misses[
                                strategy
                            ]
                        )
                        < 10
                    ):
                        example_misses[
                            strategy
                        ].append(
                            (
                                case_index,
                                n,
                                p,
                                q,
                                k,
                            )
                        )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = "
        f"{CASES}"
    )

    print(
        f"(N,k) cases                  = "
        f"{total_k_cases}"
    )

    print()

    for strategy in strategies:
        print()
        print(
            f"[{strategy.upper()}]"
        )

        for budget in PROBE_BUDGETS:
            hits = hit_counts[
                strategy
            ][budget]

            rate = (
                100.0
                * hits
                / total_k_cases
                if total_k_cases
                else 0.0
            )

            print(
                f"probes={budget:2d} "
                f"hits={hits:4d} "
                f"rate={rate:6.2f}%"
            )

    print()
    print("-" * 60)
    print("BEST STRATEGY BY BUDGET")
    print("-" * 60)

    for budget in PROBE_BUDGETS:
        best_strategy = None
        best_hits = -1

        for strategy in strategies:
            hits = hit_counts[
                strategy
            ][budget]

            if hits > best_hits:
                best_hits = hits
                best_strategy = strategy

        rate = (
            100.0
            * best_hits
            / total_k_cases
            if total_k_cases
            else 0.0
        )

        print(
            f"{budget:2d} probes -> "
            f"{best_strategy} "
            f"{best_hits}/{total_k_cases} "
            f"({rate:.2f}%)"
        )

    print()
    print("-" * 60)
    print("EXAMPLE REMAINING MISSES")
    print("-" * 60)

    for strategy in strategies:
        print()
        print(
            f"[{strategy.upper()}]"
        )

        misses = example_misses[
            strategy
        ]

        if not misses:
            print("NONE")
            continue

        for (
            case_index,
            n,
            p,
            q,
            k,
        ) in misses:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k}"
            )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
