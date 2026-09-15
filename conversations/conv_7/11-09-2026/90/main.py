import math
import random
from sympy import isprime


EXPERIMENT = 111

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

        v_k = -C(k-2,d),  d < k-1

    and v_k = 0 otherwise.
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
        p = random.randrange(P_MIN, P_MAX + 1)

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


def random_probe_sequence(k, count, rng):
    """
    Random unique z values from the complete nontrivial range:

        1 <= z <= floor(k/2)

    No linear scan is performed.
    """
    max_z = min(
        k // 2,
        MAX_Z,
    )

    if max_z <= 0:
        return []

    count = min(count, max_z)

    return rng.sample(
        range(1, max_z + 1),
        count,
    )


def run_random_search(n, p, q, k, rng):
    """
    Build one random probe sequence and evaluate its
    cumulative success at each budget.
    """

    max_budget = max(PROBE_BUDGETS)

    probes = random_probe_sequence(
        k,
        max_budget,
        rng,
    )

    results = {}

    factor_found = False
    factor_z = None
    factor_gcd = None
    factor_state = None

    budget_set = set(PROBE_BUDGETS)

    for index, z in enumerate(probes, start=1):
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
            not factor_found
            and (state == "p" or state == "q")
        ):
            factor_found = True
            factor_z = z
            factor_gcd = g
            factor_state = state

        if index in budget_set:
            results[index] = {
                "found": factor_found,
                "z": factor_z,
                "gcd": factor_gcd,
                "state": factor_state,
            }

    return results


def make_case_seed(case_index, n, k):
    """
    Convert the case identifiers into a deterministic integer seed.

    random.Random only accepts supported scalar seed types.
    """
    return (
        (case_index + 1) * 1000003
        + n * 1009
        + k * 9176
    )


def run_experiment():
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
    print(" ".join(str(x) for x in PROBE_BUDGETS))
    print()

    random.seed(111)

    total_k_cases = 0

    hit_counts = {
        budget: 0
        for budget in PROBE_BUDGETS
    }

    first_hit_counts = {
        budget: 0
        for budget in PROBE_BUDGETS
    }

    exact_first_hit_z = []

    miss_examples = []
    hit_examples = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            seed = make_case_seed(
                case_index,
                n,
                k,
            )

            rng = random.Random(seed)

            results = run_random_search(
                n,
                p,
                q,
                k,
                rng,
            )

            already_found = False

            for budget in PROBE_BUDGETS:
                result = results[budget]

                if result["found"]:
                    hit_counts[budget] += 1

                    if not already_found:
                        first_hit_counts[budget] += 1
                        already_found = True

                        exact_first_hit_z.append(
                            result["z"]
                        )

                        if len(hit_examples) < 15:
                            hit_examples.append(
                                (
                                    case_index,
                                    n,
                                    p,
                                    q,
                                    k,
                                    budget,
                                    result["z"],
                                    result["gcd"],
                                )
                            )

            if not results[
                max(PROBE_BUDGETS)
            ]["found"]:
                if len(miss_examples) < 20:
                    miss_examples.append(
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

    for budget in PROBE_BUDGETS:
        hits = hit_counts[budget]

        rate = (
            100.0 * hits / total_k_cases
            if total_k_cases
            else 0.0
        )

        print(
            f"random probes={budget:2d} "
            f"hits={hits:4d} "
            f"rate={rate:6.2f}%"
        )

    print()
    print("-" * 60)
    print("FIRST-HIT DISTRIBUTION")
    print("-" * 60)

    previous = 0

    for budget in PROBE_BUDGETS:
        cumulative = hit_counts[budget]
        new_hits = cumulative - previous

        print(
            f"<= {budget:2d} probes: "
            f"{cumulative:4d} "
            f"(+{new_hits:4d})"
        )

        previous = cumulative

    print()
    print("-" * 60)
    print("AVERAGE PROBE COUNT TO RECOVER FACTOR")
    print("-" * 60)

    weighted_hits = 0
    weighted_count = 0

    previous_hits = 0

    for budget in PROBE_BUDGETS:
        new_hits = (
            hit_counts[budget]
            - previous_hits
        )

        weighted_hits += new_hits * budget
        weighted_count += new_hits

        previous_hits = hit_counts[budget]

    final_misses = (
        total_k_cases
        - hit_counts[max(PROBE_BUDGETS)]
    )

    if weighted_count:
        print(
            f"successful cases              = "
            f"{weighted_count}"
        )

        print(
            f"approx average probes        = "
            f"{weighted_hits / weighted_count:.2f}"
        )

    print(
        f"misses after 64 probes        = "
        f"{final_misses}"
    )

    if exact_first_hit_z:
        print()
        print("-" * 60)
        print("RANDOM HIT z RANGE")
        print("-" * 60)

        print(
            f"minimum z                    = "
            f"{min(exact_first_hit_z)}"
        )

        print(
            f"maximum z                    = "
            f"{max(exact_first_hit_z)}"
        )

    print()
    print("-" * 60)
    print("EXAMPLE HITS")
    print("-" * 60)

    for (
        case_index,
        n,
        p,
        q,
        k,
        budget,
        z,
        g,
    ) in hit_examples:
        print(
            f"CASE={case_index} "
            f"N={n} "
            f"p={p} "
            f"q={q} "
            f"k={k} "
            f"budget={budget} "
            f"z={z} "
            f"gcd={g}"
        )

    print()
    print("-" * 60)
    print("EXAMPLE 64-PROBE MISSES")
    print("-" * 60)

    if miss_examples:
        for (
            case_index,
            n,
            p,
            q,
            k,
        ) in miss_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k}"
            )
    else:
        print("NONE")

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()