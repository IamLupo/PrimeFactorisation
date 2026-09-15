import math
import random
from sympy import isprime


EXPERIMENT = 110

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384

FRACTIONS = (
    2,
    4,
    8,
    16,
    32,
    64,
    128,
)

NEIGHBOR_RADIUS = 4


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
    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d)
    value %= n
    value = (-value) % n

    return math.gcd(value, n)


def base_p_digits(value, p):
    digits = []

    while value > 0:
        digits.append(value % p)
        value //= p

    if not digits:
        digits.append(0)

    return digits


def smallest_lucas_zero_d(m, p):
    """
    Construct the smallest d >= 1 such that

        C(m,d) == 0 (mod p)

    from Lucas' theorem.
    """

    digits = base_p_digits(m, p)

    candidates = []

    power = 1

    for digit in digits:
        if digit + 1 < p:
            candidate = (digit + 1) * power
            candidates.append(candidate)

        power *= p

    if not candidates:
        return None

    return min(candidates)


def exact_factor_z(p, k):
    """
    Construct the smallest odd d giving p | v_k,
    then convert d = 2z-1 into z.
    """

    d = smallest_lucas_zero_d(
        k - 2,
        p,
    )

    if d is None:
        return None

    if d % 2 == 0:
        d += 1

    if d >= k - 1:
        return None

    return (d + 1) // 2


def generate_candidates(k):
    """
    Candidate positions based only on k.

    No p or q is used.
    """

    candidates = set()

    for denominator in FRACTIONS:
        center = k // denominator

        for offset in range(
            -NEIGHBOR_RADIUS,
            NEIGHBOR_RADIUS + 1,
        ):
            z = center + offset

            if 1 <= z <= MAX_Z:
                candidates.add(z)

    return sorted(candidates)


def candidate_search(n, p, q, k):
    """
    Search the factor-blind k-based candidates.
    """

    candidates = generate_candidates(k)

    evaluations = 0

    for z in candidates:
        d = 2 * z - 1

        if d >= k - 1:
            continue

        evaluations += 1

        g = value_gcd(n, k, z)
        state = gcd_state(
            g,
            n,
            p,
            q,
        )

        if state == "p" or state == "q":
            return {
                "found": True,
                "z": z,
                "gcd": g,
                "state": state,
                "evaluations": evaluations,
            }

    return {
        "found": False,
        "z": None,
        "gcd": None,
        "state": None,
        "evaluations": evaluations,
    }


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


def normalized_bucket(ratio):
    if ratio < 0.03125:
        return "<1/32"

    if ratio < 0.0625:
        return "1/32..1/16"

    if ratio < 0.125:
        return "1/16..1/8"

    if ratio < 0.25:
        return "1/8..1/4"

    if ratio < 0.5:
        return "1/4..1/2"

    return ">=1/2"


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

    print("FRACTIONS:")
    print(" ".join(f"1/{x}" for x in FRACTIONS))
    print()

    random.seed(110)

    total_k_cases = 0

    candidate_hits = 0
    candidate_misses = 0

    total_evaluations = 0
    maximum_evaluations = 0

    exact_z_sum_ratio = 0.0
    exact_z_min_ratio = 1.0
    exact_z_max_ratio = 0.0

    buckets = {
        "<1/32": 0,
        "1/32..1/16": 0,
        "1/16..1/8": 0,
        "1/8..1/4": 0,
        "1/4..1/2": 0,
        ">=1/2": 0,
    }

    improvement_examples = []
    miss_examples = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            exact_z = exact_factor_z(
                p,
                k,
            )

            if exact_z is not None:
                ratio = exact_z / k

                exact_z_sum_ratio += ratio
                exact_z_min_ratio = min(
                    exact_z_min_ratio,
                    ratio,
                )
                exact_z_max_ratio = max(
                    exact_z_max_ratio,
                    ratio,
                )

                buckets[
                    normalized_bucket(ratio)
                ] += 1

            result = candidate_search(
                n,
                p,
                q,
                k,
            )

            total_evaluations += (
                result["evaluations"]
            )

            maximum_evaluations = max(
                maximum_evaluations,
                result["evaluations"],
            )

            if result["found"]:
                candidate_hits += 1

                if (
                    len(improvement_examples) < 20
                    and exact_z is not None
                ):
                    candidate_z = result["z"]

                    if candidate_z != exact_z:
                        improvement_examples.append(
                            (
                                case_index,
                                n,
                                p,
                                q,
                                k,
                                exact_z,
                                candidate_z,
                                result["gcd"],
                            )
                        )
            else:
                candidate_misses += 1

                if len(miss_examples) < 20:
                    miss_examples.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            exact_z,
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

    print(
        f"k-based candidate hits       = "
        f"{candidate_hits}"
    )

    print(
        f"k-based candidate misses     = "
        f"{candidate_misses}"
    )

    print(
        f"total gcd evaluations        = "
        f"{total_evaluations}"
    )

    print(
        f"maximum gcd evaluations      = "
        f"{maximum_evaluations}"
    )

    if total_k_cases:
        print(
            f"candidate hit rate           = "
            f"{100.0 * candidate_hits / total_k_cases:.2f}%"
        )

        average_ratio = (
            exact_z_sum_ratio
            / total_k_cases
        )

        print(
            f"average exact z/k            = "
            f"{average_ratio:.6f}"
        )

        print(
            f"minimum exact z/k            = "
            f"{exact_z_min_ratio:.6f}"
        )

        print(
            f"maximum exact z/k            = "
            f"{exact_z_max_ratio:.6f}"
        )

    print()
    print("-" * 60)
    print("EXACT z/k DISTRIBUTION")
    print("-" * 60)

    for bucket, count in buckets.items():
        print(
            f"{bucket:12s} "
            f"count={count}"
        )

    print()
    print("-" * 60)
    print("NON-EXACT CANDIDATE HITS")
    print("-" * 60)

    if improvement_examples:
        for (
            case_index,
            n,
            p,
            q,
            k,
            exact_z,
            candidate_z,
            g,
        ) in improvement_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"exact_z={exact_z} "
                f"candidate_z={candidate_z} "
                f"gcd={g}"
            )
    else:
        print("NONE")

    print()
    print("-" * 60)
    print("MISS EXAMPLES")
    print("-" * 60)

    for (
        case_index,
        n,
        p,
        q,
        k,
        exact_z,
    ) in miss_examples:
        print(
            f"CASE={case_index} "
            f"N={n} "
            f"p={p} "
            f"q={q} "
            f"k={k} "
            f"exact_z={exact_z}"
        )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
