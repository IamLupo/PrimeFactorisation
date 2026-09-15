import math
import random
from sympy import isprime


EXPERIMENT = 108

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384

# Test z = c * 2^j.
# Increase this if the miss rate remains significant.
MULTIPLIERS = (
    1,
    3,
    5,
    7,
    9,
    11,
    13,
    15,
    17,
    19,
    21,
    23,
    25,
    27,
    29,
    31,
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
    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d)
    value %= n
    value = (-value) % n

    return math.gcd(value, n)


def generate_candidates():
    """
    Generate candidates of the form

        z = c * 2^j

    without duplicates and without exceeding MAX_Z.
    """

    candidates = set()

    for c in MULTIPLIERS:
        z = c

        while z <= MAX_Z:
            candidates.add(z)
            z *= 2

    return sorted(candidates)


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


def search_candidates(n, p, q, k, candidates):
    """
    Search the fixed candidate family.

    The candidates themselves do not depend on p or q.
    """

    evaluations = 0

    first_factor = None

    for z in candidates:
        d = 2 * z - 1

        if d >= k - 1:
            continue

        evaluations += 1

        g = value_gcd(n, k, z)
        state = gcd_state(g, n, p, q)

        if state == "p" or state == "q":
            first_factor = (
                z,
                g,
                state,
            )
            break

    return first_factor, evaluations


def identify_multiplier(z):
    """
    Return the smallest odd c such that

        z = c * 2^j

    """

    value = z

    while value % 2 == 0:
        value //= 2

    return value


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

    print("MULTIPLIERS:")
    print(" ".join(str(x) for x in MULTIPLIERS))
    print()

    random.seed(108)

    candidates = generate_candidates()

    print(
        f"candidate count               = "
        f"{len(candidates)}"
    )

    total_k_cases = 0

    dyadic_hits = 0
    extended_hits = 0
    complete_misses = 0

    additional_recoveries = 0

    total_evaluations = 0
    max_evaluations = 0

    multiplier_counts = {
        c: 0
        for c in MULTIPLIERS
    }

    miss_examples = []
    recovery_examples = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            # First, reproduce the original dyadic search.
            dyadic_candidates = []

            z = 1

            while z <= MAX_Z:
                dyadic_candidates.append(z)
                z *= 2

            dyadic_hit = None

            for z in dyadic_candidates:
                d = 2 * z - 1

                if d >= k - 1:
                    break

                g = value_gcd(n, k, z)

                if g == p or g == q:
                    dyadic_hit = (
                        z,
                        g,
                    )
                    break

            if dyadic_hit is not None:
                dyadic_hits += 1

            result, evaluations = search_candidates(
                n,
                p,
                q,
                k,
                candidates,
            )

            total_evaluations += evaluations

            if evaluations > max_evaluations:
                max_evaluations = evaluations

            if result is not None:
                extended_hits += 1

                z_hit, g_hit, state = result

                multiplier = identify_multiplier(
                    z_hit
                )

                if multiplier in multiplier_counts:
                    multiplier_counts[multiplier] += 1

                if (
                    dyadic_hit is None
                    and len(recovery_examples) < 20
                ):
                    additional_recoveries += 1

                    recovery_examples.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            z_hit,
                            multiplier,
                            g_hit,
                        )
                    )

            else:
                complete_misses += 1

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

    print()
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
        f"candidate count               = "
        f"{len(candidates)}"
    )

    print(
        f"original dyadic hits         = "
        f"{dyadic_hits}"
    )

    print(
        f"extended-family hits         = "
        f"{extended_hits}"
    )

    print(
        f"new recoveries               = "
        f"{additional_recoveries}"
    )

    print(
        f"complete misses              = "
        f"{complete_misses}"
    )

    print(
        f"total gcd evaluations        = "
        f"{total_evaluations}"
    )

    print(
        f"maximum gcd evaluations      = "
        f"{max_evaluations}"
    )

    if total_k_cases:
        print(
            f"original dyadic rate         = "
            f"{100.0 * dyadic_hits / total_k_cases:.2f}%"
        )

        print(
            f"extended-family rate         = "
            f"{100.0 * extended_hits / total_k_cases:.2f}%"
        )

        print(
            f"remaining miss rate          = "
            f"{100.0 * complete_misses / total_k_cases:.2f}%"
        )

    print()
    print("-" * 60)
    print("RECOVERIES BY MULTIPLIER")
    print("-" * 60)

    for c in MULTIPLIERS:
        print(
            f"c={c:2d} "
            f"hits={multiplier_counts[c]}"
        )

    print()
    print("-" * 60)
    print("NEW RECOVERY EXAMPLES")
    print("-" * 60)

    if recovery_examples:
        for (
            case_index,
            n,
            p,
            q,
            k,
            z,
            multiplier,
            g,
        ) in recovery_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"z={z} "
                f"c={multiplier} "
                f"gcd={g}"
            )
    else:
        print("NONE")

    print()
    print("-" * 60)
    print("REMAINING MISS EXAMPLES")
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
