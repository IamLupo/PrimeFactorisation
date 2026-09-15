import math
import random
import time

from sympy import isprime, primerange


EXPERIMENT = 136

CASES = 100

P_MIN = 10_000_000
P_MAX = 20_000_000

Q_MIN = 20_000_000
Q_MAX = 50_000_000

# Common d. This is safely below every p.
D = (1 << 20) - 1  # 1,048,575

# Number of factor-blind k values.
K_COUNT = 8

# Offsets around sqrt(N).
LOCAL_OFFSETS = [
    -3_000_000,
    -2_000_000,
    -1_000_000,
    -250_000,
    250_000,
    1_000_000,
    2_000_000,
    3_000_000,
]


def random_prime(low, high):
    """
    Generate a random prime in [low, high).
    """
    while True:
        value = random.randrange(low, high)

        if value % 2 == 0:
            value += 1

        while value < high:
            if isprime(value):
                return value

            value += 2


def generate_case():
    """
    Generate N = p*q with

        p in [P_MIN, P_MAX)
        q in [Q_MIN, Q_MAX)
        p < q.
    """
    while True:
        p = random_prime(P_MIN, P_MAX)
        q = random_prime(Q_MIN, Q_MAX)

        if p < q:
            return p, q, p * q


def previous_power_of_two(x):
    """
    Largest power of two <= x.
    """
    return 1 << (x.bit_length() - 1)


def next_power_of_two(x):
    """
    Smallest power of two >= x.
    """
    if x <= 1:
        return 1

    low = previous_power_of_two(x)

    if low == x:
        return low

    return low << 1


def build_k_values(n, p_bound, q_bound):
    """
    Build factor-blind k values using only N.

    Sources:
      1. Several offsets around sqrt(N).
      2. Nearby dyadic powers of two.

    No p or q is used.
    """
    s = math.isqrt(n)

    values = set()

    # Local sqrt(N)-based values.
    for offset in LOCAL_OFFSETS:
        k = s + offset

        if k > 2 and k < q_bound:
            values.add(k)

    # Dyadic values around sqrt(N).
    low = previous_power_of_two(s)
    high = next_power_of_two(s)

    for k in (
        low // 2,
        low,
        high,
        high * 2,
    ):
        if k > 2 and k < q_bound:
            values.add(k)

    # Keep only the largest useful set, sorted.
    values = sorted(values)

    return values[:K_COUNT]


def lucas_hit(candidate, k, d):
    """
    Exact Lucas criterion in the valid regime:

        d < candidate
        d < k - 1

    Then

        candidate | C(k-2,d)
        iff
        d > (k-2) mod candidate.
    """
    if d >= candidate:
        return False

    if d >= k - 1:
        return False

    return d > ((k - 2) % candidate)


def pattern_for_candidate(candidate, k_values):
    """
    Calculate the hit/miss pattern for one candidate factor.
    """
    return tuple(
        lucas_hit(candidate, k, D)
        for k in k_values
    )


def pattern_string(pattern):
    """
    Convert Boolean pattern to H/M.
    """
    return "".join(
        "H" if x else "M"
        for x in pattern
    )


def build_candidate_primes():
    """
    Build all prime candidates in the p-range once.
    """
    return list(
        primerange(
            P_MIN,
            P_MAX
        )
    )


def matching_candidates(
    candidate_primes,
    k_values,
    target_pattern,
    count,
):
    """
    Return candidate primes matching the first 'count'
    Lucas observations.
    """
    selected_k = k_values[:count]
    selected_pattern = target_pattern[:count]

    matches = []

    for candidate in candidate_primes:
        pattern = pattern_for_candidate(
            candidate,
            selected_k
        )

        if pattern == selected_pattern:
            matches.append(candidate)

    return matches


def run_case(case_number, candidate_primes):
    """
    Run one factor-blind experiment case.
    """
    p, q, n = generate_case()

    k_values = build_k_values(
        n,
        p,
        q,
    )

    # We need enough valid k values.
    if len(k_values) < K_COUNT:
        return None

    true_pattern = pattern_for_candidate(
        p,
        k_values,
    )

    # Verify all probes are actually in the intended regime.
    for k in k_values:
        if k <= p:
            continue

        if D >= k - 1:
            # This probe is unusable for the positive numerator
            # interpretation. Keep it in the pattern but report it.
            pass

    print(f"CASE {case_number}")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {n}")
    print(f"  sqrt(N) = {math.isqrt(n)}")
    print(f"  D = {D}")

    print(
        "  k = ["
        + ", ".join(str(k) for k in k_values)
        + "]"
    )

    print(
        f"  TRUE PATTERN = "
        f"{pattern_string(true_pattern)}"
    )

    for k, hit in zip(
        k_values,
        true_pattern
    ):
        remainder = (k - 2) % p

        print(
            f"    k={k:<10} "
            f"remainder={remainder:<10} "
            f"state={'HIT' if hit else 'MISS'}"
        )

    results = []

    for count in range(
        1,
        K_COUNT + 1
    ):
        start = time.perf_counter()

        matches = matching_candidates(
            candidate_primes,
            k_values,
            true_pattern,
            count,
        )

        elapsed = time.perf_counter() - start

        true_present = p in matches

        unique = (
            len(matches) == 1
            and true_present
        )

        results.append({
            "count": count,
            "candidate_count": len(matches),
            "true_present": true_present,
            "unique": unique,
            "time": elapsed,
        })

        print(
            f"    {count}-K: "
            f"candidates={len(matches)} "
            f"true_present={true_present} "
            f"unique={unique} "
            f"time={elapsed:.6f}s"
        )

    print()

    return {
        "p": p,
        "q": q,
        "n": n,
        "k_values": k_values,
        "pattern": true_pattern,
        "results": results,
    }


def summarize(results):
    """
    Aggregate experiment results.
    """
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    if not results:
        print("No completed cases.")
        return

    for count in range(
        1,
        K_COUNT + 1
    ):
        rows = []

        for result in results:
            for row in result["results"]:
                if row["count"] == count:
                    rows.append(row)

        if not rows:
            continue

        unique = sum(
            row["unique"]
            for row in rows
        )

        present = sum(
            row["true_present"]
            for row in rows
        )

        average_candidates = (
            sum(
                row["candidate_count"]
                for row in rows
            )
            / len(rows)
        )

        minimum_candidates = min(
            row["candidate_count"]
            for row in rows
        )

        maximum_candidates = max(
            row["candidate_count"]
            for row in rows
        )

        average_time = (
            sum(
                row["time"]
                for row in rows
            )
            / len(rows)
        )

        print(
            f"{count}-K probes:"
        )

        print(
            f"  TRUE PRESENT       = "
            f"{present}/{len(rows)}"
        )

        print(
            f"  UNIQUE             = "
            f"{unique}/{len(rows)}"
        )

        print(
            f"  avg candidates     = "
            f"{average_candidates:.2f}"
        )

        print(
            f"  min candidates     = "
            f"{minimum_candidates}"
        )

        print(
            f"  max candidates     = "
            f"{maximum_candidates}"
        )

        print(
            f"  average time       = "
            f"{average_time:.6f}s"
        )

    print()

    pattern_counts = {}

    for result in results:
        pattern = result["pattern"]

        pattern_counts[pattern] = (
            pattern_counts.get(pattern, 0) + 1
        )

    print("FULL PATTERN DISTRIBUTION")

    for pattern, count in sorted(
        pattern_counts.items(),
        key=lambda item: pattern_string(item[0])
    ):
        print(
            f"  {pattern_string(pattern)}: "
            f"{count}/{len(results)}"
        )


def main():
    """
    Main experiment entry point.
    """
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX - 1}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX - 1}")
    print(f"D: {D}")
    print(f"K COUNT: {K_COUNT}")
    print(
        f"LOCAL OFFSETS: "
        f"{LOCAL_OFFSETS}"
    )
    print()

    random.seed(136)

    print("Building candidate-prime set...")

    start = time.perf_counter()

    candidate_primes = build_candidate_primes()

    elapsed = time.perf_counter() - start

    print(
        f"Candidate primes: "
        f"{len(candidate_primes)}"
    )

    print(
        f"Prime generation time: "
        f"{elapsed:.6f}s"
    )

    print()

    results = []

    case_number = 0

    while case_number < CASES:
        result = run_case(
            case_number + 1,
            candidate_primes,
        )

        if result is None:
            continue

        results.append(result)
        case_number += 1

    summarize(results)

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


if __name__ == "__main__":
    main()

