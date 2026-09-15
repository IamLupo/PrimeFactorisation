import math
import random
import time

from sympy import isprime


EXPERIMENT = 135

CASES = 100

P_MIN = 10_000_000
P_MAX = 20_000_000

Q_MIN = 20_000_000
Q_MAX = 50_000_000

# One d for all probes.
# d < p for every generated case.
D = (1 << 20) - 1  # 1,048,575

# Multipliers used to construct k around multiples of the hidden p.
# k = a*p + offset
MULTIPLIERS = [1, 2, 3, 4, 5, 7, 8, 11]

# Small offsets prevent all k values from being exact multiples of p.
OFFSETS = [
    1,
    17,
    257,
    4097,
    65537,
    131071,
    262143,
    524287,
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

        P_MIN <= p < P_MAX
        Q_MIN <= q < Q_MAX
        p < q
    """
    while True:
        p = random_prime(P_MIN, P_MAX)
        q = random_prime(Q_MIN, Q_MAX)

        if p < q:
            return p, q, p * q


def construct_k_values(p):
    """
    Construct k values from different multiples of the hidden p.

    We only keep k values satisfying

        p < k < q

    so that the standard ArithSeg interpretation remains meaningful.
    """
    values = []

    for multiplier, offset in zip(MULTIPLIERS, OFFSETS):
        k = multiplier * p + offset

        if p < k:
            values.append(k)

    return values


def lucas_hit(p, k, d):
    """
    Exact Lucas condition in the valid regime:

        d < p
        d < k-1

    Then

        p | C(k-2,d)
        iff
        d > (k-2) mod p.
    """
    if d >= p:
        raise ValueError(
            f"d={d} must be < p={p}"
        )

    if d >= k - 1:
        return False

    return d > ((k - 2) % p)


def pattern_for_candidate(candidate, k_values, d):
    """
    Compute the hit/miss pattern predicted for a candidate factor.
    """
    pattern = []

    for k in k_values:
        # A candidate must also satisfy the Lucas regime.
        if d >= candidate:
            pattern.append(False)
            continue

        if d >= k - 1:
            pattern.append(False)
            continue

        pattern.append(
            lucas_hit(
                candidate,
                k,
                d
            )
        )

    return tuple(pattern)


def pattern_string(pattern):
    """
    Convert Boolean pattern to H/M.
    """
    return "".join(
        "H" if x else "M"
        for x in pattern
    )


def candidate_prime_in_range(candidate):
    """
    Candidate p values must be prime and within the
    experimental p range.
    """
    if candidate < P_MIN or candidate >= P_MAX:
        return False

    return isprime(candidate)


def generate_candidate_primes():
    """
    Generate all candidate primes in the p range.

    This is done once. The range contains roughly 600k primes,
    which is manageable but not repeated per case.
    """
    candidates = []

    for candidate in range(
        P_MIN | 1,
        P_MAX,
        2
    ):
        if isprime(candidate):
            candidates.append(candidate)

    return candidates


def find_matching_candidates(
    candidate_primes,
    k_values,
    d,
    target_pattern,
):
    """
    Find all candidate primes that generate exactly
    the same hit/miss pattern.
    """
    matches = []

    for candidate in candidate_primes:
        pattern = pattern_for_candidate(
            candidate,
            k_values,
            d
        )

        if pattern == target_pattern:
            matches.append(candidate)

    return matches


def case_information(
    candidate_primes,
    k_values,
    d,
    true_p,
):
    """
    Measure how candidate count changes as k values are added.
    """
    true_pattern = pattern_for_candidate(
        true_p,
        k_values,
        d
    )

    results = []

    for count in range(
        1,
        len(k_values) + 1
    ):
        selected_k = k_values[:count]

        selected_pattern = true_pattern[:count]

        start = time.perf_counter()

        matches = find_matching_candidates(
            candidate_primes,
            selected_k,
            d,
            selected_pattern
        )

        elapsed = time.perf_counter() - start

        true_present = true_p in matches

        results.append({
            "count": count,
            "candidate_count": len(matches),
            "true_present": true_present,
            "unique": (
                len(matches) == 1
                and true_present
            ),
            "time": elapsed,
            "pattern": selected_pattern,
        })

    return true_pattern, results


def run_case(
    case_number,
    candidate_primes,
):
    """
    Run one generated semiprime.
    """
    p, q, n = generate_case()

    k_values = construct_k_values(p)

    # Keep only k values satisfying k < q.
    k_values = [
        k
        for k in k_values
        if k < q
    ]

    if len(k_values) < 4:
        return None

    # The common d must remain below p.
    if D >= p:
        raise RuntimeError(
            f"Invalid d >= p: d={D}, p={p}"
        )

    true_pattern, results = case_information(
        candidate_primes,
        k_values,
        D,
        p,
    )

    print(f"CASE {case_number}")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {n}")
    print(f"  D = {D}")
    print(f"  k = {k_values}")

    print(
        f"  pattern = "
        f"{pattern_string(true_pattern)}"
    )

    for k, bit in zip(
        k_values,
        true_pattern
    ):
        remainder = (k - 2) % p

        print(
            f"    k={k:<12} "
            f"remainder={remainder:<10} "
            f"state={'HIT' if bit else 'MISS'}"
        )

    for row in results:
        print(
            f"    {row['count']}-K: "
            f"candidates={row['candidate_count']} "
            f"true_present={row['true_present']} "
            f"unique={row['unique']} "
            f"time={row['time']:.6f}s"
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


def summarize(all_results):
    """
    Aggregate experiment results.
    """
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    if not all_results:
        print("No completed cases.")
        return

    max_k_count = max(
        len(result["k_values"])
        for result in all_results
    )

    for count in range(
        1,
        max_k_count + 1
    ):
        rows = []

        for result in all_results:
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

        avg_candidates = (
            sum(
                row["candidate_count"]
                for row in rows
            )
            / len(rows)
        )

        avg_time = (
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
            f"  average candidates = "
            f"{avg_candidates:.2f}"
        )

        print(
            f"  average time       = "
            f"{avg_time:.6f}s"
        )

    pattern_counts = {}

    for result in all_results:
        pattern = result["pattern"]

        pattern_counts[pattern] = (
            pattern_counts.get(pattern, 0) + 1
        )

    print()
    print("FULL PATTERN DISTRIBUTION")

    for pattern, count in sorted(
        pattern_counts.items(),
        key=lambda item: pattern_string(item[0])
    ):
        print(
            f"  {pattern_string(pattern)}: "
            f"{count}/{len(all_results)}"
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
    print(f"MULTIPLIERS: {MULTIPLIERS}")
    print(f"OFFSETS: {OFFSETS}")
    print()

    random.seed(135)

    print("Building candidate-prime set...")

    start = time.perf_counter()

    candidate_primes = generate_candidate_primes()

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

