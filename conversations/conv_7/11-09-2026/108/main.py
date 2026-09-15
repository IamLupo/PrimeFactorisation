import math
import random
import time

from sympy import primerange


EXPERIMENT = 133

CASES = 100

# Keep the same scale as the previous dataset.
P_MIN = 10_000_000
P_MAX = 20_000_000

Q_MIN = 20_000_000
Q_MAX = 50_000_000

K_COUNT = 4

# Candidate primes are precomputed once.
# We only count candidates in the same p-range.
CANDIDATE_PRIME_LIMIT = P_MAX

# We want d < p for every generated case.
# 2^22 - 1 = 4,194,303, safely below 10^7.
D = (1 << 22) - 1

# k offsets around a power-of-two anchor.
# We deliberately change k rather than z.
K_OFFSETS = [-4096, -2048, 2048, 4096]


def generate_prime(low, high):
    """
    Generate a random prime in [low, high).
    """
    while True:
        value = random.randrange(low, high)

        if value % 2 == 0:
            value += 1

        # trial through a small amount of odd progression
        for candidate in range(value, high, 2):
            if is_prime_fast(candidate):
                return candidate

        # Retry if we reached the upper boundary.
        continue


def is_prime_fast(n):
    """
    Small deterministic primality test for the generated range.

    SymPy primerange is used for candidate-space construction,
    while this keeps random generation simple and independent.
    """
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    small_primes = (
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True

        if n % p == 0:
            return False

    limit = math.isqrt(n)

    divisor = 41

    while divisor <= limit:
        if n % divisor == 0:
            return False

        divisor += 2

    return True


def generate_case():
    """
    Generate N = p*q with p in the requested p-range
    and q in the requested q-range.
    """
    while True:
        p = generate_prime(P_MIN, P_MAX)
        q = generate_prime(Q_MIN, Q_MAX)

        if p < q:
            n = p * q
            return p, q, n


def largest_power_two_near(value):
    """
    Return the power of two nearest to value.
    """
    lower = 1 << (value.bit_length() - 1)
    upper = lower << 1

    if value - lower <= upper - value:
        return lower

    return upper


def select_k_values(p, q, n):
    """
    Select four genuinely different k values near sqrt(N).

    Every k must satisfy

        p < k < q
        d < k - 1
    """
    s = math.isqrt(n)

    anchor = largest_power_two_near(s)

    candidates = []

    for offset in K_OFFSETS:
        k = anchor + offset

        if p < k < q:
            candidates.append(k)

    # Remove duplicates and sort.
    candidates = sorted(set(candidates))

    if len(candidates) != K_COUNT:
        return []

    for k in candidates:
        if D >= p:
            return []

        if D >= k - 1:
            return []

    return candidates


def lucas_hit(candidate, k, d):
    """
    Exact Lucas criterion for d < candidate:

        candidate | C(k-2, d)
        iff
        d > (k-2) mod candidate.
    """
    if d >= candidate:
        return False

    if d >= k - 1:
        return False

    return d > ((k - 2) % candidate)


def direct_lucas_pattern(p, k_values, d):
    """
    Compute the observable hit/miss pattern for the true factor.
    """
    return tuple(
        lucas_hit(p, k, d)
        for k in k_values
    )


def pattern_string(pattern):
    """
    Format a Boolean pattern as H/M.
    """
    return "".join(
        "H" if hit else "M"
        for hit in pattern
    )


def build_candidate_primes():
    """
    Build all candidate primes in [P_MIN, P_MAX).

    This is done once for the whole experiment.
    """
    return list(
        primerange(
            P_MIN,
            P_MAX
        )
    )


def matching_candidates(candidate_primes, k_values, d, pattern):
    """
    Return candidate primes producing exactly the same
    Lucas pattern.
    """
    matches = []

    for candidate in candidate_primes:
        candidate_pattern = direct_lucas_pattern(
            candidate,
            k_values,
            d
        )

        if candidate_pattern == pattern:
            matches.append(candidate)

    return matches


def summarize_candidate_set(matches, p):
    """
    Return useful information about the candidate set.
    """
    count = len(matches)

    if count == 0:
        classification = "ZERO"
    elif count == 1 and matches[0] == p:
        classification = "UNIQUE"
    elif p in matches:
        classification = "AMBIGUOUS"
    else:
        classification = "MISSING_TRUE_P"

    return count, classification


def run_case(case_number, candidate_primes):
    """
    Run one case and progressively evaluate the information
    provided by 1, 2, 3 and 4 k-values.
    """
    p, q, n = generate_case()

    k_values = select_k_values(
        p,
        q,
        n
    )

    if not k_values:
        return None

    pattern = direct_lucas_pattern(
        p,
        k_values,
        D
    )

    print(f"CASE {case_number}")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {n}")
    print(f"  sqrt(N) = {math.isqrt(n)}")
    print(f"  D = {D}")
    print(f"  k = {k_values}")
    print(
        f"  pattern = "
        f"{pattern_string(pattern)}"
    )

    results = []

    for count in range(1, K_COUNT + 1):
        selected_k = k_values[:count]
        selected_pattern = pattern[:count]

        start = time.perf_counter()

        matches = matching_candidates(
            candidate_primes,
            selected_k,
            D,
            selected_pattern
        )

        elapsed = time.perf_counter() - start

        candidate_count, classification = (
            summarize_candidate_set(
                matches,
                p
            )
        )

        results.append({
            "count": count,
            "candidate_count": candidate_count,
            "classification": classification,
            "time": elapsed,
        })

        preview = matches[:10]

        print(
            f"    {count}-K: "
            f"candidates={candidate_count} "
            f"classification={classification} "
            f"time={elapsed:.6f}s"
        )

        if preview:
            print(
                f"          preview={preview}"
            )

    print()

    return {
        "p": p,
        "q": q,
        "n": n,
        "k_values": k_values,
        "pattern": pattern,
        "results": results,
    }


def summarize(all_results):
    """
    Print aggregate results.
    """
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    for count in range(1, K_COUNT + 1):
        rows = []

        for result in all_results:
            for row in result["results"]:
                if row["count"] == count:
                    rows.append(row)

        if not rows:
            continue

        unique = sum(
            row["classification"] == "UNIQUE"
            for row in rows
        )

        ambiguous = sum(
            row["classification"] == "AMBIGUOUS"
            for row in rows
        )

        zero = sum(
            row["classification"] == "ZERO"
            for row in rows
        )

        missing = sum(
            row["classification"] == "MISSING_TRUE_P"
            for row in rows
        )

        average_candidates = (
            sum(
                row["candidate_count"]
                for row in rows
            )
            / len(rows)
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
            f"  UNIQUE             = "
            f"{unique}/{len(rows)}"
        )

        print(
            f"  AMBIGUOUS          = "
            f"{ambiguous}/{len(rows)}"
        )

        print(
            f"  ZERO               = "
            f"{zero}/{len(rows)}"
        )

        print(
            f"  MISSING TRUE p     = "
            f"{missing}/{len(rows)}"
        )

        print(
            f"  average candidates = "
            f"{average_candidates:.2f}"
        )

        print(
            f"  average scan time  = "
            f"{average_time:.6f}s"
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
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX}")
    print(f"D: {D}")
    print(f"K COUNT: {K_COUNT}")
    print(f"K OFFSETS: {K_OFFSETS}")
    print()

    random.seed(133)

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
            candidate_primes
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
