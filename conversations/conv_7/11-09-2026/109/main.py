import math
import random

from sympy import isprime


EXPERIMENT = 134

CASES = 100

P_MIN = 10_000_000
P_MAX = 20_000_000

Q_MIN = 20_000_000
Q_MAX = 50_000_000

# Widely separated k values.
K_VALUES = [
    1 << 21,  # 2,097,152
    1 << 22,  # 4,194,304
    1 << 23,  # 8,388,608
    1 << 24,  # 16,777,216
]

# d < p for every p in the experiment.
D = (1 << 20) - 1  # 1,048,575


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
    Generate p, q satisfying

        P_MIN <= p < P_MAX
        Q_MIN <= q < Q_MAX
        p < q
    """
    while True:
        p = random_prime(P_MIN, P_MAX)
        q = random_prime(Q_MIN, Q_MAX)

        if p < q:
            return p, q, p * q


def lucas_hit(p, k, d):
    """
    Exact Lucas criterion for d < p and d < k-1:

        p | C(k-2,d)
        iff
        d > (k-2) mod p.
    """
    if d >= p:
        raise ValueError("d must be < p")

    if d >= k - 1:
        return False

    remainder = (k - 2) % p

    return d > remainder


def pattern_for_p(p):
    """
    Return the four-bit Lucas hit/miss pattern.
    """
    return tuple(
        lucas_hit(p, k, D)
        for k in K_VALUES
    )


def pattern_string(pattern):
    """
    Convert Boolean pattern to HHMM-style text.
    """
    return "".join(
        "H" if value else "M"
        for value in pattern
    )


def candidate_p_intervals(k_values, d, p_min, p_max):
    """
    Analytically determine the p intervals producing each
    hit/miss pattern.

    For each k we require

        d > (k-2) mod p.

    We generate all quotient intervals of

        k-2 = a*p + r.

    Within one quotient region,

        r = k-2-a*p,

    so the hit condition becomes a simple linear inequality.

    The resulting intervals are then intersected across k.
    """

    # Candidate boundaries where floor((k-2)/p) changes.
    boundaries = {p_min, p_max}

    for k in k_values:
        m = k - 2

        # m / p is small in this regime.
        a = m // p_min

        while a >= 1:
            left = (m // (a + 1)) + 1
            right = m // a

            if p_min <= left < p_max:
                boundaries.add(left)

            if p_min < right < p_max:
                boundaries.add(right + 1)

            if left <= p_min:
                break

            a -= 1

    boundaries = sorted(boundaries)

    atomic_intervals = []

    for i in range(len(boundaries) - 1):
        left = boundaries[i]
        right = boundaries[i + 1] - 1

        if left > right:
            continue

        midpoint = (left + right) // 2

        pattern = pattern_for_p(midpoint)

        atomic_intervals.append(
            (left, right, pattern)
        )

    # Merge adjacent intervals with the same pattern.
    merged = []

    for left, right, pattern in atomic_intervals:
        if merged and merged[-1][2] == pattern:
            old_left, old_right, old_pattern = merged[-1]

            if old_right + 1 == left:
                merged[-1] = (
                    old_left,
                    right,
                    old_pattern,
                )
                continue

        merged.append(
            (left, right, pattern)
        )

    return merged


def count_primes_in_interval(left, right):
    """
    Count primes in [left, right].

    This is only used for the summary. The main structural
    analysis does not require prime enumeration.
    """
    if right < left:
        return 0

    count = 0

    # Lightweight sampling-free exact count for manageable ranges.
    for candidate in range(left | 1, right + 1, 2):
        if isprime(candidate):
            count += 1

    return count


def analyze_pattern_intervals(intervals):
    """
    Print structural pattern intervals.
    """
    print("PATTERN INTERVALS")

    for left, right, pattern in intervals:
        width = right - left + 1

        print(
            f"  {left:>10} .. {right:<10} "
            f"width={width:<10} "
            f"pattern={pattern_string(pattern)}"
        )

    print()


def summarize_observed_cases(cases):
    """
    Summarize the patterns observed among random cases.
    """
    counts = {}

    for case in cases:
        pattern = case["pattern"]

        counts[pattern] = counts.get(pattern, 0) + 1

    print("OBSERVED PATTERNS")

    total = len(cases)

    for pattern, count in sorted(
        counts.items(),
        key=lambda item: pattern_string(item[0])
    ):
        fraction = count / total

        print(
            f"  {pattern_string(pattern)}: "
            f"{count}/{total} "
            f"({fraction:.2%})"
        )

    print()


def random_case_analysis(intervals, p):
    """
    Find the structural interval containing p.
    """
    for left, right, pattern in intervals:
        if left <= p <= right:
            return left, right, pattern

    return None


def run_case(case_number, intervals):
    """
    Generate and analyze one random semiprime.
    """
    p, q, n = generate_case()

    pattern = pattern_for_p(p)

    located = random_case_analysis(
        intervals,
        p
    )

    print(f"CASE {case_number}")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {n}")
    print(f"  pattern = {pattern_string(pattern)}")

    if located is None:
        raise RuntimeError(
            f"p={p} did not appear in an interval"
        )

    left, right, interval_pattern = located

    print(
        f"  candidate interval = "
        f"[{left}, {right}]"
    )

    print(
        f"  interval width = "
        f"{right - left + 1}"
    )

    print(
        f"  contains true p = "
        f"{left <= p <= right}"
    )

    print()

    return {
        "p": p,
        "q": q,
        "n": n,
        "pattern": pattern,
        "interval_left": left,
        "interval_right": right,
    }


def main():
    """
    Main experiment.
    """
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX}")
    print(f"D: {D}")
    print(f"K VALUES: {K_VALUES}")
    print()

    random.seed(134)

    print("Computing exact p-pattern intervals...")

    intervals = candidate_p_intervals(
        K_VALUES,
        D,
        P_MIN,
        P_MAX,
    )

    analyze_pattern_intervals(intervals)

    results = []

    for case_number in range(1, CASES + 1):
        result = run_case(
            case_number,
            intervals,
        )

        results.append(result)

    summarize_observed_cases(results)

    print("INTERVAL INFORMATION")

    widths = []

    for result in results:
        widths.append(
            result["interval_right"]
            - result["interval_left"]
            + 1
        )

    print(
        f"  minimum interval width = "
        f"{min(widths)}"
    )

    print(
        f"  maximum interval width = "
        f"{max(widths)}"
    )

    print(
        f"  average interval width = "
        f"{sum(widths) / len(widths):.2f}"
    )

    print()

    # Theoretical information from the number of patterns.
    pattern_set = {
        result["pattern"]
        for result in results
    }

    print(
        f"Patterns observed in random sample = "
        f"{len(pattern_set)}"
    )

    print(
        f"Maximum possible information from "
        f"{len(pattern_set)} states = "
        f"{math.log2(len(pattern_set)):.4f} bits"
        if pattern_set
        else "No patterns observed"
    )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


if __name__ == "__main__":
    main()

