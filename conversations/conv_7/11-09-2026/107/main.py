import math
import random
import time

from sympy import nextprime


EXPERIMENT = 132

CASES = 6

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

D = 65_535

MAX_PROBES = 4

CASE_TIME_LIMIT = 5.0


def generate_prime(low, high):
    """
    Generate a random prime in [low, high).
    """
    value = random.randrange(low, high)

    if value % 2 == 0:
        value += 1

    return int(nextprime(value))


def generate_semiprime():
    """
    Generate an odd semiprime N = p*q with
    p < q and the requested size ranges.
    """
    while True:
        p = generate_prime(P_MIN, P_MAX)
        q = generate_prime(Q_MIN, Q_MAX)

        if p < q:
            return p, q, p * q


def largest_power_of_two_leq(n):
    """
    Return the largest power of two <= n.
    """
    return 1 << (n.bit_length() - 1)


def select_k_values(s):
    """
    Select up to MAX_PROBES dyadic k-values at or below sqrt(N).

    All selected k satisfy k <= s < q, so a successful numerator
    interval cannot contain q. This isolates the p-side behavior.
    """
    base = largest_power_of_two_leq(s)

    values = []

    k = base

    while k >= 2 and len(values) < MAX_PROBES:
        values.append(k)
        k //= 2

    values.reverse()

    return values


def lucas_hit(p, k, d):
    """
    Exact Lucas criterion in the regime d < p:

        p | C(k-2, d)  <=>  d > (k-2) mod p
    """
    if d >= p:
        raise ValueError("lucas_hit requires d < p")

    remainder = (k - 2) % p

    return remainder < d


def linear_probe(N, k, d):
    """
    Directly evaluate the numerator product

        prod_{j=0}^{d-1} (k-2-j)

    modulo N and stop immediately once a proper factor appears.

    Because d < p and k <= sqrt(N) < q, a successful probe
    isolates p rather than q.
    """
    m = k - 2

    product = 1

    for j in range(d):
        value = m - j

        product = (product * value) % N

        g = math.gcd(product, N)

        if 1 < g < N:
            return {
                "hit": True,
                "factor": g,
                "iterations": j + 1,
            }

    return {
        "hit": False,
        "factor": None,
        "iterations": d,
    }


def observe_pattern(N, p, k_values, d):
    """
    Obtain the actual observable hit/miss pattern using the
    factorization evaluator.

    The returned pattern contains only hit/miss information.
    """
    pattern = []
    timings = []

    for k in k_values:
        start = time.perf_counter()

        result = linear_probe(N, k, d)

        elapsed = time.perf_counter() - start

        hit = result["hit"]

        expected = lucas_hit(p, k, d)

        if hit != expected:
            raise RuntimeError(
                f"Lucas mismatch: p={p}, k={k}, d={d}, "
                f"observed={hit}, expected={expected}"
            )

        pattern.append(hit)
        timings.append({
            "k": k,
            "hit": hit,
            "factor": result["factor"],
            "iterations": result["iterations"],
            "time": elapsed,
        })

    return pattern, timings


def candidate_matches_pattern(candidate, k_values, d, pattern):
    """
    Check whether candidate p would produce exactly the observed
    Lucas hit/miss pattern.
    """
    for k, observed in zip(k_values, pattern):
        predicted = lucas_hit(candidate, k, d)

        if predicted != observed:
            return False

    return True


def count_candidates(k_values, d, pattern, upper_bound):
    """
    Enumerate every integer candidate c in [P_MIN, upper_bound]
    satisfying the observed Lucas pattern.

    Since d < P_MIN, the Lucas condition is valid for every
    candidate considered.
    """
    candidates = []

    for candidate in range(P_MIN, upper_bound + 1):
        if candidate_matches_pattern(candidate, k_values, d, pattern):
            candidates.append(candidate)

    return candidates


def classify_candidate_set(candidates, p):
    """
    Classify whether the observed pattern uniquely identifies p
    among the scanned candidate integers.
    """
    if len(candidates) == 1 and candidates[0] == p:
        return "UNIQUE"

    if p not in candidates:
        return "MISS"

    return "AMBIGUOUS"


def print_pattern(pattern):
    """
    Format a hit/miss pattern.
    """
    return "".join("H" if value else "." for value in pattern)


def run_case(case_number):
    """
    Run one complete experiment case.
    """
    p, q, N = generate_semiprime()

    s = math.isqrt(N)

    k_values = select_k_values(s)

    if len(k_values) < MAX_PROBES:
        print("  SKIP: insufficient k-values")
        return None

    if D >= p:
        raise RuntimeError(
            f"Invalid case: D={D} is not < p={p}"
        )

    pattern, timings = observe_pattern(
        N,
        p,
        k_values,
        D,
    )

    print(
        f"CASE {case_number}: "
        f"p={p} q={q} N={N} sqrt={s}"
    )

    print(
        f"  D={D} "
        f"k_values={k_values} "
        f"pattern={print_pattern(pattern)}"
    )

    total_probe_time = sum(item["time"] for item in timings)

    for item in timings:
        print(
            f"    k={item['k']:>8} "
            f"state={'HIT' if item['hit'] else 'MISS':>4} "
            f"iterations={item['iterations']:>6} "
            f"time={item['time']:.6f}s"
        )

    print(
        f"  probe_time={total_probe_time:.6f}s"
    )

    cumulative_results = []

    for count in range(1, MAX_PROBES + 1):
        selected_k = k_values[:count]
        selected_pattern = pattern[:count]

        start = time.perf_counter()

        candidates = count_candidates(
            selected_k,
            D,
            selected_pattern,
            s,
        )

        elapsed = time.perf_counter() - start

        classification = classify_candidate_set(
            candidates,
            p,
        )

        cumulative_results.append({
            "count": count,
            "candidate_count": len(candidates),
            "classification": classification,
            "time": elapsed,
        })

        preview = candidates[:10]

        if len(candidates) > 10:
            preview_text = f"{preview} ..."
        else:
            preview_text = str(preview)

        print(
            f"  {count}-K: "
            f"candidates={len(candidates)} "
            f"classification={classification} "
            f"time={elapsed:.6f}s"
        )

        print(
            f"        preview={preview_text}"
        )

        if elapsed > CASE_TIME_LIMIT:
            print(
                f"        WARNING: candidate scan exceeded "
                f"{CASE_TIME_LIMIT:.1f}s"
            )

            break

    print()

    return {
        "p": p,
        "q": q,
        "N": N,
        "sqrt": s,
        "pattern": pattern,
        "probe_time": total_probe_time,
        "results": cumulative_results,
    }


def summarize(results):
    """
    Print aggregate experiment statistics.
    """
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    for count in range(1, MAX_PROBES + 1):
        rows = []

        for result in results:
            for item in result["results"]:
                if item["count"] == count:
                    rows.append(item)

        if not rows:
            continue

        unique = sum(
            1
            for item in rows
            if item["classification"] == "UNIQUE"
        )

        ambiguous = sum(
            1
            for item in rows
            if item["classification"] == "AMBIGUOUS"
        )

        misses = sum(
            1
            for item in rows
            if item["classification"] == "MISS"
        )

        average_candidates = (
            sum(item["candidate_count"] for item in rows)
            / len(rows)
        )

        average_time = (
            sum(item["time"] for item in rows)
            / len(rows)
        )

        print(
            f"{count}-K probes:"
        )

        print(
            f"  UNIQUE            = {unique}/{len(rows)}"
        )

        print(
            f"  AMBIGUOUS         = {ambiguous}/{len(rows)}"
        )

        print(
            f"  MISS              = {misses}/{len(rows)}"
        )

        print(
            f"  average candidates = {average_candidates:.2f}"
        )

        print(
            f"  average scan time  = {average_time:.6f}s"
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
    print(f"MAX PROBES: {MAX_PROBES}")
    print(f"CASE TIME LIMIT: {CASE_TIME_LIMIT:.1f}s")
    print()

    random.seed(132)

    results = []

    for case_number in range(1, CASES + 1):
        result = run_case(case_number)

        if result is not None:
            results.append(result)

    summarize(results)

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


if __name__ == "__main__":
    main()
