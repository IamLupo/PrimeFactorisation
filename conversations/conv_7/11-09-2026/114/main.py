import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 140
CASES = 100

WIDTH_RATIOS = [
    0.10,
    0.25,
    0.40,
    0.49,
    0.50,
    0.60,
    0.75,
    0.90,
    1.00,
]


def generate_semiprime():
    while True:
        p = int(nextprime(random.randint(2000, 5000)))

        # Keep q well above p so q cannot accidentally enter
        # the numerator interval.
        q = int(nextprime(random.randint(4 * p, 6 * p)))

        if q > p:
            return p, q, p * q


def make_hidden_interval(p, width_ratio):
    """
    Construct an interval of approximately width_ratio*p containing p.

    The position of p is randomized inside the interval.

    We explicitly avoid putting p exactly at the midpoint because
    that would leak the answer to the binary search.
    """

    width = max(4, int(round(width_ratio * p)))

    max_left = width - 1

    left = random.randint(1, max_left)

    lower = p - left
    upper = lower + width

    midpoint = (lower + upper) // 2

    if midpoint == p:
        if left < max_left:
            left += 1
        else:
            left -= 1

        lower = p - left
        upper = lower + width

    midpoint = (lower + upper) // 2

    assert lower < p < upper
    assert midpoint != p

    return lower, upper


def make_binary_probe(lower, upper):
    midpoint = (lower + upper) // 2

    m = upper
    d = upper - midpoint
    k = m + 2

    return k, d, midpoint


def evaluate_probe(n, k, d):
    m = k - 2

    if d < 0:
        raise ValueError("negative d")

    if d >= k - 1:
        return n

    return math.gcd(math.comb(m, d), n)


def classify_signal(g, n):
    if g == 1:
        return "MISS"

    if g == n:
        return "N"

    return "HIT"


def check_invariant(lower, upper, d, q):
    """
    Conditions for the clean Lucas threshold:

        d < p
        floor(upper/p) = 1
        upper < q
    """

    if d >= lower:
        return False, "d>=lower"

    if upper >= 2 * lower:
        return False, "quotient-not-constant"

    if upper >= q:
        return False, "upper>=q"

    return True, "OK"


def run_binary_search(n, q, lower, upper, true_p):
    history = []
    invariant_failures = 0

    while lower < upper:
        k, d, midpoint = make_binary_probe(
            lower,
            upper,
        )

        valid, reason = check_invariant(
            lower,
            upper,
            d,
            q,
        )

        if not valid:
            invariant_failures += 1

        g = evaluate_probe(
            n,
            k,
            d,
        )

        signal = classify_signal(
            g,
            n,
        )

        history.append(
            {
                "lower": lower,
                "upper": upper,
                "midpoint": midpoint,
                "k": k,
                "d": d,
                "gcd": g,
                "signal": signal,
                "invariant": reason,
                "mid_is_p": midpoint == true_p,
            }
        )

        if signal == "HIT":
            lower = midpoint + 1

        elif signal == "MISS":
            upper = midpoint

        else:
            return {
                "success": False,
                "candidate": None,
                "history": history,
                "invariant_failures": invariant_failures,
            }

    candidate = lower

    return {
        "success": candidate == true_p,
        "candidate": candidate,
        "history": history,
        "invariant_failures": invariant_failures,
    }


def summarize(results):
    cases = len(results)

    successes = sum(
        result["success"]
        for result in results
    )

    clean = sum(
        result["invariant_failures"] == 0
        for result in results
    )

    total_probes = sum(
        len(result["history"])
        for result in results
    )

    midpoint_hits = sum(
        any(item["mid_is_p"] for item in result["history"])
        for result in results
    )

    return {
        "cases": cases,
        "successes": successes,
        "failures": cases - successes,
        "success_rate": successes / cases,
        "invariant_clean": clean,
        "average_probes": total_probes / cases,
        "max_probes": max(
            len(result["history"])
            for result in results
        ),
        "cases_with_midpoint_p": midpoint_hits,
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    all_results = {
        ratio: []
        for ratio in WIDTH_RATIOS
    }

    first_cases = {}

    start = time.perf_counter()

    for case in range(1, CASES + 1):
        p, q, n = generate_semiprime()

        for ratio in WIDTH_RATIOS:
            lower, upper = make_hidden_interval(
                p,
                ratio,
            )

            result = run_binary_search(
                n,
                q,
                lower,
                upper,
                p,
            )

            all_results[ratio].append(
                result
            )

            if case == 1:
                first_cases[ratio] = (
                    p,
                    q,
                    n,
                    lower,
                    upper,
                    result,
                )

    elapsed = time.perf_counter() - start

    print("SUMMARY")
    print()

    print(
        "ratio\t"
        "cases\t"
        "successes\t"
        "failures\t"
        "success_rate\t"
        "invariant_clean\t"
        "avg_probes\t"
        "max_probes\t"
        "midpoint_p"
    )

    for ratio in WIDTH_RATIOS:
        summary = summarize(
            all_results[ratio]
        )

        print(
            f"{ratio:.2f}\t"
            f"{summary['cases']}\t"
            f"{summary['successes']}\t"
            f"{summary['failures']}\t"
            f"{summary['success_rate']:.4f}\t"
            f"{summary['invariant_clean']}\t"
            f"{summary['average_probes']:.3f}\t"
            f"{summary['max_probes']}\t"
            f"{summary['cases_with_midpoint_p']}"
        )

    print()
    print("FIRST CASE TRACE")
    print()

    for ratio in WIDTH_RATIOS:
        p, q, n, lower, upper, result = first_cases[ratio]

        print(
            f"WIDTH {ratio:.2f} "
            f"p={p} "
            f"q={q} "
            f"interval=[{lower},{upper}] "
            f"candidate={result['candidate']} "
            f"success={result['success']}"
        )

        for i, item in enumerate(
            result["history"],
            1,
        ):
            print(
                f"  {i:02d}: "
                f"[{item['lower']},{item['upper']}] "
                f"mid={item['midpoint']} "
                f"k={item['k']} "
                f"d={item['d']} "
                f"gcd={item['gcd']} "
                f"{item['signal']} "
                f"{item['invariant']}"
            )

        print()

    print(f"runtime_seconds={elapsed:.6f}")
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_ID}")


if __name__ == "__main__":
    main()
