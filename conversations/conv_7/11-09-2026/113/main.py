import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 139

CASES = 100

# Width of the candidate interval as a fraction of p.
WIDTH_RATIOS = [
    0.25,
    0.50,
    0.75,
    0.90,
    0.99,
    1.00,
    1.01,
    1.10,
    1.25,
]


def generate_semiprime():
    """
    Generate p < q with q sufficiently larger than p so that
    the synthetic test interval around p never reaches q.

    q/p is approximately 4..6.
    """

    while True:
        p = int(nextprime(random.randint(2000, 5000)))

        ratio = random.uniform(4.0, 6.0)
        q = int(nextprime(int(p * ratio)))

        if q > 3 * p:
            return p, q, p * q


def make_interval(p, width_ratio):
    """
    Construct an interval [L,R] centered approximately around p.

    The interval width is approximately:

        width_ratio * p
    """

    width = max(1, int(round(width_ratio * p)))

    left = width // 2
    right = width - left

    L = max(2, p - left)
    R = p + right

    return L, R


def make_binary_probe(lower, upper):
    """
    Same binary construction as Experiment 138.

        m = R
        midpoint = floor((L+R)/2)
        d = R-midpoint
        k = m+2

    If floor(R/p)=1 throughout the candidate interval and d<p,
    then:

        HIT  <=>  p > midpoint.
    """

    if lower >= upper:
        return None

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


def check_probe_conditions(lower, upper, k, d, q):
    """
    Check the mathematical conditions required by the binary predicate.

    Condition 1:
        d < every possible p
        => d < lower

    Condition 2:
        quotient floor(R/p)=1 for every possible p
        => upper < 2*lower

    Condition 3:
        q is outside the interval and hence cannot be a numerator
        multiple.
    """

    m = k - 2

    if d >= lower:
        return False, "d>=lower"

    if m >= 2 * lower:
        return False, "quotient-not-constant"

    if m >= q:
        return False, "m>=q"

    return True, "OK"


def run_binary_search(n, q, lower, upper):
    """
    Run the adaptive binary search.

    The interval is deliberately supplied externally so this experiment
    measures the width boundary rather than the initial localization step.
    """

    history = []

    invariant_failures = []

    while lower < upper:
        probe = make_binary_probe(lower, upper)

        if probe is None:
            break

        k, d, midpoint = probe

        valid, reason = check_probe_conditions(
            lower,
            upper,
            k,
            d,
            q,
        )

        if not valid:
            invariant_failures.append(
                {
                    "lower": lower,
                    "upper": upper,
                    "k": k,
                    "d": d,
                    "reason": reason,
                }
            )

        g = evaluate_probe(n, k, d)
        signal = classify_signal(g, n)

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
            }
        )

        if signal == "HIT":
            # Intended interpretation:
            #
            #     p > midpoint
            #
            lower = midpoint + 1

        elif signal == "MISS":
            # Intended interpretation:
            #
            #     p <= midpoint
            #
            upper = midpoint

        else:
            # N or some unexpected divisor means the binary predicate
            # is no longer behaving as intended.
            return {
                "success": False,
                "candidate": None,
                "history": history,
                "invariant_failures": invariant_failures,
            }

    candidate = lower if lower == upper else None

    return {
        "success": candidate is not None,
        "candidate": candidate,
        "history": history,
        "invariant_failures": invariant_failures,
    }


def summarize_ratio(results):
    total = len(results)

    successes = sum(
        1 for result in results
        if result["success"]
    )

    invariant_clean = sum(
        1 for result in results
        if len(result["invariant_failures"]) == 0
    )

    total_probes = sum(
        len(result["history"])
        for result in results
    )

    max_probes = max(
        len(result["history"])
        for result in results
    )

    return {
        "cases": total,
        "successes": successes,
        "failures": total - successes,
        "success_rate": successes / total,
        "invariant_clean": invariant_clean,
        "average_probes": total_probes / total,
        "max_probes": max_probes,
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    all_results = {
        ratio: []
        for ratio in WIDTH_RATIOS
    }

    case_data = []

    start = time.perf_counter()

    for case in range(1, CASES + 1):
        p, q, n = generate_semiprime()

        case_data.append(
            {
                "case": case,
                "p": p,
                "q": q,
                "n": n,
            }
        )

        for ratio in WIDTH_RATIOS:
            lower, upper = make_interval(
                p,
                ratio,
            )

            result = run_binary_search(
                n,
                q,
                lower,
                upper,
            )

            all_results[ratio].append(result)

    elapsed = time.perf_counter() - start

    print("SUMMARY BY WIDTH")
    print()

    print(
        "ratio"
        "\tcases"
        "\tsuccesses"
        "\tfailures"
        "\tsuccess_rate"
        "\tinvariant_clean"
        "\tavg_probes"
        "\tmax_probes"
    )

    for ratio in WIDTH_RATIOS:
        summary = summarize_ratio(
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
            f"{summary['max_probes']}"
        )

    print()
    print("BOUNDARY CASE EXAMPLES")
    print()

    # Show the first failure for each width, if one exists.
    for ratio in WIDTH_RATIOS:
        results = all_results[ratio]

        failure_index = None

        for index, result in enumerate(results):
            if not result["success"]:
                failure_index = index
                break

        if failure_index is None:
            print(
                f"ratio={ratio:.2f} "
                f"first_failure=NONE"
            )
            continue

        data = case_data[failure_index]
        result = results[failure_index]

        print(
            f"ratio={ratio:.2f} "
            f"case={data['case']} "
            f"p={data['p']} "
            f"q={data['q']} "
            f"candidate={result['candidate']}"
        )

        if result["history"]:
            first = result["history"][0]

            print(
                f"  initial=["
                f"{first['lower']},"
                f"{first['upper']}"
                f"] "
                f"k={first['k']} "
                f"d={first['d']} "
                f"gcd={first['gcd']} "
                f"signal={first['signal']} "
                f"invariant={first['invariant']}"
            )

    print()
    print("FIRST CASE TRACE")
    print()

    # Detailed trace for one case at every width.
    first_case = case_data[0]

    print(
        f"p={first_case['p']} "
        f"q={first_case['q']} "
        f"N={first_case['n']}"
    )
    print()

    for ratio in WIDTH_RATIOS:
        lower, upper = make_interval(
            first_case["p"],
            ratio,
        )

        result = all_results[ratio][0]

        print(
            f"WIDTH {ratio:.2f}"
            f"  interval=[{lower},{upper}]"
            f"  success={result['success']}"
            f"  candidate={result['candidate']}"
        )

        for index, item in enumerate(
            result["history"],
            1,
        ):
            print(
                f"  {index:02d}: "
                f"[{item['lower']},{item['upper']}] "
                f"mid={item['midpoint']} "
                f"k={item['k']} "
                f"d={item['d']} "
                f"gcd={item['gcd']} "
                f"{item['signal']} "
                f"{item['invariant']}"
            )

        print()

    print(
        f"runtime_seconds={elapsed:.6f}"
    )
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_ID}")


if __name__ == "__main__":
    main()
