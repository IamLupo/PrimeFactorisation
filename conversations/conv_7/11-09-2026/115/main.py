import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 141

CASES_PER_RANGE = 100

ASSUMED_RATIOS = [
    1.5,
    2.0,
    3.0,
    4.0,
    5.0,
    6.0,
    8.0,
    10.0,
    16.0,
]

ACTUAL_RATIO_RANGES = [
    (1.05, 1.50),
    (1.50, 2.00),
    (2.00, 3.00),
    (3.00, 4.00),
    (4.00, 6.00),
    (6.00, 10.00),
]


def generate_semiprime(min_ratio, max_ratio):
    while True:
        p = int(nextprime(random.randint(2000, 5000)))

        ratio = random.uniform(
            min_ratio,
            max_ratio,
        )

        q = int(nextprime(int(p * ratio)))

        actual_ratio = q / p

        if (
            q > p
            and min_ratio <= actual_ratio <= max_ratio
        ):
            return p, q, p * q


def make_n_only_interval(n, assumed_ratio):
    """
    Construct the initial p-interval using only N.

    If q/p <= assumed_ratio, then:

        p >= sqrt(N / assumed_ratio)

    and always:

        p <= sqrt(N).
    """

    lower = math.ceil(
        math.sqrt(n / assumed_ratio)
    )

    upper = math.isqrt(n)

    return lower, upper


def make_binary_probe(lower, upper):
    """
    Experiment 138 construction.

        m = R
        d = R - midpoint
        k = m + 2
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

    return math.gcd(
        math.comb(m, d),
        n,
    )


def classify_signal(g, n):
    if g == 1:
        return "MISS"

    if g == n:
        return "N"

    return "HIT"


def probe_invariant(lower, upper, d, q):
    """
    Check the clean mathematical regime.

        d < p
        upper < 2p
        upper < q

    Since p is deliberately hidden from the algorithm, these checks
    are diagnostics only.
    """

    reasons = []

    # These are written as symbolic interval checks where possible.
    if d >= lower:
        reasons.append("d>=L")

    if upper >= 2 * lower:
        reasons.append("quotient-not-guaranteed")

    if upper >= q:
        reasons.append("upper>=q")

    if reasons:
        return False, ",".join(reasons)

    return True, "OK"


def run_binary_search(n, q, true_p, lower, upper):
    history = []

    while lower < upper:
        probe = make_binary_probe(
            lower,
            upper,
        )

        if probe is None:
            break

        k, d, midpoint = probe

        invariant, invariant_reason = probe_invariant(
            lower,
            upper,
            d,
            q,
        )

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
                "invariant": invariant_reason,
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
                "reason": "N",
            }

    candidate = lower

    if candidate == true_p:
        return {
            "success": True,
            "candidate": candidate,
            "history": history,
            "reason": "OK",
        }

    return {
        "success": False,
        "candidate": candidate,
        "history": history,
        "reason": "wrong-factor",
    }


def run_experiment(actual_ratio_min, actual_ratio_max):
    results = []

    for _ in range(CASES_PER_RANGE):
        p, q, n = generate_semiprime(
            actual_ratio_min,
            actual_ratio_max,
        )

        for assumed_ratio in ASSUMED_RATIOS:
            lower, upper = make_n_only_interval(
                n,
                assumed_ratio,
            )

            contains_p = lower <= p <= upper

            if not contains_p:
                result = {
                    "success": False,
                    "candidate": None,
                    "history": [],
                    "reason": "p-outside-initial-interval",
                }
            else:
                result = run_binary_search(
                    n,
                    q,
                    p,
                    lower,
                    upper,
                )

            results.append(
                {
                    "p": p,
                    "q": q,
                    "n": n,
                    "actual_ratio": q / p,
                    "assumed_ratio": assumed_ratio,
                    "lower": lower,
                    "upper": upper,
                    "contains_p": contains_p,
                    "result": result,
                }
            )

    return results


def summarize(results):
    cases = len(results)

    contained = sum(
        item["contains_p"]
        for item in results
    )

    successes = sum(
        item["result"]["success"]
        for item in results
    )

    wrong_interval = sum(
        item["result"]["reason"]
        == "p-outside-initial-interval"
        for item in results
    )

    n_failures = sum(
        item["result"]["reason"] == "N"
        for item in results
    )

    wrong_factor = sum(
        item["result"]["reason"]
        == "wrong-factor"
        for item in results
    )

    probes = [
        len(item["result"]["history"])
        for item in results
    ]

    return {
        "cases": cases,
        "contained": contained,
        "containment_rate": contained / cases,
        "successes": successes,
        "failures": cases - successes,
        "success_rate": successes / cases,
        "interval_misses": wrong_interval,
        "N_failures": n_failures,
        "wrong_factor": wrong_factor,
        "average_probes": (
            sum(probes) / cases
        ),
        "max_probes": max(probes),
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    start = time.perf_counter()

    all_groups = []

    for actual_min, actual_max in ACTUAL_RATIO_RANGES:
        results = run_experiment(
            actual_min,
            actual_max,
        )

        all_groups.append(
            (
                actual_min,
                actual_max,
                results,
            )
        )

    print("SUMMARY")
    print()

    print(
        "actual_ratio"
        "\tassumed_C"
        "\tcases"
        "\tcontained"
        "\tsuccesses"
        "\tfailures"
        "\tcontainment_rate"
        "\tsuccess_rate"
        "\tinterval_misses"
        "\tN_failures"
        "\twrong_factor"
        "\tavg_probes"
        "\tmax_probes"
    )

    for actual_min, actual_max, results in all_groups:
        for assumed_ratio in ASSUMED_RATIOS:
            subset = [
                item
                for item in results
                if item["assumed_ratio"]
                == assumed_ratio
            ]

            summary = summarize(subset)

            print(
                f"{actual_min:.2f}-{actual_max:.2f}\t"
                f"{assumed_ratio:.1f}\t"
                f"{summary['cases']}\t"
                f"{summary['contained']}\t"
                f"{summary['successes']}\t"
                f"{summary['failures']}\t"
                f"{summary['containment_rate']:.4f}\t"
                f"{summary['success_rate']:.4f}\t"
                f"{summary['interval_misses']}\t"
                f"{summary['N_failures']}\t"
                f"{summary['wrong_factor']}\t"
                f"{summary['average_probes']:.3f}\t"
                f"{summary['max_probes']}"
            )

    print()
    print("FIRST SUCCESSFUL N-ONLY TRACE")
    print()

    found = False

    for actual_min, actual_max, results in all_groups:
        for item in results:
            if (
                item["result"]["success"]
                and item["assumed_ratio"] == 2.0
            ):
                print(
                    f"actual_ratio={item['actual_ratio']:.6f}"
                )
                print(
                    f"p={item['p']}"
                )
                print(
                    f"q={item['q']}"
                )
                print(
                    f"N={item['n']}"
                )
                print(
                    f"initial=["
                    f"{item['lower']},"
                    f"{item['upper']}"
                    f"]"
                )

                for i, step in enumerate(
                    item["result"]["history"],
                    1,
                ):
                    print(
                        f"  {i:02d}: "
                        f"[{step['lower']},"
                        f"{step['upper']}] "
                        f"mid={step['midpoint']} "
                        f"k={step['k']} "
                        f"d={step['d']} "
                        f"gcd={step['gcd']} "
                        f"{step['signal']} "
                        f"{step['invariant']}"
                    )

                found = True
                break

        if found:
            break

    if not found:
        print("NONE")

    elapsed = time.perf_counter() - start

    print()
    print(
        f"runtime_seconds={elapsed:.6f}"
    )
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_ID}")


if __name__ == "__main__":
    main()
