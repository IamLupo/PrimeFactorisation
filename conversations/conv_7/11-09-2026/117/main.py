import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 143

CASES_PER_RANGE = 200

ACTUAL_RATIO_RANGES = [
    (1.05, 1.50),
    (1.50, 2.00),
    (2.00, 3.00),
    (3.00, 4.00),
    (4.00, 5.00),
    (5.00, 6.00),
    (6.00, 8.00),
    (8.00, 10.00),
    (10.00, 16.00),
]

ASSUMED_RATIO = 16.0

MAX_A = 16


def generate_semiprime(min_ratio, max_ratio):
    while True:
        p = int(nextprime(random.randint(2000, 5000)))

        ratio = random.uniform(
            min_ratio,
            max_ratio,
        )

        q = int(nextprime(int(p * ratio)))

        actual_ratio = q / p

        if min_ratio <= actual_ratio <= max_ratio:
            return p, q, p * q


def make_n_only_interval(n, assumed_ratio):
    lower = math.ceil(
        math.sqrt(n / assumed_ratio)
    )

    upper = math.isqrt(n)

    return lower, upper


def evaluate_binomial_gcd(n, m, d):
    if d < 0 or d > m:
        return n

    return math.gcd(
        math.comb(m, d),
        n,
    )


def classify_gcd(g, n, p, q):
    if g == 1:
        return "MISS"

    if g == p:
        return "P"

    if g == q:
        return "Q"

    if g == n:
        return "N"

    return "OTHER"


def quotient_of_p(m, p):
    return m // p


def run_sweep(n, p, q, lower, upper):
    """
    Use the midpoint of the N-only interval as the threshold.

    For each possible Lucas quotient a:

        d = m - a*t

    where:

        m = upper
        t = midpoint

    We only test d in the useful range:

        1 <= d < lower

    so that d is guaranteed smaller than every candidate p.
    """

    m = upper
    threshold = (lower + upper) // 2

    actual_a = quotient_of_p(
        m,
        p,
    )

    history = []

    for a in range(1, MAX_A + 1):
        d = m - a * threshold

        if d <= 0:
            continue

        if d >= lower:
            continue

        g = evaluate_binomial_gcd(
            n,
            m,
            d,
        )

        signal = classify_gcd(
            g,
            n,
            p,
            q,
        )

        history.append(
            {
                "a": a,
                "m": m,
                "threshold": threshold,
                "d": d,
                "gcd": g,
                "signal": signal,
                "actual_a": actual_a,
            }
        )

        if signal == "P":
            return {
                "success": True,
                "candidate": p,
                "reason": "direct-p",
                "actual_a": actual_a,
                "history": history,
            }

        if signal == "Q":
            return {
                "success": True,
                "candidate": q,
                "reason": "direct-q",
                "actual_a": actual_a,
                "history": history,
            }

    return {
        "success": False,
        "candidate": None,
        "reason": "no-factor",
        "actual_a": actual_a,
        "history": history,
    }


def summarize(results):
    cases = len(results)

    successes = sum(
        result["success"]
        for result in results
    )

    direct_p = sum(
        result["reason"] == "direct-p"
        for result in results
    )

    direct_q = sum(
        result["reason"] == "direct-q"
        for result in results
    )

    a_match = sum(
        result["reason"] == "direct-p"
        and any(
            step["a"] == result["actual_a"]
            and step["signal"] == "P"
            for step in result["history"]
        )
        for result in results
    )

    tested = [
        len(result["history"])
        for result in results
    ]

    return {
        "cases": cases,
        "successes": successes,
        "failures": cases - successes,
        "success_rate": successes / cases,
        "direct_p": direct_p,
        "direct_q": direct_q,
        "actual_a_hits": a_match,
        "avg_probes": (
            sum(tested) / len(tested)
            if tested
            else 0.0
        ),
        "max_probes": (
            max(tested)
            if tested
            else 0
        ),
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    start = time.perf_counter()

    all_groups = []

    for min_ratio, max_ratio in ACTUAL_RATIO_RANGES:
        results = []

        for _ in range(CASES_PER_RANGE):
            p, q, n = generate_semiprime(
                min_ratio,
                max_ratio,
            )

            lower, upper = make_n_only_interval(
                n,
                ASSUMED_RATIO,
            )

            if not (lower <= p <= upper):
                results.append(
                    {
                        "success": False,
                        "candidate": None,
                        "reason": "p-outside",
                        "actual_a": upper // p,
                        "history": [],
                    }
                )

                continue

            result = run_sweep(
                n,
                p,
                q,
                lower,
                upper,
            )

            results.append(result)

        all_groups.append(
            (
                min_ratio,
                max_ratio,
                results,
            )
        )

    print("SUMMARY")
    print()

    print(
        "actual_ratio\t"
        "cases\t"
        "successes\t"
        "failures\t"
        "success_rate\t"
        "direct_p\t"
        "direct_q\t"
        "actual_a_hits\t"
        "avg_probes\t"
        "max_probes"
    )

    for min_ratio, max_ratio, results in all_groups:
        summary = summarize(results)

        print(
            f"{min_ratio:.2f}-{max_ratio:.2f}\t"
            f"{summary['cases']}\t"
            f"{summary['successes']}\t"
            f"{summary['failures']}\t"
            f"{summary['success_rate']:.4f}\t"
            f"{summary['direct_p']}\t"
            f"{summary['direct_q']}\t"
            f"{summary['actual_a_hits']}\t"
            f"{summary['avg_probes']:.3f}\t"
            f"{summary['max_probes']}"
        )

    print()
    print("SUCCESS RATE BY ACTUAL QUOTIENT")
    print()

    quotient_results = {}

    for _, _, results in all_groups:
        for result in results:
            a = result["actual_a"]

            if a not in quotient_results:
                quotient_results[a] = []

            quotient_results[a].append(result)

    print(
        "a\t"
        "cases\t"
        "successes\t"
        "success_rate"
    )

    for a in sorted(quotient_results):
        results = quotient_results[a]

        successes = sum(
            result["success"]
            for result in results
        )

        print(
            f"{a}\t"
            f"{len(results)}\t"
            f"{successes}\t"
            f"{successes / len(results):.4f}"
        )

    print()
    print("FIRST SUCCESS EXAMPLES")
    print()

    shown = 0

    for min_ratio, max_ratio, results in all_groups:
        for result in results:
            if not result["success"]:
                continue

            print(
                f"actual_ratio={min_ratio:.2f}-{max_ratio:.2f} "
                f"actual_a={result['actual_a']} "
                f"reason={result['reason']}"
            )

            for index, step in enumerate(
                result["history"],
                1,
            ):
                print(
                    f"  {index:02d}: "
                    f"a={step['a']} "
                    f"m={step['m']} "
                    f"t={step['threshold']} "
                    f"d={step['d']} "
                    f"gcd={step['gcd']} "
                    f"{step['signal']}"
                )

            print()

            shown += 1

            if shown >= 10:
                break

        if shown >= 10:
            break

    elapsed = time.perf_counter() - start

    print(
        f"runtime_seconds={elapsed:.6f}"
    )
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_ID}")


if __name__ == "__main__":
    main()
