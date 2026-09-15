import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 144

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

MAX_A = 32


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


def find_feasible_m(lower, upper, a):
    """
    We need:

        a*upper <= m < (a+1)*lower.

    Therefore the smallest possible m is:

        m = a*upper.

    This is feasible if:

        a*upper < (a+1)*lower.
    """

    minimum_m = a * upper
    maximum_m = (a + 1) * lower - 1

    if minimum_m > maximum_m:
        return None

    return minimum_m


def choose_threshold(lower, upper):
    return (lower + upper) // 2


def make_probe(lower, upper, a):
    m = find_feasible_m(
        lower,
        upper,
        a,
    )

    if m is None:
        return None

    threshold = choose_threshold(
        lower,
        upper,
    )

    d = m - a * threshold

    if d <= 0:
        return None

    if d >= lower:
        return None

    k = m + 2

    return {
        "a": a,
        "m": m,
        "k": k,
        "d": d,
        "threshold": threshold,
    }


def evaluate_probe(n, m, d):
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


def actual_quotient(m, p):
    return m // p


def verify_probe(probe, lower, upper, p, q):
    a = probe["a"]
    m = probe["m"]
    d = probe["d"]

    actual_a = actual_quotient(
        m,
        p,
    )

    conditions = {
        "quotient_matches": actual_a == a,
        "d_lt_p": d < p,
        "quotient_all_interval": (
            a * upper <= m
            and m < (a + 1) * lower
        ),
        "q_outside": m < q,
    }

    return conditions


def run_case(n, p, q, lower, upper):
    probes = []
    feasible_count = 0

    # First determine which quotient constructions are
    # mathematically possible from N alone.
    for a in range(1, MAX_A + 1):
        m = find_feasible_m(
            lower,
            upper,
            a,
        )

        if m is not None:
            feasible_count += 1

        probe = make_probe(
            lower,
            upper,
            a,
        )

        if probe is None:
            continue

        conditions = verify_probe(
            probe,
            lower,
            upper,
            p,
            q,
        )

        g = evaluate_probe(
            n,
            probe["m"],
            probe["d"],
        )

        signal = classify_gcd(
            g,
            n,
            p,
            q,
        )

        probes.append(
            {
                **probe,
                "gcd": g,
                "signal": signal,
                "conditions": conditions,
            }
        )

        if signal == "P":
            return {
                "success": True,
                "reason": "direct-p",
                "candidate": p,
                "feasible_count": feasible_count,
                "probes": probes,
            }

        if signal == "Q":
            return {
                "success": True,
                "reason": "direct-q",
                "candidate": q,
                "feasible_count": feasible_count,
                "probes": probes,
            }

    return {
        "success": False,
        "reason": "no-factor",
        "candidate": None,
        "feasible_count": feasible_count,
        "probes": probes,
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

    total_feasible = sum(
        result["feasible_count"]
        for result in results
    )

    actual_a_matches = 0
    valid_constructive_hits = 0

    for result in results:
        for probe in result["probes"]:
            if (
                probe["conditions"]["quotient_matches"]
                and probe["conditions"]["d_lt_p"]
            ):
                actual_a_matches += 1

                if probe["signal"] == "P":
                    valid_constructive_hits += 1

    return {
        "cases": cases,
        "successes": successes,
        "failures": cases - successes,
        "success_rate": successes / cases,
        "direct_p": direct_p,
        "direct_q": direct_q,
        "avg_feasible": (
            total_feasible / cases
        ),
        "constructive_matches": actual_a_matches,
        "constructive_hits": valid_constructive_hits,
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    start = time.perf_counter()

    groups = []

    for actual_min, actual_max in ACTUAL_RATIO_RANGES:
        results = []

        for _ in range(CASES_PER_RANGE):
            p, q, n = generate_semiprime(
                actual_min,
                actual_max,
            )

            lower, upper = make_n_only_interval(
                n,
                ASSUMED_RATIO,
            )

            if not (lower <= p <= upper):
                results.append(
                    {
                        "success": False,
                        "reason": "p-outside",
                        "candidate": None,
                        "feasible_count": 0,
                        "probes": [],
                    }
                )

                continue

            result = run_case(
                n,
                p,
                q,
                lower,
                upper,
            )

            results.append(result)

        groups.append(
            (
                actual_min,
                actual_max,
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
        "avg_feasible_a\t"
        "constructive_matches\t"
        "constructive_hits"
    )

    for actual_min, actual_max, results in groups:
        summary = summarize(results)

        print(
            f"{actual_min:.2f}-{actual_max:.2f}\t"
            f"{summary['cases']}\t"
            f"{summary['successes']}\t"
            f"{summary['failures']}\t"
            f"{summary['success_rate']:.4f}\t"
            f"{summary['direct_p']}\t"
            f"{summary['direct_q']}\t"
            f"{summary['avg_feasible']:.3f}\t"
            f"{summary['constructive_matches']}\t"
            f"{summary['constructive_hits']}"
        )

    print()
    print("FEASIBILITY BY INTERVAL RATIO")
    print()

    print(
        "actual_ratio\t"
        "cases\t"
        "with_feasible_a\t"
        "average_feasible_a"
    )

    for actual_min, actual_max, results in groups:
        usable = sum(
            result["feasible_count"] > 0
            for result in results
        )

        average = (
            sum(
                result["feasible_count"]
                for result in results
            )
            / len(results)
        )

        print(
            f"{actual_min:.2f}-{actual_max:.2f}\t"
            f"{len(results)}\t"
            f"{usable}\t"
            f"{average:.3f}"
        )

    print()
    print("EXAMPLE SUCCESS TRACES")
    print()

    shown = 0

    for actual_min, actual_max, results in groups:
        for result in results:
            if not result["success"]:
                continue

            print(
                f"actual={actual_min:.2f}-"
                f"{actual_max:.2f} "
                f"reason={result['reason']} "
                f"candidate={result['candidate']} "
                f"feasible_a={result['feasible_count']}"
            )

            for i, probe in enumerate(
                result["probes"],
                1,
            ):
                print(
                    f"  {i:02d}: "
                    f"a={probe['a']} "
                    f"m={probe['m']} "
                    f"t={probe['threshold']} "
                    f"d={probe['d']} "
                    f"gcd={probe['gcd']} "
                    f"{probe['signal']} "
                    f"actual_a="
                    f"{probe['conditions']['quotient_matches']} "
                    f"d<p="
                    f"{probe['conditions']['d_lt_p']}"
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
