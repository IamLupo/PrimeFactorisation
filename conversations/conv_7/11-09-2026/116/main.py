import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 142

CASES_PER_RANGE = 200

ACTUAL_RATIO_RANGES = [
    (3.50, 3.75),
    (3.75, 4.00),
    (4.00, 4.10),
    (4.10, 4.25),
    (4.25, 4.50),
    (4.50, 4.75),
    (4.75, 5.00),
]

ASSUMED_RATIOS = [
    3.50,
    3.75,
    3.90,
    3.99,
    4.00,
    4.01,
    4.10,
    4.20,
    4.33,
    4.50,
    4.75,
    5.00,
]


def generate_semiprime(min_ratio, max_ratio):
    while True:
        p = int(nextprime(random.randint(2000, 5000)))

        ratio = random.uniform(min_ratio, max_ratio)

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


def make_probe(lower, upper):
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


def classify_gcd(g, n, p, q):
    if g == 1:
        return "MISS"

    if g == n:
        return "N"

    if g == p:
        return "P"

    if g == q:
        return "Q"

    return "OTHER"


def theoretical_state(lower, upper, p, q, d):
    """
    Diagnostic only.

    Safe binary regime requires:

        d < p
        upper < 2p
        upper < q
    """

    if d >= p:
        return "D_FAIL"

    if upper >= 2 * p:
        return "QUOTIENT_FAIL"

    if upper >= q:
        return "Q_IN_RANGE"

    return "SAFE"


def run_search(n, p, q, lower, upper):
    history = []

    while lower < upper:
        k, d, midpoint = make_probe(
            lower,
            upper,
        )

        g = evaluate_probe(
            n,
            k,
            d,
        )

        signal = classify_gcd(
            g,
            n,
            p,
            q,
        )

        theory = theoretical_state(
            lower,
            upper,
            p,
            q,
            d,
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
                "theory": theory,
            }
        )

        # Immediate factor recovery.
        if signal == "P":
            return {
                "success": True,
                "candidate": p,
                "reason": "direct-p",
                "history": history,
            }

        if signal == "Q":
            return {
                "success": True,
                "candidate": q,
                "reason": "direct-q",
                "history": history,
            }

        if signal != "MISS":
            return {
                "success": False,
                "candidate": None,
                "reason": signal,
                "history": history,
            }

        # IMPORTANT:
        #
        # We only trust MISS as a binary decision when the probe is
        # theoretically safe.
        #
        # Once outside the safe regime we stop and classify the
        # bootstrap attempt instead of pretending the binary predicate
        # is valid.
        if theory != "SAFE":
            return {
                "success": False,
                "candidate": None,
                "reason": "unsafe-miss",
                "history": history,
            }

        upper = midpoint

    if lower == p:
        return {
            "success": True,
            "candidate": lower,
            "reason": "binary-p",
            "history": history,
        }

    if lower == q:
        return {
            "success": True,
            "candidate": lower,
            "reason": "binary-q",
            "history": history,
        }

    return {
        "success": False,
        "candidate": lower,
        "reason": "wrong-candidate",
        "history": history,
    }


def run_group(actual_min, actual_max):
    results = []

    for _ in range(CASES_PER_RANGE):
        p, q, n = generate_semiprime(
            actual_min,
            actual_max,
        )

        for assumed in ASSUMED_RATIOS:
            lower, upper = make_n_only_interval(
                n,
                assumed,
            )

            contains_p = lower <= p <= upper

            if not contains_p:
                results.append(
                    {
                        "p": p,
                        "q": q,
                        "n": n,
                        "actual_ratio": q / p,
                        "assumed": assumed,
                        "contains_p": False,
                        "result": None,
                    }
                )

                continue

            result = run_search(
                n,
                p,
                q,
                lower,
                upper,
            )

            results.append(
                {
                    "p": p,
                    "q": q,
                    "n": n,
                    "actual_ratio": q / p,
                    "assumed": assumed,
                    "contains_p": True,
                    "result": result,
                }
            )

    return results


def summarize(items):
    total = len(items)

    contained = sum(
        item["contains_p"]
        for item in items
    )

    direct_p = 0
    direct_q = 0
    binary_success = 0
    unsafe_miss = 0
    wrong = 0

    first_gcd_nontrivial = 0
    first_gcd_n = 0

    probes = []

    for item in items:
        result = item["result"]

        if result is None:
            continue

        if result["success"]:
            if result["reason"] == "direct-p":
                direct_p += 1
            elif result["reason"] == "direct-q":
                direct_q += 1
            else:
                binary_success += 1

        elif result["reason"] == "unsafe-miss":
            unsafe_miss += 1
        else:
            wrong += 1

        if result["history"]:
            first = result["history"][0]

            if first["signal"] in ("P", "Q", "OTHER"):
                first_gcd_nontrivial += 1

            if first["signal"] == "N":
                first_gcd_n += 1

            probes.append(len(result["history"]))

    return {
        "cases": total,
        "contained": contained,
        "containment_rate": contained / total,
        "direct_p": direct_p,
        "direct_q": direct_q,
        "binary_success": binary_success,
        "successes": direct_p + direct_q + binary_success,
        "unsafe_miss": unsafe_miss,
        "wrong": wrong,
        "first_nontrivial": first_gcd_nontrivial,
        "first_N": first_gcd_n,
        "avg_probes": (
            sum(probes) / len(probes)
            if probes
            else 0.0
        ),
        "max_probes": (
            max(probes)
            if probes
            else 0
        ),
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    start = time.perf_counter()

    groups = []

    for actual_min, actual_max in ACTUAL_RATIO_RANGES:
        results = run_group(
            actual_min,
            actual_max,
        )

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
        "C\t"
        "cases\t"
        "contained\t"
        "successes\t"
        "direct_p\t"
        "direct_q\t"
        "binary\t"
        "unsafe_miss\t"
        "wrong\t"
        "first_nontrivial\t"
        "first_N\t"
        "avg_probes\t"
        "max_probes"
    )

    for actual_min, actual_max, results in groups:
        for assumed in ASSUMED_RATIOS:
            subset = [
                item
                for item in results
                if item["assumed"] == assumed
            ]

            s = summarize(subset)

            print(
                f"{actual_min:.2f}-{actual_max:.2f}\t"
                f"{assumed:.2f}\t"
                f"{s['cases']}\t"
                f"{s['contained']}\t"
                f"{s['successes']}\t"
                f"{s['direct_p']}\t"
                f"{s['direct_q']}\t"
                f"{s['binary_success']}\t"
                f"{s['unsafe_miss']}\t"
                f"{s['wrong']}\t"
                f"{s['first_nontrivial']}\t"
                f"{s['first_N']}\t"
                f"{s['avg_probes']:.3f}\t"
                f"{s['max_probes']}"
            )

    print()
    print("FIRST INVALID-PROBE TRACES")
    print()

    shown = 0

    for actual_min, actual_max, results in groups:
        for item in results:
            result = item["result"]

            if result is None:
                continue

            has_unsafe = any(
                step["theory"] != "SAFE"
                for step in result["history"]
            )

            if not has_unsafe:
                continue

            print(
                f"actual_ratio={item['actual_ratio']:.6f} "
                f"C={item['assumed']:.2f} "
                f"p={item['p']} "
                f"q={item['q']} "
                f"N={item['n']} "
                f"result={result['reason']}"
            )

            for i, step in enumerate(
                result["history"],
                1,
            ):
                print(
                    f"  {i:02d}: "
                    f"[{step['lower']},"
                    f"{step['upper']}] "
                    f"mid={step['midpoint']} "
                    f"d={step['d']} "
                    f"gcd={step['gcd']} "
                    f"{step['signal']} "
                    f"{step['theory']}"
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
