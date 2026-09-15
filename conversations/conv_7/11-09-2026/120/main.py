import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 146

CASES_PER_TIER = 10

# Target magnitude of N.
N_TIERS = [
    10**10,
    10**12,
    10**14,
    10**16,
    10**18,
]

# Keep the factors reasonably balanced so the experiment measures
# size scaling rather than deliberately extreme factor gaps.
RATIO_MIN = 1.05
RATIO_MAX = 1.50

# Same safe multiplicative segmentation used in Experiment 145.
ALPHA = 1.50


def generate_semiprime_near(target_n):
    """
    Generate p*q approximately target_n.

    q/p is constrained to [RATIO_MIN, RATIO_MAX].
    """

    target = math.sqrt(target_n)

    while True:
        ratio = random.uniform(
            RATIO_MIN,
            RATIO_MAX,
        )

        p_estimate = int(
            target / math.sqrt(ratio)
        )

        # Small random displacement prevents every case from
        # having nearly identical p.
        displacement = max(
            100,
            p_estimate // 50,
        )

        p_start = max(
            2,
            p_estimate - displacement,
        )

        p = int(
            nextprime(
                random.randint(
                    p_start,
                    p_estimate + displacement,
                )
            )
        )

        q_estimate = int(
            target_n / p
        )

        q = int(
            nextprime(q_estimate)
        )

        if q <= p:
            continue

        actual_n = p * q
        actual_ratio = q / p

        # Keep the generated N near the requested tier.
        if not (
            0.75 * target_n
            <= actual_n
            <= 1.50 * target_n
        ):
            continue

        if not (
            RATIO_MIN
            <= actual_ratio
            <= RATIO_MAX
        ):
            continue

        return p, q, actual_n


def make_initial_interval(n):
    """
    Use only N.

    For the true ratio range q/p <= 1.5:

        p >= sqrt(N / 1.5)

    and always:

        p <= sqrt(N).
    """

    lower = math.ceil(
        math.sqrt(n / RATIO_MAX)
    )

    upper = math.isqrt(n)

    return lower, upper


def subdivide_interval(lower, upper):
    """
    Build safe intervals satisfying:

        U <= ALPHA * L

    with ALPHA = 1.5 < 2.
    """

    intervals = []

    current = lower

    while current <= upper:
        segment_upper = min(
            upper,
            math.floor(ALPHA * current),
        )

        if segment_upper <= current:
            segment_upper = current + 1

        intervals.append(
            (
                current,
                segment_upper,
            )
        )

        if segment_upper >= upper:
            break

        current = segment_upper + 1

    return intervals


def make_probe(lower, upper):
    midpoint = (lower + upper) // 2

    m = upper
    d = upper - midpoint

    return {
        "m": m,
        "d": d,
        "k": m + 2,
        "midpoint": midpoint,
    }


def binomial_gcd_product(n, m, d):
    """
    Compute gcd(C(m,d), n) without constructing C(m,d).

    Since d < p and N=p*q:

        gcd(C(m,d), N)
        =
        gcd(product(m-d+1 ... m), N).

    The product is maintained modulo N.

    We also stop immediately when a non-trivial gcd appears.
    """

    if d < 0 or d > m:
        return n, 0

    if d == 0:
        return 1, 0

    accumulator = 1

    start = m - d + 1

    for i in range(start, m + 1):
        accumulator = (
            accumulator * i
        ) % n

        g = math.gcd(
            accumulator,
            n,
        )

        if 1 < g < n:
            return g, i - start + 1

        if g == n:
            return n, i - start + 1

    return math.gcd(
        accumulator,
        n,
    ), d


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


def binary_search_segment(
    n,
    p,
    q,
    lower,
    upper,
):
    """
    Search one mathematically safe segment.
    """

    probes = 0
    total_product_steps = 0
    max_d = 0

    while lower < upper:
        probe = make_probe(
            lower,
            upper,
        )

        m = probe["m"]
        d = probe["d"]
        midpoint = probe["midpoint"]

        # Safe-segment invariants.
        if not (
            d < lower
            and m < 2 * lower
            and m < q
        ):
            return {
                "success": False,
                "candidate": None,
                "reason": "unsafe-segment",
                "probes": probes,
                "product_steps": total_product_steps,
                "max_d": max_d,
            }

        probes += 1
        max_d = max(
            max_d,
            d,
        )

        g, used_steps = binomial_gcd_product(
            n,
            m,
            d,
        )

        total_product_steps += used_steps

        signal = classify_gcd(
            g,
            n,
            p,
            q,
        )

        if signal == "P":
            return {
                "success": True,
                "candidate": p,
                "reason": "P",
                "probes": probes,
                "product_steps": total_product_steps,
                "max_d": max_d,
            }

        if signal == "Q":
            return {
                "success": True,
                "candidate": q,
                "reason": "Q",
                "probes": probes,
                "product_steps": total_product_steps,
                "max_d": max_d,
            }

        if signal == "N":
            return {
                "success": False,
                "candidate": None,
                "reason": "N",
                "probes": probes,
                "product_steps": total_product_steps,
                "max_d": max_d,
            }

        if signal == "OTHER":
            return {
                "success": False,
                "candidate": None,
                "reason": "OTHER",
                "probes": probes,
                "product_steps": total_product_steps,
                "max_d": max_d,
            }

        # MISS:
        #
        # p <= midpoint.
        upper = midpoint

    candidate = lower

    if candidate == p:
        return {
            "success": True,
            "candidate": candidate,
            "reason": "binary",
            "probes": probes,
            "product_steps": total_product_steps,
            "max_d": max_d,
        }

    return {
        "success": False,
        "candidate": candidate,
        "reason": "wrong-segment",
        "probes": probes,
        "product_steps": total_product_steps,
        "max_d": max_d,
    }


def factor_from_n_only(n, p, q):
    """
    Run the full Experiment 145 algorithm.
    """

    lower, upper = make_initial_interval(
        n
    )

    segments = subdivide_interval(
        lower,
        upper,
    )

    total_probes = 0
    total_product_steps = 0
    max_d = 0
    segments_tested = 0

    start = time.perf_counter()

    for segment_lower, segment_upper in segments:
        segments_tested += 1

        result = binary_search_segment(
            n,
            p,
            q,
            segment_lower,
            segment_upper,
        )

        total_probes += result["probes"]
        total_product_steps += (
            result["product_steps"]
        )

        max_d = max(
            max_d,
            result["max_d"],
        )

        if result["success"]:
            elapsed = (
                time.perf_counter()
                - start
            )

            return {
                "success": True,
                "candidate": result["candidate"],
                "segments": segments_tested,
                "probes": total_probes,
                "product_steps": total_product_steps,
                "max_d": max_d,
                "runtime": elapsed,
                "reason": result["reason"],
                "initial_lower": lower,
                "initial_upper": upper,
                "total_segments": len(segments),
            }

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "success": False,
        "candidate": None,
        "segments": segments_tested,
        "probes": total_probes,
        "product_steps": total_product_steps,
        "max_d": max_d,
        "runtime": elapsed,
        "reason": "no-factor",
        "initial_lower": lower,
        "initial_upper": upper,
        "total_segments": len(segments),
    }


def run_tier(target_n):
    results = []

    for case in range(
        1,
        CASES_PER_TIER + 1,
    ):
        p, q, n = (
            generate_semiprime_near(
                target_n
            )
        )

        result = factor_from_n_only(
            n,
            p,
            q,
        )

        result["case"] = case
        result["p"] = p
        result["q"] = q
        result["n"] = n
        result["target_n"] = target_n

        result["ratio"] = q / p

        results.append(result)

    return results


def summarize(results):
    successes = sum(
        result["success"]
        for result in results
    )

    runtimes = [
        result["runtime"]
        for result in results
    ]

    probes = [
        result["probes"]
        for result in results
    ]

    product_steps = [
        result["product_steps"]
        for result in results
    ]

    max_ds = [
        result["max_d"]
        for result in results
    ]

    segments = [
        result["segments"]
        for result in results
    ]

    return {
        "cases": len(results),
        "successes": successes,
        "failures": (
            len(results) - successes
        ),
        "success_rate": (
            successes / len(results)
        ),
        "avg_runtime": (
            sum(runtimes)
            / len(runtimes)
        ),
        "max_runtime": max(runtimes),
        "avg_probes": (
            sum(probes)
            / len(probes)
        ),
        "max_probes": max(probes),
        "avg_product_steps": (
            sum(product_steps)
            / len(product_steps)
        ),
        "max_product_steps": (
            max(product_steps)
        ),
        "avg_max_d": (
            sum(max_ds)
            / len(max_ds)
        ),
        "max_d": max(max_ds),
        "avg_segments": (
            sum(segments)
            / len(segments)
        ),
        "max_segments": max(segments),
    }


def print_slowest(results):
    slowest = max(
        results,
        key=lambda result: result["runtime"],
    )

    print(
        f"slowest_case="
        f"{slowest['case']}"
    )
    print(
        f"  N={slowest['n']}"
    )
    print(
        f"  p={slowest['p']}"
    )
    print(
        f"  q={slowest['q']}"
    )
    print(
        f"  ratio={slowest['ratio']:.6f}"
    )
    print(
        f"  runtime={slowest['runtime']:.6f}s"
    )
    print(
        f"  segments={slowest['segments']}"
    )
    print(
        f"  probes={slowest['probes']}"
    )
    print(
        f"  product_steps="
        f"{slowest['product_steps']}"
    )
    print(
        f"  max_d={slowest['max_d']}"
    )
    print()


def main():
    print(
        f"START EXPERIMENT {EXPERIMENT_ID}"
    )
    print()

    random.seed(EXPERIMENT_ID)

    overall_start = time.perf_counter()

    for target_n in N_TIERS:
        print(
            f"TIER N~{target_n}"
        )
        print()

        results = run_tier(
            target_n
        )

        summary = summarize(
            results
        )

        print(
            "cases\t"
            "successes\t"
            "failures\t"
            "success_rate\t"
            "avg_runtime\t"
            "max_runtime\t"
            "avg_segments\t"
            "max_segments\t"
            "avg_probes\t"
            "max_probes\t"
            "avg_product_steps\t"
            "max_product_steps\t"
            "avg_max_d\t"
            "max_d"
        )

        print(
            f"{summary['cases']}\t"
            f"{summary['successes']}\t"
            f"{summary['failures']}\t"
            f"{summary['success_rate']:.4f}\t"
            f"{summary['avg_runtime']:.6f}\t"
            f"{summary['max_runtime']:.6f}\t"
            f"{summary['avg_segments']:.3f}\t"
            f"{summary['max_segments']}\t"
            f"{summary['avg_probes']:.3f}\t"
            f"{summary['max_probes']}\t"
            f"{summary['avg_product_steps']:.1f}\t"
            f"{summary['max_product_steps']}\t"
            f"{summary['avg_max_d']:.1f}\t"
            f"{summary['max_d']}"
        )

        print()

        print_slowest(
            results
        )

    elapsed = (
        time.perf_counter()
        - overall_start
    )

    print(
        f"TOTAL_RUNTIME={elapsed:.6f}s"
    )
    print()

    print(
        f"FINISHED EXPERIMENT {EXPERIMENT_ID}"
    )


if __name__ == "__main__":
    main()
