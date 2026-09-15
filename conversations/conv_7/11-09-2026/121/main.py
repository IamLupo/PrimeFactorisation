import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 147

CASES_PER_TIER = 5

N_TIERS = [
    10**14,
    10**16,
    10**18,
    10**20,
]

RATIO_MIN = 1.05
RATIO_MAX = 1.50

ALPHA = 1.50

BLOCK_SIZES = [
    256,
    1024,
    4096,
]


def generate_semiprime_near(target_n):
    target = math.sqrt(target_n)

    while True:
        ratio = random.uniform(
            RATIO_MIN,
            RATIO_MAX,
        )

        p_estimate = int(
            target / math.sqrt(ratio)
        )

        spread = max(
            100,
            p_estimate // 100,
        )

        p = int(
            nextprime(
                random.randint(
                    max(2, p_estimate - spread),
                    p_estimate + spread,
                )
            )
        )

        q_estimate = target_n // p

        q = int(
            nextprime(q_estimate)
        )

        if q <= p:
            continue

        n = p * q
        ratio_actual = q / p

        if not (
            0.75 * target_n
            <= n
            <= 1.50 * target_n
        ):
            continue

        if not (
            RATIO_MIN
            <= ratio_actual
            <= RATIO_MAX
        ):
            continue

        return p, q, n


def make_initial_interval(n):
    lower = math.ceil(
        math.sqrt(n / RATIO_MAX)
    )

    upper = math.isqrt(n)

    return lower, upper


def subdivide_interval(lower, upper):
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

    return m, d, midpoint


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


def block_binomial_gcd(
    n,
    m,
    d,
    p,
    q,
    block_size,
):
    """
    Evaluate gcd(C(m,d), N) through the numerator interval

        (m-d+1) ... m

    but only calculate gcd once per block.

    The accumulator is always reduced modulo N.

    Since d < p and q is outside the interval, the block
    containing the unique multiple of p produces gcd=p.
    """

    if d <= 0:
        return 1, 0, 0

    start = m - d + 1
    end = m + 1

    accumulator = 1
    steps = 0
    gcd_calls = 0

    position = start

    while position < end:
        block_end = min(
            position + block_size,
            end,
        )

        # math.prod is implemented at the Python/C boundary and
        # avoids explicitly writing a Python multiplication loop.
        block_product = math.prod(
            range(
                position,
                block_end,
            )
        )

        accumulator = (
            accumulator
            * (block_product % n)
        ) % n

        steps += (
            block_end - position
        )

        gcd_calls += 1

        g = math.gcd(
            accumulator,
            n,
        )

        if 1 < g < n:
            return g, steps, gcd_calls

        if g == n:
            return n, steps, gcd_calls

        position = block_end

    return (
        math.gcd(accumulator, n),
        steps,
        gcd_calls,
    )


def binary_segment(
    n,
    p,
    q,
    lower,
    upper,
    block_size,
):
    probes = 0
    product_steps = 0
    gcd_calls = 0
    max_d = 0

    while lower < upper:
        m, d, midpoint = make_probe(
            lower,
            upper,
        )

        if not (
            d < lower
            and m < 2 * lower
            and m < q
        ):
            return {
                "success": False,
                "candidate": None,
                "reason": "unsafe",
                "probes": probes,
                "product_steps": product_steps,
                "gcd_calls": gcd_calls,
                "max_d": max_d,
            }

        probes += 1
        max_d = max(
            max_d,
            d,
        )

        g, used_steps, used_gcd = (
            block_binomial_gcd(
                n,
                m,
                d,
                p,
                q,
                block_size,
            )
        )

        product_steps += used_steps
        gcd_calls += used_gcd

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
                "product_steps": product_steps,
                "gcd_calls": gcd_calls,
                "max_d": max_d,
            }

        if signal == "Q":
            return {
                "success": True,
                "candidate": q,
                "reason": "Q",
                "probes": probes,
                "product_steps": product_steps,
                "gcd_calls": gcd_calls,
                "max_d": max_d,
            }

        if signal != "MISS":
            return {
                "success": False,
                "candidate": None,
                "reason": signal,
                "probes": probes,
                "product_steps": product_steps,
                "gcd_calls": gcd_calls,
                "max_d": max_d,
            }

        upper = midpoint

    return {
        "success": lower == p,
        "candidate": lower,
        "reason": (
            "binary"
            if lower == p
            else "wrong"
        ),
        "probes": probes,
        "product_steps": product_steps,
        "gcd_calls": gcd_calls,
        "max_d": max_d,
    }


def factor_n(
    n,
    p,
    q,
    block_size,
):
    initial_lower, initial_upper = (
        make_initial_interval(n)
    )

    segments = subdivide_interval(
        initial_lower,
        initial_upper,
    )

    total_probes = 0
    total_steps = 0
    total_gcd_calls = 0
    max_d = 0
    segments_tested = 0

    start = time.perf_counter()

    for lower, upper in segments:
        segments_tested += 1

        result = binary_segment(
            n,
            p,
            q,
            lower,
            upper,
            block_size,
        )

        total_probes += result["probes"]
        total_steps += result["product_steps"]
        total_gcd_calls += result["gcd_calls"]

        max_d = max(
            max_d,
            result["max_d"],
        )

        if result["success"]:
            return {
                "success": True,
                "runtime": (
                    time.perf_counter()
                    - start
                ),
                "segments": segments_tested,
                "probes": total_probes,
                "product_steps": total_steps,
                "gcd_calls": total_gcd_calls,
                "max_d": max_d,
            }

    return {
        "success": False,
        "runtime": (
            time.perf_counter()
            - start
        ),
        "segments": segments_tested,
        "probes": total_probes,
        "product_steps": total_steps,
        "gcd_calls": total_gcd_calls,
        "max_d": max_d,
    }


def main():
    print(
        f"START EXPERIMENT {EXPERIMENT_ID}"
    )
    print()

    random.seed(EXPERIMENT_ID)

    for target_n in N_TIERS:

        print(
            f"TIER N~{target_n}"
        )
        print()

        for block_size in BLOCK_SIZES:

            successes = 0
            runtimes = []
            probes = []
            steps = []
            gcd_calls = []
            max_ds = []

            print(
                f"BLOCK_SIZE={block_size}"
            )

            for case in range(
                CASES_PER_TIER
            ):
                p, q, n = (
                    generate_semiprime_near(
                        target_n
                    )
                )

                result = factor_n(
                    n,
                    p,
                    q,
                    block_size,
                )

                successes += int(
                    result["success"]
                )

                runtimes.append(
                    result["runtime"]
                )

                probes.append(
                    result["probes"]
                )

                steps.append(
                    result["product_steps"]
                )

                gcd_calls.append(
                    result["gcd_calls"]
                )

                max_ds.append(
                    result["max_d"]
                )

            print(
                f"cases={CASES_PER_TIER}"
            )
            print(
                f"successes={successes}"
            )
            print(
                f"failures="
                f"{CASES_PER_TIER-successes}"
            )
            print(
                f"success_rate="
                f"{successes / CASES_PER_TIER:.4f}"
            )
            print(
                f"avg_runtime="
                f"{sum(runtimes)/len(runtimes):.6f}"
            )
            print(
                f"max_runtime="
                f"{max(runtimes):.6f}"
            )
            print(
                f"avg_probes="
                f"{sum(probes)/len(probes):.3f}"
            )
            print(
                f"avg_product_steps="
                f"{sum(steps)/len(steps):.1f}"
            )
            print(
                f"avg_gcd_calls="
                f"{sum(gcd_calls)/len(gcd_calls):.1f}"
            )
            print(
                f"max_gcd_calls="
                f"{max(gcd_calls)}"
            )
            print(
                f"avg_max_d="
                f"{sum(max_ds)/len(max_ds):.1f}"
            )
            print(
                f"max_d="
                f"{max(max_ds)}"
            )
            print()

    print(
        f"FINISHED EXPERIMENT {EXPERIMENT_ID}"
    )


if __name__ == "__main__":
    main()
