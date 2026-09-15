import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 149

TARGET_N = 10**18

RATIO_MIN = 1.05
RATIO_MAX = 1.50

ALPHA = 1.50

BLOCK_SIZE = 256


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


def binomial_gcd_chunks(
    n,
    m,
    d,
    block_size,
):
    """
    Compute gcd(C(m,d), N) using chunks of the numerator interval.

    Because d < p, the denominator d! is coprime to N, so:

        gcd(C(m,d), N)
        =
        gcd(product(m-d+1..m), N)

    We reduce modulo N after each 256-element chunk.
    """

    if d <= 0:
        return 1, 0, 0

    start = m - d + 1
    end = m + 1

    accumulator = 1

    position = start
    product_steps = 0
    gcd_calls = 0

    while position < end:
        block_end = min(
            position + block_size,
            end,
        )

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

        product_steps += (
            block_end - position
        )

        gcd_calls += 1

        g = math.gcd(
            accumulator,
            n,
        )

        if 1 < g < n:
            return (
                g,
                product_steps,
                gcd_calls,
            )

        if g == n:
            return (
                n,
                product_steps,
                gcd_calls,
            )

        position = block_end

    return (
        math.gcd(
            accumulator,
            n,
        ),
        product_steps,
        gcd_calls,
    )


def binary_search_segment(
    n,
    p,
    q,
    lower,
    upper,
):
    probes = 0
    product_steps = 0
    gcd_calls = 0

    factor_probe = None
    factor_step = None

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
                "factor_probe": factor_probe,
                "factor_step": factor_step,
            }

        probes += 1

        g, used_steps, used_gcd = (
            binomial_gcd_chunks(
                n,
                m,
                d,
                BLOCK_SIZE,
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
            factor_probe = probes
            factor_step = used_steps

            return {
                "success": True,
                "candidate": p,
                "reason": "P",
                "probes": probes,
                "product_steps": product_steps,
                "gcd_calls": gcd_calls,
                "factor_probe": factor_probe,
                "factor_step": factor_step,
            }

        if signal == "Q":
            factor_probe = probes
            factor_step = used_steps

            return {
                "success": True,
                "candidate": q,
                "reason": "Q",
                "probes": probes,
                "product_steps": product_steps,
                "gcd_calls": gcd_calls,
                "factor_probe": factor_probe,
                "factor_step": factor_step,
            }

        if signal != "MISS":
            return {
                "success": False,
                "candidate": None,
                "reason": signal,
                "probes": probes,
                "product_steps": product_steps,
                "gcd_calls": gcd_calls,
                "factor_probe": factor_probe,
                "factor_step": factor_step,
            }

        upper = midpoint

    return {
        "success": lower == p,
        "candidate": lower,
        "reason": "binary" if lower == p else "wrong",
        "probes": probes,
        "product_steps": product_steps,
        "gcd_calls": gcd_calls,
        "factor_probe": factor_probe,
        "factor_step": factor_step,
    }


def factor_n(n, p, q):
    lower, upper = make_initial_interval(n)

    segments = subdivide_interval(
        lower,
        upper,
    )

    total_probes = 0
    total_product_steps = 0
    total_gcd_calls = 0

    start = time.perf_counter()

    for index, (segment_lower, segment_upper) in enumerate(
        segments
    ):
        result = binary_search_segment(
            n,
            p,
            q,
            segment_lower,
            segment_upper,
        )

        total_probes += result["probes"]
        total_product_steps += result["product_steps"]
        total_gcd_calls += result["gcd_calls"]

        if result["success"]:
            return {
                "success": True,
                "candidate": result["candidate"],
                "reason": result["reason"],
                "runtime": (
                    time.perf_counter()
                    - start
                ),
                "segments": index + 1,
                "total_segments": len(segments),
                "probes": total_probes,
                "product_steps": total_product_steps,
                "gcd_calls": total_gcd_calls,
                "factor_probe": result["factor_probe"],
                "factor_step": result["factor_step"],
            }

    return {
        "success": False,
        "candidate": None,
        "reason": "no-factor",
        "runtime": (
            time.perf_counter()
            - start
        ),
        "segments": len(segments),
        "total_segments": len(segments),
        "probes": total_probes,
        "product_steps": total_product_steps,
        "gcd_calls": total_gcd_calls,
        "factor_probe": None,
        "factor_step": None,
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    print("GENERATING CASE")

    p, q, n = generate_semiprime_near(
        TARGET_N
    )

    print(f"p={p}")
    print(f"q={q}")
    print(f"N={n}")
    print(f"ratio={q / p:.9f}")
    print()

    lower, upper = make_initial_interval(n)

    print(
        f"initial_interval=[{lower},{upper}]"
    )

    segments = subdivide_interval(
        lower,
        upper,
    )

    print(
        f"segments={len(segments)}"
    )
    print(
        f"block_size={BLOCK_SIZE}"
    )
    print()

    start = time.perf_counter()

    result = factor_n(
        n,
        p,
        q,
    )

    total_runtime = (
        time.perf_counter()
        - start
    )

    print("RESULT")
    print()

    print(
        f"success={result['success']}"
    )
    print(
        f"candidate={result['candidate']}"
    )
    print(
        f"correct="
        f"{result['candidate'] in (p, q)}"
    )
    print(
        f"reason={result['reason']}"
    )
    print(
        f"runtime_seconds="
        f"{result['runtime']:.6f}"
    )
    print(
        f"total_runtime_seconds="
        f"{total_runtime:.6f}"
    )
    print(
        f"segments_tested="
        f"{result['segments']}"
    )
    print(
        f"total_segments="
        f"{result['total_segments']}"
    )
    print(
        f"probes="
        f"{result['probes']}"
    )
    print(
        f"product_steps="
        f"{result['product_steps']}"
    )
    print(
        f"gcd_calls="
        f"{result['gcd_calls']}"
    )
    print(
        f"factor_probe="
        f"{result['factor_probe']}"
    )
    print(
        f"factor_step="
        f"{result['factor_step']}"
    )

    print()
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT_ID}"
    )


if __name__ == "__main__":
    main()
