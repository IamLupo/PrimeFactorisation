import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 145

CASES_PER_RANGE = 100

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

# Use intervals [L,U] with U <= alpha*L.
#
# alpha < 2 is required.
ALPHA = 1.50


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


def make_initial_interval(n):
    lower = math.ceil(
        math.sqrt(n / ASSUMED_RATIO)
    )

    upper = math.isqrt(n)

    return lower, upper


def subdivide_interval(lower, upper):
    """
    Construct adjacent multiplicative intervals.

    Each interval satisfies:

        U <= ALPHA * L

    so for ALPHA < 2:

        U < 2L.

    The final interval is clipped to upper.
    """

    intervals = []

    current = lower

    while current <= upper:
        candidate_upper = math.floor(
            ALPHA * current
        )

        if candidate_upper <= current:
            candidate_upper = current + 1

        candidate_upper = min(
            candidate_upper,
            upper,
        )

        intervals.append(
            (
                current,
                candidate_upper,
            )
        )

        if candidate_upper == upper:
            break

        current = candidate_upper + 1

    return intervals


def make_probe(lower, upper):
    """
    Same valid binary construction as Experiment 138.

        m = U
        t = floor((L+U)/2)
        d = U-t
        k = m+2
    """

    if lower >= upper:
        return None

    midpoint = (lower + upper) // 2

    m = upper
    d = upper - midpoint
    k = m + 2

    return {
        "k": k,
        "m": m,
        "d": d,
        "midpoint": midpoint,
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


def binary_search_segment(n, p, q, lower, upper):
    """
    Search one safe segment.

    The segment is guaranteed to satisfy U < 2L,
    therefore floor(U/p)=1 for every p in [L,U].
    """

    history = []

    while lower < upper:
        probe = make_probe(
            lower,
            upper,
        )

        if probe is None:
            break

        m = probe["m"]
        d = probe["d"]
        midpoint = probe["midpoint"]
        k = probe["k"]

        # Mathematical invariant.
        safe = (
            d < lower
            and m < 2 * lower
            and m < q
        )

        if not safe:
            return {
                "success": False,
                "candidate": None,
                "reason": "unsafe-segment",
                "history": history,
            }

        g = evaluate_probe(
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
                "lower": lower,
                "upper": upper,
                "midpoint": midpoint,
                "m": m,
                "k": k,
                "d": d,
                "gcd": g,
                "signal": signal,
            }
        )

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

        if signal == "N" or signal == "OTHER":
            return {
                "success": False,
                "candidate": None,
                "reason": signal,
                "history": history,
            }

        # MISS:
        #
        # p <= midpoint.
        upper = midpoint

    if lower == p:
        return {
            "success": True,
            "candidate": p,
            "reason": "binary-p",
            "history": history,
        }

    return {
        "success": False,
        "candidate": lower,
        "reason": "wrong-segment",
        "history": history,
    }


def search_all_segments(n, p, q, segments):
    all_history = []

    for index, (lower, upper) in enumerate(
        segments
    ):
        # We need p to actually lie in this segment.
        # This is only diagnostic information; the algorithm below
        # does NOT use p to choose the segment.
        result = binary_search_segment(
            n,
            p,
            q,
            lower,
            upper,
        )

        all_history.append(
            {
                "segment_index": index,
                "lower": lower,
                "upper": upper,
                "result": result,
                "contains_p": lower <= p <= upper,
            }
        )

        if result["success"]:
            return {
                "success": True,
                "candidate": result["candidate"],
                "reason": result["reason"],
                "history": all_history,
            }

    return {
        "success": False,
        "candidate": None,
        "reason": "no-segment-hit",
        "history": all_history,
    }


def summarize(result):
    segments = result["history"]

    tested_segments = len(segments)

    binary_probes = sum(
        len(
            segment["result"]["history"]
        )
        for segment in segments
    )

    p_segment_index = None

    for segment in segments:
        if segment["contains_p"]:
            p_segment_index = (
                segment["segment_index"]
            )
            break

    successful_segment_index = None

    for segment in segments:
        if segment["result"]["success"]:
            successful_segment_index = (
                segment["segment_index"]
            )
            break

    return {
        "success": result["success"],
        "segments": tested_segments,
        "probes": binary_probes,
        "p_segment": p_segment_index,
        "successful_segment": successful_segment_index,
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    start = time.perf_counter()

    print(
        "actual_ratio\t"
        "cases\t"
        "successes\t"
        "failures\t"
        "success_rate\t"
        "avg_segments\t"
        "max_segments\t"
        "avg_binary_probes\t"
        "max_binary_probes"
    )

    for actual_min, actual_max in ACTUAL_RATIO_RANGES:
        successes = 0
        total_segments = 0
        total_probes = 0
        max_segments = 0
        max_probes = 0

        for _ in range(CASES_PER_RANGE):
            p, q, n = generate_semiprime(
                actual_min,
                actual_max,
            )

            lower, upper = make_initial_interval(
                n
            )

            segments = subdivide_interval(
                lower,
                upper,
            )

            result = search_all_segments(
                n,
                p,
                q,
                segments,
            )

            summary = summarize(result)

            if summary["success"]:
                successes += 1

            total_segments += summary["segments"]
            total_probes += summary["probes"]

            max_segments = max(
                max_segments,
                summary["segments"],
            )

            max_probes = max(
                max_probes,
                summary["probes"],
            )

        print(
            f"{actual_min:.2f}-{actual_max:.2f}\t"
            f"{CASES_PER_RANGE}\t"
            f"{successes}\t"
            f"{CASES_PER_RANGE-successes}\t"
            f"{successes / CASES_PER_RANGE:.4f}\t"
            f"{total_segments / CASES_PER_RANGE:.3f}\t"
            f"{max_segments}\t"
            f"{total_probes / CASES_PER_RANGE:.3f}\t"
            f"{max_probes}"
        )

    print()
    print("DETAILED EXAMPLE")
    print()

    p, q, n = generate_semiprime(
        8.0,
        10.0,
    )

    lower, upper = make_initial_interval(n)

    segments = subdivide_interval(
        lower,
        upper,
    )

    result = search_all_segments(
        n,
        p,
        q,
        segments,
    )

    print(f"p={p}")
    print(f"q={q}")
    print(f"N={n}")
    print(
        f"initial_interval=[{lower},{upper}]"
    )
    print(
        f"segments={len(segments)}"
    )
    print(
        f"candidate={result['candidate']}"
    )
    print(
        f"success={result['success']}"
    )
    print(
        f"reason={result['reason']}"
    )
    print()

    for segment in result["history"]:
        index = segment["segment_index"]
        sl = segment["lower"]
        su = segment["upper"]

        print(
            f"SEGMENT {index}: "
            f"[{sl},{su}] "
            f"contains_p="
            f"{segment['contains_p']}"
        )

        for step_index, step in enumerate(
            segment["result"]["history"],
            1,
        ):
            print(
                f"  {step_index:02d}: "
                f"mid={step['midpoint']} "
                f"d={step['d']} "
                f"gcd={step['gcd']} "
                f"{step['signal']}"
            )

        print()

    elapsed = time.perf_counter() - start

    print(
        f"runtime_seconds={elapsed:.6f}"
    )
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_ID}")


if __name__ == "__main__":
    main()
