import math
import random
import sys
import time
import tracemalloc
from array import array

from sympy import isprime


EXPERIMENT = 127

CASES = 1

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

TEST_D = (
    16_383,
    32_767,
    65_535,
    131_071,
    262_143,
    524_287,
    1_048_575,
)

TEST_TIME_LIMIT = 5.0


def generate_semiprime():
    while True:
        p = random.randrange(
            P_MIN,
            P_MAX + 1,
        )

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, Q_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        if p == q:
            continue

        return p, q


def reference_comb_mod(m, d, n):
    start = time.perf_counter()

    value = math.comb(m, d)
    value %= n

    elapsed = (
        time.perf_counter()
        - start
    )

    return value, elapsed


def batch_list(m, d, n):
    """
    Current implementation.

    Prefix values are Python integers stored in a list.
    """

    start = time.perf_counter()

    prefix = [1] * (d + 1)
    product = 1

    for j in range(1, d + 1):
        product = (
            product * j
        ) % n

        prefix[j] = product

        if (
            time.perf_counter()
            - start
            > TEST_TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "time": (
                    time.perf_counter()
                    - start
                ),
                "peak": None,
            }

    factor = math.gcd(
        product,
        n,
    )

    if factor != 1:
        return {
            "status": "factor",
            "factor": factor,
            "time": (
                time.perf_counter()
                - start
            ),
            "peak": None,
        }

    tracemalloc.start()

    inverse_product = pow(
        product,
        -1,
        n,
    )

    result = 1

    for j in range(d, 0, -1):
        inverse_j = (
            inverse_product
            * prefix[j - 1]
        ) % n

        numerator = (
            m - d + j
        )

        result *= (
            numerator % n
        )
        result %= n

        result *= inverse_j
        result %= n

        inverse_product = (
            inverse_product * j
        ) % n

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "status": "ok",
        "value": result,
        "time": elapsed,
        "peak": peak,
        "prefix_bytes_estimate": (
            (d + 1) * sys.getsizeof(1)
        ),
    }


def batch_array(m, d, n):
    """
    Same algorithm, but the prefix table is stored as
    unsigned 64-bit integers instead of Python int objects.
    """

    if n >= 2**64:
        return {
            "status": "unsupported",
            "time": 0.0,
            "peak": None,
        }

    start = time.perf_counter()

    prefix = array(
        "Q",
        [1],
    )

    product = 1

    for j in range(1, d + 1):
        product = (
            product * j
        ) % n

        prefix.append(product)

        if (
            time.perf_counter()
            - start
            > TEST_TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "time": (
                    time.perf_counter()
                    - start
                ),
                "peak": None,
            }

    factor = math.gcd(
        product,
        n,
    )

    if factor != 1:
        return {
            "status": "factor",
            "factor": factor,
            "time": (
                time.perf_counter()
                - start
            ),
            "peak": None,
        }

    tracemalloc.start()

    inverse_product = pow(
        product,
        -1,
        n,
    )

    result = 1

    for j in range(d, 0, -1):
        inverse_j = (
            inverse_product
            * prefix[j - 1]
        ) % n

        numerator = (
            m - d + j
        )

        result *= (
            numerator % n
        )
        result %= n

        result *= inverse_j
        result %= n

        inverse_product = (
            inverse_product * j
        ) % n

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "status": "ok",
        "value": result,
        "time": elapsed,
        "peak": peak,
        "prefix_bytes_estimate": len(prefix),
    }


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    random.seed(EXPERIMENT)

    p, q = generate_semiprime()
    n = p * q

    print(f"N={n}")
    print(f"p={p}")
    print(f"q={q}")
    print()

    print(
        "Each test uses m=k-2 with "
        "d approximately as requested."
    )
    print()

    for d in TEST_D:
        # Choose a nearby k so that d = 2z-1.
        k = 2 * d + 2
        m = k - 2

        print("-" * 60)
        print(
            f"d={d} "
            f"m={m}"
        )

        # -----------------------------------------------
        # EXACT REFERENCE
        # -----------------------------------------------

        if d <= 131_071:
            try:
                (
                    reference_value,
                    reference_time,
                ) = reference_comb_mod(
                    m,
                    d,
                    n,
                )

                print(
                    f"reference comb              "
                    f"time={reference_time:.6f}s"
                )
            except MemoryError:
                reference_value = None
                print(
                    "reference comb              "
                    "MEMORY ERROR"
                )
        else:
            reference_value = None

            print(
                "reference comb              "
                "SKIPPED"
            )

        # -----------------------------------------------
        # LIST
        # -----------------------------------------------

        list_result = batch_list(
            m,
            d,
            n,
        )

        if list_result["status"] == "ok":
            matches = (
                reference_value is None
                or list_result["value"]
                == reference_value
            )

            print(
                f"batch list                  "
                f"time={list_result['time']:.6f}s "
                f"status=ok "
                f"matches={matches}"
            )

            print(
                f"list prefix estimate        "
                f"~{list_result['prefix_bytes_estimate'] / (1024 * 1024):.2f} MiB"
            )

        else:
            print(
                f"batch list                  "
                f"time={list_result['time']:.6f}s "
                f"status={list_result['status']}"
            )

        # -----------------------------------------------
        # ARRAY
        # -----------------------------------------------

        array_result = batch_array(
            m,
            d,
            n,
        )

        if array_result["status"] == "ok":
            matches = (
                reference_value is None
                or array_result["value"]
                == reference_value
            )

            print(
                f"batch array                 "
                f"time={array_result['time']:.6f}s "
                f"status=ok "
                f"matches={matches}"
            )

            print(
                f"array prefix size           "
                f"~{(d + 1) * 8 / (1024 * 1024):.2f} MiB"
            )

        else:
            print(
                f"batch array                 "
                f"time={array_result['time']:.6f}s "
                f"status={array_result['status']}"
            )

        # -----------------------------------------------
        # DECISION
        # -----------------------------------------------

        if (
            list_result["status"] == "ok"
            and array_result["status"] == "ok"
        ):
            speedup = (
                list_result["time"]
                / max(
                    array_result["time"],
                    1e-12,
                )
            )

            print(
                f"array/list speedup          "
                f"{speedup:.2f}x"
            )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
