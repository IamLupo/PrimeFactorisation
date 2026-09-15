import math
import random
import time

from sympy import isprime


EXPERIMENT = 131

CASES = 10

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

MAX_PROBES = 2

TIME_LIMIT = 5.0


def generate_semiprime():
    while True:
        p = random.randrange(P_MIN, P_MAX + 1)

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, Q_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        return p, q


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def bit_reverse(value, bits):
    result = 0

    for _ in range(bits):
        result <<= 1
        result |= value & 1
        value >>= 1

    return result


def van_der_corput_sequence(max_z, count):
    count = min(count, max_z)

    if count <= 0:
        return []

    bits = max(
        1,
        (count - 1).bit_length(),
    )

    result = []
    used = set()

    index = 0

    while len(result) < count:
        reversed_value = bit_reverse(index, bits)

        z = (
            reversed_value * max_z
            // (1 << bits)
        ) + 1

        z = max(1, min(z, max_z))

        if z not in used:
            used.add(z)
            result.append(z)

        index += 1

    return result


def numerator_linear(n, left, right):
    """
    Reference O(d) numerator scan.
    """

    start = time.perf_counter()

    product = 1

    for x in range(left, right + 1):
        product = (
            product * x
        ) % n

        g = math.gcd(product, n)

        if 1 < g < n:
            return {
                "status": "factor",
                "factor": g,
                "iterations": x - left + 1,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

        if (
            time.perf_counter()
            - start
            > TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "factor": None,
                "iterations": x - left + 1,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

    return {
        "status": "none",
        "factor": None,
        "iterations": right - left + 1,
        "time": (
            time.perf_counter()
            - start
        ),
    }


def build_prefix(n, limit):
    """
    Build prefix factorials modulo n:

        P[x] = x! mod n

    while checking whether a prefix exposes a factor.
    """

    start = time.perf_counter()

    prefix = [1] * (limit + 1)

    product = 1

    for x in range(1, limit + 1):
        product = (
            product * x
        ) % n

        prefix[x] = product

        g = math.gcd(product, n)

        if 1 < g < n:
            return {
                "status": "factor",
                "factor": g,
                "prefix": prefix,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "limit": x,
            }

        if (
            time.perf_counter()
            - start
            > TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "factor": None,
                "prefix": prefix,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "limit": x,
            }

    return {
        "status": "ok",
        "factor": None,
        "prefix": prefix,
        "time": (
            time.perf_counter()
            - start
        ),
        "limit": limit,
    }


def interval_product_from_prefix(
    n,
    prefix,
    left,
    right,
):
    """
    Attempt:

        right! / (left-1)!

    modulo n.

    If the prefix denominator is not invertible,
    expose that gcd instead.
    """

    if left > right:
        return {
            "status": "empty"
        }

    denominator = prefix[left - 1]

    g = math.gcd(
        denominator,
        n,
    )

    if 1 < g < n:
        return {
            "status": "factor",
            "factor": g,
        }

    if g != 1:
        return {
            "status": "bad"
        }

    inverse = pow(
        denominator,
        -1,
        n,
    )

    value = (
        prefix[right]
        * inverse
    ) % n

    return {
        "status": "ok",
        "value": value,
    }


def binary_interval_search(
    n,
    prefix,
    left,
    right,
):
    """
    Binary search over the numerator interval.

    At each node:

        product(left..mid)
        product(mid+1..right)

    and inspect gcds.

    This is experimental; it is not assumed to be
    asymptotically superior to the linear scan yet.
    """

    start = time.perf_counter()

    nodes = 0

    stack = [
        (left, right)
    ]

    while stack:
        a, b = stack.pop()

        if a > b:
            continue

        nodes += 1

        if a == b:
            g = math.gcd(a, n)

            if 1 < g < n:
                return {
                    "status": "factor",
                    "factor": g,
                    "nodes": nodes,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                }

            continue

        mid = (a + b) // 2

        left_result = (
            interval_product_from_prefix(
                n,
                prefix,
                a,
                mid,
            )
        )

        if left_result["status"] == "factor":
            return {
                "status": "factor",
                "factor": left_result["factor"],
                "nodes": nodes,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

        if left_result["status"] == "ok":
            g = math.gcd(
                left_result["value"],
                n,
            )

            if 1 < g < n:
                stack.append(
                    (a, mid)
                )
            else:
                stack.append(
                    (a, mid)
                )

        right_result = (
            interval_product_from_prefix(
                n,
                prefix,
                mid + 1,
                b,
            )
        )

        if right_result["status"] == "factor":
            return {
                "status": "factor",
                "factor": right_result["factor"],
                "nodes": nodes,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

        if right_result["status"] == "ok":
            g = math.gcd(
                right_result["value"],
                n,
            )

            if 1 < g < n:
                stack.append(
                    (mid + 1, b)
                )
            else:
                stack.append(
                    (mid + 1, b)
                )

    return {
        "status": "none",
        "factor": None,
        "nodes": nodes,
        "time": (
            time.perf_counter()
            - start
        ),
    }


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX}")
    print(f"MAX PROBES: {MAX_PROBES}")
    print(f"TIME LIMIT: {TIME_LIMIT}s")
    print()

    random.seed(EXPERIMENT)

    linear_total = 0.0
    binary_total = 0.0
    prefix_total = 0.0

    linear_hits = 0
    binary_hits = 0

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        k = powers_between(
            p,
            q,
        )[0]

        probes = van_der_corput_sequence(
            k // 2,
            MAX_PROBES,
        )

        # Use the first probe that actually produces
        # a nontrivial Lucas numerator interval.
        z = probes[-1]

        d = 2 * z - 1

        if d >= k - 1:
            continue

        left = k - 1 - d
        right = k - 2

        print()
        print(
            f"CASE {case_index}: "
            f"p={p} "
            f"q={q} "
            f"k={k} "
            f"z={z} "
            f"d={d}"
        )

        # -------------------------------------------------
        # LINEAR
        # -------------------------------------------------

        linear = numerator_linear(
            n,
            left,
            right,
        )

        linear_total += linear["time"]

        if linear["status"] == "factor":
            linear_hits += 1

        print(
            f"  LINEAR: "
            f"status={linear['status']} "
            f"factor={linear['factor']} "
            f"iterations={linear['iterations']} "
            f"time={linear['time']:.6f}s"
        )

        if linear["status"] == "timeout":
            print(
                "  Skipping binary benchmark "
                "because linear scan timed out."
            )
            continue

        # -------------------------------------------------
        # PREFIX
        # -------------------------------------------------

        prefix = build_prefix(
            n,
            right,
        )

        prefix_total += prefix["time"]

        print(
            f"  PREFIX: "
            f"status={prefix['status']} "
            f"limit={prefix['limit']} "
            f"time={prefix['time']:.6f}s"
        )

        if prefix["status"] != "ok":
            print(
                "  Prefix construction already "
                "exposed a factor."
            )
            continue

        # -------------------------------------------------
        # BINARY
        # -------------------------------------------------

        binary = binary_interval_search(
            n,
            prefix["prefix"],
            left,
            right,
        )

        binary_total += binary["time"]

        if binary["status"] == "factor":
            binary_hits += 1

        print(
            f"  BINARY: "
            f"status={binary['status']} "
            f"factor={binary['factor']} "
            f"nodes={binary['nodes']} "
            f"time={binary['time']:.6f}s"
        )

    print()
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"linear hits                = "
        f"{linear_hits}"
    )

    print(
        f"binary hits                = "
        f"{binary_hits}"
    )

    print(
        f"linear total time          = "
        f"{linear_total:.6f}s"
    )

    print(
        f"prefix total time          = "
        f"{prefix_total:.6f}s"
    )

    print(
        f"binary total time          = "
        f"{binary_total:.6f}s"
    )

    if binary_total > 0:
        print(
            f"linear/binary speed ratio  = "
            f"{linear_total / binary_total:.2f}x"
        )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
