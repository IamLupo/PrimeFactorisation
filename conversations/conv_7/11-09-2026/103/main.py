import math
import random
import time

from sympy import isprime


EXPERIMENT = 128

CASES = 8

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

MAX_PROBES = 2

METHOD_TIME_LIMIT = 8.0


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
        reversed_value = bit_reverse(
            index,
            bits,
        )

        z = (
            reversed_value * max_z
            // (1 << bits)
        ) + 1

        z = max(
            1,
            min(z, max_z),
        )

        if z not in used:
            used.add(z)
            result.append(z)

        index += 1

    return result


def exact_comb(m, d, n):
    start = time.perf_counter()

    value = math.comb(
        m,
        d,
    ) % n

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "status": "ok",
        "value": value,
        "factor": None,
        "time": elapsed,
    }


def batch_modular(m, d, n):
    """
    Baseline batch modular evaluator.

    It builds the entire prefix table before checking
    whether the denominator contains a factor.
    """

    start = time.perf_counter()

    prefix = [1] * (d + 1)

    product = 1

    for j in range(
        1,
        d + 1,
    ):
        product = (
            product * j
        ) % n

        prefix[j] = product

        if (
            time.perf_counter()
            - start
            > METHOD_TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "time": (
                    time.perf_counter()
                    - start
                ),
                "iterations": j,
            }

    factor = math.gcd(
        product,
        n,
    )

    if 1 < factor < n:
        return {
            "status": "factor",
            "factor": factor,
            "time": (
                time.perf_counter()
                - start
            ),
            "iterations": d,
        }

    inverse_product = pow(
        product,
        -1,
        n,
    )

    result = 1

    for j in range(
        d,
        0,
        -1,
    ):
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

    return {
        "status": "ok",
        "value": result,
        "factor": None,
        "time": (
            time.perf_counter()
            - start
        ),
        "iterations": d,
    }


def batch_early_factor(m, d, n):
    """
    Batch modular evaluator with early denominator
    factor detection.

    During construction of d! mod n, test:

        gcd(j, n)

    before continuing.

    If a nontrivial gcd appears, return it immediately.
    """

    start = time.perf_counter()

    prefix = [1] * (d + 1)

    product = 1

    for j in range(
        1,
        d + 1,
    ):
        g = math.gcd(
            j,
            n,
        )

        if 1 < g < n:
            return {
                "status": "factor",
                "factor": g,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "iterations": j,
            }

        product = (
            product * j
        ) % n

        prefix[j] = product

        if (
            time.perf_counter()
            - start
            > METHOD_TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "time": (
                    time.perf_counter()
                    - start
                ),
                "iterations": j,
            }

    inverse_product = pow(
        product,
        -1,
        n,
    )

    result = 1

    for j in range(
        d,
        0,
        -1,
    ):
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

    return {
        "status": "ok",
        "value": result,
        "factor": None,
        "time": (
            time.perf_counter()
            - start
        ),
        "iterations": d,
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
    print(f"TIME LIMIT: {METHOD_TIME_LIMIT}s")
    print()

    random.seed(EXPERIMENT)

    methods = (
        "exact",
        "batch",
        "early",
    )

    total_time = {
        method: 0.0
        for method in methods
    }

    evaluations = {
        method: 0
        for method in methods
    }

    factor_hits = {
        method: 0
        for method in methods
    }

    timeouts = {
        method: 0
        for method in methods
    }

    validation_matches = {
        "batch": 0,
        "early": 0,
    }

    validation_mismatches = {
        "batch": 0,
        "early": 0,
    }

    early_iterations = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()

        n = p * q

        k_values = powers_between(
            p,
            q,
        )

        if not k_values:
            continue

        # Use the first useful power-of-two k.
        k = k_values[0]

        probes = van_der_corput_sequence(
            k // 2,
            MAX_PROBES,
        )

        print(
            f"CASE {case_index}: "
            f"N={n} "
            f"p={p} "
            f"q={q} "
            f"k={k}"
        )

        for probe_index, z in enumerate(
            probes,
            start=1,
        ):
            d = 2 * z - 1
            m = k - 2

            print(
                f"  PROBE {probe_index}: "
                f"z={z} "
                f"d={d} "
                f"d/p={d / p:.4f}"
            )

            # ---------------------------------------------
            # EXACT
            # ---------------------------------------------

            exact = exact_comb(
                m,
                d,
                n,
            )

            total_time["exact"] += (
                exact["time"]
            )

            evaluations["exact"] += 1

            if exact["status"] == "timeout":
                timeouts["exact"] += 1

            # ---------------------------------------------
            # BATCH
            # ---------------------------------------------

            batch = batch_modular(
                m,
                d,
                n,
            )

            total_time["batch"] += (
                batch["time"]
            )

            evaluations["batch"] += 1

            if batch["status"] == "timeout":
                timeouts["batch"] += 1

            elif batch["status"] == "factor":
                factor_hits["batch"] += 1

            # ---------------------------------------------
            # EARLY
            # ---------------------------------------------

            early = batch_early_factor(
                m,
                d,
                n,
            )

            total_time["early"] += (
                early["time"]
            )

            evaluations["early"] += 1

            if early["status"] == "timeout":
                timeouts["early"] += 1

            elif early["status"] == "factor":
                factor_hits["early"] += 1

                early_iterations.append(
                    early["iterations"]
                )

            # ---------------------------------------------
            # VALIDATION
            # ---------------------------------------------

            if (
                exact["status"] == "ok"
                and batch["status"] == "ok"
            ):
                if (
                    batch["value"]
                    == exact["value"]
                ):
                    validation_matches[
                        "batch"
                    ] += 1
                else:
                    validation_mismatches[
                        "batch"
                    ] += 1

            if (
                exact["status"] == "ok"
                and early["status"] == "ok"
            ):
                if (
                    early["value"]
                    == exact["value"]
                ):
                    validation_matches[
                        "early"
                    ] += 1
                else:
                    validation_mismatches[
                        "early"
                    ] += 1

    print()
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    for method in methods:
        print()
        print(
            f"[{method.upper()}]"
        )

        print(
            f"evaluations                 = "
            f"{evaluations[method]}"
        )

        print(
            f"total time                  = "
            f"{total_time[method]:.6f}s"
        )

        print(
            f"average time                = "
            f"{total_time[method] / max(1, evaluations[method]):.6f}s"
        )

        print(
            f"timeouts                    = "
            f"{timeouts[method]}"
        )

        print(
            f"factor shortcuts            = "
            f"{factor_hits[method]}"
        )

    print()
    print(
        f"early vs batch speedup       = "
        f"{total_time['batch'] / max(total_time['early'], 1e-12):.2f}x"
    )

    print(
        f"early vs exact speedup       = "
        f"{total_time['exact'] / max(total_time['early'], 1e-12):.2f}x"
    )

    print()
    print("-" * 60)
    print("VALIDATION")
    print("-" * 60)

    print(
        f"batch matches                = "
        f"{validation_matches['batch']}"
    )

    print(
        f"batch mismatches             = "
        f"{validation_mismatches['batch']}"
    )

    print(
        f"early matches                = "
        f"{validation_matches['early']}"
    )

    print(
        f"early mismatches             = "
        f"{validation_mismatches['early']}"
    )

    if early_iterations:
        print()
        print(
            f"average early factor j       = "
            f"{sum(early_iterations) / len(early_iterations):.2f}"
        )

        print(
            f"maximum early factor j       = "
            f"{max(early_iterations)}"
        )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
