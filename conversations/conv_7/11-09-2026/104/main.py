import math
import random
import time

from sympy import isprime


EXPERIMENT = 129

CASES = 5

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

MAX_PROBES = 2

TIME_LIMIT = 5.0


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


def exact_binomial(n, m, d):
    start = time.perf_counter()

    value = math.comb(
        m,
        d,
    ) % n

    elapsed = (
        time.perf_counter()
        - start
    )

    return value, elapsed


def numerator_product(n, m, d):
    """
    Computes only the numerator:

        (m-d+1) * ... * m   mod N

    No denominator and no modular inverses.

    This is valid for the gcd test when d < p,
    because d! is then coprime to N.
    """

    start = time.perf_counter()

    result = 1

    first = m - d + 1

    for x in range(
        first,
        m + 1,
    ):
        result *= x % n
        result %= n

        # Keep the benchmark bounded.
        if (
            time.perf_counter()
            - start
            > TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "value": None,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "iterations": (
                    x - first + 1
                ),
            }

    return {
        "status": "ok",
        "value": result,
        "time": (
            time.perf_counter()
            - start
        ),
        "iterations": d,
    }


def numerator_product_with_gcd(n, m, d):
    """
    Same numerator product, but periodically check
    gcd(result, N).

    This is experimental: if the accumulated numerator
    contains a factor of N, we can return it immediately.
    """

    start = time.perf_counter()

    result = 1

    first = m - d + 1

    for index, x in enumerate(
        range(first, m + 1),
        start=1,
    ):
        result *= x % n
        result %= n

        g = math.gcd(
            result,
            n,
        )

        if 1 < g < n:
            return {
                "status": "factor",
                "factor": g,
                "value": result,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "iterations": index,
            }

        if (
            time.perf_counter()
            - start
            > TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "value": None,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "iterations": index,
            }

    return {
        "status": "ok",
        "value": result,
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
    print(f"TIME LIMIT: {TIME_LIMIT}s")
    print()

    random.seed(EXPERIMENT)

    exact_total = 0.0
    numerator_total = 0.0
    numerator_gcd_total = 0.0

    exact_count = 0
    numerator_count = 0
    numerator_gcd_count = 0

    matches = 0
    mismatches = 0

    factor_shortcuts = 0
    timeouts = 0

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
                f"d/p={d / p:.6f}"
            )

            # -------------------------------------------------
            # EXACT
            # -------------------------------------------------

            exact_value, exact_time = (
                exact_binomial(
                    n,
                    m,
                    d,
                )
            )

            exact_total += exact_time
            exact_count += 1

            exact_gcd = math.gcd(
                exact_value,
                n,
            )

            print(
                f"    exact: "
                f"time={exact_time:.6f}s "
                f"gcd={exact_gcd}"
            )

            # -------------------------------------------------
            # NUMERATOR ONLY
            # -------------------------------------------------

            numerator = numerator_product(
                n,
                m,
                d,
            )

            numerator_total += (
                numerator["time"]
            )

            numerator_count += 1

            if numerator["status"] == "timeout":
                timeouts += 1

                print(
                    f"    numerator: "
                    f"TIMEOUT "
                    f"iterations={numerator['iterations']}"
                )
            else:
                numerator_gcd = math.gcd(
                    numerator["value"],
                    n,
                )

                expected = (
                    exact_gcd
                )

                if (
                    numerator_gcd
                    == expected
                ):
                    matches += 1
                else:
                    mismatches += 1

                print(
                    f"    numerator: "
                    f"time={numerator['time']:.6f}s "
                    f"gcd={numerator_gcd} "
                    f"matches={numerator_gcd == expected}"
                )

            # -------------------------------------------------
            # NUMERATOR + GCD
            # -------------------------------------------------

            early = numerator_product_with_gcd(
                n,
                m,
                d,
            )

            numerator_gcd_total += (
                early["time"]
            )

            numerator_gcd_count += 1

            if early["status"] == "factor":
                factor_shortcuts += 1

                print(
                    f"    num+gcd: "
                    f"time={early['time']:.6f}s "
                    f"factor={early['factor']} "
                    f"iterations={early['iterations']}"
                )

            elif early["status"] == "timeout":
                print(
                    f"    num+gcd: "
                    f"TIMEOUT "
                    f"iterations={early['iterations']}"
                )

            else:
                early_gcd = math.gcd(
                    early["value"],
                    n,
                )

                print(
                    f"    num+gcd: "
                    f"time={early['time']:.6f}s "
                    f"gcd={early_gcd}"
                )

    print()
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"exact evaluations            = "
        f"{exact_count}"
    )

    print(
        f"numerator evaluations        = "
        f"{numerator_count}"
    )

    print(
        f"num+gcd evaluations          = "
        f"{numerator_gcd_count}"
    )

    print()

    print(
        f"exact total time             = "
        f"{exact_total:.6f}s"
    )

    print(
        f"numerator total time         = "
        f"{numerator_total:.6f}s"
    )

    print(
        f"num+gcd total time           = "
        f"{numerator_gcd_total:.6f}s"
    )

    print()

    print(
        f"exact average                = "
        f"{exact_total / max(1, exact_count):.6f}s"
    )

    print(
        f"numerator average            = "
        f"{numerator_total / max(1, numerator_count):.6f}s"
    )

    print(
        f"num+gcd average              = "
        f"{numerator_gcd_total / max(1, numerator_gcd_count):.6f}s"
    )

    print()

    print(
        f"numerator speedup            = "
        f"{exact_total / max(numerator_total, 1e-12):.2f}x"
    )

    print(
        f"num+gcd speedup              = "
        f"{exact_total / max(numerator_gcd_total, 1e-12):.2f}x"
    )

    print()

    print(
        f"matches                      = "
        f"{matches}"
    )

    print(
        f"mismatches                   = "
        f"{mismatches}"
    )

    print(
        f"factor shortcuts             = "
        f"{factor_shortcuts}"
    )

    print(
        f"timeouts                     = "
        f"{timeouts}"
    )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
