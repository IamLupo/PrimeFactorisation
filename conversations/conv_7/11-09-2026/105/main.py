import math
import random
import time

from sympy import isprime


EXPERIMENT = 130

CASES_PER_TIER = 5

MAX_PROBES = 2

TIER_TIME_LIMIT = 15.0

# Scale the smaller prime.
TIERS = (
    {
        "name": "10^5",
        "p_min": 100_000,
        "p_max": 1_000_000,
        "q_min": 1_000_000,
        "q_max": 10_000_000,
    },
    {
        "name": "10^6",
        "p_min": 1_000_000,
        "p_max": 10_000_000,
        "q_min": 10_000_000,
        "q_max": 100_000_000,
    },
    {
        "name": "10^7",
        "p_min": 10_000_000,
        "p_max": 100_000_000,
        "q_min": 100_000_000,
        "q_max": 1_000_000_000,
    },
)


def generate_semiprime(
    p_min,
    p_max,
    q_min,
    q_max,
):
    while True:
        p = random.randrange(
            p_min,
            p_max + 1,
        )

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, q_min),
            q_max + 1,
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


def van_der_corput_sequence(
    max_z,
    count,
):
    count = min(
        count,
        max_z,
    )

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


def numerator_gcd(
    n,
    k,
    z,
):
    """
    Compute the gcd using only the numerator:

        prod_{x=k-1-d}^{k-2} x

    where

        d = 2z-1.

    Stop immediately when a nontrivial gcd appears.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return {
            "status": "trivial_n",
            "factor": n,
            "iterations": 0,
            "time": 0.0,
        }

    first = k - 1 - d
    last = k - 2

    start = time.perf_counter()

    product = 1

    for iteration, x in enumerate(
        range(first, last + 1),
        start=1,
    ):
        product *= x
        product %= n

        g = math.gcd(
            product,
            n,
        )

        if 1 < g < n:
            return {
                "status": "factor",
                "factor": g,
                "iterations": iteration,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

        if (
            time.perf_counter()
            - start
            >= TIER_TIME_LIMIT
        ):
            return {
                "status": "timeout",
                "factor": None,
                "iterations": iteration,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

    return {
        "status": "no_factor",
        "factor": None,
        "iterations": d,
        "time": (
            time.perf_counter()
            - start
        ),
    }


def test_case(p, q):
    n = p * q

    ks = powers_between(
        p,
        q,
    )

    if not ks:
        return {
            "found": False,
            "stage": "no_k",
            "time": 0.0,
            "probes": 0,
            "iterations": 0,
        }

    start = time.perf_counter()

    # Only test the first power-of-two k.
    k = ks[0]

    probes = van_der_corput_sequence(
        k // 2,
        MAX_PROBES,
    )

    total_iterations = 0

    for probe_index, z in enumerate(
        probes,
        start=1,
    ):
        result = numerator_gcd(
            n,
            k,
            z,
        )

        total_iterations += (
            result["iterations"]
        )

        if result["status"] == "factor":
            factor = result["factor"]

            if n % factor == 0:
                return {
                    "found": True,
                    "stage": "numerator",
                    "factor": factor,
                    "k": k,
                    "z": z,
                    "probe": probe_index,
                    "iterations": total_iterations,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                }

        elif result["status"] == "timeout":
            return {
                "found": False,
                "stage": "timeout",
                "k": k,
                "z": z,
                "probe": probe_index,
                "iterations": total_iterations,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

    return {
        "found": False,
        "stage": "miss",
        "k": k,
        "probe": len(probes),
        "iterations": total_iterations,
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

    print(
        f"CASES PER TIER: "
        f"{CASES_PER_TIER}"
    )

    print(
        f"MAX PROBES: "
        f"{MAX_PROBES}"
    )

    print(
        f"TIER TIME LIMIT: "
        f"{TIER_TIME_LIMIT}s"
    )

    print()

    random.seed(EXPERIMENT)

    for tier in TIERS:
        print("-" * 60)
        print(
            f"TIER {tier['name']}"
        )
        print("-" * 60)

        tier_start = time.perf_counter()

        completed = 0
        hits = 0
        timeouts = 0
        misses = 0

        total_time = 0.0
        total_iterations = 0

        maximum_time = 0.0
        maximum_iterations = 0

        slowest = None

        aborted = False

        for case_index in range(
            1,
            CASES_PER_TIER + 1,
        ):
            p, q = generate_semiprime(
                tier["p_min"],
                tier["p_max"],
                tier["q_min"],
                tier["q_max"],
            )

            n = p * q

            result = test_case(
                p,
                q,
            )

            completed += 1

            total_time += result[
                "time"
            ]

            total_iterations += result[
                "iterations"
            ]

            maximum_time = max(
                maximum_time,
                result["time"],
            )

            maximum_iterations = max(
                maximum_iterations,
                result["iterations"],
            )

            if result["found"]:
                hits += 1

            elif result["stage"] == "timeout":
                timeouts += 1

            else:
                misses += 1

            if (
                result["time"]
                >= maximum_time
            ):
                slowest = (
                    case_index,
                    p,
                    q,
                    n,
                    result,
                )

            elapsed = (
                time.perf_counter()
                - tier_start
            )

            if elapsed >= TIER_TIME_LIMIT:
                aborted = True
                break

        tier_elapsed = (
            time.perf_counter()
            - tier_start
        )

        print(
            f"completed cases             = "
            f"{completed}"
        )

        print(
            f"tier time                   = "
            f"{tier_elapsed:.6f}s"
        )

        print(
            f"aborted                     = "
            f"{aborted}"
        )

        print(
            f"hits                        = "
            f"{hits}"
        )

        print(
            f"misses                      = "
            f"{misses}"
        )

        print(
            f"timeouts                    = "
            f"{timeouts}"
        )

        if completed:
            print(
                f"average case time          = "
                f"{total_time / completed:.6f}s"
            )

            print(
                f"average iterations         = "
                f"{total_iterations / completed:.2f}"
            )

        print(
            f"maximum case time          = "
            f"{maximum_time:.6f}s"
        )

        print(
            f"maximum iterations         = "
            f"{maximum_iterations}"
        )

        if slowest is not None:
            (
                case_index,
                p,
                q,
                n,
                result,
            ) = slowest

            print()
            print(
                "SLOWEST CASE"
            )

            print(
                f"  CASE={case_index}"
            )

            print(
                f"  N={n}"
            )

            print(
                f"  p={p}"
            )

            print(
                f"  q={q}"
            )

            print(
                f"  k={result.get('k')}"
            )

            print(
                f"  z={result.get('z')}"
            )

            print(
                f"  probe={result.get('probe')}"
            )

            print(
                f"  stage={result['stage']}"
            )

            print(
                f"  iterations={result['iterations']}"
            )

            print(
                f"  time={result['time']:.6f}s"
            )

        print()

        if aborted:
            print(
                "STOPPING: tier exceeded "
                "time limit."
            )
            break

    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
