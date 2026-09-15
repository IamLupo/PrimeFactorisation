import math
import random
import time

from sympy import isprime, prevprime


EXPERIMENT = 120

# Cases per scale tier.
CASES_PER_TIER = 25

MAX_Z_PROBES = 16

# Abort a tier if its total runtime exceeds this.
TIER_TIME_LIMIT = 30.0


# ------------------------------------------------------------
# SCALE TIERS
# ------------------------------------------------------------

TIERS = (
    {
        "name": "10^4",
        "p_min": 10_000,
        "p_max": 100_000,
        "q_min": 100_000,
        "q_max": 1_000_000,
    },
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
)


def gcd_state(g, n, p, q):
    if g == 1:
        return "1"

    if g == n:
        return "n"

    if g == p:
        return "p"

    if g == q:
        return "q"

    return "other"


def bit_reverse(value, bits):
    result = 0

    for _ in range(bits):
        result <<= 1
        result |= value & 1
        value >>= 1

    return result


def van_der_corput_sequence(max_z, count):
    """
    Deterministic low-discrepancy probe sequence.
    """

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

        z = max(1, min(z, max_z))

        if z not in used:
            used.add(z)
            result.append(z)

        index += 1

    return result


def value_gcd_timed(n, k, z):
    """
    Time the two expensive parts separately:

        math.comb()
        gcd()
    """

    d = 2 * z - 1

    if d >= k - 1:
        return (
            n,
            0.0,
            0.0,
        )

    start = time.perf_counter()

    value = math.comb(
        k - 2,
        d,
    )

    value %= n
    value = (-value) % n

    comb_time = (
        time.perf_counter()
        - start
    )

    start = time.perf_counter()

    g = math.gcd(
        value,
        n,
    )

    gcd_time = (
        time.perf_counter()
        - start
    )

    return (
        g,
        comb_time,
        gcd_time,
    )


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def search_k_timed(n, k):
    """
    Run the normal VdC z search while collecting timings.
    """

    max_z = k // 2

    if max_z <= 0:
        return {
            "found": False,
            "probes": 0,
            "comb_time": 0.0,
            "gcd_time": 0.0,
            "k": k,
            "z": None,
            "gcd": None,
        }

    probes = van_der_corput_sequence(
        max_z,
        MAX_Z_PROBES,
    )

    total_comb_time = 0.0
    total_gcd_time = 0.0

    for index, z in enumerate(
        probes,
        start=1,
    ):
        (
            g,
            comb_time,
            gcd_time,
        ) = value_gcd_timed(
            n,
            k,
            z,
        )

        total_comb_time += comb_time
        total_gcd_time += gcd_time

        if 1 < g < n:
            return {
                "found": True,
                "probes": index,
                "comb_time": total_comb_time,
                "gcd_time": total_gcd_time,
                "k": k,
                "z": z,
                "gcd": g,
            }

    return {
        "found": False,
        "probes": len(probes),
        "comb_time": total_comb_time,
        "gcd_time": total_gcd_time,
        "k": k,
        "z": None,
        "gcd": None,
    }


def dyadic_stage_timed(n, p, q):
    start = time.perf_counter()

    k_tests = 0
    probes = 0
    comb_time = 0.0
    gcd_time = 0.0

    for k in powers_between(p, q):
        k_tests += 1

        result = search_k_timed(
            n,
            k,
        )

        probes += result["probes"]
        comb_time += result["comb_time"]
        gcd_time += result["gcd_time"]

        if result["found"]:
            result["stage_time"] = (
                time.perf_counter()
                - start
            )

            result["k_tests"] = k_tests

            return result

    return {
        "found": False,
        "stage_time": (
            time.perf_counter()
            - start
        ),
        "k_tests": k_tests,
        "probes": probes,
        "comb_time": comb_time,
        "gcd_time": gcd_time,
        "k": None,
        "z": None,
        "gcd": None,
    }


def sqrt_prime_stage_timed(n):
    start = time.perf_counter()

    sqrt_n = math.isqrt(n)

    r = int(
        prevprime(sqrt_n)
    )

    prevprime_time = (
        time.perf_counter()
        - start
    )

    result = search_k_timed(
        n,
        r,
    )

    result["r"] = r
    result["prevprime_time"] = (
        prevprime_time
    )
    result["stage_time"] = (
        time.perf_counter()
        - start
    )

    return result


def generate_semiprime(
    p_min,
    p_max,
    q_min,
    q_max,
):
    start = time.perf_counter()

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

        if p == q:
            continue

        elapsed = (
            time.perf_counter()
            - start
        )

        return p, q, elapsed


def run_case(
    p_min,
    p_max,
    q_min,
    q_max,
):
    generation_start = time.perf_counter()

    (
        p,
        q,
        generation_time,
    ) = generate_semiprime(
        p_min,
        p_max,
        q_min,
        q_max,
    )

    n = p * q

    case_start = time.perf_counter()

    stage1 = dyadic_stage_timed(
        n,
        p,
        q,
    )

    if stage1["found"]:
        return {
            "p": p,
            "q": q,
            "n": n,
            "stage": "dyadic",
            "generation_time": generation_time,
            "total_time": (
                time.perf_counter()
                - case_start
            ),
            "stage_time": stage1["stage_time"],
            "k_tests": stage1["k_tests"],
            "probes": stage1["probes"],
            "comb_time": stage1["comb_time"],
            "gcd_time": stage1["gcd_time"],
            "prevprime_time": 0.0,
        }

    stage2 = sqrt_prime_stage_timed(
        n,
    )

    if stage2["found"]:
        return {
            "p": p,
            "q": q,
            "n": n,
            "stage": "sqrt_prime",
            "generation_time": generation_time,
            "total_time": (
                time.perf_counter()
                - case_start
            ),
            "stage_time": (
                stage1["stage_time"]
                + stage2["stage_time"]
            ),
            "k_tests": (
                stage1["k_tests"]
            ),
            "probes": (
                stage1["probes"]
                + stage2["probes"]
            ),
            "comb_time": (
                stage1["comb_time"]
                + stage2["comb_time"]
            ),
            "gcd_time": (
                stage1["gcd_time"]
                + stage2["gcd_time"]
            ),
            "prevprime_time": (
                stage2["prevprime_time"]
            ),
        }

    return {
        "p": p,
        "q": q,
        "n": n,
        "stage": "unrecovered",
        "generation_time": generation_time,
        "total_time": (
            time.perf_counter()
            - case_start
        ),
        "stage_time": (
            stage1["stage_time"]
            + stage2["stage_time"]
        ),
        "k_tests": (
            stage1["k_tests"]
        ),
        "probes": (
            stage1["probes"]
            + stage2["probes"]
        ),
        "comb_time": (
            stage1["comb_time"]
            + stage2["comb_time"]
        ),
        "gcd_time": (
            stage1["gcd_time"]
            + stage2["gcd_time"]
        ),
        "prevprime_time": (
            stage2["prevprime_time"]
        ),
    }


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES PER TIER: {CASES_PER_TIER}")
    print(f"MAX Z PROBES: {MAX_Z_PROBES}")
    print(f"TIER TIME LIMIT: {TIER_TIME_LIMIT}s")
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
        dyadic_hits = 0
        sqrt_prime_hits = 0
        failures = 0

        total_generation = 0.0
        total_time = 0.0
        total_stage_time = 0.0
        total_comb_time = 0.0
        total_gcd_time = 0.0
        total_prevprime_time = 0.0

        total_k_tests = 0
        total_probes = 0

        maximum_case_time = 0.0
        maximum_comb_time = 0.0
        slowest_case = None

        aborted = False

        for case_index in range(
            1,
            CASES_PER_TIER + 1,
        ):
            result = run_case(
                tier["p_min"],
                tier["p_max"],
                tier["q_min"],
                tier["q_max"],
            )

            completed += 1

            total_generation += (
                result["generation_time"]
            )

            total_time += (
                result["total_time"]
            )

            total_stage_time += (
                result["stage_time"]
            )

            total_comb_time += (
                result["comb_time"]
            )

            total_gcd_time += (
                result["gcd_time"]
            )

            total_prevprime_time += (
                result["prevprime_time"]
            )

            total_k_tests += (
                result["k_tests"]
            )

            total_probes += (
                result["probes"]
            )

            if result["stage"] == "dyadic":
                dyadic_hits += 1

            elif result["stage"] == "sqrt_prime":
                sqrt_prime_hits += 1

            else:
                failures += 1

            if (
                result["total_time"]
                > maximum_case_time
            ):
                maximum_case_time = (
                    result["total_time"]
                )

                slowest_case = result

            maximum_comb_time = max(
                maximum_comb_time,
                result["comb_time"],
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
            f"{tier_elapsed:.4f}s"
        )

        print(
            f"aborted                     = "
            f"{aborted}"
        )

        if completed:
            print(
                f"average case time          = "
                f"{total_time / completed:.6f}s"
            )

            print(
                f"average generation time    = "
                f"{total_generation / completed:.6f}s"
            )

            print(
                f"average stage time         = "
                f"{total_stage_time / completed:.6f}s"
            )

            print(
                f"average comb time          = "
                f"{total_comb_time / completed:.6f}s"
            )

            print(
                f"average gcd time           = "
                f"{total_gcd_time / completed:.6f}s"
            )

        print(
            f"dyadic hits                = "
            f"{dyadic_hits}"
        )

        print(
            f"sqrt-prime hits            = "
            f"{sqrt_prime_hits}"
        )

        print(
            f"failures                   = "
            f"{failures}"
        )

        print(
            f"total k tests              = "
            f"{total_k_tests}"
        )

        print(
            f"total z probes             = "
            f"{total_probes}"
        )

        print(
            f"maximum case time          = "
            f"{maximum_case_time:.6f}s"
        )

        print(
            f"maximum comb time          = "
            f"{maximum_comb_time:.6f}s"
        )

        if slowest_case is not None:
            print()
            print("SLOWEST CASE")

            print(
                f"  N={slowest_case['n']}"
            )

            print(
                f"  p={slowest_case['p']}"
            )

            print(
                f"  q={slowest_case['q']}"
            )

            print(
                f"  stage={slowest_case['stage']}"
            )

            print(
                f"  time={slowest_case['total_time']:.6f}s"
            )

            print(
                f"  probes={slowest_case['probes']}"
            )

        print()

        if aborted:
            print(
                "STOPPING EXPERIMENT BECAUSE "
                "THIS TIER HIT THE TIME LIMIT."
            )
            break

    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
