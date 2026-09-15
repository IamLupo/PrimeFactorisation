import math
import random
import time

from sympy import isprime, prevprime


EXPERIMENT = 121

CASES = 50

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

MAX_PROBES = 16


def value_gcd_timed(n, k, z):
    """
    Compute gcd(v_k(z), n) and time only the expensive
    binomial computation.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return n, 0.0, 0

    start = time.perf_counter()

    value = math.comb(
        k - 2,
        d,
    )

    comb_time = (
        time.perf_counter()
        - start
    )

    value %= n
    value = (-value) % n

    g = math.gcd(
        value,
        n,
    )

    return g, comb_time, d


def is_factor(g, n):
    return 1 < g < n


def vdc_sequence(max_z, count):
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

    def reverse_bits(value):
        result = 0

        for _ in range(bits):
            result <<= 1
            result |= value & 1
            value >>= 1

        return result

    result = []
    used = set()
    index = 0

    while len(result) < count:
        r = reverse_bits(index)

        z = (
            r * max_z
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


def power_sequence(max_z, count):
    """
    1,2,4,8,...
    """

    result = []

    z = 1

    while (
        len(result) < count
        and z <= max_z
    ):
        result.append(z)
        z *= 2

    return result


def geometric_sequence(max_z, count):
    """
    Small-z geometric sequence.

    Starts densely, then gradually expands.
    """

    result = []

    candidates = (
        1,
        2,
        3,
        4,
        6,
        8,
        12,
        16,
        24,
        32,
        48,
        64,
        96,
        128,
        192,
        256,
        384,
        512,
        768,
        1024,
        1536,
        2048,
        3072,
        4096,
        6144,
        8192,
        12288,
        16384,
        24576,
        32768,
    )

    for z in candidates:
        if z > max_z:
            break

        result.append(z)

        if len(result) >= count:
            break

    return result


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def search_sequence(
    n,
    p,
    q,
    k,
    sequence,
):
    total_comb_time = 0.0
    total_probes = 0
    maximum_d = 0

    for z in sequence:
        if z > k // 2:
            continue

        (
            g,
            comb_time,
            d,
        ) = value_gcd_timed(
            n,
            k,
            z,
        )

        total_comb_time += comb_time
        total_probes += 1

        maximum_d = max(
            maximum_d,
            d,
        )

        if is_factor(g, n):
            return {
                "found": True,
                "z": z,
                "gcd": g,
                "probes": total_probes,
                "comb_time": total_comb_time,
                "max_d": maximum_d,
            }

    return {
        "found": False,
        "z": None,
        "gcd": None,
        "probes": total_probes,
        "comb_time": total_comb_time,
        "max_d": maximum_d,
    }


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


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX}")
    print(f"MAX PROBES: {MAX_PROBES}")
    print()

    random.seed(EXPERIMENT)

    strategies = (
        "VDC",
        "POWER",
        "GEOMETRIC",
    )

    stats = {
        strategy: {
            "hits": 0,
            "probes": 0,
            "comb_time": 0.0,
            "max_d": 0,
        }
        for strategy in strategies
    }

    examples = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        # We test the first useful power-of-two k only.
        ks = powers_between(p, q)

        if not ks:
            continue

        k = ks[0]

        max_z = k // 2

        sequences = {
            "VDC": vdc_sequence(
                max_z,
                MAX_PROBES,
            ),
            "POWER": power_sequence(
                max_z,
                MAX_PROBES,
            ),
            "GEOMETRIC": geometric_sequence(
                max_z,
                MAX_PROBES,
            ),
        }

        results = {}

        for strategy in strategies:
            result = search_sequence(
                n,
                p,
                q,
                k,
                sequences[strategy],
            )

            results[strategy] = result

            stats[strategy]["probes"] += (
                result["probes"]
            )

            stats[strategy]["comb_time"] += (
                result["comb_time"]
            )

            stats[strategy]["max_d"] = max(
                stats[strategy]["max_d"],
                result["max_d"],
            )

            if result["found"]:
                stats[strategy]["hits"] += 1

        if len(examples) < 15:
            examples.append(
                (
                    case_index,
                    n,
                    p,
                    q,
                    k,
                    results,
                )
            )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    for strategy in strategies:
        item = stats[strategy]

        print()
        print(f"[{strategy}]")

        print(
            f"hits                        = "
            f"{item['hits']}/{CASES}"
        )

        print(
            f"hit rate                    = "
            f"{100.0 * item['hits'] / CASES:.2f}%"
        )

        print(
            f"average probes              = "
            f"{item['probes'] / CASES:.2f}"
        )

        print(
            f"total comb time             = "
            f"{item['comb_time']:.4f}s"
        )

        print(
            f"average comb time           = "
            f"{item['comb_time'] / CASES:.6f}s"
        )

        print(
            f"maximum d evaluated         = "
            f"{item['max_d']}"
        )

    print()
    print("-" * 60)
    print("EXAMPLES")
    print("-" * 60)

    for (
        case_index,
        n,
        p,
        q,
        k,
        results,
    ) in examples:
        print()
        print(
            f"CASE={case_index} "
            f"N={n} "
            f"p={p} "
            f"q={q} "
            f"k={k}"
        )

        for strategy in strategies:
            result = results[strategy]

            print(
                f"  {strategy:9s} "
                f"found={result['found']} "
                f"z={result['z']} "
                f"gcd={result['gcd']} "
                f"probes={result['probes']} "
                f"d={result['max_d']} "
                f"comb={result['comb_time']:.6f}s"
            )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
