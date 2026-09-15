import math
import random
import time

from sympy import isprime


EXPERIMENT = 124

CASES = 1000

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

MAX_PROBES = 16


def lucas_zero(m, d, p):
    """
    Lucas theorem:

        C(m,d) == 0 (mod p)

    iff some base-p digit of d is greater than
    the corresponding base-p digit of m.
    """

    while m > 0 or d > 0:
        m_digit = m % p
        d_digit = d % p

        if d_digit > m_digit:
            return True

        m //= p
        d //= p

    return False


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


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def predict_factor(p, k, z):
    """
    Predict whether p divides v_k(z).
    """

    d = 2 * z - 1

    if d >= k - 1:
        return True

    return lucas_zero(
        k - 2,
        d,
        p,
    )


def run_search(p, q, k):
    """
    Search using only Lucas arithmetic.

    No giant binomial.
    No modular inverses.
    No gcd.
    """

    max_z = k // 2

    probes = van_der_corput_sequence(
        max_z,
        MAX_PROBES,
    )

    for index, z in enumerate(
        probes,
        start=1,
    ):
        if predict_factor(
            p,
            k,
            z,
        ):
            return {
                "found": True,
                "z": z,
                "probes": index,
            }

    return {
        "found": False,
        "z": None,
        "probes": len(probes),
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
    print()

    random.seed(EXPERIMENT)

    total_cases = 0
    total_k_cases = 0

    hits = 0
    misses = 0

    probes_total = 0

    total_time = 0.0

    hit_by_probe = {
        i: 0
        for i in (
            1,
            2,
            4,
            8,
            16,
        )
    }

    hard_examples = []

    experiment_start = time.perf_counter()

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()

        n = p * q

        total_cases += 1

        for k in powers_between(
            p,
            q,
        ):
            total_k_cases += 1

            start = time.perf_counter()

            result = run_search(
                p,
                q,
                k,
            )

            elapsed = (
                time.perf_counter()
                - start
            )

            total_time += elapsed

            probes_total += result[
                "probes"
            ]

            if result["found"]:
                hits += 1

                for budget in hit_by_probe:
                    if (
                        result["probes"]
                        <= budget
                    ):
                        hit_by_probe[
                            budget
                        ] += 1
            else:
                misses += 1

                if len(hard_examples) < 20:
                    hard_examples.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                        )
                    )

    total_elapsed = (
        time.perf_counter()
        - experiment_start
    )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = "
        f"{total_cases}"
    )

    print(
        f"(N,k) cases                  = "
        f"{total_k_cases}"
    )

    print(
        f"hits                         = "
        f"{hits}"
    )

    print(
        f"misses                       = "
        f"{misses}"
    )

    print(
        f"hit rate                     = "
        f"{100.0 * hits / total_k_cases:.2f}%"
    )

    print(
        f"average probes               = "
        f"{probes_total / total_k_cases:.2f}"
    )

    print(
        f"total search time            = "
        f"{total_time:.6f}s"
    )

    print(
        f"total wall time              = "
        f"{total_elapsed:.6f}s"
    )

    print(
        f"average case search time     = "
        f"{total_time / total_k_cases:.9f}s"
    )

    print()

    print("-" * 60)
    print("CUMULATIVE HIT RATE")
    print("-" * 60)

    for budget in (
        1,
        2,
        4,
        8,
        16,
    ):
        count = hit_by_probe[
            budget
        ]

        print(
            f"probes <= {budget:2d}: "
            f"{count:5d} "
            f"({100.0 * count / total_k_cases:.2f}%)"
        )

    print()
    print("-" * 60)
    print("HARD CASES")
    print("-" * 60)

    if hard_examples:
        for (
            case_index,
            n,
            p,
            q,
            k,
        ) in hard_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k}"
            )
    else:
        print("NONE")

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()