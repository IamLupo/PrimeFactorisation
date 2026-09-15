import math
import random
from sympy import isprime


EXPERIMENT = 115

CASES = 100

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_PROBES = 16

# Only the immediate neighbors.
K_OFFSETS = (
    -2,
    -1,
    0,
    1,
    2,
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


def value_gcd(n, k, z):
    """
    v_k = -C(k-2, 2z-1)
    for 2z-1 < k-1.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d) % n
    value = (-value) % n

    return math.gcd(value, n)


def van_der_corput_sequence(max_z, count):
    """
    Deterministic low-discrepancy sequence.
    """

    count = min(count, max_z)

    bits = max(
        1,
        (count - 1).bit_length(),
    )

    result = []
    used = set()

    def reverse_bits(value):
        output = 0

        for _ in range(bits):
            output <<= 1
            output |= value & 1
            value >>= 1

        return output

    index = 0

    while len(result) < count:
        reversed_value = reverse_bits(index)

        z = (
            reversed_value * max_z
            // (1 << bits)
        ) + 1

        if z < 1:
            z = 1

        if z > max_z:
            z = max_z

        if z not in used:
            used.add(z)
            result.append(z)

        index += 1

    return result


def search_k(n, p, q, k):
    """
    Search one k using at most MAX_PROBES z evaluations.
    """

    max_z = k // 2

    if max_z <= 0:
        return None

    probes = van_der_corput_sequence(
        max_z,
        MAX_PROBES,
    )

    for index, z in enumerate(
        probes,
        start=1,
    ):
        g = value_gcd(
            n,
            k,
            z,
        )

        state = gcd_state(
            g,
            n,
            p,
            q,
        )

        if state == "p" or state == "q":
            return {
                "z": z,
                "gcd": g,
                "state": state,
                "probes": index,
            }

    return None


def generate_semiprime():
    while True:
        p = random.randrange(
            P_MIN,
            P_MAX + 1,
        )

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, P_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        return p, q


def powers_between(p, q):
    powers = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        powers.append(k)
        k *= 2

    return powers


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q MAX: {Q_MAX}")
    print(f"MAX PROBES: {MAX_PROBES}")

    print()
    print("K OFFSETS:")
    print(
        " ".join(
            f"{x:+d}"
            for x in K_OFFSETS
        )
    )

    print()

    random.seed(115)

    total_power_cases = 0
    original_hits = 0
    original_misses = 0

    rescued = 0
    still_failed = 0

    rescue_examples = []
    failure_examples = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        powers = powers_between(p, q)

        for power_k in powers:
            total_power_cases += 1

            original = search_k(
                n,
                p,
                q,
                power_k,
            )

            if original is not None:
                original_hits += 1
                continue

            original_misses += 1

            found = None

            for offset in K_OFFSETS:
                if offset == 0:
                    continue

                k = power_k + offset

                if k < 2:
                    continue

                result = search_k(
                    n,
                    p,
                    q,
                    k,
                )

                if result is not None:
                    found = (
                        k,
                        result,
                    )
                    break

            if found is not None:
                rescued += 1

                if len(rescue_examples) < 20:
                    rescue_examples.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            power_k,
                            found[0],
                            found[1]["z"],
                            found[1]["gcd"],
                            found[1]["probes"],
                        )
                    )
            else:
                still_failed += 1

                if len(failure_examples) < 20:
                    failure_examples.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            power_k,
                        )
                    )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = {CASES}"
    )

    print(
        f"power-of-two k cases         = "
        f"{total_power_cases}"
    )

    print(
        f"original hits                = "
        f"{original_hits}"
    )

    print(
        f"original misses              = "
        f"{original_misses}"
    )

    print(
        f"rescued by nearby k          = "
        f"{rescued}"
    )

    print(
        f"still failed                 = "
        f"{still_failed}"
    )

    if total_power_cases:
        print(
            f"original hit rate            = "
            f"{100.0 * original_hits / total_power_cases:.2f}%"
        )

        print(
            f"rescue rate among misses     = "
            f"{100.0 * rescued / original_misses:.2f}%"
            if original_misses
            else "rescue rate among misses     = N/A"
        )

        print(
            f"final recovered rate         = "
            f"{100.0 * (original_hits + rescued) / total_power_cases:.2f}%"
        )

    print()
    print("-" * 60)
    print("RESCUE EXAMPLES")
    print("-" * 60)

    if rescue_examples:
        for (
            case_index,
            n,
            p,
            q,
            power_k,
            rescue_k,
            z,
            g,
            probes,
        ) in rescue_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"power_k={power_k} "
                f"rescue_k={rescue_k} "
                f"z={z} "
                f"gcd={g} "
                f"probes={probes}"
            )
    else:
        print("NONE")

    print()
    print("-" * 60)
    print("STILL FAILED")
    print("-" * 60)

    if failure_examples:
        for (
            case_index,
            n,
            p,
            q,
            power_k,
        ) in failure_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"power_k={power_k}"
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