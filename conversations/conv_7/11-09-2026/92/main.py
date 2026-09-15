import math
import random
from sympy import isprime


EXPERIMENT = 113

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384

PROBE_BUDGETS = (
    1,
    2,
    4,
    8,
    16,
    32,
    64,
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
    For

        d = 2z - 1

    we have

        v_k = -C(k-2,d)

    whenever d < k-1.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d)
    value %= n
    value = (-value) % n

    return math.gcd(value, n)


def lucas_zero(m, d, p):
    """
    Lucas theorem:

        C(m,d) == 0 (mod p)

    iff some base-p digit of d exceeds
    the corresponding digit of m.
    """

    while m > 0 or d > 0:
        md = m % p
        dd = d % p

        if dd > md:
            return True

        m //= p
        d //= p

    return False


def count_p_hits(k, p):
    """
    Count every z for which p divides v_k(z).

    This uses Lucas directly rather than gcd evaluation.
    """

    max_z = min(
        k // 2,
        MAX_Z,
    )

    count = 0

    for z in range(1, max_z + 1):
        d = 2 * z - 1

        if d >= k - 1:
            break

        if lucas_zero(k - 2, d, p):
            count += 1

    return count


def p_hit_intervals(k, p):
    """
    Determine contiguous intervals of z where p divides v_k.

    This is used only to characterize the geometry.
    """

    max_z = min(
        k // 2,
        MAX_Z,
    )

    intervals = []

    start = None
    previous = None

    for z in range(1, max_z + 1):
        d = 2 * z - 1

        if d >= k - 1:
            break

        hit = lucas_zero(
            k - 2,
            d,
            p,
        )

        if hit:
            if start is None:
                start = z

            previous = z

        else:
            if start is not None:
                intervals.append(
                    (
                        start,
                        previous,
                    )
                )

                start = None
                previous = None

    if start is not None:
        intervals.append(
            (
                start,
                previous,
            )
        )

    return intervals


def base_p_digits(value, p):
    digits = []

    while value > 0:
        digits.append(value % p)
        value //= p

    if not digits:
        digits.append(0)

    return digits


def smallest_lucas_zero_d(m, p):
    """
    Construct the smallest positive d satisfying Lucas zero.
    """

    digits = base_p_digits(m, p)

    candidates = []

    power = 1

    for digit in digits:
        if digit + 1 < p:
            candidates.append(
                (digit + 1) * power
            )

        power *= p

    if not candidates:
        return None

    return min(candidates)


def exact_first_z(p, k):
    d = smallest_lucas_zero_d(
        k - 2,
        p,
    )

    if d is None:
        return None

    if d % 2 == 0:
        d += 1

    if d >= k - 1:
        return None

    return (d + 1) // 2


def bit_reverse(value, bits):
    result = 0

    for _ in range(bits):
        result <<= 1
        result |= value & 1
        value >>= 1

    return result


def van_der_corput_sequence(max_z, count):
    """
    Deterministic low-discrepancy sequence.
    """

    count = min(
        count,
        max_z,
    )

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

        if z < 1:
            z = 1

        if z > max_z:
            z = max_z

        if z not in used:
            used.add(z)
            result.append(z)

        index += 1

    return result


def random_sequence(max_z, count, seed):
    rng = random.Random(seed)

    count = min(
        count,
        max_z,
    )

    return rng.sample(
        range(1, max_z + 1),
        count,
    )


def evaluate_sequence(
    n,
    p,
    q,
    k,
    sequence,
):
    found = False

    first_hit_z = None

    for z in sequence:
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
            if not found:
                found = True
                first_hit_z = z

    return found, first_hit_z


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


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q MAX: {Q_MAX}")
    print(f"MAX Z: {MAX_Z}")

    print()
    print("PROBE BUDGETS:")
    print(
        " ".join(
            str(x)
            for x in PROBE_BUDGETS
        )
    )

    print()

    random.seed(113)

    total_k_cases = 0

    density_sum = 0.0

    minimum_density = 1.0
    maximum_density = 0.0

    density_below_1pct = 0
    density_below_5pct = 0
    density_below_10pct = 0
    density_above_25pct = 0

    vdc_hits = {
        budget: 0
        for budget in PROBE_BUDGETS
    }

    random_hits = {
        budget: 0
        for budget in PROBE_BUDGETS
    }

    hard_cases = []

    first_hit_ratios = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            max_z = min(
                k // 2,
                MAX_Z,
            )

            hit_count = count_p_hits(
                k,
                p,
            )

            density = (
                hit_count / max_z
                if max_z
                else 0.0
            )

            density_sum += density

            minimum_density = min(
                minimum_density,
                density,
            )

            maximum_density = max(
                maximum_density,
                density,
            )

            if density < 0.01:
                density_below_1pct += 1

            if density < 0.05:
                density_below_5pct += 1

            if density < 0.10:
                density_below_10pct += 1

            if density > 0.25:
                density_above_25pct += 1

            exact_z = exact_first_z(
                p,
                k,
            )

            # -------------------------------------------------
            # VDC
            # -------------------------------------------------

            vdc_sequence = van_der_corput_sequence(
                max_z,
                max(PROBE_BUDGETS),
            )

            # -------------------------------------------------
            # RANDOM
            # -------------------------------------------------

            seed = (
                case_index * 1000003
                + n * 1009
                + k * 9176
            )

            random_sequence_values = (
                random_sequence(
                    max_z,
                    max(PROBE_BUDGETS),
                    seed,
                )
            )

            vdc_found_at = None
            random_found_at = None

            for budget in PROBE_BUDGETS:
                vdc_found, vdc_z = (
                    evaluate_sequence(
                        n,
                        p,
                        q,
                        k,
                        vdc_sequence[:budget],
                    )
                )

                if vdc_found:
                    vdc_hits[budget] += 1

                    if vdc_found_at is None:
                        vdc_found_at = (
                            budget,
                            vdc_z,
                        )

                random_found, random_z = (
                    evaluate_sequence(
                        n,
                        p,
                        q,
                        k,
                        random_sequence_values[
                            :budget
                        ],
                    )
                )

                if random_found:
                    random_hits[budget] += 1

                    if random_found_at is None:
                        random_found_at = (
                            budget,
                            random_z,
                        )

            if vdc_found_at is None:
                hard_cases.append(
                    {
                        "case": case_index,
                        "n": n,
                        "p": p,
                        "q": q,
                        "k": k,
                        "max_z": max_z,
                        "hit_count": hit_count,
                        "density": density,
                        "exact_z": exact_z,
                    }
                )

            if vdc_found_at is not None:
                budget, hit_z = vdc_found_at

                first_hit_ratios.append(
                    hit_z / max_z
                )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = "
        f"{CASES}"
    )

    print(
        f"(N,k) cases                  = "
        f"{total_k_cases}"
    )

    print()

    print(
        f"average p-hit density        = "
        f"{100.0 * density_sum / total_k_cases:.4f}%"
    )

    print(
        f"minimum p-hit density        = "
        f"{100.0 * minimum_density:.6f}%"
    )

    print(
        f"maximum p-hit density        = "
        f"{100.0 * maximum_density:.4f}%"
    )

    print()

    print(
        f"density < 1% cases           = "
        f"{density_below_1pct}"
    )

    print(
        f"density < 5% cases           = "
        f"{density_below_5pct}"
    )

    print(
        f"density < 10% cases          = "
        f"{density_below_10pct}"
    )

    print(
        f"density > 25% cases          = "
        f"{density_above_25pct}"
    )

    print()
    print("-" * 60)
    print("VAN DER CORPUT")
    print("-" * 60)

    for budget in PROBE_BUDGETS:
        hits = vdc_hits[budget]

        print(
            f"probes={budget:2d} "
            f"hits={hits:4d} "
            f"rate={100.0 * hits / total_k_cases:6.2f}%"
        )

    print()
    print("-" * 60)
    print("RANDOM")
    print("-" * 60)

    for budget in PROBE_BUDGETS:
        hits = random_hits[budget]

        print(
            f"probes={budget:2d} "
            f"hits={hits:4d} "
            f"rate={100.0 * hits / total_k_cases:6.2f}%"
        )

    print()
    print("-" * 60)
    print("VDC HARD CASES")
    print("-" * 60)

    print(
        f"VDC misses at 64 probes      = "
        f"{len(hard_cases)}"
    )

    for item in hard_cases[:20]:
        print(
            f"CASE={item['case']} "
            f"N={item['n']} "
            f"p={item['p']} "
            f"q={item['q']} "
            f"k={item['k']} "
            f"max_z={item['max_z']} "
            f"hits={item['hit_count']} "
            f"density={100.0 * item['density']:.6f}% "
            f"exact_z={item['exact_z']}"
        )

    print()
    print("-" * 60)
    print("VDC HIT POSITION")
    print("-" * 60)

    if first_hit_ratios:
        print(
            f"minimum hit z/max_z        = "
            f"{min(first_hit_ratios):.6f}"
        )

        print(
            f"maximum hit z/max_z        = "
            f"{max(first_hit_ratios):.6f}"
        )

        print(
            f"average hit z/max_z        = "
            f"{sum(first_hit_ratios) / len(first_hit_ratios):.6f}"
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
