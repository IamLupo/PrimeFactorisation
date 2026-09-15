import math
import random
from sympy import isprime, prevprime


EXPERIMENT = 117

CASES = 300

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_PROBES = 16


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
    For d = 2z - 1:

        v_k = -C(k-2, d)

    whenever d < k-1.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d) % n
    value = (-value) % n

    return math.gcd(value, n)


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


def search_k(n, p, q, k, max_probes):
    """
    Search one k using at most max_probes deterministic z values.
    """

    max_z = k // 2

    if max_z <= 0:
        return None

    probes = van_der_corput_sequence(
        max_z,
        max_probes,
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
                "k": k,
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


def run_dyadic_stage(n, p, q):
    """
    Original strategy:
    test every power-of-two k with p < k < q.
    """

    for k in powers_between(p, q):
        result = search_k(
            n,
            p,
            q,
            k,
            MAX_PROBES,
        )

        if result is not None:
            return result

    return None


def run_fallback_stage(n, p, q):
    """
    Proposed fallback:

        r = previous_prime(floor(sqrt(N)))
        k = r

    No modification of N is performed.
    """

    sqrt_n = math.isqrt(n)

    if sqrt_n < 3:
        return {
            "usable": False,
            "reason": "sqrt_too_small",
        }

    r = int(prevprime(sqrt_n))

    # The useful regime requires p < r < q.
    if r <= p:
        return {
            "usable": False,
            "reason": "r_le_p",
            "r": r,
        }

    if r >= q:
        return {
            "usable": False,
            "reason": "r_ge_q",
            "r": r,
        }

    result = search_k(
        n,
        p,
        q,
        r,
        MAX_PROBES,
    )

    if result is not None:
        result["r"] = r
        return {
            "usable": True,
            "found": True,
            **result,
        }

    return {
        "usable": True,
        "found": False,
        "r": r,
    }


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

    print(
        "FALLBACK: "
        "r = previous_prime(floor(sqrt(N))), k = r"
    )
    print()

    random.seed(EXPERIMENT)

    total_cases = 0

    dyadic_hits = 0
    dyadic_misses = 0

    fallback_attempts = 0
    fallback_hits = 0
    fallback_misses = 0

    r_le_p = 0
    r_ge_q = 0

    recovered_total = 0

    rescue_examples = []
    remaining_examples = []

    r_ratios = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        total_cases += 1

        # -------------------------------------------------
        # STAGE 1: ORIGINAL DYADIC SEARCH
        # -------------------------------------------------

        dyadic = run_dyadic_stage(
            n,
            p,
            q,
        )

        if dyadic is not None:
            dyadic_hits += 1
            recovered_total += 1
            continue

        dyadic_misses += 1

        # -------------------------------------------------
        # STAGE 2: PREVIOUS PRIME SQRT FALLBACK
        # -------------------------------------------------

        fallback = run_fallback_stage(
            n,
            p,
            q,
        )

        if not fallback["usable"]:
            reason = fallback["reason"]

            if reason == "r_le_p":
                r_le_p += 1

            elif reason == "r_ge_q":
                r_ge_q += 1

            if len(remaining_examples) < 20:
                remaining_examples.append(
                    (
                        case_index,
                        n,
                        p,
                        q,
                        reason,
                        fallback.get("r"),
                    )
                )

            continue

        fallback_attempts += 1

        r = fallback["r"]

        r_ratios.append(
            r / math.sqrt(n)
        )

        if fallback["found"]:
            fallback_hits += 1
            recovered_total += 1

            if len(rescue_examples) < 20:
                rescue_examples.append(
                    (
                        case_index,
                        n,
                        p,
                        q,
                        r,
                        fallback["k"],
                        fallback["z"],
                        fallback["gcd"],
                        fallback["probes"],
                    )
                )
        else:
            fallback_misses += 1

            if len(remaining_examples) < 20:
                remaining_examples.append(
                    (
                        case_index,
                        n,
                        p,
                        q,
                        "fallback_miss",
                        r,
                    )
                )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = "
        f"{total_cases}"
    )

    print(
        f"original dyadic hits         = "
        f"{dyadic_hits}"
    )

    print(
        f"original dyadic misses       = "
        f"{dyadic_misses}"
    )

    print(
        f"fallback attempts            = "
        f"{fallback_attempts}"
    )

    print(
        f"fallback hits                = "
        f"{fallback_hits}"
    )

    print(
        f"fallback misses              = "
        f"{fallback_misses}"
    )

    print(
        f"r <= p cases                 = "
        f"{r_le_p}"
    )

    print(
        f"r >= q cases                 = "
        f"{r_ge_q}"
    )

    print(
        f"total recovered              = "
        f"{recovered_total}"
    )

    print()

    print(
        f"original recovery rate       = "
        f"{100.0 * dyadic_hits / total_cases:.2f}%"
    )

    print(
        f"fallback rescue rate         = "
        f"{100.0 * fallback_hits / dyadic_misses:.2f}%"
        if dyadic_misses
        else
        "fallback rescue rate         = N/A"
    )

    print(
        f"final recovery rate          = "
        f"{100.0 * recovered_total / total_cases:.2f}%"
    )

    if r_ratios:
        print()

        print(
            f"average r/sqrt(N)            = "
            f"{sum(r_ratios) / len(r_ratios):.8f}"
        )

        print(
            f"minimum r/sqrt(N)            = "
            f"{min(r_ratios):.8f}"
        )

        print(
            f"maximum r/sqrt(N)            = "
            f"{max(r_ratios):.8f}"
        )

    print()
    print("-" * 60)
    print("FALLBACK RESCUE EXAMPLES")
    print("-" * 60)

    if rescue_examples:
        for (
            case_index,
            n,
            p,
            q,
            r,
            k,
            z,
            g,
            probes,
        ) in rescue_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"r={r} "
                f"k={k} "
                f"z={z} "
                f"gcd={g} "
                f"probes={probes}"
            )
    else:
        print("NONE")

    print()
    print("-" * 60)
    print("REMAINING FAILURES")
    print("-" * 60)

    if remaining_examples:
        for item in remaining_examples:
            (
                case_index,
                n,
                p,
                q,
                reason,
                r,
            ) = item

            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"reason={reason} "
                f"r={r}"
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
