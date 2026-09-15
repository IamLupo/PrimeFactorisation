import math
import random
from sympy import isprime, prevprime


EXPERIMENT = 118

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z_PROBES = 16


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

        v_k = -C(k-2,d)

    when d < k-1.
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


def search_k(n, p, q, k):
    """
    Search one k with a fixed maximum number of z probes.
    """
    max_z = k // 2

    if max_z <= 0:
        return None

    probes = van_der_corput_sequence(
        max_z,
        MAX_Z_PROBES,
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


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def dyadic_stage(n, p, q):
    """
    Try every power of two satisfying p < k < q.
    """
    evaluations = 0

    for k in powers_between(p, q):
        result = search_k(
            n,
            p,
            q,
            k,
        )

        evaluations += 1

        if result is not None:
            result["stage"] = "dyadic"
            result["k_evaluations"] = evaluations
            return result

    return {
        "stage": "dyadic_miss",
        "k_evaluations": evaluations,
    }


def previous_prime_sqrt_stage(n, p, q):
    """
    r = previous_prime(floor(sqrt(N)))
    k = r
    """
    sqrt_n = math.isqrt(n)
    r = int(prevprime(sqrt_n))

    if r <= p:
        return {
            "stage": "sqrt_prime_unusable",
            "r": r,
        }

    if r >= q:
        return {
            "stage": "sqrt_prime_unusable",
            "r": r,
        }

    result = search_k(
        n,
        p,
        q,
        r,
    )

    if result is None:
        return {
            "stage": "sqrt_prime_miss",
            "r": r,
        }

    result["stage"] = "sqrt_prime"
    result["r"] = r

    return result


def fermat_factor(n, p, q):
    """
    Fermat factorization:

        N = a^2 - b^2
          = (a-b)(a+b)

    Start at ceil(sqrt(N)).
    """
    a = math.isqrt(n)

    if a * a < n:
        a += 1

    iterations = 0

    while True:
        iterations += 1

        b2 = a * a - n
        b = math.isqrt(b2)

        if b * b == b2:
            factor1 = a - b
            factor2 = a + b

            if (
                factor1 > 1
                and factor2 > 1
                and factor1 * factor2 == n
            ):
                return {
                    "found": True,
                    "p": factor1,
                    "q": factor2,
                    "iterations": iterations,
                    "a": a,
                    "b": b,
                }

        a += 1


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


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q MAX: {Q_MAX}")
    print(f"MAX Z PROBES: {MAX_Z_PROBES}")
    print()

    random.seed(EXPERIMENT)

    total = 0

    dyadic_hits = 0
    sqrt_prime_hits = 0
    sqrt_prime_misses = 0

    fermat_attempts = 0
    fermat_hits = 0

    total_recovered = 0

    fermat_iterations = []

    close_factor_cases = 0
    sqrt_prime_unusable = 0

    rescue_examples = []
    fermat_examples = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        total += 1

        # -------------------------------------------------
        # STAGE 1
        # -------------------------------------------------

        stage1 = dyadic_stage(
            n,
            p,
            q,
        )

        if stage1["stage"] == "dyadic":
            dyadic_hits += 1
            total_recovered += 1
            continue

        # -------------------------------------------------
        # STAGE 2
        # -------------------------------------------------

        stage2 = previous_prime_sqrt_stage(
            n,
            p,
            q,
        )

        if stage2["stage"] == "sqrt_prime":
            sqrt_prime_hits += 1
            total_recovered += 1

            if len(rescue_examples) < 20:
                rescue_examples.append(
                    (
                        case_index,
                        n,
                        p,
                        q,
                        stage2["r"],
                        stage2["k"],
                        stage2["z"],
                        stage2["gcd"],
                        stage2["probes"],
                    )
                )

            continue

        if stage2["stage"] == "sqrt_prime_unusable":
            sqrt_prime_unusable += 1

        else:
            sqrt_prime_misses += 1

        # -------------------------------------------------
        # STAGE 3: FERMAT
        # -------------------------------------------------

        fermat_attempts += 1

        gap = q - p

        if gap <= 100:
            close_factor_cases += 1

        result = fermat_factor(
            n,
            p,
            q,
        )

        if result["found"]:
            fermat_hits += 1
            total_recovered += 1

            fermat_iterations.append(
                result["iterations"]
            )

            if len(fermat_examples) < 20:
                fermat_examples.append(
                    (
                        case_index,
                        n,
                        p,
                        q,
                        q - p,
                        result["iterations"],
                    )
                )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = {total}"
    )

    print(
        f"dyadic hits                  = {dyadic_hits}"
    )

    print(
        f"previous-prime-sqrt hits     = "
        f"{sqrt_prime_hits}"
    )

    print(
        f"previous-prime-sqrt misses   = "
        f"{sqrt_prime_misses}"
    )

    print(
        f"sqrt-prime unusable          = "
        f"{sqrt_prime_unusable}"
    )

    print(
        f"Fermat attempts              = "
        f"{fermat_attempts}"
    )

    print(
        f"Fermat hits                  = "
        f"{fermat_hits}"
    )

    print(
        f"total recovered              = "
        f"{total_recovered}"
    )

    print(
        f"close-factor cases (q-p<=100) = "
        f"{close_factor_cases}"
    )

    print()

    print(
        f"dyadic recovery rate          = "
        f"{100.0 * dyadic_hits / total:.2f}%"
    )

    print(
        f"sqrt-prime rescue rate        = "
        f"{100.0 * sqrt_prime_hits / max(1, total - dyadic_hits):.2f}%"
    )

    print(
        f"final recovery rate           = "
        f"{100.0 * total_recovered / total:.2f}%"
    )

    if fermat_iterations:
        print()

        print(
            f"minimum Fermat iterations     = "
            f"{min(fermat_iterations)}"
        )

        print(
            f"maximum Fermat iterations     = "
            f"{max(fermat_iterations)}"
        )

        print(
            f"average Fermat iterations     = "
            f"{sum(fermat_iterations) / len(fermat_iterations):.2f}"
        )

    print()
    print("-" * 60)
    print("SQRT-PRIME RESCUES")
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
    print("FERMAT CASES")
    print("-" * 60)

    if fermat_examples:
        for (
            case_index,
            n,
            p,
            q,
            gap,
            iterations,
        ) in fermat_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"gap={gap} "
                f"iterations={iterations}"
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
