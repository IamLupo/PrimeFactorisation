import math
import random
from sympy import isprime, prevprime


EXPERIMENT = 119

CASES = 300

# Larger range than Experiment 118.
P_MIN = 10_000
P_MAX = 100_000

Q_MIN = 100_000
Q_MAX = 250_000

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


def search_k(n, k):
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

        if 1 < g < n:
            return {
                "k": k,
                "z": z,
                "gcd": g,
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
    k_tests = 0
    z_probes = 0

    for k in powers_between(p, q):
        k_tests += 1

        result = search_k(
            n,
            k,
        )

        if result is not None:
            z_probes += result["probes"]

            return {
                "found": True,
                "factor": result["gcd"],
                "k": result["k"],
                "z": result["z"],
                "k_tests": k_tests,
                "z_probes": z_probes,
            }

        z_probes += MAX_Z_PROBES

    return {
        "found": False,
        "k_tests": k_tests,
        "z_probes": z_probes,
    }


def sqrt_prime_stage(n):
    sqrt_n = math.isqrt(n)
    r = int(prevprime(sqrt_n))

    return r


def sqrt_prime_search(n, p, q):
    r = sqrt_prime_stage(n)

    if not (p < r < q):
        return {
            "usable": False,
            "r": r,
        }

    result = search_k(
        n,
        r,
    )

    if result is None:
        return {
            "usable": True,
            "found": False,
            "r": r,
            "z_probes": MAX_Z_PROBES,
        }

    return {
        "usable": True,
        "found": True,
        "r": r,
        "factor": result["gcd"],
        "k": result["k"],
        "z": result["z"],
        "z_probes": result["probes"],
    }


def fermat_factor(n):
    a = math.isqrt(n)

    if a * a < n:
        a += 1

    iterations = 0

    while True:
        iterations += 1

        b2 = a * a - n
        b = math.isqrt(b2)

        if b * b == b2:
            f1 = a - b
            f2 = a + b

            if (
                f1 > 1
                and f2 > 1
                and f1 * f2 == n
            ):
                return {
                    "found": True,
                    "factor1": f1,
                    "factor2": f2,
                    "iterations": iterations,
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
            max(p + 1, Q_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        if p == q:
            continue

        return p, q


def verify_factor(n, factor):
    return (
        1 < factor < n
        and n % factor == 0
    )


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX}")
    print(f"MAX Z PROBES: {MAX_Z_PROBES}")
    print()

    random.seed(EXPERIMENT)

    dyadic_hits = 0
    sqrt_prime_hits = 0
    fermat_hits = 0
    total_recovered = 0
    failures = 0

    total_dyadic_k_tests = 0
    total_dyadic_z_probes = 0
    total_fallback_probes = 0

    fermat_iterations = []

    largest_n = 0
    largest_n_case = None

    failure_examples = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        if n > largest_n:
            largest_n = n
            largest_n_case = (
                case_index,
                n,
                p,
                q,
            )

        recovered = False

        # -------------------------------------------------
        # STAGE 1
        # -------------------------------------------------

        stage1 = dyadic_stage(
            n,
            p,
            q,
        )

        total_dyadic_k_tests += stage1[
            "k_tests"
        ]

        total_dyadic_z_probes += stage1[
            "z_probes"
        ]

        if stage1["found"]:
            factor = stage1["factor"]

            if verify_factor(
                n,
                factor,
            ):
                dyadic_hits += 1
                total_recovered += 1
                recovered = True

        if recovered:
            continue

        # -------------------------------------------------
        # STAGE 2
        # -------------------------------------------------

        stage2 = sqrt_prime_search(
            n,
            p,
            q,
        )

        if stage2.get("usable", False):
            total_fallback_probes += stage2.get(
                "z_probes",
                0,
            )

        if stage2.get("found", False):
            factor = stage2["factor"]

            if verify_factor(
                n,
                factor,
            ):
                sqrt_prime_hits += 1
                total_recovered += 1
                recovered = True

        if recovered:
            continue

        # -------------------------------------------------
        # STAGE 3
        # -------------------------------------------------

        fermat = fermat_factor(n)

        if fermat["found"]:
            f1 = fermat["factor1"]
            f2 = fermat["factor2"]

            if (
                f1 * f2 == n
                and verify_factor(n, f1)
            ):
                fermat_hits += 1
                total_recovered += 1
                recovered = True

                fermat_iterations.append(
                    fermat["iterations"]
                )

        if not recovered:
            failures += 1

            if len(failure_examples) < 20:
                failure_examples.append(
                    (
                        case_index,
                        n,
                        p,
                        q,
                    )
                )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = "
        f"{CASES}"
    )

    print(
        f"dyadic hits                  = "
        f"{dyadic_hits}"
    )

    print(
        f"sqrt-prime hits              = "
        f"{sqrt_prime_hits}"
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
        f"failures                     = "
        f"{failures}"
    )

    print()

    print(
        f"final recovery rate          = "
        f"{100.0 * total_recovered / CASES:.2f}%"
    )

    print()

    print(
        f"total dyadic k tests         = "
        f"{total_dyadic_k_tests}"
    )

    print(
        f"total dyadic z probes        = "
        f"{total_dyadic_z_probes}"
    )

    print(
        f"total sqrt-prime z probes    = "
        f"{total_fallback_probes}"
    )

    print()

    print(
        f"largest N                    = "
        f"{largest_n}"
    )

    if largest_n_case is not None:
        (
            case_index,
            n,
            p,
            q,
        ) = largest_n_case

        print(
            f"largest-N case              = "
            f"{case_index}"
        )

        print(
            f"largest-N factors           = "
            f"{p} * {q}"
        )

    if fermat_iterations:
        print()

        print(
            f"minimum Fermat iterations   = "
            f"{min(fermat_iterations)}"
        )

        print(
            f"maximum Fermat iterations   = "
            f"{max(fermat_iterations)}"
        )

        print(
            f"average Fermat iterations   = "
            f"{sum(fermat_iterations) / len(fermat_iterations):.2f}"
        )

    print()
    print("-" * 60)
    print("FAILURES")
    print("-" * 60)

    if failure_examples:
        for (
            case_index,
            n,
            p,
            q,
        ) in failure_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q}"
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
