import math
import random

from sympy import isprime


EXPERIMENT = 98


def generate_semiprime(min_prime=100, max_prime=10000):
    while True:
        p = random.randint(min_prime, max_prime)
        q = random.randint(min_prime, max_prime)

        if p >= q:
            continue

        if not isprime(p):
            continue

        if not isprime(q):
            continue

        return p * q, p, q


def power_of_two_values(p, q):
    values = []
    k = 1

    while k < q:
        if p < k < q:
            values.append(k)

        k *= 2

    return values


def arthseq_value_mod_n(n, z, k):
    """
    x = n - 1
    y = 1

    Coefficients:

        [1, n-1, 1, n-1, ...]

    repeated z times.

    Degree:

        d = 2z - 1

    Newton-form evaluation:

        v_k =
            sum_j a_j * C(k-1, j)
    """

    degree = 2 * z - 1

    if k <= degree:
        return 0

    value = 0

    for j in range(degree + 1):
        if j % 2 == 0:
            coefficient = 1
        else:
            coefficient = n - 1

        value += (
            coefficient
            * math.comb(k - 1, j)
        )

        value %= n

    return value


def gcd_value(n, z, k):
    value = arthseq_value_mod_n(
        n,
        z,
        k
    )

    return math.gcd(value, n)


def classify_gcd(g, p, q, n):
    if g == 1:
        return "1"

    if g == p:
        return "p"

    if g == q:
        return "q"

    if g == n:
        return "n"

    return "other"


def linear_first_hit(n, p, q, k, max_z):
    evaluations = 0

    for z in range(1, max_z + 1):
        g = gcd_value(
            n,
            z,
            k
        )

        evaluations += 1

        if g != 1:
            return {
                "z": z,
                "gcd": g,
                "class": classify_gcd(
                    g,
                    p,
                    q,
                    n
                ),
                "evaluations": evaluations
            }

    return {
        "z": None,
        "gcd": 1,
        "class": "1",
        "evaluations": evaluations
    }


def exponential_bracket(n, p, q, k, max_z):
    evaluations = 0

    low = 0
    z = 1

    while z <= max_z:
        g = gcd_value(
            n,
            z,
            k
        )

        evaluations += 1

        if g != 1:
            return {
                "low": low + 1,
                "high": z,
                "gcd": g,
                "evaluations": evaluations
            }

        low = z
        z *= 2

    return {
        "low": None,
        "high": None,
        "gcd": 1,
        "evaluations": evaluations
    }


def bracket_scan(n, p, q, k, low, high):
    evaluations = 0

    for z in range(low, high + 1):
        g = gcd_value(
            n,
            z,
            k
        )

        evaluations += 1

        if g != 1:
            return {
                "z": z,
                "gcd": g,
                "class": classify_gcd(
                    g,
                    p,
                    q,
                    n
                ),
                "evaluations": evaluations
            }

    return {
        "z": None,
        "gcd": 1,
        "class": "1",
        "evaluations": evaluations
    }


def analyze_power_of_two(n, p, q, k, max_z):
    linear = linear_first_hit(
        n,
        p,
        q,
        k,
        max_z
    )

    bracket = exponential_bracket(
        n,
        p,
        q,
        k,
        max_z
    )

    if bracket["low"] is None:
        optimized = {
            "z": None,
            "gcd": 1,
            "class": "1",
            "evaluations": 0
        }
    else:
        optimized = bracket_scan(
            n,
            p,
            q,
            k,
            bracket["low"],
            bracket["high"]
        )

    total_optimized_evaluations = (
        bracket["evaluations"]
        + optimized["evaluations"]
    )

    return {
        "linear": linear,
        "bracket": bracket,
        "optimized": optimized,
        "optimized_evaluations":
            total_optimized_evaluations
    }


def run_experiment(cases=500, max_z=10000):
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {cases}")
    print(f"MAX Z: {max_z}")
    print()

    tested = 0
    successful = 0

    factor_p_hits = 0
    factor_q_hits = 0
    factor_n_hits = 0
    wrong_hits = 0

    exact_matches = 0

    total_linear_evaluations = 0
    total_exponential_evaluations = 0
    total_bracket_evaluations = 0
    total_optimized_evaluations = 0

    examples_printed = 0

    for _ in range(cases):
        n, p, q = generate_semiprime()

        powers = power_of_two_values(
            p,
            q
        )

        for k in powers:
            tested += 1

            result = analyze_power_of_two(
                n,
                p,
                q,
                k,
                max_z
            )

            linear = result["linear"]
            bracket = result["bracket"]
            optimized = result["optimized"]

            total_linear_evaluations += (
                linear["evaluations"]
            )

            total_exponential_evaluations += (
                bracket["evaluations"]
            )

            total_bracket_evaluations += (
                optimized["evaluations"]
            )

            total_optimized_evaluations += (
                result["optimized_evaluations"]
            )

            if linear["z"] is None:
                continue

            successful += 1

            if linear["class"] == "p":
                factor_p_hits += 1

            elif linear["class"] == "q":
                factor_q_hits += 1

            elif linear["class"] == "n":
                factor_n_hits += 1

            else:
                wrong_hits += 1

            if (
                optimized["z"] == linear["z"]
                and optimized["gcd"] == linear["gcd"]
            ):
                exact_matches += 1

            if examples_printed < 25:
                examples_printed += 1

                print(f"CASE {examples_printed}")
                print(f"  N={n}")
                print(f"  p={p}")
                print(f"  q={q}")
                print(f"  k={k}")

                print(
                    "  linear: "
                    f"z={linear['z']} "
                    f"gcd={linear['gcd']} "
                    f"class={linear['class']}"
                )

                print(
                    "  exponential bracket: "
                    f"[{bracket['low']}, "
                    f"{bracket['high']}]"
                )

                print(
                    "  optimized: "
                    f"z={optimized['z']} "
                    f"gcd={optimized['gcd']} "
                    f"class={optimized['class']}"
                )

                print(
                    "  evaluations: "
                    f"linear={linear['evaluations']} "
                    f"exp={bracket['evaluations']} "
                    f"bracket={optimized['evaluations']} "
                    f"total={result['optimized_evaluations']}"
                )

                print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(f"POWER-OF-TWO CASES: {tested}")
    print(f"SUCCESSFUL FIRST HITS: {successful}")

    print(
        f"FIRST HIT = p: "
        f"{factor_p_hits}"
    )

    print(
        f"FIRST HIT = q: "
        f"{factor_q_hits}"
    )

    print(
        f"FIRST HIT = n: "
        f"{factor_n_hits}"
    )

    print(
        f"OTHER / WRONG: "
        f"{wrong_hits}"
    )

    print(
        f"OPTIMIZED EXACT MATCHES: "
        f"{exact_matches}"
    )

    if tested > 0:
        avg_linear = (
            total_linear_evaluations / tested
        )

        avg_exponential = (
            total_exponential_evaluations / tested
        )

        avg_bracket = (
            total_bracket_evaluations / tested
        )

        avg_optimized = (
            total_optimized_evaluations / tested
        )

        speedup = (
            total_linear_evaluations
            / total_optimized_evaluations
            if total_optimized_evaluations > 0
            else 0.0
        )

        print()
        print(
            f"AVERAGE LINEAR EVALUATIONS: "
            f"{avg_linear:.6f}"
        )

        print(
            f"AVERAGE EXPONENTIAL EVALUATIONS: "
            f"{avg_exponential:.6f}"
        )

        print(
            f"AVERAGE BRACKET-SCAN EVALUATIONS: "
            f"{avg_bracket:.6f}"
        )

        print(
            f"AVERAGE OPTIMIZED EVALUATIONS: "
            f"{avg_optimized:.6f}"
        )

        print(
            f"LINEAR/OPTIMIZED SPEEDUP: "
            f"{speedup:.6f}x"
        )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    random.seed(980011)

    run_experiment(
        cases=500,
        max_z=10000
    )


if __name__ == "__main__":
    main()