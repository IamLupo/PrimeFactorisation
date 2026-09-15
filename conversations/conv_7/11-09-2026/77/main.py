import math
import random


EXPERIMENT = 97


def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    limit = math.isqrt(n)

    d = 3

    while d <= limit:
        if n % d == 0:
            return True

        d += 2

    return True


def generate_semiprime(
    min_prime=100,
    max_prime=10000
):
    while True:
        p = random.randrange(
            min_prime,
            max_prime + 1
        )

        q = random.randrange(
            min_prime,
            max_prime + 1
        )

        if p >= q:
            continue

        if not is_prime(p):
            continue

        if not is_prime(q):
            continue

        return p * q, p, q


def gcd_v(n, z, k):
    degree = 2 * z - 1

    if k == 1:
        value = 1 % n
        return math.gcd(value, n)

    if k - 2 < degree:
        return n

    value = (
        -math.comb(
            k - 2,
            degree
        )
    ) % n

    return math.gcd(
        value,
        n
    )


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


def linear_search(
    n,
    p,
    q,
    k,
    max_z
):
    evaluations = 0

    for z in range(
        1,
        max_z + 1
    ):
        g = gcd_v(
            n,
            z,
            k
        )

        evaluations += 1

        if g != 1:
            return {
                "z": z,
                "g": g,
                "evaluations": evaluations
            }

    return {
        "z": None,
        "g": 1,
        "evaluations": evaluations
    }


def exponential_search(
    n,
    p,
    q,
    k,
    max_z
):
    evaluations = 0

    z = 1
    previous = 0

    while z <= max_z:
        g = gcd_v(
            n,
            z,
            k
        )

        evaluations += 1

        if g != 1:
            return {
                "low": previous + 1,
                "high": z,
                "g": g,
                "evaluations": evaluations
            }

        previous = z
        z *= 2

    return {
        "low": None,
        "high": None,
        "g": 1,
        "evaluations": evaluations
    }


def binary_search_first_hit(
    n,
    p,
    q,
    k,
    low,
    high
):
    evaluations = 0

    answer_z = high
    answer_g = gcd_v(
        n,
        high,
        k
    )

    evaluations += 1

    while low < high:
        mid = (
            low + high
        ) // 2

        g = gcd_v(
            n,
            mid,
            k
        )

        evaluations += 1

        if g != 1:
            answer_z = mid
            answer_g = g
            high = mid
        else:
            low = mid + 1

    return {
        "z": answer_z,
        "g": answer_g,
        "evaluations": evaluations
    }


def powers_between(p, q):
    result = []

    k = 1

    while k < q:
        if p < k < q:
            result.append(k)

        k *= 2

    return result


def analyze_case(
    n,
    p,
    q,
    max_z
):
    results = []

    for k in powers_between(
        p,
        q
    ):
        linear = linear_search(
            n,
            p,
            q,
            k,
            max_z
        )

        exponential = exponential_search(
            n,
            p,
            q,
            k,
            max_z
        )

        binary = None

        if (
            exponential["low"]
            is not None
        ):
            binary = binary_search_first_hit(
                n,
                p,
                q,
                k,
                exponential["low"],
                exponential["high"]
            )

        results.append(
            {
                "k": k,
                "linear": linear,
                "exponential": exponential,
                "binary": binary
            }
        )

    return results


def run_experiment(
    cases=500,
    max_z=10000
):
    print("=" * 60)
    print(
        f"START EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)
    print()

    print(
        f"CASES: {cases}"
    )

    print(
        f"MAX Z: {max_z}"
    )

    print()

    tested = 0
    linear_hits = 0
    binary_hits = 0
    binary_matches = 0
    factor_p_hits = 0
    factor_n_hits = 0

    total_linear_evaluations = 0
    total_exponential_evaluations = 0
    total_binary_evaluations = 0

    examples = 0

    for case in range(
        1,
        cases + 1
    ):
        n, p, q = (
            generate_semiprime()
        )

        results = analyze_case(
            n,
            p,
            q,
            max_z
        )

        for result in results:
            tested += 1

            k = result["k"]
            linear = result["linear"]
            exponential = result["exponential"]
            binary = result["binary"]

            if linear["z"] is not None:
                linear_hits += 1

                total_linear_evaluations += (
                    linear["evaluations"]
                )

                linear_class = classify_gcd(
                    linear["g"],
                    p,
                    q,
                    n
                )

                if linear_class == "p":
                    factor_p_hits += 1

                elif linear_class == "n":
                    factor_n_hits += 1

            total_exponential_evaluations += (
                exponential["evaluations"]
            )

            if binary is not None:
                binary_hits += 1

                total_binary_evaluations += (
                    binary["evaluations"]
                )

                if (
                    linear["z"]
                    == binary["z"]
                    and linear["g"]
                    == binary["g"]
                ):
                    binary_matches += 1

            if examples < 20:
                examples += 1

                print(
                    f"CASE {examples}"
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
                    f"  k={k}"
                )

                print(
                    f"  linear z="
                    f"{linear['z']} "
                    f"gcd="
                    f"{linear['g']} "
                    f"class="
                    f"{classify_gcd(linear['g'],p,q,n)}"
                )

                print(
                    f"  exponential range="
                    f"[{exponential['low']},"
                    f"{exponential['high']}]"
                )

                if binary is None:
                    print(
                        "  binary: NONE"
                    )
                else:
                    print(
                        f"  binary z="
                        f"{binary['z']} "
                        f"gcd="
                        f"{binary['g']}"
                    )

                print(
                    f"  evaluations: "
                    f"linear={linear['evaluations']} "
                    f"exp={exponential['evaluations']} "
                    f"binary="
                    f"{0 if binary is None else binary['evaluations']}"
                )

                print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(
        f"POWER-OF-TWO CASES: "
        f"{tested}"
    )

    print(
        f"LINEAR HITS: "
        f"{linear_hits}"
    )

    print(
        f"EXPONENTIAL SEARCH + RANGE FOUND: "
        f"{binary_hits}"
    )

    print(
        f"EXACT LINEAR/BINARY MATCHES: "
        f"{binary_matches}"
    )

    print(
        f"FIRST HIT = p: "
        f"{factor_p_hits}"
    )

    print(
        f"FIRST HIT = n: "
        f"{factor_n_hits}"
    )

    if tested:
        print()

        print(
            f"AVERAGE LINEAR EVALUATIONS: "
            f"{total_linear_evaluations / tested:.4f}"
        )

        print(
            f"AVERAGE EXPONENTIAL EVALUATIONS: "
            f"{total_exponential_evaluations / tested:.4f}"
        )

        print(
            f"AVERAGE BINARY EVALUATIONS: "
            f"{total_binary_evaluations / max(binary_hits,1):.4f}"
        )

    print()

    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    random.seed(970011)

    run_experiment(
        cases=500,
        max_z=10000
    )


if __name__ == "__main__":
    main()
