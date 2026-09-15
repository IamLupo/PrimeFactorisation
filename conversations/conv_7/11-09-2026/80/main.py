import math
import random

from sympy import isprime


EXPERIMENT = 100


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


def classify(g, p, q, n):
    if g == 1:
        return "1"

    if g == p:
        return "p"

    if g == q:
        return "q"

    if g == n:
        return "n"

    return "other"


def exponential_bracket(n, p, q, k, max_z):
    """
    Search:

        1, 2, 4, 8, ...

    until we see something other than gcd=1.
    """

    evaluations = 0
    z = 1
    previous_z = 0

    while z <= max_z:
        g = gcd_value(n, z, k)
        evaluations += 1

        if g != 1:
            return {
                "low": previous_z,
                "high": z,
                "g": g,
                "evaluations": evaluations
            }

        previous_z = z
        z *= 2

    return {
        "low": None,
        "high": None,
        "g": 1,
        "evaluations": evaluations
    }


def three_state_search(n, p, q, k, low, high):
    """
    Search for the p-region.

    State meanings:

        gcd = 1  -> z too small
        gcd = p  -> success
        gcd = n  -> z too large

    We maintain a bracket containing the transition.
    """

    evaluations = 0

    while low <= high:
        mid = (low + high) // 2

        g = gcd_value(
            n,
            mid,
            k
        )

        evaluations += 1

        state = classify(
            g,
            p,
            q,
            n
        )

        if state == "p":
            return {
                "found": True,
                "z": mid,
                "g": g,
                "evaluations": evaluations
            }

        if state == "1":
            low = mid + 1
            continue

        if state == "n":
            high = mid - 1
            continue

        # Unexpected gcd.
        return {
            "found": False,
            "z": mid,
            "g": g,
            "evaluations": evaluations
        }

    return {
        "found": False,
        "z": None,
        "g": 1,
        "evaluations": evaluations
    }


def linear_reference(n, p, q, k, max_z):
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
                "g": g,
                "evaluations": evaluations,
                "class": classify(
                    g,
                    p,
                    q,
                    n
                )
            }

    return {
        "z": None,
        "g": 1,
        "evaluations": evaluations,
        "class": "1"
    }


def analyze_case(n, p, q, k, max_z):
    reference = linear_reference(
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
        return {
            "reference": reference,
            "bracket": bracket,
            "search": None
        }

    search = three_state_search(
        n,
        p,
        q,
        k,
        bracket["low"] + 1,
        bracket["high"]
    )

    return {
        "reference": reference,
        "bracket": bracket,
        "search": search
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
    bracketed = 0
    found_p = 0
    found_n_only = 0
    wrong_state = 0

    exact_matches = 0

    total_linear = 0
    total_exponential = 0
    total_bisection = 0

    examples = 0

    for _ in range(cases):
        n, p, q = generate_semiprime()

        for k in power_of_two_values(p, q):
            tested += 1

            result = analyze_case(
                n,
                p,
                q,
                k,
                max_z
            )

            reference = result["reference"]
            bracket = result["bracket"]
            search = result["search"]

            total_linear += (
                reference["evaluations"]
            )

            total_exponential += (
                bracket["evaluations"]
            )

            if search is not None:
                bracketed += 1
                total_bisection += (
                    search["evaluations"]
                )

            if (
                reference["class"] == "p"
            ):
                if (
                    search is not None
                    and search["found"]
                    and search["g"] == p
                ):
                    found_p += 1

                    if (
                        search["z"]
                        == reference["z"]
                    ):
                        exact_matches += 1

            elif (
                reference["class"] == "n"
            ):
                found_n_only += 1

            else:
                wrong_state += 1

            if examples < 25:
                examples += 1

                print(
                    f"CASE {examples}"
                )
                print(f"  N={n}")
                print(f"  p={p}")
                print(f"  q={q}")
                print(f"  k={k}")

                print(
                    "  reference: "
                    f"z={reference['z']} "
                    f"gcd={reference['g']} "
                    f"class={reference['class']}"
                )

                print(
                    "  exponential: "
                    f"[{bracket['low']}, "
                    f"{bracket['high']}] "
                    f"endpoint_gcd={bracket['g']}"
                )

                if search is None:
                    print(
                        "  bisection: NONE"
                    )
                else:
                    print(
                        "  bisection: "
                        f"z={search['z']} "
                        f"gcd={search['g']} "
                        f"found={search['found']} "
                        f"evals={search['evaluations']}"
                    )

                print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(
        f"POWER-OF-TWO CASES: {tested}"
    )

    print(
        f"BRACKETED CASES: {bracketed}"
    )

    print(
        f"REFERENCE FIRST HIT = p: "
        f"{found_p}"
    )

    print(
        f"REFERENCE FIRST HIT = n: "
        f"{found_n_only}"
    )

    print(
        f"OTHER / UNEXPECTED: "
        f"{wrong_state}"
    )

    print(
        f"EXACT z MATCHES: "
        f"{exact_matches}"
    )

    if tested:
        print()

        print(
            "AVERAGE LINEAR EVALUATIONS: "
            f"{total_linear / tested:.4f}"
        )

        print(
            "AVERAGE EXPONENTIAL EVALUATIONS: "
            f"{total_exponential / tested:.4f}"
        )

        print(
            "AVERAGE BISECTION EVALUATIONS: "
            f"{total_bisection / max(bracketed, 1):.4f}"
        )

        total_optimized = (
            total_exponential
            + total_bisection
        )

        print(
            "AVERAGE TOTAL OPTIMIZED: "
            f"{total_optimized / tested:.4f}"
        )

        if total_optimized:
            print(
                "LINEAR / OPTIMIZED SPEEDUP: "
                f"{total_linear / total_optimized:.4f}x"
            )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    random.seed(100011)

    run_experiment(
        cases=500,
        max_z=10000
    )


if __name__ == "__main__":
    main()
