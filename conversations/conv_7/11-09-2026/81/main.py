import math
import random

from sympy import isprime


EXPERIMENT = 101


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


def exponential_search(n, p, q, k, max_z):
    """
    STRICTLY tests:

        z = 1, 2, 4, 8, 16, ...

    Never scans z linearly.
    """

    z = 1
    previous_z = 0
    evaluations = 0

    while z <= max_z:
        g = gcd_value(
            n,
            z,
            k
        )

        evaluations += 1

        state = classify(
            g,
            p,
            q,
            n
        )

        if state != "1":
            return {
                "found": True,
                "state": state,
                "z": z,
                "gcd": g,
                "previous_z": previous_z,
                "evaluations": evaluations,
            }

        previous_z = z
        z *= 2

    return {
        "found": False,
        "state": "1",
        "z": None,
        "gcd": 1,
        "previous_z": previous_z,
        "evaluations": evaluations,
    }


def refine_interval(n, p, q, k, low, high):
    """
    Refine [low, high] without a linear scan.

    We always test the midpoint.

    State:
        1     = z too small
        p/q   = factor found
        n     = z too large

    The search stops immediately when any proper
    factor is found.

    This is deliberately NOT a search for the first z.
    We only want a usable factor.
    """

    evaluations = 0

    intervals = [
        (low, high)
    ]

    visited = set()

    while intervals:
        left, right = intervals.pop()

        if left > right:
            continue

        mid = (
            left + right
        ) // 2

        if mid in visited:
            continue

        visited.add(mid)

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

        if state in (
            "p",
            "q",
            "other"
        ):
            return {
                "found": True,
                "state": state,
                "z": mid,
                "gcd": g,
                "evaluations": evaluations,
            }

        if state == "1":
            intervals.append(
                (
                    mid + 1,
                    right
                )
            )

        elif state == "n":
            intervals.append(
                (
                    left,
                    mid - 1
                )
            )

            # Also keep the lower side available because
            # the p-window can have already been crossed.
            if (
                mid - 1 >= left
                and mid + 1 <= right
            ):
                intervals.append(
                    (
                        mid + 1,
                        right
                    )
                )

    return {
        "found": False,
        "state": "1",
        "z": None,
        "gcd": 1,
        "evaluations": evaluations,
    }


def search(n, p, q, k, max_z):
    exponential = exponential_search(
        n,
        p,
        q,
        k,
        max_z
    )

    if not exponential["found"]:
        return {
            "found": False,
            "mode": "exponential_exhausted",
            "z": None,
            "gcd": 1,
            "state": "1",
            "evaluations":
                exponential["evaluations"],
        }

    if exponential["state"] in (
        "p",
        "q",
        "other"
    ):
        return {
            "found": True,
            "mode": "direct_power_hit",
            "z": exponential["z"],
            "gcd": exponential["gcd"],
            "state": exponential["state"],
            "evaluations":
                exponential["evaluations"],
        }

    # The only remaining nontrivial case is n.
    #
    # We know:
    #
    # previous_z -> gcd 1
    # current_z  -> gcd n
    #
    # Search the interval using midpoint subdivision.

    refined = refine_interval(
        n,
        p,
        q,
        k,
        exponential["previous_z"] + 1,
        exponential["z"] - 1
    )

    return {
        "found": refined["found"],
        "mode": "refined_after_n",
        "z": refined["z"],
        "gcd": refined["gcd"],
        "state": refined["state"],
        "evaluations":
            exponential["evaluations"]
            + refined["evaluations"],
    }


def run_experiment(
    cases=500,
    max_z=16384
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
    found = 0
    found_p = 0
    found_q = 0
    found_n = 0
    failed = 0

    total_evaluations = 0

    examples_printed = 0

    for _ in range(cases):
        n, p, q = generate_semiprime()

        for k in power_of_two_values(p, q):
            tested += 1

            result = search(
                n,
                p,
                q,
                k,
                max_z
            )

            total_evaluations += (
                result["evaluations"]
            )

            if not result["found"]:
                failed += 1
                continue

            found += 1

            if result["state"] == "p":
                found_p += 1

            elif result["state"] == "q":
                found_q += 1

            elif result["state"] == "n":
                found_n += 1

            if examples_printed < 25:
                examples_printed += 1

                print(
                    f"CASE {examples_printed}"
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
                    f"  result_z={result['z']}"
                )

                print(
                    f"  gcd={result['gcd']}"
                )

                print(
                    f"  state={result['state']}"
                )

                print(
                    f"  mode={result['mode']}"
                )

                print(
                    f"  gcd evaluations="
                    f"{result['evaluations']}"
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
        f"FOUND: "
        f"{found}"
    )

    print(
        f"FOUND p: "
        f"{found_p}"
    )

    print(
        f"FOUND q: "
        f"{found_q}"
    )

    print(
        f"FOUND n: "
        f"{found_n}"
    )

    print(
        f"FAILED: "
        f"{failed}"
    )

    if tested:
        print(
            f"AVERAGE gcd EVALUATIONS: "
            f"{total_evaluations / tested:.6f}"
        )

    print()
    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    random.seed(101001)

    run_experiment(
        cases=500,
        max_z=16384
    )


if __name__ == "__main__":
    main()