import math
import random
from sympy import isprime


EXPERIMENT = 104

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384

# Maximum number of gcd evaluations allowed for refinement.
REFINEMENT_BUDGET = 128


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
    For x=n-1, y=1:

        d = 2z - 1

    and, when d < k-1,

        v_k = -C(k-2, d).

    Therefore only C(k-2,d) is required.
    """

    degree = 2 * z - 1

    if degree >= k - 1:
        return n

    value = math.comb(k - 2, degree)
    value %= n

    value = (-value) % n

    return math.gcd(value, n)


def dyadic_search(n, p, q, k):
    """
    Strictly evaluates:

        z = 1, 2, 4, 8, ...

    Returns either a direct factor hit or the interval
    between the last gcd=1 power and the first gcd=n power.
    """

    evaluations = 0

    z = 1
    previous_z = 0

    while z <= MAX_Z:
        evaluations += 1

        g = value_gcd(n, k, z)
        state = gcd_state(g, n, p, q)

        if state == "p" or state == "q":
            return {
                "mode": "direct",
                "z": z,
                "gcd": g,
                "evaluations": evaluations,
            }

        if state == "n":
            return {
                "mode": "bracket",
                "low": previous_z + 1,
                "high": z - 1,
                "evaluations": evaluations,
            }

        previous_z = z
        z *= 2

    return {
        "mode": "limit",
        "evaluations": evaluations,
    }


def recursive_refine(n, p, q, k, low, high):
    """
    Recursively subdivide the whole interval.

    IMPORTANT:

    This is deliberately NOT ordinary binary search.

    If a midpoint has gcd=1, BOTH sides remain possible because
    the p-window can disappear again before the upper boundary.

    The search therefore explores:

        left half
        right half

    until a proper factor is found or the evaluation budget is hit.

    No linear z loop is used.
    """

    evaluations = 0

    stack = [(low, high)]

    while stack and evaluations < REFINEMENT_BUDGET:
        left, right = stack.pop()

        if left > right:
            continue

        mid = (left + right) // 2

        evaluations += 1

        g = value_gcd(n, k, mid)
        state = gcd_state(g, n, p, q)

        if state == "p" or state == "q":
            return {
                "found": True,
                "z": mid,
                "gcd": g,
                "state": state,
                "evaluations": evaluations,
            }

        if left == right:
            continue

        if state == "1":
            left_high = mid - 1
            right_low = mid + 1

            # Explore both halves.
            #
            # Push the right side first so the left side is
            # processed next.

            if right_low <= right:
                stack.append((right_low, right))

            if left <= left_high:
                stack.append((left, left_high))

            continue

        if state == "n":
            # n is normally reached at the trivial degree boundary.
            #
            # There can still be a useful region below it,
            # so continue searching the lower half.

            if left <= mid - 1:
                stack.append((left, mid - 1))

            continue

        # Unknown gcd state.
        #
        # Conservatively continue both sides.

        if right_low := mid + 1:
            if right_low <= right:
                stack.append((right_low, right))

        if left <= mid - 1:
            stack.append((left, mid - 1))

    return {
        "found": False,
        "z": None,
        "gcd": None,
        "state": None,
        "evaluations": evaluations,
    }


def generate_semiprime():
    while True:
        p = random.randrange(P_MIN, P_MAX + 1)

        if not isprime(p):
            continue

        q = random.randrange(max(p + 1, P_MIN), Q_MAX + 1)

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
    print(f"REFINEMENT BUDGET: {REFINEMENT_BUDGET}")
    print()

    random.seed(104)

    total_k_cases = 0

    direct_hits = 0
    bracket_cases = 0

    recursive_hits = 0
    recursive_failures = 0

    p_hits = 0
    q_hits = 0

    total_evaluations = 0
    maximum_evaluations = 0

    examples = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            dyadic = dyadic_search(
                n,
                p,
                q,
                k,
            )

            if dyadic["mode"] == "direct":
                direct_hits += 1

                evaluations = dyadic["evaluations"]

                total_evaluations += evaluations
                maximum_evaluations = max(
                    maximum_evaluations,
                    evaluations,
                )

                state = gcd_state(
                    dyadic["gcd"],
                    n,
                    p,
                    q,
                )

                if state == "p":
                    p_hits += 1

                elif state == "q":
                    q_hits += 1

                if len(examples) < 10:
                    examples.append(
                        (
                            "DIRECT",
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            dyadic["z"],
                            dyadic["gcd"],
                            evaluations,
                        )
                    )

                continue

            if dyadic["mode"] != "bracket":
                total_evaluations += dyadic["evaluations"]

                maximum_evaluations = max(
                    maximum_evaluations,
                    dyadic["evaluations"],
                )

                continue

            bracket_cases += 1

            low = dyadic["low"]
            high = dyadic["high"]

            refined = recursive_refine(
                n,
                p,
                q,
                k,
                low,
                high,
            )

            evaluations = (
                dyadic["evaluations"]
                + refined["evaluations"]
            )

            total_evaluations += evaluations

            maximum_evaluations = max(
                maximum_evaluations,
                evaluations,
            )

            if refined["found"]:
                recursive_hits += 1

                if refined["state"] == "p":
                    p_hits += 1

                elif refined["state"] == "q":
                    q_hits += 1

                if len(examples) < 20:
                    examples.append(
                        (
                            "RECURSIVE",
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            low,
                            high,
                            refined["z"],
                            refined["gcd"],
                            evaluations,
                        )
                    )

            else:
                recursive_failures += 1

                if len(examples) < 20:
                    examples.append(
                        (
                            "FAIL",
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            low,
                            high,
                            evaluations,
                        )
                    )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(f"semiprime cases              = {CASES}")
    print(f"(N,k) test cases             = {total_k_cases}")
    print(f"direct dyadic hits           = {direct_hits}")
    print(f"dyadic brackets              = {bracket_cases}")
    print(f"recursive factor hits        = {recursive_hits}")
    print(f"recursive failures           = {recursive_failures}")
    print(f"p hits                       = {p_hits}")
    print(f"q hits                       = {q_hits}")
    print(f"total gcd evaluations        = {total_evaluations}")
    print(f"maximum gcd evaluations      = {maximum_evaluations}")

    recovered = direct_hits + recursive_hits

    if total_k_cases:
        print(
            f"direct hit rate              = "
            f"{100.0 * direct_hits / total_k_cases:.2f}%"
        )

        print(
            f"bracket rate                 = "
            f"{100.0 * bracket_cases / total_k_cases:.2f}%"
        )

        print(
            f"recursive recovery rate     = "
            f"{100.0 * recursive_hits / bracket_cases:.2f}%"
            if bracket_cases
            else "recursive recovery rate     = N/A"
        )

        print(
            f"total recovered rate         = "
            f"{100.0 * recovered / total_k_cases:.2f}%"
        )

    print()
    print("-" * 60)
    print("EXAMPLES")
    print("-" * 60)

    for example in examples:
        if example[0] == "DIRECT":
            (
                _,
                case_index,
                n,
                p,
                q,
                k,
                z,
                g,
                evaluations,
            ) = example

            print(
                f"DIRECT     "
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"z={z} "
                f"gcd={g} "
                f"evals={evaluations}"
            )

        elif example[0] == "RECURSIVE":
            (
                _,
                case_index,
                n,
                p,
                q,
                k,
                low,
                high,
                z,
                g,
                evaluations,
            ) = example

            print(
                f"RECURSIVE  "
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"bracket=[{low},{high}] "
                f"z={z} "
                f"gcd={g} "
                f"evals={evaluations}"
            )

        else:
            (
                _,
                case_index,
                n,
                p,
                q,
                k,
                low,
                high,
                evaluations,
            ) = example

            print(
                f"FAIL       "
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"bracket=[{low},{high}] "
                f"evals={evaluations}"
            )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
