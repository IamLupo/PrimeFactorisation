import math
import random
from sympy import isprime


EXPERIMENT = 103

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384


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
    For x = n-1, y = 1 and degree d = 2z-1:

        v_k = -C(k-2, d),     d < k-1
        v_k = 0,              d >= k-1

    Only the gcd with n is required.
    """
    degree = 2 * z - 1

    if degree >= k - 1:
        return n

    value = math.comb(k - 2, degree)
    value %= n
    value = (-value) % n

    return math.gcd(value, n)


def dyadic_bracket(n, p, q, k):
    """
    Evaluate only:

        z = 1, 2, 4, 8, ...

    Returns:

        ("factor", z, gcd, evaluations)
            if a factor is found directly.

        ("bracket", low, high, evaluations)
            if we get:

                gcd(low)  = 1
                gcd(high) = n
    """

    evaluations = 0

    z = 1
    previous_z = 0

    while z <= MAX_Z:
        evaluations += 1

        g = value_gcd(n, k, z)
        state = gcd_state(g, n, p, q)

        if state == "p" or state == "q":
            return (
                "factor",
                z,
                g,
                evaluations,
            )

        if state == "n":
            return (
                "bracket",
                previous_z + 1,
                z - 1,
                evaluations,
            )

        if state != "1":
            return (
                "other",
                z,
                g,
                evaluations,
            )

        previous_z = z
        z *= 2

    return (
        "limit",
        None,
        None,
        evaluations,
    )


def binary_refine(n, p, q, k, low, high):
    """
    Search the interval [low, high] by binary subdivision.

    No linear scan is performed.

    Expected structure:

        1 ... 1 | p | p | ... | n ... n

    If midpoint is:

        1 -> search right
        p/q -> success
        n -> search left
    """

    evaluations = 0

    best_z = None
    best_g = None
    best_state = None

    while low <= high:
        mid = (low + high) // 2

        evaluations += 1

        g = value_gcd(n, k, mid)
        state = gcd_state(g, n, p, q)

        if state == "p" or state == "q":
            best_z = mid
            best_g = g
            best_state = state

            return (
                True,
                best_z,
                best_g,
                best_state,
                evaluations,
            )

        if state == "1":
            low = mid + 1
            continue

        if state == "n":
            high = mid - 1
            continue

        return (
            False,
            mid,
            g,
            state,
            evaluations,
        )

    return (
        False,
        None,
        None,
        None,
        evaluations,
    )


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
    print()

    random.seed(103)

    total_k_cases = 0

    direct_factor_hits = 0
    dyadic_n_jumps = 0

    refinement_attempts = 0
    refinement_factor_hits = 0
    refinement_failures = 0

    p_hits = 0
    q_hits = 0

    total_evaluations = 0
    max_evaluations = 0

    examples = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            result = dyadic_bracket(n, p, q, k)

            mode = result[0]

            if mode == "factor":
                _, z, g, evaluations = result

                direct_factor_hits += 1
                total_evaluations += evaluations

                if evaluations > max_evaluations:
                    max_evaluations = evaluations

                state = gcd_state(g, n, p, q)

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
                            z,
                            g,
                            evaluations,
                        )
                    )

            elif mode == "bracket":
                _, low, high, dyadic_evaluations = result

                dyadic_n_jumps += 1
                refinement_attempts += 1

                (
                    found,
                    z,
                    g,
                    state,
                    refine_evaluations,
                ) = binary_refine(
                    n,
                    p,
                    q,
                    k,
                    low,
                    high,
                )

                evaluations = (
                    dyadic_evaluations +
                    refine_evaluations
                )

                total_evaluations += evaluations

                if evaluations > max_evaluations:
                    max_evaluations = evaluations

                if found:
                    refinement_factor_hits += 1

                    if state == "p":
                        p_hits += 1
                    elif state == "q":
                        q_hits += 1

                    if len(examples) < 20:
                        examples.append(
                            (
                                "REFINED",
                                case_index,
                                n,
                                p,
                                q,
                                k,
                                z,
                                g,
                                evaluations,
                                low,
                                high,
                            )
                        )
                else:
                    refinement_failures += 1

                    if len(examples) < 20:
                        examples.append(
                            (
                                "FAIL",
                                case_index,
                                n,
                                p,
                                q,
                                k,
                                z,
                                g,
                                evaluations,
                                low,
                                high,
                            )
                        )

            else:
                total_evaluations += result[-1]

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(f"semiprime cases              = {CASES}")
    print(f"(N,k) test cases             = {total_k_cases}")
    print(f"direct dyadic factor hits    = {direct_factor_hits}")
    print(f"dyadic 1 -> n brackets       = {dyadic_n_jumps}")
    print(f"binary refinements attempted = {refinement_attempts}")
    print(f"binary factor hits           = {refinement_factor_hits}")
    print(f"binary refinement failures   = {refinement_failures}")
    print(f"p hits                       = {p_hits}")
    print(f"q hits                       = {q_hits}")
    print(f"total gcd evaluations        = {total_evaluations}")
    print(f"maximum gcd evaluations      = {max_evaluations}")

    if total_k_cases:
        direct_rate = (
            100.0 * direct_factor_hits / total_k_cases
        )

        bracket_rate = (
            100.0 * dyadic_n_jumps / total_k_cases
        )

        recovered_rate = (
            100.0
            * (
                direct_factor_hits
                + refinement_factor_hits
            )
            / total_k_cases
        )

        print(f"direct factor rate           = {direct_rate:.2f}%")
        print(f"1->n bracket rate            = {bracket_rate:.2f}%")
        print(f"total recovered factor rate  = {recovered_rate:.2f}%")

    print()
    print("-" * 60)
    print("EXAMPLES")
    print("-" * 60)

    for example in examples:
        mode = example[0]

        if mode == "DIRECT":
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
                f"DIRECT  "
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"z={z} "
                f"gcd={g} "
                f"evals={evaluations}"
            )

        elif mode == "REFINED":
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
                low,
                high,
            ) = example

            print(
                f"REFINED "
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
                z,
                g,
                evaluations,
                low,
                high,
            ) = example

            print(
                f"FAIL    "
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

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
