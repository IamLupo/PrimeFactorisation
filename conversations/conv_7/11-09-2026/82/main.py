import math
import random
from sympy import isprime


EXPERIMENT = 102
CASES = 500
P_MIN = 100
P_MAX = 12000
Q_MAX = 16000


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


def next_power_of_two(x):
    z = 1
    while z < x:
        z *= 2
    return z


def dyadic_search(n, p, q, k):
    """
    Search only z = 1, 2, 4, 8, ...

    For x=n-1, y=1 and degree d=2z-1:

        v_k = -C(k-2, 2z-1)     when 2z-1 < k-1
        v_k = 0                 otherwise

    Therefore gcd(v_k, n) can be evaluated directly.
    """

    evaluations = 0
    z = 1

    while True:
        evaluations += 1
        degree = 2 * z - 1

        if degree >= k - 1:
            g = n
        else:
            value_mod_n = math.comb(k - 2, degree) % n
            value_mod_n = (-value_mod_n) % n
            g = math.gcd(value_mod_n, n)

        state = gcd_state(g, n, p, q)

        if state in ("p", "q"):
            return {
                "found": True,
                "z": z,
                "gcd": g,
                "state": state,
                "evaluations": evaluations,
            }

        if state == "n":
            return {
                "found": False,
                "z": z,
                "gcd": g,
                "state": state,
                "evaluations": evaluations,
            }

        z *= 2


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
    """
    Return all powers of two k satisfying p < k < q.
    """
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
    print()

    random.seed(102)

    total_k_cases = 0
    factor_hits = 0
    n_jumps = 0
    other_results = 0

    p_hits = 0
    q_hits = 0

    worst_evaluations = 0
    worst_case = None

    first_misses = []
    first_hits = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        ks = powers_between(p, q)

        for k in ks:
            total_k_cases += 1

            result = dyadic_search(n, p, q, k)

            if result["state"] == "p":
                factor_hits += 1
                p_hits += 1

                if len(first_hits) < 10:
                    first_hits.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            result["z"],
                            result["gcd"],
                            result["evaluations"],
                        )
                    )

            elif result["state"] == "q":
                factor_hits += 1
                q_hits += 1

            elif result["state"] == "n":
                n_jumps += 1

                if len(first_misses) < 20:
                    first_misses.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            result["z"],
                            result["evaluations"],
                        )
                    )

            else:
                other_results += 1

            if result["evaluations"] > worst_evaluations:
                worst_evaluations = result["evaluations"]
                worst_case = (
                    case_index,
                    n,
                    p,
                    q,
                    k,
                    result["z"],
                    result["state"],
                )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(f"semiprime cases              = {CASES}")
    print(f"(N,k) test cases             = {total_k_cases}")
    print(f"proper-factor hits           = {factor_hits}")
    print(f"  p hits                     = {p_hits}")
    print(f"  q hits                     = {q_hits}")
    print(f"direct 1 -> n jumps          = {n_jumps}")
    print(f"other results                = {other_results}")

    if total_k_cases:
        hit_rate = 100.0 * factor_hits / total_k_cases
        jump_rate = 100.0 * n_jumps / total_k_cases

        print(f"factor-hit rate              = {hit_rate:.2f}%")
        print(f"direct 1->n jump rate        = {jump_rate:.2f}%")

    print(f"worst gcd evaluations        = {worst_evaluations}")

    if worst_case is not None:
        print()
        print("WORST CASE")
        print(
            f"  case={worst_case[0]} "
            f"N={worst_case[1]} "
            f"p={worst_case[2]} "
            f"q={worst_case[3]} "
            f"k={worst_case[4]} "
            f"z={worst_case[5]} "
            f"state={worst_case[6]}"
        )

    if first_hits:
        print()
        print("-" * 60)
        print("EXAMPLE FACTOR HITS")
        print("-" * 60)

        for (
            case_index,
            n,
            p,
            q,
            k,
            z,
            g,
            evaluations,
        ) in first_hits:
            print(
                f"CASE {case_index}: "
                f"N={n} p={p} q={q} "
                f"k={k} z={z} "
                f"gcd={g} evals={evaluations}"
            )

    if first_misses:
        print()
        print("-" * 60)
        print("EXAMPLE DYADIC MISSES")
        print("-" * 60)

        for (
            case_index,
            n,
            p,
            q,
            k,
            z,
            evaluations,
        ) in first_misses:
            print(
                f"CASE {case_index}: "
                f"N={n} p={p} q={q} "
                f"k={k} "
                f"first_nontrivial_dyadic_z={z} "
                f"gcd=N "
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
