import math
import random

from sympy import isprime


EXPERIMENT = 99


def generate_semiprime(
    min_prime=100,
    max_prime=10000
):
    while True:
        p = random.randint(
            min_prime,
            max_prime
        )

        q = random.randint(
            min_prime,
            max_prime
        )

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
        coefficient = (
            1
            if j % 2 == 0
            else n - 1
        )

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


def simple_z_prediction(p, k):
    """
    Valid in the regime where the relevant
    degree is below p.

    r = (k-2) mod p

    Need the smallest odd d > r.

    d = 2z - 1
    """

    r = (
        k - 2
    ) % p

    if r % 2 == 0:
        d = r + 1
    else:
        d = r + 2

    z = (
        d + 1
    ) // 2

    return z, d, r


def verify_prediction(
    n,
    p,
    q,
    k
):
    predicted_z, predicted_d, residue = (
        simple_z_prediction(
            p,
            k
        )
    )

    g = gcd_value(
        n,
        predicted_z,
        k
    )

    return {
        "k": k,
        "z": predicted_z,
        "d": predicted_d,
        "residue": residue,
        "gcd": g,
        "class": classify_gcd(
            g,
            p,
            q,
            n
        )
    }


def run_experiment(
    cases=500
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

    print()

    tested = 0
    successful = 0

    p_hits = 0
    other_hits = 0
    n_hits = 0
    one_hits = 0

    examples_printed = 0

    for _ in range(cases):
        n, p, q = (
            generate_semiprime()
        )

        for k in power_of_two_values(
            p,
            q
        ):
            tested += 1

            result = verify_prediction(
                n,
                p,
                q,
                k
            )

            g = result["gcd"]
            classification = result["class"]

            if classification == "p":
                p_hits += 1

            elif classification == "n":
                n_hits += 1

            elif classification == "1":
                one_hits += 1

            else:
                other_hits += 1

            if g != 1:
                successful += 1

            if examples_printed < 30:
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
                    f"  (k-2) mod p="
                    f"{result['residue']}"
                )

                print(
                    f"  predicted d="
                    f"{result['d']}"
                )

                print(
                    f"  predicted z="
                    f"{result['z']}"
                )

                print(
                    f"  gcd="
                    f"{g}"
                )

                print(
                    f"  class="
                    f"{classification}"
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
        f"NONTRIVIAL HITS: "
        f"{successful}"
    )

    print(
        f"gcd = p: "
        f"{p_hits}"
    )

    print(
        f"gcd = n: "
        f"{n_hits}"
    )

    print(
        f"gcd = 1: "
        f"{one_hits}"
    )

    print(
        f"OTHER: "
        f"{other_hits}"
    )

    print()

    print(
        "IMPORTANT:"
    )

    print(
        "This experiment does NOT scan z."
    )

    print(
        "It directly predicts z from"
    )

    print(
        "(k - 2) mod p."
    )

    print()

    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    random.seed(990011)

    run_experiment(
        cases=500
    )


if __name__ == "__main__":
    main()
