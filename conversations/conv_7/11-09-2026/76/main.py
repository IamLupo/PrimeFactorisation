import math
import random


EXPERIMENT = 96


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
            return False

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


def next_power_of_two_below_or_equal(n):
    return 1 << (n.bit_length() - 1)


def arthseq_value_mod_n(
    n,
    z,
    k
):
    """
    x = n-1
    y = 1

    coefficient sequence:

        1, n-1, 1, n-1, ...

    with 2z coefficients.

    Degree = 2z-1.

    Evaluate the Newton-form ArthSeq at k.
    We only need the value modulo n.
    """

    degree = 2 * z - 1

    if k - 1 < degree:
        return 0

    total = 0

    for j in range(
        degree + 1
    ):
        if j % 2 == 0:
            coefficient = 1
        else:
            coefficient = n - 1

        total += (
            coefficient
            * math.comb(
                k - 1,
                j
            )
        )

        total %= n

    return total


def gcd_v(
    n,
    z,
    k
):
    value = arthseq_value_mod_n(
        n,
        z,
        k
    )

    return math.gcd(
        value,
        n
    )


def classify_gcd(
    g,
    p,
    q,
    n
):
    if g == 1:
        return "1"

    if g == p:
        return "p"

    if g == q:
        return "q"

    if g == n:
        return "n"

    return "other"


def predicted_min_z(
    k,
    p
):
    if k <= p + 1:
        return None

    return (
        (k - p + 1) // 2
    )


def find_first_factor_hit(
    n,
    p,
    q,
    k,
    max_z
):
    for z in range(
        1,
        max_z + 1
    ):
        g = gcd_v(
            n,
            z,
            k
        )

        if g != 1:
            return z, g

    return None, 1


def analyze_power(
    n,
    p,
    q,
    k,
    max_z
):
    if not (
        p < k < q
    ):
        return None

    predicted = predicted_min_z(
        k,
        p
    )

    actual_z, actual_gcd = (
        find_first_factor_hit(
            n,
            p,
            q,
            k,
            max_z
        )
    )

    return {
        "k": k,
        "predicted_z": predicted,
        "actual_z": actual_z,
        "gcd": actual_gcd,
        "class": classify_gcd(
            actual_gcd,
            p,
            q,
            n
        )
    }


def run_case(
    n,
    p,
    q,
    max_z=100
):
    s = math.isqrt(n)

    k0 = next_power_of_two_below_or_equal(
        s
    )

    k1 = 2 * k0

    print(
        f"N={n} "
        f"p={p} "
        f"q={q} "
        f"sqrtN={s}"
    )

    print(
        f"k0={k0}"
    )

    print(
        f"k1={k1}"
    )

    print()

    for k in (
        k0,
        k1
    ):
        result = analyze_power(
            n,
            p,
            q,
            k,
            max_z
        )

        if result is None:
            print(
                f"k={k}: "
                f"NOT BETWEEN p AND q"
            )

            continue

        print(
            f"k={k}: "
            f"predicted_z="
            f"{result['predicted_z']} "
            f"actual_z="
            f"{result['actual_z']} "
            f"gcd="
            f"{result['gcd']} "
            f"class="
            f"{result['class']}"
        )

    print()


def run_experiment(
    cases=500,
    max_z=100
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

    successful = 0

    predicted_matches = 0

    wrong_factor = 0

    no_factor = 0

    k0_hits = 0
    k1_hits = 0

    examples_printed = 0

    for case in range(
        1,
        cases + 1
    ):
        n, p, q = (
            generate_semiprime()
        )

        s = math.isqrt(n)

        k0 = (
            next_power_of_two_below_or_equal(
                s
            )
        )

        k1 = 2 * k0

        for k_index, k in enumerate(
            [k0, k1]
        ):
            if not (
                p < k < q
            ):
                continue

            tested += 1

            predicted = (
                predicted_min_z(
                    k,
                    p
                )
            )

            actual_z, g = (
                find_first_factor_hit(
                    n,
                    p,
                    q,
                    k,
                    max_z
                )
            )

            if actual_z is None:
                no_factor += 1
                continue

            successful += 1

            if actual_z == predicted:
                predicted_matches += 1

            if g not in (
                p,
                q,
                n
            ):
                wrong_factor += 1

            if k_index == 0:
                k0_hits += 1
            else:
                k1_hits += 1

            if examples_printed < 20:
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
                    f"  sqrtN={s}"
                )

                print(
                    f"  k={k}"
                )

                print(
                    f"  predicted z="
                    f"{predicted}"
                )

                print(
                    f"  actual z="
                    f"{actual_z}"
                )

                print(
                    f"  gcd="
                    f"{g}"
                )

                print(
                    f"  class="
                    f"{classify_gcd(g,p,q,n)}"
                )

                print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(
        f"POWER-OF-TWO k CASES TESTED: "
        f"{tested}"
    )

    print(
        f"SUCCESSFUL FACTOR HITS: "
        f"{successful}"
    )

    print(
        f"NO HIT WITHIN MAX Z: "
        f"{no_factor}"
    )

    print(
        f"EXACT PREDICTED-Z MATCHES: "
        f"{predicted_matches}"
    )

    print(
        f"WRONG FACTOR CLASSIFICATIONS: "
        f"{wrong_factor}"
    )

    print(
        f"k0 HITS: "
        f"{k0_hits}"
    )

    print(
        f"k1 HITS: "
        f"{k1_hits}"
    )

    if tested:
        print(
            f"HIT RATE: "
            f"{successful / tested:.8f}"
        )

        print(
            f"PREDICTED-Z MATCH RATE: "
            f"{predicted_matches / tested:.8f}"
        )

    print()

    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    random.seed(960011)

    run_experiment(
        cases=500,
        max_z=100
    )


if __name__ == "__main__":
    main()
