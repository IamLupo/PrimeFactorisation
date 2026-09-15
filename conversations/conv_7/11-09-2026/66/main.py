import math
import random


EXPERIMENT = 84


def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    r = math.isqrt(n)
    d = 3

    while d <= r:
        if n % d == 0:
            return False
        d += 2

    return True


def generate_semiprime(min_value=1000, max_value=50000):
    while True:
        p = random.randrange(min_value, max_value)
        q = random.randrange(min_value, max_value)

        if p > q:
            p, q = q, p

        if p == q:
            continue

        if not is_prime(p):
            continue

        if not is_prime(q):
            continue

        return p, q


def construct_parameters(p, q):
    n = p * q
    s = math.isqrt(n)
    d = n - s * s

    x_p = s + 1 - p
    y = p + q - 2 * s - 1

    return n, s, d, x_p, y


def mod_inverse(a, m):
    return pow(a % m, -1, m)


def orbit_displacement(x, s, d, m):
    denominator = (x - s - 1) % m

    if math.gcd(denominator, m) != 1:
        return None

    numerator = (
        -(x * x - x + d - s)
    ) % m

    y = (
        numerator
        * mod_inverse(denominator, m)
    ) % m

    x2 = (1 - y - x) % m

    displacement = (
        x2 - (x % m)
    ) % m

    return displacement


def crt(values, moduli):
    x = values[0]
    modulus = moduli[0]

    for i in range(1, len(values)):
        a = values[i]
        m = moduli[i]

        t = (
            (a - x)
            * mod_inverse(modulus, m)
        ) % m

        x += modulus * t
        modulus *= m

        x %= modulus

    return x, modulus


def true_crt_displacement(x_p, s, d, moduli):
    residues = []

    for m in moduli:
        displacement = orbit_displacement(
            x_p,
            s,
            d,
            m
        )

        if displacement is None:
            return None

        residues.append(displacement)

    return crt(
        residues,
        moduli
    )


def centered(value, modulus):
    if value > modulus // 2:
        return value - modulus

    return value


def admissible_residues(s, m):
    result = []

    for x in range(m):
        denominator = (
            x - s - 1
        ) % m

        if math.gcd(
            denominator,
            m
        ) == 1:
            result.append(x)

    return result


def random_crt_displacement(
    s,
    d,
    moduli,
    admissible
):
    residues = []

    for m in moduli:
        candidates = admissible[m]

        x = random.choice(candidates)

        displacement = orbit_displacement(
            x,
            s,
            d,
            m
        )

        if displacement is None:
            return None

        residues.append(displacement)

    return crt(
        residues,
        moduli
    )


def run_experiment(num_cases=500):
    moduli = [
        3,
        5,
        7,
        11,
        13,
        17,
        19,
        23
    ]

    _, crt_modulus = crt(
        [0] * len(moduli),
        moduli
    )

    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"MODULI: {moduli}")
    print(f"CRT MODULUS: {crt_modulus}")
    print(f"CASES: {num_cases}")
    print()

    admissible_cache = {}

    for m in moduli:
        admissible_cache[m] = None

    valid_cases = 0
    exact_modular = 0
    exact_signed = 0

    random_trials = 0
    random_negative_hits = 0
    random_positive_hits = 0

    for case in range(
        1,
        num_cases + 1
    ):
        p, q = generate_semiprime()

        (
            n,
            s,
            d,
            x_p,
            y
        ) = construct_parameters(
            p,
            q
        )

        result = true_crt_displacement(
            x_p,
            s,
            d,
            moduli
        )

        if result is None:
            continue

        displacement, modulus = result

        signed_displacement = centered(
            displacement,
            modulus
        )

        true_displacement = p - q
        true_gap = q - p

        valid_cases += 1

        if (
            displacement % modulus
            == true_displacement % modulus
        ):
            exact_modular += 1

        if signed_displacement == true_displacement:
            exact_signed += 1

        if case <= 10:
            print(f"CASE {case}")
            print(f"  N = {n}")
            print(f"  p = {p}")
            print(f"  q = {q}")
            print(f"  s = {s}")
            print(f"  d = {d}")
            print(f"  x_p = {x_p}")
            print(f"  y = {y}")
            print(f"  q-p = {true_gap}")
            print(f"  p-q = {true_displacement}")
            print(
                f"  CRT displacement = "
                f"{displacement}"
            )
            print(
                f"  Signed displacement = "
                f"{signed_displacement}"
            )
            print(
                f"  Exact signed match = "
                f"{signed_displacement == true_displacement}"
            )
            print()

        # Random control.
        if admissible_cache[moduli[0]] is None:
            for m in moduli:
                admissible_cache[m] = (
                    admissible_residues(
                        s,
                        m
                    )
                )

        # Rebuild cache for each N, because s changes.
        admissible = {}

        for m in moduli:
            admissible[m] = admissible_residues(
                s,
                m
            )

        for _ in range(100):
            random_result = random_crt_displacement(
                s,
                d,
                moduli,
                admissible
            )

            if random_result is None:
                continue

            random_value, random_modulus = (
                random_result
            )

            random_signed = centered(
                random_value,
                random_modulus
            )

            random_trials += 1

            if random_signed == -true_gap:
                random_negative_hits += 1

            if random_signed == true_gap:
                random_positive_hits += 1

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(
        f"VALID CASES: "
        f"{valid_cases}/{num_cases}"
    )

    print(
        "EXACT MODULAR DISPLACEMENT: "
        f"{exact_modular}/{valid_cases}"
    )

    print(
        "EXACT SIGNED DISPLACEMENT: "
        f"{exact_signed}/{valid_cases}"
    )

    print(
        "RANDOM -(q-p) HITS: "
        f"{random_negative_hits}/{random_trials}"
    )

    print(
        "RANDOM +(q-p) HITS: "
        f"{random_positive_hits}/{random_trials}"
    )

    if random_trials > 0:
        print(
            "RANDOM -(q-p) RATE: "
            f"{random_negative_hits / random_trials:.8f}"
        )

        print(
            "RANDOM +(q-p) RATE: "
            f"{random_positive_hits / random_trials:.8f}"
        )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    random.seed(840011)

    run_experiment(
        num_cases=500
    )


if __name__ == "__main__":
    main()