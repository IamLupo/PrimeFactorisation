import math
import random


EXPERIMENT = 58


def is_prime(n):
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    ]

    for p in small_primes:
        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    r = 0

    while d % 2 == 0:
        d //= 2
        r += 1

    for a in [2, 3, 5, 7, 11, 13, 17]:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(r - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(bits):
    while True:
        n = random.getrandbits(bits)

        n |= 1
        n |= 1 << (bits - 1)

        if is_prime(n):
            return n


def generate_case(bits):
    while True:
        p = random_prime(bits)
        q = random_prime(bits)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q
        s = math.isqrt(n)

        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "y": y,
        }


def delta_mod(case, y, modulus):
    s = case["s"] % modulus
    n = case["N"] % modulus

    z = (
        2 * s
        + 1
        + y
    ) % modulus

    return (
        z * z
        - 4 * n
    ) % modulus


def square_residues(modulus):
    result = set()

    for x in range(modulus):
        result.add(
            (x * x) % modulus
        )

    return result


def local_allowed_residues(case, modulus):
    squares = square_residues(
        modulus
    )

    allowed = []

    for y in range(modulus):
        if (
            delta_mod(
                case,
                y,
                modulus,
            )
            in squares
        ):
            allowed.append(y)

    return allowed


def build_wheel(case, moduli):
    """
    Build allowed y residues modulo the product of moduli.

    The implementation progressively combines residue classes
    instead of scanning the full CRT product.
    """

    residues = [0]
    current_modulus = 1

    for modulus in moduli:
        local = local_allowed_residues(
            case,
            modulus,
        )

        new_residues = []

        for base in residues:
            for r in local:
                # Find x = base + current_modulus*t
                # satisfying x == r mod modulus.
                #
                # All moduli here are pairwise coprime.
                inv = pow(
                    current_modulus,
                    -1,
                    modulus,
                )

                t = (
                    (r - base)
                    * inv
                ) % modulus

                value = (
                    base
                    + current_modulus * t
                )

                new_residues.append(
                    value
                )

        current_modulus *= modulus

        residues = sorted(
            set(
                x % current_modulus
                for x in new_residues
            )
        )

    return current_modulus, residues


def wheel_survives(y, residues, modulus):
    r = y % modulus

    return r in residues


def sieve_count_until(case, target_y, residues, modulus):
    count = 0

    for y in range(target_y + 1):
        if wheel_survives(
            y,
            residues,
            modulus,
        ):
            count += 1

    return count


def rank_with_wheel(case, residues, modulus):
    target = case["y"]

    rank = sieve_count_until(
        case,
        target,
        residues,
        modulus,
    )

    return rank


def naive_rank(case):
    return case["y"] + 1


def density(modulus, residues):
    return (
        len(residues) / modulus
    )


def test_wheel_density(cases, moduli):
    print("============================================================")
    print("TEST 1: CRT WHEEL DENSITY")
    print("============================================================")

    densities = []

    for case in cases[:20]:
        modulus, residues = build_wheel(
            case,
            moduli,
        )

        d = density(
            modulus,
            residues,
        )

        densities.append(d)

    print(
        f"moduli = {moduli}"
    )

    print(
        f"wheel modulus = "
        f"{build_wheel(cases[0], moduli)[0]}"
    )

    print(
        f"min density  = "
        f"{min(densities):.10f}"
    )

    print(
        f"max density  = "
        f"{max(densities):.10f}"
    )

    print(
        f"mean density = "
        f"{sum(densities) / len(densities):.10f}"
    )

    print()


def test_true_y(cases, moduli):
    print("============================================================")
    print("TEST 2: TRUE y ALWAYS SURVIVES")
    print("============================================================")

    failures = 0

    for case in cases:
        modulus, residues = build_wheel(
            case,
            moduli,
        )

        if not wheel_survives(
            case["y"],
            residues,
            modulus,
        ):
            failures += 1

    print(
        f"failures = {failures}/{len(cases)}"
    )

    print()


def test_rank_reduction(cases, moduli):
    print("============================================================")
    print("TEST 3: TRUE y RANK")
    print("============================================================")

    ratios = []
    ranks = []

    for case in cases:
        modulus, residues = build_wheel(
            case,
            moduli,
        )

        raw = naive_rank(case)

        rank = rank_with_wheel(
            case,
            residues,
            modulus,
        )

        ranks.append(rank)

        if raw > 0:
            ratios.append(
                rank / raw
            )

    print(
        f"minimum sieve rank = "
        f"{min(ranks)}"
    )

    print(
        f"maximum sieve rank = "
        f"{max(ranks)}"
    )

    print(
        f"mean sieve rank    = "
        f"{sum(ranks) / len(ranks):.2f}"
    )

    print()

    print(
        f"rank/raw minimum = "
        f"{min(ratios):.10f}"
    )

    print(
        f"rank/raw maximum = "
        f"{max(ratios):.10f}"
    )

    print(
        f"rank/raw mean    = "
        f"{sum(ratios) / len(ratios):.10f}"
    )

    print()


def test_extended_wheels(cases):
    print("============================================================")
    print("TEST 4: WHEEL GROWTH")
    print("============================================================")

    configurations = [
        [8],
        [16],
        [16, 3],
        [16, 3, 5],
        [16, 3, 5, 7],
        [16, 3, 5, 7, 11],
        [16, 3, 5, 7, 11, 13],
        [16, 3, 5, 7, 11, 13, 17],
        [16, 3, 5, 7, 11, 13, 17, 19],
    ]

    for moduli in configurations:
        densities = []
        rank_ratios = []

        for case in cases[:20]:
            modulus, residues = build_wheel(
                case,
                moduli,
            )

            d = density(
                modulus,
                residues,
            )

            rank = rank_with_wheel(
                case,
                residues,
                modulus,
            )

            raw = naive_rank(case)

            densities.append(d)

            rank_ratios.append(
                rank / raw
            )

        print(
            f"moduli={str(moduli):45s}  "
            f"density={sum(densities) / len(densities):.10f}  "
            f"rank_ratio={sum(rank_ratios) / len(rank_ratios):.10f}"
        )

    print()


def test_gap_structure(cases):
    print("============================================================")
    print("TEST 5: q-p VS SIEVE RANK")
    print("============================================================")

    values = []

    for case in cases:
        gap = case["q"] - case["p"]
        y = case["y"]

        values.append(
            (
                y,
                gap,
                y + 1,
            )
        )

    values.sort()

    print(
        "First 20 sorted cases:"
    )

    for y, gap, yp1 in values[:20]:
        print(
            f"y={y:8d}  "
            f"q-p={gap:8d}  "
            f"y+1={yp1:8d}"
        )

    print()


def print_examples(cases, moduli):
    print("============================================================")
    print("TEST 6: EXAMPLES")
    print("============================================================")

    for case in cases[:10]:
        modulus, residues = build_wheel(
            case,
            moduli,
        )

        rank = rank_with_wheel(
            case,
            residues,
            modulus,
        )

        raw = naive_rank(case)

        print()
        print(
            f"N={case['N']}"
        )

        print(
            f"s={case['s']}"
        )

        print(
            f"true y={case['y']}"
        )

        print(
            f"raw rank={raw}"
        )

        print(
            f"wheel modulus={modulus}"
        )

        print(
            f"wheel residues="
            f"{len(residues)}"
        )

        print(
            f"sieve rank={rank}"
        )

        print(
            f"compression="
            f"{1 - rank / raw:.8f}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(585858)

    CASE_COUNT = 100
    PRIME_BITS = 18

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = [
        generate_case(PRIME_BITS)
        for _ in range(CASE_COUNT)
    ]

    main_moduli = [
        16,
        3,
        5,
        7,
        11,
        13,
        17,
        19,
        23,
    ]

    test_wheel_density(
        cases,
        main_moduli,
    )

    test_true_y(
        cases,
        main_moduli,
    )

    test_rank_reduction(
        cases,
        main_moduli,
    )

    test_extended_wheels(
        cases,
    )

    test_gap_structure(
        cases,
    )

    print_examples(
        cases,
        main_moduli,
    )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
