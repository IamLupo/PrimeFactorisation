import math
import random


EXPERIMENT = 57


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
            "gap": q - p,
        }


def delta(case, y):
    s = case["s"]
    n = case["N"]

    z = 2 * s + 1 + y

    return z * z - 4 * n


def is_square(n):
    if n < 0:
        return False

    r = math.isqrt(n)

    return r * r == n


def square_residues(modulus):
    """
    All quadratic residues modulo modulus.
    """

    residues = set()

    for x in range(modulus):
        residues.add(
            (x * x) % modulus
        )

    return residues


def allowed_y_residues(case, modulus):
    """
    Since

        Delta(y) = (2s+1+y)^2 - 4N

    must be a square modulo modulus, determine all y
    residues modulo modulus that survive.
    """

    residues = square_residues(modulus)

    s = case["s"] % modulus
    n = case["N"] % modulus

    allowed = []

    for y in range(modulus):
        z = (
            2 * s
            + 1
            + y
        ) % modulus

        value = (
            z * z
            - 4 * n
        ) % modulus

        if value in residues:
            allowed.append(y)

    return set(allowed)


def combined_allowed_residues(case, moduli):
    """
    Return y residues modulo the product of pairwise
    coprime moduli that survive every local square test.
    """

    modulus = math.prod(moduli)

    local_sets = {
        m: allowed_y_residues(case, m)
        for m in moduli
    }

    allowed = []

    for y in range(modulus):
        valid = True

        for m in moduli:
            if y % m not in local_sets[m]:
                valid = False
                break

        if valid:
            allowed.append(y)

    return modulus, set(allowed)


def build_sieve(moduli):
    """
    Precompute square residues for each modulus.
    """

    return {
        m: square_residues(m)
        for m in moduli
    }


def survives_sieve(case, y, sieve):
    """
    Check Delta(y) against all precomputed local square
    residue sets.
    """

    s = case["s"]
    n = case["N"]

    for modulus, residues in sieve.items():
        z = (
            2 * s
            + 1
            + y
        ) % modulus

        value = (
            z * z
            - 4 * n
        ) % modulus

        if value not in residues:
            return False

    return True


def raw_y_search(case):
    """
    Ordinary y search starting at y=0.
    """

    tested = 0

    while True:
        value = delta(
            case,
            tested,
        )

        if is_square(value):
            return tested, tested + 1

        tested += 1


def sieved_y_search(case, sieve):
    """
    y search with modular square-residue filtering.
    """

    y = 0
    full_tests = 0

    while True:
        if survives_sieve(
            case,
            y,
            sieve,
        ):
            full_tests += 1

            value = delta(
                case,
                y,
            )

            if is_square(value):
                return y, full_tests

        y += 1


def test_local_density(cases):
    print("============================================================")
    print("TEST 1: LOCAL SQUARE-RESIDUE DENSITY")
    print("============================================================")

    moduli = [
        4,
        8,
        16,
        3,
        5,
        7,
        11,
        13,
    ]

    for modulus in moduli:
        density_values = []

        for case in cases:
            allowed = allowed_y_residues(
                case,
                modulus,
            )

            density = (
                len(allowed)
                / modulus
            )

            density_values.append(
                density
            )

        print(
            f"mod={modulus:3d}  "
            f"min={min(density_values):.6f}  "
            f"max={max(density_values):.6f}  "
            f"mean={sum(density_values) / len(density_values):.6f}"
        )

    print()


def test_combined_density(cases):
    print("============================================================")
    print("TEST 2: COMBINED SIEVE DENSITY")
    print("============================================================")

    configurations = [
        [4, 8],
        [8, 3],
        [8, 3, 5],
        [8, 3, 5, 7],
        [16, 3, 5, 7],
        [16, 3, 5, 7, 11],
        [16, 3, 5, 7, 11, 13],
    ]

    for moduli in configurations:
        densities = []

        for case in cases[:100]:
            modulus, allowed = (
                combined_allowed_residues(
                    case,
                    moduli,
                )
            )

            density = (
                len(allowed)
                / modulus
            )

            densities.append(
                density
            )

        print(
            f"moduli={str(moduli):28s}  "
            f"mean density="
            f"{sum(densities) / len(densities):.8f}"
        )

    print()


def test_search_reduction(cases):
    print("============================================================")
    print("TEST 3: RAW SEARCH VS SIEVED SEARCH")
    print("============================================================")

    configurations = [
        [8],
        [8, 3],
        [8, 3, 5],
        [8, 3, 5, 7],
        [8, 3, 5, 7, 11],
        [16, 3, 5, 7, 11, 13],
    ]

    for moduli in configurations:
        sieve = build_sieve(
            moduli
        )

        raw_total = 0
        sieve_total = 0

        for case in cases:
            raw_y, _ = raw_y_search(
                case
            )

            sieved_y, full_tests = (
                sieved_y_search(
                    case,
                    sieve,
                )
            )

            if raw_y != sieved_y:
                print(
                    "ERROR: sieve changed answer"
                )
                return

            raw_total += raw_y + 1
            sieve_total += full_tests

        print(
            f"moduli={str(moduli):28s}  "
            f"raw tests={raw_total:10d}  "
            f"square tests={sieve_total:10d}  "
            f"reduction="
            f"{1 - sieve_total / raw_total:.6f}"
        )

    print()


def test_true_y_survival(cases):
    print("============================================================")
    print("TEST 4: TRUE y SURVIVAL")
    print("============================================================")

    configurations = [
        [4],
        [8],
        [16],
        [8, 3, 5],
        [16, 3, 5, 7, 11, 13],
    ]

    for moduli in configurations:
        sieve = build_sieve(
            moduli
        )

        failures = 0

        for case in cases:
            if not survives_sieve(
                case,
                case["y"],
                sieve,
            ):
                failures += 1

        print(
            f"moduli={str(moduli):28s}  "
            f"true-y failures={failures}"
        )

    print()


def test_y_residue_patterns(cases):
    print("============================================================")
    print("TEST 5: y MODULAR PATTERNS")
    print("============================================================")

    moduli = [
        8,
        16,
        3,
        5,
        7,
        11,
        13,
    ]

    for modulus in moduli:
        counts = {}

        for case in cases:
            residue = case["y"] % modulus

            counts[residue] = (
                counts.get(residue, 0) + 1
            )

        print()
        print(
            f"mod={modulus}  "
            f"distinct y residues="
            f"{len(counts)}"
        )

        for residue in sorted(counts):
            print(
                f"  y={residue:3d}  "
                f"count={counts[residue]}"
            )

    print()


def test_v2_gap_vs_y_mod8(cases):
    print("============================================================")
    print("TEST 6: v2(q-p) VS y MOD 8")
    print("============================================================")

    groups = {}

    for case in cases:
        key = case["y"] % 8
        value = (
            case["gap"]
        )

        t = 999999 if value == 0 else (
            value & -value
        ).bit_length() - 1

        groups.setdefault(
            key,
            set(),
        ).add(t)

    for residue in sorted(groups):
        print(
            f"y mod 8={residue}  "
            f"v2(q-p) values="
            f"{sorted(groups[residue])}"
        )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 7: EXAMPLES")
    print("============================================================")

    moduli = [
        16,
        3,
        5,
        7,
        11,
        13,
    ]

    sieve = build_sieve(
        moduli
    )

    for case in cases[:10]:
        raw_y, _ = raw_y_search(
            case
        )

        sieved_y, full_tests = (
            sieved_y_search(
                case,
                sieve,
            )
        )

        print()
        print(
            f"N={case['N']} "
            f"s={case['s']}"
        )

        print(
            f"true y={case['y']}"
        )

        print(
            f"raw y={raw_y}"
        )

        print(
            f"sieved y={sieved_y}"
        )

        print(
            f"full square tests="
            f"{full_tests}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(575757)

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

    test_local_density(cases)
    test_combined_density(cases)
    test_search_reduction(cases)
    test_true_y_survival(cases)
    test_y_residue_patterns(cases)
    test_v2_gap_vs_y_mod8(cases)
    print_examples(cases)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
