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
    s = case["s"]
    n = case["N"]

    z = (
        2 * s
        + 1
        + y
    ) % modulus

    n_mod = n % modulus

    return (
        z * z
        - 4 * n_mod
    ) % modulus


def build_square_tables(moduli):
    tables = {}

    for modulus in moduli:
        table = bytearray(modulus)

        for x in range(modulus):
            table[
                (x * x) % modulus
            ] = 1

        tables[modulus] = table

    return tables


def survives_sieve(case, y, tables):
    for modulus, table in tables.items():
        value = delta_mod(
            case,
            y,
            modulus,
        )

        if table[value] == 0:
            return False

    return True


def count_survivors_until(case, tables):
    """
    Count surviving y values from 0 through true y.
    """

    target = case["y"]

    survivors = 0
    full_square_tests = 0

    for y in range(target + 1):
        if survives_sieve(
            case,
            y,
            tables,
        ):
            survivors += 1
            full_square_tests += 1

            # The actual square test is only done for survivors.
            delta = (
                2 * case["s"]
                + 1
                + y
            ) ** 2 - 4 * case["N"]

            if delta < 0:
                raise RuntimeError(
                    "Negative discriminant encountered."
                )

            root = math.isqrt(delta)

            if root * root == delta:
                if y != target:
                    raise RuntimeError(
                        "Found an earlier square than true y."
                    )

    return survivors, full_square_tests


def actual_square_tests_raw(case):
    return case["y"] + 1


def test_true_y_survival(cases, tables):
    print("============================================================")
    print("TEST 1: TRUE y SURVIVAL")
    print("============================================================")

    failures = 0

    for case in cases:
        if not survives_sieve(
            case,
            case["y"],
            tables,
        ):
            failures += 1

    print(
        f"failures = {failures}/{len(cases)}"
    )

    print()


def test_density(cases, tables):
    print("============================================================")
    print("TEST 2: EMPIRICAL SIEVE DENSITY")
    print("============================================================")

    densities = []

    for case in cases:
        target = case["y"]

        if target == 0:
            continue

        survivors = 0

        for y in range(target + 1):
            if survives_sieve(
                case,
                y,
                tables,
            ):
                survivors += 1

        densities.append(
            survivors / (target + 1)
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


def test_rank(cases, tables):
    print("============================================================")
    print("TEST 3: TRUE y RANK AFTER SIEVE")
    print("============================================================")

    raw_total = 0
    sieve_total = 0

    ranks = []

    for case in cases:
        raw = actual_square_tests_raw(
            case
        )

        survivors, full_tests = (
            count_survivors_until(
                case,
                tables,
            )
        )

        raw_total += raw
        sieve_total += survivors
        ranks.append(survivors)

    print(
        f"raw candidates     = "
        f"{raw_total}"
    )

    print(
        f"sieved candidates  = "
        f"{sieve_total}"
    )

    print(
        f"reduction          = "
        f"{1 - sieve_total / raw_total:.8f}"
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


def compare_wheels(cases):
    print("============================================================")
    print("TEST 4: WHEEL STRENGTH")
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
        [16, 3, 5, 7, 11, 13, 17, 19, 23],
    ]

    for moduli in configurations:
        tables = build_square_tables(
            moduli
        )

        total_raw = 0
        total_sieved = 0

        for case in cases:
            raw = case["y"] + 1

            survivors = 0

            for y in range(case["y"] + 1):
                if survives_sieve(
                    case,
                    y,
                    tables,
                ):
                    survivors += 1

            total_raw += raw
            total_sieved += survivors

        reduction = (
            1
            - total_sieved / total_raw
        )

        print(
            f"{str(moduli):52s} "
            f"raw={total_raw:9d}  "
            f"sieved={total_sieved:8d}  "
            f"reduction={reduction:.8f}"
        )

    print()


def test_first_survivors(cases, tables):
    print("============================================================")
    print("TEST 5: FIRST SURVIVING y VALUES")
    print("============================================================")

    for index, case in enumerate(cases[:10]):
        survivors = []

        for y in range(case["y"] + 1):
            if survives_sieve(
                case,
                y,
                tables,
            ):
                survivors.append(y)

        print()
        print(
            f"CASE {index + 1}"
        )

        print(
            f"N       = {case['N']}"
        )

        print(
            f"s       = {case['s']}"
        )

        print(
            f"true y  = {case['y']}"
        )

        print(
            f"rank    = {len(survivors)}"
        )

        print(
            f"first survivors = "
            f"{survivors[:20]}"
        )

        if survivors:
            print(
                f"last survivor = "
                f"{survivors[-1]}"
            )

    print()


def test_y_structure(cases):
    print("============================================================")
    print("TEST 6: y STRUCTURE")
    print("============================================================")

    for modulus in [8, 16, 32]:
        residues = set(
            case["y"] % modulus
            for case in cases
        )

        print(
            f"mod={modulus:2d}  "
            f"distinct residues="
            f"{len(residues)}"
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

    moduli = [
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

    tables = build_square_tables(
        moduli
    )

    test_true_y_survival(
        cases,
        tables,
    )

    test_density(
        cases,
        tables,
    )

    test_rank(
        cases,
        tables,
    )

    compare_wheels(
        cases,
    )

    test_first_survivors(
        cases,
        tables,
    )

    test_y_structure(
        cases,
    )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
