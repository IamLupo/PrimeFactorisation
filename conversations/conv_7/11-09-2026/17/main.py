import math
import random


EXPERIMENT = 39


def is_prime(n):
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
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

    witnesses = [2, 3, 5, 7, 11, 13, 17]

    for a in witnesses:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        composite = True

        for _ in range(r - 1):
            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
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

        if s <= 4096:
            continue

        d = n - s * s

        x_p = s + 1 - p
        x_q = s + 1 - q

        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "x_p": x_p,
            "x_q": x_q,
            "y": y,
        }


def P(n, s, d, x, y):
    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def residue_solutions(n_mod, s_mod):
    """
    Enumerate every (x,y) mod 4 satisfying P == 0 mod 4.
    """

    d_mod = (n_mod - s_mod * s_mod) % 4

    solutions = []

    for x in range(4):
        for y in range(4):
            value = (
                x * x
                + (y - 1) * x
                + (d_mod - s_mod)
                - (s_mod + 1) * y
            ) % 4

            if value == 0:
                solutions.append((x, y))

    return solutions


def x_classes_from_s_n(n_mod, s_mod):
    """
    Project the full (x,y) solution set onto x mod 4.
    """

    solutions = residue_solutions(
        n_mod,
        s_mod,
    )

    return tuple(sorted({
        x
        for x, y in solutions
    }))


def actual_x_classes(case):
    return tuple(sorted({
        case["x_p"] % 4,
        case["x_q"] % 4,
    }))


def print_residue_table():
    print("============================================================")
    print("COMPLETE MOD-4 ALGEBRAIC TABLE")
    print("============================================================")

    for s_mod in range(4):
        # Since N = pq with odd p,q, N is always odd.
        for n_mod in [1, 3]:

            solutions = residue_solutions(
                n_mod,
                s_mod,
            )

            x_classes = tuple(sorted({
                x
                for x, y in solutions
            }))

            print(
                f"s={s_mod}  "
                f"N={n_mod}  "
                f"x-classes={x_classes}  "
                f"solutions={solutions}"
            )

    print()


def test_random_cases(cases):
    print("============================================================")
    print("TEST 1: MOD-4 ROOT CLASSES")
    print("============================================================")

    correct = 0
    spurious = 0
    impossible = 0

    by_signature = {}

    for case in cases:
        n_mod = case["N"] % 4
        s_mod = case["s"] % 4

        predicted = x_classes_from_s_n(
            n_mod,
            s_mod,
        )

        actual = actual_x_classes(case)

        key = (n_mod, s_mod)

        if key not in by_signature:
            by_signature[key] = {
                "predicted": predicted,
                "actual_sets": set(),
            }

        by_signature[key]["actual_sets"].add(actual)

        if set(actual).issubset(set(predicted)):
            if set(actual) == set(predicted):
                correct += 1
            else:
                spurious += 1
        else:
            impossible += 1

    print(f"exact match       = {correct}")
    print(f"contains spurious = {spurious}")
    print(f"missing true root = {impossible}")

    print()

    print("Observed actual root sets by (N mod 4, s mod 4):")

    for key in sorted(by_signature):
        predicted = by_signature[key]["predicted"]
        actual_sets = sorted(
            by_signature[key]["actual_sets"]
        )

        print(
            f"{key}: "
            f"predicted={predicted}  "
            f"actual={actual_sets}"
        )

    print()


def test_y_structure(cases):
    print("============================================================")
    print("TEST 2: WHERE DOES THE SPURIOUS X CLASS COME FROM?")
    print("============================================================")

    structure = {}

    for case in cases:
        n_mod = case["N"] % 4
        s_mod = case["s"] % 4

        predicted = x_classes_from_s_n(
            n_mod,
            s_mod,
        )

        actual = actual_x_classes(case)

        extra = tuple(
            x
            for x in predicted
            if x not in actual
        )

        key = (
            n_mod,
            s_mod,
            actual,
            predicted,
            extra,
        )

        structure[key] = structure.get(key, 0) + 1

    for key, count in sorted(
        structure.items(),
        key=lambda item: (
            item[0],
            item[1],
        ),
    ):
        print(
            f"count={count:5d}  "
            f"N={key[0]}  "
            f"s={key[1]}  "
            f"actual={key[2]}  "
            f"predicted={key[3]}  "
            f"extra={key[4]}"
        )

    print()


def test_factor_parity_relationship(cases):
    print("============================================================")
    print("TEST 3: FACTOR RESIDUE RELATIONSHIPS")
    print("============================================================")

    counts = {}

    for case in cases:
        p4 = case["p"] % 4
        q4 = case["q"] % 4
        s4 = case["s"] % 4
        n4 = case["N"] % 4

        xp4 = case["x_p"] % 4
        xq4 = case["x_q"] % 4

        key = (
            p4,
            q4,
            s4,
            n4,
            xp4,
            xq4,
        )

        counts[key] = counts.get(key, 0) + 1

    for key, count in sorted(counts.items()):
        print(
            f"count={count:4d}  "
            f"p={key[0]} q={key[1]} "
            f"s={key[2]} N={key[3]} "
            f"x_p={key[4]} x_q={key[5]}"
        )

    print()


def test_universal_spurious_class(cases):
    print("============================================================")
    print("TEST 4: UNIVERSAL EXTRA CLASS?")
    print("============================================================")

    signatures = {}

    for case in cases:
        n4 = case["N"] % 4
        s4 = case["s"] % 4

        predicted = x_classes_from_s_n(
            n4,
            s4,
        )

        actual = actual_x_classes(case)

        extra = tuple(
            x
            for x in predicted
            if x not in actual
        )

        key = (
            n4,
            s4,
            actual,
        )

        if key not in signatures:
            signatures[key] = set()

        signatures[key].add(extra)

    for key in sorted(signatures):
        print(
            f"N={key[0]} "
            f"s={key[1]} "
            f"actual={key[2]} "
            f"extra_classes={sorted(signatures[key])}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(39039)

    CASE_COUNT = 6000
    PRIME_BITS = 30

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    print_residue_table()

    cases = []

    for _ in range(CASE_COUNT):
        cases.append(
            generate_case(PRIME_BITS)
        )

    test_random_cases(cases)
    test_y_structure(cases)
    test_factor_parity_relationship(cases)
    test_universal_spurious_class(cases)

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
