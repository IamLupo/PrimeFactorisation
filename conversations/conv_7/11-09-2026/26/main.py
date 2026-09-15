import math
import random


EXPERIMENT = 47


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

        if s <= 0:
            continue

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": n - s * s,
            "gap": q - p,
            "gap_v2": v2(q - p),
        }


def v2(n):
    if n == 0:
        return 999999

    n = abs(n)

    return (n & -n).bit_length() - 1


def fingerprint(case, k):
    modulus = 1 << k

    return (
        case["N"] % modulus,
        case["s"] % modulus,
    )


def find_collision(cases, k):
    """
    Search for identical (N,s) mod 2^k
    but different v2(q-p).
    """

    groups = {}

    for case in cases:
        key = fingerprint(case, k)

        target = case["gap_v2"]

        if key not in groups:
            groups[key] = case
            continue

        old = groups[key]

        if old["gap_v2"] != target:
            return old, case, key

    return None, None, None


def count_conflicting_groups(cases, k):
    groups = {}

    for case in cases:
        key = fingerprint(case, k)

        if key not in groups:
            groups[key] = {case["gap_v2"]}
        else:
            groups[key].add(case["gap_v2"])

    conflicts = sum(
        1
        for values in groups.values()
        if len(values) > 1
    )

    return len(groups), conflicts


def print_case(case):
    print(f"N       = {case['N']}")
    print(f"p       = {case['p']}")
    print(f"q       = {case['q']}")
    print(f"s       = {case['s']}")
    print(f"D       = {case['D']}")
    print(f"q-p     = {case['gap']}")
    print(f"v2(q-p) = {case['gap_v2']}")


def run_scale(bits, count):
    print("============================================================")
    print(f"PRIME SIZE = {bits} BITS")
    print(f"CASES      = {count}")
    print("============================================================")

    cases = []

    for _ in range(count):
        cases.append(
            generate_case(bits)
        )

    print()
    print("Collision statistics:")

    first_collision = {}

    for k in range(10, 19):
        groups, conflicts = (
            count_conflicting_groups(
                cases,
                k,
            )
        )

        print(
            f"k={k:2d}  "
            f"groups={groups:7d}  "
            f"conflicting={conflicts:6d}"
        )

        if conflicts:
            old, new, key = find_collision(
                cases,
                k,
            )

            first_collision[k] = (
                old,
                new,
                key,
            )

    print()

    if first_collision:
        print("First observed collision at each k:")

        for k in sorted(first_collision):
            old, new, key = first_collision[k]

            print()
            print(f"k={k}")
            print(f"fingerprint={key}")

            print()
            print("CASE A")
            print_case(old)

            print()
            print("CASE B")
            print_case(new)

    else:
        print("No conflicting fingerprints found.")

    print()

    return cases


def targeted_collision_search(bits, k, attempts):
    """
    Generate cases until either a collision is found
    or the attempt limit is reached.
    """

    groups = {}

    for i in range(attempts):
        case = generate_case(bits)

        key = fingerprint(case, k)
        target = case["gap_v2"]

        if key in groups:
            old = groups[key]

            if old["gap_v2"] != target:
                return i + 1, old, case, key

        else:
            groups[key] = case

    return attempts, None, None, None


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(474747)

    print()
    print("============================================================")
    print("TARGETED COLLISION SEARCH")
    print("============================================================")
    print()

    configurations = [
        (16, 100000),
        (20, 100000),
        (24, 100000),
        (28, 100000),
        (30, 100000),
    ]

    for bits, count in configurations:
        cases = run_scale(
            bits,
            count,
        )

        del cases

    print()
    print("============================================================")
    print("TARGETED SEARCH AT k=13")
    print("============================================================")

    for bits in [16, 20, 24, 28, 30, 34]:
        attempts, old, new, key = (
            targeted_collision_search(
                bits,
                13,
                250000,
            )
        )

        print()
        print(
            f"bits={bits:2d}  "
            f"attempts={attempts}"
        )

        if old is None:
            print("No collision found.")
        else:
            print("COLLISION FOUND")
            print(f"fingerprint={key}")

            print()
            print("CASE A")
            print_case(old)

            print()
            print("CASE B")
            print_case(new)

    print()
    print("============================================================")
    print("TARGETED SEARCH AT k=14")
    print("============================================================")

    for bits in [20, 24, 28, 30, 34]:
        attempts, old, new, key = (
            targeted_collision_search(
                bits,
                14,
                250000,
            )
        )

        print()
        print(
            f"bits={bits:2d}  "
            f"attempts={attempts}"
        )

        if old is None:
            print("No collision found.")
        else:
            print("COLLISION FOUND")
            print(f"fingerprint={key}")

            print()
            print("CASE A")
            print_case(old)

            print()
            print("CASE B")
            print_case(new)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
