import math
import random


EXPERIMENT = 37


def is_prime(n):
    if n < 2:
        return False

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    r = 0

    while d % 2 == 0:
        r += 1
        d //= 2

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
        n |= (1 << (bits - 1))
        n |= 1

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

        # Keep the large-s regime used in the previous experiments.
        if s <= 4096:
            continue

        x_p = s + 1 - p
        y = p + q - 2 * s - 1
        d = n - s * s

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "x_p": x_p,
            "y": y,
        }


def H_value(n, s, x):
    d = n - s * s

    return 4 * x * x - (s + 4) * x + 3 * d - 3 * s


def fingerprint_basic(n, s, k):
    modulus = 1 << k

    return (
        n % modulus,
        s % modulus,
    )


def fingerprint_with_H(n, s, k, h_points):
    modulus = 1 << k

    values = [
        H_value(n, s, x) % modulus
        for x in h_points
    ]

    return (
        n % modulus,
        s % modulus,
        *values,
    )


def check_determinism(cases, fingerprint_function, k):
    groups = {}

    for case in cases:
        key = fingerprint_function(
            case["N"],
            case["s"],
            k,
        )

        target = case["x_p"] % 4

        if key not in groups:
            groups[key] = target
        elif groups[key] != target:
            return False, key, groups[key], target, case

    return True, None, None, None, None


def count_collisions(cases, fingerprint_function, k):
    groups = {}
    conflicting_groups = 0

    for case in cases:
        key = fingerprint_function(
            case["N"],
            case["s"],
            k,
        )

        target = case["x_p"] % 4

        if key not in groups:
            groups[key] = {target}
        else:
            groups[key].add(target)

    for values in groups.values():
        if len(values) > 1:
            conflicting_groups += 1

    return conflicting_groups


def print_example(case):
    print("Example collision:")
    print(f"N      = {case['N']}")
    print(f"p      = {case['p']}")
    print(f"q      = {case['q']}")
    print(f"s      = {case['s']}")
    print(f"D      = {case['D']}")
    print(f"x_p    = {case['x_p']}")
    print(f"x_p%4  = {case['x_p'] % 4}")
    print(f"y      = {case['y']}")


def find_collision(cases, fingerprint_function, k):
    groups = {}

    for case in cases:
        key = fingerprint_function(
            case["N"],
            case["s"],
            k,
        )

        target = case["x_p"] % 4

        if key not in groups:
            groups[key] = (target, case)
        else:
            old_target, old_case = groups[key]

            if old_target != target:
                return old_case, case

    return None, None


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(37037)

    CASE_COUNT = 12000
    PRIME_BITS = 30

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = []

    for _ in range(CASE_COUNT):
        cases.append(generate_case(PRIME_BITS))

    print("============================================================")
    print("TEST 1: BASIC FINGERPRINT")
    print("Fingerprint = (N mod 2^k, s mod 2^k)")
    print("Target      = x_p mod 4")
    print("============================================================")

    first_basic_failure = None

    for k in range(2, 21):
        conflict_count = count_collisions(
            cases,
            fingerprint_basic,
            k,
        )

        deterministic = conflict_count == 0

        print(
            f"k={k:2d}  "
            f"mod={1 << k:8d}  "
            f"conflicting_groups={conflict_count:6d}  "
            f"deterministic={'YES' if deterministic else 'NO'}"
        )

        if not deterministic and first_basic_failure is None:
            first_basic_failure = k

    print()

    if first_basic_failure is not None:
        old_case, new_case = find_collision(
            cases,
            fingerprint_basic,
            first_basic_failure,
        )

        print(
            f"First basic collision for k={first_basic_failure}"
        )

        print()
        print("CASE A")
        print_example(old_case)

        print()
        print("CASE B")
        print_example(new_case)

        print()

    print("============================================================")
    print("TEST 2: ADD H(0), H(1), H(2), H(3)")
    print("Fingerprint = (N, s, H(0), H(1), H(2), H(3)) mod 2^k")
    print("Target      = x_p mod 4")
    print("============================================================")

    h_points = [0, 1, 2, 3]

    first_H_failure = None

    for k in range(2, 21):
        conflict_count = count_collisions(
            cases,
            lambda n, s, kk: fingerprint_with_H(
                n,
                s,
                kk,
                h_points,
            ),
            k,
        )

        deterministic = conflict_count == 0

        print(
            f"k={k:2d}  "
            f"mod={1 << k:8d}  "
            f"conflicting_groups={conflict_count:6d}  "
            f"deterministic={'YES' if deterministic else 'NO'}"
        )

        if not deterministic and first_H_failure is None:
            first_H_failure = k

    print()

    if first_H_failure is not None:
        old_case, new_case = find_collision(
            cases,
            lambda n, s, kk: fingerprint_with_H(
                n,
                s,
                kk,
                h_points,
            ),
            first_H_failure,
        )

        print(
            f"First H-extended collision for k={first_H_failure}"
        )

        print()
        print("CASE A")
        print_example(old_case)

        print()
        print("CASE B")
        print_example(new_case)

        print()

    print("============================================================")
    print("TEST 3: VERIFY THE INFORMATION CONTENT OF H")
    print("============================================================")

    for k in [4, 6, 8, 10, 12, 14, 16]:
        modulus = 1 << k
        same_count = 0

        for case in cases:
            n = case["N"]
            s = case["s"]

            basic = fingerprint_basic(n, s, k)

            extended = fingerprint_with_H(
                n,
                s,
                k,
                h_points,
            )

            # H is mathematically determined by N and s.
            # This check confirms that explicitly on the data.
            reconstructed = (
                basic[0],
                basic[1],
                *[
                    H_value(n, s, x) % modulus
                    for x in h_points
                ],
            )

            if extended == reconstructed:
                same_count += 1

        print(
            f"k={k:2d}  "
            f"H-extension exactly reconstructed from (N,s): "
            f"{same_count}/{len(cases)}"
        )

    print()

    print("============================================================")
    print("TEST 4: DIRECT x_p MOD 4 DISTRIBUTION")
    print("============================================================")

    counts = [0, 0, 0, 0]

    for case in cases:
        counts[case["x_p"] % 4] += 1

    for residue in range(4):
        print(
            f"x_p mod 4 = {residue}: "
            f"{counts[residue]}"
        )

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
