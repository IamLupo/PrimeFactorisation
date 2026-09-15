import math
import random


EXPERIMENT = 46


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

        if s <= 4096:
            continue

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": n - s * s,
            "gap": q - p,
        }


def v2(n):
    if n == 0:
        return 999999

    n = abs(n)

    return (n & -n).bit_length() - 1


def fingerprint_ns(case, k):
    modulus = 1 << k

    return (
        case["N"] % modulus,
        case["s"] % modulus,
    )


def fingerprint_ns_d(case, k):
    modulus = 1 << k

    return (
        case["N"] % modulus,
        case["s"] % modulus,
        case["D"] % modulus,
    )


def fingerprint_extended(case, k):
    """
    Add simple combinations which are still entirely
    computable from N and s.
    """

    modulus = 1 << k

    n = case["N"] % modulus
    s = case["s"] % modulus

    d = (n - s * s) % modulus

    return (
        n,
        s,
        d,
        (n - s) % modulus,
        (n + s) % modulus,
        (2 * s + 1) % modulus,
        (s * s) % modulus,
    )


def count_conflicting_groups(cases, k, fingerprint_function):
    groups = {}

    for case in cases:
        key = fingerprint_function(
            case,
            k,
        )

        target = v2(case["gap"])

        if key not in groups:
            groups[key] = {target}
        else:
            groups[key].add(target)

    conflicting = 0
    total_groups = len(groups)

    for values in groups.values():
        if len(values) > 1:
            conflicting += 1

    return total_groups, conflicting


def find_collision(cases, k, fingerprint_function):
    groups = {}

    for case in cases:
        key = fingerprint_function(
            case,
            k,
        )

        target = v2(case["gap"])

        if key not in groups:
            groups[key] = (target, case)
            continue

        old_target, old_case = groups[key]

        if old_target != target:
            return old_case, case, key

    return None, None, None


def distribution(cases):
    counts = {}

    for case in cases:
        t = v2(case["gap"])

        counts[t] = counts.get(t, 0) + 1

    return counts


def print_case(case):
    print(f"N       = {case['N']}")
    print(f"p       = {case['p']}")
    print(f"q       = {case['q']}")
    print(f"s       = {case['s']}")
    print(f"D       = {case['D']}")
    print(f"q-p     = {case['gap']}")
    print(f"v2(q-p) = {v2(case['gap'])}")


def test_formula_candidates(cases, k):
    """
    Check whether simple observable 2-adic quantities
    correlate perfectly with v2(q-p).
    """

    modulus = 1 << k

    exact = 0

    for case in cases:
        n = case["N"] % modulus
        s = case["s"] % modulus
        d = (n - s * s) % modulus

        target = v2(case["gap"])

        candidates = {
            "v2(D)": v2(d),
            "v2(N-s)": v2((n - s) % modulus),
            "v2(N+s)": v2((n + s) % modulus),
            "v2(2s+1)": v2((2 * s + 1) % modulus),
        }

        # Count whether any individual candidate happens
        # to equal the target.
        if target in candidates.values():
            exact += 1

    return exact


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(464646)

    CASE_COUNT = 12000
    PRIME_BITS = 30

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = [
        generate_case(PRIME_BITS)
        for _ in range(CASE_COUNT)
    ]

    print("============================================================")
    print("TEST 1: v2(q-p) DISTRIBUTION")
    print("============================================================")

    counts = distribution(cases)

    for t in sorted(counts):
        print(
            f"v2(q-p)={t:2d}  "
            f"count={counts[t]:5d}"
        )

    print()

    print("============================================================")
    print("TEST 2: BASIC (N,s) FINGERPRINT")
    print("============================================================")

    for k in range(4, 17):
        groups, conflicting = (
            count_conflicting_groups(
                cases,
                k,
                fingerprint_ns,
            )
        )

        print(
            f"k={k:2d}  "
            f"groups={groups:6d}  "
            f"conflicting={conflicting:6d}  "
            f"deterministic="
            f"{'YES' if conflicting == 0 else 'NO'}"
        )

    print()

    print("============================================================")
    print("TEST 3: EXTENDED COMPUTABLE FINGERPRINT")
    print("============================================================")

    for k in range(4, 17):
        groups, conflicting = (
            count_conflicting_groups(
                cases,
                k,
                fingerprint_extended,
            )
        )

        print(
            f"k={k:2d}  "
            f"groups={groups:6d}  "
            f"conflicting={conflicting:6d}  "
            f"deterministic="
            f"{'YES' if conflicting == 0 else 'NO'}"
        )

    print()

    print("============================================================")
    print("TEST 4: SIMPLE v2 CANDIDATES")
    print("============================================================")

    for k in [8, 12, 16]:
        exact = test_formula_candidates(
            cases,
            k,
        )

        print(
            f"k={k:2d}  "
            f"at least one simple v2 candidate "
            f"matches target: "
            f"{exact}/{CASE_COUNT}"
        )

    print()

    print("============================================================")
    print("TEST 5: FIRST COLLISION")
    print("============================================================")

    for k in [8, 10, 12, 14, 16]:
        case_a, case_b, key = find_collision(
            cases,
            k,
            fingerprint_ns,
        )

        print()
        print(f"k={k}")

        if case_a is None:
            print("No collision found.")
        else:
            print("Fingerprint:")
            print(key)

            print()
            print("CASE A")
            print_case(case_a)

            print()
            print("CASE B")
            print_case(case_b)

    print()

    print("============================================================")
    print("TEST 6: DOES D ADD INFORMATION?")
    print("============================================================")

    for k in [8, 12, 16]:
        basic = 0
        extended = 0

        for case in cases:
            if fingerprint_ns(case, k) == (
                case["N"] % (1 << k),
                case["s"] % (1 << k),
            ):
                basic += 1

            if fingerprint_ns_d(case, k) == (
                case["N"] % (1 << k),
                case["s"] % (1 << k),
                case["D"] % (1 << k),
            ):
                extended += 1

        print(
            f"k={k:2d}  "
            f"basic={basic}  "
            f"with D={extended}"
        )

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
