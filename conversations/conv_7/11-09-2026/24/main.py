import math
import random


EXPERIMENT = 45


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

        if n % 8 != 1:
            continue

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": math.isqrt(n),
        }


def v2(n):
    if n == 0:
        return 999999

    n = abs(n)

    return (n & -n).bit_length() - 1


def v2_capped(n, k):
    """
    v2 modulo precision 2^k.
    Values >= k are reported as k.
    """

    if n == 0:
        return k

    value = v2(n)

    return min(value, k)


def lift_square_roots(n, k):
    if k == 1:
        return [
            r
            for r in range(2)
            if (r * r - n) % 2 == 0
        ]

    if k == 2:
        return [
            r
            for r in range(4)
            if (r * r - n) % 4 == 0
        ]

    if n % 8 != 1:
        return []

    roots = [1, 3, 5, 7]
    modulus = 8

    for _ in range(3, k):
        new_modulus = modulus * 2
        new_roots = []

        for r in roots:
            r0 = r
            r1 = r + modulus

            if (r0 * r0 - n) % new_modulus == 0:
                new_roots.append(r0)

            if (r1 * r1 - n) % new_modulus == 0:
                new_roots.append(r1)

        roots = sorted(set(new_roots))
        modulus = new_modulus

    return roots


def analyze_root(case, r, k):
    p = case["p"]
    q = case["q"]

    target = v2_capped(q - p, k)

    rp = v2_capped(r - p, k)
    r_plus_p = v2_capped(r + p, k)

    rq = v2_capped(r - q, k)
    r_plus_q = v2_capped(r + q, k)

    sum_p = min(rp + r_plus_p, k)
    sum_q = min(rq + r_plus_q, k)

    return {
        "v2_pq": target,

        "v2_r-p": rp,
        "v2_r+p": r_plus_p,
        "sum_p": sum_p,

        "v2_r-q": rq,
        "v2_r+q": r_plus_q,
        "sum_q": sum_q,
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(454545)

    CASE_COUNT = 1000
    PRIME_BITS = 20

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
    print("TEST 1: PRODUCT IDENTITY")
    print("============================================================")

    for k in [8, 12, 16]:
        total_roots = 0

        p_failures = 0
        q_failures = 0

        for case in cases:
            roots = lift_square_roots(
                case["N"],
                k,
            )

            for r in roots:
                total_roots += 1

                result = analyze_root(
                    case,
                    r,
                    k,
                )

                if result["sum_p"] != result["v2_pq"]:
                    p_failures += 1

                if result["sum_q"] != result["v2_pq"]:
                    q_failures += 1

        print(
            f"k={k:2d}  "
            f"roots={total_roots:5d}  "
            f"p failures={p_failures:4d}  "
            f"q failures={q_failures:4d}"
        )

    print()

    print("============================================================")
    print("TEST 2: BRANCH DISTRIBUTION")
    print("============================================================")

    for k in [8, 12, 16]:
        distribution = {}

        for case in cases:
            target = v2_capped(
                case["q"] - case["p"],
                k,
            )

            roots = lift_square_roots(
                case["N"],
                k,
            )

            for r in roots:
                result = analyze_root(
                    case,
                    r,
                    k,
                )

                key = (
                    target,
                    result["v2_r-p"],
                    result["v2_r+p"],
                )

                distribution[key] = (
                    distribution.get(key, 0) + 1
                )

        print()
        print(f"k={k}")

        for key, count in sorted(
            distribution.items()
        ):
            print(
                f"v2(q-p)={key[0]:2d}  "
                f"v2(r-p)={key[1]:2d}  "
                f"v2(r+p)={key[2]:2d}  "
                f"count={count:5d}"
            )

    print()

    print("============================================================")
    print("TEST 3: DOES EVERY ROOT SPLIT THE VALUATION?")
    print("============================================================")

    for k in [8, 12, 16]:
        split_counts = {}

        for case in cases:
            target = v2_capped(
                case["q"] - case["p"],
                k,
            )

            roots = lift_square_roots(
                case["N"],
                k,
            )

            for r in roots:
                result = analyze_root(
                    case,
                    r,
                    k,
                )

                pair = (
                    result["v2_r-p"],
                    result["v2_r+p"],
                )

                split_counts[pair] = (
                    split_counts.get(pair, 0) + 1
                )

        print()
        print(f"k={k}")

        for pair, count in sorted(
            split_counts.items()
        ):
            print(
                f"(v2(r-p), v2(r+p))="
                f"{pair}  "
                f"count={count}"
            )

    print()

    print("============================================================")
    print("TEST 4: SAME TEST FOR q")
    print("============================================================")

    for k in [8, 12, 16]:
        failures = 0
        total = 0

        for case in cases:
            roots = lift_square_roots(
                case["N"],
                k,
            )

            target = v2_capped(
                case["q"] - case["p"],
                k,
            )

            for r in roots:
                result = analyze_root(
                    case,
                    r,
                    k,
                )

                total += 1

                if result["sum_q"] != target:
                    failures += 1

        print(
            f"k={k:2d}  "
            f"checked={total:5d}  "
            f"failures={failures}"
        )

    print()

    print("============================================================")
    print("TEST 5: DETAILED EXAMPLE")
    print("============================================================")

    case = cases[0]

    print(f"N = {case['N']}")
    print(f"p = {case['p']}")
    print(f"q = {case['q']}")

    print(
        f"v2(q-p) = "
        f"{v2(case['q'] - case['p'])}"
    )

    print()

    k = 16

    roots = lift_square_roots(
        case["N"],
        k,
    )

    for index, r in enumerate(roots):
        result = analyze_root(
            case,
            r,
            k,
        )

        print(f"ROOT {index + 1}")
        print(f"r = {r}")

        print(
            f"v2(r-p) = {result['v2_r-p']}"
        )

        print(
            f"v2(r+p) = {result['v2_r+p']}"
        )

        print(
            f"sum      = {result['sum_p']}"
        )

        print(
            f"v2(r-q) = {result['v2_r-q']}"
        )

        print(
            f"v2(r+q) = {result['v2_r+q']}"
        )

        print(
            f"q sum    = {result['sum_q']}"
        )

        print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
