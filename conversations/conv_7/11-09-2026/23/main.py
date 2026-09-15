import math
import random


EXPERIMENT = 44


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

        s = math.isqrt(n)

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": n - s * s,
            "x_p": s + 1 - p,
            "x_q": s + 1 - q,
            "y": p + q - 2 * s - 1,
        }


def lift_square_roots(n, k):
    if k < 1:
        return [0]

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
            for candidate in (
                r,
                r + modulus,
            ):
                if (
                    candidate * candidate - n
                ) % new_modulus == 0:
                    new_roots.append(candidate)

        roots = sorted(set(new_roots))
        modulus = new_modulus

    return roots


def v2(n):
    if n == 0:
        return 999999

    n = abs(n)

    return (n & -n).bit_length() - 1


def distance_mod(a, b, modulus):
    d = (a - b) % modulus

    return min(
        d,
        modulus - d,
    )


def analyse_root(case, r, k):
    modulus = 1 << k

    s = case["s"]
    p = case["p"]
    q = case["q"]
    n = case["N"]

    A = s + 1 + r
    B = s + 1 - r

    quantities = {
        "r": r,
        "A": A % modulus,
        "B": B % modulus,
        "A-s": (A - s) % modulus,
        "B-s": (B - s) % modulus,
        "p": p % modulus,
        "q": q % modulus,
        "p-q": (p - q) % modulus,
        "p+q": (p + q) % modulus,
        "2s+1": (2 * s + 1) % modulus,
    }

    gcds = {
        "gcd(N,r)": math.gcd(n, r),
        "gcd(N,A)": math.gcd(n, A),
        "gcd(N,B)": math.gcd(n, B),
        "gcd(N,A-s)": math.gcd(n, A - s),
        "gcd(N,B-s)": math.gcd(n, B - s),
    }

    distances = {
        "A-p": distance_mod(A, p, modulus),
        "A-q": distance_mod(A, q, modulus),
        "B-p": distance_mod(B, p, modulus),
        "B-q": distance_mod(B, q, modulus),
        "r-(p-q)": distance_mod(r, p - q, modulus),
        "r-(p+q)": distance_mod(r, p + q, modulus),
        "A-(p+q)": distance_mod(A, p + q, modulus),
        "B-(p-q)": distance_mod(B, p - q, modulus),
    }

    valuations = {
        "v2(A-p)": v2(A - p),
        "v2(A-q)": v2(A - q),
        "v2(B-p)": v2(B - p),
        "v2(B-q)": v2(B - q),
        "v2(A-(p+q))": v2(A - (p + q)),
        "v2(B-(p-q))": v2(B - (p - q)),
    }

    return quantities, gcds, distances, valuations


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(444444)

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

    for k in [8, 12, 16]:
        print("============================================================")
        print(f"TEST k={k}")
        print("============================================================")

        modulus = 1 << k

        total_roots = 0

        gcd_hits = {
            "gcd(N,r)": 0,
            "gcd(N,A)": 0,
            "gcd(N,B)": 0,
            "gcd(N,A-s)": 0,
            "gcd(N,B-s)": 0,
        }

        distance_zero = {
            "A-p": 0,
            "A-q": 0,
            "B-p": 0,
            "B-q": 0,
            "r-(p-q)": 0,
            "r-(p+q)": 0,
            "A-(p+q)": 0,
            "B-(p-q)": 0,
        }

        max_v2 = {
            key: 0
            for key in [
                "v2(A-p)",
                "v2(A-q)",
                "v2(B-p)",
                "v2(B-q)",
                "v2(A-(p+q))",
                "v2(B-(p-q))",
            ]
        }

        for case in cases:
            roots = lift_square_roots(
                case["N"],
                k,
            )

            total_roots += len(roots)

            for r in roots:
                _, gcds, distances, valuations = (
                    analyse_root(
                        case,
                        r,
                        k,
                    )
                )

                for key, value in gcds.items():
                    if 1 < value < case["N"]:
                        gcd_hits[key] += 1

                for key, value in distances.items():
                    if value == 0:
                        distance_zero[key] += 1

                for key, value in valuations.items():
                    max_v2[key] = max(
                        max_v2[key],
                        value,
                    )

        print(f"total roots = {total_roots}")

        print()
        print("Nontrivial gcd hits:")

        for key, value in gcd_hits.items():
            print(
                f"{key:15s} = {value}"
            )

        print()
        print("Exact modular matches:")

        for key, value in distance_zero.items():
            print(
                f"{key:15s} = {value}"
            )

        print()
        print("Maximum v2 differences:")

        for key, value in max_v2.items():
            print(
                f"{key:22s} = {value}"
            )

        print()

    print("============================================================")
    print("TEST 2: ROOT SIGNATURES")
    print("============================================================")

    k = 16
    modulus = 1 << k

    signatures = {}

    for case in cases:
        roots = lift_square_roots(
            case["N"],
            k,
        )

        for r in roots:
            _, _, distances, _ = analyse_root(
                case,
                r,
                k,
            )

            nearest = min(
                distances.items(),
                key=lambda item: item[1],
            )

            key = nearest[0]

            signatures[key] = (
                signatures.get(key, 0) + 1
            )

    for key, count in sorted(
        signatures.items(),
        key=lambda item: -item[1],
    ):
        print(
            f"{key:15s} = {count}"
        )

    print()

    print("============================================================")
    print("TEST 3: EXAMPLE")
    print("============================================================")

    case = cases[0]

    print(f"N = {case['N']}")
    print(f"p = {case['p']}")
    print(f"q = {case['q']}")
    print(f"s = {case['s']}")
    print()

    roots = lift_square_roots(
        case["N"],
        16,
    )

    for index, r in enumerate(roots):
        quantities, gcds, distances, valuations = (
            analyse_root(
                case,
                r,
                16,
            )
        )

        print(f"ROOT {index + 1}")
        print(f"r = {r}")

        print()
        print("Distances:")
        for key, value in distances.items():
            print(
                f"  {key:15s} = {value}"
            )

        print()
        print("v2:")
        for key, value in valuations.items():
            print(
                f"  {key:22s} = {value}"
            )

        print()
        print("GCDs:")
        for key, value in gcds.items():
            print(
                f"  {key:15s} = {value}"
            )

        print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
