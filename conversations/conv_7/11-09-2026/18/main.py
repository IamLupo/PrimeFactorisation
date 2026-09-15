import math
import random


EXPERIMENT = 40


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


def P_mod(n, s, d, x, y, modulus):
    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    ) % modulus


def solve_y_exists(n, s, d, x, modulus):
    """
    P(x,y) = A + B*y.

    A solution modulo m exists iff
        gcd(B,m) | A.
    """

    a = (
        x * x
        - x
        + d
        - s
    )

    b = x - s - 1

    g = math.gcd(abs(b), modulus)

    return a % g == 0


def centered_residue(x, modulus):
    x %= modulus

    half = modulus // 2

    if x >= half:
        return x - modulus

    return x


def lifted_x_candidates(case, k):
    modulus = 1 << k

    residues = []

    for x in range(modulus):
        if solve_y_exists(
            case["N"],
            case["s"],
            case["D"],
            x,
            modulus,
        ):
            residues.append(x)

    centered = sorted(
        centered_residue(x, modulus)
        for x in residues
    )

    return centered


def plausible_positive(candidates, s):
    return sorted(
        x for x in candidates
        if 1 <= x <= s
    )


def plausible_negative(candidates, s):
    return sorted(
        x for x in candidates
        if -s <= x <= -1
    )


def contains_true_positive(candidates, case):
    return case["x_p"] in candidates


def contains_true_negative(candidates, case):
    return case["x_q"] in candidates


def find_threshold(case, max_k):
    """
    Find the first k where exactly one candidate lies in [1,s].
    """

    for k in range(2, max_k + 1):
        candidates = lifted_x_candidates(case, k)

        positive = plausible_positive(
            candidates,
            case["s"],
        )

        if len(positive) == 1:
            return k, positive[0], candidates

    return None, None, None


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(40040)

    CASE_COUNT = 1500
    PRIME_BITS = 20

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = []

    for _ in range(CASE_COUNT):
        cases.append(
            generate_case(PRIME_BITS)
        )

    print("============================================================")
    print("TEST 1: SIGNED ROOTS")
    print("============================================================")

    for k in range(2, 18):
        unique_positive = 0
        true_positive = 0

        for case in cases:
            candidates = lifted_x_candidates(
                case,
                k,
            )

            positive = plausible_positive(
                candidates,
                case["s"],
            )

            if len(positive) == 1:
                unique_positive += 1

                if positive[0] == case["x_p"]:
                    true_positive += 1

        print(
            f"k={k:2d}  "
            f"unique positive={unique_positive:5d}  "
            f"correct={true_positive:5d}"
        )

    print()

    print("============================================================")
    print("TEST 2: BOTH TRUE ROOTS PRESENT")
    print("============================================================")

    for k in range(2, 18):
        both_present = 0

        for case in cases:
            candidates = lifted_x_candidates(
                case,
                k,
            )

            if (
                case["x_p"] in candidates
                and case["x_q"] in candidates
            ):
                both_present += 1

        print(
            f"k={k:2d}  "
            f"x_p and x_q present="
            f"{both_present:5d}/{CASE_COUNT}"
        )

    print()

    print("============================================================")
    print("TEST 3: FIRST UNIQUE POSITIVE ROOT")
    print("============================================================")

    thresholds = {}

    failures = []

    for case in cases:
        result = find_threshold(
            case,
            20,
        )

        k, root, candidates = result

        if k is None:
            failures.append(case)
        else:
            thresholds[k] = thresholds.get(k, 0) + 1

    for k in sorted(thresholds):
        print(
            f"k={k:2d}  "
            f"cases={thresholds[k]:5d}"
        )

    print(
        f"never unique by k=20: "
        f"{len(failures)}"
    )

    print()

    print("============================================================")
    print("TEST 4: EXAMPLES")
    print("============================================================")

    shown = 0

    for case in cases:
        k, root, candidates = find_threshold(
            case,
            20,
        )

        if k is None:
            continue

        if (
            root == case["x_p"]
            and case["x_q"] in candidates
        ):
            print()
            print(f"k={k}")

            print(f"N     = {case['N']}")
            print(f"p     = {case['p']}")
            print(f"q     = {case['q']}")
            print(f"s     = {case['s']}")
            print(f"x_p   = {case['x_p']}")
            print(f"x_q   = {case['x_q']}")
            print(f"found = {root}")

            print(
                "nearby centered candidates:"
            )

            nearby = [
                x
                for x in candidates
                if abs(x - case["x_p"]) <= 2 * case["s"]
            ]

            print(nearby[:20])

            shown += 1

        if shown >= 5:
            break

    print()

    print("============================================================")
    print("TEST 5: COMPARE THRESHOLD WITH log2(2s)")
    print("============================================================")

    differences = []

    for case in cases:
        k, _, _ = find_threshold(
            case,
            20,
        )

        if k is not None:
            theoretical = math.ceil(
                math.log2(2 * case["s"])
            )

            differences.append(
                k - theoretical
            )

    if differences:
        print(
            f"min difference = "
            f"{min(differences)}"
        )

        print(
            f"max difference = "
            f"{max(differences)}"
        )

        print(
            f"mean difference = "
            f"{sum(differences) / len(differences):.4f}"
        )

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
