import math
import random


EXPERIMENT = 56


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


def v2(n):
    if n == 0:
        return 999999

    n = abs(n)

    return (n & -n).bit_length() - 1


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
        d = n - s * s

        a = s - p
        b = q - s
        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "a": a,
            "b": b,
            "y": y,
        }


def M_value(case, y):
    """
    M(y) = ab = s(y+1)-D.
    """

    return (
        case["s"] * (y + 1)
        - case["D"]
    )


def discriminant_from_M(case, y):
    """
    If

        b-a = y+1
        ab = M

    then

        (a+b)^2
        = (b-a)^2 + 4ab
        = (y+1)^2 + 4M.
    """

    M = M_value(
        case,
        y,
    )

    return (
        (y + 1) ** 2
        + 4 * M
    )


def recover_ab(case, y):
    """
    Recover a,b from the discriminant.
    """

    delta = discriminant_from_M(
        case,
        y,
    )

    if delta < 0:
        return None

    root = math.isqrt(delta)

    if root * root != delta:
        return None

    d = y + 1

    if (root - d) % 2 != 0:
        return None

    a = (root - d) // 2
    b = (root + d) // 2

    if a < 0:
        return None

    if b <= a:
        return None

    return a, b


def test_exact_identity(cases):
    print("============================================================")
    print("TEST 1: EXACT ab IDENTITY")
    print("============================================================")

    failures = 0

    for case in cases:
        M = M_value(
            case,
            case["y"],
        )

        expected = (
            case["a"]
            * case["b"]
        )

        if M != expected:
            failures += 1

    print(
        f"identity failures = "
        f"{failures}/{len(cases)}"
    )

    print()


def test_exact_recovery(cases):
    print("============================================================")
    print("TEST 2: TRUE y -> EXACT a,b")
    print("============================================================")

    failures = 0

    for case in cases:
        result = recover_ab(
            case,
            case["y"],
        )

        if result is None:
            failures += 1
            continue

        a, b = result

        if (
            a != case["a"]
            or b != case["b"]
        ):
            failures += 1

    print(
        f"recovery failures = "
        f"{failures}/{len(cases)}"
    )

    print()


def test_candidate_y(cases, max_multiplier):
    print("============================================================")
    print("TEST 3: CANDIDATE y FILTERS")
    print("============================================================")

    total_candidates = 0
    nonnegative = 0
    parity_ok = 0
    square = 0
    true_survives = 0

    for case in cases:
        true_y = case["y"]

        limit = min(
            true_y * max_multiplier + 20,
            4 * case["s"],
        )

        for y in range(0, limit + 1):
            total_candidates += 1

            M = M_value(
                case,
                y,
            )

            if M < 0:
                continue

            nonnegative += 1

            delta = (
                (y + 1) ** 2
                + 4 * M
            )

            root = math.isqrt(delta)

            if (
                (root - (y + 1)) % 2
                != 0
            ):
                continue

            parity_ok += 1

            if root * root != delta:
                continue

            square += 1

            if y == true_y:
                true_survives += 1

    print(
        f"total y candidates = "
        f"{total_candidates}"
    )

    print(
        f"M(y)>=0            = "
        f"{nonnegative}"
    )

    print(
        f"parity-valid        = "
        f"{parity_ok}"
    )

    print(
        f"perfect-square      = "
        f"{square}"
    )

    print(
        f"true y survived     = "
        f"{true_survives}/{len(cases)}"
    )

    print()


def test_simple_lower_bound(cases):
    print("============================================================")
    print("TEST 4: LOWER BOUND FROM ab >= 0")
    print("============================================================")

    failures = 0

    gaps = []

    for case in cases:
        # M = s(y+1)-D >= 0
        #
        # therefore y+1 >= D/s.
        #
        # Smallest integer y satisfying this.
        lower = math.ceil(
            case["D"] / case["s"]
        ) - 1

        if case["y"] < lower:
            failures += 1

        gaps.append(
            case["y"] - lower
        )

    print(
        f"violations = "
        f"{failures}/{len(cases)}"
    )

    print(
        f"minimum y-lower = "
        f"{min(gaps)}"
    )

    print(
        f"maximum y-lower = "
        f"{max(gaps)}"
    )

    print(
        f"mean y-lower    = "
        f"{sum(gaps) / len(gaps):.4f}"
    )

    print()


def test_M_factor_structure(cases):
    print("============================================================")
    print("TEST 5: v2 STRUCTURE OF ab")
    print("============================================================")

    differences = []

    for case in cases:
        ab = (
            case["a"]
            * case["b"]
        )

        differences.append(
            v2(ab)
        )

    counts = {}

    for value in differences:
        counts[value] = (
            counts.get(value, 0) + 1
        )

    for value in sorted(counts):
        print(
            f"v2(ab)={value:2d}  "
            f"count={counts[value]:5d}"
        )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 6: EXAMPLES")
    print("============================================================")

    for case in cases[:10]:
        M = M_value(
            case,
            case["y"],
        )

        delta = discriminant_from_M(
            case,
            case["y"],
        )

        print()
        print(
            f"N={case['N']}"
        )

        print(
            f"p={case['p']} "
            f"q={case['q']} "
            f"s={case['s']}"
        )

        print(
            f"a={case['a']} "
            f"b={case['b']}"
        )

        print(
            f"y={case['y']}"
        )

        print(
            f"D={case['D']}"
        )

        print(
            f"M=ab={M}"
        )

        print(
            f"actual ab="
            f"{case['a'] * case['b']}"
        )

        print(
            f"delta={delta}"
        )

        print(
            f"(q-p)^2="
            f"{(case['q'] - case['p']) ** 2}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(565656)

    CASE_COUNT = 300
    PRIME_BITS = 16

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = [
        generate_case(PRIME_BITS)
        for _ in range(CASE_COUNT)
    ]

    test_exact_identity(cases)
    test_exact_recovery(cases)

    test_candidate_y(
        cases,
        2,
    )

    test_simple_lower_bound(cases)
    test_M_factor_structure(cases)
    print_examples(cases)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
