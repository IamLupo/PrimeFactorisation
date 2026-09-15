import math
import random


EXPERIMENT = 48


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

        y = p + q - 2 * s - 1
        gap = q - p

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "y": y,
            "gap": gap,
            "gap_v2": v2(gap),
        }


def gap_square(case, y):
    """
    Exact identity:

        (q - p)^2
        =
        (2s + 1 + y)^2
        - 4s^2
        - 4D

    because

        p + q = 2s + 1 + y
        pq     = s^2 + D.
    """

    s = case["s"]
    d = case["D"]

    total_sum = 2 * s + 1 + y

    return (
        total_sum * total_sum
        - 4 * s * s
        - 4 * d
    )


def gap_square_mod(case, y, k):
    modulus = 1 << k

    s = case["s"] % modulus
    d = case["D"] % modulus
    y = y % modulus

    total_sum = (
        2 * s + 1 + y
    ) % modulus

    return (
        total_sum * total_sum
        - 4 * s * s
        - 4 * d
    ) % modulus


def recover_gap(case, y):
    value = gap_square(
        case,
        y,
    )

    if value < 0:
        return None

    root = math.isqrt(value)

    if root * root != value:
        return None

    return root


def possible_gap_valuations(g2, k):
    """
    Given g^2 mod 2^k, determine the possible v2(g).

    If g^2 != 0 mod 2^k:

        v2(g^2) = 2 v2(g)

    and therefore v2(g) is uniquely determined.

    If g^2 == 0 mod 2^k, only

        v2(g) >= ceil(k/2)

    is known.
    """

    modulus = 1 << k

    g2 %= modulus

    if g2 == 0:
        minimum = (k + 1) // 2

        return list(
            range(
                minimum,
                k + 1,
            )
        )

    valuation = v2(g2)

    if valuation >= k:
        minimum = (k + 1) // 2

        return list(
            range(
                minimum,
                k + 1,
            )
        )

    if valuation % 2 != 0:
        return []

    return [valuation // 2]


def test_identity(cases):
    print("============================================================")
    print("TEST 1: EXACT GAP IDENTITY")
    print("============================================================")

    failures = 0

    for case in cases:
        calculated = gap_square(
            case,
            case["y"],
        )

        expected = case["gap"] ** 2

        if calculated != expected:
            failures += 1

    print(
        f"identity failures = "
        f"{failures}/{len(cases)}"
    )

    print()


def test_recovery(cases):
    print("============================================================")
    print("TEST 2: y -> EXACT GAP")
    print("============================================================")

    recovered = 0
    failed = 0

    for case in cases:
        gap = recover_gap(
            case,
            case["y"],
        )

        if gap == case["gap"]:
            recovered += 1
        else:
            failed += 1

    print(
        f"exact gap recovered = "
        f"{recovered}/{len(cases)}"
    )

    print(
        f"recovery failures   = "
        f"{failed}/{len(cases)}"
    )

    print()


def test_y_residue(cases, k):
    """
    For each case, compare all y residues modulo 2^k
    according to the resulting g^2 residue.

    This is intentionally brute-force only for small k.
    """

    modulus = 1 << k

    unique = 0
    ambiguous = 0
    impossible = 0

    for case in cases:
        true_y = case["y"] % modulus

        true_g2 = (
            case["gap"] ** 2
        ) % modulus

        matches = []

        for y in range(modulus):
            value = gap_square_mod(
                case,
                y,
                k,
            )

            if value == true_g2:
                matches.append(y)

        if true_y not in matches:
            impossible += 1
        elif len(matches) == 1:
            unique += 1
        else:
            ambiguous += 1

    return unique, ambiguous, impossible


def test_g2_valuation(cases, k):
    """
    Determine whether g^2 mod 2^k uniquely determines v2(g).
    """

    exact = 0
    ambiguous = 0
    impossible = 0

    modulus = 1 << k

    for case in cases:
        g = case["gap"]

        g2 = (
            g * g
        ) % modulus

        possible = possible_gap_valuations(
            g2,
            k,
        )

        target = v2(g)

        if target not in possible:
            impossible += 1
        elif len(possible) == 1:
            exact += 1
        else:
            ambiguous += 1

    return exact, ambiguous, impossible


def test_true_y_valuation(cases, k):
    """
    Compute g^2 using only the true y residue and determine
    whether its 2-adic valuation reveals v2(g).
    """

    exact = 0
    ambiguous = 0
    impossible = 0

    for case in cases:
        g2 = gap_square_mod(
            case,
            case["y"],
            k,
        )

        possible = possible_gap_valuations(
            g2,
            k,
        )

        target = case["gap_v2"]

        if target not in possible:
            impossible += 1
        elif len(possible) == 1:
            exact += 1
        else:
            ambiguous += 1

    return exact, ambiguous, impossible


def print_example(case):
    print(f"N       = {case['N']}")
    print(f"p       = {case['p']}")
    print(f"q       = {case['q']}")
    print(f"s       = {case['s']}")
    print(f"D       = {case['D']}")
    print(f"y       = {case['y']}")
    print(f"q-p     = {case['gap']}")
    print(f"v2(q-p) = {case['gap_v2']}")

    total_sum = (
        2 * case["s"]
        + 1
        + case["y"]
    )

    print()
    print(
        f"p+q from y = {total_sum}"
    )

    print(
        f"actual p+q = "
        f"{case['p'] + case['q']}"
    )

    print()

    calculated = gap_square(
        case,
        case["y"],
    )

    expected = case["gap"] ** 2

    print(
        f"computed gap^2 = {calculated}"
    )

    print(
        f"actual gap^2   = {expected}"
    )

    print(
        f"equal          = "
        f"{calculated == expected}"
    )

    print()

    recovered = recover_gap(
        case,
        case["y"],
    )

    print(
        f"recovered gap  = {recovered}"
    )


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(484848)

    CASE_COUNT = 200
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

    test_identity(cases)

    test_recovery(cases)

    print("============================================================")
    print("TEST 3: y-RESIDUE AMBIGUITY")
    print("============================================================")

    for k in [4, 6, 8, 10, 12]:
        unique, ambiguous, impossible = (
            test_y_residue(
                cases,
                k,
            )
        )

        print(
            f"k={k:2d}  "
            f"unique={unique:4d}  "
            f"ambiguous={ambiguous:4d}  "
            f"impossible={impossible:4d}"
        )

    print()

    print("============================================================")
    print("TEST 4: g^2 MOD 2^k -> v2(g)")
    print("============================================================")

    for k in [4, 6, 8, 10, 12, 14]:
        exact, ambiguous, impossible = (
            test_g2_valuation(
                cases,
                k,
            )
        )

        print(
            f"k={k:2d}  "
            f"exact={exact:4d}  "
            f"ambiguous={ambiguous:4d}  "
            f"impossible={impossible:4d}"
        )

    print()

    print("============================================================")
    print("TEST 5: TRUE y -> g^2 INFORMATION")
    print("============================================================")

    for k in [4, 6, 8, 10, 12, 14]:
        exact, ambiguous, impossible = (
            test_true_y_valuation(
                cases,
                k,
            )
        )

        print(
            f"k={k:2d}  "
            f"exact={exact:4d}  "
            f"ambiguous={ambiguous:4d}  "
            f"impossible={impossible:4d}"
        )

    print()

    print("============================================================")
    print("TEST 6: EXAMPLES")
    print("============================================================")

    for case in cases[:5]:
        print()
        print_example(case)

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
