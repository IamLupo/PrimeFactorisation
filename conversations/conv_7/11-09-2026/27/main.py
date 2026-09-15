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

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "gap": q - p,
            "gap_v2": v2(q - p),
            "y": p + q - 2 * s - 1,
        }


def gap_identity_rhs(case, y):
    s = case["s"]
    d = case["D"]

    return (
        4 * d
        - 4 * s * (1 + y)
        + (2 * s + 1 + y) ** 2
    )


def gap_identity_residual(case, y, g):
    return gap_identity_rhs(case, y) - g * g


def solve_gap_square(case, y):
    g2 = gap_identity_rhs(
        case,
        y,
    )

    if g2 < 0:
        return None

    g = math.isqrt(g2)

    if g * g != g2:
        return None

    return g


def test_exact_identity(cases):
    print("============================================================")
    print("TEST 1: EXACT GAP IDENTITY")
    print("============================================================")

    failures = 0

    for case in cases:
        residual = gap_identity_residual(
            case,
            case["y"],
            case["gap"],
        )

        if residual != 0:
            failures += 1

    print(
        f"identity failures = "
        f"{failures}/{len(cases)}"
    )

    print()


def test_true_y_recovers_gap(cases):
    print("============================================================")
    print("TEST 2: y -> GAP")
    print("============================================================")

    recovered = 0

    for case in cases:
        g = solve_gap_square(
            case,
            case["y"],
        )

        if g == case["gap"]:
            recovered += 1

    print(
        f"exact gap recovered from true y = "
        f"{recovered}/{len(cases)}"
    )

    print()


def test_y_mod_2k(cases, k):
    """
    For every y modulo 2^k, calculate the implied g^2.
    Determine how many y residues produce the true g^2 residue.
    """

    modulus = 1 << k

    unique = 0
    ambiguous = 0
    impossible = 0

    for case in cases:
        true_y = case["y"] % modulus

        true_g2 = (
            case["gap"] * case["gap"]
        ) % modulus

        matching = []

        s = case["s"] % modulus
        d = case["D"] % modulus

        for y in range(modulus):
            g2 = (
                4 * d
                - 4 * s * (1 + y)
                + (2 * s + 1 + y) ** 2
            ) % modulus

            if g2 == true_g2:
                matching.append(y)

        if true_y not in matching:
            impossible += 1
        elif len(matching) == 1:
            unique += 1
        else:
            ambiguous += 1

    return unique, ambiguous, impossible


def possible_gap_valuations(g2, k):
    """
    Given g^2 mod 2^k, determine which valuations t
    are compatible with it.

    For t < ceil(k/2):

        v2(g^2) = 2t.

    If 2t >= k, then g^2 == 0 mod 2^k,
    so the only conclusion is

        v2(g) >= ceil(k/2).

    We represent the latter by k//2.
    """

    if g2 == 0:
        return list(
            range(
                (k + 1) // 2,
                k + 1,
            )
        )

    valuation = v2(g2)

    if valuation >= k:
        return list(
            range(
                (k + 1) // 2,
                k + 1,
            )
        )

    if valuation % 2 != 0:
        return []

    t = valuation // 2

    if 2 * t >= k:
        return list(
            range(
                (k + 1) // 2,
                k + 1,
            )
        )

    return [t]


def test_gap_valuation_via_square(cases, k):
    """
    Test how much information g^2 mod 2^k gives
    about v2(g).
    """

    exact = 0
    ambiguous = 0
    impossible = 0

    for case in cases:
        g = case["gap"]

        g2 = (
            g * g
        ) % (1 << k)

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


def test_y_to_gap_valuation(cases, k):
    """
    At the true y residue, calculate g^2 and determine
    the information available about v2(g).
    """

    exact = 0
    ambiguous = 0

    for case in cases:
        modulus = 1 << k

        y = case["y"] % modulus

        s = case["s"] % modulus
        d = case["D"] % modulus

        g2 = (
            4 * d
            - 4 * s * (1 + y)
            + (2 * s + 1 + y) ** 2
        ) % modulus

        possible = possible_gap_valuations(
            g2,
            k,
        )

        if len(possible) == 1:
            exact += 1
        else:
            ambiguous += 1

    return exact, ambiguous


def print_example(case):
    print(f"N      = {case['N']}")
    print(f"p      = {case['p']}")
    print(f"q      = {case['q']}")
    print(f"s      = {case['s']}")
    print(f"D      = {case['D']}")
    print(f"y      = {case['y']}")
    print(f"gap    = {case['gap']}")
    print(f"v2gap  = {case['gap_v2']}")

    print()

    lhs = 4 * case["D"]

    rhs = gap_identity_rhs(
        case,
        case["y"],
    )

    print(f"4D  = {lhs}")
    print(f"rhs = {rhs}")
    print(f"equal = {lhs == rhs}")


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

    test_exact_identity(cases)

    test_true_y_recovers_gap(cases)

    print("============================================================")
    print("TEST 3: y-RESIDUE AMBIGUITY")
    print("============================================================")

    for k in [4, 6, 8, 10, 12]:
        unique, ambiguous, impossible = (
            test_y_mod_2k(
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
            test_gap_valuation_via_square(
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
        exact, ambiguous = (
            test_y_to_gap_valuation(
                cases,
                k,
            )
        )

        print(
            f"k={k:2d}  "
            f"exact v2(g)={exact:4d}  "
            f"ambiguous={ambiguous:4d}"
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