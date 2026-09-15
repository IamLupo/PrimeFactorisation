import math
import random


EXPERIMENT = 42


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


def centered_residue(value, modulus):
    value %= modulus

    half = modulus // 2

    if value >= half:
        value -= modulus

    return value


def solve_y_mod(n, s, d, x, modulus):
    a = (
        x * x
        - x
        + d
        - s
    )

    b = x - s - 1

    if (b & 1) == 0:
        raise ValueError(
            "x does not satisfy the 2-adic parity condition"
        )

    y = (
        -a * pow(b, -1, modulus)
    ) % modulus

    return y


def admissible_x_residues(s, k):
    modulus = 1 << k

    # x == s (mod 2)
    return range(s % 2, modulus, 2)


def y_map(case, k):
    modulus = 1 << k

    result = {}

    for x in admissible_x_residues(
        case["s"],
        k,
    ):
        y = solve_y_mod(
            case["N"],
            case["s"],
            case["D"],
            x,
            modulus,
        )

        result[x] = y

    return result


def involution_partner(x, y, modulus):
    """
    If P(x,y)=0, the other quadratic root is

        x' = 1 - y - x.
    """

    return (1 - y - x) % modulus


def test_involution(case, k):
    modulus = 1 << k
    mapping = y_map(case, k)

    failures = 0

    for x, y in mapping.items():
        partner = involution_partner(
            x,
            y,
            modulus,
        )

        if partner not in mapping:
            failures += 1
            continue

        partner_y = mapping[partner]

        if partner_y != y:
            failures += 1

    return failures


def true_pair_mod(case, k):
    modulus = 1 << k

    return (
        case["x_p"] % modulus,
        case["x_q"] % modulus,
    )


def test_true_pair_relation(cases, k):
    modulus = 1 << k

    correct = 0

    for case in cases:
        xp = case["x_p"] % modulus
        xq = case["x_q"] % modulus
        y = case["y"] % modulus

        if (
            (xp + xq - (1 - y)) % modulus
            == 0
        ):
            correct += 1

    return correct


def count_orbits(case, k):
    modulus = 1 << k
    mapping = y_map(case, k)

    visited = set()
    sizes = []

    for x, y in mapping.items():
        if x in visited:
            continue

        partner = involution_partner(
            x,
            y,
            modulus,
        )

        visited.add(x)
        visited.add(partner)

        if partner == x:
            sizes.append(1)
        else:
            sizes.append(2)

    return sizes


def fixed_points(case, k):
    modulus = 1 << k
    mapping = y_map(case, k)

    result = []

    for x, y in mapping.items():
        partner = involution_partner(
            x,
            y,
            modulus,
        )

        if partner == x:
            result.append((x, y))

    return result


def pair_invariant(case, k, x, y):
    modulus = 1 << k

    partner = involution_partner(
        x,
        y,
        modulus,
    )

    # Several natural symmetric quantities.
    return {
        "x_sum": (x + partner) % modulus,
        "x_product": (x * partner) % modulus,
        "y": y,
        "p_sum": (
            (case["s"] + 1 - x)
            + (case["s"] + 1 - partner)
        ) % modulus,
        "q_sum": (
            (case["s"] + x + y)
            + (case["s"] + partner + y)
        ) % modulus,
    }


def analyze_true_pair_orbit(case, k):
    modulus = 1 << k

    xp = case["x_p"] % modulus
    xq = case["x_q"] % modulus
    y = case["y"] % modulus

    mapping = y_map(case, k)

    partner = involution_partner(
        xp,
        y,
        modulus,
    )

    return {
        "xp": xp,
        "xq": xq,
        "partner": partner,
        "partner_matches_xq": partner == xq,
        "y": y,
        "true_y_map": mapping.get(xp),
    }


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(424242)

    CASE_COUNT = 100
    PRIME_BITS = 14

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
    print("TEST 1: INVOLUTION PROPERTY")
    print("============================================================")

    for k in [6, 8, 10, 12, 14, 16]:
        total_failures = 0

        for case in cases:
            total_failures += test_involution(
                case,
                k,
            )

        print(
            f"k={k:2d}  "
            f"total involution failures="
            f"{total_failures}"
        )

    print()

    print("============================================================")
    print("TEST 2: TRUE ROOT PAIR")
    print("============================================================")

    for k in [4, 6, 8, 10, 12, 14, 16]:
        correct = test_true_pair_relation(
            cases,
            k,
        )

        print(
            f"k={k:2d}  "
            f"x_p+x_q == 1-y mod 2^k: "
            f"{correct}/{CASE_COUNT}"
        )

    print()

    print("============================================================")
    print("TEST 3: ORBIT STRUCTURE")
    print("============================================================")

    for k in [6, 8, 10, 12, 14]:
        orbit1 = 0
        orbit2 = 0
        other = 0

        for case in cases:
            sizes = count_orbits(
                case,
                k,
            )

            for size in sizes:
                if size == 1:
                    orbit1 += 1
                elif size == 2:
                    orbit2 += 1
                else:
                    other += 1

        print(
            f"k={k:2d}  "
            f"fixed-point orbits={orbit1}  "
            f"2-cycles={orbit2}  "
            f"other={other}"
        )

    print()

    print("============================================================")
    print("TEST 4: FIXED POINTS")
    print("============================================================")

    for k in [6, 8, 10, 12, 14, 16]:
        total_fixed = 0

        for case in cases:
            total_fixed += len(
                fixed_points(
                    case,
                    k,
                )
            )

        print(
            f"k={k:2d}  "
            f"total fixed points={total_fixed}"
        )

    print()

    print("============================================================")
    print("TEST 5: TRUE FACTOR PAIR")
    print("============================================================")

    for k in [8, 10, 12, 14, 16]:
        matches = 0

        for case in cases:
            result = analyze_true_pair_orbit(
                case,
                k,
            )

            if result["partner_matches_xq"]:
                matches += 1

        print(
            f"k={k:2d}  "
            f"true partner is x_q: "
            f"{matches}/{CASE_COUNT}"
        )

    print()

    print("============================================================")
    print("TEST 6: DETAILED EXAMPLES")
    print("============================================================")

    for index, case in enumerate(cases[:5]):
        k = 16

        result = analyze_true_pair_orbit(
            case,
            k,
        )

        print()
        print(f"CASE {index + 1}")
        print(f"N        = {case['N']}")
        print(f"p        = {case['p']}")
        print(f"q        = {case['q']}")
        print(f"s        = {case['s']}")
        print(f"x_p      = {case['x_p']}")
        print(f"x_q      = {case['x_q']}")
        print(f"y        = {case['y']}")
        print(f"x_p+x_q  = {case['x_p'] + case['x_q']}")
        print(f"1-y      = {1 - case['y']}")
        print(f"modulus  = {1 << k}")
        print(f"xp mod M = {result['xp']}")
        print(f"xq mod M = {result['xq']}")
        print(f"partner  = {result['partner']}")
        print(
            f"matches x_q = "
            f"{result['partner_matches_xq']}"
        )

    print()

    print("============================================================")
    print("TEST 7: SYMMETRIC INVARIANT")
    print("============================================================")

    for case in cases[:3]:
        k = 12
        modulus = 1 << k

        xp = case["x_p"] % modulus
        yp = case["y"] % modulus

        invariant = pair_invariant(
            case,
            k,
            xp,
            yp,
        )

        print()
        print(f"N = {case['N']}")
        print(f"x_p = {xp}")
        print(f"y   = {yp}")

        for key, value in invariant.items():
            print(
                f"{key:12s}= {value}"
            )

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
