import math
import random


EXPERIMENT = 43


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


def P(n, s, d, x, y):
    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def centered_residue(value, modulus):
    value %= modulus

    half = modulus // 2

    if value >= half:
        value -= modulus

    return value


def lift_square_roots(n, k):
    """
    Find all square roots of odd n modulo 2^k.

    For k >= 3, an odd n has roots modulo 2^k
    iff n == 1 mod 8.

    Starting from the four roots modulo 8,
    lift each root one bit at a time.
    """

    if k == 0:
        return [0]

    if k == 1:
        roots = []

        for r in range(2):
            if (r * r - n) % 2 == 0:
                roots.append(r)

        return roots

    if k == 2:
        roots = []

        for r in range(4):
            if (r * r - n) % 4 == 0:
                roots.append(r)

        return roots

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


def fixed_points_from_sqrt(n, s, k):
    modulus = 1 << k

    roots = lift_square_roots(
        n,
        k,
    )

    points = []

    for r in roots:
        x = (s + 1 + r) % modulus
        y = (1 - 2 * x) % modulus

        points.append(
            (
                x,
                y,
                r,
            )
        )

    return sorted(points)


def direct_fixed_points(n, s, d, k):
    """
    Check the fixed-point condition directly:

        partner = 1 - y - x
        partner == x

    We use the fact that a fixed point must satisfy
    x == s (mod 2), so only half the residues need testing.
    """

    modulus = 1 << k
    points = []

    for x in range(s % 2, modulus, 2):
        b = x - s - 1

        y = (
            -(
                x * x
                - x
                + d
                - s
            )
            * pow(b, -1, modulus)
        ) % modulus

        partner = (
            1 - y - x
        ) % modulus

        if partner == x:
            points.append(
                (x, y)
            )

    return sorted(points)


def verify_exact_identity(cases):
    print("============================================================")
    print("TEST 1: EXACT FIXED-POINT IDENTITY")
    print("============================================================")

    failures = 0
    tested = 0

    for case in cases:
        n = case["N"]
        s = case["s"]
        d = case["D"]

        for _ in range(20):
            x = random.randint(
                -2 * s,
                2 * s,
            )

            left = P(
                n,
                s,
                d,
                x,
                1 - 2 * x,
            )

            right = (
                n
                - (x - s - 1) ** 2
            )

            tested += 1

            if left != right:
                failures += 1

    print(f"tested  = {tested}")
    print(f"failures = {failures}")
    print()


def compare_fixed_points(cases):
    print("============================================================")
    print("TEST 2: DIRECT FIXED POINTS VS 2-ADIC SQUARE ROOTS")
    print("============================================================")

    for k in [4, 6, 8, 10, 12, 14, 16]:
        total_cases = len(cases)
        exact = 0
        total_fixed = 0
        total_sqrt = 0

        for case in cases:
            n = case["N"]
            s = case["s"]
            d = case["D"]

            direct = direct_fixed_points(
                n,
                s,
                d,
                k,
            )

            sqrt_points = fixed_points_from_sqrt(
                n,
                s,
                k,
            )

            direct_set = set(direct)

            sqrt_set = set(
                (x, y)
                for x, y, _ in sqrt_points
            )

            total_fixed += len(direct)
            total_sqrt += len(sqrt_points)

            if direct_set == sqrt_set:
                exact += 1

        print(
            f"k={k:2d}  "
            f"matching cases={exact:4d}/{total_cases}  "
            f"direct fixed={total_fixed:5d}  "
            f"sqrt roots={total_sqrt:5d}"
        )

    print()


def classify_mod8(cases):
    print("============================================================")
    print("TEST 3: EXISTENCE VS N MOD 8")
    print("============================================================")

    groups = {}

    for case in cases:
        n8 = case["N"] % 8
        groups.setdefault(
            n8,
            {
                "cases": 0,
                "fixed": 0,
            },
        )

        groups[n8]["cases"] += 1

        roots = lift_square_roots(
            case["N"],
            10,
        )

        groups[n8]["fixed"] += len(roots)

    for n8 in sorted(groups):
        data = groups[n8]

        print(
            f"N mod 8={n8}  "
            f"cases={data['cases']:4d}  "
            f"sqrt roots={data['fixed']:5d}"
        )

    print()


def verify_four_roots(cases):
    print("============================================================")
    print("TEST 4: FOUR ROOT STRUCTURE")
    print("============================================================")

    for k in [6, 8, 10, 12, 14, 16]:
        residue_counts = {}

        for case in cases:
            roots = lift_square_roots(
                case["N"],
                k,
            )

            count = len(roots)

            residue_counts[count] = (
                residue_counts.get(count, 0)
                + 1
            )

        print(
            f"k={k:2d}  "
            + "  ".join(
                f"{count} roots={number}"
                for count, number
                in sorted(residue_counts.items())
            )
        )

    print()


def test_gcd_factor_signal(cases):
    print("============================================================")
    print("TEST 5: DOES A MODULAR SQRT DIRECTLY REVEAL A FACTOR?")
    print("============================================================")

    for k in [8, 12, 16]:
        modulus = 1 << k

        nontrivial = 0
        total_roots = 0

        examples = []

        for case in cases:
            roots = lift_square_roots(
                case["N"],
                k,
            )

            for r in roots:
                total_roots += 1

                x = (
                    case["s"]
                    + 1
                    + r
                ) % modulus

                quantities = [
                    r,
                    r - 1,
                    r + 1,
                    x,
                    x - 1,
                    x + 1,
                ]

                found = False

                for value in quantities:
                    g = math.gcd(
                        case["N"],
                        value,
                    )

                    if (
                        1 < g < case["N"]
                    ):
                        found = True
                        break

                if found:
                    nontrivial += 1

                    if len(examples) < 3:
                        examples.append(
                            (
                                case,
                                r,
                                g,
                            )
                        )

        print(
            f"k={k:2d}  "
            f"roots={total_roots:5d}  "
            f"nontrivial gcd={nontrivial:5d}"
        )

    print()


def detailed_example(case):
    print("============================================================")
    print("TEST 6: DETAILED EXAMPLE")
    print("============================================================")

    k = 16
    modulus = 1 << k

    roots = lift_square_roots(
        case["N"],
        k,
    )

    print(f"N       = {case['N']}")
    print(f"p       = {case['p']}")
    print(f"q       = {case['q']}")
    print(f"s       = {case['s']}")
    print(f"x_p     = {case['x_p']}")
    print(f"x_q     = {case['x_q']}")
    print(f"y       = {case['y']}")
    print(f"N mod 8 = {case['N'] % 8}")
    print()

    print("Square roots of N modulo 2^k:")

    for r in roots:
        x = (
            case["s"]
            + 1
            + r
        ) % modulus

        y = (
            1 - 2 * x
        ) % modulus

        identity = (
            P(
                case["N"],
                case["s"],
                case["D"],
                x,
                y,
            )
            % modulus
        )

        print(
            f"r={r:6d}  "
            f"x={x:6d}  "
            f"y={y:6d}  "
            f"P mod M={identity}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(434343)

    CASE_COUNT = 500
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

    verify_exact_identity(cases)
    compare_fixed_points(cases)
    classify_mod8(cases)
    verify_four_roots(cases)
    test_gcd_factor_signal(cases)

    detailed_example(cases[0])

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
