import math
import random


EXPERIMENT = 78


def is_prime(n):
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3
    while d * d <= n:
        if n % d == 0:
            return False

        d += 2

    return True


def random_prime(bits):
    low = 1 << (bits - 1)
    high = (1 << bits) - 1

    while True:
        n = random.randint(low, high)
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

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
        }


def crt_two(a, p, b, q):
    t = (
        (b - a)
        * pow(p, -1, q)
    ) % q

    return (a + p * t) % (p * q)


def get_multipliers(case):
    p = case["p"]
    q = case["q"]
    n = case["N"]
    s = case["s"]

    xp = (s + 1) % p
    xq = (s + 1) % q

    roots_p = [
        xp,
        (1 - xp) % p,
    ]

    roots_q = [
        xq,
        (1 - xq) % q,
    ]

    roots = []

    for a in roots_p:
        for b in roots_q:
            roots.append(
                crt_two(a, p, b, q)
            )

    roots = sorted(set(roots))

    known_left = s + 1
    known_right = n - s

    unknown = [
        value
        for value in roots
        if value != known_left
        and value != known_right
    ]

    unknown.sort()

    A = unknown[0] - known_left

    g1 = math.gcd(A, n)
    g2 = math.gcd(
        A + 2 * s + 1,
        n,
    )

    r1 = A // g1
    r2 = (
        A + 2 * s + 1
    ) // g2

    if g1 == p:
        r = r1
        t = r2
    else:
        r = r2
        t = r1

    return r, t


def continued_fraction(num, den):
    result = []

    while den != 0:
        a = num // den
        result.append(a)

        num, den = den, num - a * den

    return result


def common_prefix(a, b):
    length = 0

    for x, y in zip(a, b):
        if x != y:
            break

        length += 1

    return length


def nearest_numerator(p, q, r):
    numerator = p * r

    # Exact integer rounding.
    a = (2 * numerator + q) // (2 * q)

    return a


def fraction_error(p, q, a, b):
    return abs(
        p * b - a * q
    ) / (q * b)


def main():
    random.seed(78001)

    bits = 18
    cases_count = 500

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"cases      = {cases_count}")
    print(f"prime bits = {bits}")
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: SAME-DENOMINATOR APPROXIMATION")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        a = nearest_numerator(
            p,
            q,
            r,
        )

        # Best same-denominator approximation.
        if abs(
            p * r - a * q
        ) > q // 2:
            failures += 1

        # t/r still has the exact CRT relation.
        if abs(
            p * r - t * q
        ) != 2 * case["s"] + 1:
            failures += 1

    print(
        f"construction failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: CF PREFIX COMPARISON")
    print("=" * 60)

    true_prefixes = []
    nearest_prefixes = []
    prefix_differences = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        a = nearest_numerator(
            p,
            q,
            r,
        )

        cf_factor = continued_fraction(
            p,
            q,
        )

        cf_true = continued_fraction(
            t,
            r,
        )

        cf_nearest = continued_fraction(
            a,
            r,
        )

        true_prefix = common_prefix(
            cf_factor,
            cf_true,
        )

        nearest_prefix = common_prefix(
            cf_factor,
            cf_nearest,
        )

        true_prefixes.append(
            true_prefix
        )

        nearest_prefixes.append(
            nearest_prefix
        )

        prefix_differences.append(
            true_prefix - nearest_prefix
        )

    print(
        "TRUE t/r:"
    )

    print(
        f"  minimum = "
        f"{min(true_prefixes)}"
    )

    print(
        f"  maximum = "
        f"{max(true_prefixes)}"
    )

    print(
        f"  mean    = "
        f"{sum(true_prefixes) / len(true_prefixes):.6f}"
    )

    print(
        f"  median  = "
        f"{sorted(true_prefixes)[len(true_prefixes)//2]}"
    )

    print()

    print(
        "NEAREST a/r:"
    )

    print(
        f"  minimum = "
        f"{min(nearest_prefixes)}"
    )

    print(
        f"  maximum = "
        f"{max(nearest_prefixes)}"
    )

    print(
        f"  mean    = "
        f"{sum(nearest_prefixes) / len(nearest_prefixes):.6f}"
    )

    print(
        f"  median  = "
        f"{sorted(nearest_prefixes)[len(nearest_prefixes)//2]}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: PREFIX ADVANTAGE OVER NEAREST RATIONAL")
    print("=" * 60)

    print(
        f"minimum difference = "
        f"{min(prefix_differences)}"
    )

    print(
        f"maximum difference = "
        f"{max(prefix_differences)}"
    )

    print(
        f"mean difference    = "
        f"{sum(prefix_differences) / len(prefix_differences):.6f}"
    )

    positive = sum(
        1
        for x in prefix_differences
        if x > 0
    )

    equal = sum(
        1
        for x in prefix_differences
        if x == 0
    )

    negative = sum(
        1
        for x in prefix_differences
        if x < 0
    )

    print(
        f"true > nearest = "
        f"{positive}/{cases_count}"
    )

    print(
        f"equal          = "
        f"{equal}/{cases_count}"
    )

    print(
        f"true < nearest = "
        f"{negative}/{cases_count}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: APPROXIMATION ERRORS")
    print("=" * 60)

    true_errors = []
    nearest_errors = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        a = nearest_numerator(
            p,
            q,
            r,
        )

        true_errors.append(
            fraction_error(
                p,
                q,
                t,
                r,
            )
        )

        nearest_errors.append(
            fraction_error(
                p,
                q,
                a,
                r,
            )
        )

    print(
        f"true mean error = "
        f"{sum(true_errors) / len(true_errors):.12e}"
    )

    print(
        f"nearest mean error = "
        f"{sum(nearest_errors) / len(nearest_errors):.12e}"
    )

    print()

    print(
        f"true median error = "
        f"{sorted(true_errors)[len(true_errors)//2]:.12e}"
    )

    print(
        f"nearest median error = "
        f"{sorted(nearest_errors)[len(nearest_errors)//2]:.12e}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: DETERMINANT COMPARISON")
    print("=" * 60)

    true_determinants = []
    nearest_determinants = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        a = nearest_numerator(
            p,
            q,
            r,
        )

        true_determinants.append(
            abs(
                p * r
                - q * t
            )
        )

        nearest_determinants.append(
            abs(
                p * r
                - q * a
            )
        )

    print(
        f"true determinant minimum = "
        f"{min(true_determinants)}"
    )

    print(
        f"true determinant maximum = "
        f"{max(true_determinants)}"
    )

    print(
        f"nearest determinant minimum = "
        f"{min(nearest_determinants)}"
    )

    print(
        f"nearest determinant maximum = "
        f"{max(nearest_determinants)}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: CAN DETERMINANT EXPLAIN PREFIX?")
    print("=" * 60)

    advantage_by_case = []

    for i in range(cases_count):
        d_true = true_determinants[i]
        d_near = nearest_determinants[i]

        prefix_true = true_prefixes[i]
        prefix_near = nearest_prefixes[i]

        advantage_by_case.append(
            (
                prefix_true - prefix_near,
                d_true / d_near
                if d_near != 0
                else float("inf"),
            )
        )

    positive_det = [
        ratio
        for diff, ratio in advantage_by_case
        if diff > 0
    ]

    zero_det = [
        ratio
        for diff, ratio in advantage_by_case
        if diff == 0
    ]

    negative_det = [
        ratio
        for diff, ratio in advantage_by_case
        if diff < 0
    ]

    print(
        f"cases true-prefix > nearest = "
        f"{len(positive_det)}"
    )

    if positive_det:
        print(
            f"mean determinant ratio there = "
            f"{sum(positive_det) / len(positive_det):.6f}"
        )

    print()

    print(
        f"cases equal = "
        f"{len(zero_det)}"
    )

    print(
        f"cases true-prefix < nearest = "
        f"{len(negative_det)}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: EXPLICIT CASES")
    print("=" * 60)

    for number, case in enumerate(
        cases[:10],
        start=1,
    ):
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        a = nearest_numerator(
            p,
            q,
            r,
        )

        cf_factor = continued_fraction(
            p,
            q,
        )

        cf_true = continued_fraction(
            t,
            r,
        )

        cf_nearest = continued_fraction(
            a,
            r,
        )

        print()
        print(f"CASE {number}")
        print(
            "p/q =",
            f"{p}/{q}",
        )

        print(
            "t/r =",
            f"{t}/{r}",
        )

        print(
            "nearest =",
            f"{a}/{r}",
        )

        print(
            "true determinant =",
            abs(p * r - q * t),
        )

        print(
            "nearest determinant =",
            abs(p * r - q * a),
        )

        print(
            "true prefix =",
            common_prefix(
                cf_factor,
                cf_true,
            ),
        )

        print(
            "nearest prefix =",
            common_prefix(
                cf_factor,
                cf_nearest,
            ),
        )

        print(
            "CF true =",
            cf_true,
        )

        print(
            "CF nearest =",
            cf_nearest,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
