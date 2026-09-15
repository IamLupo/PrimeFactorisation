import math
import random


EXPERIMENT = 76


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
        r
        for r in roots
        if r != known_left
        and r != known_right
    ]

    unknown.sort()

    A = unknown[0] - known_left

    g1 = math.gcd(A, n)
    g2 = math.gcd(
        A + 2 * s + 1,
        n,
    )

    r = A // g1
    t = (
        A + 2 * s + 1
    ) // g2

    # Determine orientation.
    if g1 == p:
        rp = r
        tq = t
    else:
        rp = t
        tq = r

    return rp, tq


def continued_fraction(num, den):
    result = []

    while den:
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


def random_ratio_like(num_bits):
    low = 1 << (num_bits - 1)
    high = (1 << num_bits) - 1

    a = random.randint(low, high)
    b = random.randint(low, high)

    while math.gcd(a, b) != 1:
        a = random.randint(low, high)
        b = random.randint(low, high)

    return a, b


def main():
    random.seed(76001)

    bits = 18
    cases_count = 500
    random_controls_per_case = 20

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"cases                 = {cases_count}")
    print(f"prime bits            = {bits}")
    print(
        f"random controls/case = "
        f"{random_controls_per_case}"
    )
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    # ---------------------------------------------------------
    # TEST 1
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: CORRECT SIGNED RELATION")
    print("=" * 60)

    failures = 0
    positive = 0
    negative = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]

        r, t = get_multipliers(case)

        value = t * q - r * p
        target = 2 * s + 1

        if value == target:
            positive += 1
        elif value == -target:
            negative += 1
        else:
            failures += 1

    print(f"+ relation = {positive}")
    print(f"- relation = {negative}")
    print(f"failures   = {failures}")
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: APPROXIMATION ERROR")
    print("=" * 60)

    errors = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        numerator = abs(
            p * r - q * t
        )

        denominator = q * r

        errors.append(
            numerator / denominator
        )

    print(
        f"minimum error = "
        f"{min(errors):.12e}"
    )

    print(
        f"maximum error = "
        f"{max(errors):.12e}"
    )

    print(
        f"mean error    = "
        f"{sum(errors) / len(errors):.12e}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: TRUE CF PREFIX LENGTH")
    print("=" * 60)

    true_prefixes = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        cf_factor = continued_fraction(
            p,
            q,
        )

        cf_multiplier = continued_fraction(
            t,
            r,
        )

        true_prefixes.append(
            common_prefix(
                cf_factor,
                cf_multiplier,
            )
        )

    print(
        f"minimum = {min(true_prefixes)}"
    )

    print(
        f"maximum = {max(true_prefixes)}"
    )

    print(
        f"mean    = "
        f"{sum(true_prefixes) / len(true_prefixes):.6f}"
    )

    print(
        f"median  = "
        f"{sorted(true_prefixes)[len(true_prefixes)//2]}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    #
    # Random rational controls with comparable
    # numerator/denominator sizes.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: RANDOM RATIONAL CONTROL")
    print("=" * 60)

    random_prefixes = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        cf_factor = continued_fraction(
            p,
            q,
        )

        for _ in range(
            random_controls_per_case
        ):
            a, b = random_ratio_like(bits)

            cf_random = continued_fraction(
                a,
                b,
            )

            random_prefixes.append(
                common_prefix(
                    cf_factor,
                    cf_random,
                )
            )

    print(
        f"controls = "
        f"{len(random_prefixes)}"
    )

    print(
        f"minimum = "
        f"{min(random_prefixes)}"
    )

    print(
        f"maximum = "
        f"{max(random_prefixes)}"
    )

    print(
        f"mean    = "
        f"{sum(random_prefixes) / len(random_prefixes):.6f}"
    )

    print(
        f"median  = "
        f"{sorted(random_prefixes)[len(random_prefixes)//2]}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: PREFIX ADVANTAGE")
    print("=" * 60)

    true_mean = (
        sum(true_prefixes)
        / len(true_prefixes)
    )

    random_mean = (
        sum(random_prefixes)
        / len(random_prefixes)
    )

    print(
        f"true mean   = {true_mean:.6f}"
    )

    print(
        f"random mean = {random_mean:.6f}"
    )

    print(
        f"difference  = "
        f"{true_mean - random_mean:.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: PREFIX >= k")
    print("=" * 60)

    for k in range(1, 11):
        true_count = sum(
            1
            for x in true_prefixes
            if x >= k
        )

        random_count = sum(
            1
            for x in random_prefixes
            if x >= k
        )

        print(
            f"k={k:2d}  "
            f"true={true_count:4d}/{cases_count}  "
            f"random={random_count:5d}/"
            f"{len(random_prefixes)}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: TRUE PREFIX CASES")
    print("=" * 60)

    for number, case in enumerate(
        cases[:10],
        start=1,
    ):
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        cf_factor = continued_fraction(
            p,
            q,
        )

        cf_multiplier = continued_fraction(
            t,
            r,
        )

        prefix = common_prefix(
            cf_factor,
            cf_multiplier,
        )

        print()
        print(f"CASE {number}")
        print(f"p/q = {p}/{q}")
        print(f"t/r = {t}/{r}")
        print(f"prefix length = {prefix}")
        print(
            "CF(p/q) =",
            cf_factor,
        )
        print(
            "CF(t/r) =",
            cf_multiplier,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
