import math
import random

EXPERIMENT = 66


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
        d = n - s * s

        xp = s + 1 - p
        xq = s + 1 - q

        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "xp": xp,
            "xq": xq,
            "y": y,
        }


def P0(case, x):
    return (
        x * x
        - x
        + case["D"]
        - case["s"]
    )


def crt_two(a, p, b, q):
    t = (
        (b - a)
        * pow(p, -1, q)
    ) % q

    return (a + p * t) % (p * q)


def local_roots(case):
    p = case["p"]
    q = case["q"]

    xp = case["xp"] % p
    xq = case["xq"] % q

    roots_p = sorted({
        xp,
        (1 - xp) % p,
    })

    roots_q = sorted({
        xq,
        (1 - xq) % q,
    })

    return roots_p, roots_q


def crt_roots(case):
    p = case["p"]
    q = case["q"]

    roots_p, roots_q = local_roots(case)

    roots_n = []

    for a in roots_p:
        for b in roots_q:
            roots_n.append(
                crt_two(
                    a,
                    p,
                    b,
                    q,
                )
            )

    return sorted(set(roots_n))


def circular_gap(a, b, n):
    return (b - a) % n


def circular_distance(a, b, n):
    d = (a - b) % n

    return min(
        d,
        n - d,
    )


def mean(values):
    if not values:
        return 0.0

    return sum(values) / len(values)


def median(values):
    if not values:
        return 0.0

    values = sorted(values)
    n = len(values)

    if n % 2:
        return float(values[n // 2])

    return (
        values[n // 2 - 1]
        + values[n // 2]
    ) / 2.0


def main():
    random.seed(66001)

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
    print("TEST 1: LOCAL ROOT STRUCTURE")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        xp = case["xp"]
        xq = case["xq"]

        roots_p, roots_q = local_roots(case)

        if xp % p not in roots_p:
            failures += 1

        if (1 - xp) % p not in roots_p:
            failures += 1

        if xq % q not in roots_q:
            failures += 1

        if (1 - xq) % q not in roots_q:
            failures += 1

        for r in roots_p:
            if P0(case, r) % p != 0:
                failures += 1

        for r in roots_q:
            if P0(case, r) % q != 0:
                failures += 1

    print(
        f"local-root failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: FOUR CRT ROOTS")
    print("=" * 60)

    count_failures = 0
    verification_failures = 0

    for case in cases:
        roots_n = crt_roots(case)
        n = case["N"]

        if len(roots_n) != 4:
            count_failures += 1

        for r in roots_n:
            if P0(case, r) % n != 0:
                verification_failures += 1

    print(
        f"root-count failures = "
        f"{count_failures}"
    )

    print(
        f"direct verification failures = "
        f"{verification_failures}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: IS x_p A CRT ROOT?")
    print("=" * 60)

    xp_hits = 0
    xq_hits = 0

    for case in cases:
        roots = crt_roots(case)
        n = case["N"]

        if case["xp"] % n in roots:
            xp_hits += 1

        if case["xq"] % n in roots:
            xq_hits += 1

    print(
        f"x_p root mod N = "
        f"{xp_hits}/{cases_count}"
    )

    print(
        f"x_q root mod N = "
        f"{xq_hits}/{cases_count}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: ROOT PAIR DIFFERENCES")
    print("=" * 60)

    gcd_p = 0
    gcd_q = 0
    gcd_trivial = 0
    gcd_other = 0

    for case in cases:
        roots = crt_roots(case)
        n = case["N"]
        p = case["p"]
        q = case["q"]

        for i in range(4):
            for j in range(i + 1, 4):

                d = abs(
                    roots[i] - roots[j]
                )

                g = math.gcd(d, n)

                if g == p:
                    gcd_p += 1
                elif g == q:
                    gcd_q += 1
                elif g == 1 or g == n:
                    gcd_trivial += 1
                else:
                    gcd_other += 1

    print(
        f"gcd = p : {gcd_p}"
    )

    print(
        f"gcd = q : {gcd_q}"
    )

    print(
        f"trivial  : {gcd_trivial}"
    )

    print(
        f"other    : {gcd_other}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: CRT ROOT SPACINGS")
    print("=" * 60)

    all_gaps = []

    for case in cases:
        roots = crt_roots(case)
        n = case["N"]

        gaps = []

        for i in range(4):
            gaps.append(
                circular_gap(
                    roots[i],
                    roots[(i + 1) % 4],
                    n,
                )
            )

        all_gaps.extend(gaps)

    print(
        f"gap count = {len(all_gaps)}"
    )

    print(
        f"minimum gap = {min(all_gaps)}"
    )

    print(
        f"maximum gap = {max(all_gaps)}"
    )

    print(
        f"mean gap = {mean(all_gaps):.6f}"
    )

    print(
        f"median gap = {median(all_gaps):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: ROOT SPACING NORMALIZATION")
    print("=" * 60)

    normalized = []

    for case in cases:
        roots = crt_roots(case)
        n = case["N"]

        for i in range(4):
            gap = circular_gap(
                roots[i],
                roots[(i + 1) % 4],
                n,
            )

            normalized.append(
                gap / n
            )

    print(
        f"minimum normalized gap = "
        f"{min(normalized):.10f}"
    )

    print(
        f"maximum normalized gap = "
        f"{max(normalized):.10f}"
    )

    print(
        f"mean normalized gap = "
        f"{mean(normalized):.10f}"
    )

    print(
        f"median normalized gap = "
        f"{median(normalized):.10f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: RELATION TO x_p AND x_q")
    print("=" * 60)

    distances_p = []
    distances_q = []

    for case in cases:
        roots = crt_roots(case)
        n = case["N"]

        for r in roots:
            distances_p.append(
                circular_distance(
                    r,
                    case["xp"],
                    n,
                )
            )

            distances_q.append(
                circular_distance(
                    r,
                    case["xq"],
                    n,
                )
            )

    print(
        "distance to x_p:"
    )

    print(
        f"  min    = "
        f"{min(distances_p)}"
    )

    print(
        f"  median = "
        f"{median(distances_p)}"
    )

    print()

    print(
        "distance to x_q:"
    )

    print(
        f"  min    = "
        f"{min(distances_q)}"
    )

    print(
        f"  median = "
        f"{median(distances_q)}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 8
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 8: ROOT SUM AND PRODUCT")
    print("=" * 60)

    sum_failures = 0

    product_values = []

    for case in cases:
        roots = crt_roots(case)
        n = case["N"]

        root_sum = sum(roots) % n

        # Because each prime sees two roots summing
        # to 1, all four CRT roots should sum to 2
        # modulo each prime, hence modulo N.
        if root_sum != 2 % n:
            sum_failures += 1

        product = 1

        for r in roots:
            product = (
                product * r
            ) % n

        product_values.append(product)

    print(
        f"sum == 2 mod N failures = "
        f"{sum_failures}"
    )

    print(
        "first product values =",
        product_values[:10],
    )

    print()

    # ---------------------------------------------------------
    # TEST 9
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 9: EXPLICIT CASES")
    print("=" * 60)

    for number, case in enumerate(
        cases[:10],
        start=1,
    ):
        roots_p, roots_q = local_roots(case)
        roots_n = crt_roots(case)

        print()
        print(f"CASE {number}")
        print(f"N       = {case['N']}")
        print(f"p       = {case['p']}")
        print(f"q       = {case['q']}")
        print(f"s       = {case['s']}")
        print(f"x_p     = {case['xp']}")
        print(f"x_q     = {case['xq']}")
        print()

        print(
            "roots mod p =",
            roots_p,
        )

        print(
            "roots mod q =",
            roots_q,
        )

        print(
            "roots mod N =",
            roots_n,
        )

        print(
            "sum mod N =",
            sum(roots_n) % case["N"],
        )

        product = 1

        for r in roots_n:
            product = (
                product * r
            ) % case["N"]

        print(
            "product mod N =",
            product,
        )

        print(
            "circular gaps =",
            [
                circular_gap(
                    roots_n[i],
                    roots_n[(i + 1) % 4],
                    case["N"],
                )
                for i in range(4)
            ],
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
