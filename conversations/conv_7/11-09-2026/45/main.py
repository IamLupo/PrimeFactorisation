import math
import random

EXPERIMENT = 63


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


def P(case, x, y):
    s = case["s"]
    d = case["D"]

    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def H(case, a, x):
    s = case["s"]
    d = case["D"]

    return (
        (a + 1) * x * x
        + (a * a - a - s - 1) * x
        + (a + 1) * (d - s)
    )


def direct_transform(case, a, x):
    total = 0

    for y in range(x + 1):
        weight = (
            math.comb(x, y)
            * (a ** (x - y))
        )

        total += P(case, x, y) * weight

    return total


def closed_transform(case, a, x):
    return (
        (a + 1) ** (x - 1)
        * H(case, a, x)
    )


def gcd_many(values):
    g = 0

    for value in values:
        g = math.gcd(g, abs(value))

    return g


def factor_recovery(g, case):
    n = case["N"]

    d = math.gcd(g, n)

    if d == 1 or d == n:
        return 0

    return d


def main():
    random.seed(63001)

    bits = 18
    cases_count = 500

    # Different arithmetic-sequence bases.
    bases = [
        1,
        2,
        3,
        4,
        5,
        7,
        11,
    ]

    print(f"START EXPERIMENT {EXPERIMENT}")
    print()
    print("Generating cases...")
    print(f"cases      = {cases_count}")
    print(f"prime bits = {bits}")
    print(f"bases      = {bases}")
    print()

    cases = [
        generate_case(bits)
        for _ in range(cases_count)
    ]

    # ---------------------------------------------------------
    # TEST 1
    # Verify generalized transform identity.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 1: GENERALIZED TRANSFORM IDENTITY")
    print("=" * 60)

    failures = 0

    for case in cases[:100]:
        # Keep x small enough for the direct transform.
        max_x = min(
            case["xp"],
            case["s"] + 1,
            40,
        )

        for x in range(max_x + 1):
            for a in bases:
                direct = direct_transform(
                    case,
                    a,
                    x,
                )

                closed = closed_transform(
                    case,
                    a,
                    x,
                )

                if direct != closed:
                    failures += 1

    print(
        f"identity failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # Evaluate H_a at the true factor coordinate.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: TRUE x_p VALUES")
    print("=" * 60)

    for a in bases:
        values = []

        for case in cases:
            values.append(
                H(case, a, case["xp"])
            )

        gcd_value = gcd_many(
            [
                values[i]
                for i in range(len(values))
            ]
        )

        print(
            f"a={a:2d}  global gcd = {gcd_value}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # Pairwise gcd at the true x_p.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: PAIRWISE H_a GCD AT TRUE x_p")
    print("=" * 60)

    pair_results = {}

    for i in range(len(bases)):
        for j in range(i + 1, len(bases)):
            a = bases[i]
            b = bases[j]

            recovered_p = 0
            recovered_q = 0
            trivial = 0

            for case in cases:
                hp = H(case, a, case["xp"])
                hq = H(case, b, case["xp"])

                g = math.gcd(
                    math.gcd(hp, hq),
                    case["N"],
                )

                if g == case["p"]:
                    recovered_p += 1
                elif g == case["q"]:
                    recovered_q += 1
                elif g == 1 or g == case["N"]:
                    trivial += 1

            pair_results[(a, b)] = (
                recovered_p,
                recovered_q,
                trivial,
            )

            print(
                f"a={a:2d}, b={b:2d}  "
                f"p={recovered_p:3d}  "
                f"q={recovered_q:3d}  "
                f"trivial={trivial:3d}"
            )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # Search x and see whether multi-base gcd detects x_p.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: MULTI-BASE x SEARCH")
    print("=" * 60)

    hits = 0
    false_hits = 0
    xp_hits = 0

    hit_distances = []

    for case in cases:
        n = case["N"]

        found = []

        # Search the physically relevant region.
        for x in range(0, case["s"] + 2):

            values = [
                H(case, a, x)
                for a in bases
            ]

            g = gcd_many(values)
            d = math.gcd(g, n)

            if d != 1 and d != n:
                found.append((x, d))

        for x, d in found:
            if d == case["p"] or d == case["q"]:
                hits += 1

                if x == case["xp"]:
                    xp_hits += 1
            else:
                false_hits += 1

        if found:
            nearest = min(
                abs(x - case["xp"])
                for x, _ in found
            )

            hit_distances.append(nearest)

    print(
        f"cases with factor gcd hit = "
        f"{hits}"
    )

    print(
        f"cases with exact x_p hit = "
        f"{xp_hits}"
    )

    print(
        f"false/non-factor hits = "
        f"{false_hits}"
    )

    if hit_distances:
        print(
            f"minimum distance to x_p = "
            f"{min(hit_distances)}"
        )

        print(
            f"maximum distance to x_p = "
            f"{max(hit_distances)}"
        )

        print(
            f"mean distance to x_p    = "
            f"{sum(hit_distances) / len(hit_distances):.6f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # Compare individual bases against combined bases.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: SINGLE vs MULTI-BASE FACTOR RECOVERY")
    print("=" * 60)

    for a in bases:
        recovered = 0

        for case in cases:
            value = H(
                case,
                a,
                case["xp"],
            )

            d = math.gcd(
                value,
                case["N"],
            )

            if d == case["p"] or d == case["q"]:
                recovered += 1

        print(
            f"a={a:2d}  "
            f"true-x factor gcd = "
            f"{recovered}/{cases_count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # Explicit examples.
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: EXPLICIT CASES")
    print("=" * 60)

    for index, case in enumerate(
        cases[:10],
        start=1,
    ):
        print()
        print(f"CASE {index}")
        print(f"N       = {case['N']}")
        print(f"p       = {case['p']}")
        print(f"q       = {case['q']}")
        print(f"s       = {case['s']}")
        print(f"x_p     = {case['xp']}")
        print(f"x_q     = {case['xq']}")
        print()

        for a in bases:
            value = H(
                case,
                a,
                case["xp"],
            )

            d = math.gcd(
                value,
                case["N"],
            )

            print(
                f"a={a:2d}  "
                f"H_a(x_p)={value}  "
                f"gcd={d}"
            )

        print()

        combined_values = [
            H(
                case,
                a,
                case["xp"],
            )
            for a in bases
        ]

        combined_gcd = math.gcd(
            case["N"],
            gcd_many(combined_values),
        )

        print(
            "combined gcd =",
            combined_gcd,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
