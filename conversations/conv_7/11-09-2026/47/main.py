import math
import random

EXPERIMENT = 65


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


def mod_inverse(a, m):
    return pow(a, -1, m)


def crt_two(a, p, b, q):
    """
    Solve:
        x = a (mod p)
        x = b (mod q)

    with p,q coprime.
    """

    # x = a + p*t
    # a + p*t = b (mod q)
    #
    # p*t = b-a (mod q)

    t = (
        (b - a)
        * mod_inverse(p, q)
    ) % q

    return a + p * t


def modular_roots_prime(case, prime):
    """
    Brute-force only over the two prime roots.
    This is deliberately used for verification and does
    not assume the quadratic formula.
    """

    roots = []

    for x in range(prime):
        if P0(case, x) % prime == 0:
            roots.append(x)

    return roots


def polynomial_root_count(case, prime):
    return len(
        modular_roots_prime(
            case,
            prime,
        )
    )


def normalize_difference(a, b, n):
    d = (a - b) % n

    return min(
        d,
        n - d,
    )


def main():
    random.seed(65001)

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
    print("TEST 1: TRUE LOCAL ROOTS")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        xp = case["xp"]
        xq = case["xq"]

        if P0(case, xp) % p != 0:
            failures += 1

        if P0(case, xq) % q != 0:
            failures += 1

    print(
        f"local root failures = {failures}/{cases_count * 2}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: NUMBER OF ROOTS MOD p AND MOD q")
    print("=" * 60)

    root_count_failures = 0

    counts = {}

    for case in cases:
        p = case["p"]
        q = case["q"]

        rp = polynomial_root_count(
            case,
            p,
        )

        rq = polynomial_root_count(
            case,
            q,
        )

        counts[(rp, rq)] = (
            counts.get((rp, rq), 0) + 1
        )

        if rp != 2 or rq != 2:
            root_count_failures += 1

    for key in sorted(counts):
        print(
            f"roots mod p = {key[0]}, "
            f"roots mod q = {key[1]} : "
            f"{counts[key]}"
        )

    print(
        f"unexpected root counts = "
        f"{root_count_failures}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: SECOND LOCAL ROOT")
    print("=" * 60)

    second_root_failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]

        roots_p = modular_roots_prime(
            case,
            p,
        )

        roots_q = modular_roots_prime(
            case,
            q,
        )

        xp_mod_p = case["xp"] % p
        xq_mod_q = case["xq"] % q

        if xp_mod_p not in roots_p:
            second_root_failures += 1

        if xq_mod_q not in roots_q:
            second_root_failures += 1

        # The other quadratic root is 1-r.
        other_p = (1 - xp_mod_p) % p
        other_q = (1 - xq_mod_q) % q

        if other_p not in roots_p:
            second_root_failures += 1

        if other_q not in roots_q:
            second_root_failures += 1

    print(
        f"second-root identity failures = "
        f"{second_root_failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: FOUR CRT ROOTS MOD N")
    print("=" * 60)

    crt_count_failures = 0
    unique_root_failures = 0

    all_spacing = []

    for case in cases:
        p = case["p"]
        q = case["q"]
        n = case["N"]

        roots_p = modular_roots_prime(
            case,
            p,
        )

        roots_q = modular_roots_prime(
            case,
            q,
        )

        roots_n = []

        for a in roots_p:
            for b in roots_q:
                x = crt_two(
                    a,
                    p,
                    b,
                    q,
                )

                x %= n

                if x not in roots_n:
                    roots_n.append(x)

        roots_n.sort()

        if len(roots_n) != 4:
            crt_count_failures += 1
            continue

        # Verify each CRT root directly.
        for x in roots_n:
            if P0(case, x) % n != 0:
                unique_root_failures += 1

        # Circular spacings.
        for i in range(4):
            a = roots_n[i]
            b = roots_n[(i + 1) % 4]

            if i == 3:
                gap = n - a + b
            else:
                gap = b - a

            all_spacing.append(gap)

    print(
        f"CRT root-count failures = "
        f"{crt_count_failures}"
    )

    print(
        f"direct root verification failures = "
        f"{unique_root_failures}"
    )

    print(
        f"total circular spacings = "
        f"{len(all_spacing)}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: RELATION OF CRT ROOTS TO x_p / x_q")
    print("=" * 60)

    xp_hits = 0
    xq_hits = 0

    distance_values = []

    for case in cases:
        p = case["p"]
        q = case["q"]
        n = case["N"]

        roots_p = modular_roots_prime(
            case,
            p,
        )

        roots_q = modular_roots_prime(
            case,
            q,
        )

        roots_n = []

        for a in roots_p:
            for b in roots_q:
                roots_n.append(
                    crt_two(
                        a,
                        p,
                        b,
                        q,
                    ) % n
                )

        roots_n = sorted(set(roots_n))

        xp = case["xp"] % n
        xq = case["xq"] % n

        if xp in roots_n:
            xp_hits += 1

        if xq in roots_n:
            xq_hits += 1

        for root in roots_n:
            distance_values.append(
                normalize_difference(
                    root,
                    xp,
                    n,
                )
            )

    print(
        f"x_p is root mod N = "
        f"{xp_hits}/{cases_count}"
    )

    print(
        f"x_q is root mod N = "
        f"{xq_hits}/{cases_count}"
    )

    print(
        "minimum circular distance "
        "from CRT root to x_p =",
        min(distance_values),
    )

    print(
        "median circular distance "
        "from CRT root to x_p =",
        sorted(distance_values)[
            len(distance_values) // 2
        ],
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: ROOT DIFFERENCE gcds")
    print("=" * 60)

    recovered_p = 0
    recovered_q = 0
    trivial = 0
    other = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        n = case["N"]

        roots_p = modular_roots_prime(
            case,
            p,
        )

        roots_q = modular_roots_prime(
            case,
            q,
        )

        roots_n = []

        for a in roots_p:
            for b in roots_q:
                roots_n.append(
                    crt_two(
                        a,
                        p,
                        b,
                        q,
                    ) % n
                )

        roots_n = sorted(set(roots_n))

        # Examine all pairwise root differences.
        case_gcds = []

        for i in range(len(roots_n)):
            for j in range(i + 1, len(roots_n)):
                difference = abs(
                    roots_n[i] - roots_n[j]
                )

                g = math.gcd(
                    difference,
                    n,
                )

                case_gcds.append(g)

                if g == p:
                    recovered_p += 1
                elif g == q:
                    recovered_q += 1
                elif g == 1 or g == n:
                    trivial += 1
                else:
                    other += 1

    print(
        f"pairwise gcd = p : {recovered_p}"
    )

    print(
        f"pairwise gcd = q : {recovered_q}"
    )

    print(
        f"pairwise trivial : {trivial}"
    )

    print(
        f"pairwise other   : {other}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: ROOT SUM RELATIONS")
    print("=" * 60)

    failures = {
        "local_sum_p": 0,
        "local_sum_q": 0,
        "crt_sum": 0,
        "crt_product": 0,
    }

    for case in cases:
        p = case["p"]
        q = case["q"]
        n = case["N"]

        roots_p = modular_roots_prime(
            case,
            p,
        )

        roots_q = modular_roots_prime(
            case,
            q,
        )

        if len(roots_p) == 2:
            if (
                (roots_p[0] + roots_p[1]) % p
                != 1 % p
            ):
                failures["local_sum_p"] += 1

        if len(roots_q) == 2:
            if (
                (roots_q[0] + roots_q[1]) % q
                != 1 % q
            ):
                failures["local_sum_q"] += 1

        roots_n = []

        for a in roots_p:
            for b in roots_q:
                roots_n.append(
                    crt_two(
                        a,
                        p,
                        b,
                        q,
                    ) % n
                )

        # For x^2-x+C, the four CRT roots
        # are generated by selecting one local
        # root from each prime.
        #
        # Check whether the sum of all four
        # roots is congruent to 2 mod N.
        if sum(roots_n) % n != 2 % n:
            failures["crt_sum"] += 1

        # Product of all four roots.
        product = 1

        for root in roots_n:
            product = (
                product * root
            ) % n

        # We deliberately report the quantity,
        # rather than assuming a target identity.
        # This lets us inspect whether it simplifies.
        if product == 0:
            failures["crt_product"] += 1

    for name, value in failures.items():
        print(
            f"{name:16s} = {value}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 8
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 8: EXPLICIT CASES")
    print("=" * 60)

    for number, case in enumerate(
        cases[:10],
        start=1,
    ):
        p = case["p"]
        q = case["q"]
        n = case["N"]

        roots_p = modular_roots_prime(
            case,
            p,
        )

        roots_q = modular_roots_prime(
            case,
            q,
        )

        roots_n = []

        for a in roots_p:
            for b in roots_q:
                roots_n.append(
                    crt_two(
                        a,
                        p,
                        b,
                        q,
                    ) % n
                )

        roots_n.sort()

        print()
        print(f"CASE {number}")
        print(f"N       = {n}")
        print(f"p       = {p}")
        print(f"q       = {q}")
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
            "P0(x_p) mod N =",
            P0(case, case["xp"]) % n,
        )

        print(
            "P0(x_q) mod N =",
            P0(case, case["xq"]) % n,
        )

        print(
            "root spacings =",
            [
                (
                    roots_n[(i + 1) % 4]
                    - roots_n[i]
                ) % n
                for i in range(4)
            ],
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
