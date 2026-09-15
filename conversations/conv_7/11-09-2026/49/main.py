import math
import random

EXPERIMENT = 67


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

        xp = s + 1 - p
        xq = s + 1 - q

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "xp": xp,
            "xq": xq,
        }


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

    roots = []

    for a in roots_p:
        for b in roots_q:
            roots.append(
                crt_two(
                    a,
                    p,
                    b,
                    q,
                )
            )

    return sorted(roots)


def circular_gaps(roots, n):
    gaps = []

    for i in range(len(roots)):
        a = roots[i]
        b = roots[(i + 1) % len(roots)]

        gaps.append((b - a) % n)

    return gaps


def main():
    random.seed(67001)

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
    print("TEST 1: MINIMUM GAP = 2s+1")
    print("=" * 60)

    failures = 0

    for case in cases:
        roots = crt_roots(case)
        gaps = circular_gaps(
            roots,
            case["N"],
        )

        target = 2 * case["s"] + 1

        if target not in gaps:
            failures += 1

    print(
        f"failures = {failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: OTHER GAP EXPRESSIONS")
    print("=" * 60)

    failures = {
        "2s+1": 0,
        "2p": 0,
        "2q": 0,
        "q-p": 0,
        "2(q-p)": 0,
        "2s+1-(q-p)": 0,
        "2s+1+(q-p)": 0,
    }

    for case in cases:
        roots = crt_roots(case)
        gaps = circular_gaps(
            roots,
            case["N"],
        )

        p = case["p"]
        q = case["q"]
        s = case["s"]

        expressions = [
            2 * s + 1,
            2 * p,
            2 * q,
            q - p,
            2 * (q - p),
            2 * s + 1 - (q - p),
            2 * s + 1 + (q - p),
        ]

        names = list(failures.keys())

        for name, value in zip(
            names,
            expressions,
        ):
            if value not in gaps:
                failures[name] += 1

    for name, count in failures.items():
        print(
            f"{name:24s} failures = {count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: GAP COMPLEMENT")
    print("=" * 60)

    failures = 0

    for case in cases:
        roots = crt_roots(case)
        gaps = circular_gaps(
            roots,
            case["N"],
        )

        n = case["N"]
        target = 2 * case["s"] + 1

        # There are four gaps whose sum is N.
        #
        # Test whether the sum of the other
        # three gaps has the expected complement.
        other = [
            g for g in gaps
            if g != target
        ]

        if len(other) != 3:
            failures += 1
            continue

        if sum(other) + target != n:
            failures += 1

    print(
        f"failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: DUPLICATE LARGE GAPS")
    print("=" * 60)

    duplicate_failures = 0

    for case in cases:
        roots = crt_roots(case)
        gaps = circular_gaps(
            roots,
            case["N"],
        )

        large = sorted(gaps)

        # Test whether the largest two are equal
        # after sorting.
        if large[-1] != large[-2]:
            duplicate_failures += 1

    print(
        f"largest-gap equality failures = "
        f"{duplicate_failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: GAP GCDs WITH N")
    print("=" * 60)

    counts = {
        "gcd(p)": 0,
        "gcd(q)": 0,
        "gcd(1)": 0,
        "other": 0,
    }

    for case in cases:
        gaps = circular_gaps(
            crt_roots(case),
            case["N"],
        )

        n = case["N"]

        for gap in gaps:
            g = math.gcd(
                gap,
                n,
            )

            if g == case["p"]:
                counts["gcd(p)"] += 1

            elif g == case["q"]:
                counts["gcd(q)"] += 1

            elif g == 1:
                counts["gcd(1)"] += 1

            else:
                counts["other"] += 1

    for name, count in counts.items():
        print(
            f"{name:10s} = {count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: RECONSTRUCT s FROM CRT ROOT GAP")
    print("=" * 60)

    failures = 0

    for case in cases:
        roots = crt_roots(case)
        gaps = circular_gaps(
            roots,
            case["N"],
        )

        recovered_s = (
            min(gaps) - 1
        ) // 2

        if recovered_s != case["s"]:
            failures += 1

    print(
        f"failures = {failures}/{cases_count}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: DOES GAP KNOW p OR q?")
    print("=" * 60)

    hit_p = 0
    hit_q = 0
    hit_both = 0
    hit_neither = 0

    for case in cases:
        gaps = circular_gaps(
            crt_roots(case),
            case["N"],
        )

        g = math.gcd(
            case["N"],
            min(gaps),
        )

        if g == case["p"]:
            hit_p += 1

        elif g == case["q"]:
            hit_q += 1

        elif g == case["N"]:
            hit_both += 1

        else:
            hit_neither += 1

    print(
        f"min-gap gcd = p : {hit_p}"
    )

    print(
        f"min-gap gcd = q : {hit_q}"
    )

    print(
        f"min-gap gcd = N : {hit_both}"
    )

    print(
        f"min-gap other   : {hit_neither}"
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
        roots = crt_roots(case)

        gaps = circular_gaps(
            roots,
            case["N"],
        )

        print()
        print(f"CASE {number}")
        print(f"N = {case['N']}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {case['s']}")
        print(
            "2s+1 =",
            2 * case["s"] + 1,
        )
        print(
            "roots =",
            roots,
        )
        print(
            "gaps =",
            gaps,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
