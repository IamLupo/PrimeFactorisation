import math
import random

EXPERIMENT = 72


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

        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "xp": xp,
            "xq": xq,
            "y": y,
        }


def crt_two(a, p, b, q):
    t = (
        (b - a)
        * pow(p, -1, q)
    ) % q

    return (a + p * t) % (p * q)


def crt_roots(case):
    p = case["p"]
    q = case["q"]

    xp = case["xp"] % p
    xq = case["xq"] % q

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
                crt_two(
                    a,
                    p,
                    b,
                    q,
                )
            )

    return sorted(set(roots))


def circular_gaps(roots, n):
    gaps = []

    for i in range(4):
        gaps.append(
            (roots[(i + 1) % 4] - roots[i]) % n
        )

    return gaps


def classify_geometry(case):
    roots = crt_roots(case)
    n = case["N"]
    s = case["s"]

    # The two universal roots.
    r0 = s + 1
    r3 = n - s

    # Remove them.
    unknown = [
        r for r in roots
        if r != r0 and r != r3
    ]

    unknown.sort()

    if len(unknown) != 2:
        return roots, None

    u = unknown[0]
    v = unknown[1]

    A = u - r0
    B = v - u
    C = r3 - v
    D = n - r3 + r0

    return roots, {
        "u": u,
        "v": v,
        "A": A,
        "B": B,
        "C": C,
        "D": D,
    }


def factor_status(g, case):
    if g == case["p"]:
        return "p"

    if g == case["q"]:
        return "q"

    if g == 1:
        return "1"

    if g == case["N"]:
        return "N"

    return "other"


def main():
    random.seed(72001)

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
    print("TEST 1: ROOT SYMMETRY")
    print("=" * 60)

    failures = 0

    for case in cases:
        roots, geometry = classify_geometry(case)

        if geometry is None:
            failures += 1
            continue

        u = geometry["u"]
        v = geometry["v"]
        n = case["N"]

        if (u + v) % n != 1:
            failures += 1

    print(
        f"u+v == 1 mod N failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: REPEATED SIDE GAP")
    print("=" * 60)

    failures = 0

    for case in cases:
        _, geometry = classify_geometry(case)

        if geometry is None:
            failures += 1
            continue

        A = geometry["A"]
        C = geometry["C"]

        if A != C:
            failures += 1

    print(
        f"A == C failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: KNOWN GAP")
    print("=" * 60)

    failures = 0

    for case in cases:
        _, geometry = classify_geometry(case)

        if geometry is None:
            failures += 1
            continue

        D = geometry["D"]

        if D != 2 * case["s"] + 1:
            failures += 1

    print(
        f"D == 2s+1 failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: GAP SUM")
    print("=" * 60)

    failures = 0

    for case in cases:
        _, geometry = classify_geometry(case)

        if geometry is None:
            failures += 1
            continue

        A = geometry["A"]
        B = geometry["B"]
        D = geometry["D"]

        if 2 * A + B + D != case["N"]:
            failures += 1

    print(
        f"2A+B+D=N failures = "
        f"{failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: GAP A gcd WITH N")
    print("=" * 60)

    counts = {
        "p": 0,
        "q": 0,
        "1": 0,
        "N": 0,
        "other": 0,
    }

    A_values = []

    for case in cases:
        _, geometry = classify_geometry(case)

        A = geometry["A"]
        A_values.append(A)

        g = math.gcd(
            A,
            case["N"],
        )

        counts[
            factor_status(g, case)
        ] += 1

    for name in [
        "p",
        "q",
        "1",
        "N",
        "other",
    ]:
        print(
            f"gcd(A,N)={name:5s} : "
            f"{counts[name]}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: GAP A / FACTOR")
    print("=" * 60)

    p_multipliers = []
    q_multipliers = []

    for case in cases:
        _, geometry = classify_geometry(case)

        A = geometry["A"]

        gp = math.gcd(
            A,
            case["N"],
        )

        if gp == case["p"]:
            p_multipliers.append(
                A // case["p"]
            )

        elif gp == case["q"]:
            q_multipliers.append(
                A // case["q"]
            )

    print(
        f"A/p cases = "
        f"{len(p_multipliers)}"
    )

    if p_multipliers:
        print(
            f"A/p minimum = "
            f"{min(p_multipliers)}"
        )

        print(
            f"A/p maximum = "
            f"{max(p_multipliers)}"
        )

        print(
            f"A/p mean = "
            f"{sum(p_multipliers) / len(p_multipliers):.6f}"
        )

    print()

    print(
        f"A/q cases = "
        f"{len(q_multipliers)}"
    )

    if q_multipliers:
        print(
            f"A/q minimum = "
            f"{min(q_multipliers)}"
        )

        print(
            f"A/q maximum = "
            f"{max(q_multipliers)}"
        )

        print(
            f"A/q mean = "
            f"{sum(q_multipliers) / len(q_multipliers):.6f}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 7
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 7: A AGAINST y AND q-p")
    print("=" * 60)

    exact_matches = {
        "A=y": 0,
        "A=(q-p)": 0,
        "A=(q-p)^2": 0,
        "A=y*p": 0,
        "A=y*q": 0,
        "A=(q-p)*p": 0,
        "A=(q-p)*q": 0,
    }

    for case in cases:
        _, geometry = classify_geometry(case)

        A = geometry["A"]
        p = case["p"]
        q = case["q"]
        y = case["y"]
        d = q - p

        tests = {
            "A=y": y,
            "A=(q-p)": d,
            "A=(q-p)^2": d * d,
            "A=y*p": y * p,
            "A=y*q": y * q,
            "A=(q-p)*p": d * p,
            "A=(q-p)*q": d * q,
        }

        for name, value in tests.items():
            if A == value:
                exact_matches[name] += 1

    for name, count in exact_matches.items():
        print(
            f"{name:18s} = {count}"
        )

    print()

    # ---------------------------------------------------------
    # TEST 8
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 8: CAN A BE RECOVERED FROM N,s?")
    print("=" * 60)

    # A is not directly determined by N and s if
    # different semiprime decompositions share N,
    # but for our sampled semiprimes we test simple
    # expressions of N,s.

    failures = {
        "N mod (2s+1)": 0,
        "D": 0,
        "N//(2s+1)": 0,
        "D//(2s+1)": 0,
    }

    for case in cases:
        _, geometry = classify_geometry(case)

        A = geometry["A"]
        n = case["N"]
        s = case["s"]
        d = n - s * s
        m = 2 * s + 1

        candidates = {
            "N mod (2s+1)": n % m,
            "D": d,
            "N//(2s+1)": n // m,
            "D//(2s+1)": d // m,
        }

        for name, value in candidates.items():
            if A == value:
                failures[name] += 1

    for name, count in failures.items():
        print(
            f"{name:20s} = {count}"
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
        roots, geometry = classify_geometry(case)

        print()
        print(f"CASE {number}")
        print(f"N = {case['N']}")
        print(f"p = {case['p']}")
        print(f"q = {case['q']}")
        print(f"s = {case['s']}")
        print(f"y = {case['y']}")
        print()

        print(
            "roots =",
            roots,
        )

        print(
            "unknown roots =",
            [
                geometry["u"],
                geometry["v"],
            ],
        )

        print(
            "gaps =",
            [
                geometry["A"],
                geometry["B"],
                geometry["C"],
                geometry["D"],
            ],
        )

        print(
            "gcd(A,N) =",
            math.gcd(
                geometry["A"],
                case["N"],
            ),
        )

        print(
            "A =",
            geometry["A"],
        )

        print(
            "B =",
            geometry["B"],
        )

        print(
            "A+B =",
            geometry["A"]
            + geometry["B"],
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
