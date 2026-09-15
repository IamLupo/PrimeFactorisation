import math
import random


EXPERIMENT = 77


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
        x
        for x in roots
        if x != known_left
        and x != known_right
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

    if g1 == p:
        rp = r
        tq = t
    else:
        rp = t
        tq = r

    return rp, tq


def continued_fraction(num, den):
    terms = []

    while den != 0:
        a = num // den
        terms.append(a)

        num, den = den, num - a * den

    return terms


def convergents(cf):
    result = []

    p_m2 = 0
    p_m1 = 1

    q_m2 = 1
    q_m1 = 0

    for a in cf:
        p = a * p_m1 + p_m2
        q = a * q_m1 + q_m2

        result.append(
            (p, q)
        )

        p_m2, p_m1 = p_m1, p
        q_m2, q_m1 = q_m1, q

    return result


def common_prefix(a, b):
    length = 0

    for x, y in zip(a, b):
        if x != y:
            break

        length += 1

    return length


def determinant(a, b, c, d):
    return a * d - b * c


def main():
    random.seed(77001)

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
    print("TEST 1: EXACT DETERMINANT")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]

        r, t = get_multipliers(case)

        d = abs(
            p * r
            - q * t
        )

        if d != 2 * s + 1:
            failures += 1

    print(
        f"determinant failures = {failures}"
    )
    print()

    # ---------------------------------------------------------
    # TEST 2
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 2: FIRST CF DIVERGENCE")
    print("=" * 60)

    divergence_indices = []
    previous_convergent_determinants = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        cf1 = continued_fraction(
            p,
            q,
        )

        cf2 = continued_fraction(
            t,
            r,
        )

        prefix = common_prefix(
            cf1,
            cf2,
        )

        divergence_indices.append(
            prefix
        )

        conv = convergents(cf1)

        if prefix == 0:
            previous_convergent_determinants.append(
                None
            )
            continue

        # The common prefix of length L means
        # the convergent at index L-1 is shared.
        a, b = conv[prefix - 1]

        d = abs(
            determinant(
                p,
                q,
                a,
                b,
            )
        )

        previous_convergent_determinants.append(
            d
        )

    print(
        f"minimum divergence index = "
        f"{min(divergence_indices)}"
    )

    print(
        f"maximum divergence index = "
        f"{max(divergence_indices)}"
    )

    print(
        f"mean divergence index = "
        f"{sum(divergence_indices) / len(divergence_indices):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: SHARED CONVERGENT DENOMINATORS")
    print("=" * 60)

    denominator_ratios = []
    determinant_ratios = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        cf1 = continued_fraction(
            p,
            q,
        )

        cf2 = continued_fraction(
            t,
            r,
        )

        prefix = common_prefix(
            cf1,
            cf2,
        )

        if prefix == 0:
            continue

        conv = convergents(cf1)

        a, b = conv[prefix - 1]

        denominator_ratios.append(
            b / r
        )

        det = abs(
            p * b
            - q * a
        )

        determinant_ratios.append(
            det
            / (2 * case["s"] + 1)
        )

    print(
        f"usable cases = "
        f"{len(denominator_ratios)}"
    )

    print(
        f"b/r minimum = "
        f"{min(denominator_ratios):.8f}"
    )

    print(
        f"b/r maximum = "
        f"{max(denominator_ratios):.8f}"
    )

    print(
        f"b/r mean = "
        f"{sum(denominator_ratios) / len(denominator_ratios):.8f}"
    )

    print()

    print(
        "determinant(convergent,p/q) / (2s+1):"
    )

    print(
        f"minimum = "
        f"{min(determinant_ratios):.8f}"
    )

    print(
        f"maximum = "
        f"{max(determinant_ratios):.8f}"
    )

    print(
        f"mean = "
        f"{sum(determinant_ratios) / len(determinant_ratios):.8f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: WHERE DOES t/r LIE?")
    print("=" * 60)

    below = 0
    above = 0
    exact = 0

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        cf = continued_fraction(
            p,
            q,
        )

        conv = convergents(cf)

        # Find the last convergent below/above t/r.
        value = t / r

        lower = None
        upper = None

        for a, b in conv:
            v = a / b

            if v < value:
                lower = (a, b)
            elif v > value and upper is None:
                upper = (a, b)

        if lower is None:
            continue

        la, lb = lower

        if (
            determinant(
                p,
                q,
                t,
                r,
            )
            == 0
        ):
            exact += 1
        elif (
            t * q
            < p * r
        ):
            below += 1
        else:
            above += 1

    print(
        f"t/r below p/q = {below}"
    )

    print(
        f"t/r above p/q = {above}"
    )

    print(
        f"exact          = {exact}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: SEMICONVERGENT TEST")
    print("=" * 60)

    semiconvergent_hits = 0

    for case in cases:
        p = case["p"]
        q = case["q"]

        r, t = get_multipliers(case)

        cf = continued_fraction(
            p,
            q,
        )

        # Generate all intermediate convergents.
        #
        # h_i,k_i are standard convergents.
        # Between consecutive convergents:
        #
        # (h_{i-2}+m h_{i-1}) /
        # (k_{i-2}+m k_{i-1})
        #
        # for 1 <= m < a_i.

        h_m2, h_m1 = 0, 1
        k_m2, k_m1 = 1, 0

        found = False

        for a_i in cf:
            for m in range(1, a_i):
                num = (
                    h_m2
                    + m * h_m1
                )

                den = (
                    k_m2
                    + m * k_m1
                )

                if num == t and den == r:
                    found = True
                    break

            if found:
                break

            h = (
                a_i * h_m1
                + h_m2
            )

            k = (
                a_i * k_m1
                + k_m2
            )

            h_m2, h_m1 = h_m1, h
            k_m2, k_m1 = k_m1, k

        if found:
            semiconvergent_hits += 1

    print(
        f"semiconvergent hits = "
        f"{semiconvergent_hits}/{cases_count}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: DETERMINANT VS CONVERGENT")
    print("=" * 60)

    ratios = []

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]

        r, t = get_multipliers(case)

        cf = continued_fraction(
            p,
            q,
        )

        conv = convergents(cf)

        prefix = common_prefix(
            cf,
            continued_fraction(
                t,
                r,
            ),
        )

        if prefix == 0:
            continue

        a, b = conv[prefix - 1]

        d = abs(
            p * b
            - q * a
        )

        ratios.append(
            d * r
            / (
                2 * s + 1
            )
        )

    print(
        f"minimum = {min(ratios):.8f}"
    )

    print(
        f"maximum = {max(ratios):.8f}"
    )

    print(
        f"mean = "
        f"{sum(ratios) / len(ratios):.8f}"
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

        cf1 = continued_fraction(
            p,
            q,
        )

        cf2 = continued_fraction(
            t,
            r,
        )

        prefix = common_prefix(
            cf1,
            cf2,
        )

        print()
        print(f"CASE {number}")
        print(f"p/q = {p}/{q}")
        print(f"t/r = {t}/{r}")
        print(
            "determinant =",
            abs(p * r - q * t),
        )
        print(
            "2s+1 =",
            2 * case["s"] + 1,
        )
        print(
            "prefix =",
            prefix,
        )
        print(
            "CF(p/q) =",
            cf1,
        )
        print(
            "CF(t/r) =",
            cf2,
        )

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
