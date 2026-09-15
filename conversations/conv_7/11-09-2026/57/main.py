import math
import random

EXPERIMENT = 75


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


def get_gap_and_multipliers(case):
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

    # Depending on orientation:
    # either g1=p,g2=q or g1=q,g2=p.
    if g1 == p:
        multiplier_p = r
        multiplier_q = t
    else:
        multiplier_p = t
        multiplier_q = r

    return (
        A,
        multiplier_p,
        multiplier_q,
    )


def continued_fraction(num, den):
    terms = []

    while den != 0:
        q = num // den
        terms.append(q)

        num, den = den, num - q * den

    return terms


def convergents(cf):
    result = []

    p_minus_2 = 0
    p_minus_1 = 1

    q_minus_2 = 1
    q_minus_1 = 0

    for a in cf:
        p = a * p_minus_1 + p_minus_2
        q = a * q_minus_1 + q_minus_2

        result.append((p, q))

        p_minus_2, p_minus_1 = p_minus_1, p
        q_minus_2, q_minus_1 = q_minus_1, q

    return result


def approximation_error(
    a,
    b,
    c,
    d,
):
    """
    |a/b - c/d| as an exact rational.
    Returns numerator, denominator.
    """
    numerator = abs(a * d - c * b)
    denominator = b * d

    return numerator, denominator


def main():
    random.seed(75001)

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
    print("TEST 1: EXACT LINEAR RELATION")
    print("=" * 60)

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]

        _, rp, rq = get_gap_and_multipliers(case)

        if (
            rq * q
            - rp * p
            != 2 * s + 1
        ):
            failures += 1

    print(
        f"tq-rp = 2s+1 failures = "
        f"{failures}"
    )
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
        s = case["s"]

        _, rp, rq = get_gap_and_multipliers(case)

        num, den = approximation_error(
            p,
            q,
            rq,
            rp,
        )

        # Error = (2s+1)/(q*r)
        assert num == 2 * s + 1
        assert den == q * rp

        errors.append(
            num / den
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
        f"mean error = "
        f"{sum(errors) / len(errors):.12e}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 3
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 3: IS t/r A CONTINUED-FRACTION CONVERGENT?")
    print("=" * 60)

    exact_convergent_hits = 0
    semiconvergent_possible = 0

    for case in cases:
        p = case["p"]
        q = case["q"]

        _, rp, rq = get_gap_and_multipliers(case)

        cf = continued_fraction(
            p,
            q,
        )

        conv = convergents(cf)

        if (rq, rp) in conv:
            exact_convergent_hits += 1

        # A loose semiconvergent test:
        #
        # q*r? We check whether denominator rp
        # occurs among intermediate Euclidean
        # denominators.

        denom_set = {
            denominator
            for _, denominator in conv
        }

        if rp in denom_set:
            semiconvergent_possible += 1

    print(
        f"exact convergent hits = "
        f"{exact_convergent_hits}/{cases_count}"
    )

    print(
        f"same denominator as convergent = "
        f"{semiconvergent_possible}/{cases_count}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 4
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 4: CONTINUED-FRACTION QUALITY")
    print("=" * 60)

    quality = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        _, rp, rq = get_gap_and_multipliers(case)

        numerator = abs(
            p * rp
            - q * rq
        )

        # For p/q vs t/r the numerator is |p*r-q*t|.
        # We know from tq-rp=m:
        #
        # p*r - q*t = -m.
        #
        quality_value = (
            numerator
            * q
            * rp
        )

        quality.append(
            quality_value
        )

    print(
        "For a true convergent a/b to p/q, "
        "the quantity"
    )

    print(
        "|p*b-q*a|*q*b"
    )

    print(
        "is typically very small."
    )

    print(
        f"minimum = {min(quality)}"
    )

    print(
        f"maximum = {max(quality)}"
    )

    print(
        f"mean    = "
        f"{sum(quality) / len(quality):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 5
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 5: EUCLIDEAN QUOTIENT STRUCTURE")
    print("=" * 60)

    same_cf_length = 0
    common_prefix = []

    for case in cases:
        p = case["p"]
        q = case["q"]

        _, rp, rq = get_gap_and_multipliers(case)

        cf_factor = continued_fraction(
            p,
            q,
        )

        cf_multiplier = continued_fraction(
            rq,
            rp,
        )

        if len(cf_factor) == len(
            cf_multiplier
        ):
            same_cf_length += 1

        prefix = 0

        for a, b in zip(
            cf_factor,
            cf_multiplier,
        ):
            if a != b:
                break

            prefix += 1

        common_prefix.append(prefix)

    print(
        f"same CF length = "
        f"{same_cf_length}/{cases_count}"
    )

    print(
        f"minimum common prefix = "
        f"{min(common_prefix)}"
    )

    print(
        f"maximum common prefix = "
        f"{max(common_prefix)}"
    )

    print(
        f"mean common prefix = "
        f"{sum(common_prefix) / len(common_prefix):.6f}"
    )

    print()

    # ---------------------------------------------------------
    # TEST 6
    # ---------------------------------------------------------

    print("=" * 60)
    print("TEST 6: MULTIPLIER PRODUCT VS FACTOR GAP")
    print("=" * 60)

    matches = {
        "r*t = (q-p)^2": 0,
        "r*t = (q-p)^2/2": 0,
        "r*t = y^2": 0,
        "r*t = y*(q-p)": 0,
        "r*t = s^2": 0,
        "r*t = D": 0,
    }

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]
        n = case["N"]

        y = p + q - 2 * s - 1
        D = n - s * s

        _, rp, rq = get_gap_and_multipliers(case)

        product = rp * rq
        d = q - p

        tests = {
            "r*t = (q-p)^2": product == d * d,
            "r*t = (q-p)^2/2": (
                2 * product == d * d
            ),
            "r*t = y^2": product == y * y,
            "r*t = y*(q-p)": (
                product == y * d
            ),
            "r*t = s^2": product == s * s,
            "r*t = D": product == D,
        }

        for name, value in tests.items():
            if value:
                matches[name] += 1

    for name, count in matches.items():
        print(
            f"{name:24s} = {count}"
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
        s = case["s"]

        A, rp, rq = get_gap_and_multipliers(
            case
        )

        cf_factor = continued_fraction(
            p,
            q,
        )

        cf_multiplier = continued_fraction(
            rq,
            rp,
        )

        print()
        print(f"CASE {number}")
        print(f"p = {p}")
        print(f"q = {q}")
        print(f"s = {s}")
        print(f"A = {A}")
        print(f"r = {rp}")
        print(f"t = {rq}")
        print(
            "2s+1 =",
            2 * s + 1,
        )

        print(
            "t*q-r*p =",
            rq * q - rp * p,
        )

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
