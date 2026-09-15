import math
import random


EXPERIMENT = 55


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


def v2(n):
    if n == 0:
        return 999999

    n = abs(n)

    return (n & -n).bit_length() - 1


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

        a = s - p
        b = q - s

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "a": a,
            "b": b,
            "radius": max(a, b),
            "sum": a + b,
            "difference": b - a,
        }


def test_exact_identities(cases):
    print("============================================================")
    print("TEST 1: EXACT a,b IDENTITIES")
    print("============================================================")

    failures = 0

    for case in cases:
        s = case["s"]
        d = case["D"]
        a = case["a"]
        b = case["b"]

        # D = s(b-a) - ab
        lhs1 = d
        rhs1 = s * (b - a) - a * b

        # q-p = a+b
        lhs2 = case["q"] - case["p"]
        rhs2 = a + b

        # max(a,b) = (a+b+|b-a|)/2
        lhs3 = case["radius"]
        rhs3 = (
            (a + b) + abs(b - a)
        ) // 2

        if lhs1 != rhs1:
            failures += 1

        if lhs2 != rhs2:
            failures += 1

        if lhs3 != rhs3:
            failures += 1

    print(
        f"identity failures = {failures}"
    )

    print()


def test_obvious_bounds(cases):
    print("============================================================")
    print("TEST 2: SIMPLE D-BASED BOUNDS")
    print("============================================================")

    bound_names = [
        "sqrt(D)",
        "D/s",
        "2sqrt(D)",
    ]

    covered = {
        name: 0
        for name in bound_names
    }

    violations = {
        name: 0
        for name in bound_names
    }

    for case in cases:
        d = abs(case["D"])
        s = case["s"]
        radius = case["radius"]

        candidates = {
            "sqrt(D)": math.isqrt(d),
            "D/s": d // s if s != 0 else 0,
            "2sqrt(D)": 2 * math.isqrt(d),
        }

        for name, bound in candidates.items():
            if radius <= bound:
                covered[name] += 1
            else:
                violations[name] += 1

    for name in bound_names:
        print(
            f"{name:12s}  "
            f"covers={covered[name]:5d}/{len(cases)}  "
            f"violations={violations[name]:5d}"
        )

    print()


def derive_radius_bound(case):
    """
    From

        |D| <= sR + R^2

    derive the positive solution of

        R^2 + sR - |D| >= 0.

    This gives a lower bound on R.
    """

    s = case["s"]
    d = abs(case["D"])

    root = math.sqrt(
        s * s + 4 * d
    )

    return (
        root - s
    ) / 2


def test_derived_bound(cases):
    print("============================================================")
    print("TEST 3: DERIVED RADIUS BOUND")
    print("============================================================")

    ratios = []
    violations = 0

    for case in cases:
        bound = derive_radius_bound(case)
        radius = case["radius"]

        if radius < bound:
            violations += 1

        if radius > 0:
            ratios.append(
                bound / radius
            )

    print(
        f"violations = {violations}"
    )

    print(
        f"bound/radius min  = "
        f"{min(ratios):.6f}"
    )

    print(
        f"bound/radius max  = "
        f"{max(ratios):.6f}"
    )

    print(
        f"bound/radius mean = "
        f"{sum(ratios) / len(ratios):.6f}"
    )

    print()


def test_D_sign(cases):
    print("============================================================")
    print("TEST 4: SIGN OF D VS b-a")
    print("============================================================")

    same_sign = 0
    opposite_sign = 0
    zero_difference = 0

    for case in cases:
        d = case["D"]
        difference = case["difference"]

        if difference == 0:
            zero_difference += 1
        elif d * difference > 0:
            same_sign += 1
        else:
            opposite_sign += 1

    print(
        f"D and (b-a) same sign     = "
        f"{same_sign}"
    )

    print(
        f"D and (b-a) opposite sign = "
        f"{opposite_sign}"
    )

    print(
        f"b-a = 0                    = "
        f"{zero_difference}"
    )

    print()


def test_ratio_to_D(cases):
    print("============================================================")
    print("TEST 5: NORMALIZED D")
    print("============================================================")

    values = []

    for case in cases:
        s = case["s"]
        d = abs(case["D"])
        radius = case["radius"]

        if s == 0 or radius == 0:
            continue

        values.append(
            d / (s * radius)
        )

    print(
        f"min |D|/(sR)  = "
        f"{min(values):.8f}"
    )

    print(
        f"max |D|/(sR)  = "
        f"{max(values):.8f}"
    )

    print(
        f"mean |D|/(sR) = "
        f"{sum(values) / len(values):.8f}"
    )

    print()


def test_exact_radius_from_a_b(cases):
    print("============================================================")
    print("TEST 6: EXACT RADIUS STRUCTURE")
    print("============================================================")

    a_equals_b = 0
    a_less_b = 0
    a_greater_b = 0

    for case in cases:
        a = case["a"]
        b = case["b"]

        if a == b:
            a_equals_b += 1
        elif a < b:
            a_less_b += 1
        else:
            a_greater_b += 1

    print(
        f"a=b  = {a_equals_b}"
    )

    print(
        f"a<b  = {a_less_b}"
    )

    print(
        f"a>b  = {a_greater_b}"
    )

    print()

    print("Largest-distance side:")

    print(
        f"s-p larger = {a_greater_b}"
    )

    print(
        f"q-s larger = {a_less_b}"
    )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 7: EXAMPLES")
    print("============================================================")

    for case in cases[:10]:
        print()

        print(
            f"N={case['N']}"
        )

        print(
            f"p={case['p']} "
            f"q={case['q']} "
            f"s={case['s']}"
        )

        print(
            f"a=s-p={case['a']}"
        )

        print(
            f"b=q-s={case['b']}"
        )

        print(
            f"D={case['D']}"
        )

        print(
            f"q-p={case['sum']}"
        )

        print(
            f"b-a={case['difference']}"
        )

        print(
            f"radius={case['radius']}"
        )

        difference = case["difference"]

        if difference != 0:
            ratio = (
                case["D"]
                / (
                    case["s"]
                    * difference
                )
            )

            print(
                f"D/[s(b-a)]={ratio:.8f}"
            )
        else:
            print(
                "D/[s(b-a)]=undefined"
            )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(555555)

    CASE_COUNT = 5000
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

    test_exact_identities(cases)
    test_obvious_bounds(cases)
    test_derived_bound(cases)
    test_D_sign(cases)
    test_ratio_to_D(cases)
    test_exact_radius_from_a_b(cases)
    print_examples(cases)

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()