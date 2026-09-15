import math
import random


EXPERIMENT = 50


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

        x_p = s + 1 - p
        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "x_p": x_p,
            "y": y,
        }


def H_value(case, x):
    s = case["s"]
    d = case["D"]

    return (
        4 * x * x
        - (s + 4) * x
        + 3 * d
        - 3 * s
    )


def H_derivative_twice_root(case):
    """
    H(x) = 4x^2 - (s+4)x + C

    Vertex:
        x = (s+4)/8
    """

    return (case["s"] + 4) / 8.0


def integer_root_estimates(case):
    """
    Real roots of H(x), when discriminant >= 0.
    """

    s = case["s"]
    d = case["D"]

    a = 4
    b = -(s + 4)
    c = 3 * d - 3 * s

    disc = b * b - 4 * a * c

    if disc < 0:
        return []

    root = math.sqrt(disc)

    r1 = (-b - root) / (2 * a)
    r2 = (-b + root) / (2 * a)

    return [r1, r2]


def nearest_integer(value):
    return min(
        [math.floor(value), math.ceil(value)],
        key=lambda x: abs(x - value),
    )


def distances_to_H_roots(case):
    xp = case["x_p"]

    roots = integer_root_estimates(case)

    if not roots:
        return None

    return min(
        abs(xp - root)
        for root in roots
    )


def distance_to_vertex(case):
    return abs(
        case["x_p"]
        - H_derivative_twice_root(case)
    )


def finite_difference_H(case, x):
    return (
        H_value(case, x + 1)
        - H_value(case, x)
    )


def second_difference_H(case):
    """
    Δ²H = 8 identically.
    """

    return (
        H_value(case, 2)
        - 2 * H_value(case, 1)
        + H_value(case, 0)
    )


def H_root_residual_at_true_x(case):
    """
    H(x_p) is known to contain p algebraically:

        H(x_p) = p(4p+3q-7s-4)

    This is only used for validation.
    """

    return H_value(
        case,
        case["x_p"],
    )


def scan_local_H(case, radius):
    """
    Search around the vertex and return the x which
    minimizes |H(x)|.
    """

    center = nearest_integer(
        H_derivative_twice_root(case)
    )

    candidates = range(
        max(0, center - radius),
        center + radius + 1,
    )

    best = min(
        candidates,
        key=lambda x: abs(H_value(case, x))
    )

    return best, abs(H_value(case, best))


def test_H_root_proximity(cases):
    print("============================================================")
    print("TEST 1: IS x_p NEAR A ROOT OF H?")
    print("============================================================")

    distances = []

    no_real_roots = 0

    for case in cases:
        value = distances_to_H_roots(case)

        if value is None:
            no_real_roots += 1
        else:
            distances.append(value)

    print(
        f"cases with real H roots = "
        f"{len(distances)}/{len(cases)}"
    )

    print(
        f"no real H roots          = "
        f"{no_real_roots}/{len(cases)}"
    )

    if distances:
        print(
            f"minimum distance = "
            f"{min(distances):.6f}"
        )

        print(
            f"maximum distance = "
            f"{max(distances):.6f}"
        )

        print(
            f"mean distance    = "
            f"{sum(distances) / len(distances):.6f}"
        )

    print()


def test_vertex_proximity(cases):
    print("============================================================")
    print("TEST 2: IS x_p NEAR THE H VERTEX?")
    print("============================================================")

    distances = [
        distance_to_vertex(case)
        for case in cases
    ]

    print(
        f"minimum distance = "
        f"{min(distances):.6f}"
    )

    print(
        f"maximum distance = "
        f"{max(distances):.6f}"
    )

    print(
        f"mean distance    = "
        f"{sum(distances) / len(distances):.6f}"
    )

    print()


def test_local_minimum(cases):
    print("============================================================")
    print("TEST 3: LOCAL MINIMUM OF |H|")
    print("============================================================")

    for radius in [1, 2, 4, 8, 16, 32]:
        hits = 0

        for case in cases:
            x, _ = scan_local_H(
                case,
                radius,
            )

            if x == case["x_p"]:
                hits += 1

        print(
            f"radius={radius:2d}  "
            f"x_p found={hits}/{len(cases)}"
        )

    print()


def test_H_xp_values(cases):
    print("============================================================")
    print("TEST 4: H(x_p) STRUCTURE")
    print("============================================================")

    failures = 0

    for case in cases:
        p = case["p"]
        q = case["q"]
        s = case["s"]

        expected = (
            p
            * (
                4 * p
                + 3 * q
                - 7 * s
                - 4
            )
        )

        actual = H_value(
            case,
            case["x_p"],
        )

        if actual != expected:
            failures += 1

    print(
        f"H(x_p) identity failures = "
        f"{failures}/{len(cases)}"
    )

    print()


def test_small_y_prediction(cases):
    print("============================================================")
    print("TEST 5: HOW SPECIAL IS y?")
    print("============================================================")

    counts = {
        "y<=1": 0,
        "y<=3": 0,
        "y<=7": 0,
        "y<=15": 0,
        "y<=31": 0,
        "y<=63": 0,
    }

    for case in cases:
        y = case["y"]

        if y <= 1:
            counts["y<=1"] += 1

        if y <= 3:
            counts["y<=3"] += 1

        if y <= 7:
            counts["y<=7"] += 1

        if y <= 15:
            counts["y<=15"] += 1

        if y <= 31:
            counts["y<=31"] += 1

        if y <= 63:
            counts["y<=63"] += 1

    for key, value in counts.items():
        print(
            f"{key:8s} = {value:5d}"
        )

    print()


def test_xp_prediction_from_H(cases):
    print("============================================================")
    print("TEST 6: CAN SIMPLE H FEATURES PREDICT x_p?")
    print("============================================================")

    exact = 0

    for case in cases:
        s = case["s"]

        vertex = H_derivative_twice_root(
            case
        )

        root_estimates = (
            integer_root_estimates(case)
        )

        candidates = [
            nearest_integer(vertex),
            math.floor(vertex),
            math.ceil(vertex),
        ]

        for root in root_estimates:
            candidates.append(
                nearest_integer(root)
            )

        if case["x_p"] in candidates:
            exact += 1

    print(
        f"x_p matched a simple H feature = "
        f"{exact}/{len(cases)}"
    )

    print()


def print_examples(cases):
    print("============================================================")
    print("TEST 7: EXAMPLES")
    print("============================================================")

    for case in cases[:10]:
        roots = integer_root_estimates(
            case
        )

        vertex = H_derivative_twice_root(
            case
        )

        print()
        print(
            f"N={case['N']}"
        )
        print(
            f"p={case['p']} q={case['q']}"
        )
        print(
            f"s={case['s']}"
        )
        print(
            f"x_p={case['x_p']}"
        )
        print(
            f"y={case['y']}"
        )
        print(
            f"H(x_p)={H_value(case, case['x_p'])}"
        )
        print(
            f"H vertex={vertex:.6f}"
        )
        print(
            f"H roots={roots}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(505050)

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

    test_H_root_proximity(cases)
    test_vertex_proximity(cases)
    test_local_minimum(cases)
    test_H_xp_values(cases)
    test_small_y_prediction(cases)
    test_xp_prediction_from_H(cases)
    print_examples(cases)

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
