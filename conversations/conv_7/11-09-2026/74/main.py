import math


EXPERIMENT = 93


def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    r = math.isqrt(n)

    d = 3

    while d <= r:
        if n % d == 0:
            return False

        d += 2

    return True


def smallest_semiprime_below_200():
    p = 23
    q = 37
    return p * q, p, q
 
    for n in range(4, 200):
        for p in range(2, math.isqrt(n) + 1):
            if n % p != 0:
                continue

            q = n // p

            if (
                is_prime(p)
                and is_prime(q)
            ):
                return n, p, q

    return None


def arthseq_value(coefficients, k):
    """
    Evaluate ArthSeq([a1, a2, ..., ad])
    at position k using the Newton/binomial form:

        v_k =
            a1*C(k-1,0)
          + a2*C(k-1,1)
          + a3*C(k-1,2)
          + ...

    """

    value = 0

    for j, coefficient in enumerate(coefficients):
        value += (
            coefficient
            * math.comb(k - 1, j)
        )

    return value


def arthseq_values(coefficients, count):
    return [
        arthseq_value(
            coefficients,
            k
        )
        for k in range(
            1,
            count + 1
        )
    ]


def repeated_coefficients(x, y, repetitions):
    coefficients = []

    for _ in range(repetitions):
        coefficients.append(y)
        coefficients.append(x)

    return coefficients


def gcd_pattern(values, n):
    return [
        math.gcd(
            value,
            n
        )
        for value in values
    ]


def classify_gcd(value, p, q, n):
    if value == 1:
        return "1"

    if value == p:
        return "p"

    if value == q:
        return "q"

    if value == n:
        return "n"

    return "?"


def print_pattern(
    coefficients,
    n,
    p,
    q,
    count=100
):
    values = arthseq_values(
        coefficients,
        count
    )

    gcds = gcd_pattern(
        values,
        n
    )

    labels = [
        classify_gcd(
            gcd_value,
            p,
            q,
            n
        )
        for gcd_value in gcds
    ]

    print(
        f"COEFFICIENTS = {coefficients}"
    )

    print(
        "VALUES:"
    )

    print(values)

    print(
        "GCD VALUES:"
    )

    print(gcds)

    print(
        "GCD CLASSES:"
    )

    print(labels)

    print()


def count_classes(
    values,
    n,
    p,
    q
):
    counts = {
        "1": 0,
        "p": 0,
        "q": 0,
        "n": 0,
        "?": 0,
    }

    for value in values:
        gcd_value = math.gcd(
            value,
            n
        )

        category = classify_gcd(
            gcd_value,
            p,
            q,
            n
        )

        counts[category] += 1

    return counts


def first_hit_positions(
    values,
    n,
    p,
    q
):
    positions = {
        "p": [],
        "q": [],
        "n": [],
    }

    for index, value in enumerate(
        values,
        start=1
    ):
        gcd_value = math.gcd(
            value,
            n
        )

        if gcd_value == p:
            positions["p"].append(
                index
            )

        elif gcd_value == q:
            positions["q"].append(
                index
            )

        elif gcd_value == n:
            positions["n"].append(
                index
            )

    return positions


def print_summary(
    repetitions,
    coefficients,
    n,
    p,
    q,
    count=100
):
    values = arthseq_values(
        coefficients,
        count
    )

    counts = count_classes(
        values,
        n,
        p,
        q
    )

    positions = first_hit_positions(
        values,
        n,
        p,
        q
    )

    print(
        f"REPETITIONS = {repetitions}"
    )

    print(
        f"DEGREE = {len(coefficients) - 1}"
    )

    print(
        f"COUNTS = {counts}"
    )

    print(
        f"p POSITIONS = "
        f"{positions['p']}"
    )

    print(
        f"q POSITIONS = "
        f"{positions['q']}"
    )

    print(
        f"n POSITIONS = "
        f"{positions['n']}"
    )

    print()


def run_experiment():
    print("=" * 60)
    print(
        f"START EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)
    print()

    result = smallest_semiprime_below_200()

    if result is None:
        print(
            "NO SEMIPRIME FOUND"
        )
        return

    n, p, q = result

    x = n - 1
    y = 1

    print(
        f"N = {n}"
    )

    print(
        f"P = {p}"
    )

    print(
        f"Q = {q}"
    )

    print(
        f"X = {x}"
    )

    print(
        f"Y = {y}"
    )

    print()

    print(
        "1. FULL FIRST-100 PATTERNS"
    )
    print("-" * 60)

    for repetitions in range(
        1,
        6
    ):
        coefficients = repeated_coefficients(
            x,
            y,
            repetitions
        )

        print(
            f"REPETITIONS = "
            f"{repetitions}"
        )

        print_pattern(
            coefficients,
            n,
            p,
            q,
            count=100
        )

    print(
        "2. HIT-COUNT GROWTH"
    )
    print("-" * 60)

    for repetitions in range(
        1,
        11
    ):
        coefficients = repeated_coefficients(
            x,
            y,
            repetitions
        )

        print_summary(
            repetitions,
            coefficients,
            n,
            p,
            q,
            count=100
        )

    print(
        "3. FIRST HIT POSITIONS"
    )
    print("-" * 60)

    for repetitions in range(
        1,
        16
    ):
        coefficients = repeated_coefficients(
            x,
            y,
            repetitions
        )

        values = arthseq_values(
            coefficients,
            100
        )

        positions = first_hit_positions(
            values,
            n,
            p,
            q
        )

        print(
            f"L={repetitions:2d} "
            f"p={positions['p']} "
            f"q={positions['q']} "
            f"n={positions['n']}"
        )

    print()

    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
