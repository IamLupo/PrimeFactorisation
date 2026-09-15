import math
import random
from sympy import isprime


EXPERIMENT = 109

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384


def gcd_state(g, n, p, q):
    if g == 1:
        return "1"

    if g == n:
        return "n"

    if g == p:
        return "p"

    if g == q:
        return "q"

    return "other"


def value_gcd(n, k, z):
    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d)
    value %= n
    value = (-value) % n

    return math.gcd(value, n)


def base_p_digits(value, p):
    digits = []

    while value > 0:
        digits.append(value % p)
        value //= p

    if not digits:
        digits.append(0)

    return digits


def smallest_lucas_zero_d(m, p):
    """
    Find the smallest positive d such that

        C(m,d) == 0 mod p

    using Lucas' theorem.

    We construct candidates by making one base-p digit
    of d exceed the corresponding digit of m while
    keeping lower and higher digits minimal.
    """

    digits = base_p_digits(m, p)

    candidates = []

    power = 1

    for digit in digits:
        if digit + 1 < p:
            candidate = (digit + 1) * power
            candidates.append(candidate)

        power *= p

    if not candidates:
        return None

    return min(candidates)


def smallest_valid_z(p, k):
    """
    Find the smallest z for which

        C(k-2, 2z-1) == 0 mod p.

    The required d must be odd.
    """

    m = k - 2

    d = smallest_lucas_zero_d(m, p)

    if d is None:
        return None

    # Make d odd.
    if d % 2 == 0:
        d += 1

    if d >= k - 1:
        return None

    return (d + 1) // 2


def odd_part(value):
    """
    Remove all factors of 2.
    """

    while value % 2 == 0:
        value //= 2

    return value


def v2(value):
    """
    2-adic valuation.
    """

    count = 0

    while value % 2 == 0:
        value //= 2
        count += 1

    return count


def generate_semiprime():
    while True:
        p = random.randrange(
            P_MIN,
            P_MAX + 1,
        )

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, P_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        if p == q:
            continue

        return p, q


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def original_dyadic_hit(n, p, q, k):
    """
    Original search:

        z = 1,2,4,8,...

    """

    z = 1

    while z <= MAX_Z:
        d = 2 * z - 1

        if d >= k - 1:
            return None

        g = value_gcd(
            n,
            k,
            z,
        )

        state = gcd_state(
            g,
            n,
            p,
            q,
        )

        if state == "p" or state == "q":
            return z

        if state == "n":
            return None

        z *= 2

    return None


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q MAX: {Q_MAX}")
    print(f"MAX Z: {MAX_Z}")
    print()

    random.seed(109)

    total_k_cases = 0

    dyadic_hits = 0
    dyadic_misses = 0

    constructed_hits = 0
    construction_failures = 0

    max_c = 0
    max_c_example = None

    c_histogram = {}

    miss_c_values = []

    total_exact_z = 0

    examples = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            z_exact = smallest_valid_z(
                p,
                k,
            )

            if z_exact is None:
                construction_failures += 1
                continue

            total_exact_z += 1

            g = value_gcd(
                n,
                k,
                z_exact,
            )

            state = gcd_state(
                g,
                n,
                p,
                q,
            )

            if state == "p":
                constructed_hits += 1
            else:
                construction_failures += 1

            c = odd_part(z_exact)

            exponent = v2(z_exact)

            c_histogram[c] = (
                c_histogram.get(c, 0) + 1
            )

            if c > max_c:
                max_c = c

                max_c_example = (
                    case_index,
                    n,
                    p,
                    q,
                    k,
                    z_exact,
                    c,
                    exponent,
                    g,
                    state,
                )

            z_dyadic = original_dyadic_hit(
                n,
                p,
                q,
                k,
            )

            if z_dyadic is None:
                dyadic_misses += 1

                miss_c_values.append(c)

                if len(examples) < 25:
                    examples.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            z_exact,
                            c,
                            exponent,
                        )
                    )
            else:
                dyadic_hits += 1

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"semiprime cases              = "
        f"{CASES}"
    )

    print(
        f"(N,k) cases                  = "
        f"{total_k_cases}"
    )

    print(
        f"exact Lucas constructions    = "
        f"{total_exact_z}"
    )

    print(
        f"constructed p hits           = "
        f"{constructed_hits}"
    )

    print(
        f"construction failures        = "
        f"{construction_failures}"
    )

    print(
        f"original dyadic hits         = "
        f"{dyadic_hits}"
    )

    print(
        f"original dyadic misses       = "
        f"{dyadic_misses}"
    )

    if total_exact_z:
        print(
            f"construction success rate   = "
            f"{100.0 * constructed_hits / total_exact_z:.2f}%"
        )

    if total_k_cases:
        print(
            f"dyadic success rate          = "
            f"{100.0 * dyadic_hits / total_k_cases:.2f}%"
        )

    print()

    print(
        f"maximum required odd c       = "
        f"{max_c}"
    )

    if miss_c_values:
        print(
            f"maximum c among dyadic misses = "
            f"{max(miss_c_values)}"
        )

        print(
            f"minimum c among dyadic misses = "
            f"{min(miss_c_values)}"
        )

    print()
    print("-" * 60)
    print("MULTIPLIER HISTOGRAM")
    print("-" * 60)

    for c in sorted(c_histogram):
        print(
            f"c={c:5d} "
            f"count={c_histogram[c]}"
        )

    if max_c_example is not None:
        print()
        print("-" * 60)
        print("MAXIMUM c EXAMPLE")
        print("-" * 60)

        (
            case_index,
            n,
            p,
            q,
            k,
            z,
            c,
            exponent,
            g,
            state,
        ) = max_c_example

        print(
            f"CASE={case_index}"
        )

        print(
            f"N={n}"
        )

        print(
            f"p={p}"
        )

        print(
            f"q={q}"
        )

        print(
            f"k={k}"
        )

        print(
            f"exact_z={z}"
        )

        print(
            f"odd_multiplier_c={c}"
        )

        print(
            f"power_of_two_exponent={exponent}"
        )

        print(
            f"gcd={g}"
        )

        print(
            f"state={state}"
        )

    print()
    print("-" * 60)
    print("DYADIC MISS / EXACT-z COMPARISON")
    print("-" * 60)

    if examples:
        for (
            case_index,
            n,
            p,
            q,
            k,
            z,
            c,
            exponent,
        ) in examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"exact_z={z} "
                f"c={c} "
                f"2^j={2 ** exponent}"
            )
    else:
        print("NONE")

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
