import math
import random
from sympy import isprime


EXPERIMENT = 107

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

    value = math.comb(k - 2, d) % n
    value = (-value) % n

    return math.gcd(value, n)


def lucas_zero(m, d, p):
    """
    Lucas theorem:

        C(m,d) == 0 (mod p)

    iff at least one base-p digit of d exceeds
    the corresponding digit of m.
    """

    while m > 0 or d > 0:
        m_digit = m % p
        d_digit = d % p

        if d_digit > m_digit:
            return True

        m //= p
        d //= p

    return False


def smallest_lucas_zero_d(m, p):
    """
    Find the smallest d >= 1 such that

        C(m,d) == 0 mod p.

    This is done directly from the base-p digits rather
    than scanning all d.

    For d < p this is simply:

        d > (m mod p).

    For the general case we construct the smallest d
    whose first violating digit is as cheap as possible.
    """

    digits = []

    value = m

    while value > 0:
        digits.append(value % p)
        value //= p

    if not digits:
        return None

    # Try each digit position as the first violating
    # position. All lower digits are chosen minimally.
    candidates = []

    for i in range(len(digits)):
        current = digits[i]

        # Need d_i > m_i.
        if current + 1 >= p:
            continue

        # Keep higher digits zero and lower digits zero.
        candidate = (current + 1) * (p ** i)

        candidates.append(candidate)

    if not candidates:
        return None

    return min(candidates)


def smallest_admissible_z(p, k):
    """
    Find the smallest z such that

        d = 2z-1

    satisfies the Lucas zero condition.
    """

    m = k - 2

    # We need d odd.
    #
    # Starting from the smallest Lucas-zero d, move
    # upward by at most one if necessary to obtain odd d.

    d0 = smallest_lucas_zero_d(m, p)

    if d0 is None:
        return None

    if d0 % 2 == 0:
        d0 += 1

    if d0 >= k - 1:
        return None

    return (d0 + 1) // 2


def dyadic_first_hit(n, p, q, k):
    """
    Actual dyadic search.
    """

    z = 1

    while z <= MAX_Z:
        g = value_gcd(n, k, z)
        state = gcd_state(g, n, p, q)

        if state in ("p", "q"):
            return z

        if state == "n":
            return None

        z *= 2

    return None


def generate_semiprime():
    while True:
        p = random.randrange(P_MIN, P_MAX + 1)

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, P_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        if p != q:
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

    random.seed(107)

    total_k_cases = 0

    lucas_constructed = 0
    actual_constructed_hits = 0

    dyadic_hits = 0
    dyadic_misses = 0

    improvement_cases = 0

    construction_mismatches = 0

    largest_dyadic_ratio = 0.0
    largest_ratio_example = None

    miss_examples = []

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            z_lucas = smallest_admissible_z(p, k)

            if z_lucas is not None:
                lucas_constructed += 1

                g = value_gcd(
                    n,
                    k,
                    z_lucas,
                )

                state = gcd_state(
                    g,
                    n,
                    p,
                    q,
                )

                if state == "p":
                    actual_constructed_hits += 1
                else:
                    construction_mismatches += 1

            z_dyadic = dyadic_first_hit(
                n,
                p,
                q,
                k,
            )

            if z_dyadic is None:
                dyadic_misses += 1

                if z_lucas is not None:
                    improvement_cases += 1

                    if len(miss_examples) < 20:
                        miss_examples.append(
                            (
                                case_index,
                                n,
                                p,
                                q,
                                k,
                                z_lucas,
                            )
                        )
            else:
                dyadic_hits += 1

                if (
                    z_lucas is not None
                    and z_lucas < z_dyadic
                ):
                    improvement_cases += 1

                    ratio = (
                        z_dyadic
                        / z_lucas
                    )

                    if ratio > largest_dyadic_ratio:
                        largest_dyadic_ratio = ratio

                        largest_ratio_example = (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            z_lucas,
                            z_dyadic,
                            ratio,
                        )

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
        f"Lucas-construction cases     = "
        f"{lucas_constructed}"
    )

    print(
        f"constructed actual p hits    = "
        f"{actual_constructed_hits}"
    )

    print(
        f"construction mismatches      = "
        f"{construction_mismatches}"
    )

    print(
        f"dyadic hits                  = "
        f"{dyadic_hits}"
    )

    print(
        f"dyadic misses                = "
        f"{dyadic_misses}"
    )

    print(
        f"cases improved               = "
        f"{improvement_cases}"
    )

    if lucas_constructed:
        print(
            f"construction success rate   = "
            f"{100.0 * actual_constructed_hits / lucas_constructed:.2f}%"
        )

    if total_k_cases:
        print(
            f"dyadic success rate          = "
            f"{100.0 * dyadic_hits / total_k_cases:.2f}%"
        )

    print(
        f"largest dyadic/Lucas ratio   = "
        f"{largest_dyadic_ratio:.2f}"
    )

    if largest_ratio_example is not None:
        (
            case_index,
            n,
            p,
            q,
            k,
            z_lucas,
            z_dyadic,
            ratio,
        ) = largest_ratio_example

        print()
        print("LARGEST IMPROVEMENT")
        print(
            f"  CASE={case_index}"
        )
        print(
            f"  N={n}"
        )
        print(
            f"  p={p}"
        )
        print(
            f"  q={q}"
        )
        print(
            f"  k={k}"
        )
        print(
            f"  smallest Lucas z={z_lucas}"
        )
        print(
            f"  dyadic z={z_dyadic}"
        )
        print(
            f"  ratio={ratio:.2f}"
        )

    print()
    print("-" * 60)
    print("DYADIC MISS EXAMPLES")
    print("-" * 60)

    for (
        case_index,
        n,
        p,
        q,
        k,
        z_lucas,
    ) in miss_examples:
        print(
            f"CASE={case_index} "
            f"N={n} "
            f"p={p} "
            f"q={q} "
            f"k={k} "
            f"Lucas_z={z_lucas}"
        )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
