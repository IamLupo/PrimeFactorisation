import math
import random
import time

from sympy import isprime


EXPERIMENT = 126

CASES = 3

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

MAX_PROBES = 2


def generate_semiprime():
    while True:
        p = random.randrange(
            P_MIN,
            P_MAX + 1,
        )

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, Q_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
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


def vdc_z(k):
    """
    The second VdC probe is the important one in
    the previous experiments.
    """
    max_z = k // 2

    # For the power-of-two ranges we are using,
    # the second VdC point is approximately k/4.
    return (max_z // 2) + 1


def exact_comb(m, d, n):
    start = time.perf_counter()

    value = math.comb(
        m,
        d,
    ) % n

    elapsed = time.perf_counter() - start

    return value, elapsed


def batch_comb_mod(m, d, n):
    """
    Batch inversion.

    This assumes gcd(j,n)=1 for all 1 <= j <= d.

    Instead of computing j^(-1) individually:

        P_j = j!

    is built modulo n.

    Then P_d is inverted ONCE.

    A backward pass recovers every inverse.
    """

    start = time.perf_counter()

    prefix = [1] * (d + 1)

    product = 1

    for j in range(1, d + 1):
        product = (
            product * j
        ) % n

        prefix[j] = product

    factor = math.gcd(
        product,
        n,
    )

    if factor != 1:
        return None, (
            time.perf_counter() - start
        ), factor

    inverse_product = pow(
        product,
        -1,
        n,
    )

    result = 1

    for j in range(d, 0, -1):
        inverse_j = (
            inverse_product
            * prefix[j - 1]
        ) % n

        numerator = (
            m - d + j
        )

        result *= (
            numerator % n
        )
        result %= n

        result *= inverse_j
        result %= n

        inverse_product = (
            inverse_product * j
        ) % n

    elapsed = (
        time.perf_counter()
        - start
    )

    return result, elapsed, None


def recurrence_with_single_batch_inverse(
    m,
    d,
    n,
):
    """
    Same mathematics as the recurrence, but instead of
    computing pow(j,-1,n) d times, calculate all inverses
    with one batch inversion.
    """

    # This is deliberately kept separate from the above
    # implementation so we can benchmark the exact method.
    return batch_comb_mod(
        m,
        d,
        n,
    )


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX}")
    print()

    random.seed(EXPERIMENT)

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        ks = powers_between(
            p,
            q,
        )

        if not ks:
            continue

        k = ks[0]

        z = vdc_z(k)

        d = 2 * z - 1
        m = k - 2

        print("-" * 60)
        print(
            f"CASE {case_index}"
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
            f"z={z}"
        )

        print(
            f"d={d}"
        )

        print(
            f"d/p={d / p:.6f}"
        )

        # -------------------------------------------------
        # EXACT
        # -------------------------------------------------

        exact_value, exact_time = (
            exact_comb(
                m,
                d,
                n,
            )
        )

        exact_gcd = math.gcd(
            exact_value,
            n,
        )

        print()
        print("EXACT math.comb")

        print(
            f"  time={exact_time:.6f}s"
        )

        print(
            f"  gcd={exact_gcd}"
        )

        # -------------------------------------------------
        # BATCH
        # -------------------------------------------------

        (
            batch_value,
            batch_time,
            batch_factor,
        ) = batch_comb_mod(
            m,
            d,
            n,
        )

        print()
        print("BATCH MODULAR")

        print(
            f"  time={batch_time:.6f}s"
        )

        if batch_factor is not None:
            print(
                f"  factor={batch_factor}"
            )

        else:
            batch_gcd = math.gcd(
                batch_value,
                n,
            )

            print(
                f"  gcd={batch_gcd}"
            )

            print(
                f"  matches_exact="
                f"{batch_value == exact_value}"
            )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()