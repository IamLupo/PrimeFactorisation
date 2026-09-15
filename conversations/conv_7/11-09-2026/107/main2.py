import math
import random

from sympy import isprime


EXPERIMENT = 132

CASES = 100

N_MIN =   100_000_000_000_000
N_MAX = 1_000_000_000_000_000

# Number of k values to generate for each semiprime.
K_COUNT = 4


def random_prime(low, high):
    """
    Generate a random prime in [low, high].
    """
    while True:
        value = random.randrange(low, high)

        if value < 2:
            continue

        if value % 2 == 0:
            value += 1

        if isprime(value):
            return value


def generate_case():
    """
    Generate p, q with

        p < q
        10^7 <= p*q < 10^8
    """
    p = random_prime(math.floor(math.sqrt(N_MIN)), math.floor(math.isqrt(N_MAX)))
    q = random_prime(math.ceil(math.sqrt(N_MIN)), math.floor(math.isqrt(N_MAX)))

    if p >= q:
        p, q = q, p

    n = p * q

    return p, q, n


def dyadic_k_values(p, q):
    """
    Return dyadic k values satisfying

        p < k < q.

    Only powers of two are used.
    """
    values = []

    k = 1

    while k < q:
        if p < k < q:
            values.append(k)

        k *= 2

    return values


def generate_z_values(k):
    """
    Generate the standard dyadic z values

        z = 1, 2, 4, 8, ...

    subject to

        d = 2z - 1 < k - 1.

    This keeps the numerator interval entirely positive.
    """
    values = []

    z = 1

    while True:
        d = 2 * z - 1

        if d >= k - 1:
            break

        values.append(z)

        z *= 2

    return values


def select_z_values(k):
    """
    Select the largest few valid dyadic z values.

    These are the probes closest to the central part of
    the ArithSeg interval while remaining valid.
    """
    values = generate_z_values(k)

    if len(values) <= K_COUNT:
        return values

    return values[-K_COUNT:]


def main():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"N RANGE: {N_MIN}..{N_MAX - 1}")
    print(f"K COUNT: {K_COUNT}")
    print()

    #random.seed(132)

    case_number = 0

    while case_number < CASES:
        p, q, n = generate_case()

        k_values = dyadic_k_values(p, q)

        if not k_values:
            continue

        case_number += 1

        print(f"CASE {case_number}")
        print(f"  p = {p}")
        print(f"  q = {q}")
        print(f"  N = {n}")
        print(f"  sqrt(N) = {math.isqrt(n)}")
        print()

        print("  ARITHSEG DATA")

        # Use the first K_COUNT usable dyadic k values.
        selected_k = k_values[:K_COUNT]

        for k in selected_k:
            z_values = select_z_values(k)

            print(f"    k = {k}")

            for z in z_values:
                d = 2 * z - 1

                # ArithSeg numerator interval.
                m = k - 2
                left = m - d + 1
                right = m

                valid = (
                    d < p
                    and d < k - 1
                    and left > 0
                )

                print(
                    f"      z={z:<8} "
                    f"d={d:<8} "
                    f"interval=[{left},{right}] "
                    f"d<p={d < p} "
                    f"valid={valid}"
                )

            print()

        print("-" * 60)

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


if __name__ == "__main__":
    main()

