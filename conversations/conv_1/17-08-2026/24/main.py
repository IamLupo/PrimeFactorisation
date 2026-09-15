import math
import random


def sieve_primes(lo: int, hi: int) -> list[int]:
    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, math.isqrt(hi) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [x for x in range(lo, hi + 1) if sieve[x]]


def power_sum(n: int, s: int, k: int) -> int:
    """Return p^k + q^k from n=pq and s=p+q."""
    if k == 0:
        return 2
    if k == 1:
        return s

    r0 = 2
    r1 = s

    for _ in range(2, k + 1):
        r0, r1 = r1, s * r1 - n * r0

    return r1


def sigma_semiprime(n: int, s: int, k: int) -> int:
    """
    For n=pq:
        sigma_k(n) = 1 + p^k + q^k + n^k
    """
    return 1 + power_sum(n, s, k) + n**k


def h6(n: int, s: int) -> int:
    """
    H6 = E1 / 6

    E1 = (n^2 - n + 1) sigma_1(n) - sigma_3(n)
    """
    numerator = (
        (n * n - n + 1) * sigma_semiprime(n, s, 1)
        - sigma_semiprime(n, s, 3)
    )

    if numerator % 6 != 0:
        raise ArithmeticError("H6 integrality failure")

    return numerator // 6


def h8(n: int, s: int) -> int:
    """
    24 H8 =
        -n^2 sigma_1(n)
        +(n^2+1) sigma_3(n)
        -sigma_5(n)
    """
    numerator = (
        -n * n * sigma_semiprime(n, s, 1)
        + (n * n + 1) * sigma_semiprime(n, s, 3)
        - sigma_semiprime(n, s, 5)
    )

    if numerator % 24 != 0:
        raise ArithmeticError("H8 integrality failure")

    return numerator // 24


def generate_cases(
    count: int = 20,
    prime_min: int = 2_000_000,
    prime_max: int = 4_200_000,
    seed: int = 12345,
) -> None:

    rng = random.Random(seed)

    primes = sieve_primes(
        prime_min,
        prime_max,
    )

    seen = set()

    generated = 0

    while generated < count:

        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        n = p * q
        s = p + q

        H6 = h6(n, s)
        H8 = h8(n, s)

        # Independent identity checks.
        A = s * ((n + 1) ** 2 - s ** 2)

        assert 6 * H6 == A
        assert 24 * H8 == A * (s ** 2 - 3 * n)
        assert 4 * H8 == H6 * (s ** 2 - 3 * n)

        print(f"{p} {q} {n} {H6} {H8} {4 * H8 // H6}")

        generated += 1


if __name__ == "__main__":
    generate_cases(
        count=1000,
        prime_min=2_000_000,
        prime_max=4_200_000,
        seed=12345,
    )

