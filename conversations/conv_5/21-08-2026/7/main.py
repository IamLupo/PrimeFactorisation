#!/usr/bin/env python3

from math import isqrt


MAX_N = 500_000


def prime_sieve(limit):
    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if not sieve[p]:
            continue

        start = p * p
        count = (limit - start) // p + 1
        sieve[start:limit + 1:p] = b"\x00" * count

    return sieve


def main():

    sieve = prime_sieve(MAX_N)

    primes = [
        p
        for p in range(3, MAX_N + 1, 2)
        if sieve[p]
    ]

    for i, p in enumerate(primes):

        if p * p > MAX_N:
            break

        max_q = MAX_N // p

        for q in primes[i:]:

            if q > max_q:
                break

            n = p * q

            # Only n = 1 mod 8.
            if n % 8 != 7:
                continue

            # 2020 x,y construction.
            x_num = q - p + 6
            y_num = p + q

            # Require integer x,y.
            #if x_num % 4 != 0:
            #    continue

            if y_num % 4 != 0:
                continue

            x = x_num // 4
            y = y_num // 4

            print(
                f"{n}	{p}	{q}	{x}	{y}	{x+y}"
            )


if __name__ == "__main__":
    main()