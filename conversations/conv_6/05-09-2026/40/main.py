#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 123
# Pollard-rho control + carry verification
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

R00 = R1[0] * R2[0]
R11 = R1[1] * R2[1]
M = math.lcm(R00, R11)

RANDOM_SEED = 123


# ============================================================
# PRIME TEST
# ============================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    )

    for p in small:

        if n % p == 0:
            return n == p

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in (
        2, 3, 5, 7,
        11, 13, 17
    ):

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


def random_prime(bits):

    while True:

        p = random.getrandbits(bits)

        p |= 1 << (bits - 1)
        p |= 1

        if is_probable_prime(p):
            return p


def make_semiprime(bits):

    fb = bits // 2

    while True:

        p = random_prime(fb)
        q = random_prime(fb)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n.bit_length() >= bits - 1:
            return n, p, q


# ============================================================
# CARRY
# ============================================================

def carry_E(r1, r2, k, l, a, b):

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    return c1 + c2 + c3


# ============================================================
# FULL CARRY VERIFICATION
# ============================================================

def verify_factorization(n, p, q):

    if p > q:
        p, q = q, p

    if p * q != n:
        return False

    for r1 in R1:

        k = p // r1
        a = p % r1

        for r2 in R2:

            l = q // r2
            b = q % r2

            Q = n // (r1 * r2)

            E1 = Q - k * l

            E2 = carry_E(
                r1,
                r2,
                k,
                l,
                a,
                b
            )

            if E1 != E2:
                return False

    return True


# ============================================================
# POLLARD RHO
# ============================================================

def pollard_rho(n):

    if n % 2 == 0:
        return 2

    if n % 3 == 0:
        return 3

    while True:

        c = random.randrange(
            1,
            n - 1
        )

        x = random.randrange(
            2,
            n - 1
        )

        y = x

        d = 1

        iterations = 0

        while d == 1:

            x = (
                x * x
                + c
            ) % n

            y = (
                y * y
                + c
            ) % n

            y = (
                y * y
                + c
            ) % n

            d = math.gcd(
                abs(x - y),
                n
            )

            iterations += 1

            if iterations > 10_000_000:
                break

        if d > 1 and d < n:
            return d


# ============================================================
# POLLARD RHO WITH COUNTER
# ============================================================

def rho_factor(n):

    start = time.perf_counter()

    if is_probable_prime(n):
        return {
            "factor": n,
            "iterations": 0,
            "time":
                time.perf_counter()
                - start,
        }

    while True:

        if n % 2 == 0:

            return {
                "factor": 2,
                "iterations": 0,
                "time":
                    time.perf_counter()
                    - start,
            }

        c = random.randrange(
            1,
            n - 1
        )

        x = random.randrange(
            2,
            n - 1
        )

        y = x

        d = 1
        iterations = 0

        while d == 1:

            x = (
                x * x
                + c
            ) % n

            y = (
                y * y
                + c
            ) % n

            y = (
                y * y
                + c
            ) % n

            d = math.gcd(
                abs(x - y),
                n
            )

            iterations += 1

        if d != n:

            return {
                "factor": d,
                "iterations": iterations,
                "time":
                    time.perf_counter()
                    - start,
            }


# ============================================================
# DIRECT O(sqrt(n)) CONTROL
# ============================================================

def direct_factor(n):

    limit = math.isqrt(n)

    tested = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        tested += 1

        if n % p == 0:

            q = n // p

            return {
                "p": p,
                "q": q,
                "tested": tested,
                "time":
                    time.perf_counter()
                    - start,
            }

    return {
        "p": None,
        "q": None,
        "tested": tested,
        "time":
            time.perf_counter()
            - start,
    }


# ============================================================
# MODULAR-CARRY SEARCH
# ============================================================

def carry_modular_factor(n):

    """
    Experiment-113/121 style control.

    For every p < sqrt(n):

        pq == n (mod 15)
        pq == n (mod 899)

    gives q modulo 13485.

    The second quotient bucket then limits q.
    """

    limit = math.isqrt(n)

    p_tested = 0
    q_candidates = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        p_tested += 1

        if math.gcd(
            p,
            M
        ) != 1:

            continue

        q_res = (
            n
            * pow(
                p,
                -1,
                M
            )
        ) % M

        Q11 = n // R11

        q_low = (
            R11 * Q11
            + p - 1
        ) // p

        q_high = (
            R11 * (Q11 + 1)
            - 1
        ) // p

        q = (
            q_low
            + (q_res - q_low) % M
        )

        while q <= q_high:

            q_candidates += 1

            if p <= q:

                if p * q == n:

                    return {
                        "p": p,
                        "q": q,
                        "p_tested":
                            p_tested,
                        "q_candidates":
                            q_candidates,
                        "time":
                            time.perf_counter()
                            - start,
                    }

            q += M

    return {
        "p": None,
        "q": None,
        "p_tested":
            p_tested,
        "q_candidates":
            q_candidates,
        "time":
            time.perf_counter()
            - start,
    }


# ============================================================
# HYBRID: POLLARD RHO + CARRY
# ============================================================

def hybrid_factor(n):

    result = rho_factor(n)

    factor = result["factor"]

    if factor == n:
        return {
            **result,
            "verified": False,
        }

    other = n // factor

    if factor > other:
        factor, other = other, factor

    verified = verify_factorization(
        n,
        factor,
        other
    )

    return {
        **result,
        "p": factor,
        "q": other,
        "verified": verified,
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    direct,
    modular,
    hybrid,
):

    print()
    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print()
    print(
        "n      =",
        n
    )

    print(
        "true p =",
        true_p
    )

    print(
        "true q =",
        true_q
    )

    print()
    print("TRUE GAP")

    print(
        "q-p =",
        true_q - true_p
    )

    print()
    print("DIRECT sqrt(n) SCAN")

    print(
        "p tested =",
        direct["tested"]
    )

    print(
        "runtime  = %.6f s"
        % direct["time"]
    )

    print()
    print("TWO-CELL MODULAR SCAN")

    print(
        "p tested      =",
        modular["p_tested"]
    )

    print(
        "q candidates  =",
        modular["q_candidates"]
    )

    print(
        "runtime       = %.6f s"
        % modular["time"]
    )

    print()
    print("POLLARD RHO")

    print(
        "iterations =",
        hybrid["iterations"]
    )

    print(
        "runtime    = %.6f s"
        % hybrid["time"]
    )

    print(
        "factor     =",
        hybrid["factor"]
    )

    print(
        "p           =",
        hybrid["p"]
    )

    print(
        "q           =",
        hybrid["q"]
    )

    print(
        "carry grid verified =",
        hybrid["verified"]
    )

    print()
    print("TIMING RATIOS")

    if hybrid["time"] > 0:

        print(
            "direct / rho = %.3fx"
            % (
                direct["time"]
                / hybrid["time"]
            )
        )

        print(
            "modular / rho = %.3fx"
            % (
                modular["time"]
                / hybrid["time"]
            )
        )

    print()
    print(
        "TRUE FACTORIZATION FOUND =",
        (
            hybrid["p"] == true_p
            and
            hybrid["q"] == true_q
        )
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(
        RANDOM_SEED
    )

    print("=" * 72)
    print("START EXPERIMENT 123")
    print("Pollard-rho control + carry verification")
    print("=" * 72)

    print()
    print("R00 =", R00)
    print("R11 =", R11)
    print("M   =", M)

    for bits in (
        30,
        36,
        42,
        48,
        54,
    ):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(
            bits
        )

        direct = direct_factor(
            n
        )

        modular = carry_modular_factor(
            n
        )

        hybrid = hybrid_factor(
            n
        )

        report(
            n,
            p,
            q,
            direct,
            modular,
            hybrid
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 123")
    print("=" * 72)


if __name__ == "__main__":
    main()
