import sympy
import math


# ============================================================
# SETTINGS
# ============================================================

MODULI = [
    3, 5, 7, 11, 13,
    17, 19, 23, 29, 31
]

TEST_RANGES = [
    (10**3, 10**4),
    (10**4, 10**5),
    (10**5, 10**6),
]


# ============================================================
# SEMIPRIME GENERATION
# ============================================================

def generate_semiprime(low, high):
    p = sympy.randprime(low, high)
    q = sympy.randprime(low, high)

    return p, q, p * q


# ============================================================
# FAST FINGERPRINT FACTOR TEST
# ============================================================

def test_moduli(p, q, n, moduli, prime_set, factor_limit):

    M = math.prod(moduli)

    candidate_count = 0
    pairs = []

    for candidate_p in prime_set:

        # p must have an inverse modulo M.
        if math.gcd(candidate_p, M) != 1:
            continue

        inv = pow(candidate_p, -1, M)

        # Required q residue.
        rq = (n * inv) % M

        # If rq is zero, q cannot be a normal prime.
        if rq < 2:
            continue

        # If M > factor_limit, rq is already the only
        # possible positive q <= factor_limit.
        if M > factor_limit:

            if rq > factor_limit:
                continue

            if not sympy.isprime(rq):
                continue

            candidate_count += 1

            a, b = sorted((candidate_p, rq))

            pairs.append((a, b))

        else:

            # M is small enough that multiple q values can
            # exist in the range.
            q_candidate = rq

            if q_candidate < 2:
                q_candidate += (
                    (2 - q_candidate + M - 1) // M
                ) * M

            while q_candidate <= factor_limit:

                if sympy.isprime(q_candidate):

                    candidate_count += 1

                    a, b = sorted(
                        (candidate_p, q_candidate)
                    )

                    pairs.append((a, b))

                q_candidate += M

    pairs = sorted(set(pairs))

    actual = tuple(sorted((p, q)))

    return M, candidate_count, pairs, actual in pairs


# ============================================================
# MAIN
# ============================================================

print("=" * 90)
print("FAST SCALING TEST")
print("=" * 90)

print("Moduli:")
print(MODULI)
print()


for low, high in TEST_RANGES:

    p, q, n = generate_semiprime(low, high)

    print()
    print("=" * 90)
    print(f"FACTOR RANGE: {low:,} - {high:,}")
    print("=" * 90)

    print(f"p = {p:,}")
    print(f"q = {q:,}")
    print(f"n = {n:,}")

    print()

    # We only need primes up to the smaller factor,
    # because p <= sqrt(n).
    limit = math.isqrt(n)

    primes = list(sympy.primerange(2, limit + 1))

    print(
        f"sqrt(n) = {limit:,}"
    )

    print(
        f"Prime candidates = {len(primes):,}"
    )

    print()

    # --------------------------------------------------------
    # Incrementally add moduli
    # --------------------------------------------------------

    for i in range(1, len(MODULI) + 1):

        current = MODULI[:i]

        M, count, pairs, found = test_moduli(
            p,
            q,
            n,
            current,
            primes,
            high
        )

        print(
            f"[{i:2}] "
            f"M={M:,} "
            f"candidate_pairs={len(pairs):,} "
            f"actual={found}"
        )

        if len(pairs) <= 10:

            for a, b in pairs:

                marker = ""

                if (a, b) == tuple(sorted((p, q))):
                    marker = " <-- ACTUAL"

                print(
                    f"     {a:,} * {b:,}"
                    f"{marker}"
                )

        # Once unique, no need to add more moduli.
        if len(pairs) == 1:

            print(
                "     UNIQUE FACTOR PAIR FOUND"
            )

            break
