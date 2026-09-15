import sympy
import math
import random


# ============================================================
# EXPERIMENT SETTINGS
# ============================================================

# Add moduli one by one.
MODULI = [
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
]

# Different factor-size experiments.
# Each tuple is (low, high) for p and q.
TEST_RANGES = [
    (10**3, 10**4),
    (10**4, 10**5),
    (10**5, 10**6),
    (10**6, 10**7),
]


# ============================================================
# HELPERS
# ============================================================

def make_semiprime(low, high):
    p = sympy.randprime(low, high)
    q = sympy.randprime(low, high)

    return p, q, p * q


def fingerprint(x, moduli):
    return tuple(x % m for m in moduli)


def valid_factor_residue(p_candidate, n, M):
    """
    Return q residue satisfying

        p_candidate * q = n (mod M)

    or None if p_candidate is not invertible modulo M.
    """

    if math.gcd(p_candidate, M) != 1:
        return None

    inv = pow(p_candidate, -1, M)

    return (n * inv) % M


# ============================================================
# ONE SCALE TEST
# ============================================================

def run_test(p, q, n, limit):

    print()
    print("=" * 90)
    print(f"TEST: factors approximately {limit:,}")
    print("=" * 90)

    print(f"p = {p}")
    print(f"q = {q}")
    print(f"n = {n}")
    print(f"sqrt(n) = {math.isqrt(n)}")
    print()

    # All primes up to sqrt(n).
    candidate_primes = list(
        sympy.primerange(2, math.isqrt(n) + 1)
    )

    print(
        f"Prime candidates p <= sqrt(n): "
        f"{len(candidate_primes):,}"
    )

    print()

    # --------------------------------------------------------
    # Incrementally add moduli
    # --------------------------------------------------------

    for count in range(1, len(MODULI) + 1):

        moduli = MODULI[:count]
        M = math.prod(moduli)

        possible_p = 0
        possible_pairs = []

        # ----------------------------------------------------
        # For each candidate p:
        #
        # q = n * p^-1 (mod M)
        #
        # Then q must be an actual prime factor in
        # the allowed range.
        # ----------------------------------------------------

        for candidate_p in candidate_primes:

            rq = valid_factor_residue(
                candidate_p,
                n,
                M
            )

            if rq is None:
                continue

            # We only care about q values in the same
            # factor-size range used for this experiment.
            #
            # q = rq + k*M

            if rq == 0:
                k = 0
            else:
                k = 0

            # Find the smallest k giving q >= 2.
            if rq < 2:
                k = (2 - rq + M - 1) // M

            candidate_q = rq + k * M

            # If M is smaller than the factors, there can
            # be many q values in the class. Search them.
            while candidate_q <= limit:

                if candidate_q >= 2:

                    if sympy.isprime(candidate_q):

                        possible_p += 1

                        possible_pairs.append(
                            (
                                candidate_p,
                                candidate_q
                            )
                        )

                candidate_q += M

        # ----------------------------------------------------
        # Remove duplicate ordered pairs.
        # ----------------------------------------------------

        possible_pairs = sorted(
            set(
                tuple(sorted(pair))
                for pair in possible_pairs
            )
        )

        actual_pair = tuple(sorted((p, q)))

        actual_found = actual_pair in possible_pairs

        print(
            f"[{count:2}] "
            f"moduli={moduli} "
            f"M={M:,}"
        )

        print(
            f"     prime-compatible p values = "
            f"{possible_p:,}"
        )

        print(
            f"     unordered prime pairs      = "
            f"{len(possible_pairs):,}"
        )

        print(
            f"     actual factor pair found   = "
            f"{actual_found}"
        )

        # Print candidates only when the list becomes small.
        if len(possible_pairs) <= 20:

            print("     candidates:")

            for a, b in possible_pairs:

                marker = ""

                if (a, b) == actual_pair:
                    marker = " <-- ACTUAL"

                print(
                    f"         {a:,} * {b:,}"
                    f"{marker}"
                )

        print()


# ============================================================
# MAIN
# ============================================================

print("=" * 90)
print("SCALING TEST FOR MULTIPLE MODULAR FINGERPRINTS")
print("=" * 90)

print("Moduli:")
print(MODULI)
print()

for low, high in TEST_RANGES:

    p, q, n = make_semiprime(low, high)

    # The search for q is bounded by the larger factor range.
    run_test(
        p,
        q,
        n,
        high
    )
