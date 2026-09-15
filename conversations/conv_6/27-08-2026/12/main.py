import sympy
import math


# ============================================================
# SETTINGS
# ============================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

# Only the final two layers
MODULI9 = [3, 5, 7, 11, 13, 17, 19, 23, 29]
MODULI10 = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


# ============================================================
# GENERATE SEMIPRIME
# ============================================================

p = sympy.randprime(FACTOR_LOW, FACTOR_HIGH)
q = sympy.randprime(FACTOR_LOW, FACTOR_HIGH)

n = p * q

actual_pair = tuple(sorted((p, q)))

sqrt_n = math.isqrt(n)

primes = list(sympy.primerange(2, sqrt_n + 1))


# ============================================================
# FAST CANDIDATE TEST
# ============================================================

def find_candidates(moduli):

    M = math.prod(moduli)

    candidates = []

    for candidate_p in primes:

        # Need modular inverse
        if math.gcd(candidate_p, M) != 1:
            continue

        # Required q residue:
        #
        # candidate_p * q = n (mod M)
        #
        # q = n / candidate_p (mod M)

        rq = (n * pow(candidate_p, -1, M)) % M

        # If M is larger than the allowed q range,
        # rq is the only possible q.
        if M > FACTOR_HIGH:

            if rq < 2:
                continue

            if rq > FACTOR_HIGH:
                continue

            if not sympy.isprime(rq):
                continue

            a, b = sorted((candidate_p, rq))
            candidates.append((a, b))

        else:

            # General case: q = rq + k*M
            q_candidate = rq

            if q_candidate < 2:
                q_candidate += (
                    (2 - q_candidate + M - 1) // M
                ) * M

            while q_candidate <= FACTOR_HIGH:

                if sympy.isprime(q_candidate):

                    a, b = sorted(
                        (candidate_p, q_candidate)
                    )

                    candidates.append((a, b))

                q_candidate += M

    return M, sorted(set(candidates))


# ============================================================
# RUN ONLY LAST TWO LAYERS
# ============================================================

print("=" * 80)
print("FINAL TWO FINGERPRINT LAYERS")
print("=" * 80)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print(f"sqrt(n) = {sqrt_n}")
print(f"Prime candidates <= sqrt(n) = {len(primes)}")
print()


for label, moduli in [
    ("LAYER 9", MODULI9),
    ("LAYER 10", MODULI10),
]:

    M, candidates = find_candidates(moduli)

    print("=" * 80)
    print(label)
    print("=" * 80)

    print(f"MODULI = {moduli}")
    print(f"M      = {M:,}")
    print()

    print(f"Candidate pairs = {len(candidates)}")

    actual_found = actual_pair in candidates

    print(f"Actual pair found = {actual_found}")
    print()

    if len(candidates) <= 30:

        for a, b in candidates:

            marker = ""

            if (a, b) == actual_pair:
                marker = " <-- ACTUAL"

            print(
                f"{a:>10,} * {b:>10,}"
                f"{marker}"
            )

    print()
