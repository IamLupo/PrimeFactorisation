import sympy
import math
import statistics


# ============================================================
# SETTINGS
# ============================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

TRIALS = 100

MODULI9 = [3, 5, 7, 11, 13, 17, 19, 23, 29]
MODULI10 = MODULI9 + [31]


# ============================================================
# CANDIDATE SEARCH
# ============================================================

def count_candidates(p, q, n, primes, moduli):

    M = math.prod(moduli)

    actual_pair = tuple(sorted((p, q)))

    candidates = set()

    for candidate_p in primes:

        if math.gcd(candidate_p, M) != 1:
            continue

        rq = (n * pow(candidate_p, -1, M)) % M

        # In our experiment M >> FACTOR_HIGH, so rq is the
        # only possible q in the tested range.

        if rq < FACTOR_LOW or rq > FACTOR_HIGH:
            continue

        if not sympy.isprime(rq):
            continue

        pair = tuple(sorted((candidate_p, rq)))

        candidates.add(pair)

    return M, candidates, actual_pair in candidates


# ============================================================
# MAIN
# ============================================================

results = []


print("=" * 100)
print("STATISTICAL TEST — FINAL TWO FINGERPRINT LAYERS")
print("=" * 100)

print(f"Factor range : {FACTOR_LOW:,} - {FACTOR_HIGH:,}")
print(f"Trials       : {TRIALS}")
print()

for trial in range(1, TRIALS + 1):

    # Generate factors
    p = sympy.randprime(FACTOR_LOW, FACTOR_HIGH)
    q = sympy.randprime(FACTOR_LOW, FACTOR_HIGH)

    n = p * q

    # Only primes <= sqrt(n)
    limit = math.isqrt(n)

    primes = list(sympy.primerange(2, limit + 1))

    # --------------------------------------------------------
    # Layer 9
    # --------------------------------------------------------

    M9, candidates9, found9 = count_candidates(
        p, q, n, primes, MODULI9
    )

    # --------------------------------------------------------
    # Layer 10
    # --------------------------------------------------------

    M10, candidates10, found10 = count_candidates(
        p, q, n, primes, MODULI10
    )

    results.append(
        (
            len(candidates9),
            len(candidates10),
            found9,
            found10,
        )
    )

    print(
        f"{trial:3}: "
        f"p={p:6} "
        f"q={q:6} "
        f"C9={len(candidates9):3} "
        f"C10={len(candidates10):3} "
        f"actual={found9 and found10}"
    )


# ============================================================
# SUMMARY
# ============================================================

c9 = [x[0] for x in results]
c10 = [x[1] for x in results]

unique9 = sum(x == 1 for x in c9)
unique10 = sum(x == 1 for x in c10)

zero9 = sum(x == 0 for x in c9)
zero10 = sum(x == 0 for x in c10)

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(f"Trials                    : {TRIALS}")

print()

print("LAYER 9")
print(f"  Minimum candidates      : {min(c9)}")
print(f"  Maximum candidates      : {max(c9)}")
print(f"  Mean candidates         : {statistics.mean(c9):.3f}")
print(f"  Median candidates       : {statistics.median(c9)}")
print(f"  Exactly 1 candidate     : {unique9}")
print(f"  Zero candidates         : {zero9}")

print()

print("LAYER 10")
print(f"  Minimum candidates      : {min(c10)}")
print(f"  Maximum candidates      : {max(c10)}")
print(f"  Mean candidates         : {statistics.mean(c10):.3f}")
print(f"  Median candidates       : {statistics.median(c10)}")
print(f"  Exactly 1 candidate     : {unique10}")
print(f"  Zero candidates         : {zero10}")

print()

print("=" * 100)
print("HISTOGRAM")
print("=" * 100)

print("C9 distribution:")

from collections import Counter

for count, amount in sorted(Counter(c9).items()):
    print(f"  {count:4} candidates : {amount:4} trials")

print()

print("C10 distribution:")

for count, amount in sorted(Counter(c10).items()):
    print(f"  {count:4} candidates : {amount:4} trials")
