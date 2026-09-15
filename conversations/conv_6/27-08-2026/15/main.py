import sympy
import math
from collections import Counter


# ============================================================
# SETTINGS
# ============================================================

TRIALS = 100

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

MODULI = [
    3, 5, 7, 11, 13,
    17, 19, 23, 29, 31
]


# ============================================================
# PRECOMPUTE LAYERS
# ============================================================

LAYERS = []

M = 1

for m in MODULI:
    M *= m

    LAYERS.append({
        "modulus": M,
        "moduli": MODULI[:len(LAYERS) + 1]
    })


# ============================================================
# FIND NUMBER OF SURVIVING PAIRS
# ============================================================

def count_candidates(n, primes, M, factor_high):

    count = 0

    for p_candidate in primes:

        if math.gcd(p_candidate, M) != 1:
            continue

        # p*q = n (mod M)
        #
        # q = n * p^-1 (mod M)

        q_residue = (
            n * pow(p_candidate, -1, M)
        ) % M

        # M is already much larger than the factor range
        # for the interesting layers.
        if q_residue < FACTOR_LOW:
            continue

        if q_residue > factor_high:
            continue

        if not sympy.isprime(q_residue):
            continue

        count += 1

    return count


# ============================================================
# MAIN
# ============================================================

print("=" * 90)
print("THRESHOLD EXPERIMENT")
print("=" * 90)

print(
    f"Factor range: "
    f"{FACTOR_LOW:,} - {FACTOR_HIGH:,}"
)

print(f"Trials: {TRIALS}")
print()


first_unique_layers = []

for trial in range(1, TRIALS + 1):

    # --------------------------------------------------------
    # Generate semiprime
    # --------------------------------------------------------

    p = sympy.randprime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    q = sympy.randprime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    n = p * q

    sqrt_n = math.isqrt(n)

    primes = list(
        sympy.primerange(
            2,
            sqrt_n + 1
        )
    )

    actual_pair = tuple(sorted((p, q)))

    # --------------------------------------------------------
    # Find FIRST unique layer
    # --------------------------------------------------------

    first_unique = None

    counts = []

    for layer_number, layer in enumerate(LAYERS, start=1):

        M = layer["modulus"]

        count = count_candidates(
            n,
            primes,
            M,
            FACTOR_HIGH
        )

        counts.append(count)

        if count == 1:

            first_unique = layer_number

            break

    first_unique_layers.append(first_unique)

    # --------------------------------------------------------
    # Print one compact line
    # --------------------------------------------------------

    print(
        f"{trial:3}: "
        f"p={p:6,} "
        f"q={q:6,} "
        f"first_unique={first_unique} "
        f"counts={counts}"
    )


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 90)
print("SUMMARY")
print("=" * 90)

distribution = Counter(first_unique_layers)

for layer in sorted(distribution):

    amount = distribution[layer]

    M = LAYERS[layer - 1]["modulus"]

    print(
        f"Layer {layer:2}: "
        f"M={M:,} "
        f"{amount:3}/{TRIALS}"
    )


# ============================================================
# IMPORTANT RATIO
# ============================================================

print()
print("=" * 90)
print("MODULUS / FACTOR-SIZE RATIOS")
print("=" * 90)

for i, layer in enumerate(LAYERS, start=1):

    M = layer["modulus"]

    ratio = M / FACTOR_HIGH

    print(
        f"Layer {i:2}: "
        f"M={M:,} "
        f"M / max_factor={ratio:,.2f}"
    )
