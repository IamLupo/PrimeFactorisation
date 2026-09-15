import sympy
import math
from collections import Counter


# ============================================================
# SETTINGS
# ============================================================

MODULI_LAYERS = [
    [3, 5, 7, 11, 13],
    [3, 5, 7, 11, 13, 17],
    [3, 5, 7, 11, 13, 17, 19],
    [3, 5, 7, 11, 13, 17, 19, 23],
    [3, 5, 7, 11, 13, 17, 19, 23, 29],
    [3, 5, 7, 11, 13, 17, 19, 23, 29, 31],
]

TEST_RANGES = [
    (10_000, 100_000),
    (100_000, 1_000_000),
    (1_000_000, 10_000_000),
]

TRIALS = 100


# ============================================================
# PRECOMPUTE MODULUS PRODUCTS
# ============================================================

LAYER_DATA = []

for moduli in MODULI_LAYERS:
    LAYER_DATA.append(
        (
            moduli,
            math.prod(moduli)
        )
    )


# ============================================================
# COUNT SURVIVING PRIME PAIRS
# ============================================================

def find_candidates(n, primes, modulus, factor_high):
    """
    For every prime p <= sqrt(n), derive

        q = n * p^(-1) mod modulus

    and test whether the resulting q is a prime inside the
    tested factor range.

    No p*q == n test is used to construct the candidates.
    """

    candidates = set()

    for p_candidate in primes:

        if math.gcd(p_candidate, modulus) != 1:
            continue

        # p*q = n (mod M)
        # q = n*p^-1 (mod M)
        q_residue = (
            n * pow(p_candidate, -1, modulus)
        ) % modulus

        # If M > factor_high, there can be at most one
        # positive q in [2, factor_high] for this residue.
        if modulus > factor_high:

            q_candidate = q_residue

            if not (2 <= q_candidate <= factor_high):
                continue

            if not sympy.isprime(q_candidate):
                continue

            candidates.add(
                tuple(sorted((p_candidate, q_candidate)))
            )

        else:

            # General case: there can be multiple q values
            # in the same residue class.
            if q_residue < 2:
                k = (2 - q_residue + modulus - 1) // modulus
            else:
                k = 0

            q_candidate = q_residue + k * modulus

            while q_candidate <= factor_high:

                if q_candidate >= 2 and sympy.isprime(q_candidate):

                    candidates.add(
                        tuple(sorted(
                            (p_candidate, q_candidate)
                        ))
                    )

                q_candidate += modulus

    return candidates


# ============================================================
# MAIN
# ============================================================

print("=" * 100)
print("MODULAR FINGERPRINT SCALING EXPERIMENT")
print("=" * 100)
print()

print("Layers:")
for i, (moduli, modulus) in enumerate(LAYER_DATA, 1):
    print(
        f"  Layer {i:2}: "
        f"M = {modulus:,} "
        f"moduli = {moduli}"
    )

print()


# ============================================================
# SCALE TESTS
# ============================================================

for factor_low, factor_high in TEST_RANGES:

    print()
    print("=" * 100)
    print(
        f"FACTOR RANGE: "
        f"{factor_low:,} - {factor_high:,}"
    )
    print("=" * 100)
    print()

    # One list per layer.
    layer_counts = [
        []
        for _ in LAYER_DATA
    ]

    first_unique = []

    for trial in range(1, TRIALS + 1):

        # ----------------------------------------------------
        # Generate semiprime
        # ----------------------------------------------------

        p = sympy.randprime(
            factor_low,
            factor_high
        )

        q = sympy.randprime(
            factor_low,
            factor_high
        )

        n = p * q

        actual_pair = tuple(sorted((p, q)))

        sqrt_n = math.isqrt(n)

        primes = list(
            sympy.primerange(
                2,
                sqrt_n + 1
            )
        )

        counts_this_trial = []
        unique_at = None

        # ----------------------------------------------------
        # Test each layer
        # ----------------------------------------------------

        for layer_index, (moduli, modulus) in enumerate(
            LAYER_DATA,
            start=1
        ):

            candidates = find_candidates(
                n,
                primes,
                modulus,
                factor_high
            )

            count = len(candidates)

            # Store result
            layer_counts[layer_index - 1].append(count)
            counts_this_trial.append(count)

            # Sanity check
            if actual_pair not in candidates:

                print()
                print("=" * 100)
                print("ERROR")
                print("=" * 100)
                print(f"p = {p}")
                print(f"q = {q}")
                print(f"n = {n}")
                print(f"Layer = {layer_index}")
                print(f"M = {modulus}")
                print(
                    "Actual factor pair disappeared."
                )
                print()

                raise RuntimeError(
                    "Actual factor pair missing"
                )

            if unique_at is None and count == 1:
                unique_at = layer_index

        first_unique.append(unique_at)

        # ----------------------------------------------------
        # Keep output small
        # ----------------------------------------------------

        if trial % 10 == 0:

            print(
                f"trial {trial:3}: "
                f"counts={counts_this_trial} "
                f"first_unique={unique_at}"
            )

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("-" * 100)
    print("SUMMARY")
    print("-" * 100)
    print()

    for layer_index, ((moduli, modulus), counts) in enumerate(
        zip(LAYER_DATA, layer_counts),
        start=1
    ):

        unique_count = sum(
            count == 1
            for count in counts
        )

        multi_count = sum(
            count > 1
            for count in counts
        )

        zero_count = sum(
            count == 0
            for count in counts
        )

        print(
            f"Layer {layer_index:2}: "
            f"M={modulus:,}"
        )

        print(
            f"    exactly 1 : {unique_count:3}/{TRIALS}"
        )

        print(
            f"    >1       : {multi_count:3}/{TRIALS}"
        )

        print(
            f"    0        : {zero_count:3}/{TRIALS}"
        )

        print(
            f"    min      : {min(counts)}"
        )

        print(
            f"    max      : {max(counts)}"
        )

        print(
            f"    average  : "
            f"{sum(counts) / len(counts):.4f}"
        )

        print()

    # ========================================================
    # FIRST UNIQUE DISTRIBUTION
    # ========================================================

    print("-" * 100)
    print("FIRST UNIQUE LAYER")
    print("-" * 100)

    distribution = Counter(first_unique)

    for layer in sorted(
        distribution,
        key=lambda x: (
            999 if x is None else x
        )
    ):

        amount = distribution[layer]

        if layer is None:
            label = "never unique"
        else:
            label = f"layer {layer}"

        print(
            f"{label:16} : "
            f"{amount:3}/{TRIALS}"
        )

    print()