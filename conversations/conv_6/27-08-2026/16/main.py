import sympy
import math


# ============================================================
# SETTINGS
# ============================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

TRIALS = 100

MODULI_LAYERS = [
    (
        8,
        [3, 5, 7, 11, 13, 17, 19, 23]
    ),
    (
        9,
        [3, 5, 7, 11, 13, 17, 19, 23, 29]
    ),
    (
        10,
        [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    ),
]


# ============================================================
# FIND CANDIDATES
# ============================================================

def find_candidates(n, primes, modulus):
    """
    Find unordered prime pairs (p,q) in the tested range
    satisfying

        p*q == n (mod modulus)

    The actual equality p*q == n is NOT used to construct
    the candidate set.
    """

    candidates = set()

    for p_candidate in primes:

        # Need inverse modulo M.
        if math.gcd(p_candidate, modulus) != 1:
            continue

        # p*q = n (mod M)
        #
        # q = n * p^(-1) (mod M)

        q_candidate = (
            n * pow(p_candidate, -1, modulus)
        ) % modulus

        # Because M is much larger than the factor range
        # in this experiment, q_candidate is the only possible
        # representative below FACTOR_HIGH.
        if not (
            FACTOR_LOW
            <= q_candidate
            <= FACTOR_HIGH
        ):
            continue

        if not sympy.isprime(q_candidate):
            continue

        pair = tuple(
            sorted(
                (p_candidate, q_candidate)
            )
        )

        candidates.add(pair)

    return sorted(candidates)


# ============================================================
# MAIN
# ============================================================

print("=" * 100)
print("COLLISION STRUCTURE EXPERIMENT")
print("=" * 100)

print(
    f"Factor range : "
    f"{FACTOR_LOW:,} - {FACTOR_HIGH:,}"
)

print(
    f"Trials       : {TRIALS}"
)

print()

# Store results as a list instead of a dictionary,
# avoiding the previous layer-index bug.

all_results = []


for trial in range(1, TRIALS + 1):

    # --------------------------------------------------------
    # Generate random semiprime
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

    actual_pair = tuple(
        sorted((p, q))
    )

    # Only p <= sqrt(n) needs to be considered.
    sqrt_n = math.isqrt(n)

    primes = list(
        sympy.primerange(
            2,
            sqrt_n + 1
        )
    )

    print()
    print("-" * 100)
    print(
        f"TRIAL {trial:3}"
    )
    print("-" * 100)

    print(
        f"p = {p:,}"
    )

    print(
        f"q = {q:,}"
    )

    print(
        f"n = {n:,}"
    )

    print(
        f"sqrt(n) = {sqrt_n:,}"
    )

    print()

    trial_results = []

    # --------------------------------------------------------
    # Layers 8, 9, 10
    # --------------------------------------------------------

    for layer_number, moduli in MODULI_LAYERS:

        M = math.prod(moduli)

        candidates = find_candidates(
            n,
            primes,
            M
        )

        # Sanity check
        if actual_pair not in candidates:

            print()
            print("ERROR: actual pair disappeared")
            print(
                f"Layer = {layer_number}"
            )
            print(
                f"M = {M:,}"
            )
            print(
                f"Candidates = {len(candidates)}"
            )

            raise RuntimeError(
                "Actual factor pair missing."
            )

        false_candidates = [
            pair
            for pair in candidates
            if pair != actual_pair
        ]

        # ----------------------------------------------------
        # Print layer summary
        # ----------------------------------------------------

        print(
            f"LAYER {layer_number}"
        )

        print(
            f"  M = {M:,}"
        )

        print(
            f"  candidates = {len(candidates)}"
        )

        print(
            f"  false collisions = "
            f"{len(false_candidates)}"
        )

        # ----------------------------------------------------
        # Collision details
        # ----------------------------------------------------

        collision_data = []

        for a, b in false_candidates:

            product = a * b

            difference = product - n

            # Since this candidate has the same fingerprint:
            #
            #     product = n (mod M)
            #
            # therefore:
            #
            #     difference = k*M

            if difference % M != 0:

                raise RuntimeError(
                    "Fingerprint collision did not produce "
                    "an exact multiple of M."
                )

            k = difference // M

            collision_data.append(
                (
                    a,
                    b,
                    product,
                    difference,
                    k
                )
            )

        # ----------------------------------------------------
        # Print collisions
        # ----------------------------------------------------

        for (
            a,
            b,
            product,
            difference,
            k
        ) in collision_data:

            print(
                f"    FALSE: "
                f"{a:,} * {b:,} = {product:,}"
            )

            print(
                f"            "
                f"product - n = "
                f"{difference:+,}"
            )

            print(
                f"            "
                f"k = {k:+,}"
            )

        trial_results.append(
            {
                "layer": layer_number,
                "M": M,
                "candidate_count": len(candidates),
                "false_count": len(false_candidates),
                "collisions": collision_data,
            }
        )

    all_results.append(
        {
            "trial": trial,
            "p": p,
            "q": q,
            "n": n,
            "layers": trial_results,
        }
    )


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

for layer_number, moduli in MODULI_LAYERS:

    M = math.prod(moduli)

    values = []

    for result in all_results:

        layer_result = next(
            item
            for item in result["layers"]
            if item["layer"] == layer_number
        )

        values.append(
            layer_result["candidate_count"]
        )

    unique = sum(
        value == 1
        for value in values
    )

    multi = sum(
        value > 1
        for value in values
    )

    print()
    print(
        f"LAYER {layer_number}"
    )

    print(
        f"  M                  = {M:,}"
    )

    print(
        f"  min candidates     = {min(values)}"
    )

    print(
        f"  max candidates     = {max(values)}"
    )

    print(
        f"  average candidates = "
        f"{sum(values) / len(values):.3f}"
    )

    print(
        f"  unique             = "
        f"{unique}/{TRIALS}"
    )

    print(
        f"  collisions         = "
        f"{multi}/{TRIALS}"
    )


# ============================================================
# COLLISION k SUMMARY
# ============================================================

print()
print("=" * 100)
print("FALSE-COLLISION k VALUES")
print("=" * 100)

for layer_number, moduli in MODULI_LAYERS:

    k_values = []

    for result in all_results:

        layer_result = next(
            item
            for item in result["layers"]
            if item["layer"] == layer_number
        )

        for collision in layer_result["collisions"]:

            k = collision[4]

            k_values.append(k)

    print()
    print(
        f"LAYER {layer_number}"
    )

    if not k_values:

        print(
            "  No false collisions."
        )

        continue

    print(
        f"  collisions = {len(k_values)}"
    )

    print(
        f"  min k = {min(k_values):+,}"
    )

    print(
        f"  max k = {max(k_values):+,}"
    )

    print(
        f"  unique k values = "
        f"{sorted(set(k_values))}"
    )