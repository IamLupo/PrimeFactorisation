from math import isqrt
from collections import defaultdict


# ============================================================
# SETTINGS
# ============================================================

R1 = 11
R2 = 13

WINDOW = 25

TEST_CASES = [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
    (37, 71),
    (41, 73),
    (43, 79),
    (47, 83),
    (53, 89),
    (59, 97),
    (61, 101),
    (67, 103),
    (71, 107),
]


# ============================================================
# FACTORIZATION
# ============================================================

def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def factor_pairs(n):
    """
    Return all factor pairs p <= q with p*q=n.

    This is used ONLY for the nearby numbers n+x
    and for final ground-truth evaluation of n.
    """
    result = []

    d = 2

    while d * d <= n:

        if n % d == 0:
            q = n // d

            result.append((d, q))

        d += 1

    return result


# ============================================================
# COORDINATE TRANSFORMATION
# ============================================================

def state_from_factors(p, q):
    """
    Given p*q=n:

        p = a + k*r1
        q = b + l*r2

    Return the complete coordinate state.
    """

    a = p % R1
    k = p // R1

    b = q % R2
    l = q // R2

    K = k * l

    T = (p * q) // (R1 * R2)

    E = T - K

    S = p + q

    return {
        "p": p,
        "q": q,
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "K": K,
        "T": T,
        "E": E,
        "S": S,
    }


# ============================================================
# NEIGHBOR SIGNATURE
# ============================================================

def neighbor_signature(n):
    """
    Construct a signature using ONLY factorizations of
    neighboring integers n+x, x != 0.

    We deliberately do not factor n here.
    """

    signature = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        pairs = factor_pairs(nx)

        #
        # Only semiprime neighbors are useful here:
        # exactly one factor pair with both factors prime.
        #
        semiprime_states = []

        for p, q in pairs:

            if not is_prime(p):
                continue

            if not is_prime(q):
                continue

            state = state_from_factors(p, q)

            semiprime_states.append(state)

        if not semiprime_states:
            continue

        #
        # Preserve the complete nearby state.
        #
        signature.append(
            (
                x,
                tuple(
                    (
                        state["a"],
                        state["b"],
                        state["k"],
                        state["l"],
                        state["K"],
                        state["T"],
                        state["E"],
                    )
                    for state in semiprime_states
                )
            )
        )

    return tuple(signature)


# ============================================================
# EXTRACT NEIGHBOR FEATURES
# ============================================================

def extract_features(signature):
    """
    Convert the neighborhood into several simple
    feature sets.

    These are the quantities used to predict K of n.
    """

    K_values = []
    E_values = []
    T_values = []

    KT_values = []
    KE_values = []
    TE_values = []

    for x, states in signature:

        for a, b, k, l, K, T, E in states:

            K_values.append(K)
            E_values.append(E)
            T_values.append(T)

            KT_values.append((K, T))
            KE_values.append((K, E))
            TE_values.append((T, E))

    return {
        "K": sorted(set(K_values)),
        "E": sorted(set(E_values)),
        "T": sorted(set(T_values)),
        "KT": sorted(set(KT_values)),
        "KE": sorted(set(KE_values)),
        "TE": sorted(set(TE_values)),
    }


# ============================================================
# PREDICTION RULES
# ============================================================

def candidate_K_from_neighbor_Ks(features):
    """
    Baseline:

    Any K observed in a nearby semiprime is considered
    a candidate for the target K.
    """

    return set(features["K"])


def candidate_K_from_E_T(features, n):
    """
    From

        E_x = T_x - K_x

    recover every observed K_x.

    This is intentionally just a consistency baseline.
    """

    candidates = set()

    for T, E in features["TE"]:
        candidates.add(T - E)

    return candidates


# ============================================================
# SIGNATURE HAMMING DISTANCE
# ============================================================

def signature_distance(sig_a, sig_b):
    """
    Compare two neighborhood signatures.

    We construct sets of exact observations:

        (offset, a, b, k, l, K, T, E)

    and count symmetric differences.
    """

    A = set()

    for x, states in sig_a:
        for state in states:
            A.add((x,) + state)

    B = set()

    for x, states in sig_b:
        for state in states:
            B.add((x,) + state)

    return len(A ^ B)


# ============================================================
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 13")
print()

print(
    "n\tp\tq\t"
    "true_K\ttrue_T\ttrue_E\t"
    "neighbor_count\t"
    "unique_neighbor_K\t"
    "true_K_seen\t"
    "K_rank\t"
    "candidate_reduction_pct\t"
    "same_signature_candidates"
)

#
# Build signatures without using the factors of the target.
#
signatures = {}

for p, q in TEST_CASES:

    n = p * q

    signatures[n] = neighbor_signature(n)


for p, q in TEST_CASES:

    n = p * q

    #
    # Ground truth is used ONLY here for evaluation.
    #
    true_state = state_from_factors(p, q)

    true_K = true_state["K"]
    true_T = true_state["T"]
    true_E = true_state["E"]

    sig = signatures[n]

    features = extract_features(sig)

    candidates_1 = candidate_K_from_neighbor_Ks(
        features
    )

    candidates_2 = candidate_K_from_E_T(
        features,
        n
    )

    candidates = (
        candidates_1
        | candidates_2
    )

    sorted_candidates = sorted(candidates)

    if true_K in sorted_candidates:
        K_rank = (
            sorted_candidates.index(true_K) + 1
        )
        true_K_seen = True
    else:
        K_rank = -1
        true_K_seen = False

    #
    # Candidate-space reduction relative to a
    # deliberately broad range [0, T].
    #
    #
    # The true K satisfies 0 <= K <= T.
    #
    broad_space = true_T + 1

    if broad_space > 0:
        reduction = (
            100.0
            * (
                1.0
                - len(candidates) / broad_space
            )
        )
    else:
        reduction = 0.0

    #
    # Number of nearby semiprime observations.
    #
    neighbor_count = sum(
        len(states)
        for _, states in sig
    )

    #
    # Find other tested numbers with identical
    # neighborhood signatures.
    #
    same_signature_candidates = []

    for other_n, other_sig in signatures.items():

        if other_n == n:
            continue

        if other_sig == sig:
            same_signature_candidates.append(
                other_n
            )

    print(
        f"{n}\t"
        f"{p}\t"
        f"{q}\t"
        f"{true_K}\t"
        f"{true_T}\t"
        f"{true_E}\t"
        f"{neighbor_count}\t"
        f"{len(candidates)}\t"
        f"{true_K_seen}\t"
        f"{K_rank}\t"
        f"{reduction:.2f}\t"
        f"{len(same_signature_candidates)}"
    )


print()
print("NEIGHBOR SIGNATURE DETAILS")
print("-" * 120)

for p, q in TEST_CASES:

    n = p * q

    true_state = state_from_factors(p, q)

    sig = signatures[n]

    features = extract_features(sig)

    print()
    print(
        f"n={n}, "
        f"true K={true_state['K']}, "
        f"T={true_state['T']}, "
        f"E={true_state['E']}"
    )

    print(
        "Observed neighboring K values:",
        features["K"]
    )

    print(
        "Observed neighboring E values:",
        features["E"]
    )

    print(
        "Observed neighboring T values:",
        features["T"]
    )


print()
print("PAIRWISE NEIGHBOR SIGNATURE DISTANCES")
print("-" * 120)

for i in range(len(TEST_CASES)):

    p1, q1 = TEST_CASES[i]
    n1 = p1 * q1

    for j in range(i + 1, len(TEST_CASES)):

        p2, q2 = TEST_CASES[j]
        n2 = p2 * q2

        distance = signature_distance(
            signatures[n1],
            signatures[n2]
        )

        print(
            f"{n1}\t"
            f"{n2}\t"
            f"{distance}"
        )


print()
print("FINISHED EXPERIMENT 13")

