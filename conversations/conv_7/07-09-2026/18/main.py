from math import isqrt
from collections import defaultdict


R1 = 11
R2 = 13

WINDOW = 50

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
# BASIC NUMBER THEORY
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


def prime_factor_pairs(n):
    """
    Return all factor pairs (p,q), p <= q,
    where both p and q are prime.
    """
    result = []

    d = 2

    while d * d <= n:

        if n % d == 0:

            q = n // d

            if is_prime(d) and is_prime(q):
                result.append((d, q))

        d += 1

    return result


# ============================================================
# MODULAR STATE
# ============================================================

def make_state(p, q, n):
    """
    State:

        p = a + k*R1
        q = b + l*R2

        K = k*l
        T = floor(n/(R1*R2))
        E = T-K
    """

    a = p % R1
    k = p // R1

    b = q % R2
    l = q // R2

    K = k * l

    T = n // (R1 * R2)

    E = T - K

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
    }


# ============================================================
# BUILD NEIGHBOR STATES
# ============================================================

def neighbor_states(n):
    """
    Return all semiprime states for n+x,
    excluding x=0.

    The target factorization of n is NOT used here.
    """

    states = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        pairs = prime_factor_pairs(nx)

        for p, q in pairs:

            s = make_state(p, q, nx)

            states.append(
                {
                    "x": x,
                    **s,
                }
            )

    return states


# ============================================================
# PAIRWISE DIFFERENCE SIGNATURE
# ============================================================

def pairwise_signatures(states):
    """
    For every pair of nearby factorizations, compute:

        dK = Kx - Ky
        dT = Tx - Ty
        dE = Ex - Ey

    and verify the exact identity

        dE = dT - dK
    """

    result = []

    for i in range(len(states)):

        a = states[i]

        for j in range(i + 1, len(states)):

            b = states[j]

            dK = a["K"] - b["K"]
            dT = a["T"] - b["T"]
            dE = a["E"] - b["E"]

            identity_ok = (
                dE == dT - dK
            )

            result.append(
                {
                    "x1": a["x"],
                    "x2": b["x"],
                    "K1": a["K"],
                    "K2": b["K"],
                    "T1": a["T"],
                    "T2": b["T"],
                    "E1": a["E"],
                    "E2": b["E"],
                    "dK": dK,
                    "dT": dT,
                    "dE": dE,
                    "identity_ok": identity_ok,
                }
            )

    return result


# ============================================================
# TRANSLATION CANDIDATES
# ============================================================

def translation_candidates(states):
    """
    Treat each nearby state as a hypothetical translation
    of the target K.

    From

        E_x = T_x - K_x

    and

        E_x - E = (T_x-T) - (K_x-K)

    we get

        K = K_x + E_x - E - (T_x-T).

    Since E is unknown, every nearby state gives:

        K - E = K_x - E_x - (T_x-T).

    The important quantity is therefore:

        C_x = K_x - E_x - (T_x-T).

    If several neighbors produce the same C_x,
    they agree on the same candidate value of K-E.
    """

    candidates = defaultdict(list)

    for s in states:

        x = s["x"]

        T_target = None

        #
        # T for the target is floor(n/(R1*R2)).
        # This function receives states only, so infer the target
        # T by removing x from the neighbor n stored below.
        #
        # We attach target_T later.
        #
        T_target = s["target_T"]

        C = (
            s["K"]
            - s["E"]
            - (s["T"] - T_target)
        )

        candidates[C].append(x)

    return candidates


# ============================================================
# MAIN EXPERIMENT
# ============================================================

print("START EXPERIMENT 15")
print()

print(
    "n\tp\tq\t"
    "true_K\ttrue_T\ttrue_E\t"
    "neighbors\t"
    "pair_count\t"
    "identity_failures\t"
    "unique_C\t"
    "max_C_support\t"
    "true_K_minus_E_seen"
)

for p, q in TEST_CASES:

    n = p * q

    #
    # Ground truth ONLY for evaluation.
    #
    target = make_state(p, q, n)

    states = neighbor_states(n)

    target_T = target["T"]

    #
    # Attach target_T so the translation calculation
    # remains independent of target K/E.
    #
    for s in states:
        s["target_T"] = target_T

    pairs = pairwise_signatures(states)

    identity_failures = sum(
        not x["identity_ok"]
        for x in pairs
    )

    #
    # Compute the quantity

    #   C_x =
    #       K_x - E_x - (T_x-T)
    #
    # which, from the exact descent relation, satisfies
    #
    #   C_x = K-E
    #
    # whenever the nearby state is on the same translated
    # branch.
    #

    C_values = defaultdict(list)

    for s in states:

        C = (
            s["K"]
            - s["E"]
            - (s["T"] - target_T)
        )

        C_values[C].append(
            s["x"]
        )

    true_K_minus_E = (
        target["K"]
        - target["E"]
    )

    true_seen = (
        true_K_minus_E in C_values
    )

    max_support = (
        max(
            (len(v) for v in C_values.values()),
            default=0
        )
    )

    print(
        f"{n}\t"
        f"{p}\t"
        f"{q}\t"
        f"{target['K']}\t"
        f"{target['T']}\t"
        f"{target['E']}\t"
        f"{len(states)}\t"
        f"{len(pairs)}\t"
        f"{identity_failures}\t"
        f"{len(C_values)}\t"
        f"{max_support}\t"
        f"{true_seen}"
    )


print()
print("C-VALUE DETAILS")
print("-" * 110)

for p, q in TEST_CASES:

    n = p * q

    target = make_state(p, q, n)

    states = neighbor_states(n)

    target_T = target["T"]

    C_values = defaultdict(list)

    for s in states:

        C = (
            s["K"]
            - s["E"]
            - (s["T"] - target_T)
        )

        C_values[C].append(
            s["x"]
        )

    print()
    print(
        f"n={n} "
        f"true K={target['K']} "
        f"true E={target['E']} "
        f"true K-E={target['K'] - target['E']}"
    )

    for C in sorted(C_values):

        offsets = C_values[C]

        print(
            f"  C={C:4d} "
            f"support={len(offsets):2d} "
            f"offsets={offsets}"
        )


print()
print("PAIRWISE IDENTITY CHECK")
print("-" * 110)

total_pairs = 0
total_failures = 0

for p, q in TEST_CASES:

    n = p * q

    states = neighbor_states(n)

    pairs = pairwise_signatures(states)

    total_pairs += len(pairs)

    failures = [
        x for x in pairs
        if not x["identity_ok"]
    ]

    total_failures += len(failures)

    print(
        f"n={n}: "
        f"pairs={len(pairs)} "
        f"identity_failures={len(failures)}"
    )


print()
print(
    f"TOTAL PAIRS={total_pairs}"
)

print(
    f"TOTAL IDENTITY FAILURES={total_failures}"
)

print()
print("FINISHED EXPERIMENT 15")
