import math
import random
import time

import sympy


# ============================================================
# START EXPERIMENT 84
# ============================================================
#
# MULTI-LEVEL CANDIDATE PROPAGATION
#
# Hypothesis:
#
#   Small r:
#       a,b are cheap to enumerate.
#
#   Larger r:
#       a,b become expensive,
#       but an old candidate (k,l,a,b)
#       determines p,q and therefore the new
#       (k',l',a',b') exactly.
#
# We intentionally DO NOT require p*q == n
# during the propagation stages.
#
# We only use:
#
#     floor(n / (r1*r2)) = K + E
#
# as the level filter.
#
# At the final stage, we verify p*q == n.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

CASE_SIZES = [
    10**9,
    10**12,
    10**16,
]

INITIAL_R1 = 11
INITIAL_R2 = 13

# Search window around sqrt(n)/r for k and l
K_WINDOW = 8

# Number of levels
LEVELS = 7

# Limit printed candidate examples
PRINT_EXAMPLES = 8


# ------------------------------------------------------------
# Prime helper
# ------------------------------------------------------------

def next_prime_2x(r):
    """
    Pick the first prime >= 2*r.
    """
    return int(sympy.nextprime(2 * r - 1))


# ------------------------------------------------------------
# Carry calculation
# ------------------------------------------------------------

def compute_state(k, l, a, b, r1, r2):
    """
    Compute K and E from one candidate state.
    """

    K = k * l

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    R = r1 * r2

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    E = c1 + c2 + c3

    return K, E, c1, c2, c3


# ------------------------------------------------------------
# Generate a candidate cloud at the first level
# ------------------------------------------------------------

def initial_candidates(n, r1, r2, window):
    """
    Search a small neighborhood around the expected values

        k ~ sqrt(n)/r1
        l ~ sqrt(n)/r2

    and enumerate every possible a,b.

    We KEEP candidates satisfying only

        floor(n/R) == K + E

    We deliberately do NOT test p*q == n here.
    """

    R = r1 * r2
    Q = n // R

    sqrt_n = math.isqrt(n)

    k_center = max(1, sqrt_n // r1)
    l_center = max(1, sqrt_n // r2)

    candidates = []

    raw = 0

    k_min = max(1, k_center - window)
    k_max = k_center + window

    l_min = max(1, l_center - window)
    l_max = l_center + window

    for k in range(k_min, k_max + 1):

        for l in range(l_min, l_max + 1):

            K = k * l

            for a in range(r1):

                for b in range(r2):

                    raw += 1

                    K2, E, c1, c2, c3 = compute_state(
                        k, l, a, b, r1, r2
                    )

                    if K2 + E != Q:
                        continue

                    candidates.append({
                        "k": k,
                        "l": l,
                        "a": a,
                        "b": b,
                        "K": K2,
                        "E": E,
                        "c1": c1,
                        "c2": c2,
                        "c3": c3,
                    })

    return candidates, raw


# ------------------------------------------------------------
# Propagate candidates to a larger r1,r2
# ------------------------------------------------------------

def propagate_candidates(
    candidates,
    n,
    old_r1,
    old_r2,
    new_r1,
    new_r2,
):
    """
    Transform every candidate through p and q:

        p = old_r1*k + a
        q = old_r2*l + b

    Then compute the exact residues at the new level:

        k' = p // new_r1
        a' = p %  new_r1

        l' = q // new_r2
        b' = q %  new_r2

    No factorization is performed.

    We then apply only:

        floor(n/(new_r1*new_r2)) == K' + E'
    """

    R = new_r1 * new_r2
    Q = n // R

    out = {}
    total = len(candidates)

    for candidate in candidates:

        k = candidate["k"]
        l = candidate["l"]
        a = candidate["a"]
        b = candidate["b"]

        # Recover the candidate factors implied by the state.
        p = old_r1 * k + a
        q = old_r2 * l + b

        # New level representation.
        new_k, new_a = divmod(p, new_r1)
        new_l, new_b = divmod(q, new_r2)

        if new_k <= 0 or new_l <= 0:
            continue

        K, E, c1, c2, c3 = compute_state(
            new_k,
            new_l,
            new_a,
            new_b,
            new_r1,
            new_r2,
        )

        # The only filter.
        if K + E != Q:
            continue

        key = (
            new_k,
            new_l,
            new_a,
            new_b,
        )

        out[key] = {
            "k": new_k,
            "l": new_l,
            "a": new_a,
            "b": new_b,
            "K": K,
            "E": E,
            "c1": c1,
            "c2": c2,
            "c3": c3,
        }

    return list(out.values()), total


# ------------------------------------------------------------
# Final verification
# ------------------------------------------------------------

def verify_candidates(
    candidates,
    r1,
    r2,
    n,
):
    """
    Finally check p*q == n.
    """

    valid = []

    for candidate in candidates:

        p = r1 * candidate["k"] + candidate["a"]
        q = r2 * candidate["l"] + candidate["b"]

        if p * q == n:

            item = dict(candidate)

            item["p"] = p
            item["q"] = q

            valid.append(item)

    return valid


# ------------------------------------------------------------
# Generate a balanced semiprime
# ------------------------------------------------------------

def generate_semiprime(target):
    """
    Generate two random primes near sqrt(target).
    """

    root = math.isqrt(target)

    low = max(3, int(root * 0.85))
    high = max(low + 10, int(root * 1.15))

    p = int(sympy.randprime(low, high))

    q_low = max(3, target // p)

    # Search around the corresponding q region.
    q_low = int(q_low * 0.90)
    q_high = int(q_low * 1.25) + 100

    q = int(sympy.randprime(max(3, q_low), q_high))

    return p, q, p * q


# ------------------------------------------------------------
# Print a few candidates
# ------------------------------------------------------------

def print_examples(candidates, r1, r2):

    if not candidates:
        print("    examples: none")
        return

    print("    examples:")

    for candidate in candidates[:PRINT_EXAMPLES]:

        p = r1 * candidate["k"] + candidate["a"]
        q = r2 * candidate["l"] + candidate["b"]

        print(
            "      "
            f"k={candidate['k']} "
            f"l={candidate['l']} "
            f"a={candidate['a']} "
            f"b={candidate['b']} "
            f"K={candidate['K']} "
            f"E={candidate['E']} "
            f"p*q==n? {p*q == CURRENT_N}"
        )


# ------------------------------------------------------------
# One complete test
# ------------------------------------------------------------

def run_case(target):

    global CURRENT_N

    p, q, n = generate_semiprime(target)
    CURRENT_N = n

    print()
    print("=" * 78)
    print(f"TARGET ~ {target:,}")
    print(f"n      = {n}")
    print(f"p      = {p}")
    print(f"q      = {q}")
    print(f"|p-q|  = {abs(p-q)}")
    print("=" * 78)

    # --------------------------------------------------------
    # LEVEL 0
    # --------------------------------------------------------

    r1 = INITIAL_R1
    r2 = INITIAL_R2

    print()
    print("LEVEL 0")
    print(f"r1={r1}")
    print(f"r2={r2}")

    start = time.perf_counter()

    candidates, raw = initial_candidates(
        n,
        r1,
        r2,
        K_WINDOW,
    )

    elapsed = time.perf_counter() - start

    true_count = len(
        verify_candidates(
            candidates,
            r1,
            r2,
            n,
        )
    )

    print(f"raw states       = {raw:,}")
    print(f"surviving        = {len(candidates):,}")
    print(f"true factorizers = {true_count:,}")
    print(f"time             = {elapsed:.6f} s")

    print_examples(candidates, r1, r2)

    # --------------------------------------------------------
    # PROPAGATION
    # --------------------------------------------------------

    for level in range(1, LEVELS + 1):

        new_r1 = next_prime_2x(r1)
        new_r2 = next_prime_2x(r2)

        old_count = len(candidates)

        start = time.perf_counter()

        candidates, transformed = propagate_candidates(
            candidates,
            n,
            r1,
            r2,
            new_r1,
            new_r2,
        )

        elapsed = time.perf_counter() - start

        true_count = len(
            verify_candidates(
                candidates,
                new_r1,
                new_r2,
                n,
            )
        )

        print()
        print(f"LEVEL {level}")
        print(f"r1={new_r1}")
        print(f"r2={new_r2}")
        print()
        print(f"previous candidates = {old_count:,}")
        print(f"transformed         = {transformed:,}")
        print(f"surviving           = {len(candidates):,}")
        print(f"true factorizers    = {true_count:,}")
        print(f"reduction           = {old_count / max(1, len(candidates)):.3f}x")
        print(f"time                = {elapsed:.6f} s")

        print_examples(candidates, new_r1, new_r2)

        r1 = new_r1
        r2 = new_r2

        if not candidates:
            print()
            print("ALL CANDIDATES DIED.")
            print("The propagation hypothesis failed for this case.")
            break

    # --------------------------------------------------------
    # FINAL RESULT
    # --------------------------------------------------------

    print()
    print("-" * 78)

    final_valid = verify_candidates(
        candidates,
        r1,
        r2,
        n,
    )

    if final_valid:

        print("FINAL RESULT: TRUE FACTOR FOUND")

        seen = set()

        for item in final_valid:

            pair = tuple(sorted((item["p"], item["q"])))

            if pair in seen:
                continue

            seen.add(pair)

            print(
                f"p={pair[0]}"
            )
            print(
                f"q={pair[1]}"
            )

    else:

        print("FINAL RESULT: no exact factorization among survivors")

    print("-" * 78)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("MULTI-LEVEL CANDIDATE PROPAGATION TEST")
    print("=" * 78)

    for target in CASE_SIZES:
        run_case(target)

    print()
    print("Experiment completed.")


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 84
# ============================================================
