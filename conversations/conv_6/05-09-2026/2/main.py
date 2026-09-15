import math
import time
import sympy


# ============================================================
# START EXPERIMENT 85
# ============================================================
#
# MULTI-LEVEL CANDIDATE PROPAGATION
#
# Corrected experiment.
#
# Stage 0:
#
#   - r1,r2 are SMALL
#   - a,b are exhaustively searched
#   - K is generated from
#
#         K = floor(n/(r1*r2)) - E
#
#   - divisor pairs (k,l) of K are tested
#
# Stage 1+:
#
#   - candidates are propagated to larger r1,r2
#   - no new a,b search is performed
#   - each old candidate gives the new representation exactly
#
# Final stage:
#
#   - p*q == n is checked
#
# Purpose:
#
#   Does a cheap small-r candidate cloud collapse as r grows?
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

LEVELS = 7

# Search this many E values.
#
# E is normally O(k+l).
# For balanced factors:
#
#     k,l ~ sqrt(K)
#
# so E is roughly O(sqrt(K)).
#
# This is deliberately a TEST parameter, not assumed
# to be an algorithmic solution.
#
E_LIMIT_FACTOR = 3.0

PRINT_EXAMPLES = 10


# ------------------------------------------------------------
# Generate semiprime
# ------------------------------------------------------------

def generate_semiprime(target):

    root = math.isqrt(target)

    p = int(
        sympy.randprime(
            max(3, int(root * 0.85)),
            int(root * 1.15),
        )
    )

    q_est = target // p

    q = int(
        sympy.randprime(
            max(3, int(q_est * 0.90)),
            int(q_est * 1.10) + 100,
        )
    )

    return p, q, p * q


# ------------------------------------------------------------
# Carry state
# ------------------------------------------------------------

def compute_state(k, l, a, b, r1, r2):

    R = r1 * r2

    K = k * l

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    E = c1 + c2 + c3

    return K, E, c1, c2, c3


# ------------------------------------------------------------
# Initial candidate generation
# ------------------------------------------------------------

def initial_candidates(n, r1, r2):

    R = r1 * r2
    Q = n // R

    # For balanced k,l:
    #
    #     E < k + l + 3
    #
    # and approximately k,l ~ sqrt(Q).
    #
    # Give ourselves a generous experimental range.
    e_limit = int(
        E_LIMIT_FACTOR * math.isqrt(Q)
    ) + 20

    candidates = {}

    tested_E = 0
    tested_K = 0
    tested_divisor_pairs = 0
    tested_ab = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        tested_E += 1
        tested_K += 1

        # Factor K to get all k*l = K possibilities.
        #
        # This is acceptable for the SMALL-r experiment.
        factorization = sympy.factorint(K)

        divisors = sympy.divisors(K)

        for k in divisors:

            l = K // k

            tested_divisor_pairs += 1

            for a in range(r1):

                for b in range(r2):

                    tested_ab += 1

                    K2, E2, c1, c2, c3 = compute_state(
                        k,
                        l,
                        a,
                        b,
                        r1,
                        r2,
                    )

                    if K2 != K:
                        continue

                    if E2 != E:
                        continue

                    key = (
                        k,
                        l,
                        a,
                        b,
                    )

                    candidates[key] = {
                        "k": k,
                        "l": l,
                        "a": a,
                        "b": b,
                        "K": K2,
                        "E": E2,
                        "c1": c1,
                        "c2": c2,
                        "c3": c3,
                    }

    elapsed = time.perf_counter() - start

    stats = {
        "e_limit": e_limit,
        "tested_E": tested_E,
        "tested_K": tested_K,
        "tested_divisor_pairs": tested_divisor_pairs,
        "tested_ab": tested_ab,
        "time": elapsed,
    }

    return list(candidates.values()), stats


# ------------------------------------------------------------
# Propagate
# ------------------------------------------------------------

def propagate(
    candidates,
    n,
    old_r1,
    old_r2,
    new_r1,
    new_r2,
):

    Q = n // (new_r1 * new_r2)

    out = {}

    for candidate in candidates:

        k = candidate["k"]
        l = candidate["l"]
        a = candidate["a"]
        b = candidate["b"]

        # Candidate factors at the OLD level.
        p = old_r1 * k + a
        q = old_r2 * l + b

        # New representation.
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

    return list(out.values())


# ------------------------------------------------------------
# Exact verification
# ------------------------------------------------------------

def verify(candidates, r1, r2, n):

    valid = []

    for candidate in candidates:

        p = (
            r1 * candidate["k"]
            + candidate["a"]
        )

        q = (
            r2 * candidate["l"]
            + candidate["b"]
        )

        if p * q == n:

            item = dict(candidate)

            item["p"] = p
            item["q"] = q

            valid.append(item)

    return valid


# ------------------------------------------------------------
# Is the true state present?
# ------------------------------------------------------------

def true_state_present(
    candidates,
    p,
    q,
    r1,
    r2,
):

    k, a = divmod(p, r1)
    l, b = divmod(q, r2)

    for c in candidates:

        if (
            c["k"] == k
            and c["l"] == l
            and c["a"] == a
            and c["b"] == b
        ):
            return True

    return False


# ------------------------------------------------------------
# Print examples
# ------------------------------------------------------------

def print_examples(candidates, r1, r2):

    if not candidates:

        print("    examples: none")
        return

    print("    examples:")

    for c in candidates[:PRINT_EXAMPLES]:

        p = r1 * c["k"] + c["a"]
        q = r2 * c["l"] + c["b"]

        print(
            "      "
            f"k={c['k']} "
            f"l={c['l']} "
            f"a={c['a']} "
            f"b={c['b']} "
            f"K={c['K']} "
            f"E={c['E']} "
            f"p={p} "
            f"q={q}"
        )


# ------------------------------------------------------------
# Generate next prime approximately 2x
# ------------------------------------------------------------

def next_prime_2x(r):

    return int(
        sympy.nextprime(
            2 * r - 1
        )
    )


# ------------------------------------------------------------
# Run case
# ------------------------------------------------------------

def run_case(target):

    p, q, n = generate_semiprime(target)

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

    candidates, stats = initial_candidates(
        n,
        r1,
        r2,
    )

    true_alive = true_state_present(
        candidates,
        p,
        q,
        r1,
        r2,
    )

    print()
    print(f"E limit            = {stats['e_limit']:,}")
    print(f"E values tested    = {stats['tested_E']:,}")
    print(f"K values tested    = {stats['tested_K']:,}")
    print(
        "divisor pairs      = "
        f"{stats['tested_divisor_pairs']:,}"
    )
    print(
        "a,b combinations   = "
        f"{stats['tested_ab']:,}"
    )
    print(f"surviving          = {len(candidates):,}")
    print(f"TRUE candidate alive? {true_alive}")
    print(f"time               = {stats['time']:.6f} s")

    print_examples(
        candidates,
        r1,
        r2,
    )

    # --------------------------------------------------------
    # PROPAGATION
    # --------------------------------------------------------

    for level in range(1, LEVELS + 1):

        new_r1 = next_prime_2x(r1)
        new_r2 = next_prime_2x(r2)

        old_count = len(candidates)

        start = time.perf_counter()

        candidates = propagate(
            candidates,
            n,
            r1,
            r2,
            new_r1,
            new_r2,
        )

        elapsed = time.perf_counter() - start

        true_alive = true_state_present(
            candidates,
            p,
            q,
            new_r1,
            new_r2,
        )

        print()
        print(f"LEVEL {level}")
        print(f"r1={new_r1}")
        print(f"r2={new_r2}")
        print()
        print(
            f"previous candidates = "
            f"{old_count:,}"
        )
        print(
            f"surviving           = "
            f"{len(candidates):,}"
        )
        print(
            f"reduction           = "
            f"{old_count / max(1, len(candidates)):.3f}x"
        )
        print(
            f"TRUE candidate alive? "
            f"{true_alive}"
        )
        print(
            f"time                = "
            f"{elapsed:.6f} s"
        )

        print_examples(
            candidates,
            new_r1,
            new_r2,
        )

        if not true_alive:

            print()
            print(
                "WARNING: TRUE CANDIDATE DIED "
                "AT THIS LEVEL."
            )

        if not candidates:

            print()
            print("ALL CANDIDATES DIED.")
            break

        r1 = new_r1
        r2 = new_r2

    # --------------------------------------------------------
    # FINAL EXACT TEST
    # --------------------------------------------------------

    print()
    print("-" * 78)

    valid = verify(
        candidates,
        r1,
        r2,
        n,
    )

    if valid:

        print(
            "FINAL RESULT: "
            "TRUE FACTORIZATION FOUND"
        )

        seen = set()

        for item in valid:

            pair = tuple(
                sorted(
                    (
                        item["p"],
                        item["q"],
                    )
                )
            )

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

        print(
            "FINAL RESULT: "
            "no exact factorization among survivors"
        )

    print("-" * 78)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("MULTI-LEVEL CANDIDATE PROPAGATION")
    print("EXPERIMENT 85")
    print("=" * 78)

    for target in CASE_SIZES:

        run_case(target)

    print()
    print("Experiment completed.")


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 85
# ============================================================
