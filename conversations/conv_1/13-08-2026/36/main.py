#!/usr/bin/env python3

import math
import random
import time
from sympy import isprime

# ============================================================
# KAPPA EXPERIMENT 36
# N-ONLY CRT RESIDUE RECOVERY
#
# GOAL:
#   Given ONLY n, derive all possible factor-residue pairs
#   modulo the auxiliary moduli, combine them with generalized
#   CRT, and attempt to lift the resulting p-residue classes
#   back into the prime interval.
#
# NO CSV OUTPUT
# ============================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Prefixes to test.
PREFIXES = [3, 4, 5, 6, 7]

TARGETS = 12

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

# Keep this experiment compact.
PRINT_RESIDUE_CLASSES = 5


# ============================================================
# MODULI
# ============================================================

MODULI = [r * r + 3 for r in R_VALUES]


# ============================================================
# PRIME GENERATION
# ============================================================

def random_prime():
    while True:
        x = random.randrange(PRIME_LOW, PRIME_HIGH)
        if x % 2 == 0:
            x += 1

        if isprime(x):
            return x


def make_targets(count):
    targets = []

    while len(targets) < count:
        p = random_prime()
        q = random_prime()

        if p == q:
            continue

        n = p * q

        targets.append((p, q, n))

    return targets


# ============================================================
# GENERALIZED CRT
# ============================================================

def crt_pair(r1, m1, r2, m2):
    """
    Solve:

        x = r1 (mod m1)
        x = r2 (mod m2)

    Returns (r, lcm(m1,m2)), or None if incompatible.
    """

    g = math.gcd(m1, m2)

    if (r2 - r1) % g != 0:
        return None

    m1g = m1 // g
    m2g = m2 // g

    # Solve:
    #
    # r1 + m1*t = r2 (mod m2)
    #
    rhs = (r2 - r1) // g

    inv = pow(m1g, -1, m2g)

    t = (rhs * inv) % m2g

    modulus = m1 * m2g
    residue = (r1 + m1 * t) % modulus

    return residue, modulus


def merge_factor_pair_state(state, local_pairs, modulus):
    """
    state:
        set of (p_residue, q_residue)

    local_pairs:
        possible (p_residue, q_residue) modulo `modulus`

    Returns:
        (new_state, new_modulus)
    """

    if state is None:
        return set(local_pairs), modulus

    old_modulus = state[1]
    old_pairs = state[0]

    new_modulus = math.lcm(old_modulus, modulus)

    merged = set()

    for rp, rq in old_pairs:
        for lp, lq in local_pairs:

            p_merge = crt_pair(
                rp,
                old_modulus,
                lp,
                modulus,
            )

            if p_merge is None:
                continue

            q_merge = crt_pair(
                rq,
                old_modulus,
                lq,
                modulus,
            )

            if q_merge is None:
                continue

            pr, pm = p_merge
            qr, qm = q_merge

            if pm != new_modulus or qm != new_modulus:
                raise RuntimeError("CRT modulus inconsistency")

            merged.add((pr, qr))

    return merged, new_modulus


# ============================================================
# LOCAL N-ONLY FACTOR RESIDUE ENUMERATION
# ============================================================

def local_factor_pairs_from_n(n, m):
    """
    IMPORTANT:
    This function receives only n and m.

    It does NOT use p or q.

    Enumerate all residue pairs (a,b) satisfying:

        a*b = n (mod m)

    No invertibility assumption is made.
    """

    n_mod = n % m

    pairs = []

    for a in range(m):
        for b in range(m):
            if (a * b) % m == n_mod:
                pairs.append((a, b))

    return pairs


# ============================================================
# INTEGER LIFT
# ============================================================

def lifts_in_interval(residue, modulus):
    """
    Return all integers x in [PRIME_LOW, PRIME_HIGH]
    satisfying:

        x = residue (mod modulus)
    """

    if modulus <= 0:
        return []

    residue %= modulus

    if residue == 0:
        first = (
            (PRIME_LOW + modulus - 1)
            // modulus
        ) * modulus
    else:
        k = (
            PRIME_LOW - residue + modulus - 1
        ) // modulus

        first = residue + k * modulus

    if first > PRIME_HIGH:
        return []

    count = (
        (PRIME_HIGH - first) // modulus
    ) + 1

    return [
        first + i * modulus
        for i in range(count)
    ]


def interval_lift_count(residue, modulus):
    """
    Fast count without constructing the list.
    """

    residue %= modulus

    if residue == 0:
        first = (
            (PRIME_LOW + modulus - 1)
            // modulus
        ) * modulus
    else:
        first = residue + (
            (PRIME_LOW - residue + modulus - 1)
            // modulus
        ) * modulus

    if first > PRIME_HIGH:
        return 0

    return (
        (PRIME_HIGH - first) // modulus
    ) + 1


# ============================================================
# N-ONLY RECOVERY
# ============================================================

def recover_from_n(n, prefix):
    """
    Recover candidate factor pairs using ONLY n.

    Steps:

      1. For each modulus:
             a*b = n mod m

      2. Combine all possible residue pairs by generalized CRT.

      3. Once the combined modulus becomes large enough,
         lift candidate p residues into the prime interval.

      4. Check exact divisibility:
             n % p == 0
         and verify q = n//p.

    Returns a diagnostic dictionary.
    """

    state = None

    local_counts = []
    combined_counts = []

    for i in range(prefix):

        m = MODULI[i]

        local_pairs = local_factor_pairs_from_n(n, m)

        local_counts.append(
            (m, len(local_pairs))
        )

        state = merge_factor_pair_state(
            state,
            local_pairs,
            m,
        )

        pair_set, combined_modulus = state

        combined_counts.append(
            (combined_modulus, len(pair_set))
        )

    pair_set, combined_modulus = state

    # --------------------------------------------------------
    # Lift only p-residue classes.
    #
    # q is not independently brute-forced.
    # Once p is known:
    #
    #       q = n // p
    #
    # --------------------------------------------------------

    exact_factor_candidates = set()

    prime_lift_candidates = 0

    for rp, rq in pair_set:

        p_lifts = lifts_in_interval(
            rp,
            combined_modulus,
        )

        prime_lift_candidates += len(p_lifts)

        for p in p_lifts:

            if not isprime(p):
                continue

            if n % p != 0:
                continue

            q = n // p

            if q < PRIME_LOW or q > PRIME_HIGH:
                continue

            if not isprime(q):
                continue

            # Check BOTH CRT residue vectors.
            if q % combined_modulus != rq:
                continue

            if p * q != n:
                continue

            pair = tuple(sorted((p, q)))

            exact_factor_candidates.add(pair)

    return {
        "local_counts": local_counts,
        "combined_counts": combined_counts,
        "combined_modulus": combined_modulus,
        "crt_pair_count": len(pair_set),
        "prime_lift_candidates": prime_lift_candidates,
        "factor_candidates": sorted(
            exact_factor_candidates
        ),
        "residue_pairs": pair_set,
    }


# ============================================================
# CONTROL: TARGET SIGNATURE FROM KNOWN FACTORS
# ============================================================

def factor_residue_vector(p, prefix):
    """
    CONTROL ONLY.

    Used after recovery to confirm that the true factors
    belong to one of the n-derived CRT classes.
    """

    vector = []

    for i in range(prefix):
        m = MODULI[i]

        vector.append(
            (
                p % m,
            )
        )

    return tuple(vector)


def verify_target_in_crt_classes(p, q, residue_pairs, modulus):
    rp = p % modulus
    rq = q % modulus

    if (rp, rq) in residue_pairs:
        return True

    if (rq, rp) in residue_pairs:
        return True

    return False


# ============================================================
# FORMAT
# ============================================================

def print_small_list(values, limit=PRINT_RESIDUE_CLASSES):

    if not values:
        print("    none")
        return

    shown = values[:limit]

    for value in shown:
        print(f"    {value}")

    if len(values) > limit:
        print(
            f"    ... {len(values) - limit} more"
        )


# ============================================================
# MAIN
# ============================================================

def main():

    start_total = time.time()

    print("=" * 78)
    print("KAPPA EXPERIMENT 36")
    print("N-ONLY CRT RESIDUE RECOVERY")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime interval   = [{PRIME_LOW:,}, {PRIME_HIGH:,}]")
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")
    print()

    # --------------------------------------------------------
    # Modulus inventory
    # --------------------------------------------------------

    print("-" * 78)
    print("MODULUS INVENTORY")
    print("-" * 78)

    for r, m in zip(R_VALUES, MODULI):

        units = sum(
            1
            for x in range(1, m)
            if math.gcd(x, m) == 1
        )

        print(
            f"r={r:3d} "
            f"m={m:7d} "
            f"units={units:7d}"
        )

    print()

    # --------------------------------------------------------
    # Targets
    # --------------------------------------------------------

    print("-" * 78)
    print("GENERATING TARGETS")
    print("-" * 78)

    targets = make_targets(TARGETS)

    for i, (p, q, n) in enumerate(
        targets,
        start=1,
    ):
        print(
            f"target {i:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    print()

    # --------------------------------------------------------
    # Main experiment
    # --------------------------------------------------------

    global_summary = {
        prefix: {
            "unique": 0,
            "correct": 0,
            "not_found": 0,
            "ambiguous": 0,
        }
        for prefix in PREFIXES
    }

    for target_index, (true_p, true_q, n) in enumerate(
        targets,
        start=1,
    ):

        print("=" * 78)
        print(
            f"TARGET {target_index}: "
            f"n={n}"
        )
        print("=" * 78)

        print(
            f"true factors (CONTROL ONLY): "
            f"{true_p}, {true_q}"
        )

        print()

        for prefix in PREFIXES:

            t0 = time.time()

            result = recover_from_n(
                n,
                prefix,
            )

            elapsed = time.time() - t0

            combined_modulus = result[
                "combined_modulus"
            ]

            crt_pair_count = result[
                "crt_pair_count"
            ]

            lift_count = result[
                "prime_lift_candidates"
            ]

            factors = result[
                "factor_candidates"
            ]

            contains_true = (
                tuple(sorted((true_p, true_q)))
                in factors
            )

            if len(factors) == 0:
                status = "NO FACTOR FOUND"

                global_summary[prefix][
                    "not_found"
                ] += 1

            elif len(factors) == 1:

                if contains_true:
                    status = "UNIQUE / CORRECT"

                    global_summary[prefix][
                        "unique"
                    ] += 1

                    global_summary[prefix][
                        "correct"
                    ] += 1
                else:
                    status = "UNIQUE / WRONG"

                    global_summary[prefix][
                        "unique"
                    ] += 1
            else:
                status = (
                    f"AMBIGUOUS ({len(factors)})"
                )

                global_summary[prefix][
                    "ambiguous"
                ] += 1

            print()
            print(
                f"PREFIX {prefix:2d}"
            )
            print("-" * 78)

            print(
                f"last modulus             = "
                f"{MODULI[prefix-1]}"
            )

            print(
                f"combined CRT modulus     = "
                f"{combined_modulus:,}"
            )

            print(
                f"CRT residue-pair classes = "
                f"{crt_pair_count:,}"
            )

            print(
                f"prime interval lifts     = "
                f"{lift_count:,}"
            )

            print(
                f"exact factor candidates  = "
                f"{len(factors)}"
            )

            print(
                f"status                   = "
                f"{status}"
            )

            print(
                f"true pair in candidates  = "
                f"{contains_true}"
            )

            # The local counts are useful because they show
            # how much information comes from each modulus.
            print()
            print("LOCAL n-ONLY SOLUTION COUNTS")

            for m, count in result[
                "local_counts"
            ]:
                print(
                    f"  m={m:6d} "
                    f"solutions={count:6d}"
                )

            print()
            print("CRT GROWTH")

            for j, (
                cm,
                cc,
            ) in enumerate(
                result["combined_counts"],
                start=1,
            ):

                print(
                    f"  through modulus "
                    f"{MODULI[j-1]:6d}: "
                    f"M={cm:9,} "
                    f"classes={cc:9,}"
                )

            if len(factors) > 0:
                print()
                print(
                    "EXACT FACTOR CANDIDATES"
                )

                for pair in factors[:10]:
                    print(
                        f"  {pair}"
                    )

            print(
                f"\nelapsed = {elapsed:.3f}s"
            )

    # --------------------------------------------------------
    # Global summary
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("GLOBAL SUMMARY")
    print("=" * 78)

    print()
    print(
        f"{'prefix':>8s}"
        f"{'unique':>10s}"
        f"{'correct':>10s}"
        f"{'ambiguous':>12s}"
        f"{'none':>10s}"
    )

    print("-" * 78)

    for prefix in PREFIXES:

        row = global_summary[prefix]

        print(
            f"{prefix:8d}"
            f"{row['unique']:10d}"
            f"{row['correct']:10d}"
            f"{row['ambiguous']:12d}"
            f"{row['not_found']:10d}"
        )

    print()
    print("=" * 78)
    print("INTERPRETATION")
    print("=" * 78)

    print(
        """
This experiment is deliberately different from the previous
signature-collision experiments.

The recovery stage receives ONLY:

    n
    and the auxiliary moduli.

For each modulus m it solves:

    a*b = n (mod m)

for every residue pair (a,b).

Those residue pairs are then combined with generalized CRT.

No target p or q is used to construct the candidate classes.

Once the combined modulus is large enough, the experiment lifts
the possible p residues into the allowed prime interval.

For every lifted p it computes:

    q = n / p

and accepts the pair only when:

    p*q = n
    p and q are prime
    q is in the interval
    q has the required CRT residue.

This produces four useful outcomes:

    UNIQUE / CORRECT
        n-only residue information recovers the actual factor pair
        uniquely inside the tested interval.

    AMBIGUOUS
        several exact factor candidates survive.

    NO FACTOR FOUND
        the tested residue system does not produce a factor in the
        chosen interval.

    UNIQUE / WRONG
        a diagnostic failure showing that the residue constraints
        selected a different pair.

The critical result is whether increasing the CRT prefix changes:

    many residue classes
        ->
    few residue classes
        ->
    one exact factor pair.

If this occurs repeatedly, the experiment has crossed from
"signature collision analysis" into an actual n-only candidate
recovery mechanism.

IMPORTANT:
A successful UNIQUE / CORRECT result is still conditional on the
finite prime interval tested here. It is not by itself a proof of
a general factoring algorithm.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.time() - start_total:.3f}s"
    )

    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    import time

    main()

