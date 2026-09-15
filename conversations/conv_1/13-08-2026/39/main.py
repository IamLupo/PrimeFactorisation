#!/usr/bin/env python3

import math
import random
import time
from collections import defaultdict

from sympy import isprime


# ==============================================================================
# KAPPA EXPERIMENT 39
# N-ONLY CRT CANDIDATE RECOVERY / EXACT FACTOR TEST
# NO CSV OUTPUT
# ==============================================================================

SEED = 20260814
random.seed(SEED)

TARGETS = 12

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

PREFIXES = [3, 4, 5, 6, 7]

# Only print the full candidate list when it is small.
PRINT_CANDIDATES_AT_MOST = 20

# Rank mode:
#   "smallest-first" checks candidates in ascending a,b order.
#   "closest" checks candidates whose two factors are closest together first.
#   "random" checks candidates in random order.
SEARCH_ORDER = "closest"


# ==============================================================================
# BASIC
# ==============================================================================

def generate_primes(low, high):
    primes = []

    # 2
    if low <= 2 < high:
        primes.append(2)

    start = max(3, low | 1)

    for n in range(start, high, 2):
        if isprime(n):
            primes.append(n)

    return primes


def factor_pair_from_population(primes):
    while True:
        p = random.choice(primes)
        q = random.choice(primes)

        if p != q:
            return p, q


# ==============================================================================
# MODULI
# ==============================================================================

# Experiment-37 / 38 family.
MODULI = [r * r + 3 for r in R_VALUES]


def lcm_many(values):
    out = 1

    for x in values:
        out = math.lcm(out, x)

    return out


PREFIX_MODULI = {}
PREFIX_M = {}

for k in PREFIXES:
    mods = MODULI[:k]

    PREFIX_MODULI[k] = mods
    PREFIX_M[k] = lcm_many(mods)


# ==============================================================================
# CRT / RESIDUE REPRESENTATION
# ==============================================================================

def residue_key(x, mods):
    return tuple(x % m for m in mods)


def crt_value_from_residues(residues, mods):
    """
    Recover the unique CRT residue modulo lcm(mods).

    The moduli are not necessarily pairwise coprime, so we use
    generalized CRT consistency explicitly.
    """

    x = 0
    current_modulus = 1

    for a, m in zip(residues, mods):

        g = math.gcd(current_modulus, m)

        if (a - x) % g != 0:
            return None

        # Solve:
        #
        # x + current_modulus * t == a (mod m)
        #
        m1 = current_modulus // g
        m2 = m // g

        if m2 == 1:
            t = 0
        else:
            rhs = (a - x) // g
            inv = pow(m1, -1, m2)
            t = (rhs * inv) % m2

        x += current_modulus * t
        current_modulus = math.lcm(current_modulus, m)

        x %= current_modulus

    return x


# ==============================================================================
# PRIME RESIDUE INDEX
# ==============================================================================

def build_prime_index(primes, mods):
    """
    key -> list of primes having that complete residue vector.
    """

    index = defaultdict(list)

    for p in primes:
        index[residue_key(p, mods)].append(p)

    return index


# ==============================================================================
# N-ONLY COMPATIBLE RESIDUE CLASSES
# ==============================================================================

def compatible_second_residue_vectors(n, mods, first_key):
    """
    Given the residue vector of p, compute the possible residue vector(s)
    for q satisfying

        p*q == n mod m

    for every modulus.

    If p is a unit modulo m, q is unique modulo m.
    """

    out = []

    for a, m in zip(first_key, mods):

        n_mod = n % m

        if math.gcd(a, m) != 1:
            return None

        q_mod = (n_mod * pow(a, -1, m)) % m

        out.append(q_mod)

    return tuple(out)


# ==============================================================================
# EXACT N-ONLY CANDIDATE GENERATION
# ==============================================================================

def build_candidate_pairs(n, primes, index, mods):
    """
    Return all unordered prime pairs (a,b) from the population satisfying
    the complete N-only CRT condition.

    The search is over residue classes, not over prime pairs.
    """

    candidates = set()

    for a_key, a_values in index.items():

        b_key = compatible_second_residue_vectors(
            n,
            mods,
            a_key
        )

        if b_key is None:
            continue

        b_values = index.get(b_key)

        if not b_values:
            continue

        for a in a_values:
            for b in b_values:

                if a == b:
                    continue

                pair = (a, b) if a < b else (b, a)
                candidates.add(pair)

    return sorted(candidates)


# ==============================================================================
# SEARCH ORDER
# ==============================================================================

def order_candidates(candidates):
    if SEARCH_ORDER == "smallest-first":
        return sorted(candidates)

    if SEARCH_ORDER == "random":
        result = list(candidates)
        random.shuffle(result)
        return result

    if SEARCH_ORDER == "closest":
        return sorted(
            candidates,
            key=lambda pair: (
                abs(pair[0] - pair[1]),
                pair[0],
                pair[1],
            )
        )

    raise ValueError(
        f"Unknown SEARCH_ORDER={SEARCH_ORDER!r}"
    )


# ==============================================================================
# EXACT RECOVERY
# ==============================================================================

def exact_recovery(n, candidates):
    """
    Test candidate products until a*b == n.

    Returns:
        found_pair
        number_checked
        index_of_true_pair_if_present
    """

    ordered = order_candidates(candidates)

    for i, (a, b) in enumerate(ordered, start=1):
        if a * b == n:
            return (a, b), i, len(ordered)

    return None, len(ordered), len(ordered)


# ==============================================================================
# TARGET GENERATION
# ==============================================================================

def generate_targets(primes, count):
    targets = []

    seen = set()

    while len(targets) < count:

        p, q = factor_pair_from_population(primes)

        key = tuple(sorted((p, q)))

        if key in seen:
            continue

        seen.add(key)

        n = p * q

        targets.append((p, q, n))

    return targets


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 39")
    print("N-ONLY CRT CANDIDATE RECOVERY / EXACT FACTOR TEST")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print()

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime interval   = [{PRIME_LOW:,}, {PRIME_HIGH:,})")
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")
    print(f"search order     = {SEARCH_ORDER}")
    print()

    # ------------------------------------------------------------------
    # MODULUS INVENTORY
    # ------------------------------------------------------------------

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r, m in zip(R_VALUES, MODULI):
        print(
            f"r={r:3d} "
            f"m={m:8d}"
        )

    print()

    # ------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------

    print("=" * 78)
    print("2. PRIME POPULATION")
    print("=" * 78)

    gen_start = time.perf_counter()

    primes = generate_primes(
        PRIME_LOW,
        PRIME_HIGH
    )

    gen_time = time.perf_counter() - gen_start

    prime_set = set(primes)

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {gen_time:.3f}s")
    print(
        f"unordered prime pairs = "
        f"{len(primes) * (len(primes) - 1) // 2:,}"
    )
    print()

    # ------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------

    print("=" * 78)
    print("3. TARGETS")
    print("=" * 78)

    targets = generate_targets(
        primes,
        TARGETS
    )

    for i, (p, q, n) in enumerate(targets, start=1):
        print(
            f"target {i:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    print()

    # ------------------------------------------------------------------
    # INDEXES
    # ------------------------------------------------------------------

    print("=" * 78)
    print("4. BUILDING PRIME RESIDUE INDEXES")
    print("=" * 78)

    indexes = {}

    for prefix in PREFIXES:

        mods = PREFIX_MODULI[prefix]
        M = PREFIX_M[prefix]

        start = time.perf_counter()

        index = build_prime_index(
            primes,
            mods
        )

        indexes[prefix] = index

        elapsed = time.perf_counter() - start

        singleton_classes = sum(
            1
            for values in index.values()
            if len(values) == 1
        )

        max_class = max(
            len(values)
            for values in index.values()
        )

        print(
            f"prefix={prefix:2d} "
            f"M={M:12,d} "
            f"classes={len(index):8,d} "
            f"singleton={singleton_classes:8,d} "
            f"max_class={max_class:5d} "
            f"time={elapsed:.3f}s"
        )

    print()

    # ------------------------------------------------------------------
    # EXPERIMENT
    # ------------------------------------------------------------------

    print("=" * 78)
    print("5. N-ONLY CRT CANDIDATE RECOVERY")
    print("=" * 78)

    print()
    print(
        "The CRT condition is used to construct the candidate population"
    )
    print(
        "directly from n. Every candidate is then tested exactly against n."
    )
    print()

    all_stats = {
        prefix: []
        for prefix in PREFIXES
    }

    for target_index, (p, q, n) in enumerate(targets, start=1):

        true_pair = tuple(sorted((p, q)))

        print("-" * 78)
        print(
            f"TARGET {target_index:2d} "
            f"(control factors {true_pair[0]}, {true_pair[1]})"
        )
        print(f"n = {n}")
        print("-" * 78)

        for prefix in PREFIXES:

            mods = PREFIX_MODULI[prefix]
            M = PREFIX_M[prefix]
            index = indexes[prefix]

            start = time.perf_counter()

            candidates = build_candidate_pairs(
                n,
                primes,
                index,
                mods
            )

            construct_time = time.perf_counter() - start

            found, checked, total = exact_recovery(
                n,
                candidates
            )

            recover_time = time.perf_counter() - start

            true_present = true_pair in candidates

            # The ideal unrestricted unit-class reference.
            # Here we compare against 1/M just as Experiment 38 did.
            fraction = (
                len(candidates)
                /
                (len(primes) * (len(primes) - 1) // 2)
            )

            reference = 1.0 / M

            ratio = (
                fraction / reference
                if reference > 0
                else float("inf")
            )

            all_stats[prefix].append(
                (
                    len(candidates),
                    fraction,
                    ratio,
                    checked,
                    found == true_pair,
                    true_present,
                )
            )

            print()
            print(
                f"PREFIX {prefix:2d} "
                f"M={M:,}"
            )
            print(
                f"  candidates             = "
                f"{len(candidates):,}"
            )
            print(
                f"  candidate fraction     = "
                f"{fraction:.12g}"
            )
            print(
                f"  1/M                    = "
                f"{reference:.12g}"
            )
            print(
                f"  ratio / (1/M)          = "
                f"{ratio:.8g}"
            )
            print(
                f"  true pair present      = "
                f"{true_present}"
            )
            print(
                f"  exact recovery         = "
                f"{found == true_pair}"
            )
            print(
                f"  candidates checked     = "
                f"{checked:,}"
            )
            print(
                f"  construction time      = "
                f"{construct_time:.4f}s"
            )
            print(
                f"  total recovery time    = "
                f"{recover_time:.4f}s"
            )

            if found is not None:
                print(
                    f"  recovered pair         = "
                    f"({found[0]}, {found[1]})"
                )

            if (
                len(candidates) <= PRINT_CANDIDATES_AT_MOST
            ):
                print("  candidate list:")
                for pair in order_candidates(candidates):
                    marker = "  <-- TRUE" if pair == true_pair else ""
                    print(
                        f"    {pair}{marker}"
                    )

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. GLOBAL RECOVERY SUMMARY")
    print("=" * 78)
    print()

    print(
        f"{'prefix':>6s}"
        f"{'M':>14s}"
        f"{'median candidates':>22s}"
        f"{'median fraction':>20s}"
        f"{'median ratio':>16s}"
        f"{'recovered':>12s}"
    )
    print("-" * 78)

    for prefix in PREFIXES:

        rows = all_stats[prefix]

        candidate_values = sorted(
            row[0] for row in rows
        )

        fraction_values = sorted(
            row[1] for row in rows
        )

        ratio_values = sorted(
            row[2] for row in rows
        )

        def median(values):
            n = len(values)

            if n % 2:
                return values[n // 2]

            return (
                values[n // 2 - 1]
                + values[n // 2]
            ) / 2.0

        median_candidates = median(candidate_values)
        median_fraction = median(fraction_values)
        median_ratio = median(ratio_values)

        recovered = sum(
            1
            for row in rows
            if row[4]
        )

        print(
            f"{prefix:6d}"
            f"{PREFIX_M[prefix]:14,d}"
            f"{median_candidates:22,.1f}"
            f"{median_fraction:20.12g}"
            f"{median_ratio:16.8g}"
            f"{recovered:8d}/{TARGETS:<3d}"
        )

    # ------------------------------------------------------------------
    # RECOVERY DEPTH
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. EXACT RECOVERY DEPTH")
    print("=" * 78)
    print()

    print(
        "For each target, 'checked' is the number of CRT-compatible"
    )
    print(
        "candidate pairs actually multiplied before n was recovered."
    )
    print()

    for prefix in PREFIXES:

        rows = all_stats[prefix]

        checks = [
            row[3]
            for row in rows
        ]

        successes = sum(
            1
            for row in rows
            if row[4]
        )

        print(
            f"prefix {prefix:2d}: "
            f"median checked = {sorted(checks)[len(checks)//2]:,} "
            f"success = {successes}/{TARGETS}"
        )

    # ------------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. DIAGNOSTIC")
    print("=" * 78)
    print()

    prefix7 = all_stats[7]

    median_candidates_7 = sorted(
        row[0] for row in prefix7
    )[len(prefix7) // 2]

    print(
        f"PREFIX 7 median candidate population = "
        f"{median_candidates_7:,.0f}"
    )

    if median_candidates_7 <= 1000:
        print(
            "RESULT: the CRT filter reduces this finite population "
            "to a very small exact-search problem."
        )
    elif median_candidates_7 <= 10000:
        print(
            "RESULT: the CRT filter produces a manageable but nontrivial "
            "candidate-search problem."
        )
    else:
        print(
            "RESULT: the CRT filter still leaves a large candidate "
            "population for this interval."
        )

    print()
    print(
        "IMPORTANT:"
    )
    print(
        "This experiment proves only exact recovery within the specified"
    )
    print(
        "prime interval. It does not establish a general factoring method"
    )
    print(
        "for an unknown factor interval."
    )

    total_time = time.perf_counter() - t0

    print()
    print("=" * 78)
    print("EXPERIMENT 39 COMPLETE")
    print("=" * 78)
    print(
        f"total runtime = {total_time:.3f}s "
        f"({total_time / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

