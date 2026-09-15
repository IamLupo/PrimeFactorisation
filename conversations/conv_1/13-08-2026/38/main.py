#!/usr/bin/env python3

import math
import random
import time
from collections import Counter

from sympy import isprime


# ============================================================================
# KAPPA EXPERIMENT 38
# EXACT CRT RESIDUE-CLASS CANDIDATE COUNT
# NO CSV OUTPUT
#
# Purpose:
#
# Experiment 37 tried to enumerate compatible prime candidates directly.
# That is unnecessary.
#
# If
#
#     p*q = n (mod M)
#
# and p is a unit modulo M, then
#
#     q = n * p^(-1) (mod M).
#
# We therefore work entirely at the residue-class level.
#
# For every possible residue a:
#
#     b = n * a^(-1) mod M
#
# We count how many primes occupy residue a and residue b.
#
# This gives the EXACT number of prime-pair candidates surviving the
# N-only CRT constraint, without enumerating the prime-pair combinations.
#
# Only the prime population is enumerated.
# No candidate-pair enumeration is performed.
# ============================================================================


SEED = 20260814
random.seed(SEED)

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

PREFIXES = [3, 4, 5, 6, 7]

TARGETS = 12

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000


# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

def modulus_from_r(r):
    return r * r + 3


def gcd(a, b):
    return math.gcd(a, b)


def lcm(a, b):
    return abs(a // gcd(a, b) * b)


def combined_modulus(moduli):
    M = 1

    for m in moduli:
        M = lcm(M, m)

    return M


def fmt(x):
    return f"{x:,}"


def ratio(a, b):
    if b == 0:
        return 0.0

    return a / b


# ============================================================================
# PRIME GENERATION
# ============================================================================

def generate_prime_population():
    start = time.perf_counter()

    primes = []

    x = PRIME_LOW

    if x <= 2:
        x = 2
    elif x % 2 == 0:
        x += 1

    while x < PRIME_HIGH:

        if isprime(x):
            primes.append(x)

        x += 2

    elapsed = time.perf_counter() - start

    return primes, elapsed


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(primes):
    targets = []

    while len(targets) < TARGETS:

        p = random.choice(primes)
        q = random.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        # Keep examples reasonably balanced.
        if q > 2.2 * p:
            continue

        pair = (p, q)

        if pair in targets:
            continue

        targets.append(pair)

    return targets


# ============================================================================
# PRIME RESIDUE COUNTS
# ============================================================================

def build_residue_counts(primes, M):

    counts = Counter()

    unit_count = 0
    nonunit_count = 0

    for p in primes:

        if gcd(p, M) == 1:
            counts[p % M] += 1
            unit_count += 1

        else:
            nonunit_count += 1

    return counts, unit_count, nonunit_count


# ============================================================================
# N-ONLY CRT RESIDUE PAIR CLASSES
# ============================================================================

def compatible_residue_pairs(n, M, residue_counts):
    """
    Return all residue-pair classes

        (a,b)

    satisfying

        a*b = n mod M

    with both a,b units.

    The construction uses n only.
    """

    pairs = []

    n_mod = n % M

    for a in residue_counts:

        if gcd(a, M) != 1:
            continue

        inv_a = pow(a, -1, M)
        b = (n_mod * inv_a) % M

        if b not in residue_counts:
            continue

        pairs.append((a, b))

    return pairs


# ============================================================================
# EXACT PRIME-PAIR COUNT FROM RESIDUE CLASSES
# ============================================================================

def exact_unordered_candidate_count(
    residue_pairs,
    residue_counts,
):
    """
    Convert compatible residue classes into an exact number of
    unordered prime pairs.

    For a != b:

        count[a] * count[b]

    For a == b:

        count[a] choose 2

    Each unordered residue-class combination is counted once.
    """

    total = 0

    seen = set()

    for a, b in residue_pairs:

        key = tuple(sorted((a, b)))

        if key in seen:
            continue

        seen.add(key)

        ca = residue_counts.get(a, 0)
        cb = residue_counts.get(b, 0)

        if ca == 0 or cb == 0:
            continue

        if a == b:
            total += ca * (ca - 1) // 2

        else:
            total += ca * cb

    return total


# ============================================================================
# EXACT ORDERED COUNT
# ============================================================================

def exact_ordered_candidate_count(
    residue_pairs,
    residue_counts,
):
    total = 0

    for a, b in residue_pairs:

        ca = residue_counts.get(a, 0)
        cb = residue_counts.get(b, 0)

        total += ca * cb

        if a == b:
            total -= ca

    return total


# ============================================================================
# RESIDUE-CLASS BREAKDOWN
# ============================================================================

def class_breakdown(
    residue_pairs,
    residue_counts,
    n,
    M,
):
    """
    Produce a compact residue-class breakdown.

    Only the largest classes are printed.
    """

    rows = []
    seen = set()

    for a, b in residue_pairs:

        key = tuple(sorted((a, b)))

        if key in seen:
            continue

        seen.add(key)

        ca = residue_counts[a]
        cb = residue_counts[b]

        if a == b:
            contribution = ca * (ca - 1) // 2

        else:
            contribution = ca * cb

        rows.append(
            (
                contribution,
                a,
                b,
                ca,
                cb,
            )
        )

    rows.sort(reverse=True)

    return rows


# ============================================================================
# EXACT TARGET VERIFICATION
# ============================================================================

def true_residue_pair_present(
    p,
    q,
    residue_counts,
    M,
    compatible_pairs,
):
    a = p % M
    b = q % M

    direct = (a, b) in compatible_pairs
    reverse = (b, a) in compatible_pairs

    return direct or reverse


# ============================================================================
# PRIME PAIR BASELINE
# ============================================================================

def total_unordered_prime_pairs(prime_count):
    return prime_count * (prime_count - 1) // 2


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 38")
    print("EXACT CRT RESIDUE-CLASS CANDIDATE COUNT")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print()
    print(f"random seed    = {SEED}")
    print(f"targets        = {TARGETS}")
    print(
        f"prime interval = "
        f"[{PRIME_LOW:,}, {PRIME_HIGH:,}]"
    )
    print(f"prefixes       = {PREFIXES}")
    print(f"R values       = {R_VALUES}")

    # ---------------------------------------------------------------------
    # PRIME POPULATION
    # ---------------------------------------------------------------------

    print()
    print("-" * 78)
    print("1. PRIME POPULATION")
    print("-" * 78)

    primes, prime_time = generate_prime_population()

    print(
        f"prime population = {fmt(len(primes))}"
    )

    print(
        f"generation time  = {prime_time:.3f}s"
    )

    total_pairs = total_unordered_prime_pairs(len(primes))

    print(
        f"unordered prime pairs = {fmt(total_pairs)}"
    )

    # ---------------------------------------------------------------------
    # TARGETS
    # ---------------------------------------------------------------------

    print()
    print("-" * 78)
    print("2. TARGETS")
    print("-" * 78)

    targets = generate_targets(primes)

    for i, (p, q) in enumerate(targets, 1):

        n = p * q

        print(
            f"target {i:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    # ---------------------------------------------------------------------
    # PREFIX EXPERIMENT
    # ---------------------------------------------------------------------

    global_results = []

    for prefix in PREFIXES:

        print()
        print("=" * 78)
        print(f"3. PREFIX {prefix}")
        print("=" * 78)

        moduli = [
            modulus_from_r(r)
            for r in R_VALUES[:prefix]
        ]

        M = combined_modulus(moduli)

        print(f"moduli               = {moduli}")
        print(f"combined modulus     = {fmt(M)}")
        print(f"modulus bits         = {M.bit_length()}")

        # ---------------------------------------------------------------
        # Build residue counts exactly once for this M.
        # ---------------------------------------------------------------

        t0 = time.perf_counter()

        residue_counts, unit_count, nonunit_count = (
            build_residue_counts(
                primes,
                M
            )
        )

        index_time = time.perf_counter() - t0

        print()
        print("RESIDUE POPULATION")
        print("-" * 78)

        print(
            f"unit primes          = {fmt(unit_count)}"
        )

        print(
            f"non-unit primes      = {fmt(nonunit_count)}"
        )

        print(
            f"occupied unit classes = "
            f"{fmt(len(residue_counts))}"
        )

        print(
            f"index time            = "
            f"{index_time:.3f}s"
        )

        # ---------------------------------------------------------------
        # Target loop.
        # ---------------------------------------------------------------

        prefix_candidates = []

        for target_num, (p, q) in enumerate(
            targets,
            1
        ):

            n = p * q

            t1 = time.perf_counter()

            residue_pairs = compatible_residue_pairs(
                n=n,
                M=M,
                residue_counts=residue_counts,
            )

            unordered_candidates = (
                exact_unordered_candidate_count(
                    residue_pairs,
                    residue_counts,
                )
            )

            ordered_candidates = (
                exact_ordered_candidate_count(
                    residue_pairs,
                    residue_counts,
                )
            )

            target_present = (
                true_residue_pair_present(
                    p,
                    q,
                    residue_counts,
                    M,
                    residue_pairs,
                )
            )

            elapsed = time.perf_counter() - t1

            all_pair_fraction = ratio(
                unordered_candidates,
                total_pairs
            )

            unit_pair_total = (
                unit_count *
                (unit_count - 1) //
                2
            )

            unit_pair_fraction = ratio(
                unordered_candidates,
                unit_pair_total
            )

            # -----------------------------------------------------------
            # Ideal random-product comparison.
            #
            # If residues were distributed uniformly, p*q=n mod M
            # would leave roughly 1/M of ordered pairs.
            #
            # For unordered pairs this is only a rough comparison,
            # so we explicitly label it as a heuristic.
            # -----------------------------------------------------------

            ideal_fraction = 1 / M

            print()
            print(
                f"TARGET {target_num:2d}: "
                f"n={n}"
            )

            print("-" * 78)

            print(
                f"compatible residue-pair classes = "
                f"{fmt(len(residue_pairs))}"
            )

            print(
                f"ordered prime candidates         = "
                f"{fmt(ordered_candidates)}"
            )

            print(
                f"unordered prime candidates       = "
                f"{fmt(unordered_candidates)}"
            )

            print(
                f"candidate / all prime pairs      = "
                f"{all_pair_fraction:.12g}"
            )

            print(
                f"candidate / unit prime pairs     = "
                f"{unit_pair_fraction:.12g}"
            )

            print(
                f"ideal 1/M reference              = "
                f"{ideal_fraction:.12g}"
            )

            print(
                f"observed / (1/M)                 = "
                f"{all_pair_fraction / ideal_fraction:.6g}"
            )

            print(
                f"true factor pair residue-valid  = "
                f"{target_present}"
            )

            print(
                f"computation time                 = "
                f"{elapsed:.6f}s"
            )

            # -----------------------------------------------------------
            # Largest residue-class contributions.
            # -----------------------------------------------------------

            rows = class_breakdown(
                residue_pairs,
                residue_counts,
                n,
                M,
            )

            print()
            print("LARGEST RESIDUE-CLASS CONTRIBUTIONS")
            print("-" * 78)

            print(
                f"{'a':>8s}"
                f"{'b':>8s}"
                f"{'count(a)':>14s}"
                f"{'count(b)':>14s}"
                f"{'pairs':>18s}"
            )

            for contribution, a, b, ca, cb in rows[:8]:

                print(
                    f"{a:8d}"
                    f"{b:8d}"
                    f"{ca:14,d}"
                    f"{cb:14,d}"
                    f"{contribution:18,d}"
                )

            prefix_candidates.append(
                unordered_candidates
            )

        # ---------------------------------------------------------------
        # Prefix median.
        # ---------------------------------------------------------------

        ordered_values = sorted(
            prefix_candidates
        )

        middle = len(ordered_values) // 2

        if len(ordered_values) % 2:

            median_candidates = ordered_values[middle]

        else:

            median_candidates = (
                ordered_values[middle - 1]
                + ordered_values[middle]
            ) / 2

        median_fraction = (
            median_candidates /
            total_pairs
        )

        global_results.append(
            (
                prefix,
                M,
                median_candidates,
                median_fraction,
            )
        )

    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================

    print()
    print("=" * 78)
    print("4. GLOBAL SUMMARY")
    print("=" * 78)

    print()

    print(
        f"{'prefix':>8s}"
        f"{'M':>12s}"
        f"{'median candidates':>22s}"
        f"{'median fraction':>20s}"
        f"{'ratio to 1/M':>18s}"
    )

    print("-" * 78)

    for prefix, M, median_candidates, median_fraction in (
        global_results
    ):

        reference = 1 / M

        ratio_to_reference = (
            median_fraction /
            reference
        )

        print(
            f"{prefix:8d}"
            f"{M:12,d}"
            f"{median_candidates:22,.1f}"
            f"{median_fraction:20.12g}"
            f"{ratio_to_reference:18.6g}"
        )

    # =========================================================================
    # IMPORTANT INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("5. INTERPRETATION")
    print("=" * 78)

    print("""
This experiment removes the expensive prime-pair enumeration completely.

The only large operation is generating the prime population.

For each prefix and target, the N-only CRT condition gives:

    q = n * p^(-1) mod M

at the residue level.

Instead of looping over every possible prime p and then every compatible
prime q, we count how many primes occupy each residue class.

The resulting candidate count is therefore exact for the finite prime
population.

The key statistic is:

    candidate fraction
        =
    CRT-compatible prime pairs
    --------------------------------
    all unordered prime pairs

Compare this with:

    1 / M

If the measured fraction stays near 1/M, then the CRT relation is acting
mainly as an ordinary modular-product filter.

If it is dramatically smaller, that would indicate extra structure in
the prime distribution or in the chosen moduli.

But even a very small candidate fraction is not yet a factorization
algorithm.

The practical threshold is:

    number of CRT-compatible prime pairs
        small enough to enumerate directly.

The experiment also reports whether the true factor pair is present,
which is a control only. The candidate population itself is constructed
from n and the CRT constraint.

No CSV files are produced.
""")

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

