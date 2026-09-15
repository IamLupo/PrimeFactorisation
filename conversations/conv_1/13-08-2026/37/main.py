#!/usr/bin/env python3

import math
import random
import time
from collections import defaultdict

from sympy import isprime


# ============================================================================
# KAPPA EXPERIMENT 37
# N-ONLY CRT PRIME-CANDIDATE DENSITY
# NO CSV OUTPUT
#
# Goal:
#   Measure whether the CRT constraints derived from n meaningfully reduce
#   the prime-pair search space BEFORE exact divisibility testing.
#
# Important:
#   p, q are used only as controls for verification.
#   The CRT candidate sets themselves are constructed from n only.
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

# Number of primes used for the density population.
# None = use every prime in the interval.
#
# Keeping this as None gives an exact interval experiment.
MAX_PRIME_POPULATION = None


# ============================================================================
# BASIC HELPERS
# ============================================================================

def gcd(a, b):
    return math.gcd(a, b)


def lcm(a, b):
    return abs(a // gcd(a, b) * b)


def combined_modulus(moduli):
    M = 1
    for m in moduli:
        M = lcm(M, m)
    return M


def choose_moduli(prefix):
    return [r * r + 3 for r in R_VALUES[:prefix]]


def fmt_int(x):
    return f"{x:,}"


def fmt_ratio(num, den):
    if den == 0:
        return "0"
    return f"{num / den:.8g}"


# ============================================================================
# PRIME GENERATION
# ============================================================================

def generate_primes():
    print("-" * 78)
    print("GENERATING PRIME POPULATION")
    print("-" * 78)

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

            if (
                MAX_PRIME_POPULATION is not None
                and len(primes) >= MAX_PRIME_POPULATION
            ):
                break

        x += 2

    elapsed = time.perf_counter() - start

    print(f"prime population = {fmt_int(len(primes))}")
    print(f"generation time  = {elapsed:.3f}s")
    print()

    return primes


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(primes):
    print("-" * 78)
    print("GENERATING TARGETS")
    print("-" * 78)

    targets = []

    while len(targets) < TARGETS:
        p = random.choice(primes)
        q = random.choice(primes)

        if p == q:
            continue

        n = p * q

        # Keep targets reasonably balanced.
        ratio = max(p, q) / min(p, q)

        if ratio > 2.2:
            continue

        pair = tuple(sorted((p, q)))

        if any(t[0:2] == pair for t in targets):
            continue

        targets.append((p, q, n))

    for i, (p, q, n) in enumerate(targets, 1):
        print(
            f"target {i:2d}: "
            f"p={p} q={q} n={n}"
        )

    print()

    return targets


# ============================================================================
# PRIME RESIDUE INDEX
# ============================================================================

def build_residue_index(primes, M):
    """
    Index primes by residue modulo M.

    Only primes coprime to M are included because the N-only inversion

        q = n * p^{-1} mod M

    is only valid for units.
    """

    index = defaultdict(list)

    unit_primes = 0
    nonunit_primes = 0

    for p in primes:
        if gcd(p, M) == 1:
            index[p % M].append(p)
            unit_primes += 1
        else:
            nonunit_primes += 1

    return index, unit_primes, nonunit_primes


# ============================================================================
# N-ONLY CRT CANDIDATE ANALYSIS
# ============================================================================

def analyze_target(
    p_true,
    q_true,
    n,
    primes,
    prefix,
    residue_index,
    M,
):
    """
    Entire candidate analysis is n-only.

    The true factors are used only afterwards for control.

    For every prime p with gcd(p,M)=1:

        q = n * p^{-1} mod M

    is the unique compatible q residue modulo M.

    We then look up all primes in that residue class.
    """

    target_residue_p = p_true % M
    target_residue_q = q_true % M

    # Counts before exact divisibility.
    candidate_prime_pairs = 0
    candidate_ordered_pairs = 0

    # Number of distinct p values which have at least one compatible q.
    compatible_p_count = 0

    # Number of p values for which the compatible q residue exists.
    unique_q_residue_hits = 0

    # Prime-pair multiplicity distribution.
    multiplicities = defaultdict(int)

    # Actual control pair membership.
    true_pair_found = False

    # Exact factor verification controls only.
    exact_divisor_hits = 0

    # Candidate pairs are counted once with p < q.
    for p in primes:

        if gcd(p, M) != 1:
            continue

        inv_p = pow(p, -1, M)
        q_residue = (n * inv_p) % M

        q_list = residue_index.get(q_residue)

        if not q_list:
            continue

        unique_q_residue_hits += 1

        compatible_p_count += 1
        multiplicities[len(q_list)] += 1

        for q in q_list:

            if p == q:
                continue

            candidate_ordered_pairs += 1

            if p < q:
                candidate_prime_pairs += 1

                if {p, q} == {p_true, q_true}:
                    true_pair_found = True

                # Exact divisibility is deliberately performed only
                # as a final control.
                if n % p == 0 and p * q == n:
                    exact_divisor_hits += 1

    return {
        "candidate_prime_pairs": candidate_prime_pairs,
        "candidate_ordered_pairs": candidate_ordered_pairs,
        "compatible_p_count": compatible_p_count,
        "unique_q_residue_hits": unique_q_residue_hits,
        "multiplicities": multiplicities,
        "true_pair_found": true_pair_found,
        "exact_divisor_hits": exact_divisor_hits,
        "target_residue_p": target_residue_p,
        "target_residue_q": target_residue_q,
    }


# ============================================================================
# INTEGER SEARCH BASELINES
# ============================================================================

def integer_interval_size():
    return PRIME_HIGH - PRIME_LOW


def estimate_prime_pair_baseline(prime_count):
    """
    Number of unordered prime pairs in the searchable population.
    """

    if prime_count < 2:
        return 0

    return prime_count * (prime_count - 1) // 2


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 37")
    print("N-ONLY CRT PRIME-CANDIDATE DENSITY")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print()

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime interval   = [{PRIME_LOW:,}, {PRIME_HIGH:,}]")
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")
    print(f"max prime pool   = {MAX_PRIME_POPULATION}")
    print()

    # ----------------------------------------------------------------------
    # Modulus inventory
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r in R_VALUES:
        m = r * r + 3

        # Count units by direct inspection.
        unit_count = sum(
            1 for x in range(m)
            if gcd(x, m) == 1
        )

        print(
            f"r={r:3d} "
            f"m={m:7d} "
            f"units={unit_count:7d}"
        )

    print()

    # ----------------------------------------------------------------------
    # Prime population
    # ----------------------------------------------------------------------

    primes = generate_primes()

    if len(primes) < 2:
        raise RuntimeError("Prime population is too small.")

    # ----------------------------------------------------------------------
    # Targets
    # ----------------------------------------------------------------------

    targets = generate_targets(primes)

    # ----------------------------------------------------------------------
    # Baselines
    # ----------------------------------------------------------------------

    total_prime_pairs = estimate_prime_pair_baseline(len(primes))

    print("=" * 78)
    print("2. PRIME-PAIR BASELINE")
    print("=" * 78)

    print(
        f"prime population              = {fmt_int(len(primes))}"
    )
    print(
        f"unordered prime pairs         = {fmt_int(total_prime_pairs)}"
    )

    print()

    # ----------------------------------------------------------------------
    # Global summary containers
    # ----------------------------------------------------------------------

    global_results = {
        prefix: []
        for prefix in PREFIXES
    }

    # ----------------------------------------------------------------------
    # Prefix loop
    # ----------------------------------------------------------------------

    for prefix in PREFIXES:

        print("=" * 78)
        print(f"3. PREFIX {prefix}")
        print("=" * 78)

        moduli = choose_moduli(prefix)
        M = combined_modulus(moduli)

        print(f"moduli                       = {moduli}")
        print(f"combined CRT modulus         = {fmt_int(M)}")
        print(f"modulus bits                 = {M.bit_length()}")

        if M >= PRIME_HIGH - PRIME_LOW:
            print(
                "NOTE: CRT modulus is now comparable to or larger than "
                "the search interval."
            )

        print()

        # --------------------------------------------------------------
        # Build one index for this M.
        # --------------------------------------------------------------

        index_start = time.perf_counter()

        residue_index, unit_prime_count, nonunit_prime_count = (
            build_residue_index(primes, M)
        )

        index_elapsed = time.perf_counter() - index_start

        print("PRIME RESIDUE INDEX")
        print("-" * 78)
        print(
            f"unit primes                     = "
            f"{fmt_int(unit_prime_count)}"
        )
        print(
            f"non-unit primes                 = "
            f"{fmt_int(nonunit_prime_count)}"
        )
        print(
            f"distinct residue classes        = "
            f"{fmt_int(len(residue_index))}"
        )
        print(
            f"index time                      = "
            f"{index_elapsed:.3f}s"
        )
        print()

        # --------------------------------------------------------------
        # Target analyses
        # --------------------------------------------------------------

        for target_num, (p, q, n) in enumerate(targets, 1):

            print(
                f"TARGET {target_num:2d}: "
                f"n={n}"
            )

            print(
                f"true factors CONTROL ONLY: "
                f"{p}, {q}"
            )

            # Check whether target factors are units.
            target_unit_p = gcd(p, M) == 1
            target_unit_q = gcd(q, M) == 1

            if not (target_unit_p and target_unit_q):
                print(
                    "WARNING: target factor is non-unit modulo M; "
                    "unit-only inversion does not cover this target."
                )

                print()
                continue

            start = time.perf_counter()

            result = analyze_target(
                p_true=p,
                q_true=q,
                n=n,
                primes=primes,
                prefix=prefix,
                residue_index=residue_index,
                M=M,
            )

            elapsed = time.perf_counter() - start

            candidates = result["candidate_prime_pairs"]

            # Density relative to all unordered prime pairs.
            pair_fraction = (
                candidates / total_prime_pairs
                if total_prime_pairs
                else 0.0
            )

            # Density relative to unit-prime pairs.
            unit_prime_pairs = (
                unit_prime_count * (unit_prime_count - 1) // 2
            )

            unit_fraction = (
                candidates / unit_prime_pairs
                if unit_prime_pairs
                else 0.0
            )

            # How many candidates per possible p?
            avg_candidates_per_p = (
                result["candidate_ordered_pairs"]
                / unit_prime_count
                if unit_prime_count
                else 0.0
            )

            print()
            print("N-ONLY CRT PRIME FILTER")
            print("-" * 78)

            print(
                f"candidate prime pairs          = "
                f"{fmt_int(candidates)}"
            )

            print(
                f"candidate ordered pairs        = "
                f"{fmt_int(result['candidate_ordered_pairs'])}"
            )

            print(
                f"compatible prime p values      = "
                f"{fmt_int(result['compatible_p_count'])}"
            )

            print(
                f"avg q candidates per p         = "
                f"{avg_candidates_per_p:.6f}"
            )

            print(
                f"candidate / all-prime-pairs    = "
                f"{pair_fraction:.10g}"
            )

            print(
                f"candidate / unit-prime-pairs   = "
                f"{unit_fraction:.10g}"
            )

            print(
                f"true pair in CRT candidates     = "
                f"{result['true_pair_found']}"
            )

            print(
                f"exact divisors among candidates = "
                f"{result['exact_divisor_hits']}"
            )

            print(
                f"analysis time                  = "
                f"{elapsed:.3f}s"
            )

            print()

            # ----------------------------------------------------------
            # Multiplicity summary
            # ----------------------------------------------------------

            mult = result["multiplicities"]

            if mult:
                max_mult = max(mult)
                weighted_pairs = sum(
                    k * count
                    for k, count in mult.items()
                )

                print("Q-RESIDUE MULTIPLICITY")
                print("-" * 78)
                print(
                    f"distinct multiplicities          = "
                    f"{len(mult)}"
                )
                print(
                    f"maximum q-list multiplicity       = "
                    f"{max_mult}"
                )
                print(
                    f"weighted q candidates             = "
                    f"{fmt_int(weighted_pairs)}"
                )

                # Show only a compact histogram.
                top = sorted(
                    mult.items(),
                    key=lambda kv: kv[0]
                )[:12]

                print("first multiplicities:")
                for k, count in top:
                    print(
                        f"  {k:5d} q-values : "
                        f"{fmt_int(count):>10s} p-values"
                    )

            else:
                print(
                    "No compatible prime residue classes found."
                )

            print()

            global_results[prefix].append({
                "candidate_pairs": candidates,
                "pair_fraction": pair_fraction,
                "unit_fraction": unit_fraction,
                "true_pair": result["true_pair_found"],
                "exact_hits": result["exact_divisor_hits"],
                "unit_primes": unit_prime_count,
                "nonunit_primes": nonunit_prime_count,
                "M": M,
            })

    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================

    print("=" * 78)
    print("4. GLOBAL SUMMARY")
    print("=" * 78)

    print()
    print(
        f"{'prefix':>6s} "
        f"{'M':>12s} "
        f"{'median candidates':>20s} "
        f"{'median fraction':>18s} "
        f"{'true pair found':>17s}"
    )

    for prefix in PREFIXES:

        rows = global_results[prefix]

        if not rows:
            continue

        candidates = sorted(
            row["candidate_pairs"]
            for row in rows
        )

        fractions = sorted(
            row["pair_fraction"]
            for row in rows
        )

        middle = len(rows) // 2

        if len(rows) % 2:
            median_candidates = candidates[middle]
            median_fraction = fractions[middle]
        else:
            median_candidates = (
                candidates[middle - 1] + candidates[middle]
            ) / 2

            median_fraction = (
                fractions[middle - 1] + fractions[middle]
            ) / 2

        true_found = sum(
            1 for row in rows
            if row["true_pair"]
        )

        M = rows[0]["M"]

        print(
            f"{prefix:6d} "
            f"{M:12,d} "
            f"{median_candidates:20,.1f} "
            f"{median_fraction:18.10g} "
            f"{true_found:12d}/{len(rows)}"
        )

    print()

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print("=" * 78)
    print("INTERPRETATION")
    print("=" * 78)

    print("""
This experiment isolates the question left open by the previous CRT tests.

For a target n and a combined modulus M:

    p*q = n (mod M)

and, when gcd(p,M)=1,

    q = n * p^(-1) (mod M).

Therefore the experiment can construct the compatible prime-pair
population using n alone.

No exact divisibility test is needed to construct the population.

The measurements are:

    1. number of compatible prime pairs
    2. fraction of all prime pairs surviving
    3. average number of q candidates per compatible p
    4. whether the true pair is present
    5. only afterwards, as a control, whether exact division succeeds

The important comparison is:

    CRT candidate fraction
        versus
    1 / M

If the CRT behaved like an almost-random product constraint,
the surviving pair density should decrease roughly with the
modulus scale.

But that is not sufficient for a factoring breakthrough.

The actual breakthrough criterion is much stronger:

    CRT constraints
        ->
    very small prime candidate set
        ->
    small enough to test directly

If the candidate population remains large even when M grows,
then the CRT structure is functioning primarily as a modular
filter rather than as a practical factorization mechanism.

IMPORTANT:

This experiment only treats the unit case gcd(p,M)=1.
Non-unit primes are reported separately rather than silently
discarded.

No CSV files are produced.
""")

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

