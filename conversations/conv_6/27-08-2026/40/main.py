#!/usr/bin/env python3

"""
====================================================================================================
MULTI-MODULUS PRIME-COMBINATION FINGERPRINT EXPERIMENT
====================================================================================================

Original idea:

    n mod r

where r is chosen from a set of primes.

For a factor pair:

    n = p*q

we examine the fingerprint

    (p*q mod r1, p*q mod r2, ...)

for different combinations of prime moduli.

This experiment tests all 255 non-empty subsets of:

    [3, 5, 7, 11, 13, 17, 19, 23]

The important distinction is:

    A subset of moduli creates a smaller fingerprint.

    The complete set creates:

        n mod M

    where

        M = 3*5*7*11*13*17*19*23
          = 111,546,435

For each actual semiprime anchor n=p*q:

    C = p*q mod M

we construct its EXACT C-class without building a giant product table.

For each prime a:

    a*b = C (mod M)

gives

    b = C * inverse(a,M) (mod M)

Because M > 100,000, there is at most one candidate b in the allowed
factor range.

This makes the experiment memory-efficient.

The experiment measures:

    1. Number of candidates surviving each modulus combination.
    2. Fraction of false candidates eliminated.
    3. Number of combinations required to isolate the actual pair.
    4. Whether the actual anchor behaves differently from a random anchor
       after conditioning on the exact C-class.
    5. Which modulus combinations provide the strongest discrimination.

IMPORTANT:

A candidate in the exact C-class already satisfies

    a*b == p*q (mod M)

and therefore also satisfies

    a*b == p*q (mod r)

for every prime r dividing M.

Therefore:

    The purpose here is NOT to distinguish members of the exact C-class
    using those same divisors.

Instead, we compare fingerprints at progressively smaller moduli BEFORE
conditioning on the complete M residue.

We therefore perform two related analyses:

    A. GLOBAL FINGERPRINT ANALYSIS
       Random prime pairs vs actual anchor pairs.

    B. C-CLASS STRUCTURE
       Exact candidates sharing the complete M fingerprint.

The key question is how the fingerprint grows as more moduli are combined.

====================================================================================================
"""

import math
import random
import statistics
from itertools import combinations
from collections import Counter, defaultdict


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
RANDOM_CONTROLS = 300

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

T_MIN = -25
T_MAX = 25

SEED = 1_511_464_998

# Number of random prime pairs used for the global fingerprint baseline.
GLOBAL_RANDOM_PAIRS = 50_000

# Number of examples to print.
EXAMPLES = 15


# ================================================================================================
# PRIME GENERATION
# ================================================================================================

def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)

    is_prime[0] = 0
    is_prime[1] = 0

    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            is_prime[start:limit + 1:p] = b"\x00" * count

    return [x for x in range(limit + 1) if is_prime[x]]


ALL_PRIMES = sieve(FACTOR_MAX)

FACTOR_PRIMES = [
    p for p in ALL_PRIMES
    if FACTOR_MIN <= p <= FACTOR_MAX
]

PRIME_SET = set(FACTOR_PRIMES)


# ================================================================================================
# MODULUS COMBINATIONS
# ================================================================================================

COMBINATIONS = []

for size in range(1, len(MODULI) + 1):
    for combo in combinations(MODULI, size):
        modulus = math.prod(combo)

        COMBINATIONS.append({
            "size": size,
            "moduli": combo,
            "modulus": modulus,
        })


# ================================================================================================
# BASIC HELPERS
# ================================================================================================

def random_prime_pair(rng):
    while True:
        a = rng.choice(FACTOR_PRIMES)
        b = rng.choice(FACTOR_PRIMES)

        if a != b:
            return tuple(sorted((a, b)))


def generate_unique_semiprimes(count, rng):
    seen = set()
    result = []

    while len(result) < count:

        p, q = random_prime_pair(rng)
        n = p * q

        if n not in seen:
            seen.add(n)
            result.append((p, q))

    return result


# ================================================================================================
# EXACT C-CLASS WITHOUT MASSIVE MEMORY USAGE
# ================================================================================================

def exact_c_class(p, q):
    """
    Return all prime pairs (a,b) satisfying

        a*b == p*q (mod M)

    without constructing all ~35 million prime products.

    For each a:

        b = C * inverse(a,M) mod M

    Since M > FACTOR_MAX, b has at most one representative in the
    allowed factor interval.
    """

    C = (p * q) % M

    result = set()

    for a in FACTOR_PRIMES:

        # Every prime factor in our range is coprime to M because
        # the factors of M are all below 23.
        inv_a = pow(a, -1, M)

        b = (C * inv_a) % M

        if FACTOR_MIN <= b <= FACTOR_MAX:
            if b in PRIME_SET and a != b:

                pair = tuple(sorted((a, b)))
                result.add(pair)

    return sorted(result)


# ================================================================================================
# FINGERPRINT
# ================================================================================================

def fingerprint(product, moduli):
    """
    Return the tuple

        (product mod r1, product mod r2, ...)

    """

    return tuple(product % r for r in moduli)


def fingerprint_pair(a, b, moduli):
    return fingerprint(a * b, moduli)


# ================================================================================================
# SHIFTED EVENTS
# ================================================================================================

def shifted_events(p, q):
    """
    Exact shifted products:

        a*b = p*q + 2*M*t

    for the requested t range.

    We do not search all pairs.

    Instead, for each t, factor the exact target using a bounded
    trial through the available prime factors.

    This remains small because only 51 targets are tested.
    """

    base = p * q
    events = []

    for t in range(T_MIN, T_MAX + 1):

        if t == 0:
            continue

        target = base + 2 * M * t

        if target <= 0:
            continue

        # Search divisors only among allowed primes.
        limit = math.isqrt(target)

        for a in FACTOR_PRIMES:

            if a > limit:
                break

            if target % a == 0:

                b = target // a

                if (
                    b in PRIME_SET
                    and a != b
                    and FACTOR_MIN <= b <= FACTOR_MAX
                ):

                    pair = tuple(sorted((a, b)))

                    events.append({
                        "t": t,
                        "pair": pair,
                        "product": target,
                    })

                    break

    return events


# ================================================================================================
# GLOBAL RANDOM PAIR SAMPLE
# ================================================================================================

def generate_random_pairs(count, rng):
    return [
        random_prime_pair(rng)
        for _ in range(count)
    ]


# ================================================================================================
# GLOBAL FINGERPRINT ANALYSIS
# ================================================================================================

def analyze_global_fingerprints(actual_anchors, random_pairs):

    """
    For each subset of the eight moduli:

        - count unique fingerprints
        - measure collision rate
        - measure average fingerprint class size

    We do this for:

        actual anchor products

    and

        random prime-pair products
    """

    actual_products = [
        p * q for p, q in actual_anchors
    ]

    random_products = [
        p * q for p, q in random_pairs
    ]

    rows = []

    for combo in COMBINATIONS:

        moduli = combo["moduli"]

        actual_counts = Counter(
            fingerprint(x, moduli)
            for x in actual_products
        )

        random_counts = Counter(
            fingerprint(x, moduli)
            for x in random_products
        )

        actual_unique = len(actual_counts)
        random_unique = len(random_counts)

        actual_max = max(actual_counts.values())
        random_max = max(random_counts.values())

        actual_entropy = entropy_from_counter(actual_counts)
        random_entropy = entropy_from_counter(random_counts)

        rows.append({
            "size": combo["size"],
            "moduli": moduli,
            "modulus": combo["modulus"],

            "actual_unique": actual_unique,
            "random_unique": random_unique,

            "actual_entropy": actual_entropy,
            "random_entropy": random_entropy,

            "actual_max_collision": actual_max,
            "random_max_collision": random_max,

            "actual_collision_pairs":
                sum(v * (v - 1) // 2 for v in actual_counts.values()),

            "random_collision_pairs":
                sum(v * (v - 1) // 2 for v in random_counts.values()),
        })

    return rows


# ================================================================================================
# ENTROPY
# ================================================================================================

def entropy_from_counter(counter):
    total = sum(counter.values())

    if total == 0:
        return 0.0

    h = 0.0

    for count in counter.values():

        p = count / total

        h -= p * math.log2(p)

    return h


# ================================================================================================
# ANCHOR CLASS ANALYSIS
# ================================================================================================

def analyze_anchor(p, q):

    base = p * q

    candidates = exact_c_class(p, q)

    false_candidates = [
        pair
        for pair in candidates
        if pair != tuple(sorted((p, q)))
    ]

    return {
        "anchor": tuple(sorted((p, q))),
        "product": base,
        "C": base % M,
        "class_size": len(candidates),
        "false_candidates": false_candidates,
        "events": shifted_events(p, q),
    }


# ================================================================================================
# C-CLASS AGGREGATION
# ================================================================================================

def aggregate_classes(results):

    class_sizes = []
    false_counts = []

    total_events = 0

    for result in results:

        class_sizes.append(result["class_size"])
        false_counts.append(
            len(result["false_candidates"])
        )

        total_events += len(result["events"])

    return {
        "class_sizes": class_sizes,
        "false_counts": false_counts,
        "total_events": total_events,
    }


# ================================================================================================
# COMBINATION SURVIVAL
# ================================================================================================

def candidate_survival_for_combination(
    anchor_product,
    candidates,
    moduli
):
    """
    Number of candidates sharing the anchor's fingerprint under
    the specified modulus combination.
    """

    anchor_fp = fingerprint(
        anchor_product,
        moduli
    )

    count = 0

    for a, b in candidates:

        if fingerprint_pair(a, b, moduli) == anchor_fp:
            count += 1

    return count


def analyze_combination_survival(results):

    rows = []

    for combo in COMBINATIONS:

        moduli = combo["moduli"]

        survivors = []

        for result in results:

            count = candidate_survival_for_combination(
                result["product"],
                result["false_candidates"],
                moduli
            )

            survivors.append(count)

        rows.append({
            "size": combo["size"],
            "moduli": moduli,
            "modulus": combo["modulus"],

            "mean_false_survivors":
                statistics.mean(survivors)
                if survivors else 0.0,

            "max_false_survivors":
                max(survivors)
                if survivors else 0,

            "zero_survivor_anchors":
                sum(x == 0 for x in survivors),

            "anchors_with_survivors":
                sum(x > 0 for x in survivors),
        })

    return rows


# ================================================================================================
# SHIFTED EVENT MODULAR SIGNATURES
# ================================================================================================

def analyze_shifted_events(results):

    """
    For every shifted collision event, calculate its fingerprint under
    every modulus combination.

    Because the shifted event differs from the anchor by

        2*M*t

    its residue modulo every divisor of M is identical.

    This explicitly demonstrates what information the prime-modulus
    fingerprint can and cannot distinguish.
    """

    total_events = 0
    all_identity_failures = 0

    per_modulus = {
        r: {
            "events": 0,
            "failures": 0,
        }
        for r in MODULI
    }

    for result in results:

        base = result["product"]

        for event in result["events"]:

            total_events += 1

            a, b = event["pair"]

            shifted = a * b

            if shifted != base + 2 * M * event["t"]:
                all_identity_failures += 1

            for r in MODULI:

                if shifted % r != base % r:
                    per_modulus[r]["failures"] += 1

                per_modulus[r]["events"] += 1

    return {
        "events": total_events,
        "identity_failures": all_identity_failures,
        "per_modulus": per_modulus,
    }


# ================================================================================================
# PRINTING
# ================================================================================================

def print_header():

    print("=" * 100)
    print("MULTI-MODULUS PRIME-COMBINATION FINGERPRINT EXPERIMENT")
    print("=" * 100)

    print(f"M                    = {M:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"T range              = [{T_MIN}, {T_MAX}]")
    print(f"actual anchors       = {TRIALS}")
    print(f"random controls      = {RANDOM_CONTROLS}")
    print(f"global random pairs  = {GLOBAL_RANDOM_PAIRS:,}")
    print(f"prime moduli         = {MODULI}")
    print(f"number of combinations = {len(COMBINATIONS)}")
    print(f"factor primes        = {len(FACTOR_PRIMES):,}")
    print(f"seed                 = {SEED:,}")
    print()


def print_global_by_size(rows):

    print("=" * 100)
    print("GLOBAL FINGERPRINT INFORMATION BY NUMBER OF MODULI")
    print("=" * 100)

    print(
        f"{'K':>3} "
        f"{'ACT UNIQUE':>12} "
        f"{'RND UNIQUE':>12} "
        f"{'ACT H':>12} "
        f"{'RND H':>12} "
        f"{'BEST ACT COLL':>15} "
        f"{'BEST RND COLL':>15}"
    )

    print("-" * 100)

    for k in range(1, len(MODULI) + 1):

        subset = [
            x for x in rows
            if x["size"] == k
        ]

        actual_unique = statistics.mean(
            x["actual_unique"] for x in subset
        )

        random_unique = statistics.mean(
            x["random_unique"] for x in subset
        )

        actual_h = statistics.mean(
            x["actual_entropy"] for x in subset
        )

        random_h = statistics.mean(
            x["random_entropy"] for x in subset
        )

        best_actual = max(
            x["actual_max_collision"] for x in subset
        )

        best_random = max(
            x["random_max_collision"] for x in subset
        )

        print(
            f"{k:3d} "
            f"{actual_unique:12.2f} "
            f"{random_unique:12.2f} "
            f"{actual_h:12.5f} "
            f"{random_h:12.5f} "
            f"{best_actual:15d} "
            f"{best_random:15d}"
        )

    print()


def print_best_combinations(rows):

    print("=" * 100)
    print("MOST COLLISION-RESISTANT MODULUS COMBINATIONS")
    print("=" * 100)

    for k in range(1, len(MODULI) + 1):

        subset = [
            x for x in rows
            if x["size"] == k
        ]

        subset.sort(
            key=lambda x: (
                -x["actual_unique"],
                -x["actual_entropy"]
            )
        )

        best = subset[0]

        print(
            f"K={k} "
            f"moduli={best['moduli']} "
            f"product={best['modulus']:,} "
            f"actual_unique={best['actual_unique']} "
            f"actual_entropy={best['actual_entropy']:.6f}"
        )

    print()


def print_c_class_summary(results):

    agg = aggregate_classes(results)

    print("=" * 100)
    print("EXACT C-CLASS STRUCTURE")
    print("=" * 100)

    print(
        f"mean class size      = "
        f"{statistics.mean(agg['class_sizes']):.6f}"
    )

    print(
        f"minimum class size   = "
        f"{min(agg['class_sizes'])}"
    )

    print(
        f"maximum class size   = "
        f"{max(agg['class_sizes'])}"
    )

    print(
        f"mean false candidates = "
        f"{statistics.mean(agg['false_counts']):.6f}"
    )

    print(
        f"anchors with false candidates = "
        f"{sum(x > 0 for x in agg['false_counts'])}"
    )

    print()


def print_survival_summary(rows):

    print("=" * 100)
    print("EXACT C-CLASS SURVIVAL")
    print("=" * 100)

    print(
        f"{'K':>3} "
        f"{'COMBINATIONS':>12} "
        f"{'MEAN FALSE SURVIVORS':>22} "
        f"{'MAX FALSE':>12} "
        f"{'ZERO':>8}"
    )

    print("-" * 100)

    for k in range(1, len(MODULI) + 1):

        subset = [
            x for x in rows
            if x["size"] == k
        ]

        mean_survivors = statistics.mean(
            x["mean_false_survivors"]
            for x in subset
        )

        max_survivors = max(
            x["max_false_survivors"]
            for x in subset
        )

        zero = sum(
            x["zero_survivor_anchors"]
            for x in subset
        )

        print(
            f"{k:3d} "
            f"{len(subset):12d} "
            f"{mean_survivors:22.6f} "
            f"{max_survivors:12d} "
            f"{zero:8d}"
        )

    print()


def print_shifted_analysis(shift):

    print("=" * 100)
    print("SHIFTED COLLISION MODULAR TEST")
    print("=" * 100)

    print(
        f"shifted collision events = "
        f"{shift['events']}"
    )

    print(
        f"identity failures       = "
        f"{shift['identity_failures']}"
    )

    print()

    print(
        f"{'r':>5} "
        f"{'EVENTS':>12} "
        f"{'RESIDUE FAILURES':>18}"
    )

    print("-" * 100)

    for r in MODULI:

        row = shift["per_modulus"][r]

        print(
            f"{r:5d} "
            f"{row['events']:12d} "
            f"{row['failures']:18d}"
        )

    print()


def print_examples(results):

    print("=" * 100)
    print("EXACT C-CLASS EXAMPLES")
    print("=" * 100)

    shown = 0

    for result in results:

        if shown >= EXAMPLES:
            break

        if not result["false_candidates"]:
            continue

        p, q = result["anchor"]

        print()
        print(
            f"anchor=({p:,},{q:,}) "
            f"C={result['C']:,} "
            f"class={result['class_size']}"
        )

        for a, b in result["false_candidates"]:

            print(
                f"    false=({a:,},{b:,}) "
                f"product={a*b:,}"
            )

        shown += 1

    print()


# ================================================================================================
# MAIN
# ================================================================================================

def main():

    print_header()

    rng = random.Random(SEED)

    # --------------------------------------------------------------------------------------------
    # Actual anchors
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("GENERATING ACTUAL ANCHORS")
    print("=" * 100)

    actual_anchors = generate_unique_semiprimes(
        TRIALS,
        rng
    )

    print(
        f"anchors = {len(actual_anchors)}"
    )

    print()

    # --------------------------------------------------------------------------------------------
    # Actual exact C-classes
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING EXACT C-CLASSES")
    print("=" * 100)

    actual_results = []

    for i, (p, q) in enumerate(
        actual_anchors,
        start=1
    ):

        result = analyze_anchor(p, q)

        actual_results.append(result)

        if i % 25 == 0 or i == TRIALS:

            print(
                f"anchor {i:3d}/{TRIALS}"
            )

    print()

    # --------------------------------------------------------------------------------------------
    # Random global controls
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("GENERATING GLOBAL RANDOM PRIME-PAIR CONTROL")
    print("=" * 100)

    random_pairs = generate_random_pairs(
        GLOBAL_RANDOM_PAIRS,
        rng
    )

    print(
        f"random pairs = {len(random_pairs):,}"
    )

    print()

    # --------------------------------------------------------------------------------------------
    # Global fingerprint experiment
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("GLOBAL FINGERPRINT ANALYSIS")
    print("=" * 100)

    global_rows = analyze_global_fingerprints(
        actual_anchors,
        random_pairs
    )

    print_global_by_size(
        global_rows
    )

    print_best_combinations(
        global_rows
    )

    # --------------------------------------------------------------------------------------------
    # C-class analysis
    # --------------------------------------------------------------------------------------------

    print_c_class_summary(
        actual_results
    )

    survival_rows = analyze_combination_survival(
        actual_results
    )

    print_survival_summary(
        survival_rows
    )

    # --------------------------------------------------------------------------------------------
    # Shifted events
    # --------------------------------------------------------------------------------------------

    shift = analyze_shifted_events(
        actual_results
    )

    print_shifted_analysis(
        shift
    )

    # --------------------------------------------------------------------------------------------
    # Sanity examples
    # --------------------------------------------------------------------------------------------

    print_examples(
        actual_results
    )

    # --------------------------------------------------------------------------------------------
    # Final interpretation
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
The experiment separates two different questions.

QUESTION 1
----------

How much information does a collection of prime moduli provide?

For a set S:

    fingerprint_S(n)
        =
    (n mod r) for r in S

Adding moduli makes the fingerprint more specific.

However, once all prime factors of M are included,

    M = 3*5*7*11*13*17*19*23,

the complete fingerprint is equivalent to n mod M by CRT.

Therefore a larger collection does NOT create information beyond n mod M.

It creates a finer representation of the same residue information.

QUESTION 2
----------

Does the actual semiprime population have a different fingerprint collision
structure from generic random prime products?

That is what the GLOBAL section tests.

The C-CLASS section asks a different question:

    Once C = p*q mod M is fixed,

    are there multiple prime pairs on that exact modular hyperbola?

Because every member of the same C-class has the same residue modulo every
prime divisor of M, no combination of those same divisors can distinguish
two members of that exact class.

This is mathematically necessary.

Therefore:

    more divisors of M
        -> stronger fingerprint
        -> fewer collisions among unrestricted products

but:

    conditioning on the complete M residue
        -> all remaining candidates have identical fingerprints
           under every divisor of M.

For the shifted experiment:

    a*b = p*q + 2*M*t

implies

    a*b == p*q (mod r)

for every r dividing M.

Consequently every one of the eight prime-modulus tests must agree for
every shifted collision.

That is not experimental evidence; it follows algebraically.

The useful next step is therefore NOT simply adding more divisors of M.

The useful next step is to test additional moduli that are NOT divisors
of M, or to study the factor residues themselves rather than only the
product residue.

"""
    )

    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()