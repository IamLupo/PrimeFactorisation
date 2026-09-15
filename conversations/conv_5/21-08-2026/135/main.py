#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 696
==========================================================================================

ALL-DIVISOR H -> p CONGRUENCE RECOVERY

Experiment 695 established the exact identity for Frame A:

    C = q + 3

    n + c = p(C-3) + c
          = pC + (c - 3p)

Therefore, for

    H_c = gcd(C, n+c)

we always have:

    H_c | C
    H_c | (c - 3p)

and hence:

    3p ≡ c (mod H_c).

The critical flaw in Experiment 695:

    TEST 6 and TEST 7 used the TRUE hidden H_c.

That does not test recovery.

This experiment removes that information.

For each shift c:

    M_c = n+c
    factor(M_c)
    enumerate EVERY divisor H of M_c

For every divisor H we solve:

    3p ≡ c (mod H)

giving a possible p residue class.

We then apply only public structural bounds:

    p > MIN_P
    p < q
    q > p
    p <= sqrt(n)

but we DO NOT initially use:

    q = n/p

and we DO NOT select H using the true C.

The main measurements are:

    1. how many p values survive one shift?
    2. how many survive multiple shifts?
    3. does [3,81,137] actually intersect the all-divisor
       congruence channels?
    4. does the true p remain present?
    5. how often does the resulting p set become unique?
    6. after the congruence stage, how much exact factor
       validation is still required?

An exact n % p check is performed ONLY as a final validation
stage, never during candidate generation.

==========================================================================================
"""

from math import gcd, isqrt
from collections import Counter
from sympy import primerange, factorint, divisors


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 400
MIN_P = 137

SHIFTS = [3, 81, 137]

MAX_EXAMPLES = 20


# =============================================================================
# BASIC FUNCTIONS
# =============================================================================

def generate_semiprimes(limit):
    primes = list(primerange(3, limit))

    for i, p in enumerate(primes):
        for q in primes[i:]:
            yield p * q, p, q


def ordered_states(states):
    """
    We deliberately restrict to:

        q > p > MIN_P

    so that the structural ordering is part of the experiment.
    """

    return [
        (n, p, q)
        for n, p, q in states
        if q > p > MIN_P
    ]


def solve_linear_congruence_3p(c, H):
    """
    Solve:

        3p ≡ c (mod H)

    Returns:

        (residue, modulus)

    describing:

        p ≡ residue (mod modulus)

    or None if no solution exists.
    """

    g = gcd(3, H)

    if c % g != 0:
        return None

    reduced_H = H // g
    reduced_c = c // g
    reduced_a = 3 // g

    if reduced_H == 1:
        return 0, 1

    inverse = pow(reduced_a, -1, reduced_H)

    residue = (
        reduced_c * inverse
    ) % reduced_H

    return residue, reduced_H


def values_from_residue(
    residue,
    modulus,
    lower,
    upper,
):
    """
    Return all integers x in [lower, upper] satisfying

        x ≡ residue (mod modulus).
    """

    if lower > upper:
        return set()

    if modulus == 1:
        return set(range(lower, upper + 1))

    first = lower + (
        (residue - lower) % modulus
    )

    if first > upper:
        return set()

    return set(
        range(
            first,
            upper + 1,
            modulus,
        )
    )


def all_divisor_p_candidates(n, c):
    """
    Generate p candidates using ONLY:

        H | (n+c)
        3p ≡ c (mod H)
        MIN_P < p <= sqrt(n)

    No true H is used.

    No q=n/p check is used.

    Returns:

        candidates:
            union of all p values admitted by at least one divisor H

        witnesses:
            H -> residue/modulus information
    """

    M = n + c

    if M <= 0:
        return set(), {}

    p_low = MIN_P + 1
    p_high = isqrt(n)

    candidates = set()
    witnesses = {}

    for H in divisors(M):

        solved = solve_linear_congruence_3p(c, H)

        if solved is None:
            continue

        residue, modulus = solved

        values = values_from_residue(
            residue,
            modulus,
            p_low,
            p_high,
        )

        if not values:
            continue

        witnesses[H] = (
            residue,
            modulus,
            values,
        )

        candidates.update(values)

    return candidates, witnesses


# =============================================================================
# START
# =============================================================================

print("=" * 90)
print("EXPERIMENT 696 START")
print("=" * 90)
print()

print(f"prime limit={PRIME_LIMIT}")
print(f"minimum p={MIN_P}")
print(f"shifts={SHIFTS}")

states = list(generate_semiprimes(PRIME_LIMIT))

print(f"semiprimes={len(states)}")

ordered = ordered_states(states)

print(f"ordered states={len(ordered)}")
print()


# =============================================================================
# TEST 0
# =============================================================================

print("=" * 90)
print("TEST 0: BASELINE")
print("=" * 90)

failures = 0

for n, p, q in states:

    if p * q != n:
        failures += 1

print(f"checked={len(states)}")
print(f"failures={failures}")
print()


# =============================================================================
# TEST 1:
# ALL-DIVISOR p CANDIDATES FOR EACH SHIFT
# =============================================================================

print("=" * 90)
print("TEST 1: ALL DIVISORS H -> p CANDIDATES")
print("=" * 90)
print()

single_shift_sets = {}

for c in SHIFTS:

    unique = 0
    ambiguous = 0
    missing = 0

    candidate_distribution = Counter()

    examples = []

    for n, p, q in ordered:

        candidates, witnesses = all_divisor_p_candidates(
            n,
            c,
        )

        single_shift_sets[(n, p, q, c)] = candidates

        count = len(candidates)

        candidate_distribution[count] += 1

        if p not in candidates:
            missing += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append(
                    (
                        n,
                        p,
                        q,
                        sorted(candidates),
                    )
                )

        elif count == 1:
            unique += 1

        else:
            ambiguous += 1

    total = len(ordered)

    print(
        f"c={c:3d}"
    )

    print(
        f"    true p present="
        f"{total-missing}/{total}"
    )

    print(
        f"    unique={unique}"
        f" ambiguous={ambiguous}"
        f" missing={missing}"
    )

    print(
        f"    unique ratio="
        f"{unique/total:.6f}"
    )

    print(
        f"    mean candidates="
        f"{sum(k*v for k,v in candidate_distribution.items())/total:.4f}"
    )

    print(
        f"    minimum candidates="
        f"{min(candidate_distribution)}"
    )

    print(
        f"    maximum candidates="
        f"{max(candidate_distribution)}"
    )

    if examples:

        print(
            "    missing examples:"
        )

        for (
            n,
            p,
            q,
            values,
        ) in examples:

            print(
                f"        n={n} "
                f"true p={p} q={q}"
            )

            print(
                f"            candidates={values}"
            )

    print()


# =============================================================================
# TEST 2:
# MULTI-SHIFT UNION
# =============================================================================

print("=" * 90)
print("TEST 2: MULTI-SHIFT UNION")
print("=" * 90)
print()

for r in range(1, len(SHIFTS) + 1):

    active = SHIFTS[:r]

    counts = Counter()

    true_present = 0
    unique = 0
    missing = 0

    examples = []

    for n, p, q in ordered:

        union = set()

        for c in active:

            candidates = single_shift_sets[
                (n, p, q, c)
            ]

            union.update(candidates)

        count = len(union)

        counts[count] += 1

        if p not in union:

            missing += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append(
                    (
                        n,
                        p,
                        q,
                        sorted(union),
                    )
                )

        else:

            true_present += 1

            if count == 1:
                unique += 1

    print(
        f"shifts={active}"
    )

    print(
        f"    true p present="
        f"{true_present}/{len(ordered)}"
    )

    print(
        f"    unique={unique}"
        f" ambiguous={true_present-unique}"
        f" missing={missing}"
    )

    print(
        f"    mean candidates="
        f"{sum(k*v for k,v in counts.items())/len(ordered):.4f}"
    )

    print()


# =============================================================================
# TEST 3:
# MULTI-SHIFT INTERSECTION
#
# This is the critical experiment.
# =============================================================================

print("=" * 90)
print("TEST 3: MULTI-SHIFT INTERSECTION")
print("=" * 90)
print()

intersection_results = {}

for r in range(1, len(SHIFTS) + 1):

    active = SHIFTS[:r]

    counts = Counter()

    true_present = 0
    unique = 0
    missing = 0

    examples = []

    for n, p, q in ordered:

        possible = None

        for c in active:

            candidates = single_shift_sets[
                (n, p, q, c)
            ]

            if possible is None:
                possible = set(candidates)
            else:
                possible &= candidates

        possible = possible or set()

        intersection_results[
            (n, p, q, tuple(active))
        ] = possible

        count = len(possible)

        counts[count] += 1

        if p not in possible:

            missing += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append(
                    (
                        n,
                        p,
                        q,
                        sorted(possible),
                    )
                )

        else:

            true_present += 1

            if count == 1:
                unique += 1

    total = len(ordered)

    print(
        f"shifts={active}"
    )

    print(
        f"    true p present="
        f"{true_present}/{total}"
    )

    print(
        f"    unique={unique}"
        f" ambiguous={true_present-unique}"
        f" missing={missing}"
    )

    print(
        f"    unique ratio="
        f"{unique/total:.6f}"
    )

    print(
        f"    mean candidates="
        f"{sum(k*v for k,v in counts.items())/total:.4f}"
    )

    print(
        f"    minimum candidates="
        f"{min(counts)}"
    )

    print(
        f"    maximum candidates="
        f"{max(counts)}"
    )

    print(
        "    candidate-count distribution:"
    )

    # Keep this readable for large runs.
    for count, frequency in sorted(counts.items()):

        if count <= 20 or frequency <= 3:

            print(
                f"        {count:6d} -> {frequency}"
            )

    if examples:

        print(
            "    first missing examples:"
        )

        for (
            n,
            p,
            q,
            values,
        ) in examples:

            print(
                f"        n={n} "
                f"true p={p} "
                f"q={q}"
            )

            print(
                f"            intersection={values}"
            )

    print()


# =============================================================================
# TEST 4:
# HOW MUCH INFORMATION COMES FROM NONTRIVIAL H?
#
# Repeat the intersection while removing H <= 2.
# =============================================================================

print("=" * 90)
print("TEST 4: INTERSECTION USING ONLY H > 2")
print("=" * 90)
print()

for r in range(1, len(SHIFTS) + 1):

    active = SHIFTS[:r]

    unique = 0
    ambiguous = 0
    missing = 0

    candidate_counts = Counter()

    for n, p, q in ordered:

        possible = None

        for c in active:

            M = n + c

            candidates = set()

            for H in divisors(M):

                if H <= 2:
                    continue

                solved = solve_linear_congruence_3p(
                    c,
                    H,
                )

                if solved is None:
                    continue

                residue, modulus = solved

                candidates.update(
                    values_from_residue(
                        residue,
                        modulus,
                        MIN_P + 1,
                        isqrt(n),
                    )
                )

            if possible is None:
                possible = candidates
            else:
                possible &= candidates

        possible = possible or set()

        count = len(possible)

        candidate_counts[count] += 1

        if p not in possible:
            missing += 1
        elif count == 1:
            unique += 1
        else:
            ambiguous += 1

    print(
        f"shifts={active}"
    )

    print(
        f"    true p present="
        f"{len(ordered)-missing}/{len(ordered)}"
    )

    print(
        f"    unique={unique}"
        f" ambiguous={ambiguous}"
        f" missing={missing}"
    )

    print(
        f"    unique ratio="
        f"{unique/len(ordered):.6f}"
    )

    print()


# =============================================================================
# TEST 5:
# TRUE H VS ALL H
#
# This explicitly compares the old hidden-H experiment against the new
# all-divisor experiment.
# =============================================================================

print("=" * 90)
print("TEST 5: TRUE-H CONTROL VS ALL-DIVISOR SEARCH")
print("=" * 90)
print()

for c in SHIFTS:

    true_H_unique = 0
    all_H_unique = 0

    true_H_counts = Counter()
    all_H_counts = Counter()

    for n, p, q in ordered:

        C = q + 3

        true_H = gcd(C, n + c)

        solved = solve_linear_congruence_3p(
            c,
            true_H,
        )

        if solved is None:
            true_values = set()
        else:
            residue, modulus = solved

            true_values = values_from_residue(
                residue,
                modulus,
                MIN_P + 1,
                isqrt(n),
            )

        all_values = single_shift_sets[
            (n, p, q, c)
        ]

        true_H_counts[len(true_values)] += 1
        all_H_counts[len(all_values)] += 1

        if len(true_values) == 1:
            true_H_unique += 1

        if len(all_values) == 1:
            all_H_unique += 1

    print(
        f"c={c:3d}"
    )

    print(
        f"    true-H unique="
        f"{true_H_unique}/{len(ordered)}"
    )

    print(
        f"    all-H unique="
        f"{all_H_unique}/{len(ordered)}"
    )

print()


# =============================================================================
# TEST 6:
# REPRESENTATIVE STATES
# =============================================================================

print("=" * 90)
print("TEST 6: REPRESENTATIVE ALL-DIVISOR CHANNELS")
print("=" * 90)
print()

shown = 0

for n, p, q in ordered:

    if shown >= MAX_EXAMPLES:
        break

    print(
        f"n={n} p={p} q={q} C={q+3}"
    )

    for c in SHIFTS:

        M = n + c
        fac = factorint(M)

        candidates, witnesses = all_divisor_p_candidates(
            n,
            c,
        )

        true_C = q + 3
        true_H = gcd(true_C, M)

        true_solution = solve_linear_congruence_3p(
            c,
            true_H,
        )

        if true_solution is None:
            true_residue = None
        else:
            true_residue = true_solution

        print(
            f"    c={c:3d}"
        )

        print(
            f"        M={M}"
        )

        print(
            f"        factor(M)={fac}"
        )

        print(
            f"        number of divisors="
            f"{len(divisors(M))}"
        )

        print(
            f"        true H={true_H}"
        )

        print(
            f"        true H congruence="
            f"{true_residue}"
        )

        print(
            f"        all-H p candidates="
            f"{len(candidates)}"
        )

        # Show only the first useful witnesses.
        witness_rows = []

        for H, (
            residue,
            modulus,
            values,
        ) in sorted(witnesses.items()):

            if len(witness_rows) >= 8:
                break

            witness_rows.append(
                (
                    H,
                    residue,
                    modulus,
                    len(values),
                )
            )

        for (
            H,
            residue,
            modulus,
            count,
        ) in witness_rows:

            print(
                f"            H={H:6d} "
                f"p≡{residue} "
                f"(mod {modulus}) "
                f"count={count}"
            )

    # Full intersection.
    possible = None

    for c in SHIFTS:

        values = single_shift_sets[
            (n, p, q, c)
        ]

        if possible is None:
            possible = set(values)
        else:
            possible &= values

    possible = possible or set()

    print(
        f"    final intersection="
        f"{sorted(possible)}"
    )

    print(
        f"    true p={p}"
    )

    print("-" * 70)

    shown += 1


# =============================================================================
# TEST 7:
# FINAL EXACT VALIDATION
#
# Only now do we use n % p and primality.
#
# This test answers:
#
#     If all-divisor congruence filtering is performed first,
#     how much exact factor validation remains?
# =============================================================================

print("=" * 90)
print("TEST 7: FINAL EXACT VALIDATION")
print("=" * 90)
print()

active = SHIFTS

unique_before_validation = 0
ambiguous_before_validation = 0

exact_unique_after_validation = 0
exact_ambiguous_after_validation = 0
exact_missing = 0

final_candidate_distribution = Counter()

examples = []

for n, p, q in ordered:

    possible = None

    for c in active:

        values = single_shift_sets[
            (n, p, q, c)
        ]

        if possible is None:
            possible = set(values)
        else:
            possible &= values

    possible = possible or set()

    if len(possible) == 1:
        unique_before_validation += 1
    else:
        ambiguous_before_validation += 1

    validated = []

    for candidate_p in sorted(possible):

        # This is intentionally the FIRST place where the exact
        # factor relation n=pq is tested.

        if not factorint(candidate_p):
            continue

        if not (
            len(factorint(candidate_p)) == 1
            and list(
                factorint(candidate_p).values()
            )[0] == 1
        ):
            continue

        if candidate_p <= MIN_P:
            continue

        if n % candidate_p != 0:
            continue

        candidate_q = n // candidate_p

        if candidate_q <= candidate_p:
            continue

        if not (
            len(factorint(candidate_q)) == 1
            and list(
                factorint(candidate_q).values()
            )[0] == 1
        ):
            continue

        validated.append(
            (
                candidate_p,
                candidate_q,
            )
        )

    final_candidate_distribution[
        len(validated)
    ] += 1

    true_pair = (p, q)

    if true_pair not in validated:

        exact_missing += 1

        if len(examples) < MAX_EXAMPLES:

            examples.append(
                (
                    n,
                    p,
                    q,
                    sorted(validated),
                    sorted(possible),
                )
            )

    elif len(validated) == 1:

        exact_unique_after_validation += 1

    else:

        exact_ambiguous_after_validation += 1

print(
    f"unique before exact validation="
    f"{unique_before_validation}"
)

print(
    f"ambiguous before exact validation="
    f"{ambiguous_before_validation}"
)

print(
    f"exact unique after validation="
    f"{exact_unique_after_validation}"
)

print(
    f"exact ambiguous after validation="
    f"{exact_ambiguous_after_validation}"
)

print(
    f"exact missing="
    f"{exact_missing}"
)

print(
    "final validated candidate distribution:"
)

for count, frequency in sorted(
    final_candidate_distribution.items()
):

    print(
        f"    {count:4d} -> {frequency}"
    )

if examples:

    print(
        "validation failures:"
    )

    for (
        n,
        p,
        q,
        validated,
        before,
    ) in examples:

        print(
            f"    n={n} true=({p},{q})"
        )

        print(
            f"        before={before}"
        )

        print(
            f"        validated={validated}"
        )

print()


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("=" * 90)
print("FINAL STRUCTURAL SUMMARY")
print("=" * 90)

print(
"""
Experiment 696 removes the hidden-H assumption from Experiment 695.

For each shift:

    M_c = n+c

we factor M_c and enumerate every divisor:

    H | M_c.

For every such H we solve:

    3p ≡ c (mod H).

Thus the observable channel is:

    factor(n+c)
        |
        v
    all H | (n+c)
        |
        v
    p ≡ r_H (mod m_H)
        |
        v
    p > 137
    p <= sqrt(n)
        |
        v
    candidate p set.

For multiple shifts:

    P_3
    P_81
    P_137

we compute:

    P = P_3 ∩ P_81 ∩ P_137.

This is the critical test.

Previously:

    true H
        |
        v
    congruence
        |
        v
    p

was tested.

Now:

    all observable H
        |
        v
    all possible congruences
        |
        v
    intersection
        |
        v
    p

is tested.

Therefore, if the intersection becomes small while the true p
remains present, the divisor lattices themselves are supplying
nontrivial information about the hidden factor.

The final exact factor check is deliberately isolated at the end.

The strongest possible outcome would be:

    all-divisor H channels
        ->
    unique p
        ->
    q=n/p.

A weaker outcome is:

    all-divisor H channels
        ->
    small p candidate set
        ->
    one final exact divisibility test.

A negative outcome is:

    all-divisor H channels
        ->
    huge p candidate set.

That would show that the congruence identity is algebraically true
but does not provide enough information for factor recovery by itself.

The essential comparison is:

    TRUE-H candidate count
        versus
    ALL-H candidate count.

Only the second one measures the actual observable attack/recovery
channel.
"""
)

print()
print(f"BASELINE FAILURES={failures}")
print("=" * 90)
print("EXPERIMENT 696 FINISHED")
print("=" * 90)

n = 2140324650240744961264423072839333563008614715144755017797754920881418023447140136643345519095804679610992851872470914587687396261921557363047454770520805119056493106687691590019759405693457452230589325976697471681738069364894699871578494975937497937
p = 64135289477071580278790190170577389084825014742943447208116859632024532344630238623598752668347708737661925585694639798853367
q = 33372027594978156556226010605355114227940760344767554666784520987023841729210037080257448673296881877565718986258036932062711

for c in range(3, 100, 2):
    print(f"{c} {gcd(q + 3, n + c)} {gcd(q + 3, c - 3*p)}")
