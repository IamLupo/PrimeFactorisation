#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 695
==========================================================================================

DIVISOR-LATTICE -> C RECOVERY WITHOUT q|n SEARCH

Previous result:

    H_c | (q-c)

was FALSE.

The exact identity for Frame A is:

    C = q + 3

    H_c = gcd(C, n+c)

and therefore:

    H_c | C
    H_c | (c - 3p)

because

    n+c = pq+c
         = p(C-3)+c
         = pC + (c-3).

Hence:

    H_c | (c - 3p).

This experiment tests the corrected channel.

For every divisor H of n+c we generate possible C values:

    C = kH

and impose the frame relation:

    p = (c - tH) / 3

for some integer t,

because:

    H | (c - 3p).

Then:

    q = C - 3
    n = p*q

The crucial point:

    NO scan over q
    NO q | n search over an interval

The only exact factor validation is performed after the candidate
C has been constructed from the divisor lattice.

The experiment measures:

    1. how many H divisors are useful;
    2. how many C candidates each shift generates;
    3. how much [3,81,137] reduces the C space;
    4. whether the true C is always generated;
    5. how often the C candidate becomes unique before exact n
       divisibility is checked;
    6. whether p > 137 materially helps.

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

def prime_semiprimes(limit):
    primes = list(primerange(3, limit))

    for i, p in enumerate(primes):
        for q in primes[i:]:
            yield p * q, p, q


def frame_a_C(p, q):
    return q + 3


def valid_prime_factorization(n, p, q):
    if p <= MIN_P:
        return False

    if p >= q:
        return False

    if p * q != n:
        return False

    fp = factorint(p)
    fq = factorint(q)

    return (
        len(fp) == 1
        and list(fp.values())[0] == 1
        and len(fq) == 1
        and list(fq.values())[0] == 1
    )


# =============================================================================
# CANDIDATE GENERATION
# =============================================================================

def generate_C_from_H(n, c):
    """
    Generate C candidates from the complete divisor lattice of M=n+c.

    We use:

        H | C

    and:

        H | (c - 3p)

    so:

        c - 3p = H*t

    and therefore:

        p = (c - H*t)/3.

    Since C=q+3 and p*q=n:

        q = C-3
        p = n/(C-3)

    Candidate generation starts from H and does NOT scan q.

    Instead, C is generated as a multiple of H in the frame range.

    For each candidate C we retain only the algebraically compatible
    p-congruence class.

    The exact n divisibility test is deliberately separated from
    the first candidate-generation stage.
    """

    M = n + c

    if M <= 0:
        return set(), {}

    H_values = divisors(M)

    # q > p > MIN_P
    #
    # q > sqrt(n)
    #
    # therefore:
    #
    # C = q+3 > sqrt(n)+3
    #
    # and:
    #
    # q < n/MIN_P
    #
    # therefore:
    #
    # C < n/MIN_P + 3

    C_min = isqrt(n) + 4
    C_max = n // MIN_P + 3

    candidates = set()
    witnesses = {}

    for H in H_values:

        if H <= 0:
            continue

        # C must be a positive multiple of H.
        k_min = (C_min + H - 1) // H
        k_max = C_max // H

        for k in range(k_min, k_max + 1):

            C = k * H
            q = C - 3

            if q <= 0:
                continue

            # -------------------------------------------------------------
            # Algebraic condition:
            #
            # H | (c - 3p)
            #
            # but p = n/q is NOT used yet.
            #
            # We therefore record H as a divisor candidate for C.
            #
            # The next condition tests whether some integer p satisfying
            #
            #     c - 3p = H*t
            #
            # exists.
            #
            # This means:
            #
            #     3p == c (mod H)
            #
            # -------------------------------------------------------------

            g = gcd(3, H)

            if c % g != 0:
                continue

            # Solve:
            #
            #     3p = c (mod H)
            #
            # Reduce:
            #
            #     (3/g)p = c/g (mod H/g)
            #
            mod = H // g
            a = 3 // g
            b = c // g

            if mod == 1:
                residue = 0
            else:
                inv = pow(a, -1, mod)
                residue = (b * inv) % mod

            candidates.add(C)

            witnesses.setdefault(
                C,
                []
            ).append(
                {
                    "H": H,
                    "p_residue": residue,
                    "p_modulus": mod,
                }
            )

    return candidates, witnesses


# =============================================================================
# START
# =============================================================================

print("=" * 90)
print("EXPERIMENT 695 START")
print("=" * 90)
print()

print(f"prime limit={PRIME_LIMIT}")
print(f"minimum p={MIN_P}")
print(f"shifts={SHIFTS}")

states = list(prime_semiprimes(PRIME_LIMIT))

print(f"semiprimes={len(states)}")
print()


# =============================================================================
# ORDERED STATE SUBSET
# =============================================================================

ordered_states = [
    (n, p, q)
    for n, p, q in states
    if q > p > MIN_P
]

print(f"ordered states={len(ordered_states)}")
print()


# =============================================================================
# TEST 0: BASELINE
# =============================================================================

print("=" * 90)
print("TEST 0: BASELINE")
print("=" * 90)

baseline_failures = 0

for n, p, q in states:

    if p * q != n:
        baseline_failures += 1

print(f"checked={len(states)}")
print(f"failures={baseline_failures}")
print()


# =============================================================================
# TEST 1:
# VERIFY THE CORRECT ALGEBRAIC CHANNEL
#
# H | C
# H | (c-3p)
# =============================================================================

print("=" * 90)
print("TEST 1: CORRECT ALGEBRAIC CHANNEL")
print("=" * 90)
print()

for c in SHIFTS:

    total = 0
    failures = 0

    for n, p, q in ordered_states:

        C = q + 3
        H = gcd(C, n + c)

        ok1 = (C % H == 0)
        ok2 = ((c - 3 * p) % H == 0)

        total += 1

        if not (ok1 and ok2):
            failures += 1

    print(
        f"c={c:3d}: "
        f"H|C AND H|(c-3p): "
        f"{total-failures}/{total}"
    )

print()


# =============================================================================
# TEST 2:
# CANDIDATE GENERATION FROM ONE SHIFT
# =============================================================================

print("=" * 90)
print("TEST 2: ONE SHIFT -> C CANDIDATES")
print("=" * 90)
print()

single_shift_stats = {}

for c in SHIFTS:

    total = len(ordered_states)

    true_present = 0
    unique = 0
    ambiguous = 0
    missing = 0

    count_distribution = Counter()

    examples = []

    for n, p, q in ordered_states:

        true_C = q + 3

        candidates, witnesses = generate_C_from_H(
            n=n,
            c=c,
        )

        count = len(candidates)
        count_distribution[count] += 1

        if true_C not in candidates:
            missing += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append(
                    (
                        n,
                        p,
                        q,
                        true_C,
                        sorted(candidates),
                    )
                )

        else:

            true_present += 1

            if count == 1:
                unique += 1
            else:
                ambiguous += 1

    single_shift_stats[c] = {
        "unique": unique,
        "ambiguous": ambiguous,
        "missing": missing,
        "distribution": count_distribution,
    }

    print(f"c={c}")
    print(f"    true C present={true_present}/{total}")
    print(f"    unique={unique}")
    print(f"    ambiguous={ambiguous}")
    print(f"    missing={missing}")

    if total:
        print(
            f"    unique ratio={unique/total:.6f}"
        )

    print("    candidate-count distribution:")

    for count, frequency in sorted(
        count_distribution.items()
    ):
        print(
            f"        {count:6d} -> {frequency}"
        )

    if examples:

        print("    missing examples:")

        for n, p, q, C, candidates in examples:

            print(
                f"        n={n} p={p} q={q} "
                f"C={C}"
            )

            print(
                f"            candidates={candidates}"
            )

    print()


# =============================================================================
# TEST 3:
# MULTI-SHIFT INTERSECTION
# =============================================================================

print("=" * 90)
print("TEST 3: MULTI-SHIFT INTERSECTION")
print("=" * 90)
print()

for r in range(1, len(SHIFTS) + 1):

    active_shifts = SHIFTS[:r]

    total = len(ordered_states)

    unique = 0
    ambiguous = 0
    missing = 0

    count_distribution = Counter()

    examples = []

    for n, p, q in ordered_states:

        true_C = q + 3

        possible = None

        for c in active_shifts:

            candidates, _ = generate_C_from_H(
                n=n,
                c=c,
            )

            if possible is None:
                possible = set(candidates)

            else:
                possible &= candidates

        possible = possible or set()

        count_distribution[len(possible)] += 1

        if true_C not in possible:

            missing += 1

            if len(examples) < MAX_EXAMPLES:

                examples.append(
                    (
                        n,
                        p,
                        q,
                        true_C,
                        sorted(possible),
                    )
                )

        elif len(possible) == 1:

            unique += 1

        else:

            ambiguous += 1

    print(
        f"shifts={active_shifts}"
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

    print("    candidate-count distribution:")

    for count, frequency in sorted(
        count_distribution.items()
    ):
        print(
            f"        {count:6d} -> {frequency}"
        )

    if examples:

        print("    first unresolved examples:")

        for n, p, q, C, candidates in examples:

            print(
                f"        n={n} p={p} q={q} "
                f"C={C}"
            )

            print(
                f"            candidates={candidates}"
            )

    print()


# =============================================================================
# TEST 4:
# REMOVE TRIVIAL H=1 AND H=2 CHANNELS
#
# This is important because H=1 gives no information and H=2 is
# automatically compatible with every odd C/p structure.
# =============================================================================

print("=" * 90)
print("TEST 4: NONTRIVIAL H CHANNEL")
print("=" * 90)
print()

for c in SHIFTS:

    total = 0
    usable = 0
    unique = 0
    ambiguous = 0
    missing = 0

    for n, p, q in ordered_states:

        true_C = q + 3

        M = n + c

        H_values = [
            H
            for H in divisors(M)
            if H > 2
        ]

        candidates = set()

        C_min = isqrt(n) + 4
        C_max = n // MIN_P + 3

        for H in H_values:

            k_min = (C_min + H - 1) // H
            k_max = C_max // H

            for k in range(k_min, k_max + 1):

                C = k * H

                g = gcd(3, H)

                if c % g != 0:
                    continue

                candidates.add(C)

        total += 1

        if true_C in candidates:

            usable += 1

            if len(candidates) == 1:
                unique += 1
            else:
                ambiguous += 1

        else:
            missing += 1

    print(
        f"c={c:3d}: "
        f"true C present={usable}/{total}, "
        f"unique={unique}, "
        f"ambiguous={ambiguous}, "
        f"missing={missing}"
    )

print()


# =============================================================================
# TEST 5:
# TRUE C'S ACTUAL H VALUES
#
# Show how much of C is exposed by H.
# =============================================================================

print("=" * 90)
print("TEST 5: PRIME-POWER EXPOSURE OF TRUE C")
print("=" * 90)
print()

for c in SHIFTS:

    exposure = Counter()

    examples = []

    for n, p, q in ordered_states:

        C = q + 3
        H = gcd(C, n + c)

        if H == C:
            kind = "full"
        elif H == 1:
            kind = "none"
        else:
            kind = "partial"

        exposure[kind] += 1

        if len(examples) < MAX_EXAMPLES and kind != "full":

            examples.append(
                (
                    n,
                    p,
                    q,
                    C,
                    H,
                    C // H,
                    kind,
                )
            )

    print(f"c={c}")
    print(
        f"    full    ={exposure['full']}"
    )
    print(
        f"    partial ={exposure['partial']}"
    )
    print(
        f"    none    ={exposure['none']}"
    )

    if examples:

        print("    examples:")

        for (
            n,
            p,
            q,
            C,
            H,
            residual,
            kind,
        ) in examples:

            print(
                f"        n={n} p={p} q={q} "
                f"C={C} H={H} "
                f"C/H={residual} "
                f"{kind}"
            )

    print()


# =============================================================================
# TEST 6:
# CAN THE DIVISOR H DIRECTLY RECONSTRUCT p?
#
# We now use:
#
#     H | (c - 3p)
#
# and the known bound:
#
#     MIN_P < p <= sqrt(n)
#
# to generate p congruence classes.
#
# No n%q search is performed here.
# =============================================================================

print("=" * 90)
print("TEST 6: H -> p CONGRUENCE CHANNEL")
print("=" * 90)
print()

for c in SHIFTS:

    total = 0
    true_residue_found = 0
    unique_p_residue = 0

    for n, p, q in ordered_states:

        C = q + 3
        H = gcd(C, n + c)

        g = gcd(3, H)

        if c % g != 0:
            total += 1
            continue

        mod = H // g

        if mod == 1:
            residue = 0
        else:
            a = 3 // g
            b = c // g

            residue = (
                b * pow(a, -1, mod)
            ) % mod

        total += 1

        p_min = MIN_P + 1
        p_max = isqrt(n)

        values = []

        first = p_min

        if first % mod != residue:
            first += (
                residue - first
            ) % mod

        for candidate_p in range(
            first,
            p_max + 1,
            mod
        ):
            values.append(candidate_p)

        if p in values:
            true_residue_found += 1

        if len(values) == 1:
            unique_p_residue += 1

    print(
        f"c={c:3d}: "
        f"true p residue present="
        f"{true_residue_found}/{total}"
    )

    print(
        f"        unique p in congruence "
        f"interval={unique_p_residue}"
    )

print()


# =============================================================================
# TEST 7:
# MULTI-SHIFT p CONGRUENCE INTERSECTION
#
# Here the useful object is:
#
#     p ≡ r_c (mod m_c)
#
# for several c.
# =============================================================================

print("=" * 90)
print("TEST 7: MULTI-SHIFT p CONGRUENCE INTERSECTION")
print("=" * 90)
print()

total = len(ordered_states)

unique_p = 0
ambiguous_p = 0
missing_p = 0

examples = []

for n, p, q in ordered_states:

    possible = set(
        range(
            MIN_P + 1,
            isqrt(n) + 1
        )
    )

    for c in SHIFTS:

        C = q + 3
        H = gcd(C, n + c)

        g = gcd(3, H)

        if c % g != 0:
            possible.clear()
            break

        mod = H // g

        if mod == 1:
            continue

        a = 3 // g
        b = c // g

        residue = (
            b * pow(a, -1, mod)
        ) % mod

        possible = {
            candidate_p
            for candidate_p in possible
            if candidate_p % mod == residue
        }

        if not possible:
            break

    # Primality and exact factor relation only AFTER the congruence stage.
    valid = set()

    for candidate_p in possible:

        fac = factorint(candidate_p)

        if (
            len(fac) == 1
            and list(fac.values())[0] == 1
        ):
            valid.add(candidate_p)

    if p not in valid:

        missing_p += 1

        if len(examples) < MAX_EXAMPLES:
            examples.append(
                (
                    n,
                    p,
                    q,
                    sorted(valid),
                )
            )

    elif len(valid) == 1:

        unique_p += 1

    else:

        ambiguous_p += 1

print(f"states={total}")
print(f"true p recovered={total-missing_p}")
print(f"missing={missing_p}")
print(f"unique p={unique_p}")
print(f"ambiguous p={ambiguous_p}")

if total:
    print(
        f"unique ratio={unique_p/total:.6f}"
    )

if examples:

    print("missing examples:")

    for n, p, q, values in examples:

        print(
            f"    n={n} p={p} q={q}"
            f" valid_p={values}"
        )

print()


# =============================================================================
# TEST 8:
# REPRESENTATIVE STATES
# =============================================================================

print("=" * 90)
print("TEST 8: REPRESENTATIVE DIVISOR CHANNELS")
print("=" * 90)
print()

printed = 0

for n, p, q in ordered_states:

    if printed >= MAX_EXAMPLES:
        break

    C = q + 3

    print(
        f"n={n} p={p} q={q} C={C}"
    )

    for c in SHIFTS:

        M = n + c
        fac = factorint(M)

        H = gcd(C, M)

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
            f"        H={H}"
        )
        print(
            f"        H|C={C % H == 0}"
        )
        print(
            f"        H|(c-3p)={(c-3*p) % H == 0}"
        )

        g = gcd(3, H)

        if c % g == 0:

            mod = H // g

            if mod == 1:
                residue = 0
            else:
                a = 3 // g
                b = c // g
                residue = (
                    b * pow(a, -1, mod)
                ) % mod

            print(
                f"        p ≡ {residue} "
                f"(mod {mod})"
            )

    print("-" * 60)

    printed += 1


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("=" * 90)
print("FINAL STRUCTURAL SUMMARY")
print("=" * 90)

print(
    """
The previous experiment established an important correction:

    H_c | (q-c)

is not generally true.

The exact Frame-A identity is:

    C = q+3

    n+c = p(C-3)+c
        = pC + (c-3p)

therefore:

    H_c = gcd(C,n+c)

implies:

    H_c | C
    H_c | (c-3p).

This experiment removes the q-congruence assumption.

The new channel is:

    factor(n+c)
         |
         v
    all H | (n+c)
         |
         v
    H | C
    H | (c-3p)
         |
         v
    C = kH
         |
         v
    p ≡ solution of 3p ≡ c (mod H)
         |
         v
    p > 137
    p <= sqrt(n)
         |
         v
    candidate p / C states

The key measurement is no longer:

    "Can I find q by testing q | n?"

because that can hide the actual source of the factorization.

Instead we measure:

    1. How many C values are generated from H-divisors?
    2. Does the true C always occur?
    3. How many C candidates survive one shift?
    4. How many survive [3,81,137]?
    5. How much does the p-congruence shrink the search?
    6. Can the p-congruence become unique before exact n-factor
       validation?

A particularly important outcome would be:

    divisor lattice
        ->
    small p-congruence class
        ->
    unique prime p
        ->
    q=n/p.

That would be a substantially cleaner information chain than the
previous experiment because it avoids using q|n as the main search
mechanism.
"""
)

print()
print(f"BASELINE FAILURES={baseline_failures}")
print("=" * 90)
print("EXPERIMENT 695 FINISHED")
print("=" * 90)
