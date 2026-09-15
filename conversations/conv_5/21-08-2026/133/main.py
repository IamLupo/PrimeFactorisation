#!/usr/bin/env python3
"""
==========================================================================================
EXPERIMENT 694
==========================================================================================

H_c -> q-c DIVISIBILITY + ORDERED FACTOR RECOVERY

Hypothesis under test:

    Frame A:
        C = q + 3

    and for shifted channels:

        H_3   | (q - 3)
        H_81  | (q - 81)
        H_137 | (q - 137)

where:

        H_c = gcd(C, n+c)

The experiment tests four things:

    1. Whether the proposed divisibility identities are actually true.
    2. Whether p > 137 is sufficient to make the q-search useful.
    3. Whether H-divisors of n+c can generate candidate q values
       without using the known factorization.
    4. Whether combining c=3,81,137 produces a unique q.

The known p,q are used ONLY for validation.

Candidate generation is NOT allowed to use:
    p
    q
    C
    true H

==========================================================================================
"""

from math import gcd, isqrt
from collections import Counter, defaultdict
from sympy import primerange, factorint, divisors


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 400

SHIFTS = [3, 81, 137]

# Ordering restriction proposed in the research:
#     q > p > 137
MIN_P = 137

# Maximum number of printed examples per category
MAX_EXAMPLES = 20


# =============================================================================
# BASIC FUNCTIONS
# =============================================================================

def v2(x: int) -> int:
    """Return v_2(x) for nonzero positive integer x."""
    if x == 0:
        return 10**9

    x = abs(x)
    r = 0

    while x % 2 == 0:
        x //= 2
        r += 1

    return r


def prime_semiprimes(limit: int):
    """
    Generate ordered semiprimes

        n = p*q
        p <= q

    from primes < limit.
    """
    primes = list(primerange(2, limit))

    for i, p in enumerate(primes):
        for q in primes[i:]:
            yield p * q, p, q


def frame_a_C(p: int, q: int) -> int:
    """
    Frame-A canonical C from the current hypothesis.
    """
    return q + 3


def candidate_q_values_from_H_divisors(n: int, c: int, p_min: int):
    """
    Candidate generation using ONLY:

        M = n+c
        factor(M)
        divisors(M)

    and the hypothesis:

        H_c | (q-c)

    Hence:

        q = c + k*H_c.

    We additionally impose:

        q > sqrt(n)
        q < n/p_min

    because q > p > p_min.

    No hidden p or q is used.
    """

    M = n + c

    if M <= 0:
        return set(), {}

    fac = factorint(M)
    H_values = divisors(M)

    q_min = isqrt(n) + 1
    q_max = (n - 1) // p_min

    candidates = set()

    witness = {}

    for H in H_values:
        if H <= 0:
            continue

        # q = c mod H
        # q = c + kH

        first = c

        if first < q_min:
            k = (q_min - first + H - 1) // H
            first += k * H

        for q in range(first, q_max + 1, H):
            if q <= p_min:
                continue

            # q itself must be prime.
            # factorint(q) is deliberately avoided here:
            # trial divisibility through the configured prime list
            # is sufficient for these experiment sizes.
            qfac = factorint(q)

            if len(qfac) != 1 or list(qfac.values())[0] != 1:
                continue

            if n % q != 0:
                continue

            candidates.add(q)

            if q not in witness:
                witness[q] = H

    return candidates, {
        "factorization": fac,
        "divisors": H_values,
        "witness": witness,
        "q_min": q_min,
        "q_max": q_max,
    }


# =============================================================================
# START
# =============================================================================

print("=" * 90)
print("EXPERIMENT 694 START")
print("=" * 90)
print()
print(f"prime limit={PRIME_LIMIT}")
print(f"shifts={SHIFTS}")
print(f"minimum p={MIN_P}")
print()

primes = list(primerange(2, PRIME_LIMIT))
print(f"odd primes={len([p for p in primes if p > 2])}")

states = list(prime_semiprimes(PRIME_LIMIT))
print(f"semiprimes={len(states)}")

print()


# =============================================================================
# TEST 0: BASELINE
# =============================================================================

print("=" * 90)
print("TEST 0: BASELINE")
print("=" * 90)

baseline_failures = 0

for n, p, q in states:
    if n != p * q:
        baseline_failures += 1

print(f"checked={len(states)}")
print(f"failures={baseline_failures}")
print()


# =============================================================================
# TEST 1:
# ACTUAL H_c VS HYPOTHESIZED q-c DIVISIBILITY
# =============================================================================

print("=" * 90)
print("TEST 1: DOES H_c DIVIDE (q-c)?")
print("=" * 90)
print()

relation_stats = {}

for c in SHIFTS:
    ok = 0
    fail = 0
    examples_fail = []

    for n, p, q in states:

        if p <= MIN_P:
            continue

        C = q + 3
        H = gcd(C, n + c)

        if (q - c) % H == 0:
            ok += 1
        else:
            fail += 1

            if len(examples_fail) < MAX_EXAMPLES:
                examples_fail.append(
                    (n, p, q, C, c, H, q - c)
                )

    relation_stats[c] = (ok, fail)

    print(
        f"c={c:3d}: "
        f"H_c | (q-c) = {ok}/{ok+fail} "
        f"failures={fail}"
    )

    if examples_fail:
        print("  counterexamples:")
        for n, p, q, C, c, H, target in examples_fail:
            print(
                f"      n={n} p={p} q={q} "
                f"C={C} H={H} q-c={target}"
            )

print()


# =============================================================================
# TEST 2:
# ALGEBRAIC RELATION CHANNEL
# =============================================================================

print("=" * 90)
print("TEST 2: ALGEBRAIC GCD RELATION")
print("=" * 90)
print()

for c in SHIFTS:

    exact = 0
    failures = 0
    examples = []

    for n, p, q in states:

        if p <= MIN_P:
            continue

        C = q + 3

        H = gcd(C, n + c)

        # Since C=q+3:
        #
        # n+c = pq+c
        #
        # modulo C:
        #
        # q = -3
        #
        # n+c = -3p+c
        #
        # therefore:
        #
        # H | (c - 3p).

        predicted = abs(c - 3 * p)

        # gcd(C,n+c) must equal gcd(C,c-3p).
        rhs = gcd(C, predicted)

        if H == rhs:
            exact += 1
        else:
            failures += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append(
                    (n, p, q, C, c, H, rhs)
                )

    total = exact + failures

    print(
        f"c={c:3d}: "
        f"H = gcd(C,c-3p): "
        f"{exact}/{total} exact"
    )

    if examples:
        print("  failures:")
        for item in examples:
            print(
                f"      n={item[0]} p={item[1]} q={item[2]} "
                f"C={item[3]} c={item[4]} "
                f"H={item[5]} rhs={item[6]}"
            )

print()


# =============================================================================
# TEST 3:
# ORDERING FILTER q > p > 137
# =============================================================================

print("=" * 90)
print("TEST 3: ORDERING FILTER q > p > 137")
print("=" * 90)
print()

ordered_states = [
    (n, p, q)
    for n, p, q in states
    if q > p > MIN_P
]

print(f"ordered states={len(ordered_states)}")

if ordered_states:
    print(
        f"smallest p={min(p for _, p, _ in ordered_states)}"
    )
    print(
        f"largest p={max(p for _, p, _ in ordered_states)}"
    )
    print(
        f"smallest q={min(q for _, _, q in ordered_states)}"
    )
    print(
        f"largest q={max(q for _, _, q in ordered_states)}"
    )

print()


# =============================================================================
# TEST 4:
# SINGLE SHIFT -> q RECOVERY USING H DIVISORS
# =============================================================================

print("=" * 90)
print("TEST 4: SINGLE-SHIFT H-LATTICE -> q RECOVERY")
print("=" * 90)
print()

for c in SHIFTS:

    total = 0
    recovered = 0
    unique = 0
    ambiguous = 0
    missing = 0

    candidate_counts = Counter()
    examples = []

    for n, p, q in ordered_states:

        total += 1

        candidates, info = candidate_q_values_from_H_divisors(
            n=n,
            c=c,
            p_min=MIN_P,
        )

        candidate_counts[len(candidates)] += 1

        if q in candidates:
            recovered += 1

            if len(candidates) == 1:
                unique += 1
            else:
                ambiguous += 1

        else:
            missing += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append(
                    (n, p, q, candidates)
                )

    print(f"c={c}")
    print(f"  states={total}")
    print(f"  true q recovered={recovered}")
    print(f"  missing={missing}")
    print(f"  unique={unique}")
    print(f"  ambiguous={ambiguous}")

    if total:
        print(f"  unique ratio={unique/total:.6f}")

    print("  candidate-count distribution:")
    for count, frequency in sorted(candidate_counts.items()):
        print(f"      {count:4d} -> {frequency}")

    if examples:
        print("  missing examples:")
        for n, p, q, candidates in examples:
            print(
                f"      n={n} p={p} q={q} "
                f"candidates={sorted(candidates)}"
            )

    print()


# =============================================================================
# TEST 5:
# MULTI-SHIFT INTERSECTION
# =============================================================================

print("=" * 90)
print("TEST 5: MULTI-SHIFT INTERSECTION -> q")
print("=" * 90)
print()

multi_total = 0
multi_recovered = 0
multi_unique = 0
multi_ambiguous = 0
multi_missing = 0

candidate_count_distribution = Counter()
multi_examples = []

for n, p, q in ordered_states:

    multi_total += 1

    possible = None
    all_single = {}

    for c in SHIFTS:

        candidates, info = candidate_q_values_from_H_divisors(
            n=n,
            c=c,
            p_min=MIN_P,
        )

        all_single[c] = candidates

        if possible is None:
            possible = set(candidates)
        else:
            possible &= candidates

    possible = possible or set()

    candidate_count_distribution[len(possible)] += 1

    if q in possible:
        multi_recovered += 1

        if len(possible) == 1:
            multi_unique += 1
        else:
            multi_ambiguous += 1
    else:
        multi_missing += 1

        if len(multi_examples) < MAX_EXAMPLES:
            multi_examples.append(
                (n, p, q, all_single, possible)
            )

print(f"states={multi_total}")
print(f"true q recovered={multi_recovered}")
print(f"missing={multi_missing}")
print(f"unique={multi_unique}")
print(f"ambiguous={multi_ambiguous}")

if multi_total:
    print(f"unique ratio={multi_unique/multi_total:.6f}")

print("candidate-count distribution:")
for count, frequency in sorted(
    candidate_count_distribution.items()
):
    print(f"    {count:4d} -> {frequency}")

if multi_examples:
    print()
    print("missing examples:")

    for n, p, q, all_single, possible in multi_examples:
        print(
            f"    n={n} p={p} q={q} "
            f"intersection={sorted(possible)}"
        )

        for c in SHIFTS:
            print(
                f"        c={c:3d}: "
                f"{sorted(all_single[c])}"
            )

print()


# =============================================================================
# TEST 6:
# DIRECT q-CONGRUENCE CONSISTENCY
# =============================================================================

print("=" * 90)
print("TEST 6: H-VECTOR q-CONGRUENCE CONSISTENCY")
print("=" * 90)
print()

signature_counts = Counter()
signature_examples = defaultdict(list)

for n, p, q in ordered_states:

    hs = []

    for c in SHIFTS:
        C = q + 3
        H = gcd(C, n + c)

        hs.append(H)

    signature = tuple(hs)

    signature_counts[signature] += 1

    if len(signature_examples[signature]) < 3:
        signature_examples[signature].append(
            (n, p, q)
        )

print(f"distinct H-vectors={len(signature_counts)}")

multiplicity = Counter(signature_counts.values())

print("H-vector multiplicity:")
for m, count in sorted(multiplicity.items()):
    print(f"    {m:4d} -> {count}")

print()


# =============================================================================
# TEST 7:
# CAN THE H-VECTOR GIVE q DIRECTLY?
#
# For each actual H_c, test:
#
#     q = c (mod H_c)
#
# and combine all congruences using brute-force CRT-compatible
# enumeration inside the ordered q interval.
# =============================================================================

print("=" * 90)
print("TEST 7: H-VECTOR -> ORDERED q INTERVAL")
print("=" * 90)
print()

direct_total = 0
direct_unique = 0
direct_ambiguous = 0
direct_missing = 0

direct_counts = Counter()
direct_examples = []

for n, p, q in ordered_states:

    direct_total += 1

    q_min = isqrt(n) + 1
    q_max = (n - 1) // MIN_P

    # Start with entire ordered interval.
    possible = set(range(q_min, q_max + 1))

    for c in SHIFTS:

        C = q + 3
        H = gcd(C, n + c)

        if H <= 1:
            continue

        possible = {
            x
            for x in possible
            if (x - c) % H == 0
        }

        if not possible:
            break

    # Prime + factor test is only used after the congruence filtering.
    valid = set()

    for candidate_q in possible:

        if candidate_q <= MIN_P:
            continue

        fac = factorint(candidate_q)

        if len(fac) != 1 or list(fac.values())[0] != 1:
            continue

        if n % candidate_q != 0:
            continue

        candidate_p = n // candidate_q

        if not (candidate_p > MIN_P):
            continue

        if candidate_p >= candidate_q:
            continue

        valid.add(candidate_q)

    direct_counts[len(valid)] += 1

    if q in valid:

        if len(valid) == 1:
            direct_unique += 1
        else:
            direct_ambiguous += 1

    else:
        direct_missing += 1

        if len(direct_examples) < MAX_EXAMPLES:
            direct_examples.append(
                (n, p, q, sorted(valid), sorted(possible))
            )

print(f"states={direct_total}")
print(f"true q recovered={direct_total-direct_missing}")
print(f"missing={direct_missing}")
print(f"unique={direct_unique}")
print(f"ambiguous={direct_ambiguous}")

if direct_total:
    print(f"unique ratio={direct_unique/direct_total:.6f}")

print("candidate-count distribution:")
for count, frequency in sorted(direct_counts.items()):
    print(f"    {count:4d} -> {frequency}")

print()


# =============================================================================
# TEST 8:
# REPRESENTATIVE EXAMPLES
# =============================================================================

print("=" * 90)
print("TEST 8: REPRESENTATIVE ORDERED STATES")
print("=" * 90)
print()

printed = 0

for n, p, q in ordered_states:

    if printed >= MAX_EXAMPLES:
        break

    C = q + 3

    print(f"n={n} p={p} q={q} C=q+3={C}")
    print(f"sqrt(n)={isqrt(n)}")
    print()

    for c in SHIFTS:

        M = n + c
        fac = factorint(M)
        divs = divisors(M)

        H = gcd(C, M)

        relation = ((q - c) % H == 0)

        print(
            f"    c={c:3d}"
            f"  M={M}"
            f"  factor(M)={fac}"
            f"  H={H}"
            f"  H|(q-c)={relation}"
            f"  q-c={q-c}"
        )

    print()

    # Candidate generation using the proposed rule.
    for c in SHIFTS:
        candidates, _ = candidate_q_values_from_H_divisors(
            n=n,
            c=c,
            p_min=MIN_P,
        )

        print(
            f"    candidate q from c={c}: "
            f"{sorted(candidates)}"
        )

    possible = None

    for c in SHIFTS:
        candidates, _ = candidate_q_values_from_H_divisors(
            n=n,
            c=c,
            p_min=MIN_P,
        )

        if possible is None:
            possible = set(candidates)
        else:
            possible &= candidates

    print(
        f"    multi-shift intersection: "
        f"{sorted(possible or set())}"
    )

    print("=" * 60)

    printed += 1


# =============================================================================
# TEST 9:
# IMPORTANT COUNTERCHECK
#
# The user's proposed relation:
#
#     H_c | (q-c)
#
# is NOT automatically implied by:
#
#     H_c = gcd(q+3,n+c).
#
# This test looks for the first failure in the ordered domain.
# =============================================================================

print("=" * 90)
print("TEST 9: FIRST COUNTEREXAMPLE TO H_c | (q-c)")
print("=" * 90)
print()

for c in SHIFTS:

    found = False

    for n, p, q in ordered_states:

        C = q + 3
        H = gcd(C, n + c)

        if (q - c) % H != 0:

            print(
                f"c={c}:"
            )
            print(
                f"    n={n}"
            )
            print(
                f"    p={p}"
            )
            print(
                f"    q={q}"
            )
            print(
                f"    C=q+3={C}"
            )
            print(
                f"    H=gcd(C,n+c)={H}"
            )
            print(
                f"    q-c={q-c}"
            )
            print(
                f"    (q-c) mod H={(q-c)%H}"
            )

            found = True
            break

    if not found:
        print(
            f"c={c}: no counterexample "
            f"in q>p>{MIN_P}"
        )

print()


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("=" * 90)
print("FINAL STRUCTURAL SUMMARY")
print("=" * 90)
print()

print(
    """
The experiment tests the proposed Frame-A recovery channel:

    C = q + 3

and:

    H_c = gcd(C, n+c)

with candidate divisibility:

    H_c | (q-c)

for:

    c = 3, 81, 137.

The ordering restriction is:

    q > p > 137.

Therefore:

    sqrt(n) < q < n/137.

The candidate-generation mechanism is:

    factor(n+c)
        |
        v
    all divisors H of n+c
        |
        v
    assume H | (q-c)
        |
        v
    q = c + kH
        |
        v
    enforce sqrt(n) < q < n/137
        |
        v
    q prime
        |
        v
    q | n
        |
        v
    p=n/q > 137
        |
        v
    q > p

The critical distinction is that this experiment does NOT assume
the true H during candidate generation.

It independently tests whether the proposed divisibility law is
valid and whether the resulting congruence channels can recover q.

If:

    H_c | (q-c)

survives universally, then the three channels provide:

    q ≡ 3   (mod H_3)
    q ≡ 81  (mod H_81)
    q ≡ 137 (mod H_137)

and their intersection can be used as a progressively narrower
q-congruence class.

If the relation fails, the algebraic identity that is always true
for C=q+3 is:

    H_c = gcd(q+3, c-3p).

That failure would mean that H_c is constraining the p-side affine
quantity rather than directly the q-c difference.

The experiment therefore distinguishes these two possibilities
instead of assuming them equivalent.
"""
)

print()
print(f"BASELINE FAILURES={baseline_failures}")
print("=" * 90)
print("EXPERIMENT 694 FINISHED")
print("=" * 90)
