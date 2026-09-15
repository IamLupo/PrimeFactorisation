#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 693
DIVISOR-LATTICE -> C RECOVERY WITH EXPLICIT FACTOR-ORDERING
==========================================================================================

Goal
----
Experiment 692 showed that:

    gcd(C, n+c) = H_c

is tautological when H_c is computed from the guessed C.

This experiment reverses the direction.

For each shift c:

    M = n + c
    factor(M)
    D = all divisors of M

Every admissible H must satisfy:

    H | M

We then construct C candidates directly from the divisor lattice:

    C = H * k

subject to the frame-specific bounds.

For every constructed C we independently test:

    Frame A:
        q = C - 3
        p = n / q

    Frame B:
        p = C - 1
        q = n / p

and require:

    p, q are integers
    p and q are prime
    p <= q

The experiment measures whether the divisor lattice itself can
select the canonical C.

IMPORTANT
---------
No hidden p, q, C, or H is used to construct candidates.

The known factorization is used ONLY as a validation oracle.

Tests
-----
0. Baseline
1. Single-shift divisor-derived C candidates
2. Multi-shift UNION of divisor-derived C candidates
3. Multi-shift INTERSECTION of candidate C sets
4. Exact gcd consistency after candidate generation
5. Factor-ordering reduction
6. Candidate-count distribution
7. Compare H=1/2/etc. pruning
8. Representative examples
9. Controls proving the true C is never discarded

==========================================================================================
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import gcd, isqrt, lcm
from sympy import primerange, factorint
from sympy.ntheory import divisors, isprime


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 200

# Use the shifts discovered in the previous experiments.
SHIFTS = [3, 9, 81, 137]

# Number of examples to print.
MAX_EXAMPLES = 20

# Candidate H filters.
# The first is the complete divisor lattice.
H_FILTERS = {
    "all": lambda h: True,
    "H>1": lambda h: h > 1,
    "H>=4": lambda h: h >= 4,
    "odd(H)>1": lambda h: (h % 2 == 1 and h > 1),
}


# =============================================================================
# BASIC HELPERS
# =============================================================================

def v2(x: int) -> int:
    """Return v_2(x) for x != 0."""
    x = abs(x)
    if x == 0:
        return 10**9

    r = 0
    while (x & 1) == 0:
        x >>= 1
        r += 1
    return r


def all_semiprimes(primes: list[int]) -> list[tuple[int, int, int]]:
    """
    Return (n, p, q) for all p <= q.
    """
    out = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            out.append((p * q, p, q))

    return out


def factor_divisor_lattice(n: int) -> tuple[dict[int, int], list[int]]:
    """
    Return factorization and complete positive divisor lattice.
    """
    fac = factorint(n)
    ds = divisors(n)
    return fac, ds


# =============================================================================
# FRAME CANDIDATE GENERATION
# =============================================================================

def frame_a_candidates_from_H(n: int, H: int) -> list[int]:
    """
    Frame A:
        C = q + 3
        q = C - 3
        p = n / q

    Since p <= q:
        q >= sqrt(n)

    Since p >= 2:
        q <= n/2

    Thus:

        sqrt(n) + 3 <= C <= n/2 + 3

    We generate only multiples of H.
    """
    lo = isqrt(n)
    if lo * lo < n:
        lo += 1

    c_min = lo + 3
    c_max = n // 2 + 3

    first = ((c_min + H - 1) // H) * H

    return list(range(first, c_max + 1, H))


def frame_b_candidates_from_H(n: int, H: int) -> list[int]:
    """
    Frame B:
        C = p + 1
        p = C - 1
        q = n / p

    Since p <= q:
        p <= sqrt(n)

    Therefore:

        4 <= C <= sqrt(n) + 1

    Again only multiples of H are generated.
    """
    c_min = 4
    c_max = isqrt(n) + 1

    first = ((c_min + H - 1) // H) * H

    return list(range(first, c_max + 1, H))


def validate_frame_a(n: int, C: int) -> tuple[int, int] | None:
    """
    Validate:

        C = q + 3
        q = C - 3
        p = n/q

    Require p,q prime and p <= q.
    """
    q = C - 3

    if q <= 1:
        return None

    if n % q != 0:
        return None

    p = n // q

    if p > q:
        return None

    if not isprime(p) or not isprime(q):
        return None

    return p, q


def validate_frame_b(n: int, C: int) -> tuple[int, int] | None:
    """
    Validate:

        C = p + 1
        p = C - 1
        q = n/p

    Require p,q prime and p <= q.
    """
    p = C - 1

    if p <= 1:
        return None

    if n % p != 0:
        return None

    q = n // p

    if p > q:
        return None

    if not isprime(p) or not isprime(q):
        return None

    return p, q


def true_C_for_frame(frame: str, p: int, q: int) -> int:
    if frame == "A":
        return q + 3
    if frame == "B":
        return p + 1
    raise ValueError(frame)


def frame_candidates_from_H(
    n: int,
    H: int,
    frame: str,
) -> list[int]:

    if frame == "A":
        raw = frame_a_candidates_from_H(n, H)
    elif frame == "B":
        raw = frame_b_candidates_from_H(n, H)
    else:
        raise ValueError(frame)

    return raw


# =============================================================================
# INDEPENDENT DIVISOR-LATTICE SEARCH
# =============================================================================

def search_candidates_for_shift(
    n: int,
    c: int,
    frame: str,
    h_filter,
) -> tuple[set[int], dict[int, set[int]]]:
    """
    IMPORTANT:
    H is generated from the divisor lattice of M=n+c.

    We DO NOT compute H from C.

    For every lattice divisor H:
        C = H*k
    within the frame bounds.

    After generation, candidate C values are validated.
    """
    M = n + c

    if M <= 0:
        return set(), {}

    fac, lattice = factor_divisor_lattice(M)

    all_C = set()
    source_H = defaultdict(set)

    for H in lattice:
        if not h_filter(H):
            continue

        candidates = frame_candidates_from_H(n, H, frame)

        for C in candidates:
            # Independent frame algebra check.
            if frame == "A":
                pair = validate_frame_a(n, C)
            else:
                pair = validate_frame_b(n, C)

            if pair is not None:
                all_C.add(C)
                source_H[C].add(H)

    return all_C, source_H


# =============================================================================
# MULTI-SHIFT SEARCH
# =============================================================================

def search_all_shifts(
    n: int,
    frame: str,
    shifts: list[int],
    h_filter,
) -> dict[int, set[int]]:
    """
    Return candidate-C sets separately for every shift.
    """
    out = {}

    for c in shifts:
        candidates, _ = search_candidates_for_shift(
            n=n,
            c=c,
            frame=frame,
            h_filter=h_filter,
        )
        out[c] = candidates

    return out


def intersection_candidates(
    per_shift: dict[int, set[int]]
) -> set[int]:
    """
    Intersection of candidate C sets across shifts.
    """
    sets = list(per_shift.values())

    if not sets:
        return set()

    result = set(sets[0])

    for s in sets[1:]:
        result &= s

    return result


def union_candidates(
    per_shift: dict[int, set[int]]
) -> set[int]:
    result = set()

    for s in per_shift.values():
        result |= s

    return result


# =============================================================================
# EXACT GCD INFORMATION FOR VALIDATION ONLY
# =============================================================================

def exact_H_vector(C: int, n: int, shifts: list[int]) -> tuple[int, ...]:
    """
    This is calculated ONLY after a C candidate has been generated.
    """
    return tuple(gcd(C, n + c) for c in shifts)


def candidate_matches_shift_divisor_lattice(
    C: int,
    n: int,
    shifts: list[int],
) -> bool:
    """
    Every gcd(C,n+c) must actually be one of the divisors of n+c.

    This is mathematically automatic, but explicitly checking it
    documents the information channel.
    """
    for c in shifts:
        M = n + c
        H = gcd(C, M)

        if M % H != 0:
            return False

    return True


# =============================================================================
# TEST 0
# =============================================================================

print("=" * 90)
print("EXPERIMENT 693 START")
print("=" * 90)

print()
print(f"prime limit={PRIME_LIMIT}")
print(f"shifts={SHIFTS}")

primes = list(primerange(3, PRIME_LIMIT + 1))
states = all_semiprimes(primes)

print(f"odd primes={len(primes)}")
print(f"semiprimes={len(states)}")


# =============================================================================
# TEST 0: BASELINE
# =============================================================================

baseline_failures = 0

for n, p, q in states:
    if p * q != n:
        baseline_failures += 1

print()
print("=" * 90)
print("TEST 0: BASELINE")
print("=" * 90)
print(f"checked={len(states)}")
print(f"failures={baseline_failures}")


# =============================================================================
# TEST 1:
# TRUE C MUST BE GENERATED FROM THE DIVISOR LATTICE
# =============================================================================

print()
print("=" * 90)
print("TEST 1: TRUE C IS GENERATED WITHOUT USING TRUE H")
print("=" * 90)

for frame in ("A", "B"):
    checked = 0
    missing = 0
    examples = []

    for n, p, q in states:
        C_true = true_C_for_frame(frame, p, q)

        found = False

        # Try every shift, using only divisor lattices.
        for c in SHIFTS:
            candidates, _ = search_candidates_for_shift(
                n=n,
                c=c,
                frame=frame,
                h_filter=H_FILTERS["all"],
            )

            if C_true in candidates:
                found = True
                break

        checked += 1

        if not found:
            missing += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append((n, p, q, C_true))

    print()
    print(f"FRAME {frame}")
    print("-" * 90)
    print(f"checked={checked}")
    print(f"true C missing={missing}")
    print(f"true C generated={checked - missing}")

    if examples:
        print("missing examples:")
        for row in examples:
            print("   ", row)


# =============================================================================
# TEST 2:
# SINGLE SHIFT CANDIDATE COUNTS
# =============================================================================

print()
print("=" * 90)
print("TEST 2: SINGLE-SHIFT DIVISOR-LATTICE -> C")
print("=" * 90)

for frame in ("A", "B"):
    print()
    print(f"FRAME {frame}")
    print("-" * 90)

    for c in SHIFTS:
        counts = Counter()
        unique = 0
        ambiguous = 0
        missing = 0

        for n, p, q in states:
            C_true = true_C_for_frame(frame, p, q)

            candidates, _ = search_candidates_for_shift(
                n=n,
                c=c,
                frame=frame,
                h_filter=H_FILTERS["all"],
            )

            k = len(candidates)
            counts[k] += 1

            if C_true not in candidates:
                missing += 1
            elif k == 1:
                unique += 1
            else:
                ambiguous += 1

        total = len(states)

        print(
            f"c={c:4d} "
            f"unique={unique:4d} "
            f"ambiguous={ambiguous:4d} "
            f"missing={missing:4d} "
            f"mean_candidates="
            f"{sum(k*v for k,v in counts.items())/total:.4f}"
        )


# =============================================================================
# TEST 3:
# MULTI-SHIFT UNION / INTERSECTION
# =============================================================================

print()
print("=" * 90)
print("TEST 3: MULTI-SHIFT UNION / INTERSECTION")
print("=" * 90)

for frame in ("A", "B"):
    print()
    print(f"FRAME {frame}")
    print("-" * 90)

    for prefix_len in range(1, len(SHIFTS) + 1):
        active_shifts = SHIFTS[:prefix_len]

        union_total = 0
        intersection_total = 0

        union_unique = 0
        intersection_unique = 0

        union_missing = 0
        intersection_missing = 0

        for n, p, q in states:
            C_true = true_C_for_frame(frame, p, q)

            per_shift = search_all_shifts(
                n=n,
                frame=frame,
                shifts=active_shifts,
                h_filter=H_FILTERS["all"],
            )

            U = union_candidates(per_shift)
            I = intersection_candidates(per_shift)

            union_total += len(U)
            intersection_total += len(I)

            if len(U) == 1:
                union_unique += 1

            if len(I) == 1:
                intersection_unique += 1

            if C_true not in U:
                union_missing += 1

            if C_true not in I:
                intersection_missing += 1

        total = len(states)

        print(
            f"shifts={active_shifts}"
        )
        print(
            f"    union: "
            f"unique={union_unique:4d} "
            f"missing={union_missing:4d} "
            f"mean_candidates={union_total/total:.4f}"
        )
        print(
            f"    intersection: "
            f"unique={intersection_unique:4d} "
            f"missing={intersection_missing:4d} "
            f"mean_candidates={intersection_total/total:.4f}"
        )


# =============================================================================
# TEST 4:
# H-FILTER INFORMATION BOUNDARY
# =============================================================================

print()
print("=" * 90)
print("TEST 4: WHICH DIVISORS H ACTUALLY CARRY USEFUL INFORMATION?")
print("=" * 90)

for frame in ("A", "B"):
    print()
    print(f"FRAME {frame}")
    print("-" * 90)

    for filter_name, h_filter in H_FILTERS.items():
        total_candidates = 0
        unique = 0
        ambiguous = 0
        missing = 0

        for n, p, q in states:
            C_true = true_C_for_frame(frame, p, q)

            per_shift = search_all_shifts(
                n=n,
                frame=frame,
                shifts=SHIFTS,
                h_filter=h_filter,
            )

            candidates = intersection_candidates(per_shift)

            total_candidates += len(candidates)

            if C_true not in candidates:
                missing += 1
            elif len(candidates) == 1:
                unique += 1
            else:
                ambiguous += 1

        total = len(states)

        print(
            f"{filter_name:10s} "
            f"unique={unique:4d} "
            f"ambiguous={ambiguous:4d} "
            f"missing={missing:4d} "
            f"mean_candidates={total_candidates/total:.4f}"
        )


# =============================================================================
# TEST 5:
# FACTOR-ORDERING TEST
# =============================================================================

print()
print("=" * 90)
print("TEST 5: FACTOR-ORDERING REDUCTION")
print("=" * 90)

for frame in ("A", "B"):
    before = Counter()
    after = Counter()

    for n, p, q in states:
        C_true = true_C_for_frame(frame, p, q)

        # Generate candidates WITHOUT p <= q.
        raw = set()

        for c in SHIFTS:
            M = n + c
            _, lattice = factor_divisor_lattice(M)

            for H in lattice:
                for C in frame_candidates_from_H(n, H, frame):
                    if frame == "A":
                        q2 = C - 3
                        if q2 > 1 and n % q2 == 0:
                            p2 = n // q2
                            if isprime(p2) and isprime(q2):
                                raw.add(C)

                    else:
                        p2 = C - 1
                        if p2 > 1 and n % p2 == 0:
                            q2 = n // p2
                            if isprime(p2) and isprime(q2):
                                raw.add(C)

        raw_count = len(raw)

        ordered = set()

        for C in raw:
            if frame == "A":
                pair = validate_frame_a(n, C)
            else:
                pair = validate_frame_b(n, C)

            if pair is not None:
                ordered.add(C)

        ordered_count = len(ordered)

        before[raw_count] += 1
        after[ordered_count] += 1

    print()
    print(f"FRAME {frame}")
    print(f"raw candidates: total states={len(states)}")
    print(f"ordered candidates: total states={len(states)}")

    print("largest raw candidate multiplicities:")
    for k, v in sorted(before.items(), reverse=True)[:10]:
        print(f"    {k:5d} -> {v}")

    print("largest ordered candidate multiplicities:")
    for k, v in sorted(after.items(), reverse=True)[:10]:
        print(f"    {k:5d} -> {v}")


# =============================================================================
# TEST 6:
# REPRESENTATIVE STATES
# =============================================================================

print()
print("=" * 90)
print("TEST 6: REPRESENTATIVE STATES")
print("=" * 90)

interesting = []

for n, p, q in states:
    for frame in ("A", "B"):
        C_true = true_C_for_frame(frame, p, q)

        per_shift = search_all_shifts(
            n=n,
            frame=frame,
            shifts=SHIFTS,
            h_filter=H_FILTERS["all"],
        )

        I = intersection_candidates(per_shift)

        if len(I) <= 5:
            interesting.append(
                (
                    n,
                    p,
                    q,
                    frame,
                    C_true,
                    I,
                )
            )

interesting = interesting[:MAX_EXAMPLES]

for n, p, q, frame, C_true, candidates in interesting:
    print()
    print(
        f"frame={frame} n={n} p={p} q={q} "
        f"true_C={C_true}"
    )
    print(f"intersection candidates={sorted(candidates)}")

    for c in SHIFTS:
        M = n + c
        fac, lattice = factor_divisor_lattice(M)
        H_true = gcd(C_true, M)

        print(
            f"    c={c:4d} "
            f"M={M:6d} "
            f"factor(M)={fac} "
            f"H_true={H_true:6d} "
            f"H_is_divisor={H_true in lattice}"
        )


# =============================================================================
# TEST 7:
# TRUE-C CONTROL
# =============================================================================

print()
print("=" * 90)
print("TEST 7: TRUE-C CONTROL")
print("=" * 90)

control_failures = 0

for n, p, q in states:
    for frame in ("A", "B"):
        C_true = true_C_for_frame(frame, p, q)

        for c in SHIFTS:
            M = n + c
            fac = factorint(M)

            if M % gcd(C_true, M) != 0:
                control_failures += 1

print(f"checked={len(states) * 2 * len(SHIFTS)}")
print(f"control failures={control_failures}")


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print()
print("=" * 90)
print("FINAL STRUCTURAL SUMMARY")
print("=" * 90)

print(
    """
The experiment reverses the direction of the previous tests.

Previous approach:
    guess C
       |
       v
    compute H = gcd(C,n+c)

That makes H a consequence of C.

Current approach:
    factor(n+c)
       |
       v
    enumerate H | (n+c)
       |
       v
    construct C = H*k
       |
       v
    test frame equation
       |
       v
    test p,q prime
       |
       v
    enforce p <= q

The important distinction is:

    H is now generated BEFORE C is known.

The experiment measures whether this independent divisor lattice
is sufficient to recover C.

The critical quantities are:

    single-shift candidate count
    multi-shift union
    multi-shift intersection
    unique-C rate
    missing-C rate

A successful result would be:

    true C present
    AND
    candidate count approaches 1

A particularly important control is factor ordering:

    the apparent two-factor solutions may simply correspond to

        (p,q)

    versus

        (q,p)

which should collapse once p <= q is enforced.

No hidden p, q, C, or H is used for candidate generation.
The known factorization is used only for validation.
"""
)

print()
print(f"BASELINE FAILURES={baseline_failures}")
print(f"CONTROL FAILURES={control_failures}")

print("=" * 90)
print("EXPERIMENT 693 FINISHED")
print("=" * 90)
