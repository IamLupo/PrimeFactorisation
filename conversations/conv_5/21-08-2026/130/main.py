#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 691 START
==========================================================================================

OBSERVABLE MULTI-SHIFT DIVISOR-LATTICE -> C RECOVERY

Core idea
---------

For a fixed canonical C:

    H_c = gcd(C, n+c)

Therefore:

    H_c | (n+c)
    H_c | C

So for several shifts:

    L = lcm(H_c1, H_c2, ..., H_cr)

satisfies:

    L | C.

The previous experiment computed H_c directly from the TRUE C.
That made eventual recovery of C partly tautological.

This experiment removes that information.

We factor only:

    n+c

and generate ALL divisors of each n+c.

We then consider possible H values from those divisor lattices.

For every tuple:

    h_1 in Div(n+c_1)
    h_2 in Div(n+c_2)
    ...

we compute

    L = lcm(h_1,h_2,...).

Every genuine tuple must satisfy:

    L | C.

Therefore the candidate C values are generated from:

    C = L * k.

Then we test the frame equations:

FRAME A:
    C = q+3
    q = C-3

FRAME B:
    C = p+1
    p = C-1

The experimental harness knows p,q ONLY for validation.
The candidate-generation process itself does NOT use factorint(n).

Main questions
--------------

1. Does the divisor lattice contain the true H_c for every shift?
2. Does some divisor combination produce a strong L | C constraint?
3. How many C candidates remain after the LCM constraint?
4. How many remain after the elementary frame equation?
5. Does adding more shifts reduce the candidate set?
6. Can the true C become UNIQUE using only factorizations of n+c?

This is the first experiment that attempts to make the
multi-shift idea observable rather than using hidden H values.

==========================================================================================
"""

from __future__ import annotations

from math import gcd, lcm, isqrt
from collections import Counter

from sympy import primerange, factorint, isprime


# ----------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------

PRIME_LIMIT = 200

SHIFT_SETS = [
    [3],
    [3, 9],
    [3, 9, 81],
    [3, 9, 81, 137],
]

EXAMPLE_LIMIT = 20

# Safety bound for number of divisor combinations.
MAX_L_COMBINATIONS = 100000

# Maximum number of C candidates retained per state.
MAX_C_CANDIDATES = 100000


# ----------------------------------------------------------------------------------------
# BASIC HELPERS
# ----------------------------------------------------------------------------------------

def divisors_from_factorization(factors: dict[int, int]) -> list[int]:
    divisors = [1]

    for p, e in factors.items():
        old = list(divisors)
        divisors = []

        power = 1

        for _ in range(e + 1):
            for d in old:
                divisors.append(d * power)

            power *= p

    return sorted(divisors)


def all_divisors(n: int) -> list[int]:
    return divisors_from_factorization(factorint(n))


def canonical_C(frame: str, p: int, q: int) -> int:
    if frame == "A":
        return q + 3
    return p + 1


def canonical_H(frame: str, p: int, q: int, n: int, c: int) -> int:
    C = canonical_C(frame, p, q)
    return gcd(C, n + c)


def frame_candidate_valid(n: int, frame: str, C: int) -> bool:
    """
    Candidate test WITHOUT factorint(n).

    Frame A:
        q = C - 3
        p = n/q

    Frame B:
        p = C - 1
        q = n/p
    """

    if C <= 0:
        return False

    if frame == "A":
        q = C - 3

        if q < 3 or not isprime(q):
            return False

        if n % q != 0:
            return False

        p = n // q

        return p >= 3 and isprime(p)

    if frame == "B":
        p = C - 1

        if p < 3 or not isprime(p):
            return False

        if n % p != 0:
            return False

        q = n // p

        return q >= 3 and isprime(q)

    raise ValueError(frame)


# ----------------------------------------------------------------------------------------
# VALIDATION
# ----------------------------------------------------------------------------------------

def baseline(states):
    failures = 0

    for n, p, q in states:
        if n != p * q:
            failures += 1

    print("=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ----------------------------------------------------------------------------------------
# TEST 1
# ----------------------------------------------------------------------------------------

def test_divisor_lattice_contains_true_H(states, shifts):
    print("=" * 90)
    print(f"TEST 1: TRUE H IS IN EVERY DIVISOR LATTICE shifts={shifts}")
    print("=" * 90)

    failures = 0

    for frame in ("A", "B"):

        for n, p, q in states:

            for c in shifts:

                M = n + c
                divs = all_divisors(M)

                H = canonical_H(frame, p, q, n, c)

                if H not in divs:
                    failures += 1

                    print(
                        f"missing H "
                        f"frame={frame} n={n} c={c} "
                        f"H={H}"
                    )

                    if failures >= 20:
                        break

            if failures >= 20:
                break

        if failures >= 20:
            break

    print(f"failures={failures}")
    print()

    return failures


# ----------------------------------------------------------------------------------------
# LCM GENERATION
# ----------------------------------------------------------------------------------------

def generate_l_values(divisor_lists: list[list[int]]) -> set[int]:
    """
    Generate distinct LCM values from divisor choices.

        h_i in Div(n+c_i)

        L = lcm(h_1,...,h_r)

    We deduplicate aggressively.

    Prune values that are already larger than the maximum possible C
    for the given frame outside this function if desired.
    """

    current = {1}

    for divs in divisor_lists:

        next_values = set()

        for L in current:
            for d in divs:
                v = lcm(L, d)
                next_values.add(v)

                if len(next_values) > MAX_L_COMBINATIONS:
                    break

            if len(next_values) > MAX_L_COMBINATIONS:
                break

        current = next_values

        if len(current) > MAX_L_COMBINATIONS:
            break

    return current


# ----------------------------------------------------------------------------------------
# CANDIDATE GENERATION
# ----------------------------------------------------------------------------------------

def generate_C_candidates(
    n: int,
    frame: str,
    L_values: set[int],
) -> set[int]:

    candidates = set()

    sqrt_n = isqrt(n)

    if frame == "A":
        # q >= p and pq=n
        #
        # q >= sqrt(n)
        # q <= n/3 for p >= 3
        #
        # C=q+3
        min_C = sqrt_n + 3
        max_C = n // 3 + 3

    else:
        # p <= q and pq=n
        #
        # p <= sqrt(n)
        #
        # C=p+1
        min_C = 4
        max_C = sqrt_n + 1

    for L in L_values:

        if L <= 0:
            continue

        first_k = (min_C + L - 1) // L

        last_k = max_C // L

        for k in range(first_k, last_k + 1):

            C = L * k

            if C < min_C or C > max_C:
                continue

            candidates.add(C)

            if len(candidates) >= MAX_C_CANDIDATES:
                return candidates

    return candidates


# ----------------------------------------------------------------------------------------
# TEST 2
# ----------------------------------------------------------------------------------------

def test_observable_candidate_reduction(states, shifts):
    print("=" * 90)
    print(
        "TEST 2: OBSERVABLE LCM -> C CANDIDATE REDUCTION "
        f"shifts={shifts}"
    )
    print("=" * 90)

    stats = {
        "A": Counter(),
        "B": Counter(),
    }

    examples = {
        "A": [],
        "B": [],
    }

    for frame in ("A", "B"):

        for n, p, q in states:

            # Build divisor lattices ONLY from n+c.
            divisor_lists = [
                all_divisors(n + c)
                for c in shifts
            ]

            L_values = generate_l_values(divisor_lists)

            C_candidates = generate_C_candidates(
                n,
                frame,
                L_values,
            )

            true_C = canonical_C(frame, p, q)

            if true_C not in C_candidates:
                stats[frame]["TRUE_C_MISSING"] += 1

                if len(examples[frame]) < EXAMPLE_LIMIT:
                    examples[frame].append(
                        (
                            n,
                            p,
                            q,
                            true_C,
                            len(L_values),
                            len(C_candidates),
                            [],
                        )
                    )

                continue

            stats[frame][len(C_candidates)] += 1

            if len(examples[frame]) < EXAMPLE_LIMIT and len(C_candidates) <= 10:

                examples[frame].append(
                    (
                        n,
                        p,
                        q,
                        true_C,
                        len(L_values),
                        len(C_candidates),
                        sorted(C_candidates)[:20],
                    )
                )

        print(f"FRAME {frame}")
        print("-" * 90)

        total = sum(
            v for k, v in stats[frame].items()
            if k != "TRUE_C_MISSING"
        )

        missing = stats[frame]["TRUE_C_MISSING"]

        print(f"states={len(states)}")
        print(f"true C missing={missing}")

        if total:
            unique = stats[frame][1]
            ambiguous = total - unique

            print(f"unique C={unique}")
            print(f"ambiguous C={ambiguous}")
            print(
                f"unique ratio={unique/total:.6f}"
            )

        print("candidate-count distribution:")
        for key in sorted(
            (k for k in stats[frame] if isinstance(k, int))
        ):
            print(
                f"    {key:6d} -> {stats[frame][key]}"
            )

        if examples[frame]:
            print("examples:")

            for (
                n,
                p,
                q,
                true_C,
                lcount,
                ccount,
                candidates,
            ) in examples[frame]:

                print(
                    f"    n={n} p={p} q={q} "
                    f"true_C={true_C}"
                )
                print(
                    f"        L-values={lcount} "
                    f"C-candidates={ccount}"
                )
                print(
                    f"        candidates={candidates}"
                )

        print()


# ----------------------------------------------------------------------------------------
# TEST 3
# ----------------------------------------------------------------------------------------

def test_true_H_lcm_position(states, shifts):
    """
    This does use the hidden true H ONLY as a diagnostic.

    It determines where the true LCM sits among all observable L-values.

    This measures how much information is theoretically retained by
    the divisor lattice.
    """

    print("=" * 90)
    print(
        "TEST 3: TRUE LCM POSITION INSIDE OBSERVABLE LATTICE "
        f"shifts={shifts}"
    )
    print("=" * 90)

    for frame in ("A", "B"):

        missing = 0
        exact_L_equals_C = 0
        smaller_L = 0

        ratio_buckets = Counter()

        for n, p, q in states:

            divisor_lists = [
                all_divisors(n + c)
                for c in shifts
            ]

            L_values = generate_l_values(divisor_lists)

            true_L = 1

            for c in shifts:
                H = canonical_H(
                    frame,
                    p,
                    q,
                    n,
                    c,
                )

                true_L = lcm(true_L, H)

            C = canonical_C(frame, p, q)

            if true_L not in L_values:
                missing += 1
                continue

            if true_L == C:
                exact_L_equals_C += 1
            else:
                smaller_L += 1

            ratio = true_L / C

            if ratio == 1:
                ratio_buckets["1"] += 1
            elif ratio >= 0.5:
                ratio_buckets[">=1/2"] += 1
            elif ratio >= 0.25:
                ratio_buckets[">=1/4"] += 1
            else:
                ratio_buckets["<1/4"] += 1

        print(f"FRAME {frame}")
        print(f"true L missing={missing}")
        print(f"true L == C={exact_L_equals_C}")
        print(f"true L < C={smaller_L}")

        print("true L / C:")
        for k in ("1", ">=1/2", ">=1/4", "<1/4"):
            print(
                f"    {k:8s} -> {ratio_buckets[k]}"
            )

        print()


# ----------------------------------------------------------------------------------------
# TEST 4
# ----------------------------------------------------------------------------------------

def test_shift_progression(states):
    print("=" * 90)
    print("TEST 4: SHIFT PROGRESSION")
    print("=" * 90)

    for shifts in SHIFT_SETS:

        print(
            f"shifts={shifts}"
        )

        # Run only the candidate-reduction summary.
        for frame in ("A", "B"):

            unique = 0
            ambiguous = 0
            missing = 0
            total_candidates = 0

            for n, p, q in states:

                divisor_lists = [
                    all_divisors(n + c)
                    for c in shifts
                ]

                L_values = generate_l_values(divisor_lists)

                C_candidates = generate_C_candidates(
                    n,
                    frame,
                    L_values,
                )

                true_C = canonical_C(frame, p, q)

                if true_C not in C_candidates:
                    missing += 1
                    continue

                if len(C_candidates) == 1:
                    unique += 1
                else:
                    ambiguous += 1

                total_candidates += len(C_candidates)

            usable = unique + ambiguous

            mean_candidates = (
                total_candidates / usable
                if usable
                else 0
            )

            print(
                f"    frame={frame} "
                f"unique={unique} "
                f"ambiguous={ambiguous} "
                f"missing={missing} "
                f"mean_candidates={mean_candidates:.3f}"
            )

        print()


# ----------------------------------------------------------------------------------------
# TEST 5
# ----------------------------------------------------------------------------------------

def test_what_extra_constraint_is_needed(states):
    """
    For ambiguous candidate sets, inspect whether the correct C can be
    selected by checking additional observable quantities:

        gcd(C, n+c)

for an extra shift that was NOT used to generate the candidate.

The candidate C is accepted only if it reproduces an H value that is
present in the divisor lattice of the new n+c.

IMPORTANT:
this still does not tell us WHICH divisor is the true H.
It only tests whether an additional shift can eliminate candidates.
"""

    print("=" * 90)
    print("TEST 5: EXTRA-SHIFT COMPATIBILITY FILTER")
    print("=" * 90)

    base_shifts = [3, 9]
    extra_shifts = [81, 137, 201, 251]

    for frame in ("A", "B"):

        print(f"FRAME {frame}")

        total_ambiguous = 0

        for n, p, q in states:

            divisor_lists = [
                all_divisors(n + c)
                for c in base_shifts
            ]

            L_values = generate_L_values = generate_l_values(divisor_lists)

            candidates = generate_C_candidates(
                n,
                frame,
                L_values,
            )

            if len(candidates) <= 1:
                continue

            total_ambiguous += 1

        print(
            f"    base-shift ambiguous states={total_ambiguous}"
        )

        for extra_c in extra_shifts:

            remaining = 0
            eliminated = 0

            for n, p, q in states:

                divisor_lists = [
                    all_divisors(n + c)
                    for c in base_shifts
                ]

                L_values = generate_l_values(divisor_lists)

                candidates = generate_C_candidates(
                    n,
                    frame,
                    L_values,
                )

                if len(candidates) <= 1:
                    continue

                M = n + extra_c
                divs = set(all_divisors(M))

                filtered = []

                for C in candidates:
                    H = gcd(C, M)

                    if H in divs:
                        filtered.append(C)

                if len(filtered) > 1:
                    remaining += 1
                elif len(filtered) == 1:
                    eliminated += 1

            print(
                f"    extra c={extra_c:3d}: "
                f"remaining_ambiguous={remaining} "
                f"resolved={eliminated}"
            )

        print()


# ----------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------

def main():
    print("=" * 90)
    print("EXPERIMENT 691 START")
    print("=" * 90)

    print(f"prime limit={PRIME_LIMIT}")

    primes = list(primerange(3, PRIME_LIMIT + 1))

    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append((p * q, p, q))

    print(f"odd primes={len(primes)}")
    print(f"semiprimes={len(states)}")
    print()

    failures = 0

    failures += baseline(states)

    for shifts in SHIFT_SETS:

        failures += test_divisor_lattice_contains_true_H(
            states,
            shifts,
        )

        test_observable_candidate_reduction(
            states,
            shifts,
        )

        test_true_H_lcm_position(
            states,
            shifts,
        )

    test_shift_progression(states)

    test_what_extra_constraint_is_needed(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
The crucial correction to Experiment 690 is:

    DO NOT compute H_c from the hidden true C.

Instead:

    M_c = n+c
    factor(M_c)
    D_c = complete divisor lattice of M_c

and therefore:

    H_c belongs to D_c.

For several shifts:

    H_c1 | C
    H_c2 | C
    ...

so:

    L = lcm(H_c1,H_c2,...) | C.

The experiment therefore generates possible L values directly
from the observable divisor lattices.

Every resulting candidate has:

    C = L * k.

The frame equations then restrict the allowed C:

    Frame A:
        q = C-3
        pq=n

    Frame B:
        p = C-1
        pq=n

The key information boundary is now measurable:

    divisor lattices
          |
          v
    possible H_c
          |
          v
    possible L = lcm(H_c)
          |
          v
    possible C = L*k
          |
          v
    canonical frame candidates

A strong result would be:

    very few C candidates remain,
    preferably one,

using ONLY factorizations of n+c.

A weak result would be:

    true C is present,
    but many C candidates survive.

That would mean the shifted-factorization channel contains
the correct information but does not yet select it uniquely.

TOTAL BASELINE FAILURES=%d
==========================================================================================
EXPERIMENT 691 FINISHED
==========================================================================================
"""
        % failures
    )


if __name__ == "__main__":
    main()
