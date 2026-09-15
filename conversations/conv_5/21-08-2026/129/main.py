#!/usr/bin/env python3
"""
==========================================================================================
EXPERIMENT 690
==========================================================================================

MULTI-SHIFT H-LATTICE ACCUMULATION

Idea under test
---------------

For a FIXED canonical C, define

    H_c = gcd(C, n + c)

for several odd shifts c.

Since every H_c divides C,

    lcm(H_c1, H_c2, ...) | C.

Therefore multiple different H observations may progressively reveal
different prime-power components of the same C.

This experiment does NOT change C when c changes.

Canonical constructions under test:

    FRAME A:
        C_A = q + 3

    FRAME B:
        C_B = p + 1

For each fixed C, we test:

    H_c = gcd(C, n+c)

and accumulate

    L_r = lcm(H_c1, ..., H_cr)

Questions
---------

1. Does L_r frequently converge exactly to C?
2. How many shifts are typically required?
3. Do the known shifts c=3, 9, 81, 137 help?
4. Which prime-power components of C remain hidden?
5. If L < C, is C/L constrained to a small set?
6. Can a candidate C be recovered from:
       n
       frame
       L = lcm(H_c)
   without knowing the original factorization?

Important:
----------

The factorization p*q is used ONLY by the validation harness to know
the true C and H_c. The final candidate search does not use factorint(n).

==========================================================================================
"""

from __future__ import annotations

from math import gcd, lcm
from collections import defaultdict, Counter
from sympy import primerange, factorint, isprime


# ----------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------

PRIME_LIMIT = 200

# Small set containing the shifts already suggested by the research.
CORE_SHIFTS = [3, 9, 81, 137]

# Broad odd-shift family.
# All are odd, therefore n+c is even for odd n.
BROAD_SHIFTS = list(range(1, 202, 2))

# Candidate search limit:
# We test multiples of L in a finite interval around plausible C values.
#
# For Frame A, C=q+3 can be larger than sqrt(n).
# For Frame B, C=p+1 is at most sqrt(n)+1.
#
# This search is deliberately separated from the hidden-factor validation.
SEARCH_MULTIPLIER_LIMIT = 1000

EXAMPLE_LIMIT = 20


# ----------------------------------------------------------------------------------------
# BASIC HELPERS
# ----------------------------------------------------------------------------------------

def v2(x: int) -> int:
    if x == 0:
        return 10**9
    x = abs(x)
    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    x = abs(x)
    return x >> v2(x)


def canonical_data(p: int, q: int):
    """
    Return both mathematical affine frames.

    A:
        C = q+3

    B:
        C = p+1
    """
    return {
        "A": {
            "C": q + 3,
            "factor": q,
            "anchor": 3,
        },
        "B": {
            "C": p + 1,
            "factor": p,
            "anchor": -1,
        },
    }


def candidate_C_valid(n: int, frame: str, C: int) -> bool:
    """
    Check whether C could actually be the canonical C for this n/frame.

    Frame A:
        q = C-3
        p = n/(C-3)

    Frame B:
        p = C-1
        q = n/(C-1)

    We do NOT factor n here.
    We merely test the arithmetic consequences of a proposed C.
    """
    if C <= 0:
        return False

    if frame == "A":
        q = C - 3
        if q <= 1 or not isprime(q):
            return False
        if n % q != 0:
            return False
        p = n // q
        return p == 1 or isprime(p)

    if frame == "B":
        p = C - 1
        if p <= 1 or not isprime(p):
            return False
        if n % p != 0:
            return False
        q = n // p
        return q == 1 or isprime(q)

    raise ValueError(frame)


def collect_divisors(n: int) -> list[int]:
    """
    Generate all divisors of n using SymPy factorint.
    """
    fac = factorint(n)
    divisors = [1]

    for p, e in fac.items():
        old = list(divisors)
        divisors = []

        power = 1
        for _ in range(e + 1):
            divisors.extend(d * power for d in old)
            power *= p

    return sorted(divisors)


# ----------------------------------------------------------------------------------------
# TEST 0
# ----------------------------------------------------------------------------------------

def test_baseline(states):
    failures = 0

    for n, p, q in states:
        assert n == p * q

        for frame, data in canonical_data(p, q).items():
            C = data["C"]

            if frame == "A":
                expected = q + 3
            else:
                expected = p + 1

            if C != expected:
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

def test_multi_shift_lcm(states):
    print("=" * 90)
    print("TEST 1: MULTI-SHIFT LCM GROWTH")
    print("=" * 90)

    summary = {}
    examples = []

    for frame in ("A", "B"):
        print(f"FRAME {frame}")
        print("-" * 90)

        exact_by_r = Counter()
        ratio_buckets = Counter()
        shift_needed = []

        total = 0

        for n, p, q in states:
            C = canonical_data(p, q)[frame]["C"]

            L = 1
            exact_at = None

            for r, c in enumerate(BROAD_SHIFTS, start=1):
                H = gcd(C, n + c)
                L = lcm(L, H)

                if L == C and exact_at is None:
                    exact_at = r

            total += 1

            if exact_at is not None:
                exact_by_r[exact_at] += 1
                shift_needed.append(exact_at)
            else:
                exact_by_r[0] += 1

            ratio = L / C

            if ratio == 1:
                ratio_buckets["exact"] += 1
            elif ratio >= 0.75:
                ratio_buckets[">=3/4"] += 1
            elif ratio >= 0.5:
                ratio_buckets[">=1/2"] += 1
            elif ratio >= 0.25:
                ratio_buckets[">=1/4"] += 1
            else:
                ratio_buckets["<1/4"] += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    (
                        frame,
                        n,
                        p,
                        q,
                        C,
                        exact_at,
                        L,
                    )
                )

        exact_total = sum(v for k, v in exact_by_r.items() if k != 0)

        print(f"states={total}")
        print(f"exactly recovered by full shift set={exact_total}/{total}")

        if shift_needed:
            print(f"maximum shifts required={max(shift_needed)}")
            print(f"mean shifts required={sum(shift_needed)/len(shift_needed):.4f}")

        print("exact-at-shift distribution:")
        for k in sorted(exact_by_r):
            if k == 0:
                print(f"    never -> {exact_by_r[k]}")
            else:
                print(f"    {k:3d} -> {exact_by_r[k]}")

        print("final L/C distribution:")
        for key in ("exact", ">=3/4", ">=1/2", ">=1/4", "<1/4"):
            print(f"    {key:8s} -> {ratio_buckets[key]}")

        print()

        summary[frame] = (exact_total, total)

    print("EXAMPLES")
    print("-" * 90)

    shown = 0
    for frame, n, p, q, C, exact_at, L in examples:
        if shown >= EXAMPLE_LIMIT:
            break

        print(
            f"frame={frame} n={n} p={p} q={q} "
            f"C={C} exact_shift_index={exact_at} final_L={L}"
        )

        for c in CORE_SHIFTS:
            H = gcd(C, n + c)
            print(
                f"    c={c:3d} "
                f"H={H:6d} "
                f"v2(H)={v2(H):2d}"
            )

        shown += 1

    print()

    return summary


# ----------------------------------------------------------------------------------------
# TEST 2
# ----------------------------------------------------------------------------------------

def test_core_shift_set(states):
    print("=" * 90)
    print("TEST 2: CORE SHIFTS [3,9,81,137]")
    print("=" * 90)

    for frame in ("A", "B"):
        exact = 0
        total = 0
        l_values = []
        missing_prime_power_examples = []

        for n, p, q in states:
            C = canonical_data(p, q)[frame]["C"]

            L = 1

            for c in CORE_SHIFTS:
                H = gcd(C, n + c)
                L = lcm(L, H)

            total += 1

            if L == C:
                exact += 1
            else:
                # Determine which prime powers of C are still missing.
                missing = []

                for prime, exponent in factorint(C).items():
                    required = prime ** exponent

                    if L % required != 0:
                        missing.append(required)

                if len(missing_prime_power_examples) < 15:
                    missing_prime_power_examples.append(
                        (n, C, L, missing)
                    )

            l_values.append(L)

        print(
            f"frame={frame} "
            f"exact={exact}/{total} "
            f"ratio={exact/total:.6f}"
        )

        if missing_prime_power_examples:
            print("first unresolved prime-power examples:")
            for n, C, L, missing in missing_prime_power_examples:
                print(
                    f"    n={n} C={C} L={L} "
                    f"missing={missing}"
                )

        print()

# ----------------------------------------------------------------------------------------
# TEST 3
# ----------------------------------------------------------------------------------------

def test_incremental_shift_growth(states):
    print("=" * 90)
    print("TEST 3: INCREMENTAL SHIFT ACCUMULATION")
    print("=" * 90)

    for frame in ("A", "B"):
        print(f"FRAME {frame}")
        print("-" * 90)

        L = {0: 1}

        for c in BROAD_SHIFTS:
            exact = 0
            total = 0

            for n, p, q in states:
                C = canonical_data(p, q)[frame]["C"]

                current = 1

                for previous_c in BROAD_SHIFTS:
                    if previous_c > c:
                        break

                    current = lcm(current, gcd(C, n + previous_c))

                total += 1

                if current == C:
                    exact += 1

            L[c] = exact

            print(
                f"    through c={c:3d}: "
                f"exact={exact:6d}/{total}"
            )

        print()

# ----------------------------------------------------------------------------------------
# TEST 4
# ----------------------------------------------------------------------------------------

def test_candidate_search(states):
    """
    Given L = lcm(H_c), find all C values that:

        1. are multiples of L
        2. satisfy the frame's canonical factor equation
        3. lie in a finite search range.

    This is NOT used as a proof of a factorization algorithm.
    It measures whether the accumulated H information makes C
    unique once the elementary frame equation is enforced.
    """

    print("=" * 90)
    print("TEST 4: LCM -> CANONICAL C CANDIDATE REDUCTION")
    print("=" * 90)

    for frame in ("A", "B"):
        total = 0
        none = 0
        unique = 0
        ambiguous = 0

        candidate_count_distribution = Counter()
        examples = []

        for n, p, q in states:
            C_true = canonical_data(p, q)[frame]["C"]

            # Use the four core shifts as the first accumulated constraint.
            L = 1
            for c in CORE_SHIFTS:
                L = lcm(L, gcd(C_true, n + c))

            if L <= 0:
                continue

            total += 1

            # Search finite multiples of L.
            #
            # Include a broad range because the two frames have different
            # natural sizes for C.
            max_C = n + 10

            candidates = []

            for C in range(L, max_C + 1, L):
                if candidate_C_valid(n, frame, C):
                    candidates.append(C)

                    # Once we know the true one and candidate list is huge,
                    # don't artificially truncate the classification.
                    if len(candidates) > 1000:
                        break

            if C_true not in candidates:
                # It should always be present because C_true is a multiple
                # of every H_c, hence of L.
                print(
                    f"ERROR: true C missing "
                    f"frame={frame} n={n} C={C_true} L={L}"
                )
                none += 1
                continue

            candidate_count_distribution[len(candidates)] += 1

            if len(candidates) == 0:
                none += 1
            elif len(candidates) == 1:
                unique += 1
            else:
                ambiguous += 1

            if len(examples) < EXAMPLE_LIMIT and len(candidates) > 1:
                examples.append(
                    (n, p, q, C_true, L, candidates[:20])
                )

        print(f"frame={frame}")
        print(f"total={total}")
        print(f"none={none}")
        print(f"unique={unique}")
        print(f"ambiguous={ambiguous}")

        print("candidate-count distribution:")
        for count in sorted(candidate_count_distribution):
            print(
                f"    {count:4d} -> "
                f"{candidate_count_distribution[count]}"
            )

        if examples:
            print("ambiguous examples:")
            for n, p, q, C_true, L, candidates in examples:
                print(
                    f"    n={n} "
                    f"true_C={C_true} "
                    f"L={L}"
                )
                print(
                    f"        candidates={candidates}"
                )

        print()


# ----------------------------------------------------------------------------------------
# TEST 5
# ----------------------------------------------------------------------------------------

def test_lcm_prime_power_recovery(states):
    """
    For each C, inspect every prime-power component.

    A component is visible if it divides at least one H_c.
    We report how many shifts are needed before each component
    appears in the LCM.
    """

    print("=" * 90)
    print("TEST 5: PRIME-POWER VISIBILITY")
    print("=" * 90)

    for frame in ("A", "B"):
        visibility = Counter()
        unresolved = []

        for n, p, q in states:
            C = canonical_data(p, q)[frame]["C"]
            components = [
                prime ** exponent
                for prime, exponent in factorint(C).items()
            ]

            L = 1
            found_at = {component: None for component in components}

            for index, c in enumerate(BROAD_SHIFTS, start=1):
                H = gcd(C, n + c)
                L = lcm(L, H)

                for component in components:
                    if found_at[component] is None and L % component == 0:
                        found_at[component] = index

            for component, index in found_at.items():
                if index is None:
                    visibility["never"] += 1
                    if len(unresolved) < 20:
                        unresolved.append(
                            (n, C, component)
                        )
                else:
                    visibility[index] += 1

        print(f"FRAME {frame}")
        print("component recovery index distribution:")

        # Only print the first useful part.
        for index in sorted(visibility):
            print(f"    {index} -> {visibility[index]}")

        if unresolved:
            print("unresolved components:")
            for n, C, component in unresolved:
                print(
                    f"    n={n} C={C} "
                    f"missing_component={component}"
                )

        print()


# ----------------------------------------------------------------------------------------
# TEST 6
# ----------------------------------------------------------------------------------------

def test_known_shift_comparison(states):
    """
    Compare:

        only c=3
        only c=9
        c in [3,9]
        c in [3,9,81]
        c in [3,9,81,137]

    This directly tests whether multiple H values are useful
    compared with a single H.
    """

    print("=" * 90)
    print("TEST 6: SHIFT-SET COMPARISON")
    print("=" * 90)

    shift_sets = [
        [3],
        [9],
        [3, 9],
        [3, 9, 81],
        [3, 9, 81, 137],
    ]

    for frame in ("A", "B"):
        print(f"FRAME {frame}")
        print("-" * 90)

        for shifts in shift_sets:
            exact = 0

            for n, p, q in states:
                C = canonical_data(p, q)[frame]["C"]

                L = 1
                for c in shifts:
                    L = lcm(L, gcd(C, n + c))

                if L == C:
                    exact += 1

            print(
                f"    shifts={str(shifts):20s} "
                f"exact={exact:6d}/{len(states)} "
                f"ratio={exact/len(states):.6f}"
            )

        print()


# ----------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------

def main():
    print("=" * 90)
    print("EXPERIMENT 690 START")
    print("=" * 90)
    print(f"prime limit={PRIME_LIMIT}")

    primes = list(primerange(3, PRIME_LIMIT + 1))
    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            # Keep this experiment within the same class of positive odd semiprimes.
            states.append((n, p, q))

    print(f"odd primes={len(primes)}")
    print(f"semiprimes={len(states)}")
    print()

    failures = 0

    failures += test_baseline(states)

    test_multi_shift_lcm(states)
    test_core_shift_set(states)
    test_incremental_shift_growth(states)
    test_candidate_search(states)
    test_lcm_prime_power_recovery(states)
    test_known_shift_comparison(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
For a FIXED canonical C:

    H_c = gcd(C, n+c)

for every tested odd shift c.

Therefore:

    H_c | C

and hence:

    lcm(H_c1, H_c2, ..., H_cr) | C.

The experiment tests whether different shifts expose different
prime-power components of the SAME C.

This is fundamentally different from changing C with the shift.

The main quantities are:

    H_1, H_2, ..., H_r

    L_r = lcm(H_1,...,H_r)

and the target relation:

    L_r = C.

If L_r < C, then:

    C = L_r * K

for some integer K.

The important questions are:

    1. How often does L_r = C?
    2. How many shifts are required?
    3. Which prime powers of C remain hidden?
    4. How much does [3,9,81,137] improve over one shift?
    5. Once L_r is known, how many canonical C candidates remain?

A positive result would mean that multiple shifted gcd channels
progressively reveal independent pieces of the SAME canonical C.

That would be much stronger evidence for the multi-shift recovery
mechanism than Experiment 689, because Experiment 689 changed C
when it changed the anchor a.
"""
    )

    print(f"TOTAL BASELINE FAILURES={failures}")
    print("=" * 90)
    print("EXPERIMENT 690 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
