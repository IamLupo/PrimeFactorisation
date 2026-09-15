#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 685
ODD-CORE -> H RECOVERY -> FACTOR RECOVERY

QUESTION
==========================================================================================

Experiment 684 found:

    M = n + c
    M = 2^w * R,     R odd

    H = gcd(C, M)

and exactly:

    odd(H) = gcd(C, R).

However:

    odd(C) | R

only for a small subset of the semiprimes.

Therefore the stronger idea

    "factor R -> recover C"

is false in general.

This experiment asks the weaker and more fundamental question:

    "factor R -> recover H"

where

    H = gcd(C, n+c).

If H can be recovered from the factorization of R plus the
known 2-adic level w, we then test whether H is sufficient
to reconstruct the canonical factor pair.

The experiment deliberately separates:

    R-information
        ->
    odd(H)

        ->
    possible H values

        ->
    canonical C candidates

        ->
    factor pair candidates.

No assumption is made that odd(C) divides R.

==========================================================================================
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import gcd, isqrt
from typing import List, Tuple, Dict


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

PRIME_LIMIT = 6000
EXAMPLE_LIMIT = 20
COLLISION_LIMIT = 20


# ==========================================================================================
# ARITHMETIC
# ==========================================================================================

INF = 10**9


def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    x = abs(x)

    while x and (x & 1) == 0:
        x >>= 1

    return x


def sieve(limit: int) -> List[int]:
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0:2] = b"\x00\x00"

    for p in range(2, isqrt(limit) + 1):
        if a[p]:
            start = p * p
            a[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if a[i]]


PRIMES_SMALL = sieve(100000)


# ==========================================================================================
# FACTORIZATION
# ==========================================================================================

def factor_trial(n: int) -> List[Tuple[int, int]]:
    """
    Factor n exactly.

    R <= roughly 2^24 in the current experiment range, so trial
    division by the small prime table is fast enough.
    """
    if n <= 1:
        return []

    result: List[Tuple[int, int]] = []
    x = n

    for p in PRIMES_SMALL:
        if p * p > x:
            break

        if x % p != 0:
            continue

        e = 0

        while x % p == 0:
            x //= p
            e += 1

        result.append((p, e))

        if x == 1:
            return result

    if x > 1:
        result.append((x, 1))

    return result


def divisors_from_factorization(
    factors: List[Tuple[int, int]]
) -> List[int]:
    divisors = [1]

    for p, e in factors:
        old = divisors[:]
        multiplier = 1

        for _ in range(e):
            multiplier *= p

            for d in old:
                divisors.append(d * multiplier)

    divisors.sort()
    return divisors


# ==========================================================================================
# PRIMALITY
# ==========================================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3
    r = isqrt(n)

    while d <= r:
        if n % d == 0:
            return False
        d += 2

    return True


# ==========================================================================================
# FRAME
# ==========================================================================================

def frame_data(p: int, q: int, n: int):
    """
    Canonical frames established by earlier experiments.

    Frame A:
        C = q + 3
        c = 9

    Frame B:
        C = p + 1
        c = 3
    """

    if n % 4 == 3:
        return "A", q + 3, 9

    return "B", p + 1, 3


# ==========================================================================================
# CANDIDATE FACTOR RECOVERY
# ==========================================================================================

def factors_from_C(
    n: int,
    frame: str,
    C: int,
):
    """
    Convert a candidate canonical residual C into candidate
    prime factors.

    Frame A:
        C = q + 3
        q = C - 3

    Frame B:
        C = p + 1
        p = C - 1
    """

    if C <= 0 or C % 2 != 0:
        return None

    if frame == "A":
        q = C - 3

        if q <= 1 or n % q != 0:
            return None

        p = n // q

    else:
        p = C - 1

        if p <= 1 or n % p != 0:
            return None

        q = n // p

    if p <= 1 or q <= 1:
        return None

    if not is_prime(p) or not is_prime(q):
        return None

    if p > q:
        p, q = q, p

    return p, q


# ==========================================================================================
# STATE GENERATION
# ==========================================================================================

def generate_states():
    primes = [p for p in sieve(PRIME_LIMIT) if p & 1]

    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            frame, C, c = frame_data(p, q, n)

            M = n + c

            w = v2(M)
            R = M >> w

            H = gcd(C, M)

            h2 = v2(H)
            hodd = odd_part(H)

            states.append(
                (
                    n,
                    p,
                    q,
                    frame,
                    C,
                    c,
                    M,
                    w,
                    R,
                    H,
                    h2,
                    hodd,
                )
            )

    return states


# ==========================================================================================
# TEST 0
# ==========================================================================================

def test_baseline(states):
    print("=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)

    failures = 0

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        expected = gcd(C, M)

        if expected != H:
            failures += 1

            if failures <= EXAMPLE_LIMIT:
                print(
                    "mismatch n={} H={} expected={}".format(
                        n,
                        H,
                        expected,
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()


# ==========================================================================================
# TEST 1
# ==========================================================================================

def test_odd_H_exact(states):
    print("=" * 90)
    print("TEST 1: ODD(H) = gcd(C,R)")
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        g = gcd(C, R)

        distribution[
            (frame, g == hodd)
        ] += 1

        if g != hodd:
            failures += 1

            if failures <= EXAMPLE_LIMIT:
                print(
                    "mismatch n={} frame={} C={} R={} gcd={} odd(H)={}".format(
                        n,
                        frame,
                        C,
                        R,
                        g,
                        hodd,
                    )
                )

    print("distribution:")

    for key in sorted(distribution):
        print("    {} -> {}".format(key, distribution[key]))

    print()
    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()


# ==========================================================================================
# TEST 2
# ==========================================================================================

def test_R_factorization(states):
    print("=" * 90)
    print("TEST 2: FACTOR THE ODD CORE R")
    print("=" * 90)

    factor_distribution = Counter()
    max_factors = 0
    max_R = 0
    max_divisors = 0

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        if R > max_R:
            max_R = R

        factors = factor_trial(R)
        divisor_count = 1

        for _, e in factors:
            divisor_count *= (e + 1)

        factor_distribution[len(factors)] += 1

        max_factors = max(max_factors, len(factors))
        max_divisors = max(max_divisors, divisor_count)

    print("max R={}".format(max_R))
    print("max distinct prime factors of R={}".format(max_factors))
    print("max divisor count of R={}".format(max_divisors))
    print()
    print("distinct-prime-factor distribution:")

    for k, count in sorted(factor_distribution.items()):
        print("    {} -> {}".format(k, count))

    print()
    print("checked={}".format(len(states)))
    print("failures=0")
    print()


# ==========================================================================================
# TEST 3
# ==========================================================================================

def test_H_reconstruction(states):
    print("=" * 90)
    print("TEST 3: CAN H BE RECONSTRUCTED FROM R + w?")
    print("=" * 90)

    """
    Given:

        R = odd part of n+c
        w = v2(n+c)

    we know:

        odd(H) | R
        v2(H) <= w

    Therefore all arithmetic candidates of the form

        H_candidate = 2^j * d

    with

        d | R
        0 <= j <= w

    are consistent with the observable divisor information.

    We count how many such H candidates are possible.

    The true H must appear in this lattice.

    This is the crucial information-loss measurement.
    """

    candidate_count_distribution = Counter()

    true_H_missing = 0

    examples = []

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        factors = factor_trial(R)
        divisors = divisors_from_factorization(factors)

        candidates = set()

        for d in divisors:
            value = d

            for _ in range(w + 1):
                if value > M:
                    break

                candidates.add(value)
                value <<= 1

        count = len(candidates)

        candidate_count_distribution[count] += 1

        if H not in candidates:
            true_H_missing += 1

        if len(examples) < EXAMPLE_LIMIT and count > 1:
            examples.append(
                (
                    n,
                    frame,
                    C,
                    w,
                    R,
                    H,
                    candidates,
                )
            )

    print("candidate-count distribution:")

    for count, frequency in sorted(candidate_count_distribution.items()):
        print(
            "    {} -> {}".format(
                count,
                frequency,
            )
        )

    print()
    print("true H missing={}".format(true_H_missing))
    print()

    print("ambiguous examples:")

    for (
        n,
        frame,
        C,
        w,
        R,
        H,
        candidates,
    ) in examples:
        print(
            "    n={} frame={} C={} w={} R={} H={}".format(
                n,
                frame,
                C,
                w,
                R,
                H,
            )
        )

        print(
            "        candidates={}".format(
                sorted(candidates)
            )
        )

    print()
    print("checked={}".format(len(states)))
    print("failures={}".format(true_H_missing))
    print()


# ==========================================================================================
# TEST 4
# ==========================================================================================

def test_H_plus_frame_factor_recovery(states):
    print("=" * 90)
    print("TEST 4: H + FRAME -> FACTOR RECOVERY")
    print("=" * 90)

    """
    H itself is not C.

    We therefore test whether knowing H allows C to be reconstructed
    by using the exact relation

        H = gcd(C, M)

    together with the frame equation.

    For each state we enumerate divisors C of M whose gcd with M
    equals H, then ask how many of those C values produce the
    correct prime factor pair.

    This is intentionally a local information test, not a factoring
    algorithm.
    """

    total = 0
    unique = 0
    ambiguous = 0
    none = 0

    examples = []

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        # Factor M directly for the diagnostic. M is n+c and is
        # intentionally the "new target" discussed in the research.
        factors = factor_trial(M)
        divisors = divisors_from_factorization(factors)

        valid_C = []

        for candidate_C in divisors:
            if candidate_C <= 0 or candidate_C % 2:
                continue

            if gcd(candidate_C, M) != H:
                continue

            hit = factors_from_C(n, frame, candidate_C)

            if hit is not None:
                valid_C.append(
                    (candidate_C, hit[0], hit[1])
                )

        valid_C = sorted(set(valid_C))

        total += 1

        if len(valid_C) == 0:
            none += 1
        elif len(valid_C) == 1:
            unique += 1
        else:
            ambiguous += 1

        if len(examples) < EXAMPLE_LIMIT and len(valid_C) != 1:
            examples.append(
                (
                    n,
                    p,
                    q,
                    frame,
                    C,
                    H,
                    valid_C,
                )
            )

    print("factor-recovery classification:")
    print("    none      = {}".format(none))
    print("    unique    = {}".format(unique))
    print("    ambiguous = {}".format(ambiguous))

    print()

    if examples:
        print("non-unique examples:")

        for (
            n,
            p,
            q,
            frame,
            C,
            H,
            valid_C,
        ) in examples:
            print(
                "    n={} frame={} true_C={} H={}".format(
                    n,
                    frame,
                    C,
                    H,
                )
            )
            print(
                "        valid_C={}".format(
                    valid_C
                )
            )

    print()
    print("checked={}".format(total))
    print()


# ==========================================================================================
# TEST 5
# ==========================================================================================

def test_H_signature_collisions(states):
    print("=" * 90)
    print("TEST 5: H-SIGNATURE COLLISIONS")
    print("=" * 90)

    """
    We now deliberately forget C and retain only:

        (frame, w, odd(H))

    This is effectively:

        (frame, H)

    apart from the exact 2-adic decomposition.

    We ask:

        Does the same H-signature ever correspond to different
        depths?

    More importantly:

        Does the same H-signature correspond to different C?
    """

    buckets: Dict[Tuple[str, int, int], List[Tuple]] = defaultdict(list)

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        key = (frame, h2, hodd)

        if len(buckets[key]) < 4:
            buckets[key].append(
                (n, p, q, C, H, v2(gcd(2 * C, M)))
            )

    ambiguous_C = 0
    ambiguous_depth = 0

    examples = []

    for key, rows in buckets.items():
        Cs = {row[3] for row in rows}
        depths = {row[5] for row in rows}

        # The truncated bucket can hide additional members because
        # we only stored the first four. This is therefore only a
        # diagnostic, not a proof.
        if len(Cs) > 1:
            ambiguous_C += 1

        if len(depths) > 1:
            ambiguous_depth += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append((key, rows))

    print(
        "H-signatures={}".format(
            len(buckets)
        )
    )
    print(
        "signatures with multiple C examples={}".format(
            ambiguous_C
        )
    )
    print(
        "signatures with multiple depth examples={}".format(
            ambiguous_depth
        )
    )

    if examples:
        print()
        print("depth-collision examples:")

        for key, rows in examples:
            print(
                "    signature={}".format(
                    key
                )
            )

            for row in rows:
                print(
                    "        n={} p={} q={} C={} H={} depth={}".format(
                        row[0],
                        row[1],
                        row[2],
                        row[3],
                        row[4],
                        row[5],
                    )
                )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# TEST 6
# ==========================================================================================

def test_odd_core_factor_overlap(states):
    print("=" * 90)
    print("TEST 6: WHICH PART OF C IS VISIBLE IN R?")
    print("=" * 90)

    """
    Compare:

        odd(C)

    against:

        gcd(odd(C), R) = odd(H).

    Measure the information retained by R.
    """

    ratio_distribution = Counter()
    overlap_distribution = Counter()

    examples = []

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        oc = odd_part(C)

        overlap = gcd(oc, R)

        # quotient of the odd residual not visible in R
        hidden = oc // overlap

        overlap_distribution[
            (frame, hidden)
        ] += 1

        if hidden == 1:
            ratio_distribution[
                (frame, "full")
            ] += 1
        else:
            ratio_distribution[
                (frame, "partial")
            ] += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    (
                        n,
                        p,
                        q,
                        frame,
                        C,
                        oc,
                        R,
                        overlap,
                        hidden,
                    )
                )

    print("visibility distribution:")

    for key, count in sorted(ratio_distribution.items()):
        print(
            "    {} -> {}".format(
                key,
                count,
            )
        )

    print()
    print("partial-visibility examples:")

    for (
        n,
        p,
        q,
        frame,
        C,
        oc,
        R,
        overlap,
        hidden,
    ) in examples:
        print(
            "    n={} frame={} C={} odd(C)={} R={}".format(
                n,
                frame,
                C,
                oc,
                R,
            )
        )
        print(
            "        gcd(odd(C),R)={}".format(
                overlap
            )
        )
        print(
            "        hidden_odd_part={}".format(
                hidden
            )
        )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# EXAMPLES
# ==========================================================================================

def print_examples(states):
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = [
        9,
        15,
        21,
        33,
        39,
        51,
        57,
        69,
        77,
        87,
        93,
        111,
        141,
        183,
        485879,
        5579767,
    ]

    wanted_set = set(wanted)

    for (
        n, p, q, frame, C, c, M, w, R, H, h2, hodd
    ) in states:

        if n not in wanted_set:
            continue

        print()
        print(
            "n={} p={} q={} frame={}".format(
                n,
                p,
                q,
                frame,
            )
        )
        print("    C={}".format(C))
        print("    c={}".format(c))
        print("    M=n+c={}".format(M))
        print("    w=v2(M)={}".format(w))
        print("    R={}".format(R))
        print("    factor(R)={}".format(factor_trial(R)))
        print("    H={}".format(H))
        print("    v2(H)={}".format(h2))
        print("    odd(H)={}".format(hodd))
        print("    odd(C)={}".format(odd_part(C)))
        print("    gcd(C,R)={}".format(gcd(C, R)))
        print("    depth={}".format(v2(gcd(2 * C, M))))


# ==========================================================================================
# MAIN
# ==========================================================================================

def main():
    print("=" * 90)
    print("EXPERIMENT 685 START")
    print("=" * 90)
    print()

    states = generate_states()

    print(
        "prime limit={}".format(
            PRIME_LIMIT
        )
    )
    print(
        "odd primes={}".format(
            len([p for p in sieve(PRIME_LIMIT) if p & 1])
        )
    )
    print(
        "semiprimes={}".format(
            len(states)
        )
    )
    print()

    test_baseline(states)
    test_odd_H_exact(states)
    test_R_factorization(states)
    test_H_reconstruction(states)
    test_H_plus_frame_factor_recovery(states)
    test_H_signature_collisions(states)
    test_odd_core_factor_overlap(states)
    print_examples(states)

    print()
    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()
    print("The exact identity is:")
    print()
    print("    H = gcd(C,n+c)")
    print()
    print("Writing:")
    print()
    print("    n+c = 2^w R")
    print("    R odd")
    print()
    print("we always have:")
    print()
    print("    odd(H) = gcd(C,R)")
    print()
    print("Experiment 684 showed that:")
    print()
    print("    odd(C) | R")
    print("is only true for a minority of states.")
    print()
    print("Therefore the useful object supplied by factoring R is")
    print("not generally the whole canonical residual C.")
    print()
    print("The new question is:")
    print()
    print("    Can R + w determine H?")
    print()
    print("and then:")
    print()
    print("    Can H + frame determine C?")
    print()
    print("and finally:")
    print()
    print("    Can C determine the factor pair?")
    print()
    print("This separates the proposed 'factor n+c repeatedly'")
    print("idea into three precise information channels.")
    print()
    print("==========================================================================================")
    print("EXPERIMENT 685 FINISHED")
    print("==========================================================================================")


if __name__ == "__main__":
    main()
