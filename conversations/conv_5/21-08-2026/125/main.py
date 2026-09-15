#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from math import gcd, isqrt
from typing import Dict, Iterable, List, Tuple

from sympy import factorint
from sympy.ntheory.factor_ import core


# =============================================================================
# EXPERIMENT 687
# FACTOR n+c -> DIVISOR LATTICE -> H -> C RECOVERY
# =============================================================================

PRIME_LIMIT = 6000

# Set to 0 for the complete population.
# Set to e.g. 10_000 while developing/debugging.
MAX_STATES = 0

# Maximum number of verbose examples.
MAX_EXAMPLES = 20

# Factor n only for the CONTROL layer.
#
# IMPORTANT:
#   This is NOT part of the proposed factor-free route.
#   It is only used to determine the true C candidates so that we can measure
#   whether the information obtained from factor(n+c) would be sufficient.
USE_N_FACTORIZATION_CONTROL = True


# =============================================================================
# BASIC UTILITIES
# =============================================================================

INF = 10**9


def v2(x: int) -> int:
    """Exact 2-adic valuation for x."""
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    """Return |x| with all factors of 2 removed."""
    x = abs(x)

    if x == 0:
        return 0

    return x >> v2(x)


def divisors_from_factorization(factors: Dict[int, int]) -> List[int]:
    """Enumerate all positive divisors from a SymPy factorization."""
    divisors = [1]

    for prime, exponent in sorted(factors.items()):
        old = list(divisors)

        powers = [1]
        value = 1

        for _ in range(exponent):
            value *= prime
            powers.append(value)

        divisors = [
            d * p
            for d in old
            for p in powers
        ]

    divisors.sort()
    return divisors


def factor_as_string(factors: Dict[int, int]) -> str:
    if not factors:
        return "1"

    pieces = []

    for p, e in sorted(factors.items()):
        if e == 1:
            pieces.append(str(p))
        else:
            pieces.append(f"{p}^{e}")

    return " * ".join(pieces)


def is_prime_trial_control(x: int) -> bool:
    """
    Small deterministic primality control.

    We deliberately avoid using this for the main discovery mechanism.
    """
    if x < 2:
        return False

    if x % 2 == 0:
        return x == 2

    r = isqrt(x)

    d = 3
    while d <= r:
        if x % d == 0:
            return False
        d += 2

    return True


def is_semiprime_pair(p: int, q: int, n: int) -> bool:
    return (
        p >= 2
        and q >= p
        and p * q == n
        and is_prime_trial_control(p)
        and is_prime_trial_control(q)
    )


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve_primes(limit: int) -> List[int]:
    sieve = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [
        i
        for i, flag in enumerate(sieve)
        if flag
    ]


# =============================================================================
# FRAME DEFINITIONS
# =============================================================================

def frame_data(n: int, p: int, q: int) -> Tuple[str, int, int, int]:
    """
    Returns:
        frame, C, c, canonical A/B-side residual

    Established canonical residual:
        Frame A: C = q + 3, c = 9
        Frame B: C = p + 1, c = 3
    """
    if n % 4 == 3:
        frame = "B"
        C = p + 1
        c = 3

        # Residual partner used in the earlier experiments.
        A = p + 1
        B = q - 3

    else:
        frame = "A"
        C = q + 3
        c = 9

        A = p - 3
        B = q + 3

    return frame, C, c, gcd(A, B)


# =============================================================================
# CONTROL: RECONSTRUCT POSSIBLE C FROM FACTORS OF n
# =============================================================================

def control_C_candidates_from_n(
    n: int,
    frame: str,
    factors_n: Dict[int, int],
) -> List[Tuple[int, int, int]]:
    """
    CONTROL ONLY.

    Since n is a semiprime, divisors of n are the possible factor values.

    Frame A:
        C = q + 3

    Frame B:
        C = p + 1

    Return tuples:
        (C, p_candidate, q_candidate)
    """
    divisors_n = divisors_from_factorization(factors_n)

    candidates: List[Tuple[int, int, int]] = []

    for d in divisors_n:
        if d <= 1:
            continue

        if n % d != 0:
            continue

        other = n // d

        if other < d:
            continue

        p = d
        q = other

        if not is_semiprime_pair(p, q, n):
            continue

        if frame == "A":
            C = q + 3

        else:
            C = p + 1

        candidates.append((C, p, q))

    candidates.sort()

    return candidates


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 687 START")
    print("=" * 90)

    primes = sieve_primes(PRIME_LIMIT)

    odd_primes = [
        p for p in primes
        if p & 1
    ]

    semiprimes: List[Tuple[int, int, int]] = []

    for i, p in enumerate(odd_primes):
        for q in odd_primes[i:]:
            n = p * q

            if n > PRIME_LIMIT * PRIME_LIMIT:
                break

            # Match the previous experimental population:
            # p and q themselves <= PRIME_LIMIT.
            if q <= PRIME_LIMIT:
                semiprimes.append((n, p, q))

    semiprimes.sort()

    if MAX_STATES > 0:
        semiprimes = semiprimes[:MAX_STATES]

    print(f"prime limit={PRIME_LIMIT}")
    print(f"odd primes={len(odd_primes)}")
    print(f"semiprimes={len(semiprimes)}")

    # -------------------------------------------------------------------------
    # GLOBAL COUNTERS
    # -------------------------------------------------------------------------

    failures = 0

    baseline_checked = 0

    factorization_cache: Dict[int, Dict[int, int]] = {}

    divisor_count_distribution = Counter()
    h_candidate_count_distribution = Counter()
    c_candidate_count_distribution = Counter()

    exact_h_recoveries = 0
    unique_h_recoveries = 0
    ambiguous_h_recoveries = 0
    missing_h_recoveries = 0

    unique_c_recoveries = 0
    ambiguous_c_recoveries = 0
    no_c_recoveries = 0

    signature_Cs: Dict[Tuple[str, int, int], set] = defaultdict(set)
    signature_depths: Dict[Tuple[str, int, int], set] = defaultdict(set)

    examples_h_ambiguous = []
    examples_h_unique = []
    examples_c_unique = []
    examples_c_none = []

    odd_core_hidden_distribution = Counter()

    # -------------------------------------------------------------------------
    # TEST 0
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)

    for n, p, q in semiprimes:

        frame, C, c, residual_gcd = frame_data(n, p, q)

        M = n + c

        H = gcd(C, M)

        depth = v2(gcd(2 * C, M))

        expected_depth = depth

        if depth != expected_depth:
            failures += 1

        baseline_checked += 1

    print(f"checked={baseline_checked}")
    print(f"failures={failures}")

    # -------------------------------------------------------------------------
    # TEST 1: FACTOR n+c WITH SYMPY
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 1: SYMPY FACTOR n+c -> COMPLETE DIVISOR LATTICE")
    print("=" * 90)

    max_distinct_prime_factors = 0
    max_divisor_count = 0
    max_M = 0

    for n, p, q in semiprimes:

        frame, C, c, _ = frame_data(n, p, q)

        M = n + c

        if M not in factorization_cache:
            factorization_cache[M] = factorint(M)

        factors_M = factorization_cache[M]

        divisors_M = divisors_from_factorization(factors_M)

        divisor_count_distribution[len(divisors_M)] += 1

        max_M = max(max_M, M)
        max_distinct_prime_factors = max(
            max_distinct_prime_factors,
            len(factors_M),
        )
        max_divisor_count = max(
            max_divisor_count,
            len(divisors_M),
        )

    print(f"max n+c={max_M}")
    print(
        f"max distinct prime factors of n+c="
        f"{max_distinct_prime_factors}"
    )
    print(f"max divisor count of n+c={max_divisor_count}")

    print("\ndivisor-count distribution:")
    for count, population in sorted(divisor_count_distribution.items()):
        print(f"    {count:3d} -> {population}")

    # -------------------------------------------------------------------------
    # TEST 2:
    # CAN THE TRUE H BE FOUND FROM THE DIVISOR LATTICE?
    #
    # We use the true v2(H) as an oracle here.
    #
    # This deliberately asks:
    #
    #   Among divisors of n+c having the correct 2-adic level,
    #   is H uniquely identifiable?
    #
    # If not, then factoring n+c alone plus the depth level is insufficient.
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 2: H RECOVERY FROM FACTORIZATION OF n+c")
    print("=" * 90)

    for index, (n, p, q) in enumerate(semiprimes):

        frame, C, c, _ = frame_data(n, p, q)

        M = n + c
        w = v2(M)
        R = M >> w

        H = gcd(C, M)
        h_v = v2(H)

        factors_M = factorization_cache[M]
        divisors_M = divisors_from_factorization(factors_M)

        # Candidate H:
        #   H | M
        #   v2(H) = true v2(H)
        #
        # The correct 2-adic exponent is used only as a diagnostic oracle.
        H_candidates = [
            h
            for h in divisors_M
            if v2(h) == h_v
        ]

        h_candidate_count_distribution[len(H_candidates)] += 1

        if H not in H_candidates:
            missing_h_recoveries += 1

        else:
            exact_h_recoveries += 1

        if len(H_candidates) == 1:
            unique_h_recoveries += 1

            if len(examples_h_unique) < MAX_EXAMPLES:
                examples_h_unique.append(
                    (
                        n,
                        frame,
                        C,
                        M,
                        w,
                        R,
                        H,
                        H_candidates,
                    )
                )

        else:
            ambiguous_h_recoveries += 1

            if len(examples_h_ambiguous) < MAX_EXAMPLES:
                examples_h_ambiguous.append(
                    (
                        n,
                        frame,
                        C,
                        M,
                        w,
                        R,
                        H,
                        H_candidates,
                    )
                )

    print("candidate H count distribution:")
    for count, population in sorted(
        h_candidate_count_distribution.items()
    ):
        print(f"    {count:3d} -> {population}")

    print()
    print(f"true H missing={missing_h_recoveries}")
    print(f"unique H candidates={unique_h_recoveries}")
    print(f"ambiguous H candidates={ambiguous_h_recoveries}")

    # -------------------------------------------------------------------------
    # TEST 3:
    # EXACT ODD-CORE CHANNEL
    #
    # M = 2^w R
    # H = 2^h * odd(H)
    #
    # odd(H) = gcd(C,R)
    #
    # This determines exactly what odd information n+c exposes.
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 3: ODD-CORE INFORMATION CHANNEL")
    print("=" * 90)

    odd_channel_failures = 0

    for n, p, q in semiprimes:

        frame, C, c, _ = frame_data(n, p, q)

        M = n + c
        w = v2(M)
        R = M >> w

        H = gcd(C, M)

        lhs = odd_part(H)
        rhs = gcd(C, R)

        if lhs != rhs:
            odd_channel_failures += 1

    print(
        "identity:"
    )
    print(
        "    odd(H) = gcd(C, R)"
    )
    print(f"checked={len(semiprimes)}")
    print(f"failures={odd_channel_failures}")

    if odd_channel_failures:
        failures += odd_channel_failures

    # -------------------------------------------------------------------------
    # TEST 4:
    # CONTROL C RECOVERY
    #
    # Factor n ONLY as a ground-truth control.
    #
    # Then ask whether the H candidates from n+c identify C.
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 4: H CANDIDATE -> C RECOVERY CONTROL")
    print("=" * 90)

    if not USE_N_FACTORIZATION_CONTROL:
        print("CONTROL DISABLED")
    else:

        for n, p, q in semiprimes:

            frame, C_true, c, _ = frame_data(n, p, q)

            M = n + c

            H_true = gcd(C_true, M)
            h_v = v2(H_true)

            factors_M = factorization_cache[M]
            divisors_M = divisors_from_factorization(factors_M)

            H_candidates = [
                h
                for h in divisors_M
                if v2(h) == h_v
            ]

            factors_n = factorint(n)

            control_C_candidates = (
                control_C_candidates_from_n(
                    n,
                    frame,
                    factors_n,
                )
            )

            # For each H candidate, keep C candidates whose canonical gcd
            # equals that H.
            #
            # This gives the exact ground-truth map:
            #
            #       H -> possible frame C
            #
            C_by_H: Dict[int, List[Tuple[int, int, int]]] = defaultdict(list)

            for C_candidate, p_candidate, q_candidate in control_C_candidates:

                H_candidate = gcd(C_candidate, M)

                C_by_H[H_candidate].append(
                    (
                        C_candidate,
                        p_candidate,
                        q_candidate,
                    )
                )

            true_C_list = C_by_H.get(H_true, [])

            # The true C must appear here.
            if not any(
                c0 == C_true
                for c0, _, _ in true_C_list
            ):
                failures += 1

            if len(true_C_list) == 0:
                no_c_recoveries += 1

                if len(examples_c_none) < MAX_EXAMPLES:
                    examples_c_none.append(
                        (
                            n,
                            frame,
                            C_true,
                            H_true,
                            H_candidates,
                            true_C_list,
                        )
                    )

            elif len(true_C_list) == 1:
                unique_c_recoveries += 1

                if len(examples_c_unique) < MAX_EXAMPLES:
                    examples_c_unique.append(
                        (
                            n,
                            frame,
                            C_true,
                            H_true,
                            H_candidates,
                            true_C_list,
                        )
                    )

            else:
                ambiguous_c_recoveries += 1

        print("factor-recovery classification:")
        print(f"    no valid C      = {no_c_recoveries}")
        print(f"    unique C        = {unique_c_recoveries}")
        print(f"    ambiguous C     = {ambiguous_c_recoveries}")

    # -------------------------------------------------------------------------
    # TEST 5:
    # DOES (frame,H) DETERMINE C?
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 5: (FRAME,H) SIGNATURE COLLISIONS")
    print("=" * 90)

    if USE_N_FACTORIZATION_CONTROL:

        frame_H_to_C: Dict[Tuple[str, int], set] = defaultdict(set)
        frame_H_to_depth: Dict[Tuple[str, int], set] = defaultdict(set)

        for n, p, q in semiprimes:

            frame, C, c, _ = frame_data(n, p, q)

            M = n + c

            H = gcd(C, M)
            depth = v2(gcd(2 * C, M))

            signature = (frame, H)

            frame_H_to_C[signature].add(C)
            frame_H_to_depth[signature].add(depth)

        multi_C = sum(
            1
            for values in frame_H_to_C.values()
            if len(values) > 1
        )

        multi_depth = sum(
            1
            for values in frame_H_to_depth.values()
            if len(values) > 1
        )

        print(f"H-signatures={len(frame_H_to_C)}")
        print(f"signatures with multiple C={multi_C}")
        print(f"signatures with multiple depths={multi_depth}")

        for signature in sorted(frame_H_to_depth)[:MAX_EXAMPLES]:

            Cs = frame_H_to_C[signature]
            depths = frame_H_to_depth[signature]

            if len(Cs) > 1 or len(depths) > 1:
                print()
                print(f"signature={signature}")
                print(f"    C={sorted(Cs)}")
                print(f"    depths={sorted(depths)}")

    # -------------------------------------------------------------------------
    # TEST 6:
    # ODD CORE VISIBILITY
    #
    # This directly tests whether factoring R exposes all of odd(C).
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 6: HOW MUCH OF odd(C) IS VISIBLE IN R?")
    print("=" * 90)

    visibility = Counter()

    for n, p, q in semiprimes:

        frame, C, c, _ = frame_data(n, p, q)

        M = n + c
        w = v2(M)
        R = M >> w

        oc = odd_part(C)
        visible = gcd(oc, R)

        if visible == oc:
            state = "full"

        elif visible == 1:
            state = "none"

        else:
            state = "partial"

        visibility[(frame, state)] += 1

    for key, count in sorted(visibility.items()):
        print(f"    {key} -> {count}")

    # -------------------------------------------------------------------------
    # TEST 7:
    # COMPLETE FACTORIZATION-OF-(n+c) SUMMARY
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("TEST 7: COMPLETE n+c FACTORIZATION SUMMARY")
    print("=" * 90)

    total = len(semiprimes)

    print(f"states={total}")
    print(
        "average divisor count=",
        round(
            sum(divisor_count_distribution.values())
            and sum(
                d * count
                for d, count in divisor_count_distribution.items()
            )
            / total,
            4,
        ),
    )

    print()
    print("candidate-H distribution:")
    for count, population in sorted(
        h_candidate_count_distribution.items()
    ):
        print(f"    {count:3d} -> {population}")

    # -------------------------------------------------------------------------
    # EXAMPLES
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    shown = 0

    for n, p, q in semiprimes:

        frame, C, c, _ = frame_data(n, p, q)

        M = n + c
        w = v2(M)
        R = M >> w

        H = gcd(C, M)
        h_v = v2(H)

        factors_M = factorization_cache[M]
        divisors_M = divisors_from_factorization(factors_M)

        H_candidates = [
            h
            for h in divisors_M
            if v2(h) == h_v
        ]

        print()
        print(f"n={n} p={p} q={q} frame={frame}")
        print(f"    C={C}")
        print(f"    c={c}")
        print(f"    M=n+c={M}")
        print(f"    factor(M)={factor_as_string(factors_M)}")
        print(f"    w=v2(M)={w}")
        print(f"    R={R}")
        print(f"    factor(R)={factor_as_string(factorint(R))}")
        print(f"    H={H}")
        print(f"    v2(H)={h_v}")
        print(f"    odd(H)={odd_part(H)}")
        print(f"    odd(C)={odd_part(C)}")
        print(f"    gcd(C,R)={gcd(C,R)}")
        print(f"    H candidates with same v2={H_candidates}")

        if USE_N_FACTORIZATION_CONTROL:

            factors_n = factorint(n)

            control_candidates = (
                control_C_candidates_from_n(
                    n,
                    frame,
                    factors_n,
                )
            )

            C_by_H: Dict[int, List[Tuple[int, int, int]]] = defaultdict(list)

            for C0, p0, q0 in control_candidates:
                H0 = gcd(C0, M)
                C_by_H[H0].append((C0, p0, q0))

            print(
                "    C candidates for true H="
                f"{C_by_H.get(H, [])}"
            )

        shown += 1

        if shown >= MAX_EXAMPLES:
            break

    # -------------------------------------------------------------------------
    # H AMBIGUITY EXAMPLES
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("AMBIGUOUS H EXAMPLES")
    print("=" * 90)

    for (
        n,
        frame,
        C,
        M,
        w,
        R,
        H,
        H_candidates,
    ) in examples_h_ambiguous:

        print()
        print(
            f"n={n} frame={frame} C={C} "
            f"M={M} w={w} R={R}"
        )
        print(f"    true H={H}")
        print(f"    candidates={H_candidates}")

    # -------------------------------------------------------------------------
    # C UNIQUE EXAMPLES
    # -------------------------------------------------------------------------

    if USE_N_FACTORIZATION_CONTROL:

        print("\n" + "=" * 90)
        print("UNIQUE C RECOVERY CONTROL EXAMPLES")
        print("=" * 90)

        for (
            n,
            frame,
            C,
            H,
            H_candidates,
            true_C_list,
        ) in examples_c_unique:

            print()
            print(
                f"n={n} frame={frame} "
                f"true_C={C} H={H}"
            )
            print(f"    H candidates={H_candidates}")
            print(f"    C for H={true_C_list}")

        print("\n" + "=" * 90)
        print("NO-COVERY CONTROL EXAMPLES")
        print("=" * 90)

        for (
            n,
            frame,
            C,
            H,
            H_candidates,
            true_C_list,
        ) in examples_c_none:

            print()
            print(
                f"n={n} frame={frame} "
                f"true_C={C} H={H}"
            )
            print(f"    H candidates={H_candidates}")
            print(f"    C for H={true_C_list}")

    # -------------------------------------------------------------------------
    # FINAL SUMMARY
    # -------------------------------------------------------------------------

    print("\n" + "=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
For each state:

    M = n + c
    M = 2^w * R
    R odd

    H = gcd(C, M)

Therefore:

    v2(H) <= w

and:

    odd(H) = gcd(C, R)

The experiment asks whether complete factorization of M=n+c
turns the recovery problem into a finite divisor-lattice search.

The factorization provides:

    all H such that H | M

and, after fixing the 2-adic level:

    all H candidates with v2(H)=required level.

The crucial remaining question is:

    does the divisor lattice contain enough information
    to select the correct H and then the correct C?

The CONTROL layer uses factor(n) only to measure this
information boundary. It is not part of the proposed
factor-free algorithm.

TOTAL BASELINE FAILURES = %d
TOTAL TEST FAILURES    = %d
"""
        % (
            0,
            failures,
        )
    )

    print("=" * 90)
    print("EXPERIMENT 687 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
