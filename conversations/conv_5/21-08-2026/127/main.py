#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from math import gcd
from statistics import mean
from sympy import factorint
from sympy.ntheory import divisors


# =============================================================================
# EXPERIMENT 688
#
# SHIFTED-TARGET HALVING + COMPLETE H-DIVISOR LATTICE
#
# Main hypothesis:
#
#     M = n + c
#     M = 2^w * R
#
# Since n and c are odd:
#
#     n + c is always even
#
# Define:
#
#     n1 = (n + c) / 2
#
# Then n1 < n.
#
# We test:
#
#   1. How much does this actually reduce the target?
#   2. Can the reduced target be recursively halved?
#   3. Does factoring n+c give every possible H divisor?
#   4. Is the true H always in that divisor lattice?
#   5. How many H candidates remain?
#   6. Can H candidates be connected to C?
#   7. Does factoring the smaller target n1 reveal C/H information?
#
# No factorization of n is used by the proposed channels.
# Known p,q are used ONLY as validation controls.
# =============================================================================


PRIME_LIMIT = 1000
SHOW_EXAMPLES = 20
MAX_RECURSION = 16


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def odd_primes(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start::p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(3, limit + 1, 2) if sieve[i]]


def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return 10**9

    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    x = abs(x)

    if x == 0:
        return 0

    return x >> v2(x)


def frame_data(p: int, q: int, n: int):
    """
    Frame selection used throughout the experiments.

    n % 4 == 3 -> Frame A
        C = q + 3
        c = 9

    n % 4 == 1 -> Frame B
        C = p + 1
        c = 3
    """

    if n % 4 == 3:
        frame = "A"
        C = q + 3
        c = 9
        p0 = 3
    else:
        frame = "B"
        C = p + 1
        c = 3
        p0 = -1

    return frame, C, c, p0


def canonical_values(p: int, q: int, n: int):
    frame, C, c, p0 = frame_data(p, q, n)

    M = n + c
    H = gcd(C, M)

    return {
        "frame": frame,
        "C": C,
        "c": c,
        "p0": p0,
        "M": M,
        "w": v2(M),
        "R": odd_part(M),
        "H": H,
        "v2H": v2(H),
        "depth": v2(gcd(2 * C, M)),
    }


# =============================================================================
# SHIFTED TARGET RECURSION
# =============================================================================

def shifted_halving_chain(n: int, c: int, max_depth: int = MAX_RECURSION):
    """
    Start from n and repeatedly apply

        x -> (x + c) / 2

    whenever x + c is even.

    The important point is that after the first step the target
    can become even, in which case another halving is impossible
    with the same odd c.

    Returns a list including the initial target.
    """

    chain = [n]
    x = n

    for _ in range(max_depth):
        if (x + c) % 2 != 0:
            break

        x = (x + c) // 2
        chain.append(x)

    return chain


# =============================================================================
# DIVISOR HELPERS
# =============================================================================

def all_divisors_from_factorization(factors: dict[int, int]) -> list[int]:
    """
    Generate the complete positive divisor lattice from factorint output.
    """
    out = [1]

    for prime, exponent in factors.items():
        previous = list(out)
        powers = [prime ** e for e in range(1, exponent + 1)]

        out = []
        for d in previous:
            for power in [1] + powers:
                out.append(d * power)

    return sorted(set(out))


def h_candidates(M: int, required_v2: int | None = None):
    """
    Every H candidate is a divisor of M.

    Optionally restrict to a known v2(H) as a CONTROL condition.
    The production-style factor-free search would not know this value,
    so both unrestricted and valuation-filtered counts are reported.
    """

    fac = factorint(M)
    ds = all_divisors_from_factorization(fac)

    if required_v2 is None:
        return ds

    return [d for d in ds if v2(d) == required_v2]


# =============================================================================
# FACTOR-OF-REDUCED-TARGET ANALYSIS
# =============================================================================

def target_factor_summary(x: int):
    """
    Factor a shifted target and return useful information.
    """
    if x <= 1:
        return {
            "x": x,
            "factorization": {},
            "odd_core": x,
            "divisors": [1],
        }

    fac = factorint(x)
    ds = all_divisors_from_factorization(fac)

    return {
        "x": x,
        "factorization": fac,
        "odd_core": odd_part(x),
        "divisors": ds,
    }


def contains_value_or_multiple(ds: list[int], value: int) -> bool:
    """
    Check whether the target divisor lattice contains value,
    or a divisor of it.
    """
    if value == 0:
        return False

    for d in ds:
        if value % d == 0:
            return True

    return False


# =============================================================================
# DATA GENERATION
# =============================================================================

def build_states():
    primes = odd_primes(PRIME_LIMIT)

    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            data = canonical_values(p, q, n)
            data["p"] = p
            data["q"] = q
            data["n"] = n

            states.append(data)

    return states


# =============================================================================
# OUTPUT HELPERS
# =============================================================================

def print_header(title: str):
    print()
    print("=" * 90)
    print(title)
    print("=" * 90)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print_header("EXPERIMENT 688 START")

    states = build_states()

    print(f"prime limit={PRIME_LIMIT}")
    print(f"odd primes={len(odd_primes(PRIME_LIMIT))}")
    print(f"semiprimes={len(states)}")

    failures = 0

    # -------------------------------------------------------------------------
    # TEST 0
    # -------------------------------------------------------------------------

    print_header("TEST 0: BASELINE")

    for s in states:
        p = s["p"]
        q = s["q"]
        n = s["n"]

        if p * q != n:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    # -------------------------------------------------------------------------
    # TEST 1
    #
    # Shifted target:
    #
    #     n1 = (n+c)/2
    #
    # Measure actual reduction.
    # -------------------------------------------------------------------------

    print_header("TEST 1: FIRST SHIFTED TARGET REDUCTION")

    reductions = []
    target_values = []
    strict_reductions = 0

    for s in states:
        n = s["n"]
        c = s["c"]

        M = n + c

        if M % 2 != 0:
            print(f"ERROR: n+c not even for n={n}")
            failures += 1
            continue

        n1 = M // 2

        if not (n1 < n):
            failures += 1
            continue

        reduction = n / n1

        reductions.append(reduction)
        target_values.append(n1)

        if n1 < n:
            strict_reductions += 1

    print(f"checked={len(states)}")
    print(f"strict reductions={strict_reductions}")
    print(f"minimum n1={min(target_values)}")
    print(f"maximum n1={max(target_values)}")
    print(f"mean reduction ratio n/n1={mean(reductions):.6f}")

    # -------------------------------------------------------------------------
    # TEST 2
    #
    # Can the same shifted-halving transformation continue?
    # -------------------------------------------------------------------------

    print_header("TEST 2: RECURSIVE SHIFTED HALVING")

    chain_lengths = Counter()
    continuation = Counter()
    max_chain = 0

    chain_examples = []

    for s in states:
        n = s["n"]
        c = s["c"]

        chain = shifted_halving_chain(n, c)

        length = len(chain) - 1
        chain_lengths[length] += 1
        max_chain = max(max_chain, length)

        if len(chain) >= 3:
            continuation[length] += 1

        if len(chain_examples) < SHOW_EXAMPLES and length >= 1:
            chain_examples.append((s, chain))

    print("halving-step distribution:")
    for depth, count in sorted(chain_lengths.items()):
        print(f"    steps={depth:2d} -> {count}")

    print(f"maximum recursive halving steps={max_chain}")

    if chain_examples:
        print()
        print("examples:")

        for s, chain in chain_examples:
            print(
                f"    n={s['n']} frame={s['frame']} c={s['c']} "
                f"chain={chain}"
            )

    # -------------------------------------------------------------------------
    # TEST 3
    #
    # Factor M=n+c and generate EVERY H divisor.
    # -------------------------------------------------------------------------

    print_header("TEST 3: COMPLETE H DIVISOR LATTICE")

    true_H_missing = 0
    all_H_candidate_counts = Counter()
    v2_filtered_counts = Counter()

    lattice_examples = []

    factor_cache: dict[int, dict[int, int]] = {}

    for s in states:
        M = s["M"]
        true_H = s["H"]

        if M not in factor_cache:
            factor_cache[M] = factorint(M)

        factors = factor_cache[M]
        ds = all_divisors_from_factorization(factors)

        all_H_candidate_counts[len(ds)] += 1

        filtered = [d for d in ds if v2(d) == s["v2H"]]
        v2_filtered_counts[len(filtered)] += 1

        if true_H not in ds:
            true_H_missing += 1

        if len(lattice_examples) < SHOW_EXAMPLES:
            lattice_examples.append(
                (
                    s,
                    factors,
                    ds,
                    filtered,
                )
            )

    print(f"unique M factorizations={len(factor_cache)}")
    print(f"true H missing={true_H_missing}")
    print("all-divisor count distribution:")

    for count, freq in sorted(all_H_candidate_counts.items()):
        print(f"    {count:3d} -> {freq}")

    print()
    print("H candidates after true-v2 control filter:")

    for count, freq in sorted(v2_filtered_counts.items()):
        print(f"    {count:3d} -> {freq}")

    # -------------------------------------------------------------------------
    # TEST 4
    #
    # Show that every H is literally one of the divisors.
    # -------------------------------------------------------------------------

    print_header("TEST 4: H IS ALWAYS A DIVISOR OF n+c")

    checked = 0

    for s in states:
        M = s["M"]
        H = s["H"]

        factors = factor_cache[M]
        ds = all_divisors_from_factorization(factors)

        if M % H != 0:
            print(
                f"FAIL H∤M: n={s['n']} M={M} H={H}"
            )
            failures += 1

        if H not in ds:
            print(
                f"FAIL H not generated: n={s['n']} M={M} H={H}"
            )
            failures += 1

        checked += 1

    print(f"checked={checked}")
    print(f"failures={failures}")

    # -------------------------------------------------------------------------
    # TEST 5
    #
    # Compare n+c factorization with n1=(n+c)/2.
    #
    # The reduced target differs only by removing one factor of 2:
    #
    #     M = 2*n1
    #
    # We test whether factorizing n1 gives useful overlap with H/C.
    # -------------------------------------------------------------------------

    print_header("TEST 5: FACTORIZATION OF THE REDUCED TARGET n1")

    reduced_factor_overlap = Counter()
    reduced_factor_contains_true_H = 0
    reduced_factor_contains_oddH = 0

    reduced_examples = []

    for s in states:
        M = s["M"]
        n1 = M // 2

        fac1 = factorint(n1)
        d1 = all_divisors_from_factorization(fac1)

        H = s["H"]
        oddH = odd_part(H)
        C = s["C"]

        if H in d1:
            reduced_factor_contains_true_H += 1

        if oddH in d1:
            reduced_factor_contains_oddH += 1

        overlap = set(fac1.keys()) & set(factor_cache[M].keys())
        reduced_factor_overlap[len(overlap)] += 1

        if len(reduced_examples) < SHOW_EXAMPLES:
            reduced_examples.append(
                (s, n1, fac1, d1)
            )

    print(
        f"n1 contains true H as divisor="
        f"{reduced_factor_contains_true_H}"
    )

    print(
        f"n1 contains odd(H) as divisor="
        f"{reduced_factor_contains_oddH}"
    )

    print("number of shared prime factors between n+c and n1:")

    for count, freq in sorted(reduced_factor_overlap.items()):
        print(f"    {count} -> {freq}")

    # -------------------------------------------------------------------------
    # TEST 6
    #
    # Does the smaller target reveal C directly?
    #
    # We test several natural channels:
    #
    #   C | n1
    #   odd(C) | n1
    #   C | some divisor of n1
    #   gcd(C,n1)
    #
    # These are intentionally diagnostic rather than assumptions.
    # -------------------------------------------------------------------------

    print_header("TEST 6: DOES THE REDUCED TARGET REVEAL C?")

    c_divides_n1 = 0
    odd_c_divides_n1 = 0
    gcd_c_n1_equals_c = 0
    gcd_c_n1_values = Counter()

    for s in states:
        n1 = s["M"] // 2
        C = s["C"]

        g = gcd(C, n1)

        gcd_c_n1_values[g] += 1

        if n1 % C == 0:
            c_divides_n1 += 1

        oC = odd_part(C)

        if oC != 0 and n1 % oC == 0:
            odd_c_divides_n1 += 1

        if g == C:
            gcd_c_n1_equals_c += 1

    print(f"C | n1={c_divides_n1}")
    print(f"odd(C) | n1={odd_c_divides_n1}")
    print(f"gcd(C,n1)=C={gcd_c_n1_equals_c}")

    # -------------------------------------------------------------------------
    # TEST 7
    #
    # Does reduced factorization make the H lattice smaller?
    # -------------------------------------------------------------------------

    print_header("TEST 7: H CANDIDATE REDUCTION AFTER HALVING")

    before_counts = []
    after_counts = []

    reduction_distribution = Counter()

    for s in states:
        M = s["M"]
        n1 = M // 2
        true_v = s["v2H"]

        before = factor_cache[M]
        ds_before = all_divisors_from_factorization(before)

        f1 = factorint(n1)
        ds_after = all_divisors_from_factorization(f1)

        before_h = [d for d in ds_before if v2(d) == true_v]
        after_h = [d for d in ds_after if v2(d) == true_v]

        before_counts.append(len(before_h))
        after_counts.append(len(after_h))

        if len(before_h) > 0:
            ratio_num = len(after_h)
            ratio_den = len(before_h)

            if ratio_num == 0:
                reduction_distribution["zero"] += 1
            elif ratio_num < ratio_den:
                reduction_distribution["smaller"] += 1
            elif ratio_num == ratio_den:
                reduction_distribution["same"] += 1
            else:
                reduction_distribution["larger"] += 1

    print(f"mean H candidates before={mean(before_counts):.4f}")
    print(f"mean H candidates after ={mean(after_counts):.4f}")

    print("candidate-count relation:")
    for label, count in reduction_distribution.items():
        print(f"    {label} -> {count}")

    # -------------------------------------------------------------------------
    # TEST 8
    #
    # Factor-free conceptual recovery experiment:
    #
    # We know only:
    #
    #   frame
    #   M=n+c factorization
    #   all H|M
    #
    # For every H, test simple C multiples near the actual scale using
    # divisibility conditions derived from the frame.
    #
    # This is a CONTROL search and intentionally does NOT use factor(n).
    #
    # To keep runtime bounded, C=H*t is searched up to C <= n.
    # We additionally restrict to C values satisfying:
    #
    #   Frame A: C-3 divides n
    #   Frame B: C-1 divides n
    #
    # That final test is an integer divisibility test, NOT factorint(n).
    # =============================================================================

    print_header("TEST 8: DIRECT H -> C RECOVERY WITHOUT factorint(n)")

    recover_none = 0
    recover_unique = 0
    recover_ambiguous = 0

    direct_recovery_examples = []

    for s in states:
        n = s["n"]
        frame = s["frame"]
        C_true = s["C"]
        H_true = s["H"]

        M = s["M"]
        ds = all_divisors_from_factorization(factor_cache[M])

        candidates = []

        # Restrict H candidates to correct 2-adic level as a CONTROL
        # because previous experiments establish that this level is
        # recoverable from the canonical depth law.
        Hs = [h for h in ds if v2(h) == s["v2H"]]

        # Search multiples H*t, but only over the frame-valid factor side.
        #
        # Instead of iterating to n/H for every H, use the fact that
        # C corresponds to a factor:
        #
        # Frame A: q = C - 3
        # Frame B: p = C - 1
        #
        # The divisibility test is exact.
        #
        # We cap t by n//H, which is safe but potentially large.
        # To avoid pathological runtime, we search only up to sqrt(n)
        # multiples and separately test the complementary large-C branch.
        limit = int(n ** 0.5) + 2

        seen = set()

        for H in Hs:
            max_t = min(n // H, limit)

            for t in range(1, max_t + 1):
                C = H * t

                if C <= 0:
                    continue

                if frame == "A":
                    q_candidate = C - 3

                    if q_candidate <= 1:
                        continue

                    if n % q_candidate == 0:
                        p_candidate = n // q_candidate

                        if p_candidate * q_candidate == n:
                            item = (
                                C,
                                p_candidate,
                                q_candidate,
                            )
                            seen.add(item)

                else:
                    p_candidate = C - 1

                    if p_candidate <= 1:
                        continue

                    if n % p_candidate == 0:
                        q_candidate = n // p_candidate

                        if p_candidate * q_candidate == n:
                            item = (
                                C,
                                p_candidate,
                                q_candidate,
                            )
                            seen.add(item)

        candidates = sorted(seen)

        if not candidates:
            recover_none += 1
        elif len(candidates) == 1:
            recover_unique += 1
        else:
            recover_ambiguous += 1

        if (
            len(direct_recovery_examples) < SHOW_EXAMPLES
            and candidates
        ):
            direct_recovery_examples.append(
                (
                    s,
                    candidates,
                )
            )

    print(f"none      = {recover_none}")
    print(f"unique    = {recover_unique}")
    print(f"ambiguous = {recover_ambiguous}")

    if direct_recovery_examples:
        print()
        print("recovery examples:")

        for s, candidates in direct_recovery_examples:
            print(
                f"    n={s['n']} frame={s['frame']} "
                f"true_C={s['C']} true_H={s['H']}"
            )
            print(f"        candidates={candidates}")

    # -------------------------------------------------------------------------
    # TEST 9
    #
    # Representative examples.
    # -------------------------------------------------------------------------

    print_header("TEST 9: REPRESENTATIVE EXAMPLES")

    for s in states[:SHOW_EXAMPLES]:
        n = s["n"]
        c = s["c"]
        M = s["M"]
        H = s["H"]
        C = s["C"]

        chain = shifted_halving_chain(n, c)

        facM = factor_cache[M]
        dsM = all_divisors_from_factorization(facM)

        print()
        print(
            f"n={n} p={s['p']} q={s['q']} frame={s['frame']}"
        )
        print(f"    C={C}")
        print(f"    c={c}")
        print(f"    M=n+c={M}")
        print(f"    factor(M)={facM}")
        print(f"    w=v2(M)={s['w']}")
        print(f"    R={s['R']}")
        print(f"    H={H}")
        print(f"    H | M={M % H == 0}")
        print(f"    number of H divisors={len(dsM)}")
        print(f"    H candidates at true v2={len([d for d in dsM if v2(d) == s['v2H']])}")
        print(f"    n halving chain={chain}")

    # -------------------------------------------------------------------------
    # FINAL SUMMARY
    # -------------------------------------------------------------------------

    print_header("FINAL STRUCTURAL SUMMARY")

    print(
        """
For every state:

    M = n + c
    M = 2^w R
    R odd

Because n and c are odd:

    M is always even

so:

    n1 = (n+c)/2 < n

and:

    M = 2*n1.

The exact gcd object remains:

    H = gcd(C,M)

Therefore:

    H | M

and complete factorization of M generates the COMPLETE
divisor lattice containing the true H.

The experiment therefore separates two questions:

    1. EXISTENCE:
         Is the true H available from factor(n+c)?

       Yes, necessarily.

    2. SELECTION:
         Can we identify which divisor H is the canonical one
         without already factoring n?

The new shifted-target channel is:

    n
      -> (n+c)/2
      -> factor the smaller target
      -> inspect its divisors / prime factors
      -> compare against H, odd(H), C, odd(C)

The recursive step is only available when the new target plus
the same odd shift is again even. The experiment measures whether
that produces a genuine recursive descent or whether the process
usually terminates after the first halving.

The critical next question is therefore:

    Does factoring the reduced target n1 provide information
    that selects the correct H among the divisors of n+c?

"""
    )

    print(f"TOTAL FAILURES = {failures}")
    print("EXPERIMENT 688 FINISHED")


if __name__ == "__main__":
    main()
