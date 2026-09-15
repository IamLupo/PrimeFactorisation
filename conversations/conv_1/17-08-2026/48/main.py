#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 115
CONSTRUCTIVE CRT COLLISION / FORCED N-RESIDUE EQUIVALENCE
SAME N MOD M, DIFFERENT S / T / PAPER QUOTIENTS
CRT COLLISION -> SQUARE-DIFFERENCE BRIDGE
NO FACTOR-PAIR SEARCH IN INVERSE STAGE
NO CSV
NO SKLEARN
==============================================================================

GOAL
----

Experiment 114 showed that random residue fingerprints eventually become
unique simply because the fingerprint space becomes large.

This experiment removes that ambiguity.

We deliberately construct many DIFFERENT semiprimes

    N_i = p_i q_i

such that

    N_i == N_0 (mod M)

for exactly the same modulus M.

We do this by fixing

    p_i == a (mod M)
    q_i == b (mod M)

so that

    p_i q_i == a*b (mod M)

for every target.

We then test whether the following quantities vary despite identical N mod M:

    S = p + q
    T = S^2 - 3N
    T-N = (p-q)^2
    Q_(1,3)
    Q_(1,5)
    Q_(1,7)
    Q_(3,5)
    Q_(3,7)
    Q_(5,7)

The square-root bridge is also tested through

    (p-q)^2 == T-N (mod M).

The important question is:

    Same N mod M
        +
    different S / Q residues
        ->
    bounded N-residue information is insufficient.

This is a CONSTRUCTIVE collision experiment, not a random-collision experiment.

No factorization is performed from N.
p,q are used only to construct and validate the oracle family.
==============================================================================

"""

from __future__ import annotations

import math
import random
import sys
import time

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


# =============================================================================
# PARAMETERS
# =============================================================================

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

# Collision moduli.
MODULI = (
    5 * 7,                 # 35
    5 * 7 * 11,            # 385
    5 * 7 * 11 * 13,       # 5005
)

# Number of pairs requested per residue-class experiment.
PAIRS_PER_CLASS = 20

# Number of distinct residue-class constructions.
CLASS_TRIALS = 12

# Need enough primes in each residue class.
MIN_PRIMES_PER_CLASS = 8

DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

RNG_SEED = 115


# =============================================================================
# DATA TYPES
# =============================================================================

@dataclass(frozen=True)
class PrimePair:
    p: int
    q: int

    @property
    def N(self) -> int:
        return self.p * self.q

    @property
    def S(self) -> int:
        return self.p + self.q

    @property
    def diff(self) -> int:
        return self.p - self.q


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi < lo:
        return []

    root = math.isqrt(hi)

    base = bytearray(b"\x01") * (root + 1)
    base[:2] = b"\x00\x00"

    for p in range(2, math.isqrt(root) + 1):
        if base[p]:
            start = p * p
            base[start:root + 1:p] = (
                b"\x00" *
                (((root - start) // p) + 1)
            )

    base_primes = [
        p for p in range(2, root + 1)
        if base[p]
    ]

    segment_size = 250_000
    result: List[int] = []

    for left in range(lo, hi + 1, segment_size):
        right = min(left + segment_size - 1, hi)
        segment = bytearray(b"\x01") * (right - left + 1)

        for p in base_primes:
            first = max(
                p * p,
                ((left + p - 1) // p) * p,
            )

            if first > right:
                continue

            segment[
                first - left:
                right - left + 1:
                p
            ] = (
                b"\x00" *
                (((right - first) // p) + 1)
            )

        result.extend(
            left + i
            for i, flag in enumerate(segment)
            if flag
        )

    return result


# =============================================================================
# PAPER QUOTIENTS
# =============================================================================

def raw_detector(
    k: int,
    ell: int,
    p: int,
    q: int,
) -> int:
    return (
        p**k * (q + 1)**ell
        + q**k * (p + 1)**ell
        - p**ell * (q + 1)**k
        - q**ell * (p + 1)**k
    )


def quotient_value(
    k: int,
    ell: int,
    p: int,
    q: int,
) -> int:
    numerator = raw_detector(k, ell, p, q)
    divisor = p + q + 1

    if numerator % divisor != 0:
        raise ArithmeticError(
            f"Detector ({k},{ell}) not divisible for "
            f"p={p}, q={q}"
        )

    return numerator // divisor


# =============================================================================
# RESIDUE CLASS CONSTRUCTION
# =============================================================================

def valid_residues(modulus: int) -> List[int]:
    """
    Prime residues must be invertible mod M because p,q > M and prime.
    """
    return [
        r
        for r in range(1, modulus)
        if math.gcd(r, modulus) == 1
    ]


def prime_classes(
    primes: Sequence[int],
    modulus: int,
) -> Dict[int, List[int]]:
    classes: Dict[int, List[int]] = defaultdict(list)

    for p in primes:
        classes[p % modulus].append(p)

    return dict(classes)


def choose_residue_pair(
    classes: Dict[int, List[int]],
    modulus: int,
    rng: random.Random,
) -> Tuple[int, int]:
    """
    Choose a,b such that both residue classes contain enough primes.

    No factor search is involved.
    """
    candidates = [
        r
        for r, values in classes.items()
        if math.gcd(r, modulus) == 1
        and len(values) >= MIN_PRIMES_PER_CLASS
    ]

    if len(candidates) < 2:
        raise RuntimeError(
            f"Insufficient populated residue classes mod {modulus}."
        )

    for _ in range(10_000):
        a = rng.choice(candidates)
        b = rng.choice(candidates)

        if a != b:
            return a, b

    raise RuntimeError(
        f"Could not select residue pair mod {modulus}."
    )


def build_forced_collision_family(
    primes: Sequence[int],
    modulus: int,
    a: int,
    b: int,
    count: int,
) -> List[PrimePair]:
    classes = prime_classes(
        primes,
        modulus,
    )

    P = classes.get(a, [])
    Q = classes.get(b, [])

    if len(P) < count or len(Q) < count:
        raise RuntimeError(
            f"Not enough primes in residue classes "
            f"a={a}, b={b} mod {modulus}: "
            f"|P|={len(P)}, |Q|={len(Q)}"
        )

    # Pair different positions so that N and S vary.
    pairs: List[PrimePair] = []

    for i in range(count):
        p = P[i]

        # Reverse one side to avoid a simple synchronized walk.
        q = Q[count - 1 - i]

        if p == q:
            continue

        pairs.append(
            PrimePair(p=p, q=q)
        )

    if len(pairs) < count // 2:
        raise RuntimeError(
            "Constructed collision family too small."
        )

    return pairs


# =============================================================================
# EXACT COLLISION VALIDATION
# =============================================================================

def validate_same_N_residue(
    pairs: Sequence[PrimePair],
    modulus: int,
) -> None:
    if not pairs:
        raise ValueError("Empty pair family.")

    target_residue = pairs[0].N % modulus

    for pair in pairs:
        if pair.N % modulus != target_residue:
            raise ArithmeticError(
                "CRT collision construction failed."
            )

        if pair.p % modulus != pairs[0].p % modulus:
            raise ArithmeticError(
                "p residue construction failed."
            )

        if pair.q % modulus != pairs[0].q % modulus:
            raise ArithmeticError(
                "q residue construction failed."
            )


def distinct_count(
    values: Sequence[int],
    modulus: int,
) -> int:
    return len({
        value % modulus
        for value in values
    })


# =============================================================================
# COLLISION ANALYSIS
# =============================================================================

def analyze_collision_family(
    pairs: Sequence[PrimePair],
    modulus: int,
) -> None:
    print()
    print(f"  M = {modulus}")
    print("  " + "-" * 74)

    validate_same_N_residue(
        pairs,
        modulus,
    )

    N0 = pairs[0].N % modulus

    S_values = [p.S for p in pairs]
    T_values = [
        p.S * p.S - 3 * p.N
        for p in pairs
    ]
    D2_values = [
        p.diff * p.diff
        for p in pairs
    ]

    print(
        f"  forced N residue      = {N0}"
    )
    print(
        f"  number of pairs       = {len(pairs)}"
    )
    print(
        f"  distinct N residues   = "
        f"{distinct_count([p.N for p in pairs], modulus)}"
    )
    print(
        f"  distinct S residues   = "
        f"{distinct_count(S_values, modulus)}"
    )
    print(
        f"  distinct T residues   = "
        f"{distinct_count(T_values, modulus)}"
    )
    print(
        f"  distinct (p-q)^2      = "
        f"{distinct_count(D2_values, modulus)}"
    )

    # Exact square-difference bridge.
    bridge_failures = 0

    for pair in pairs:
        T = pair.S * pair.S - 3 * pair.N

        lhs = (T - pair.N) % modulus
        rhs = (pair.diff * pair.diff) % modulus

        if lhs != rhs:
            bridge_failures += 1

    print(
        f"  T-N=(p-q)^2 failures  = "
        f"{bridge_failures}/{len(pairs)}"
    )

    if bridge_failures:
        raise ArithmeticError(
            "Square-difference bridge failed."
        )

    # Detector variation.
    for detector in DETECTORS:
        k, ell = detector

        Q_values = [
            quotient_value(
                k,
                ell,
                pair.p,
                pair.q,
            )
            for pair in pairs
        ]

        unique_Q = distinct_count(
            Q_values,
            modulus,
        )

        print(
            f"  Q_({k},{ell}) distinct mod M = "
            f"{unique_Q}"
        )


# =============================================================================
# CRT CONFLICT TEST
# =============================================================================

def conflict_test(
    pairs: Sequence[PrimePair],
    modulus: int,
) -> None:
    """
    Strongest test:

        N is identical modulo M

    while output residues are allowed to vary.

    Therefore no function

        f(N mod M)

    can recover the varying output residue.
    """
    print()
    print("  DIRECT CONFLICT TEST")
    print("  " + "-" * 74)

    N_residue = pairs[0].N % modulus

    assert all(
        pair.N % modulus == N_residue
        for pair in pairs
    )

    S_residues = {
        pair.S % modulus
        for pair in pairs
    }

    print(
        f"  common N mod M       = {N_residue}"
    )
    print(
        f"  number of S mod M    = {len(S_residues)}"
    )

    conflict_count = 0

    for detector in DETECTORS:
        residues = {
            quotient_value(
                detector[0],
                detector[1],
                pair.p,
                pair.q,
            ) % modulus
            for pair in pairs
        }

        print(
            f"  Q_{detector} values mod M = "
            f"{len(residues)}"
        )

        if len(residues) > 1:
            conflict_count += 1

    print(
        f"  detectors with genuine "
        f"N-mod-M conflict = {conflict_count}/{len(DETECTORS)}"
    )


# =============================================================================
# CHAINED MODULUS EXPERIMENT
# =============================================================================

def modulus_experiment(
    primes: Sequence[int],
    modulus: int,
    rng: random.Random,
) -> None:
    print()
    print("=" * 78)
    print(f"CONSTRUCTIVE CRT FAMILY: M={modulus}")
    print("=" * 78)

    classes = prime_classes(
        primes,
        modulus,
    )

    populated = [
        (r, len(values))
        for r, values in classes.items()
        if len(values) >= MIN_PRIMES_PER_CLASS
        and math.gcd(r, modulus) == 1
    ]

    populated.sort(
        key=lambda item: item[1],
        reverse=True,
    )

    print(
        f"populated invertible residue classes = "
        f"{len(populated)}"
    )

    if len(populated) < 2:
        raise RuntimeError(
            f"Too few residue classes for M={modulus}"
        )

    success = 0

    for trial in range(CLASS_TRIALS):
        a, b = choose_residue_pair(
            classes,
            modulus,
            rng,
        )

        available = min(
            len(classes[a]),
            len(classes[b]),
        )

        count = min(
            PAIRS_PER_CLASS,
            available,
        )

        if count < 4:
            continue

        pairs = build_forced_collision_family(
            primes,
            modulus,
            a,
            b,
            count,
        )

        if len(pairs) < 4:
            continue

        print()
        print(
            f"TRIAL {trial + 1}: "
            f"a={a}, b={b}, "
            f"a*b mod M={(a*b) % modulus}, "
            f"pairs={len(pairs)}"
        )

        analyze_collision_family(
            pairs,
            modulus,
        )

        conflict_test(
            pairs,
            modulus,
        )

        success += 1

    print()
    print(
        f"completed collision trials = "
        f"{success}/{CLASS_TRIALS}"
    )


# =============================================================================
# TWO-STAGE SQUARE-ROOT BRIDGE
# =============================================================================

def square_root_bridge_experiment(
    primes: Sequence[int],
    modulus: int,
    rng: random.Random,
) -> None:
    print()
    print("=" * 78)
    print(f"SQUARE-DIFFERENCE BRIDGE: M={modulus}")
    print("=" * 78)

    classes = prime_classes(
        primes,
        modulus,
    )

    a, b = choose_residue_pair(
        classes,
        modulus,
        rng,
    )

    count = min(
        PAIRS_PER_CLASS,
        len(classes[a]),
        len(classes[b]),
    )

    pairs = build_forced_collision_family(
        primes,
        modulus,
        a,
        b,
        count,
    )

    Nres = pairs[0].N % modulus

    print(
        f"forced N residue = {Nres}"
    )
    print(
        f"a={a}, b={b}"
    )

    d2_residues = {
        ((p.p - p.q) ** 2) % modulus
        for p in pairs
    }

    t_minus_n = {
        (
            (p.S * p.S - 3 * p.N - p.N)
            % modulus
        )
        for p in pairs
    }

    print(
        f"distinct (p-q)^2 mod M = "
        f"{len(d2_residues)}"
    )
    print(
        f"distinct T-N mod M      = "
        f"{len(t_minus_n)}"
    )

    if d2_residues != t_minus_n:
        raise ArithmeticError(
            "Square-difference sets disagree."
        )

    # Show explicit distinct roots of the same quadratic relation where
    # available. Here the roots are represented by ±(p-q) mod M.
    root_residues = {
        (p.p - p.q) % modulus
        for p in pairs
    }

    print(
        f"distinct (p-q) residues = "
        f"{len(root_residues)}"
    )

    print(
        "bridge identity status = PASS"
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 115")
    print("CONSTRUCTIVE CRT COLLISION / FORCED N-RESIDUE EQUIVALENCE")
    print("SAME N MOD M, DIFFERENT S / T / PAPER QUOTIENTS")
    print("CRT COLLISION -> SQUARE-DIFFERENCE BRIDGE")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    rng = random.Random(
        RNG_SEED
    )

    # -------------------------------------------------------------------------
    # 1. PRIME POPULATION
    # -------------------------------------------------------------------------
    print()
    print("1. PRIME POPULATION")
    print("-" * 78)

    start = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    elapsed = time.perf_counter() - start

    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = {elapsed:.6f}s"
    )

    if len(primes) < 100_000:
        raise ArithmeticError(
            "Unexpectedly small prime population."
        )

    # -------------------------------------------------------------------------
    # 2. CONSTRUCTIVE MODULUS TESTS
    # -------------------------------------------------------------------------
    print()
    print("2. CONSTRUCTIVE CRT COLLISION TESTS")
    print("-" * 78)

    for modulus in MODULI:
        modulus_experiment(
            primes,
            modulus,
            rng,
        )

    # -------------------------------------------------------------------------
    # 3. SQUARE DIFFERENCE BRIDGE
    # -------------------------------------------------------------------------
    print()
    print("3. SQUARE-DIFFERENCE / ROOT BRIDGE")
    print("-" * 78)

    square_root_bridge_experiment(
        primes,
        385,
        rng,
    )

    # -------------------------------------------------------------------------
    # 4. FINAL SUMMARY
    # -------------------------------------------------------------------------
    elapsed_total = (
        time.perf_counter()
        - total_start
    )

    print()
    print("=" * 78)
    print("4. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The decisive construction is:"
    )
    print()
    print(
        "    p_i == a (mod M)"
    )
    print(
        "    q_i == b (mod M)"
    )
    print()
    print(
        "which forces"
    )
    print()
    print(
        "    N_i = p_i q_i == a*b (mod M)"
    )
    print()
    print(
        "for every constructed target."
    )

    print()
    print(
        "If S, T, or Q_(k,l) still take multiple residues modulo M,"
    )
    print(
        "then no function of the single bounded residue N mod M can"
    )
    print(
        "recover those quantities."
    )

    print()
    print(
        "The square-root bridge checks independently that"
    )
    print()
    print(
        "    T-N = (p-q)^2"
    )
    print()
    print(
        "so the unresolved information is exactly a hidden modular"
    )
    print(
        "square-root / sign-selection problem."
    )

    print()
    print(
        "This is the key reason to compare the detector programme with"
    )
    print(
        "a^2 == b^2 (mod N): both routes ultimately require access to"
    )
    print(
        "information distinguishing the two quadratic roots."
    )

    print()
    print(
        "A particularly interesting outcome is:"
    )
    print(
        "    same N residue"
    )
    print(
        "    + different Q residue"
    )
    print(
        "    + different (p-q)^2 residue"
    )
    print(
        "    => detector quotient carries genuinely more information"
    )
    print(
        "       than bounded N residue alone."
    )

    print()
    print(
        f"total runtime = {elapsed_total:.6f}s"
    )
    print("=" * 78)
    print("EXPERIMENT 115 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(
            f"\nFATAL: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise
