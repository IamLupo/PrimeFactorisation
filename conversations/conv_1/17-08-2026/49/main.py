#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 116R
MEMORY-SAFE MULTIPLICATIVE CRT COLLISION SEARCH
QUADRATIC TRACE / NORM BRANCH CONSTRUCTION
SAME N MOD M, DIFFERENT S MOD M
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================

This version fixes Experiment 116's memory problem.

We do NOT build all O(phi(M)^2) residue pairs.

Instead, for a chosen product P and candidate traces S, solve

    x^2 - S*x + P = 0 (mod M)

using CRT over the prime-power factors of M.

Two different traces give two different decompositions

    (a,b), (c,d)

with

    ab == cd == P (mod M)

but

    a+b != c+d (mod M).

The resulting residue branches are then populated with actual primes.
"""

from __future__ import annotations

import math
import random
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple


# =============================================================================
# PARAMETERS
# =============================================================================

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

MODULI = (
    35,
    385,
    5005,
    55055,
)

TRIALS_PER_MODULUS = 20
PRIMES_PER_BRANCH = 10

RNG_SEED = 116


DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)


# =============================================================================
# DATA
# =============================================================================

@dataclass(frozen=True)
class PairRecord:
    p: int
    q: int

    @property
    def N(self) -> int:
        return self.p * self.q

    @property
    def S(self) -> int:
        return self.p + self.q

    @property
    def D(self) -> int:
        return self.p - self.q

    @property
    def T(self) -> int:
        return self.S * self.S - 3 * self.N


# =============================================================================
# SIEVE
# =============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if lo > hi:
        return []

    root = math.isqrt(hi)

    base = bytearray(b"\x01") * (root + 1)
    if root >= 0:
        base[0] = 0
    if root >= 1:
        base[1] = 0

    for p in range(2, math.isqrt(root) + 1):
        if base[p]:
            start = p * p
            for x in range(start, root + 1, p):
                base[x] = 0

    base_primes = [
        p for p in range(2, root + 1)
        if base[p]
    ]

    out: List[int] = []
    segment = 250_000

    for left in range(lo, hi + 1, segment):
        right = min(hi, left + segment - 1)
        flags = bytearray(b"\x01") * (right - left + 1)

        for p in base_primes:
            first = max(
                p * p,
                ((left + p - 1) // p) * p,
            )
            if first > right:
                continue

            for x in range(first, right + 1, p):
                flags[x - left] = 0

        out.extend(
            left + i
            for i, flag in enumerate(flags)
            if flag
        )

    return out


# =============================================================================
# DETECTORS
# =============================================================================

def raw_detector(k: int, ell: int, p: int, q: int) -> int:
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
    denominator = p + q + 1

    remainder = numerator % denominator

    if remainder != 0:
        raise ArithmeticError(
            f"Q_({k},{ell}) quotient failure: "
            f"remainder={remainder}"
        )

    return numerator // denominator


# =============================================================================
# MODULAR UTILITIES
# =============================================================================

def factor_prime_powers(n: int) -> List[int]:
    out: List[int] = []
    x = n
    p = 2

    while p * p <= x:
        if x % p == 0:
            power = 1
            while x % p == 0:
                x //= p
                power *= p
            out.append(power)

        p += 1 if p == 2 else 2

    if x > 1:
        out.append(x)

    return out


def crt_pair(a1: int, m1: int, a2: int, m2: int) -> int:
    """
    m1,m2 coprime.
    """
    inv = pow(m1, -1, m2)
    t = ((a2 - a1) * inv) % m2
    return (a1 + m1 * t) % (m1 * m2)


def crt_many(values: Sequence[int], moduli: Sequence[int]) -> int:
    if not values:
        raise ValueError("Empty CRT input.")

    x = values[0]
    m = moduli[0]

    for a, q in zip(values[1:], moduli[1:]):
        x = crt_pair(x, m, a, q)
        m *= q

    return x % m


# =============================================================================
# MODULAR SQUARE ROOTS
# =============================================================================

def tonelli_shanks(a: int, p: int) -> Optional[int]:
    """
    Return one square root of a mod odd prime p, or None.
    """
    a %= p

    if a == 0:
        return 0

    if p == 2:
        return a

    if pow(a, (p - 1) // 2, p) != 1:
        return None

    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    q = p - 1
    s = 0

    while q % 2 == 0:
        s += 1
        q //= 2

    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    c = pow(z, q, p)
    x = pow(a, (q + 1) // 2, p)
    t = pow(a, q, p)
    m = s

    while t != 1:
        i = 1
        tt = (t * t) % p

        while tt != 1:
            tt = (tt * tt) % p
            i += 1

            if i >= m:
                return None

        b = pow(c, 1 << (m - i - 1), p)

        x = (x * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return x


def sqrt_mod_prime_power(
    a: int,
    prime_power: int,
) -> List[int]:
    """
    Modular square roots for the prime powers occurring here.

    For this experiment all prime powers are odd and squarefree-modulus
    instances are the main case, but this routine also handles p^e
    by brute lifting from the prime root.
    """
    if prime_power <= 1:
        return [0]

    # Factor p^e.
    p = 2
    while p * p <= prime_power:
        if prime_power % p == 0:
            break
        p += 1

    if prime_power % p != 0:
        p = prime_power

    e = 0
    q = prime_power
    power = 1

    while q % p == 0:
        q //= p
        e += 1
        power *= p

    if q != 1:
        raise ArithmeticError(
            f"Unexpected prime-power decomposition: {prime_power}"
        )

    # For the moduli used here e=1. Keep a generic lifting path.
    if e == 1:
        root = tonelli_shanks(a, p)
        if root is None:
            return []

        r = root % p
        return sorted({r, (-r) % p})

    # Direct enumeration of roots is acceptable only for a small p^e.
    # Never invoke this for the large global modulus directly.
    if prime_power > 100_000:
        raise ArithmeticError(
            f"Large non-prime prime-power unsupported: {prime_power}"
        )

    return [
        x
        for x in range(prime_power)
        if (x * x - a) % prime_power == 0
    ]


# =============================================================================
# QUADRATIC BRANCH SOLVER
# =============================================================================

def quadratic_roots_mod(
    trace: int,
    product: int,
    modulus: int,
) -> List[Tuple[int, int]]:
    """
    Solve

        x^2 - trace*x + product == 0 (mod modulus)

    via discriminant

        Delta = trace^2 - 4*product.

    For each square root y of Delta:

        x = (trace +/- y) / 2.

    All moduli used are odd, so 2 is invertible.
    """
    trace %= modulus
    product %= modulus

    delta = (
        trace * trace
        - 4 * product
    ) % modulus

    pp = factor_prime_powers(modulus)

    root_lists: List[List[int]] = []

    for m in pp:
        local_delta = delta % m
        roots = sqrt_mod_prime_power(
            local_delta,
            m,
        )

        if not roots:
            return []

        root_lists.append(roots)

    # Combine local roots by CRT.
    roots_mod_M = [0]

    current_mod = 1

    for local_mod, local_roots in zip(pp, root_lists):
        new_roots: List[int] = []

        for x0 in roots_mod_M:
            for r in local_roots:
                if current_mod == 1:
                    x = r % local_mod
                else:
                    x = crt_pair(
                        x0,
                        current_mod,
                        r,
                        local_mod,
                    )

                new_roots.append(x)

        roots_mod_M = sorted(
            set(new_roots)
        )

        current_mod *= local_mod

    inv2 = pow(2, -1, modulus)

    result = set()

    for root in roots_mod_M:
        x1 = ((trace + root) * inv2) % modulus
        x2 = ((trace - root) * inv2) % modulus

        if (
            x1 * x1
            - trace * x1
            + product
        ) % modulus != 0:
            continue

        if (
            x2 * x2
            - trace * x2
            + product
        ) % modulus != 0:
            continue

        result.add(
            tuple(sorted((x1, x2)))
        )

    return sorted(result)


# =============================================================================
# SEARCH WITHOUT LARGE MEMORY
# =============================================================================

def find_collision(
    modulus: int,
    rng: random.Random,
    attempts: int = 30_000,
) -> Optional[
    Tuple[int, int, int, int, int, int]
]:
    """
    Return

        (P, S1, a, b, S2, c, d)

    conceptually, but packed as a 6-tuple:

        (P, S1, a, b, S2, c, d)

    where

        ab == cd == P mod M
        a+b == S1 mod M
        c+d == S2 mod M
        S1 != S2.

    We search one product and several traces rather than all pairs.
    """
    residues = [
        r
        for r in range(1, modulus)
        if math.gcd(r, modulus) == 1
    ]

    if len(residues) < 8:
        return None

    for _ in range(attempts):
        P = rng.choice(residues)

        # Pick traces randomly. We do not build a product->pair table.
        trace_seen: Dict[int, Tuple[int, int]] = {}

        for _ in range(50):
            S = rng.randrange(modulus)

            if S in trace_seen:
                continue

            roots = quadratic_roots_mod(
                S,
                P,
                modulus,
            )

            # A valid factorization requires a nontrivial pair.
            for a, b in roots:
                if (
                    math.gcd(a, modulus) != 1
                    or math.gcd(b, modulus) != 1
                ):
                    continue

                if a == b:
                    continue

                trace_seen[S] = (a, b)
                break

        if len(trace_seen) < 2:
            continue

        traces = list(trace_seen.items())
        rng.shuffle(traces)

        for i in range(len(traces)):
            S1, pair1 = traces[i]

            for j in range(i + 1, len(traces)):
                S2, pair2 = traces[j]

                if S1 == S2:
                    continue

                a, b = pair1
                c, d = pair2

                disc1 = (
                    (a - b) * (a - b)
                ) % modulus

                disc2 = (
                    (c - d) * (c - d)
                ) % modulus

                if disc1 == disc2:
                    continue

                return (
                    P,
                    S1,
                    a,
                    b,
                    S2,
                    c,
                    d,
                )

    return None


# =============================================================================
# PRIME RESIDUE CLASSES
# =============================================================================

def build_prime_classes(
    primes: Sequence[int],
    modulus: int,
) -> Dict[int, List[int]]:
    classes: Dict[int, List[int]] = defaultdict(list)

    for p in primes:
        classes[p % modulus].append(p)

    return dict(classes)


def select_primes(
    classes: Dict[int, List[int]],
    residue: int,
    count: int,
) -> Optional[Tuple[int, ...]]:
    values = classes.get(residue)

    if values is None or len(values) < count:
        return None

    return tuple(values[:count])


# =============================================================================
# CONSTRUCT BRANCHES
# =============================================================================

def construct_pairs(
    classes: Dict[int, List[int]],
    a: int,
    b: int,
    c: int,
    d: int,
) -> List[Tuple[PairRecord, PairRecord]]:
    pa = select_primes(
        classes,
        a,
        PRIMES_PER_BRANCH,
    )
    pb = select_primes(
        classes,
        b,
        PRIMES_PER_BRANCH,
    )
    pc = select_primes(
        classes,
        c,
        PRIMES_PER_BRANCH,
    )
    pd = select_primes(
        classes,
        d,
        PRIMES_PER_BRANCH,
    )

    if None in (pa, pb, pc, pd):
        return []

    assert pa is not None
    assert pb is not None
    assert pc is not None
    assert pd is not None

    result = []

    for i in range(PRIMES_PER_BRANCH):
        left = PairRecord(
            pa[i],
            pb[-1 - i],
        )

        right = PairRecord(
            pc[i],
            pd[-1 - i],
        )

        if left.p == left.q:
            continue

        if right.p == right.q:
            continue

        result.append(
            (left, right)
        )

    return result


# =============================================================================
# VALIDATION
# =============================================================================

def verify_pair_collision(
    pairs: Sequence[Tuple[PairRecord, PairRecord]],
    modulus: int,
) -> Dict[str, int]:
    stats = {
        "N_same": 0,
        "S_different": 0,
        "D2_different": 0,
        "T_different": 0,
    }

    for left, right in pairs:
        if left.N % modulus == right.N % modulus:
            stats["N_same"] += 1

        if left.S % modulus != right.S % modulus:
            stats["S_different"] += 1

        if (
            (left.D * left.D) % modulus
            != (right.D * right.D) % modulus
        ):
            stats["D2_different"] += 1

        if left.T % modulus != right.T % modulus:
            stats["T_different"] += 1

    return stats


def detector_profile(
    pairs: Sequence[Tuple[PairRecord, PairRecord]],
    modulus: int,
) -> Dict[Tuple[int, int], Tuple[int, int]]:
    result = {}

    for detector in DETECTORS:
        k, ell = detector

        same = 0
        diff = 0

        for left, right in pairs:
            ql = (
                quotient_value(k, ell, left.p, left.q)
                % modulus
            )
            qr = (
                quotient_value(k, ell, right.p, right.q)
                % modulus
            )

            if ql == qr:
                same += 1
            else:
                diff += 1

        result[detector] = (same, diff)

    return result


# =============================================================================
# SINGLE MODULUS
# =============================================================================

def run_modulus(
    primes: Sequence[int],
    modulus: int,
    rng: random.Random,
) -> None:
    print()
    print("=" * 78)
    print(f"MODULUS M={modulus}")
    print("=" * 78)

    classes = build_prime_classes(
        primes,
        modulus,
    )

    populated = sum(
        len(v) >= PRIMES_PER_BRANCH
        for v in classes.values()
    )

    print(
        f"populated residue classes = {populated}"
    )

    successes = 0
    failed = 0

    start = time.perf_counter()

    for trial in range(1, TRIALS_PER_MODULUS + 1):
        collision = find_collision(
            modulus,
            rng,
        )

        if collision is None:
            print(
                f"trial {trial}: no collision found"
            )
            failed += 1
            continue

        P, S1, a, b, S2, c, d = collision

        pairs = construct_pairs(
            classes,
            a,
            b,
            c,
            d,
        )

        if len(pairs) < 4:
            failed += 1
            continue

        stats = verify_pair_collision(
            pairs,
            modulus,
        )

        if stats["N_same"] != len(pairs):
            raise ArithmeticError(
                "Constructed branches do not have identical N residue."
            )

        if stats["S_different"] != len(pairs):
            raise ArithmeticError(
                "Constructed branches do not have distinct S residues."
            )

        if stats["D2_different"] != len(pairs):
            raise ArithmeticError(
                "Constructed branches do not have distinct "
                "(p-q)^2 residues."
            )

        if stats["T_different"] != len(pairs):
            raise ArithmeticError(
                "Constructed branches do not have distinct T residues."
            )

        profile = detector_profile(
            pairs,
            modulus,
        )

        separated = sum(
            diff > 0
            for same, diff in profile.values()
        )

        print()
        print(
            f"TRIAL {trial}: "
            f"(a,b)=({a},{b}) "
            f"(c,d)=({c},{d})"
        )

        print(
            f"  product P        = {P}"
        )
        print(
            f"  trace S1         = {S1}"
        )
        print(
            f"  trace S2         = {S2}"
        )
        print(
            f"  pair count       = {len(pairs)}"
        )

        print(
            f"  N same           = {stats['N_same']}/{len(pairs)}"
        )
        print(
            f"  S different      = "
            f"{stats['S_different']}/{len(pairs)}"
        )
        print(
            f"  D^2 different    = "
            f"{stats['D2_different']}/{len(pairs)}"
        )
        print(
            f"  T different      = "
            f"{stats['T_different']}/{len(pairs)}"
        )

        print()
        print("  DETECTOR SEPARATION")
        print("  " + "-" * 72)

        for detector in DETECTORS:
            same, diff = profile[detector]

            print(
                f"  Q_{detector}: "
                f"same={same:2d} "
                f"different={diff:2d}"
            )

        if separated > 0:
            successes += 1
        else:
            failed += 1

        # One explicit witness.
        left, right = pairs[0]

        print()
        print("  WITNESS")
        print("  " + "-" * 72)

        print(
            f"  left : p={left.p} q={left.q} "
            f"N={left.N % modulus} "
            f"S={left.S % modulus} "
            f"D2={(left.D*left.D) % modulus}"
        )

        print(
            f"  right: p={right.p} q={right.q} "
            f"N={right.N % modulus} "
            f"S={right.S % modulus} "
            f"D2={(right.D*right.D) % modulus}"
        )

        for detector in DETECTORS:
            k, ell = detector

            ql = (
                quotient_value(k, ell, left.p, left.q)
                % modulus
            )

            qr = (
                quotient_value(k, ell, right.p, right.q)
                % modulus
            )

            print(
                f"  Q_{detector} = {ql} vs {qr}"
            )

    elapsed = time.perf_counter() - start

    print()
    print(
        f"successful detector-separation trials = "
        f"{successes}"
    )
    print(
        f"failed trials = {failed}"
    )
    print(
        f"elapsed = {elapsed:.6f}s"
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 116R")
    print("MEMORY-SAFE MULTIPLICATIVE CRT COLLISION")
    print("SAME N MOD M / DIFFERENT TRACE BRANCHES")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    rng = random.Random(RNG_SEED)

    # -------------------------------------------------------------------------
    # PRIME POPULATION
    # -------------------------------------------------------------------------
    print()
    print("1. PRIME POPULATION")
    print("-" * 78)

    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # -------------------------------------------------------------------------
    # MODULUS TESTS
    # -------------------------------------------------------------------------
    print()
    print("2. MEMORY-SAFE COLLISION SEARCH")
    print("-" * 78)

    for modulus in MODULI:
        run_modulus(
            primes,
            modulus,
            rng,
        )

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("3. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "Experiment 116R constructs multiplicative collisions "
        "without enumerating the full residue-pair space."
    )

    print()
    print(
        "The target configuration is:"
    )
    print(
        "    N1 == N2 (mod M)"
    )
    print(
        "    S1 != S2 (mod M)"
    )
    print(
        "    (p1-q1)^2 != (p2-q2)^2 (mod M)"
    )

    print()
    print(
        "The critical test is whether the paper quotients "
        "separate these branches."
    )

    print()
    print(
        "This is the cleaner continuation of Experiment 115:"
    )
    print(
        "we fix the norm/product residue without fixing the "
        "individual factor residues."
    )

    print()
    print(
        "The search uses only modular quadratic equations and CRT."
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - total_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 116R COMPLETE")
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