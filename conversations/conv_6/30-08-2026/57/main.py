#!/usr/bin/env python3

"""
============================================================================================
DETERMINISTIC SMOOTH n+x -> COFACTOR -> K RECURSION EXPERIMENT
============================================================================================

For odd n choose x deterministically so that:

    n + x = S * m

where S is smooth.

The experiment explicitly treats S as already factored.

Only the remaining cofactor m must be factored.

From the complete factorization of n+x we generate all divisor
pairs:

    n+x = p_x * q_x

and for each pair calculate:

    k_x = floor(p_x / r1)
    l_x = floor(q_x / r2)

    K_x = k_x * l_x

    T_x = floor((n+x)/(r1*r2))
    E_x = T_x - K_x

For the original n:

    n = p*q

    k = floor(p/r1)
    l = floor(q/r2)

    K = k*l

    T = floor(n/(r1*r2))
    E = T-K

The experiment is specifically interested in whether a
deterministically constructed n+x produces quotient information
that can constrain the original K.

IMPORTANT:

This script does NOT search over x.

x is always:

    x = (-n) mod S

No original p interval is enumerated.

The auxiliary number n+x may have more than two prime factors.
That is intentional.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

N_ANCHORS = 30

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

# Smooth numbers.
#
# Because n is odd and S is even, x is odd.
#
SMOOTH_S = [
    30,
    210,
    2310,
    2730,
    30030,
]

PRINT_EVERY = 5


# ==========================================================================================
# DATA STRUCTURES
# ==========================================================================================

@dataclass
class Anchor:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q


@dataclass
class AuxPair:
    S: int
    x: int
    nx: int
    cofactor: int

    px: int
    qx: int

    kx: int
    lx: int

    Kx: int
    Tx: int
    Ex: int

    factors: dict[int, int]


# ==========================================================================================
# PRIME / FACTOR UTILITIES
# ==========================================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False

        d += 2

    return True


def prime_pool(lo: int, hi: int) -> list[int]:
    out = []

    for n in range(lo, hi + 1):
        if is_prime(n):
            out.append(n)

    return out


def factor_integer(n: int) -> dict[int, int]:
    """
    Exact trial-division factorization for the relatively small
    experimental values.
    """
    n = int(n)

    if n < 2:
        return {}

    factors: dict[int, int] = {}

    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2

    d = 3

    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d

        d += 2

    if n > 1:
        factors[n] = factors.get(n, 0) + 1

    return factors


def multiply_factorization(
    factors: dict[int, int],
) -> int:

    value = 1

    for p, e in factors.items():
        value *= p ** e

    return value


def divisor_pairs_from_factors(
    factors: dict[int, int],
) -> list[tuple[int, int]]:

    value = multiply_factorization(factors)

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(divisors)

        power = 1

        for _ in range(exponent):

            power *= prime

            for d in old:
                divisors.append(d * power)

    divisors = sorted(set(divisors))

    pairs = []

    for d in divisors:

        if d * d > value:
            break

        if value % d == 0:
            pairs.append(
                (d, value // d)
            )

    return pairs


# ==========================================================================================
# ANCHORS
# ==========================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    anchors: list[Anchor] = []
    seen: set[int] = set()

    while len(anchors) < count:

        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        anchors.append(
            Anchor(
                p=p,
                q=q,
            )
        )

    anchors.sort(
        key=lambda a: a.n
    )

    return anchors


# ==========================================================================================
# CLOSE MODULUS PAIRS
# ==========================================================================================

def build_close_pairs(
    mod_primes: list[int],
    ratio: float,
) -> list[tuple[int, int]]:

    pairs = []

    for i, r1 in enumerate(mod_primes):

        for r2 in mod_primes[i + 1:]:

            if (r2 - r1) <= ratio * r1:
                pairs.append(
                    (r1, r2)
                )

    return pairs


# ==========================================================================================
# DETERMINISTIC x
# ==========================================================================================

def deterministic_x(
    n: int,
    S: int,
) -> int:
    """
    Smallest positive x such that:

        n+x == 0 (mod S)
    """
    x = (-n) % S

    if x == 0:
        x = S

    return x


# ==========================================================================================
# ORIGINAL COORDINATES
# ==========================================================================================

def original_coordinates(
    p: int,
    q: int,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int, int]:

    R = r1 * r2

    k = p // r1
    l = q // r2

    K = k * l

    T = n_div = (p * q) // R

    E = T - K

    return k, l, K, T, E


# ==========================================================================================
# AUXILIARY FACTORIZATION
# ==========================================================================================

def build_auxiliary_factorization(
    n: int,
    S: int,
) -> tuple[int, int, dict[int, int]]:

    x = deterministic_x(
        n,
        S,
    )

    nx = n + x

    if nx % S != 0:
        raise AssertionError(
            "Construction failed: n+x is not divisible by S."
        )

    cofactor = nx // S

    # Factor only the unknown cofactor.
    cofactor_factors = factor_integer(
        cofactor
    )

    # Factorization of S is known.
    S_factors = factor_integer(
        S
    )

    full_factors = dict(
        S_factors
    )

    for prime, exponent in cofactor_factors.items():
        full_factors[prime] = (
            full_factors.get(prime, 0)
            + exponent
        )

    reconstructed = multiply_factorization(
        full_factors
    )

    if reconstructed != nx:
        raise AssertionError(
            "Full factorization does not reconstruct n+x."
        )

    return x, cofactor, full_factors


# ==========================================================================================
# GENERATE AUXILIARY QUOTIENT PAIRS
# ==========================================================================================

def auxiliary_quotient_pairs(
    n: int,
    r1: int,
    r2: int,
    S: int,
) -> list[AuxPair]:

    x, cofactor, factors = (
        build_auxiliary_factorization(
            n,
            S,
        )
    )

    nx = n + x

    R = r1 * r2

    pairs = divisor_pairs_from_factors(
        factors
    )

    out = []

    for px, qx in pairs:

        # Only use factors in the experimental
        # factor range.
        if px < FACTOR_MIN:
            continue

        if px > FACTOR_MAX:
            continue

        if qx < FACTOR_MIN:
            continue

        if qx > FACTOR_MAX:
            continue

        if px > qx:
            continue

        if px * qx != nx:
            continue

        kx = px // r1
        lx = qx // r2

        Kx = kx * lx

        Tx = nx // R

        Ex = Tx - Kx

        out.append(
            AuxPair(
                S=S,
                x=x,
                nx=nx,
                cofactor=cofactor,
                px=px,
                qx=qx,
                kx=kx,
                lx=lx,
                Kx=Kx,
                Tx=Tx,
                Ex=Ex,
                factors=factors,
            )
        )

    return out


# ==========================================================================================
# K FACTORIZATION
# ==========================================================================================

def K_factorization(
    K: int,
) -> tuple[dict[int, int], list[tuple[int, int]]]:

    if K <= 0:
        return {}, []

    factors = factor_integer(K)

    pairs = divisor_pairs_from_factors(
        factors
    )

    return factors, pairs


# ==========================================================================================
# TRAJECTORY IDENTITY
# ==========================================================================================

def verify_trajectory_identity(
    n: int,
    x: int,
    r1: int,
    r2: int,
    K: int,
    E: int,
    Kx: int,
    Ex: int,
) -> None:

    R = r1 * r2

    T = n // R

    Tx = (n + x) // R

    left = Kx - K

    right = (
        (Tx - T)
        - (Ex - E)
    )

    if left != right:
        raise AssertionError(
            "K/E trajectory identity failed."
        )


# ==========================================================================================
# RUN ONE ANCHOR
# ==========================================================================================

def run_anchor(
    anchor: Anchor,
    modulus_pair: tuple[int, int],
) -> dict:

    p = anchor.p
    q = anchor.q
    n = anchor.n

    r1, r2 = modulus_pair

    k, l, K, T, E = original_coordinates(
        p,
        q,
        r1,
        r2,
    )

    result = {
        "n": n,
        "p": p,
        "q": q,
        "r1": r1,
        "r2": r2,
        "k": k,
        "l": l,
        "K": K,
        "T": T,
        "E": E,
        "aux": [],
    }

    for S in SMOOTH_S:

        aux_pairs = auxiliary_quotient_pairs(
            n=n,
            r1=r1,
            r2=r2,
            S=S,
        )

        # A single n+x can have multiple relevant
        # factor pairs.
        for aux in aux_pairs:

            verify_trajectory_identity(
                n=n,
                x=aux.x,
                r1=r1,
                r2=r2,
                K=K,
                E=E,
                Kx=aux.Kx,
                Ex=aux.Ex,
            )

            Kx_factors, Kx_pairs = (
                K_factorization(
                    aux.Kx
                )
            )

            result["aux"].append(
                {
                    "S": aux.S,
                    "x": aux.x,
                    "nx": aux.nx,
                    "cofactor": aux.cofactor,
                    "full_factors": aux.factors,
                    "px": aux.px,
                    "qx": aux.qx,
                    "kx": aux.kx,
                    "lx": aux.lx,
                    "Kx": aux.Kx,
                    "Tx": aux.Tx,
                    "Ex": aux.Ex,
                    "Kx_factors": Kx_factors,
                    "Kx_pairs": Kx_pairs,
                }
            )

    return result


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    start = time.perf_counter()

    print("=" * 100)
    print(
        "DETERMINISTIC SMOOTH n+x -> COFACTOR -> "
        "K RECURSION EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"N anchors                 = {N_ANCHORS}"
    )

    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )

    print(
        f"modulus range             = "
        f"{MOD_MIN:,} - {MOD_MAX:,}"
    )

    print(
        f"close ratio               = "
        f"{CLOSE_RATIO:.2%}"
    )

    print(
        f"seed                      = "
        f"{SEED:,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # PRIME POOLS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "BUILDING PRIME POOLS"
    )
    print("=" * 100)

    factor_primes = prime_pool(
        FACTOR_MIN,
        FACTOR_MAX,
    )

    modulus_primes = prime_pool(
        MOD_MIN,
        MOD_MAX,
    )

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # MODULUS PAIRS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "BUILDING CLOSE MODULUS PAIRS"
    )
    print("=" * 100)

    close_pairs = build_close_pairs(
        modulus_primes,
        CLOSE_RATIO,
    )

    print(
        f"close modulus pairs       = "
        f"{len(close_pairs):,}"
    )

    print(
        f"smooth S values           = "
        f"{SMOOTH_S}"
    )

    print()

    if not close_pairs:
        raise RuntimeError(
            "No close modulus pairs found."
        )

    # --------------------------------------------------------------------------------------
    # ANCHORS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "BUILDING ANCHORS"
    )
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print(
        f"actual anchors            = "
        f"{len(anchors):,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # MODULUS ASSIGNMENT
    # --------------------------------------------------------------------------------------

    rng = random.Random(
        SEED
    )

    selected_pairs = [
        rng.choice(close_pairs)
        for _ in anchors
    ]

    # --------------------------------------------------------------------------------------
    # EXPERIMENT
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "RUNNING DETERMINISTIC AUXILIARY FACTORIZATIONS"
    )
    print("=" * 100)

    results = []

    for index, (anchor, pair) in enumerate(
        zip(
            anchors,
            selected_pairs,
        ),
        start=1,
    ):

        result = run_anchor(
            anchor,
            pair,
        )

        results.append(result)

        if (
            index == 1
            or index == len(anchors)
            or index % PRINT_EVERY == 0
        ):
            print(
                f"anchor {index:3d}/"
                f"{len(anchors):3d}"
            )

    print()

    # --------------------------------------------------------------------------------------
    # RAW RESULTS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "RAW RECURSIVE TRAJECTORIES"
    )
    print("=" * 100)

    for result in results:

        print(
            f"n={result['n']:,} "
            f"p={result['p']:,} "
            f"q={result['q']:,} "
            f"mods=("
            f"{result['r1']},"
            f"{result['r2']})"
        )

        print(
            f"    original: "
            f"k={result['k']} "
            f"l={result['l']} "
            f"K={result['K']} "
            f"T={result['T']} "
            f"E={result['E']}"
        )

        if not result["aux"]:
            print(
                "    no valid auxiliary factor pair "
                "in the selected factor range"
            )

        for aux in result["aux"]:

            K = result["K"]
            E = result["E"]
            T = result["T"]

            delta_K = aux["Kx"] - K
            delta_T = aux["Tx"] - T
            delta_E = aux["Ex"] - E

            print()

            print(
                f"    S={aux['S']:,} "
                f"x={aux['x']:,}"
            )

            print(
                f"        n+x={aux['nx']:,}"
            )

            print(
                f"        (n+x)/S="
                f"{aux['cofactor']:,}"
            )

            print(
                f"        full factorization="
                f"{aux['full_factors']}"
            )

            print(
                f"        auxiliary factors="
                f"({aux['px']:,},"
                f"{aux['qx']:,})"
            )

            print(
                f"        (kx,lx)=("
                f"{aux['kx']},"
                f"{aux['lx']})"
            )

            print(
                f"        Kx={aux['Kx']}"
            )

            print(
                f"        Tx={aux['Tx']}"
            )

            print(
                f"        Ex={aux['Ex']}"
            )

            print(
                f"        factor(Kx)="
                f"{aux['Kx_factors']}"
            )

            print(
                f"        divisor pairs(Kx)="
                f"{aux['Kx_pairs']}"
            )

            print(
                "        trajectory:"
            )

            print(
                f"            delta_T = "
                f"{delta_T}"
            )

            print(
                f"            delta_E = "
                f"{delta_E}"
            )

            print(
                f"            delta_K = "
                f"{delta_K}"
            )

            print(
                f"            "
                f"(delta_T-delta_E) = "
                f"{delta_T - delta_E}"
            )

            print(
                f"            identity = "
                f"{delta_K == delta_T - delta_E}"
            )

        print()

    # --------------------------------------------------------------------------------------
    # DETERMINISTIC CONSTRUCTION VALIDATION
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "CONSTRUCTION VALIDATION"
    )
    print("=" * 100)

    for result in results:

        n = result["n"]

        for aux in result["aux"]:

            S = aux["S"]
            x = aux["x"]
            nx = aux["nx"]

            if nx != n + x:
                raise AssertionError(
                    "n+x mismatch."
                )

            if nx % S != 0:
                raise AssertionError(
                    "n+x is not divisible by S."
                )

            if aux["cofactor"] != nx // S:
                raise AssertionError(
                    "Cofactor mismatch."
                )

    print(
        "all deterministic x constructions verified"
    )

    print(
        "all n+x factorizations reconstruct exactly"
    )

    print(
        "all K/E trajectory identities verified"
    )

    print(
        "no x search performed"
    )

    print(
        "no original p interval enumeration performed"
    )

    # --------------------------------------------------------------------------------------
    # RECURSIVE DESCRIPTION
    # --------------------------------------------------------------------------------------

    print()

    print("=" * 100)
    print(
        "RECURSIVE STRUCTURE"
    )
    print("=" * 100)

    print(
        "Level 0:"
    )
    print(
        "    n"
    )

    print(
        "Level 1:"
    )
    print(
        "    x = (-n) mod S"
    )
    print(
        "    n+x = S*m"
    )

    print(
        "Level 2:"
    )
    print(
        "    factor(m)"
    )
    print(
        "    reconstruct divisors of n+x"
    )

    print(
        "Level 3:"
    )
    print(
        "    p_x = a_x + k_x*r1"
    )
    print(
        "    q_x = b_x + l_x*r2"
    )

    print(
        "Level 4:"
    )
    print(
        "    K_x = k_x*l_x"
    )

    print(
        "Level 5:"
    )
    print(
        "    factor(K_x)"
    )

    print()

    print(
        "The unresolved bridge remains:"
    )

    print(
        "    K_x"
        "  --->"
        "  information about original K"
    )

    print()

    # --------------------------------------------------------------------------------------
    # TIMING
    # --------------------------------------------------------------------------------------

    elapsed = time.perf_counter() - start

    print("=" * 100)
    print(
        "TIMING"
    )
    print("=" * 100)

    print(
        f"total runtime              = "
        f"{elapsed:.3f} s"
    )

    print()

    print("=" * 100)
    print(
        "EXPERIMENT COMPLETE"
    )
    print("=" * 100)


if __name__ == "__main__":
    run()
