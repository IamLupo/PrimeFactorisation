#!/usr/bin/env python3

"""
============================================================================================
COMPACT DETERMINISTIC n+x -> Kx -> (k,l) RECURSIVE EXPERIMENT
============================================================================================

For each original semiprime:

    n = p*q

construct deterministic smooth neighbours:

    x = (-n) mod S

so that:

    n+x = S*m.

The known smooth part S is removed first.

The remaining cofactor m is factored.

Every valid factor pair of n+x gives:

    kx = floor(px/r1)
    lx = floor(qx/r2)

    Kx = kx*lx

    Tx = floor((n+x)/(r1*r2))
    Ex = Tx-Kx.

Original quantities:

    k = floor(p/r1)
    l = floor(q/r2)

    K = k*l

    T = floor(n/(r1*r2))
    E = T-K.

The experiment tests the recursive bridge:

    n+x
      ->
    factor(n+x)
      ->
    (kx,lx)
      ->
    Kx
      ->
    factor(Kx).

The output is intentionally compact.

No search over x is performed.
No original p interval is enumerated.
"""

from __future__ import annotations

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

# Use only a few S values so the experiment remains quick.
SMOOTH_S = (
    210,
    2310,
    30030,
)

# Number of compact examples printed at start/end.
SHOW_START = 3
SHOW_END = 3


# ==========================================================================================
# DATA
# ==========================================================================================

@dataclass
class Anchor:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q


@dataclass
class Aux:
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

    number_of_Kx_pairs: int


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
    return [
        n
        for n in range(lo, hi + 1)
        if is_prime(n)
    ]


def factor_integer(n: int) -> dict[int, int]:
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


def divisor_count_from_factors(
    factors: dict[int, int],
) -> int:
    count = 1

    for exponent in factors.values():
        count *= exponent + 1

    return count


def divisor_pairs(
    n: int,
) -> list[tuple[int, int]]:

    if n <= 0:
        return []

    factors = factor_integer(n)

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
        if d * d > n:
            break

        if n % d == 0:
            pairs.append((d, n // d))

    return pairs


# ==========================================================================================
# ANCHORS
# ==========================================================================================

def build_anchors(
    primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    seen: set[int] = set()
    out: list[Anchor] = []

    while len(out) < count:

        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        out.append(
            Anchor(p, q)
        )

    out.sort(key=lambda a: a.n)

    return out


# ==========================================================================================
# CLOSE MODULUS PAIRS
# ==========================================================================================

def build_close_pairs(
    primes: list[int],
    ratio: float,
) -> list[tuple[int, int]]:

    pairs = []

    for i, r1 in enumerate(primes):

        for r2 in primes[i + 1:]:

            if r2 - r1 <= ratio * r1:
                pairs.append((r1, r2))

    return pairs


# ==========================================================================================
# ORIGINAL COORDINATES
# ==========================================================================================

def original_state(
    anchor: Anchor,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int, int]:

    n = anchor.n
    R = r1 * r2

    k = anchor.p // r1
    l = anchor.q // r2

    K = k * l

    T = n // R

    E = T - K

    return k, l, K, T, E


# ==========================================================================================
# DETERMINISTIC x
# ==========================================================================================

def deterministic_x(
    n: int,
    S: int,
) -> int:

    x = (-n) % S

    if x == 0:
        x = S

    return x


# ==========================================================================================
# AUXILIARY FACTORIZATION
# ==========================================================================================

def build_auxiliaries(
    n: int,
    r1: int,
    r2: int,
    S: int,
) -> list[Aux]:

    x = deterministic_x(n, S)
    nx = n + x

    if nx % S != 0:
        raise AssertionError(
            "n+x is not divisible by S"
        )

    cofactor = nx // S

    # Factorization of S is known.
    sf = factor_integer(S)

    # Factor only the remaining cofactor.
    mf = factor_integer(cofactor)

    full = dict(sf)

    for p, e in mf.items():
        full[p] = full.get(p, 0) + e

    reconstructed = 1

    for p, e in full.items():
        reconstructed *= p ** e

    if reconstructed != nx:
        raise AssertionError(
            "Factorization reconstruction failure"
        )

    factor_pairs = divisor_pairs(nx)

    out = []

    R = r1 * r2
    Tx = nx // R

    for px, qx in factor_pairs:

        if px < FACTOR_MIN:
            continue

        if px > FACTOR_MAX:
            continue

        if qx < FACTOR_MIN:
            continue

        if qx > FACTOR_MAX:
            continue

        if px * qx != nx:
            continue

        kx = px // r1
        lx = qx // r2

        Kx = kx * lx
        Ex = Tx - Kx

        Kx_pairs = divisor_pairs(Kx)

        out.append(
            Aux(
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
                number_of_Kx_pairs=len(Kx_pairs),
            )
        )

    return out


# ==========================================================================================
# IDENTITY
# ==========================================================================================

def verify_identity(
    n: int,
    r1: int,
    r2: int,
    K: int,
    E: int,
    aux: Aux,
) -> None:

    R = r1 * r2

    T = n // R

    delta_k = aux.Kx - K

    delta_t = aux.Tx - T

    delta_e = aux.Ex - E

    if delta_k != delta_t - delta_e:
        raise AssertionError(
            "K trajectory identity failed"
        )


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    start = time.perf_counter()

    print("=" * 100)
    print(
        "COMPACT DETERMINISTIC n+x -> Kx -> (k,l) "
        "RECURSIVE EXPERIMENT"
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
    print(
        f"S values                  = "
        f"{list(SMOOTH_S)}"
    )
    print()

    # --------------------------------------------------------------------------------------
    # PRIME POOLS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOLS")
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
    # MODULI
    # --------------------------------------------------------------------------------------

    close_pairs = build_close_pairs(
        modulus_primes,
        CLOSE_RATIO,
    )

    print("=" * 100)
    print("BUILDING CLOSE MODULUS PAIRS")
    print("=" * 100)

    print(
        f"close modulus pairs       = "
        f"{len(close_pairs):,}"
    )

    print()

    if not close_pairs:
        raise RuntimeError(
            "No close modulus pairs."
        )

    # --------------------------------------------------------------------------------------
    # ANCHORS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING ANCHORS")
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

    rng = random.Random(SEED)

    selected_pairs = [
        rng.choice(close_pairs)
        for _ in anchors
    ]

    # --------------------------------------------------------------------------------------
    # RUN
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("RUNNING RECURSIVE n+x CONSTRUCTION")
    print("=" * 100)

    all_records: list[dict] = []

    successful_aux = 0

    identity_checks = 0

    for index, (anchor, modulus_pair) in enumerate(
        zip(anchors, selected_pairs),
        start=1,
    ):

        r1, r2 = modulus_pair

        k, l, K, T, E = original_state(
            anchor,
            r1,
            r2,
        )

        auxiliaries = []

        for S in SMOOTH_S:

            aux_list = build_auxiliaries(
                anchor.n,
                r1,
                r2,
                S,
            )

            for aux in aux_list:

                verify_identity(
                    anchor.n,
                    r1,
                    r2,
                    K,
                    E,
                    aux,
                )

                identity_checks += 1

                auxiliaries.append(aux)

                successful_aux += 1

        all_records.append(
            {
                "anchor": anchor,
                "r1": r1,
                "r2": r2,
                "k": k,
                "l": l,
                "K": K,
                "T": T,
                "E": E,
                "aux": auxiliaries,
            }
        )

        if (
            index == 1
            or index == N_ANCHORS
            or index % 5 == 0
        ):
            print(
                f"anchor {index:3d}/{N_ANCHORS:3d}"
            )

    # --------------------------------------------------------------------------------------
    # COMPACT START RECORDS
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        f"FIRST {min(SHOW_START, len(all_records))} "
        "ANCHORS"
    )
    print("=" * 100)

    for record in all_records[:SHOW_START]:

        anchor = record["anchor"]

        print()
        print(
            f"n={anchor.n:,}"
        )

        print(
            f"    original "
            f"(p,q)=({anchor.p:,},{anchor.q:,}) "
            f"mods=({record['r1']},{record['r2']})"
        )

        print(
            f"    (k,l)=({record['k']},{record['l']}) "
            f"K={record['K']} "
            f"T={record['T']} "
            f"E={record['E']}"
        )

        for aux in record["aux"][:4]:

            print(
                f"    S={aux.S:,} "
                f"x={aux.x:,} "
                f"n+x={aux.nx:,} "
                f"factors=({aux.px:,},{aux.qx:,})"
            )

            print(
                f"        "
                f"(kx,lx)=({aux.kx},{aux.lx}) "
                f"Kx={aux.Kx} "
                f"Tx={aux.Tx} "
                f"Ex={aux.Ex} "
                f"KxPairs={aux.number_of_Kx_pairs}"
            )

            print(
                f"        "
                f"dK={aux.Kx-record['K']} "
                f"dE={aux.Ex-record['E']} "
                f"dT={aux.Tx-record['T']}"
            )

    # --------------------------------------------------------------------------------------
    # COMPACT END RECORDS
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        f"LAST {min(SHOW_END, len(all_records))} "
        "ANCHORS"
    )
    print("=" * 100)

    for record in all_records[-SHOW_END:]:

        anchor = record["anchor"]

        print()
        print(
            f"n={anchor.n:,}"
        )

        print(
            f"    original "
            f"(p,q)=({anchor.p:,},{anchor.q:,}) "
            f"mods=({record['r1']},{record['r2']})"
        )

        print(
            f"    (k,l)=({record['k']},{record['l']}) "
            f"K={record['K']} "
            f"T={record['T']} "
            f"E={record['E']}"
        )

        for aux in record["aux"][:4]:

            print(
                f"    S={aux.S:,} "
                f"x={aux.x:,} "
                f"n+x={aux.nx:,} "
                f"factors=({aux.px:,},{aux.qx:,})"
            )

            print(
                f"        "
                f"(kx,lx)=({aux.kx},{aux.lx}) "
                f"Kx={aux.Kx} "
                f"Tx={aux.Tx} "
                f"Ex={aux.Ex} "
                f"KxPairs={aux.number_of_Kx_pairs}"
            )

            print(
                f"        "
                f"dK={aux.Kx-record['K']} "
                f"dE={aux.Ex-record['E']} "
                f"dT={aux.Tx-record['T']}"
            )

    # --------------------------------------------------------------------------------------
    # COMPACT SUMMARY
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"anchors                    = "
        f"{len(all_records)}"
    )

    print(
        f"successful auxiliary pairs = "
        f"{successful_aux}"
    )

    print(
        f"identity checks             = "
        f"{identity_checks}"
    )

    print()

    print(
        "The recursive objects produced are:"
    )

    print(
        "    n+x"
    )

    print(
        "      -> factorization of n+x"
    )

    print(
        "      -> (kx,lx)"
    )

    print(
        "      -> Kx=kx*lx"
    )

    print(
        "      -> factorization of Kx"
    )

    print()

    print(
        "The unresolved bridge remains:"
    )

    print(
        "    Kx  --->  original K"
    )

    print(
        "No assumption Kx=K is made."
    )

    # --------------------------------------------------------------------------------------
    # TIMING
    # --------------------------------------------------------------------------------------

    elapsed = time.perf_counter() - start

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime              = "
        f"{elapsed:.3f} s"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
