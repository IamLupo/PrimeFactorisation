#!/usr/bin/env python3

"""
============================================================================================
DETERMINISTIC SMOOTH-n+x / K-TRAJECTORY EXPERIMENT
============================================================================================

For an odd target n, construct x deterministically so that

    n + x == 0 (mod S)

for selected smooth S.

Then:

    n+x = S*m

is factored in the experiment.

From the auxiliary factors:

    p_x = a_x + k_x*r1
    q_x = b_x + l_x*r2

we obtain:

    k_x = floor(p_x/r1)
    l_x = floor(q_x/r2)

and:

    K_x = k_x*l_x

while:

    T_x = floor((n+x)/(r1*r2))
    E_x = T_x-K_x.

The original target has:

    K = k*l
    T = floor(n/(r1*r2))
    E = T-K.

The purpose is to expose the recursive relationship between:

    n
      ->
    n+x
      ->
    (k_x,l_x)
      ->
    K_x.

No search over x is performed.
No enumeration of the original p interval is performed.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from typing import Optional


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

# Smooth multipliers.
#
# All contain 2, so for odd n the resulting deterministic x is odd.
#
SMOOTH_S = [
    30,          # 2*3*5
    210,         # 2*3*5*7
    2310,        # 2*3*5*7*11
    2730,        # 2*3*5*7*13
    30030,       # 2*3*5*7*11*13
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
class Auxiliary:
    S: int
    x: int
    nx: int
    cofactor: int

    p_x: int
    q_x: int

    r1: int
    r2: int

    k_x: int
    l_x: int

    K_x: int
    T_x: int
    E_x: int


# ==========================================================================================
# NUMBER THEORY
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
    """
    Trial factorization for the small auxiliary/K values used here.
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


def divisor_pairs(n: int) -> list[tuple[int, int]]:
    """
    Return unordered positive divisor pairs (a,b), a <= b.
    """
    if n <= 0:
        return []

    factors = factor_integer(n)

    divisors = [1]

    for prime, exponent in factors.items():
        current = list(divisors)
        power = 1

        for _ in range(exponent):
            power *= prime

            for d in current:
                divisors.append(d * power)

    divisors = sorted(set(divisors))

    pairs = []

    for d in divisors:
        if d * d > n:
            break

        if n % d == 0:
            pairs.append((d, n // d))

    return pairs


def factor_semiprime_exact(n: int) -> Optional[tuple[int, int]]:
    """
    Exact semiprime factorization for the deliberately small
    experimental numbers.

    Returns (p,q) with p <= q when both are prime.
    """
    if n < 4:
        return None

    if n % 2 == 0:
        q = n // 2

        if is_prime(q):
            return 2, q

    d = 3

    while d * d <= n:
        if n % d == 0:
            q = n // d

            if is_prime(d) and is_prime(q):
                return d, q

        d += 2

    return None


# ==========================================================================================
# ANCHOR GENERATION
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

    anchors.sort(key=lambda a: a.n)

    return anchors


# ==========================================================================================
# CLOSE MODULUS PAIRS
# ==========================================================================================

def build_close_pairs(
    mod_primes: list[int],
    ratio: float,
) -> list[tuple[int, int]]:

    pairs: list[tuple[int, int]] = []

    for i, r1 in enumerate(mod_primes):

        for r2 in mod_primes[i + 1:]:

            gap = r2 - r1

            if gap <= ratio * r1:
                pairs.append((r1, r2))

    return pairs


# ==========================================================================================
# DETERMINISTIC x
# ==========================================================================================

def deterministic_x(
    n: int,
    S: int,
) -> int:
    """
    Smallest positive x satisfying

        n+x == 0 (mod S).
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

    T = (p * q) // R

    E = T - K

    return k, l, K, T, E


# ==========================================================================================
# AUXILIARY FACTORIZATION
# ==========================================================================================

def build_auxiliary(
    n: int,
    r1: int,
    r2: int,
    S: int,
) -> Optional[Auxiliary]:

    x = deterministic_x(n, S)

    nx = n + x

    if nx % S != 0:
        raise AssertionError(
            "Construction error: n+x is not divisible by S."
        )

    cofactor = nx // S

    aux_factor = factor_semiprime_exact(nx)

    if aux_factor is None:
        return None

    p_x, q_x = aux_factor

    if p_x > q_x:
        p_x, q_x = q_x, p_x

    k_x = p_x // r1
    l_x = q_x // r2

    K_x = k_x * l_x

    R = r1 * r2

    T_x = nx // R

    E_x = T_x - K_x

    return Auxiliary(
        S=S,
        x=x,
        nx=nx,
        cofactor=cofactor,
        p_x=p_x,
        q_x=q_x,
        r1=r1,
        r2=r2,
        k_x=k_x,
        l_x=l_x,
        K_x=K_x,
        T_x=T_x,
        E_x=E_x,
    )


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
    K_x: int,
    E_x: int,
) -> None:

    R = r1 * r2

    T = n // R
    T_x = (n + x) // R

    lhs = K_x - K
    rhs = (T_x - T) - (E_x - E)

    if lhs != rhs:
        raise AssertionError(
            "Trajectory identity failed."
        )


# ==========================================================================================
# Kx FACTORIZATION
# ==========================================================================================

def Kx_factor_structure(
    K_x: int,
    r1: int,
    r2: int,
) -> tuple[dict[int, int], list[tuple[int, int]]]:

    if K_x <= 0:
        return {}, []

    factors = factor_integer(K_x)

    pairs = divisor_pairs(K_x)

    admissible = []

    for k, l in pairs:

        p_lo = k * r1
        p_hi = (k + 1) * r1 - 1

        q_lo = l * r2
        q_hi = (l + 1) * r2 - 1

        if p_hi < FACTOR_MIN:
            continue

        if p_lo > FACTOR_MAX:
            continue

        if q_hi < FACTOR_MIN:
            continue

        if q_lo > FACTOR_MAX:
            continue

        admissible.append((k, l))

    return factors, admissible


# ==========================================================================================
# SINGLE ANCHOR
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

        aux = build_auxiliary(
            n=n,
            r1=r1,
            r2=r2,
            S=S,
        )

        if aux is None:
            continue

        verify_trajectory_identity(
            n=n,
            x=aux.x,
            r1=r1,
            r2=r2,
            K=K,
            E=E,
            K_x=aux.K_x,
            E_x=aux.E_x,
        )

        factors, pairs = Kx_factor_structure(
            aux.K_x,
            r1,
            r2,
        )

        result["aux"].append(
            {
                "S": aux.S,
                "x": aux.x,
                "nx": aux.nx,
                "cofactor": aux.cofactor,
                "p_x": aux.p_x,
                "q_x": aux.q_x,
                "kx": aux.k_x,
                "lx": aux.l_x,
                "Kx": aux.K_x,
                "Tx": aux.T_x,
                "Ex": aux.E_x,
                "Kx_factors": factors,
                "Kx_pairs": pairs,
            }
        )

    return result


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    start = time.perf_counter()

    print("=" * 100)
    print("DETERMINISTIC SMOOTH-n+x / K-TRAJECTORY EXPERIMENT")
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
        f"seed                      = {SEED:,}"
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
    # MODULUS PAIRS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING CLOSE MODULUS PAIRS")
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

    # --------------------------------------------------------------------------------------
    # SELECT PAIRS
    # --------------------------------------------------------------------------------------

    rng = random.Random(SEED)

    selected_pairs = [
        rng.choice(close_pairs)
        for _ in anchors
    ]

    # --------------------------------------------------------------------------------------
    # RUN
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("RUNNING DETERMINISTIC n+x TRAJECTORIES")
    print("=" * 100)

    results = []

    for index, (anchor, pair) in enumerate(
        zip(anchors, selected_pairs),
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
                f"anchor {index:3d}/{len(anchors):3d}"
            )

    print()

    # --------------------------------------------------------------------------------------
    # RAW RECORDS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("RAW TRAJECTORY RECORDS")
    print("=" * 100)

    for result in results:

        print(
            f"n={result['n']:,} "
            f"p={result['p']:,} "
            f"q={result['q']:,} "
            f"mods=({result['r1']},{result['r2']})"
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
                "    no auxiliary semiprime found "
                "for selected S values"
            )

        for aux in result["aux"]:

            Kx = aux["Kx"]
            K = result["K"]

            Tx = aux["Tx"]
            T = result["T"]

            Ex = aux["Ex"]
            E = result["E"]

            delta_K = Kx - K
            delta_T = Tx - T
            delta_E = Ex - E

            print()

            print(
                f"    S={aux['S']:,}"
            )

            print(
                f"        x={aux['x']:,}"
            )

            print(
                f"        n+x={aux['nx']:,}"
            )

            print(
                f"        cofactor=(n+x)/S="
                f"{aux['cofactor']:,}"
            )

            print(
                f"        auxiliary factors="
                f"({aux['p_x']:,},{aux['q_x']:,})"
            )

            print(
                f"        (kx,lx)="
                f"({aux['kx']},{aux['lx']})"
            )

            print(
                f"        Kx={Kx}"
            )

            print(
                f"        Tx={Tx}"
            )

            print(
                f"        Ex={Ex}"
            )

            print(
                f"        factor(Kx)="
                f"{aux['Kx_factors']}"
            )

            print(
                f"        divisor pairs(Kx)="
                f"{aux['Kx_pairs']}"
            )

            print()

            print(
                "        trajectory identity:"
            )

            print(
                f"            Kx-K = {delta_K}"
            )

            print(
                f"            (Tx-T)-(Ex-E) = "
                f"{delta_T - delta_E}"
            )

            print(
                f"            delta_T = {delta_T}"
            )

            print(
                f"            delta_E = {delta_E}"
            )

            print(
                f"            identity_holds = "
                f"{delta_K == delta_T - delta_E}"
            )

        print()

    # --------------------------------------------------------------------------------------
    # CONSTRUCTION CHECK
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("DETERMINISTIC CONSTRUCTION CHECK")
    print("=" * 100)

    for result in results:

        n = result["n"]

        for aux in result["aux"]:

            S = aux["S"]
            x = aux["x"]

            if x <= 0:
                raise AssertionError(
                    "x must be positive."
                )

            if (n + x) % S != 0:
                raise AssertionError(
                    "n+x is not divisible by S."
                )

    print(
        "all auxiliary x values satisfy:"
    )

    print(
        "    n+x == 0 (mod S)"
    )

    print(
        "all trajectory identities verified"
    )

    print(
        "no x search was performed"
    )

    print(
        "no original p enumeration was performed"
    )

    # --------------------------------------------------------------------------------------
    # RECURSIVE K OBJECT
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RECURSIVE K OBJECT")
    print("=" * 100)

    print(
        "For every successful auxiliary factorization:"
    )

    print()
    print(
        "    n+x"
    )
    print(
        "      -> (p_x,q_x)"
    )
    print(
        "      -> (k_x,l_x)"
    )
    print(
        "      -> K_x = k_x*l_x"
    )
    print(
        "      -> factor(K_x)"
    )
    print(
        "      -> divisor pairs of K_x"
    )

    print()
    print(
        "The unresolved bridge is:"
    )

    print()
    print(
        "    known auxiliary K_x"
    )
    print(
        "             |"
    )
    print(
        "             v"
    )
    print(
        "    constraint on original K"
    )

    print()
    print(
        "The experiment deliberately does not assume"
    )
    print(
        "that K_x == K."
    )

    print(
        "It exposes the quantities required to test"
    )
    print(
        "such a relationship mathematically."
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
        f"total runtime              = {elapsed:.3f} s"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()