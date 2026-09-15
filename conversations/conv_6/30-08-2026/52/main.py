#!/usr/bin/env python3

"""
============================================================================================
NEARBY-FACTORIZATION BRIDGE / KNOWN (Tx,Ex,kx,lx) -> ORIGINAL FACTOR EXPERIMENT
============================================================================================

Research question:

    Suppose n is hard to factor, but a nearby number

        n_x = n + x

    has an easy factorization:

        n_x = p_x * q_x.

    From that easy factorization we know:

        T_x
        E_x
        k_x
        l_x
        K_x = k_x*l_x.

Can this information be transported back to the original

        n = p*q

to recover p and q efficiently?

IMPORTANT:

    The experiment NEVER gives the original p,q to the reconstruction
    algorithm.

    p,q are used only for measuring recovery quality afterward.

The experiment tests:

    A. K_x == K ?
    B. |K_x-K| is small?
    C. k_x,l_x predict k,l?
    D. simple coordinate-difference bounds recover p,q?
    E. enumerating only a SMALL neighborhood around (kx,lx) works?
    F. combining multiple nearby factorizations improves recovery?

============================================================================================
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

SEED = 1_511_464_998

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

# Nearby factor perturbation.
# n+x is constructed from:
#
#     p_x = p + dp
#     q_x = q + dq
#
# without the algorithm being told dp,dq.
PERTURB_RADIUS = 50

# Reconstruction neighborhood around kx,lx.
K_RADIUS = 20
L_RADIUS = 20

# Product neighborhood around Kx.
K_PRODUCT_RADIUS = 500

# Number of independent nearby factorizations per anchor.
NEARBY_SAMPLES = 8

# A few relatively small modulus pairs are enough here.
CLOSE_RATIO = 0.20

PROGRESS_EVERY = 10


# ==========================================================================================
# DATA STRUCTURES
# ==========================================================================================

@dataclass
class Anchor:
    p: int
    q: int
    n: int


@dataclass
class NearbyFactor:
    x: int
    nx: int

    px: int
    qx: int

    r1: int
    r2: int

    ax: int
    bx: int
    kx: int
    lx: int

    Tx: int
    Ex: int
    Kx: int


@dataclass
class ReconstructionResult:
    exact: bool
    method: str
    candidates: int
    found_p: Optional[int]
    found_q: Optional[int]


# ==========================================================================================
# PRIME GENERATION
# ==========================================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            step = p
            start = p * p
            is_prime[start : limit + 1 : step] = b"\x00" * (
                ((limit - start) // step) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


# ==========================================================================================
# CLOSE MODULUS PAIRS
# ==========================================================================================

def build_close_pairs(primes: list[int]) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []

    for i, r1 in enumerate(primes):
        for r2 in primes[i + 1:]:
            if r2 - r1 <= int(r1 * CLOSE_RATIO):
                pairs.append((r1, r2))

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

    possible = [
        (p, q)
        for p in factor_primes
        for q in factor_primes
        if p <= q
        and FACTOR_MIN <= p <= FACTOR_MAX
        and FACTOR_MIN <= q <= FACTOR_MAX
    ]

    rng.shuffle(possible)

    chosen: list[Anchor] = []
    seen_n: set[int] = set()

    for p, q in possible:
        n = p * q

        if n in seen_n:
            continue

        chosen.append(Anchor(p=p, q=q, n=n))
        seen_n.add(n)

        if len(chosen) >= count:
            break

    if len(chosen) < count:
        raise RuntimeError(
            f"Could only build {len(chosen)} anchors."
        )

    return chosen


# ==========================================================================================
# TRUE COORDINATES
# ==========================================================================================

def coordinates(p: int, q: int, r1: int, r2: int):
    a = p % r1
    b = q % r2

    k = (p - a) // r1
    l = (q - b) // r2

    R = r1 * r2
    T = p * q // R
    E = T - k * l
    K = k * l

    return a, b, k, l, T, E, K


# ==========================================================================================
# BUILD A NEARBY FACTORIZATION
# ==========================================================================================

def make_nearby_factorization(
    anchor: Anchor,
    rng: random.Random,
) -> tuple[int, int, int]:
    """
    Construct:

        px = p + dp
        qx = q + dq

    with small perturbations.

    The reconstruction code is only given px,qx indirectly through
    the factorization of nx.

    We avoid zero perturbation so that nx != n.
    """

    p = anchor.p
    q = anchor.q

    for _ in range(1000):
        dp = rng.randint(-PERTURB_RADIUS, PERTURB_RADIUS)
        dq = rng.randint(-PERTURB_RADIUS, PERTURB_RADIUS)

        if dp == 0 and dq == 0:
            continue

        px = p + dp
        qx = q + dq

        if px < FACTOR_MIN or px > FACTOR_MAX:
            continue

        if qx < FACTOR_MIN or qx > FACTOR_MAX:
            continue

        # Ensure both remain prime.
        if px <= 1 or qx <= 1:
            continue

        # The anchor factor pool is prime, but perturbations do not preserve
        # primality. Caller checks primality using a set.
        return px, qx, dp * q + dq * p + dp * dq

    raise RuntimeError("Unable to construct nearby factorization.")


# ==========================================================================================
# BUILD NEARBY SAMPLE
# ==========================================================================================

def make_nearby_sample(
    anchor: Anchor,
    r1: int,
    r2: int,
    factor_set: set[int],
    rng: random.Random,
) -> Optional[NearbyFactor]:
    for _ in range(5000):
        px, qx, x = make_nearby_factorization(anchor, rng)

        if px not in factor_set or qx not in factor_set:
            continue

        nx = px * qx

        if nx == anchor.n:
            continue

        ax, bx, kx, lx, Tx, Ex, Kx = coordinates(
            px, qx, r1, r2
        )

        return NearbyFactor(
            x=x,
            nx=nx,
            px=px,
            qx=qx,
            r1=r1,
            r2=r2,
            ax=ax,
            bx=bx,
            kx=kx,
            lx=lx,
            Tx=Tx,
            Ex=Ex,
            Kx=Kx,
        )

    return None


# ==========================================================================================
# EXACT VERIFICATION
# ==========================================================================================

def exact_factor_from_candidates(
    n: int,
    candidates: list[tuple[int, int]],
    true_p: int,
    true_q: int,
    method: str,
) -> ReconstructionResult:

    tested = 0

    seen: set[tuple[int, int]] = set()

    for p, q in candidates:
        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))
        tested += 1

        if p * q == n:
            return ReconstructionResult(
                exact=True,
                method=method,
                candidates=tested,
                found_p=p,
                found_q=q,
            )

    return ReconstructionResult(
        exact=False,
        method=method,
        candidates=tested,
        found_p=None,
        found_q=None,
    )


# ==========================================================================================
# RECONSTRUCTION 1
# ==========================================================================================

def reconstruct_from_kx_lx(
    anchor: Anchor,
    sample: NearbyFactor,
    k_radius: int,
    l_radius: int,
) -> ReconstructionResult:

    """
    Candidate generation:

        k = kx + dk
        l = lx + dl

    Then derive the possible p,q intervals from the quotient coordinates.

    Because residues are unknown, every coordinate block is represented
    by its full interval.

    This is deliberately a fairly weak reconstruction rule.
    """

    n = anchor.n
    r1 = sample.r1
    r2 = sample.r2

    candidates: list[tuple[int, int]] = []

    k_lo = max(0, sample.kx - k_radius)
    k_hi = sample.kx + k_radius

    l_lo = max(0, sample.lx - l_radius)
    l_hi = sample.lx + l_radius

    for k in range(k_lo, k_hi + 1):
        p_lo = k * r1
        p_hi = (k + 1) * r1 - 1

        for l in range(l_lo, l_hi + 1):
            q_lo = l * r2
            q_hi = (l + 1) * r2 - 1

            # We still have to search the block, but this is much smaller
            # than the original complete p range.
            for p in range(p_lo, p_hi + 1):
                if p < FACTOR_MIN or p > FACTOR_MAX:
                    continue

                if p > q_hi:
                    continue

                q = n // p

                if q < q_lo or q > q_hi:
                    continue

                candidates.append((p, q))

    return exact_factor_from_candidates(
        n,
        candidates,
        anchor.p,
        anchor.q,
        "kx/lx neighborhood",
    )


# ==========================================================================================
# RECONSTRUCTION 2
# ==========================================================================================

def reconstruct_from_Kx(
    anchor: Anchor,
    sample: NearbyFactor,
    radius: int,
) -> ReconstructionResult:

    """
    Treat:

        K = k*l

    as the hidden quantity.

    We test:

        K in [Kx-radius, Kx+radius]

    and factor every candidate K.

    """

    n = anchor.n
    r1 = sample.r1
    r2 = sample.r2

    candidates: list[tuple[int, int]] = []

    for K in range(
        max(0, sample.Kx - radius),
        sample.Kx + radius + 1,
    ):
        if K == 0:
            continue

        root = int(math.isqrt(K))

        for k in range(1, root + 1):
            if K % k != 0:
                continue

            l = K // k

            # Both orientations.
            quotient_pairs = ((k, l), (l, k))

            for kk, ll in quotient_pairs:
                p_lo = kk * r1
                p_hi = (kk + 1) * r1 - 1

                q_lo = ll * r2
                q_hi = (ll + 1) * r2 - 1

                # p*q=n means q=floor(n/p).
                for p in range(p_lo, p_hi + 1):
                    if p < FACTOR_MIN or p > FACTOR_MAX:
                        continue

                    if p > q_hi:
                        continue

                    q = n // p

                    if q_lo <= q <= q_hi:
                        candidates.append((p, q))

    return exact_factor_from_candidates(
        n,
        candidates,
        anchor.p,
        anchor.q,
        "Kx neighborhood",
    )


# ==========================================================================================
# RECONSTRUCTION 3
# ==========================================================================================

def reconstruct_using_E_transform(
    anchor: Anchor,
    sample: NearbyFactor,
    e_radius: int,
) -> ReconstructionResult:

    """
    Use:

        K = T - E

    for the original n.

    We know Tx and Ex for n+x.

    Test a local neighborhood around Ex while correcting for the
    floor change:

        T = floor(n/R)

        Tx = floor((n+x)/R).

    Thus:

        Ex + deltaE

    produces:

        K = T - E.

    """

    n = anchor.n
    r1 = sample.r1
    r2 = sample.r2
    R = r1 * r2

    T = n // R
    dT = sample.Tx - T

    candidates: list[tuple[int, int]] = []

    for dE in range(-e_radius, e_radius + 1):
        E = sample.Ex + dE - dT
        K = T - E

        if K <= 0:
            continue

        root = int(math.isqrt(K))

        for k in range(1, root + 1):
            if K % k:
                continue

            l = K // k

            for kk, ll in ((k, l), (l, k)):
                p_lo = kk * r1
                p_hi = (kk + 1) * r1 - 1

                q_lo = ll * r2
                q_hi = (ll + 1) * r2 - 1

                for p in range(p_lo, p_hi + 1):
                    if p < FACTOR_MIN or p > FACTOR_MAX:
                        continue

                    if p > q_hi:
                        continue

                    q = n // p

                    if q_lo <= q <= q_hi:
                        candidates.append((p, q))

    return exact_factor_from_candidates(
        n,
        candidates,
        anchor.p,
        anchor.q,
        "Ex-neighborhood",
    )


# ==========================================================================================
# MAIN EXPERIMENT
# ==========================================================================================

def run() -> None:

    start_total = time.perf_counter()

    print("=" * 92)
    print("NEARBY-FACTORIZATION BRIDGE / KNOWN (Tx,Ex,kx,lx) -> ORIGINAL FACTOR")
    print("=" * 92)
    print(f"N anchors                 = {N_ANCHORS}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range             = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"perturbation radius       = {PERTURB_RADIUS}")
    print(f"nearby samples / anchor   = {NEARBY_SAMPLES}")
    print(f"k neighborhood radius    = {K_RADIUS}")
    print(f"l neighborhood radius    = {L_RADIUS}")
    print(f"K neighborhood radius    = {K_PRODUCT_RADIUS}")
    print(f"seed                      = {SEED}")
    print()

    # ------------------------------------------------------------------
    # Prime pools
    # ------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING PRIME POOLS")
    print("=" * 92)

    factor_primes = sieve(FACTOR_MAX)
    factor_primes = [
        p for p in factor_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]
    factor_set = set(factor_primes)

    modulus_primes = sieve(MOD_MAX)
    modulus_primes = [
        p for p in modulus_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    close_pairs = build_close_pairs(modulus_primes)

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print(f"close modulus pairs       = {len(close_pairs):,}")
    print()

    # ------------------------------------------------------------------
    # Anchors
    # ------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING ANCHORS")
    print("=" * 92)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print(f"actual anchors            = {len(anchors)}")
    print()

    rng = random.Random(SEED)

    # ------------------------------------------------------------------
    # Results
    # ------------------------------------------------------------------

    total_samples = 0

    same_K = 0
    same_E = 0

    sum_abs_dK = 0
    sum_abs_dE = 0

    max_abs_dK = 0
    max_abs_dE = 0

    kx_success = 0
    Kx_success = 0
    Ex_success = 0

    kx_candidates = 0
    Kx_candidates = 0
    Ex_candidates = 0

    interesting: list[tuple] = []

    # ------------------------------------------------------------------
    # Search
    # ------------------------------------------------------------------

    print("=" * 92)
    print("RUNNING NEARBY-FACTORIZATION BRIDGE")
    print("=" * 92)

    for idx, anchor in enumerate(anchors, 1):

        # Select one modulus pair for this anchor.
        r1, r2 = rng.choice(close_pairs)

        for _ in range(NEARBY_SAMPLES):

            sample = make_nearby_sample(
                anchor,
                r1,
                r2,
                factor_set,
                rng,
            )

            if sample is None:
                continue

            total_samples += 1

            # Original coordinates, only for evaluation.
            _, _, k, l, T, E, K = coordinates(
                anchor.p,
                anchor.q,
                r1,
                r2,
            )

            dK = sample.Kx - K
            dE = sample.Ex - E

            if dK == 0:
                same_K += 1

            if dE == 0:
                same_E += 1

            sum_abs_dK += abs(dK)
            sum_abs_dE += abs(dE)

            max_abs_dK = max(max_abs_dK, abs(dK))
            max_abs_dE = max(max_abs_dE, abs(dE))

            # ----------------------------------------------------------
            # Reconstruction tests
            # ----------------------------------------------------------

            result_kx = reconstruct_from_kx_lx(
                anchor,
                sample,
                K_RADIUS,
                L_RADIUS,
            )

            result_Kx = reconstruct_from_Kx(
                anchor,
                sample,
                K_PRODUCT_RADIUS,
            )

            result_Ex = reconstruct_using_E_transform(
                anchor,
                sample,
                e_radius=K_PRODUCT_RADIUS,
            )

            if result_kx.exact:
                kx_success += 1
                kx_candidates += result_kx.candidates

            if result_Kx.exact:
                Kx_success += 1
                Kx_candidates += result_Kx.candidates

            if result_Ex.exact:
                Ex_success += 1
                Ex_candidates += result_Ex.candidates

            # ----------------------------------------------------------
            # Interesting samples:
            # nearby factorization changes, but K/E survive.
            # ----------------------------------------------------------

            if (
                sample.px != anchor.p
                or sample.qx != anchor.q
            ):
                if sample.Kx == K or sample.Ex == E:
                    interesting.append(
                        (
                            anchor,
                            sample,
                            k,
                            l,
                            K,
                            E,
                        )
                    )

        if idx % PROGRESS_EVERY == 0 or idx == len(anchors):
            print(f"anchor {idx:4d}/{len(anchors)}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    print()
    print("=" * 92)
    print("SUMMARY")
    print("=" * 92)

    print(f"anchors analyzed           = {len(anchors):,}")
    print(f"nearby factorizations      = {total_samples:,}")
    print()

    if total_samples:
        print("NEARBY INVARIANT TEST")
        print(f"    Kx = K                  = {same_K:,}/{total_samples:,}")
        print(
            f"    fraction               = "
            f"{same_K / total_samples:.6f}"
        )
        print()

        print("    Ex = E                  =", same_E,
              "/", total_samples)
        print(
            f"    fraction               = "
            f"{same_E / total_samples:.6f}"
        )
        print()

        print("K DIFFERENCE")
        print(
            f"    mean |Kx-K|            = "
            f"{sum_abs_dK / total_samples:.6f}"
        )
        print(
            f"    max  |Kx-K|            = "
            f"{max_abs_dK:,}"
        )
        print()

        print("E DIFFERENCE")
        print(
            f"    mean |Ex-E|            = "
            f"{sum_abs_dE / total_samples:.6f}"
        )
        print(
            f"    max  |Ex-E|            = "
            f"{max_abs_dE:,}"
        )
        print()

    print("=" * 92)
    print("RECONSTRUCTION FROM KNOWN NEARBY FACTORIZATION")
    print("=" * 92)

    print("METHOD 1: kx/lx neighborhood")
    print(
        f"    exact recovery         = "
        f"{kx_success}/{total_samples}"
    )
    print(
        f"    recovery rate          = "
        f"{kx_success / total_samples if total_samples else 0:.6f}"
    )
    print(
        f"    total candidate tests  = "
        f"{kx_candidates:,}"
    )
    print()

    print("METHOD 2: Kx neighborhood")
    print(
        f"    exact recovery         = "
        f"{Kx_success}/{total_samples}"
    )
    print(
        f"    recovery rate          = "
        f"{Kx_success / total_samples if total_samples else 0:.6f}"
    )
    print(
        f"    total candidate tests  = "
        f"{Kx_candidates:,}"
    )
    print()

    print("METHOD 3: Ex neighborhood")
    print(
        f"    exact recovery         = "
        f"{Ex_success}/{total_samples}"
    )
    print(
        f"    recovery rate          = "
        f"{Ex_success / total_samples if total_samples else 0:.6f}"
    )
    print(
        f"    total candidate tests  = "
        f"{Ex_candidates:,}"
    )
    print()

    # ------------------------------------------------------------------
    # Interesting cases
    # ------------------------------------------------------------------

    print("=" * 92)
    print("NONTRIVIAL E/K STABILITY CASES")
    print("=" * 92)

    shown = 0

    for anchor, sample, k, l, K, E in interesting[:20]:

        print(
            f"n={anchor.n:,} x={sample.x:+,} "
            f"orig=({anchor.p:,},{anchor.q:,}) "
            f"near=({sample.px:,},{sample.qx:,})"
        )

        print(
            f"    mods=({sample.r1},{sample.r2}) "
            f"K={K}->{sample.Kx} "
            f"E={E}->{sample.Ex}"
        )

        print(
            f"    coords=(k,l)=({k},{l}) "
            f"-> ({sample.kx},{sample.lx})"
        )

        print()

        shown += 1

    if shown == 0:
        print("No nontrivial invariant cases found.")
        print()

    # ------------------------------------------------------------------
    # Explicit examples
    # ------------------------------------------------------------------

    print("=" * 92)
    print("EXAMPLE DATA")
    print("=" * 92)

    example_count = 0

    for anchor in anchors:

        if example_count >= 20:
            break

        r1, r2 = rng.choice(close_pairs)

        sample = make_nearby_sample(
            anchor,
            r1,
            r2,
            factor_set,
            rng,
        )

        if sample is None:
            continue

        _, _, k, l, T, E, K = coordinates(
            anchor.p,
            anchor.q,
            r1,
            r2,
        )

        print(
            f"n={anchor.n:,}"
            f"  n+x={sample.nx:,}"
            f"  x={sample.x:+,}"
        )

        print(
            f"    original hidden: "
            f"(p,q)=({anchor.p:,},{anchor.q:,}) "
            f"(k,l)=({k},{l}) "
            f"K={K} E={E} T={T}"
        )

        print(
            f"    known nearby:   "
            f"(kx,lx)=({sample.kx},{sample.lx}) "
            f"Kx={sample.Kx} "
            f"Ex={sample.Ex} "
            f"Tx={sample.Tx}"
        )

        print(
            f"    dK={sample.Kx-K:+d} "
            f"dE={sample.Ex-E:+d}"
        )

        print()

        example_count += 1

    # ------------------------------------------------------------------
    # Algebra
    # ------------------------------------------------------------------

    print("=" * 92)
    print("ALGEBRAIC RELATION")
    print("=" * 92)

    print(
        """
For every modulus pair:

    R  = r1*r2
    T  = floor(n/R)
    E  = T-k*l

For the nearby number:

    nx = n+x
    Tx = floor(nx/R)
    Ex = Tx-kx*lx

Therefore:

    K  = k*l  = T - E
    Kx = kx*lx = Tx - Ex

and exactly:

    Kx-K = (Tx-T) - (Ex-E).

Thus if the nearby factorization supplies:

    Tx, Ex, kx, lx

we know Kx exactly.

The experiment asks whether Kx is predictive of K.

But knowing Kx is not sufficient by itself.

We need some exploitable relation such as:

    K = Kx
    K approximately Kx
    K = F(Kx,x)
    E = F(Ex,x)
    (k,l) = F(kx,lx,x)

The reconstruction tests deliberately distinguish these
possibilities.
"""
    )

    print("=" * 92)
    print("IMPORTANT INTERPRETATION")
    print("=" * 92)

    print(
        """
A successful result here would be materially different from the
earlier modular experiments.

Earlier:

    n
      |
      +--> congruence filter
      |
      +--> enumerate candidates.

Here:

    n+x
      |
      +--> easy factorization
      |
      +--> known Tx,Ex,kx,lx
      |
      +--> predict information about n
      |
      +--> generate a SMALL candidate set
      |
      +--> exact n factorization.

The crucial metric is not simply:

    Kx == K

It is:

    Can the factorization of n+x reduce the search for n?

The strongest useful outcome would be:

    nearby factorization
        ->
    narrow K / (k,l) range
        ->
    only a handful of candidates for n
        ->
    exact recovery.

If the required neighborhood grows proportionally with the
original factor size, then the nearby factorization has not
removed the underlying search.

"""
    )

    total_time = time.perf_counter() - start_total

    print("=" * 92)
    print("TIMING")
    print("=" * 92)
    print(f"total runtime               = {total_time:.3f} seconds")

    print()
    print("=" * 92)
    print("EXPERIMENT COMPLETE")
    print("=" * 92)


if __name__ == "__main__":
    run()
