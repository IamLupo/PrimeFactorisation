#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 26
HASH-INDEXED CONSTRUCTIVE SIGNATURE COLLISION DEPTH
NO CSV OUTPUT
==============================================================================

Purpose
-------
Continue Experiment 25, but extend the constructive search from Fpair to
all tested signatures.

The search never enumerates all P^2 prime pairs.

For algebraic signatures, each candidate first factor determines the required
signature of the second factor, which is then resolved through a hash index.

Signatures
----------
Fpair
Fsorted
Fsum
Fprod
Fdiff
cube_sum
order_pair
order_sorted

The experiment searches prefixes 3..7 and reports explicit wrong pairs.

Important
---------
NONE FOUND means only that the finite prime pool contains no collision.
It does NOT prove mathematical uniqueness.

No CSV output.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

TARGETS = 12

PRIME_POOL_SIZE = 5_000
PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

# Search at most this many candidate first primes for each target/signature.
# None means the complete pool.
MAX_CANDIDATES_PER_SEARCH = None

# If True, print the first explicit collision only.
FIRST_COLLISION_ONLY = True

# Signature names.
SIGNATURES = [
    "Fpair",
    "Fsorted",
    "Fsum",
    "Fprod",
    "Fdiff",
    "cube_sum",
    "order_pair",
    "order_sorted",
]


# ============================================================================
# FIXED TARGETS FROM EXPERIMENT 25
# ============================================================================

TARGET_PAIRS = [
    (2280911, 3904289),
    (2451809, 2094101),
    (2013329, 2228711),
    (2276231, 4125593),
    (2452433, 2120863),
    (2585953, 4093603),
    (3611719, 2242469),
    (3455059, 2706413),
    (2120101, 2922551),
    (3244361, 3267037),
    (2708999, 3641857),
    (3353047, 2857303),
]


# ============================================================================
# MODULUS INFORMATION
# ============================================================================

def factor_integer(n: int) -> Dict[int, int]:
    """Small integer factorization by trial division."""
    factors: Dict[int, int] = {}
    d = 2

    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d

        d = 3 if d == 2 else d + 2

    if n > 1:
        factors[n] = factors.get(n, 0) + 1

    return factors


def phi(n: int) -> int:
    result = n
    for p in factor_integer(n):
        result -= result // p
    return result


MODULI = [2 * r * r - 1 for r in R_VALUES]


def modulus_inventory() -> List[Tuple[int, int, int]]:
    rows = []

    for r, m in zip(R_VALUES, MODULI):
        rows.append((r, m, phi(m)))

    return rows


# ============================================================================
# PRIME GENERATION
# ============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

    for p in small:
        if n % p == 0:
            return n == p

    d = 41
    step = 2

    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step

    return True


def generate_prime_pool(
    size: int,
    low: int,
    high: int,
    seed: int,
) -> List[int]:

    rng = random.Random(seed)

    # Start with the fixed targets so they are guaranteed to be available.
    pool_set = set()

    for p, q in TARGET_PAIRS:
        if low <= p <= high and is_prime(p):
            pool_set.add(p)
        if low <= q <= high and is_prime(q):
            pool_set.add(q)

    while len(pool_set) < size:
        x = rng.randrange(low | 1, high + 1, 2)

        if is_prime(x):
            pool_set.add(x)

    return sorted(pool_set)


# ============================================================================
# RESIDUE CACHE
# ============================================================================

def build_residue_cache(
    pool: Sequence[int],
) -> Dict[int, Tuple[int, ...]]:

    cache: Dict[int, Tuple[int, ...]] = {}

    for p in pool:
        cache[p] = tuple(p % m for m in MODULI)

    return cache


# ============================================================================
# MULTIPLICATIVE ORDERS
# ============================================================================

def multiplicative_order(a: int, m: int, phi_m: int, phi_factors) -> int:
    """
    Compute ord_m(a), assuming gcd(a,m)=1.

    We reduce the candidate order from phi(m) using its prime factors.
    """

    if math.gcd(a, m) != 1:
        return 0

    order = phi_m

    for prime_factor in phi_factors:
        while order % prime_factor == 0:
            candidate = order // prime_factor

            if pow(a, candidate, m) == 1:
                order = candidate
            else:
                break

    return order


PHI_VALUES = [phi(m) for m in MODULI]
PHI_FACTORS = [tuple(factor_integer(x).keys()) for x in PHI_VALUES]


def build_order_cache(
    pool: Sequence[int],
) -> Dict[int, Tuple[int, ...]]:

    cache: Dict[int, Tuple[int, ...]] = {}

    for p in pool:
        values = []

        for m, phi_m, phi_f in zip(
            MODULI,
            PHI_VALUES,
            PHI_FACTORS,
        ):
            values.append(
                multiplicative_order(
                    p % m,
                    m,
                    phi_m,
                    phi_f,
                )
            )

        cache[p] = tuple(values)

    return cache


# ============================================================================
# PREFIX SIGNATURES
# ============================================================================

def fpair_signature(
    p: int,
    q: int,
    prefix: int,
    residues,
) -> Tuple[int, ...]:

    rp = residues[p][:prefix]
    rq = residues[q][:prefix]

    result = []

    for a, b in zip(rp, rq):
        result.extend((a, b))

    return tuple(result)


def fsorted_signature(
    p: int,
    q: int,
    prefix: int,
    residues,
) -> Tuple[int, ...]:

    rp = residues[p][:prefix]
    rq = residues[q][:prefix]

    result = []

    for a, b in zip(rp, rq):
        if a <= b:
            result.extend((a, b))
        else:
            result.extend((b, a))

    return tuple(result)


def fsum_signature(
    p: int,
    q: int,
    prefix: int,
    residues,
) -> Tuple[int, ...]:

    return tuple(
        (a + b) % m
        for a, b, m in zip(
            residues[p][:prefix],
            residues[q][:prefix],
            MODULI[:prefix],
        )
    )


def fprod_signature(
    p: int,
    q: int,
    prefix: int,
    residues,
) -> Tuple[int, ...]:

    return tuple(
        (a * b) % m
        for a, b, m in zip(
            residues[p][:prefix],
            residues[q][:prefix],
            MODULI[:prefix],
        )
    )


def fdiff_signature(
    p: int,
    q: int,
    prefix: int,
    residues,
) -> Tuple[int, ...]:

    return tuple(
        (a - b) % m
        for a, b, m in zip(
            residues[p][:prefix],
            residues[q][:prefix],
            MODULI[:prefix],
        )
    )


def cube_sum_signature(
    p: int,
    q: int,
    prefix: int,
    residues,
) -> Tuple[int, ...]:

    return tuple(
        (pow(a, 3, m) + pow(b, 3, m)) % m
        for a, b, m in zip(
            residues[p][:prefix],
            residues[q][:prefix],
            MODULI[:prefix],
        )
    )


def order_pair_signature(
    p: int,
    q: int,
    prefix: int,
    orders,
) -> Tuple[int, ...]:

    op = orders[p][:prefix]
    oq = orders[q][:prefix]

    result = []

    for a, b in zip(op, oq):
        result.extend((a, b))

    return tuple(result)


def order_sorted_signature(
    p: int,
    q: int,
    prefix: int,
    orders,
) -> Tuple[int, ...]:

    op = orders[p][:prefix]
    oq = orders[q][:prefix]

    result = []

    for a, b in zip(op, oq):
        if a <= b:
            result.extend((a, b))
        else:
            result.extend((b, a))

    return tuple(result)


# ============================================================================
# GENERIC SIGNATURE DISPATCH
# ============================================================================

def signature(
    name: str,
    p: int,
    q: int,
    prefix: int,
    residues,
    orders,
):
    if name == "Fpair":
        return fpair_signature(p, q, prefix, residues)

    if name == "Fsorted":
        return fsorted_signature(p, q, prefix, residues)

    if name == "Fsum":
        return fsum_signature(p, q, prefix, residues)

    if name == "Fprod":
        return fprod_signature(p, q, prefix, residues)

    if name == "Fdiff":
        return fdiff_signature(p, q, prefix, residues)

    if name == "cube_sum":
        return cube_sum_signature(p, q, prefix, residues)

    if name == "order_pair":
        return order_pair_signature(p, q, prefix, orders)

    if name == "order_sorted":
        return order_sorted_signature(p, q, prefix, orders)

    raise ValueError(f"Unknown signature: {name}")


# ============================================================================
# INDEX CONSTRUCTION
# ============================================================================

def build_pair_index(
    pool: Sequence[int],
    name: str,
    prefix: int,
    residues,
    orders,
):
    """
    Build an index of prime-pair signatures.

    Important optimization:
    We do NOT build all P^2 pairs.

    Instead, for each first prime p, we store the signature contribution
    needed to solve for q where possible.

    The implementation below uses algebraic complement indexes.
    """

    index = defaultdict(list)

    if name in ("Fpair", "Fsorted", "order_pair", "order_sorted"):
        # These signatures can be handled using per-prime signature vectors.
        #
        # For Fpair and order_pair, a collision requires each factor's
        # individual vector to match the corresponding target vector.
        #
        # Fsorted/order_sorted additionally permit swapping.
        for p in pool:
            if name == "Fpair":
                key = residues[p][:prefix]

            elif name == "order_pair":
                key = orders[p][:prefix]

            else:
                # The individual vector is still useful as a lookup index.
                key = residues[p][:prefix] if name == "Fsorted" else orders[p][:prefix]

            index[key].append(p)

        return index

    # For scalar signatures, index each prime by its individual residue
    # vector. This permits algebraic complement lookup.
    for p in pool:
        index[residues[p][:prefix]].append(p)

    return index


# ============================================================================
# TARGET / PAIR VALIDATION
# ============================================================================

def pair_is_valid(
    p: int,
    q: int,
    target_p: int,
    target_q: int,
) -> bool:

    if p == q:
        return False

    if {p, q} == {target_p, target_q}:
        return False

    return True


def verify_collision(
    name: str,
    p: int,
    q: int,
    target_p: int,
    target_q: int,
    prefix: int,
    residues,
    orders,
) -> bool:

    if not pair_is_valid(p, q, target_p, target_q):
        return False

    target_sig = signature(
        name,
        target_p,
        target_q,
        prefix,
        residues,
        orders,
    )

    candidate_sig = signature(
        name,
        p,
        q,
        prefix,
        residues,
        orders,
    )

    return target_sig == candidate_sig


# ============================================================================
# CONSTRUCTIVE SEARCH
# ============================================================================

def find_collision(
    name: str,
    prefix: int,
    target_p: int,
    target_q: int,
    pool: Sequence[int],
    residues,
    orders,
) -> Optional[Tuple[int, int]]:

    target_sig = signature(
        name,
        target_p,
        target_q,
        prefix,
        residues,
        orders,
    )

    # ------------------------------------------------------------------------
    # FPAIR
    # ------------------------------------------------------------------------

    if name == "Fpair":
        target_p_vec = residues[target_p][:prefix]
        target_q_vec = residues[target_q][:prefix]

        class_p = [
            x for x in pool
            if residues[x][:prefix] == target_p_vec
        ]

        class_q = [
            x for x in pool
            if residues[x][:prefix] == target_q_vec
        ]

        for p in class_p:
            for q in class_q:
                if pair_is_valid(p, q, target_p, target_q):
                    if verify_collision(
                        name,
                        p,
                        q,
                        target_p,
                        target_q,
                        prefix,
                        residues,
                        orders,
                    ):
                        return p, q

        return None

    # ------------------------------------------------------------------------
    # FSORTED
    # ------------------------------------------------------------------------

    if name == "Fsorted":
        target_rp = residues[target_p][:prefix]
        target_rq = residues[target_q][:prefix]

        # Try both orientations.
        orientations = [
            (target_rp, target_rq),
            (target_rq, target_rp),
        ]

        by_vector = defaultdict(list)

        for x in pool:
            by_vector[residues[x][:prefix]].append(x)

        for vec_p, vec_q in orientations:
            for p in by_vector.get(vec_p, ()):
                for q in by_vector.get(vec_q, ()):

                    if not pair_is_valid(p, q, target_p, target_q):
                        continue

                    if verify_collision(
                        name,
                        p,
                        q,
                        target_p,
                        target_q,
                        prefix,
                        residues,
                        orders,
                    ):
                        return p, q

        return None

    # ------------------------------------------------------------------------
    # ORDER PAIR
    # ------------------------------------------------------------------------

    if name == "order_pair":
        target_op = orders[target_p][:prefix]
        target_oq = orders[target_q][:prefix]

        by_order = defaultdict(list)

        for x in pool:
            by_order[orders[x][:prefix]].append(x)

        for p in by_order.get(target_op, ()):
            for q in by_order.get(target_oq, ()):

                if not pair_is_valid(p, q, target_p, target_q):
                    continue

                if verify_collision(
                    name,
                    p,
                    q,
                    target_p,
                    target_q,
                    prefix,
                    residues,
                    orders,
                ):
                    return p, q

        return None

    # ------------------------------------------------------------------------
    # ORDER SORTED
    # ------------------------------------------------------------------------

    if name == "order_sorted":
        target_op = orders[target_p][:prefix]
        target_oq = orders[target_q][:prefix]

        by_order = defaultdict(list)

        for x in pool:
            by_order[orders[x][:prefix]].append(x)

        orientations = [
            (target_op, target_oq),
            (target_oq, target_op),
        ]

        for vec_p, vec_q in orientations:
            for p in by_order.get(vec_p, ()):
                for q in by_order.get(vec_q, ()):

                    if not pair_is_valid(p, q, target_p, target_q):
                        continue

                    if verify_collision(
                        name,
                        p,
                        q,
                        target_p,
                        target_q,
                        prefix,
                        residues,
                        orders,
                    ):
                        return p, q

        return None

    # ------------------------------------------------------------------------
    # ALGEBRAIC SIGNATURES
    # ------------------------------------------------------------------------

    # Build a residue-vector -> primes index.
    by_residue = defaultdict(list)

    for x in pool:
        by_residue[residues[x][:prefix]].append(x)

    target_res = residues[target_p][:prefix]
    target_q_res = residues[target_q][:prefix]

    for p in pool:

        rp = residues[p][:prefix]

        required = []

        possible = True

        for i, m in enumerate(MODULI[:prefix]):

            a = rp[i]

            if name == "Fsum":
                b = (target_res[i] + target_q_res[i] - a) % m

            elif name == "Fdiff":
                # Target is p-q.
                b = (a - (
                    target_res[i] - target_q_res[i]
                )) % m

            elif name == "Fprod":
                target_product = (
                    target_res[i] * target_q_res[i]
                ) % m

                if math.gcd(a, m) != 1:
                    possible = False
                    break

                inv = pow(a, -1, m)
                b = (target_product * inv) % m

            elif name == "cube_sum":
                target_cube = (
                    pow(target_res[i], 3, m)
                    + pow(target_q_res[i], 3, m)
                ) % m

                # We cannot directly invert the cube map in general.
                # Instead, build a cube-residue index below.
                required = None
                break

            else:
                raise ValueError(name)

            required.append(b)

        # Cube-sum gets a separate index.
        if name == "cube_sum":
            continue

        required_key = tuple(required)

        for q in by_residue.get(required_key, ()):

            if not pair_is_valid(p, q, target_p, target_q):
                continue

            if verify_collision(
                name,
                p,
                q,
                target_p,
                target_q,
                prefix,
                residues,
                orders,
            ):
                return p, q

    # ------------------------------------------------------------------------
    # CUBE SUM
    # ------------------------------------------------------------------------

    if name == "cube_sum":

        cube_index = defaultdict(list)

        for x in pool:
            key = tuple(
                pow(a, 3, m)
                for a, m in zip(
                    residues[x][:prefix],
                    MODULI[:prefix],
                )
            )

            cube_index[key].append(x)

        target_cube = tuple(
            (
                pow(a, 3, m)
                + pow(b, 3, m)
            ) % m
            for a, b, m in zip(
                target_res,
                target_q_res,
                MODULI[:prefix],
            )
        )

        for p in pool:

            pcube = tuple(
                pow(a, 3, m)
                for a, m in zip(
                    residues[p][:prefix],
                    MODULI[:prefix],
                )
            )

            required = tuple(
                (tc - pc) % m
                for tc, pc, m in zip(
                    target_cube,
                    pcube,
                    MODULI[:prefix],
                )
            )

            for q in cube_index.get(required, ()):

                if not pair_is_valid(p, q, target_p, target_q):
                    continue

                if verify_collision(
                    name,
                    p,
                    q,
                    target_p,
                    target_q,
                    prefix,
                    residues,
                    orders,
                ):
                    return p, q

    return None


# ============================================================================
# THEORETICAL FPAIR BASELINE
# ============================================================================

def theoretical_log10_baseline(prefix: int) -> float:
    total = 0.0

    for m in MODULI[:prefix]:
        units = phi(m)
        total -= math.log10(units)

    return total


# ============================================================================
# TARGET INFORMATION
# ============================================================================

def print_target_info(
    target_number: int,
    p: int,
    q: int,
):

    n = p * q
    s = p + q
    delta = (p - q) ** 2

    print(f"TARGET {target_number}")
    print("=" * 78)
    print(f"p       = {p}")
    print(f"q       = {q}")
    print(f"n       = {n}")
    print(f"s       = {s}")
    print(f"Delta   = {delta}")
    print(f"n bits  = {n.bit_length()}")
    print()


# ============================================================================
# MAIN
# ============================================================================

def main():

    random.seed(SEED)

    experiment_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 26")
    print("HASH-INDEXED CONSTRUCTIVE SIGNATURE COLLISION DEPTH")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print(f"random seed      = {SEED}")
    print(f"R values         = {R_VALUES}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {PRIME_POOL_SIZE:,}")
    print(f"prime interval   = [{PRIME_LOW:,}, {PRIME_HIGH:,}]")
    print(f"prefixes         = {PREFIXES}")
    print()

    # ------------------------------------------------------------------------
    # 1. MODULUS INVENTORY
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r, m, units in modulus_inventory():
        print(
            f"r={r:3d} "
            f"m={m:6d} "
            f"units={units:7d}"
        )

    print()

    # ------------------------------------------------------------------------
    # 2. THEORETICAL BASELINES
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("2. THEORETICAL FPAIR PREFIX BASELINES")
    print("=" * 78)

    print(
        f"{'prefix':>8} "
        f"{'last m':>10} "
        f"{'log10 generic freq':>24}"
    )

    for prefix in PREFIXES:
        print(
            f"{prefix:8d} "
            f"{MODULI[prefix - 1]:10d} "
            f"{theoretical_log10_baseline(prefix):24.6f}"
        )

    print()

    # ------------------------------------------------------------------------
    # 3. PRIME POOL
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. PRIME POOL GENERATION")
    print("=" * 78)

    t0 = time.perf_counter()

    print("Generating unique prime pool...")

    pool = generate_prime_pool(
        PRIME_POOL_SIZE,
        PRIME_LOW,
        PRIME_HIGH,
        SEED,
    )

    print(f"generated primes = {len(pool):,}")
    print(
        f"pool generation time = "
        f"{time.perf_counter() - t0:.2f}s"
    )
    print()

    # ------------------------------------------------------------------------
    # 4. RESIDUE CACHE
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("4. PRIME RESIDUE CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    residues = build_residue_cache(pool)

    print("residue cache complete")
    print(
        f"cache time = "
        f"{time.perf_counter() - t0:.2f}s"
    )
    print()

    # ------------------------------------------------------------------------
    # 5. ORDER CACHE
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. MULTIPLICATIVE ORDER CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    orders = build_order_cache(pool)

    print("order cache complete")
    print(
        f"cache time = "
        f"{time.perf_counter() - t0:.2f}s"
    )
    print()

    # ------------------------------------------------------------------------
    # 6. CONSTRUCTIVE SEARCH
    # ------------------------------------------------------------------------

    totals = {
        name: {
            prefix: 0
            for prefix in PREFIXES
        }
        for name in SIGNATURES
    }

    survival = {
        name: {
            prefix: 0
            for prefix in PREFIXES
        }
        for name in SIGNATURES
    }

    for target_number, (target_p, target_q) in enumerate(
        TARGET_PAIRS[:TARGETS],
        start=1,
    ):

        print_target_info(
            target_number,
            target_p,
            target_q,
        )

        print("CONSTRUCTIVE COLLISION SEARCH")
        print("-" * 78)

        for prefix in PREFIXES:

            print()
            print(
                f"PREFIX {prefix} "
                f"(m={MODULI[prefix - 1]})"
            )

            for name in SIGNATURES:

                t0 = time.perf_counter()

                candidate = find_collision(
                    name,
                    prefix,
                    target_p,
                    target_q,
                    pool,
                    residues,
                    orders,
                )

                elapsed = time.perf_counter() - t0

                if candidate is None:
                    result = "NONE FOUND"
                else:
                    result = f"{candidate}"

                    totals[name][prefix] += 1
                    survival[name][prefix] += 1

                print(
                    f"  {name:14s} "
                    f"{result:32s} "
                    f"[{elapsed:.3f}s]"
                )

        print()

    # ------------------------------------------------------------------------
    # 7. SUMMARY
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("7. COLLISION SURVIVAL SUMMARY")
    print("=" * 78)

    print()
    print(
        f"{'signature':16s}",
        end="",
    )

    for prefix in PREFIXES:
        print(
            f"{prefix:>10d}",
            end="",
        )

    print()

    print("-" * 78)

    for name in SIGNATURES:

        print(f"{name:16s}", end="")

        for prefix in PREFIXES:
            print(
                f"{survival[name][prefix]:>7d}/{TARGETS:<2d}",
                end="",
            )

        print()

    print()

    # ------------------------------------------------------------------------
    # 8. FIRST DISAPPEARANCE
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("8. FIRST OBSERVED COLLISION DISAPPEARANCE")
    print("=" * 78)

    for name in SIGNATURES:

        first_zero = None

        for prefix in PREFIXES:
            if survival[name][prefix] == 0:
                first_zero = prefix
                break

        if first_zero is None:
            text = "not reached"
        else:
            text = str(first_zero)

        print(
            f"{name:16s} "
            f"first zero-survival prefix = {text}"
        )

    print()

    # ------------------------------------------------------------------------
    # 9. INTERPRETATION
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("9. INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 26 extends Experiment 25 from Fpair-only constructive
search to the complete signature family.

For each prefix, the program attempts to produce an explicit wrong
prime pair having exactly the same signature as the target.

The important distinction is:

    random collision
        versus
    explicit constructive collision

A constructive collision gives an actual pair (a,b), rather than
just estimating how often collisions occur.

The search is hash-indexed and avoids enumerating all P^2 prime pairs.

For Fsum, Fprod, Fdiff and cube_sum, the residue vector of the first
prime determines the required residue vector of the second prime.

For Fpair and the ordering signatures, the individual residue/order
classes are indexed directly.

Therefore the expensive part is approximately proportional to the
prime-pool size times the number of prefixes, rather than to the
number of possible prime pairs.

Interpretation of NONE FOUND:

    NONE FOUND != proof of uniqueness.

It means that this finite prime pool contains no explicit collision
for that signature and prefix.

The key empirical quantity is the survival curve:

    prefix 3 -> prefix 4 -> prefix 5 -> prefix 6 -> prefix 7

If Fpair loses all constructive collisions at prefix 6 while a
compressed signature still has surviving collisions, that would be
evidence that Fpair retains more discriminating information.

If several signatures disappear at the same prefix, that suggests
the modulus-prefix itself is dominating the collision depth.

The experiment is therefore designed to distinguish:

    signature information
from
    modulus-prefix information.
"""
    )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)

    elapsed_total = time.perf_counter() - experiment_start

    print(
        f"total runtime = {elapsed_total:.2f}s "
        f"({elapsed_total / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()