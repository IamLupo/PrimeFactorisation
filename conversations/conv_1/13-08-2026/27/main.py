#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 27
CROSS-SIGNATURE COLLISION SEPARATION
FIXED ORDER CACHE

NO CSV OUTPUT
==============================================================================

Purpose
-------
Experiment 27 asks a different question from Experiments 25/26.

For each target factor pair (p,q), and for each modulus prefix, we find
constructive collisions for the individual signatures:

    Fpair
    Fsorted
    Fsum
    Fprod
    Fdiff
    cube_sum
    order_pair
    order_sorted

Then every found collision is tested against ALL other signatures.

This distinguishes:

    collision in signature A
        from
    collision simultaneously surviving signature A + B

The latter is the important "cross-signature" collision.

Important implementation detail
--------------------------------
The multiplicative-order cache is lazy.

The previous Experiment 27 crashed because order_cache contained only
pool primes while target primes were queried as well:

    KeyError: 2280911

This version computes an order on demand and caches it, so targets and
pool elements are handled uniformly.

==============================================================================
"""

from __future__ import annotations

import math
import random
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 20260814

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47
]

TARGETS = 12

POOL_SIZE = 5_000
PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

# Experiment 27 modulus family:
#
#     m = r^2 + r + 1
#
# This gives:
#
#     2 -> 7
#     3 -> 13
#     5 -> 31
#     7 -> 57
#     ...
#
MODULUS_FORMULA = lambda r: r * r + r + 1


SIGNATURE_NAMES = [
    "Fpair",
    "Fsorted",
    "Fsum",
    "Fprod",
    "Fdiff",
    "cube_sum",
    "order_pair",
    "order_sorted",
]


# These are the exact target factors used in Experiments 25/26.
#
# Keeping them fixed makes cross-experiment comparison meaningful.
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


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q

    @property
    def delta(self) -> int:
        return (self.p - self.q) ** 2


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin for the relevant integer range."""

    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    if n in small_primes:
        return True

    for p in small_primes:
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # More than enough for < 2^64.
    bases = (
        2, 325, 9375, 28178,
        450775, 9780504, 1795265022
    )

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = pow(x, 2, n)

            if x == n - 1:
                break
        else:
            return False

    return True


def generate_prime_pool(
    size: int,
    lo: int,
    hi: int,
    seed: int,
) -> List[int]:
    """
    Generate a deterministic unique prime pool.

    Candidates are sampled randomly from the interval. This preserves the
    experiment's random-seed behaviour without constructing a huge sieve.
    """

    rng = random.Random(seed)

    result = set()

    while len(result) < size:
        x = rng.randint(lo, hi)

        if x in result:
            continue

        if is_prime(x):
            result.add(x)

    return sorted(result)


def factor_integer(n: int) -> Dict[int, int]:
    """Small trial-division factorisation; moduli here are tiny."""

    factors: Dict[int, int] = {}

    d = 2

    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d

        if d == 2:
            d = 3
        else:
            d += 2

    if n > 1:
        factors[n] = factors.get(n, 0) + 1

    return factors


def phi(n: int) -> int:
    result = n

    for p in factor_integer(n):
        result -= result // p

    return result


def multiplicative_order(
    a: int,
    modulus: int,
    phi_cache: Dict[int, int],
    factor_cache: Dict[int, Dict[int, int]],
) -> Optional[int]:
    """
    Multiplicative order of a modulo modulus.

    Returns None when gcd(a, modulus) != 1.
    """

    a %= modulus

    if math.gcd(a, modulus) != 1:
        return None

    if modulus not in phi_cache:
        phi_cache[modulus] = phi(modulus)

    order = phi_cache[modulus]

    if modulus not in factor_cache:
        factor_cache[modulus] = factor_integer(order)

    for prime_factor in factor_cache[modulus]:
        while order % prime_factor == 0:
            candidate = order // prime_factor

            if pow(a, candidate, modulus) == 1:
                order = candidate
            else:
                break

    return order


# =============================================================================
# LAZY ORDER CACHE
# =============================================================================

class OrderCache:
    """
    Lazy multiplicative-order cache.

    Crucially, this accepts ANY prime/integer, including target factors.

    This fixes the Experiment 27 KeyError.
    """

    def __init__(self, moduli: Sequence[int]):
        self.moduli = tuple(moduli)

        self.orders: Dict[int, Dict[int, Optional[int]]] = {}

        self.phi_cache: Dict[int, int] = {}
        self.factor_cache: Dict[int, Dict[int, int]] = {}

    def get(
        self,
        value: int,
        r_values: Sequence[int],
    ) -> Dict[int, Optional[int]]:
        if value not in self.orders:
            self.orders[value] = {}

        row = self.orders[value]

        for r, modulus in zip(r_values, self.moduli):
            if r not in row:
                row[r] = multiplicative_order(
                    value,
                    modulus,
                    self.phi_cache,
                    self.factor_cache,
                )

        return row

    def get_one(
        self,
        value: int,
        r: int,
        modulus: int,
    ) -> Optional[int]:

        if value not in self.orders:
            self.orders[value] = {}

        row = self.orders[value]

        if r not in row:
            row[r] = multiplicative_order(
                value,
                modulus,
                self.phi_cache,
                self.factor_cache,
            )

        return row[r]


# =============================================================================
# RESIDUE CACHE
# =============================================================================

class ResidueCache:
    def __init__(
        self,
        pool: Sequence[int],
        moduli: Sequence[int],
        r_values: Sequence[int],
    ):
        self.data: Dict[int, Dict[int, int]] = {}

        for p in pool:
            self.data[p] = {
                r: p % m
                for r, m in zip(r_values, moduli)
            }

    def get(
        self,
        value: int,
        r_values: Sequence[int],
        moduli: Sequence[int],
    ) -> Dict[int, int]:

        if value not in self.data:
            self.data[value] = {
                r: value % m
                for r, m in zip(r_values, moduli)
            }

        return self.data[value]


# =============================================================================
# SIGNATURE FUNCTIONS
# =============================================================================

def fpair(
    p: int,
    q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(p, r_values, moduli)
    rq = residue_cache.get(q, r_values, moduli)

    return tuple(
        (rp[r], rq[r])
        for r in prefix
    )


def fsorted(
    p: int,
    q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(p, r_values, moduli)
    rq = residue_cache.get(q, r_values, moduli)

    return tuple(
        tuple(sorted((rp[r], rq[r])))
        for r in prefix
    )


def fsum(
    p: int,
    q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(p, r_values, moduli)
    rq = residue_cache.get(q, r_values, moduli)

    return tuple(
        (rp[r] + rq[r]) % m
        for r, m in zip(prefix, moduli)
    )


def fprod(
    p: int,
    q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(p, r_values, moduli)
    rq = residue_cache.get(q, r_values, moduli)

    return tuple(
        (rp[r] * rq[r]) % m
        for r, m in zip(prefix, moduli)
    )


def fdiff(
    p: int,
    q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(p, r_values, moduli)
    rq = residue_cache.get(q, r_values, moduli)

    return tuple(
        (rp[r] - rq[r]) % m
        for r, m in zip(prefix, moduli)
    )


def cube_sum(
    p: int,
    q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(p, r_values, moduli)
    rq = residue_cache.get(q, r_values, moduli)

    return tuple(
        (pow(rp[r], 3, m) + pow(rq[r], 3, m)) % m
        for r, m in zip(prefix, moduli)
    )


def order_pair(
    p: int,
    q: int,
    prefix: Sequence[int],
    order_cache: OrderCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    return tuple(
        (
            order_cache.get_one(p, r, m),
            order_cache.get_one(q, r, m),
        )
        for r, m in zip(prefix, moduli)
    )


def order_sorted(
    p: int,
    q: int,
    prefix: Sequence[int],
    order_cache: OrderCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    return tuple(
        tuple(sorted((
            order_cache.get_one(p, r, m),
            order_cache.get_one(q, r, m),
        ), key=lambda x: (-1 if x is None else x)))
        for r, m in zip(prefix, moduli)
    )


# =============================================================================
# GENERIC SIGNATURE DISPATCH
# =============================================================================

def signature(
    name: str,
    p: int,
    q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    order_cache: OrderCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    if name == "Fpair":
        return fpair(
            p, q, prefix,
            residue_cache, moduli, r_values
        )

    if name == "Fsorted":
        return fsorted(
            p, q, prefix,
            residue_cache, moduli, r_values
        )

    if name == "Fsum":
        return fsum(
            p, q, prefix,
            residue_cache, moduli, r_values
        )

    if name == "Fprod":
        return fprod(
            p, q, prefix,
            residue_cache, moduli, r_values
        )

    if name == "Fdiff":
        return fdiff(
            p, q, prefix,
            residue_cache, moduli, r_values
        )

    if name == "cube_sum":
        return cube_sum(
            p, q, prefix,
            residue_cache, moduli, r_values
        )

    if name == "order_pair":
        return order_pair(
            p, q, prefix,
            order_cache, moduli, r_values
        )

    if name == "order_sorted":
        return order_sorted(
            p, q, prefix,
            order_cache, moduli, r_values
        )

    raise ValueError(f"Unknown signature: {name}")


# =============================================================================
# HASH INDEXES
# =============================================================================

def build_pair_index(
    pool: Sequence[int],
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    """
    Index individual residue vectors.

    Used by Fpair/Fsorted.
    """

    index = defaultdict(list)

    for p in pool:
        rp = residue_cache.get(
            p, r_values, moduli
        )

        key = tuple(rp[r] for r in prefix)

        index[key].append(p)

    return index


def build_single_residue_index(
    pool: Sequence[int],
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    """
    Same physical structure as pair index, retained separately for clarity.
    """

    return build_pair_index(
        pool,
        prefix,
        residue_cache,
        moduli,
        r_values,
    )


# =============================================================================
# CONSTRUCTIVE COLLISION SEARCH
# =============================================================================

def find_fpair_collision(
    target_p: int,
    target_q: int,
    pool: Sequence[int],
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    """
    Find another pair with exactly the same ordered residue vectors.
    """

    target_p_vec = tuple(
        residue_cache.get(
            target_p, r_values, moduli
        )[r]
        for r in prefix
    )

    target_q_vec = tuple(
        residue_cache.get(
            target_q, r_values, moduli
        )[r]
        for r in prefix
    )

    by_vec = defaultdict(list)

    for x in pool:
        vec = tuple(
            residue_cache.get(
                x, r_values, moduli
            )[r]
            for r in prefix
        )

        by_vec[vec].append(x)

    p_candidates = by_vec.get(target_p_vec, [])
    q_candidates = by_vec.get(target_q_vec, [])

    for a in p_candidates:
        for b in q_candidates:
            if a == b:
                continue

            if {a, b} == {target_p, target_q}:
                continue

            return (a, b)

    return None


def find_fsorted_collision(
    target_p: int,
    target_q: int,
    pool: Sequence[int],
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    """
    Find another pair whose per-modulus residue pairs match as unordered pairs.
    """

    target_key = fsorted(
        target_p,
        target_q,
        prefix,
        residue_cache,
        moduli,
        r_values,
    )

    # Build signatures of pool elements lazily through a vector index.
    # For only 5000 primes this remains inexpensive.
    vectors = {}

    for x in pool:
        vectors[x] = tuple(
            residue_cache.get(
                x, r_values, moduli
            )[r]
            for r in prefix
        )

    # We need pairs, but avoid P^2 by grouping elements by their complete
    # residue vector. Candidate residue-vector pairs are inferred from target.
    target_p_vec = vectors.get(
        target_p,
        tuple(
            residue_cache.get(
                target_p, r_values, moduli
            )[r]
            for r in prefix
        ),
    )

    target_q_vec = vectors.get(
        target_q,
        tuple(
            residue_cache.get(
                target_q, r_values, moduli
            )[r]
            for r in prefix
        ),
    )

    by_vec = defaultdict(list)

    for x, vec in vectors.items():
        by_vec[vec].append(x)

    candidates = []

    for a in by_vec.get(target_p_vec, []):
        for b in by_vec.get(target_q_vec, []):
            if a == b:
                continue

            if {a, b} == {target_p, target_q}:
                continue

            candidates.append((a, b))

    for a in by_vec.get(target_q_vec, []):
        for b in by_vec.get(target_p_vec, []):
            if a == b:
                continue

            if {a, b} == {target_p, target_q}:
                continue

            candidates.append((a, b))

    for pair in candidates:
        if fsorted(
            pair[0],
            pair[1],
            prefix,
            residue_cache,
            moduli,
            r_values,
        ) == target_key:
            return pair

    return None


def derive_sum_pair_key(
    target_p: int,
    target_q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(target_p, r_values, moduli)
    rq = residue_cache.get(target_q, r_values, moduli)

    return tuple(
        (rp[r] + rq[r]) % m
        for r, m in zip(prefix, moduli)
    )


def derive_prod_pair_key(
    target_p: int,
    target_q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(target_p, r_values, moduli)
    rq = residue_cache.get(target_q, r_values, moduli)

    return tuple(
        (rp[r] * rq[r]) % m
        for r, m in zip(prefix, moduli)
    )


def derive_diff_pair_key(
    target_p: int,
    target_q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(target_p, r_values, moduli)
    rq = residue_cache.get(target_q, r_values, moduli)

    return tuple(
        (rp[r] - rq[r]) % m
        for r, m in zip(prefix, moduli)
    )


def derive_cube_sum_key(
    target_p: int,
    target_q: int,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    rp = residue_cache.get(target_p, r_values, moduli)
    rq = residue_cache.get(target_q, r_values, moduli)

    return tuple(
        (
            pow(rp[r], 3, m)
            + pow(rq[r], 3, m)
        ) % m
        for r, m in zip(prefix, moduli)
    )


def find_derived_collision(
    name: str,
    target_p: int,
    target_q: int,
    pool: Sequence[int],
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    order_cache: OrderCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    """
    Generic constructive search for Fsum/Fprod/Fdiff/cube_sum.

    We index the first prime and derive the required residue of the second
    prime for every modulus in the prefix.

    For a small pool, direct verification is cheap.
    """

    target_sig = signature(
        name,
        target_p,
        target_q,
        prefix,
        residue_cache,
        order_cache,
        moduli,
        r_values,
    )

    index = defaultdict(list)

    for x in pool:
        sx = signature(
            name,
            x,
            target_q,
            prefix,
            residue_cache,
            order_cache,
            moduli,
            r_values,
        )

        index[sx].append(x)

    # First try replacing p while retaining q.
    for a in index.get(target_sig, []):
        if a == target_p:
            continue

        if a == target_q:
            continue

        return (a, target_q)

    # General fallback: construct a signature map for pool elements.
    #
    # We deliberately avoid P^2. For each possible first prime we only need
    # to look up a complementary second-prime signature.

    element_keys = {}

    for x in pool:
        rx = residue_cache.get(
            x, r_values, moduli
        )

        if name == "Fsum":
            key = tuple(
                rx[r]
                for r in prefix
            )

        elif name == "Fprod":
            key = tuple(
                rx[r]
                for r in prefix
            )

        elif name == "Fdiff":
            key = tuple(
                rx[r]
                for r in prefix
            )

        elif name == "cube_sum":
            key = tuple(
                rx[r]
                for r in prefix
            )

        else:
            key = None

        element_keys[x] = key

    # Use a simple pair-index fallback.  This remains O(P) candidate checks
    # because the actual signature equality is evaluated only for candidates
    # from compatible residue classes.

    for a in pool:

        for b in pool:

            # Do NOT perform a quadratic scan.
            # This branch is intentionally disabled.
            break

        break

    # The robust route for these signatures is an indexed lookup using the
    # complete signature itself. Build pair signatures only among buckets
    # implied by the target. Since each residue vector is usually unique at
    # higher prefixes, this is very small.

    target_rp = residue_cache.get(
        target_p, r_values, moduli
    )
    target_rq = residue_cache.get(
        target_q, r_values, moduli
    )

    by_vector = defaultdict(list)

    for x in pool:
        rx = residue_cache.get(
            x, r_values, moduli
        )

        vec = tuple(
            rx[r]
            for r in prefix
        )

        by_vector[vec].append(x)

    for a in pool:
        ra = residue_cache.get(
            a, r_values, moduli
        )

        avec = tuple(
            ra[r]
            for r in prefix
        )

        if name == "Fsum":
            required = tuple(
                (target_rp[r] + target_rq[r] - ra[r]) % m
                for r, m in zip(prefix, moduli)
            )

        elif name == "Fprod":
            required = []

            possible = True

            for r, m in zip(prefix, moduli):
                x = ra[r]
                target = (
                    target_rp[r] * target_rq[r]
                ) % m

                # Need x*y = target mod m.
                # If x is not invertible, enumerate compatible residues.
                if math.gcd(x, m) != 1:
                    possible = False
                    break

                inv = pow(x, -1, m)
                required.append((target * inv) % m)

            if not possible:
                continue

            required = tuple(required)

        elif name == "Fdiff":
            required = tuple(
                (target_rp[r] - target_rq[r] + ra[r]) % m
                for r, m in zip(prefix, moduli)
            )

        elif name == "cube_sum":
            required = []

            possible = True

            for r, m in zip(prefix, moduli):
                target = (
                    pow(target_rp[r], 3, m)
                    + pow(target_rq[r], 3, m)
                ) % m

                # Need y^3 = target - x^3 mod m.
                # There can be multiple cube roots, so use the direct
                # compatibility scan inside the residue bucket only.
                required.append(
                    (target - pow(ra[r], 3, m)) % m
                )

            required = tuple(required)

            # Cube roots cannot be inferred by inversion in general.
            # Handle below using a global cube signature index.

        else:
            continue

        for b in by_vector.get(required, []):
            if a == b:
                continue

            if {a, b} == {target_p, target_q}:
                continue

            if signature(
                name,
                a,
                b,
                prefix,
                residue_cache,
                order_cache,
                moduli,
                r_values,
            ) == target_sig:
                return (a, b)

    # Special indexed fallback for cube_sum.
    if name == "cube_sum":

        cube_index = defaultdict(list)

        for x in pool:
            rx = residue_cache.get(
                x, r_values, moduli
            )

            key = tuple(
                pow(rx[r], 3, m)
                for r, m in zip(prefix, moduli)
            )

            cube_index[key].append(x)

        target_key = cube_sum(
            target_p,
            target_q,
            prefix,
            residue_cache,
            moduli,
            r_values,
        )

        for a in pool:
            ra = residue_cache.get(
                a, r_values, moduli
            )

            avec = tuple(
                pow(ra[r], 3, m)
                for r, m in zip(prefix, moduli)
            )

            required = tuple(
                (target - av) % m
                for target, av, m in zip(
                    target_key,
                    avec,
                    moduli,
                )
            )

            # Convert required sum-of-cubes residue to matching bucket.
            # The cube index is keyed by the actual cube residues.
            for b in cube_index.get(required, []):
                if a == b:
                    continue

                if {a, b} == {target_p, target_q}:
                    continue

                if signature(
                    name,
                    a,
                    b,
                    prefix,
                    residue_cache,
                    order_cache,
                    moduli,
                    r_values,
                ) == target_sig:
                    return (a, b)

    return None


def find_order_collision(
    name: str,
    target_p: int,
    target_q: int,
    pool: Sequence[int],
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    order_cache: OrderCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    target_sig = signature(
        name,
        target_p,
        target_q,
        prefix,
        residue_cache,
        order_cache,
        moduli,
        r_values,
    )

    # Index individual order vectors.
    index = defaultdict(list)

    for x in pool:
        if name == "order_pair":
            key = tuple(
                order_cache.get_one(
                    x,
                    r,
                    m,
                )
                for r, m in zip(prefix, moduli)
            )
        else:
            key = tuple(
                order_cache.get_one(
                    x,
                    r,
                    m,
                )
                for r, m in zip(prefix, moduli)
            )

        index[key].append(x)

    target_p_key = tuple(
        order_cache.get_one(
            target_p,
            r,
            m,
        )
        for r, m in zip(prefix, moduli)
    )

    target_q_key = tuple(
        order_cache.get_one(
            target_q,
            r,
            m,
        )
        for r, m in zip(prefix, moduli)
    )

    candidate_vectors = []

    if name == "order_pair":
        candidate_vectors = [
            (target_p_key, target_q_key),
        ]
    else:
        candidate_vectors = [
            (target_p_key, target_q_key),
            (target_q_key, target_p_key),
        ]

    for va, vb in candidate_vectors:

        for a in index.get(va, []):
            for b in index.get(vb, []):

                if a == b:
                    continue

                if {a, b} == {target_p, target_q}:
                    continue

                if signature(
                    name,
                    a,
                    b,
                    prefix,
                    residue_cache,
                    order_cache,
                    moduli,
                    r_values,
                ) == target_sig:
                    return (a, b)

    return None


def find_collision(
    name: str,
    target_p: int,
    target_q: int,
    pool: Sequence[int],
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    order_cache: OrderCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    if name == "Fpair":
        return find_fpair_collision(
            target_p,
            target_q,
            pool,
            prefix,
            residue_cache,
            moduli,
            r_values,
        )

    if name == "Fsorted":
        return find_fsorted_collision(
            target_p,
            target_q,
            pool,
            prefix,
            residue_cache,
            moduli,
            r_values,
        )

    if name in (
        "Fsum",
        "Fprod",
        "Fdiff",
        "cube_sum",
    ):
        return find_derived_collision(
            name,
            target_p,
            target_q,
            pool,
            prefix,
            residue_cache,
            order_cache,
            moduli,
            r_values,
        )

    if name in (
        "order_pair",
        "order_sorted",
    ):
        return find_order_collision(
            name,
            target_p,
            target_q,
            pool,
            prefix,
            residue_cache,
            order_cache,
            moduli,
            r_values,
        )

    raise ValueError(name)


# =============================================================================
# CROSS-SIGNATURE TEST
# =============================================================================

def cross_signature_matches(
    candidate: Tuple[int, int],
    target: Target,
    prefix: Sequence[int],
    residue_cache: ResidueCache,
    order_cache: OrderCache,
    moduli: Sequence[int],
    r_values: Sequence[int],
):
    """
    Compare one candidate against the target under every signature.

    Returns:
        dict signature -> bool
    """

    a, b = candidate

    result = {}

    for name in SIGNATURE_NAMES:

        target_sig = signature(
            name,
            target.p,
            target.q,
            prefix,
            residue_cache,
            order_cache,
            moduli,
            r_values,
        )

        candidate_sig = signature(
            name,
            a,
            b,
            prefix,
            residue_cache,
            order_cache,
            moduli,
            r_values,
        )

        result[name] = (
            target_sig == candidate_sig
        )

    return result


# =============================================================================
# FORMATTING
# =============================================================================

def pair_text(pair):
    if pair is None:
        return "NONE FOUND"

    return f"({pair[0]}, {pair[1]})"


def target_from_pair(pair) -> Target:
    return Target(pair[0], pair[1])


# =============================================================================
# MODULUS INVENTORY
# =============================================================================

def print_modulus_inventory():
    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r in R_VALUES:
        m = MODULUS_FORMULA(r)
        units = sum(
            1
            for x in range(1, m)
            if math.gcd(x, m) == 1
        )

        print(
            f"r={r:3d} "
            f"m={m:7d} "
            f"units={units:7d}"
        )

    print()


def print_theoretical_baselines():
    print("=" * 78)
    print("2. THEORETICAL FPAIR PREFIX BASELINES")
    print("=" * 78)

    # For an ordered pair of generic residues, approximate collision
    # probability as product of 1/m over the prefix.
    #
    # This is only a baseline, not a theorem about prime-restricted pairs.

    cumulative = 0.0

    for i, r in enumerate(R_VALUES[:7], start=1):

        m = MODULUS_FORMULA(r)

        cumulative += math.log10(m)

        print(
            f"{i:8d}"
            f"{m:12d}"
            f"{-cumulative:25.6f}"
        )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    total_start = time.perf_counter()

    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 27")
    print("CROSS-SIGNATURE COLLISION SEPARATION")
    print("FIXED ORDER CACHE")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print(f"random seed      = {SEED}")
    print(f"R values         = {R_VALUES}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {POOL_SIZE:,}")
    print(
        f"prime interval   = "
        f"[{PRIME_LO:,}, {PRIME_HI:,}]"
    )
    print(f"prefixes         = {PREFIXES}")
    print()

    # -------------------------------------------------------------------------
    # MODULI
    # -------------------------------------------------------------------------

    moduli_all = [
        MODULUS_FORMULA(r)
        for r in R_VALUES
    ]

    print_modulus_inventory()

    # -------------------------------------------------------------------------
    # BASELINES
    # -------------------------------------------------------------------------

    print_theoretical_baselines()

    # -------------------------------------------------------------------------
    # PRIME POOL
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("3. PRIME POOL GENERATION")
    print("=" * 78)

    t0 = time.perf_counter()

    print("Generating unique prime pool...")

    pool = generate_prime_pool(
        POOL_SIZE,
        PRIME_LO,
        PRIME_HI,
        SEED,
    )

    print(
        f"generated primes = {len(pool):,}"
    )

    print(
        f"pool generation time = "
        f"{time.perf_counter() - t0:.2f}s"
    )

    print()

    pool_set = set(pool)

    # -------------------------------------------------------------------------
    # TARGET VALIDATION
    # -------------------------------------------------------------------------

    targets = [
        target_from_pair(pair)
        for pair in TARGET_PAIRS[:TARGETS]
    ]

    missing_targets = []

    for t in targets:
        if t.p not in pool_set:
            missing_targets.append(t.p)

        if t.q not in pool_set:
            missing_targets.append(t.q)

    # The experiment deliberately uses fixed targets. If the randomly generated
    # pool does not contain one, we add the target factors to the searchable
    # population. This also guarantees that the target is always available.
    #
    # This is NOT a collision: target factors are excluded when evaluating
    # candidate pairs.

    if missing_targets:
        print(
            "NOTE: some fixed target factors were not in the generated pool."
        )
        print(
            "Adding target factors to the searchable population."
        )
        print()

        pool = sorted(
            set(pool).union(
                missing_targets
            )
        )

    # -------------------------------------------------------------------------
    # RESIDUE CACHE
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("4. PRIME RESIDUE CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    residue_cache = ResidueCache(
        pool,
        moduli_all,
        R_VALUES,
    )

    # IMPORTANT:
    # The cache can now also service target factors.
    for t in targets:
        residue_cache.get(
            t.p,
            R_VALUES,
            moduli_all,
        )

        residue_cache.get(
            t.q,
            R_VALUES,
            moduli_all,
        )

    print("residue cache complete")
    print(
        f"cache time = "
        f"{time.perf_counter() - t0:.2f}s"
    )
    print()

    # -------------------------------------------------------------------------
    # ORDER CACHE
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("5. MULTIPLICATIVE ORDER CACHE")
    print("=" * 78)

    print(
        "using lazy order cache "
        "(pool + target primes)"
    )

    t0 = time.perf_counter()

    order_cache = OrderCache(
        moduli_all
    )

    # Warm the cache for the pool. This is still fast enough for the current
    # experiment, but unlike the old version it also explicitly warms targets.
    for x in pool:
        order_cache.get(
            x,
            R_VALUES,
        )

    for t in targets:
        order_cache.get(
            t.p,
            R_VALUES,
        )
        order_cache.get(
            t.q,
            R_VALUES,
        )

    print("order cache complete")
    print(
        f"cache time = "
        f"{time.perf_counter() - t0:.2f}s"
    )
    print()

    # -------------------------------------------------------------------------
    # HASH INDEX INFORMATION
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("6. HASH INDEX CONSTRUCTION")
    print("=" * 78)

    for prefix_len in PREFIXES:

        prefix = R_VALUES[:prefix_len]
        moduli = moduli_all[:prefix_len]

        t0 = time.perf_counter()

        index = build_pair_index(
            pool,
            prefix,
            residue_cache,
            moduli,
            R_VALUES,
        )

        print(
            f"prefix {prefix_len}: "
            f"{len(index):,} residue classes "
            f"[{time.perf_counter() - t0:.3f}s]"
        )

    print()

    # -------------------------------------------------------------------------
    # CROSS-SIGNATURE EXPERIMENT
    # -------------------------------------------------------------------------

    survival = {
        name: {
            prefix_len: 0
            for prefix_len in PREFIXES
        }
        for name in SIGNATURE_NAMES
    }

    cross_survival = {
        name: {
            prefix_len: 0
            for prefix_len in PREFIXES
        }
        for name in SIGNATURE_NAMES
    }

    # Number of collisions that fool at least one OTHER signature.
    any_cross = {
        name: {
            prefix_len: 0
            for prefix_len in PREFIXES
        }
        for name in SIGNATURE_NAMES
    }

    # Number of collisions that fool ALL signatures.
    all_cross = {
        name: {
            prefix_len: 0
            for prefix_len in PREFIXES
        }
        for name in SIGNATURE_NAMES
    }

    # For each source signature, record the other signatures it accidentally
    # satisfies.
    cross_matrix = {
        source: {
            other: {
                prefix_len: 0
                for prefix_len in PREFIXES
            }
            for other in SIGNATURE_NAMES
            if other != source
        }
        for source in SIGNATURE_NAMES
    }

    # Store examples so the output contains concrete evidence.
    examples = {
        prefix_len: {}
        for prefix_len in PREFIXES
    }

    print("=" * 78)
    print("7. CROSS-SIGNATURE COLLISION SEARCH")
    print("=" * 78)
    print()

    for target_idx, target in enumerate(targets, start=1):

        print("TARGET", target_idx)
        print("=" * 78)

        print(f"p       = {target.p}")
        print(f"q       = {target.q}")
        print(f"n       = {target.n}")
        print(f"s       = {target.s}")
        print(f"Delta   = {target.delta}")
        print(
            f"n bits  = "
            f"{target.n.bit_length()}"
        )

        print()
        print("CROSS-SIGNATURE SEARCH")
        print("-" * 78)

        for prefix_len in PREFIXES:

            prefix = R_VALUES[:prefix_len]
            moduli = moduli_all[:prefix_len]

            print()
            print(
                f"PREFIX {prefix_len} "
                f"(m={moduli[-1]})"
            )

            examples[prefix_len].setdefault(
                target_idx,
                {}
            )

            found_for_target = {}

            for source_name in SIGNATURE_NAMES:

                t0 = time.perf_counter()

                candidate = find_collision(
                    source_name,
                    target.p,
                    target.q,
                    pool,
                    prefix,
                    residue_cache,
                    order_cache,
                    moduli,
                    R_VALUES,
                )

                elapsed = (
                    time.perf_counter()
                    - t0
                )

                if candidate is None:
                    print(
                        f"  {source_name:14s} "
                        f"NONE FOUND"
                        f"{'':>24}"
                        f"[{elapsed:.3f}s]"
                    )

                    found_for_target[
                        source_name
                    ] = None

                    continue

                survival[source_name][
                    prefix_len
                ] += 1

                matches = cross_signature_matches(
                    candidate,
                    target,
                    prefix,
                    residue_cache,
                    order_cache,
                    moduli,
                    R_VALUES,
                )

                found_for_target[
                    source_name
                ] = (
                    candidate,
                    matches,
                )

                matching_other = [
                    name
                    for name, matched in matches.items()
                    if name != source_name
                    and matched
                ]

                all_match = all(
                    matches.values()
                )

                if matching_other:
                    cross_survival[
                        source_name
                    ][prefix_len] += 1

                    any_cross[
                        source_name
                    ][prefix_len] += 1

                    for other in matching_other:
                        cross_matrix[
                            source_name
                        ][other][prefix_len] += 1

                if all_match:
                    all_cross[
                        source_name
                    ][prefix_len] += 1

                print(
                    f"  {source_name:14s} "
                    f"{pair_text(candidate):32s}"
                    f"[{elapsed:.3f}s]"
                )

                if matching_other:
                    print(
                        "      ALSO MATCHES: "
                        + ", ".join(
                            matching_other
                        )
                    )
                else:
                    print(
                        "      ALSO MATCHES: "
                        "none"
                    )

                if all_match:
                    print(
                        "      *** ALL SIGNATURES MATCH ***"
                    )

            examples[prefix_len][
                target_idx
            ] = found_for_target

        print()

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("8. SELF-COLLISION SURVIVAL")
    print("=" * 78)

    header = (
        f"{'signature':16s}"
        + "".join(
            f"{p:>12d}"
            for p in PREFIXES
        )
    )

    print(header)
    print("-" * 78)

    for name in SIGNATURE_NAMES:
        row = (
            f"{name:16s}"
            + "".join(
                f"{survival[name][p]:>12d}/{TARGETS}"
                for p in PREFIXES
            )
        )
        print(row)

    print()

    # -------------------------------------------------------------------------
    # CROSS COLLISION SUMMARY
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("9. COLLISIONS THAT ALSO FOOL ANOTHER SIGNATURE")
    print("=" * 78)

    print(header)
    print("-" * 78)

    for name in SIGNATURE_NAMES:
        row = (
            f"{name:16s}"
            + "".join(
                f"{any_cross[name][p]:>12d}/{TARGETS}"
                for p in PREFIXES
            )
        )
        print(row)

    print()

    # -------------------------------------------------------------------------
    # ALL-SIGNATURE SURVIVAL
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("10. COLLISIONS THAT FOOL ALL SIGNATURES")
    print("=" * 78)

    print(header)
    print("-" * 78)

    for name in SIGNATURE_NAMES:
        row = (
            f"{name:16s}"
            + "".join(
                f"{all_cross[name][p]:>12d}/{TARGETS}"
                for p in PREFIXES
            )
        )
        print(row)

    print()

    # -------------------------------------------------------------------------
    # CROSS MATRIX
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("11. CROSS-SIGNATURE SEPARATION MATRIX")
    print("=" * 78)

    print(
        "Each cell = number of source collisions that also satisfy"
    )
    print(
        "the destination signature."
    )
    print()

    for prefix_len in PREFIXES:

        print(
            f"PREFIX {prefix_len}"
        )
        print("-" * 78)

        print(
            f"{'source':16s}"
            + "".join(
                f"{other[:12]:>14s}"
                for other in SIGNATURE_NAMES
            )
        )

        for source in SIGNATURE_NAMES:

            row = f"{source:16s}"

            for other in SIGNATURE_NAMES:

                if source == other:
                    value = survival[
                        source
                    ][prefix_len]

                else:
                    value = cross_matrix[
                        source
                    ][other][prefix_len]

                row += f"{value:>14d}"

            print(row)

        print()

    # -------------------------------------------------------------------------
    # EXPLICIT EXAMPLES
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("12. EXPLICIT CROSS-COLLISION EXAMPLES")
    print("=" * 78)

    for prefix_len in PREFIXES:

        printed = False

        for target_idx in range(
            1,
            TARGETS + 1,
        ):

            data = examples[
                prefix_len
            ].get(target_idx, {})

            for source_name, item in data.items():

                if item is None:
                    continue

                candidate, matches = item

                others = [
                    name
                    for name, matched in matches.items()
                    if name != source_name
                    and matched
                ]

                if not others:
                    continue

                print(
                    f"prefix {prefix_len}, "
                    f"target {target_idx}, "
                    f"source={source_name}"
                )

                print(
                    f"  target    = "
                    f"({targets[target_idx-1].p}, "
                    f"{targets[target_idx-1].q})"
                )

                print(
                    f"  candidate = "
                    f"{candidate}"
                )

                print(
                    f"  also matches = "
                    f"{', '.join(others)}"
                )

                printed = True

                # Keep output manageable.
                if printed:
                    break

            if printed:
                break

        if not printed:
            print(
                f"prefix {prefix_len}: "
                "no cross-signature example found"
            )

        print()

    # -------------------------------------------------------------------------
    # FIRST ZERO SURVIVAL
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("13. FIRST ZERO SELF-COLLISION SURVIVAL")
    print("=" * 78)

    for name in SIGNATURE_NAMES:

        first_zero = None

        for p in PREFIXES:
            if survival[name][p] == 0:
                first_zero = p
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

    # -------------------------------------------------------------------------
    # INTERPRETATION
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("14. INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 27 measures CROSS-SIGNATURE separation.

Experiment 26 asked:

    Does signature X have a constructive collision?

Experiment 27 asks the stronger question:

    If signature X has a collision, does that same wrong pair
    also fool signature Y?

A useful distinction is:

    SELF COLLISION
        candidate matches the source signature only

    CROSS COLLISION
        candidate matches source signature and at least one
        different signature

    ALL-SIGNATURE COLLISION
        candidate matches every tested signature

The last two quantities are especially important.

If a source signature has many self-collisions but very few
cross-collisions, it is providing information that is not redundant
with the other signatures.

If almost every collision also fools several other signatures, then
those signatures are highly correlated for this modulus prefix.

NONE FOUND still does not prove uniqueness.

The experiment is finite and pool-dependent.

The most interesting comparison is therefore:

    self-survival
        versus
    cross-survival

across prefixes 3, 4, 5, 6 and 7.
"""
    )

    # -------------------------------------------------------------------------
    # RUNTIME
    # -------------------------------------------------------------------------

    total_runtime = (
        time.perf_counter()
        - total_start
    )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)
    print(
        f"total runtime = "
        f"{total_runtime:.2f}s "
        f"({total_runtime / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()