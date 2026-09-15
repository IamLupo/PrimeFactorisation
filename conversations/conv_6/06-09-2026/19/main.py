#!/usr/bin/env python3

"""
START EXPERIMENT 176

CROSS-CRT EXACT-PRODUCT BRIDGE

Core idea
---------

For an A-side residue:

    p == a (mod r)

and:

    p*q == n,

if gcd(a,r)=1 then:

    q == n*a^(-1) (mod r).

Likewise:

    q == b (mod s)

implies:

    p == n*b^(-1) (mod s).

Experiments 168-175 used the induced congruence mainly as a way
to construct t-residue conditions.

Experiment 176 uses it differently.

Suppose we have a partial A assignment:

    p == a0 (mod r0)
    p == a1 (mod r1)

Then we immediately know:

    q == c0 (mod r0)
    q == c1 (mod r1).

Those are genuine q congruences.

The B-side independently supplies:

    q == b0 (mod s0)
    q == b1 (mod s1)
    ...

Because every R1 modulus is coprime to every R2 modulus, we can
combine these into one CRT system for q.

If the combined modulus M exceeds sqrt(n), then there is at most
one q <= sqrt(n) represented by that residue class.

Therefore we can test q directly:

    q_candidate
        |
        v
    n % q_candidate
        |
        v
      exact factor

The same procedure is done in the other direction.

This is a genuine bidirectional use of:

    p*q = n

inside the residue search.

The search is deliberately NOT a generic 10-variable CSP.

It incrementally chooses only a few A residues and then uses
the cheapest B-domain CRT bridge.

Measurements:

    A partial states
    B partial states
    bridge CRT combinations
    candidate q/p values
    divisibility tests
    exact recoveries

FINISHED EXPERIMENT 176
"""

from __future__ import annotations

import itertools
import math
import random
import time

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from sympy import randprime


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS_LIST = [30, 36, 42, 48, 54]

SEED = 1762026

MAX_PARTIAL_STATES = 2_000_000
MAX_BRIDGE_COMBINATIONS = 2_000_000
MAX_EXACT_TESTS = 2_000_000


# ============================================================
# CRT
# ============================================================

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> Tuple[int, int]:

    g = math.gcd(
        m1,
        m2,
    )

    if g != 1:
        raise ValueError(
            f"CRT moduli not coprime: {m1}, {m2}"
        )

    t = (
        (a2 - a1)
        * pow(m1, -1, m2)
    ) % m2

    x = a1 + m1 * t
    M = m1 * m2

    return (
        x % M,
        M,
    )


def crt_many(
    residues: List[int],
    moduli: List[int],
) -> Tuple[int, int]:

    x = 0
    M = 1

    for a, m in zip(
        residues,
        moduli,
    ):

        x, M = crt_pair(
            x,
            M,
            a,
            m,
        )

    return x, M


# ============================================================
# C
# ============================================================

def exact_C(
    n: int,
    r: int,
    s: int,
    a: int,
    b: int,
) -> int:

    z = a * b

    beta = (
        (n - z)
        * pow(r, -1, s)
    ) % s

    alpha = (
        (n - z)
        * pow(s, -1, r)
    ) % r

    return (
        r * beta
        + s * alpha
        + z
    ) // (r * s)


# ============================================================
# C GRID
# ============================================================

def make_C_grid(
    n: int,
    p: int,
    q: int,
) -> List[List[int]]:

    A = [
        p % r
        for r in R1
    ]

    B = [
        q % s
        for s in R2
    ]

    return [
        [
            exact_C(
                n,
                r,
                s,
                A[i],
                B[j],
            )
            for j, s in enumerate(R2)
        ]
        for i, r in enumerate(R1)
    ]


# ============================================================
# RELATIONS
# ============================================================

@dataclass
class Relation:

    r: int
    s: int

    target_c: int

    pairs: Set[Tuple[int, int]]

    by_a: Dict[int, Set[int]]
    by_b: Dict[int, Set[int]]


def build_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
) -> Relation:

    pairs = set()

    by_a = {}
    by_b = {}

    # 0 is excluded because a prime factor is larger than
    # every radix used in this experiment.
    for a in range(1, r):

        for b in range(1, s):

            if exact_C(
                n,
                r,
                s,
                a,
                b,
            ) != target_c:

                continue

            pairs.add(
                (a, b)
            )

            by_a.setdefault(
                a,
                set(),
            ).add(b)

            by_b.setdefault(
                b,
                set(),
            ).add(a)

    return Relation(
        r=r,
        s=s,
        target_c=target_c,
        pairs=pairs,
        by_a=by_a,
        by_b=by_b,
    )


def build_relations(
    n: int,
    c_grid: List[List[int]],
) -> Dict[
    Tuple[int, int],
    Relation,
]:

    result = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            result[(i, j)] = build_relation(
                n,
                r,
                s,
                c_grid[i][j],
            )

    return result


# ============================================================
# INITIAL DOMAINS
# ============================================================

def make_A_domains() -> List[Set[int]]:

    return [
        set(range(1, r))
        for r in R1
    ]


def make_B_domains() -> List[Set[int]]:

    return [
        set(range(1, s))
        for s in R2
    ]


# ============================================================
# C PROPAGATION
# ============================================================

def propagate(
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> bool:

    changed = True

    while changed:

        changed = False

        for i in range(len(R1)):

            for j in range(len(R2)):

                rel = relations[
                    (i, j)
                ]

                new_A = {
                    a
                    for a in A[i]
                    if (
                        rel.by_a.get(
                            a,
                            set(),
                        )
                        & B[j]
                    )
                }

                if new_A != A[i]:

                    A[i] = new_A
                    changed = True

                    if not new_A:
                        return False

                new_B = {
                    b
                    for b in B[j]
                    if (
                        rel.by_b.get(
                            b,
                            set(),
                        )
                        & A[i]
                    )
                }

                if new_B != B[j]:

                    B[j] = new_B
                    changed = True

                    if not new_B:
                        return False

    return True


# ============================================================
# DIRECT ASSIGNMENT CHECK
# ============================================================

def A_value_supported(
    i: int,
    a: int,
    B: List[Set[int]],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> bool:

    for j in range(len(R2)):

        if not (
            relations[(i, j)]
            .by_a
            .get(a, set())
            & B[j]
        ):

            return False

    return True


def B_value_supported(
    j: int,
    b: int,
    A: List[Set[int]],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> bool:

    for i in range(len(R1)):

        if not (
            relations[(i, j)]
            .by_b
            .get(b, set())
            & A[i]
        ):

            return False

    return True


# ============================================================
# FACTOR-INDUCED CONGRUENCES
# ============================================================

def implied_q_residue(
    n: int,
    a: int,
    r: int,
) -> Optional[int]:

    a %= r

    if a == 0:
        return None

    return (
        n
        * pow(a, -1, r)
    ) % r


def implied_p_residue(
    n: int,
    b: int,
    s: int,
) -> Optional[int]:

    b %= s

    if b == 0:
        return None

    return (
        n
        * pow(b, -1, s)
    ) % s


# ============================================================
# CHEAPEST OPPOSITE DOMAIN SUBSETS
# ============================================================

def choose_bridge_subsets(
    fixed_moduli: List[int],
    opposite_domains: List[Set[int]],
    opposite_moduli: List[int],
    sqrt_n: int,
) -> List[
    Tuple[
        Tuple[int, ...],
        int,
        int,
    ]
]:

    """
    Return useful subsets of opposite domains.

    Each subset is:

        (indices, residue-combination count, total modulus)

    where the total modulus together with fixed_moduli exceeds sqrt(n).
    """

    result = []

    fixed_product = math.prod(
        fixed_moduli
    )

    if fixed_product > sqrt_n:

        result.append(
            (
                (),
                1,
                fixed_product,
            )
        )

    for mask in range(
        1,
        1 << len(opposite_moduli),
    ):

        indices = tuple(
            i
            for i in range(
                len(opposite_moduli)
            )
            if mask & (1 << i)
        )

        combo_count = 1

        modulus = fixed_product

        for i in indices:

            combo_count *= len(
                opposite_domains[i]
            )

            modulus *= (
                opposite_moduli[i]
            )

        if modulus <= sqrt_n:
            continue

        if combo_count > MAX_BRIDGE_COMBINATIONS:
            continue

        result.append(
            (
                indices,
                combo_count,
                modulus,
            )
        )

    result.sort(
        key=lambda x: (
            x[1],
            len(x[0]),
            x[2],
        )
    )

    return result


# ============================================================
# BRIDGE SEARCH
# ============================================================

@dataclass
class BridgeStats:

    partial_states: int = 0
    bridge_combinations: int = 0
    q_candidates: int = 0
    p_candidates: int = 0
    divisions: int = 0
    exact_hits: int = 0


def bridge_q_from_A_state(
    *,
    n: int,
    sqrt_n: int,
    A_indices: Tuple[int, ...],
    A_values: Tuple[int, ...],
    B_domains: List[Set[int]],
    stats: BridgeStats,
    c_grid: List[List[int]],
) -> Optional[Tuple[int, int]]:

    """
    A partial state gives q congruences modulo R1.

    Combine those with a small subset of B residues.

    Once modulus > sqrt(n), each CRT class has at most one q
    in [1,sqrt(n)].
    """

    q_residues = []
    q_moduli = []

    for i, a in zip(
        A_indices,
        A_values,
    ):

        r = R1[i]

        qr = implied_q_residue(
            n,
            a,
            r,
        )

        if qr is None:
            return None

        q_residues.append(qr)
        q_moduli.append(r)

    subsets = choose_bridge_subsets(
        q_moduli,
        B_domains,
        R2,
        sqrt_n,
    )

    for B_indices, _, bridge_modulus in subsets:

        if bridge_modulus <= sqrt_n:
            continue

        lists = [
            sorted(
                B_domains[j]
            )
            for j in B_indices
        ]

        for B_values in itertools.product(
            *lists
        ):

            stats.bridge_combinations += 1

            residues = list(
                q_residues
            )

            moduli = list(
                q_moduli
            )

            residues.extend(
                B_values
            )

            moduli.extend(
                R2[j]
                for j in B_indices
            )

            q0, M = crt_many(
                residues,
                moduli,
            )

            if M <= sqrt_n:
                continue

            if not (
                1
                < q0
                <= sqrt_n
            ):
                continue

            stats.q_candidates += 1
            stats.divisions += 1

            if n % q0 != 0:
                continue

            p0 = n // q0

            stats.exact_hits += 1

            if verify_factor(
                n,
                p0,
                q0,
                c_grid,
            ):

                return (
                    p0,
                    q0,
                )

    return None


def bridge_p_from_B_state(
    *,
    n: int,
    sqrt_n: int,
    B_indices: Tuple[int, ...],
    B_values: Tuple[int, ...],
    A_domains: List[Set[int]],
    stats: BridgeStats,
    c_grid: List[List[int]],
) -> Optional[Tuple[int, int]]:

    """
    Symmetric B -> p bridge.
    """

    p_residues = []
    p_moduli = []

    for j, b in zip(
        B_indices,
        B_values,
    ):

        s = R2[j]

        pr = implied_p_residue(
            n,
            b,
            s,
        )

        if pr is None:
            return None

        p_residues.append(pr)
        p_moduli.append(s)

    subsets = choose_bridge_subsets(
        p_moduli,
        A_domains,
        R1,
        sqrt_n,
    )

    for A_indices, _, bridge_modulus in subsets:

        if bridge_modulus <= sqrt_n:
            continue

        lists = [
            sorted(
                A_domains[i]
            )
            for i in A_indices
        ]

        for A_values in itertools.product(
            *lists
        ):

            stats.bridge_combinations += 1

            residues = list(
                p_residues
            )

            moduli = list(
                p_moduli
            )

            residues.extend(
                A_values
            )

            moduli.extend(
                R1[i]
                for i in A_indices
            )

            p0, M = crt_many(
                residues,
                moduli,
            )

            if M <= sqrt_n:
                continue

            if not (
                1
                < p0
                <= sqrt_n
            ):
                continue

            stats.p_candidates += 1
            stats.divisions += 1

            if n % p0 != 0:
                continue

            q0 = n // p0

            stats.exact_hits += 1

            if verify_factor(
                n,
                p0,
                q0,
                c_grid,
            ):

                return (
                    p0,
                    q0,
                )

    return None


# ============================================================
# ORDERING
# ============================================================

def order_variables(
    domains: List[Set[int]],
) -> List[int]:

    return sorted(
        range(len(domains)),
        key=lambda i: (
            len(domains[i]),
            -(
                R1[i]
            ),
            i,
        ),
    )


# ============================================================
# SEMIPRIME
# ============================================================

def make_semiprime(
    bits: int,
    rng: random.Random,
) -> Tuple[int, int, int]:

    lo = 1 << (
        bits // 2 - 1
    )

    hi = 1 << (
        bits // 2 + 1
    )

    p = int(
        randprime(
            lo,
            hi,
        )
    )

    q = int(
        randprime(
            lo,
            hi,
        )
    )

    while p == q:

        q = int(
            randprime(
                lo,
                hi,
            )
        )

    return (
        p,
        q,
        p * q,
    )


# ============================================================
# SAMPLE
# ============================================================

def run_sample(
    bits: int,
    sample_id: int,
    rng: random.Random,
) -> dict:

    p, q, n = make_semiprime(
        bits,
        rng,
    )

    sqrt_n = math.isqrt(n)

    true_A = [
        p % r
        for r in R1
    ]

    true_B = [
        q % s
        for s in R2
    ]

    c_grid = make_C_grid(
        n,
        p,
        q,
    )

    relations = build_relations(
        n,
        c_grid,
    )

    A = make_A_domains()
    B = make_B_domains()

    # Initial C propagation.
    start_prop = time.perf_counter()

    initial_ok = propagate(
        A,
        B,
        relations,
    )

    propagation_time = (
        time.perf_counter()
        - start_prop
    )

    initial_A = [
        set(x)
        for x in A
    ]

    initial_B = [
        set(x)
        for x in B
    ]

    A_order = order_variables(
        A
    )

    B_order = order_variables(
        B
    )

    true_survives = (
        initial_ok
        and all(
            true_A[i]
            in A[i]
            for i in range(len(R1))
        )
        and all(
            true_B[j]
            in B[j]
            for j in range(len(R2))
        )
    )

    # --------------------------------------------------------
    # Search state.
    # --------------------------------------------------------

    A_stats = BridgeStats()
    B_stats = BridgeStats()

    exact = None

    search_start = time.perf_counter()

    # --------------------------------------------------------
    # A-side partial CRT search.
    #
    # We deliberately search only a few A residues before
    # making bridge attempts.
    # --------------------------------------------------------

    A_prefix_limit = len(A_order)

    for depth in range(
        1,
        A_prefix_limit + 1,
    ):

        prefix = tuple(
            A_order[:depth]
        )

        A_lists = [
            sorted(
                A[i]
            )
            for i in prefix
        ]

        prefix_count = math.prod(
            len(x)
            for x in A_lists
        )

        if (
            A_stats.partial_states
            + prefix_count
            > MAX_PARTIAL_STATES
        ):
            break

        for values in itertools.product(
            *A_lists
        ):

            A_stats.partial_states += 1

            # Quick C support.
            good = True

            for i, a in zip(
                prefix,
                values,
            ):

                if not A_value_supported(
                    i,
                    a,
                    B,
                    relations,
                ):

                    good = False
                    break

            if not good:
                continue

            result = bridge_q_from_A_state(
                n=n,
                sqrt_n=sqrt_n,
                A_indices=prefix,
                A_values=values,
                B_domains=B,
                stats=A_stats,
                c_grid=c_grid,
            )

            if result is not None:

                exact = (
                    result[0],
                    result[1],
                    "A->q-bridge",
                )

                break

        if exact is not None:
            break

    # --------------------------------------------------------
    # Symmetric B-side search if needed.
    # --------------------------------------------------------

    if exact is None:

        for depth in range(
            1,
            len(B_order) + 1,
        ):

            prefix = tuple(
                B_order[:depth]
            )

            B_lists = [
                sorted(
                    B[j]
                )
                for j in prefix
            ]

            prefix_count = math.prod(
                len(x)
                for x in B_lists
            )

            if (
                B_stats.partial_states
                + prefix_count
                > MAX_PARTIAL_STATES
            ):
                break

            for values in itertools.product(
                *B_lists
            ):

                B_stats.partial_states += 1

                good = True

                for j, b in zip(
                    prefix,
                    values,
                ):

                    if not B_value_supported(
                        j,
                        b,
                        A,
                        relations,
                    ):

                        good = False
                        break

                if not good:
                    continue

                result = bridge_p_from_B_state(
                    n=n,
                    sqrt_n=sqrt_n,
                    B_indices=prefix,
                    B_values=values,
                    A_domains=A,
                    stats=B_stats,
                    c_grid=c_grid,
                )

                if result is not None:

                    exact = (
                        result[0],
                        result[1],
                        "B->p-bridge",
                    )

                    break

            if exact is not None:
                break

    search_time = (
        time.perf_counter()
        - search_start
    )

    total_time = (
        propagation_time
        + search_time
    )

    recovered = False

    if exact is not None:

        ep, eq, direction = exact

        recovered = (
            ep * eq == n
            and {
                ep,
                eq,
            }
            == {
                p,
                q,
            }
        )

    # --------------------------------------------------------
    # OUTPUT
    # --------------------------------------------------------

    print()
    print("=" * 90)
    print(
        f"EXPERIMENT 176 SAMPLE {sample_id}"
    )
    print("=" * 90)

    print(
        f"bits              = {bits}"
    )

    print(
        f"p                 = {p}"
    )

    print(
        f"q                 = {q}"
    )

    print(
        f"n                 = {n}"
    )

    print(
        f"sqrt(n)           = {sqrt_n}"
    )

    print()
    print(
        f"true A            = {true_A}"
    )

    print(
        f"true B            = {true_B}"
    )

    print()
    print(
        "C-grid:"
    )

    for row in c_grid:

        print(
            "   ",
            row,
        )

    print()
    print(
        "--- INITIAL C PROPAGATION ---"
    )

    print(
        f"initial OK                  = "
        f"{initial_ok}"
    )

    print(
        f"A widths                    = "
        f"{[len(x) for x in initial_A]}"
    )

    print(
        f"B widths                    = "
        f"{[len(x) for x in initial_B]}"
    )

    print(
        f"true residues survive       = "
        f"{true_survives}"
    )

    print(
        f"propagation time            = "
        f"{propagation_time:.6f}s"
    )

    print()
    print(
        "--- A -> Q BRIDGE ---"
    )

    print(
        f"partial states              = "
        f"{A_stats.partial_states}"
    )

    print(
        f"bridge combinations         = "
        f"{A_stats.bridge_combinations}"
    )

    print(
        f"q candidates                = "
        f"{A_stats.q_candidates}"
    )

    print(
        f"divisibility tests          = "
        f"{A_stats.divisions}"
    )

    print(
        f"exact hits                  = "
        f"{A_stats.exact_hits}"
    )

    print()
    print(
        "--- B -> P BRIDGE ---"
    )

    print(
        f"partial states              = "
        f"{B_stats.partial_states}"
    )

    print(
        f"bridge combinations         = "
        f"{B_stats.bridge_combinations}"
    )

    print(
        f"p candidates                = "
        f"{B_stats.p_candidates}"
    )

    print(
        f"divisibility tests          = "
        f"{B_stats.divisions}"
    )

    print(
        f"exact hits                  = "
        f"{B_stats.exact_hits}"
    )

    print()
    print(
        f"exact                       = "
        f"{exact}"
    )

    print(
        f"recovered                   = "
        f"{recovered}"
    )

    print()
    print(
        f"search time                 = "
        f"{search_time:.6f}s"
    )

    print(
        f"total time                  = "
        f"{total_time:.6f}s"
    )

    print("=" * 90)

    return {
        "bits": bits,
        "initial_A_widths": [
            len(x)
            for x in initial_A
        ],
        "initial_B_widths": [
            len(x)
            for x in initial_B
        ],
        "true_survives": true_survives,
        "A_partial_states": A_stats.partial_states,
        "A_bridge_combinations": A_stats.bridge_combinations,
        "A_candidates": A_stats.q_candidates,
        "A_exact_hits": A_stats.exact_hits,
        "B_partial_states": B_stats.partial_states,
        "B_bridge_combinations": B_stats.bridge_combinations,
        "B_candidates": B_stats.p_candidates,
        "B_exact_hits": B_stats.exact_hits,
        "exact": exact,
        "recovered": recovered,
        "total_time": total_time,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    rng = random.Random(
        SEED
    )

    print(
        "START EXPERIMENT 176"
    )

    print()

    print(
        "R1 =",
        R1,
    )

    print(
        "R2 =",
        R2,
    )

    print(
        "BITS =",
        BITS_LIST,
    )

    print(
        "MAX_PARTIAL_STATES =",
        MAX_PARTIAL_STATES,
    )

    print(
        "MAX_BRIDGE_COMBINATIONS =",
        MAX_BRIDGE_COMBINATIONS,
    )

    print()

    results = []

    for sample_id, bits in enumerate(
        BITS_LIST,
        start=1,
    ):

        results.append(
            run_sample(
                bits,
                sample_id,
                rng,
            )
        )

    print()
    print("=" * 125)
    print(
        "EXPERIMENT 176 SUMMARY"
    )
    print("=" * 125)

    print()

    print(
        "bits | A-partial | A-bridge | A-cand | "
        "B-partial | B-bridge | B-cand | "
        "exact | recovered | time"
    )

    print("-" * 125)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{r['A_partial_states']:>9} | "
            f"{r['A_bridge_combinations']:>8} | "
            f"{r['A_candidates']:>6} | "
            f"{r['B_partial_states']:>9} | "
            f"{r['B_bridge_combinations']:>8} | "
            f"{r['B_candidates']:>6} | "
            f"{str(r['exact'] is not None):>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['total_time']:.6f}s"
        )

    print()
    print(
        "FINISHED EXPERIMENT 176"
    )


if __name__ == "__main__":
    main()
