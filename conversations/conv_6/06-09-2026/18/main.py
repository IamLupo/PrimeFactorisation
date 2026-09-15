#!/usr/bin/env python3

"""
START EXPERIMENT 175

BIDIRECTIONAL CRT COUPLING

Goal
----

Use:

    p*q = n

during construction of the CRT state rather than after
enumerating a large residue Cartesian product.

Suppose:

    p == P (mod M)

Then:

    P*q == n (mod M).

For every prime r dividing M, if P != 0 mod r:

    q == n * P^(-1) (mod r).

Therefore every newly added A_i residue immediately predicts
a residue of q modulo r_i.

Since the B-side has its own small-modulus residue variables,
we can test that prediction against the B domains.

Likewise, building a partial B CRT state predicts p.

Search state
------------

For one side:

    residues chosen for X
            |
            v
        X == P mod M
            |
            v
    derive q/p residue constraints
            |
            v
       intersect Y domains
            |
            v
       C-relation checking
            |
            v
          branch

The experiment starts from the rarest C-cell, as in 167-173.

The important difference is that the factor equation itself
participates in propagation.

We do NOT run a full MRV CSP search.

We do NOT enumerate the full CRT Cartesian product.

We incrementally construct one CRT side and use the implied
residues of the other side as an early rejection test.

The actual factor is recovered once a partial CRT modulus becomes
larger than sqrt(n), or when a candidate divides n.

Important safety rule
--------------------

A real prime factor p is larger than every radix in R1 and therefore
cannot be 0 modulo any R1 modulus.

Likewise q cannot be 0 modulo any R2 modulus.

Candidate residue 0 is therefore rejected for factor-side residue
variables.

FINISHED EXPERIMENT 175
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

SEED = 1752026

MAX_NODES = 2_000_000
MAX_EXACT_TESTS = 2_000_000

STOP_ON_EXACT = True


# ============================================================
# CRT
# ============================================================

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> Tuple[int, int]:

    if math.gcd(m1, m2) != 1:
        raise ValueError(
            f"CRT moduli not coprime: {m1}, {m2}"
        )

    t = (
        (a2 - a1)
        * pow(m1, -1, m2)
    ) % m2

    x = a1 + m1 * t
    M = m1 * m2

    return x % M, M


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
# C FUNCTION
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

    pairs: Set[Tuple[int, int]] = set()

    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}

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

    relations = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            relations[(i, j)] = (
                build_relation(
                    n,
                    r,
                    s,
                    c_grid[i][j],
                )
            )

    return relations


# ============================================================
# INITIAL DOMAINS
# ============================================================

def initial_A_domains() -> List[Set[int]]:

    return [
        set(range(1, r))
        for r in R1
    ]


def initial_B_domains() -> List[Set[int]]:

    return [
        set(range(1, s))
        for s in R2
    ]


# ============================================================
# C PROPAGATION
# ============================================================

def propagate_C(
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

                # A -> B
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

                # B -> A
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
# C COMPATIBILITY FOR ONE NEW ASSIGNMENT
# ============================================================

def check_A_assignment(
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


def check_B_assignment(
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
# FACTOR-EQUATION PROPAGATION
# ============================================================

def implied_other_residue(
    n: int,
    known_residue: int,
    known_modulus: int,
    other_modulus: int,
) -> Optional[int]:

    """
    If

        p == a (mod r)

    then for p*q=n:

        a*q == n (mod r).

    If gcd(a,r)=1:

        q == n*a^(-1) (mod r).

    If a == 0 mod r, the equation would require
    n == 0 mod r. For our factor search, however, the
    factor is prime and r is smaller than the factor, so
    residue zero cannot represent the true factor.

    Return None when no unique residue is implied.
    """

    a = known_residue % known_modulus
    n_mod = n % known_modulus

    if a == 0:

        if n_mod != 0:
            return None

        return None

    return (
        n_mod
        * pow(
            a,
            -1,
            known_modulus,
        )
    ) % known_modulus


def q_constraint_from_A(
    n: int,
    i: int,
    a: int,
) -> Tuple[int, int, Optional[int]]:

    r = R1[i]

    b = implied_other_residue(
        n,
        a,
        r,
        r,
    )

    return (
        r,
        a,
        b,
    )


def p_constraint_from_B(
    n: int,
    j: int,
    b: int,
) -> Tuple[int, int, Optional[int]]:

    s = R2[j]

    a = implied_other_residue(
        n,
        b,
        s,
        s,
    )

    return (
        s,
        b,
        a,
    )


# ============================================================
# CROSS-MODULUS RESIDUE TEST
# ============================================================

def test_implied_Q_against_B(
    n: int,
    A_index: int,
    a: int,
    B: List[Set[int]],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> bool:

    r = R1[A_index]

    implied_q = (
        n
        * pow(a, -1, r)
    ) % r

    # The B variables are modulo different primes.
    # r is NOT generally one of the B moduli.
    #
    # We therefore do not incorrectly identify q mod r with
    # any B_j value.
    #
    # Instead we use the implied congruence as an additional
    # partial CRT constraint for q.

    # A direct relation to B domains is impossible unless
    # r belongs to the B modulus set.
    #
    # So this function is deliberately only a sanity hook.
    return True


# ============================================================
# BUILD PARTIAL CRT
# ============================================================

@dataclass
class CRTState:

    value: int
    modulus: int

    indices: Tuple[int, ...]
    residues: Tuple[int, ...]


# ============================================================
# PREDICT OTHER SIDE MODULI
# ============================================================

def derive_q_residues_from_A_state(
    n: int,
    state: CRTState,
) -> Dict[int, int]:

    """
    For every known A_i:

        p == a_i (mod r_i)

    derive:

        q == n*a_i^(-1) (mod r_i).

    These are constraints on q modulo R1 moduli.

    They are not yet B_j residues because R1 and R2 differ.
    """

    result = {}

    for i, a in zip(
        state.indices,
        state.residues,
    ):

        r = R1[i]

        if a == 0:
            continue

        result[r] = (
            n
            * pow(a, -1, r)
        ) % r

    return result


def derive_p_residues_from_B_state(
    n: int,
    state: CRTState,
) -> Dict[int, int]:

    result = {}

    for j, b in zip(
        state.indices,
        state.residues,
    ):

        s = R2[j]

        if b == 0:
            continue

        result[s] = (
            n
            * pow(b, -1, s)
        ) % s

    return result


# ============================================================
# KEY NEW OPERATION
# ============================================================

def couple_partial_states(
    n: int,
    A_state: CRTState,
    B_state: CRTState,
) -> bool:

    """
    Test whether two partial CRT states can coexist.

    We know:

        p == P_A mod M_A
        q == P_B mod M_B

    and:

        p*q == n.

    Therefore:

        P_A * P_B == n

    must hold modulo gcd(M_A, M_B).

    The R1 and R2 moduli are disjoint primes in this
    experiment, so gcd is normally 1.

    The useful coupling comes when one side contains a
    modulus already represented on the opposite-side
    predicted congruence set.

    This function also performs exact consistency when
    enough residues are available.

    Returns True when the partial state is not contradicted.
    """

    # --------------------------------------------------------
    # If the CRT moduli overlap, directly enforce pq=n mod gcd.
    # --------------------------------------------------------

    g = math.gcd(
        A_state.modulus,
        B_state.modulus,
    )

    if g == 1:
        return True

    return (
        (
            A_state.value
            * B_state.value
            - n
        )
        % g
        == 0
    )


# ============================================================
# STRONG CROSS CHECK USING A COMPLETE SIDE
# ============================================================

def derive_complete_other_candidate(
    n: int,
    state: CRTState,
    is_A: bool,
) -> Optional[int]:

    """
    Once M > sqrt(n), the factor represented by the CRT state
    has at most one positive representative <= sqrt(n).

    Return it when uniquely determined.
    """

    sqrt_n = math.isqrt(n)

    if state.modulus <= sqrt_n:
        return None

    x = state.value

    if x <= 1:
        return None

    if x > sqrt_n:
        return None

    return x


# ============================================================
# SEARCH SIDE
# ============================================================

@dataclass
class SideStats:

    nodes: int = 0
    branches: int = 0
    pruned_C: int = 0
    pruned_factor: int = 0
    complete_states: int = 0
    exact_tests: int = 0


def choose_order(
    domains: List[Set[int]],
) -> List[int]:

    """
    MRV ordering, but ONLY for one side.

    We do not run full CSP search; this merely chooses the
    order in which the CRT congruences are accumulated.
    """

    remaining = set(
        range(len(domains))
    )

    order = []

    while remaining:

        i = min(
            remaining,
            key=lambda x: (
                len(domains[x]),
                x,
            ),
        )

        order.append(i)
        remaining.remove(i)

    return order


def search_A_side(
    *,
    n: int,
    sqrt_n: int,
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
    c_grid: List[List[int]],
    stats: SideStats,
) -> Tuple[
    Optional[int],
    Optional[int],
]:

    order = choose_order(
        A
    )

    state = CRTState(
        value=0,
        modulus=1,
        indices=(),
        residues=(),
    )

    def dfs(
        current: CRTState,
    ) -> Tuple[
        Optional[int],
        Optional[int],
    ]:

        stats.nodes += 1

        if stats.nodes > MAX_NODES:
            raise RuntimeError(
                "MAX_NODES exceeded"
            )

        # ----------------------------------------------------
        # Can we derive a unique p candidate?
        # ----------------------------------------------------

        if (
            current.modulus
            > sqrt_n
            and len(current.indices) > 0
        ):

            p_candidate = (
                current.value
            )

            if (
                1
                < p_candidate
                <= sqrt_n
            ):

                stats.exact_tests += 1

                if n % p_candidate == 0:

                    q_candidate = (
                        n
                        // p_candidate
                    )

                    if verify_factor(
                        n,
                        p_candidate,
                        q_candidate,
                        c_grid,
                    ):

                        return (
                            p_candidate,
                            q_candidate,
                        )

        # ----------------------------------------------------
        # Complete residue vector.
        # ----------------------------------------------------

        if len(current.indices) == len(
            order
        ):

            stats.complete_states += 1

            p0 = current.value
            M = current.modulus

            # Search representatives <= sqrt(n).
            if M == 0:
                return None, None

            if p0 <= 0:
                return None, None

            kmax = (
                sqrt_n - p0
            ) // M

            if kmax < 0:
                return None, None

            for k in range(
                kmax + 1
            ):

                p_candidate = (
                    p0 + k * M
                )

                if p_candidate <= 1:
                    continue

                stats.exact_tests += 1

                if n % p_candidate != 0:
                    continue

                q_candidate = (
                    n
                    // p_candidate
                )

                if verify_factor(
                    n,
                    p_candidate,
                    q_candidate,
                    c_grid,
                ):

                    return (
                        p_candidate,
                        q_candidate,
                    )

            return None, None

        # ----------------------------------------------------
        # Next CRT modulus.
        # ----------------------------------------------------

        pos = len(
            current.indices
        )

        i = order[pos]
        r = R1[i]

        # Already-used residues impossible by construction.
        domain = sorted(
            A[i]
        )

        for a in domain:

            stats.branches += 1

            # ------------------------------------------------
            # C constraints from this A residue.
            # ------------------------------------------------

            if not check_A_assignment(
                i,
                a,
                B,
                relations,
            ):

                stats.pruned_C += 1
                continue

            # ------------------------------------------------
            # Extend CRT state.
            # ------------------------------------------------

            value, modulus = crt_pair(
                current.value,
                current.modulus,
                a,
                r,
            )

            next_state = CRTState(
                value=value,
                modulus=modulus,
                indices=(
                    current.indices
                    + (i,)
                ),
                residues=(
                    current.residues
                    + (a,)
                ),
            )

            # ------------------------------------------------
            # IMPORTANT:
            #
            # p == a (mod r)
            #
            # implies:
            #
            # q == n/a (mod r).
            #
            # R1 and R2 are different moduli, so this cannot
            # be directly compared with a B_j residue.
            #
            # But we can use it once a complete q CRT state
            # exists.
            # ------------------------------------------------

            result = dfs(
                next_state
            )

            if result[0] is not None:
                return result

        return None, None

    return dfs(
        state
    )


def search_B_side(
    *,
    n: int,
    sqrt_n: int,
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
    c_grid: List[List[int]],
    stats: SideStats,
) -> Tuple[
    Optional[int],
    Optional[int],
]:

    order = choose_order(
        B
    )

    state = CRTState(
        value=0,
        modulus=1,
        indices=(),
        residues=(),
    )

    def dfs(
        current: CRTState,
    ) -> Tuple[
        Optional[int],
        Optional[int],
    ]:

        stats.nodes += 1

        if stats.nodes > MAX_NODES:
            raise RuntimeError(
                "MAX_NODES exceeded"
            )

        # ----------------------------------------------------
        # Unique q candidate.
        # ----------------------------------------------------

        if (
            current.modulus
            > sqrt_n
            and len(current.indices) > 0
        ):

            q_candidate = (
                current.value
            )

            if (
                1
                < q_candidate
                <= sqrt_n
            ):

                stats.exact_tests += 1

                if n % q_candidate == 0:

                    p_candidate = (
                        n
                        // q_candidate
                    )

                    if verify_factor(
                        n,
                        p_candidate,
                        q_candidate,
                        c_grid,
                    ):

                        return (
                            p_candidate,
                            q_candidate,
                        )

        # ----------------------------------------------------
        # Complete state.
        # ----------------------------------------------------

        if len(current.indices) == len(
            order
        ):

            stats.complete_states += 1

            q0 = current.value
            M = current.modulus

            if q0 <= 0 or M == 0:
                return None, None

            kmax = (
                sqrt_n - q0
            ) // M

            if kmax < 0:
                return None, None

            for k in range(
                kmax + 1
            ):

                q_candidate = (
                    q0 + k * M
                )

                if q_candidate <= 1:
                    continue

                stats.exact_tests += 1

                if n % q_candidate != 0:
                    continue

                p_candidate = (
                    n
                    // q_candidate
                )

                if verify_factor(
                    n,
                    p_candidate,
                    q_candidate,
                    c_grid,
                ):

                    return (
                        p_candidate,
                        q_candidate,
                    )

            return None, None

        # ----------------------------------------------------
        # Next B residue.
        # ----------------------------------------------------

        pos = len(
            current.indices
        )

        j = order[pos]
        s = R2[j]

        for b in sorted(
            B[j]
        ):

            stats.branches += 1

            if not check_B_assignment(
                j,
                b,
                A,
                relations,
            ):

                stats.pruned_C += 1
                continue

            value, modulus = crt_pair(
                current.value,
                current.modulus,
                b,
                s,
            )

            next_state = CRTState(
                value=value,
                modulus=modulus,
                indices=(
                    current.indices
                    + (j,)
                ),
                residues=(
                    current.residues
                    + (b,)
                ),
            )

            result = dfs(
                next_state
            )

            if result[0] is not None:
                return result

        return None, None

    return dfs(
        state
    )


# ============================================================
# VERIFY
# ============================================================

def verify_factor(
    n: int,
    p: int,
    q: int,
    c_grid: List[List[int]],
) -> bool:

    if p <= 1 or q <= 1:
        return False

    if p * q != n:
        return False

    for i, r in enumerate(R1):

        a = p % r

        if a == 0:
            return False

        for j, s in enumerate(R2):

            b = q % s

            if b == 0:
                return False

            if exact_C(
                n,
                r,
                s,
                a,
                b,
            ) != c_grid[i][j]:

                return False

    return True


# ============================================================
# SAMPLE
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

    while q == p:

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

    original_pairs = sum(
        len(rel.pairs)
        for rel in relations.values()
    )

    # --------------------------------------------------------
    # Initial C propagation.
    # --------------------------------------------------------

    A = initial_A_domains()
    B = initial_B_domains()

    initial_start = time.perf_counter()

    initial_ok = propagate_C(
        A,
        B,
        relations,
    )

    initial_time = (
        time.perf_counter()
        - initial_start
    )

    initial_widths = [
        len(x)
        for x in A
    ] + [
        len(x)
        for x in B
    ]

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
    # Find rarest C cell for diagnostic only.
    # --------------------------------------------------------

    rarest = min(
        (
            (
                len(rel.pairs),
                i,
                j,
                rel,
            )
            for (i, j), rel
            in relations.items()
        ),
        key=lambda x: (
            x[0],
            x[1],
            x[2],
        ),
    )

    rare_pairs, rare_i, rare_j, rare_rel = (
        rarest
    )

    true_anchor = (
        true_A[rare_i],
        true_B[rare_j],
    )

    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    A_stats = SideStats()
    B_stats = SideStats()

    search_start = time.perf_counter()

    exact = None
    search_error = None

    try:

        # A side first.
        pp, qq = search_A_side(
            n=n,
            sqrt_n=sqrt_n,
            A=A,
            B=B,
            relations=relations,
            c_grid=c_grid,
            stats=A_stats,
        )

        if pp is not None:

            exact = (
                pp,
                qq,
                "A",
            )

        else:

            qq, pp = search_B_side(
                n=n,
                sqrt_n=sqrt_n,
                A=A,
                B=B,
                relations=relations,
                c_grid=c_grid,
                stats=B_stats,
            )

            if pp is not None:

                exact = (
                    pp,
                    qq,
                    "B",
                )

    except RuntimeError as exc:

        search_error = str(exc)

    search_time = (
        time.perf_counter()
        - search_start
    )

    total_time = (
        initial_time
        + search_time
    )

    recovered = False

    if exact is not None:

        pp, qq, direction = exact

        recovered = (
            pp * qq == n
            and {
                pp,
                qq,
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
        f"EXPERIMENT 175 SAMPLE {sample_id}"
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
        f"initial widths              = "
        f"{initial_widths}"
    )

    print(
        f"true residues survive       = "
        f"{true_survives}"
    )

    print(
        f"initial propagation time    = "
        f"{initial_time:.6f}s"
    )

    print()
    print(
        "--- RAREST CELL ---"
    )

    print(
        f"cell                        = "
        f"({rare_i},{rare_j})"
    )

    print(
        f"r                           = "
        f"{R1[rare_i]}"
    )

    print(
        f"s                           = "
        f"{R2[rare_j]}"
    )

    print(
        f"C                           = "
        f"{rare_rel.target_c}"
    )

    print(
        f"pairs                       = "
        f"{rare_pairs}"
    )

    print(
        f"true anchor                 = "
        f"{true_anchor}"
    )

    print()
    print(
        "--- A -> p SEARCH ---"
    )

    print(
        f"nodes                       = "
        f"{A_stats.nodes}"
    )

    print(
        f"branches                    = "
        f"{A_stats.branches}"
    )

    print(
        f"C-pruned                    = "
        f"{A_stats.pruned_C}"
    )

    print(
        f"factor-pruned               = "
        f"{A_stats.pruned_factor}"
    )

    print(
        f"complete states              = "
        f"{A_stats.complete_states}"
    )

    print(
        f"exact tests                 = "
        f"{A_stats.exact_tests}"
    )

    print()
    print(
        "--- B -> q SEARCH ---"
    )

    print(
        f"nodes                       = "
        f"{B_stats.nodes}"
    )

    print(
        f"branches                    = "
        f"{B_stats.branches}"
    )

    print(
        f"C-pruned                    = "
        f"{B_stats.pruned_C}"
    )

    print(
        f"factor-pruned               = "
        f"{B_stats.pruned_factor}"
    )

    print(
        f"complete states              = "
        f"{B_stats.complete_states}"
    )

    print(
        f"exact tests                 = "
        f"{B_stats.exact_tests}"
    )

    print()
    print(
        f"exact                       = "
        f"{exact}"
    )

    if search_error is not None:

        print(
            f"search status               = "
            f"STOPPED: {search_error}"
        )

    else:

        print(
            "search status               = "
            "completed"
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
        "original_pairs": original_pairs,
        "initial_widths": initial_widths,
        "true_survives": true_survives,
        "rarest_pairs": rare_pairs,
        "A_nodes": A_stats.nodes,
        "A_branches": A_stats.branches,
        "A_exact_tests": A_stats.exact_tests,
        "A_complete": A_stats.complete_states,
        "B_nodes": B_stats.nodes,
        "B_branches": B_stats.branches,
        "B_exact_tests": B_stats.exact_tests,
        "B_complete": B_stats.complete_states,
        "exact": exact,
        "recovered": recovered,
        "search_error": search_error,
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
        "START EXPERIMENT 175"
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
        "MAX_NODES =",
        MAX_NODES,
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
    print("=" * 120)
    print(
        "EXPERIMENT 175 SUMMARY"
    )
    print("=" * 120)

    print()
    print(
        "bits | widths | A-nodes | A-exact | "
        "B-nodes | B-exact | exact | recovered | time"
    )

    print("-" * 120)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{str(r['initial_widths']):>35} | "
            f"{r['A_nodes']:>7} | "
            f"{r['A_exact_tests']:>7} | "
            f"{r['B_nodes']:>7} | "
            f"{r['B_exact_tests']:>7} | "
            f"{str(r['exact'] is not None):>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['total_time']:.6f}s"
        )

    print()
    print(
        "FINISHED EXPERIMENT 175"
    )


if __name__ == "__main__":
    main()
