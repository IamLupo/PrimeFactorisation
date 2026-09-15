#!/usr/bin/env python3
"""
EXPERIMENT 671

SHIFTED SYMMETRIC BIT-BOUND THEOREM

Goal:
    Determine whether the exact bit depth needed to recover v2(C)
    from the symmetric pair (S,D) is controlled by the maximum
    observed v2(C).

Previous result:
    FRAME A: 2C = S - D + 6
    FRAME B: 2C = S + D + 2

This experiment avoids O(N^2) pair separation.

It tests:

    1. Exact baseline canonical depth law.

    2. Exact reconstruction of C from shifted symmetric numerators.

    3. For each k:
           (frame, S mod 2^k, D mod 2^k)
       -> truncated v2(C).

    4. The first k where the signature is exact for the finite domain.

    5. Whether:
           k <= max(v2(C)) + constant
       follows empirically.

    6. Whether the one-bit shift is visible:
           raw SD modulus 2^k
       versus
           shifted numerator modulus 2^k.

    7. A direct algebraic bound:
           knowing S,D mod 2^(k+1)
           determines C mod 2^k.

    8. Prime-range scaling of the minimal exact k.

No O(N^2) comparisons are used.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import gcd
from typing import DefaultDict, Dict, Iterable, List, Optional, Set, Tuple


INF_V2 = 10**9


# =============================================================================
# DATA
# =============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str

    A: int
    B: int
    C: int

    S: int
    D: int

    c: int

    depth: int
    vc: int
    wn: int


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def sieve_primes(limit: int) -> List[int]:
    """Return all odd primes <= limit."""
    if limit < 3:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    root = int(limit ** 0.5)

    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [p for p in range(3, limit + 1, 2) if sieve[p]]


def v2(x: int) -> int:
    """2-adic valuation, with v2(0)=INF_V2."""
    if x == 0:
        return INF_V2

    x = abs(x)
    return (x & -x).bit_length() - 1


def v2_gcd(*values: int) -> int:
    g = 0
    for value in values:
        g = gcd(g, abs(value))
    return v2(g)


# =============================================================================
# FRAME / RESIDUAL MAP
# =============================================================================

def build_state(p: int, q: int) -> Optional[State]:
    """
    Construct the state using the established frame convention.

    FRAME A:
        n == 3 (mod 4)
        A = p - 3
        B = q + 3
        C = B = q + 3
        c = 9

    FRAME B:
        n == 1 (mod 4)
        A = p + 1
        B = q - 3
        C = A = p + 1
        c = 3
    """
    n = p * q

    r = n & 3

    if r == 3:
        frame = "A"

        A = p - 3
        B = q + 3
        C = B

        c = 9

    elif r == 1:
        frame = "B"

        A = p + 1
        B = q - 3
        C = A

        c = 3

    else:
        return None

    S = p + q
    D = p - q

    # Canonical exact law:
    depth = min(v2(C) + 1, v2(n + c))

    vc = v2(C)
    wn = v2(n + c)

    return State(
        n=n,
        p=p,
        q=q,
        frame=frame,
        A=A,
        B=B,
        C=C,
        S=S,
        D=D,
        c=c,
        depth=depth,
        vc=vc,
        wn=wn,
    )


def build_states(prime_limit: int) -> List[State]:
    primes = sieve_primes(prime_limit)

    states: List[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            state = build_state(p, q)

            if state is not None:
                states.append(state)

    return states


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: List[State]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = min(v2(s.C) + 1, v2(s.n + s.c))

        if predicted != s.depth:
            failures += 1

            if failures <= 20:
                print(
                    "mismatch "
                    f"n={s.n} frame={s.frame} "
                    f"C={s.C} vc={s.vc} wn={s.wn} "
                    f"predicted={predicted} actual={s.depth}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_shifted_numerator_identity(states: List[State]) -> int:
    print("=" * 90)
    print("TEST 1: EXACT SHIFTED SYMMETRIC NUMERATOR")
    print("=" * 90)

    failures = 0

    for s in states:
        if s.frame == "A":
            lhs = 2 * s.C
            rhs = s.S - s.D + 6
        else:
            lhs = 2 * s.C
            rhs = s.S + s.D + 2

        if lhs != rhs:
            failures += 1

            if failures <= 20:
                print(
                    "mismatch "
                    f"n={s.n} frame={s.frame} "
                    f"2C={lhs} RHS={rhs}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# SIGNATURE HELPERS
# =============================================================================

def truncated_v2(x: int, k: int) -> int:
    """
    Return min(v2(x), k).
    """
    return min(v2(x), k)


def signature_sd(s: State, k: int) -> Tuple[str, int, int]:
    modulus = 1 << k
    return (
        s.frame,
        s.S % modulus,
        s.D % modulus,
    )


def signature_sd_shifted(s: State, k: int) -> Tuple[str, int, int]:
    """
    Signature carrying one extra bit of symmetric information.

    Since:
        2C = S +/- D + constant

    knowing S,D mod 2^(k+1) determines C mod 2^k.
    """
    modulus = 1 << (k + 1)

    return (
        s.frame,
        s.S % modulus,
        s.D % modulus,
    )


def signature_shifted_numerator(s: State, k: int) -> Tuple[str, int]:
    """
    Direct shifted numerator signature.

    This stores 2C mod 2^(k+1), which determines C mod 2^k.
    """
    modulus = 1 << (k + 1)

    if s.frame == "A":
        numerator = s.S - s.D + 6
    else:
        numerator = s.S + s.D + 2

    return (
        s.frame,
        numerator % modulus,
    )


# =============================================================================
# GENERIC EXACTNESS TEST
# =============================================================================

def ambiguity_count(
    states: List[State],
    signature_fn,
    value_fn,
) -> Tuple[int, int, Set[int]]:
    """
    Return:
        number of signatures
        ambiguous signature count
        observed value set
    """
    buckets: DefaultDict[Tuple, Set[int]] = defaultdict(set)
    all_values: Set[int] = set()

    for s in states:
        signature = signature_fn(s)
        value = value_fn(s)

        buckets[signature].add(value)
        all_values.add(value)

    ambiguous = sum(1 for values in buckets.values() if len(values) > 1)

    return len(buckets), ambiguous, all_values


# =============================================================================
# TEST 2
# =============================================================================

def test_raw_sd_scaling(states: List[State], max_k: int = 16) -> int:
    print("=" * 90)
    print("TEST 2: RAW (S,D) -> TRUNCATED v2(C)")
    print("=" * 90)

    first_exact: Optional[int] = None

    max_vc = max(s.vc for s in states if s.vc < INF_V2)

    print(f"observed max finite v2(C)={max_vc}")
    print()

    for k in range(1, max_k + 1):
        signatures, ambiguous, values = ambiguity_count(
            states,
            lambda s, kk=k: signature_sd(s, kk),
            lambda s, kk=k: truncated_v2(s.C, kk),
        )

        status = "EXACT" if ambiguous == 0 else "AMBIGUOUS"

        print(
            f"k={k:2d} "
            f"signatures={signatures:7d} "
            f"ambiguous={ambiguous:6d} "
            f"values={sorted(values)} "
            f"status={status}"
        )

        if ambiguous == 0 and first_exact is None:
            first_exact = k

    print()

    if first_exact is None:
        print("first exact raw SD k=None")
    else:
        print(f"first exact raw SD k={first_exact}")

    if first_exact is not None:
        print(
            "comparison to max v2(C): "
            f"first_exact={first_exact}, "
            f"max_v2C={max_vc}, "
            f"difference={first_exact - max_vc}"
        )

    print()

    return 0


# =============================================================================
# TEST 3
# =============================================================================

def test_shifted_sd_scaling(states: List[State], max_k: int = 16) -> int:
    print("=" * 90)
    print("TEST 3: SHIFTED (S,D) -> TRUNCATED v2(C)")
    print("=" * 90)

    first_exact: Optional[int] = None

    max_vc = max(s.vc for s in states if s.vc < INF_V2)

    for k in range(1, max_k + 1):
        signatures, ambiguous, values = ambiguity_count(
            states,
            lambda s, kk=k: signature_sd_shifted(s, kk),
            lambda s, kk=k: truncated_v2(s.C, kk),
        )

        status = "EXACT" if ambiguous == 0 else "AMBIGUOUS"

        print(
            f"k={k:2d} "
            f"signatures={signatures:7d} "
            f"ambiguous={ambiguous:6d} "
            f"values={sorted(values)} "
            f"status={status}"
        )

        if ambiguous == 0 and first_exact is None:
            first_exact = k

    print()

    if first_exact is None:
        print("first exact shifted SD k=None")
    else:
        print(f"first exact shifted SD k={first_exact}")
        print(
            "max v2(C)="
            f"{max_vc}, "
            "first_exact-max_v2C="
            f"{first_exact - max_vc}"
        )

    print()

    return 0


# =============================================================================
# TEST 4
# =============================================================================

def test_direct_shifted_numerator(states: List[State], max_k: int = 16) -> int:
    print("=" * 90)
    print("TEST 4: SHIFTED NUMERATOR -> TRUNCATED v2(C)")
    print("=" * 90)

    first_exact: Optional[int] = None

    for k in range(1, max_k + 1):
        signatures, ambiguous, values = ambiguity_count(
            states,
            lambda s, kk=k: signature_shifted_numerator(s, kk),
            lambda s, kk=k: truncated_v2(s.C, kk),
        )

        status = "EXACT" if ambiguous == 0 else "AMBIGUOUS"

        print(
            f"k={k:2d} "
            f"signatures={signatures:7d} "
            f"ambiguous={ambiguous:6d} "
            f"values={sorted(values)} "
            f"status={status}"
        )

        if ambiguous == 0 and first_exact is None:
            first_exact = k

    print()

    if first_exact is None:
        print("first exact shifted numerator k=None")
    else:
        print(
            f"first exact shifted numerator k={first_exact}"
        )

    print()

    return 0


# =============================================================================
# TEST 5
# =============================================================================

def test_bound(states: List[State], max_k: int = 16) -> int:
    print("=" * 90)
    print("TEST 5: BIT-BOUND CHECK")
    print("=" * 90)

    failures = 0

    max_vc = max(s.vc for s in states if s.vc < INF_V2)

    candidate_bounds = [
        ("max_v2C", max_vc),
        ("max_v2C+1", max_vc + 1),
        ("max_v2C+2", max_vc + 2),
    ]

    for name, k in candidate_bounds:
        signatures, ambiguous, _ = ambiguity_count(
            states,
            lambda s, kk=k: signature_sd(s, kk),
            lambda s, kk=k: truncated_v2(s.C, kk),
        )

        exact = ambiguous == 0

        print(
            f"{name:12s} "
            f"k={k:2d} "
            f"signatures={signatures:7d} "
            f"ambiguous={ambiguous:6d} "
            f"exact={exact}"
        )

        if not exact:
            failures += 1

    print()
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_algebraic_bit_loss(states: List[State], max_k: int = 16) -> int:
    print("=" * 90)
    print("TEST 6: ONE-BIT INFORMATION SHIFT")
    print("=" * 90)

    failures = 0

    for k in range(1, max_k + 1):
        modulus = 1 << (k + 1)

        for s in states:
            if s.frame == "A":
                expected = (s.S - s.D + 6) % modulus
            else:
                expected = (s.S + s.D + 2) % modulus

            actual = (2 * s.C) % modulus

            if actual != expected:
                failures += 1

                if failures <= 20:
                    print(
                        "mismatch "
                        f"n={s.n} frame={s.frame} k={k} "
                        f"actual={actual} expected={expected}"
                    )

    checks = len(states) * max_k

    print(f"checks={checks}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_examples(states: List[State]) -> int:
    print("=" * 90)
    print("TEST 7: EXAMPLES")
    print("=" * 90)

    targets = {
        9,
        15,
        21,
        33,
        39,
        57,
        69,
        77,
        87,
        93,
        111,
        141,
        183,
        213,
        485879,
        5579767,
    }

    selected = [s for s in states if s.n in targets]

    selected.sort(key=lambda s: s.n)

    for s in selected:
        print(
            f"n={s.n} p={s.p} q={s.q} frame={s.frame}"
        )
        print(
            f"    S={s.S} D={s.D} C={s.C}"
        )
        print(
            f"    v2(C)={s.vc}"
        )
        print(
            f"    v2(n+c)={s.wn}"
        )
        print(
            f"    depth={s.depth}"
        )

        if s.frame == "A":
            shifted = s.S - s.D + 6
        else:
            shifted = s.S + s.D + 2

        print(
            f"    shifted numerator={shifted}"
        )

        print()

    return 0


# =============================================================================
# PRIME-RANGE SCALING
# =============================================================================

def scaling_experiment(
    limits: Iterable[int],
    max_k: int = 16,
) -> None:
    print("=" * 90)
    print("TEST 8: PRIME-RANGE SCALING")
    print("=" * 90)

    print(
        "limit    states      max_v2C    "
        "raw_exact_k   shifted_exact_k"
    )

    for limit in limits:
        states = build_states(limit)

        finite_vc = [
            s.vc for s in states
            if s.vc < INF_V2
        ]

        max_vc = max(finite_vc) if finite_vc else 0

        raw_exact: Optional[int] = None
        shifted_exact: Optional[int] = None

        for k in range(1, max_k + 1):
            _, raw_ambiguous, _ = ambiguity_count(
                states,
                lambda s, kk=k: signature_sd(s, kk),
                lambda s, kk=k: truncated_v2(s.C, kk),
            )

            if raw_ambiguous == 0:
                raw_exact = k
                break

        for k in range(1, max_k + 1):
            _, shifted_ambiguous, _ = ambiguity_count(
                states,
                lambda s, kk=k: signature_sd_shifted(s, kk),
                lambda s, kk=k: truncated_v2(s.C, kk),
            )

            if shifted_ambiguous == 0:
                shifted_exact = k
                break

        print(
            f"{limit:5d} "
            f"{len(states):10d} "
            f"{max_vc:10d} "
            f"{str(raw_exact):12s} "
            f"{str(shifted_exact):16s}"
        )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    PRIME_LIMIT = 6000
    MAX_K = 16

    print("=" * 90)
    print("EXPERIMENT 671 START")
    print("=" * 90)
    print()
    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve_primes(PRIME_LIMIT)
    print(f"odd primes={len(primes)}")

    states = build_states(PRIME_LIMIT)
    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_shifted_numerator_identity(states)

    test_raw_sd_scaling(
        states,
        max_k=MAX_K,
    )

    test_shifted_sd_scaling(
        states,
        max_k=MAX_K,
    )

    test_direct_shifted_numerator(
        states,
        max_k=MAX_K,
    )

    total_failures += test_bound(states)

    total_failures += test_algebraic_bit_loss(
        states,
        max_k=MAX_K,
    )

    test_examples(states)

    scaling_experiment(
        limits=(100, 300, 1000, 3000, 6000),
        max_k=MAX_K,
    )

    print("=" * 90)
    print("EXPERIMENT 671 FINISHED")
    print("=" * 90)

    if total_failures == 0:
        print("TOTAL FAILURES=0")
        print("STATUS=ALL CORE TESTS PASSED")
    else:
        print(f"TOTAL FAILURES={total_failures}")
        print("STATUS=COUNTEREXAMPLES FOUND")


if __name__ == "__main__":
    main()
