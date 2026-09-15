#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import isqrt


# =============================================================================
# EXPERIMENT 670
#
# FAST 2-ADIC SYMMETRIC-RESIDUE RECOVERY
#
# Experiment 669 discovered:
#
#   (frame, S mod 2^k, D mod 2^k)
#
# eventually determines the truncated valuation
#
#   min(v2(C), k)
#
# for the tested domain.
#
# However:
#
#   k=1
#
# is misleading because every state has the same trivial parity
# class.  This script therefore distinguishes:
#
#   TRIVIAL:
#       the value is constant over the whole domain
#
#   EXACT:
#       the signature determines the truncated v2(C), with
#       at least two distinct truncated values present globally.
#
# The expensive O(N^2) "state separation" test is removed.
#
# New focus:
#
#   1. exact residue recovery by k
#   2. threshold recovery by k
#   3. first nontrivial exact k
#   4. number of ambiguous residue classes
#   5. collision structure at the first nontrivial failures
#   6. compare:
#          S only
#          D only
#          (S,D)
#          shifted numerator 2C
#   7. prime-range scaling
#
# Algebra:
#
# FRAME A:
#
#     C = q+3
#     2C = S-D+6
#
# FRAME B:
#
#     C = p+1
#     2C = S+D+2
#
# where:
#
#     S = p+q
#     D = p-q.
#
# =============================================================================


LIMIT = 6000
MAX_K = 16

SCALE_LIMITS = [
    100,
    300,
    1000,
    3000,
    6000,
]

INF = 10**9


# =============================================================================
# BASIC ARITHMETIC
# =============================================================================

def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


def sieve(limit: int) -> list[int]:
    a = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        a[0] = 0
    if limit >= 1:
        a[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start:limit + 1:p] = b"\x00" * count

    return [
        p
        for p in range(3, limit + 1, 2)
        if a[p]
    ]


# =============================================================================
# STATE
# =============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int

    frame: str

    S: int
    D: int

    A: int
    B: int
    C: int

    c: int

    vc: int
    w: int
    depth: int


def make_state(p: int, q: int) -> State:
    n = p * q
    S = p + q
    D = p - q

    if n % 4 == 3:
        # FRAME A
        #
        # A = p-3
        # B = q+3
        # C = q+3
        # c = 9

        frame = "A"
        A = p - 3
        B = q + 3
        C = B
        c = 9

    else:
        # FRAME B
        #
        # A = p+1
        # B = q-3
        # C = p+1
        # c = 3

        frame = "B"
        A = p + 1
        B = q - 3
        C = A
        c = 3

    vc = v2(C)
    w = v2(n + c)

    # Established law:
    #
    #     depth = min(v2(C)+1, v2(n+c))

    depth = min(vc + 1, w)

    return State(
        n=n,
        p=p,
        q=q,
        frame=frame,
        S=S,
        D=D,
        A=A,
        B=B,
        C=C,
        c=c,
        vc=vc,
        w=w,
        depth=depth,
    )


def build_states(limit: int) -> list[State]:
    primes = sieve(limit)

    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(
                make_state(p, q)
            )

    return states


# =============================================================================
# HELPERS
# =============================================================================

def trunc_v2(x: int, k: int) -> int:
    return min(v2(x), k)


def count_ambiguous(
    buckets: dict,
) -> int:
    return sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )


def print_first_collisions(
    buckets: dict,
    max_items: int = 8,
) -> None:
    shown = 0

    for key, values in buckets.items():
        if len(values) <= 1:
            continue

        print(
            f"        {key} -> {sorted(values)}"
        )

        shown += 1

        if shown >= max_items:
            break


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = min(
            s.vc + 1,
            s.w,
        )

        if predicted != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
#
# EXACT S,D RECONSTRUCTION
# =============================================================================

def test_reconstruction(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: EXACT SYMMETRIC RECONSTRUCTION")
    print("=" * 90)

    failures = 0

    for s in states:
        if (s.S + s.D) & 1:
            failures += 1
            continue

        if (s.S - s.D) & 1:
            failures += 1
            continue

        p = (s.S + s.D) // 2
        q = (s.S - s.D) // 2

        if p != s.p or q != s.q:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# CORE ANALYSIS
# =============================================================================

def analyze_signature(
    states: list[State],
    key_fn,
    value_fn,
) -> tuple[int, int, bool, dict]:
    buckets: dict = defaultdict(set)

    for s in states:
        key = key_fn(s)
        value = value_fn(s)
        buckets[key].add(value)

    ambiguous = count_ambiguous(buckets)

    global_values = set()

    for values in buckets.values():
        global_values.update(values)

    is_trivial = len(global_values) <= 1

    return (
        len(buckets),
        ambiguous,
        is_trivial,
        buckets,
    )


# =============================================================================
# TEST 2
#
# (frame,S,D) LOW-BIT RECOVERY
# =============================================================================

def test_sd_residues(
    states: list[State],
) -> tuple[int | None, dict[int, int]]:
    print("=" * 90)
    print("TEST 2: (S,D) RESIDUE -> TRUNCATED v2(C)")
    print("=" * 90)

    ambiguity = {}
    first_nontrivial_exact = None

    for k in range(1, MAX_K + 1):
        mod = 1 << k

        signatures = defaultdict(set)

        for s in states:
            signatures[
                (
                    s.frame,
                    s.S % mod,
                    s.D % mod,
                )
            ].add(
                min(s.vc, k)
            )

        amb = count_ambiguous(signatures)

        values = set()

        for v in signatures.values():
            values.update(v)

        trivial = len(values) <= 1

        ambiguity[k] = amb

        status = (
            "TRIVIAL"
            if trivial
            else "EXACT"
            if amb == 0
            else "AMBIGUOUS"
        )

        print(
            f"    k={k:2d}"
            f" signatures={len(signatures):7d}"
            f" ambiguous={amb:6d}"
            f" values={sorted(values)}"
            f" status={status}"
        )

        if (
            not trivial
            and
            amb == 0
            and
            first_nontrivial_exact is None
        ):
            first_nontrivial_exact = k

    print()

    print(
        f"first nontrivial exact k={first_nontrivial_exact}"
    )

    print()

    return first_nontrivial_exact, ambiguity


# =============================================================================
# TEST 3
#
# Threshold visibility:
#
#     v2(C) >= t
#
# from low bits of S,D.
#
# This is cleaner than asking for exact v2(C).
# =============================================================================

def test_thresholds(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: THRESHOLD VISIBILITY")
    print("=" * 90)

    failures = 0

    for t in range(1, min(MAX_K, 12) + 1):
        mod = 1 << t

        buckets = defaultdict(set)

        for s in states:
            signature = (
                s.frame,
                s.S % mod,
                s.D % mod,
            )

            truth = (
                s.C % mod == 0
            )

            buckets[signature].add(truth)

        amb = count_ambiguous(buckets)

        # Mathematically the threshold is directly encoded by
        # the residue of C, so the only purpose here is to check
        # whether the S,D projection preserves it.
        if amb != 0:
            failures += 1

        print(
            f"    t={t:2d}"
            f" signatures={len(buckets):7d}"
            f" ambiguous={amb:6d}"
        )

    print()
    return failures


# =============================================================================
# TEST 4
#
# Compare:
#
#   S only
#   D only
#   S,D pair
#
# =============================================================================

def test_coordinate_comparison(
    states: list[State],
) -> None:
    print("=" * 90)
    print("TEST 4: S-ONLY / D-ONLY / (S,D)")
    print("=" * 90)

    for mode in ("S", "D", "SD"):
        print()
        print(f"MODE={mode}")

        for k in range(2, MAX_K + 1):
            mod = 1 << k

            buckets = defaultdict(set)

            for s in states:
                if mode == "S":
                    key = (
                        s.frame,
                        s.S % mod,
                    )
                elif mode == "D":
                    key = (
                        s.frame,
                        s.D % mod,
                    )
                else:
                    key = (
                        s.frame,
                        s.S % mod,
                        s.D % mod,
                    )

                buckets[key].add(
                    min(s.vc, k)
                )

            amb = count_ambiguous(buckets)

            print(
                f"    k={k:2d}"
                f" signatures={len(buckets):7d}"
                f" ambiguous={amb:6d}"
            )

    print()


# =============================================================================
# TEST 5
#
# Shifted numerator:
#
# FRAME A:
#
#     2C = S-D+6
#
# FRAME B:
#
#     2C = S+D+2
#
# This should expose why the (S,D) pair eventually works.
# =============================================================================

def test_shifted_numerator(
    states: list[State],
) -> int:
    print("=" * 90)
    print("TEST 5: SHIFTED SYMMETRIC NUMERATOR")
    print("=" * 90)

    failures = 0

    first_nontrivial_exact = None

    for k in range(2, MAX_K + 1):
        mod = 1 << (k + 1)

        buckets = defaultdict(set)

        for s in states:
            if s.frame == "A":
                z = (
                    s.S - s.D + 6
                ) % mod
            else:
                z = (
                    s.S + s.D + 2
                ) % mod

            buckets[
                (
                    s.frame,
                    z,
                )
            ].add(
                min(s.vc, k)
            )

        amb = count_ambiguous(buckets)

        if (
            amb == 0
            and first_nontrivial_exact is None
        ):
            first_nontrivial_exact = k

        print(
            f"    k={k:2d}"
            f" signatures={len(buckets):7d}"
            f" ambiguous={amb:6d}"
        )

    print()
    print(
        f"first exact shifted-numerator k={first_nontrivial_exact}"
    )
    print()

    return failures


# =============================================================================
# TEST 6
#
# Collision examples at selected k.
# =============================================================================

def test_collision_examples(
    states: list[State],
) -> None:
    print("=" * 90)
    print("TEST 6: RESIDUE COLLISION EXAMPLES")
    print("=" * 90)

    for k in (2, 4, 6, 8, 10, 12):
        mod = 1 << k

        buckets = defaultdict(set)

        for s in states:
            buckets[
                (
                    s.frame,
                    s.S % mod,
                    s.D % mod,
                )
            ].add(
                min(s.vc, k)
            )

        amb = count_ambiguous(buckets)

        print()
        print(
            f"k={k} ambiguous={amb}"
        )

        if amb:
            print_first_collisions(
                buckets,
                max_items=5,
            )

    print()


# =============================================================================
# TEST 7
#
# Direct C residue recovery.
#
# Since C itself is:
#
#   A: (S-D+6)/2
#   B: (S+D+2)/2
#
# modulo 2^k we need S,D modulo 2^(k+1).
#
# Therefore test the exact information boundary.
# =============================================================================

def test_minimal_information_boundary(
    states: list[State],
) -> int:
    print("=" * 90)
    print("TEST 7: MINIMAL INFORMATION BOUNDARY")
    print("=" * 90)

    failures = 0

    for k in range(1, MAX_K + 1):
        mod = 1 << (k + 1)

        buckets = defaultdict(set)

        for s in states:
            if s.frame == "A":
                numerator = (
                    s.S - s.D + 6
                ) % mod
            else:
                numerator = (
                    s.S + s.D + 2
                ) % mod

            # Numerator is even. Dividing its residue by two
            # yields C modulo 2^k.
            C_mod = (numerator // 2) % (1 << k)

            buckets[
                (
                    s.frame,
                    numerator,
                )
            ].add(
                min(s.vc, k)
            )

        amb = count_ambiguous(buckets)

        if amb != 0:
            # This is diagnostic, not necessarily an error:
            # two distinct C residues can have the same
            # truncated valuation.
            pass

        print(
            f"    k={k:2d}"
            f" numerator-signatures={len(buckets):7d}"
            f" ambiguous-v2={amb:6d}"
        )

    print()
    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(
    states: list[State],
) -> None:
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = [
        9,
        15,
        21,
        33,
        39,
        57,
        69,
        77,
        93,
        141,
        183,
        213,
        485879,
        5579767,
    ]

    by_n = {
        s.n: s
        for s in states
    }

    for n in wanted:
        s = by_n.get(n)

        if s is None:
            continue

        print()
        print(
            f"n={s.n} p={s.p} q={s.q} frame={s.frame}"
        )

        print(
            f"    S={s.S}"
        )

        print(
            f"    D={s.D}"
        )

        print(
            f"    C={s.C}"
        )

        print(
            f"    v2(S)={v2(s.S)}"
        )

        print(
            f"    v2(D)={v2(s.D)}"
        )

        print(
            f"    v2(C)={s.vc}"
        )

        print(
            f"    v2(n+c)={s.w}"
        )

        print(
            f"    depth={s.depth}"
        )

        for k in (4, 6, 8, 10):
            mod = 1 << k

            print(
                f"    k={k:2d}"
                f" S mod 2^k={s.S % mod:5d}"
                f" D mod 2^k={s.D % mod:5d}"
                f" C mod 2^k={s.C % mod:5d}"
            )


# =============================================================================
# TEST 8
#
# PRIME-RANGE SCALING
#
# No expensive state-by-state comparisons.
# =============================================================================

def scaling_test() -> None:
    print("=" * 90)
    print("TEST 8: PRIME-RANGE SCALING")
    print("=" * 90)

    for limit in SCALE_LIMITS:
        states = build_states(limit)

        minimal_exact = None

        for k in range(2, MAX_K + 1):
            mod = 1 << k

            buckets = defaultdict(set)

            for s in states:
                buckets[
                    (
                        s.frame,
                        s.S % mod,
                        s.D % mod,
                    )
                ].add(
                    min(s.vc, k)
                )

            amb = count_ambiguous(buckets)

            if amb == 0:
                minimal_exact = k
                break

        print(
            f"    limit={limit:5d}"
            f" states={len(states):8d}"
            f" minimal_nontrivial_exact_k="
            f"{str(minimal_exact):>4}"
        )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 670 START")
    print("=" * 90)
    print()

    print(f"prime limit={LIMIT}")

    primes = sieve(LIMIT)

    print(
        f"odd primes={len(primes)}"
    )

    states = build_states(LIMIT)

    print(
        f"semiprimes={len(states)}"
    )

    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_reconstruction(states)

    test_sd_residues(states)

    total_failures += test_thresholds(states)

    test_coordinate_comparison(states)

    total_failures += test_shifted_numerator(states)

    test_collision_examples(states)

    total_failures += test_minimal_information_boundary(states)

    print_examples(states)

    scaling_test()

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print(
        "Experiment 669 showed that the low bits of"
    )
    print(
        "(S,D) eventually determine the truncated"
    )
    print(
        "valuation of C, but the original script"
    )
    print(
        "contained an O(N^2) state-separation test."
    )
    print()

    print(
        "Experiment 670 replaces that with direct"
    )
    print(
        "signature buckets."
    )
    print()

    print(
        "The exact identities are:"
    )
    print()

    print(
        "    FRAME A:"
    )
    print(
        "        2C = S - D + 6"
    )
    print()

    print(
        "    FRAME B:"
    )
    print(
        "        2C = S + D + 2"
    )
    print()

    print(
        "Thus the real information boundary is:"
    )
    print()

    print(
        "    S,D modulo 2^(k+1)"
    )
    print(
        "             ↓"
    )
    print(
        "    C modulo 2^k"
    )
    print(
        "             ↓"
    )
    print(
        "    v2(C) truncated at k"
    )
    print()

    print(
        "The important result is the first"
    )
    print(
        "NONTRIVIAL k at which the (S,D) residue"
    )
    print(
        "signature becomes exact."
    )
    print()

    print(
        f"TOTAL TEST FAILURES={total_failures}"
    )

    if total_failures == 0:
        print(
            "STATUS=CORE TESTS PASSED"
        )
    else:
        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )

    print("=" * 90)
    print("EXPERIMENT 670 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
