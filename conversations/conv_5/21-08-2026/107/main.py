#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# =============================================================================
# EXPERIMENT 669
#
# 2-ADIC SYMMETRIC-COORDINATE RECOVERY OF v2(C)
#
# Goal:
#
# Experiment 668 showed that:
#
#     v2(S)
#     v2(D)
#     (v2(S),v2(D))
#
# do NOT determine v2(C).
#
# Here:
#
#     S = p + q
#     D = p - q
#
# and C is the canonical residual:
#
#     FRAME A: C = q + 3
#     FRAME B: C = p + 1
#
# We now retain the actual residue information:
#
#     S mod 2^k
#     D mod 2^k
#
# and ask:
#
#     Does (frame, S mod 2^k, D mod 2^k)
#     determine min(v2(C), k)?
#
# Since:
#
#     p = (S + D)/2
#     q = (S - D)/2,
#
# the pair (S,D) contains the factor coordinates.
#
# The interesting question is whether LOW 2-ADIC BITS
# of the symmetric coordinates already contain the
# canonical branch valuation.
#
# Tests:
#
#   0. baseline
#   1. exact reconstruction of p,q from S,D
#   2. residue-signature ambiguity for k=2..16
#   3. first k that determines truncated v2(C)
#   4. scaling of minimal k by prime limit
#   5. compare S-only, D-only, and (S,D)
#   6. threshold formulation
#   7. explicit collision examples
#
# =============================================================================


PRIME_LIMITS = [100, 300, 1000, 3000, 6000]

MAIN_LIMIT = 6000
MAX_K = 16
INF = 10**9


# =============================================================================
# ARITHMETIC
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
        p for p in range(3, limit + 1)
        if a[p] and (p & 1)
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
        # C = B
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
        # C = A
        # c = 3
        frame = "B"
        A = p + 1
        B = q - 3
        C = A
        c = 3

    vc = v2(C)
    w = v2(n + c)
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
            states.append(make_state(p, q))

    return states


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        if s.depth != min(s.vc + 1, s.w):
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
#
# S,D -> p,q exactly.
# =============================================================================

def test_symmetric_reconstruction(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: EXACT SYMMETRIC RECONSTRUCTION")
    print("=" * 90)

    failures = 0

    for s in states:
        if (s.S + s.D) % 2 != 0:
            failures += 1
            continue

        if (s.S - s.D) % 2 != 0:
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
# TRUNCATED VALUATION
# =============================================================================

def truncated_v2(x: int, k: int) -> int:
    return min(v2(x), k)


# =============================================================================
# GENERIC SIGNATURE ANALYSIS
# =============================================================================

def analyze_signature(
    states: list[State],
    key_fn,
    value_fn,
) -> tuple[int, int, dict]:
    buckets: dict[tuple, set[int]] = defaultdict(set)

    for s in states:
        key = key_fn(s)
        value = value_fn(s)
        buckets[key].add(value)

    ambiguous = {
        key: values
        for key, values in buckets.items()
        if len(values) > 1
    }

    return (
        len(buckets),
        len(ambiguous),
        ambiguous,
    )


# =============================================================================
# TEST 2
#
# FULL (S,D) RESIDUE SIGNATURE
#
# Check whether:
#
#     (frame, S mod 2^k, D mod 2^k)
#
# determines:
#
#     min(v2(C),k)
# =============================================================================

def test_sd_residue_levels(
    states: list[State],
    max_k: int,
) -> dict[int, int]:
    print("=" * 90)
    print("TEST 2: (S,D) RESIDUE -> TRUNCATED v2(C)")
    print("=" * 90)

    minimal_exact = None
    ambiguity_by_k: dict[int, int] = {}

    for k in range(1, max_k + 1):
        modulus = 1 << k

        buckets: dict[tuple, set[int]] = defaultdict(set)

        for s in states:
            signature = (
                s.frame,
                s.S % modulus,
                s.D % modulus,
            )

            buckets[signature].add(
                truncated_v2(s.C, k)
            )

        ambiguous = sum(
            1
            for values in buckets.values()
            if len(values) > 1
        )

        ambiguity_by_k[k] = ambiguous

        print(
            f"    k={k:2d}"
            f" signatures={len(buckets):6d}"
            f" ambiguous={ambiguous:6d}"
        )

        if ambiguous == 0 and minimal_exact is None:
            minimal_exact = k

    print()

    if minimal_exact is None:
        print("minimal exact k=None")
    else:
        print(f"minimal exact k={minimal_exact}")

    print()

    return ambiguity_by_k


# =============================================================================
# TEST 3
#
# S-only and D-only comparison.
# =============================================================================

def test_single_coordinate_levels(
    states: list[State],
    max_k: int,
) -> None:
    print("=" * 90)
    print("TEST 3: S-ONLY VS D-ONLY VS (S,D)")
    print("=" * 90)

    for mode in ("S", "D"):
        print()
        print(f"MODE={mode}")

        for k in range(1, max_k + 1):
            modulus = 1 << k

            buckets: dict[tuple, set[int]] = defaultdict(set)

            for s in states:
                if mode == "S":
                    residue = s.S % modulus
                else:
                    residue = s.D % modulus

                buckets[
                    (s.frame, residue)
                ].add(
                    truncated_v2(s.C, k)
                )

            ambiguous = sum(
                1
                for values in buckets.values()
                if len(values) > 1
            )

            print(
                f"    k={k:2d}"
                f" signatures={len(buckets):6d}"
                f" ambiguous={ambiguous:6d}"
            )

    print()


# =============================================================================
# TEST 4
#
# Does a FIXED k determine FULL v2(C) over the tested domain?
#
# Important distinction:
#
#   truncated v2(C) may be determined
#   while full v2(C) exceeds k.
#
# We therefore specifically search for collisions where:
#
#   same signature
#   different actual vc
# =============================================================================

def test_full_vc_collisions(
    states: list[State],
    max_k: int,
) -> None:
    print("=" * 90)
    print("TEST 4: FULL v2(C) COLLISIONS")
    print("=" * 90)

    for k in range(1, max_k + 1):
        modulus = 1 << k

        buckets: dict[tuple, dict[int, list[State]]] = defaultdict(
            lambda: defaultdict(list)
        )

        for s in states:
            signature = (
                s.frame,
                s.S % modulus,
                s.D % modulus,
            )

            buckets[signature][s.vc].append(s)

        collisions = []

        for signature, by_vc in buckets.items():
            if len(by_vc) > 1:
                collisions.append(
                    (signature, by_vc)
                )

        print(
            f"    k={k:2d}"
            f" full-vc ambiguous={len(collisions):6d}"
        )

        if collisions and k <= 8:
            for signature, by_vc in collisions[:3]:
                values = sorted(by_vc)
                print(
                    f"        {signature} -> {values}"
                )

    print()


# =============================================================================
# TEST 5
#
# First k needed for EACH state to distinguish its
# truncated C valuation from all states sharing the
# same residue signature.
# =============================================================================

def test_state_separation(
    states: list[State],
    max_k: int,
) -> None:
    print("=" * 90)
    print("TEST 5: STATE SEPARATION LEVEL")
    print("=" * 90)

    first_exact: list[int | None] = []

    for idx, s in enumerate(states):
        found = None

        for k in range(1, max_k + 1):
            modulus = 1 << k

            target = truncated_v2(s.C, k)

            valid = True

            for t in states:
                if t.frame != s.frame:
                    continue

                if (
                    t.S % modulus == s.S % modulus
                    and
                    t.D % modulus == s.D % modulus
                    and
                    truncated_v2(t.C, k) != target
                ):
                    valid = False
                    break

            if valid:
                found = k
                break

        first_exact.append(found)

    counter = defaultdict(int)

    for value in first_exact:
        counter[value] += 1

    print("first-exact-k distribution:")

    for key in sorted(counter, key=lambda x: (x is None, x)):
        print(
            f"    k={str(key):>4}"
            f" count={counter[key]}"
        )

    print()


# =============================================================================
# TEST 6
#
# Threshold version:
#
#     v2(C) >= t
#
# iff:
#
#     C == 0 mod 2^t
#
# and therefore test whether this condition is visible
# from the low bits of S,D.
# =============================================================================

def test_threshold_visibility(
    states: list[State],
    max_t: int,
) -> int:
    print("=" * 90)
    print("TEST 6: THRESHOLD VISIBILITY OF C")
    print("=" * 90)

    failures = 0

    for t in range(1, max_t + 1):
        modulus = 1 << t

        buckets: dict[tuple, set[bool]] = defaultdict(set)

        for s in states:
            signature = (
                s.frame,
                s.S % modulus,
                s.D % modulus,
            )

            truth = (s.C % modulus == 0)

            buckets[signature].add(truth)

        ambiguous = sum(
            1
            for values in buckets.values()
            if len(values) > 1
        )

        print(
            f"    t={t:2d}"
            f" threshold-ambiguous={ambiguous:6d}"
        )

    print()

    return failures


# =============================================================================
# TEST 7
#
# Direct algebraic relation:
#
#     FRAME A:
#         C = q+3 = (S-D)/2 + 3
#
#         2C = S-D+6
#
#     FRAME B:
#         C = p+1 = (S+D)/2 + 1
#
#         2C = S+D+2
#
# This means C valuation can be tested directly against
# residue arithmetic on S,D.
# =============================================================================

def test_algebraic_c_forms(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 7: ALGEBRAIC C REPRESENTATION")
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

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(states: list[State]) -> None:
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

    by_n = {s.n: s for s in states}

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
            f"    v2(C)={s.vc}"
        )

        print(
            f"    v2(S)={v2(s.S)}"
        )

        print(
            f"    v2(D)={v2(s.D)}"
        )

        print(
            f"    v2(n+c)={s.w}"
        )

        print(
            f"    depth={s.depth}"
        )

        for k in (3, 4, 6, 8, 10):
            modulus = 1 << k

            print(
                f"    k={k:2d}:"
                f" S mod 2^k={s.S % modulus:>5d}"
                f" D mod 2^k={s.D % modulus:>5d}"
                f" C mod 2^k={s.C % modulus:>5d}"
            )


# =============================================================================
# TEST 8
#
# Compare direct C reconstruction from S,D modulo 2^k.
#
# Since 2 is not invertible modulo 2^k, work modulo 2^(k+1):
#
#     FRAME A:
#         2C = S-D+6
#
#     FRAME B:
#         2C = S+D+2
#
# Then determine whether the exact truncated v2(C)
# follows from the right side.
# =============================================================================

def test_shifted_sum_difference(
    states: list[State],
    max_k: int,
) -> None:
    print("=" * 90)
    print("TEST 8: SHIFTED SYMMETRIC NUMERATOR")
    print("=" * 90)

    for k in range(1, max_k + 1):
        modulus = 1 << (k + 1)

        buckets: dict[tuple, set[int]] = defaultdict(set)

        for s in states:
            if s.frame == "A":
                z = (s.S - s.D + 6) % modulus
            else:
                z = (s.S + s.D + 2) % modulus

            buckets[
                (s.frame, z)
            ].add(
                truncated_v2(s.C, k)
            )

        ambiguous = sum(
            1
            for values in buckets.values()
            if len(values) > 1
        )

        print(
            f"    k={k:2d}"
            f" signatures={len(buckets):6d}"
            f" ambiguous={ambiguous:6d}"
        )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 669 START")
    print("=" * 90)
    print()

    print(f"main prime limit={MAIN_LIMIT}")

    primes = sieve(MAIN_LIMIT)

    print(f"odd primes={len(primes)}")

    states = build_states(MAIN_LIMIT)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_symmetric_reconstruction(states)

    test_sd_residue_levels(
        states,
        MAX_K,
    )

    test_single_coordinate_levels(
        states,
        MAX_K,
    )

    test_full_vc_collisions(
        states,
        MAX_K,
    )

    test_state_separation(
        states,
        MAX_K,
    )

    total_failures += test_threshold_visibility(
        states,
        min(MAX_K, 12),
    )

    total_failures += test_algebraic_c_forms(states)

    test_shifted_sum_difference(
        states,
        MAX_K,
    )

    print_examples(states)

    # =====================================================================
    # SCALING TEST
    # =====================================================================

    print("=" * 90)
    print("TEST 9: PRIME-RANGE SCALING OF (S,D) RESIDUE RECOVERY")
    print("=" * 90)

    for limit in PRIME_LIMITS:
        local_states = build_states(limit)

        exact_k = None

        for k in range(1, MAX_K + 1):
            modulus = 1 << k

            buckets: dict[tuple, set[int]] = defaultdict(set)

            for s in local_states:
                buckets[
                    (
                        s.frame,
                        s.S % modulus,
                        s.D % modulus,
                    )
                ].add(
                    truncated_v2(s.C, k)
                )

            ambiguous = any(
                len(values) > 1
                for values in buckets.values()
            )

            if not ambiguous:
                exact_k = k
                break

        print(
            f"    limit={limit:5d}"
            f" states={len(local_states):8d}"
            f" minimal_exact_k={str(exact_k):>4}"
        )

    print()

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print(
        "Experiment 668 established that"
    )
    print(
        "v2(p+q), v2(p-q), and their pair do not"
    )
    print(
        "determine v2(C)."
    )
    print()

    print(
        "Experiment 669 asks whether the missing"
    )
    print(
        "information is actually present in the LOW"
    )
    print(
        "2-ADIC BITS of the symmetric coordinates:"
    )
    print()

    print(
        "    (S mod 2^k, D mod 2^k)"
    )
    print()

    print(
        "The exact algebraic identities are:"
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
        "Therefore this experiment distinguishes two"
    )
    print(
        "possibilities:"
    )
    print()

    print(
        "1. A finite number of low symmetric bits"
    )
    print(
        "   determines the truncated C valuation."
    )
    print()
    print(
        "2. Even the full low-bit symmetric projection"
    )
    print(
        "   exhibits collisions, meaning the canonical"
    )
    print(
        "   branch valuation remains genuinely branch-local."
    )
    print()

    print(
        f"TOTAL TEST FAILURES={total_failures}"
    )

    if total_failures == 0:
        print("STATUS=CORE ALGEBRAIC TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 669 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
