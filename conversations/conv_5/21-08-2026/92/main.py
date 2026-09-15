#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isqrt


# =============================================================================
# EXPERIMENT 654
# ONE-BIT TERMINAL BRANCH THEOREM
# =============================================================================

INF = 10**9
PRIME_LIMIT = 6000


# =============================================================================
# NUMBER THEORY
# =============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def sieve_odd_primes(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [p for p in range(3, limit + 1, 2) if sieve[p]]


# =============================================================================
# STATE
# =============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str

    A: int
    B: int

    p0: int
    q0: int
    c: int


def build_states(primes: list[int]) -> list[State]:
    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            if n % 4 == 3:
                # FRAME A
                frame = "A"
                A = p - 3
                B = q + 3
                p0 = 3
                q0 = -3
                c = 9

            else:
                # FRAME B
                frame = "B"
                A = p + 1
                B = q - 3
                p0 = -1
                q0 = 3
                c = 3

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    A=A,
                    B=B,
                    p0=p0,
                    q0=q0,
                    c=c,
                )
            )

    return states


# =============================================================================
# CORE QUANTITIES
# =============================================================================

def transformed_numerators(s: State) -> tuple[int, int]:
    if s.frame == "A":
        return (
            s.B - s.A,
            s.A + s.B,
        )

    return (
        3 * s.A - s.B,
        3 * s.A + s.B,
    )


def direct_depth(s: State) -> int:
    u, v = transformed_numerators(s)
    return min(v2(u), v2(v))


def tp(s: State) -> int:
    return v2(s.p - s.p0)


def w(s: State) -> int:
    return v2(s.n + s.c)


def law_depth(s: State) -> int:
    return min(tp(s) + 1, w(s))


# =============================================================================
# PREFIX TESTS
# =============================================================================

def prefix_equal_p(s: State, bits: int) -> bool:
    """
    p == p0 (mod 2^bits)
    """
    modulus = 1 << bits
    return (s.p - s.p0) % modulus == 0


def prefix_equal_n(s: State, bits: int) -> bool:
    """
    n == -c (mod 2^bits)
    """
    modulus = 1 << bits
    return (s.n + s.c) % modulus == 0


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)

    failures = 0

    for s in states:
        d1 = direct_depth(s)
        d2 = law_depth(s)

        if d1 != d2:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s.n} frame={s.frame} "
                    f"p={s.p} q={s.q} "
                    f"direct={d1} law={d2}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_branch_first_mismatch(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 1: BRANCH-FIRST-MISMATCH THEOREM")
    print("=" * 90)

    failures = 0

    for s in states:
        d = direct_depth(s)
        W = w(s)
        T = tp(s)

        # Only relevant when branch, rather than n, determines depth.
        if T < W:
            predicted = T + 1

            if predicted != d:
                failures += 1

            # At depth d:
            #   p must still match through d-1 bits.
            # At depth d+1:
            #   p must fail.
            if not prefix_equal_p(s, d - 1):
                failures += 1

            if prefix_equal_p(s, d):
                failures += 1

            # n must survive the actual depth.
            if not prefix_equal_n(s, d):
                failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_n_bound_branch(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 2: N-BOUND TERMINAL THEOREM")
    print("=" * 90)

    failures = 0

    for s in states:
        d = direct_depth(s)
        W = w(s)
        T = tp(s)

        if W <= T + 1:
            if d != W:
                failures += 1

            # To reach the n-bound, p must survive through W-1.
            if W >= 2:
                if not prefix_equal_p(s, W - 1):
                    failures += 1

            # n survives exactly through W and fails at W+1.
            if not prefix_equal_n(s, W):
                failures += 1

            if prefix_equal_n(s, W + 1):
                failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_single_bit_partition(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 3: SINGLE-BIT TERMINAL PARTITION")
    print("=" * 90)

    """
    For each possible depth d:

        d < w:
            p matches anchor through d-1
            and first fails at d.

        d == w:
            p matches anchor through w-1.
    """

    failures = 0

    for s in states:
        d = direct_depth(s)
        W = w(s)

        if d < W:
            # Branch stopped first.
            lhs = (
                prefix_equal_p(s, d - 1)
                and not prefix_equal_p(s, d)
            )

            if not lhs:
                failures += 1

        else:
            # N stopped first or simultaneously.
            if not prefix_equal_p(s, d - 1):
                failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_depth_from_one_failed_bit(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 4: DEPTH FROM FIRST FAILED BRANCH BIT")
    print("=" * 90)

    failures = 0

    for s in states:
        W = w(s)
        actual = direct_depth(s)

        branch_limit = min(W, 32)

        first_failure = None

        for bit in range(branch_limit):
            # bit corresponds to modulus 2^(bit+1)
            bits = bit + 1

            if not prefix_equal_p(s, bits):
                first_failure = bits
                break

        if first_failure is None:
            predicted = W
        else:
            predicted = min(first_failure + 1, W)

        if predicted != actual:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s.n} frame={s.frame} "
                    f"W={W} first_failure={first_failure} "
                    f"predicted={predicted} actual={actual}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_branch_code_is_exact(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 5: FIRST-DIFFERING-BIT CODE")
    print("=" * 90)

    """
    Encode the branch using:

        code = tp

    but also construct it only from the first failed prefix.
    """

    failures = 0
    distribution = Counter()

    for s in states:
        T = tp(s)
        W = w(s)

        # If p exactly equals p0, there is no finite branch failure.
        if T >= INF // 2:
            code = None
        else:
            code = T

        observed = None

        for k in range(1, min(W, 32) + 1):
            if not prefix_equal_p(s, k):
                observed = k
                break

        if observed is None:
            reconstructed_tp = INF
        else:
            reconstructed_tp = observed - 1

        # We only need the truncated value up to W.
        expected = min(T, W)

        if min(reconstructed_tp, W) != expected:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s.n} frame={s.frame} "
                    f"T={T} W={W} "
                    f"reconstructed={reconstructed_tp}"
                )

        if code is None:
            distribution[("no-break", W)] += 1
        else:
            distribution[("break", min(T, W))] += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    print("\n    branch-code distribution:")
    for key, count in sorted(distribution.items()):
        print(f"        {key}: {count}")

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_threshold_sets(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 6: THRESHOLD SET NESTING")
    print("=" * 90)

    """
    Verify:

        S_(d+1) subset S_d

    where:

        S_d = states satisfying
              p == p0 mod 2^(d-1)
              and
              n == -c mod 2^d.
    """

    failures = 0

    max_d = 12

    for s in states:
        for d in range(1, max_d):
            current = (
                prefix_equal_p(s, d)
                and prefix_equal_n(s, d)
            )

            nxt = (
                prefix_equal_p(s, d + 1)
                and prefix_equal_n(s, d + 1)
            )

            if nxt and not current:
                failures += 1

                if failures <= 20:
                    print(
                        f"    nesting violation n={s.n} "
                        f"frame={s.frame} d={d}"
                    )

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_threshold_cardinalities(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 7: THRESHOLD CARDINALITIES")
    print("=" * 90)

    for frame in ("A", "B"):
        print(f"\nFRAME {frame}")

        for d in range(1, 13):
            total = 0

            for s in states:
                if s.frame != frame:
                    continue

                if (
                    prefix_equal_p(s, d)
                    and prefix_equal_n(s, d)
                ):
                    total += 1

            print(
                f"    d={d:2d} surviving={total}"
            )

    return 0


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(states: list[State]) -> None:
    print("\n" + "=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = {
        9, 15, 21, 33, 39, 57, 69,
        77, 87, 93, 111, 141, 183,
        213, 485879, 5579767
    }

    shown = 0

    for s in states:
        if s.n not in wanted:
            continue

        T = tp(s)
        W = w(s)
        D = direct_depth(s)

        print(f"\nn={s.n} p={s.p} q={s.q} frame={s.frame}")
        print(f"    p0={s.p0} c={s.c}")
        print(f"    tp={T}")
        print(f"    w={W}")
        print(f"    depth={D}")

        for d in range(1, min(D + 2, 10)):
            p_ok = prefix_equal_p(s, d - 1)
            p_fail = not prefix_equal_p(s, d)
            n_ok = prefix_equal_n(s, d)

            print(
                f"    d={d}: "
                f"p-prefix(d-1)={p_ok} "
                f"p-first-fails(d)={p_fail} "
                f"n-prefix={n_ok}"
            )

        shown += 1

        if shown >= 16:
            break


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 654 START")
    print("=" * 90)

    print(f"\nprime limit={PRIME_LIMIT}")

    primes = sieve_odd_primes(PRIME_LIMIT)
    states = build_states(primes)

    print(f"odd primes={len(primes)}")
    print(f"semiprimes={len(states)}")

    failures = 0

    failures += test_baseline(states)
    failures += test_branch_first_mismatch(states)
    failures += test_n_bound_branch(states)
    failures += test_single_bit_partition(states)
    failures += test_depth_from_one_failed_bit(states)
    failures += test_branch_code_is_exact(states)
    failures += test_threshold_sets(states)

    test_threshold_cardinalities(states)
    print_examples(states)

    print("\n" + "=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
Experiment 653 established the exact threshold theorem:

    depth >= d

        iff

    p == p0 (mod 2^(d-1))
    AND
    n == -c (mod 2^d).

Experiment 654 rewrites this as a one-bit stopping rule.

For a finite:

    tp = v2(p-p0)

the branch survives exactly through the prefix

    p == p0 (mod 2^tp)

and fails at

    2^(tp+1).

Therefore:

    branch depth = tp + 1.

The n-side survives exactly through:

    w = v2(n+c).

Hence:

    depth = min(tp+1, w).

The important reduction is that the factor-side state is
not a general 2-adic value anymore.

It is simply:

    FIRST DIFFERING BIT OF p FROM p0.

Equivalently, the branch is a path

    p0
     |
     +-- matching bit
     |
     +-- matching bit
     |
     +-- first differing bit

and the depth is the first stopping event between that
branch and the n-side congruence.

The experiment therefore tests whether the complete
integer depth can be represented as a single terminal
branch bit rather than a stored valuation.
"""
    )

    print(f"\nTOTAL FAILURES={failures}")

    if failures == 0:
        print("STATUS=ALL TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 654 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
