#!/usr/bin/env python3
# =============================================================================
# EXPERIMENT 672
# =============================================================================
#
# N-ONLY GCD / POLYNOMIAL RECONSTRUCTION SEARCH
#
# Established:
#
#   FRAME A:
#       C = q + 3
#       c = 9
#
#       depth = v2(gcd(2C, n+9))
#
#   FRAME B:
#       C = p + 1
#       c = 3
#
#       depth = v2(gcd(2C, n+3))
#
# Equivalent:
#
#       depth = min(v2(C)+1, v2(n+c))
#
# This experiment searches whether the C-dependent part
# can be replaced by an n-only arithmetic expression.
#
# Tested families:
#
#   1. gcd(n+c, n+r)
#   2. gcd(n+c, a*n+b)
#   3. gcd(n+c, n^2+a*n+b)
#   4. Small n-only valuation signatures
#   5. Direct same-v2(n+c) ambiguity
#
# No O(N^2) state separation is used.
#
# =============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import gcd


PRIME_LIMIT = 6000
INF = 10**9


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def v2(x: int) -> int:
    """Return v2(x), with v2(0)=INF."""
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def sieve(limit: int) -> list[int]:
    """Return odd primes <= limit."""
    if limit < 3:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0] = 0
    is_prime[1] = 0

    p = 2
    while p * p <= limit:
        if is_prime[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            is_prime[start : limit + 1 : p] = b"\x00" * count
        p += 1

    return [p for p in range(3, limit + 1, 2) if is_prime[p]]


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
    C: int

    c: int

    depth: int
    w: int
    vc: int


# =============================================================================
# STATE CONSTRUCTION
# =============================================================================

def make_state(p: int, q: int) -> State:
    n = p * q

    # Established frame convention:
    #
    #   A: n == 3 (mod 4)
    #   B: n == 1 (mod 4)

    if n % 4 == 3:
        frame = "A"

        A = p - 3
        B = q + 3
        C = B
        c = 9

    else:
        frame = "B"

        A = p + 1
        B = q - 3
        C = A
        c = 3

    depth = v2(gcd(2 * C, n + c))
    w = v2(n + c)
    vc = v2(C)

    return State(
        n=n,
        p=p,
        q=q,
        frame=frame,
        A=A,
        B=B,
        C=C,
        c=c,
        depth=depth,
        w=w,
        vc=vc,
    )


def build_states(primes: list[int]) -> list[State]:
    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(make_state(p, q))

    return states


# =============================================================================
# BASELINE
# =============================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = v2(gcd(2 * s.C, s.n + s.c))

        if predicted != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# GENERIC CANDIDATE
# =============================================================================

def predict_from_candidate(s: State, candidate: int) -> int:
    """
    Candidate n-only reconstruction:

        v2(gcd(n+c, candidate))

    This intentionally does NOT include the explicit factor 2.
    The search is testing whether the n-only candidate can
    reproduce the final depth directly.
    """
    return v2(gcd(s.n + s.c, candidate))


# =============================================================================
# TEST 1
# =============================================================================

def search_shift_family(
    states: list[State],
    r_min: int,
    r_max: int,
) -> None:
    print("=" * 90)
    print("TEST 1: SHIFT FAMILY")
    print("=" * 90)

    exact: list[int] = []

    for r in range(r_min, r_max + 1):
        failures = 0
        first = None

        for s in states:
            candidate = s.n + r
            predicted = predict_from_candidate(s, candidate)

            if predicted != s.depth:
                failures += 1

                if first is None:
                    first = (s, predicted)

                # Early exit:
                break

        if failures == 0:
            exact.append(r)

        print(
            f"r={r:4d} "
            f"failures={failures:1d} "
            f"exact={'YES' if failures == 0 else 'NO'}"
        )

        if first is not None and failures == 1:
            s, predicted = first

            print(
                f"    first mismatch: "
                f"n={s.n} p={s.p} q={s.q} "
                f"frame={s.frame} "
                f"depth={s.depth} "
                f"predicted={predicted}"
            )

    print()
    print(f"exact shifts={exact}")
    print()


# =============================================================================
# TEST 2
# =============================================================================

def search_affine_family(
    states: list[State],
    a_min: int,
    a_max: int,
    b_min: int,
    b_max: int,
) -> None:
    print("=" * 90)
    print("TEST 2: AFFINE n-ONLY FAMILY")
    print("=" * 90)

    tested = 0
    exact: list[tuple[int, int]] = []

    for a in range(a_min, a_max + 1):
        for b in range(b_min, b_max + 1):
            tested += 1

            exact_here = True

            for s in states:
                candidate = a * s.n + b
                predicted = predict_from_candidate(s, candidate)

                if predicted != s.depth:
                    exact_here = False
                    break

            if exact_here:
                exact.append((a, b))

    print(f"tested affine expressions={tested}")
    print(f"exact affine expressions={len(exact)}")

    if exact:
        print("exact expressions:")
        for a, b in exact:
            print(f"    F(n) = ({a})*n + ({b})")

    print()


# =============================================================================
# TEST 3
# =============================================================================

def search_quadratic_family(
    states: list[State],
    a_min: int,
    a_max: int,
    b_min: int,
    b_max: int,
) -> None:
    print("=" * 90)
    print("TEST 3: QUADRATIC n-ONLY FAMILY")
    print("=" * 90)

    tested = 0
    exact: list[tuple[int, int]] = []

    for a in range(a_min, a_max + 1):
        for b in range(b_min, b_max + 1):
            tested += 1

            exact_here = True

            for s in states:
                candidate = s.n * s.n + a * s.n + b
                predicted = predict_from_candidate(s, candidate)

                if predicted != s.depth:
                    exact_here = False
                    break

            if exact_here:
                exact.append((a, b))

    print(f"tested quadratic expressions={tested}")
    print(f"exact quadratic expressions={len(exact)}")

    if exact:
        print("exact expressions:")
        for a, b in exact:
            print(f"    F(n) = n^2 + ({a})*n + ({b})")

    print()


# =============================================================================
# N-ONLY SIGNATURE
# =============================================================================

def n_only_signature(s: State) -> tuple[int, ...]:
    """
    Nearby 2-adic information around n+c.

    The offsets are chosen to remain small and entirely
    n-only.
    """
    base = s.n + s.c

    offsets = (
        0,
        2,
        4,
        6,
        8,
        10,
        12,
        16,
        24,
        32,
    )

    return tuple(v2(base + off) for off in offsets)


def test_n_only_signature(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: N-ONLY MULTI-VALUATION SIGNATURE")
    print("=" * 90)

    buckets: dict[
        tuple[str, tuple[int, ...]],
        set[int],
    ] = defaultdict(set)

    witnesses: dict[
        tuple[str, tuple[int, ...]],
        State,
    ] = {}

    for s in states:
        signature = (s.frame, n_only_signature(s))

        buckets[signature].add(s.depth)
        witnesses.setdefault(signature, s)

    ambiguous = {
        sig: depths
        for sig, depths in buckets.items()
        if len(depths) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print()
        print("first counterexamples:")

        for index, (sig, depths) in enumerate(ambiguous.items()):
            if index >= 15:
                break

            frame, signature = sig
            s = witnesses[sig]

            print(
                f"    frame={frame} "
                f"depths={sorted(depths)}"
            )

            print(
                f"        n={s.n} "
                f"p={s.p} "
                f"q={s.q}"
            )

            print(
                f"        signature={signature}"
            )

    print()

    return len(ambiguous)


# =============================================================================
# TEST 5
# =============================================================================

def test_same_w_different_depth(states: list[State]) -> None:
    print("=" * 90)
    print("TEST 5: SAME n-SIDE VALUATION, DIFFERENT DEPTH")
    print("=" * 90)

    buckets: dict[tuple[str, int], list[State]] = defaultdict(list)

    for s in states:
        buckets[(s.frame, s.w)].append(s)

    ambiguous_count = 0

    for (frame, w), group in buckets.items():
        depths = {s.depth for s in group}

        if len(depths) <= 1:
            continue

        ambiguous_count += 1

        if ambiguous_count <= 15:
            print(
                f"    frame={frame} "
                f"w={w} "
                f"depths={sorted(depths)}"
            )

            shown = 0
            seen_depths: set[int] = set()

            for s in group:
                if s.depth in seen_depths:
                    continue

                seen_depths.add(s.depth)

                print(
                    f"        n={s.n:<10} "
                    f"p={s.p:<5} "
                    f"q={s.q:<5} "
                    f"depth={s.depth:<3} "
                    f"v2(C)={s.vc}"
                )

                shown += 1

                if shown >= 4:
                    break

    print()
    print(f"ambiguous w-buckets={ambiguous_count}")
    print()


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
        87,
        93,
        111,
        141,
        183,
        213,
        485879,
        5579767,
    ]

    by_n = {s.n: s for s in states}

    shown = 0

    for n in wanted:
        s = by_n.get(n)

        if s is None:
            continue

        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(f"    A={s.A}")
        print(f"    B={s.B}")
        print(f"    C={s.C}")
        print(f"    c={s.c}")
        print(f"    n+c={s.n + s.c}")
        print(f"    v2(C)={s.vc}")
        print(f"    v2(n+c)={s.w}")
        print(f"    depth={s.depth}")

        canonical = gcd(2 * s.C, s.n + s.c)

        print(f"    canonical gcd={canonical}")
        print(f"    v2(canonical gcd)={v2(canonical)}")

        print()

        shown += 1

    print(f"examples shown={shown}")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 672 START")
    print("=" * 90)
    print()

    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve(PRIME_LIMIT)

    print(f"odd primes={len(primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    # -------------------------------------------------------------------------
    # Baseline
    # -------------------------------------------------------------------------

    total_failures += test_baseline(states)

    # -------------------------------------------------------------------------
    # Shift family
    # -------------------------------------------------------------------------

    search_shift_family(
        states,
        r_min=-64,
        r_max=64,
    )

    # -------------------------------------------------------------------------
    # Affine family
    # -------------------------------------------------------------------------

    search_affine_family(
        states,
        a_min=-3,
        a_max=3,
        b_min=-32,
        b_max=32,
    )

    # -------------------------------------------------------------------------
    # Quadratic family
    # -------------------------------------------------------------------------

    search_quadratic_family(
        states,
        a_min=-8,
        a_max=8,
        b_min=-32,
        b_max=32,
    )

    # -------------------------------------------------------------------------
    # Multi-valuation signature
    # -------------------------------------------------------------------------

    ambiguous_signature_count = test_n_only_signature(states)

    # -------------------------------------------------------------------------
    # Direct information-loss demonstration
    # -------------------------------------------------------------------------

    test_same_w_different_depth(states)

    # -------------------------------------------------------------------------
    # Examples
    # -------------------------------------------------------------------------

    print_examples(states)

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Established canonical laws:")
    print()

    print("    FRAME A:")
    print("        C = q + 3")
    print("        c = 9")
    print("        depth = v2(gcd(2C, n+9))")
    print()

    print("    FRAME B:")
    print("        C = p + 1")
    print("        c = 3")
    print("        depth = v2(gcd(2C, n+3))")
    print()

    print("Equivalent valuation form:")
    print()
    print("    depth = min(v2(C)+1, v2(n+c))")
    print()

    print("Remaining question:")
    print()
    print("    Can the factor-dependent C be replaced")
    print("    by an n-only arithmetic construction?")
    print()

    print("This experiment searches:")
    print()
    print("    fixed shifts")
    print("    affine polynomials")
    print("    quadratic polynomials")
    print("    multi-valuation n-only signatures")
    print()

    print(
        "An ambiguous n-only signature means the tested "
        "information is insufficient."
    )

    print(
        "It does NOT prove that no arbitrary n-only "
        "function could exist."
    )

    print()

    print(
        f"n-only signature ambiguous buckets="
        f"{ambiguous_signature_count}"
    )

    print()

    print(f"TOTAL TEST FAILURES={total_failures}")

    print("=" * 90)
    print("EXPERIMENT 672 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()