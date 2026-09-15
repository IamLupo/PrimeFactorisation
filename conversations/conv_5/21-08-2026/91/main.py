#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 653
==========================================================================================

ALGEBRAIC 2-ADIC THRESHOLD THEOREM

Goal:
    Stop searching for empirical signatures.

    Directly verify, for every semiprime state and every relevant
    level d, the exact equivalence:

        depth >= d

        iff

        p == p0 (mod 2^(d-1))
        AND
        n == -c (mod 2^d).

Frames:

    FRAME A:
        A = p - 3
        B = q + 3
        p0 = 3
        q0 = -3
        c  = 9

    FRAME B:
        A = p + 1
        B = q - 3
        p0 = -1
        q0 = 3
        c  = 3

Also tests:

    1. direct transformed-numerator depth
    2. depth = min(v2(p-p0)+1, v2(n+c))
    3. threshold equivalence
    4. q-prefix is redundant once p-prefix+n-prefix are known
    5. n-prefix alone is insufficient
    6. first failed congruence exactly identifies the depth
    7. threshold sets are nested
    8. constructive modular factor lift
    9. symbolic-style integer identities

This is intended as a proof-oriented experiment rather than
another brute-force formula search.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


INF = 10**9


# =========================================================================================
# BASIC NUMBER THEORY
# =========================================================================================

def v2(x: int) -> int:
    """2-adic valuation of nonzero integer x."""
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def odd_primes(limit: int) -> list[int]:
    """Return odd primes <= limit."""
    if limit < 3:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [p for p in range(3, limit + 1, 2) if sieve[p]]


# =========================================================================================
# STATE
# =========================================================================================

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


def build_states(limit: int) -> list[State]:
    primes = odd_primes(limit)

    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            # Historical frame convention established by previous experiments.
            if n % 4 == 3:
                frame = "A"

                A = p - 3
                B = q + 3

                p0 = 3
                q0 = -3
                c = 9

            else:
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


# =========================================================================================
# DIRECT STRUCTURAL QUANTITIES
# =========================================================================================

def transformed_numerators(s: State) -> tuple[int, int]:
    """
    Return (2X, 2Y).

    FRAME A:
        2X = B-A
        2Y = A+B

    FRAME B:
        2X = 3A-B
        2Y = 3A+B
    """
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
    """
    Exact depth from the transformed numerator theorem.
    """
    u, v = transformed_numerators(s)
    return min(v2(u), v2(v))


def tp_value(s: State) -> int:
    """
    Branch valuation:

        tp = v2(p-p0)
    """
    return v2(s.p - s.p0)


def tq_value(s: State) -> int:
    """
    Other factor branch valuation:

        tq = v2(q-q0)
    """
    return v2(s.q - s.q0)


def w_value(s: State) -> int:
    """
    n-only valuation:

        w = v2(n+c)
    """
    return v2(s.n + s.c)


def predicted_depth(s: State) -> int:
    """
    Experiment 651 law:

        depth = min(tp+1, w)
    """
    return min(tp_value(s) + 1, w_value(s))


# =========================================================================================
# THRESHOLD PREDICATES
# =========================================================================================

def p_prefix(s: State, d: int) -> bool:
    """
    Prime branch survives to threshold d iff

        p == p0 (mod 2^(d-1)).
    """
    modulus = 1 << (d - 1)
    return (s.p - s.p0) % modulus == 0


def q_prefix(s: State, d: int) -> bool:
    """
    Other factor satisfies its corresponding prefix:

        q == q0 (mod 2^(d-1)).
    """
    modulus = 1 << (d - 1)
    return (s.q - s.q0) % modulus == 0


def n_prefix(s: State, d: int) -> bool:
    """
    Observable-side prefix:

        n == -c (mod 2^d).
    """
    modulus = 1 << d
    return (s.n + s.c) % modulus == 0


def threshold_rhs(s: State, d: int) -> bool:
    """
    Candidate threshold theorem.
    """
    return p_prefix(s, d) and n_prefix(s, d)


# =========================================================================================
# TEST 0
# =========================================================================================

def test_odd_domain(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for s in states:
        if s.p % 2 == 0 or s.q % 2 == 0 or s.n % 2 == 0:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =========================================================================================
# TEST 1
# =========================================================================================

def test_direct_depth(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 1: DIRECT DEPTH VS EXPERIMENT 651 LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        d1 = direct_depth(s)
        d2 = predicted_depth(s)

        if d1 != d2:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s.n} frame={s.frame} "
                    f"p={s.p} q={s.q} "
                    f"direct={d1} predicted={d2}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =========================================================================================
# TEST 2
# =========================================================================================

def test_threshold_theorem(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 2: EXACT THRESHOLD THEOREM")
    print("=" * 90)

    failures = 0
    checks = 0

    for s in states:
        depth = direct_depth(s)

        # Test one level below, at depth, and above.
        test_max = min(depth + 2, 32)

        for d in range(1, test_max + 1):
            checks += 1

            lhs = depth >= d
            rhs = threshold_rhs(s, d)

            if lhs != rhs:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch n={s.n} frame={s.frame} "
                        f"d={d} depth={depth} "
                        f"lhs={lhs} rhs={rhs}"
                    )

    print(f"checks={checks}")
    print(f"failures={failures}")

    return failures


# =========================================================================================
# TEST 3
# =========================================================================================

def test_q_prefix_redundancy(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 3: q-PREFIX REDUNDANCY")
    print("=" * 90)

    """
    The threshold theorem claims:

        p-prefix + n-prefix
            =>
        q-prefix.

    This means the second prime factor is not an independent
    branch variable once p and n are known.

    We test:

        p == p0 (mod 2^(d-1))
        AND
        n == -c (mod 2^d)

        =>
        q == q0 (mod 2^(d-1)).
    """

    failures = 0
    checks = 0

    for s in states:
        depth = direct_depth(s)

        for d in range(1, min(depth + 2, 32) + 1):
            checks += 1

            lhs = p_prefix(s, d) and n_prefix(s, d)
            rhs = q_prefix(s, d)

            if lhs and not rhs:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch n={s.n} frame={s.frame} d={d} "
                        f"p={s.p} q={s.q}"
                    )

    print(f"checks={checks}")
    print(f"failures={failures}")

    return failures


# =========================================================================================
# TEST 4
# =========================================================================================

def test_n_prefix_not_sufficient(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 4: n-PREFIX ALONE IS NOT SUFFICIENT")
    print("=" * 90)

    """
    Group states by:

        (frame, d, n mod 2^d)

    and look for different actual depths.

    This demonstrates precisely why the factor-prefix is
    necessary.
    """

    buckets: dict[tuple[str, int, int], set[int]] = defaultdict(set)

    for s in states:
        depth = direct_depth(s)

        # Only test interesting levels.
        for d in range(2, min(depth + 2, 16) + 1):
            residue = s.n % (1 << d)
            buckets[(s.frame, d, residue)].add(depth)

    ambiguous = {
        key: vals
        for key, vals in buckets.items()
        if len(vals) > 1
    }

    print(f"buckets={len(buckets)}")
    print(f"ambiguous n-only buckets={len(ambiguous)}")

    shown = 0

    for (frame, d, residue), depths in sorted(ambiguous.items()):
        print(
            f"    frame={frame} d={d} "
            f"n mod 2^d={residue} depths={sorted(depths)}"
        )

        shown += 1
        if shown >= 12:
            break

    return 0


# =========================================================================================
# TEST 5
# =========================================================================================

def test_first_failed_congruence(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 5: DEPTH AS FIRST FAILED CONGRUENCE")
    print("=" * 90)

    failures = 0

    for s in states:
        depth = direct_depth(s)

        # The theorem predicts that d <= depth survives and
        # d = depth+1 fails.
        for d in range(1, min(depth + 2, 32) + 1):
            survives = threshold_rhs(s, d)

            expected = d <= depth

            if survives != expected:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch n={s.n} frame={s.frame} "
                        f"d={d} depth={depth} "
                        f"survives={survives} expected={expected}"
                    )

        # Explicitly inspect the first failed threshold.
        if depth < 32:
            if threshold_rhs(s, depth + 1):
                failures += 1

                if failures <= 20:
                    print(
                        f"    first-failure mismatch n={s.n} "
                        f"frame={s.frame} depth={depth}"
                    )

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =========================================================================================
# TEST 6
# =========================================================================================

def test_constructive_inverse(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 6: CONSTRUCTIVE MODULAR FACTOR LIFT")
    print("=" * 90)

    """
    For each state and each d <= w:

        FRAME A:
            p0 = 3
            q0 = -3

        FRAME B:
            p0 = -1
            q0 = 3

    We verify directly that the threshold conditions imply
    the required factor residue.

    Then we construct a compatible q residue from n and p0:

        q ≡ n * p0^{-1} (mod 2^d)

    and verify it equals q0.
    """

    failures = 0
    checks = 0

    examples = []

    for s in states:
        w = w_value(s)

        for d in range(1, min(w, 16) + 1):
            checks += 1

            modulus = 1 << d

            # p0 must be invertible modulo 2^d.
            p0_mod = s.p0 % modulus

            try:
                p0_inv = pow(p0_mod, -1, modulus)
            except ValueError:
                failures += 1
                continue

            constructed_q = (s.n * p0_inv) % modulus
            expected_q = s.q0 % modulus

            if constructed_q != expected_q:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch n={s.n} frame={s.frame} d={d} "
                        f"constructed_q={constructed_q} "
                        f"expected_q={expected_q}"
                    )

            if len(examples) < 10:
                examples.append(
                    (
                        s.frame,
                        d,
                        s.n,
                        modulus,
                        s.p0,
                        constructed_q,
                        expected_q,
                    )
                )

    print(f"checks={checks}")
    print(f"failures={failures}")

    print("\n    sample constructive lifts:")
    for item in examples:
        frame, d, n, modulus, p0, q_constructed, q_expected = item

        print(
            f"        frame={frame} d={d} "
            f"n={n} mod={modulus} "
            f"p0={p0} "
            f"q={q_constructed} "
            f"expected={q_expected}"
        )

    return failures


# =========================================================================================
# TEST 7
# =========================================================================================

def test_symbolic_integer_identities(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 7: INTEGER IDENTITIES BEHIND THE THRESHOLD")
    print("=" * 90)

    """
    Frame A:

        A = p-3
        B = q+3

        A+B = p+q
        B-A = q-p+6

    Frame B:

        A = p+1
        B = q-3

        3A+B = 3p+q
        3A-B = 3p-q+6

    These identities make the threshold reduction explicit.
    """

    failures = 0

    for s in states:
        if s.frame == "A":
            u = s.A + s.B
            v = s.B - s.A

            expected_u = s.p + s.q
            expected_v = s.q - s.p + 6

        else:
            u = 3 * s.A + s.B
            v = 3 * s.A - s.B

            expected_u = 3 * s.p + s.q
            expected_v = 3 * s.p - s.q + 6

        if u != expected_u or v != expected_v:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s.n} frame={s.frame} "
                    f"u={u} expected_u={expected_u} "
                    f"v={v} expected_v={expected_v}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    return failures


# =========================================================================================
# TEST 8
# =========================================================================================

def test_threshold_counts(states: list[State]) -> int:
    print("\n" + "=" * 90)
    print("TEST 8: THRESHOLD SURVIVAL COUNTS")
    print("=" * 90)

    """
    Count how many actual states survive each threshold.

    This gives the empirical distribution of the two stopping
    mechanisms without doing any signature search.
    """

    counts: dict[tuple[str, int], int] = defaultdict(int)

    max_d = 12

    for s in states:
        depth = direct_depth(s)

        for d in range(1, max_d + 1):
            if depth >= d:
                counts[(s.frame, d)] += 1

    for frame in ("A", "B"):
        print(f"\nFRAME {frame}")

        for d in range(1, max_d + 1):
            print(
                f"    d={d:2d} "
                f"surviving={counts[(frame, d)]}"
            )

    return 0


# =========================================================================================
# EXAMPLES
# =========================================================================================

def print_examples(states: list[State]) -> None:
    print("\n" + "=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = {9, 15, 21, 33, 39, 57, 69, 77, 87, 93, 111, 141, 183, 213}

    shown = 0

    for s in states:
        if s.n not in wanted:
            continue

        depth = direct_depth(s)
        tp = tp_value(s)
        w = w_value(s)

        print(f"\nn={s.n} p={s.p} q={s.q} frame={s.frame}")
        print(
            f"    p0={s.p0} q0={s.q0} c={s.c}"
        )
        print(
            f"    tp=v2(p-p0)={tp} "
            f"w=v2(n+c)={w}"
        )
        print(
            f"    depth={depth} "
            f"min(tp+1,w)={min(tp + 1, w)}"
        )

        limit = min(depth + 1, 8)

        for d in range(1, limit + 1):
            pp = p_prefix(s, d)
            qp = q_prefix(s, d)
            np = n_prefix(s, d)
            th = threshold_rhs(s, d)

            print(
                f"    d={d}: "
                f"p-prefix={pp} "
                f"q-prefix={qp} "
                f"n-prefix={np} "
                f"threshold={th}"
            )

        shown += 1

        if shown >= 14:
            break


# =========================================================================================
# MAIN
# =========================================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 653 START")
    print("=" * 90)

    PRIME_LIMIT = 6000

    print(f"\nprime limit={PRIME_LIMIT}")

    states = build_states(PRIME_LIMIT)

    primes = odd_primes(PRIME_LIMIT)

    print(f"odd primes={len(primes)}")
    print(f"semiprimes={len(states)}")

    total_failures = 0

    total_failures += test_odd_domain(states)
    total_failures += test_direct_depth(states)
    total_failures += test_threshold_theorem(states)
    total_failures += test_q_prefix_redundancy(states)

    # Diagnostic only; this intentionally expects ambiguity.
    test_n_prefix_not_sufficient(states)

    total_failures += test_first_failed_congruence(states)
    total_failures += test_constructive_inverse(states)
    total_failures += test_symbolic_integer_identities(states)

    test_threshold_counts(states)
    print_examples(states)

    print("\n" + "=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
The experiment tests the exact threshold form:

    depth >= d

        iff

    p == p0 (mod 2^(d-1))

    and

    n == -c (mod 2^d).

FRAME A:

    p0 =  3
    q0 = -3
    c  =  9

FRAME B:

    p0 = -1
    q0 =  3
    c  =  3.

Equivalently:

    depth =
        min(
            v2(p-p0) + 1,
            v2(n+c)
        ).

The important algebraic point is that q does not need to
be stored independently.

Because p is odd and therefore invertible modulo 2^d:

    q = n * p^(-1) (mod 2^d).

Once:

    p == p0 (mod 2^(d-1))

and

    n == -c (mod 2^d),

the required q-prefix follows automatically.

Therefore the surviving factor information is exactly one
2-adic branch:

    p -> p0.

The remaining question is no longer whether the threshold
law is true.  It is whether the first differing bit of p
from p0 has some independent number-theoretic description.
"""
    )

    print(f"TOTAL FAILURES={total_failures}")

    if total_failures == 0:
        print("STATUS=ALL THEOREM TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 653 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
