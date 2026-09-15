#!/usr/bin/env python3

from collections import defaultdict, Counter


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 6000
MAX_LEVEL = 24
INF = 10**9


# =============================================================================
# HEADER
# =============================================================================

def banner(title="EXPERIMENT 657 START"):
    print("=" * 90)
    print(title)
    print("=" * 90)


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def sieve(limit: int) -> list[int]:
    composite = bytearray(limit + 1)
    primes = []

    for p in range(2, limit + 1):
        if not composite[p]:
            primes.append(p)

            if p * p <= limit:
                composite[p * p : limit + 1 : p] = b"\x01" * (
                    ((limit - p * p) // p) + 1
                )

    return primes


def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


# =============================================================================
# FRAME DEFINITIONS
# =============================================================================

def frame_of(n: int) -> str:
    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(f"Unexpected odd n={n}")


def frame_parameters(frame: str):
    if frame == "A":
        # A = p - 3
        # B = q + 3
        # n + 9
        return 3, -3, 9

    if frame == "B":
        # A = p + 1
        # B = q - 3
        # n + 3
        return -1, 3, 3

    raise ValueError(frame)


# =============================================================================
# STATE GENERATION
# =============================================================================

def build_states(primes: list[int]):
    odd_primes = [p for p in primes if p & 1]

    states = []

    for i, p in enumerate(odd_primes):
        for q in odd_primes[i:]:
            n = p * q

            frame = frame_of(n)
            p0, q0, c = frame_parameters(frame)

            A = p - p0
            B = q - q0

            if frame == "A":
                # 2X = B-A
                # 2Y = A+B
                two_x = B - A
                two_y = A + B
            else:
                # 2X = 3A-B
                # 2Y = 3A+B
                two_x = 3 * A - B
                two_y = 3 * A + B

            depth = min(v2(two_x), v2(two_y))

            tp = v2(p - p0)
            tq = v2(q - q0)
            w = v2(n + c)

            states.append(
                {
                    "p": p,
                    "q": q,
                    "n": n,
                    "frame": frame,
                    "p0": p0,
                    "q0": q0,
                    "c": c,
                    "A": A,
                    "B": B,
                    "tp": tp,
                    "tq": tq,
                    "w": w,
                    "depth": depth,
                }
            )

    return states


# =============================================================================
# PREFIX PREDICATES
# =============================================================================

def p_prefix(state, d: int) -> bool:
    """
    Threshold theorem uses:

        p == p0 (mod 2^(d-1))
    """
    if d <= 1:
        return True

    modulus = 1 << (d - 1)

    return (state["p"] - state["p0"]) % modulus == 0


def q_prefix(state, d: int) -> bool:
    """
    q-prefix at level d:

        q == q0 (mod 2^(d-1))
    """
    if d <= 1:
        return True

    modulus = 1 << (d - 1)

    return (state["q"] - state["q0"]) % modulus == 0


def n_prefix(state, d: int) -> bool:
    """
    n-prefix at level d:

        n == -c (mod 2^d)
    """
    if d <= 0:
        return True

    modulus = 1 << d

    return (state["n"] + state["c"]) % modulus == 0


def threshold(state, d: int) -> bool:
    return p_prefix(state, d) and n_prefix(state, d)


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states):
    print("=" * 90)
    print("TEST 0: EXPERIMENT 656 BASELINE")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = min(s["tp"] + 1, s["w"])

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch "
                    f"n={s['n']} "
                    f"p={s['p']} "
                    f"q={s['q']} "
                    f"tp={s['tp']} "
                    f"w={s['w']} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_threshold_theorem(states):
    print("=" * 90)
    print("TEST 1: EXACT THRESHOLD THEOREM")
    print("=" * 90)

    failures = 0
    checks = 0

    for s in states:

        # Test a little beyond the true depth so the first failed
        # threshold is also checked.
        max_d = min(MAX_LEVEL, s["depth"] + 2)

        for d in range(1, max_d + 1):
            predicted = threshold(s, d)
            actual = s["depth"] >= d

            checks += 1

            if predicted != actual:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch "
                        f"n={s['n']} "
                        f"d={d} "
                        f"p-prefix={p_prefix(s, d)} "
                        f"n-prefix={n_prefix(s, d)} "
                        f"predicted={predicted} "
                        f"actual={actual}"
                    )

    print(f"checks={checks}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_q_redundancy(states):
    print("=" * 90)
    print("TEST 2: Q-PREFIX CONSEQUENCES")
    print("=" * 90)

    """
    Correct distinction:

        p == p0 mod 2^(d-1)
        n == -c mod 2^d

    implies:

        q == q0 mod 2^(d-1)

    because q is determined by n * p^{-1} modulo 2^(d-1).

    But p-prefix + q-prefix alone generally gives only
    n == -c mod 2^(d-1), NOT modulo 2^d.

    So we test the correct one-directional implication.
    """

    failures = 0
    checks = 0

    for s in states:

        max_d = min(MAX_LEVEL, s["depth"] + 2)

        for d in range(1, max_d + 1):

            pp = p_prefix(s, d)
            np = n_prefix(s, d)
            qp = q_prefix(s, d)

            # Correct implication:
            # p-prefix(d) AND n-prefix(d) -> q-prefix(d)
            predicted_q = pp and np

            checks += 1

            if predicted_q and not qp:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch "
                        f"n={s['n']} "
                        f"d={d} "
                        f"p={s['p']} "
                        f"q={s['q']} "
                        f"p-prefix={pp} "
                        f"n-prefix={np} "
                        f"q-prefix={qp}"
                    )

    print(f"checks={checks}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_local_bit_automaton(states):
    print("=" * 90)
    print("TEST 3: TWO-CHANNEL BIT AUTOMATON")
    print("=" * 90)

    """
    Maintain only two states:

        p_alive
        n_alive

    At every level:

        p_alive:
            has p matched p0 at every required previous bit?

        n_alive:
            has n matched -c at every required previous bit?

    The threshold is simply:

        p_alive AND n_alive
    """

    failures = 0
    checks = 0

    for s in states:

        p_alive = True
        n_alive = True

        max_d = min(MAX_LEVEL, s["depth"] + 2)

        for d in range(1, max_d + 1):

            # -------------------------------------------------------------
            # New p bit
            #
            # p-prefix(d) examines the first d-1 bits.
            # Thus the newly introduced p bit at level d is bit d-2.
            # -------------------------------------------------------------
            if d >= 2:
                bit_position = d - 2
                mask = 1 << bit_position

                p_bit = (s["p"] & mask) != 0
                p0_bit = (s["p0"] & mask) != 0

                if p_bit != p0_bit:
                    p_alive = False

            # -------------------------------------------------------------
            # New n bit
            #
            # n-prefix(d) examines d bits.
            # Thus the newly introduced n bit is bit d-1.
            # -------------------------------------------------------------
            bit_position = d - 1
            mask = 1 << bit_position

            n_bit = (s["n"] & mask) != 0
            target_bit = ((-s["c"]) & mask) != 0

            if n_bit != target_bit:
                n_alive = False

            predicted = p_alive and n_alive
            actual = s["depth"] >= d

            checks += 1

            if predicted != actual:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch "
                        f"n={s['n']} "
                        f"d={d} "
                        f"p_alive={p_alive} "
                        f"n_alive={n_alive} "
                        f"predicted={predicted} "
                        f"actual={actual}"
                    )

    print(f"checks={checks}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_level_signature(states):
    print("=" * 90)
    print("TEST 4: LEVEL-INDEXED LOCAL SIGNATURE")
    print("=" * 90)

    """
    This is where the previous Experiment 656 style of test
    could generate artificial ambiguity.

    The level d is part of the state.

    Signature:

        frame
        d
        p_alive
        n_alive
        current p bit
        anchor p bit
        current n bit
        target n bit

    We check whether that local state uniquely determines the
    threshold bit.
    """

    buckets = defaultdict(set)

    for s in states:

        p_alive = True
        n_alive = True

        max_d = min(MAX_LEVEL, s["depth"] + 2)

        for d in range(1, max_d + 1):

            if d >= 2:
                p_mask = 1 << (d - 2)

                p_bit = int(bool(s["p"] & p_mask))
                p0_bit = int(bool(s["p0"] & p_mask))

                if p_bit != p0_bit:
                    p_alive = False

            else:
                p_bit = 0
                p0_bit = 0

            n_mask = 1 << (d - 1)

            n_bit = int(bool(s["n"] & n_mask))
            target_n_bit = int(bool((-s["c"]) & n_mask))

            if n_bit != target_n_bit:
                n_alive = False

            signature = (
                s["frame"],
                d,
                int(p_alive),
                int(n_alive),
                p_bit,
                p0_bit,
                n_bit,
                target_n_bit,
            )

            buckets[signature].add(p_alive and n_alive)

    ambiguous = [
        (signature, values)
        for signature, values in buckets.items()
        if len(values) > 1
    ]

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    for signature, values in ambiguous[:20]:
        print(f"    {signature} -> {sorted(values)}")

    print()

    return len(ambiguous)


# =============================================================================
# TEST 5
# =============================================================================

def test_first_failure_encoding(states):
    print("=" * 90)
    print("TEST 5: FIRST-FAILURE ENCODING")
    print("=" * 90)

    failures = 0
    distribution = Counter()

    for s in states:

        if s["tp"] >= INF // 2:
            first_failure = None
            code = "no-break"
        else:
            first_failure = s["tp"] + 1
            code = f"break@{first_failure}"

        if first_failure is None:
            predicted = s["w"]
        else:
            predicted = min(first_failure, s["w"])

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch "
                    f"n={s['n']} "
                    f"first_failure={first_failure} "
                    f"w={s['w']} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

        distribution[(s["frame"], code)] += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    print()
    print("    branch-code distribution:")

    for key in sorted(distribution, key=str):
        print(f"        {key}: {distribution[key]}")

    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_threshold_tree(states):
    print("=" * 90)
    print("TEST 6: THRESHOLD SET NESTING")
    print("=" * 90)

    failures = 0
    checks = 0

    for s in states:

        max_d = min(MAX_LEVEL, s["depth"] + 2)

        previous = True

        for d in range(1, max_d + 1):

            current = threshold(s, d)

            checks += 1

            # Once false, it must remain false.
            if previous and not current:
                previous = False

            elif not previous and current:
                failures += 1

                if failures <= 20:
                    print(
                        f"    non-monotone threshold "
                        f"n={s['n']} d={d}"
                    )

    print(f"checked={checks}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_exact_depth_partition(states):
    print("=" * 90)
    print("TEST 7: EXACT TERMINAL EVENT PARTITION")
    print("=" * 90)

    counts = Counter()
    failures = 0

    for s in states:

        tp = s["tp"]
        w = s["w"]

        if tp >= INF // 2:
            branch_first = INF
        else:
            branch_first = tp + 1

        if branch_first < w:
            event = "branch-first"
        elif branch_first > w:
            event = "n-first"
        else:
            event = "simultaneous"

        predicted = min(branch_first, w)

        if predicted != s["depth"]:
            failures += 1

        counts[(s["frame"], event)] += 1

    for key in sorted(counts, key=str):
        print(f"    {key}: {counts[key]}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(states):
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

    by_n = {s["n"]: s for s in states}

    shown = 0

    for n in wanted:

        s = by_n.get(n)

        if s is None:
            continue

        print()
        print(
            f"n={s['n']} "
            f"p={s['p']} "
            f"q={s['q']} "
            f"frame={s['frame']}"
        )

        print(
            f"    p0={s['p0']} "
            f"q0={s['q0']} "
            f"c={s['c']}"
        )

        print(
            f"    tp={s['tp']} "
            f"tq={s['tq']} "
            f"w={s['w']} "
            f"depth={s['depth']}"
        )

        if s["tp"] >= INF // 2:
            first_failure = None
        else:
            first_failure = s["tp"] + 1

        print(
            f"    first_failed_branch_level={first_failure}"
        )

        if first_failure is None:
            predicted = s["w"]
        else:
            predicted = min(first_failure, s["w"])

        print(
            f"    predicted depth="
            f"min(first_failure,w)={predicted}"
        )

        max_d = min(MAX_LEVEL, s["depth"] + 1)

        for d in range(1, max_d + 1):
            print(
                f"    d={d}: "
                f"p-prefix={p_prefix(s,d)} "
                f"q-prefix={q_prefix(s,d)} "
                f"n-prefix={n_prefix(s,d)} "
                f"threshold={threshold(s,d)}"
            )

        shown += 1

    print()
    print(f"examples shown={shown}")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    banner()

    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve(PRIME_LIMIT)
    odd = [p for p in primes if p & 1]

    print(f"odd primes={len(odd)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_threshold_theorem(states)
    total_failures += test_q_redundancy(states)
    total_failures += test_local_bit_automaton(states)
    total_failures += test_level_signature(states)
    total_failures += test_first_failure_encoding(states)
    total_failures += test_threshold_tree(states)
    total_failures += test_exact_depth_partition(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
The corrected level-indexed theorem is:

    depth >= d

        iff

    p == p0 (mod 2^(d-1))

    AND

    n == -c (mod 2^d).

with:

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

The important correction relative to the previous script is
that q-prefix redundancy is tested only in the direction that
is actually implied:

    p-prefix(d) AND n-prefix(d)
        =>
    q-prefix(d).

The converse does not generally provide the extra n bit.

The experiment also represents the process as two binary
survival channels:

    P-channel:
        Does p still match p0?

    N-channel:
        Does n still match -c?

The depth is the first level at which either channel fails.

If the local-state test is exact, the hierarchy has an
explicit finite transition description even though the
underlying first-difference position can grow arbitrarily.
"""
    )

    print()
    print(f"TOTAL FAILURES={total_failures}")

    if total_failures == 0:
        print("STATUS=ALL TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 657 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()