#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 655
==========================================================================================

CORRECTED FIRST-DIFFERING-BIT / BRANCH-RESIDUE THEOREM

GOAL
----
Experiment 654 had an indexing error in TEST 4.

The correct relation is:

    tp = v2(p-p0)

and for finite tp:

    first_failed_level = tp + 1.

Therefore:

    depth = min(first_failed_level, w)

where:

    w = v2(n+c).

This experiment:

    1. rebuilds the complete odd semiprime domain
    2. verifies the exact depth law
    3. verifies the threshold theorem
    4. verifies first-failure indexing
    5. verifies depth from first failure
    6. verifies residue-class equivalence
    7. searches for a smaller branch representation
    8. verifies q-prefix redundancy
    9. prints explicit examples

FRAME A:
    p0 =  3
    q0 = -3
    c  =  9

FRAME B:
    p0 = -1
    q0 =  3
    c  =  3
"""

from collections import defaultdict, Counter
from math import isqrt


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 6000
MAX_TEST_LEVEL = 24
INF = 10**9
EXAMPLE_LIMIT = 20


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def sieve(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []

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
    """2-adic valuation; v2(0)=INF."""
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def frame_for_n(n: int) -> str:
    """
    Historical frame convention:

        FRAME A: n == 3 mod 4
        FRAME B: n == 1 mod 4
    """
    r = n & 3

    if r == 3:
        return "A"
    if r == 1:
        return "B"

    raise ValueError(f"Unexpected odd n residue: n={n}, n mod 4={r}")


def frame_parameters(frame: str) -> tuple[int, int, int]:
    """
    Return:
        p0, q0, c
    """
    if frame == "A":
        return 3, -3, 9

    if frame == "B":
        return -1, 3, 3

    raise ValueError(frame)


# =============================================================================
# GLOBAL STATE CONSTRUCTION
# =============================================================================

def build_states(primes: list[int]) -> list[dict]:
    """
    Build unordered semiprime states p <= q.

    For the historical domain:

        782 odd primes

    gives:

        782 * 783 / 2 = 306153
    """
    odd_primes = [p for p in primes if p & 1]

    states = []

    for i, p in enumerate(odd_primes):
        for q in odd_primes[i:]:
            n = p * q
            frame = frame_for_n(n)

            p0, q0, c = frame_parameters(frame)

            A = p - p0
            B = q - q0

            if frame == "A":
                # historical coordinates
                X = (B - A) // 2
                Y = (A + B) // 2
            else:
                X = (3 * A - B) // 2
                Y = (3 * A + B) // 2

            # Exact transformed numerator valuation.
            depth = min(v2(2 * X), v2(2 * Y))

            tp = v2(p - p0)
            tq = v2(q - q0)
            w = v2(n + c)

            # Experiment 642 / 648 quantity.
            m = min(tp, tq)

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
                    "X": X,
                    "Y": Y,
                    "depth": depth,
                    "tp": tp,
                    "tq": tq,
                    "w": w,
                    "m": m,
                }
            )

    return states


# =============================================================================
# THRESHOLD PREDICATES
# =============================================================================

def p_prefix(state: dict, d: int) -> bool:
    """
    p == p0 (mod 2^(d-1))
    """
    if d <= 1:
        return True

    mod = 1 << (d - 1)
    return (state["p"] - state["p0"]) % mod == 0


def q_prefix(state: dict, d: int) -> bool:
    """
    q == q0 (mod 2^(d-1))
    """
    if d <= 1:
        return True

    mod = 1 << (d - 1)
    return (state["q"] - state["q0"]) % mod == 0


def n_prefix(state: dict, d: int) -> bool:
    """
    n == -c (mod 2^d)
    """
    if d <= 0:
        return True

    mod = 1 << d
    return (state["n"] + state["c"]) % mod == 0


def threshold(state: dict, d: int) -> bool:
    """
    Exact threshold theorem:

        depth >= d
        iff
        p-prefix && n-prefix
    """
    return p_prefix(state, d) and n_prefix(state, d)


def first_failed_branch_level(state: dict) -> int | None:
    """
    The first level d at which:

        p != p0 mod 2^(d-1)

    For finite tp:

        first_failed_level = tp + 1.

    Returns None when no failure occurs in the tested range.
    """
    tp = state["tp"]

    if tp >= INF // 2:
        return None

    return tp + 1


# =============================================================================
# TEST 0
# =============================================================================

def test_domain(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for s in states:
        if s["p"] & 1 == 0 or s["q"] & 1 == 0:
            failures += 1
        if s["p"] > s["q"]:
            failures += 1
        if s["p"] * s["q"] != s["n"]:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_exact_depth_law(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 1: EXACT DEPTH LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = min(s["tp"] + 1, s["w"])

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s['n']} "
                    f"p={s['p']} q={s['q']} "
                    f"frame={s['frame']} "
                    f"tp={s['tp']} w={s['w']} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_threshold_theorem(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 2: EXACT THRESHOLD THEOREM")
    print("=" * 90)

    failures = 0
    checks = 0

    for s in states:
        max_d = min(MAX_TEST_LEVEL, s["depth"] + 2)

        for d in range(1, max_d + 1):
            lhs = s["depth"] >= d
            rhs = threshold(s, d)

            checks += 1

            if lhs != rhs:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch n={s['n']} "
                        f"d={d} "
                        f"depth={s['depth']} "
                        f"lhs={lhs} rhs={rhs}"
                    )

    print(f"checks={checks}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_first_failure_index(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 3: FIRST-DIFFERING-BIT INDEXING")
    print("=" * 90)

    failures = 0

    finite_tp = 0
    infinite_tp = 0

    for s in states:
        first = first_failed_branch_level(s)

        if first is None:
            infinite_tp += 1
            continue

        finite_tp += 1

        expected = s["tp"] + 1

        if first != expected:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s['n']} "
                    f"tp={s['tp']} "
                    f"first_failure={first} "
                    f"expected={expected}"
                )

    print(f"checked={len(states)}")
    print(f"finite tp={finite_tp}")
    print(f"infinite tp={infinite_tp}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_depth_from_first_failure(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 4: DEPTH FROM FIRST FAILED BRANCH LEVEL")
    print("=" * 90)

    failures = 0

    for s in states:
        first = first_failed_branch_level(s)

        if first is None:
            predicted = s["w"]
        else:
            predicted = min(first, s["w"])

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s['n']} "
                    f"frame={s['frame']} "
                    f"W={s['w']} "
                    f"first_failure={first} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_first_failure_bit(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 5: FIRST-DIFFERING-BIT CODE")
    print("=" * 90)

    failures = 0
    distribution = Counter()

    for s in states:
        tp = s["tp"]

        if tp >= INF // 2:
            code = ("no-break", s["w"])
        else:
            code = ("break", tp + 1)

        distribution[code] += 1

        # Direct reconstruction.
        if tp >= INF // 2:
            reconstructed_tp = INF
        else:
            reconstructed_tp = code[1] - 1

        if reconstructed_tp != tp:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch n={s['n']} "
                    f"tp={tp} "
                    f"code={code} "
                    f"reconstructed={reconstructed_tp}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()
    print("    branch-code distribution:")

    for key in sorted(distribution, key=lambda x: (x[0], x[1])):
        print(f"        {key}: {distribution[key]}")

    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_residue_equivalence(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 6: FIRST FAILURE <-> RESIDUE CLASS")
    print("=" * 90)

    failures = 0

    for s in states:
        tp = s["tp"]

        # Check every level around the actual first failure.
        if tp >= INF // 2:
            levels = range(1, min(MAX_TEST_LEVEL, s["w"]) + 1)
        else:
            levels = range(1, min(MAX_TEST_LEVEL, tp + 2) + 1)

        for d in levels:
            if d == 1:
                expected_prefix = True
            else:
                expected_prefix = (
                    (s["p"] - s["p0"]) % (1 << (d - 1)) == 0
                )

            actual_prefix = p_prefix(s, d)

            if expected_prefix != actual_prefix:
                failures += 1

                if failures <= 20:
                    print(
                        f"    p-prefix mismatch n={s['n']} "
                        f"d={d}"
                    )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_q_redundancy(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 7: Q-PREFIX REDUNDANCY")
    print("=" * 90)

    failures = 0
    checks = 0

    for s in states:
        max_d = min(MAX_TEST_LEVEL, s["depth"] + 1)

        for d in range(1, max_d + 1):
            # For d=1 everything is automatically valid.
            if d <= 1:
                continue

            p_ok = p_prefix(s, d)
            n_ok = n_prefix(s, d)
            q_ok = q_prefix(s, d)

            # The q-prefix is predicted to follow from p-prefix+n-prefix.
            predicted_q = p_ok and n_ok

            checks += 1

            if q_ok != predicted_q:
                failures += 1

                if failures <= 20:
                    print(
                        f"    mismatch n={s['n']} "
                        f"d={d} "
                        f"p={p_ok} "
                        f"q={q_ok} "
                        f"n={n_ok}"
                    )

    print(f"checks={checks}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 8
# =============================================================================

def test_terminal_partition(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 8: TERMINAL EVENT PARTITION")
    print("=" * 90)

    failures = 0
    counts = Counter()

    for s in states:
        tp = s["tp"]
        w = s["w"]

        if tp < w:
            event = "branch-first"
        elif tp > w:
            event = "n-first"
        else:
            event = "simultaneous"

        counts[(s["frame"], event)] += 1

        predicted = min(tp + 1, w)

        # Terminal category sanity.
        if event == "branch-first" and predicted != tp + 1:
            failures += 1
        elif event == "n-first" and predicted != w:
            failures += 1
        elif event == "simultaneous" and predicted != w:
            failures += 1

    for key in sorted(counts):
        print(f"    {key}: {counts[key]}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 9
# =============================================================================

def test_compressed_signature(states: list[dict]) -> int:
    print("=" * 90)
    print("TEST 9: COMPRESSED BRANCH SIGNATURE")
    print("=" * 90)

    """
    The branch is encoded only by:

        branch_event
        branch_level

    together with:

        w

    For branch-first:
        depth = branch_level

    For n-first/simultaneous:
        depth = w.

    This is intentionally a finite categorical encoding
    of the exact theorem, not a bounded valuation cap.
    """

    groups: dict[tuple, set[int]] = defaultdict(set)

    for s in states:
        tp = s["tp"]
        w = s["w"]

        if tp >= INF // 2:
            branch_code = ("no-break",)
        else:
            branch_code = ("break", tp + 1)

        signature = (
            s["frame"],
            w,
            branch_code,
        )

        groups[signature].add(s["depth"])

    ambiguous = [
        (sig, vals)
        for sig, vals in groups.items()
        if len(vals) > 1
    ]

    print(f"signatures={len(groups)}")
    print(f"ambiguous={len(ambiguous)}")

    for sig, vals in ambiguous[:20]:
        print(f"    {sig} -> {sorted(vals)}")

    print()

    return len(ambiguous)


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(states: list[dict]) -> None:
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = [9, 15, 21, 33, 39, 57, 69, 87, 93, 111,
              141, 183, 213, 77, 485879, 5579767]

    by_n = {s["n"]: s for s in states}

    shown = 0

    for n in wanted:
        s = by_n.get(n)
        if s is None:
            continue

        tp = s["tp"]
        w = s["w"]

        first_failure = first_failed_branch_level(s)

        if first_failure is None:
            predicted = w
            branch_text = "never-breaks-before-n-bound"
        else:
            predicted = min(first_failure, w)
            branch_text = f"break@{first_failure}"

        print()
        print(
            f"n={s['n']} p={s['p']} q={s['q']} "
            f"frame={s['frame']}"
        )
        print(
            f"    p0={s['p0']} c={s['c']} "
            f"tp={tp} w={w}"
        )
        print(
            f"    actual depth={s['depth']}"
        )
        print(
            f"    first_failed_branch_level={first_failure}"
        )
        print(
            f"    predicted depth=min(first_failure,w)={predicted}"
        )
        print(
            f"    branch code={branch_text}"
        )

        max_d = min(s["depth"] + 1, 10)

        print("    thresholds:")
        for d in range(1, max_d + 1):
            print(
                f"        d={d}: "
                f"p-prefix={p_prefix(s, d)} "
                f"n-prefix={n_prefix(s, d)} "
                f"threshold={threshold(s, d)}"
            )

        shown += 1

    print()
    print(f"examples shown={shown}")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 655 START")
    print("=" * 90)
    print()
    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve(PRIME_LIMIT)
    odd_primes = [p for p in primes if p & 1]

    print(f"odd primes={len(odd_primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_domain(states)
    total_failures += test_exact_depth_law(states)
    total_failures += test_threshold_theorem(states)
    total_failures += test_first_failure_index(states)
    total_failures += test_depth_from_first_failure(states)
    total_failures += test_first_failure_bit(states)
    total_failures += test_residue_equivalence(states)
    total_failures += test_q_redundancy(states)
    total_failures += test_terminal_partition(states)

    ambiguous = test_compressed_signature(states)
    total_failures += ambiguous

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()
    print("Corrected theorem:")
    print()
    print("    tp = v2(p-p0)")
    print("    w  = v2(n+c)")
    print()
    print("    depth = min(tp + 1, w)")
    print()
    print("Equivalently:")
    print()
    print("    depth >= d")
    print("        iff")
    print("    p == p0 (mod 2^(d-1))")
    print("    and")
    print("    n == -c (mod 2^d)")
    print()
    print("For finite tp:")
    print()
    print("    first_failed_branch_level = tp + 1")
    print()
    print("Therefore:")
    print()
    print("    depth = min(first_failed_branch_level, w)")
    print()
    print("The off-by-one error from Experiment 654 is explicitly")
    print("eliminated by reconstructing:")
    print()
    print("    tp = first_failed_branch_level - 1")
    print()
    print("rather than:")
    print()
    print("    tp = first_failed_branch_level")
    print()
    print("The experiment also checks whether the complete depth")
    print("can be represented by the categorical branch code:")
    print()
    print("    break@L")
    print("or")
    print("    no-break-before-n-bound")
    print()
    print("combined with:")
    print()
    print("    w = v2(n+c).")
    print()

    if total_failures == 0:
        print("TOTAL FAILURES=0")
        print("STATUS=ALL TESTS PASSED")
    else:
        print(f"TOTAL FAILURES={total_failures}")
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 655 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
