#!/usr/bin/env python3
"""
==========================================================================================
EXPERIMENT 689
==========================================================================================

AFFINE ANCHOR / SHIFTED CONSTANT FAMILY SEARCH

Known canonical cases:

    FRAME A:
        C = q + 3
        c = 9

    FRAME B:
        C = p + 1
        c = 3

New hypothesis:

    C = factor + a
    c = a^2

We search:

    Frame A:
        C = q + a

    Frame B:
        C = p + a

and ask whether some fixed c produces the same exact depth law:

    depth = v2(gcd(2*C, n+c))

The experiment separately checks:

    1. depth equality
    2. H equality
    3. whether C itself equals the canonical C
    4. whether simple relations between a and c emerge

All arithmetic is exact integer arithmetic.
"""

from __future__ import annotations

from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 250

A_MIN = -15
A_MAX = 15

C_MIN = -200
C_MAX = 200

MAX_EXAMPLES = 20


# =============================================================================
# BASIC UTILITIES
# =============================================================================

def primes_upto(limit: int) -> list[int]:
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [p for p in range(3, limit + 1, 2) if sieve[p]]


def v2(x: int) -> int:
    if x == 0:
        return 10**9

    x = abs(x)
    return (x & -x).bit_length() - 1


def frame_for(p: int, q: int, n: int) -> str:
    """
    Existing experimentally established frame selection.
    """
    return "A" if n % 4 == 3 else "B"


def canonical_C_c(p: int, q: int, n: int) -> tuple[int, int]:
    frame = frame_for(p, q, n)

    if frame == "A":
        return q + 3, 9

    return p + 1, 3


def canonical_H(p: int, q: int, n: int) -> int:
    C, c = canonical_C_c(p, q, n)
    return gcd_abs(C, n + c)


def canonical_depth(p: int, q: int, n: int) -> int:
    C, c = canonical_C_c(p, q, n)
    return candidate_depth(C, n, c)


def gcd_abs(a: int, b: int) -> int:
    """
    Local integer gcd implementation so the experiment remains explicit.
    """
    a = abs(a)
    b = abs(b)

    while b:
        a, b = b, a % b

    return a


def candidate_C(frame: str, p: int, q: int, a: int) -> int:
    if frame == "A":
        return q + a

    return p + a


def candidate_H(C: int, n: int, c: int) -> int:
    M = n + c

    if C == 0 or M == 0:
        return 0

    return gcd_abs(C, M)


def candidate_depth(C: int, n: int, c: int) -> int:
    M = n + c

    if C == 0 or M == 0:
        return 10**9

    return v2(gcd_abs(2 * C, M))


# =============================================================================
# GENERATE DATASET
# =============================================================================

primes = primes_upto(PRIME_LIMIT)

states = []

for i, p in enumerate(primes):
    for q in primes[i:]:
        n = p * q

        frame = frame_for(p, q, n)

        C0, c0 = canonical_C_c(p, q, n)
        H0 = canonical_H(p, q, n)
        d0 = canonical_depth(p, q, n)

        states.append(
            {
                "p": p,
                "q": q,
                "n": n,
                "frame": frame,
                "C0": C0,
                "c0": c0,
                "H0": H0,
                "depth": d0,
            }
        )


print("=" * 90)
print("EXPERIMENT 689 START")
print("=" * 90)
print()
print(f"prime limit={PRIME_LIMIT}")
print(f"odd primes={len(primes)}")
print(f"semiprimes={len(states)}")
print()


# =============================================================================
# TEST 0: BASELINE
# =============================================================================

print("=" * 90)
print("TEST 0: BASELINE")
print("=" * 90)

baseline_failures = 0

for s in states:
    C = s["C0"]
    c = s["c0"]

    H = candidate_H(C, s["n"], c)
    d = candidate_depth(C, s["n"], c)

    if H != s["H0"] or d != s["depth"]:
        baseline_failures += 1

print(f"checked={len(states)}")
print(f"failures={baseline_failures}")
print()


# =============================================================================
# TEST 1: c = a^2
# =============================================================================

print("=" * 90)
print("TEST 1: SPECIAL FAMILY c = a^2")
print("=" * 90)

square_results = []

for a in range(A_MIN, A_MAX + 1):

    if a == 0:
        continue

    c = a * a

    for frame in ("A", "B"):

        frame_states = [s for s in states if s["frame"] == frame]

        failures = 0
        H_matches = 0
        depth_matches = 0
        first_mismatch = None

        for s in frame_states:

            C = candidate_C(frame, s["p"], s["q"], a)

            if C == 0 or s["n"] + c == 0:
                failures += 1

                if first_mismatch is None:
                    first_mismatch = s

                continue

            H = candidate_H(C, s["n"], c)
            d = candidate_depth(C, s["n"], c)

            if H == s["H0"]:
                H_matches += 1

            if d == s["depth"]:
                depth_matches += 1
            else:
                failures += 1

                if first_mismatch is None:
                    first_mismatch = s

        count = len(frame_states)

        square_results.append(
            {
                "frame": frame,
                "a": a,
                "c": c,
                "count": count,
                "failures": failures,
                "H_matches": H_matches,
                "depth_matches": depth_matches,
                "first_mismatch": first_mismatch,
            }
        )

        status = "EXACT" if failures == 0 else "NO"

        print(
            f"frame={frame} "
            f"a={a:3d} "
            f"c=a^2={c:5d} "
            f"depth={depth_matches:6d}/{count:<6d} "
            f"H={H_matches:6d}/{count:<6d} "
            f"status={status}"
        )

print()


# =============================================================================
# TEST 2: SEARCH ALL c FOR EACH a
# =============================================================================

print("=" * 90)
print("TEST 2: FIXED AFFINE ANCHOR a -> SEARCH ALL c")
print("=" * 90)

exact_pairs = []
best_pairs = []

for frame in ("A", "B"):

    frame_states = [s for s in states if s["frame"] == frame]
    frame_count = len(frame_states)

    print()
    print(f"FRAME {frame}")
    print("-" * 90)

    for a in range(A_MIN, A_MAX + 1):

        if a == 0:
            continue

        best_score = -1
        best_c = None
        best_h_matches = 0

        exact_for_a = []

        for c in range(C_MIN, C_MAX + 1):

            depth_matches = 0
            H_matches = 0

            for s in frame_states:

                C = candidate_C(frame, s["p"], s["q"], a)

                if C == 0 or s["n"] + c == 0:
                    continue

                H = candidate_H(C, s["n"], c)
                d = candidate_depth(C, s["n"], c)

                if d == s["depth"]:
                    depth_matches += 1

                if H == s["H0"]:
                    H_matches += 1

            if depth_matches > best_score:
                best_score = depth_matches
                best_c = c
                best_h_matches = H_matches

            if depth_matches == frame_count:
                exact_for_a.append((c, H_matches))

        if exact_for_a:

            for c, h_matches in exact_for_a:
                exact_pairs.append(
                    {
                        "frame": frame,
                        "a": a,
                        "c": c,
                        "H_matches": h_matches,
                    }
                )

            print(
                f"a={a:3d} "
                f"EXACT c={exact_for_a} "
                f"square_c={a * a}"
            )

        else:

            best_pairs.append(
                {
                    "frame": frame,
                    "a": a,
                    "c": best_c,
                    "depth_matches": best_score,
                    "frame_count": frame_count,
                    "H_matches": best_h_matches,
                }
            )

            print(
                f"a={a:3d} "
                f"best_c={best_c:5d} "
                f"depth={best_score:6d}/{frame_count:<6d} "
                f"H={best_h_matches:6d}/{frame_count:<6d} "
                f"square_c={a * a}"
            )

print()


# =============================================================================
# TEST 3: RELATIONSHIPS BETWEEN a AND c
# =============================================================================

print("=" * 90)
print("TEST 3: RELATIONSHIPS BETWEEN a AND c")
print("=" * 90)

if not exact_pairs:
    print("NO exact anchor/shift pairs found.")
else:
    for item in exact_pairs:

        a = item["a"]
        c = item["c"]

        relations = []

        if c == a * a:
            relations.append("c=a^2")

        if c == a:
            relations.append("c=a")

        if c == 2 * a:
            relations.append("c=2a")

        if c == a * a + 1:
            relations.append("c=a^2+1")

        if c == a * a - 1:
            relations.append("c=a^2-1")

        if c == a * (a + 1):
            relations.append("c=a(a+1)")

        if c == a * (a - 1):
            relations.append("c=a(a-1)")

        relation_text = ", ".join(relations)
        if not relation_text:
            relation_text = "no-simple-relation"

        print(
            f"frame={item['frame']} "
            f"a={a:3d} "
            f"c={c:5d} "
            f"H-matches={item['H_matches']:6d} "
            f"relation={relation_text}"
        )

print()


# =============================================================================
# TEST 4: DOES THE CANDIDATE C EQUAL THE CANONICAL C?
# =============================================================================

print("=" * 90)
print("TEST 4: CANDIDATE C VS CANONICAL C")
print("=" * 90)

for item in exact_pairs:

    frame = item["frame"]
    a = item["a"]
    c = item["c"]

    frame_states = [s for s in states if s["frame"] == frame]

    C_equal = 0

    for s in frame_states:

        C = candidate_C(frame, s["p"], s["q"], a)

        if C == s["C0"]:
            C_equal += 1

    print(
        f"frame={frame} "
        f"a={a:3d} "
        f"c={c:5d} "
        f"C_exact={C_equal}/{len(frame_states)}"
    )

print()


# =============================================================================
# TEST 5: SEARCH FOR THE KNOWN ANCHORS
# =============================================================================

print("=" * 90)
print("TEST 5: KNOWN CANONICAL ANCHORS CONTROL")
print("=" * 90)

known_cases = [
    ("A", 3, 9),
    ("B", 1, 3),
]

for frame, a, c in known_cases:

    frame_states = [s for s in states if s["frame"] == frame]

    H_matches = 0
    depth_matches = 0
    C_matches = 0

    for s in frame_states:

        C = candidate_C(frame, s["p"], s["q"], a)
        H = candidate_H(C, s["n"], c)
        d = candidate_depth(C, s["n"], c)

        if C == s["C0"]:
            C_matches += 1

        if H == s["H0"]:
            H_matches += 1

        if d == s["depth"]:
            depth_matches += 1

    print(
        f"frame={frame} "
        f"a={a} "
        f"c={c} "
        f"C={C_matches}/{len(frame_states)} "
        f"H={H_matches}/{len(frame_states)} "
        f"depth={depth_matches}/{len(frame_states)}"
    )

print()


# =============================================================================
# TEST 6: BEST NON-EXACT PAIRS
# =============================================================================

print("=" * 90)
print("TEST 6: BEST NON-EXACT ANCHOR/SHIFT PAIRS")
print("=" * 90)

best_pairs_sorted = sorted(
    best_pairs,
    key=lambda x: (
        x["depth_matches"] / x["frame_count"]
        if x["frame_count"]
        else 0.0
    ),
    reverse=True,
)

for item in best_pairs_sorted[:20]:

    ratio = (
        item["depth_matches"] / item["frame_count"]
        if item["frame_count"]
        else 0.0
    )

    print(
        f"frame={item['frame']} "
        f"a={item['a']:3d} "
        f"best_c={item['c']:5d} "
        f"depth={item['depth_matches']:6d}/{item['frame_count']:<6d} "
        f"ratio={ratio:.6f} "
        f"H={item['H_matches']:6d}"
    )

print()


# =============================================================================
# TEST 7: REPRESENTATIVE SQUARE-LAW FAILURES
# =============================================================================

print("=" * 90)
print("TEST 7: REPRESENTATIVE c=a^2 COUNTEREXAMPLES")
print("=" * 90)

shown = 0

for item in square_results:

    if shown >= MAX_EXAMPLES:
        break

    a = item["a"]

    if a not in (1, 3, 5, 7, 9):
        continue

    mismatch = item["first_mismatch"]

    if mismatch is None:
        continue

    frame = item["frame"]
    c = item["c"]

    C = candidate_C(frame, mismatch["p"], mismatch["q"], a)
    H = candidate_H(C, mismatch["n"], c)
    d = candidate_depth(C, mismatch["n"], c)

    print(
        f"frame={frame} a={a} c={c}\n"
        f"    n={mismatch['n']} "
        f"p={mismatch['p']} "
        f"q={mismatch['q']}\n"
        f"    candidate C={C}\n"
        f"    candidate H={H}\n"
        f"    candidate depth={d}\n"
        f"    true C={mismatch['C0']}\n"
        f"    true H={mismatch['H0']}\n"
        f"    true depth={mismatch['depth']}\n"
    )

    shown += 1


# =============================================================================
# TEST 8: CHECK WHETHER c=a^2 IS AT LEAST CONSISTENT WITH ANCHOR STRUCTURE
# =============================================================================

print("=" * 90)
print("TEST 8: ANCHOR/SQUARE CONSISTENCY")
print("=" * 90)

for a in (1, 3, 5, 7, 9):

    print(
        f"a={a:3d} "
        f"hypothesized c={a*a:5d}"
    )

print()
print(
    "The key result is whether any a != known anchors "
    "survives TEST 2 exactly."
)
print()


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("=" * 90)
print("FINAL STRUCTURAL SUMMARY")
print("=" * 90)

print(
    f"""
KNOWN EXACT CASES:

    Frame A:
        C = q + 3
        c = 9

    Frame B:
        C = p + 1
        c = 3

NEW HYPOTHESIS:

    C = factor + a
    c = a^2

SEARCH RANGE:

    a = [{A_MIN}, {A_MAX}]
    c = [{C_MIN}, {C_MAX}]

TESTS:

    depth:
        v2(gcd(2C, n+c))

    H:
        gcd(C, n+c)

The important distinction is:

    depth-exact
        does NOT necessarily mean
    C-exact.

Therefore the experiment separately checks:

    1. depth equality
    2. H equality
    3. C equality
    4. algebraic relationship between a and c

RESULT:

    baseline failures = {baseline_failures}
    exact anchor/shift pairs = {len(exact_pairs)}

"""

)

if exact_pairs:
    print("EXACT PAIRS FOUND:")
    for item in exact_pairs:
        print(
            f"    frame={item['frame']} "
            f"a={item['a']} "
            f"c={item['c']} "
            f"H-matches={item['H_matches']}"
        )
else:
    print("NO NEW EXACT ANCHOR/SHIFT PAIRS FOUND.")

print()
print("=" * 90)
print("EXPERIMENT 689 FINISHED")
print("=" * 90)