#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 686
RECURSIVE SHIFTED-CORE DESCENT

MOTIVATION
==========================================================================================

Experiments 684-685 established:

    M = n + c
    M = 2^w * R
    R odd

    H = gcd(C, M)

    odd(H) = gcd(C, R)

But:

    odd(C) | R

is false for most states.

Therefore factoring the FIRST odd core R does not generally expose C.

The new idea is:

    n+c
       ↓ strip powers of 2
    R_0
       ↓ apply small odd shifts
    R_1
       ↓ strip powers of 2
    R_2
       ↓ ...
    
Question:

    Does the hidden odd part of C
    eventually appear somewhere in this
    recursive shifted-core tree?

More strongly:

    Can the recursive tree determine
        odd(C)
    or
        v2(C)
    or
        depth
    without using p or q?

We test both shifts:

    x -> odd(x + 3)
    x -> odd(x + 9)

because those are the two canonical offsets
already present in the discovered frames.

IMPORTANT
==========================================================================================

This experiment does NOT assume that the recursive
target remains a semiprime.

It only studies arithmetic information flow.

For every original state we know the true:

    C
    odd(C)
    v2(C)
    H
    depth

and compare them against the recursive odd-core tree.

==========================================================================================
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import gcd, isqrt
from typing import Dict, List, Tuple, Set


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

PRIME_LIMIT = 6000

MAX_DEPTH = 8

SHIFTS = (3, 9)

EXAMPLE_LIMIT = 20


# ==========================================================================================
# BASIC ARITHMETIC
# ==========================================================================================

INF = 10**9


def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    x = abs(x)

    while x and (x & 1) == 0:
        x >>= 1

    return x


def sieve(limit: int) -> List[int]:
    a = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        a[0] = 0

    if limit >= 1:
        a[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if not a[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        a[start : limit + 1 : p] = b"\x00" * count

    return [i for i in range(2, limit + 1) if a[i]]


def build_states():
    primes = [p for p in sieve(PRIME_LIMIT) if p & 1]

    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            if n % 4 == 3:
                frame = "A"
                C = q + 3
                c = 9
            else:
                frame = "B"
                C = p + 1
                c = 3

            M = n + c
            w = v2(M)
            R = M >> w

            H = gcd(C, M)

            states.append(
                {
                    "n": n,
                    "p": p,
                    "q": q,
                    "frame": frame,
                    "C": C,
                    "c": c,
                    "M": M,
                    "w": w,
                    "R0": R,
                    "H": H,
                    "odd_H": odd_part(H),
                    "odd_C": odd_part(C),
                    "v2_C": v2(C),
                    "depth": v2(gcd(2 * C, M)),
                }
            )

    return states


# ==========================================================================================
# RECURSIVE CORE TREE
# ==========================================================================================

def next_core(x: int, shift: int) -> int:
    """
    x is expected to be odd.

    Return odd(x + shift).
    """
    y = x + shift

    if y == 0:
        return 0

    return odd_part(y)


def recursive_core_tree(
    start: int,
    max_depth: int = MAX_DEPTH,
) -> Dict[int, Set[Tuple[int, int]]]:
    """
    Breadth-first recursive tree.

    Return:

        depth -> {(core, shift_used)}

    We retain the shift used to reach each node so that
    later collisions can be analyzed.

    Duplicate cores are allowed across different depths
    but the first occurrence is enough for discovery.
    """

    layers: Dict[int, Set[Tuple[int, int]]] = defaultdict(set)

    current = {start}
    seen = {start}

    for depth in range(max_depth + 1):
        for x in current:
            layers[depth].add((x, 0))

        if depth == max_depth:
            break

        nxt = set()

        for x in current:
            for shift in SHIFTS:
                y = next_core(x, shift)

                if y <= 0:
                    continue

                if y not in seen:
                    seen.add(y)
                    nxt.add(y)

                layers[depth + 1].add((y, shift))

        current = nxt

        if not current:
            break

    return layers


def flatten_tree(
    layers: Dict[int, Set[Tuple[int, int]]]
) -> Dict[int, Tuple[int, int]]:
    """
    Map core -> (first_depth, shift).

    """
    out = {}

    for depth in sorted(layers):
        for core, shift in layers[depth]:
            if core not in out:
                out[core] = (depth, shift)

    return out


# ==========================================================================================
# TEST 0
# ==========================================================================================

def test_baseline(states):
    print("=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)

    failures = 0

    for s in states:
        expected = gcd(s["C"], s["M"])

        if expected != s["H"]:
            failures += 1

            if failures <= EXAMPLE_LIMIT:
                print(
                    "mismatch n={} expected_H={} actual_H={}".format(
                        s["n"],
                        expected,
                        s["H"],
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()


# ==========================================================================================
# TEST 1
# ==========================================================================================

def test_recursive_odd_C_visibility(states):
    print("=" * 90)
    print("TEST 1: DOES RECURSIVE SHIFTED DESCENT REVEAL odd(C)?")
    print("=" * 90)

    found = 0
    missing = 0

    first_hits = Counter()

    examples = []

    for s in states:
        tree = recursive_core_tree(
            s["R0"],
            MAX_DEPTH,
        )

        flat = flatten_tree(tree)

        target = s["odd_C"]

        if target in flat:
            found += 1

            d, shift = flat[target]
            first_hits[d] += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    (
                        s,
                        d,
                        shift,
                        target,
                    )
                )
        else:
            missing += 1

    print("odd(C) found in recursive tree={}".format(found))
    print("odd(C) never found={}".format(missing))
    print()

    print("first-hit depth distribution:")

    for d, count in sorted(first_hits.items()):
        print(
            "    depth={} -> {}".format(
                d,
                count,
            )
        )

    print()

    if examples:
        print("examples:")

        for s, d, shift, target in examples:
            print(
                "    n={} frame={} C={} odd(C)={} R0={}".format(
                    s["n"],
                    s["frame"],
                    s["C"],
                    target,
                    s["R0"],
                )
            )

            print(
                "        hit_depth={} shift={}".format(
                    d,
                    shift,
                )
            )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# TEST 2
# ==========================================================================================

def test_recursive_v2C_reconstruction(states):
    print("=" * 90)
    print("TEST 2: RECURSIVE CORE -> v2(C)")
    print("=" * 90)

    """
    Search for a very simple rule:

        v2(C)

    being equal to:

        v2(core + shift)

    or:

        v2(core - target)

    for some node in the recursive tree.

    This is exploratory and records the earliest node
    that reproduces the true v2(C).
    """

    found = 0
    missing = 0

    hit_depth_distribution = Counter()

    examples = []

    for s in states:
        target_v = s["v2_C"]

        layers = recursive_core_tree(
            s["R0"],
            MAX_DEPTH,
        )

        hit = None

        for depth in sorted(layers):
            for core, _shift in layers[depth]:

                # Candidate local valuations.
                vals = {
                    v2(core + 3),
                    v2(core + 9),
                }

                if target_v in vals:
                    hit = depth
                    break

            if hit is not None:
                break

        if hit is None:
            missing += 1
        else:
            found += 1
            hit_depth_distribution[hit] += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    (
                        s,
                        hit,
                    )
                )

    print("states with matching recursive v2={}".format(found))
    print("states without matching recursive v2={}".format(missing))
    print()

    print("first matching depth:")

    for d, count in sorted(hit_depth_distribution.items()):
        print(
            "    depth={} -> {}".format(
                d,
                count,
            )
        )

    print()

    for s, hit in examples:
        print(
            "    n={} frame={} C={} v2(C)={} R0={}".format(
                s["n"],
                s["frame"],
                s["C"],
                s["v2_C"],
                s["R0"],
            )
        )
        print(
            "        first_match_depth={}".format(
                hit
            )
        )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# TEST 3
# ==========================================================================================

def test_recursive_depth_recovery(states):
    print("=" * 90)
    print("TEST 3: CAN THE RECURSIVE TREE RECOVER DEPTH?")
    print("=" * 90)

    """
    Search for arithmetic signatures that reproduce the actual depth.

    We deliberately use only the recursive odd-core tree and w.

    For every state collect:

        signature = (
            w,
            sorted(
                v2(core+3),
                v2(core+9)
            )
            over all reachable nodes
        )

    Then ask whether identical signatures always produce
    identical depth.

    This is an information-compression test.
    """

    buckets: Dict[Tuple, Set[int]] = defaultdict(set)

    examples_by_sig = {}

    for s in states:
        layers = recursive_core_tree(
            s["R0"],
            MAX_DEPTH,
        )

        vals = set()

        for depth in sorted(layers):
            for core, _shift in layers[depth]:
                vals.add(v2(core + 3))
                vals.add(v2(core + 9))

        signature = (
            s["frame"],
            s["w"],
            tuple(sorted(vals)),
        )

        buckets[signature].add(s["depth"])

        if signature not in examples_by_sig:
            examples_by_sig[signature] = s

    ambiguous = 0
    exact = 0

    examples = []

    for sig, depths in buckets.items():

        if len(depths) > 1:
            ambiguous += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    (
                        sig,
                        depths,
                        examples_by_sig[sig],
                    )
                )
        else:
            exact += 1

    print("signatures={}".format(len(buckets)))
    print("exact signatures={}".format(exact))
    print("ambiguous signatures={}".format(ambiguous))

    print()

    if examples:
        print("first ambiguous signatures:")

        for sig, depths, s in examples:
            print(
                "    signature={} depths={}".format(
                    sig,
                    sorted(depths),
                )
            )

            print(
                "        n={} p={} q={} frame={} R0={} C={} depth={}".format(
                    s["n"],
                    s["p"],
                    s["q"],
                    s["frame"],
                    s["R0"],
                    s["C"],
                    s["depth"],
                )
            )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# TEST 4
# ==========================================================================================

def test_recursive_odd_H_visibility(states):
    print("=" * 90)
    print("TEST 4: DOES RECURSION RECOVER odd(H)?")
    print("=" * 90)

    found = 0
    missing = 0

    first_hits = Counter()

    examples = []

    for s in states:
        target = s["odd_H"]

        layers = recursive_core_tree(
            s["R0"],
            MAX_DEPTH,
        )

        flat = flatten_tree(layers)

        if target in flat:
            found += 1

            d, shift = flat[target]
            first_hits[d] += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    (
                        s,
                        d,
                        shift,
                    )
                )

        else:
            missing += 1

    print("odd(H) found={}".format(found))
    print("odd(H) missing={}".format(missing))
    print()

    print("first-hit depth distribution:")

    for d, count in sorted(first_hits.items()):
        print(
            "    depth={} -> {}".format(
                d,
                count,
            )
        )

    print()

    for s, d, shift in examples:
        print(
            "    n={} frame={} R0={} odd(H)={} H={}".format(
                s["n"],
                s["frame"],
                s["R0"],
                s["odd_H"],
                s["H"],
            )
        )

        print(
            "        hit_depth={} shift={}".format(
                d,
                shift,
            )
        )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# TEST 5
# ==========================================================================================

def test_core_shrink(states):
    print("=" * 90)
    print("TEST 5: CORE SIZE REDUCTION")
    print("=" * 90)

    """
    Measure how aggressively recursive shifted reduction
    shrinks the target.

    For each state we follow every discovered core and record:

        min(core) / R0
        max core depth
        number of unique nodes
    """

    ratio_buckets = Counter()

    max_nodes = 0
    max_depth_seen = 0

    examples = []

    for s in states:
        layers = recursive_core_tree(
            s["R0"],
            MAX_DEPTH,
        )

        flat = flatten_tree(layers)

        cores = list(flat.keys())

        if not cores:
            continue

        min_core = min(cores)

        ratio_num = min_core
        ratio_den = max(1, s["R0"])

        # Logarithmic-ish classification without floating point.
        if min_core == 1:
            bucket = "1"
        elif min_core * 2 <= s["R0"]:
            bucket = "<=R/2"
        elif min_core * 4 <= s["R0"] * 3:
            bucket = "<=3R/4"
        else:
            bucket = ">3R/4"

        ratio_buckets[bucket] += 1

        max_nodes = max(
            max_nodes,
            len(flat),
        )

        if len(layers) - 1 > max_depth_seen:
            max_depth_seen = len(layers) - 1

        if len(examples) < EXAMPLE_LIMIT:
            examples.append(
                (
                    s,
                    min_core,
                    len(flat),
                )
            )

    print("minimum-core ratio buckets:")

    for bucket, count in sorted(ratio_buckets.items()):
        print(
            "    {} -> {}".format(
                bucket,
                count,
            )
        )

    print()
    print("maximum unique recursive nodes={}".format(max_nodes))
    print("maximum recursive depth={}".format(max_depth_seen))

    print()
    print("examples:")

    for s, minimum, count in examples:
        print(
            "    n={} R0={} min_core={} nodes={}".format(
                s["n"],
                s["R0"],
                minimum,
                count,
            )
        )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# TEST 6
# ==========================================================================================

def test_simple_recursive_signature(states):
    print("=" * 90)
    print("TEST 6: SIMPLE RECURSIVE SIGNATURE")
    print("=" * 90)

    """
    Try a deliberately small signature:

        (frame, w, v2(R0+3), v2(R0+9))

    This asks whether the first recursive split already
    contains enough information to determine depth.
    """

    buckets: Dict[Tuple, Set[int]] = defaultdict(set)

    for s in states:
        sig = (
            s["frame"],
            s["w"],
            v2(s["R0"] + 3),
            v2(s["R0"] + 9),
        )

        buckets[sig].add(s["depth"])

    ambiguous = 0
    examples = []

    for sig, depths in buckets.items():
        if len(depths) > 1:
            ambiguous += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    (
                        sig,
                        depths,
                    )
                )

    print("signatures={}".format(len(buckets)))
    print("ambiguous={}".format(ambiguous))

    print()

    for sig, depths in examples:
        print(
            "    signature={} depths={}".format(
                sig,
                sorted(depths),
            )
        )

    print()
    print("checked={}".format(len(states)))
    print()


# ==========================================================================================
# EXAMPLES
# ==========================================================================================

def print_examples(states):
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = {
        9,
        15,
        21,
        33,
        39,
        51,
        57,
        69,
        77,
        87,
        93,
        111,
        141,
        183,
        485879,
        5579767,
    }

    count = 0

    for s in states:
        if s["n"] not in wanted:
            continue

        count += 1

        print()
        print(
            "n={} p={} q={} frame={}".format(
                s["n"],
                s["p"],
                s["q"],
                s["frame"],
            )
        )

        print("    C={}".format(s["C"]))
        print("    c={}".format(s["c"]))
        print("    M={}".format(s["M"]))
        print("    w={}".format(s["w"]))
        print("    R0={}".format(s["R0"]))
        print("    H={}".format(s["H"]))
        print("    odd(H)={}".format(s["odd_H"]))
        print("    odd(C)={}".format(s["odd_C"]))
        print("    v2(C)={}".format(s["v2_C"]))
        print("    depth={}".format(s["depth"]))

        layers = recursive_core_tree(
            s["R0"],
            min(MAX_DEPTH, 4),
        )

        print("    recursive cores:")

        for d in sorted(layers):
            cores = sorted(
                core
                for core, _shift in layers[d]
            )

            if len(cores) > 12:
                cores = cores[:12] + ["..."]

            print(
                "        level {}: {}".format(
                    d,
                    cores,
                )
            )

        if count >= EXAMPLE_LIMIT:
            break


# ==========================================================================================
# MAIN
# ==========================================================================================

def main():
    print("=" * 90)
    print("EXPERIMENT 686 START")
    print("=" * 90)
    print()

    states = build_states()

    odd_prime_count = len(
        [p for p in sieve(PRIME_LIMIT) if p & 1]
    )

    print("prime limit={}".format(PRIME_LIMIT))
    print("odd primes={}".format(odd_prime_count))
    print("semiprimes={}".format(len(states)))
    print()

    test_baseline(states)
    test_recursive_odd_C_visibility(states)
    test_recursive_v2C_reconstruction(states)
    test_recursive_depth_recovery(states)
    test_recursive_odd_H_visibility(states)
    test_core_shrink(states)
    test_simple_recursive_signature(states)
    print_examples(states)

    print()
    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Established from the previous experiments:")
    print()
    print("    M = n+c")
    print("    M = 2^w R")
    print("    R odd")
    print()
    print("    H = gcd(C,M)")
    print("    odd(H) = gcd(C,R)")
    print()

    print("Experiment 684 showed:")
    print()
    print("    odd(C) | R")
    print()
    print("is false for most states.")
    print()

    print("Experiment 685 showed:")
    print()
    print("    R + w does not uniquely determine H")
    print()
    print("from divisor information alone.")
    print()
    print("The new experiment tests whether this information")
    print("loss can be repaired by recursively transforming")
    print("the odd core:")
    print()
    print("    R -> odd(R+3)")
    print("    R -> odd(R+9)")
    print()
    print("The target questions are:")
    print()
    print("    1. Does recursive descent reveal odd(C)?")
    print("    2. Does recursive descent reveal v2(C)?")
    print("    3. Does recursive descent determine depth?")
    print("    4. Does recursive descent recover odd(H)?")
    print()
    print("A positive result would give evidence for a genuine")
    print("recursive 'factor n+c, strip 2, retarget' mechanism.")
    print()
    print("A negative result would establish that the first")
    print("shifted-core transformation loses information that")
    print("cannot be recovered by the same small-shift recursion.")
    print()

    print("=" * 90)
    print("EXPERIMENT 686 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
