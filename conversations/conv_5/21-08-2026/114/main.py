#!/usr/bin/env python3
"""
==========================================================================================
EXPERIMENT 676
==========================================================================================

2-ADIC DISCRIMINANT ROOT DEPTH CLASSIFICATION

Motivation
----------
Experiment 675 established:

    Delta = S^2 - 4n
    D^2 = Delta

and showed that the remaining (n,S) ambiguity is completely
explained by the different 2-adic roots D of Delta.

The important new observation is visible in the examples:

    n=15:
        roots split into depth 2 and depth 3

    n=39:
        roots split into depth 2 and depth 4

    n=93:
        roots split into depth 3 and depth 5

    n=485879:
        roots split into depth 2 and depth 8

So the next question is:

    CAN THE DEPTH OF A ROOT D BE EXPRESSED DIRECTLY
    FROM A SIMPLE 2-ADIC PROPERTY OF D?

Candidates:

    1. v2(D - D0)
    2. v2(D + D0)
    3. v2(D - S)
    4. v2(D + S)
    5. v2(D^2 - Delta)
    6. the first differing bit between D and the
       actual/nearby discriminant branch
    7. a simple residue class of D modulo 2^k

We do not compare all states pairwise.

All tests use hash buckets and direct root enumeration.

Important:
----------
A root D modulo 2^k is not necessarily the actual integer D.
The experiment therefore distinguishes:

    actual D

from

    alternative 2-adic roots.

The target is a structural rule for the root branch itself.

==========================================================================================
"""

from collections import defaultdict, Counter
from math import gcd
from typing import Dict, List, Tuple


PRIME_LIMIT = 6000
ROOT_BITS = [6, 8, 10, 12]
MAX_EXAMPLES = 20
INF = 10**9


# ==========================================================================================
# BASIC NUMBER THEORY
# ==========================================================================================

def v2(x: int) -> int:
    x = abs(x)
    if x == 0:
        return INF
    return (x & -x).bit_length() - 1


def sieve(limit: int) -> List[int]:
    a = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        a[0] = 0
    if limit >= 1:
        a[1] = 0

    p = 2
    while p * p <= limit:
        if a[p]:
            start = p * p
            a[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )
        p += 1

    return [i for i in range(limit + 1) if a[i]]


def odd_primes(limit: int) -> List[int]:
    return [p for p in sieve(limit) if p & 1]


# ==========================================================================================
# STATE CONSTRUCTION
# ==========================================================================================

def make_state(p: int, q: int) -> dict:
    n = p * q
    S = p + q
    D = p - q
    Delta = D * D

    if n % 4 == 3:
        frame = "A"
        c = 9
        p0 = 3
        q0 = -3
        A = p - 3
        B = q + 3
        C = B
    else:
        frame = "B"
        c = 3
        p0 = -1
        q0 = 3
        A = p + 1
        B = q - 3
        C = A

    depth = v2(gcd(2 * C, n + c))

    return {
        "p": p,
        "q": q,
        "n": n,
        "S": S,
        "D": D,
        "Delta": Delta,
        "frame": frame,
        "c": c,
        "p0": p0,
        "q0": q0,
        "A": A,
        "B": B,
        "C": C,
        "depth": depth,
    }


def generate_states(limit: int) -> List[dict]:
    ps = odd_primes(limit)
    states = []

    for i, p in enumerate(ps):
        for q in ps[i:]:
            states.append(make_state(p, q))

    return states


# ==========================================================================================
# 2-ADIC ROOTS
# ==========================================================================================

def square_roots_mod_2k(delta: int, k: int) -> List[int]:
    """
    Enumerate x modulo 2^k satisfying:

        x^2 == delta mod 2^k

    by binary lifting.
    """
    if k <= 0:
        return [0]

    # Root modulo 1.
    roots = [0]

    # Lift one bit at a time.
    for bits in range(1, k + 1):
        mod = 1 << bits
        prev_mod = mod >> 1

        if bits == 1:
            candidates = [0, 1]
        else:
            candidates = []
            for r in roots:
                candidates.append(r)
                candidates.append(r + prev_mod)

        roots = sorted({
            x % mod
            for x in candidates
            if (x * x - delta) % mod == 0
        })

    return roots


def trunc_v2_residue(x: int, k: int) -> int:
    """
    v2 of x modulo 2^k, capped at k.
    """
    x &= (1 << k) - 1

    if x == 0:
        return k

    return (x & -x).bit_length() - 1


# ==========================================================================================
# FRAME DEPTH FROM A ROOT
# ==========================================================================================

def C_from_root(frame: str, S: int, D_res: int, k: int) -> int:
    """
    Recover C modulo 2^(k-1).

        A:
            2C = S - D + 6

        B:
            2C = S + D + 2
    """
    mod = 1 << k
    half_mod = 1 << (k - 1)

    S &= mod - 1
    D_res &= mod - 1

    if frame == "A":
        numerator = (S - D_res + 6) & (mod - 1)
    else:
        numerator = (S + D_res + 2) & (mod - 1)

    return (numerator // 2) % half_mod


def root_depth(frame: str, n: int, S: int, D_res: int, k: int) -> int:
    """
    Depth visible from the root residue.
    """
    vc = trunc_v2_residue(
        C_from_root(frame, S, D_res, k),
        k - 1,
    )

    c = 9 if frame == "A" else 3
    vn = min(v2(n + c), k - 1)

    return min(vc + 1, vn)


# ==========================================================================================
# TEST 0
# ==========================================================================================

def test_baseline(states: List[dict]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)

    failures = 0

    for st in states:
        predicted = min(
            v2(st["C"]) + 1,
            v2(st["n"] + st["c"]),
        )

        if predicted != st["depth"]:
            failures += 1
            if failures <= 10:
                print(
                    "mismatch n={} p={} q={} frame={} actual={} predicted={}".format(
                        st["n"],
                        st["p"],
                        st["q"],
                        st["frame"],
                        st["depth"],
                        predicted,
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# ==========================================================================================
# TEST 1
# ==========================================================================================

def test_root_depth_partition(
    states: List[dict],
    k: int,
) -> int:
    """
    Build root-depth classes for each (n,S).

    Question:

        Does each root naturally map to one of a small number
        of depth classes?

    Also compare actual root D with all alternatives.
    """
    print("=" * 90)
    print("TEST 1: ROOT -> DEPTH PARTITION k={}".format(k))
    print("=" * 90)

    mod = 1 << k

    # (frame, n mod 2^k, S mod 2^k)
    buckets = defaultdict(lambda: defaultdict(list))

    for st in states:
        sig = (
            st["frame"],
            st["n"] % mod,
            st["S"] % mod,
        )

        roots = square_roots_mod_2k(st["Delta"], k)

        for r in roots:
            d = root_depth(
                st["frame"],
                st["n"],
                st["S"],
                r,
                k,
            )

            buckets[sig][d].append(r)

    ambiguous_root_classes = 0
    total_root_classes = 0
    max_classes = 0
    distribution = Counter()

    examples = []

    for sig, depth_map in buckets.items():
        classes = len(depth_map)

        total_root_classes += 1
        distribution[classes] += 1
        max_classes = max(max_classes, classes)

        if classes > 1:
            ambiguous_root_classes += 1

            if len(examples) < MAX_EXAMPLES:
                examples.append(
                    (
                        sig,
                        {
                            d: sorted(rs)[:12]
                            for d, rs in depth_map.items()
                        },
                    )
                )

    print("signatures={}".format(total_root_classes))
    print("root-depth-class distribution={}".format(
        sorted(distribution.items())
    ))
    print("max depth classes per (n,S)={}".format(max_classes))
    print("multi-depth root signatures={}".format(
        ambiguous_root_classes
    ))

    if examples:
        print("examples:")
        for sig, depth_map in examples:
            print("    {}".format(sig))
            for d in sorted(depth_map):
                print("        depth={} roots={}".format(
                    d,
                    depth_map[d],
                ))

    print()

    # This test is diagnostic, not an expected-failure theorem.
    return 0


# ==========================================================================================
# TEST 2
# ==========================================================================================

def test_root_valuation_signatures(
    states: List[dict],
    k: int,
) -> int:
    """
    Test simple D-based signatures against the root-derived depth.

    Signatures tested:

        v2(D)
        v2(D-S)
        v2(D+S)
        v2(D^2-S^2)
        v2(D-D_actual)
        v2(D+D_actual)

    The last two are mainly diagnostic and are evaluated only
    within each actual state.

    A signature is successful if one state bucket maps to one
    root-derived depth.
    """
    print("=" * 90)
    print("TEST 2: ROOT VALUATION SIGNATURES k={}".format(k))
    print("=" * 90)

    mod = 1 << k
    signatures = {
        "v2(D)": defaultdict(set),
        "v2(D-S)": defaultdict(set),
        "v2(D+S)": defaultdict(set),
        "v2(D^2-S^2)": defaultdict(set),
    }

    counts = Counter()

    for st in states:
        roots = square_roots_mod_2k(st["Delta"], k)

        for r in roots:
            d = root_depth(
                st["frame"],
                st["n"],
                st["S"],
                r,
                k,
            )

            values = {
                "v2(D)": trunc_v2_residue(r, k),
                "v2(D-S)": trunc_v2_residue(r - st["S"], k),
                "v2(D+S)": trunc_v2_residue(r + st["S"], k),
                "v2(D^2-S^2)": trunc_v2_residue(
                    (r * r) - (st["S"] * st["S"]),
                    k,
                ),
            }

            for name, value in values.items():
                signatures[name][
                    (
                        st["frame"],
                        st["n"] % mod,
                        st["S"] % mod,
                        value,
                    )
                ].add(d)

    failures = 0

    for name, buckets in signatures.items():
        ambiguous = sum(
            1
            for depths in buckets.values()
            if len(depths) > 1
        )

        counts[name] = ambiguous

        print(
            "{:<18} signatures={:8d} ambiguous={:8d}".format(
                name,
                len(buckets),
                ambiguous,
            )
        )

    print()

    return failures


# ==========================================================================================
# TEST 3
# ==========================================================================================

def test_first_differing_branch_bit(
    states: List[dict],
    k: int,
) -> int:
    """
    Compare the root residue against the actual branch.

    For every exact state, classify alternative roots by:

        v2(r - D_actual)

    and ask whether the root-derived depth depends only on this
    branch displacement.

    This is the most direct continuation of Experiments 649-675.
    """
    print("=" * 90)
    print("TEST 3: ROOT BRANCH DISPLACEMENT k={}".format(k))
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(set)
    examples = []

    for st in states:
        roots = square_roots_mod_2k(st["Delta"], k)
        d_actual = st["D"] % mod

        for r in roots:
            displacement = trunc_v2_residue(
                r - d_actual,
                k,
            )

            depth = root_depth(
                st["frame"],
                st["n"],
                st["S"],
                r,
                k,
            )

            sig = (
                st["frame"],
                st["n"] % mod,
                st["S"] % mod,
                displacement,
            )

            buckets[sig].add(depth)

    ambiguous = {
        sig: depths
        for sig, depths in buckets.items()
        if len(depths) > 1
    }

    print("signatures={}".format(len(buckets)))
    print("ambiguous={}".format(len(ambiguous)))

    if ambiguous:
        print("first ambiguous signatures:")

        for sig, depths in list(ambiguous.items())[:10]:
            print("    {} -> {}".format(
                sig,
                sorted(depths),
            ))

    print()

    return len(ambiguous)


# ==========================================================================================
# TEST 4
# ==========================================================================================

def test_root_pairing_structure(
    states: List[dict],
    k: int,
) -> int:
    """
    Study the canonical paired-root structure.

    For a root r, compare:

        r
        -r
        r + 2^(k-1)
        -r + 2^(k-1)

    where valid.

    We record which depths occur in each root orbit.

    This is intentionally descriptive: it tells us whether the
    depth is invariant under the standard 2-adic square-root
    symmetries or selects one branch.
    """
    print("=" * 90)
    print("TEST 4: 2-ADIC ROOT ORBIT STRUCTURE k={}".format(k))
    print("=" * 90)

    mod = 1 << k

    orbit_distribution = Counter()
    examples = []

    for st in states[:]:
        roots = square_roots_mod_2k(st["Delta"], k)
        root_set = set(roots)

        seen = set()

        for r in roots:
            if r in seen:
                continue

            candidates = {
                r % mod,
                (-r) % mod,
                (r + (mod // 2)) % mod,
                (-r + (mod // 2)) % mod,
            }

            orbit = sorted(candidates & root_set)

            seen.update(orbit)

            depths = tuple(sorted({
                root_depth(
                    st["frame"],
                    st["n"],
                    st["S"],
                    x,
                    k,
                )
                for x in orbit
            }))

            orbit_distribution[(len(orbit), depths)] += 1

            if (
                len(depths) > 1
                and len(examples) < MAX_EXAMPLES
            ):
                examples.append(
                    (
                        st["n"],
                        st["frame"],
                        orbit,
                        depths,
                    )
                )

    print("orbit distribution:")
    for key, count in sorted(orbit_distribution.items()):
        print("    {}: {}".format(key, count))

    if examples:
        print()
        print("mixed-depth root orbits:")
        for n, frame, orbit, depths in examples:
            print(
                "    n={} frame={} orbit={} depths={}".format(
                    n,
                    frame,
                    orbit,
                    depths,
                )
            )

    print()

    return 0


# ==========================================================================================
# TEST 5
# ==========================================================================================

def show_examples(
    states: List[dict],
    k: int,
) -> None:
    print("=" * 90)
    print("TEST 5: EXAMPLES")
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

    by_n = {st["n"]: st for st in states}
    mod = 1 << k

    for n in wanted:
        st = by_n.get(n)

        if st is None:
            continue

        roots = square_roots_mod_2k(
            st["Delta"],
            k,
        )

        print(
            "n={} p={} q={} frame={}".format(
                st["n"],
                st["p"],
                st["q"],
                st["frame"],
            )
        )
        print("    S={}".format(st["S"]))
        print("    D(actual)={}".format(st["D"]))
        print("    Delta={}".format(st["Delta"]))
        print("    C={}".format(st["C"]))
        print("    v2(C)={}".format(v2(st["C"])))
        print("    v2(n+c)={}".format(v2(st["n"] + st["c"])))
        print("    depth={}".format(st["depth"]))
        print("    roots mod 2^{} count={}".format(
            k,
            len(roots),
        ))

        rows = []

        for r in roots:
            rows.append(
                (
                    r,
                    trunc_v2_residue(r, k),
                    trunc_v2_residue(r - st["S"], k),
                    trunc_v2_residue(r + st["S"], k),
                    root_depth(
                        st["frame"],
                        st["n"],
                        st["S"],
                        r,
                        k,
                    ),
                )
            )

        print(
            "    root rows: "
            "(D, v2(D), v2(D-S), v2(D+S), root-depth)"
        )

        # Avoid gigantic example output.
        for row in rows[:16]:
            print("        {}".format(row))

        if len(rows) > 16:
            print(
                "        ... {} more roots".format(
                    len(rows) - 16
                )
            )

        print()


# ==========================================================================================
# MAIN
# ==========================================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 676 START")
    print("=" * 90)
    print()
    print("prime limit={}".format(PRIME_LIMIT))

    ps = odd_primes(PRIME_LIMIT)
    states = generate_states(PRIME_LIMIT)

    print("odd primes={}".format(len(ps)))
    print("semiprimes={}".format(len(states)))
    print()

    total_failures = 0

    total_failures += test_baseline(states)

    for k in ROOT_BITS:
        total_failures += test_root_depth_partition(
            states,
            k,
        )

        total_failures += test_root_valuation_signatures(
            states,
            k,
        )

        # The displacement experiment is the most important
        # one, but only run at higher precision.
        if k >= 8:
            total_failures += test_first_differing_branch_bit(
                states,
                k,
            )

        if k >= 10:
            total_failures += test_root_pairing_structure(
                states,
                k,
            )

    show_examples(states, max(ROOT_BITS))

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()
    print("Experiment 675 established:")
    print()
    print("    Delta = S^2 - 4n")
    print("    D^2   = Delta")
    print()
    print("and demonstrated that the ambiguity in (n,S)")
    print("is completely explained by the available 2-adic")
    print("D-root branches.")
    print()
    print("Experiment 676 now asks:")
    print()
    print("    WHICH PROPERTY OF THE ROOT D")
    print("    selects the depth?")
    print()
    print("The principal candidates are:")
    print()
    print("    v2(D)")
    print("    v2(D-S)")
    print("    v2(D+S)")
    print("    v2(D^2-S^2)")
    print("    v2(D-D_actual)")
    print()
    print("The important structural distinction is:")
    print()
    print("    (n,S)")
    print("        -> several 2-adic D branches")
    print()
    print("    (n,S,D)")
    print("        -> exact branch")
    print("        -> exact depth")
    print()
    print("If a small valuation of D or a simple branch")
    print("displacement also determines the depth, then")
    print("the discriminant-root description can be compressed")
    print("again.")
    print()
    print("TOTAL FAILURES={}".format(total_failures))

    if total_failures == 0:
        print("STATUS=ROOT STRUCTURE TESTS COMPLETED")
    else:
        print("STATUS=BRANCH SIGNATURE COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 676 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
