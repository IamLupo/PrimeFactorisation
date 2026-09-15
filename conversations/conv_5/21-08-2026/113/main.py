#!/usr/bin/env python3
"""
==========================================================================================
EXPERIMENT 675
==========================================================================================

2-ADIC DISCRIMINANT ROOT BRANCH -> DEPTH

Goal
----
Experiment 674 suggested that the remaining ambiguity after
observing (n,S) is a 2-adic square-root branch of

    Delta = S^2 - 4n
    D^2 = Delta

with

    S = p+q
    D = p-q.

This experiment avoids O(N^2) state separation.

We test:

    1. D is always a root of Delta modulo 2^k.
    2. Enumerated 2-adic roots reproduce the actual truncated depth.
    3. (n,S,D) modulo 2^k determines truncated depth.
    4. The actual D branch is sufficient to select the depth.
    5. Ambiguous (n,S) buckets are explained by different D-root branches.
    6. The number of possible root branches and depth branches is compared.

Frames:

    FRAME A:
        n == 3 (mod 4)
        C = q+3
        c = 9
        2C = S-D+6

    FRAME B:
        n == 1 (mod 4)
        C = p+1
        c = 3
        2C = S+D+2

Exact depth law:

    depth = v2(gcd(2C, n+c))
          = min(v2(C)+1, v2(n+c)).

All arithmetic is integer/exact.
No O(N^2) comparisons are used.
==========================================================================================
"""

from collections import defaultdict, Counter
from math import isqrt
from typing import Dict, List, Tuple, Iterable


PRIME_LIMIT = 6000

# Root moduli tested.
ROOT_BITS = [4, 6, 8, 10, 12]

# Examples to print.
EXAMPLE_LIMIT = 16

INF_V2 = 10**9


# ==========================================================================================
# BASIC NUMBER THEORY
# ==========================================================================================

def v2(x: int) -> int:
    """Return v2(|x|), with a large sentinel for x=0."""
    x = abs(x)
    if x == 0:
        return INF_V2
    return (x & -x).bit_length() - 1


def primes_upto(limit: int) -> List[int]:
    """Standard sieve."""
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    p = 2
    while p * p <= limit:
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )
        p += 1

    return [i for i in range(2, limit + 1) if sieve[i]]


def odd_primes_upto(limit: int) -> List[int]:
    return [p for p in primes_upto(limit) if p & 1]


# ==========================================================================================
# STATE GENERATION
# ==========================================================================================

def classify_frame(n: int) -> str:
    r = n & 3
    if r == 3:
        return "A"
    if r == 1:
        return "B"
    raise ValueError("Odd semiprime must be 1 or 3 modulo 4.")


def frame_parameters(frame: str) -> Tuple[int, int, int]:
    """
    Return (p0, c, frame_constant).

    frame_constant is used for the shifted numerator:
        A: 2C = S-D+6
        B: 2C = S+D+2
    """
    if frame == "A":
        return 3, 9, 6
    if frame == "B":
        return -1, 3, 2
    raise ValueError(frame)


def make_state(p: int, q: int) -> dict:
    n = p * q
    frame = classify_frame(n)

    if frame == "A":
        A = p - 3
        B = q + 3
        C = B
        c = 9
        two_c_numerator = (p + q) - (p - q) + 6
    else:
        A = p + 1
        B = q - 3
        C = A
        c = 3
        two_c_numerator = (p + q) + (p - q) + 2

    S = p + q
    D = p - q
    delta = S * S - 4 * n

    assert delta == D * D
    assert two_c_numerator == 2 * C

    depth = v2(__import__("math").gcd(2 * C, n + c))
    direct_depth = min(v2(C) + 1, v2(n + c))

    assert depth == direct_depth

    return {
        "p": p,
        "q": q,
        "n": n,
        "frame": frame,
        "S": S,
        "D": D,
        "delta": delta,
        "A": A,
        "B": B,
        "C": C,
        "c": c,
        "depth": depth,
    }


def generate_states(limit: int) -> List[dict]:
    ps = odd_primes_upto(limit)

    states: List[dict] = []

    for i, p in enumerate(ps):
        for q in ps[i:]:
            states.append(make_state(p, q))

    return states


# ==========================================================================================
# 2-ADIC ROOT ENUMERATION
# ==========================================================================================

def roots_square_mod_power_of_two(delta: int, bits: int) -> List[int]:
    """
    Enumerate all x modulo 2^bits satisfying

        x^2 == delta (mod 2^bits)

    by binary Hensel-style lifting.

    Starting from roots mod 1:

        x = 0

    then repeatedly lift r -> r or r + 2^e.

    The root counts remain small for the data studied here.
    """
    if bits <= 0:
        return [0]

    # Mod 2 initially.
    roots = [0]
    current_bits = 1

    delta &= (1 << bits) - 1

    # Add the root 1 if it exists mod 2.
    if (delta & 1) == 1:
        roots = [1]

    while current_bits < bits:
        modulus = 1 << current_bits
        next_modulus = modulus << 1

        next_roots = []

        for r in roots:
            r0 = r
            r1 = r + modulus

            if ((r0 * r0 - delta) % next_modulus) == 0:
                next_roots.append(r0)

            if ((r1 * r1 - delta) % next_modulus) == 0:
                next_roots.append(r1)

        roots = sorted(set(next_roots))
        current_bits += 1

    return roots


def residue_v2(x: int, bits: int) -> int:
    """
    v2 of an integer residue modulo 2^bits, capped at bits.

    For example:
        residue 0 mod 2^bits -> bits
    """
    x &= (1 << bits) - 1

    if x == 0:
        return bits

    return (x & -x).bit_length() - 1


def c_from_root_residue(
    frame: str,
    S_res: int,
    D_res: int,
    bits: int,
) -> int:
    """
    Recover C modulo 2^(bits-1) from S,D modulo 2^bits.

    FRAME A:
        2C = S-D+6

    FRAME B:
        2C = S+D+2

    Since S and D are even, the numerator is divisible by 2.
    """
    modulus = 1 << bits
    c_modulus = 1 << (bits - 1)

    S_res &= modulus - 1
    D_res &= modulus - 1

    if frame == "A":
        numerator = (S_res - D_res + 6) & (modulus - 1)
    else:
        numerator = (S_res + D_res + 2) & (modulus - 1)

    assert numerator % 2 == 0

    return (numerator // 2) % c_modulus


def depth_from_root_residue(
    frame: str,
    n_value: int,
    S_res: int,
    D_res: int,
    bits: int,
) -> int:
    """
    Compute the depth visible through a D-root modulo 2^bits.

    C is known modulo 2^(bits-1), hence v2(C) is known up to
    bits-1.

    The n-side valuation is computed exactly here and then
    capped by bits-1 because that is the available C precision.
    """
    c_modulus_bits = bits - 1
    C_res = c_from_root_residue(frame, S_res, D_res, bits)

    vc = residue_v2(C_res, c_modulus_bits)

    c = 9 if frame == "A" else 3
    vn = min(v2(n_value + c), c_modulus_bits)

    return min(vc + 1, vn)


# ==========================================================================================
# TEST 0
# ==========================================================================================

def test_baseline(states: List[dict]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
    print("=" * 90)

    failures = 0

    for st in states:
        depth1 = st["depth"]

        depth2 = min(
            v2(st["C"]) + 1,
            v2(st["n"] + st["c"]),
        )

        if depth1 != depth2:
            failures += 1
            if failures <= 10:
                print(
                    "mismatch",
                    st["n"],
                    st["frame"],
                    depth1,
                    depth2,
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# ==========================================================================================
# TEST 1
# ==========================================================================================

def test_discriminant_root_identity(states: List[dict]) -> int:
    print("=" * 90)
    print("TEST 1: EXACT DISCRIMINANT ROOT IDENTITY")
    print("=" * 90)

    failures = 0

    for st in states:
        lhs = st["D"] * st["D"]
        rhs = st["S"] * st["S"] - 4 * st["n"]

        if lhs != rhs:
            failures += 1
            if failures <= 10:
                print(
                    "mismatch n={} p={} q={} D={} S={} lhs={} rhs={}".format(
                        st["n"],
                        st["p"],
                        st["q"],
                        st["D"],
                        st["S"],
                        lhs,
                        rhs,
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# ==========================================================================================
# TEST 2
# ==========================================================================================

def test_root_image_and_depth(states: List[dict], bits: int) -> Tuple[int, dict]:
    """
    For each exact state:

      delta = S^2 - 4n

    enumerate all D roots mod 2^bits.

    Verify:
      actual D is among roots

    and:
      actual truncated depth is among root-generated depths.

    Returns failures and a diagnostic record.
    """
    failures = 0

    root_count_counter = Counter()
    root_depth_count_counter = Counter()

    examples = []

    cap = bits - 1

    for st in states:
        delta = st["delta"]

        roots = roots_square_mod_power_of_two(delta, bits)
        root_set = set(roots)

        actual_D = st["D"] & ((1 << bits) - 1)

        if actual_D not in root_set:
            failures += 1
            if failures <= 10:
                print(
                    "missing actual root:",
                    st["n"],
                    st["frame"],
                    "Dmod={}".format(actual_D),
                )
            continue

        generated_depths = set()

        S_res = st["S"] & ((1 << bits) - 1)

        for root in roots:
            d = depth_from_root_residue(
                st["frame"],
                st["n"],
                S_res,
                root,
                bits,
            )
            generated_depths.add(min(d, cap))

        actual_depth = min(st["depth"], cap)

        if actual_depth not in generated_depths:
            failures += 1

            if len(examples) < EXAMPLE_LIMIT:
                examples.append(
                    {
                        "n": st["n"],
                        "p": st["p"],
                        "q": st["q"],
                        "frame": st["frame"],
                        "actual_depth": actual_depth,
                        "roots": roots,
                        "generated_depths": sorted(generated_depths),
                    }
                )

        root_count_counter[len(roots)] += 1
        root_depth_count_counter[len(generated_depths)] += 1

    print("k={}".format(bits))
    print("root-count distribution={}".format(sorted(root_count_counter.items())))
    print(
        "root-depth-count distribution={}".format(
            sorted(root_depth_count_counter.items())
        )
    )
    print("failures={}".format(failures))

    if examples:
        print("counterexamples:")
        for ex in examples:
            print(
                "    n={} p={} q={} frame={}".format(
                    ex["n"],
                    ex["p"],
                    ex["q"],
                    ex["frame"],
                )
            )
            print("        actual_depth={}".format(ex["actual_depth"]))
            print("        roots={}".format(ex["roots"]))
            print(
                "        generated_depths={}".format(
                    ex["generated_depths"]
                )
            )

    print()

    diagnostics = {
        "root_counts": root_count_counter,
        "root_depth_counts": root_depth_count_counter,
        "examples": examples,
    }

    return failures, diagnostics


# ==========================================================================================
# TEST 3
# ==========================================================================================

def test_root_branch_signature(states: List[dict], bits: int) -> int:
    """
    Test the key hypothesis:

        (frame, n mod 2^bits, S mod 2^bits, D mod 2^bits)

    determines truncated depth.

    This is a direct bucket test, not pairwise comparison.
    """
    print("=" * 90)
    print(
        "TEST 3: ROOT BRANCH SIGNATURE -> TRUNCATED DEPTH k={}".format(
            bits
        )
    )
    print("=" * 90)

    modulus = 1 << bits
    cap = bits - 1

    buckets: Dict[Tuple[str, int, int, int], set] = defaultdict(set)

    for st in states:
        sig = (
            st["frame"],
            st["n"] % modulus,
            st["S"] % modulus,
            st["D"] % modulus,
        )

        buckets[sig].add(min(st["depth"], cap))

    ambiguous = {
        sig: depths
        for sig, depths in buckets.items()
        if len(depths) > 1
    }

    print("signatures={}".format(len(buckets)))
    print("ambiguous={}".format(len(ambiguous)))

    if ambiguous:
        print("first ambiguous signatures:")
        shown = 0

        for sig, depths in ambiguous.items():
            print("    {} -> {}".format(sig, sorted(depths)))
            shown += 1
            if shown >= 10:
                break

    failures = len(ambiguous)

    print("failures={}".format(failures))
    print()

    return failures


# ==========================================================================================
# TEST 4
# ==========================================================================================

def test_n_s_ambiguity_explained_by_roots(
    states: List[dict],
    bits: int,
) -> int:
    """
    For each (frame, n, S) residue signature, compare:

        actual depth set

    to the union of depths produced by all D roots.

    The point is not that n,S should determine depth.
    The point is that its ambiguity should be fully explained
    by the D-root branches.
    """
    print("=" * 90)
    print(
        "TEST 4: (n,S)-AMBIGUITY VS 2-ADIC ROOT BRANCHES k={}".format(
            bits
        )
    )
    print("=" * 90)

    modulus = 1 << bits
    cap = bits - 1

    actual_buckets: Dict[Tuple[str, int, int], set] = defaultdict(set)
    representative: Dict[Tuple[str, int, int], dict] = {}

    for st in states:
        sig = (
            st["frame"],
            st["n"] % modulus,
            st["S"] % modulus,
        )

        actual_buckets[sig].add(min(st["depth"], cap))
        representative.setdefault(sig, st)

    mismatches = 0
    ambiguous = 0
    resolved_by_roots = 0
    examples = []

    for sig, actual_depths in actual_buckets.items():
        st = representative[sig]

        if len(actual_depths) <= 1:
            continue

        ambiguous += 1

        roots = roots_square_mod_power_of_two(
            st["delta"],
            bits,
        )

        S_res = st["S"] % modulus

        root_depths = set()

        for root in roots:
            root_depths.add(
                depth_from_root_residue(
                    st["frame"],
                    st["n"],
                    S_res,
                    root,
                    bits,
                )
            )

        root_depths = {min(x, cap) for x in root_depths}

        if actual_depths.issubset(root_depths):
            resolved_by_roots += 1
        else:
            mismatches += 1

            if len(examples) < 10:
                examples.append(
                    (
                        sig,
                        sorted(actual_depths),
                        roots,
                        sorted(root_depths),
                    )
                )

    print("ambiguous n,S signatures={}".format(ambiguous))
    print("explained by root branches={}".format(resolved_by_roots))
    print("root-image mismatches={}".format(mismatches))

    if examples:
        print("first root-image mismatches:")
        for sig, actual_depths, roots, root_depths in examples:
            print("    signature={}".format(sig))
            print("        actual_depths={}".format(actual_depths))
            print("        roots={}".format(roots))
            print("        root_depths={}".format(root_depths))

    failures = mismatches

    print("failures={}".format(failures))
    print()

    return failures


# ==========================================================================================
# TEST 5
# ==========================================================================================

def show_selected_examples(
    states: List[dict],
    bits: int,
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

    cap = bits - 1
    modulus = 1 << bits

    shown = 0

    for n in wanted:
        st = by_n.get(n)

        if st is None:
            continue

        roots = roots_square_mod_power_of_two(
            st["delta"],
            bits,
        )

        root_depth_pairs = []

        S_res = st["S"] % modulus

        for root in roots:
            d = depth_from_root_residue(
                st["frame"],
                st["n"],
                S_res,
                root,
                bits,
            )

            root_depth_pairs.append(
                (root, min(d, cap))
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
        print("    Delta={}".format(st["delta"]))
        print("    v2(C)={}".format(v2(st["C"])))
        print("    v2(n+c)={}".format(v2(st["n"] + st["c"])))
        print("    actual depth={}".format(st["depth"]))
        print("    D mod 2^{}={}".format(bits, st["D"] % modulus))
        print("    roots mod 2^{}={}".format(bits, roots))
        print("    root -> truncated-depth={}".format(root_depth_pairs))
        print()

        shown += 1
        if shown >= EXAMPLE_LIMIT:
            break


# ==========================================================================================
# MAIN
# ==========================================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 675 START")
    print("=" * 90)
    print()
    print("prime limit={}".format(PRIME_LIMIT))

    odd_primes = odd_primes_upto(PRIME_LIMIT)

    print("odd primes={}".format(len(odd_primes)))

    states = generate_states(PRIME_LIMIT)

    print("semiprimes={}".format(len(states)))
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_discriminant_root_identity(states)

    root_diagnostics = {}

    # Main root tests.
    for bits in ROOT_BITS:
        failures, diagnostics = test_root_image_and_depth(
            states,
            bits,
        )
        total_failures += failures
        root_diagnostics[bits] = diagnostics

        total_failures += test_root_branch_signature(
            states,
            bits,
        )

        total_failures += test_n_s_ambiguity_explained_by_roots(
            states,
            bits,
        )

    show_selected_examples(states, max(ROOT_BITS))

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()
    print("The tested discriminant identity is:")
    print()
    print("    Delta = S^2 - 4n")
    print("    D^2   = Delta")
    print()
    print("The factor pair is recovered from:")
    print()
    print("    p = (S+D)/2")
    print("    q = (S-D)/2")
    print()
    print("Canonical residual:")
    print()
    print("    FRAME A:")
    print("        C = q+3")
    print("        2C = S-D+6")
    print()
    print("    FRAME B:")
    print("        C = p+1")
    print("        2C = S+D+2")
    print()
    print("Canonical depth law:")
    print()
    print("    depth = v2(gcd(2C,n+c))")
    print()
    print("The experiment asks whether the ambiguity remaining")
    print("after observing (n,S) is exactly the choice of the")
    print("2-adic square-root branch D of Delta.")
    print()
    print("The strongest passing outcome is:")
    print()
    print("    (n,S,D) mod 2^k")
    print("        -> truncated depth")
    print()
    print("and:")
    print()
    print("    actual (n,S)-depth ambiguity")
    print("        is contained in")
    print("    the depth image of the 2-adic D-root branches.")
    print()
    print("TOTAL FAILURES={}".format(total_failures))

    if total_failures == 0:
        print("STATUS=ALL ROOT-BRANCH TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 675 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
