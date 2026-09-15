#!/usr/bin/env python3

from collections import defaultdict, Counter
from math import gcd


# =============================================================================
# EXPERIMENT 677
#
# DISCRIMINANT ROOT DISTANCE THEOREM
#
# Candidate theorem:
#
#     FRAME A:
#         T = S + 6
#         c = 9
#
#     FRAME B:
#         T = S + 2
#         c = 3
#
#     depth = min(v2(D - T), v2(n + c))
#
# Equivalently:
#
#     depth >= d
#
#         iff
#
#     D == T (mod 2^d)
#     and
#     n == -c (mod 2^d)
#
# This experiment checks:
#
#   1. Exact integer identity.
#   2. Exact threshold theorem.
#   3. Root residue -> depth.
#   4. Root distance classes.
#   5. Whether roots of the same distance class always
#      have the same depth.
#   6. Whether the root branches are exactly organized
#      by distance from T.
#
# Root generation is deduplicated at every lift level.
# =============================================================================


PRIME_LIMIT = 6000
ROOT_BITS = [4, 6, 8, 10, 12]
INF = 10**9
MAX_EXAMPLES = 12


# =============================================================================
# BASIC
# =============================================================================

def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


def sieve(limit: int):
    a = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        a[0] = 0

    if limit >= 1:
        a[1] = 0

    p = 2

    while p * p <= limit:
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start:limit + 1:p] = b"\x00" * count

        p += 1

    return [i for i in range(limit + 1) if a[i]]


def odd_primes(limit: int):
    return [p for p in sieve(limit) if p & 1]


# =============================================================================
# STATE
# =============================================================================

def make_state(p: int, q: int):
    n = p * q
    S = p + q
    D = p - q

    if n % 4 == 3:
        frame = "A"
        c = 9

        p0 = 3
        q0 = -3

        A = p - 3
        B = q + 3

        C = B
        T = S + 6

    else:
        frame = "B"
        c = 3

        p0 = -1
        q0 = 3

        A = p + 1
        B = q - 3

        C = A
        T = S + 2

    depth = v2(gcd(2 * C, n + c))

    return {
        "p": p,
        "q": q,
        "n": n,
        "S": S,
        "D": D,
        "Delta": D * D,
        "frame": frame,
        "c": c,
        "p0": p0,
        "q0": q0,
        "A": A,
        "B": B,
        "C": C,
        "T": T,
        "depth": depth,
    }


def generate_states(limit: int):
    ps = odd_primes(limit)

    states = []

    for i, p in enumerate(ps):
        for q in ps[i:]:
            states.append(make_state(p, q))

    return ps, states


# =============================================================================
# ROOT GENERATION
# =============================================================================

def square_roots_mod_2k(delta: int, k: int):
    """
    Deduplicated binary Hensel-style lifting.

    Returns all roots r in [0, 2^k) satisfying

        r^2 == delta mod 2^k
    """

    roots = {0}

    for bits in range(1, k + 1):
        half = 1 << (bits - 1)
        mod = 1 << bits

        nxt = set()

        for r in roots:
            a = r
            b = r + half

            if (a * a - delta) % mod == 0:
                nxt.add(a)

            if (b * b - delta) % mod == 0:
                nxt.add(b)

        roots = nxt

    return sorted(roots)


def trunc_v2(x: int, k: int):
    """
    v2 modulo 2^k, capped at k.
    """
    x &= (1 << k) - 1

    if x == 0:
        return k

    return (x & -x).bit_length() - 1


# =============================================================================
# THEORETICAL ROOT DEPTH
# =============================================================================

def predicted_root_depth(frame, n, S, D, k):
    if frame == "A":
        T = S + 6
        c = 9
    else:
        T = S + 2
        c = 3

    root_distance = trunc_v2(D - T, k)
    n_bound = min(v2(n + c), k)

    return min(root_distance, n_bound)


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states):
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
                    "mismatch n={} p={} q={} frame={} "
                    "actual={} predicted={}".format(
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


# =============================================================================
# TEST 1
# =============================================================================

def test_exact_root_distance_identity(states):
    print("=" * 90)
    print("TEST 1: EXACT ROOT-DISTANCE IDENTITY")
    print("=" * 90)

    failures = 0

    for st in states:
        lhs = st["depth"]

        rhs = min(
            v2(st["D"] - st["T"]),
            v2(st["n"] + st["c"]),
        )

        if lhs != rhs:
            failures += 1

            if failures <= 20:
                print(
                    "mismatch n={} frame={} D={} T={} "
                    "v2(D-T)={} w={} depth={} predicted={}".format(
                        st["n"],
                        st["frame"],
                        st["D"],
                        st["T"],
                        v2(st["D"] - st["T"]),
                        v2(st["n"] + st["c"]),
                        lhs,
                        rhs,
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_threshold_theorem(states, max_d=12):
    print("=" * 90)
    print("TEST 2: ROOT-DISTANCE THRESHOLD THEOREM")
    print("=" * 90)

    checks = 0
    failures = 0

    for st in states:
        for d in range(1, max_d + 1):

            lhs = st["depth"] >= d

            mod = 1 << d

            if st["frame"] == "A":
                T = st["S"] + 6
                c = 9
            else:
                T = st["S"] + 2
                c = 3

            rhs = (
                (st["D"] - T) % mod == 0
                and
                (st["n"] + c) % mod == 0
            )

            checks += 1

            if lhs != rhs:
                failures += 1

                if failures <= 20:
                    print(
                        "mismatch n={} frame={} d={} "
                        "depth={} D={} T={} n+c={} "
                        "lhs={} rhs={}".format(
                            st["n"],
                            st["frame"],
                            d,
                            st["depth"],
                            st["D"],
                            T,
                            st["n"] + c,
                            lhs,
                            rhs,
                        )
                    )

    print("checks={}".format(checks))
    print("failures={}".format(failures))
    print()

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_root_residue_depth(states, k):
    print("=" * 90)
    print(
        "TEST 3: ROOT RESIDUE -> TRUNCATED DEPTH k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(set)

    for st in states:
        roots = square_roots_mod_2k(
            st["Delta"],
            k,
        )

        for r in roots:
            depth = predicted_root_depth(
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
                r,
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
            print(
                "    {} -> {}".format(
                    sig,
                    sorted(depths),
                )
            )

    print()

    return len(ambiguous)


# =============================================================================
# TEST 4
# =============================================================================

def test_distance_class_structure(states, k):
    print("=" * 90)
    print(
        "TEST 4: ROOT DISTANCE CLASSIFICATION k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(set)

    for st in states:
        roots = square_roots_mod_2k(
            st["Delta"],
            k,
        )

        for r in roots:

            dist = trunc_v2(
                r - st["T"],
                k,
            )

            depth = predicted_root_depth(
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
                dist,
            )

            buckets[sig].add(depth)

    ambiguous = {
        sig: values
        for sig, values in buckets.items()
        if len(values) > 1
    }

    print("distance signatures={}".format(len(buckets)))
    print("ambiguous distance signatures={}".format(
        len(ambiguous)
    ))

    if ambiguous:
        print("first ambiguous:")
        for sig, values in list(ambiguous.items())[:10]:
            print(
                "    {} -> {}".format(
                    sig,
                    sorted(values),
                )
            )

    print()

    return len(ambiguous)


# =============================================================================
# TEST 5
# =============================================================================

def test_root_distance_distributions(states, k):
    print("=" * 90)
    print(
        "TEST 5: ROOT DISTANCE / DEPTH DISTRIBUTION k={}".format(k)
    )
    print("=" * 90)

    distribution = Counter()
    examples = []

    for st in states:
        roots = square_roots_mod_2k(
            st["Delta"],
            k,
        )

        classes = defaultdict(list)

        for r in roots:

            dist = trunc_v2(
                r - st["T"],
                k,
            )

            depth = predicted_root_depth(
                st["frame"],
                st["n"],
                st["S"],
                r,
                k,
            )

            classes[depth].append(r)

        distribution[
            (st["frame"], tuple(sorted(classes.keys())))
        ] += 1

        if (
            len(classes) > 1
            and len(examples) < MAX_EXAMPLES
        ):
            examples.append(
                (
                    st,
                    {
                        depth: sorted(rs)
                        for depth, rs in classes.items()
                    },
                )
            )

    print("depth-class distribution:")

    for key, count in sorted(distribution.items()):
        print("    {}: {}".format(key, count))

    print()

    if examples:
        print("examples:")
        for st, classes in examples:
            print(
                "    n={} frame={} S={}".format(
                    st["n"],
                    st["frame"],
                    st["S"],
                )
            )

            for depth in sorted(classes):
                roots = classes[depth]

                compact = roots[:16]

                print(
                    "        depth={} roots={}{}".format(
                        depth,
                        compact,
                        " ..." if len(roots) > 16 else "",
                    )
                )

                pairs = []

                for r in compact[:8]:
                    pairs.append(
                        (
                            r,
                            trunc_v2(
                                r - st["T"],
                                k,
                            ),
                        )
                    )

                print(
                    "            (root,distance)={}".format(
                        pairs
                    )
                )

    print()

    return 0


# =============================================================================
# TEST 6
# =============================================================================

def test_actual_root(st, k):
    """
    Diagnostic for one actual factorization.
    """

    roots = square_roots_mod_2k(
        st["Delta"],
        k,
    )

    actual = st["D"] % (1 << k)

    target = st["T"] % (1 << k)

    actual_distance = trunc_v2(
        actual - target,
        k,
    )

    actual_depth = predicted_root_depth(
        st["frame"],
        st["n"],
        st["S"],
        actual,
        k,
    )

    return (
        actual,
        target,
        actual_distance,
        actual_depth,
        roots,
    )


# =============================================================================
# EXAMPLES
# =============================================================================

def show_examples(states, k):
    print("=" * 90)
    print("TEST 6: EXAMPLES")
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

    by_n = {
        st["n"]: st
        for st in states
    }

    for n in wanted:
        st = by_n.get(n)

        if st is None:
            continue

        actual, target, distance, pred_depth, roots = \
            test_actual_root(st, k)

        print(
            "n={} p={} q={} frame={}".format(
                st["n"],
                st["p"],
                st["q"],
                st["frame"],
            )
        )

        print("    S={}".format(st["S"]))
        print("    Delta={}".format(st["Delta"]))
        print("    D(actual)={}".format(st["D"]))
        print("    T={}".format(st["T"]))
        print("    actual root mod 2^k={}".format(actual))
        print("    target T mod 2^k={}".format(target))
        print("    v2(D-T)={} ".format(distance))
        print("    v2(n+c)={}".format(
            v2(st["n"] + st["c"])
        ))
        print("    actual depth={}".format(
            st["depth"]
        ))
        print("    predicted depth={}".format(
            pred_depth
        ))
        print("    root count={}".format(
            len(roots)
        ))

        rows = []

        for r in roots[:24]:
            dist = trunc_v2(
                r - st["T"],
                k,
            )

            depth = predicted_root_depth(
                st["frame"],
                st["n"],
                st["S"],
                r,
                k,
            )

            rows.append(
                (
                    r,
                    dist,
                    depth,
                )
            )

        print(
            "    root -> (v2(D-T),depth):"
        )

        print("        {}".format(rows))

        print()


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 90)
    print("EXPERIMENT 677 START")
    print("=" * 90)
    print()

    print("prime limit={}".format(PRIME_LIMIT))

    primes, states = generate_states(
        PRIME_LIMIT
    )

    print("odd primes={}".format(
        len(primes)
    ))

    print("semiprimes={}".format(
        len(states)
    ))

    print()

    total_failures = 0

    total_failures += test_baseline(
        states
    )

    total_failures += test_exact_root_distance_identity(
        states
    )

    total_failures += test_threshold_theorem(
        states,
        max_d=12,
    )

    for k in [6, 8, 10, 12]:

        total_failures += test_root_residue_depth(
            states,
            k,
        )

        total_failures += test_distance_class_structure(
            states,
            k,
        )

        test_root_distance_distributions(
            states,
            k,
        )

    show_examples(
        states,
        12,
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("The new candidate identity is:")

    print()
    print("    FRAME A:")
    print("        T = S + 6")
    print("        c = 9")
    print()

    print("    FRAME B:")
    print("        T = S + 2")
    print("        c = 3")
    print()

    print("In both frames:")
    print()
    print("    depth")
    print("      = min(v2(D-T), v2(n+c))")
    print()

    print("Therefore:")
    print()
    print("    depth >= d")
    print()
    print("        iff")
    print()
    print("    D == T (mod 2^d)")
    print("    and")
    print("    n == -c (mod 2^d)")
    print()

    print("This converts the discriminant-root branch problem")
    print("into a direct 2-adic distance problem:")
    print()
    print("    root D")
    print("       ->")
    print("    distance from target T")
    print("       ->")
    print("    first branch failure")
    print("       ->")
    print("    depth.")
    print()

    print("TOTAL FAILURES={}".format(
        total_failures
    ))

    if total_failures == 0:
        print("STATUS=ALL ROOT-DISTANCE TESTS PASSED")
    else:
        print("STATUS=ROOT-DISTANCE COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 677 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
