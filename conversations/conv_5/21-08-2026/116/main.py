#!/usr/bin/env python3

from collections import defaultdict, Counter
from math import gcd


# =============================================================================
# EXPERIMENT 678
#
# CORRECTED DISCRIMINANT-ROOT TARGET THEOREM
#
# Frame A:
#
#     2C = S - D + 6
#     D target = S + 6
#     c = 9
#
# Frame B:
#
#     2C = S + D + 2
#     D target = -(S + 2)
#     c = 3
#
# Unified:
#
#     T =
#         S + 6     (A)
#        -S - 2     (B)
#
#     depth = min(v2(D-T), v2(n+c))
#
# Threshold:
#
#     depth >= d
#
#       iff
#
#     D == T (mod 2^d)
#     and
#     n == -c (mod 2^d)
#
# This version intentionally avoids enumerating all roots modulo
# 2^10 / 2^12 for every state.
# =============================================================================


PRIME_LIMIT = 6000
THRESHOLD_MAX = 16

INF = 10**9


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

        A = p - 3
        B = q + 3

        C = B

        # 2C = S - D + 6
        T = S + 6

    else:
        frame = "B"

        c = 3

        A = p + 1
        B = q - 3

        C = A

        # 2C = S + D + 2
        #
        # Therefore:
        # D - (-(S+2)) = D + S + 2
        T = -S - 2

    depth = v2(gcd(2 * C, n + c))

    return {
        "p": p,
        "q": q,
        "n": n,
        "S": S,
        "D": D,
        "frame": frame,
        "c": c,
        "A": A,
        "B": B,
        "C": C,
        "T": T,
        "depth": depth,
    }


def generate_states(limit: int):
    primes = odd_primes(limit)

    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(make_state(p, q))

    return primes, states


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
                    "depth={} predicted={}".format(
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

def test_corrected_root_distance(states):
    print("=" * 90)
    print("TEST 1: CORRECTED ROOT-DISTANCE IDENTITY")
    print("=" * 90)

    failures = 0

    for st in states:
        root_distance = v2(st["D"] - st["T"])
        n_bound = v2(st["n"] + st["c"])

        predicted = min(
            root_distance,
            n_bound,
        )

        if predicted != st["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    "mismatch n={} frame={} "
                    "D={} T={} "
                    "v2(D-T)={} w={} "
                    "depth={} predicted={}".format(
                        st["n"],
                        st["frame"],
                        st["D"],
                        st["T"],
                        root_distance,
                        n_bound,
                        st["depth"],
                        predicted,
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_threshold(states):
    print("=" * 90)
    print("TEST 2: CORRECTED THRESHOLD THEOREM")
    print("=" * 90)

    checks = 0
    failures = 0

    for st in states:

        for d in range(1, THRESHOLD_MAX + 1):
            mod = 1 << d

            lhs = st["depth"] >= d

            rhs = (
                (st["D"] - st["T"]) % mod == 0
                and
                (st["n"] + st["c"]) % mod == 0
            )

            checks += 1

            if lhs != rhs:
                failures += 1

                if failures <= 20:
                    print(
                        "mismatch n={} frame={} d={} "
                        "D={} T={} n+c={} "
                        "depth={} lhs={} rhs={}".format(
                            st["n"],
                            st["frame"],
                            d,
                            st["D"],
                            st["T"],
                            st["n"] + st["c"],
                            st["depth"],
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

def test_shift_identity(states):
    print("=" * 90)
    print("TEST 3: SHIFTED FACTOR IDENTITY")
    print("=" * 90)

    failures = 0

    for st in states:

        if st["frame"] == "A":
            lhs = st["D"] - (st["S"] + 6)

        else:
            lhs = st["D"] + st["S"] + 2

        rhs = 2 * st["C"]

        if v2(lhs) != v2(rhs):
            failures += 1

            if failures <= 20:
                print(
                    "mismatch n={} frame={} "
                    "lhs={} rhs={} v2(lhs)={} v2(rhs)={}".format(
                        st["n"],
                        st["frame"],
                        lhs,
                        rhs,
                        v2(lhs),
                        v2(rhs),
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_threshold_counts(states):
    print("=" * 90)
    print("TEST 4: SURVIVAL COUNTS")
    print("=" * 90)

    counts = Counter()

    for st in states:
        counts[(st["frame"], st["depth"])] += 1

    for frame in ("A", "B"):
        print("FRAME {}".format(frame))

        for d in range(1, THRESHOLD_MAX + 1):
            surviving = sum(
                1
                for st in states
                if st["frame"] == frame
                and st["depth"] >= d
            )

            print(
                "    d={:2d} surviving={}".format(
                    d,
                    surviving,
                )
            )

        print()

    return 0


# =============================================================================
# TEST 5
# =============================================================================

def test_root_target_signatures(states, k=8):
    print("=" * 90)
    print(
        "TEST 5: ROOT TARGET RESIDUE -> TRUNCATED DEPTH k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k
    buckets = defaultdict(set)

    for st in states:

        D = st["D"] % mod
        T = st["T"] % mod

        root_distance = trunc_v2(
            D - T,
            k,
        )

        n_bound = min(
            v2(st["n"] + st["c"]),
            k,
        )

        depth = min(
            root_distance,
            n_bound,
        )

        signature = (
            st["frame"],
            st["n"] % mod,
            st["S"] % mod,
            D,
            T,
        )

        buckets[signature].add(depth)

    ambiguous = sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )

    print("signatures={}".format(len(buckets)))
    print("ambiguous={}".format(ambiguous))

    if ambiguous:
        print("first ambiguous signatures:")

        shown = 0

        for sig, values in buckets.items():

            if len(values) <= 1:
                continue

            print(
                "    {} -> {}".format(
                    sig,
                    sorted(values),
                )
            )

            shown += 1

            if shown >= 10:
                break

    print()

    return ambiguous


# =============================================================================
# TEST 6
# =============================================================================

def test_distance_signature(states, k=8):
    print("=" * 90)
    print(
        "TEST 6: ROOT 2-ADIC DISTANCE -> TRUNCATED DEPTH k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(set)

    for st in states:

        D = st["D"] % mod
        T = st["T"] % mod

        distance = trunc_v2(
            D - T,
            k,
        )

        n_bound = min(
            v2(st["n"] + st["c"]),
            k,
        )

        depth = min(
            distance,
            n_bound,
        )

        signature = (
            st["frame"],
            st["n"] % mod,
            st["S"] % mod,
            distance,
        )

        buckets[signature].add(depth)

    ambiguous = sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )

    print("distance signatures={}".format(len(buckets)))
    print("ambiguous={}".format(ambiguous))

    if ambiguous:
        print("first ambiguous:")

        shown = 0

        for sig, values in buckets.items():

            if len(values) <= 1:
                continue

            print(
                "    {} -> {}".format(
                    sig,
                    sorted(values),
                )
            )

            shown += 1

            if shown >= 10:
                break

    print()

    return ambiguous


# =============================================================================
# TEST 7
# =============================================================================

def test_root_lift_sanity(states):
    """
    Fast sanity check: use only the actual discriminant root D.
    No exhaustive root lifting.
    """

    print("=" * 90)
    print("TEST 7: ACTUAL DISCRIMINANT ROOT SANITY")
    print("=" * 90)

    failures = 0

    for st in states:

        delta = st["S"] * st["S"] - 4 * st["n"]

        if delta != st["D"] * st["D"]:
            failures += 1

            if failures <= 10:
                print(
                    "discriminant mismatch n={} D={} delta={}".format(
                        st["n"],
                        st["D"],
                        delta,
                    )
                )
            continue

        reconstructed_depth = min(
            v2(st["D"] - st["T"]),
            v2(st["n"] + st["c"]),
        )

        if reconstructed_depth != st["depth"]:
            failures += 1

            if failures <= 10:
                print(
                    "root mismatch n={} frame={} "
                    "D={} T={} actual={} predicted={}".format(
                        st["n"],
                        st["frame"],
                        st["D"],
                        st["T"],
                        st["depth"],
                        reconstructed_depth,
                    )
                )

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# =============================================================================
# TRUNCATED v2
# =============================================================================

def trunc_v2(x, k):
    x &= (1 << k) - 1

    if x == 0:
        return k

    return (x & -x).bit_length() - 1


# =============================================================================
# EXAMPLES
# =============================================================================

def show_examples(states):
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

    by_n = {
        st["n"]: st
        for st in states
    }

    for n in wanted:
        st = by_n.get(n)

        if st is None:
            continue

        distance = v2(
            st["D"] - st["T"]
        )

        w = v2(
            st["n"] + st["c"]
        )

        predicted = min(
            distance,
            w,
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
        print("    D={}".format(st["D"]))

        print("    T={}".format(
            st["T"]
        ))

        print("    D-T={}".format(
            st["D"] - st["T"]
        ))

        print("    v2(D-T)={}".format(
            distance
        ))

        print("    n+c={}".format(
            st["n"] + st["c"]
        ))

        print("    v2(n+c)={}".format(
            w
        ))

        print("    depth={}".format(
            st["depth"]
        ))

        print("    predicted={}".format(
            predicted
        ))

        print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 678 START")
    print("=" * 90)
    print()

    print("prime limit={}".format(
        PRIME_LIMIT
    ))

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

    total_failures += test_shift_identity(
        states
    )

    total_failures += test_corrected_root_distance(
        states
    )

    total_failures += test_threshold(
        states
    )

    test_threshold_counts(
        states
    )

    total_failures += test_root_target_signatures(
        states,
        6,
    )

    total_failures += test_distance_signature(
        states,
        6,
    )

    # Only one higher-resolution check.
    # This avoids the expensive exhaustive root enumeration
    # which caused the previous experiment to run for too long.

    total_failures += test_root_target_signatures(
        states,
        10,
    )

    total_failures += test_distance_signature(
        states,
        10,
    )

    total_failures += test_root_lift_sanity(
        states
    )

    show_examples(
        states
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Corrected frame targets:")
    print()

    print("    FRAME A:")
    print("        T = S + 6")
    print("        c = 9")
    print()

    print("    FRAME B:")
    print("        T = -S - 2")
    print("        c = 3")
    print()

    print("Unified candidate:")
    print()
    print("    depth")
    print("      = min(v2(D-T), v2(n+c))")
    print()

    print("Threshold form:")
    print()
    print("    depth >= d")
    print("        iff")
    print("    D == T (mod 2^d)")
    print("    and")
    print("    n == -c (mod 2^d)")
    print()

    print("The previous Frame-B mismatch came from using")
    print("T = S + 2 instead of")
    print("T = -(S + 2).")
    print()

    print("TOTAL FAILURES={}".format(
        total_failures
    ))

    if total_failures == 0:
        print("STATUS=ALL CORRECTED ROOT-DISTANCE TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 678 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
