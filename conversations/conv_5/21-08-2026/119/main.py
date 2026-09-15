#!/usr/bin/env python3

from collections import defaultdict
from math import gcd


# =============================================================================
# EXPERIMENT 681
#
# INFORMATION TRADEOFF:
#
#     CAN EXTRA n-BITS REPLACE THE MISSING S-BITS?
#
# Experiment 679:
#
#     (frame, n mod 2^k, D mod 2^k)
#
#     is not sufficient.
#
# Experiment 680:
#
#     (frame,
#      n mod 2^k,
#      D mod 2^k,
#      S mod 2^r)
#
# was tested.
#
# Minimal r:
#
#     k=4  -> 4
#     k=6  -> 6
#     k=8  -> 8
#     k=10 -> 10
#     k=12 -> 10
#
# This experiment asks the dual question:
#
#     Instead of supplying extra S bits,
#     how many EXTRA n bits are required?
#
# We test:
#
#     (frame,
#      D mod 2^k,
#      n mod 2^t)
#
# for t = 0 ... LIMIT.
#
# We deliberately separate:
#
#     t < k
#     t = k
#     t > k
#
# because the discriminant equation
#
#     S^2 = D^2 + 4n
#
# may make higher n precision relevant to the
# missing square-root branch.
#
# A second test measures the information tradeoff:
#
#     (extra n bits, extra S bits)
#
# to find minimal combinations that determine
# truncated depth.
#
# No O(N^2) state comparisons are used.
# Everything is bucket-based.
# =============================================================================


PRIME_LIMIT = 6000

K_VALUES = (4, 6, 8, 10, 12)

# Extra precision explored for n.
MAX_EXTRA_N_BITS = 16


# =============================================================================
# BASIC ARITHMETIC
# =============================================================================

INF = 10 ** 9


def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


# =============================================================================
# SIEVE
# =============================================================================

def sieve(limit: int):
    values = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        values[0] = 0

    if limit >= 1:
        values[1] = 0

    p = 2

    while p * p <= limit:

        if values[p]:

            start = p * p

            count = ((limit - start) // p) + 1

            values[start:limit + 1:p] = b"\x00" * count

        p += 1

    return [
        i
        for i in range(limit + 1)
        if values[i]
    ]


def odd_primes(limit: int):
    return [
        p
        for p in sieve(limit)
        if p & 1
    ]


# =============================================================================
# STATE CONSTRUCTION
# =============================================================================

def make_state(p: int, q: int):

    n = p * q
    S = p + q
    D = p - q

    # Same frame convention used in Experiments 665-680.
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
        T = -S - 2

    depth = v2(
        gcd(
            2 * C,
            n + c,
        )
    )

    return {
        "p": p,
        "q": q,
        "n": n,
        "S": S,
        "D": D,
        "A": A,
        "B": B,
        "C": C,
        "T": T,
        "c": c,
        "frame": frame,
        "depth": depth,
    }


def generate_states(limit: int):

    primes = odd_primes(limit)

    states = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            states.append(
                make_state(p, q)
            )

    return primes, states


# =============================================================================
# BASELINE
# =============================================================================

def test_baseline(states):

    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
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

    print(
        "checked={}".format(
            len(states)
        )
    )

    print(
        "failures={}".format(
            failures
        )
    )

    print()

    return failures


# =============================================================================
# GENERIC BUCKET TEST
# =============================================================================

def bucket_signature(
    st,
    k,
    t,
):
    """
    Signature:

        frame
        D mod 2^k
        n mod 2^t

    Target:

        depth truncated at k.

    t=0 means no n information.
    """

    d_mod = 1 << k

    if t == 0:
        n_mod = 1
        n_value = 0
    else:
        n_mod = 1 << t
        n_value = st["n"] % n_mod

    return (
        st["frame"],
        st["D"] % d_mod,
        n_value,
    )


def ambiguity_count_nd(
    states,
    k,
    t,
):
    buckets = defaultdict(set)

    for st in states:

        signature = bucket_signature(
            st,
            k,
            t,
        )

        target = min(
            st["depth"],
            k,
        )

        buckets[signature].add(
            target
        )

    ambiguous = sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )

    return (
        len(buckets),
        ambiguous,
    )


# =============================================================================
# TEST 1
# =============================================================================

def test_extra_n_bits(states, k):

    print("=" * 90)
    print(
        "TEST 1: MINIMUM n-BITS AFTER (D MOD 2^k) k={}".format(
            k
        )
    )
    print("=" * 90)

    exact_t = None

    max_t = k + MAX_EXTRA_N_BITS

    for t in range(0, max_t + 1):

        signatures, ambiguous = ambiguity_count_nd(
            states,
            k,
            t,
        )

        if t == 0:
            status_prefix = "NO n INFO"
        elif t < k:
            status_prefix = "LESS THAN k"
        elif t == k:
            status_prefix = "BASELINE n PRECISION"
        else:
            status_prefix = "EXTRA n BITS"

        status = (
            "EXACT"
            if ambiguous == 0
            else "AMBIGUOUS"
        )

        print(
            "    t={:2d} signatures={:8d} ambiguous={:7d} "
            "{} {}".format(
                t,
                signatures,
                ambiguous,
                status_prefix,
                status,
            )
        )

        if ambiguous == 0 and exact_t is None:

            exact_t = t

    print()

    print(
        "minimal exact t={}".format(
            exact_t
        )
    )

    print()

    return exact_t


# =============================================================================
# TEST 2
# =============================================================================

def test_one_extra_n_bit(states, k):

    print("=" * 90)
    print(
        "TEST 2: DOES ONE EXTRA n BIT SUFFICE? k={}".format(
            k
        )
    )
    print("=" * 90)

    t = k + 1

    signatures, ambiguous = ambiguity_count_nd(
        states,
        k,
        t,
    )

    print(
        "signature=("
        "frame,"
        "D mod 2^k,"
        "n mod 2^(k+1)"
        ")"
    )

    print(
        "signatures={}".format(
            signatures
        )
    )

    print(
        "ambiguous={}".format(
            ambiguous
        )
    )

    print()

    return ambiguous


# =============================================================================
# TEST 3
# =============================================================================

def test_n_precision_monotonicity(states, k):

    print("=" * 90)
    print(
        "TEST 3: n-BIT INFORMATION MONOTONICITY k={}".format(
            k
        )
    )
    print("=" * 90)

    previous = None
    failures = 0

    for t in range(0, k + MAX_EXTRA_N_BITS + 1):

        _, ambiguous = ambiguity_count_nd(
            states,
            k,
            t,
        )

        if previous is not None:

            if ambiguous > previous:

                failures += 1

                print(
                    "monotonicity failure "
                    "t={} previous={} current={}".format(
                        t,
                        previous,
                        ambiguous,
                    )
                )

        previous = ambiguous

    print(
        "checked={} n-bit levels".format(
            k + MAX_EXTRA_N_BITS + 1
        )
    )

    print(
        "failures={}".format(
            failures
        )
    )

    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_n_s_tradeoff(states, k):

    print("=" * 90)
    print(
        "TEST 4: n-BITS / S-BITS INFORMATION TRADEOFF k={}".format(
            k
        )
    )
    print("=" * 90)

    # We search all small combinations.
    #
    # Signature:
    #
    #     frame
    #     D mod 2^k
    #     n mod 2^t
    #     S mod 2^r
    #
    # and look for exact combinations.
    #
    # This does not attempt to prove minimality globally.
    # It simply identifies the observed boundary.

    results = []

    max_t = k + 6
    max_r = k

    for t in range(0, max_t + 1):

        for r in range(0, max_r + 1):

            buckets = defaultdict(set)

            d_mod = 1 << k

            if t == 0:
                n_mask = None
            else:
                n_mask = (1 << t) - 1

            if r == 0:
                s_mask = None
            else:
                s_mask = (1 << r) - 1

            for st in states:

                if n_mask is None:
                    n_value = 0
                else:
                    n_value = st["n"] & n_mask

                if s_mask is None:
                    s_value = 0
                else:
                    s_value = st["S"] & s_mask

                signature = (
                    st["frame"],
                    st["D"] % d_mod,
                    n_value,
                    s_value,
                )

                buckets[signature].add(
                    min(
                        st["depth"],
                        k,
                    )
                )

            ambiguous = sum(
                1
                for values in buckets.values()
                if len(values) > 1
            )

            if ambiguous == 0:

                results.append(
                    (t, r)
                )

    if not results:

        print(
            "no exact (t,r) combinations found"
        )

        print()

        return []

    # Pareto frontier:
    #
    # A pair is dominated if another pair uses
    # no more n bits and no more S bits.

    frontier = []

    for t, r in sorted(results):

        dominated = False

        for t2, r2 in results:

            if (
                t2 <= t
                and r2 <= r
                and (t2 < t or r2 < r)
            ):

                dominated = True
                break

        if not dominated:

            frontier.append(
                (t, r)
            )

    print(
        "exact combinations={}".format(
            len(results)
        )
    )

    print(
        "Pareto frontier:"
    )

    for t, r in frontier:

        print(
            "    n bits={:2d}  S bits={:2d}".format(
                t,
                r,
            )
        )

    print()

    return frontier


# =============================================================================
# TEST 5
# =============================================================================

def test_discriminant_reconstruction(states):

    print("=" * 90)
    print(
        "TEST 5: DISCRIMINANT RECONSTRUCTION"
    )
    print("=" * 90)

    failures = 0

    for st in states:

        lhs = st["S"] * st["S"]

        rhs = (
            st["D"] * st["D"]
            + 4 * st["n"]
        )

        if lhs != rhs:

            failures += 1

            if failures <= 10:

                print(
                    "mismatch n={} S={} D={}".format(
                        st["n"],
                        st["S"],
                        st["D"],
                    )
                )

    print(
        "checked={}".format(
            len(states)
        )
    )

    print(
        "failures={}".format(
            failures
        )
    )

    print()

    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def show_collision_examples(
    states,
    k,
    t,
):
    print("=" * 90)
    print(
        "EXAMPLES: (D,n)-COLLISIONS k={} t={}".format(
            k,
            t,
        )
    )
    print("=" * 90)

    buckets = defaultdict(list)

    for st in states:

        signature = bucket_signature(
            st,
            k,
            t,
        )

        buckets[signature].append(
            st
        )

    shown = 0

    for signature, items in buckets.items():

        depth_groups = defaultdict(list)

        for st in items:

            d = min(
                st["depth"],
                k,
            )

            depth_groups[d].append(
                st
            )

        if len(depth_groups) <= 1:
            continue

        print(
            "signature={} depths={}".format(
                signature,
                sorted(depth_groups),
            )
        )

        for depth in sorted(depth_groups):

            st = depth_groups[depth][0]

            print(
                "    depth={:2d} "
                "n={} S={} D={} p={} q={}".format(
                    depth,
                    st["n"],
                    st["S"],
                    st["D"],
                    st["p"],
                    st["q"],
                )
            )

        print()

        shown += 1

        if shown >= 10:
            break

    if shown == 0:
        print(
            "no collision examples"
        )

    print()


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

        print(
            "n={} p={} q={} frame={}".format(
                st["n"],
                st["p"],
                st["q"],
                st["frame"],
            )
        )

        print(
            "    S={}".format(
                st["S"]
            )
        )

        print(
            "    D={}".format(
                st["D"]
            )
        )

        print(
            "    n+c={}".format(
                st["n"] + st["c"]
            )
        )

        print(
            "    C={}".format(
                st["C"]
            )
        )

        print(
            "    v2(C)={}".format(
                v2(st["C"])
            )
        )

        print(
            "    v2(n+c)={}".format(
                v2(st["n"] + st["c"])
            )
        )

        print(
            "    depth={}".format(
                st["depth"]
            )
        )

        print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 681 START")
    print("=" * 90)
    print()

    print(
        "prime limit={}".format(
            PRIME_LIMIT
        )
    )

    primes, states = generate_states(
        PRIME_LIMIT
    )

    print(
        "odd primes={}".format(
            len(primes)
        )
    )

    print(
        "semiprimes={}".format(
            len(states)
        )
    )

    print()

    total_failures = 0

    exact_results = {}

    # -------------------------------------------------------------------------
    # BASELINE
    # -------------------------------------------------------------------------

    total_failures += test_baseline(
        states
    )

    # -------------------------------------------------------------------------
    # MAIN TEST
    # -------------------------------------------------------------------------

    for k in K_VALUES:

        exact_t = test_extra_n_bits(
            states,
            k,
        )

        exact_results[k] = exact_t

        test_one_extra_n_bit(
            states,
            k,
        )

        total_failures += test_n_precision_monotonicity(
            states,
            k,
        )

        test_n_s_tradeoff(
            states,
            k,
        )

        show_collision_examples(
            states,
            k,
            k,
        )

    # -------------------------------------------------------------------------
    # CONTROL
    # -------------------------------------------------------------------------

    total_failures += test_discriminant_reconstruction(
        states
    )

    # -------------------------------------------------------------------------
    # EXAMPLES
    # -------------------------------------------------------------------------

    show_examples(
        states
    )

    # -------------------------------------------------------------------------
    # FINAL
    # -------------------------------------------------------------------------

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print(
        "Experiment 680 found:"
    )

    print()

    print(
        "    minimal S bits:"
    )

    for k in K_VALUES:

        print(
            "        k={:2d} -> r={}".format(
                k,
                exact_results[k],
            )
        )

    print()

    print(
        "Experiment 681 asks the dual question:"
    )

    print()

    print(
        "    Can extra precision in n"
    )

    print(
        "    replace the missing S information?"
    )

    print()

    print(
        "Tested signature:"
    )

    print()

    print(
        "    (frame,"
    )

    print(
        "     D mod 2^k,"
    )

    print(
        "     n mod 2^t)"
    )

    print()

    print(
        "The important comparison is:"
    )

    print()

    print(
        "    minimal t"
    )

    print(
        "versus"
    )

    print(
        "    k."
    )

    print()

    print(
        "If t ~= k,"
    )

    print(
        "then n precision does not compress the"
    )

    print(
        "missing branch information."
    )

    print()

    print(
        "If t < k,"
    )

    print(
        "then additional n bits are more efficient"
    )

    print(
        "than additional S bits."
    )

    print()

    print(
        "If t > k,"
    )

    print(
        "then the discriminant equation requires"
    )

    print(
        "higher n precision to resolve the"
    )

    print(
        "2-adic square-root branch."
    )

    print()

    print(
        "TOTAL FAILURES={}".format(
            total_failures
        )
    )

    if total_failures == 0:

        print(
            "STATUS=ALL CORE TESTS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )

    print("=" * 90)
    print("EXPERIMENT 681 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
