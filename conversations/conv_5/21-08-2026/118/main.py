#!/usr/bin/env python3

from collections import defaultdict
from math import gcd


# =============================================================================
# EXPERIMENT 680
#
# MINIMAL EXTRA INFORMATION AFTER (n,D)
#
# Experiment 679 established:
#
#     (frame, n mod 2^k, D mod 2^k)
#
# does NOT generally determine depth.
#
# But:
#
#     (frame, n mod 2^k, S mod 2^k, D mod 2^k)
#
# does determine it.
#
# The new question is:
#
#     HOW MANY LOW BITS OF S
#     are actually missing?
#
# We test:
#
#     (frame,
#      n mod 2^k,
#      D mod 2^k,
#      S mod 2^r)
#
# for increasing r.
#
# For each k we report the smallest r for which the signature
# determines the truncated depth exactly.
#
# This gives an empirical "information boundary" between D and S.
#
# We also test whether only ONE extra S bit is sufficient.
#
# =============================================================================


PRIME_LIMIT = 6000

K_VALUES = (4, 6, 8, 10, 12)

INF = 10**9


# =============================================================================
# BASIC ARITHMETIC
# =============================================================================

def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


# =============================================================================
# PRIME SIEVE
# =============================================================================

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
            a[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

        p += 1

    return [i for i in range(limit + 1) if a[i]]


def odd_primes(limit: int):
    return [p for p in sieve(limit) if p & 1]


# =============================================================================
# STATE GENERATION
# =============================================================================

def make_state(p: int, q: int):
    n = p * q
    S = p + q
    D = p - q

    # Frame convention used in all previous experiments.
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
# BUCKET TEST
# =============================================================================

def ambiguity_count(
    states,
    k,
    r,
):
    """
    Signature:

        (frame,
         n mod 2^k,
         D mod 2^k,
         S mod 2^r)

    Target:

        min(depth, k)
    """

    nmod = 1 << k
    smod = 1 << r
    dmod = nmod

    buckets = defaultdict(set)

    for st in states:

        signature = (
            st["frame"],
            st["n"] % nmod,
            st["D"] % dmod,
            st["S"] % smod,
        )

        value = min(
            st["depth"],
            k,
        )

        buckets[signature].add(
            value
        )

    ambiguous = sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )

    return len(buckets), ambiguous


# =============================================================================
# TEST 1
# =============================================================================

def test_minimal_s_bits(states, k):
    print("=" * 90)
    print(
        "TEST 1: MINIMUM EXTRA S BITS AFTER (n,D) k={}".format(k)
    )
    print("=" * 90)

    print(
        "r= 0 means NO S information"
    )

    exact_r = None

    # r=0 ... k
    for r in range(0, k + 1):

        signatures, ambiguous = ambiguity_count(
            states,
            k,
            r,
        )

        status = (
            "EXACT"
            if ambiguous == 0
            else "AMBIGUOUS"
        )

        print(
            "    r={:2d} signatures={:7d} ambiguous={:6d} {}".format(
                r,
                signatures,
                ambiguous,
                status,
            )
        )

        if ambiguous == 0 and exact_r is None:
            exact_r = r

    print()

    if exact_r is None:
        print(
            "minimal exact r=None"
        )
    else:
        print(
            "minimal exact r={}".format(
                exact_r
            )
        )

    print()

    return exact_r


# =============================================================================
# TEST 2
# =============================================================================

def test_one_bit_s(states, k):
    print("=" * 90)
    print(
        "TEST 2: DOES ONE EXTRA S BIT SUFFICE? k={}".format(k)
    )
    print("=" * 90)

    signatures, ambiguous = ambiguity_count(
        states,
        k,
        1,
    )

    print(
        "signature=(frame,n mod 2^k,D mod 2^k,S mod 2)"
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

def test_same_nd_different_s(states, k):
    """
    Show concrete cases where identical (n,D) information produces
    different depths, together with the S values that separate them.
    """

    print("=" * 90)
    print(
        "TEST 3: CONCRETE (n,D) COLLISIONS k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(list)

    for st in states:

        signature = (
            st["frame"],
            st["n"] % mod,
            st["D"] % mod,
        )

        buckets[signature].append(
            st
        )

    shown = 0

    for signature, items in buckets.items():

        depths = sorted(
            {
                min(
                    st["depth"],
                    k,
                )
                for st in items
            }
        )

        if len(depths) <= 1:
            continue

        print(
            "signature={} depths={}".format(
                signature,
                depths,
            )
        )

        # Show at most one state for each depth.
        representatives = {}

        for st in items:

            d = min(
                st["depth"],
                k,
            )

            if d not in representatives:
                representatives[d] = st

        for d in sorted(representatives):

            st = representatives[d]

            print(
                "    depth={:2d} n={} S={} D={} p={} q={}".format(
                    d,
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
            "no ambiguous signatures found"
        )

    print()

    return 0


# =============================================================================
# TEST 4
# =============================================================================

def test_d_plus_s_minus_n_identity(states):
    """
    Algebraic sanity check:

        D^2 = S^2 - 4n

    and therefore:

        S^2 = D^2 + 4n

    We verify that every state satisfies the identity and then
    investigate the residual sign/branch ambiguity.
    """

    print("=" * 90)
    print(
        "TEST 4: DISCRIMINANT CONSISTENCY"
    )
    print("=" * 90)

    failures = 0

    for st in states:

        lhs = st["D"] * st["D"]
        rhs = st["S"] * st["S"] - 4 * st["n"]

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
# TEST 5
# =============================================================================

def test_shifted_s_bit_progression(states, k):
    """
    Directly tests whether the extra S information helps
    monotonically.

    We compare:

        r bits of S

    against:

        r+1 bits of S.

    We expect ambiguity to be non-increasing.
    """

    print("=" * 90)
    print(
        "TEST 5: S-BIT INFORMATION MONOTONICITY k={}".format(k)
    )
    print("=" * 90)

    previous = None
    failures = 0

    for r in range(0, k + 1):

        _, ambiguous = ambiguity_count(
            states,
            k,
            r,
        )

        if previous is not None and ambiguous > previous:

            failures += 1

            print(
                "monotonicity failure r={} previous={} current={}".format(
                    r,
                    previous,
                    ambiguous,
                )
            )

        previous = ambiguous

    print(
        "checked={} S-bit levels".format(
            k + 1
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
# TEST 6
# =============================================================================

def test_compare_to_full_s(states, k):
    """
    Control:

        full S mod 2^k

    should resolve the depth exactly.
    """

    print("=" * 90)
    print(
        "TEST 6: FULL S CONTROL k={}".format(k)
    )
    print("=" * 90)

    signatures, ambiguous = ambiguity_count(
        states,
        k,
        k,
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
            "    S mod 2^4={}".format(
                st["S"] % 16
            )
        )

        print(
            "    S mod 2^8={}".format(
                st["S"] % 256
            )
        )

        print(
            "    D mod 2^8={}".format(
                st["D"] % 256
            )
        )

        print(
            "    n mod 2^8={}".format(
                st["n"] % 256
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
    print("EXPERIMENT 680 START")
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
    # MAIN INFORMATION-BOUNDARY TEST
    # -------------------------------------------------------------------------

    for k in K_VALUES:

        exact_r = test_minimal_s_bits(
            states,
            k,
        )

        exact_results[k] = exact_r

        # Specifically test one extra S bit.
        test_one_bit_s(
            states,
            k,
        )

        # Show concrete collisions.
        test_same_nd_different_s(
            states,
            k,
        )

        # Monotonicity.
        total_failures += test_shifted_s_bit_progression(
            states,
            k,
        )

        # Full-S control.
        total_failures += test_compare_to_full_s(
            states,
            k,
        )

    # -------------------------------------------------------------------------
    # ALGEBRAIC CONTROL
    # -------------------------------------------------------------------------

    total_failures += test_d_plus_s_minus_n_identity(
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
        "Experiment 679 established:"
    )
    print()

    print(
        "    (n,D) mod 2^k"
    )
    print(
        "        is generally insufficient."
    )
    print()

    print(
        "Experiment 680 asks:"
    )
    print()

    print(
        "    How many bits of S are actually missing?"
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
        "     n mod 2^k,"
    )
    print(
        "     D mod 2^k,"
    )
    print(
        "     S mod 2^r)"
    )
    print()

    print(
        "The key quantity is minimal r."
    )
    print()

    for k in K_VALUES:

        r = exact_results[k]

        print(
            "    k={} -> minimal S bits={}".format(
                k,
                r,
            )
        )

    print()

    print(
        "If r is small compared with k,"
    )
    print(
        "then only a small amount of the symmetric"
    )
    print(
        "coordinate is needed to resolve the"
    )
    print(
        "discriminant-root ambiguity."
    )

    print()

    print(
        "If r grows roughly with k,"
    )
    print(
        "then S carries genuinely growing"
    )
    print(
        "information not present in (n,D)."
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
    print("EXPERIMENT 680 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
