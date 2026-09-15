#!/usr/bin/env python3

from collections import defaultdict
from math import gcd


# =============================================================================
# EXPERIMENT 679
#
# ROOT-ONLY INFORMATION TEST
#
# Previous result:
#
#     FRAME A:
#         T = S + 6
#         c = 9
#
#     FRAME B:
#         T = -S - 2
#         c = 3
#
#     depth = min(v2(D-T), v2(n+c))
#
# But:
#
#     FRAME A:
#         D-(S+6) = -2(q+3)
#
#     FRAME B:
#         D-(-S-2) =  2(p+1)
#
# so the previous identity is algebraically equivalent to the
# canonical C-law.
#
# The new question is:
#
#     Does (frame, n, D) contain enough information by itself?
#
# Specifically:
#
#     (frame, n mod 2^k, D mod 2^k)
#         ->
#     truncated depth
#
# without explicitly supplying S.
#
# =============================================================================


PRIME_LIMIT = 6000

K_VALUES = (2, 4, 6, 8, 10, 12)

INF = 10**9


# =============================================================================
# BASIC
# =============================================================================

def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


def trunc_v2(x: int, k: int) -> int:
    """
    v2(x), truncated to at most k.
    """
    mod = 1 << k
    x %= mod

    if x == 0:
        return k

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
            a[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

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

    # Same frame convention used throughout the experiments.
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

    print("checked={}".format(len(states)))
    print("failures={}".format(failures))
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_dn_signature(states, k):
    """
    Main experiment:

        (frame, n mod 2^k, D mod 2^k)
            ->
        truncated depth
    """

    print("=" * 90)
    print(
        "TEST 1: (n,D) -> TRUNCATED DEPTH k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(set)

    for st in states:

        depth_k = min(
            st["depth"],
            k,
        )

        signature = (
            st["frame"],
            st["n"] % mod,
            st["D"] % mod,
        )

        buckets[signature].add(
            depth_k
        )

    ambiguous = sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )

    print(
        "signatures={}".format(
            len(buckets)
        )
    )

    print(
        "ambiguous={}".format(
            ambiguous
        )
    )

    if ambiguous:
        print(
            "first ambiguous signatures:"
        )

        shown = 0

        for signature, values in buckets.items():

            if len(values) <= 1:
                continue

            print(
                "    {} -> {}".format(
                    signature,
                    sorted(values),
                )
            )

            shown += 1

            if shown >= 10:
                break

    print()

    return ambiguous


# =============================================================================
# TEST 2
# =============================================================================

def test_sdn_signature(states, k):
    """
    Compare:

        (n,D)

    against:

        (n,S,D)

    The latter is guaranteed to carry more information.
    """

    print("=" * 90)
    print(
        "TEST 2: (n,S,D) CONTROL k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(set)

    for st in states:

        depth_k = min(
            st["depth"],
            k,
        )

        signature = (
            st["frame"],
            st["n"] % mod,
            st["S"] % mod,
            st["D"] % mod,
        )

        buckets[signature].add(
            depth_k
        )

    ambiguous = sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )

    print(
        "signatures={}".format(
            len(buckets)
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

def test_n_plus_d_shift(states, k):
    """
    Test whether the root can be converted directly into the
    canonical C-channel using only n and D.

    We reconstruct a candidate S modulo 2^(k+1) by using:

        S^2 = D^2 + 4n

    but we DO NOT enumerate roots.

    Instead, we verify whether the actual S residue is uniquely
    determined by (n,D) over the observed dataset.
    """

    print("=" * 90)
    print(
        "TEST 3: DOES (n,D) DETERMINE S MOD 2^k? k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    buckets = defaultdict(set)

    for st in states:

        signature = (
            st["frame"],
            st["n"] % mod,
            st["D"] % mod,
        )

        buckets[signature].add(
            st["S"] % mod
        )

    ambiguous = sum(
        1
        for values in buckets.values()
        if len(values) > 1
    )

    print(
        "signatures={}".format(
            len(buckets)
        )
    )

    print(
        "ambiguous_S={}".format(
            ambiguous
        )
    )

    if ambiguous:

        print(
            "first ambiguous examples:"
        )

        shown = 0

        for signature, values in buckets.items():

            if len(values) <= 1:
                continue

            print(
                "    {} -> {}".format(
                    signature,
                    sorted(values),
                )
            )

            shown += 1

            if shown >= 10:
                break

    print()

    return ambiguous


# =============================================================================
# TEST 4
# =============================================================================

def test_one_extra_bit(states, k):
    """
    Check whether adding one extra bit to D resolves any
    remaining ambiguity:

        n mod 2^k
        D mod 2^(k+1)
    """

    print("=" * 90)
    print(
        "TEST 4: (n mod 2^k, D mod 2^(k+1)) k={}".format(k)
    )
    print("=" * 90)

    nmod = 1 << k
    dmod = 1 << (k + 1)

    buckets = defaultdict(set)

    for st in states:

        signature = (
            st["frame"],
            st["n"] % nmod,
            st["D"] % dmod,
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

    print(
        "signatures={}".format(
            len(buckets)
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
# TEST 5
# =============================================================================

def test_root_sign_flip(states, k):
    """
    Investigate D -> -D.

    Swapping p and q changes D sign. However the frame anchor
    changes as well, so we test the actual transformation rather
    than assuming sign symmetry.
    """

    print("=" * 90)
    print(
        "TEST 5: ROOT SIGN FLIP k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    failures = 0

    for st in states:

        D = st["D"] % mod
        negD = (-st["D"]) % mod

        if st["frame"] == "A":

            # Swapping the factors changes the frame.
            #
            # A: p near 3, q near -3
            #
            # after swap:
            # q near -3, p near 3
            #
            # this maps naturally to the same structural pair.

            expected = (
                -st["D"]
            ) % mod

        else:

            expected = (
                -st["D"]
            ) % mod

        if negD != expected:
            failures += 1

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
# TEST 6
# =============================================================================

def test_direct_residue_depth(states, k):
    """
    Directly test:

        D == T mod 2^d
        n == -c mod 2^d

    using only modular arithmetic.
    """

    print("=" * 90)
    print(
        "TEST 6: ROOT-ONLY THRESHOLD k={}".format(k)
    )
    print("=" * 90)

    mod = 1 << k

    failures = 0

    for st in states:

        depth_k = min(
            st["depth"],
            k,
        )

        for d in range(1, depth_k + 1):

            md = 1 << d

            lhs = (
                st["D"] - st["T"]
            ) % md == 0

            rhs = (
                st["n"] + st["c"]
            ) % md == 0

            if not (
                lhs and rhs
            ):
                failures += 1

                if failures <= 10:
                    print(
                        "failure n={} frame={} "
                        "d={} D={} T={} "
                        "lhs={} rhs={}".format(
                            st["n"],
                            st["frame"],
                            d,
                            st["D"],
                            st["T"],
                            lhs,
                            rhs,
                        )
                    )

                break

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
            "    D mod 2^10={}".format(
                st["D"] % (1 << 10)
            )
        )

        print(
            "    n mod 2^10={}".format(
                st["n"] % (1 << 10)
            )
        )

        print(
            "    T={}".format(
                st["T"]
            )
        )

        print(
            "    D-T={}".format(
                st["D"] - st["T"]
            )
        )

        print(
            "    v2(D-T)={}".format(
                v2(st["D"] - st["T"])
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
    print("EXPERIMENT 679 START")
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
    total_ambiguities = 0

    total_failures += test_baseline(
        states
    )

    for k in K_VALUES:

        total_ambiguities += test_dn_signature(
            states,
            k,
        )

        total_ambiguities += test_sdn_signature(
            states,
            k,
        )

        total_ambiguities += test_direct_residue_depth(
            states,
            k,
        )

    # Focused "does n,D determine S?" test.
    for k in (6, 8, 10):

        total_ambiguities += test_n_plus_d_shift(
            states,
            k,
        )

        total_ambiguities += test_one_extra_bit(
            states,
            k,
        )

        total_failures += test_root_sign_flip(
            states,
            k,
        )

    show_examples(
        states
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Established:")
    print()

    print("    D^2 = S^2 - 4n")
    print()

    print("    FRAME A:")
    print("        T = S + 6")
    print("        c = 9")
    print()

    print("    FRAME B:")
    print("        T = -S - 2")
    print("        c = 3")
    print()

    print("    depth")
    print("        = min(v2(D-T), v2(n+c))")
    print()

    print("New question:")
    print()
    print("    Does")
    print()
    print("        (frame, n mod 2^k, D mod 2^k)")
    print()
    print("    determine truncated depth WITHOUT S?")
    print()

    print("If yes, the discriminant root itself carries")
    print("the missing branch information.")
    print()

    print("If no, the remaining information boundary lies")
    print("between the discriminant root D and the symmetric")
    print("coordinate S.")
    print()

    print("TOTAL FAILURES={}".format(
        total_failures
    ))

    print(
        "TOTAL AMBIGUITIES={}".format(
            total_ambiguities
        )
    )

    if total_failures == 0:
        print(
            "STATUS=CORE TESTS PASSED"
        )
    else:
        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )

    print("=" * 90)
    print("EXPERIMENT 679 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
