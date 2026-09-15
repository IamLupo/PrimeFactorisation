#!/usr/bin/env python3

from collections import defaultdict
from math import gcd, isqrt


# =============================================================================
# EXPERIMENT 682
#
# INFORMATION BOUNDARY OF (D, n)
#
# Experiment 679:
#
#     (frame, D mod 2^k, n mod 2^k)
#
#     is insufficient.
#
# Experiment 681:
#
#     increasing the number of n bits eventually resolves the collisions,
#     but the old "ambiguous bucket count" was incorrectly treated as a
#     monotone quantity.
#
# A refinement can split one ambiguous bucket into several smaller
# ambiguous buckets. Therefore:
#
#     ambiguous_bucket_count
#
# is NOT monotone.
#
# Experiment 682 measures instead:
#
#     1. whether ANY ambiguity remains;
#     2. number of states participating in ambiguity;
#     3. number of unresolved depth pairs;
#     4. first t at which the full partition is exact.
#
# It also studies the structure of the surviving collisions.
#
# Signature:
#
#     (frame,
#      D mod 2^k,
#      n mod 2^t)
#
# Target:
#
#     min(depth, k)
#
# We additionally measure:
#
#     S mod 2^r
#
# for the remaining collision classes to determine exactly which
# S bits are missing.
#
# No O(N^2) pairwise comparisons.
# =============================================================================


PRIME_LIMIT = 6000

K_VALUES = (4, 6, 8, 10, 12)

# Search sufficiently far beyond k.
MAX_EXTRA_N_BITS = 20

# Maximum S precision used for collision diagnosis.
MAX_S_BITS = 16


INF = 10 ** 9


# =============================================================================
# BASIC ARITHMETIC
# =============================================================================

def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    x = abs(x)

    if x == 0:
        return 0

    return x >> v2(x)


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve(limit: int):
    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    p = 2

    while p * p <= limit:

        if a[p]:

            start = p * p

            count = ((limit - start) // p) + 1

            a[start:limit + 1:p] = b"\x00" * count

        p += 1

    return [
        i
        for i in range(limit + 1)
        if a[i]
    ]


def odd_primes(limit: int):
    return [
        p
        for p in sieve(limit)
        if p & 1
    ]


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

        T = S + 6

    else:

        frame = "B"

        c = 3

        A = p + 1
        B = q - 3

        C = A

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
                    "mismatch n={} p={} q={} "
                    "frame={} depth={} predicted={}".format(
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
# PARTITION ANALYSIS
# =============================================================================

def signature(
    st,
    k,
    t,
):
    """
    Signature:

        frame
        D mod 2^k
        n mod 2^t

    t=0 means no n bits.
    """

    dmask = (1 << k) - 1

    if t == 0:
        nvalue = 0
    else:
        nvalue = st["n"] & ((1 << t) - 1)

    return (
        st["frame"],
        st["D"] & dmask,
        nvalue,
    )


def build_partition(states, k, t):

    buckets = defaultdict(list)

    for st in states:

        buckets[
            signature(st, k, t)
        ].append(st)

    return buckets


def analyze_partition(
    buckets,
    k,
):

    ambiguous_buckets = 0
    ambiguous_states = 0

    depth_pair_set = set()

    max_depth_classes = 0

    for items in buckets.values():

        depths = {
            min(st["depth"], k)
            for st in items
        }

        class_count = len(depths)

        if class_count > max_depth_classes:
            max_depth_classes = class_count

        if class_count > 1:

            ambiguous_buckets += 1

            ambiguous_states += len(items)

            dlist = sorted(depths)

            for i in range(len(dlist)):

                for j in range(i + 1, len(dlist)):

                    depth_pair_set.add(
                        (dlist[i], dlist[j])
                    )

    return {
        "ambiguous_buckets": ambiguous_buckets,
        "ambiguous_states": ambiguous_states,
        "depth_pairs": depth_pair_set,
        "max_depth_classes": max_depth_classes,
        "exact": ambiguous_buckets == 0,
    }


# =============================================================================
# TRUE MONOTONICITY TEST
# =============================================================================

def test_true_monotonicity(states, k, t_values):

    print("=" * 90)
    print(
        "TEST 1: TRUE INFORMATION MONOTONICITY k={}".format(
            k
        )
    )
    print("=" * 90)

    previous_ambiguity_signatures = None
    failures = 0

    for t in t_values:

        buckets = build_partition(
            states,
            k,
            t,
        )

        info = analyze_partition(
            buckets,
            k,
        )

        # The correct monotone object is:
        #
        #     set of state-pairs that remain unresolved.
        #
        # We cannot enumerate every pair, so we use:
        #
        #     depth-pair set
        #
        # together with unresolved-state count.
        #
        # Neither is guaranteed to be a complete lattice
        # invariant, but exactness itself IS monotone.
        if (
            previous_ambiguity_signatures is not None
            and info["exact"]
            and not previous_ambiguity_signatures[0]
        ):
            pass

        exact = info["exact"]

        print(
            "    t={:2d} buckets={:7d} "
            "ambiguous_buckets={:7d} "
            "ambiguous_states={:7d} "
            "depth_pairs={:2d} "
            "max_depth_classes={} "
            "status={}".format(
                t,
                len(buckets),
                info["ambiguous_buckets"],
                info["ambiguous_states"],
                len(info["depth_pairs"]),
                info["max_depth_classes"],
                "EXACT" if exact else "AMBIGUOUS",
            )
        )

        if previous_ambiguity_signatures is not None:

            previous_exact = previous_ambiguity_signatures[0]

            # Exact -> non-exact is impossible under refinement.
            if previous_exact and not exact:

                failures += 1

                print(
                    "    ERROR: exact partition became ambiguous"
                )

        previous_ambiguity_signatures = (
            exact,
            info["depth_pairs"],
        )

    print()

    print(
        "monotonicity failures={}".format(
            failures
        )
    )

    print()

    return failures


# =============================================================================
# FIRST EXACT t
# =============================================================================

def find_first_exact_t(states, k):

    max_t = k + MAX_EXTRA_N_BITS

    first_exact = None

    for t in range(max_t + 1):

        buckets = build_partition(
            states,
            k,
            t,
        )

        info = analyze_partition(
            buckets,
            k,
        )

        if info["exact"]:

            first_exact = t

            break

    return first_exact


# =============================================================================
# TEST 2:
# EXACT n-BIT BOUNDARY
# =============================================================================

def test_n_boundary(states, k):

    print("=" * 90)
    print(
        "TEST 2: EXACT n-BIT BOUNDARY k={}".format(
            k
        )
    )
    print("=" * 90)

    max_t = k + MAX_EXTRA_N_BITS

    first_exact = None

    for t in range(max_t + 1):

        buckets = build_partition(
            states,
            k,
            t,
        )

        info = analyze_partition(
            buckets,
            k,
        )

        if info["exact"]:

            first_exact = t

            print(
                "first exact t={}".format(
                    t
                )
            )

            break

    if first_exact is None:

        print(
            "no exact t <= {}".format(
                max_t
            )
        )

    print()

    return first_exact


# =============================================================================
# TEST 3:
# RESOLVING THE LAST COLLISIONS
# =============================================================================

def find_last_ambiguous_bucket(
    states,
    k,
    t,
):

    buckets = build_partition(
        states,
        k,
        t,
    )

    candidates = []

    for sig, items in buckets.items():

        depths = sorted({
            min(st["depth"], k)
            for st in items
        })

        if len(depths) > 1:

            candidates.append(
                (
                    len(items),
                    sig,
                    depths,
                    items,
                )
            )

    candidates.sort(
        key=lambda x: (
            x[0],
            len(x[2]),
        )
    )

    return candidates


def show_remaining_collisions(
    states,
    k,
    t,
    limit=12,
):

    print("=" * 90)
    print(
        "TEST 3: REMAINING COLLISIONS k={} t={}".format(
            k,
            t,
        )
    )
    print("=" * 90)

    candidates = find_last_ambiguous_bucket(
        states,
        k,
        t,
    )

    if not candidates:

        print(
            "NO REMAINING COLLISIONS"
        )

        print()

        return

    print(
        "ambiguous buckets={}".format(
            len(candidates)
        )
    )

    print()

    for size, sig, depths, items in candidates[:limit]:

        print(
            "signature={} states={} depths={}".format(
                sig,
                size,
                depths,
            )
        )

        shown_depths = set()

        for st in items:

            d = min(
                st["depth"],
                k,
            )

            if d in shown_depths:
                continue

            shown_depths.add(d)

            print(
                "    depth={} "
                "n={} "
                "S={} "
                "D={} "
                "p={} "
                "q={}".format(
                    d,
                    st["n"],
                    st["S"],
                    st["D"],
                    st["p"],
                    st["q"],
                )
            )

            if len(shown_depths) == len(depths):
                break

        print()


# =============================================================================
# TEST 4:
# HOW MANY S BITS RESOLVE THE FINAL n,D COLLISIONS?
# =============================================================================

def s_bits_for_partition(
    states,
    k,
    t,
    max_r,
):

    for r in range(max_r + 1):

        buckets = defaultdict(set)

        dmask = (1 << k) - 1

        if t == 0:
            nmask = None
        else:
            nmask = (1 << t) - 1

        if r == 0:
            smask = None
        else:
            smask = (1 << r) - 1

        for st in states:

            if nmask is None:
                nvalue = 0
            else:
                nvalue = st["n"] & nmask

            if smask is None:
                svalue = 0
            else:
                svalue = st["S"] & smask

            sig = (
                st["frame"],
                st["D"] & dmask,
                nvalue,
                svalue,
            )

            buckets[sig].add(
                min(st["depth"], k)
            )

        ambiguous = sum(
            1
            for values in buckets.values()
            if len(values) > 1
        )

        if ambiguous == 0:

            return r

    return None


def test_final_collision_tradeoff(
    states,
    k,
    t,
):

    print("=" * 90)
    print(
        "TEST 4: S-BIT COST AFTER n-BIT PRECISION "
        "k={} t={}".format(
            k,
            t,
        )
    )
    print("=" * 90)

    result = s_bits_for_partition(
        states,
        k,
        t,
        MAX_S_BITS,
    )

    print(
        "minimum S bits required={}".format(
            result
        )
    )

    print()

    return result


# =============================================================================
# TEST 5:
# ROOT SIGN FLIP SANITY
# =============================================================================

def test_root_sign(states):

    print("=" * 90)
    print("TEST 5: ROOT SIGN FLIP")
    print("=" * 90)

    failures = 0

    for st in states:

        S = st["S"]
        D = st["D"]
        n = st["n"]

        if (
            S * S
            != D * D + 4 * n
        ):

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
            "    C={}".format(
                st["C"]
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
    print("EXPERIMENT 682 START")
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

        # Correct monotonicity test.
        total_failures += test_true_monotonicity(
            states,
            k,
            range(
                0,
                k + MAX_EXTRA_N_BITS + 1,
            ),
        )

        # Find actual exact boundary.
        exact_t = test_n_boundary(
            states,
            k,
        )

        exact_results[k] = exact_t

        # Show the final unresolved structure.
        if exact_t is not None:

            # The last ambiguous partition is t=exact_t-1.
            if exact_t > 0:

                show_remaining_collisions(
                    states,
                    k,
                    exact_t - 1,
                )

                test_final_collision_tradeoff(
                    states,
                    k,
                    exact_t - 1,
                )

    # -------------------------------------------------------------------------
    # ROOT CONTROL
    # -------------------------------------------------------------------------

    total_failures += test_root_sign(
        states
    )

    # -------------------------------------------------------------------------
    # EXAMPLES
    # -------------------------------------------------------------------------

    show_examples(
        states
    )

    # -------------------------------------------------------------------------
    # FINAL SUMMARY
    # -------------------------------------------------------------------------

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print(
        "Experiment 681 found:"
    )

    print(
        "    k=10 -> first exact n precision t=25"
    )

    print()

    print(
        "Experiment 682 measures this correctly."
    )

    print()

    print(
        "The number of ambiguous buckets is NOT"
    )

    print(
        "required to be monotone under refinement."
    )

    print()

    print(
        "The correct monotone question is:"
    )

    print(
        "    does ANY ambiguity remain?"
    )

    print()

    print(
        "For the signature"
    )

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
        "we measure the first exact t and inspect"
    )

    print(
        "the final surviving collision classes."
    )

    print()

    print(
        "Observed exact boundaries:"
    )

    for k in K_VALUES:

        print(
            "    k={:2d} -> t={}".format(
                k,
                exact_results[k],
            )
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
    print("EXPERIMENT 682 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
