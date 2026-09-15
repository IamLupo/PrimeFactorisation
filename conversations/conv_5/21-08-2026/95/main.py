#!/usr/bin/env python3

from math import gcd


# =============================================================================
# CONFIG
# =============================================================================

PRIME_LIMIT = 6000
INF = 10**9


# =============================================================================
# BASIC
# =============================================================================

def sieve(limit: int):
    composite = bytearray(limit + 1)
    primes = []

    for p in range(2, limit + 1):
        if not composite[p]:
            primes.append(p)

            if p * p <= limit:
                composite[p * p:limit + 1:p] = b"\x01" * (
                    ((limit - p * p) // p) + 1
                )

    return primes


def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


# =============================================================================
# FRAME
# =============================================================================

def frame_of(n: int) -> str:
    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(f"Unexpected odd n={n}")


def params(frame: str):
    if frame == "A":
        return {
            "p0": 3,
            "q0": -3,
            "c": 9,
        }

    return {
        "p0": -1,
        "q0": 3,
        "c": 3,
    }


# =============================================================================
# STATE BUILD
# =============================================================================

def build_states(primes):
    odd_primes = [p for p in primes if p & 1]

    states = []

    for i, p in enumerate(odd_primes):
        for q in odd_primes[i:]:
            n = p * q

            frame = frame_of(n)
            cfg = params(frame)

            p0 = cfg["p0"]
            q0 = cfg["q0"]
            c = cfg["c"]

            A = p - p0
            B = q - q0

            if frame == "A":
                two_x = B - A
                two_y = A + B
            else:
                two_x = 3 * A - B
                two_y = 3 * A + B

            actual_depth = min(
                v2(two_x),
                v2(two_y),
            )

            tp = v2(A)
            tq = v2(B)
            w = v2(n + c)

            states.append(
                {
                    "p": p,
                    "q": q,
                    "n": n,
                    "frame": frame,
                    "p0": p0,
                    "q0": q0,
                    "c": c,
                    "A": A,
                    "B": B,
                    "tp": tp,
                    "tq": tq,
                    "w": w,
                    "depth": actual_depth,
                }
            )

    return states


# =============================================================================
# 3-ARG GCD
# =============================================================================

def gcd3(a: int, b: int, c: int) -> int:
    return gcd(gcd(abs(a), abs(b)), abs(c))


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states):
    print("=" * 90)
    print("TEST 0: PREVIOUS EXACT DEPTH LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = min(
            s["tp"] + 1,
            s["w"],
        )

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch "
                    f"n={s['n']} "
                    f"p={s['p']} "
                    f"q={s['q']} "
                    f"tp={s['tp']} "
                    f"tq={s['tq']} "
                    f"w={s['w']} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_p_gcd_identity(states):
    print("=" * 90)
    print("TEST 1: DEPTH = v2(gcd(2(p-p0), n+c))")
    print("=" * 90)

    failures = 0

    for s in states:
        G = gcd(
            2 * abs(s["p"] - s["p0"]),
            abs(s["n"] + s["c"]),
        )

        predicted = v2(G)

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch "
                    f"n={s['n']} "
                    f"gcd={G} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_q_gcd_identity(states):
    print("=" * 90)
    print("TEST 2: DEPTH = v2(gcd(2(q-q0), n+c))")
    print("=" * 90)

    failures = 0

    for s in states:
        G = gcd(
            2 * abs(s["q"] - s["q0"]),
            abs(s["n"] + s["c"]),
        )

        predicted = v2(G)

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch "
                    f"n={s['n']} "
                    f"gcd={G} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_symmetric_gcd_identity(states):
    print("=" * 90)
    print(
        "TEST 3: "
        "DEPTH = v2(gcd(2A, 2B, n+c))"
    )
    print("=" * 90)

    failures = 0

    for s in states:
        G = gcd3(
            2 * s["A"],
            2 * s["B"],
            s["n"] + s["c"],
        )

        predicted = v2(G)

        if predicted != s["depth"]:
            failures += 1

            if failures <= 20:
                print(
                    f"    mismatch "
                    f"n={s['n']} "
                    f"A={s['A']} "
                    f"B={s['B']} "
                    f"n+c={s['n'] + s['c']} "
                    f"G={G} "
                    f"predicted={predicted} "
                    f"actual={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_gcd_equalities(states):
    print("=" * 90)
    print("TEST 4: INTEGER GCD CROSS-EQUALITY")
    print("=" * 90)

    failures = 0

    for s in states:

        gp = gcd(
            2 * abs(s["p"] - s["p0"]),
            abs(s["n"] + s["c"]),
        )

        gq = gcd(
            2 * abs(s["q"] - s["q0"]),
            abs(s["n"] + s["c"]),
        )

        gs = gcd3(
            2 * s["A"],
            2 * s["B"],
            s["n"] + s["c"],
        )

        if v2(gp) != s["depth"]:
            failures += 1

        if v2(gq) != s["depth"]:
            failures += 1

        if v2(gs) != s["depth"]:
            failures += 1

        if not (
            v2(gp)
            == v2(gq)
            == v2(gs)
            == s["depth"]
        ):
            if failures <= 20:
                print(
                    f"    mismatch "
                    f"n={s['n']} "
                    f"gp={gp} "
                    f"gq={gq} "
                    f"gs={gs} "
                    f"depth={s['depth']}"
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_odd_part(states):
    print("=" * 90)
    print("TEST 5: ODD-PART BEHAVIOR OF THE SYMMETRIC GCD")
    print("=" * 90)

    """
    We are only claiming a 2-adic identity.

    Therefore the odd part of

        gcd(2A,2B,n+c)

    is allowed to vary.

    This test measures that odd part rather than incorrectly
    asserting an integer equality.
    """

    distribution = {}

    for s in states:
        G = gcd3(
            2 * s["A"],
            2 * s["B"],
            s["n"] + s["c"],
        )

        odd = abs(G) >> v2(G)

        key = (s["frame"], v2(G))

        distribution[key] = distribution.get(key, set())
        distribution[key].add(odd)

    print("    distinct odd parts per frame/depth:")
    for key in sorted(distribution):
        values = distribution[key]

        print(
            f"        {key}: "
            f"count={len(values)} "
            f"sample={sorted(values)[:10]}"
        )

    print()
    return 0


# =============================================================================
# TEST 6
# =============================================================================

def test_branch_reduction(states):
    print("=" * 90)
    print("TEST 6: BRANCH VALUATION COLLAPSE")
    print("=" * 90)

    failures = 0

    for s in states:

        # Original:
        original = min(
            s["tp"] + 1,
            s["w"],
        )

        # Symmetric gcd:
        G = gcd3(
            2 * s["A"],
            2 * s["B"],
            s["n"] + s["c"],
        )

        collapsed = v2(G)

        if original != collapsed:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def examples(states):

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

    by_n = {s["n"]: s for s in states}

    for n in wanted:

        s = by_n.get(n)

        if s is None:
            continue

        Gp = gcd(
            2 * abs(s["p"] - s["p0"]),
            abs(s["n"] + s["c"]),
        )

        Gq = gcd(
            2 * abs(s["q"] - s["q0"]),
            abs(s["n"] + s["c"]),
        )

        Gs = gcd3(
            2 * s["A"],
            2 * s["B"],
            s["n"] + s["c"],
        )

        print()
        print(
            f"n={s['n']} "
            f"p={s['p']} "
            f"q={s['q']} "
            f"frame={s['frame']}"
        )

        print(
            f"    A={s['A']} "
            f"B={s['B']} "
            f"p0={s['p0']} "
            f"q0={s['q0']} "
            f"c={s['c']}"
        )

        print(
            f"    tp={s['tp']} "
            f"tq={s['tq']} "
            f"w={s['w']} "
            f"depth={s['depth']}"
        )

        print(
            f"    2(p-p0)={2*(s['p']-s['p0'])}"
        )

        print(
            f"    2(q-q0)={2*(s['q']-s['q0'])}"
        )

        print(
            f"    n+c={s['n'] + s['c']}"
        )

        print(
            f"    gcd_p={Gp} "
            f"v2={v2(Gp)}"
        )

        print(
            f"    gcd_q={Gq} "
            f"v2={v2(Gq)}"
        )

        print(
            f"    gcd_sym={Gs} "
            f"v2={v2(Gs)}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 658 START")
    print("=" * 90)

    print()
    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve(PRIME_LIMIT)
    odd_primes = [p for p in primes if p & 1]

    print(f"odd primes={len(odd_primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    failures = 0

    failures += test_baseline(states)
    failures += test_p_gcd_identity(states)
    failures += test_q_gcd_identity(states)
    failures += test_symmetric_gcd_identity(states)
    failures += test_gcd_equalities(states)
    test_odd_part(states)
    failures += test_branch_reduction(states)

    examples(states)

    print()
    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
The previous exact law is:

    depth =
        min(
            v2(p-p0)+1,
            v2(n+c)
        ).

Because:

    v2(gcd(x,y))
        =
    min(v2(x),v2(y)),

we obtain immediately:

    depth =
        v2(
            gcd(
                2(p-p0),
                n+c
            )
        ).

The same must hold from q:

    depth =
        v2(
            gcd(
                2(q-q0),
                n+c
            )
        ).

Since:

    A = p-p0
    B = q-q0

and:

    m = min(v2(A),v2(B)),

the symmetric form is:

    depth =
        v2(
            gcd(
                2A,
                2B,
                n+c
            )
        ).

This is the new object being tested.

The important question is whether the entire hierarchy can
therefore be expressed as the 2-adic valuation of one
integer gcd, instead of as a first-differing-bit process.

Frame A:

    A = p-3
    B = q+3
    c = 9

Frame B:

    A = p+1
    B = q-3
    c = 3.

The expected unified identity is:

    depth =
        v2(
            gcd(
                2A,
                2B,
                n+c
            )
        ).

"""
    )

    print(f"TOTAL FAILURES={failures}")

    if failures == 0:
        print("STATUS=ALL TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 658 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
