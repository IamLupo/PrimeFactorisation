#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 629
# ==============================================================================
#
# MINIMAL GLOBAL 2-ADIC STATE AFTER PRIMARY n-VALUATION
#
# Experiment 628 established:
#
#     FRAME A:
#         primary = v2(n+9)
#
#     FRAME B:
#         primary = v2(n+3)
#
# and, crucially:
#
#     (primary, v2(X), v2(Y))
#
# determines depth exactly.
#
# This experiment asks whether one of v2(X), v2(Y) is redundant.
#
# Tested signatures:
#
#     primary
#     primary + v2X
#     primary + v2Y
#     primary + min(v2X,v2Y)
#     primary + max(v2X,v2Y)
#     primary + |v2X-v2Y|
#     primary + equality(v2X,v2Y)
#     primary + min + equality
#     primary + v2X + v2Y
#
# We also test frame-conditioned variants because A/B may permit
# different reduced descriptions.
#
# No global offset search is performed.
#
# ==============================================================================


PRIME_LIMIT = 6250
INF = 10**9


# ==============================================================================
# DATA
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


# ==============================================================================
# v2
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


# ==============================================================================
# SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):

        if not sieve[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        sieve[start : limit + 1 : p] = (
            b"\x00"
        ) * count

    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIMES
# ==============================================================================

def generate_states(
    primes: list[int],
) -> list[State]:

    states = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            states.append(
                State(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return states


# ==============================================================================
# FRAME
# ==============================================================================

def frame_of(n: int) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"unexpected odd n mod 4={r}"
    )


# ==============================================================================
# RAW FACTOR RESIDUALS
# ==============================================================================

def raw_residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        return (
            p - 3,
            q + 3,
        )

    return (
        p + 1,
        q - 3,
    )


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        X = (
            q - p + 6
        ) // 2

        Y = (
            p + q
        ) // 2

    else:

        X = (
            3 * p - q + 6
        ) // 2

        Y = (
            3 * p + q
        ) // 2

    return X, Y


# ==============================================================================
# EXACT DEPTH
# ==============================================================================

def exact_depth(
    state: State,
) -> int:

    frame = frame_of(state.n)

    A, B = raw_residuals(
        frame,
        state.p,
        state.q,
    )

    alpha = v2(A)
    beta = v2(B)

    return (
        min(alpha, beta)
        + int(alpha == beta)
    )


# ==============================================================================
# PRIMARY VALUATION
# ==============================================================================

def primary_value(
    state: State,
    frame: str,
) -> int:

    if frame == "A":

        return v2(state.n + 9)

    return v2(state.n + 3)


# ==============================================================================
# BUILD RECORDS
# ==============================================================================

def build_records(
    states: list[State],
):

    records = []

    for i, state in enumerate(states):

        frame = frame_of(state.n)

        primary = primary_value(
            state,
            frame,
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        vx = v2(X)
        vy = v2(Y)

        depth = exact_depth(state)

        records.append(
            (
                i,
                frame,
                primary,
                vx,
                vy,
                depth,
            )
        )

    return records


# ==============================================================================
# PRIMARY COLLISION SETS
# ==============================================================================

def find_ambiguous_primary_buckets(
    records,
):

    buckets = defaultdict(list)

    for record in records:

        (
            i,
            frame,
            primary,
            vx,
            vy,
            depth,
        ) = record

        buckets[
            (
                frame,
                primary,
            )
        ].append(record)

    ambiguous = {}

    for key, rows in buckets.items():

        depths = {
            row[5]
            for row in rows
        }

        if len(depths) > 1:

            ambiguous[key] = rows

    return ambiguous


# ==============================================================================
# SIGNATURE CHECK
# ==============================================================================

def signature_exact(
    rows,
    signature_fn,
):

    mapping = {}

    for row in rows:

        sig = signature_fn(row)
        depth = row[5]

        old = mapping.get(sig)

        if old is None:

            mapping[sig] = depth

        elif old != depth:

            return False

    return True


def total_ambiguous(
    rows,
    signature_fn,
):

    mapping = {}
    bad = set()

    for row in rows:

        sig = signature_fn(row)
        depth = row[5]

        old = mapping.get(sig)

        if old is None:

            mapping[sig] = depth

        elif old != depth:

            bad.add(sig)

    return len(bad)


# ==============================================================================
# TEST SIGNATURE FAMILY
# ==============================================================================

def test_signatures(
    ambiguous,
):

    print()
    print("=" * 90)
    print("TEST 1: GLOBAL SIGNATURE REDUCTION")
    print("=" * 90)

    candidates = {

        "primary":
            lambda r: (
                r[2],
            ),

        "primary + v2X":
            lambda r: (
                r[2],
                r[3],
            ),

        "primary + v2Y":
            lambda r: (
                r[2],
                r[4],
            ),

        "primary + min":
            lambda r: (
                r[2],
                min(r[3], r[4]),
            ),

        "primary + max":
            lambda r: (
                r[2],
                max(r[3], r[4]),
            ),

        "primary + |diff|":
            lambda r: (
                r[2],
                abs(r[3] - r[4]),
            ),

        "primary + equality":
            lambda r: (
                r[2],
                int(r[3] == r[4]),
            ),

        "primary + min + equality":
            lambda r: (
                r[2],
                min(r[3], r[4]),
                int(r[3] == r[4]),
            ),

        "primary + v2X + v2Y":
            lambda r: (
                r[2],
                r[3],
                r[4],
            ),
    }

    exact = []

    for name, fn in candidates.items():

        ambiguous_count = 0

        for rows in ambiguous.values():

            ambiguous_count += (
                total_ambiguous(
                    rows,
                    fn,
                )
            )

        status = (
            "EXACT"
            if ambiguous_count == 0
            else "not exact"
        )

        print(
            f"    {name:<28}"
            f" ambiguous={ambiguous_count:<6}"
            f" {status}"
        )

        if ambiguous_count == 0:

            exact.append(name)

    return exact


# ==============================================================================
# FRAME-CONDITIONED TEST
# ==============================================================================

def test_frame_conditioned(
    ambiguous,
):

    print()
    print("=" * 90)
    print("TEST 2: FRAME-CONDITIONED REDUCTION")
    print("=" * 90)

    candidates = {

        "primary + v2X":
            lambda r: (
                r[2],
                r[3],
            ),

        "primary + v2Y":
            lambda r: (
                r[2],
                r[4],
            ),

        "primary + min":
            lambda r: (
                r[2],
                min(r[3], r[4]),
            ),

        "primary + equality":
            lambda r: (
                r[2],
                int(r[3] == r[4]),
            ),

        "primary + min + equality":
            lambda r: (
                r[2],
                min(r[3], r[4]),
                int(r[3] == r[4]),
            ),

        "primary + v2X + v2Y":
            lambda r: (
                r[2],
                r[3],
                r[4],
            ),
    }

    for frame in ("A", "B"):

        print()
        print(
            f"FRAME {frame}"
        )

        frame_rows = [
            rows
            for (
                f,
                _,
            ), rows in ambiguous.items()
            if f == frame
        ]

        for name, fn in candidates.items():

            ambiguous_count = 0

            for rows in frame_rows:

                ambiguous_count += (
                    total_ambiguous(
                        rows,
                        fn,
                    )
                )

            print(
                f"    {name:<28}"
                f" ambiguous={ambiguous_count}"
            )


# ==============================================================================
# RELATION PRIMARY -> min(v2X,v2Y)
# ==============================================================================

def test_primary_vs_min(
    records,
):

    print()
    print("=" * 90)
    print("TEST 3: PRIMARY VS GLOBAL MIN-VALUATION")
    print("=" * 90)

    mapping = defaultdict(set)

    for (
        i,
        frame,
        primary,
        vx,
        vy,
        depth,
    ) in records:

        m = min(vx, vy)

        mapping[
            (
                frame,
                primary,
            )
        ].add(m)

    ambiguous = []

    for key, values in mapping.items():

        if len(values) > 1:

            ambiguous.append(
                (
                    key,
                    sorted(values),
                )
            )

    print(
        f"    ambiguous primary->min "
        f"buckets={len(ambiguous)}"
    )

    for key, values in ambiguous[:20]:

        print(
            f"        {key}"
            f" -> {values}"
        )


# ==============================================================================
# RELATION PRIMARY + ONE VALUATION
# ==============================================================================

def test_one_valuation_formula(
    records,
):

    print()
    print("=" * 90)
    print("TEST 4: PRIMARY + ONE GLOBAL VALUATION")
    print("=" * 90)

    tests = {

        "primary,v2X":
            lambda primary, vx, vy: (
                primary,
                vx,
            ),

        "primary,v2Y":
            lambda primary, vx, vy: (
                primary,
                vy,
            ),

        "primary,min":
            lambda primary, vx, vy: (
                primary,
                min(vx, vy),
            ),

        "primary,max":
            lambda primary, vx, vy: (
                primary,
                max(vx, vy),
            ),

        "primary,diff":
            lambda primary, vx, vy: (
                primary,
                abs(vx - vy),
            ),

    }

    for name, fn in tests.items():

        mapping = {}
        ambiguous = 0

        for (
            i,
            frame,
            primary,
            vx,
            vy,
            depth,
        ) in records:

            sig = fn(
                primary,
                vx,
                vy,
            )

            old = mapping.get(sig)

            if old is None:

                mapping[sig] = depth

            elif old != depth:

                ambiguous += 1

        print(
            f"    {name:<20}"
            f" collisions={ambiguous}"
        )


# ==============================================================================
# LOOK AT EXACT COLLISION
# ==============================================================================

def print_key_collisions(
    ambiguous,
    states,
):

    print()
    print("=" * 90)
    print("TEST 5: KEY COLLISION EXAMPLES")
    print("=" * 90)

    shown = 0

    for key, rows in sorted(
        ambiguous.items()
    ):

        if shown >= 12:

            break

        frame, primary = key

        print()
        print(
            f"    FRAME={frame}"
            f" primary={primary}"
        )

        signatures = defaultdict(
            list
        )

        for row in rows:

            (
                i,
                frame,
                primary,
                vx,
                vy,
                depth,
            ) = row

            signatures[
                (
                    vx,
                    vy,
                )
            ].append(
                (
                    states[i],
                    depth,
                )
            )

        for sig, entries in list(
            signatures.items()
        )[:8]:

            depth_set = sorted({
                depth
                for _, depth in entries
            })

            if len(depth_set) < 2:
                continue

            print(
                f"        "
                f"(v2X,v2Y)={sig}"
                f" -> depths={depth_set}"
            )

            for state, depth in entries[:5]:

                print(
                    f"            "
                    f"n={state.n} "
                    f"p={state.p} "
                    f"q={state.q} "
                    f"depth={depth}"
                )

        shown += 1


# ==============================================================================
# DIRECT DEPTH IDENTITY
# ==============================================================================

def verify_global_identity(
    records,
):

    print()
    print("=" * 90)
    print("TEST 6: GLOBAL DEPTH IDENTITY")
    print("=" * 90)

    failures = 0

    for (
        i,
        frame,
        primary,
        vx,
        vy,
        depth,
    ) in records:

        predicted = (
            1 + min(vx, vy)
        )

        if predicted != depth:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" index={i}"
                    f" frame={frame}"
                    f" primary={primary}"
                    f" v2X={vx}"
                    f" v2Y={vy}"
                    f" depth={depth}"
                    f" predicted={predicted}"
                )

    print(
        f"checked={len(records)} "
        f"failures={failures}"
    )


# ==============================================================================
# FINAL COLLISION REDUCTION
# ==============================================================================

def final_reduction(
    ambiguous,
):

    print()
    print("=" * 90)
    print("TEST 7: MINIMAL EXACT STATE")
    print("=" * 90)

    candidates = [

        (
            "primary + v2X",
            lambda r: (
                r[2],
                r[3],
            ),
        ),

        (
            "primary + v2Y",
            lambda r: (
                r[2],
                r[4],
            ),
        ),

        (
            "primary + min",
            lambda r: (
                r[2],
                min(r[3], r[4]),
            ),
        ),

        (
            "primary + min + equality",
            lambda r: (
                r[2],
                min(r[3], r[4]),
                int(r[3] == r[4]),
            ),
        ),

        (
            "primary + v2X + v2Y",
            lambda r: (
                r[2],
                r[3],
                r[4],
            ),
        ),
    ]

    exact = []

    for name, fn in candidates:

        total = 0

        for rows in ambiguous.values():

            total += total_ambiguous(
                rows,
                fn,
            )

        if total == 0:

            exact.append(name)

        print(
            f"    {name:<32}"
            f" ambiguous={total}"
        )

    return exact


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 629 START")
    print("=" * 90)

    print()
    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )

    print()
    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )

    print()
    print("[3] BUILD GLOBAL RECORDS")

    records = build_records(
        states
    )

    print(
        f"    states={len(records)}"
    )

    print()
    print("[4] BUILD PRIMARY COLLISION BUCKETS")

    ambiguous = (
        find_ambiguous_primary_buckets(
            records
        )
    )

    print(
        f"    ambiguous buckets="
        f"{len(ambiguous)}"
    )

    collision_states = sum(
        len(rows)
        for rows in ambiguous.values()
    )

    print(
        f"    collision states="
        f"{collision_states}"
    )

    # --------------------------------------------------------------------------
    # Main reduction
    # --------------------------------------------------------------------------

    exact = test_signatures(
        ambiguous
    )

    test_frame_conditioned(
        ambiguous
    )

    test_primary_vs_min(
        records
    )

    test_one_valuation_formula(
        records
    )

    print_key_collisions(
        ambiguous,
        states,
    )

    verify_global_identity(
        records
    )

    minimal = final_reduction(
        ambiguous
    )

    # --------------------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print("FINAL SUMMARY")
    print("=" * 90)

    print()
    print(
        "Experiment 628 established that:"
    )

    print()
    print(
        "    (primary, v2X, v2Y)"
        " -> depth"
    )

    print(
        "is exact."
    )

    print()
    print(
        "Experiment 629 asks whether this"
        " can be reduced further."
    )

    print()
    print(
        "EXACT GLOBAL SIGNATURES:"
    )

    for name in exact:

        print(
            f"    {name}"
        )

    print()

    print(
        "MINIMAL EXACT CANDIDATES:"
    )

    for name in minimal:

        print(
            f"    {name}"
        )

    print()
    print(
        "If only:"
    )

    print(
        "    primary + min(v2X,v2Y)"
    )

    print(
        "is exact, then the equality branch"
        " and the larger valuation are redundant."
    )

    print()
    print(
        "If:"
    )

    print(
        "    primary + v2X"
        " or"
        " primary + v2Y"
    )

    print(
        "is exact, then one entire global"
        " coordinate valuation disappears."
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 629 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
