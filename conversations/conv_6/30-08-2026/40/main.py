#!/usr/bin/env python3

import math
import random
import statistics
import time
from dataclasses import dataclass

import z3


# ================================================================================================
# CONFIGURATION
# ================================================================================================

SEED = 1_511_464_998

FACTOR_LO = 10_000
FACTOR_HI = 100_000

MOD_MIN = 300
MOD_MAX = 3_000
CLOSE_RATIO = 0.20

SEARCH_ANCHORS = 20
TIMEOUT_MS = 5_000


# ================================================================================================
# DATA
# ================================================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int
    n: int


# ================================================================================================
# PRIME GENERATION
# ================================================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0] = 0
    a[1] = 0

    for p in range(2, math.isqrt(limit) + 1):
        if a[p]:
            start = p * p
            a[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, v in enumerate(a) if v]


def build_prime_pool(lo: int, hi: int) -> list[int]:
    return [p for p in sieve(hi) if p >= lo]


# ================================================================================================
# ANCHORS
# ================================================================================================

def build_anchors(
    primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    result = []
    seen = set()

    while len(result) < count:
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        result.append(
            Anchor(
                p=p,
                q=q,
                n=p * q,
            )
        )

    return result


# ================================================================================================
# CLOSE MODULUS TRIPLES
# ================================================================================================

def build_close_triples(
    primes: list[int],
    count: int,
    seed: int,
) -> list[tuple[int, int, int]]:

    rng = random.Random(seed)

    candidates = []

    for i, r1 in enumerate(primes):

        max_r = int(r1 * (1.0 + CLOSE_RATIO))

        close = [
            r
            for r in primes[i + 1:]
            if r <= max_r
        ]

        # Keep construction bounded.
        close = close[:12]

        for j in range(len(close)):
            for k in range(j + 1, len(close)):

                r2 = close[j]
                r3 = close[k]

                if (
                    math.gcd(r1, r2) == 1
                    and math.gcd(r1, r3) == 1
                    and math.gcd(r2, r3) == 1
                ):
                    candidates.append(
                        (r1, r2, r3)
                    )

    if not candidates:
        raise RuntimeError("No valid modulus triples.")

    rng.shuffle(candidates)

    result = []
    seen = set()

    for triple in candidates:
        if triple in seen:
            continue

        seen.add(triple)
        result.append(triple)

        if len(result) >= count:
            break

    if len(result) < count:
        raise RuntimeError(
            f"Only {len(result)} unique triples available."
        )

    return result


# ================================================================================================
# MODEL 0
#
# Pure nonlinear integer factorization.
# ================================================================================================

def build_plain_model(
    n: int,
    lo: int,
    hi: int,
):

    p = z3.Int("p")
    q = z3.Int("q")

    s = z3.Solver()
    s.set(timeout=TIMEOUT_MS)

    s.add(
        p >= lo,
        p <= hi,
        q >= lo,
        q <= hi,
        p <= q,
        p * q == n,
    )

    return s, p, q


# ================================================================================================
# MODEL 1
#
# Unknown residues:
#
#     a_i = p mod r_i
#     b_i = q mod r_i
#
# but the residues are NOT supplied.
#
# We explicitly connect them through:
#
#     p = a_i + k_i r_i
#     q = b_i + l_i r_i
#
# This is still mathematically equivalent to ordinary integer
# factorization, but tests whether the coordinate representation
# changes solver behaviour.
# ================================================================================================

def build_coordinate_model(
    n: int,
    lo: int,
    hi: int,
    mods: tuple[int, int, int],
):

    p = z3.Int("p")
    q = z3.Int("q")

    s = z3.Solver()
    s.set(timeout=TIMEOUT_MS)

    s.add(
        p >= lo,
        p <= hi,
        q >= lo,
        q <= hi,
        p <= q,
        p * q == n,
    )

    for i, r in enumerate(mods):

        a = z3.Int(f"a_{i}")
        b = z3.Int(f"b_{i}")

        k = z3.Int(f"k_{i}")
        ell = z3.Int(f"l_{i}")

        s.add(
            0 <= a,
            a < r,
            0 <= b,
            b < r,
            k >= 0,
            ell >= 0,
        )

        s.add(
            p == a + k * r,
            q == b + ell * r,
        )

    return s, p, q


# ================================================================================================
# MODEL 2
#
# Unknown residues, but explicitly enforce the multiplicative
# residue relation:
#
#     a_i b_i == n (mod r_i)
#
# IMPORTANT:
#
# z3's Int expressions support "%" directly.
#
# Do NOT use z3.Mod(...).
# ================================================================================================

def build_residue_model(
    n: int,
    lo: int,
    hi: int,
    mods: tuple[int, int, int],
):

    p = z3.Int("p")
    q = z3.Int("q")

    s = z3.Solver()
    s.set(timeout=TIMEOUT_MS)

    s.add(
        p >= lo,
        p <= hi,
        q >= lo,
        q <= hi,
        p <= q,
        p * q == n,
    )

    for i, r in enumerate(mods):

        a = z3.Int(f"ra_{i}")
        b = z3.Int(f"rb_{i}")

        k = z3.Int(f"rk_{i}")
        ell = z3.Int(f"rl_{i}")

        s.add(
            0 <= a,
            a < r,
            0 <= b,
            b < r,
            k >= 0,
            ell >= 0,
        )

        s.add(
            p == a + k * r,
            q == b + ell * r,
        )

        # Correct Z3Py syntax.
        s.add(
            (a * b - n) % r == 0
        )

    return s, p, q


# ================================================================================================
# MODEL 3
#
# Unknown residues plus CLOSE-MODULUS difference coordinates.
#
# For:
#
#     r2 = r1 + d2
#     r3 = r1 + d3
#
# the same integer p is represented as:
#
#     p = a1 + k1*r1
#       = a2 + k2*r2
#       = a3 + k3*r3
#
# We expose the relationships between the quotient coordinates
# instead of merely adding independent decompositions.
# ================================================================================================

def build_close_coordinate_model(
    n: int,
    lo: int,
    hi: int,
    mods: tuple[int, int, int],
):

    r1, r2, r3 = mods

    d2 = r2 - r1
    d3 = r3 - r1

    p = z3.Int("p")
    q = z3.Int("q")

    a1 = z3.Int("a1")
    a2 = z3.Int("a2")
    a3 = z3.Int("a3")

    b1 = z3.Int("b1")
    b2 = z3.Int("b2")
    b3 = z3.Int("b3")

    k1 = z3.Int("k1")
    k2 = z3.Int("k2")
    k3 = z3.Int("k3")

    l1 = z3.Int("l1")
    l2 = z3.Int("l2")
    l3 = z3.Int("l3")

    s = z3.Solver()
    s.set(timeout=TIMEOUT_MS)

    s.add(
        p >= lo,
        p <= hi,
        q >= lo,
        q <= hi,
        p <= q,
        p * q == n,
    )

    residues = [
        (a1, b1, k1, l1, r1),
        (a2, b2, k2, l2, r2),
        (a3, b3, k3, l3, r3),
    ]

    for a, b, k, ell, r in residues:

        s.add(
            0 <= a,
            a < r,
            0 <= b,
            b < r,
            k >= 0,
            ell >= 0,
        )

        s.add(
            p == a + k * r,
            q == b + ell * r,
        )

    # Same p represented at r1 and r2:
    #
    # a1 + k1*r1 = a2 + k2*r2
    #
    # which exposes the close-modulus relation.
    s.add(
        a1 + k1 * r1
        ==
        a2 + k2 * r2
    )

    s.add(
        a1 + k1 * r1
        ==
        a3 + k3 * r3
    )

    # Same q.
    s.add(
        b1 + l1 * r1
        ==
        b2 + l2 * r2
    )

    s.add(
        b1 + l1 * r1
        ==
        b3 + l3 * r3
    )

    # Explicit difference identities.
    #
    # a2-a1 = k1*r1-k2*r2
    #         = (k1-k2)*r1-k2*d2
    #
    # These are equivalent to the representation equalities,
    # but help expose the special role of the small gaps d2,d3.
    s.add(
        a2 - a1
        ==
        (k1 - k2) * r1 - k2 * d2
    )

    s.add(
        a3 - a1
        ==
        (k1 - k3) * r1 - k3 * d3
    )

    s.add(
        b2 - b1
        ==
        (l1 - l2) * r1 - l2 * d2
    )

    s.add(
        b3 - b1
        ==
        (l1 - l3) * r1 - l3 * d3
    )

    return s, p, q


# ================================================================================================
# RUN
# ================================================================================================

def execute(
    builder,
    n: int,
    true_p: int,
    true_q: int,
    lo: int,
    hi: int,
    mods,
):

    start = time.perf_counter()

    try:
        solver, p_var, q_var = builder(
            n,
            lo,
            hi,
            mods,
        )
    except Exception as exc:
        return {
            "status": "builder-error",
            "time_ms": 0.0,
            "p": None,
            "q": None,
            "exact": False,
            "error": repr(exc),
        }

    result = solver.check()

    elapsed = (
        time.perf_counter() - start
    ) * 1000.0

    recovered_p = None
    recovered_q = None
    exact = False

    if result == z3.sat:

        model = solver.model()

        recovered_p = model.eval(
            p_var
        ).as_long()

        recovered_q = model.eval(
            q_var
        ).as_long()

        exact = (
            recovered_p * recovered_q == n
            and recovered_p in (true_p, true_q)
            and recovered_q in (true_p, true_q)
            and recovered_p != recovered_q
        )

    return {
        "status": str(result),
        "time_ms": elapsed,
        "p": recovered_p,
        "q": recovered_q,
        "exact": exact,
        "error": None,
    }


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    total_start = time.perf_counter()

    print("=" * 100)
    print(
        "THREE-CLOSE-PRIME SMT UNKNOWN-RESIDUE / "
        "CLOSE-COORDINATE EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"factor range             = "
        f"{FACTOR_LO:,} - {FACTOR_HI:,}"
    )
    print(
        f"modulus prime range      = "
        f"{MOD_MIN:,} - {MOD_MAX:,}"
    )
    print(
        f"close ratio              = "
        f"{CLOSE_RATIO:.0%}"
    )
    print(
        f"search anchors            = "
        f"{SEARCH_ANCHORS}"
    )
    print(
        f"SMT timeout/model         = "
        f"{TIMEOUT_MS:,} ms"
    )

    modulus_primes = build_prime_pool(
        MOD_MIN,
        MOD_MAX,
    )

    factor_primes = build_prime_pool(
        FACTOR_LO,
        FACTOR_HI,
    )

    anchors = build_anchors(
        factor_primes,
        SEARCH_ANCHORS,
        SEED,
    )

    triples = build_close_triples(
        modulus_primes,
        SEARCH_ANCHORS,
        SEED ^ 0x123456,
    )

    print()
    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    print(
        f"unique modulus triples    = "
        f"{len(set(triples)):,}"
    )

    # --------------------------------------------------------------------------------------------
    # VALIDATION
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("VALIDATING THREE-MODULUS IDENTITIES")
    print("=" * 100)

    failures = 0

    for anchor, mods in zip(anchors, triples):

        for r in mods:

            if anchor.n % r != (
                anchor.p * anchor.q
            ) % r:

                failures += 1

    print(
        f"identity failures         = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Modular identity validation failed."
        )

    # --------------------------------------------------------------------------------------------
    # MODELS
    # --------------------------------------------------------------------------------------------

    builders = {
        0: None,
        1: build_coordinate_model,
        2: build_residue_model,
        3: build_close_coordinate_model,
    }

    results = {
        0: [],
        1: [],
        2: [],
        3: [],
    }

    names = {
        0: "SMT0  p*q=n",
        1: "SMT1  unknown residues",
        2: "SMT2  unknown residues + product congruence",
        3: "SMT3  unknown residues + close-coordinate relations",
    }

    # --------------------------------------------------------------------------------------------
    # SEARCH
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RUNNING SMT SEARCH")
    print("=" * 100)

    for idx, (anchor, mods) in enumerate(
        zip(anchors, triples)
    ):

        print(
            f"anchor {idx + 1:2d}/{SEARCH_ANCHORS}",
            end="\r",
            flush=True,
        )

        r0 = execute(
            build_plain_model,
            anchor.n,
            anchor.p,
            anchor.q,
            FACTOR_LO,
            FACTOR_HI,
            mods,
        )

        r1 = execute(
            build_coordinate_model,
            anchor.n,
            anchor.p,
            anchor.q,
            FACTOR_LO,
            FACTOR_HI,
            mods,
        )

        r2 = execute(
            build_residue_model,
            anchor.n,
            anchor.p,
            anchor.q,
            FACTOR_LO,
            FACTOR_HI,
            mods,
        )

        r3 = execute(
            build_close_coordinate_model,
            anchor.n,
            anchor.p,
            anchor.q,
            FACTOR_LO,
            FACTOR_HI,
            mods,
        )

        results[0].append(r0)
        results[1].append(r1)
        results[2].append(r2)
        results[3].append(r3)

    print(" " * 100, end="\r")

    # --------------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"factor-prime baseline      = "
        f"{len(factor_primes):,}"
    )

    print()

    for model_id in range(4):

        rs = results[model_id]

        sat = sum(
            x["status"] == "sat"
            for x in rs
        )

        unknown = sum(
            x["status"] == "unknown"
            for x in rs
        )

        errors = sum(
            x["status"] == "builder-error"
            for x in rs
        )

        exact = sum(
            x["exact"]
            for x in rs
        )

        times = [
            x["time_ms"]
            for x in rs
            if x["status"] in ("sat", "unsat")
        ]

        print(names[model_id])

        print(
            f"    sat                      = "
            f"{sat}/{SEARCH_ANCHORS}"
        )

        print(
            f"    unknown/timeout          = "
            f"{unknown}/{SEARCH_ANCHORS}"
        )

        print(
            f"    builder errors           = "
            f"{errors}"
        )

        print(
            f"    exact recovery           = "
            f"{exact}/{SEARCH_ANCHORS}"
        )

        if times:
            print(
                f"    mean solve time          = "
                f"{statistics.mean(times):.3f} ms"
            )

            print(
                f"    median solve time        = "
                f"{statistics.median(times):.3f} ms"
            )

            print(
                f"    maximum solve time       = "
                f"{max(times):.3f} ms"
            )

        print()

    # --------------------------------------------------------------------------------------------
    # SPEED RATIOS
    # --------------------------------------------------------------------------------------------

    means = {}

    for model_id in range(4):

        ts = [
            x["time_ms"]
            for x in results[model_id]
            if x["status"] in ("sat", "unsat")
        ]

        means[model_id] = (
            statistics.mean(ts)
            if ts
            else float("inf")
        )

    print("=" * 100)
    print("SMT SPEED COMPARISON")
    print("=" * 100)

    print(
        f"SMT0 mean                 = "
        f"{means[0]:,.3f} ms"
    )

    for model_id in (1, 2, 3):

        ratio = (
            means[model_id] / means[0]
            if math.isfinite(means[model_id])
            else float("inf")
        )

        print(
            f"SMT{model_id}/SMT0            = "
            f"{ratio:.4f}"
        )

    # --------------------------------------------------------------------------------------------
    # ANCHOR RESULTS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID       p       q    r1   r2   r3"
        "    SMT0    SMT1    SMT2    SMT3"
    )

    print("-" * 100)

    for i, anchor in enumerate(anchors):

        r1, r2, r3 = triples[i]

        vals = []

        for model_id in range(4):

            r = results[model_id][i]

            if r["status"] == "unknown":
                vals.append("TIMEOUT")

            elif r["status"] == "builder-error":
                vals.append("ERROR")

            else:
                vals.append(
                    f"{r['time_ms']:7.1f}"
                )

        print(
            f"{i+1:3d} "
            f"{anchor.p:8,d} "
            f"{anchor.q:8,d} "
            f"{r1:4d} "
            f"{r2:4d} "
            f"{r3:4d} "
            f"{vals[0]:>9} "
            f"{vals[1]:>9} "
            f"{vals[2]:>9} "
            f"{vals[3]:>9}"
        )

    # --------------------------------------------------------------------------------------------
    # MATHEMATICAL INTERPRETATION
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
This experiment removes the strongest weakness of the previous SMT
test: no true residue is given to the solver.

For every modulus r_i the solver must discover:

    a_i = p mod r_i
    b_i = q mod r_i

itself.

The four models are:

MODEL 0

    p*q = n

This is the nonlinear integer factorization baseline.

MODEL 1

    p = a_i + k_i*r_i
    q = b_i + l_i*r_i

with unknown a_i,b_i.

This tests the coordinate representation.

MODEL 2 additionally states:

    a_i*b_i == n (mod r_i)

but the residues remain unknown.

MODEL 3 exposes the close-modulus relationships:

    r2 = r1 + d2
    r3 = r1 + d3

and forces the same p and q to have consistent coordinates at all
three moduli.

The central question is now:

    does the small modulus gap create a solver-friendly constraint
    that is substantially easier than p*q=n alone?

This is still NOT evidence of a new factoring algorithm.

The decisive measurement is scaling with n.

For a useful result we would want:

    solver time grows slowly with factor size
    while ordinary factor enumeration grows rapidly.

A merely successful SAT/SMT solve on 32-bit or 40-bit examples is
not enough, because the solver may simply be performing implicit
factor search.

The next serious test is therefore to increase:

    10^5
    10^6
    10^7
    10^8
    ...

and record:

    bit length
    explicit prime candidates
    SMT0 time
    SMT1 time
    SMT2 time
    SMT3 time
    timeout rate
    exact recovery.

If SMT3 remains stable while the explicit candidate domain becomes
much larger, then the close-modulus coordinate representation has
found something worth investigating.

If all models degrade rapidly with bit length, the conclusion is
that SMT is changing the representation of the search rather than
removing the underlying complexity.

No true residue values are supplied to any model.
"""
    )

    # --------------------------------------------------------------------------------------------
    # TIMING
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime              = "
        f"{time.perf_counter() - total_start:.3f} seconds"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()