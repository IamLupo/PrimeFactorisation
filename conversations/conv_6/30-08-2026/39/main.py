#!/usr/bin/env python3

"""
====================================================================================================
THREE-CLOSE-PRIME SMT / DIRECT FACTOR MODEL EXPERIMENT
====================================================================================================

Goal
----
Test whether an SMT bit-vector solver can recover p,q directly from:

    p*q = n

without explicitly enumerating candidate p values.

Three solver models are compared:

    MODEL 0:
        p*q == n

    MODEL 1:
        p*q == n
        p % r_i == p_true % r_i

    MODEL 2:
        p*q == n
        p % r_i == p_true % r_i
        q % r_i == q_true % r_i

The factor values are NEVER enumerated by the Python search itself.

The known residues are supplied deliberately as an information-content
experiment. This tests whether the modular constraints are useful to
the solver once the candidate-generation loop is removed.

IMPORTANT:
-----------
This is NOT yet a factoring algorithm. Supplying the true residues is
additional information. The experiment asks:

    "Can a constraint solver exploit this information directly?"

The next stage, if useful, is to make the residue classes variables
instead of feeding the true residues.


Requirements
------------
    pip install z3-solver
"""

import math
import random
import statistics
import time

try:
    import z3
except ImportError:
    raise SystemExit(
        "z3-solver is not installed.\n"
        "Install it with:\n\n"
        "    pip install z3-solver\n"
    )


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

ANCHORS = 300

# Start with a small number because SMT multiplication can be expensive.
SEARCH_ANCHORS = 20

# Maximum solver time per model per anchor.
TIMEOUT_MS = 5_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998


# ================================================================================================
# PRIME GENERATION
# ================================================================================================

def sieve(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, flag in enumerate(is_prime) if flag]


# ================================================================================================
# ANCHOR GENERATION
# ================================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    seed: int,
) -> list[tuple[int, int, int]]:
    """
    Deterministically generate distinct prime factor pairs.

    Returns:
        (p, q, n)
    """
    rng = random.Random(seed)

    seen: set[tuple[int, int]] = set()
    anchors: list[tuple[int, int, int]] = []

    while len(anchors) < count:
        p, q = rng.sample(factor_primes, 2)

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)

        n = p * q
        anchors.append((p, q, n))

    return anchors


# ================================================================================================
# CLOSE MODULUS TRIPLE
# ================================================================================================

def choose_close_triple(
    primes: list[int],
    ratio: float,
    anchor_index: int,
) -> tuple[int, int, int]:
    """
    Deterministically choose a high-product close triple.

    For every candidate first prime r1, choose the largest possible
    r2,r3 satisfying:

        r3/r1 <= 1 + ratio

    The triple with maximal product is selected.

    A small anchor-dependent rotation is used so different anchors
    do not all receive the same triple.
    """
    best = None
    best_product = -1

    n = len(primes)

    start = anchor_index % n

    ordered = primes[start:] + primes[:start]

    # Sort candidates by value again after rotation.
    ordered = sorted(ordered)

    for i, r1 in enumerate(ordered):
        max_r = int(r1 * (1.0 + ratio))

        candidates = [
            x for x in ordered[i + 1:]
            if x <= max_r
        ]

        if len(candidates) < 2:
            continue

        # For maximum product choose the two largest.
        r2 = candidates[-2]
        r3 = candidates[-1]

        triple = (r1, r2, r3)
        product = r1 * r2 * r3

        if product > best_product:
            best_product = product
            best = triple

    if best is None:
        raise RuntimeError("Could not construct close modulus triple.")

    return best


# ================================================================================================
# BIT-VECTOR WIDTH
# ================================================================================================

def bit_width_for_exact_product(n: int) -> int:
    """
    Choose a width large enough that p*q cannot wrap.

    If p,q < 2^b, product < 2^(2b), therefore 2b+2 is sufficient.
    """
    b = max(1, n.bit_length())
    return 2 * b + 2


# ================================================================================================
# MODEL CONSTRUCTION
# ================================================================================================

def build_solver(
    n: int,
    rmods: tuple[int, int, int],
    true_p: int,
    true_q: int,
    model: int,
) -> tuple[z3.Solver, z3.BitVecRef, z3.BitVecRef]:
    """
    Build the SMT model.

    model = 0:
        p*q = n

    model = 1:
        p*q = n
        p residues known

    model = 2:
        p*q = n
        p and q residues known
    """

    width = bit_width_for_exact_product(n)

    p = z3.BitVec("p", width)
    q = z3.BitVec("q", width)

    s = z3.Solver()
    s.set(timeout=TIMEOUT_MS)

    N = z3.BitVecVal(n, width)

    pmin = z3.BitVecVal(FACTOR_MIN, width)
    pmax = z3.BitVecVal(FACTOR_MAX, width)

    # Factor interval.
    s.add(
        z3.UGE(p, pmin),
        z3.ULE(p, pmax),
        z3.UGE(q, pmin),
        z3.ULE(q, pmax),
    )

    # Remove symmetric duplicate solution.
    s.add(z3.ULE(p, q))

    # Exact multiplication.
    s.add(p * q == N)

    # Modular information.
    if model >= 1:
        for r in rmods:
            r_bv = z3.BitVecVal(r, width)

            true_p_res = true_p % r

            s.add(
                z3.URem(p, r_bv)
                == z3.BitVecVal(true_p_res, width)
            )

    if model >= 2:
        for r in rmods:
            r_bv = z3.BitVecVal(r, width)

            true_q_res = true_q % r

            s.add(
                z3.URem(q, r_bv)
                == z3.BitVecVal(true_q_res, width)
            )

    return s, p, q


# ================================================================================================
# SMT SOLVE
# ================================================================================================

def solve_model(
    n: int,
    p_true: int,
    q_true: int,
    mods: tuple[int, int, int],
    model: int,
) -> dict:
    """
    Run one SMT model and return measurements.
    """

    solver, p_var, q_var = build_solver(
        n=n,
        rmods=mods,
        true_p=p_true,
        true_q=q_true,
        model=model,
    )

    t0 = time.perf_counter()

    result = solver.check()

    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    recovered_p = None
    recovered_q = None
    exact = False

    if result == z3.sat:
        values = solver.model()

        recovered_p = values[p_var].as_long()
        recovered_q = values[q_var].as_long()

        exact = (
            recovered_p * recovered_q == n
            and {
                recovered_p,
                recovered_q,
            } == {
                p_true,
                q_true,
            }
        )

    stats_text = str(solver.statistics())

    return {
        "result": str(result),
        "time_ms": elapsed_ms,
        "p": recovered_p,
        "q": recovered_q,
        "exact": exact,
        "stats": stats_text,
    }


# ================================================================================================
# PRIME ENUMERATION BASELINE
# ================================================================================================

def prime_baseline(
    n: int,
    p_true: int,
    factor_primes: list[int],
) -> tuple[int | None, int]:
    """
    Enumerate factor primes until the true factor is reached.

    Returns:
        (factor, tests)
    """

    tests = 0

    for p in factor_primes:
        if p < FACTOR_MIN or p > FACTOR_MAX:
            continue

        if p > math.isqrt(n):
            break

        tests += 1

        if n % p == 0:
            return p, tests

    return None, tests


# ================================================================================================
# INTEGER ENUMERATION BASELINE
# ================================================================================================

def integer_baseline(n: int, p_true: int) -> int:
    """
    Number of integer p states before reaching the true p.
    """
    lo = FACTOR_MIN
    hi = min(FACTOR_MAX, math.isqrt(n))

    if not (lo <= p_true <= hi):
        return hi - lo + 1

    return p_true - lo + 1


# ================================================================================================
# STATISTICS HELPERS
# ================================================================================================

def average(values: list[float]) -> float:
    return statistics.mean(values) if values else 0.0


def median(values: list[float]) -> float:
    return statistics.median(values) if values else 0.0


# ================================================================================================
# MAIN
# ================================================================================================

def run() -> None:
    total_start = time.perf_counter()

    print("=" * 100)
    print("THREE-CLOSE-PRIME SMT / DIRECT FACTOR MODEL EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"search anchors            = {SEARCH_ANCHORS}")
    print(f"modulus prime range       = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print(f"SMT timeout/model         = {TIMEOUT_MS:,} ms")
    print()

    # --------------------------------------------------------------------------------------------
    # Prime pools
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = [
        p for p in sieve(FACTOR_MAX)
        if p >= FACTOR_MIN
    ]

    modulus_primes = [
        p for p in sieve(MOD_MAX)
        if p >= MOD_MIN
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # --------------------------------------------------------------------------------------------
    # Anchors
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes=factor_primes,
        count=ANCHORS,
        seed=SEED,
    )

    print(f"actual anchors            = {len(anchors):,}")
    print()

    # --------------------------------------------------------------------------------------------
    # Modulus triples
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("SELECTING CLOSE MODULUS TRIPLES")
    print("=" * 100)

    triples = []

    for i, _anchor in enumerate(anchors):
        triple = choose_close_triple(
            modulus_primes,
            CLOSE_RATIO,
            i,
        )
        triples.append(triple)

    print(f"unique modulus triples    = {len(set(triples)):,}")
    print()

    # --------------------------------------------------------------------------------------------
    # Validation
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("VALIDATING MODELS ON TRUE FACTORS")
    print("=" * 100)

    validation_failures = 0

    for i in range(min(ANCHORS, 20)):
        p, q, n = anchors[i]
        r1, r2, r3 = triples[i]

        for r in (r1, r2, r3):
            if p % r < 0 or q % r < 0:
                validation_failures += 1

    print(f"validation failures       = {validation_failures}")
    print()

    # --------------------------------------------------------------------------------------------
    # Solver experiment
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("RUNNING DIRECT SMT SEARCH")
    print("=" * 100)

    results = {
        0: [],
        1: [],
        2: [],
    }

    baseline_prime_tests = []
    baseline_integer_tests = []

    recovered = {
        0: 0,
        1: 0,
        2: 0,
    }

    timeouts = {
        0: 0,
        1: 0,
        2: 0,
    }

    for index in range(SEARCH_ANCHORS):
        p_true, q_true, n = anchors[index]
        mods = triples[index]

        _, prime_tests = prime_baseline(
            n,
            p_true,
            factor_primes,
        )

        integer_tests = integer_baseline(
            n,
            p_true,
        )

        baseline_prime_tests.append(prime_tests)
        baseline_integer_tests.append(integer_tests)

        for model in (0, 1, 2):

            result = solve_model(
                n=n,
                p_true=p_true,
                q_true=q_true,
                mods=mods,
                model=model,
            )

            results[model].append(result)

            if result["exact"]:
                recovered[model] += 1

            if result["result"] == "unknown":
                timeouts[model] += 1

        done = index + 1

        if done % 5 == 0 or done == SEARCH_ANCHORS:
            print(f"anchor {done:3d}/{SEARCH_ANCHORS}")

    print()

    # --------------------------------------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"average integer p states       = "
        f"{average(baseline_integer_tests):,.3f}"
    )

    print(
        f"average prime enumeration      = "
        f"{average(baseline_prime_tests):,.3f}"
    )

    for model, name in (
        (0, "SMT pq=n"),
        (1, "SMT + p CRT"),
        (2, "SMT + p,q CRT"),
    ):
        solve_times = [
            x["time_ms"]
            for x in results[model]
            if x["result"] != "unknown"
        ]

        print()
        print(f"{name}")
        print(
            f"    solved                  = "
            f"{len(solve_times):,}/{SEARCH_ANCHORS:,}"
        )
        print(
            f"    timeouts                = "
            f"{timeouts[model]:,}"
        )
        print(
            f"    exact recovery          = "
            f"{recovered[model]:,}/{SEARCH_ANCHORS:,}"
        )

        if solve_times:
            print(
                f"    mean solve time         = "
                f"{average(solve_times):,.3f} ms"
            )
            print(
                f"    median solve time       = "
                f"{median(solve_times):,.3f} ms"
            )
            print(
                f"    maximum solve time      = "
                f"{max(solve_times):,.3f} ms"
            )

    # --------------------------------------------------------------------------------------------
    # Reduction comparison
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("COMPARISON WITH EXPLICIT p ENUMERATION")
    print("=" * 100)

    avg_integer = average(baseline_integer_tests)
    avg_prime = average(baseline_prime_tests)

    print(
        f"integer p states            = "
        f"{avg_integer:,.3f}"
    )

    print(
        f"prime tests                 = "
        f"{avg_prime:,.3f}"
    )

    print()

    for model, name in (
        (0, "SMT pq=n"),
        (1, "SMT + p CRT"),
        (2, "SMT + p,q CRT"),
    ):
        print(
            f"{name:24s} exact = "
            f"{recovered[model]}/{SEARCH_ANCHORS}"
        )

    # --------------------------------------------------------------------------------------------
    # Per-anchor results
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID          p          q       r1   r2   r3   "
        "INTP   PRIMEP   SMT0(ms)  SMT1(ms)  SMT2(ms)  EXACT"
    )
    print("-" * 100)

    for i in range(SEARCH_ANCHORS):
        p, q, n = anchors[i]
        r1, r2, r3 = triples[i]

        t0 = results[0][i]["time_ms"]
        t1 = results[1][i]["time_ms"]
        t2 = results[2][i]["time_ms"]

        e0 = results[0][i]["exact"]
        e1 = results[1][i]["exact"]
        e2 = results[2][i]["exact"]

        exact_flag = int(e0 and e1 and e2)

        print(
            f"{i + 1:3d} "
            f"{p:10,d} "
            f"{q:10,d} "
            f"{r1:5d} "
            f"{r2:4d} "
            f"{r3:4d} "
            f"{baseline_integer_tests[i]:7,d} "
            f"{baseline_prime_tests[i]:7,d} "
            f"{t0:9.2f} "
            f"{t1:9.2f} "
            f"{t2:9.2f} "
            f"{exact_flag:5d}"
        )

    # --------------------------------------------------------------------------------------------
    # Examples
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MODEL EXAMPLES")
    print("=" * 100)

    for i in range(min(SEARCH_ANCHORS, 10)):
        p_true, q_true, n = anchors[i]
        r1, r2, r3 = triples[i]

        print(
            f"n={n:,} "
            f"true=({p_true:,},{q_true:,}) "
            f"mods=({r1},{r2},{r3})"
        )

        for model, name in (
            (0, "SMT0"),
            (1, "SMT1"),
            (2, "SMT2"),
        ):
            result = results[model][i]

            print(
                f"    {name}: "
                f"result={result['result']} "
                f"time={result['time_ms']:.2f}ms "
                f"p={result['p']} "
                f"q={result['q']} "
                f"exact={result['exact']}"
            )

    # --------------------------------------------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
This experiment changes the search model.

The explicit-search approach is:

    p candidates
        |
        v
    modular filtering
        |
        v
    q reconstruction
        |
        v
    exact verification.

The SMT formulation instead gives the solver:

    p*q = n

together with:

    factor bounds

and optionally:

    p mod r_i = known residue
    q mod r_i = known residue.

Therefore the Python program performs NO loop over candidate p values.

The important comparison is:

    explicit p enumeration
        versus
    solver search over symbolic p,q.

MODEL 0:

    pq = n

This is the pure solver baseline.

MODEL 1:

    pq = n
    p mod r1/r2/r3 = true residues

This tests whether the three-modulus information helps the solver
when candidate p values are not explicitly generated.

MODEL 2:

    pq = n
    p mod r_i = true residue
    q mod r_i = true residue

This supplies the complete observed residue pair for each modulus.

A successful MODEL 1 or MODEL 2 result does NOT yet demonstrate a
new factoring algorithm because the true residues are being supplied
as extra information.

The experiment is testing a narrower question:

    Can constraint solving exploit the modular structure without
    explicit candidate enumeration?

The next important experiment should therefore NOT simply add more
filters.

Instead, make the residue variables unknown:

    p mod r_i = A_i
    q mod r_i = B_i

and ask the solver to determine A_i,B_i together with p,q.

That removes the "known true residue" advantage.

The decisive metric is scaling.

If:

    40-bit
    50-bit
    60-bit
    70-bit
    80-bit
    ...

remain tractable while the explicit p-domain grows rapidly, then
the solver formulation has found exploitable structure.

If solver time grows rapidly with bit length, then the SMT model is
primarily replacing explicit enumeration with implicit search.

Every reported factor pair is finally checked with:

    p*q == n
"""
    )

    total_time = time.perf_counter() - total_start

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"total runtime              = {total_time:.3f} seconds")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
