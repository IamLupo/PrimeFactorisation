import math
import random
import statistics
import time
from collections import Counter

# ============================================================
# THREE-CLOSE-PRIME DIRECT K-STATE / RESIDUE-LATTICE EXPERIMENT
# ============================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ANCHORS = 300
SEARCH_ANCHORS = 300

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20
SEED = 1_511_464_998

random.seed(SEED)


# ============================================================
# PRIME GENERATION
# ============================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    r = math.isqrt(n)
    d = 3
    while d <= r:
        if n % d == 0:
            return False
        d += 2
    return True


def primes_in_range(lo: int, hi: int):
    return [x for x in range(lo, hi + 1) if is_prime(x)]


# ============================================================
# ANCHORS
# ============================================================

def build_anchors(factor_primes, count):
    """
    Construct deterministic factor pairs.

    We prefer pairs whose product is close to M modulo the same
    style of anchor selection used by the preceding experiments.
    The exact anchor-selection mechanism is deliberately explicit
    here so the experiment is reproducible.
    """
    rng = random.Random(SEED)

    pairs = []

    # Generate a large random candidate pool.
    for _ in range(count * 40):
        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        if p >= q:
            continue

        n = p * q

        # Keep products in the intended scale.
        if FACTOR_MIN <= p <= FACTOR_MAX and FACTOR_MIN <= q <= FACTOR_MAX:
            pairs.append((n, p, q))

    # Remove duplicate pairs.
    seen = set()
    unique = []

    for item in pairs:
        key = (item[1], item[2])
        if key in seen:
            continue
        seen.add(key)
        unique.append(item)

    # Sort toward products close to M * 1-ish scale.
    unique.sort(
        key=lambda x: (
            abs((x[0] / max(M, 1)) - round(x[0] / max(M, 1))),
            -x[0],
        )
    )

    if len(unique) < count:
        raise RuntimeError(
            f"Unable to construct {count} usable anchors; got {len(unique)}."
        )

    return unique[:count]


# ============================================================
# CLOSE MODULUS TRIPLES
# ============================================================

def select_close_triples(mod_primes, anchors, count):
    """
    Select an anchor-specific triple.

    We score triples by how close the three primes are, while
    avoiding repeated triples where possible.
    """
    rng = random.Random(SEED)

    triples = []

    for n, p, q in anchors:
        # Prefer moduli near the geometric scale of sqrt(n).
        target = math.isqrt(n)

        candidates = sorted(
            mod_primes,
            key=lambda r: abs(r - target)
        )

        # Restrict to a local neighbourhood.
        local = candidates[: min(40, len(candidates))]

        best = None
        best_score = None

        for i in range(len(local)):
            for j in range(i + 1, len(local)):
                for k in range(j + 1, len(local)):
                    r1, r2, r3 = local[i], local[j], local[k]

                    span = r3 - r1
                    if span > CLOSE_RATIO * max(r1, 1):
                        continue

                    if math.gcd(r1, r2) != 1:
                        continue
                    if math.gcd(r1, r3) != 1:
                        continue
                    if math.gcd(r2, r3) != 1:
                        continue

                    score = (
                        span,
                        abs(r2 - r1),
                        abs(r3 - r2),
                        -r1,
                    )

                    if best_score is None or score < best_score:
                        best_score = score
                        best = (r1, r2, r3)

        if best is None:
            raise RuntimeError(
                f"No suitable modulus triple for n={n}."
            )

        triples.append(best)

    # Make deterministic but allow naturally varying triples.
    return triples


# ============================================================
# MODULAR HELPERS
# ============================================================

def inv_mod(a: int, p: int) -> int:
    g = math.gcd(a, p)
    if g != 1:
        raise ValueError(f"{a} is not invertible modulo {p}.")
    return pow(a, -1, p)


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


# ============================================================
# VALIDATION
# ============================================================

def validate_triple(n: int, p: int, q: int, mods):
    r1, r2, r3 = mods

    if p * q != n:
        return False

    # Actual residues.
    a1 = p % r1
    a3 = p % r3

    # Recover k from the r1 representation.
    k = (p - a1) // r1

    # Direct r3 coordinate elimination:
    #
    # a1 + k*r1 == a3 (mod r3)
    #
    # => k == (a3-a1)*r1^{-1} (mod r3)
    k0 = ((a3 - a1) * inv_mod(r1, r3)) % r3

    if k % r3 != k0:
        return False

    # Reconstruct p.
    if a1 + k * r1 != p:
        return False

    return True


# ============================================================
# DIRECT K-STATE SEARCH
# ============================================================

def direct_k_search(n: int, true_p: int, true_q: int, mods):
    r1, r2, r3 = mods

    # We only search p <= q to avoid duplicate factor pairs.
    p_lo = FACTOR_MIN
    p_hi = min(FACTOR_MAX, math.isqrt(n))

    # --------------------------------------------------------
    # 1. Enumerate possible a = p mod r1.
    #
    # This is a residue domain, not a p domain.
    # --------------------------------------------------------

    a_states = 0
    residue_lattice_states = 0
    feasible_k_states = 0
    hyperbola_states = 0
    l_residue_states = 0
    final_candidates = 0
    exact_solutions = []

    # We do not enumerate p.
    #
    # For each a, determine the k range that puts p into [p_lo,p_hi].
    for a in range(r1):
        k_min = max(
            0,
            ceil_div(p_lo - a, r1),
        )

        k_max = (p_hi - a) // r1

        if k_max < k_min:
            continue

        a_states += 1

        # ----------------------------------------------------
        # Third-modulus residue states.
        #
        # We enumerate a3, but directly convert each (a,a3)
        # state into one k congruence.
        # ----------------------------------------------------

        for a3 in range(r3):
            residue_lattice_states += 1

            k0 = ((a3 - a) * inv_mod(r1, r3)) % r3

            # Find the first/last k in [k_min,k_max]
            # satisfying k == k0 mod r3.
            m_min = ceil_div(k_min - k0, r3)
            m_max = (k_max - k0) // r3

            if m_max < m_min:
                continue

            # In most realistic instances there will be at most
            # one such k because r3 is large relative to k-range.
            for m in range(m_min, m_max + 1):
                k = k0 + m * r3

                if not (k_min <= k <= k_max):
                    continue

                p = a + k * r1

                if not (p_lo <= p <= p_hi):
                    continue

                feasible_k_states += 1

                if p <= 1:
                    continue

                # ------------------------------------------------
                # Hyperbola determines the corresponding q interval.
                # ------------------------------------------------

                q_real = n / p

                # q must be >= p because p <= sqrt(n).
                q_center = int(q_real)

                q_candidates = {
                    q_center - 1,
                    q_center,
                    q_center + 1,
                }

                hyperbola_states += len(q_candidates)

                for q in q_candidates:
                    if q < FACTOR_MIN or q > FACTOR_MAX:
                        continue
                    if p * q != n:
                        continue

                    if not is_prime(q):
                        continue

                    # Third-modulus consistency for q.
                    b3 = q % r3

                    # Equivalent direct check.
                    if (q - b3) % r3 != 0:
                        continue

                    l_residue_states += 1

                    exact_solutions.append((p, q))

    exact_solutions = sorted(set(exact_solutions))

    if true_p * true_q == n:
        true_pair = tuple(sorted((true_p, true_q)))
    else:
        true_pair = None

    recovered = true_pair in exact_solutions if true_pair else False

    return {
        "a_states": a_states,
        "residue_lattice_states": residue_lattice_states,
        "feasible_k_states": feasible_k_states,
        "hyperbola_states": hyperbola_states,
        "l_residue_states": l_residue_states,
        "final_candidates": len(exact_solutions),
        "exact_solutions": len(exact_solutions),
        "recovered": recovered,
        "solutions": exact_solutions,
    }


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def run():
    t0 = time.perf_counter()

    print("=" * 100)
    print("THREE-CLOSE-PRIME DIRECT K-STATE / RESIDUE-LATTICE EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range       = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # --------------------------------------------------------
    # Prime pools
    # --------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = primes_in_range(FACTOR_MIN, FACTOR_MAX)
    modulus_primes = primes_in_range(MOD_MIN, MOD_MAX)

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # --------------------------------------------------------
    # Anchors
    # --------------------------------------------------------

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(factor_primes, ANCHORS)

    print(f"actual anchors            = {len(anchors):,}")
    print()

    # --------------------------------------------------------
    # Modulus triples
    # --------------------------------------------------------

    print("=" * 100)
    print("SELECTING ANCHOR-SPECIFIC CLOSE TRIPLES")
    print("=" * 100)

    triples = select_close_triples(
        modulus_primes,
        anchors,
        ANCHORS,
    )

    unique_triples = len(set(triples))

    for i in range(25, ANCHORS + 1, 25):
        print(f"anchor {i:3d}/{ANCHORS}")

    print()
    print("=" * 100)
    print("VALIDATING THREE-MODULUS RESIDUE IDENTITY")
    print("=" * 100)

    failures = 0

    for (n, p, q), mods in zip(anchors, triples):
        if not validate_triple(n, p, q, mods):
            failures += 1

    print(f"identity failures         = {failures}")
    print(f"unique modulus triples   = {unique_triples}")
    print()

    if failures:
        raise RuntimeError(
            "Direct k residue identity validation failed."
        )

    # --------------------------------------------------------
    # Search
    # --------------------------------------------------------

    print("=" * 100)
    print("RUNNING DIRECT K-STATE SEARCH")
    print("=" * 100)

    stats = []

    for idx, ((n, p, q), mods) in enumerate(
        zip(anchors[:SEARCH_ANCHORS], triples[:SEARCH_ANCHORS]),
        start=1,
    ):
        result = direct_k_search(
            n,
            p,
            q,
            mods,
        )

        result["n"] = n
        result["p"] = p
        result["q"] = q
        result["r1"] = mods[0]
        result["r2"] = mods[1]
        result["r3"] = mods[2]

        stats.append(result)

        if idx % 25 == 0:
            print(f"anchor {idx:3d}/{SEARCH_ANCHORS}")

    # --------------------------------------------------------
    # Aggregate
    # --------------------------------------------------------

    def avg(key):
        return statistics.mean(x[key] for x in stats)

    recovered = sum(x["recovered"] for x in stats)

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(f"anchors analyzed            = {len(stats)}")
    print(f"average a states            = {avg('a_states'):.3f}")
    print(f"average residue-lattice states = {avg('residue_lattice_states'):.3f}")
    print(f"average feasible k states   = {avg('feasible_k_states'):.3f}")
    print(f"average hyperbola states    = {avg('hyperbola_states'):.3f}")
    print(f"average l residue states    = {avg('l_residue_states'):.3f}")
    print(f"average final candidates    = {avg('final_candidates'):.3f}")
    print(f"average exact solutions     = {avg('exact_solutions'):.3f}")
    print()

    # --------------------------------------------------------
    # Key ratios
    # --------------------------------------------------------

    print("=" * 100)
    print("DIRECT SEARCH REDUCTION")
    print("=" * 100)

    avg_k = avg("feasible_k_states")
    avg_residue = avg("residue_lattice_states")
    avg_final = avg("final_candidates")

    k_domain_estimate = (
        statistics.mean(
            max(
                0,
                (min(FACTOR_MAX, math.isqrt(x["n"])) // x["r1"])
                - (FACTOR_MIN // x["r1"])
                + 2,
            )
            for x in stats
        )
    )

    print(f"estimated raw k-domain      = {k_domain_estimate:.3f}")

    if k_domain_estimate:
        print(
            f"feasible k / raw k-domain  = "
            f"{avg_k / k_domain_estimate:.9f}"
        )

    print(
        f"final / feasible k          = "
        f"{avg_final / max(avg_k, 1):.9f}"
    )

    print()

    # --------------------------------------------------------
    # Recovery
    # --------------------------------------------------------

    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)
    print(f"correctly recovered         = {recovered}/{len(stats)}")
    print(
        f"recovery rate               = "
        f"{100.0 * recovered / len(stats):.4f}%"
    )
    print()

    # --------------------------------------------------------
    # Distribution
    # --------------------------------------------------------

    print("=" * 100)
    print("FEASIBLE-k DISTRIBUTION")
    print("=" * 100)

    dist = Counter(x["feasible_k_states"] for x in stats)

    for value, count in sorted(dist.items()):
        print(
            f"k states = {value:5d} anchors = {count:4d}"
        )

    print()

    # --------------------------------------------------------
    # Strongest / weakest
    # --------------------------------------------------------

    strongest = sorted(
        stats,
        key=lambda x: x["feasible_k_states"]
    )[:20]

    weakest = sorted(
        stats,
        key=lambda x: x["feasible_k_states"],
        reverse=True
    )[:20]

    print("=" * 100)
    print("STRONGEST DIRECT-k COLLAPSES")
    print("=" * 100)

    for x in strongest:
        print(
            f"n={x['n']:,} "
            f"p={x['p']:,} q={x['q']:,} "
            f"mods=({x['r1']},{x['r2']},{x['r3']}) "
            f"a={x['a_states']} "
            f"residue={x['residue_lattice_states']} "
            f"k={x['feasible_k_states']} "
            f"hyper={x['hyperbola_states']} "
            f"final={x['final_candidates']} "
            f"recovered={x['recovered']}"
        )

    print()

    print("=" * 100)
    print("WEAKEST DIRECT-k COLLAPSES")
    print("=" * 100)

    for x in weakest:
        print(
            f"n={x['n']:,} "
            f"p={x['p']:,} q={x['q']:,} "
            f"mods=({x['r1']},{x['r2']},{x['r3']}) "
            f"a={x['a_states']} "
            f"residue={x['residue_lattice_states']} "
            f"k={x['feasible_k_states']} "
            f"hyper={x['hyperbola_states']} "
            f"final={x['final_candidates']} "
            f"recovered={x['recovered']}"
        )

    print()

    # --------------------------------------------------------
    # Anchor table
    # --------------------------------------------------------

    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)
    print(
        " ID       p       q    r1    r2    r3    "
        "A   RESIDUE       K   HYPER   FINAL EXACT"
    )
    print("-" * 100)

    for i, x in enumerate(stats, start=1):
        print(
            f"{i:3d} "
            f"{x['p']:7d} "
            f"{x['q']:7d} "
            f"{x['r1']:5d} "
            f"{x['r2']:5d} "
            f"{x['r3']:5d} "
            f"{x['a_states']:4d} "
            f"{x['residue_lattice_states']:10d} "
            f"{x['feasible_k_states']:5d} "
            f"{x['hyperbola_states']:7d} "
            f"{x['final_candidates']:5d} "
            f"{int(x['recovered']):5d}"
        )

    # --------------------------------------------------------
    # Timing
    # --------------------------------------------------------

    elapsed = time.perf_counter() - t0

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
This experiment tests the proposed transition from candidate
filtering to direct quotient-state generation.

We write:

    p = a + k*r1

with:

    0 <= a < r1.

The third modulus gives:

    p == a3 (mod r3)

and therefore:

    a + k*r1 == a3 (mod r3)

so:

    k == (a3-a)*r1^(-1) (mod r3).

For each pair (a,a3), there is therefore one arithmetic
progression of possible k values:

    k = k0 + m*r3.

The allowed p interval gives a finite interval for k.

The experiment intersects those two conditions directly.

No p values are enumerated.

After a surviving k is obtained:

    p = a + k*r1

and the hyperbola gives:

    q = n/p.

The final candidate is then checked exactly.

The important distinction is:

    OLD:
        generate candidates
        -> modular filtering

    NEW:
        residue state
        -> directly generate feasible k
        -> reconstruct p
        -> derive q
        -> exact verification.

The central measurement is therefore:

    feasible-k states

compared with the actual raw k domain.

If feasible-k states are consistently around O(1), the residue
condition is genuinely collapsing the quotient coordinate before
candidate generation.

If feasible-k states remain comparable to the complete k domain,
then the third residue is merely reordering/reparameterizing the
same search.

There is another important structural observation:

    k = floor(p/r1)

and p is restricted to 10,000..100,000.

Consequently k is intrinsically small:

    k = O(FACTOR_MAX/r1).

Any experiment reporting tens of thousands of k states for this
literal parameterization is therefore counting a different
coordinate or has a quotient-range bug.

This experiment makes that distinction explicit.

Every recovered factor pair is finally checked with:

    p*q == n.

"""
    )

    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"total runtime              = {elapsed:.3f} seconds")
    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
