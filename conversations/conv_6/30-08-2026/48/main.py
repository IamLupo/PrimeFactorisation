import math
import random
import statistics
import time


# =============================================================================
# CONFIGURATION
# =============================================================================

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS_EVERY = 10


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start::p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [
        i
        for i, flag in enumerate(sieve)
        if flag
    ]


# =============================================================================
# ANCHORS
# =============================================================================

def build_anchors(factor_primes, count, seed):
    rng = random.Random(seed)

    anchors = []

    while len(anchors) < count:
        p, q = rng.sample(factor_primes, 2)

        if p > q:
            p, q = q, p

        anchors.append((p, q))

    return anchors


# =============================================================================
# MODULUS PAIRS
# =============================================================================

def build_close_pairs(modulus_primes):
    pairs = []

    for i, r1 in enumerate(modulus_primes):
        for r2 in modulus_primes[i + 1:]:
            if r2 > r1 * (1.0 + CLOSE_RATIO):
                break

            pairs.append((r1, r2))

    return pairs


def choose_modulus_pair(n, close_pairs, rng):
    valid = [
        pair
        for pair in close_pairs
        if pair[0] * pair[1] < n
    ]

    if not valid:
        return None

    valid.sort(
        key=lambda pair: abs(n - pair[0] * pair[1])
    )

    # Randomize slightly among the closest pairs so that we do not
    # always select the exact same pair.
    top = valid[:min(40, len(valid))]

    return rng.choice(top)


# =============================================================================
# FACTOR K
# =============================================================================

def factor_pairs(k_value):
    """
    Enumerate positive divisor pairs (k,l) such that k*l = k_value.
    """
    if k_value <= 0:
        return []

    result = []

    root = math.isqrt(k_value)

    for d in range(1, root + 1):
        if k_value % d == 0:
            result.append(
                (d, k_value // d)
            )

    return result


# =============================================================================
# INTEGER CEILING DIVISION
# =============================================================================

def ceil_div(a, b):
    return -((-a) // b)


# =============================================================================
# TRUE COORDINATES
# =============================================================================

def true_coordinates(p, q, r1, r2):
    a = p % r1
    b = q % r2

    k = p // r1
    l = q // r2

    # Cross-term quotients.
    c1 = (a * l) // r1
    c2 = (b * k) // r2

    # Cross-term remainders.
    d1 = (a * l) % r1
    d2 = (b * k) % r2

    # Final carry.
    c3 = (
        d1 * r2
        + d2 * r1
        + a * b
    ) // (r1 * r2)

    E = c1 + c2 + c3

    return {
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "E": E,
    }


# =============================================================================
# CARRY -> RESIDUE INTERVALS
# =============================================================================

def carry_intervals(k, l, c1, c2, r1, r2):
    """
    From

        c1 = floor(a*l/r1)
        c2 = floor(b*k/r2)

    derive all possible integer values of a and b.
    """

    # c1 <= a*l/r1 < c1+1
    #
    # => c1*r1 <= a*l < (c1+1)*r1

    a_lo = ceil_div(
        c1 * r1,
        l,
    )

    a_hi = (
        ((c1 + 1) * r1) - 1
    ) // l

    # c2 <= b*k/r2 < c2+1

    b_lo = ceil_div(
        c2 * r2,
        k,
    )

    b_hi = (
        ((c2 + 1) * r2) - 1
    ) // k

    # Residue ranges.
    a_lo = max(0, a_lo)
    a_hi = min(r1 - 1, a_hi)

    b_lo = max(0, b_lo)
    b_hi = min(r2 - 1, b_hi)

    return (
        a_lo,
        a_hi,
        b_lo,
        b_hi,
    )


# =============================================================================
# ONE ANCHOR
# =============================================================================

def run_anchor(p_true, q_true, r1, r2):

    n = p_true * q_true

    R = r1 * r2

    T = n // R

    true_data = true_coordinates(
        p_true,
        q_true,
        r1,
        r2,
    )

    true_k = true_data["k"]
    true_l = true_data["l"]

    true_E = true_data["E"]

    true_K = true_k * true_l

    # -------------------------------------------------------------------------
    # E WINDOW
    # -------------------------------------------------------------------------

    # Since:
    #
    #     E = T - k*l
    #
    # and:
    #
    #     k <= FACTOR_MAX/r1
    #     l <= FACTOR_MAX/r2
    #
    # we use the safe upper bound:
    #
    #     E <= k+l

    k_max = FACTOR_MAX // r1
    l_max = FACTOR_MAX // r2

    E_max = min(
        T,
        k_max + l_max,
    )

    E_candidate_count = 0

    # -------------------------------------------------------------------------
    # COUNTERS
    # -------------------------------------------------------------------------

    root_pairs = 0
    admissible_pairs = 0

    carry_compatible_pairs = 0
    carry_cells = 0

    carry_a_states = 0
    carry_b_states = 0

    product_range_survivors = 0

    exact_pairs = set()

    true_E_seen = False
    true_carry_seen = False

    # -------------------------------------------------------------------------
    # E SEARCH
    # -------------------------------------------------------------------------

    for E in range(E_max + 1):

        K = T - E

        if K <= 0:
            continue

        E_candidate_count += 1

        # ---------------------------------------------------------------------
        # FACTOR K
        # ---------------------------------------------------------------------

        pairs = factor_pairs(K)

        for k, l in pairs:

            # Candidate quotient rectangle:
            #
            #     p = k*r1 + a
            #     0 <= a < r1
            #
            #     q = l*r2 + b
            #     0 <= b < r2

            p_cell_lo = k * r1
            p_cell_hi = (k + 1) * r1 - 1

            q_cell_lo = l * r2
            q_cell_hi = (l + 1) * r2 - 1

            # Intersect with requested factor range.
            if p_cell_hi < FACTOR_MIN:
                continue

            if q_cell_hi < FACTOR_MIN:
                continue

            if p_cell_lo > FACTOR_MAX:
                continue

            if q_cell_lo > FACTOR_MAX:
                continue

            root_pairs += 1

            # Actual p/q integer ranges after factor-box restriction.
            p_lo = max(
                p_cell_lo,
                FACTOR_MIN,
            )

            p_hi = min(
                p_cell_hi,
                FACTOR_MAX,
            )

            q_lo = max(
                q_cell_lo,
                FACTOR_MIN,
            )

            q_hi = min(
                q_cell_hi,
                FACTOR_MAX,
            )

            if p_lo > p_hi or q_lo > q_hi:
                continue

            admissible_pairs += 1

            # -----------------------------------------------------------------
            # CARRY COMPATIBILITY
            # -----------------------------------------------------------------

            found_carry_for_pair = False

            # c3 is tiny. The theoretical bound here is small enough that
            # 0,1,2 is sufficient for this construction.
            for c3 in (0, 1, 2):

                carry_sum = E - c3

                if carry_sum < 0:
                    continue

                # c1 + c2 = E-c3
                #
                # 0 <= c1 < l
                # 0 <= c2 < k

                c1_lo = max(
                    0,
                    carry_sum - (k - 1),
                )

                c1_hi = min(
                    l - 1,
                    carry_sum,
                )

                if c1_lo > c1_hi:
                    continue

                for c1 in range(
                    c1_lo,
                    c1_hi + 1,
                ):

                    c2 = carry_sum - c1

                    if c2 < 0 or c2 >= k:
                        continue

                    found_carry_for_pair = True

                    # ---------------------------------------------------------
                    # RESIDUE INTERVALS
                    # ---------------------------------------------------------

                    (
                        a_lo,
                        a_hi,
                        b_lo,
                        b_hi,
                    ) = carry_intervals(
                        k,
                        l,
                        c1,
                        c2,
                        r1,
                        r2,
                    )

                    if a_lo > a_hi or b_lo > b_hi:
                        continue

                    # ---------------------------------------------------------
                    # CONVERT TO ACTUAL p/q RANGES
                    # ---------------------------------------------------------

                    carry_p_lo = k * r1 + a_lo
                    carry_p_hi = k * r1 + a_hi

                    carry_q_lo = l * r2 + b_lo
                    carry_q_hi = l * r2 + b_hi

                    # Intersect with original factor bounds.
                    carry_p_lo = max(
                        carry_p_lo,
                        p_lo,
                    )

                    carry_p_hi = min(
                        carry_p_hi,
                        p_hi,
                    )

                    carry_q_lo = max(
                        carry_q_lo,
                        q_lo,
                    )

                    carry_q_hi = min(
                        carry_q_hi,
                        q_hi,
                    )

                    if (
                        carry_p_lo > carry_p_hi
                        or carry_q_lo > carry_q_hi
                    ):
                        continue

                    carry_cells += 1

                    a_count = (
                        carry_p_hi
                        - carry_p_lo
                        + 1
                    )

                    b_count = (
                        carry_q_hi
                        - carry_q_lo
                        + 1
                    )

                    carry_a_states += a_count
                    carry_b_states += b_count

                    # ---------------------------------------------------------
                    # PRODUCT RANGE TEST
                    #
                    # p*q is monotone for positive p,q, so the rectangle
                    # contains n only if:
                    #
                    #     p_min*q_min <= n <= p_max*q_max
                    # ---------------------------------------------------------

                    product_lo = (
                        carry_p_lo
                        * carry_q_lo
                    )

                    product_hi = (
                        carry_p_hi
                        * carry_q_hi
                    )

                    if (
                        product_lo <= n
                        <= product_hi
                    ):
                        product_range_survivors += 1

                        # -----------------------------------------------------
                        # EXACT VALIDATION
                        #
                        # This scan is only used to measure how much work
                        # remains AFTER carry inversion.
                        # -----------------------------------------------------

                        for p in range(
                            carry_p_lo,
                            carry_p_hi + 1,
                        ):

                            if n % p != 0:
                                continue

                            q = n // p

                            if not (
                                carry_q_lo
                                <= q
                                <= carry_q_hi
                            ):
                                continue

                            if p * q != n:
                                continue

                            exact_pairs.add(
                                (p, q)
                            )

                    # ---------------------------------------------------------
                    # TRUE-CARRY VALIDATION
                    # ---------------------------------------------------------

                    if (
                        E == true_E
                        and
                        k == true_k
                        and
                        l == true_l
                        and
                        c1 == true_data["c1"]
                        and
                        c2 == true_data["c2"]
                        and
                        c3 == true_data["c3"]
                    ):
                        true_carry_seen = True

            if found_carry_for_pair:
                carry_compatible_pairs += 1

            if (
                E == true_E
                and
                k == true_k
                and
                l == true_l
            ):
                true_E_seen = True

    # -------------------------------------------------------------------------
    # RECOVERY
    # -------------------------------------------------------------------------

    recovered = (
        (p_true, q_true) in exact_pairs
        or
        (q_true, p_true) in exact_pairs
    )

    return {
        "n": n,

        "r1": r1,
        "r2": r2,

        "T": T,

        "true_E": true_E,
        "E_max": E_max,
        "true_K": true_K,

        "true_k": true_k,
        "true_l": true_l,

        "c1": true_data["c1"],
        "c2": true_data["c2"],
        "c3": true_data["c3"],

        "E_candidates": E_candidate_count,

        "root_pairs": root_pairs,
        "admissible_pairs": admissible_pairs,

        "carry_compatible_pairs":
            carry_compatible_pairs,

        "carry_cells":
            carry_cells,

        "carry_a_states":
            carry_a_states,

        "carry_b_states":
            carry_b_states,

        "product_range_survivors":
            product_range_survivors,

        "exact_count":
            len(exact_pairs),

        "recovered":
            recovered,

        "true_E_seen":
            true_E_seen,

        "true_carry_seen":
            true_carry_seen,
    }


# =============================================================================
# CORRELATION
# =============================================================================

def correlation(xs, ys):

    if len(xs) != len(ys):
        return 0.0

    if len(xs) < 2:
        return 0.0

    mx = statistics.mean(xs)
    my = statistics.mean(ys)

    numerator = sum(
        (x - mx) * (y - my)
        for x, y in zip(xs, ys)
    )

    dx = math.sqrt(
        sum(
            (x - mx) ** 2
            for x in xs
        )
    )

    dy = math.sqrt(
        sum(
            (y - my) ** 2
            for y in ys
        )
    )

    if dx == 0 or dy == 0:
        return 0.0

    return numerator / (dx * dy)


# =============================================================================
# MAIN
# =============================================================================

def run():

    total_start = time.perf_counter()

    print("=" * 92)
    print(
        "TWO-MODULUS E CARRY INVERSION / "
        "CARRY-CELL PRUNING EXPERIMENT"
    )
    print("=" * 92)

    print(
        f"N anchors                 = {N_ANCHORS:,}"
    )

    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )

    print(
        f"modulus range             = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )

    print(
        f"close ratio               = "
        f"{CLOSE_RATIO:.0%}"
    )

    print(
        f"seed                      = "
        f"{SEED:,}"
    )

    # =========================================================================
    # PRIME POOLS
    # =========================================================================

    print()
    print("=" * 92)
    print("BUILDING PRIME POOLS")
    print("=" * 92)

    factor_primes = [
        p
        for p in sieve_primes(FACTOR_MAX)
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p
        for p in sieve_primes(MODULUS_MAX)
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    # =========================================================================
    # ANCHORS
    # =========================================================================

    print()
    print("=" * 92)
    print("BUILDING ANCHORS")
    print("=" * 92)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print(
        f"actual anchors            = "
        f"{len(anchors):,}"
    )

    # =========================================================================
    # MODULUS PAIRS
    # =========================================================================

    print()
    print("=" * 92)
    print("BUILDING CLOSE MODULUS PAIRS")
    print("=" * 92)

    close_pairs = build_close_pairs(
        modulus_primes
    )

    print(
        f"close modulus pairs       = "
        f"{len(close_pairs):,}"
    )

    # =========================================================================
    # SEARCH
    # =========================================================================

    rng = random.Random(SEED + 1)

    results = []

    identity_failures = 0

    print()
    print("=" * 92)
    print("RUNNING CARRY-CELL PRUNING SEARCH")
    print("=" * 92)

    search_start = time.perf_counter()

    for index, (p, q) in enumerate(
        anchors,
        1,
    ):

        n = p * q

        pair = choose_modulus_pair(
            n,
            close_pairs,
            rng,
        )

        if pair is None:
            continue

        r1, r2 = pair

        result = run_anchor(
            p,
            q,
            r1,
            r2,
        )

        # ---------------------------------------------------------------------
        # EXACT IDENTITY VALIDATION
        # ---------------------------------------------------------------------

        lhs = (
            result["T"]
            -
            result["true_k"]
            * result["true_l"]
        )

        rhs = (
            result["c1"]
            + result["c2"]
            + result["c3"]
        )

        if lhs != result["true_E"]:
            identity_failures += 1

        if rhs != result["true_E"]:
            identity_failures += 1

        results.append(
            (
                p,
                q,
                result,
            )
        )

        if index % PROGRESS_EVERY == 0:
            print(
                f"anchor {index:3d}/{N_ANCHORS}"
            )

    search_runtime = (
        time.perf_counter()
        - search_start
    )

    if not results:
        raise RuntimeError(
            "No valid anchors were produced."
        )

    # =========================================================================
    # HELPERS
    # =========================================================================

    def avg(key):
        return statistics.mean(
            r[2][key]
            for r in results
        )

    recovery = sum(
        1
        for _, _, r in results
        if r["recovered"]
    )

    true_E_seen = sum(
        1
        for _, _, r in results
        if r["true_E_seen"]
    )

    true_carry_seen = sum(
        1
        for _, _, r in results
        if r["true_carry_seen"]
    )

    # =========================================================================
    # SUMMARY
    # =========================================================================

    print()
    print("=" * 92)
    print("SUMMARY")
    print("=" * 92)

    print(
        f"anchors analyzed             = "
        f"{len(results):,}"
    )

    print(
        f"identity failures            = "
        f"{identity_failures:,}"
    )

    print()

    print(
        f"average E-window candidates  = "
        f"{avg('E_candidates'):,.3f}"
    )

    print(
        f"average quotient pairs       = "
        f"{avg('root_pairs'):,.3f}"
    )

    print(
        f"average admissible pairs     = "
        f"{avg('admissible_pairs'):,.3f}"
    )

    print(
        f"average carry-compatible "
        f"pairs                       = "
        f"{avg('carry_compatible_pairs'):,.3f}"
    )

    print(
        f"average carry cells          = "
        f"{avg('carry_cells'):,.3f}"
    )

    print(
        f"average a states in cells    = "
        f"{avg('carry_a_states'):,.3f}"
    )

    print(
        f"average b states in cells    = "
        f"{avg('carry_b_states'):,.3f}"
    )

    print(
        f"average product-range "
        f"survivors                   = "
        f"{avg('product_range_survivors'):,.3f}"
    )

    # =========================================================================
    # PRUNING
    # =========================================================================

    print()
    print("=" * 92)
    print("CARRY PRUNING")
    print("=" * 92)

    avg_admissible = avg(
        "admissible_pairs"
    )

    avg_carry_pairs = avg(
        "carry_compatible_pairs"
    )

    if avg_admissible > 0:

        ratio = (
            avg_carry_pairs
            / avg_admissible
        )

        print(
            f"carry-compatible / "
            f"admissible              = "
            f"{ratio:.9f}"
        )

        print(
            f"carry pruning reduction    = "
            f"{(1.0 - ratio) * 100:.6f}%"
        )

    avg_carry_cells = avg(
        "carry_cells"
    )

    avg_product_survivors = avg(
        "product_range_survivors"
    )

    if avg_carry_cells > 0:

        print(
            f"product-range / "
            f"carry-cells              = "
            f"{avg_product_survivors / avg_carry_cells:.9f}"
        )

    # =========================================================================
    # RECOVERY
    # =========================================================================

    print()
    print("=" * 92)
    print("RECOVERY")
    print("=" * 92)

    print(
        f"correctly recovered            = "
        f"{recovery}/{len(results)}"
    )

    print(
        f"recovery rate                  = "
        f"{100.0 * recovery / len(results):.4f}%"
    )

    print(
        f"true E retained                = "
        f"{true_E_seen}/{len(results)}"
    )

    print(
        f"true carry retained            = "
        f"{true_carry_seen}/{len(results)}"
    )

    # =========================================================================
    # CORRELATIONS
    # =========================================================================

    true_E_values = [
        r[2]["true_E"]
        for r in results
    ]

    c1_values = [
        r[2]["c1"]
        for r in results
    ]

    c2_values = [
        r[2]["c2"]
        for r in results
    ]

    c3_values = [
        r[2]["c3"]
        for r in results
    ]

    k_values = [
        r[2]["true_k"]
        for r in results
    ]

    l_values = [
        r[2]["true_l"]
        for r in results
    ]

    print()
    print("=" * 92)
    print("TRUE CARRY STRUCTURE")
    print("=" * 92)

    print(
        f"average true E                = "
        f"{statistics.mean(true_E_values):.3f}"
    )

    print(
        f"average c1                   = "
        f"{statistics.mean(c1_values):.3f}"
    )

    print(
        f"average c2                   = "
        f"{statistics.mean(c2_values):.3f}"
    )

    print(
        f"average c3                   = "
        f"{statistics.mean(c3_values):.3f}"
    )

    print()

    print(
        f"corr(E,c1)                   = "
        f"{correlation(true_E_values, c1_values):.6f}"
    )

    print(
        f"corr(E,c2)                   = "
        f"{correlation(true_E_values, c2_values):.6f}"
    )

    print(
        f"corr(E,c3)                   = "
        f"{correlation(true_E_values, c3_values):.6f}"
    )

    print()

    print(
        f"corr(c1,l)                   = "
        f"{correlation(c1_values, l_values):.6f}"
    )

    print(
        f"corr(c2,k)                   = "
        f"{correlation(c2_values, k_values):.6f}"
    )

    # =========================================================================
    # EXAMPLES
    # =========================================================================

    print()
    print("=" * 92)
    print("EXAMPLES")
    print("=" * 92)

    for p, q, r in results[:20]:

        print(
            f"n={r['n']:,} "
            f"p={p:,} q={q:,} "
            f"mods=({r['r1']},{r['r2']})"
        )

        print(
            f"    T={r['T']} "
            f"E={r['true_E']} "
            f"K={r['true_K']} "
            f"Emax={r['E_max']}"
        )

        print(
            f"    true (k,l)=("
            f"{r['true_k']},"
            f"{r['true_l']}) "
            f"carry=("
            f"{r['c1']},"
            f"{r['c2']},"
            f"{r['c3']})"
        )

        print(
            f"    quotientPairs="
            f"{r['admissible_pairs']} "
            f"carryPairs="
            f"{r['carry_compatible_pairs']} "
            f"carryCells="
            f"{r['carry_cells']}"
        )

        print(
            f"    aStates="
            f"{r['carry_a_states']} "
            f"bStates="
            f"{r['carry_b_states']} "
            f"productRange="
            f"{r['product_range_survivors']} "
            f"exact="
            f"{r['exact_count']} "
            f"recovered="
            f"{r['recovered']}"
        )

    # =========================================================================
    # STRONGEST
    # =========================================================================

    def carry_ratio(entry):

        r = entry[2]

        if r["admissible_pairs"] == 0:
            return 1.0

        return (
            r["carry_compatible_pairs"]
            / r["admissible_pairs"]
        )

    ranked = sorted(
        results,
        key=carry_ratio,
    )

    print()
    print("=" * 92)
    print("STRONGEST CARRY PRUNING")
    print("=" * 92)

    for p, q, r in ranked[:20]:

        ratio = carry_ratio(
            (p, q, r)
        )

        print(
            f"n={r['n']:,} "
            f"p={p:,} q={q:,} "
            f"mods=({r['r1']},{r['r2']}) "
            f"pairs={r['admissible_pairs']} "
            f"carryPairs={r['carry_compatible_pairs']} "
            f"cells={r['carry_cells']} "
            f"ratio={ratio:.6f}"
        )

    # =========================================================================
    # WEAKEST
    # =========================================================================

    print()
    print("=" * 92)
    print("WEAKEST CARRY PRUNING")
    print("=" * 92)

    for p, q, r in ranked[-20:]:

        ratio = carry_ratio(
            (p, q, r)
        )

        print(
            f"n={r['n']:,} "
            f"p={p:,} q={q:,} "
            f"mods=({r['r1']},{r['r2']}) "
            f"pairs={r['admissible_pairs']} "
            f"carryPairs={r['carry_compatible_pairs']} "
            f"cells={r['carry_cells']} "
            f"ratio={ratio:.6f}"
        )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 92)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 92)

    print(
        r"""
The experiment begins with the exact two-modulus identity:

    p = a + k*r1
    q = b + l*r2

and:

    n =
        k*l*r1*r2
        + a*l*r2
        + b*k*r1
        + a*b.

Let:

    T = floor(n/(r1*r2))

and:

    E = T-k*l.

The carry decomposition is:

    a*l = c1*r1 + d1
    b*k = c2*r2 + d2

with:

    c1 = floor(a*l/r1)
    c2 = floor(b*k/r2).

The remaining terms produce:

    E = c1 + c2 + c3.

The previous question was whether this decomposition is merely
descriptive or whether it can actually prune candidate quotient
pairs.

This experiment tests that directly.

For every observable E:

    K = T-E

is factored:

    K = k*l.

Every admissible quotient pair must then satisfy:

    E = c1+c2+c3

with:

    0 <= c1 < l
    0 <= c2 < k.

If the carry split is valid, c1 and c2 produce intervals for the
unknown residues:

    ceil(c1*r1/l) <= a
                           <= floor(((c1+1)*r1-1)/l)

and:

    ceil(c2*r2/k) <= b
                           <= floor(((c2+1)*r2-1)/k).

Therefore each quotient pair is transformed into one or more
rectangles in the hidden residue plane (a,b).

The search hierarchy is:

    E window
        |
        v
    factor K=T-E
        |
        v
    admissible (k,l)
        |
        v
    carry compatibility
        |
        v
    residue rectangle
        |
        v
    product range
        |
        v
    exact factor.

The critical statistic is:

    carry-compatible pairs
        /
    admissible pairs.

If this is close to 1, then the carry decomposition gives almost
no additional information beyond the already-known E bound.

If it is substantially below 1, then the carry decomposition is
providing a genuine additional filter.

The next level is even more important.

For a surviving carry decomposition, the widths of the residue
intervals are approximately:

    r1/l
    r2/k.

Thus large k and l can make the hidden residue intervals narrow.

That gives a possible mechanism for recursive compression:

    unknown p,q
        ->
    unknown k,l
        ->
    carry restrictions
        ->
    narrow a,b
        ->
    narrow p,q.

However, this experiment deliberately does NOT claim that this has
been achieved.

The exact integer scan remaining inside the carry rectangle is
validation work. It is not included as evidence of candidate
generation efficiency.

The important question is whether the carry rectangles themselves
become small enough that their generation is substantially cheaper
than searching the original factor interval.

The final condition remains exact:

    p*q == n.

Only exact factor pairs count as recovery.
"""
    )

    # =========================================================================
    # TIMING
    # =========================================================================

    print()
    print("=" * 92)
    print("TIMING")
    print("=" * 92)

    print(
        f"search runtime                = "
        f"{search_runtime:.3f} s"
    )

    print(
        f"total runtime                 = "
        f"{time.perf_counter() - total_start:.3f} s"
    )

    print()
    print("=" * 92)
    print("EXPERIMENT COMPLETE")
    print("=" * 92)


if __name__ == "__main__":
    run()