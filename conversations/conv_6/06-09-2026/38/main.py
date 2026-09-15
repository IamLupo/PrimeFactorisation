#!/usr/bin/env python3

import math
import time


# ============================================================================
# START EXPERIMENT 194
#
# SECOND-ORDER QR SIEVE COST MODEL
#
# Experiment 193 found that an additive model
#
#     T(S) = sum(cost[p]) + rho(S) * N * test_cost
#
# is not accurate enough.
#
# The reason is that applying a prime to an already-partially-sieved
# bytearray has a context-dependent cost.
#
# Experiment 194 therefore measures:
#
#     marginal_cost(S, p)
#
# directly:
#
#     cost(S + p) - cost(S)
#
# and builds an empirical SECOND-ORDER model.
#
# We use the measured time of actual sieve passes on a representative range.
#
# For every subset S encountered during optimization:
#
#     sieve_time(S)
#     survivor_density(S)
#
# are measured/cached.
#
# Then we use dynamic construction:
#
#     T_model(S)
#         = measured_sieve_time(S)
#         + expected_survivors(S) * fermat_test_cost
#
# The search is still over 2^16 possible subsets, but actual sieve
# measurements are only performed for promising subsets / incremental
# states rather than blindly assuming independent costs.
#
# We additionally measure the important quantity:
#
#     marginal cost of p after S
#
# to determine whether prime interactions are significant.
#
# ============================================================================


PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59
]

BLOCK_SIZE = 1 << 20

# Measurement range.
MEASURE_RANGE = 1 << 22

# Actual Fermat search range.
MAX_Y = 25_000_000

# Number of top predicted subsets to verify exactly.
TOP_VERIFY = 20

# Candidate subsets retained after each beam layer.
BEAM_WIDTH = 250

# Number of random subsets used to calibrate the model.
RANDOM_CALIBRATION = 100

# Maximum number of exact sieve measurements.
MAX_EXACT_MEASUREMENTS = 2500

# Existing instances.
INSTANCES = {
    48: (8390069, 33547589),
    54: (124517461, 144517463),
    60: (1058841403, 1088841421),
    66: (8569934017, 8609934041),
}


# ============================================================================
# LEGENDRE SYMBOL / ANALYTICAL DENSITY
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:

    a %= p

    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1

    if value == p - 1:
        return -1

    raise RuntimeError("invalid Legendre symbol")


def prime_density(n: int, p: int) -> float:

    chi = legendre_symbol(-n, p)

    if n % p == 0:
        allowed = (p + 1) // 2
    else:
        allowed = (p + chi) // 2

    return allowed / p


def subset_density(n: int, subset) -> float:

    result = 1.0

    for p in subset:
        result *= prime_density(n, p)

    return result


# ============================================================================
# QUADRATIC-RESIDUE BAD RESIDUES
# ============================================================================

def build_bad_residues(n: int, p: int):

    qr = {
        (x * x) % p
        for x in range(p)
    }

    return [
        r
        for r in range(p)
        if (n + r * r) % p not in qr
    ]


# ============================================================================
# SIEVE STATE
# ============================================================================

class PrimeData:

    def __init__(self, n: int):

        self.n = n

        self.bad = {}

        for p in PRIMES:
            self.bad[p] = build_bad_residues(
                n,
                p
            )


# ============================================================================
# MEASURE ONE SUBSET
# ============================================================================

def measure_sieve(
    prime_data: PrimeData,
    subset,
    limit: int,
    block_size: int,
):

    start = time.perf_counter()

    survivors = 0

    block_start = 0

    while block_start < limit:

        block_end = min(
            block_start + block_size,
            limit
        )

        block_len = block_end - block_start

        alive = bytearray(
            b"\x01"
        ) * block_len

        for p in subset:

            for residue in prime_data.bad[p]:

                first = (
                    block_start
                    + ((residue - block_start) % p)
                )

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // p
                ) + 1

                alive[idx:block_len:p] = (
                    b"\x00"
                ) * count

        survivors += alive.count(1)

        block_start = block_end

    elapsed = time.perf_counter() - start

    density = survivors / limit

    return elapsed, survivors, density


# ============================================================================
# FERMAT TEST COST
# ============================================================================

def measure_test_cost(n: int):

    samples = min(
        MEASURE_RANGE,
        1_000_000
    )

    start = time.perf_counter()

    hits = 0

    for y in range(samples):

        value = n + y * y

        x = math.isqrt(value)

        if x * x == value:
            hits += 1

    elapsed = time.perf_counter() - start

    return elapsed / samples, hits


# ============================================================================
# EXACT FERMAT RECOVERY
# ============================================================================

def segmented_fermat(
    n: int,
    subset,
    max_y: int,
    block_size: int,
    prime_data: PrimeData,
):

    start = time.perf_counter()

    survivors = 0

    block_start = 0

    while block_start <= max_y:

        block_end = min(
            block_start + block_size,
            max_y + 1
        )

        block_len = block_end - block_start

        alive = bytearray(
            b"\x01"
        ) * block_len

        for p in subset:

            for residue in prime_data.bad[p]:

                first = (
                    block_start
                    + ((residue - block_start) % p)
                )

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // p
                ) + 1

                alive[idx:block_len:p] = (
                    b"\x00"
                ) * count

        pos = 0

        while True:

            pos = alive.find(1, pos)

            if pos < 0:
                break

            y = block_start + pos

            survivors += 1

            value = n + y * y

            x = math.isqrt(value)

            if x * x == value:

                p = x - y
                q = x + y

                if p > 1 and q > 1 and p * q == n:

                    elapsed = time.perf_counter() - start

                    return (
                        y,
                        (p, q),
                        survivors,
                        elapsed,
                    )

            pos += 1

        block_start = block_end

    elapsed = time.perf_counter() - start

    return (
        None,
        None,
        survivors,
        elapsed,
    )


# ============================================================================
# BITMASK UTILITIES
# ============================================================================

def tuple_to_mask(subset):

    mask = 0

    for p in subset:

        i = PRIMES.index(p)

        mask |= 1 << i

    return mask


def mask_to_tuple(mask):

    return tuple(
        PRIMES[i]
        for i in range(len(PRIMES))
        if mask & (1 << i)
    )


# ============================================================================
# FIRST-ORDER ANALYTICAL MODEL
# ============================================================================

def analytical_model(
    n: int,
    subset,
    y_limit: int,
    test_cost: float,
    individual_costs,
):

    rho = subset_density(
        n,
        subset
    )

    sieve_cost = sum(
        individual_costs[p]
        for p in subset
    )

    predicted_survivors = (
        y_limit + 1
    ) * rho

    predicted_total = (
        sieve_cost
        + predicted_survivors * test_cost
    )

    return predicted_total, rho


# ============================================================================
# RANDOM SUBSET CALIBRATION
# ============================================================================

def random_mask():

    import random

    return random.randrange(
        1 << len(PRIMES)
    )


def calibrate_model(
    n: int,
    prime_data: PrimeData,
    test_cost: float,
):

    print()
    print("=" * 72)
    print("CALIBRATION")
    print("=" * 72)

    individual_costs = {}

    # ------------------------------------------------------------
    # Individual prime costs.
    # ------------------------------------------------------------

    for p in PRIMES:

        subset = (p,)

        elapsed, survivors, density = (
            measure_sieve(
                prime_data,
                subset,
                MEASURE_RANGE,
                BLOCK_SIZE,
            )
        )

        individual_costs[p] = elapsed

        print(
            f"prime={p:2d} "
            f"density={density:.12f} "
            f"time={elapsed:.8f}s"
        )

    # ------------------------------------------------------------
    # Random subset calibration.
    # ------------------------------------------------------------

    import random

    random_masks = set()

    while len(random_masks) < RANDOM_CALIBRATION:

        random_masks.add(
            random_mask()
        )

    calibration = []

    for mask in random_masks:

        subset = mask_to_tuple(mask)

        actual_time, survivors, density = (
            measure_sieve(
                prime_data,
                subset,
                MEASURE_RANGE,
                BLOCK_SIZE,
            )
        )

        additive_time = sum(
            individual_costs[p]
            for p in subset
        )

        # Interaction term:
        #
        # actual - additive
        #
        interaction = (
            actual_time
            - additive_time
        )

        calibration.append(
            (
                subset,
                actual_time,
                additive_time,
                interaction,
                density,
                survivors,
            )
        )

    interactions = [
        x[3]
        for x in calibration
    ]

    mean_interaction = (
        sum(interactions) / len(interactions)
        if interactions
        else 0.0
    )

    print()
    print("SECOND-ORDER INTERACTION STATISTICS")
    print(
        f"    samples = "
        f"{len(calibration)}"
    )
    print(
        f"    mean interaction = "
        f"{mean_interaction:.8f}s"
    )
    print(
        f"    min interaction = "
        f"{min(interactions):.8f}s"
    )
    print(
        f"    max interaction = "
        f"{max(interactions):.8f}s"
    )

    return individual_costs, calibration


# ============================================================================
# SECOND-ORDER REGRESSION
# ============================================================================

def fit_pairwise_model(
    calibration,
    n: int,
):

    """
    Fit:

        T(S) ≈
            intercept
            + sum_i c_i
            + sum_{i<j} c_ij

    using ordinary least squares.

    We solve this ourselves with normal equations and Gaussian elimination.
    """

    pairs = []

    for i in range(len(PRIMES)):

        for j in range(i + 1, len(PRIMES)):

            pairs.append(
                (i, j)
            )

    feature_count = (
        1
        + len(PRIMES)
        + len(pairs)
    )

    # X^T X
    XTX = [
        [0.0 for _ in range(feature_count)]
        for _ in range(feature_count)
    ]

    XTy = [
        0.0
        for _ in range(feature_count)
    ]

    for subset, actual_time, _, _, _, _ in calibration:

        mask = tuple_to_mask(subset)

        x = [0.0] * feature_count

        x[0] = 1.0

        for i in range(len(PRIMES)):

            if mask & (1 << i):

                x[1 + i] = 1.0

        offset = 1 + len(PRIMES)

        for k, (i, j) in enumerate(pairs):

            if (
                mask & (1 << i)
                and mask & (1 << j)
            ):

                x[offset + k] = 1.0

        for i in range(feature_count):

            XTy[i] += (
                x[i] * actual_time
            )

            for j in range(feature_count):

                XTX[i][j] += (
                    x[i] * x[j]
                )

    # ------------------------------------------------------------
    # Gaussian elimination.
    # ------------------------------------------------------------

    A = [
        XTX[i] + [XTy[i]]
        for i in range(feature_count)
    ]

    for col in range(feature_count):

        pivot = max(
            range(col, feature_count),
            key=lambda r: abs(A[r][col])
        )

        if abs(A[pivot][col]) < 1e-14:

            continue

        if pivot != col:

            A[col], A[pivot] = (
                A[pivot],
                A[col]
            )

        pivot_value = A[col][col]

        for j in range(col, feature_count + 1):

            A[col][j] /= pivot_value

        for r in range(feature_count):

            if r == col:
                continue

            factor = A[r][col]

            if factor == 0:
                continue

            for j in range(col, feature_count + 1):

                A[r][j] -= (
                    factor * A[col][j]
                )

    coefficients = [
        A[i][-1]
        for i in range(feature_count)
    ]

    return coefficients, pairs


# ============================================================================
# PREDICT WITH SECOND-ORDER MODEL
# ============================================================================

def pairwise_prediction(
    subset,
    coefficients,
    pairs,
):

    mask = tuple_to_mask(subset)

    value = coefficients[0]

    for i in range(len(PRIMES)):

        if mask & (1 << i):

            value += coefficients[1 + i]

    offset = 1 + len(PRIMES)

    for k, (i, j) in enumerate(pairs):

        if (
            mask & (1 << i)
            and mask & (1 << j)
        ):

            value += coefficients[offset + k]

    return value


# ============================================================================
# BEAM SEARCH FOR BEST SUBSETS
# ============================================================================

def beam_search(
    n: int,
    test_cost: float,
    prime_data: PrimeData,
    coefficients,
    pairs,
):

    """
    Build subsets incrementally.

    At each cardinality retain BEAM_WIDTH candidates.

    Score:

        second-order sieve-time prediction
        + expected Fermat-test work.
    """

    current = {
        tuple()
    }

    all_candidates = []

    for depth in range(
        1,
        len(PRIMES) + 1
    ):

        next_set = set()

        for subset in current:

            used = set(subset)

            for p in PRIMES:

                if p in used:
                    continue

                candidate = tuple(
                    sorted(
                        subset + (p,)
                    )
                )

                next_set.add(candidate)

        scored = []

        for subset in next_set:

            sieve_prediction = pairwise_prediction(
                subset,
                coefficients,
                pairs,
            )

            rho = subset_density(
                n,
                subset
            )

            expected_survivors = (
                MAX_Y + 1
            ) * rho

            total_prediction = (
                sieve_prediction
                + expected_survivors * test_cost
            )

            scored.append(
                (
                    total_prediction,
                    subset,
                    rho,
                )
            )

        scored.sort(
            key=lambda x: x[0]
        )

        current = {
            x[1]
            for x in scored[:BEAM_WIDTH]
        }

        all_candidates.extend(
            scored[:BEAM_WIDTH]
        )

    all_candidates.sort(
        key=lambda x: x[0]
    )

    # Remove duplicate subsets.
    seen = set()
    unique = []

    for item in all_candidates:

        subset = item[1]

        if subset in seen:
            continue

        seen.add(subset)
        unique.append(item)

    return unique


# ============================================================================
# EXACT VERIFICATION
# ============================================================================

def verify_candidates(
    n: int,
    p: int,
    q: int,
    candidates,
    prime_data: PrimeData,
):

    verified = []

    for rank, (_, subset, predicted_rho) in enumerate(
        candidates[:TOP_VERIFY],
        start=1,
    ):

        start = time.perf_counter()

        found_y, result, survivors, actual_time = (
            segmented_fermat(
                n,
                subset,
                MAX_Y,
                BLOCK_SIZE,
                prime_data,
            )
        )

        verification_time = (
            time.perf_counter()
            - start
        )

        correct = (
            result == (p, q)
            or result == (q, p)
        )

        verified.append(
            (
                actual_time,
                subset,
                survivors,
                found_y,
                correct,
                predicted_rho,
                verification_time,
            )
        )

        print()
        print(f"RANK {rank}")
        print(f"    subset = {list(subset)}")
        print(
            f"    predicted density = "
            f"{predicted_rho:.12f}"
        )
        print(
            f"    actual survivors = "
            f"{survivors}"
        )
        print(
            f"    found y = "
            f"{found_y}"
        )
        print(
            f"    result = "
            f"{result}"
        )
        print(
            f"    correct = "
            f"{correct}"
        )
        print(
            f"    actual time = "
            f"{actual_time:.8f}s"
        )
        print(
            f"    verification wall time = "
            f"{verification_time:.8f}s"
        )

    verified.sort(
        key=lambda x: x[0]
    )

    return verified


# ============================================================================
# ANALYZE INSTANCE
# ============================================================================

def run_instance(bits: int, p: int, q: int):

    n = p * q

    true_y = (q - p) // 2

    print()
    print("=" * 72)
    print(f"START INSTANCE {bits}-BIT")
    print("=" * 72)

    print(f"p = {p}")
    print(f"q = {q}")
    print(f"n = {n}")
    print(f"true y = {true_y}")

    if true_y > MAX_Y:

        print()
        print(
            "SKIPPED: true y exceeds MAX_Y"
        )

        print(
            f"FINISHED INSTANCE {bits}-BIT"
        )

        return

    prime_data = PrimeData(n)

    test_cost, accidental = (
        measure_test_cost(n)
    )

    print()
    print("FERMAT TEST COST")
    print(
        f"    cost = "
        f"{test_cost:.12e}s"
    )
    print(
        f"    accidental squares = "
        f"{accidental}"
    )

    # ------------------------------------------------------------------------
    # Calibration
    # ------------------------------------------------------------------------

    individual_costs, calibration = (
        calibrate_model(
            n,
            prime_data,
            test_cost,
        )
    )

    # ------------------------------------------------------------------------
    # Fit pairwise interaction model
    # ------------------------------------------------------------------------

    print()
    print("=" * 72)
    print("PAIRWISE MODEL")
    print("=" * 72)

    coefficients, pairs = (
        fit_pairwise_model(
            calibration,
            n,
        )
    )

    print(
        f"    coefficients = "
        f"{len(coefficients)}"
    )

    # Calculate training residual.
    residuals = []

    for subset, actual, _, _, _, _ in calibration:

        predicted = pairwise_prediction(
            subset,
            coefficients,
            pairs,
        )

        residuals.append(
            actual - predicted
        )

    mse = (
        sum(x * x for x in residuals)
        / len(residuals)
    )

    rmse = math.sqrt(mse)

    print(
        f"    training RMSE = "
        f"{rmse:.10e}s"
    )

    # ------------------------------------------------------------------------
    # Beam search
    # ------------------------------------------------------------------------

    print()
    print("=" * 72)
    print("BEAM SEARCH")
    print("=" * 72)

    candidates = beam_search(
        n,
        test_cost,
        prime_data,
        coefficients,
        pairs,
    )

    print(
        f"    generated candidate states = "
        f"{len(candidates)}"
    )

    print()
    print("TOP PREDICTED CANDIDATES")

    for rank, item in enumerate(
        candidates[:15],
        start=1,
    ):

        prediction, subset, density = item

        print()
        print(f"RANK {rank}")
        print(f"    primes = {list(subset)}")
        print(
            f"    predicted density = "
            f"{density:.12f}"
        )
        print(
            f"    predicted total = "
            f"{prediction:.8f}s"
        )

    # ------------------------------------------------------------------------
    # Exact verification
    # ------------------------------------------------------------------------

    print()
    print("=" * 72)
    print("EXACT VERIFICATION")
    print("=" * 72)

    verified = verify_candidates(
        n,
        p,
        q,
        candidates,
        prime_data,
    )

    if verified:

        best = verified[0]

        print()
        print("BEST ACTUAL CANDIDATE")
        print(
            f"    primes = "
            f"{list(best[1])}"
        )
        print(
            f"    time = "
            f"{best[0]:.8f}s"
        )
        print(
            f"    survivors = "
            f"{best[2]}"
        )
        print(
            f"    found y = "
            f"{best[3]}"
        )
        print(
            f"    correct = "
            f"{best[4]}"
        )

    print()
    print(f"FINISHED INSTANCE {bits}-BIT")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("START EXPERIMENT 194")
    print()
    print("Second-order empirical QR sieve optimization")
    print()
    print(f"PRIMES = {PRIMES}")
    print(f"BLOCK_SIZE = {BLOCK_SIZE}")
    print(f"MEASURE_RANGE = {MEASURE_RANGE}")
    print(f"MAX_Y = {MAX_Y}")
    print(f"TOP_VERIFY = {TOP_VERIFY}")
    print(f"BEAM_WIDTH = {BEAM_WIDTH}")
    print(f"RANDOM_CALIBRATION = {RANDOM_CALIBRATION}")

    # Start with 48-bit and 54-bit because these have reasonably quick
    # recovery times while still giving enough survivors for meaningful
    # model fitting.

    for bits in [48, 54]:

        p, q = INSTANCES[bits]

        run_instance(
            bits,
            p,
            q,
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 194")
    print("=" * 72)


if __name__ == "__main__":
    main()
