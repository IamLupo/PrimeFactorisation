#!/usr/bin/env python3

"""
START EXPERIMENT 185

SEGMENTED PRIME FACTOR SEARCH

Experiments 180-184 established:

    1. C-derived relations are tautological.
    2. Small-radix residue classes do not select the factor.
    3. q ≡ n*p^{-1} (mod t) gives no independent sieve condition
       when t does not divide n.
    4. Naive modular sieving in Python is slower than plain trial
       division because of interpreter overhead.

Experiment 185 therefore tests the most direct computational use of
the prime structure:

    Instead of testing every odd integer p <= sqrt(n),
    generate only prime p values with a segmented sieve.

Then test:

    n % p == 0

only for prime p.

This is still trial division. The experiment is NOT claiming a new
factorization algorithm.

The purpose is to establish a clean optimized baseline for all future
experiments.

Three methods are compared:

    A) odd integer trial division
    B) segmented prime trial division
    C) segmented prime trial division with the project radices
       explicitly removed from the candidate stream

Method C should be identical in mathematical candidate density to B,
because the R1/R2 values are already prime candidates. This tests
whether the special radix structure provides anything computationally
once prime-only enumeration is used.

No C-values are used.

No factor residues are assumed.

The hidden factors are used only for verification.

"""


import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

# Number of integers handled by one segmented-sieve block.
SEGMENT_SIZE = 1_000_000

# Safety cap.
MAX_DIVISION_TESTS = 20_000_000


# ============================================================
# PRIME TEST
# ============================================================

def is_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37,
    ]

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = 41
    step = 2

    while d * d <= n:

        if n % d == 0:
            return False

        d += step
        step = 6 - step

    return True


# ============================================================
# SEMIPRIME GENERATION
# ============================================================

def make_semiprime(bits):

    low = 1 << (
        bits // 2 - 1
    )

    high = 1 << (
        bits // 2 + 1
    )

    for p in range(
        low | 1,
        high,
        2,
    ):

        if not is_prime(p):
            continue

        target = (
            (1 << bits)
            // p
        )

        for delta in range(
            -1000,
            1001,
            2,
        ):

            q = target + delta

            if q <= p:
                continue

            if not is_prime(q):
                continue

            n = p * q

            if n.bit_length() == bits:
                return p, q, n

    raise RuntimeError(
        f"Could not generate {bits}-bit semiprime"
    )


# ============================================================
# BASE PRIME GENERATION
# ============================================================

def small_primes(limit):

    """
    Ordinary sieve for primes <= limit.

    limit is only sqrt(SEGMENT_SIZE), so this is tiny.
    """

    if limit < 2:
        return []

    sieve = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

    sieve[0:2] = (
        b"\x00\x00"
    )

    root = math.isqrt(
        limit
    )

    for p in range(
        2,
        root + 1,
    ):

        if not sieve[p]:
            continue

        start = p * p

        step_count = (
            (limit - start)
            // p
            + 1
        )

        sieve[
            start:
            limit + 1:
            p
        ] = b"\x00" * step_count

    return [
        i
        for i in range(
            2,
            limit + 1,
        )
        if sieve[i]
    ]


# ============================================================
# SEGMENTED PRIME GENERATOR
# ============================================================

def segmented_primes(
    limit,
):
    """
    Yield every prime <= limit using a segmented sieve.
    """

    if limit < 2:
        return

    root = math.isqrt(
        limit
    )

    base_primes = small_primes(
        root
    )

    # Yield 2 separately.
    yield 2

    low = 3

    if low > limit:
        return

    # Keep only odd numbers in the segment.
    while low <= limit:

        high = min(
            low + 2 * SEGMENT_SIZE - 2,
            limit,
        )

        # Number of odd integers:
        #
        # low, low+2, ..., high
        #
        length = (
            (high - low)
            // 2
            + 1
        )

        segment = bytearray(
            b"\x01"
        ) * length

        for p in base_primes:

            if p == 2:
                continue

            p2 = p * p

            if p2 > high:
                break

            # First multiple >= low.
            first = max(
                p2,
                (
                    (low + p - 1)
                    // p
                ) * p,
            )

            # Make first odd.
            if first % 2 == 0:
                first += p

            if first > high:
                continue

            index = (
                first - low
            ) // 2

            step = p

            count = (
                (length - 1 - index)
                // step
                + 1
            )

            segment[
                index:
                length:
                step
            ] = b"\x00" * count

        for i, ok in enumerate(
            segment
        ):

            if ok:

                yield (
                    low + 2 * i
                )

        low = high + 2


# ============================================================
# BASELINE ODD SCAN
# ============================================================

def odd_trial_division(n):

    limit = math.isqrt(
        n
    )

    tests = 0

    start = time.perf_counter()

    # Check 2 separately.
    if n % 2 == 0:

        return {
            "found": 2,
            "tests": 1,
            "time":
                time.perf_counter()
                - start,
            "aborted": None,
        }

    for p in range(
        3,
        limit + 1,
        2,
    ):

        tests += 1

        if (
            tests
            > MAX_DIVISION_TESTS
        ):

            return {
                "found": None,
                "tests": tests,
                "time":
                    time.perf_counter()
                    - start,
                "aborted":
                    "MAX_DIVISION_TESTS",
            }

        if n % p == 0:

            return {
                "found": p,
                "tests": tests,
                "time":
                    time.perf_counter()
                    - start,
                "aborted": None,
            }

    return {
        "found": None,
        "tests": tests,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# SEGMENTED PRIME TRIAL DIVISION
# ============================================================

def segmented_prime_trial(
    n,
    excluded_radices=None,
):

    limit = math.isqrt(
        n
    )

    if excluded_radices is None:
        excluded_radices = set()
    else:
        excluded_radices = set(
            excluded_radices
        )

    tests = 0

    generated_primes = 0

    skipped_radix_primes = 0

    start = time.perf_counter()

    for p in segmented_primes(
        limit
    ):

        generated_primes += 1

        # We don't need to skip a radix prime if it is itself a
        # possible factor. Check n first.
        if p in excluded_radices:

            if n % p == 0:

                return {
                    "found": p,
                    "tests": tests + 1,
                    "generated_primes":
                        generated_primes,
                    "skipped_radix_primes":
                        skipped_radix_primes,
                    "time":
                        time.perf_counter()
                        - start,
                    "aborted": None,
                }

            skipped_radix_primes += 1
            continue

        tests += 1

        if (
            tests
            > MAX_DIVISION_TESTS
        ):

            return {
                "found": None,
                "tests": tests,
                "generated_primes":
                    generated_primes,
                "skipped_radix_primes":
                    skipped_radix_primes,
                "time":
                    time.perf_counter()
                    - start,
                "aborted":
                    "MAX_DIVISION_TESTS",
            }

        if n % p == 0:

            return {
                "found": p,
                "tests": tests,
                "generated_primes":
                    generated_primes,
                "skipped_radix_primes":
                    skipped_radix_primes,
                "time":
                    time.perf_counter()
                    - start,
                "aborted": None,
            }

    return {
        "found": None,
        "tests": tests,
        "generated_primes":
            generated_primes,
        "skipped_radix_primes":
            skipped_radix_primes,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# PRIME COUNT ESTIMATE
# ============================================================

def prime_density_estimate(
    limit,
):

    if limit <= 1:
        return 0.0

    return (
        limit
        / math.log(limit)
    )


# ============================================================
# RUN INSTANCE
# ============================================================

def run_instance(bits):

    print()
    print("=" * 72)
    print(
        f"START INSTANCE {bits}-BIT"
    )
    print("=" * 72)

    p, q, n = make_semiprime(
        bits
    )

    limit = math.isqrt(
        n
    )

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"sqrt(n) = {limit}"
    )

    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    print()

    # --------------------------------------------------------
    # Expected candidate counts.
    # --------------------------------------------------------

    odd_count = (
        (limit - 1)
        // 2
    )

    estimate = prime_density_estimate(
        limit
    )

    print(
        "CANDIDATE COUNTS"
    )

    print(
        f"    odd integers <= sqrt(n) = "
        f"{odd_count}"
    )

    print(
        f"    rough pi(sqrt(n)) estimate = "
        f"{estimate:.0f}"
    )

    print()

    # --------------------------------------------------------
    # Baseline odd scan.
    # --------------------------------------------------------

    print(
        "METHOD A: ODD TRIAL DIVISION"
    )

    baseline = odd_trial_division(
        n
    )

    print(
        f"    tests = "
        f"{baseline['tests']}"
    )

    print(
        f"    time = "
        f"{baseline['time']:.6f}s"
    )

    print(
        f"    found = "
        f"{baseline['found']}"
    )

    if baseline["aborted"]:

        print(
            f"    status = "
            f"{baseline['aborted']}"
        )

    # --------------------------------------------------------
    # Segmented prime search.
    # --------------------------------------------------------

    print()
    print(
        "METHOD B: SEGMENTED PRIME TRIAL DIVISION"
    )

    prime_result = segmented_prime_trial(
        n
    )

    print(
        f"    generated primes = "
        f"{prime_result['generated_primes']}"
    )

    print(
        f"    n % p tests = "
        f"{prime_result['tests']}"
    )

    print(
        f"    time = "
        f"{prime_result['time']:.6f}s"
    )

    print(
        f"    found = "
        f"{prime_result['found']}"
    )

    if prime_result["aborted"]:

        print(
            f"    status = "
            f"{prime_result['aborted']}"
        )

    # --------------------------------------------------------
    # Segmented prime search with project radices removed.
    #
    # This should provide almost no benefit because there are
    # only ten distinct radix primes.
    # --------------------------------------------------------

    print()
    print(
        "METHOD C: PRIME SEARCH WITH R1/R2 RADICES REMOVED"
    )

    radix_set = set(
        R1 + R2
    )

    radix_result = segmented_prime_trial(
        n,
        excluded_radices=radix_set,
    )

    print(
        f"    generated primes = "
        f"{radix_result['generated_primes']}"
    )

    print(
        f"    n % p tests = "
        f"{radix_result['tests']}"
    )

    print(
        f"    skipped radix primes = "
        f"{radix_result['skipped_radix_primes']}"
    )

    print(
        f"    time = "
        f"{radix_result['time']:.6f}s"
    )

    print(
        f"    found = "
        f"{radix_result['found']}"
    )

    if radix_result["aborted"]:

        print(
            f"    status = "
            f"{radix_result['aborted']}"
        )

    # --------------------------------------------------------
    # Verification.
    # --------------------------------------------------------

    print()
    print(
        "VERIFICATION"
    )

    print(
        f"    hidden p prime = "
        f"{is_prime(p)}"
    )

    print(
        f"    hidden q prime = "
        f"{is_prime(q)}"
    )

    print(
        f"    Method A correct = "
        f"{baseline['found'] == p}"
    )

    print(
        f"    Method B correct = "
        f"{prime_result['found'] == p}"
    )

    print(
        f"    Method C correct = "
        f"{radix_result['found'] == p}"
    )

    if prime_result["found"] is not None:

        print()
        print(
            "SPEEDUP B/A"
        )

        if (
            prime_result["time"]
            > 0
        ):

            print(
                f"    "
                f"{baseline['time'] / prime_result['time']:.3f}x"
            )

    if radix_result["found"] is not None:

        print()
        print(
            "SPEEDUP C/A"
        )

        if (
            radix_result["time"]
            > 0
        ):

            print(
                f"    "
                f"{baseline['time'] / radix_result['time']:.3f}x"
            )

    print()
    print(
        f"FINISHED INSTANCE {bits}-BIT"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print(
        "START EXPERIMENT 185"
    )

    print()

    print(
        f"R1 = {R1}"
    )

    print(
        f"R2 = {R2}"
    )

    print(
        f"BITS = {BITS}"
    )

    print(
        f"SEGMENT_SIZE = "
        f"{SEGMENT_SIZE}"
    )

    print(
        f"MAX_DIVISION_TESTS = "
        f"{MAX_DIVISION_TESTS}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 185"
    )


if __name__ == "__main__":
    main()
