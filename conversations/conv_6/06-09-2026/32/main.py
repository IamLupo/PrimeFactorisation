#!/usr/bin/env python3

"""
START EXPERIMENT 188

FAST GAP-WHEEL FERMAT SEARCH

Experiment 187 found that modular difference-of-squares filtering can
make the Fermat gap search substantially sparser.

However, the implementation was inefficient:

    next gap = min(all current arithmetic-progressions)

which requires scanning every wheel residue class for every candidate.

This experiment fixes the iterator.

For

    n = p*q

define

    d = q-p
    s = p+q

then

    s^2-d^2 = 4n.

For every modulus m we require:

    s^2-d^2 = 4n (mod m).

For a given m, determine the allowed d residues:

    exists s:
        s^2 = d^2 + 4n (mod m).

The allowed residues from several coprime moduli are combined into a
single wheel:

    d = r (mod M)

where

    M = product(moduli).

Instead of maintaining one "current value" per residue and repeatedly
calling min(), we:

    1. sort all admissible residues in one period;
    2. enumerate that sorted list;
    3. advance by exactly M to the next period.

Because d for odd p,q is even, parity is handled directly by doubling
the wheel period.

The experiment compares:

    A) ordinary Fermat
    B) fast modular gap wheel

and records:

    * wheel density
    * number of admissible gaps before the true gap
    * number of square tests
    * wall-clock time
    * speedup

We deliberately test more primes than Experiment 187.

No C-values are used.

No factor residues are supplied.

The hidden p,q are used only for verification.

"""


import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

FILTER_PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31
]

WHEEL_DEPTHS = [
    1, 2, 3, 4, 5, 6, 7
]

MAX_GAP_TESTS = 50_000_000

MAX_WHEEL_CLASSES = 5_000_000


# ============================================================
# PRIME TEST
# ============================================================

def is_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
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
        "semiprime generation failed"
    )


# ============================================================
# GENERALIZED CRT
# ============================================================

def egcd(a, b):

    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(
        b,
        a % b,
    )

    return (
        g,
        y1,
        x1 - (a // b) * y1,
    )


def inv_mod(a, m):

    g, x, _ = egcd(
        a,
        m,
    )

    if g != 1:
        raise ValueError(
            f"{a} not invertible modulo {m}"
        )

    return x % m


def crt_pair(
    a1,
    m1,
    a2,
    m2,
):

    g = math.gcd(
        m1,
        m2,
    )

    if (
        a2 - a1
    ) % g != 0:
        return None

    m1r = m1 // g
    m2r = m2 // g

    rhs = (
        (a2 - a1)
        // g
    ) % m2r

    if m2r == 1:

        t = 0

    else:

        t = (
            rhs
            * inv_mod(
                m1r % m2r,
                m2r,
            )
        ) % m2r

    x = (
        a1
        + m1 * t
    )

    modulus = (
        m1
        * m2r
    )

    return (
        x % modulus,
        modulus,
    )


# ============================================================
# ALLOWED GAP RESIDUES MOD m
# ============================================================

def allowed_gap_residues(
    n,
    m,
):
    """
    Return d residues modulo m for which there exists s satisfying

        s^2 - d^2 = 4n mod m.
    """

    target = (
        4 * n
    ) % m

    square_set = set()

    for s in range(m):

        square_set.add(
            (s * s) % m
        )

    allowed = []

    for d in range(m):

        d2 = (
            d * d
        ) % m

        required = (
            d2 + target
        ) % m

        if required in square_set:

            allowed.append(d)

    return allowed


# ============================================================
# BUILD RAW WHEEL
# ============================================================

def build_raw_wheel(
    n,
    primes,
):
    """
    Build CRT-combined gap residues modulo product(primes).
    """

    states = {
        (0, 1)
    }

    for m in primes:

        allowed = allowed_gap_residues(
            n,
            m,
        )

        next_states = set()

        for residue, modulus in states:

            for a in allowed:

                merged = crt_pair(
                    residue,
                    modulus,
                    a,
                    m,
                )

                if merged is None:
                    continue

                next_states.add(
                    merged
                )

                if (
                    len(next_states)
                    > MAX_WHEEL_CLASSES
                ):
                    raise RuntimeError(
                        "too many wheel classes"
                    )

        states = next_states

    modulus = math.prod(
        primes
    )

    residues = sorted(
        residue
        for residue, current_modulus
        in states
        if current_modulus == modulus
    )

    return residues, modulus


# ============================================================
# BUILD FAST EVEN-GAP WHEEL
# ============================================================

def build_even_gap_wheel(
    n,
    primes,
):
    """
    Build admissible EVEN gap values.

    We search:

        d = 2*k.

    Instead of storing only parity-compatible residues modulo M,
    transform the problem to k modulo M.

    Since d=2k:

        d^2 = 4k^2.

    Thus

        s^2 - 4k^2 = 4n.

    We still construct d residues first, then convert them to
    the corresponding k residues.

    """

    residues, modulus = build_raw_wheel(
        n,
        primes,
    )

    even_residues = [
        r
        for r in residues
        if r % 2 == 0
    ]

    # We want k = d/2 modulo M.
    #
    # The d wheel has period M. Because M is odd, multiplication
    # by 2 is invertible modulo M.
    #
    # Therefore d = 2k mod M gives:
    #
    # k = d * inv(2) mod M.

    inv2 = inv_mod(
        2,
        modulus,
    )

    k_residues = sorted(
        (
            r
            * inv2
        ) % modulus
        for r in even_residues
    )

    # Remove duplicates.
    k_residues = sorted(
        set(k_residues)
    )

    return (
        k_residues,
        modulus,
    )


# ============================================================
# COUNT VALUES <= TARGET
# ============================================================

def count_before(
    target,
    residues,
    modulus,
):
    """
    Count even gaps d <= target represented by the wheel.
    """

    if target < 0:
        return 0

    max_k = (
        target // 2
    )

    count = 0

    for r in residues:

        if r > max_k:
            break

        count += 1

    if max_k >= modulus:

        full_periods = (
            max_k // modulus
        )

        remainder = (
            max_k
            % modulus
        )

        count = (
            full_periods
            * len(residues)
        )

        count += sum(
            1
            for r in residues
            if r <= remainder
        )

        return count

    return count


# ============================================================
# ORDINARY GAP SEARCH
# ============================================================

def ordinary_gap_search(
    n,
):

    # For odd p,q, d=q-p is even.
    d = 0

    tests = 0

    start = time.perf_counter()

    while (
        tests
        < MAX_GAP_TESTS
    ):

        value = (
            4 * n
            + d * d
        )

        s = math.isqrt(
            value
        )

        tests += 1

        if (
            s * s
            == value
        ):

            if (
                (s - d) % 2
                == 0
            ):

                p = (
                    s - d
                ) // 2

                q = (
                    s + d
                ) // 2

                if (
                    p > 1
                    and p * q == n
                ):

                    return {
                        "found": (
                            p,
                            q,
                        ),
                        "tests": tests,
                        "time":
                            time.perf_counter()
                            - start,
                        "aborted": None,
                    }

        d += 2

    return {
        "found": None,
        "tests": tests,
        "time":
            time.perf_counter()
            - start,
        "aborted":
            "MAX_GAP_TESTS",
    }


# ============================================================
# FAST WHEEL GAP SEARCH
# ============================================================

def wheel_gap_search(
    n,
    k_residues,
    modulus,
):
    """
    Search

        d = 2*k

    where k belongs to one of the admissible residue classes modulo M.

    We enumerate the sorted residues directly instead of repeatedly
    calling min() over all active arithmetic progressions.
    """

    # True d begins at 0, so k begins at 0.
    k = 0

    tests = 0

    start = time.perf_counter()

    period = modulus

    residue_count = len(
        k_residues
    )

    period_index = 0

    while (
        tests
        < MAX_GAP_TESTS
    ):

        base = (
            period_index
            * period
        )

        for residue in k_residues:

            k = (
                base
                + residue
            )

            d = (
                2 * k
            )

            value = (
                4 * n
                + d * d
            )

            s = math.isqrt(
                value
            )

            tests += 1

            if (
                s * s
                == value
            ):

                if (
                    (s - d) % 2
                    == 0
                ):

                    p = (
                        s - d
                    ) // 2

                    q = (
                        s + d
                    ) // 2

                    if (
                        p > 1
                        and p * q == n
                    ):

                        return {
                            "found": (
                                p,
                                q,
                            ),
                            "tests": tests,
                            "time":
                                time.perf_counter()
                                - start,
                            "aborted": None,
                        }

            # Do not continue beyond the actual true factor gap
            # once d becomes larger than a safe Fermat region.
            #
            # This is only a safety check. The exact solution cannot
            # require d larger than q-p, but we don't know q-p here.
            #
            # Therefore no artificial mathematical cutoff is used.

            if (
                tests
                >= MAX_GAP_TESTS
            ):
                break

        period_index += 1

    return {
        "found": None,
        "tests": tests,
        "time":
            time.perf_counter()
            - start,
        "aborted":
            "MAX_GAP_TESTS",
    }


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

    gap = q - p

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"sqrt(n) = "
        f"{math.isqrt(n)}"
    )

    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    print(
        f"true gap = {gap}"
    )

    # --------------------------------------------------------
    # Ordinary search.
    # --------------------------------------------------------

    print()
    print(
        "ORDINARY GAP SEARCH"
    )

    ordinary = ordinary_gap_search(
        n
    )

    print(
        f"    tests = "
        f"{ordinary['tests']}"
    )

    print(
        f"    time = "
        f"{ordinary['time']:.6f}s"
    )

    print(
        f"    found = "
        f"{ordinary['found']}"
    )

    if ordinary["aborted"]:

        print(
            f"    status = "
            f"{ordinary['aborted']}"
        )

    # --------------------------------------------------------
    # Wheel depths.
    # --------------------------------------------------------

    for depth in WHEEL_DEPTHS:

        primes = FILTER_PRIMES[
            :depth
        ]

        print()
        print(
            f"WHEEL DEPTH {depth}"
        )

        print(
            f"    primes = "
            f"{primes}"
        )

        try:

            t0 = time.perf_counter()

            residues, modulus = (
                build_even_gap_wheel(
                    n,
                    primes,
                )
            )

            build_time = (
                time.perf_counter()
                - t0
            )

        except RuntimeError as exc:

            print(
                f"    build failed: "
                f"{exc}"
            )

            continue

        density = (
            len(residues)
            / modulus
        )

        true_k = (
            gap // 2
        )

        true_residue = (
            true_k
            % modulus
        )

        survives = (
            true_residue
            in residues
        )

        # Count wheel candidates before true gap.
        before = count_before(
            gap,
            residues,
            modulus,
        )

        print(
            f"    modulus = "
            f"{modulus}"
        )

        print(
            f"    allowed k residues = "
            f"{len(residues)}"
        )

        print(
            f"    k density = "
            f"{density:.12f}"
        )

        print(
            f"    true k residue = "
            f"{true_residue}"
        )

        print(
            f"    true gap survives = "
            f"{survives}"
        )

        print(
            f"    filtered gaps <= true gap = "
            f"{before}"
        )

        print(
            f"    wheel build time = "
            f"{build_time:.6f}s"
        )

        if not survives:

            print(
                "    ERROR: true factor gap "
                "was filtered out"
            )

            continue

        result = wheel_gap_search(
            n,
            residues,
            modulus,
        )

        print(
            f"    wheel tests = "
            f"{result['tests']}"
        )

        print(
            f"    wheel search time = "
            f"{result['time']:.6f}s"
        )

        if result["found"]:

            fp, fq = result[
                "found"
            ]

            print(
                f"    FOUND = "
                f"({fp},{fq})"
            )

            print(
                f"    correct = "
                f"{fp * fq == n}"
            )

            if result["time"] > 0:

                print(
                    f"    speedup vs ordinary = "
                    f"{ordinary['time'] / result['time']:.3f}x"
                )

            # Once a factor is found, do not run larger wheels
            # for this instance unless explicitly desired.
            #
            # We continue because the point of the experiment is
            # to observe how the wheel density changes.
        else:

            if result["aborted"]:

                print(
                    f"    status = "
                    f"{result['aborted']}"
                )

    # --------------------------------------------------------
    # Exact identity.
    # --------------------------------------------------------

    x = (
        p + q
    ) // 2

    y = (
        q - p
    ) // 2

    print()
    print(
        "EXACT FERMAT IDENTITY"
    )

    print(
        f"    x = "
        f"{x}"
    )

    print(
        f"    y = "
        f"{y}"
    )

    print(
        f"    x²-y² = "
        f"{x*x-y*y}"
    )

    print(
        f"    n = "
        f"{n}"
    )

    print(
        f"    identity = "
        f"{x*x-y*y == n}"
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
        "START EXPERIMENT 188"
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
        f"FILTER_PRIMES = "
        f"{FILTER_PRIMES}"
    )

    print(
        f"WHEEL_DEPTHS = "
        f"{WHEEL_DEPTHS}"
    )

    print(
        f"MAX_GAP_TESTS = "
        f"{MAX_GAP_TESTS}"
    )

    print(
        f"MAX_WHEEL_CLASSES = "
        f"{MAX_WHEEL_CLASSES}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 188"
    )


if __name__ == "__main__":
    main()
