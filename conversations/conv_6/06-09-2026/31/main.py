#!/usr/bin/env python3

"""
START EXPERIMENT 187

FACTOR-GAP SEARCH WITH MODULAR SQUARE FILTER

Fermat factorization can be written in terms of the factor gap:

    d = q - p

For

    n = p*q

we have

    p + q = sqrt(4n + d^2).

Therefore:

    4n + d^2 = s^2

for some integer s, and

    p = (s-d)/2
    q = (s+d)/2.

Ordinary Fermat searches increasing x=(p+q)/2.

Experiment 187 instead searches the gap d directly.

The exact condition is:

    s^2 = 4n + d^2.

Equivalently:

    s^2 - d^2 = 4n.

For every small modulus m we can precompute which residues d mod m
can possibly satisfy

    s^2-d^2 = 4n mod m.

Then only those d-values are tested.

This differs from Experiment 186:

    Exp 186:
        search x=(p+q)/2

    Exp 187:
        search d=q-p

The experiment measures:

    * actual factor gap
    * number of d values below the true gap
    * modularly admissible d density
    * ordinary gap search
    * filtered gap search

The hidden p,q are used only for verification.

No C-values are used.

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
    3, 5, 7, 11, 13, 17, 19, 23
]

FILTER_DEPTHS = [
    1, 2, 3, 4, 5, 6
]

MAX_GAP_STEPS = 20_000_000

MAX_RESIDUES = 2_000_000


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

    x = a1 + m1 * t

    modulus = (
        m1
        * m2r
    )

    return (
        x % modulus,
        modulus,
    )


# ============================================================
# MODULAR GAP RESIDUES
# ============================================================

def allowed_gap_residues(
    n,
    m,
):
    """
    Return d residues modulo m for which there exists an s residue
    satisfying

        s^2 - d^2 = 4n mod m.
    """

    target = (
        4 * n
    ) % m

    squares = set()

    for s in range(m):
        squares.add(
            (s * s) % m
        )

    allowed = set()

    for d in range(m):

        d2 = (
            d * d
        ) % m

        required = (
            target
            + d2
        ) % m

        if required in squares:
            allowed.add(d)

    return allowed


# ============================================================
# BUILD GAP WHEEL
# ============================================================

def build_gap_wheel(
    n,
    primes,
):
    """
    Combine allowed d residues from all selected moduli.
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
                    > MAX_RESIDUES
                ):
                    raise RuntimeError(
                        "too many gap residue classes"
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
# ORDINARY GAP SEARCH
# ============================================================

def ordinary_gap_search(
    n,
    start_d,
):

    steps = 0

    d = start_d

    start = time.perf_counter()

    while steps < MAX_GAP_STEPS:

        value = (
            4 * n
            + d * d
        )

        s = math.isqrt(
            value
        )

        if (
            s * s
            == value
        ):

            # Need same parity.
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
                        "steps":
                            steps + 1,
                        "time":
                            time.perf_counter()
                            - start,
                        "aborted":
                            None,
                    }

        d += 2
        steps += 1

    return {
        "found": None,
        "steps": steps,
        "time":
            time.perf_counter()
            - start,
        "aborted":
            "MAX_GAP_STEPS",
    }


# ============================================================
# FILTERED GAP SEARCH
# ============================================================

def filtered_gap_search(
    n,
    residues,
    modulus,
    start_d,
):

    # Gap must have same parity as q-p.
    #
    # For odd p,q, d is even.
    #
    # So only retain even d.

    usable_residues = [
        r
        for r in residues
        if (
            r % 2
            == start_d % 2
        )
    ]

    if not usable_residues:
        return {
            "found": None,
            "steps": 0,
            "time": 0.0,
            "aborted":
                "NO_RESIDUES",
        }

    current = []

    for residue in usable_residues:

        if residue >= start_d:

            first = residue

        else:

            k = (
                start_d
                - residue
                + modulus
                - 1
            ) // modulus

            first = (
                residue
                + k * modulus
            )

        if (
            first
            % 2
            != start_d % 2
        ):
            first += 1

            if first >= start_d:
                pass

        current.append(
            first
        )

    steps = 0

    start = time.perf_counter()

    while steps < MAX_GAP_STEPS:

        idx = min(
            range(
                len(current)
            ),
            key=current.__getitem__,
        )

        d = current[idx]

        value = (
            4 * n
            + d * d
        )

        s = math.isqrt(
            value
        )

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
                        "steps":
                            steps + 1,
                        "time":
                            time.perf_counter()
                            - start,
                        "aborted":
                            None,
                    }

        current[idx] += (
            2 * modulus
        )

        steps += 1

    return {
        "found": None,
        "steps": steps,
        "time":
            time.perf_counter()
            - start,
        "aborted":
            "MAX_GAP_STEPS",
    }


# ============================================================
# COUNT FILTERED GAPS
# ============================================================

def count_filtered_gaps(
    start_d,
    end_d,
    residues,
    modulus,
):

    if end_d < start_d:
        return 0

    count = 0

    for residue in residues:

        if (
            residue
            % 2
            != start_d % 2
        ):
            continue

        if residue >= start_d:

            first = residue

        else:

            k = (
                start_d
                - residue
                + modulus
                - 1
            ) // modulus

            first = (
                residue
                + k * modulus
            )

        if (
            first
            % 2
            != start_d % 2
        ):

            first += 1

        if first > end_d:
            continue

        count += (
            (
                end_d
                - first
            )
            // (
                2 * modulus
            )
            + 1
        )

    return count


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

    sqrt_n = math.isqrt(
        n
    )

    gap = q - p

    # For odd p,q the gap is even.
    start_d = 0

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"sqrt(n) = {sqrt_n}"
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

    print()

    # --------------------------------------------------------
    # Ordinary search.
    # --------------------------------------------------------

    print(
        "ORDINARY GAP SEARCH"
    )

    ordinary = ordinary_gap_search(
        n,
        start_d,
    )

    print(
        f"    steps = "
        f"{ordinary['steps']}"
    )

    print(
        f"    time = "
        f"{ordinary['time']:.6f}s"
    )

    if ordinary["found"]:

        print(
            f"    FOUND = "
            f"{ordinary['found']}"
        )

    if ordinary["aborted"]:

        print(
            f"    status = "
            f"{ordinary['aborted']}"
        )

    # --------------------------------------------------------
    # Modular wheels.
    # --------------------------------------------------------

    for depth in FILTER_DEPTHS:

        primes = FILTER_PRIMES[
            :depth
        ]

        print()
        print(
            f"GAP WHEEL DEPTH {depth}"
        )

        print(
            f"    primes = {primes}"
        )

        try:

            t0 = time.perf_counter()

            residues, modulus = (
                build_gap_wheel(
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
                f"    wheel build failed: "
                f"{exc}"
            )

            continue

        parity_residues = [
            r
            for r in residues
            if (
                r % 2
                == start_d % 2
            )
        ]

        # Effective density among even gaps.
        if modulus % 2 == 0:

            effective_density = (
                len(parity_residues)
                / modulus
            )

        else:

            effective_density = (
                len(parity_residues)
                / modulus
            )

        true_residue = (
            gap % modulus
        )

        survives = (
            true_residue
            in residues
        )

        count_before = (
            count_filtered_gaps(
                start_d,
                gap,
                residues,
                modulus,
            )
        )

        print(
            f"    modulus = "
            f"{modulus}"
        )

        print(
            f"    total allowed residues = "
            f"{len(residues)}"
        )

        print(
            f"    parity-compatible residues = "
            f"{len(parity_residues)}"
        )

        print(
            f"    effective density = "
            f"{effective_density:.12f}"
        )

        print(
            f"    true gap residue = "
            f"{true_residue}"
        )

        print(
            f"    true gap survives = "
            f"{survives}"
        )

        print(
            f"    filtered gaps <= true gap = "
            f"{count_before}"
        )

        print(
            f"    wheel build time = "
            f"{build_time:.6f}s"
        )

        if (
            len(parity_residues)
            > 500_000
        ):

            print(
                "    filtered search = "
                "SKIP (too many classes)"
            )

            continue

        result = filtered_gap_search(
            n,
            residues,
            modulus,
            start_d,
        )

        print(
            f"    filtered iterations = "
            f"{result['steps']}"
        )

        print(
            f"    search time = "
            f"{result['time']:.6f}s"
        )

        if result["found"]:

            print(
                f"    FOUND = "
                f"{result['found']}"
            )

        if result["aborted"]:

            print(
                f"    status = "
                f"{result['aborted']}"
            )

    # --------------------------------------------------------
    # Exact identity verification.
    # --------------------------------------------------------

    s_true = (
        p + q
    )

    print()
    print(
        "TRUE GAP IDENTITY"
    )

    print(
        f"    d = q-p = "
        f"{gap}"
    )

    print(
        f"    s = p+q = "
        f"{s_true}"
    )

    print(
        f"    s²-d² = "
        f"{s_true*s_true - gap*gap}"
    )

    print(
        f"    4n = "
        f"{4*n}"
    )

    print(
        f"    exact = "
        f"{s_true*s_true - gap*gap == 4*n}"
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
        "START EXPERIMENT 187"
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
        f"FILTER_DEPTHS = "
        f"{FILTER_DEPTHS}"
    )

    print(
        f"MAX_GAP_STEPS = "
        f"{MAX_GAP_STEPS}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 187"
    )


if __name__ == "__main__":
    main()
