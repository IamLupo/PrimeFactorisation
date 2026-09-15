#!/usr/bin/env python3

"""
START EXPERIMENT 186B

MODULAR FERMAT SEARCH

We study Fermat factorization:

    n = p*q

with

    x = (p+q)/2
    y = (q-p)/2

so

    x^2 - y^2 = n.

Ordinary Fermat searches:

    x = ceil(sqrt(n)), ceil(sqrt(n))+1, ...

until

    x^2-n

is a perfect square.

This experiment asks whether modular difference-of-squares conditions
can filter the Fermat x values.

For a modulus m:

    x^2-y^2 = n (mod m)

must hold for the true x,y.

For each prime m we therefore determine the allowed residues

    x (mod m)

for which at least one y (mod m) exists satisfying the equation.

Those x residue classes are combined with CRT into a wheel.

We then compare:

    A) ordinary Fermat
    B) modularly filtered Fermat

No factor residues are supplied.

No C-values are used.

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

WHEEL_PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23
]

WHEEL_DEPTHS = [
    1, 2, 3, 4
]

MAX_FERMAT_STEPS = 20_000_000

MAX_WHEEL_CLASSES = 2_000_000


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
        f"Could not generate {bits}-bit semiprime"
    )


# ============================================================
# EXTENDED GCD
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


# ============================================================
# GENERALIZED CRT
# ============================================================

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
# MODULAR DIFFERENCE OF SQUARES
# ============================================================

def allowed_x_residues(
    n,
    m,
):
    """
    Determine all x mod m for which some y mod m exists such that

        x^2 - y^2 = n mod m.
    """

    n_mod = n % m

    allowed = set()

    # Precompute y^2.
    square_values = set()

    for y in range(m):

        square_values.add(
            (y * y) % m
        )

    for x in range(m):

        x2 = (
            x * x
        ) % m

        # Need:
        #
        # x^2-y^2 = n
        #
        # therefore:
        #
        # y^2 = x^2-n

        required = (
            x2 - n_mod
        ) % m

        if required in square_values:

            allowed.add(x)

    return allowed


# ============================================================
# BUILD WHEEL
# ============================================================

def build_x_wheel(
    n,
    primes,
):
    """
    Build CRT classes of x satisfying the modular
    difference-of-squares condition for every prime in primes.
    """

    states = {
        (0, 1)
    }

    for m in primes:

        allowed = allowed_x_residues(
            n,
            m,
        )

        new_states = set()

        for x, modulus in states:

            for residue in allowed:

                merged = crt_pair(
                    x,
                    modulus,
                    residue,
                    m,
                )

                if merged is None:
                    continue

                new_states.add(
                    merged
                )

                if (
                    len(new_states)
                    > MAX_WHEEL_CLASSES
                ):

                    raise RuntimeError(
                        "Too many wheel classes"
                    )

        states = new_states

    modulus = math.prod(
        primes
    )

    residues = sorted(
        x
        for x, current_modulus in states
        if current_modulus == modulus
    )

    return residues, modulus


# ============================================================
# ORDINARY FERMAT
# ============================================================

def fermat_search(
    n,
):

    sqrt_n = math.isqrt(
        n
    )

    x = sqrt_n

    if (
        x * x
        < n
    ):
        x += 1

    steps = 0

    start = time.perf_counter()

    while steps < MAX_FERMAT_STEPS:

        d = (
            x * x
            - n
        )

        y = math.isqrt(
            d
        )

        if (
            y * y
            == d
        ):

            p = x - y
            q = x + y

            if (
                p > 1
                and p * q == n
            ):

                return {
                    "found": (
                        p,
                        q,
                    ),
                    "steps": steps + 1,
                    "time":
                        time.perf_counter()
                        - start,
                    "aborted": None,
                }

        x += 1
        steps += 1

    return {
        "found": None,
        "steps": steps,
        "time":
            time.perf_counter()
            - start,
        "aborted":
            "MAX_FERMAT_STEPS",
    }


# ============================================================
# FILTERED FERMAT
# ============================================================

def filtered_fermat(
    n,
    residues,
    modulus,
):

    sqrt_n = math.isqrt(
        n
    )

    x0 = sqrt_n

    if (
        x0 * x0
        < n
    ):
        x0 += 1

    # --------------------------------------------------------
    # Build one current value for each residue class.
    # --------------------------------------------------------

    current = []

    for residue in residues:

        if residue >= x0:

            first = residue

        else:

            k = (
                x0
                - residue
                + modulus
                - 1
            ) // modulus

            first = (
                residue
                + k * modulus
            )

        current.append(
            first
        )

    steps = 0

    start = time.perf_counter()

    while steps < MAX_FERMAT_STEPS:

        # Find smallest currently available x.
        idx = min(
            range(len(current)),
            key=current.__getitem__,
        )

        x = current[idx]

        if x < x0:

            k = (
                x0
                - residues[idx]
                + modulus
                - 1
            ) // modulus

            x = (
                residues[idx]
                + k * modulus
            )

            current[idx] = x

        d = (
            x * x
            - n
        )

        y = math.isqrt(
            d
        )

        if (
            y * y
            == d
        ):

            p = x - y
            q = x + y

            if (
                p > 1
                and p * q == n
            ):

                return {
                    "found": (
                        p,
                        q,
                    ),
                    "steps": steps + 1,
                    "time":
                        time.perf_counter()
                        - start,
                    "aborted": None,
                }

        current[idx] += modulus

        steps += 1

    return {
        "found": None,
        "steps": steps,
        "time":
            time.perf_counter()
            - start,
        "aborted":
            "MAX_FERMAT_STEPS",
    }


# ============================================================
# COUNT FILTERED X VALUES
# ============================================================

def count_filtered_values(
    start_x,
    end_x,
    residues,
    modulus,
):

    if end_x < start_x:
        return 0

    count = 0

    for residue in residues:

        if residue >= start_x:

            first = residue

        else:

            k = (
                start_x
                - residue
                + modulus
                - 1
            ) // modulus

            first = (
                residue
                + k * modulus
            )

        if first > end_x:
            continue

        count += (
            (end_x - first)
            // modulus
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

    x_true = (
        p + q
    ) // 2

    y_true = (
        q - p
    ) // 2

    start_x = sqrt_n

    if (
        start_x * start_x
        < n
    ):
        start_x += 1

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
        f"true x = {x_true}"
    )

    print(
        f"true y = {y_true}"
    )

    print(
        f"Fermat distance = "
        f"{x_true - start_x}"
    )

    print(
        f"factor gap = "
        f"{q - p}"
    )

    # --------------------------------------------------------
    # Modular wheel tests.
    # --------------------------------------------------------

    for depth in WHEEL_DEPTHS:

        primes = WHEEL_PRIMES[
            :depth
        ]

        print()
        print(
            f"WHEEL DEPTH {depth}"
        )

        print(
            f"    primes = {primes}"
        )

        try:

            t0 = time.perf_counter()

            residues, modulus = (
                build_x_wheel(
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

        density = (
            len(residues)
            / modulus
        )

        true_residue = (
            x_true % modulus
        )

        survives = (
            true_residue
            in residues
        )

        count_before_true = (
            count_filtered_values(
                start_x,
                x_true,
                residues,
                modulus,
            )
        )

        print(
            f"    modulus = "
            f"{modulus}"
        )

        print(
            f"    allowed classes = "
            f"{len(residues)}"
        )

        print(
            f"    density = "
            f"{density:.12f}"
        )

        print(
            f"    true x residue = "
            f"{true_residue}"
        )

        print(
            f"    true x survives = "
            f"{survives}"
        )

        print(
            f"    filtered x values "
            f"to true x = "
            f"{count_before_true}"
        )

        print(
            f"    wheel build time = "
            f"{build_time:.6f}s"
        )

        # ----------------------------------------------------
        # Don't run a wheel with too many classes.
        # ----------------------------------------------------

        if (
            len(residues)
            > 500_000
        ):

            print(
                "    filtered Fermat = "
                "SKIP (too many classes)"
            )

            continue

        result = filtered_fermat(
            n,
            residues,
            modulus,
        )

        print(
            f"    Fermat iterations = "
            f"{result['steps']}"
        )

        print(
            f"    search time = "
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

        if result["aborted"]:

            print(
                f"    status = "
                f"{result['aborted']}"
            )

    # --------------------------------------------------------
    # Ordinary Fermat.
    # --------------------------------------------------------

    print()
    print(
        "ORDINARY FERMAT"
    )

    result = fermat_search(
        n
    )

    print(
        f"    iterations = "
        f"{result['steps']}"
    )

    print(
        f"    time = "
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

    if result["aborted"]:

        print(
            f"    status = "
            f"{result['aborted']}"
        )

    # --------------------------------------------------------
    # Direct modular verification.
    # --------------------------------------------------------

    print()
    print(
        "TRUE MODULAR EQUATION"
    )

    for m in WHEEL_PRIMES:

        lhs = (
            x_true * x_true
            - y_true * y_true
        ) % m

        rhs = n % m

        print(
            f"    m={m} "
            f"x²-y² mod m={lhs} "
            f"n mod m={rhs} "
            f"ok={lhs == rhs}"
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
        "START EXPERIMENT 186B"
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
        f"WHEEL_PRIMES = "
        f"{WHEEL_PRIMES}"
    )

    print(
        f"WHEEL_DEPTHS = "
        f"{WHEEL_DEPTHS}"
    )

    print(
        f"MAX_FERMAT_STEPS = "
        f"{MAX_FERMAT_STEPS}"
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
        "FINISHED EXPERIMENT 186B"
    )


if __name__ == "__main__":
    main()