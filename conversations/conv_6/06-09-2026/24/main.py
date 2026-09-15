#!/usr/bin/env python3

"""
START EXPERIMENT 180

C-MATRIX TAUTOLOGY / INFORMATION TEST

Purpose
-------

Experiments 177B-179B demonstrated that:

    observed C-values
        ->
    residue constraints
        ->
    fast CRT factor recovery.

Experiment 180 asks the crucial question:

    Does C itself provide information about the factors when C is
    calculated from n?

For arbitrary nonzero residues

    a_i mod r_i
    b_j mod s_j

define

    z_ij = a_i*b_j

and

    C_ij =
        floor(
            (r_i*beta + s_j*alpha + z_ij)
            / (r_i*s_j)
        )

where

    beta  = (n-z_ij)*r_i^{-1} mod s_j
    alpha = (n-z_ij)*s_j^{-1} mod r_i.

Then the chosen (a_i,b_j) pair is automatically valid for that
C_ij by construction.

Therefore the C-matrix may be a representation of the residue
assignment rather than an independently constraining observation.

This experiment demonstrates that explicitly.

Tests:

    1. Hidden true factor residues.
    2. All-ones residues.
    3. Random residue vectors.
    4. Random residue vectors with CRT checks.

For every assignment we construct its C-matrix and then verify that
every chosen (a_i,b_j) belongs to the cell relation for that C-value.

If arbitrary assignments pass this test, then C-derived-only CSP
cannot distinguish the factor residue vector from arbitrary residue
vectors.

The exact factor equation pq=n is then the only missing information.

"""


import math
import random
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

RANDOM_TRIALS = 1000

RANDOM_SEED = 180

SHOW_RANDOM_TRIALS = 10


# ============================================================
# NUMBER THEORY
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
            f"{a} not invertible mod {m}"
        )

    return x % m


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
# C VALUE
# ============================================================

def c_value(
    n,
    r,
    s,
    z,
):

    inv_r_s = inv_mod(
        r % s,
        s,
    )

    inv_s_r = inv_mod(
        s % r,
        r,
    )

    beta = (
        (n - z)
        * inv_r_s
    ) % s

    alpha = (
        (n - z)
        * inv_s_r
    ) % r

    return (
        r * beta
        + s * alpha
        + z
    ) // (
        r * s
    )


# ============================================================
# PRECOMPUTED CELL MAP
# ============================================================

def build_cell_map(
    n,
    r,
    s,
):
    """
    Build:

        C -> set((a,b))

    for one cell.
    """

    inv_r_s = inv_mod(
        r % s,
        s,
    )

    inv_s_r = inv_mod(
        s % r,
        r,
    )

    result = {}

    for a in range(
        1,
        r,
    ):

        for b in range(
            1,
            s,
        ):

            z = a * b

            beta = (
                (n - z)
                * inv_r_s
            ) % s

            alpha = (
                (n - z)
                * inv_s_r
            ) % r

            c = (
                r * beta
                + s * alpha
                + z
            ) // (
                r * s
            )

            if c not in result:
                result[c] = set()

            result[c].add(
                (a, b)
            )

    return result


def build_all_cell_maps(n):

    cells = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            cells[(i, j)] = (
                build_cell_map(
                    n,
                    r,
                    s,
                )
            )

    return cells


# ============================================================
# GENERATE C MATRIX
# ============================================================

def generate_c_matrix(
    n,
    A,
    B,
):

    C = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            a = A[i]
            b = B[j]

            z = a * b

            c = c_value(
                n,
                r,
                s,
                z,
            )

            C[(i, j)] = c

    return C


# ============================================================
# VERIFY ASSIGNMENT AGAINST ITS C MATRIX
# ============================================================

def verify_assignment(
    A,
    B,
    C,
    cell_maps,
):

    failures = []

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            a = A[i]
            b = B[j]
            c = C[(i, j)]

            relation = cell_maps[
                (i, j)
            ]

            pairs = relation.get(
                c,
                set(),
            )

            if (a, b) not in pairs:

                failures.append({
                    "cell": (
                        r,
                        s,
                    ),
                    "a": a,
                    "b": b,
                    "c": c,
                })

    return failures


# ============================================================
# RANDOM RESIDUE VECTORS
# ============================================================

def random_residue_vector(
    radices,
    rng,
):

    return [
        rng.randrange(
            1,
            r,
        )
        for r in radices
    ]


# ============================================================
# CRT RESIDUE
# ============================================================

def crt_vector(
    residues,
    moduli,
):

    x = 0
    M = 1

    for a, m in zip(
        residues,
        moduli,
    ):

        g = math.gcd(
            M,
            m,
        )

        if (
            x - a
        ) % g != 0:

            return None, None

        M1 = M // g
        m1 = m // g

        rhs = (
            (a - x)
            // g
        )

        if m1 == 1:

            t = 0

        else:

            inv = inv_mod(
                M1 % m1,
                m1,
            )

            t = (
                rhs * inv
            ) % m1

        x = (
            x
            + M * t
        )

        M = (
            M
            * m1
        )

        x %= M

    return x, M


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

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"p = {p}"
    )

    print(
        f"q = {q}"
    )

    print()

    # --------------------------------------------------------
    # Build C maps.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    cell_maps = build_all_cell_maps(
        n
    )

    build_time = (
        time.perf_counter()
        - t0
    )

    print(
        f"cell-map construction = "
        f"{build_time:.6f}s"
    )

    # --------------------------------------------------------
    # True residue vectors.
    # --------------------------------------------------------

    true_A = [
        p % r
        for r in R1
    ]

    true_B = [
        q % s
        for s in R2
    ]

    print()
    print(
        f"TRUE A = {true_A}"
    )

    print(
        f"TRUE B = {true_B}"
    )

    true_C = generate_c_matrix(
        n,
        true_A,
        true_B,
    )

    true_failures = (
        verify_assignment(
            true_A,
            true_B,
            true_C,
            cell_maps,
        )
    )

    print(
        f"TRUE assignment failures = "
        f"{len(true_failures)}"
    )

    # --------------------------------------------------------
    # All-ones assignment.
    # --------------------------------------------------------

    ones_A = [
        1
        for _ in R1
    ]

    ones_B = [
        1
        for _ in R2
    ]

    ones_C = generate_c_matrix(
        n,
        ones_A,
        ones_B,
    )

    ones_failures = (
        verify_assignment(
            ones_A,
            ones_B,
            ones_C,
            cell_maps,
        )
    )

    print()
    print(
        "ALL-ONES ASSIGNMENT"
    )

    print(
        f"    A = {ones_A}"
    )

    print(
        f"    B = {ones_B}"
    )

    print(
        f"    failures = "
        f"{len(ones_failures)}"
    )

    print(
        f"    exactly same C rule = "
        f"{len(ones_failures) == 0}"
    )

    # --------------------------------------------------------
    # Random assignments.
    # --------------------------------------------------------

    rng = random.Random(
        RANDOM_SEED
    )

    successful = 0

    failures = 0

    examples = []

    t0 = time.perf_counter()

    for trial in range(
        RANDOM_TRIALS
    ):

        A = random_residue_vector(
            R1,
            rng,
        )

        B = random_residue_vector(
            R2,
            rng,
        )

        C = generate_c_matrix(
            n,
            A,
            B,
        )

        bad = verify_assignment(
            A,
            B,
            C,
            cell_maps,
        )

        if bad:

            failures += 1

            if len(examples) < 5:

                examples.append(
                    (
                        A,
                        B,
                        bad,
                    )
                )

        else:

            successful += 1

            if (
                len(
                    examples
                )
                < SHOW_RANDOM_TRIALS
            ):

                examples.append(
                    (
                        A,
                        B,
                        None,
                    )
                )

        if trial < SHOW_RANDOM_TRIALS:

            print()
            print(
                f"RANDOM TRIAL {trial + 1}"
            )

            print(
                f"    A = {A}"
            )

            print(
                f"    B = {B}"
            )

            print(
                f"    C(0,0) = "
                f"{C[(0,0)]}"
            )

            print(
                f"    valid = "
                f"{len(bad) == 0}"
            )

    random_time = (
        time.perf_counter()
        - t0
    )

    print()
    print(
        "RANDOM ASSIGNMENT RESULTS"
    )

    print(
        f"    trials = "
        f"{RANDOM_TRIALS}"
    )

    print(
        f"    successful = "
        f"{successful}"
    )

    print(
        f"    failures = "
        f"{failures}"
    )

    print(
        f"    success rate = "
        f"{successful / RANDOM_TRIALS:.6f}"
    )

    print(
        f"    time = "
        f"{random_time:.6f}s"
    )

    # --------------------------------------------------------
    # Show one deliberately false residue vector and its C
    # matrix.
    # --------------------------------------------------------

    false_A = [
        1,
        2,
        3,
        4,
        5,
    ]

    false_B = [
        1,
        2,
        3,
        4,
        5,
    ]

    false_C = generate_c_matrix(
        n,
        false_A,
        false_B,
    )

    false_failures = (
        verify_assignment(
            false_A,
            false_B,
            false_C,
            cell_maps,
        )
    )

    print()
    print(
        "DELIBERATELY FALSE RESIDUES"
    )

    print(
        f"    A = {false_A}"
    )

    print(
        f"    B = {false_B}"
    )

    print(
        f"    failures = "
        f"{len(false_failures)}"
    )

    print(
        "    These residues are not the hidden "
        "factor residues, but their generated "
        "C-matrix is nevertheless perfectly "
        "self-consistent."
    )

    # --------------------------------------------------------
    # CRT reconstruction of the hidden residues.
    # --------------------------------------------------------

    p_reconstructed, p_modulus = (
        crt_vector(
            true_A,
            R1,
        )
    )

    q_reconstructed, q_modulus = (
        crt_vector(
            true_B,
            R2,
        )
    )

    print()
    print(
        "TRUE RESIDUE CRT"
    )

    print(
        f"    p0 = {p_reconstructed}"
    )

    print(
        f"    Mp = {p_modulus}"
    )

    print(
        f"    q0 = {q_reconstructed}"
    )

    print(
        f"    Mq = {q_modulus}"
    )

    print(
        f"    actual p mod Mp = "
        f"{p % p_modulus}"
    )

    print(
        f"    actual q mod Mq = "
        f"{q % q_modulus}"
    )

    # --------------------------------------------------------
    # Final interpretation.
    # --------------------------------------------------------

    print()
    print(
        "INTERPRETATION"
    )

    if (
        successful == RANDOM_TRIALS
        and len(ones_failures) == 0
        and len(false_failures) == 0
    ):

        print(
            "    C-derived cell relations accept "
            "arbitrary residue vectors."
        )

        print(
            "    Therefore C(n,r,s,a*b) is not "
            "an independent observable constraint."
        )

        print(
            "    The actual factorization information "
            "must enter somewhere else."
        )

    else:

        print(
            "    Not all arbitrary assignments passed."
        )

        print(
            "    Further C-structure investigation "
            "may still be useful."
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
        "START EXPERIMENT 180"
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
        f"RANDOM_TRIALS = {RANDOM_TRIALS}"
    )

    print(
        f"RANDOM_SEED = {RANDOM_SEED}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 180"
    )


if __name__ == "__main__":
    main()
