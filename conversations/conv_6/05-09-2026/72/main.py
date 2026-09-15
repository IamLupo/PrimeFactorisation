# ================================================================
# START EXPERIMENT 155
# Carry-pattern CSP / row-column constraint propagation
#
# Goal:
#
#   Test whether the exact C-matrix contains useful information
#   once the shared residues
#
#       a_i = p mod r_i
#       b_j = q mod s_j
#
#   are coupled across rows and columns.
#
# Core decomposition:
#
#   E_ij = U_ij + V_ij + C_ij
#
# where
#
#   U_ij = floor(k_i*b_j / s_j)
#   V_ij = floor(ell_j*a_i / r_i)
#
# and
#
#   C_ij =
#     floor(
#       (r_i*beta_ij + s_j*alpha_ij + a_i*b_j)
#       / (r_i*s_j)
#     )
#
# with
#
#   beta_ij  = (k_i*b_j) mod s_j
#   alpha_ij = (ell_j*a_i) mod r_i
#
# Thus
#
#   C_ij in {0,1,2}.
#
# For the experiment we know the true C pattern from the generated
# factorization. The question is how strongly that oracle pattern
# constrains the hidden residue vectors a_i and b_j when combined
# with n mod (r_i*s_j).
#
# We DO NOT factor n here.
#
# ================================================================

import math
import random
import time


# ================================================================
# CONFIGURATION
# ================================================================

R1 = [
    17,
    43,
    59,
    71,
    83,
]

R2 = [
    19,
    37,
    61,
    73,
    89,
]

BIT_SIZES = [
    30,
    36,
    42,
    48,
    54,
]

SAMPLES_PER_SIZE = 3


# ================================================================
# PRIMALITY / SEMIPRIME GENERATION
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37,
    ]

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in [
        2, 3, 5, 7, 11, 13, 17,
    ]:

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        witness_failed = True

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                witness_failed = False
                break

        if witness_failed:
            return False

    return True


def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= (1 << (bits - 1))

        if is_probable_prime(x):
            return x


def random_semiprime(bits):

    p_bits = bits // 2
    q_bits = bits - p_bits

    while True:

        p = random_prime(p_bits)
        q = random_prime(q_bits)

        if p != q:

            return p, q


# ================================================================
# FACTOR RESIDUE DATA
# ================================================================

def factor_data(p, q):

    k = []
    a = []

    for r in R1:

        k.append(p // r)
        a.append(p % r)

    ell = []
    b = []

    for s in R2:

        ell.append(q // s)
        b.append(q % s)

    return k, a, ell, b


# ================================================================
# EXACT C MATRIX
# ================================================================

def build_C(
    p,
    q,
):

    k, a, ell, b = factor_data(
        p,
        q,
    )

    C = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            beta = (
                (k[i] * b[j])
                % s
            )

            alpha = (
                (ell[j] * a[i])
                % r
            )

            numerator = (
                r * beta
                + s * alpha
                + a[i] * b[j]
            )

            denominator = (
                r * s
            )

            cij = (
                numerator
                // denominator
            )

            row.append(cij)

        C.append(row)

    return C, a, b


# ================================================================
# LINEAR CONGRUENCE SOLVER
#
# Solve:
#
#     A*x == B (mod M)
#
# and return all x in [0,M-1].
#
# The radices here are prime, but this routine also handles the
# non-coprime case correctly.
# ================================================================

def solve_linear_congruence(
    A,
    B,
    M,
):

    A %= M
    B %= M

    g = math.gcd(
        A,
        M,
    )

    if B % g != 0:
        return []

    A1 = A // g
    B1 = B // g
    M1 = M // g

    if M1 == 1:
        return list(range(M))

    inv = pow(
        A1,
        -1,
        M1,
    )

    x0 = (
        B1 * inv
    ) % M1

    solutions = []

    for t in range(g):

        x = (
            x0
            + t * M1
        ) % M

        solutions.append(x)

    return solutions


# ================================================================
# CELL RELATION
#
# For fixed (r,s), a,b must satisfy:
#
#   n = pq mod rs
#
# which gives
#
#   n - ab == r*x*b + s*y*a  (mod rs)
#
# where
#
#   x = k mod s
#   y = ell mod r.
#
# Equivalently:
#
#   r*b*x == n-ab (mod s)
#
#   s*a*y == n-ab (mod r)
#
# We solve those two congruences.
#
# Then we test whether at least one compatible (x,y) gives the
# observed C value.
#
# ================================================================

def build_cell_relation(
    n,
    r,
    s,
    observed_c,
):

    relation = set()

    n_mod_rs = (
        n % (r * s)
    )

    for a in range(r):

        for b in range(s):

            target = (
                n_mod_rs
                - a * b
            ) % (r * s)

            target_mod_s = (
                target % s
            )

            target_mod_r = (
                target % r
            )

            # ----------------------------------------------------
            # Solve for x = k mod s
            #
            # r*b*x == target (mod s)
            # ----------------------------------------------------

            x_solutions = solve_linear_congruence(
                (r * b) % s,
                target_mod_s,
                s,
            )

            if not x_solutions:
                continue

            # ----------------------------------------------------
            # Solve for y = ell mod r
            #
            # s*a*y == target (mod r)
            # ----------------------------------------------------

            y_solutions = solve_linear_congruence(
                (s * a) % r,
                target_mod_r,
                r,
            )

            if not y_solutions:
                continue

            # ----------------------------------------------------
            # Check whether any x,y state produces the observed C.
            # ----------------------------------------------------

            compatible = False

            for x in x_solutions:

                beta = (
                    x * b
                ) % s

                for y in y_solutions:

                    alpha = (
                        y * a
                    ) % r

                    numerator = (
                        r * beta
                        + s * alpha
                        + a * b
                    )

                    cij = (
                        numerator
                        // (r * s)
                    )

                    if cij == observed_c:

                        compatible = True
                        break

                if compatible:
                    break

            if compatible:

                relation.add(
                    (a, b)
                )

    return relation


# ================================================================
# ARC CONSISTENCY
#
# Variables:
#
#   A_i = possible a_i values
#   B_j = possible b_j values
#
# Each cell (i,j) contains a binary relation:
#
#   R_ij subset [0,r_i) x [0,s_j)
#
# We repeatedly remove values which have no supporting partner.
#
# This is classic arc-consistency, but implemented directly.
# ================================================================

def arc_consistency(
    domains_a,
    domains_b,
    relations,
):

    changed = True
    iterations = 0

    while changed:

        changed = False
        iterations += 1

        # --------------------------------------------------------
        # A -> B
        # --------------------------------------------------------

        for i in range(len(R1)):

            for j in range(len(R2)):

                relation = relations[i][j]

                old_a = domains_a[i]

                supported_a = set()

                b_domain = domains_b[j]

                for a in old_a:

                    for aa, bb in relation:

                        if aa == a and bb in b_domain:

                            supported_a.add(a)
                            break

                if supported_a != old_a:

                    domains_a[i] = supported_a
                    changed = True

        # --------------------------------------------------------
        # B -> A
        # --------------------------------------------------------

        for i in range(len(R1)):

            for j in range(len(R2)):

                relation = relations[i][j]

                old_b = domains_b[j]

                supported_b = set()

                a_domain = domains_a[i]

                for b in old_b:

                    for aa, bb in relation:

                        if (
                            bb == b
                            and aa in a_domain
                        ):

                            supported_b.add(b)
                            break

                if supported_b != old_b:

                    domains_b[j] = supported_b
                    changed = True

    return iterations


# ================================================================
# DOMAIN HELPERS
# ================================================================

def domain_product(domains):

    result = 1

    for d in domains:

        result *= max(
            1,
            len(d),
        )

    return result


def domain_summary(
    name,
    domains,
    radices,
):

    print(name)

    for i, domain in enumerate(domains):

        values = sorted(domain)

        print(
            f"  {i}: "
            f"{len(values):4d}/{radices[i]:4d}",
            values
            if len(values) <= 20
            else (
                values[:10]
                + ["..."]
                + values[-5:]
            )
        )

    print()


# ================================================================
# MAIN EXPERIMENT
# ================================================================

print("=" * 72)
print("START EXPERIMENT 155")
print("Carry-pattern CSP / row-column constraint propagation")
print("=" * 72)
print()

print("R1 =", R1)
print("R2 =", R2)
print(
    "SAMPLES PER BIT SIZE =",
    SAMPLES_PER_SIZE
)
print()

aggregate = {
    "cells": len(R1) * len(R2) * len(BIT_SIZES) * SAMPLES_PER_SIZE,
    "true_pair_survived": 0,
    "cells_with_reduction": 0,
    "cells_empty": 0,
}

for bits in BIT_SIZES:

    print("=" * 72)
    print(
        f"BIT SIZE = {bits}"
    )
    print("=" * 72)
    print()

    for sample in range(
        1,
        SAMPLES_PER_SIZE + 1
    ):

        t0 = time.perf_counter()

        p, q = random_semiprime(
            bits
        )

        n = p * q

        C, true_a, true_b = build_C(
            p,
            q,
        )

        print("-" * 72)
        print(
            f"SAMPLE {sample}/{SAMPLES_PER_SIZE}"
        )

        print(
            "n bits =",
            n.bit_length()
        )

        print(
            "n      =",
            n
        )

        print(
            "p      =",
            p
        )

        print(
            "q      =",
            q
        )

        print()

        print(
            "TRUE a =",
            true_a
        )

        print(
            "TRUE b =",
            true_b
        )

        print()

        # --------------------------------------------------------
        # Print oracle C.
        # --------------------------------------------------------

        print(
            "ORACLE C MATRIX"
        )

        for row in C:
            print(row)

        print()

        # --------------------------------------------------------
        # Build every binary cell relation.
        # --------------------------------------------------------

        relations = [
            [None for _ in R2]
            for _ in R1
        ]

        total_possible_pairs = 0
        total_mod_pairs = 0
        total_c_pairs = 0

        min_cell = None
        max_cell = None

        print(
            "CELL RELATIONS"
        )

        for i, r in enumerate(R1):

            for j, s in enumerate(R2):

                relation = build_cell_relation(
                    n,
                    r,
                    s,
                    C[i][j],
                )

                relations[i][j] = relation

                possible = (
                    r * s
                )

                total_possible_pairs += possible
                total_c_pairs += len(
                    relation
                )

                # Check the true pair.
                true_pair = (
                    true_a[i],
                    true_b[j],
                )

                survived = (
                    true_pair
                    in relation
                )

                if survived:

                    aggregate[
                        "true_pair_survived"
                    ] += 1

                else:

                    print(
                        "WARNING: TRUE PAIR "
                        "NOT IN RELATION:",
                        i,
                        j,
                        true_pair,
                    )

                if len(relation) < possible:

                    aggregate[
                        "cells_with_reduction"
                    ] += 1

                if len(relation) == 0:

                    aggregate[
                        "cells_empty"
                    ] += 1

                if (
                    min_cell is None
                    or len(relation)
                    < min_cell
                ):

                    min_cell = len(
                        relation
                    )

                if (
                    max_cell is None
                    or len(relation)
                    > max_cell
                ):

                    max_cell = len(
                        relation
                    )

                print(
                    f"  ({i},{j}) "
                    f"r={r:2d} s={s:2d} "
                    f"C={C[i][j]} "
                    f"states="
                    f"{len(relation):5d}/"
                    f"{possible:5d} "
                    f"true="
                    f"{'YES' if survived else 'NO'}"
                )

        print()

        print(
            "CELL STATE SUMMARY"
        )

        print(
            "  minimum relation size =",
            min_cell
        )

        print(
            "  maximum relation size =",
            max_cell
        )

        print(
            "  total raw pairs       =",
            total_possible_pairs
        )

        print(
            "  total C-compatible    =",
            total_c_pairs
        )

        print(
            "  raw -> C reduction    =",
            (
                total_c_pairs
                / total_possible_pairs
            )
        )

        print()

        # --------------------------------------------------------
        # Initial domains.
        # --------------------------------------------------------

        domains_a = [
            set(
                range(r)
            )
            for r in R1
        ]

        domains_b = [
            set(
                range(s)
            )
            for s in R2
        ]

        initial_a_product = domain_product(
            domains_a
        )

        initial_b_product = domain_product(
            domains_b
        )

        print(
            "INITIAL DOMAIN SIZE"
        )

        print(
            "  A product =",
            initial_a_product
        )

        print(
            "  B product =",
            initial_b_product
        )

        print()

        # --------------------------------------------------------
        # Direct row/column support before iterative propagation.
        #
        # This measures how much one entire row or column says
        # without feeding reductions back into neighbouring cells.
        # --------------------------------------------------------

        row_domains = []

        for i, r in enumerate(R1):

            allowed = set(
                range(r)
            )

            for j in range(len(R2)):

                relation = relations[i][j]

                supported = {
                    aa
                    for aa, bb in relation
                    if bb >= 0
                }

                allowed &= supported

            row_domains.append(
                allowed
            )

        col_domains = []

        for j, s in enumerate(R2):

            allowed = set(
                range(s)
            )

            for i in range(len(R1)):

                relation = relations[i][j]

                supported = {
                    bb
                    for aa, bb in relation
                    if aa >= 0
                }

                allowed &= supported

            col_domains.append(
                allowed
            )

        print(
            "ROW-SHARED DOMAINS "
            "(one pass)"
        )

        for i, domain in enumerate(
            row_domains
        ):

            print(
                f"  a[{i}] = "
                f"{len(domain):3d}/{R1[i]} "
                f"true="
                f"{true_a[i] in domain}"
            )

        print()

        print(
            "COLUMN-SHARED DOMAINS "
            "(one pass)"
        )

        for j, domain in enumerate(
            col_domains
        ):

            print(
                f"  b[{j}] = "
                f"{len(domain):3d}/{R2[j]} "
                f"true="
                f"{true_b[j] in domain}"
            )

        print()

        # --------------------------------------------------------
        # Install those one-pass domains as the starting CSP.
        # --------------------------------------------------------

        domains_a = [
            set(x)
            for x in row_domains
        ]

        domains_b = [
            set(x)
            for x in col_domains
        ]

        before_ac_a = [
            len(x)
            for x in domains_a
        ]

        before_ac_b = [
            len(x)
            for x in domains_b
        ]

        # --------------------------------------------------------
        # Arc consistency.
        # --------------------------------------------------------

        iterations = arc_consistency(
            domains_a,
            domains_b,
            relations,
        )

        print(
            "ARC CONSISTENCY"
        )

        print(
            "  iterations =",
            iterations
        )

        print()

        domain_summary(
            "FINAL A DOMAINS",
            domains_a,
            R1,
        )

        domain_summary(
            "FINAL B DOMAINS",
            domains_b,
            R2,
        )

        # --------------------------------------------------------
        # True residue verification.
        # --------------------------------------------------------

        true_a_alive = all(
            true_a[i]
            in domains_a[i]
            for i in range(len(R1))
        )

        true_b_alive = all(
            true_b[j]
            in domains_b[j]
            for j in range(len(R2))
        )

        print(
            "TRUE STATE SURVIVAL"
        )

        print(
            "  all true a_i survive =",
            true_a_alive
        )

        print(
            "  all true b_j survive =",
            true_b_alive
        )

        print()

        # --------------------------------------------------------
        # Final domain sizes.
        # --------------------------------------------------------

        final_a_product = domain_product(
            domains_a
        )

        final_b_product = domain_product(
            domains_b
        )

        print(
            "FINAL DOMAIN PRODUCTS"
        )

        print(
            "  A product =",
            final_a_product
        )

        print(
            "  B product =",
            final_b_product
        )

        print()

        # --------------------------------------------------------
        # Reduction factors.
        # --------------------------------------------------------

        raw_a = 1
        raw_b = 1

        for r in R1:
            raw_a *= r

        for s in R2:
            raw_b *= s

        row_reduction = (
            final_a_product
            / raw_a
        )

        col_reduction = (
            final_b_product
            / raw_b
        )

        print(
            "REDUCTION"
        )

        print(
            "  raw A states   =",
            raw_a
        )

        print(
            "  final A states =",
            final_a_product
        )

        print(
            "  A fraction     =",
            f"{row_reduction:.12e}"
        )

        print(
            "  raw B states   =",
            raw_b
        )

        print(
            "  final B states =",
            final_b_product
        )

        print(
            "  B fraction     =",
            f"{col_reduction:.12e}"
        )

        print()

        # --------------------------------------------------------
        # Singleton detection.
        # --------------------------------------------------------

        a_singletons = sum(
            len(x) == 1
            for x in domains_a
        )

        b_singletons = sum(
            len(x) == 1
            for x in domains_b
        )

        print(
            "SINGLETON VARIABLES"
        )

        print(
            "  A singletons =",
            a_singletons,
            "/",
            len(R1)
        )

        print(
            "  B singletons =",
            b_singletons,
            "/",
            len(R2)
        )

        print()

        # --------------------------------------------------------
        # Empty-domain detection.
        # --------------------------------------------------------

        a_empty = [
            i
            for i, x in enumerate(
                domains_a
            )
            if not x
        ]

        b_empty = [
            j
            for j, x in enumerate(
                domains_b
            )
            if not x
        ]

        print(
            "EMPTY DOMAINS"
        )

        print(
            "  A =",
            a_empty
        )

        print(
            "  B =",
            b_empty
        )

        print()

        elapsed = (
            time.perf_counter()
            - t0
        )

        print(
            "SAMPLE RUNTIME =",
            f"{elapsed:.6f} s"
        )

        print()

# ================================================================
# AGGREGATE
# ================================================================

print("=" * 72)
print("EXPERIMENT 155 AGGREGATE")
print("=" * 72)
print()

print(
    "total cells tested =",
    aggregate["cells"]
)

print(
    "cells with reduction =",
    aggregate["cells_with_reduction"]
)

print(
    "empty cell relations =",
    aggregate["cells_empty"]
)

print(
    "true cell pairs surviving =",
    aggregate["true_pair_survived"],
    "/",
    aggregate["cells"]
)

print()

if aggregate["true_pair_survived"] == \
        aggregate["cells"]:

    print(
        "TRUE PAIR SURVIVAL: PASS"
    )

else:

    print(
        "TRUE PAIR SURVIVAL: FAIL"
    )

print()

print(
    "Interpretation:"
)

print(
    "The experiment tests whether the exact ternary C pattern,"
)

print(
    "combined with n mod (r_i*s_j), becomes a strong CSP"
)

print(
    "when shared row residues a_i and column residues b_j"
)

print(
    "are propagated across the entire 5x5 grid."
)

print()

print("=" * 72)
print("FINISHED EXPERIMENT 155")
print("=" * 72)
