# ================================================================
# START EXPERIMENT 157
# Global quotient-residue consistency on strongest C cells
#
# Motivation:
#
# Experiment 156 treated each cell independently.
#
# For a cell (i,j), we allowed:
#
#     x_ij = k_i mod s_j
#     y_ij = ell_j mod r_i
#
# to be independently chosen.
#
# But this is NOT globally correct.
#
# All x_ij belonging to the same row must originate from the
# SAME integer k_i, and all y_ij belonging to the same column
# must originate from the SAME ell_j.
#
# An even stronger representation is to work directly with
#
#     p mod M
#     q mod M
#
# where M contains all radices involved in a selected subgrid.
#
# Since
#
#     p*q = n (mod M),
#
# once p is invertible modulo M,
#
#     q = n * p^{-1} (mod M).
#
# Therefore a single p-residue determines q-residue.
#
# We then calculate the ACTUAL residues
#
#     a_i = p mod r_i
#     b_j = q mod s_j
#
# and the exact C values for the selected cells.
#
# This experiment asks:
#
#     How many globally consistent p mod M states
#     reproduce the oracle C pattern?
#
# We compare:
#
#     independent cell relation count
#
# against
#
#     globally consistent subgrid count.
#
# No sqrt(n) factor search.
# No p/q search.
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
# PRIMALITY
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

        x = pow(
            a,
            d,
            n
        )

        if x == 1 or x == n - 1:
            continue

        composite = True

        for _ in range(s - 1):

            x = (
                x * x
            ) % n

            if x == n - 1:

                composite = False
                break

        if composite:

            return False

    return True


def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= (
            1 << (bits - 1)
        )

        if is_probable_prime(x):

            return x


def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p != q:

            return p, q


# ================================================================
# TRUE C
# ================================================================

def build_C(
    p,
    q,
):

    k = []
    a = []

    for r in R1:

        k.append(
            p // r
        )

        a.append(
            p % r
        )

    ell = []
    b = []

    for s in R2:

        ell.append(
            q // s
        )

        b.append(
            q % s
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

            cij = (
                numerator
                // (r * s)
            )

            row.append(cij)

        C.append(row)

    return C, a, b


# ================================================================
# PRIME PRODUCT
# ================================================================

def product(values):

    result = 1

    for x in values:

        result *= x

    return result


# ================================================================
# SELECT STRONGEST CELLS
#
# Instead of arbitrary cells, use the oracle C values and choose
# cells that have historically produced strong local reductions.
#
# We first calculate the exact global consistency score for every
# possible 2x2 row/column block.
# ================================================================

def all_2x2_blocks():

    blocks = []

    for i0 in range(len(R1)):

        for i1 in range(
            i0 + 1,
            len(R1)
        ):

            for j0 in range(len(R2)):

                for j1 in range(
                    j0 + 1,
                    len(R2)
                ):

                    cells = [
                        (i0, j0),
                        (i0, j1),
                        (i1, j0),
                        (i1, j1),
                    ]

                    blocks.append(
                        cells
                    )

    return blocks


def block_modulus(cells):

    rows = sorted(
        set(
            i for i, j in cells
        )
    )

    cols = sorted(
        set(
            j for i, j in cells
        )
    )

    row_modulus = product(
        R1[i]
        for i in rows
    )

    col_modulus = product(
        R2[j]
        for j in cols
    )

    return (
        row_modulus
        * col_modulus
    )


# ================================================================
# GLOBAL SUBGRID ENUMERATION
#
# For a selected row/column block:
#
#     M = product(r_i) * product(s_j)
#
# Enumerate p mod M.
#
# q mod M is then forced by:
#
#     q = n * p^{-1} mod M
#
# whenever gcd(p,M)=1.
#
# Then calculate C directly from p_mod_M and q_mod_M.
#
# Because all factors of M are prime and the generated factors are
# far larger than the radices, the true p and q are invertible with
# overwhelming probability; nevertheless gcd is checked explicitly.
# ================================================================

def global_subgrid_states(
    n,
    C,
    cells,
):

    rows = sorted(
        set(
            i for i, j in cells
        )
    )

    cols = sorted(
        set(
            j for i, j in cells
        )
    )

    M_r = product(
        R1[i]
        for i in rows
    )

    M_s = product(
        R2[j]
        for j in cols
    )

    M = M_r * M_s

    survivors = []

    tested = 0
    invertible = 0

    # ------------------------------------------------------------
    # We only need p residues modulo M.
    #
    # p itself is irrelevant; p mod M determines all required
    # a_i and k_i mod s_j because M contains r_i and s_j.
    # ------------------------------------------------------------

    for pmod in range(1, M):

        tested += 1

        if math.gcd(
            pmod,
            M
        ) != 1:

            continue

        invertible += 1

        qmod = (
            n
            * pow(
                pmod,
                -1,
                M
            )
        ) % M

        good = True

        for i, j in cells:

            r = R1[i]
            s = R2[j]

            # ----------------------------------------------------
            # Recover a_i and b_j.
            # ----------------------------------------------------

            a = pmod % r
            b = qmod % s

            # ----------------------------------------------------
            # Recover k_i mod s.
            #
            # Since p = r*k + a:
            #
            #     k = (p-a)/r
            #
            # We do this modulo s by obtaining p modulo r*s.
            # Because M is a multiple of r*s, pmod contains
            # exactly the needed information.
            # ----------------------------------------------------

            p_rs = pmod % (
                r * s
            )

            k_mod_s = (
                (
                    p_rs - a
                )
                // r
            ) % s

            # ----------------------------------------------------
            # Recover ell_j mod r.
            # ----------------------------------------------------

            q_rs = qmod % (
                r * s
            )

            ell_mod_r = (
                (
                    q_rs - b
                )
                // s
            ) % r

            beta = (
                k_mod_s
                * b
            ) % s

            alpha = (
                ell_mod_r
                * a
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

            if cij != C[i][j]:

                good = False
                break

        if good:

            survivors.append(
                (
                    pmod,
                    qmod
                )
            )

    return (
        M,
        tested,
        invertible,
        survivors,
    )


# ================================================================
# LOCAL CELL RELATION COUNT
#
# This deliberately reproduces the weaker Experiment 156 notion:
# each cell is checked independently while allowing arbitrary
# hidden quotient residues.
#
# This gives us a baseline.
# ================================================================

def local_cell_relation_size(
    n,
    r,
    s,
    c,
):

    count = 0

    for a in range(r):

        for b in range(s):

            nmod = (
                n
                % (r * s)
            )

            target = (
                nmod
                - a * b
            ) % (
                r * s
            )

            # Need existence of k mod s.
            possible_k = False

            for x in range(s):

                beta = (
                    x * b
                ) % s

                if (
                    r * beta
                    - target
                ) % s == 0:

                    possible_k = True
                    break

            if not possible_k:

                continue

            # Need existence of ell mod r.
            possible_ell = False

            for y in range(r):

                alpha = (
                    y * a
                ) % r

                if (
                    s * alpha
                    - target
                ) % r == 0:

                    possible_ell = True
                    break

            if not possible_ell:

                continue

            # Check whether some x,y creates C=c.
            found = False

            for x in range(s):

                beta = (
                    x * b
                ) % s

                for y in range(r):

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

                    if cij == c:

                        found = True
                        break

                if found:
                    break

            if found:
                count += 1

    return count


# ================================================================
# BLOCK INFORMATION
# ================================================================

def print_block(
    C,
    cells,
):

    print(
        "cells =",
        cells
    )

    print(
        "C values =",
        [
            C[i][j]
            for i, j in cells
        ]
    )

    rows = sorted(
        set(i for i, j in cells)
    )

    cols = sorted(
        set(j for i, j in cells)
    )

    print(
        "rows =",
        rows
    )

    print(
        "cols =",
        cols
    )

    M = block_modulus(
        cells
    )

    print(
        "M =",
        M
    )

    print()


# ================================================================
# MAIN
# ================================================================

print("=" * 72)
print("START EXPERIMENT 157")
print("Global quotient-residue consistency")
print("=" * 72)
print()

print(
    "R1 =",
    R1
)

print(
    "R2 =",
    R2
)

print()

blocks = all_2x2_blocks()

print(
    "2x2 blocks =",
    len(blocks)
)

print()

# Aggregate statistics.
total_samples = 0
total_blocks = 0
zero_survivor_blocks = 0
true_survivor_blocks = 0

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
            q
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
            "n =",
            n
        )

        print(
            "true p =",
            p
        )

        print(
            "true q =",
            q
        )

        print()

        print(
            "C MATRIX"
        )

        for row in C:

            print(row)

        print()

        # --------------------------------------------------------
        # Evaluate every 2x2 block.
        #
        # There are only 100 blocks, but M can be around
        # 500k, so this is still manageable at these dimensions.
        # --------------------------------------------------------

        block_results = []

        for cells in blocks:

            # ----------------------------------------------------
            # Fast pre-score:
            #
            # Sum of individual information values.
            #
            # This determines which blocks are most promising.
            # ----------------------------------------------------

            local_info = 0.0
            local_count = 0

            for i, j in cells:

                count = local_cell_relation_size(
                    n,
                    R1[i],
                    R2[j],
                    C[i][j]
                )

                raw = (
                    R1[i]
                    * R2[j]
                )

                local_count += count

                if count > 0:

                    local_info += (
                        math.log2(
                            raw
                            / count
                        )
                    )

            block_results.append(
                (
                    local_info,
                    cells,
                    local_count
                )
            )

        block_results.sort(
            key=lambda x: -x[0]
        )

        # --------------------------------------------------------
        # To keep runtime controlled, perform full global
        # enumeration only on the strongest 10 blocks.
        # --------------------------------------------------------

        print(
            "TOP GLOBAL-CONSISTENCY BLOCKS"
        )

        print(
            "rank  cells                     "
            "M        local-info"
        )

        print(
            "-" * 72
        )

        for rank, item in enumerate(
            block_results[:10],
            start=1
        ):

            info, cells, local_count = item

            M = block_modulus(
                cells
            )

            print(
                f"{rank:4d}  "
                f"{str(cells):25s} "
                f"{M:8d}   "
                f"{info:.6f}"
            )

        print()

        # --------------------------------------------------------
        # Run global consistency on top 10.
        # --------------------------------------------------------

        for rank, item in enumerate(
            block_results[:10],
            start=1
        ):

            info, cells, local_count = item

            print(
                "-" * 72
            )

            print(
                f"GLOBAL BLOCK {rank}/10"
            )

            print_block(
                C,
                cells
            )

            print(
                "Independent local compatible "
                "cell-state total =",
                local_count
            )

            M, tested, invertible, survivors = (
                global_subgrid_states(
                    n,
                    C,
                    cells
                )
            )

            total_samples += 1
            total_blocks += 1

            print(
                "GLOBAL ENUMERATION"
            )

            print(
                "  M =",
                M
            )

            print(
                "  p residues tested =",
                tested
            )

            print(
                "  invertible p residues =",
                invertible
            )

            print(
                "  globally consistent survivors =",
                len(survivors)
            )

            # ----------------------------------------------------
            # True state check.
            #
            # p mod M and q mod M must appear among survivors.
            # ----------------------------------------------------

            true_pmod = p % M
            true_qmod = q % M

            true_found = (
                (
                    true_pmod,
                    true_qmod
                )
                in survivors
            )

            print(
                "  true residue state survived =",
                true_found
            )

            if true_found:

                true_survivor_blocks += 1

            if len(survivors) == 0:

                zero_survivor_blocks += 1

            # ----------------------------------------------------
            # Print survivors for small sets.
            # ----------------------------------------------------

            if len(survivors) <= 20:

                print(
                    "  survivor states:"
                )

                for pp, qq in survivors:

                    marker = ""

                    if (
                        pp == true_pmod
                        and
                        qq == true_qmod
                    ):

                        marker = " <-- TRUE"

                    print(
                        "    pmod =",
                        pp,
                        "qmod =",
                        qq,
                        marker
                    )

            else:

                print(
                    "  first survivor states:"
                )

                for pp, qq in survivors[:10]:

                    marker = ""

                    if (
                        pp == true_pmod
                        and
                        qq == true_qmod
                    ):

                        marker = " <-- TRUE"

                    print(
                        "    pmod =",
                        pp,
                        "qmod =",
                        qq,
                        marker
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
print("EXPERIMENT 157 AGGREGATE")
print("=" * 72)
print()

print(
    "global blocks evaluated =",
    total_blocks
)

print(
    "zero-survivor blocks =",
    zero_survivor_blocks
)

print(
    "true residue surviving blocks =",
    true_survivor_blocks
)

print()

print(
    "Interpretation:"
)

print(
    "Experiment 156 allowed each cell to choose its hidden"
)

print(
    "quotient residues independently."
)

print(
    "Experiment 157 instead forces every selected 2x2 block"
)

print(
    "to come from one globally consistent p mod M and q mod M."
)

print()

print(
    "The critical comparison is:"
)

print(
    "    local cell compatibility"
)

print(
    "versus"
)

print(
    "    globally consistent residue states"
)

print()

print("=" * 72)
print("FINISHED EXPERIMENT 157")
print("=" * 72)
