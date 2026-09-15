# ================================================================
# START EXPERIMENT 154
# High-dimensional ternary carry-correction rank test
#
# We observed in Experiment 153:
#
#     E = U + V + C
#
# with
#
#     U_ij = floor(k_i*b_j/s_j)
#     V_ij = floor(ell_j*a_i/r_i)
#     C_ij in {0,1,2}
#
# and, surprisingly:
#
#     det(C) = 0
#
# for every 4x4 case tested.
#
# This experiment tests whether that rank deficiency persists on
# a larger 5x5 grid with independent radices.
#
# No factor search.
# No sqrt(n) search.
# No residue enumeration.
#
# ================================================================

import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

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


# Number of independent random samples per bit size.
SAMPLES_PER_SIZE = 5


# ================================================================
# Miller-Rabin
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37,
    ]

    for prime in small_primes:

        if n == prime:
            return True

        if n % prime == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in [
        2, 3, 5, 7, 11, 13, 17
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

        passed = False

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                passed = True
                break

        if not passed:
            return False

    return True


def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

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
# Matrix construction
# ================================================================

def build_data(p, q):

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

    return k, a, ell, b


def build_Q(n):

    Q = []

    for r in R1:

        row = []

        for s in R2:

            row.append(
                n // (r * s)
            )

        Q.append(row)

    return Q


def build_E(Q, k, ell):

    E = []

    for i in range(len(R1)):

        row = []

        for j in range(len(R2)):

            row.append(
                Q[i][j]
                - k[i] * ell[j]
            )

        E.append(row)

    return E


# ================================================================
# Exact U,V,C decomposition
# ================================================================

def build_UVC(
    k,
    a,
    ell,
    b,
):

    U = []
    V = []
    C = []

    size_r = len(R1)
    size_s = len(R2)

    for i in range(size_r):

        U_row = []
        V_row = []
        C_row = []

        r = R1[i]

        for j in range(size_s):

            s = R2[j]

            Uij = (
                k[i] * b[j]
            ) // s

            Vij = (
                ell[j] * a[i]
            ) // r

            numerator = (
                k[i] * b[j] * r
                + ell[j] * a[i] * s
                + a[i] * b[j]
            )

            denominator = r * s

            Eij = (
                numerator // denominator
            )

            Cij = (
                Eij
                - Uij
                - Vij
            )

            U_row.append(Uij)
            V_row.append(Vij)
            C_row.append(Cij)

        U.append(U_row)
        V.append(V_row)
        C.append(C_row)

    return U, V, C


# ================================================================
# Exact determinant
# ================================================================

def determinant(M):

    n = len(M)

    if n == 0:
        return 1

    if any(
        len(row) != n
        for row in M
    ):
        raise ValueError(
            "determinant requires square matrix"
        )

    if n == 1:
        return M[0][0]

    if n == 2:
        return (
            M[0][0] * M[1][1]
            - M[0][1] * M[1][0]
        )

    total = 0

    for col in range(n):

        sub = []

        for i in range(1, n):

            row = []

            for j in range(n):

                if j != col:
                    row.append(
                        M[i][j]
                    )

            sub.append(row)

        sign = (
            1
            if col % 2 == 0
            else -1
        )

        total += (
            sign
            * M[0][col]
            * determinant(sub)
        )

    return total


# ================================================================
# Combinations
# ================================================================

def combinations_indices(n, k):

    result = []

    def recurse(
        start,
        chosen,
    ):

        if len(chosen) == k:

            result.append(
                tuple(chosen)
            )

            return

        remaining = (
            k - len(chosen)
        )

        max_start = (
            n - remaining
        )

        for x in range(
            start,
            max_start + 1
        ):

            recurse(
                x + 1,
                chosen + [x]
            )

    recurse(0, [])

    return result


# ================================================================
# Exact rank
# ================================================================

def exact_rank(M):

    rows = len(M)

    cols = len(M[0])

    max_size = min(
        rows,
        cols,
    )

    row_sets = {}

    col_sets = {}

    for size in range(
        max_size,
        0,
        -1
    ):

        if size not in row_sets:

            row_sets[size] = \
                combinations_indices(
                    rows,
                    size
                )

        if size not in col_sets:

            col_sets[size] = \
                combinations_indices(
                    cols,
                    size
                )

        for rs in row_sets[size]:

            for cs in col_sets[size]:

                sub = []

                for i in rs:

                    row = []

                    for j in cs:

                        row.append(
                            M[i][j]
                        )

                    sub.append(row)

                if determinant(sub) != 0:
                    return size

    return 0


# ================================================================
# Minor statistics
# ================================================================

def minor_statistics(
    M,
    size,
):

    rows = len(M)
    cols = len(M[0])

    if (
        size > rows
        or size > cols
    ):
        return {
            "count": 0,
            "zero": 0,
            "max_abs": 0,
            "distinct": 0,
        }

    row_sets = combinations_indices(
        rows,
        size
    )

    col_sets = combinations_indices(
        cols,
        size
    )

    values = []

    for rs in row_sets:

        for cs in col_sets:

            sub = []

            for i in rs:

                row = []

                for j in cs:

                    row.append(
                        M[i][j]
                    )

                sub.append(row)

            values.append(
                determinant(sub)
            )

    return {
        "count": len(values),
        "zero": sum(
            x == 0
            for x in values
        ),
        "max_abs": max(
            abs(x)
            for x in values
        ),
        "distinct": len(
            set(values)
        ),
    }


# ================================================================
# Norm
# ================================================================

def matrix_norm(M):

    total = 0.0

    for row in M:

        for x in row:

            total += (
                float(x)
                * float(x)
            )

    return math.sqrt(total)


def max_abs(M):

    return max(
        abs(x)
        for row in M
        for x in row
    )


# ================================================================
# Matrix addition
# ================================================================

def add(A, B):

    return [
        [
            A[i][j] + B[i][j]
            for j in range(len(A[0]))
        ]
        for i in range(len(A))
    ]


# ================================================================
# Display
# ================================================================

def print_matrix(name, M):

    print(name)

    for row in M:

        print(
            "["
            + ", ".join(
                str(x)
                for x in row
            )
            + "]"
        )

    print()


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 154")
print("High-dimensional ternary carry-correction rank test")
print("=" * 72)
print()

print("R1 =", R1)
print("R2 =", R2)
print(
    "SAMPLES PER BIT SIZE =",
    SAMPLES_PER_SIZE
)
print()

global_rank_histogram = {}
global_det_zero_count = 0
global_det_nonzero_count = 0
global_c_zero_count = 0
global_c_one_count = 0
global_c_two_count = 0

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

        start = time.perf_counter()

        p, q = random_semiprime(bits)

        n = p * q

        print(
            "-" * 72
        )

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

        # --------------------------------------------------------
        # Build factor data.
        # --------------------------------------------------------

        k, a, ell, b = build_data(
            p,
            q
        )

        Q = build_Q(
            n
        )

        E = build_E(
            Q,
            k,
            ell
        )

        U, V, C = build_UVC(
            k,
            a,
            ell,
            b
        )

        UV = add(
            U,
            V
        )

        reconstructed_E = add(
            UV,
            C
        )

        # --------------------------------------------------------
        # Exact decomposition.
        # --------------------------------------------------------

        decomposition_ok = (
            reconstructed_E == E
        )

        print(
            "E = U + V + C:",
            "PASS"
            if decomposition_ok
            else "FAIL"
        )

        print()

        # --------------------------------------------------------
        # Correction bounds.
        # --------------------------------------------------------

        C_min = min(
            x
            for row in C
            for x in row
        )

        C_max = max(
            x
            for row in C
            for x in row
        )

        C_values = sorted(
            set(
                x
                for row in C
                for x in row
            )
        )

        print(
            "C VALUES =",
            C_values
        )

        print(
            "min(C) =",
            C_min
        )

        print(
            "max(C) =",
            C_max
        )

        print(
            "C IN {0,1,2}:",
            all(
                x in (0, 1, 2)
                for row in C
                for x in row
            )
        )

        print()

        # --------------------------------------------------------
        # Count correction symbols.
        # --------------------------------------------------------

        c0 = sum(
            x == 0
            for row in C
            for x in row
        )

        c1 = sum(
            x == 1
            for row in C
            for x in row
        )

        c2 = sum(
            x == 2
            for row in C
            for x in row
        )

        global_c_zero_count += c0
        global_c_one_count += c1
        global_c_two_count += c2

        print(
            "C COUNTS:"
        )

        print(
            "  0 =",
            c0
        )

        print(
            "  1 =",
            c1
        )

        print(
            "  2 =",
            c2
        )

        print()

        # --------------------------------------------------------
        # Rank.
        # --------------------------------------------------------

        rank_C = exact_rank(
            C
        )

        rank_U = exact_rank(
            U
        )

        rank_V = exact_rank(
            V
        )

        rank_E = exact_rank(
            E
        )

        rank_UV = exact_rank(
            UV
        )

        print(
            "EXACT RANKS"
        )

        print(
            "rank(C)   =",
            rank_C
        )

        print(
            "rank(U)   =",
            rank_U
        )

        print(
            "rank(V)   =",
            rank_V
        )

        print(
            "rank(U+V) =",
            rank_UV
        )

        print(
            "rank(E)   =",
            rank_E
        )

        print()

        global_rank_histogram[
            rank_C
        ] = (
            global_rank_histogram.get(
                rank_C,
                0
            )
            + 1
        )

        # --------------------------------------------------------
        # 2x2 minors.
        # --------------------------------------------------------

        C2 = minor_statistics(
            C,
            2
        )

        C3 = minor_statistics(
            C,
            3
        )

        C4 = minor_statistics(
            C,
            4
        )

        C5 = minor_statistics(
            C,
            5
        )

        print(
            "C MINOR STATISTICS"
        )

        for size, stats in [
            (2, C2),
            (3, C3),
            (4, C4),
            (5, C5),
        ]:

            print(
                f"{size}x{size}: "
                f"count={stats['count']}, "
                f"zero={stats['zero']}, "
                f"max_abs={stats['max_abs']}, "
                f"distinct={stats['distinct']}"
            )

        print()

        # --------------------------------------------------------
        # Full determinant.
        # --------------------------------------------------------

        det_C = determinant(
            C
        )

        print(
            "det(C) =",
            det_C
        )

        if det_C == 0:

            global_det_zero_count += 1

        else:

            global_det_nonzero_count += 1

        print()

        # --------------------------------------------------------
        # Norm.
        # --------------------------------------------------------

        norm_U = matrix_norm(
            U
        )

        norm_V = matrix_norm(
            V
        )

        norm_UV = matrix_norm(
            UV
        )

        norm_C = matrix_norm(
            C
        )

        norm_E = matrix_norm(
            E
        )

        print(
            "NORM SCALE"
        )

        print(
            "||U|| =",
            f"{norm_U:.12e}"
        )

        print(
            "||V|| =",
            f"{norm_V:.12e}"
        )

        print(
            "||U+V|| =",
            f"{norm_UV:.12e}"
        )

        print(
            "||C|| =",
            f"{norm_C:.12e}"
        )

        print(
            "||E|| =",
            f"{norm_E:.12e}"
        )

        if norm_E != 0:

            print(
                "||C|| / ||E|| =",
                f"{norm_C / norm_E:.12e}"
            )

        print()

        # --------------------------------------------------------
        # Print C matrix.
        # --------------------------------------------------------

        print_matrix(
            "C MATRIX",
            C
        )

        # --------------------------------------------------------
        # Row/column pattern.
        # --------------------------------------------------------

        unique_rows = len(
            set(
                tuple(row)
                for row in C
            )
        )

        unique_columns = len(
            set(
                tuple(
                    C[i][j]
                    for i in range(5)
                )
                for j in range(5)
            )
        )

        print(
            "C PATTERN"
        )

        print(
            "unique rows =",
            unique_rows,
            "/ 5"
        )

        print(
            "unique columns =",
            unique_columns,
            "/ 5"
        )

        print()

        # --------------------------------------------------------
        # All rows as vectors and row sums.
        # --------------------------------------------------------

        print(
            "C ROW SUMS"
        )

        for i, row in enumerate(C):

            print(
                f"row {i}:",
                sum(row),
                row
            )

        print()

        # --------------------------------------------------------
        # Check whether identical rows/columns imply the
        # determinant zero. This is diagnostic only.
        # --------------------------------------------------------

        repeated_rows = (
            unique_rows < 5
        )

        repeated_columns = (
            unique_columns < 5
        )

        print(
            "REPEATED-STRUCTURE FLAGS"
        )

        print(
            "repeated C rows =",
            repeated_rows
        )

        print(
            "repeated C columns =",
            repeated_columns
        )

        print()

        # --------------------------------------------------------
        # Exact runtime.
        # --------------------------------------------------------

        elapsed = (
            time.perf_counter()
            - start
        )

        print(
            "SAMPLE RUNTIME =",
            f"{elapsed:.6f} s"
        )

        print()

    print()


# ================================================================
# Aggregate results
# ================================================================

print("=" * 72)
print("EXPERIMENT 154 AGGREGATE")
print("=" * 72)
print()

print(
    "C rank histogram =",
    dict(
        sorted(
            global_rank_histogram.items()
        )
    )
)

print()

print(
    "det(C) = 0 cases =",
    global_det_zero_count
)

print(
    "det(C) != 0 cases =",
    global_det_nonzero_count
)

print()

total_c = (
    global_c_zero_count
    + global_c_one_count
    + global_c_two_count
)

print(
    "GLOBAL C SYMBOL COUNTS"
)

print(
    "C=0 =",
    global_c_zero_count
)

print(
    "C=1 =",
    global_c_one_count
)

print(
    "C=2 =",
    global_c_two_count
)

if total_c:

    print(
        "C=0 fraction =",
        f"{global_c_zero_count / total_c:.12e}"
    )

    print(
        "C=1 fraction =",
        f"{global_c_one_count / total_c:.12e}"
    )

    print(
        "C=2 fraction =",
        f"{global_c_two_count / total_c:.12e}"
    )

print()

# ================================================================
# Final structural conclusion
# ================================================================

if global_det_nonzero_count == 0:

    print(
        "ALL TESTED C MATRICES ARE SINGULAR."
    )

else:

    print(
        "AT LEAST ONE C MATRIX HAS NONZERO DETERMINANT."
    )

print()

print("=" * 72)
print("FINISHED EXPERIMENT 154")
print("=" * 72)
