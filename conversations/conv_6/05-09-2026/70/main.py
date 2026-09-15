# ================================================================
# START EXPERIMENT 153
# Exact carry decomposition:
#
#     E = U + V + C
#
# where
#
#     U_ij = floor(k_i * b_j / s_j)
#     V_ij = floor(ell_j * a_i / r_i)
#
# and
#
#     C_ij in {0,1,2}.
#
# Goal:
#
#   Determine whether the apparently full-rank E matrix is really
#   composed of two structured floor matrices plus a tiny
#   correction matrix.
#
# No factor search.
# No residue enumeration.
# No sqrt(n) scan.
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
]

R2 = [
    19,
    37,
    61,
    73,
]

BIT_SIZES = [
    30,
    36,
    42,
    48,
    54,
]


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

        x = pow(a, d, n)

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
# Factor decomposition
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
# Quotient/carry matrices
# ================================================================

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

    for i in range(4):

        row = []

        for j in range(4):

            row.append(
                Q[i][j]
                - k[i] * ell[j]
            )

        E.append(row)

    return E


# ================================================================
# Exact decomposition
# ================================================================

def build_U_V_C(k, a, ell, b):

    U = []
    V = []
    C = []

    # Exact rational values are computed as integer numerators,
    # avoiding floating-point arithmetic.

    for i, r in enumerate(R1):

        U_row = []
        V_row = []
        C_row = []

        for j, s in enumerate(R2):

            # ----------------------------------------------------
            # U = floor(k_i*b_j/s_j)
            # ----------------------------------------------------

            Uij = (
                k[i] * b[j]
            ) // s

            # ----------------------------------------------------
            # V = floor(ell_j*a_i/r_i)
            # ----------------------------------------------------

            Vij = (
                ell[j] * a[i]
            ) // r

            U_row.append(Uij)
            V_row.append(Vij)

            # ----------------------------------------------------
            # Full E:
            #
            # floor(
            #       kb/s
            #     + ell*a/r
            #     + ab/(rs)
            # )
            #
            # Compute it exactly using common denominator rs.
            # ----------------------------------------------------

            numerator = (
                k[i] * b[j] * r
                + ell[j] * a[i] * s
                + a[i] * b[j]
            )

            denominator = (
                r * s
            )

            Eij = (
                numerator
                // denominator
            )

            Cij = (
                Eij
                - Uij
                - Vij
            )

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

    if any(len(row) != n for row in M):
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

    def rec(start, chosen):

        if len(chosen) == k:

            result.append(
                tuple(chosen)
            )

            return

        remaining = k - len(chosen)

        max_start = n - remaining

        for x in range(
            start,
            max_start + 1
        ):

            rec(
                x + 1,
                chosen + [x]
            )

    rec(0, [])

    return result


# ================================================================
# Minor values
# ================================================================

def minor_values(M, size):

    rows = len(M)
    cols = len(M[0])

    if size > rows or size > cols:
        return []

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

    return values


# ================================================================
# Rank
# ================================================================

def exact_rank(M):

    rows = len(M)
    cols = len(M[0])

    for size in range(
        min(rows, cols),
        0,
        -1
    ):

        values = minor_values(
            M,
            size
        )

        if any(
            value != 0
            for value in values
        ):
            return size

    return 0


# ================================================================
# Norm
# ================================================================

def norm(M):

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
# Matrix arithmetic
# ================================================================

def add_matrices(A, B):

    return [
        [
            A[i][j] + B[i][j]
            for j in range(
                len(A[0])
            )
        ]
        for i in range(len(A))
    ]


def subtract_matrices(A, B):

    return [
        [
            A[i][j] - B[i][j]
            for j in range(
                len(A[0])
            )
        ]
        for i in range(len(A))
    ]


# ================================================================
# Matrix display
# ================================================================

def print_matrix(name, M):

    print(name)

    for row in M:

        print(
            "["
            + ", ".join(str(x) for x in row)
            + "]"
        )

    print()


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 153")
print("Exact carry decomposition E = U + V + C")
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

for bits in BIT_SIZES:

    start = time.perf_counter()

    print("=" * 72)
    print(
        f"GENERATING {bits}-BIT SEMIPRIME"
    )
    print("=" * 72)
    print()

    p, q = random_semiprime(bits)
    n = p * q

    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print(
        "n      =",
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

    # ------------------------------------------------------------
    # Hidden factor data.
    # ------------------------------------------------------------

    k, a, ell, b = factor_data(
        p,
        q
    )

    print(
        "k =",
        k
    )

    print(
        "a =",
        a
    )

    print(
        "ell =",
        ell
    )

    print(
        "b =",
        b
    )

    print()

    # ------------------------------------------------------------
    # Q and E.
    # ------------------------------------------------------------

    Q = build_Q(n)

    E = build_E(
        Q,
        k,
        ell
    )

    # ------------------------------------------------------------
    # Decomposition.
    # ------------------------------------------------------------

    U, V, C = build_U_V_C(
        k,
        a,
        ell,
        b
    )

    UV = add_matrices(
        U,
        V
    )

    reconstructed_E = add_matrices(
        UV,
        C
    )

    print_matrix(
        "E MATRIX",
        E
    )

    print_matrix(
        "U MATRIX = floor(k_i*b_j/s_j)",
        U
    )

    print_matrix(
        "V MATRIX = floor(ell_j*a_i/r_i)",
        V
    )

    print_matrix(
        "C MATRIX = E-U-V",
        C
    )

    # ------------------------------------------------------------
    # Exact identity check.
    # ------------------------------------------------------------

    identity_ok = (
        reconstructed_E == E
    )

    print(
        "E = U + V + C:",
        "PASS"
        if identity_ok
        else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # Ternary correction test.
    # ------------------------------------------------------------

    c_values = sorted(
        set(
            x
            for row in C
            for x in row
        )
    )

    c_min = min(
        x
        for row in C
        for x in row
    )

    c_max = max(
        x
        for row in C
        for x in row
    )

    print(
        "C VALUES =",
        c_values
    )

    print(
        "min(C) =",
        c_min
    )

    print(
        "max(C) =",
        c_max
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

    # ------------------------------------------------------------
    # Value frequencies.
    # ------------------------------------------------------------

    counts = {
        0: 0,
        1: 0,
        2: 0,
    }

    for row in C:

        for x in row:

            if x not in counts:
                counts[x] = 0

            counts[x] += 1

    print(
        "C VALUE COUNTS =",
        counts
    )

    print()

    # ------------------------------------------------------------
    # Ranks.
    # ------------------------------------------------------------

    rank_E = exact_rank(E)
    rank_U = exact_rank(U)
    rank_V = exact_rank(V)
    rank_C = exact_rank(C)
    rank_UV = exact_rank(UV)

    print(
        "EXACT RANKS"
    )

    print(
        "rank(E)  =",
        rank_E
    )

    print(
        "rank(U)  =",
        rank_U
    )

    print(
        "rank(V)  =",
        rank_V
    )

    print(
        "rank(U+V) =",
        rank_UV
    )

    print(
        "rank(C)  =",
        rank_C
    )

    print()

    # ------------------------------------------------------------
    # Minor statistics.
    # ------------------------------------------------------------

    print(
        "MINOR STATISTICS"
    )

    print()

    for name, matrix in [
        ("U", U),
        ("V", V),
        ("U+V", UV),
        ("C", C),
        ("E", E),
    ]:

        m2 = minor_values(
            matrix,
            2
        )

        m3 = minor_values(
            matrix,
            3
        )

        m4 = minor_values(
            matrix,
            4
        )

        print(
            f"{name}: "
            f"2x2 zero={sum(x == 0 for x in m2)}/{len(m2)}, "
            f"3x3 zero={sum(x == 0 for x in m3)}/{len(m3)}, "
            f"4x4 zero={sum(x == 0 for x in m4)}/{len(m4)}"
        )

    print()

    # ------------------------------------------------------------
    # Determinants.
    # ------------------------------------------------------------

    det_E = determinant(E)
    det_U = determinant(U)
    det_V = determinant(V)
    det_UV = determinant(UV)
    det_C = determinant(C)

    print(
        "DETERMINANTS"
    )

    print(
        "det(U)   =",
        det_U
    )

    print(
        "det(V)   =",
        det_V
    )

    print(
        "det(U+V) =",
        det_UV
    )

    print(
        "det(C)   =",
        det_C
    )

    print(
        "det(E)   =",
        det_E
    )

    print()

    # ------------------------------------------------------------
    # Magnitude hierarchy.
    # ------------------------------------------------------------

    nU = norm(U)
    nV = norm(V)
    nC = norm(C)
    nE = norm(E)
    nUV = norm(UV)

    print(
        "NORM HIERARCHY"
    )

    print()

    print(
        "||U|| =",
        f"{nU:.12e}"
    )

    print(
        "||V|| =",
        f"{nV:.12e}"
    )

    print(
        "||U+V|| =",
        f"{nUV:.12e}"
    )

    print(
        "||C|| =",
        f"{nC:.12e}"
    )

    print(
        "||E|| =",
        f"{nE:.12e}"
    )

    if nE != 0:

        print(
            "||C|| / ||E|| =",
            f"{nC / nE:.12e}"
        )

        print(
            "||C|| / ||U+V|| =",
            f"{nC / nUV:.12e}"
        )

    print()

    # ------------------------------------------------------------
    # Max-entry hierarchy.
    # ------------------------------------------------------------

    print(
        "MAX ENTRY HIERARCHY"
    )

    print()

    print(
        "max |U| =",
        max_abs(U)
    )

    print(
        "max |V| =",
        max_abs(V)
    )

    print(
        "max |C| =",
        max_abs(C)
    )

    print(
        "max |E| =",
        max_abs(E)
    )

    print()

    # ------------------------------------------------------------
    # Row and column structure of C.
    # ------------------------------------------------------------

    print(
        "C ROW VECTORS"
    )

    for i, row in enumerate(C):

        print(
            f"row {i}:",
            row,
            "sum=",
            sum(row)
        )

    print()

    print(
        "C COLUMN VECTORS"
    )

    for j in range(4):

        col = [
            C[i][j]
            for i in range(4)
        ]

        print(
            f"column {j}:",
            col,
            "sum=",
            sum(col)
        )

    print()

    # ------------------------------------------------------------
    # Check whether C is constant / separable / low-complexity.
    # ------------------------------------------------------------

    unique_rows = len(
        set(
            tuple(row)
            for row in C
        )
    )

    unique_cols = len(
        set(
            tuple(
                C[i][j]
                for i in range(4)
            )
            for j in range(4)
        )
    )

    print(
        "C PATTERN COMPLEXITY"
    )

    print()

    print(
        "unique C rows =",
        unique_rows,
        "/ 4"
    )

    print(
        "unique C columns =",
        unique_cols,
        "/ 4"
    )

    print()

    # ------------------------------------------------------------
    # Verify the exact fractional interpretation.
    #
    # C = floor(frac(U-term) + frac(V-term) + ab/(rs)).
    #
    # We print the numerator components to make the source of
    # the ternary correction explicit.
    # ------------------------------------------------------------

    print(
        "FRACTIONAL CARRY SOURCE"
    )

    print()

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            x_rem = (
                k[i] * b[j]
            ) % s

            y_rem = (
                ell[j] * a[i]
            ) % r

            numerator = (
                x_rem * r
                + y_rem * s
                + a[i] * b[j]
            )

            denominator = (
                r * s
            )

            predicted_C = (
                numerator
                // denominator
            )

            print(
                f"({i},{j}) "
                f"x_rem={x_rem:3d} "
                f"y_rem={y_rem:3d} "
                f"ab={a[i] * b[j]:4d} "
                f"C={C[i][j]} "
                f"predicted={predicted_C}"
            )

    print()

    # ------------------------------------------------------------
    # Runtime.
    # ------------------------------------------------------------

    elapsed = (
        time.perf_counter()
        - start
    )

    print(
        "TOTAL RUNTIME =",
        f"{elapsed:.6f} s"
    )

    print()

    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 153")
print("=" * 72)
