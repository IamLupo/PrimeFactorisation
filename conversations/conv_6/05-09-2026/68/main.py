# ================================================================
# START EXPERIMENT 151
# 4x4 carry-matrix rank / minor structure
#
# Corrected version:
#   - exact rank now supports rectangular matrices
#   - minor generation supports arbitrary m x n matrices
#   - all previous diagnostics retained
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

    for p in small_primes:

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
# Matrix construction
# ================================================================

def build_matrices(n, p, q):

    K = []
    A = []

    for r in R1:

        K.append(
            p // r
        )

        A.append(
            p % r
        )

    L = []
    B = []

    for s in R2:

        L.append(
            q // s
        )

        B.append(
            q % s
        )

    Q = []
    E = []

    for i, r in enumerate(R1):

        qrow = []
        erow = []

        for j, s in enumerate(R2):

            modulus = r * s

            qij = n // modulus

            eij = (
                qij
                - K[i] * L[j]
            )

            qrow.append(qij)
            erow.append(eij)

        Q.append(qrow)
        E.append(erow)

    return K, A, L, B, Q, E


# ================================================================
# Combinations
# ================================================================

def combinations_indices(n, k):

    out = []

    if k < 0 or k > n:
        return out

    def rec(start, chosen):

        if len(chosen) == k:

            out.append(
                tuple(chosen)
            )

            return

        remaining = k - len(chosen)

        max_value = n - remaining

        for x in range(
            start,
            max_value + 1
        ):

            rec(
                x + 1,
                chosen + [x]
            )

    rec(0, [])

    return out


# ================================================================
# Exact determinant
# ================================================================

def determinant(M):

    n = len(M)

    if n == 0:
        return 1

    if len(M[0]) != n:
        raise ValueError(
            "determinant() requires a square matrix"
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
# Matrix dimensions
# ================================================================

def matrix_shape(M):

    if not M:
        return 0, 0

    return len(M), len(M[0])


# ================================================================
# All k x k minors of a rectangular matrix
# ================================================================

def minor_values(M, k):

    rows, cols = matrix_shape(M)

    if k < 1:
        return []

    if k > rows or k > cols:
        return []

    row_sets = combinations_indices(
        rows,
        k
    )

    col_sets = combinations_indices(
        cols,
        k
    )

    values = []

    for row_indices in row_sets:

        for col_indices in col_sets:

            sub = []

            for i in row_indices:

                row = []

                for j in col_indices:

                    row.append(
                        M[i][j]
                    )

                sub.append(row)

            values.append(
                determinant(sub)
            )

    return values


# ================================================================
# Exact rank of arbitrary rectangular matrix
# ================================================================

def exact_rank(M):

    rows, cols = matrix_shape(M)

    if rows == 0 or cols == 0:
        return 0

    maximum = min(
        rows,
        cols
    )

    for k in range(
        maximum,
        0,
        -1
    ):

        values = minor_values(
            M,
            k
        )

        for value in values:

            if value != 0:
                return k

    return 0


# ================================================================
# Norm
# ================================================================

def frobenius_norm(M):

    total = 0.0

    for row in M:

        for x in row:

            total += (
                float(x)
                * float(x)
            )

    return math.sqrt(total)


def max_abs(M):

    if not M:
        return 0

    return max(
        abs(x)
        for row in M
        for x in row
    )


# ================================================================
# Row differences
# ================================================================

def row_difference(M):

    rows, cols = matrix_shape(M)

    if rows <= 1:
        return []

    out = []

    for i in range(rows - 1):

        row = []

        for j in range(cols):

            row.append(
                M[i + 1][j]
                - M[i][j]
            )

        out.append(row)

    return out


# ================================================================
# Column differences
# ================================================================

def column_difference(M):

    rows, cols = matrix_shape(M)

    if cols <= 1:
        return []

    out = []

    for i in range(rows):

        row = []

        for j in range(cols - 1):

            row.append(
                M[i][j + 1]
                - M[i][j]
            )

        out.append(row)

    return out


# ================================================================
# Matrix display
# ================================================================

def print_matrix(name, M):

    print(name)

    if not M:

        print("[]")
        print()
        return

    for row in M:

        print(
            "["
            + ", ".join(str(x) for x in row)
            + "]"
        )

    print()


# ================================================================
# Print minor statistics
# ================================================================

def print_minor_stats(name, M, k):

    rows, cols = matrix_shape(M)

    values = minor_values(
        M,
        k
    )

    if not values:

        print(
            f"{name}: no {k}x{k} minors "
            f"for {rows}x{cols} matrix"
        )

        return

    zero_count = sum(
        value == 0
        for value in values
    )

    maximum = max(
        abs(value)
        for value in values
    )

    distinct = len(
        set(values)
    )

    print(
        f"{name}: "
        f"{k}x{k} minors = {len(values)}, "
        f"zero = {zero_count}, "
        f"max_abs = {maximum}, "
        f"distinct = {distinct}"
    )


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 151")
print("4x4 carry-matrix rank / minor structure")
print("=" * 72)
print()

print("R1 =", R1)
print("R2 =", R2)
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

    K, A, L, B, Q, E = build_matrices(
        n,
        p,
        q
    )

    K_matrix = [
        [
            K[i] * L[j]
            for j in range(4)
        ]
        for i in range(4)
    ]

    print("TRUE K")
    print(K)
    print()

    print("TRUE L")
    print(L)
    print()

    print("TRUE A")
    print(A)
    print()

    print("TRUE B")
    print(B)
    print()

    print_matrix(
        "Q MATRIX",
        Q
    )

    print_matrix(
        "K MATRIX",
        K_matrix
    )

    print_matrix(
        "E MATRIX",
        E
    )

    # ------------------------------------------------------------
    # Basic exact identity Q = K + E.
    # ------------------------------------------------------------

    identity_ok = True

    for i in range(4):

        for j in range(4):

            if Q[i][j] != (
                K_matrix[i][j]
                + E[i][j]
            ):

                identity_ok = False

    print(
        "Q = K + E:",
        "PASS"
        if identity_ok
        else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # Exact ranks.
    # ------------------------------------------------------------

    rank_Q = exact_rank(Q)
    rank_K = exact_rank(K_matrix)
    rank_E = exact_rank(E)

    print("EXACT RANKS")
    print()

    print(
        "rank(Q) =",
        rank_Q
    )

    print(
        "rank(K) =",
        rank_K
    )

    print(
        "rank(E) =",
        rank_E
    )

    print()

    # ------------------------------------------------------------
    # E rank interpretation.
    # ------------------------------------------------------------

    if rank_E == 1:

        print(
            "E HAS RANK 1"
        )

    elif rank_E < 4:

        print(
            "E HAS LOWER RANK THAN 4"
        )

    else:

        print(
            "E HAS FULL 4x4 RANK"
        )

    print()

    # ------------------------------------------------------------
    # Exact minors.
    # ------------------------------------------------------------

    print("MINOR STATISTICS")
    print()

    print_minor_stats(
        "K",
        K_matrix,
        2
    )

    print_minor_stats(
        "E",
        E,
        2
    )

    print_minor_stats(
        "E",
        E,
        3
    )

    print_minor_stats(
        "E",
        E,
        4
    )

    print_minor_stats(
        "Q",
        Q,
        2
    )

    print_minor_stats(
        "Q",
        Q,
        3
    )

    print_minor_stats(
        "Q",
        Q,
        4
    )

    print()

    # ------------------------------------------------------------
    # Determinants.
    # ------------------------------------------------------------

    det_Q = determinant(Q)
    det_K = determinant(K_matrix)
    det_E = determinant(E)

    print("FULL DETERMINANTS")
    print()

    print(
        "det(Q) =",
        det_Q
    )

    print(
        "det(K) =",
        det_K
    )

    print(
        "det(E) =",
        det_E
    )

    print()

    # ------------------------------------------------------------
    # Determinant ratios.
    # ------------------------------------------------------------

    print("DETERMINANT SCALE")
    print()

    if det_E != 0:

        print(
            "det(Q) / det(E) =",
            f"{det_Q / det_E:.12e}"
        )

    else:

        print(
            "det(Q) / det(E) = undefined"
        )

    emax = max_abs(E)

    if emax > 0:

        print(
            "|det(Q)| / max(|E|)^4 =",
            f"{abs(det_Q) / (emax ** 4):.12e}"
        )

        print(
            "|det(E)| / max(|E|)^4 =",
            f"{abs(det_E) / (emax ** 4):.12e}"
        )

    print()

    # ------------------------------------------------------------
    # Difference matrices.
    # ------------------------------------------------------------

    drow_E = row_difference(E)
    dcol_E = column_difference(E)

    print_matrix(
        "ROW DIFFERENCES OF E",
        drow_E
    )

    print_matrix(
        "COLUMN DIFFERENCES OF E",
        dcol_E
    )

    # ------------------------------------------------------------
    # Rectangular ranks.
    # ------------------------------------------------------------

    rank_drow = exact_rank(
        drow_E
    )

    rank_dcol = exact_rank(
        dcol_E
    )

    print("DIFFERENCE RANKS")
    print()

    print(
        "shape(row-difference(E)) =",
        matrix_shape(drow_E)
    )

    print(
        "shape(column-difference(E)) =",
        matrix_shape(dcol_E)
    )

    print(
        "rank(row-difference(E)) =",
        rank_drow
    )

    print(
        "rank(column-difference(E)) =",
        rank_dcol
    )

    print()

    # ------------------------------------------------------------
    # Difference minors.
    # ------------------------------------------------------------

    print("DIFFERENCE MINOR STATISTICS")
    print()

    print_minor_stats(
        "row-difference(E)",
        drow_E,
        2
    )

    print_minor_stats(
        "row-difference(E)",
        drow_E,
        3
    )

    print_minor_stats(
        "column-difference(E)",
        dcol_E,
        2
    )

    print_minor_stats(
        "column-difference(E)",
        dcol_E,
        3
    )

    print()

    # ------------------------------------------------------------
    # Norm scales.
    # ------------------------------------------------------------

    norm_Q = frobenius_norm(Q)
    norm_K = frobenius_norm(K_matrix)
    norm_E = frobenius_norm(E)
    norm_drow = frobenius_norm(
        drow_E
    )
    norm_dcol = frobenius_norm(
        dcol_E
    )

    print("NORM SCALE")
    print()

    print(
        "||Q|| =",
        f"{norm_Q:.12e}"
    )

    print(
        "||K|| =",
        f"{norm_K:.12e}"
    )

    print(
        "||E|| =",
        f"{norm_E:.12e}"
    )

    if norm_Q != 0:

        print(
            "||E|| / ||Q|| =",
            f"{norm_E / norm_Q:.12e}"
        )

    if norm_E != 0:

        print(
            "||row-diff(E)|| / ||E|| =",
            f"{norm_drow / norm_E:.12e}"
        )

        print(
            "||col-diff(E)|| / ||E|| =",
            f"{norm_dcol / norm_E:.12e}"
        )

    print()

    # ------------------------------------------------------------
    # Carry scale.
    # ------------------------------------------------------------

    print("CARRY SCALE")
    print()

    print(
        "max |E| =",
        max_abs(E)
    )

    print(
        "max radix R1 =",
        max(R1)
    )

    print(
        "max radix R2 =",
        max(R2)
    )

    print(
        "max |E| / sqrt(n) =",
        f"{max_abs(E) / math.sqrt(n):.12e}"
    )

    print()

    # ------------------------------------------------------------
    # Minor diversity.
    # ------------------------------------------------------------

    E2 = minor_values(
        E,
        2
    )

    E3 = minor_values(
        E,
        3
    )

    E4 = minor_values(
        E,
        4
    )

    print("MINOR DIVERSITY")
    print()

    print(
        "distinct E 2x2 minors =",
        len(set(E2)),
        "/",
        len(E2)
    )

    print(
        "distinct E 3x3 minors =",
        len(set(E3)),
        "/",
        len(E3)
    )

    print(
        "distinct E 4x4 minors =",
        len(set(E4)),
        "/",
        len(E4)
    )

    print()

    # ------------------------------------------------------------
    # Look for unusually strong structure in E.
    # ------------------------------------------------------------

    print("STRUCTURAL FLAGS")
    print()

    print(
        "E rank < 4:",
        rank_E < 4
    )

    print(
        "row-difference rank < 3:",
        rank_drow < 3
    )

    print(
        "column-difference rank < 3:",
        rank_dcol < 3
    )

    print(
        "det(E) == 0:",
        det_E == 0
    )

    print(
        "all E 3x3 minors zero:",
        all(x == 0 for x in E3)
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
print("FINISHED EXPERIMENT 151")
print("=" * 72)