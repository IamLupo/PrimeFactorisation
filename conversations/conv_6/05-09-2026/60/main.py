# ================================================================
# START EXPERIMENT 143
# Nested-radix Q matrix = base-B digit information of n
# Exact Hankel / digit-collapse test
# ================================================================

import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

BASE = 7

R1 = [
    BASE,
    BASE**2,
    BASE**3,
]

R2 = [
    BASE,
    BASE**2,
    BASE**3,
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

    witnesses = [
        2, 3, 5, 7, 11, 13, 17
    ]

    for a in witnesses:

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

        x |= 1 << (bits - 1)
        x |= 1

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
# Base-B digits
# ================================================================

def base_digits(x, base, count=None):

    if x == 0:
        return [0]

    digits = []

    y = x

    while y:

        digits.append(y % base)
        y //= base

    if count is not None:

        while len(digits) < count:
            digits.append(0)

    return digits


def reconstruct_digits(digits, base):

    value = 0
    power = 1

    for d in digits:

        value += d * power
        power *= base

    return value


# ================================================================
# Nested quotient matrix
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


# ================================================================
# Nested carry matrix
# ================================================================

def build_K_E(n, p, q):

    K = []
    A = []

    for r in R1:

        k = p // r
        a = p % r

        K.append(k)
        A.append(a)

    L = []
    B = []

    for s in R2:

        ell = q // s
        b = q % s

        L.append(ell)
        B.append(b)

    Q = build_Q(n)

    E = []

    for i in range(len(R1)):

        row = []

        for j in range(len(R2)):

            row.append(
                Q[i][j] - K[i] * L[j]
            )

        E.append(row)

    return K, A, L, B, Q, E


# ================================================================
# Exact Hankel test
# ================================================================

def hankel_expected(n):

    """
    Since

        R1[i] = B^(i+1)
        R2[j] = B^(j+1),

    we have

        R1[i] R2[j] = B^(i+j+2).

    Therefore

        Q[i][j] = floor(n / B^(i+j+2)).

    """

    Q_expected = []

    for i in range(3):

        row = []

        for j in range(3):

            exponent = i + j + 2

            modulus = BASE ** exponent

            row.append(
                n // modulus
            )

        Q_expected.append(row)

    return Q_expected


# ================================================================
# Exact quotient-digit extraction
# ================================================================

def quotient_digit(n, position):

    """
    Base-B digit at position 'position':

        d_position
          = floor(n/B^position)
            - B floor(n/B^(position+1))
    """

    left = n // (BASE ** position)
    right = n // (BASE ** (position + 1))

    return left - BASE * right


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
print("START EXPERIMENT 143")
print("Nested-radix Q matrix = base-B digit information of n")
print("Exact Hankel / digit-collapse test")
print("=" * 72)
print()

print("BASE =", BASE)
print("R1   =", R1)
print("R2   =", R2)
print()


for bits in BIT_SIZES:

    start = time.perf_counter()

    print("=" * 72)
    print(f"GENERATING {bits}-BIT SEMIPRIME")
    print("=" * 72)
    print()

    p, q = random_semiprime(bits)
    n = p * q

    print("-" * 72)
    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)
    print()

    # ------------------------------------------------------------
    # Factor decomposition
    # ------------------------------------------------------------

    K, A, L, Bdigits, Q, E = build_K_E(
        n,
        p,
        q
    )

    print("TRUE K")
    print(K)
    print()

    print("TRUE A")
    print(A)
    print()

    print("TRUE L")
    print(L)
    print()

    print("TRUE B")
    print(Bdigits)
    print()

    print_matrix(
        "Q MATRIX",
        Q
    )

    print_matrix(
        "E MATRIX",
        E
    )

    # ------------------------------------------------------------
    # Exact Hankel reconstruction
    # ------------------------------------------------------------

    Q_expected = hankel_expected(n)

    print_matrix(
        "Q RECONSTRUCTED DIRECTLY FROM n AND B",
        Q_expected
    )

    hankel_exact = (
        Q == Q_expected
    )

    print(
        "Q = floor(n / B^(i+j+2)):",
        "PASS"
        if hankel_exact
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Explicit Hankel equalities
    # ------------------------------------------------------------

    print("HANKEL EQUALITIES")
    print()

    checks = [
        ("Q01 = Q10", Q[0][1], Q[1][0]),
        ("Q02 = Q11", Q[0][2], Q[1][1]),
        ("Q11 = Q20", Q[1][1], Q[2][0]),
        ("Q12 = Q21", Q[1][2], Q[2][1]),
    ]

    all_hankel = True

    for name, lhs, rhs in checks:

        passed = lhs == rhs

        if not passed:
            all_hankel = False

        print(
            f"{name}: "
            f"{lhs} == {rhs} -> "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()

    print(
        "ALL HANKEL EQUALITIES:",
        "PASS"
        if all_hankel
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Base-B digits of n
    # ------------------------------------------------------------

    digits = base_digits(
        n,
        BASE
    )

    print("BASE-B DIGITS OF n")
    print(digits)
    print()

    # ------------------------------------------------------------
    # Difference identities.
    #
    # For every i,j:
    #
    # Q[i][j] - B Q[i+1][j]
    #
    # equals digit i+j+2 of n,
    # whenever i+1 exists.
    # ------------------------------------------------------------

    print("VERTICAL QUOTIENT DIGIT EXTRACTION")
    print()

    vertical_checks = []

    for i in range(2):

        for j in range(3):

            position = i + j + 2

            observed = (
                Q[i][j]
                - BASE * Q[i + 1][j]
            )

            expected = quotient_digit(
                n,
                position
            )

            vertical_checks.append(
                observed == expected
            )

            print(
                f"(i,j)=({i},{j}) "
                f"position={position}: "
                f"{observed} "
                f"vs digit={expected} "
                f"[{'PASS' if observed == expected else 'FAIL'}]"
            )

    print()

    print(
        "VERTICAL DIGIT EXTRACTION:",
        "PASS"
        if all(vertical_checks)
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Horizontal extraction.
    # ------------------------------------------------------------

    print("HORIZONTAL QUOTIENT DIGIT EXTRACTION")
    print()

    horizontal_checks = []

    for i in range(3):

        for j in range(2):

            position = i + j + 2

            observed = (
                Q[i][j]
                - BASE * Q[i][j + 1]
            )

            expected = quotient_digit(
                n,
                position
            )

            horizontal_checks.append(
                observed == expected
            )

            print(
                f"(i,j)=({i},{j}) "
                f"position={position}: "
                f"{observed} "
                f"vs digit={expected} "
                f"[{'PASS' if observed == expected else 'FAIL'}]"
            )

    print()

    print(
        "HORIZONTAL DIGIT EXTRACTION:",
        "PASS"
        if all(horizontal_checks)
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Diagonal compression.
    #
    # A 3x3 Q matrix should contain only five distinct values:
    #
    # q_2, q_3, q_4, q_5, q_6.
    # ------------------------------------------------------------

    diagonal_values = {}

    for i in range(3):

        for j in range(3):

            s = i + j

            diagonal_values.setdefault(
                s,
                set()
            ).add(Q[i][j])

    print("DIAGONAL VALUE STRUCTURE")
    print()

    diagonal_exact = True

    for s in sorted(diagonal_values):

        values = diagonal_values[s]

        print(
            f"i+j={s}: "
            f"values={sorted(values)}"
        )

        if len(values) != 1:
            diagonal_exact = False

    print()

    print(
        "ONE VALUE PER DIAGONAL:",
        "PASS"
        if diagonal_exact
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Number of genuinely independent Q entries.
    # ------------------------------------------------------------

    independent_Q_values = [
        Q[0][0],
        Q[0][1],
        Q[0][2],
        Q[1][2],
        Q[2][2],
    ]

    distinct_Q_values = set(
        x
        for row in Q
        for x in row
    )

    print("Q INFORMATION COUNT")
    print()

    print(
        "total matrix entries =",
        9
    )

    print(
        "distinct Q values    =",
        len(distinct_Q_values)
    )

    print(
        "expected diagonal values =",
        5
    )

    print(
        "independent diagonal representatives =",
        len(independent_Q_values)
    )

    print()

    # ------------------------------------------------------------
    # Reconstruct all of Q from just the five diagonal values.
    # ------------------------------------------------------------

    Q_from_five = [
        [
            independent_Q_values[i + j]
            for j in range(3)
        ]
        for i in range(3)
    ]

    print_matrix(
        "Q RECONSTRUCTED FROM FIVE VALUES",
        Q_from_five
    )

    print(
        "Q RECONSTRUCTION FROM 5 SCALARS:",
        "PASS"
        if Q_from_five == Q
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Compare Q information with n alone.
    #
    # The five values themselves are simply:
    #
    # floor(n/B^2), ..., floor(n/B^6).
    # ------------------------------------------------------------

    five_from_n = [
        n // (BASE ** k)
        for k in range(2, 7)
    ]

    print("FIVE VALUES GENERATED DIRECTLY FROM n")
    print(five_from_n)
    print()

    print(
        "Q FIVE-VALUE VECTOR = n-DIVISION VECTOR:",
        "PASS"
        if five_from_n == independent_Q_values
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Factor-specific residual E.
    #
    # Q itself has collapsed to known information from n.
    #
    # Now inspect whether E contains additional structure that
    # survives the collapse.
    # ------------------------------------------------------------

    print("FACTOR-DEPENDENT RESIDUAL E")
    print()

    E_diagonal_values = {}

    for i in range(3):

        for j in range(3):

            s = i + j

            E_diagonal_values.setdefault(
                s,
                []
            ).append(E[i][j])

    for s in sorted(E_diagonal_values):

        print(
            f"E diagonal i+j={s}: "
            f"{E_diagonal_values[s]}"
        )

    print()

    # ------------------------------------------------------------
    # Test whether E itself is Hankel.
    # ------------------------------------------------------------

    E_hankel = True

    for i in range(3):

        for j in range(3):

            for u in range(3):

                for v in range(3):

                    if i + j == u + v:

                        if E[i][j] != E[u][v]:

                            E_hankel = False

    print(
        "E IS HANKEL:",
        "YES"
        if E_hankel
        else "NO"
    )
    print()

    # ------------------------------------------------------------
    # Compare nested Q with ordinary n information.
    #
    # If Q is completely reconstructible from n without p or q,
    # the nested quotient matrix is not giving a factor oracle.
    # ------------------------------------------------------------

    print("EXACT INFORMATION TEST")
    print()

    print(
        "Q depends on p,q individually:",
        "NO"
    )

    print(
        "Q reconstructible from n alone:",
        "YES"
        if hankel_exact
        else "NO"
    )

    print(
        "Q factor-specific residual information:",
        "NONE"
    )
    print()

    # ------------------------------------------------------------
    # Verify the direct formula for all Q entries.
    # ------------------------------------------------------------

    formula_checks = []

    for i in range(3):

        for j in range(3):

            direct = (
                n
                // (
                    R1[i] * R2[j]
                )
            )

            formula = (
                n
                // (
                    BASE ** (i + j + 2)
                )
            )

            formula_checks.append(
                direct == formula
            )

    print(
        "NESTED-RADIX FORMULA FOR ALL Q ENTRIES:",
        "PASS"
        if all(formula_checks)
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Runtime
    # ------------------------------------------------------------

    elapsed = time.perf_counter() - start

    print("-" * 72)
    print(
        f"{bits}-BIT CASE RUNTIME = "
        f"{elapsed:.6f} s"
    )
    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 143")
print("=" * 72)
