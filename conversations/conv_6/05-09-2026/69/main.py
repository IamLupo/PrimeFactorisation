# ================================================================
# START EXPERIMENT 152
# Exact carry-floor residue constraint system
#
# Goal:
#
#   Attack the exact formula
#
#     E_ij =
#       floor(
#           k_i*b_j/s_j
#         + ell_j*a_i/r_i
#         + a_i*b_j/(r_i*s_j)
#       )
#
#   without scanning k_i or ell_j up to sqrt(n).
#
# For a candidate residue vector
#
#     a_i = p mod r_i
#     b_j = q mod s_j
#
# the corresponding factors must satisfy
#
#     p = r_i*k_i + a_i
#     q = s_j*ell_j + b_j.
#
# Since
#
#     Q_ij = k_i*ell_j + E_ij
#
# we obtain an integer interval for k_i*ell_j.
#
# The experiment searches only the bounded residue space and
# checks whether the implied quotient products can be mutually
# consistent with a rank-1 K matrix.
#
# This is a residue-space experiment, NOT a sqrt(n) search.
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


# ------------------------------------------------
# Safety limits
# ------------------------------------------------

MAX_RESIDUE_STATES = 2_000_000


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
# Exact factor decomposition
# ================================================================

def factor_state(p, q):

    a = [
        p % r
        for r in R1
    ]

    k = [
        p // r
        for r in R1
    ]

    b = [
        q % s
        for s in R2
    ]

    ell = [
        q // s
        for s in R2
    ]

    return k, a, ell, b


# ================================================================
# Quotient matrix
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
# Exact carry formula
# ================================================================

def carry_value(
    k,
    a,
    ell,
    b,
    r,
    s,
):
    """
    Exact:

      E = floor(
            k*b/s
          + ell*a/r
          + a*b/(r*s)
          )
    """

    # Compute the whole quantity exactly using one numerator.

    numerator = (
        k * b * r
        + ell * a * s
        + a * b
    )

    denominator = (
        r * s
    )

    return numerator // denominator


# ================================================================
# Derive exact K intervals from residues
# ================================================================

def carry_bounds_from_residues(
    n,
    Q,
    a0,
    ai,
    b0,
    bj,
    i,
    j,
):
    """
    We don't know k_i and ell_j individually.

    We use the factor relation:

        p = r_i*k_i + a_i
        q = s_j*ell_j + b_j

    together with

        Q_ij = k_i*ell_j + E_ij.

    For fixed residues, E_ij is an affine-floor function of
    k_i and ell_j.

    This function returns a conservative exact interval for
    k_i*ell_j based only on the factor bounds induced by n.
    """

    r = R1[i]
    s = R2[j]

    # Crude but exact-safe balanced upper bound.
    root = math.isqrt(
        4 * n
    )

    while root * root < 4 * n:
        root += 1

    k_max = root // r + 1
    ell_max = root // s + 1

    # Since
    #
    # E = floor(
    #      kb/s + ell*a/r + ab/(rs)
    #     )
    #
    # and
    #
    # 0 <= a < r
    # 0 <= b < s,
    #
    # we have
    #
    # 0 <= E <= k + ell.
    #
    E_min = 0
    E_max = k_max + ell_max

    return (
        max(0, Q[i][j] - E_max),
        Q[i][j] - E_min
    )


# ================================================================
# Recover possible k0 intervals from residue differences
# ================================================================

def possible_k0_values_from_residue_vector(
    n,
    Q,
    a,
):
    """
    Given a candidate residue vector

        a = [a0,a1,a2,a3],

    we know

        p = r_i*k_i + a_i

    for every i.

    Therefore

        r_i*k_i + a_i
        =
        r_0*k_0 + a_0.

    Thus

        k_i =
          (r0*k0 + a0 - ai)/ri.

    Rather than scanning k0, derive the congruence conditions:

        r0*k0 = ai-a0 (mod ri).

    CRT gives a single k0 residue class.

    Return:
        residue,
        modulus.
    """

    residue = 0
    modulus = 1

    for i in range(1, len(R1)):

        r = R1[i]

        delta = (
            a[i] - a[0]
        )

        # r0*k0 = delta (mod r)
        #
        # r0 is invertible modulo r for our selected radices.

        target = (
            delta
            * pow(
                R1[0],
                -1,
                r
            )
        ) % r

        residue, modulus = crt_pair(
            residue,
            modulus,
            target,
            r
        )

    return residue, modulus


# ================================================================
# CRT
# ================================================================

def crt_pair(
    a1,
    m1,
    a2,
    m2,
):
    """
    Solve:

      x = a1 mod m1
      x = a2 mod m2
    """

    inv = pow(
        m1,
        -1,
        m2
    )

    t = (
        (a2 - a1)
        * inv
    ) % m2

    x = a1 + m1 * t

    return (
        x % (m1 * m2),
        m1 * m2,
    )


# ================================================================
# Same for ell0
# ================================================================

def possible_ell0_class(
    b,
):
    """
    Since

        s_i*ell_i + b_i
        =
        s_0*ell_0 + b_0,

    obtain a CRT class for ell0.
    """

    residue = 0
    modulus = 1

    for j in range(1, len(R2)):

        s = R2[j]

        delta = (
            b[j] - b[0]
        )

        target = (
            delta
            * pow(
                R2[0],
                -1,
                s
            )
        ) % s

        residue, modulus = crt_pair(
            residue,
            modulus,
            target,
            s
        )

    return residue, modulus


# ================================================================
# Reconstruct k-vector from k0
# ================================================================

def reconstruct_k(
    k0,
    a,
):
    """
    Return k_i implied by

        r_i*k_i + a_i
        =
        r0*k0 + a0.
    """

    p = (
        R1[0] * k0
        + a[0]
    )

    k = []

    for i, r in enumerate(R1):

        if p < a[i]:
            return None

        numerator = (
            p - a[i]
        )

        if numerator % r != 0:
            return None

        ki = numerator // r

        if ki <= 0:
            return None

        k.append(ki)

    return k


# ================================================================
# Reconstruct ell-vector
# ================================================================

def reconstruct_ell(
    ell0,
    b,
):
    """
    Return ell_i implied by

        s_i*ell_i + b_i
        =
        s0*ell0 + b0.
    """

    q = (
        R2[0] * ell0
        + b[0]
    )

    ell = []

    for j, s in enumerate(R2):

        numerator = (
            q - b[j]
        )

        if numerator < 0:
            return None

        if numerator % s != 0:
            return None

        lj = numerator // s

        if lj <= 0:
            return None

        ell.append(lj)

    return ell


# ================================================================
# Exact residue-state verification
# ================================================================

def verify_residue_state(
    n,
    Q,
    a,
    b,
):
    """
    Given complete residue vectors, derive the exact p/q residue
    CRT classes and test whether any compatible positive factors
    exist within the balanced search region.

    IMPORTANT:
      This routine does NOT scan every k0.

    It only tests the unique CRT residue class for k0 and ell0
    against the observable Q product intervals.
    """

    root4 = math.isqrt(
        4 * n
    )

    p_max = root4 + 1
    q_max = root4 + 1

    k0_res, k0_mod = \
        possible_k0_values_from_residue_vector(
            n,
            Q,
            a
        )

    ell0_res, ell0_mod = \
        possible_ell0_class(
            b
        )

    # Compute observable k0/ell0 ranges.
    k0_min = 1
    k0_max = (
        p_max
        // R1[0]
        + 1
    )

    ell0_min = 1
    ell0_max = (
        q_max
        // R2[0]
        + 1
    )

    # Jump to the first CRT-compatible values.
    k0_values = enumerate_class(
        k0_res,
        k0_mod,
        k0_min,
        k0_max
    )

    ell0_values = enumerate_class(
        ell0_res,
        ell0_mod,
        ell0_min,
        ell0_max
    )

    # We don't want this routine to become another sqrt scan.
    #
    # A residue vector should normally give sparse CRT classes.
    #
    # Only use the resulting finite CRT points.

    for k0 in k0_values:

        k = reconstruct_k(
            k0,
            a
        )

        if k is None:
            continue

        p = (
            R1[0] * k0
            + a[0]
        )

        if p <= 1 or p > p_max:
            continue

        for ell0 in ell0_values:

            ell = reconstruct_ell(
                ell0,
                b
            )

            if ell is None:
                continue

            q = (
                R2[0] * ell0
                + b[0]
            )

            if q <= 1 or q > q_max:
                continue

            # Check all nine exact quotient equations.
            exact = True

            for i, r in enumerate(R1):

                for j, s in enumerate(R2):

                    E = carry_value(
                        k[i],
                        a[i],
                        ell[j],
                        b[j],
                        r,
                        s,
                    )

                    expected = (
                        k[i] * ell[j]
                        + E
                    )

                    if expected != Q[i][j]:

                        exact = False
                        break

                if not exact:
                    break

            if not exact:
                continue

            if p * q == n:

                return (
                    p,
                    q,
                    k,
                    ell,
                )

    return None


# ================================================================
# Enumerate CRT class
# ================================================================

def enumerate_class(
    residue,
    modulus,
    lo,
    hi,
):
    """
    Enumerate x = residue mod modulus in [lo,hi].
    """

    if lo > hi:
        return []

    if residue < lo:

        x = residue + (
            (
                lo
                - residue
                + modulus
                - 1
            )
            // modulus
        ) * modulus

    else:

        x = residue

    result = []

    while x <= hi:

        result.append(x)

        x += modulus

        # Safety.
        if len(result) > MAX_RESIDUE_STATES:
            break

    return result


# ================================================================
# Exhaustive residue-vector experiment
# ================================================================

def run_residue_experiment(
    n,
    p,
    q,
    Q,
):
    """
    Enumerate a_i and b_j vectors.

    This measures whether exact carry equations reduce the
    residue-space before any large quotient search.

    WARNING:
      The complete Cartesian space is large:
        17*43*59*71
      times
        19*37*61*73.

    Therefore we exploit the first residue relation to build
    candidate p residues modulo the R1 product and q residues
    modulo the R2 product, then couple them through n.
    """

    M_p = math.prod(R1)
    M_q = math.prod(R2)

    print(
        "M_p =",
        M_p
    )

    print(
        "M_q =",
        M_q
    )

    # ------------------------------------------------------------
    # Enumerate p residue classes using CRT directly.
    # ------------------------------------------------------------

    p_residue_classes = []

    for a0 in range(R1[0]):

        for a1 in range(R1[1]):

            for a2 in range(R1[2]):

                for a3 in range(R1[3]):

                    a = [
                        a0,
                        a1,
                        a2,
                        a3,
                    ]

                    residue, modulus = \
                        possible_k0_values_from_residue_vector(
                            n,
                            Q,
                            a
                        )

                    # Corresponding p residue modulo M_p.
                    p_residue = (
                        R1[0] * residue
                        + a0
                    ) % M_p

                    p_residue_classes.append(
                        (
                            p_residue,
                            tuple(a),
                        )
                    )

    print(
        "raw p residue vectors =",
        len(p_residue_classes)
    )

    # Deduplicate.
    p_residue_classes = list(
        dict.fromkeys(
            p_residue_classes
        )
    )

    print(
        "unique p residue classes =",
        len(p_residue_classes)
    )

    print()

    # ------------------------------------------------------------
    # Same for q.
    # ------------------------------------------------------------

    q_residue_classes = []

    for b0 in range(R2[0]):

        for b1 in range(R2[1]):

            for b2 in range(R2[2]):

                for b3 in range(R2[3]):

                    b = [
                        b0,
                        b1,
                        b2,
                        b3,
                    ]

                    residue, modulus = \
                        possible_ell0_class(
                            b
                        )

                    q_residue = (
                        R2[0] * residue
                        + b0
                    ) % M_q

                    q_residue_classes.append(
                        (
                            q_residue,
                            tuple(b),
                        )
                    )

    q_residue_classes = list(
        dict.fromkeys(
            q_residue_classes
        )
    )

    print(
        "unique q residue classes =",
        len(q_residue_classes)
    )

    print()

    # ------------------------------------------------------------
    # Important:
    #
    # Do NOT Cartesian-product the two giant lists.
    #
    # Instead derive q residue from each p residue:
    #
    #     q = n * p^{-1} mod gcd/product.
    #
    # Since M_p and M_q are different, use their gcd and exact
    # shared radices only as a preliminary coupling.
    #
    # For this experiment we use the common radix subset:
    #
    #     gcd(M_p,M_q) = 1
    #
    # so there is no direct common-modulus coupling.
    #
    # Instead use the carry equations as the coupling filter.
    # ------------------------------------------------------------

    candidate_count = 0
    true_a_present = False
    true_b_present = False

    true_a = tuple(
        p % r
        for r in R1
    )

    true_b = tuple(
        q % s
        for s in R2
    )

    # ------------------------------------------------------------
    # We use only the first two columns and rows for a cheap
    # exact carry filter.
    #
    # For each candidate residue vector, estimate whether the
    # corresponding K rank-1 products can occupy Q intervals.
    # ------------------------------------------------------------

    # Build a map from p residue to vectors.
    p_map = {}

    for pres, avec in p_residue_classes:

        p_map.setdefault(
            pres,
            []
        ).append(
            avec
        )

    q_map = {}

    for qres, bvec in q_residue_classes:

        q_map.setdefault(
            qres,
            []
        ).append(
            bvec
        )

    # The direct product would still be huge. Instead, use
    # the fact that n = p*q modulo R1[i]*R2[j].
    #
    # Test residue pairs via the nine small moduli.
    #
    # We generate only q residue candidates forced by p modulo
    # each small product and intersect them.
    #
    # Start from each p vector but stop if candidate growth gets
    # excessive.

    filtered = []

    max_pairs = 2_000_000

    for pres, avecs in p_map.items():

        # Determine q modulo each R1[i]*R2[j] using each possible
        # p residue modulo R1[i].
        #
        # Since q is unknown modulo the R2 product, recover the
        # required residue modulo each s_j directly:
        #
        #     q = n / p (mod s_j)
        #
        # where p residue modulo s_j is not known.
        #
        # Therefore this route itself doesn't close.
        #
        # We instead record the fact explicitly and avoid pretending
        # it is a solved coupling.
        #
        # For the experiment we only count residue states.

        if true_a in avecs:

            true_a_present = True

        candidate_count += len(avecs)

        if candidate_count > max_pairs:
            break

    for qres, bvecs in q_map.items():

        if true_b in bvecs:
            true_b_present = True

    return (
        len(p_residue_classes),
        len(q_residue_classes),
        candidate_count,
        true_a_present,
        true_b_present,
    )


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 152")
print("Exact carry-floor residue constraint system")
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

    K, A, L, B = factor_state(
        p,
        q
    )

    Q = build_Q(
        n
    )

    print(
        "TRUE a =",
        A
    )

    print(
        "TRUE b =",
        B
    )

    print()

    # ------------------------------------------------------------
    # Verify exact carry equations.
    # ------------------------------------------------------------

    E = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            e = carry_value(
                K[i],
                A[i],
                L[j],
                B[j],
                r,
                s,
            )

            row.append(e)

        E.append(row)

    print(
        "TRUE E MATRIX"
    )

    for row in E:

        print(
            "["
            + ", ".join(str(x) for x in row)
            + "]"
        )

    print()

    exact_ok = True

    for i in range(4):

        for j in range(4):

            if (
                K[i] * L[j]
                + E[i][j]
                != Q[i][j]
            ):

                exact_ok = False

    print(
        "Q = K + E EXACTLY:",
        "PASS"
        if exact_ok
        else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # Residue state enumeration.
    # ------------------------------------------------------------

    enum_start = time.perf_counter()

    (
        p_count,
        q_count,
        partial_count,
        true_a_present,
        true_b_present,
    ) = run_residue_experiment(
        n,
        p,
        q,
        Q
    )

    enum_elapsed = (
        time.perf_counter()
        - enum_start
    )

    print(
        "RESIDUE ENUMERATION"
    )

    print()

    print(
        "p residue states =",
        p_count
    )

    print(
        "q residue states =",
        q_count
    )

    print(
        "partial states inspected =",
        partial_count
    )

    print(
        "true A present =",
        true_a_present
    )

    print(
        "true B present =",
        true_b_present
    )

    print(
        "enumeration runtime =",
        f"{enum_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # Show exact true residue CRT classes.
    # ------------------------------------------------------------

    p_modulus = math.prod(R1)
    q_modulus = math.prod(R2)

    true_p_residue = p % p_modulus
    true_q_residue = q % q_modulus

    print(
        "TRUE RESIDUE CRT CLASSES"
    )

    print()

    print(
        "p mod product(R1) =",
        true_p_residue
    )

    print(
        "q mod product(R2) =",
        true_q_residue
    )

    print(
        "p CRT modulus =",
        p_modulus
    )

    print(
        "q CRT modulus =",
        q_modulus
    )

    print()

    # ------------------------------------------------------------
    # Critical question:
    #
    # Does residue enumeration itself reduce the factor problem?
    #
    # Count of possible p residues is independent of n except
    # for coprimality.
    # ------------------------------------------------------------

    sqrt_n = math.isqrt(n)

    print(
        "COMPLEXITY COMPARISON"
    )

    print()

    print(
        "sqrt(n) =",
        sqrt_n
    )

    print(
        "R1 product =",
        p_modulus
    )

    print(
        "R2 product =",
        q_modulus
    )

    print()

    if p_modulus > sqrt_n:

        print(
            "R1 product > sqrt(n): "
            "one p residue class identifies at most one p "
            "below sqrt(n)."
        )

    else:

        print(
            "R1 product <= sqrt(n): "
            "multiple p values share a residue class."
        )

    print()

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
print("FINISHED EXPERIMENT 152")
print("=" * 72)
