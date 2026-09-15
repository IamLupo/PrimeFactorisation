# ================================================================
# START EXPERIMENT 145
# CRT-accelerated exact carry-bound ratio lifting
#
# This is Experiment 144 with the expensive integer interval scan
# removed.
#
# Core equations:
#
#   43*k1 = 17*k0 - delta1
#   59*k2 = 17*k0 - delta2
#
# Therefore:
#
#   k0 = c_k (mod 43*59)
#
# for every fixed (delta1, delta2).
#
# Likewise:
#
#   ell0 = c_l (mod 37*61)
#
# for every fixed (epsilon1, epsilon2).
#
# We only enumerate members of these exact CRT residue classes
# inside the observable interval.
# ================================================================

import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

R1 = [17, 43, 59]
R2 = [19, 37, 61]

BIT_SIZES = [30, 36, 42, 48, 54]

BALANCE_RATIO = 4


# ================================================================
# Miller-Rabin
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
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

    for a in [2, 3, 5, 7, 11, 13, 17]:

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


def random_balanced_semiprime(bits):

    p_bits = bits // 2
    q_bits = bits - p_bits

    while True:

        p = random_prime(p_bits)
        q = random_prime(q_bits)

        if p == q:
            continue

        lo = min(p, q)
        hi = max(p, q)

        if hi <= BALANCE_RATIO * lo:
            return p, q


# ================================================================
# Exact quotient data
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
# True factor decomposition
# ================================================================

def build_true_KL(p, q):

    K = []
    A = []

    for r in R1:

        K.append(p // r)
        A.append(p % r)

    L = []
    B = []

    for s in R2:

        L.append(q // s)
        B.append(q % s)

    return K, A, L, B


def build_E(n, p, q):

    K, A, L, B = build_true_KL(p, q)

    Q = build_Q(n)

    E = []

    for i in range(3):

        row = []

        for j in range(3):

            row.append(
                Q[i][j] - K[i] * L[j]
            )

        E.append(row)

    return E


# ================================================================
# Carry bound
# ================================================================

def factor_upper_bound(n):

    # max(p,q) <= sqrt(BALANCE_RATIO*n)

    return math.isqrt(
        BALANCE_RATIO * n
    ) + 1


def carry_bound(n, i, j):

    max_factor = factor_upper_bound(n)

    kmax = max_factor // R1[i]
    lmax = max_factor // R2[j]

    return kmax + lmax


def build_K_intervals(n, Q):

    intervals = []

    for i in range(3):

        row = []

        for j in range(3):

            B = carry_bound(
                n,
                i,
                j
            )

            lo = max(
                1,
                Q[i][j] - B
            )

            hi = Q[i][j]

            row.append(
                (lo, hi)
            )

        intervals.append(row)

    return intervals


# ================================================================
# Exact fraction helpers
# ================================================================

def ratio_interval(num_interval, den_interval):

    nlo, nhi = num_interval
    dlo, dhi = den_interval

    return (
        Fraction(nlo, dhi),
        Fraction(nhi, dlo)
    )


# We use Fraction only after importing it below.
from fractions import Fraction


def row_ratio_intervals(KI):

    R = {}

    for i in (1, 2):

        for j in range(3):

            R[(i, j)] = ratio_interval(
                KI[i][j],
                KI[0][j]
            )

    return R


def column_ratio_intervals(KI):

    R = {}

    for j in (1, 2):

        for i in range(3):

            R[(j, i)] = ratio_interval(
                KI[i][j],
                KI[i][0]
            )

    return R


# ================================================================
# Exact interval for k0 from
#
#   rho = r0/r - delta/(r*k0)
# ================================================================

def k0_interval_from_delta(
    delta,
    r0,
    r,
    rho_interval
):

    if delta == 0:
        return None

    lo, hi = rho_interval

    d1 = (
        Fraction(r0, 1)
        - r * lo
    )

    d2 = (
        Fraction(r0, 1)
        - r * hi
    )

    # If zero occurs inside the denominator interval,
    # this particular observable does not bound k0.
    if d1 == 0 or d2 == 0:
        return None

    x1 = Fraction(delta, 1) / d1
    x2 = Fraction(delta, 1) / d2

    return (
        min(x1, x2),
        max(x1, x2)
    )


def intersect(intervals):

    if not intervals:
        return None

    lo = intervals[0][0]
    hi = intervals[0][1]

    for a, b in intervals[1:]:

        lo = max(lo, a)
        hi = min(hi, b)

        if lo > hi:
            return None

    return lo, hi


def integer_interval(interval):

    if interval is None:
        return None

    lo, hi = interval

    lo_int = (
        (lo.numerator + lo.denominator - 1)
        // lo.denominator
    )

    hi_int = (
        hi.numerator
        // hi.denominator
    )

    if lo_int > hi_int:
        return None

    return lo_int, hi_int


# ================================================================
# CRT
# ================================================================

def crt_pair(a1, m1, a2, m2):

    """
    Solve

        x = a1 mod m1
        x = a2 mod m2

    assuming gcd(m1,m2)=1.
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

    modulus = m1 * m2

    return x % modulus, modulus


def linear_congruence_residue(
    coefficient,
    value,
    modulus
):

    """
    Solve

        coefficient*x = value (mod modulus)

    when coefficient is invertible modulo modulus.
    """

    inv = pow(
        coefficient,
        -1,
        modulus
    )

    return (
        value * inv
    ) % modulus


# ================================================================
# Generate k0 CRT class
# ================================================================

def k0_crt_class(delta1, delta2):

    a1 = linear_congruence_residue(
        R1[0],
        delta1,
        R1[1]
    )

    a2 = linear_congruence_residue(
        R1[0],
        delta2,
        R1[2]
    )

    return crt_pair(
        a1,
        R1[1],
        a2,
        R1[2]
    )


# ================================================================
# Generate ell0 CRT class
# ================================================================

def ell0_crt_class(epsilon1, epsilon2):

    a1 = linear_congruence_residue(
        R2[0],
        epsilon1,
        R2[1]
    )

    a2 = linear_congruence_residue(
        R2[0],
        epsilon2,
        R2[2]
    )

    return crt_pair(
        a1,
        R2[1],
        a2,
        R2[2]
    )


# ================================================================
# Enumerate a residue class inside an integer interval
#
# This replaces the O(interval-width) scan.
# ================================================================

def enumerate_residue_class(
    residue,
    modulus,
    lo,
    hi
):

    if lo > hi:
        return []

    # First x >= lo satisfying x = residue (mod modulus).

    if residue < lo:

        k = (
            (lo - residue + modulus - 1)
            // modulus
        )

        x = residue + k * modulus

    else:

        x = residue

    out = []

    while x <= hi:

        out.append(x)

        x += modulus

    return out


# ================================================================
# Build exact K candidates
# ================================================================

def exact_k_states(
    n,
    Q,
    ratio_intervals
):

    delta1_range = range(
        -(R1[0] - 1),
        R1[1]
    )

    delta2_range = range(
        -(R1[0] - 1),
        R1[2]
    )

    states = []

    interval_states = 0

    crt_points = 0

    M = R1[1] * R1[2]

    for d1 in delta1_range:

        # Derive k0 interval using row 1.

        parts1 = []

        for j in range(3):

            part = k0_interval_from_delta(
                d1,
                R1[0],
                R1[1],
                ratio_intervals[(1, j)]
            )

            if part is not None:
                parts1.append(part)

        interval1 = intersect(parts1)

        if interval1 is None:
            continue

        int1 = integer_interval(interval1)

        if int1 is None:
            continue

        for d2 in delta2_range:

            parts2 = []

            for j in range(3):

                part = k0_interval_from_delta(
                    d2,
                    R1[0],
                    R1[2],
                    ratio_intervals[(2, j)]
                )

                if part is not None:
                    parts2.append(part)

            interval2 = intersect(parts2)

            if interval2 is None:
                continue

            combined = intersect(
                [
                    interval1,
                    interval2,
                ]
            )

            if combined is None:
                continue

            integer_bounds = integer_interval(
                combined
            )

            if integer_bounds is None:
                continue

            interval_states += 1

            lo, hi = integer_bounds

            residue, modulus = k0_crt_class(
                d1,
                d2
            )

            points = enumerate_residue_class(
                residue,
                modulus,
                lo,
                hi
            )

            crt_points += len(points)

            for k0 in points:

                # Recover k1,k2 exactly.

                x1 = (
                    R1[0] * k0 - d1
                )

                x2 = (
                    R1[0] * k0 - d2
                )

                if x1 % R1[1] != 0:
                    continue

                if x2 % R1[2] != 0:
                    continue

                k1 = x1 // R1[1]
                k2 = x2 // R1[2]

                if k0 <= 0 or k1 <= 0 or k2 <= 0:
                    continue

                states.append(
                    (
                        k0,
                        k1,
                        k2,
                        d1,
                        d2,
                    )
                )

    return (
        sorted(set(states)),
        interval_states,
        crt_points,
        M
    )


# ================================================================
# Build exact ell candidates
# ================================================================

def exact_ell_states(
    n,
    Q,
    ratio_intervals
):

    epsilon1_range = range(
        -(R2[0] - 1),
        R2[1]
    )

    epsilon2_range = range(
        -(R2[0] - 1),
        R2[2]
    )

    states = []

    interval_states = 0
    crt_points = 0

    M = R2[1] * R2[2]

    for e1 in epsilon1_range:

        parts1 = []

        for i in range(3):

            part = k0_interval_from_delta(
                e1,
                R2[0],
                R2[1],
                ratio_intervals[(1, i)]
            )

            if part is not None:
                parts1.append(part)

        interval1 = intersect(parts1)

        if interval1 is None:
            continue

        for e2 in epsilon2_range:

            parts2 = []

            for i in range(3):

                part = k0_interval_from_delta(
                    e2,
                    R2[0],
                    R2[2],
                    ratio_intervals[(2, i)]
                )

                if part is not None:
                    parts2.append(part)

            interval2 = intersect(parts2)

            if interval2 is None:
                continue

            combined = intersect(
                [
                    interval1,
                    interval2,
                ]
            )

            if combined is None:
                continue

            integer_bounds = integer_interval(
                combined
            )

            if integer_bounds is None:
                continue

            interval_states += 1

            lo, hi = integer_bounds

            residue, modulus = \
                ell0_crt_class(
                    e1,
                    e2
                )

            points = enumerate_residue_class(
                residue,
                modulus,
                lo,
                hi
            )

            crt_points += len(points)

            for ell0 in points:

                x1 = (
                    R2[0] * ell0 - e1
                )

                x2 = (
                    R2[0] * ell0 - e2
                )

                if x1 % R2[1] != 0:
                    continue

                if x2 % R2[2] != 0:
                    continue

                ell1 = x1 // R2[1]
                ell2 = x2 // R2[2]

                if ell0 <= 0 or ell1 <= 0 or ell2 <= 0:
                    continue

                states.append(
                    (
                        ell0,
                        ell1,
                        ell2,
                        e1,
                        e2,
                    )
                )

    return (
        sorted(set(states)),
        interval_states,
        crt_points,
        M
    )


# ================================================================
# Expand quotient state into possible p
# ================================================================

def expand_k_state(state):

    k0, k1, k2, d1, d2 = state

    candidates = []

    for a0 in range(R1[0]):

        a1 = a0 + d1
        a2 = a0 + d2

        if not (0 <= a1 < R1[1]):
            continue

        if not (0 <= a2 < R1[2]):
            continue

        p = R1[0] * k0 + a0

        if p // R1[1] != k1:
            continue

        if p // R1[2] != k2:
            continue

        candidates.append(p)

    return candidates


def expand_ell_state(state):

    ell0, ell1, ell2, e1, e2 = state

    candidates = []

    for b0 in range(R2[0]):

        b1 = b0 + e1
        b2 = b0 + e2

        if not (0 <= b1 < R2[1]):
            continue

        if not (0 <= b2 < R2[2]):
            continue

        q = R2[0] * ell0 + b0

        if q // R2[1] != ell1:
            continue

        if q // R2[2] != ell2:
            continue

        candidates.append(q)

    return candidates


# ================================================================
# Exact factor verification
# ================================================================

def verify_factor(n, p):

    if p <= 1:
        return None

    if n % p != 0:
        return None

    q = n // p

    if q <= 1:
        return None

    return p, q


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 145")
print("CRT-accelerated exact carry-bound ratio lifting")
print("=" * 72)
print()

print("R1 =", R1)
print("R2 =", R2)
print("BALANCE RATIO =", BALANCE_RATIO)
print()

K_CRT_MODULUS = R1[1] * R1[2]
L_CRT_MODULUS = R2[1] * R2[2]

print(
    "K CRT MODULUS =",
    K_CRT_MODULUS
)

print(
    "ELL CRT MODULUS =",
    L_CRT_MODULUS
)

print()


for bits in BIT_SIZES:

    case_start = time.perf_counter()

    print("=" * 72)
    print(
        f"GENERATING {bits}-BIT BALANCED SEMIPRIME"
    )
    print("=" * 72)
    print()

    p, q = random_balanced_semiprime(bits)
    n = p * q

    print("-" * 72)
    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    print(
        "factor ratio =",
        f"{max(p,q) / min(p,q):.6f}"
    )

    print()

    K_true, A_true, L_true, B_true = \
        build_true_KL(
            p,
            q
        )

    E_true = build_E(
        n,
        p,
        q
    )

    delta1_true = (
        R1[0] * K_true[0]
        - R1[1] * K_true[1]
    )

    delta2_true = (
        R1[0] * K_true[0]
        - R1[2] * K_true[2]
    )

    epsilon1_true = (
        R2[0] * L_true[0]
        - R2[1] * L_true[1]
    )

    epsilon2_true = (
        R2[0] * L_true[0]
        - R2[2] * L_true[2]
    )

    print("TRUE K =", K_true)
    print("TRUE A =", A_true)
    print("TRUE L =", L_true)
    print("TRUE B =", B_true)

    print()

    print(
        "TRUE DELTAS =",
        delta1_true,
        delta2_true
    )

    print(
        "TRUE EPSILONS =",
        epsilon1_true,
        epsilon2_true
    )

    print()

    # ------------------------------------------------------------
    # Observable Q
    # ------------------------------------------------------------

    Q = build_Q(n)

    # ------------------------------------------------------------
    # Observable carry-derived intervals
    # ------------------------------------------------------------

    KI = build_K_intervals(
        n,
        Q
    )

    row_ratios = row_ratio_intervals(
        KI
    )

    column_ratios = column_ratio_intervals(
        KI
    )

    # ------------------------------------------------------------
    # K reconstruction
    # ------------------------------------------------------------

    t_k = time.perf_counter()

    (
        k_states,
        k_interval_states,
        k_crt_points,
        k_modulus
    ) = exact_k_states(
        n,
        Q,
        row_ratios
    )

    k_elapsed = time.perf_counter() - t_k

    # ------------------------------------------------------------
    # L reconstruction
    # ------------------------------------------------------------

    t_l = time.perf_counter()

    (
        ell_states,
        ell_interval_states,
        ell_crt_points,
        ell_modulus
    ) = exact_ell_states(
        n,
        Q,
        column_ratios
    )

    ell_elapsed = time.perf_counter() - t_l

    print("CRT-LIFT K RESULTS")
    print()

    print(
        "interval residue-pairs =",
        k_interval_states
    )

    print(
        "actual CRT k0 points tested =",
        k_crt_points
    )

    print(
        "K states =",
        len(k_states)
    )

    print(
        "K modulus =",
        k_modulus
    )

    print(
        "K lifting runtime =",
        f"{k_elapsed:.6f} s"
    )

    print()

    print("CRT-LIFT ELL RESULTS")
    print()

    print(
        "interval residue-pairs =",
        ell_interval_states
    )

    print(
        "actual CRT ell0 points tested =",
        ell_crt_points
    )

    print(
        "ELL states =",
        len(ell_states)
    )

    print(
        "ELL modulus =",
        ell_modulus
    )

    print(
        "ELL lifting runtime =",
        f"{ell_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # True-state membership
    # ------------------------------------------------------------

    true_K_present = False

    for state in k_states:

        if (
            state[0] == K_true[0]
            and
            state[1] == K_true[1]
            and
            state[2] == K_true[2]
        ):
            true_K_present = True
            break

    true_L_present = False

    for state in ell_states:

        if (
            state[0] == L_true[0]
            and
            state[1] == L_true[1]
            and
            state[2] == L_true[2]
        ):
            true_L_present = True
            break

    print(
        "TRUE K VECTOR PRESENT =",
        true_K_present
    )

    print(
        "TRUE ELL VECTOR PRESENT =",
        true_L_present
    )

    print()

    # ------------------------------------------------------------
    # Expand to p and q.
    # ------------------------------------------------------------

    p_candidates = set()

    for state in k_states:

        for pc in expand_k_state(state):

            p_candidates.add(pc)

    q_candidates = set()

    for state in ell_states:

        for qc in expand_ell_state(state):

            q_candidates.add(qc)

    print("FACTOR CANDIDATES")
    print()

    print(
        "p candidates =",
        len(p_candidates)
    )

    print(
        "q candidates =",
        len(q_candidates)
    )

    print()

    # ------------------------------------------------------------
    # Exact p scan only over generated candidates.
    # ------------------------------------------------------------

    verified = set()

    for pc in p_candidates:

        result = verify_factor(
            n,
            pc
        )

        if result is None:
            continue

        p_found, q_found = result

        if q_found in q_candidates:

            verified.add(
                tuple(
                    sorted(
                        (
                            p_found,
                            q_found
                        )
                    )
                )
            )

    print("EXACT VERIFICATION")
    print()

    print(
        "verified factor pairs =",
        len(verified)
    )

    for pair in sorted(verified):

        print(
            "  ",
            pair
        )

    print()

    # ------------------------------------------------------------
    # True CRT classes.
    # ------------------------------------------------------------

    true_k_residue, _ = k0_crt_class(
        delta1_true,
        delta2_true
    )

    true_l_residue, _ = ell0_crt_class(
        epsilon1_true,
        epsilon2_true
    )

    print("TRUE CRT CLASSES")
    print()

    print(
        "true k0 =",
        K_true[0]
    )

    print(
        "true k0 residue =",
        true_k_residue,
        "mod",
        K_CRT_MODULUS
    )

    print(
        "true ell0 =",
        L_true[0]
    )

    print(
        "true ell0 residue =",
        true_l_residue,
        "mod",
        L_CRT_MODULUS
    )

    print()

    # ------------------------------------------------------------
    # Candidate density
    # ------------------------------------------------------------

    print("CANDIDATE DENSITY")
    print()

    if k_interval_states:

        print(
            "K CRT points / interval states =",
            f"{k_crt_points / k_interval_states:.6f}"
        )

    if ell_interval_states:

        print(
            "ELL CRT points / interval states =",
            f"{ell_crt_points / ell_interval_states:.6f}"
        )

    print()

    sqrt_n = math.isqrt(n)

    print("COMPARISON")
    print()

    print(
        "sqrt(n) =",
        sqrt_n
    )

    print(
        "K states =",
        len(k_states)
    )

    print(
        "ELL states =",
        len(ell_states)
    )

    print(
        "p candidates =",
        len(p_candidates)
    )

    print(
        "q candidates =",
        len(q_candidates)
    )

    print(
        "verified =",
        len(verified)
    )

    print()

    elapsed = time.perf_counter() - case_start

    print(
        f"TOTAL RUNTIME = {elapsed:.6f} s"
    )

    print()
    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 145")
print("=" * 72)
