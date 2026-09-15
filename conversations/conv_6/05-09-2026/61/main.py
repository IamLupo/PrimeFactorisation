# ================================================================
# START EXPERIMENT 144
# Exact carry-bound ratio lifting
#
# Goal:
#   Recover k0 and ell0 from observable Q intervals only,
#   without scanning to sqrt(n).
#
# Core facts:
#
#   Q_ij = k_i * ell_j + E_ij
#
#   0 <= E_ij <= k_i + ell_j
#
# For a balanced semiprime p,q with q/p <= BALANCE_RATIO:
#
#   p,q <= sqrt(n * BALANCE_RATIO)
#
# Therefore a completely n-derived upper bound for E_ij exists.
#
# Radix relation:
#
#   r0*k0 - ri*ki = delta_i
#
# where
#
#   delta_i = ai - a0
#
# and
#
#   delta_i in [-(r0-1), ri-1].
#
# Likewise:
#
#   s0*ell0 - sj*ellj = epsilon_j.
#
# ================================================================

from fractions import Fraction
import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

R1 = [17, 43, 59]
R2 = [19, 37, 61]

BIT_SIZES = [30, 36, 42, 48, 54]

# We deliberately test balanced semiprimes.
#
# q/p <= 4  implies
#
# max(p,q) <= 2*sqrt(n).
#
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
# Matrix helpers
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


def build_true_K(p, q):

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

    return K, A, L, B


def build_true_E(n, p, q):

    K, A, L, B = build_true_K(p, q)

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
# Carry bound
# ================================================================

def factor_upper_bound_from_balance(n):

    # q/p <= BALANCE_RATIO.
    #
    # If pq=n and q/p <= B then
    #
    # q <= sqrt(B*n)
    #
    # and similarly for p.
    #
    # For B=4:
    #
    # max(p,q) <= 2*sqrt(n).

    x = math.isqrt(
        BALANCE_RATIO * n
    )

    # Make the bound safely inclusive.
    while x * x < BALANCE_RATIO * n:
        x += 1

    return x


def carry_upper_bound(n, i, j):

    max_factor = factor_upper_bound_from_balance(n)

    max_k = max_factor // R1[i]
    max_ell = max_factor // R2[j]

    # Exact universal inequality:
    #
    # E_ij <= k_i + ell_j
    #
    return max_k + max_ell


def build_K_intervals(n, Q):

    """
    Since

        Q_ij = K_ij + E_ij

    and

        0 <= E_ij <= B_ij,

    we obtain

        Q_ij - B_ij <= K_ij <= Q_ij.
    """

    intervals = []

    for i in range(3):

        row = []

        for j in range(3):

            B = carry_upper_bound(
                n,
                i,
                j
            )

            lower = max(
                1,
                Q[i][j] - B
            )

            upper = Q[i][j]

            row.append(
                (lower, upper)
            )

        intervals.append(row)

    return intervals


# ================================================================
# Ratio intervals
# ================================================================

def quotient_interval(num_interval, den_interval):
    """
    Positive intervals:

        numerator in [Nlo, Nhi]
        denominator in [Dlo, Dhi]

    therefore

        num/den in
        [Nlo/Dhi, Nhi/Dlo].
    """

    Nlo, Nhi = num_interval
    Dlo, Dhi = den_interval

    return (
        Fraction(Nlo, Dhi),
        Fraction(Nhi, Dlo)
    )


def row_ratio_intervals(K_intervals):

    """
    For each i>0 and column j:

        ki/k0 = Kij / K0j.

    Produce exact ratio intervals.
    """

    ratios = {}

    for i in (1, 2):

        for j in range(3):

            ratios[(i, j)] = quotient_interval(
                K_intervals[i][j],
                K_intervals[0][j]
            )

    return ratios


def column_ratio_intervals(K_intervals):

    """
    For each j>0 and row i:

        ell_j/ell0 = Kij / Ki0.
    """

    ratios = {}

    for j in (1, 2):

        for i in range(3):

            ratios[(j, i)] = quotient_interval(
                K_intervals[i][j],
                K_intervals[i][0]
            )

    return ratios


# ================================================================
# Convert ratio interval -> k0 interval
# ================================================================

def delta_to_base_interval(delta, radix0, radi, ratio_interval):

    """
    Exact relation:

        radix0*k0 - radi*ki = delta

    hence

        ki/k0
          = radix0/radi
            - delta/(radi*k0)

    Therefore

        k0
          = delta /
            (radix0 - radi*rho)

    where rho = ki/k0.

    rho is known only inside an interval.

    Return the exact resulting k0 interval.
    """

    lo, hi = ratio_interval

    if delta == 0:

        # Then

        #   ki/k0 = radix0/radi

        # and this equation alone contains no magnitude
        # information about k0.

        return None

    d_lo = (
        Fraction(radix0, 1)
        - radi * lo
    )

    d_hi = (
        Fraction(radix0, 1)
        - radi * hi
    )

    x1 = Fraction(delta, 1) / d_lo
    x2 = Fraction(delta, 1) / d_hi

    lo_x = min(x1, x2)
    hi_x = max(x1, x2)

    return lo_x, hi_x


def intersect_intervals(intervals):

    """
    Exact intersection.
    """

    if not intervals:
        return None

    lo = intervals[0][0]
    hi = intervals[0][1]

    for a, b in intervals[1:]:

        if a > lo:
            lo = a

        if b < hi:
            hi = b

        if lo > hi:
            return None

    return lo, hi


def integer_candidates_from_interval(interval):

    if interval is None:
        return []

    lo, hi = interval

    # ceil(lo)
    lo_int = (
        (lo.numerator + lo.denominator - 1)
        // lo.denominator
    )

    # floor(hi)
    hi_int = (
        hi.numerator
        // hi.denominator
    )

    if lo_int > hi_int:
        return []

    return list(
        range(lo_int, hi_int + 1)
    )


# ================================================================
# Recover k0 candidates
# ================================================================

def recover_k_candidates(
    ratio_intervals,
    deltas,
):

    """
    Enumerate only the small radix residue-difference states.

    There is NO sqrt(n)-scale k0 scan.
    """

    by_delta = {}

    for delta in deltas:

        interval_parts = []

        possible = True

        for i in (1, 2):

            # We have three independent estimates:
            #
            # column 0, 1, 2.
            #
            # Intersect all of them.

            for j in range(3):

                ratio = ratio_intervals[
                    (i, j)
                ]

                part = delta_to_base_interval(
                    delta,
                    R1[0],
                    R1[i],
                    ratio
                )

                if part is None:

                    # delta=0 carries no magnitude information.
                    continue

                interval_parts.append(part)

        if not interval_parts:

            by_delta[delta] = None
            continue

        intersection = intersect_intervals(
            interval_parts
        )

        if intersection is None:
            continue

        by_delta[delta] = intersection

    return by_delta


def combine_k_states(
    delta1_intervals,
    delta2_intervals,
):

    states = []

    for d1, interval1 in delta1_intervals.items():

        if interval1 is None:
            continue

        for d2, interval2 in delta2_intervals.items():

            if interval2 is None:
                continue

            inter = intersect_intervals(
                [
                    interval1,
                    interval2,
                ]
            )

            if inter is None:
                continue

            for k0 in integer_candidates_from_interval(
                inter
            ):

                # Recover k1,k2 exactly.
                num1 = R1[0] * k0 - d1
                num2 = R1[0] * k0 - d2

                if num1 % R1[1] != 0:
                    continue

                if num2 % R1[2] != 0:
                    continue

                k1 = num1 // R1[1]
                k2 = num2 // R1[2]

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

    # Remove duplicates.
    return sorted(set(states))


# ================================================================
# Recover ell candidates
# ================================================================

def recover_ell_candidates(
    ratio_intervals,
    epsilons,
):

    by_epsilon = {}

    for epsilon in epsilons:

        interval_parts = []

        for j in (1, 2):

            for i in range(3):

                ratio = ratio_intervals[
                    (j, i)
                ]

                part = delta_to_base_interval(
                    epsilon,
                    R2[0],
                    R2[j],
                    ratio
                )

                if part is None:
                    continue

                interval_parts.append(part)

        if not interval_parts:

            by_epsilon[epsilon] = None
            continue

        intersection = intersect_intervals(
            interval_parts
        )

        if intersection is None:
            continue

        by_epsilon[epsilon] = intersection

    return by_epsilon


def combine_ell_states(
    epsilon1_intervals,
    epsilon2_intervals,
):

    states = []

    for e1, interval1 in epsilon1_intervals.items():

        if interval1 is None:
            continue

        for e2, interval2 in epsilon2_intervals.items():

            if interval2 is None:
                continue

            inter = intersect_intervals(
                [
                    interval1,
                    interval2,
                ]
            )

            if inter is None:
                continue

            for ell0 in integer_candidates_from_interval(
                inter
            ):

                num1 = R2[0] * ell0 - e1
                num2 = R2[0] * ell0 - e2

                if num1 % R2[1] != 0:
                    continue

                if num2 % R2[2] != 0:
                    continue

                ell1 = num1 // R2[1]
                ell2 = num2 // R2[2]

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

    return sorted(set(states))


# ================================================================
# Convert quotient vectors -> actual p,q candidates
# ================================================================

def expand_k_state(state):

    k0, k1, k2, d1, d2 = state

    candidates = []

    # d_i = a_i - a0
    #
    # so
    #
    # a_i = a0 + d_i.

    for a0 in range(R1[0]):

        a1 = a0 + d1
        a2 = a0 + d2

        if not (0 <= a1 < R1[1]):
            continue

        if not (0 <= a2 < R1[2]):
            continue

        p = R1[0] * k0 + a0

        # Exact consistency.
        if p // R1[1] != k1:
            continue

        if p // R1[2] != k2:
            continue

        candidates.append(p)

    return candidates


def expand_ell_state(state):

    ell0, ell1, ell2, e1, e2 = state

    candidates = []

    # epsilon_j = b_j - b0

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
# Exact final Q-box verification
# ================================================================

def verify_candidate(n, p, q, K_intervals):

    if p <= 0 or q <= 0:
        return False

    if p * q != n:
        return False

    K, A, L, B = build_true_K(p, q)

    for i in range(3):

        for j in range(3):

            kij = K[i] * L[j]

            lo, hi = K_intervals[i][j]

            if not (lo <= kij <= hi):
                return False

    return True


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 144")
print("Exact carry-bound ratio lifting")
print("=" * 72)
print()

print("R1 =", R1)
print("R2 =", R2)
print(
    "BALANCE RATIO =",
    BALANCE_RATIO
)
print()

for bits in BIT_SIZES:

    case_start = time.perf_counter()

    print("=" * 72)
    print(f"GENERATING {bits}-BIT BALANCED SEMIPRIME")
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

    # ------------------------------------------------------------
    # True hidden structure.
    # ------------------------------------------------------------

    K_true, A_true, L_true, B_true = \
        build_true_K(
            p,
            q
        )

    E_true = build_true_E(
        n,
        p,
        q
    )

    print("TRUE k =", K_true)
    print("TRUE a =", A_true)

    print("TRUE ell =", L_true)
    print("TRUE b   =", B_true)

    print()

    # Exact radix differences.
    delta1_true = (
        R1[0] * K_true[0]
        - R1[1] * K_true[1]
    )

    delta2_true = (
        R1[0] * K_true[0]
        - R1[2] * K_true[2]
    )

    eps1_true = (
        R2[0] * L_true[0]
        - R2[1] * L_true[1]
    )

    eps2_true = (
        R2[0] * L_true[0]
        - R2[2] * L_true[2]
    )

    print(
        "TRUE DELTAS =",
        delta1_true,
        delta2_true
    )

    print(
        "TRUE EPSILONS =",
        eps1_true,
        eps2_true
    )

    print()

    # ------------------------------------------------------------
    # Observable Q.
    # ------------------------------------------------------------

    Q = build_Q(n)

    print_matrix(
        "Q MATRIX",
        Q
    )

    # ------------------------------------------------------------
    # Build rigorous K intervals from n alone.
    # ------------------------------------------------------------

    K_intervals = build_K_intervals(
        n,
        Q
    )

    print("OBSERVABLE K INTERVALS")
    print()

    for i in range(3):

        for j in range(3):

            lo, hi = K_intervals[i][j]

            Bbound = Q[i][j] - lo

            print(
                f"K[{i},{j}] in "
                f"[{lo}, {hi}] "
                f"(carry bound {Bbound})"
            )

    print()

    # ------------------------------------------------------------
    # Ratio intervals for rows.
    # ------------------------------------------------------------

    row_ratios = row_ratio_intervals(
        K_intervals
    )

    print("ROW RATIO INTERVALS")
    print()

    for key in sorted(row_ratios):

        lo, hi = row_ratios[key]

        print(
            f"k{key[0]}/k0 "
            f"from column {key[1]}: "
            f"[{float(lo):.15e}, "
            f"{float(hi):.15e}]"
        )

    print()

    # ------------------------------------------------------------
    # Ratio intervals for columns.
    # ------------------------------------------------------------

    column_ratios = column_ratio_intervals(
        K_intervals
    )

    print("COLUMN RATIO INTERVALS")
    print()

    for key in sorted(column_ratios):

        lo, hi = column_ratios[key]

        print(
            f"ell{key[0]}/ell0 "
            f"from row {key[1]}: "
            f"[{float(lo):.15e}, "
            f"{float(hi):.15e}]"
        )

    print()

    # ------------------------------------------------------------
    # Delta ranges from radix residues.
    # ------------------------------------------------------------

    delta1_values = range(
        -(R1[0] - 1),
        R1[1]
    )

    delta2_values = range(
        -(R1[0] - 1),
        R1[2]
    )

    epsilon1_values = range(
        -(R2[0] - 1),
        R2[1]
    )

    epsilon2_values = range(
        -(R2[0] - 1),
        R2[2]
    )

    print("DELTA SEARCH SIZE")
    print(
        "delta1 states =",
        len(delta1_values)
    )
    print(
        "delta2 states =",
        len(delta2_values)
    )
    print(
        "combined k residue states =",
        len(delta1_values)
        * len(delta2_values)
    )

    print()

    # ------------------------------------------------------------
    # Turn exact ratio intervals into k0 intervals.
    # ------------------------------------------------------------

    k0_from_delta1 = recover_k_candidates(
        {
            (1, 0): row_ratios[(1, 0)],
            (1, 1): row_ratios[(1, 1)],
            (1, 2): row_ratios[(1, 2)],
            (2, 0): row_ratios[(2, 0)],
            (2, 1): row_ratios[(2, 1)],
            (2, 2): row_ratios[(2, 2)],
        },
        delta1_values
    )

    # The function above intentionally treats every ratio under
    # a single delta.  For delta-specific row lifting we need the
    # row-1 ratios only.
    k_delta1 = {}

    for d in delta1_values:

        parts = []

        for j in range(3):

            part = delta_to_base_interval(
                d,
                R1[0],
                R1[1],
                row_ratios[(1, j)]
            )

            if part is not None:
                parts.append(part)

        if parts:

            inter = intersect_intervals(parts)

            if inter is not None:
                k_delta1[d] = inter

    k_delta2 = {}

    for d in delta2_values:

        parts = []

        for j in range(3):

            part = delta_to_base_interval(
                d,
                R1[0],
                R1[2],
                row_ratios[(2, j)]
            )

            if part is not None:
                parts.append(part)

        if parts:

            inter = intersect_intervals(parts)

            if inter is not None:
                k_delta2[d] = inter

    # ------------------------------------------------------------
    # Combine delta states.
    # ------------------------------------------------------------

    k_states = combine_k_states(
        k_delta1,
        k_delta2
    )

    print("K LIFT RESULTS")
    print()

    print(
        "delta1 surviving states =",
        len(k_delta1)
    )

    print(
        "delta2 surviving states =",
        len(k_delta2)
    )

    print(
        "integer K states =",
        len(k_states)
    )

    print()

    # ------------------------------------------------------------
    # Epsilon lifting.
    # ------------------------------------------------------------

    ell_epsilon1 = {}

    for e in epsilon1_values:

        parts = []

        for i in range(3):

            part = delta_to_base_interval(
                e,
                R2[0],
                R2[1],
                column_ratios[(1, i)]
            )

            if part is not None:
                parts.append(part)

        if parts:

            inter = intersect_intervals(parts)

            if inter is not None:
                ell_epsilon1[e] = inter

    ell_epsilon2 = {}

    for e in epsilon2_values:

        parts = []

        for i in range(3):

            part = delta_to_base_interval(
                e,
                R2[0],
                R2[2],
                column_ratios[(2, i)]
            )

            if part is not None:
                parts.append(part)

        if parts:

            inter = intersect_intervals(parts)

            if inter is not None:
                ell_epsilon2[e] = inter

    ell_states = combine_ell_states(
        ell_epsilon1,
        ell_epsilon2
    )

    print("ELL LIFT RESULTS")
    print()

    print(
        "epsilon1 surviving states =",
        len(ell_epsilon1)
    )

    print(
        "epsilon2 surviving states =",
        len(ell_epsilon2)
    )

    print(
        "integer ELL states =",
        len(ell_states)
    )

    print()

    # ------------------------------------------------------------
    # Expand quotient states into p/q candidates.
    # ------------------------------------------------------------

    p_candidates = set()

    for state in k_states:

        for pc in expand_k_state(state):

            p_candidates.add(pc)

    q_candidates = set()

    for state in ell_states:

        for qc in expand_ell_state(state):

            q_candidates.add(qc)

    print("FACTOR CANDIDATE EXPANSION")
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
    # Exact factor verification.
    # ------------------------------------------------------------

    verified = []

    for pc in p_candidates:

        # Since n = p*q, q is forced.
        if n % pc != 0:
            continue

        qc = n // pc

        if qc not in q_candidates:
            continue

        if verify_candidate(
            n,
            pc,
            qc,
            K_intervals
        ):
            verified.append(
                (pc, qc)
            )

    verified = sorted(
        set(verified)
    )

    print("EXACT VERIFICATION")
    print()

    print(
        "verified factor pairs =",
        len(verified)
    )

    for pair in verified[:10]:

        print(
            "  ",
            pair
        )

    print()

    # ------------------------------------------------------------
    # Hidden truth membership.
    # ------------------------------------------------------------

    true_k_state_present = False

    for state in k_states:

        if (
            state[0] == K_true[0]
            and
            state[1] == K_true[1]
            and
            state[2] == K_true[2]
        ):
            true_k_state_present = True
            break

    true_ell_state_present = False

    for state in ell_states:

        if (
            state[0] == L_true[0]
            and
            state[1] == L_true[1]
            and
            state[2] == L_true[2]
        ):
            true_ell_state_present = True
            break

    print(
        "TRUE K VECTOR PRESENT =",
        true_k_state_present
    )

    print(
        "TRUE ELL VECTOR PRESENT =",
        true_ell_state_present
    )

    print()

    # ------------------------------------------------------------
    # Interval width diagnostics for the true delta state.
    # ------------------------------------------------------------

    print("TRUE-STATE INTERVAL DIAGNOSTICS")
    print()

    if delta1_true in k_delta1:

        lo, hi = k_delta1[
            delta1_true
        ]

        print(
            "delta1 true interval:"
        )

        print(
            "  [",
            float(lo),
            ",",
            float(hi),
            "]"
        )

        print(
            "  width =",
            float(hi - lo)
        )

    else:

        print(
            "delta1 true interval: EMPTY"
        )

    if delta2_true in k_delta2:

        lo, hi = k_delta2[
            delta2_true
        ]

        print(
            "delta2 true interval:"
        )

        print(
            "  [",
            float(lo),
            ",",
            float(hi),
            "]"
        )

        print(
            "  width =",
            float(hi - lo)
        )

    else:

        print(
            "delta2 true interval: EMPTY"
        )

    print()

    if eps1_true in ell_epsilon1:

        lo, hi = ell_epsilon1[
            eps1_true
        ]

        print(
            "epsilon1 true interval:"
        )

        print(
            "  [",
            float(lo),
            ",",
            float(hi),
            "]"
        )

        print(
            "  width =",
            float(hi - lo)
        )

    else:

        print(
            "epsilon1 true interval: EMPTY"
        )

    if eps2_true in ell_epsilon2:

        lo, hi = ell_epsilon2[
            eps2_true
        ]

        print(
            "epsilon2 true interval:"
        )

        print(
            "  [",
            float(lo),
            ",",
            float(hi),
            "]"
        )

        print(
            "  width =",
            float(hi - lo)
        )

    else:

        print(
            "epsilon2 true interval: EMPTY"
        )

    print()

    # ------------------------------------------------------------
    # Square-root control.
    # ------------------------------------------------------------

    sqrt_n = math.isqrt(n)

    print("COMPARISON")
    print()

    print(
        "sqrt(n) =",
        sqrt_n
    )

    print(
        "K lift states =",
        len(k_states)
    )

    print(
        "ELL lift states =",
        len(ell_states)
    )

    print(
        "verified factor pairs =",
        len(verified)
    )

    if len(k_states) > 0:

        print(
            "sqrt(n) / K-states =",
            f"{sqrt_n / len(k_states):.6f}"
        )

    print()

    elapsed = time.perf_counter() - case_start

    print(
        f"runtime = {elapsed:.6f} s"
    )

    print()
    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 144")
print("=" * 72)
