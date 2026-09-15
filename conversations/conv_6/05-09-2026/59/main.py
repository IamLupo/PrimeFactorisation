# ================================================================
# START EXPERIMENT 142
# Nested radix / base-B digit carry propagation
# ================================================================

import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

BASE = 7

R1 = [BASE, BASE**2, BASE**3]
R2 = [BASE, BASE**2, BASE**3]

BIT_SIZES = [30, 36, 42, 48, 54]

# Maximum number of low-order base-B digits to propagate.
MAX_LOCAL_DIGITS = 8

# Safety limit for local states.
MAX_STATES = 2_000_000


# ================================================================
# Miller-Rabin
# ================================================================

def is_probable_prime(n):
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
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

    # Deterministic enough for our experimental range.
    witnesses = [2, 3, 5, 7, 11, 13, 17]

    for a in witnesses:
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

        x |= (1 << (bits - 1))
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
# Base-B utilities
# ================================================================

def base_digits(x, base, count):
    digits = []
    y = x

    for _ in range(count):
        digits.append(y % base)
        y //= base

    return digits


def reconstruct_digits(digits, base):
    value = 0
    power = 1

    for d in digits:
        value += d * power
        power *= base

    return value


# ================================================================
# Exact base-B multiplication
# ================================================================

def multiplication_raw_coefficients(p_digits, q_digits):
    """
    raw[t] = sum_{i+j=t} p_i q_j
    """

    length = len(p_digits) + len(q_digits) - 1

    raw = [0] * length

    for i, pi in enumerate(p_digits):
        for j, qj in enumerate(q_digits):
            raw[i + j] += pi * qj

    return raw


def propagate_carries(raw, base):
    """
    Given raw base-B coefficients, produce:
        digits[t]
        carries[t+1]

    such that

        raw[t] + carry[t]
          = digit[t] + base * carry[t+1].
    """

    digits = [0] * len(raw)
    carries = [0] * (len(raw) + 1)

    for t in range(len(raw)):

        value = raw[t] + carries[t]

        digits[t] = value % base
        carries[t + 1] = value // base

    return digits, carries


# ================================================================
# Nested-radix quotient/carry matrices
# ================================================================

def build_nested_matrices(n, p, q):

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

    Q = []
    E = []

    for i, r in enumerate(R1):

        qrow = []
        erow = []

        for j, s in enumerate(R2):

            modulus = r * s

            qij = n // modulus
            eij = qij - K[i] * L[j]

            qrow.append(qij)
            erow.append(eij)

        Q.append(qrow)
        E.append(erow)

    return K, A, L, B, Q, E


# ================================================================
# Correct local digit convolution
# ================================================================

def local_convolution(pd, qd, position):
    """
    Calculate

        sum_{i+j=position} pd[i] * qd[j]

    using only the currently known digits.
    """

    total = 0

    for i in range(position + 1):

        j = position - i

        if i >= len(pd):
            continue

        if j >= len(qd):
            continue

        total += pd[i] * qd[j]

    return total


# ================================================================
# Local carry-state enumeration
# ================================================================

def enumerate_local_digit_states(
    target_digits,
    base,
    max_digits,
    true_p_digits,
    true_q_digits,
):
    """
    Enumerate all low-order digit pairs satisfying the exact
    multiplication equations.

    At position t:

        convolution_t + carry_in
            = target_digit_t + base * carry_out

    Only low-order digits are considered.

    Returns:
        final states
        branching counts
        whether the true prefix survives
        whether swapped prefix survives
    """

    states = [
        {
            "pd": [],
            "qd": [],
            "carry": 0,
        }
    ]

    branch_counts = []

    positions = min(
        max_digits,
        len(target_digits)
    )

    stopped_early = False

    for pos in range(positions):

        target_digit = target_digits[pos]

        new_states = []

        for state in states:

            carry_in = state["carry"]
            old_pd = state["pd"]
            old_qd = state["qd"]

            for pd_digit in range(base):

                for qd_digit in range(base):

                    # Build the candidate prefix first.
                    candidate_pd = old_pd + [pd_digit]
                    candidate_qd = old_qd + [qd_digit]

                    conv = local_convolution(
                        candidate_pd,
                        candidate_qd,
                        pos
                    )

                    value = conv + carry_in

                    if value % base != target_digit:
                        continue

                    carry_out = value // base

                    new_states.append(
                        {
                            "pd": candidate_pd,
                            "qd": candidate_qd,
                            "carry": carry_out,
                        }
                    )

        states = new_states

        branch_counts.append(len(states))

        print(
            f"digit position {pos}: "
            f"states = {len(states)}"
        )

        if len(states) > MAX_STATES:

            print(
                "STATE EXPLOSION: "
                f"{len(states)} > {MAX_STATES}"
            )

            stopped_early = True
            break

        if len(states) == 0:

            print(
                "NO SURVIVING STATES"
            )

            break

    # ------------------------------------------------------------
    # True-prefix test
    # ------------------------------------------------------------

    effective_positions = min(
        len(states[0]["pd"]) if states else 0,
        len(true_p_digits),
        len(true_q_digits),
    )

    true_survives = False
    swapped_survives = False

    if states:

        true_prefix_p = true_p_digits[:effective_positions]
        true_prefix_q = true_q_digits[:effective_positions]

        swapped_prefix_p = true_q_digits[:effective_positions]
        swapped_prefix_q = true_p_digits[:effective_positions]

        for state in states:

            if (
                state["pd"] == true_prefix_p
                and
                state["qd"] == true_prefix_q
            ):
                true_survives = True

            if (
                state["pd"] == swapped_prefix_p
                and
                state["qd"] == swapped_prefix_q
            ):
                swapped_survives = True

    return (
        states,
        branch_counts,
        true_survives,
        swapped_survives,
        stopped_early,
    )


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
# EXPERIMENT
# ================================================================

print("=" * 72)
print("START EXPERIMENT 142")
print("Nested radix / base-B digit carry propagation")
print("=" * 72)
print()

print("BASE =", BASE)
print("R1   =", R1)
print("R2   =", R2)
print()

for bits in BIT_SIZES:

    print("=" * 72)
    print(f"GENERATING {bits}-BIT SEMIPRIME")
    print("=" * 72)

    case_start = time.perf_counter()

    p, q = random_semiprime(bits)
    n = p * q

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)
    print()

    # ------------------------------------------------------------
    # Base-B representations
    # ------------------------------------------------------------

    p_digits_count = max(
        1,
        math.ceil(math.log(p + 1, BASE))
    )

    q_digits_count = max(
        1,
        math.ceil(math.log(q + 1, BASE))
    )

    p_digits = base_digits(
        p,
        BASE,
        p_digits_count
    )

    q_digits = base_digits(
        q,
        BASE,
        q_digits_count
    )

    print("p base-B digits:")
    print(p_digits)
    print()

    print("q base-B digits:")
    print(q_digits)
    print()

    assert reconstruct_digits(
        p_digits,
        BASE
    ) == p

    assert reconstruct_digits(
        q_digits,
        BASE
    ) == q

    # ------------------------------------------------------------
    # Exact product convolution
    # ------------------------------------------------------------

    raw = multiplication_raw_coefficients(
        p_digits,
        q_digits
    )

    product_digits, true_carries = propagate_carries(
        raw,
        BASE
    )

    reconstructed_n = reconstruct_digits(
        product_digits,
        BASE
    )

    print("RAW CONVOLUTION COEFFICIENTS")
    print(raw)
    print()

    print("PRODUCT DIGITS")
    print(product_digits)
    print()

    print("TRUE BASE-B CARRIES")
    print(true_carries)
    print()

    print(
        "BASE-B PRODUCT RECONSTRUCTION:",
        "PASS"
        if reconstructed_n == n
        else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Nested radix matrices
    # ------------------------------------------------------------

    K, A, L, B, Q, E = build_nested_matrices(
        n,
        p,
        q
    )

    print("NESTED-RADIX QUOTIENT DATA")
    print()

    print("K =", K)
    print("A =", A)
    print("L =", L)
    print("B =", B)
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
    # Successive quotient differences
    # ------------------------------------------------------------

    print("SUCCESSIVE QUOTIENT DIFFERENCES")
    print()

    for j in range(3):

        d01 = Q[0][j] - BASE * Q[1][j]
        d12 = Q[1][j] - BASE * Q[2][j]

        print(
            f"column {j}: "
            f"Q0 - B*Q1 = {d01}, "
            f"Q1 - B*Q2 = {d12}"
        )

    print()

    for i in range(3):

        d01 = Q[i][0] - BASE * Q[i][1]
        d12 = Q[i][1] - BASE * Q[i][2]

        print(
            f"row {i}: "
            f"Q0 - B*Q1 = {d01}, "
            f"Q1 - B*Q2 = {d12}"
        )

    print()

    # ------------------------------------------------------------
    # Carry differences
    # ------------------------------------------------------------

    print("SUCCESSIVE CARRY-MATRIX DIFFERENCES")
    print()

    for j in range(3):

        d01 = E[0][j] - BASE * E[1][j]
        d12 = E[1][j] - BASE * E[2][j]

        print(
            f"column {j}: "
            f"E0 - B*E1 = {d01}, "
            f"E1 - B*E2 = {d12}"
        )

    print()

    for i in range(3):

        d01 = E[i][0] - BASE * E[i][1]
        d12 = E[i][1] - BASE * E[i][2]

        print(
            f"row {i}: "
            f"E0 - B*E1 = {d01}, "
            f"E1 - B*E2 = {d12}"
        )

    print()

    # ------------------------------------------------------------
    # Exact local multiplication equations
    # ------------------------------------------------------------

    print("BASE-B LOCAL MULTIPLICATION EQUATIONS")
    print()

    for t in range(len(raw)):

        raw_t = raw[t]
        carry_in = true_carries[t]
        digit = product_digits[t]
        carry_out = true_carries[t + 1]

        lhs = raw_t + carry_in
        rhs = digit + BASE * carry_out

        status = "PASS" if lhs == rhs else "FAIL"

        print(
            f"position {t}: "
            f"raw={raw_t}, "
            f"carry_in={carry_in}, "
            f"digit={digit}, "
            f"carry_out={carry_out}, "
            f"{lhs} = {rhs} [{status}]"
        )

    print()

    # ------------------------------------------------------------
    # Low-order information in n
    # ------------------------------------------------------------

    print("LOW-BASE INFORMATION")
    print()

    prefix = 0
    power = 1

    for m in range(1, 4):

        prefix += product_digits[m - 1] * power
        power *= BASE

        modulus = BASE ** m

        print(
            f"B^{m} = {modulus:<5} "
            f"n mod B^{m} = {n % modulus:<5} "
            f"digit-prefix = {prefix:<5}"
        )

    print()

    # ------------------------------------------------------------
    # Local state enumeration
    # ------------------------------------------------------------

    print("LOCAL DIGIT-CARRY STATE ENUMERATION")
    print()

    (
        states,
        branch_counts,
        true_survives,
        swapped_survives,
        stopped_early,
    ) = enumerate_local_digit_states(
        target_digits=product_digits,
        base=BASE,
        max_digits=MAX_LOCAL_DIGITS,
        true_p_digits=p_digits,
        true_q_digits=q_digits,
    )

    print()

    print("LOCAL BRANCH COUNTS")
    print(branch_counts)
    print()

    print("FINAL LOCAL STATES =", len(states))
    print()

    print(
        "TRUE PREFIX SURVIVES =",
        true_survives
    )

    print(
        "SWAPPED PREFIX SURVIVES =",
        swapped_survives
    )

    print(
        "STATE ENUMERATION STOPPED EARLY =",
        stopped_early
    )

    print()

    # ------------------------------------------------------------
    # Display a few surviving states.
    # ------------------------------------------------------------

    display_count = min(10, len(states))

    if display_count:

        print(
            f"FIRST {display_count} SURVIVING STATES"
        )

        for idx in range(display_count):

            state = states[idx]

            p_partial = reconstruct_digits(
                state["pd"],
                BASE
            )

            q_partial = reconstruct_digits(
                state["qd"],
                BASE
            )

            print(
                f"{idx:3d}: "
                f"p_prefix={p_partial:<10} "
                f"q_prefix={q_partial:<10} "
                f"carry={state['carry']:<10}"
            )

        print()

    # ------------------------------------------------------------
    # Compare against sqrt(n)
    # ------------------------------------------------------------

    sqrt_n = math.isqrt(n)

    print("SEARCH-SCALE COMPARISON")
    print()

    print("sqrt(n) =", sqrt_n)
    print(
        "local digit states =",
        len(states)
    )

    if len(states) > 0:

        print(
            "sqrt/state ratio =",
            f"{sqrt_n / len(states):.6f}"
        )

    print()

    elapsed = time.perf_counter() - case_start

    print(
        f"runtime = {elapsed:.6f} s"
    )

    print()
    print("-" * 72)
    print(f"COMPLETED {bits}-BIT CASE")
    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 142")
print("=" * 72)