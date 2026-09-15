# ================================================================
# START EXPERIMENT 150
# Incremental discriminant CRT compression
#
# Goal:
#
#   Avoid the 21-billion-state Cartesian-product explosion from
#   Experiment 149.
#
#   We incrementally combine
#
#       S^2 - d^2 = 4n (mod r)
#
#   one radix at a time, while monitoring the number of distinct
#   (S,d) CRT states.
#
#   We also simultaneously maintain the equivalent (p,q) CRT
#   representation.
#
#   The purpose is to determine whether the discriminant
#   representation provides ANY state compression before the full
#   CRT modulus is reached.
#
# Important:
#
#   We stop automatically if the state count exceeds MAX_STATES.
#   Therefore this experiment cannot consume all system memory.
#
# ================================================================

import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

RADICES = [
    17,
    43,
    59,
    71,
    83,
    97,
]

BIT_SIZES = [
    30,
    36,
    42,
    48,
    54,
]

# Hard safety limit.
MAX_STATES = 2_000_000


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

    return (
        x % (m1 * m2),
        m1 * m2,
    )


# ================================================================
# Local discriminant states
# ================================================================

def local_sd_states(n, r):

    """
    All solutions of

        S^2 - d^2 = 4n (mod r).

    Since r is prime and n is not divisible by r,
    there are exactly r-1 local states.
    """

    target = (
        4 * n
    ) % r

    states = []

    for s in range(r):

        s2 = (
            s * s
        ) % r

        for d in range(r):

            if (
                (s2 - d * d) % r
                == target
            ):

                states.append(
                    (s, d)
                )

    return states


# ================================================================
# Local factor states
# ================================================================

def local_pq_states(n, r):

    """
    All invertible p residues.

    q is forced by

        p*q = n (mod r).
    """

    states = []

    for p in range(1, r):

        if math.gcd(p, r) != 1:
            continue

        q = (
            n
            * pow(
                p,
                -1,
                r
            )
        ) % r

        states.append(
            (p, q)
        )

    return states


# ================================================================
# Incremental SD CRT combination
# ================================================================

def incremental_sd_crt(n):

    """
    Incrementally combine local (S,d) states.

    Returns:
        history
        states
        final_modulus
        stopped
    """

    current = [
        (
            0,
            0,
        )
    ]

    modulus = 1

    history = []

    stopped = False

    for r in RADICES:

        local = local_sd_states(
            n,
            r
        )

        next_states = []

        new_modulus = (
            modulus * r
        )

        for S0, d0 in current:

            for s, d in local:

                S_new, _ = crt_pair(
                    S0,
                    modulus,
                    s,
                    r
                )

                d_new, _ = crt_pair(
                    d0,
                    modulus,
                    d,
                    r
                )

                next_states.append(
                    (
                        S_new,
                        d_new,
                    )
                )

                if (
                    len(next_states)
                    > MAX_STATES
                ):

                    stopped = True
                    break

            if stopped:
                break

        state_count = len(next_states)

        history.append(
            {
                "r": r,
                "local": len(local),
                "states": state_count,
                "modulus": new_modulus,
                "stopped": stopped,
            }
        )

        if stopped:

            return (
                history,
                next_states,
                new_modulus,
                True,
            )

        current = sorted(
            set(next_states)
        )

        modulus = new_modulus

    return (
        history,
        current,
        modulus,
        False,
    )


# ================================================================
# Incremental factor CRT combination
# ================================================================

def incremental_pq_crt(n):

    """
    Directly combine (p,q) residue states.
    """

    current = [
        (
            0,
            0,
        )
    ]

    modulus = 1

    history = []

    stopped = False

    for r in RADICES:

        local = local_pq_states(
            n,
            r
        )

        next_states = []

        new_modulus = (
            modulus * r
        )

        for p0, q0 in current:

            for p, q in local:

                p_new, _ = crt_pair(
                    p0,
                    modulus,
                    p,
                    r
                )

                q_new, _ = crt_pair(
                    q0,
                    modulus,
                    q,
                    r
                )

                next_states.append(
                    (
                        p_new,
                        q_new,
                    )
                )

                if (
                    len(next_states)
                    > MAX_STATES
                ):

                    stopped = True
                    break

            if stopped:
                break

        state_count = len(next_states)

        history.append(
            {
                "r": r,
                "local": len(local),
                "states": state_count,
                "modulus": new_modulus,
                "stopped": stopped,
            }
        )

        if stopped:

            return (
                history,
                next_states,
                new_modulus,
                True,
            )

        current = sorted(
            set(next_states)
        )

        modulus = new_modulus

    return (
        history,
        current,
        modulus,
        False,
    )


# ================================================================
# Convert (S,d) -> (p,q)
# ================================================================

def sd_to_pq(S, d, modulus):

    inv2 = pow(
        2,
        -1,
        modulus
    )

    p = (
        (S + d)
        * inv2
    ) % modulus

    q = (
        (S - d)
        * inv2
    ) % modulus

    return p, q


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 150")
print("Incremental discriminant CRT compression")
print("=" * 72)
print()

print(
    "RADICES =",
    RADICES
)

print(
    "MAX_STATES =",
    MAX_STATES
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

    S_true = p + q
    d_true = p - q

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

    print(
        "true S =",
        S_true
    )

    print(
        "true d =",
        d_true
    )

    print()

    # ------------------------------------------------------------
    # Incremental SD CRT
    # ------------------------------------------------------------

    sd_start = time.perf_counter()

    (
        sd_history,
        sd_states,
        sd_modulus,
        sd_stopped
    ) = incremental_sd_crt(n)

    sd_elapsed = (
        time.perf_counter()
        - sd_start
    )

    print(
        "INCREMENTAL (S,d) CRT"
    )

    print()

    for entry in sd_history:

        print(
            f"r={entry['r']:3d} "
            f"local={entry['local']:4d} "
            f"states={entry['states']:>10d} "
            f"modulus={entry['modulus']:>15d} "
            f"stopped={entry['stopped']}"
        )

    print()

    print(
        "SD FINAL STATES =",
        len(sd_states)
    )

    print(
        "SD FINAL MODULUS =",
        sd_modulus
    )

    print(
        "SD STOPPED =",
        sd_stopped
    )

    print(
        "SD RUNTIME =",
        f"{sd_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # Direct p,q CRT
    # ------------------------------------------------------------

    pq_start = time.perf_counter()

    (
        pq_history,
        pq_states,
        pq_modulus,
        pq_stopped
    ) = incremental_pq_crt(n)

    pq_elapsed = (
        time.perf_counter()
        - pq_start
    )

    print(
        "INCREMENTAL (p,q) CRT"
    )

    print()

    for entry in pq_history:

        print(
            f"r={entry['r']:3d} "
            f"local={entry['local']:4d} "
            f"states={entry['states']:>10d} "
            f"modulus={entry['modulus']:>15d} "
            f"stopped={entry['stopped']}"
        )

    print()

    print(
        "PQ FINAL STATES =",
        len(pq_states)
    )

    print(
        "PQ FINAL MODULUS =",
        pq_modulus
    )

    print(
        "PQ STOPPED =",
        pq_stopped
    )

    print(
        "PQ RUNTIME =",
        f"{pq_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # True state membership for the largest completed modulus.
    # ------------------------------------------------------------

    if not sd_stopped:

        true_sd = (
            S_true % sd_modulus,
            d_true % sd_modulus
        )

        print(
            "TRUE SD STATE PRESENT =",
            true_sd in sd_states
        )

        print()

    if not pq_stopped:

        true_pq = (
            p % pq_modulus,
            q % pq_modulus
        )

        print(
            "TRUE PQ STATE PRESENT =",
            true_pq in pq_states
        )

        print()

    # ------------------------------------------------------------
    # If both completed, verify exact equivalence.
    # ------------------------------------------------------------

    if (
        not sd_stopped
        and
        not pq_stopped
    ):

        converted = set()

        for S, d in sd_states:

            converted.add(
                sd_to_pq(
                    S,
                    d,
                    sd_modulus
                )
            )

        pq_set = set(
            pq_states
        )

        print(
            "SD -> PQ STATE COUNT =",
            len(converted)
        )

        print(
            "DIRECT PQ STATE COUNT =",
            len(pq_set)
        )

        print(
            "EXACT STATE EQUIVALENCE =",
            converted == pq_set
        )

        print()

    # ------------------------------------------------------------
    # Theoretical local product.
    # ------------------------------------------------------------

    expected_local_product = 1

    for r in RADICES:

        expected_local_product *= (
            r - 1
        )

    print(
        "THEORETICAL FULL CRT STATE COUNT =",
        expected_local_product
    )

    print()

    # ------------------------------------------------------------
    # Compare state growth to modulus.
    # ------------------------------------------------------------

    print(
        "STATE GROWTH SUMMARY"
    )

    print()

    for entry in sd_history:

        density = (
            entry["states"]
            /
            entry["modulus"]
        )

        print(
            f"SD r={entry['r']:3d}: "
            f"states/modulus={density:.12e}"
        )

    print()

    # ------------------------------------------------------------
    # Runtime
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
print("FINISHED EXPERIMENT 150")
print("=" * 72)
