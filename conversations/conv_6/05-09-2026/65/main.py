# ================================================================
# START EXPERIMENT 148
# Multi-radix discriminant congruence intersection
#
# Goal:
#
#   Let
#
#       S = p + q
#
#   Then
#
#       S^2 - 4n = (p-q)^2.
#
#   Therefore, for every odd prime radix r:
#
#       S^2 = 4n (mod r).
#
#   Each radix therefore gives at most two possible residues
#   for S modulo r.
#
#   We combine all residue choices with CRT and then lift S only
#   inside the mathematically bounded balanced-factor interval:
#
#       2*sqrt(n) <= S <= 5/2*sqrt(n)
#
#   when max(p,q)/min(p,q) <= 4.
#
#   This is deliberately NOT a p-scan.
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

BALANCE_RATIO = 4


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
# Number-theoretic helpers
# ================================================================

def legendre_symbol(a, p):

    a %= p

    if a == 0:
        return 0

    return pow(
        a,
        (p - 1) // 2,
        p
    )


def tonelli_shanks(n, p):
    """
    Return one square root x of n modulo odd prime p.

    Return None when n is not a quadratic residue modulo p.
    """

    n %= p

    if n == 0:
        return 0

    if p == 2:
        return n

    if legendre_symbol(n, p) != 1:
        return None

    # Fast case p = 3 mod 4.
    if p % 4 == 3:

        return pow(
            n,
            (p + 1) // 4,
            p
        )

    # Write p-1 = q * 2^s with q odd.
    q = p - 1
    s = 0

    while q % 2 == 0:
        q //= 2
        s += 1

    # Find quadratic non-residue z.
    z = 2

    while legendre_symbol(z, p) != p - 1:
        z += 1

    c = pow(z, q, p)
    x = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s

    while t != 1:

        i = 1
        t2i = (t * t) % p

        while t2i != 1:

            t2i = (t2i * t2i) % p
            i += 1

            if i >= m:
                raise RuntimeError(
                    "Tonelli-Shanks failed"
                )

        b = pow(
            c,
            1 << (m - i - 1),
            p
        )

        x = (x * b) % p

        b2 = (b * b) % p

        t = (t * b2) % p
        c = b2
        m = i

    return x


def modular_square_roots(value, p):
    """
    Return the two roots of x^2 = value (mod p)
    for odd prime p.
    """

    value %= p

    if value == 0:
        return [0]

    root = tonelli_shanks(
        value,
        p
    )

    if root is None:
        return []

    other = (-root) % p

    if root == other:
        return [root]

    return sorted([
        root,
        other,
    ])


# ================================================================
# CRT
# ================================================================

def crt_pair(a1, m1, a2, m2):
    """
    Solve

        x = a1 (mod m1)
        x = a2 (mod m2)

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


def combine_residue_classes(classes, modulus):
    """
    Starting from:

        x = residue (mod modulus)

    extend with a list of residues, one prime at a time.
    """

    current = [
        (
            0,
            1,
        )
    ]

    for residues, prime in classes:

        next_current = []

        for x, m in current:

            for r in residues:

                y, new_m = crt_pair(
                    x,
                    m,
                    r,
                    prime
                )

                next_current.append(
                    (
                        y,
                        new_m,
                    )
                )

        current = next_current

    # Remove duplicate residue classes.
    return sorted(set(current))


# ================================================================
# Integer square root interval
# ================================================================

def ceil_sqrt_fraction(num, den):
    """
    ceil(sqrt(num/den)) using integer arithmetic.
    """

    x = math.isqrt(
        num // den
    )

    while x * x * den < num:
        x += 1

    return x


# ================================================================
# Balanced S interval
# ================================================================

def balanced_sum_interval(n):
    """
    For pq=n and max(p,q)/min(p,q) <= 4:

        2*sqrt(n) <= p+q <= 5/2*sqrt(n).
    """

    sqrt_n = math.isqrt(n)

    while (sqrt_n + 1) * (sqrt_n + 1) <= n:
        sqrt_n += 1

    # Lower bound:
    #
    # S >= 2 sqrt(n).
    #
    # We need ceil(2 sqrt(n)).
    #
    lower = 2 * sqrt_n

    while lower * lower < 4 * n:
        lower += 1

    # Upper bound:
    #
    # S <= (5/2) sqrt(n)
    #
    # equivalently
    #
    # 2S <= 5 sqrt(n)
    #
    # so
    #
    # 4S^2 <= 25n.
    #
    # Use exact integer correction.

    upper = (
        5 * sqrt_n
    ) // 2

    while (
        4 * upper * upper
        > 25 * n
    ):
        upper -= 1

    while (
        4 * (upper + 1) * (upper + 1)
        <= 25 * n
    ):
        upper += 1

    return lower, upper


# ================================================================
# Lift CRT residue class inside [lo, hi]
# ================================================================

def lift_residue_class(
    residue,
    modulus,
    lo,
    hi
):

    if lo > hi:
        return []

    # Find first x >= lo with x == residue mod modulus.

    if residue >= lo:

        x = residue

    else:

        jumps = (
            (lo - residue + modulus - 1)
            // modulus
        )

        x = residue + jumps * modulus

    if x > hi:
        return []

    count = (
        (hi - x)
        // modulus
    ) + 1

    return [
        x + i * modulus
        for i in range(count)
    ]


# ================================================================
# Exact factorization from S
# ================================================================

def factor_from_sum(n, S):

    D = S * S - 4 * n

    if D < 0:
        return None

    root = math.isqrt(D)

    if root * root != D:
        return None

    if (S - root) % 2 != 0:
        return None

    p = (S - root) // 2
    q = (S + root) // 2

    if p <= 1 or q <= 1:
        return None

    if p * q != n:
        return None

    return tuple(
        sorted(
            (
                p,
                q,
            )
        )
    )


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 148")
print("Multi-radix discriminant congruence intersection")
print("=" * 72)
print()

print("RADICES =", RADICES)
print("BALANCE RATIO =", BALANCE_RATIO)
print()

modulus = 1

for r in RADICES:
    modulus *= r

print(
    "COMBINED CRT MODULUS =",
    modulus
)

print()


for bits in BIT_SIZES:

    start = time.perf_counter()

    print("=" * 72)
    print(
        f"GENERATING {bits}-BIT BALANCED SEMIPRIME"
    )
    print("=" * 72)
    print()

    p, q = random_balanced_semiprime(bits)
    n = p * q

    S_true = p + q

    print("-" * 72)
    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)
    print("true S =", S_true)
    print()

    print(
        "factor ratio =",
        f"{max(p,q) / min(p,q):.6f}"
    )

    print()

    # ------------------------------------------------------------
    # Check all radices are coprime to n.
    # ------------------------------------------------------------

    bad_radices = [
        r
        for r in RADICES
        if math.gcd(n, r) != 1
    ]

    if bad_radices:

        print(
            "RADICES SHARING FACTORS WITH n =",
            bad_radices
        )

        print(
            "Skipping case."
        )

        print()

        continue

    # ------------------------------------------------------------
    # Modular discriminant roots.
    # ------------------------------------------------------------

    classes = []

    print(
        "LOCAL DISCRIMINANT CONGRUENCES"
    )

    print()

    local_root_product = 1

    local_root_counts = []

    local_checks = []

    for r in RADICES:

        value = (
            4 * n
        ) % r

        roots = modular_square_roots(
            value,
            r
        )

        local_root_counts.append(
            len(roots)
        )

        if len(roots) == 0:

            print(
                f"r={r}: NO ROOT -> FAIL"
            )

            local_checks.append(False)

            continue

        local_root_product *= len(roots)

        true_residue = S_true % r

        true_present = (
            true_residue in roots
        )

        local_checks.append(
            true_present
        )

        print(
            f"r={r:3d} "
            f"4n mod r={value:3d} "
            f"S roots={roots} "
            f"true S mod r={true_residue} "
            f"true_present={true_present}"
        )

        classes.append(
            (
                roots,
                r,
            )
        )

    print()

    print(
        "LOCAL ROOT COUNT PRODUCT =",
        local_root_product
    )

    print(
        "ALL TRUE LOCAL ROOTS PRESENT =",
        all(local_checks)
    )

    print()

    # ------------------------------------------------------------
    # Combine all roots by CRT.
    # ------------------------------------------------------------

    crt_start = time.perf_counter()

    crt_classes = combine_residue_classes(
        classes,
        1
    )

    crt_elapsed = (
        time.perf_counter()
        - crt_start
    )

    print(
        "CRT S RESIDUE CLASSES =",
        len(crt_classes)
    )

    print(
        "CRT RUNTIME =",
        f"{crt_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # Check true S.
    # ------------------------------------------------------------

    true_crt_residue = (
        S_true % modulus
    )

    true_class_present = any(
        residue == true_crt_residue
        for residue, m in crt_classes
    )

    print(
        "TRUE S MOD M =",
        true_crt_residue
    )

    print(
        "TRUE CRT CLASS PRESENT =",
        true_class_present
    )

    print()

    # ------------------------------------------------------------
    # Balanced S interval.
    # ------------------------------------------------------------

    S_lo, S_hi = balanced_sum_interval(n)

    print(
        "BALANCED S INTERVAL"
    )

    print(
        "S lower =",
        S_lo
    )

    print(
        "S upper =",
        S_hi
    )

    print(
        "interval width =",
        S_hi - S_lo + 1
    )

    print()

    # ------------------------------------------------------------
    # Lift CRT classes.
    # ------------------------------------------------------------

    total_lift_candidates = 0

    true_lift_survives = False

    factor_hits = set()

    lifted_by_class = []

    for residue, m in crt_classes:

        candidates = lift_residue_class(
            residue,
            m,
            S_lo,
            S_hi
        )

        if candidates:

            lifted_by_class.append(
                (
                    residue,
                    m,
                    candidates,
                )
            )

        total_lift_candidates += len(
            candidates
        )

        for S in candidates:

            if S == S_true:
                true_lift_survives = True

            # For odd p,q, S must be even.
            #
            # This is safe for our generated semiprimes because
            # p and q are odd.
            if S % 2 != 0:
                continue

            hit = factor_from_sum(
                n,
                S
            )

            if hit is not None:

                factor_hits.add(hit)

    print(
        "LIFTING"
    )

    print()

    print(
        "CRT classes with >=1 lift =",
        len(lifted_by_class)
    )

    print(
        "TOTAL S CANDIDATES =",
        total_lift_candidates
    )

    print(
        "TRUE S SURVIVES =",
        true_lift_survives
    )

    print()

    # ------------------------------------------------------------
    # Show lifted candidates.
    # ------------------------------------------------------------

    print(
        "FIRST CRT-LIFT CLASSES"
    )

    shown = 0

    for residue, m, candidates in lifted_by_class:

        print(
            f"residue={residue} "
            f"mod={m} "
            f"candidates={candidates[:10]}"
            + (
                " ..."
                if len(candidates) > 10
                else ""
            )
        )

        shown += 1

        if shown >= 10:
            break

    print()

    # ------------------------------------------------------------
    # Exact factor result.
    # ------------------------------------------------------------

    print(
        "EXACT FACTORIZATION"
    )

    print()

    print(
        "factor hits =",
        len(factor_hits)
    )

    for hit in sorted(factor_hits):

        print(
            "  ",
            hit
        )

    print()

    # ------------------------------------------------------------
    # Complexity comparison.
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
        "S interval width =",
        S_hi - S_lo + 1
    )

    print(
        "CRT classes =",
        len(crt_classes)
    )

    print(
        "S candidates after CRT =",
        total_lift_candidates
    )

    if total_lift_candidates > 0:

        print(
            "sqrt(n) / CRT-S-candidates =",
            f"{sqrt_n / total_lift_candidates:.6f}"
        )

    else:

        print(
            "sqrt(n) / CRT-S-candidates = INF"
        )

    print()

    # ------------------------------------------------------------
    # Modulus threshold.
    # ------------------------------------------------------------

    print(
        "MODULUS THRESHOLD"
    )

    print()

    print(
        "M =",
        modulus
    )

    print(
        "M / sqrt(n) =",
        f"{modulus / sqrt_n:.6f}"
    )

    if modulus > S_hi - S_lo + 1:

        print(
            "M exceeds the entire balanced S interval."
        )

        print(
            "Each CRT class can therefore contribute at most "
            "one S candidate."
        )

    else:

        print(
            "M does not exceed the balanced S interval."
        )

    print()

    # ------------------------------------------------------------
    # Direct verification of discriminant congruence.
    # ------------------------------------------------------------

    print(
        "TRUE DISCRIMINANT CHECK"
    )

    D_true = (
        S_true * S_true
        - 4 * n
    )

    print(
        "D = S^2 - 4n =",
        D_true
    )

    print(
        "D is perfect square =",
        math.isqrt(D_true) ** 2 == D_true
    )

    for r in RADICES:

        print(
            f"D mod {r} =",
            D_true % r
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
print("FINISHED EXPERIMENT 148")
print("=" * 72)
