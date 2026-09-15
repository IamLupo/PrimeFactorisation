# ================================================================
# START EXPERIMENT 147
# Cross-radix quotient elimination / gcd leak test
#
# Goal:
#
#   Test whether several exact quotient relations derived from
#   different radices can produce a non-trivial gcd containing p
#   or q, WITHOUT scanning candidate factors.
#
# For each radix r:
#
#     p = r*k_r + a_r
#
# and
#
#     n = p*q
#
# so
#
#     floor(n/r)
#       = k_r*q + floor(a_r*q/r).
#
# Define
#
#     Q_r = floor(n/r)
#
# and
#
#     T_r = r*Q_r - n.
#
# Then
#
#     T_r = a_r*q - r*floor(a_r*q/r).
#
# Hence T_r is a small remainder modulo q:
#
#     T_r = -(n mod r)  [up to the sign convention]
#
# but combining several radices may reveal a common q-divisible
# expression after eliminating the unknown quotient terms.
#
# This experiment constructs several natural eliminations and
# performs gcd tests.
#
# IMPORTANT:
#
#   A non-trivial gcd is interesting only if it equals p or q.
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

        witness_passed = False

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                witness_passed = True
                break

        if not witness_passed:
            return False

    return True


def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= (1 << (bits - 1))

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
# Exact radix state
# ================================================================

def radix_state(n, p, q, r):

    k = p // r
    a = p % r

    Q = n // r
    rem = n % r

    # Exact identity:
    #
    # n = r*k*q + a*q
    #
    # Q = k*q + floor(a*q/r)
    #
    carry = (a * q) // r

    assert Q == k * q + carry

    # r*Q - n
    #
    # = r*(kq + carry) - (rkq + aq)
    #
    # = r*carry - aq
    #
    T = r * Q - n

    assert T == r * carry - a * q

    return {
        "r": r,
        "k": k,
        "a": a,
        "Q": Q,
        "rem": rem,
        "carry": carry,
        "T": T,
    }


# ================================================================
# Print helper
# ================================================================

def print_state(s):

    print(
        f"r={s['r']:3d} "
        f"k={s['k']:>12d} "
        f"a={s['a']:3d} "
        f"Q={s['Q']:>18d} "
        f"rem={s['rem']:3d} "
        f"carry={s['carry']:>18d} "
        f"T={s['T']:>18d}"
    )


# ================================================================
# Pairwise elimination candidates
# ================================================================

def build_pairwise_expressions(states):

    expressions = []

    for i in range(len(states)):

        s1 = states[i]

        for j in range(i + 1, len(states)):

            s2 = states[j]

            r1 = s1["r"]
            r2 = s2["r"]

            Q1 = s1["Q"]
            Q2 = s2["Q"]

            T1 = s1["T"]
            T2 = s2["T"]

            # ----------------------------------------------------
            # Expression A
            #
            # r1*Q1 - r2*Q2
            #
            # Since:
            #
            # r*Q = n + T
            #
            # this is simply T1-T2.
            # ----------------------------------------------------

            A = (
                r1 * Q1
                - r2 * Q2
            )

            # ----------------------------------------------------
            # Expression B
            #
            # r2*T1 - r1*T2
            # ----------------------------------------------------

            B = (
                r2 * T1
                - r1 * T2
            )

            # ----------------------------------------------------
            # Expression C
            #
            # Q1*r2 - Q2*r1
            #
            # algebraically same scale as A, but kept as a
            # separate control.
            # ----------------------------------------------------

            C = (
                Q1 * r2
                - Q2 * r1
            )

            # ----------------------------------------------------
            # Expression D
            #
            # T1*r2 - T2*r1
            # ----------------------------------------------------

            D = (
                T1 * r2
                - T2 * r1
            )

            expressions.append(
                {
                    "i": i,
                    "j": j,
                    "r1": r1,
                    "r2": r2,
                    "A": A,
                    "B": B,
                    "C": C,
                    "D": D,
                }
            )

    return expressions


# ================================================================
# Triple elimination
# ================================================================

def build_triple_expressions(states):

    expressions = []

    for i in range(len(states)):

        for j in range(i + 1, len(states)):

            for k in range(j + 1, len(states)):

                s1 = states[i]
                s2 = states[j]
                s3 = states[k]

                r1 = s1["r"]
                r2 = s2["r"]
                r3 = s3["r"]

                Q1 = s1["Q"]
                Q2 = s2["Q"]
                Q3 = s3["Q"]

                T1 = s1["T"]
                T2 = s2["T"]
                T3 = s3["T"]

                # ------------------------------------------------
                # Eliminate the explicit n contribution.
                #
                # Since:
                #
                #   r_i Q_i = n + T_i
                #
                # use:
                #
                #   (r2-r3)T1
                # + (r3-r1)T2
                # + (r1-r2)T3
                #
                # ------------------------------------------------

                X = (
                    (r2 - r3) * T1
                    + (r3 - r1) * T2
                    + (r1 - r2) * T3
                )

                # ------------------------------------------------
                # Same construction directly using Q.
                # ------------------------------------------------

                Y = (
                    (r2 - r3) * (r1 * Q1)
                    + (r3 - r1) * (r2 * Q2)
                    + (r1 - r2) * (r3 * Q3)
                )

                expressions.append(
                    {
                        "i": i,
                        "j": j,
                        "k": k,
                        "X": X,
                        "Y": Y,
                    }
                )

    return expressions


# ================================================================
# GCD classification
# ================================================================

def classify_gcd(g, p, q, n):

    if g == 0:
        return "ZERO"

    if g == 1:
        return "ONE"

    if g == n:
        return "N"

    if g == p:
        return "P"

    if g == q:
        return "Q"

    if n % g == 0:
        return "NONTRIVIAL_DIVISOR"

    return "OTHER"


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 147")
print("Cross-radix quotient elimination / gcd leak test")
print("=" * 72)
print()

print("RADICES =", RADICES)
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

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    print()

    # ------------------------------------------------------------
    # Build exact radix states.
    # ------------------------------------------------------------

    states = []

    for r in RADICES:

        s = radix_state(
            n,
            p,
            q,
            r
        )

        states.append(s)

    print("RADIX STATES")
    print()

    for s in states:
        print_state(s)

    print()

    # ------------------------------------------------------------
    # Verify fundamental identity.
    # ------------------------------------------------------------

    print("FUNDAMENTAL IDENTITIES")
    print()

    identity_pass = True

    for s in states:

        r = s["r"]
        k = s["k"]
        a = s["a"]
        Q = s["Q"]
        carry = s["carry"]

        lhs = Q
        rhs = k * q + carry

        if lhs != rhs:
            identity_pass = False

        lhs2 = (
            r * Q - n
        )

        rhs2 = (
            r * carry
            - a * q
        )

        if lhs2 != rhs2:
            identity_pass = False

    print(
        "Q_r = k_r*q + floor(a_r*q/r):",
        "PASS" if identity_pass else "FAIL"
    )

    print(
        "r*Q_r - n = r*carry - a*q:",
        "PASS" if identity_pass else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # Pairwise expressions.
    # ------------------------------------------------------------

    pair_exprs = build_pairwise_expressions(
        states
    )

    print(
        "PAIRWISE GCD TESTS"
    )

    print()

    pair_hits = {
        "P": 0,
        "Q": 0,
        "NONTRIVIAL_DIVISOR": 0,
        "N": 0,
        "OTHER": 0,
        "ONE": 0,
        "ZERO": 0,
    }

    interesting_pairs = []

    for e in pair_exprs:

        A = e["A"]
        B = e["B"]
        C = e["C"]
        D = e["D"]

        gcds = [
            ("A,B", math.gcd(abs(A), abs(B))),
            ("A,C", math.gcd(abs(A), abs(C))),
            ("A,D", math.gcd(abs(A), abs(D))),
            ("B,C", math.gcd(abs(B), abs(C))),
            ("B,D", math.gcd(abs(B), abs(D))),
            ("C,D", math.gcd(abs(C), abs(D))),
        ]

        for name, g in gcds:

            cls = classify_gcd(
                g,
                p,
                q,
                n
            )

            pair_hits[cls] += 1

            if cls in (
                "P",
                "Q",
                "NONTRIVIAL_DIVISOR",
            ):

                interesting_pairs.append(
                    (
                        e["r1"],
                        e["r2"],
                        name,
                        g,
                        cls,
                    )
                )

    print(
        "pair gcd classified as P =",
        pair_hits["P"]
    )

    print(
        "pair gcd classified as Q =",
        pair_hits["Q"]
    )

    print(
        "pair gcd classified as nontrivial divisor =",
        pair_hits["NONTRIVIAL_DIVISOR"]
    )

    print(
        "pair gcd classified as N =",
        pair_hits["N"]
    )

    print()

    # ------------------------------------------------------------
    # Triple expressions.
    # ------------------------------------------------------------

    triple_exprs = build_triple_expressions(
        states
    )

    print(
        "TRIPLE ELIMINATION GCD TESTS"
    )

    print()

    triple_hits = {
        "P": 0,
        "Q": 0,
        "NONTRIVIAL_DIVISOR": 0,
        "N": 0,
        "OTHER": 0,
        "ONE": 0,
        "ZERO": 0,
    }

    interesting_triples = []

    for e in triple_exprs:

        X = e["X"]
        Y = e["Y"]

        g = math.gcd(
            abs(X),
            abs(Y)
        )

        cls = classify_gcd(
            g,
            p,
            q,
            n
        )

        triple_hits[cls] += 1

        if cls in (
            "P",
            "Q",
            "NONTRIVIAL_DIVISOR",
        ):

            interesting_triples.append(
                (
                    states[e["i"]]["r"],
                    states[e["j"]]["r"],
                    states[e["k"]]["r"],
                    g,
                    cls,
                )
            )

    print(
        "triple gcd = P:",
        triple_hits["P"]
    )

    print(
        "triple gcd = Q:",
        triple_hits["Q"]
    )

    print(
        "triple gcd = nontrivial divisor:",
        triple_hits["NONTRIVIAL_DIVISOR"]
    )

    print(
        "triple gcd = N:",
        triple_hits["N"]
    )

    print()

    # ------------------------------------------------------------
    # Explicit factor divisibility test.
    #
    # Instead of only looking at gcd(A,B), test whether p or q
    # divides each expression individually.
    #
    # This identifies whether the algebra naturally produces
    # factor-multiples even when two expressions have a larger
    # gcd.
    # ------------------------------------------------------------

    print("FACTOR DIVISIBILITY MATRIX")
    print()

    for s in states:

        r = s["r"]

        Qr = s["Q"]
        Tr = s["T"]

        expressions = {
            "Q": Qr,
            "rQ": r * Qr,
            "T": Tr,
            "rT": r * Tr,
            "Q-n": Qr - n,
            "rQ-n": r * Qr - n,
        }

        print(
            f"r={r}"
        )

        for name, value in expressions.items():

            div_p = (
                value % p == 0
            )

            div_q = (
                value % q == 0
            )

            print(
                f"  {name:6s}: "
                f"div_p={div_p} "
                f"div_q={div_q} "
                f"value={value}"
            )

        print()

    # ------------------------------------------------------------
    # Show any interesting hits.
    # ------------------------------------------------------------

    if interesting_pairs:

        print(
            "INTERESTING PAIRWISE HITS"
        )

        for hit in interesting_pairs[:20]:

            print(
                "  ",
                hit
            )

        print()

    else:

        print(
            "INTERESTING PAIRWISE HITS: NONE"
        )

        print()

    if interesting_triples:

        print(
            "INTERESTING TRIPLE HITS"
        )

        for hit in interesting_triples[:20]:

            print(
                "  ",
                hit
            )

        print()

    else:

        print(
            "INTERESTING TRIPLE HITS: NONE"
        )

        print()

    # ------------------------------------------------------------
    # Magnitude diagnostics.
    # ------------------------------------------------------------

    print("MAGNITUDE DIAGNOSTICS")
    print()

    for s in states:

        print(
            f"r={s['r']:3d}: "
            f"|Q|={abs(s['Q'])} "
            f"|T|={abs(s['T'])} "
            f"|q|={q}"
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
print("FINISHED EXPERIMENT 147")
print("=" * 72)
