#!/usr/bin/env python3

from math import gcd, isqrt
from itertools import combinations, product
import random


# ============================================================================
# KAPPA MULTI-AUXILIARY ELIMINATION / U-COLLAPSE SEARCH
# ============================================================================

F = lambda x: x * x - x + 1


# ----------------------------------------------------------------------------
# Prime generation
# ----------------------------------------------------------------------------

def primes_up_to(limit):
    sieve = [True] * (limit + 1)
    sieve[0:2] = [False, False]

    for i in range(2, isqrt(limit) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False

    return [i for i in range(2, limit + 1) if sieve[i]]


PRIMES = primes_up_to(1000)


# ----------------------------------------------------------------------------
# Basic algebra
# ----------------------------------------------------------------------------

def A_from_pq(p, q):
    return F(p) * F(q)


def A_from_nu(n, u):
    return (
        u * u
        - (n + 1) * u
        + n * n
        - n
        + 1
    )


def K_from_A_r(A, r):
    return 1 - A * F(r)


def K_sequence(A, rs):
    return [K_from_A_r(A, r) for r in rs]


def factor_pairs(n):
    out = []

    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            out.append((d, n // d))

    return out


# ----------------------------------------------------------------------------
# Semiprime generation
# ----------------------------------------------------------------------------

def random_prime(bits):
    candidates = [p for p in PRIMES if p.bit_length() >= bits]
    if not candidates:
        raise ValueError("Prime table too small")
    return random.choice(candidates)


def primes_near_bits(bits):
    lo = 1 << (bits - 1)
    hi = 1 << bits

    return [
        p for p in PRIMES
        if lo <= p < hi
    ]


def make_semiprime(total_bits):
    # Generate factors of roughly half the requested size.
    fb = total_bits // 2

    candidates = primes_near_bits(fb)

    if len(candidates) < 10:
        # Fall back to a larger generated prime list.
        raise RuntimeError(
            "Increase prime-generation limit for this target size."
        )

    p = random.choice(candidates)
    q = random.choice(candidates)

    while q == p:
        q = random.choice(candidates)

    n = p * q

    return p, q, n


# ----------------------------------------------------------------------------
# Auxiliary sequences
# ----------------------------------------------------------------------------

def build_auxiliary_sequences():
    first = tuple(PRIMES[:9])

    odd_first = tuple(
        p for p in PRIMES
        if p % 2 == 1
    )[:8]

    every_other = tuple(PRIMES[::2][:8])

    larger_first = tuple(
        p for p in PRIMES
        if p >= 11
    )[:6]

    return {
        "first": first,
        "odd_first": odd_first,
        "every_other": every_other,
        "larger_first": larger_first,
    }


# ----------------------------------------------------------------------------
# Utility formatting
# ----------------------------------------------------------------------------

def bits(x):
    x = abs(x)
    return 1 if x == 0 else x.bit_length()


def line(char="-", n=78):
    print(char * n)


# ============================================================================
# 1. AUXILIARY ALGEBRA
# ============================================================================

def test_auxiliary_algebra(A, rs):
    print()
    print("=" * 78)
    print("1. MULTI-AUXILIARY ALGEBRA")
    print("=" * 78)

    Ks = K_sequence(A, rs)

    print()
    print("Basic relation:")
    print()
    print("    K_r = 1 - A F(r)")
    print()

    print("r      F(r)              K_r")
    line()

    for r, k in zip(rs, Ks):
        print(f"{r:2d} {F(r):12d} {k:30d}")

    # Pairwise elimination.
    print()
    print("PAIRWISE ELIMINATION")
    line()

    for r1, r2 in combinations(rs[:5], 2):
        k1 = K_from_A_r(A, r1)
        k2 = K_from_A_r(A, r2)

        lhs = k1 - k2
        rhs = -A * (F(r1) - F(r2))

        print(
            f"r={r1:2d}, s={r2:2d}: "
            f"K_r-K_s={lhs}, "
            f"identity={lhs == rhs}"
        )


# ============================================================================
# 2. THREE-AUXILIARY ELIMINATION
# ============================================================================

def test_three_auxiliary_elimination(A, rs):
    print()
    print("=" * 78)
    print("2. THREE-AUXILIARY A-ELIMINATION")
    print("=" * 78)

    print()
    print(
        "We search for relations involving three K-values where A cancels."
    )

    print()
    print(
        "For r,s,t the expected identity is:"
    )
    print()
    print(
        " (K_r-K_s)(F(t)-F(r))"
    )
    print(
        " -(K_r-K_t)(F(s)-F(r)) = 0"
    )

    print()
    print("r,s,t      residual")
    line()

    discoveries = 0

    for r, s, t in combinations(rs, 3):

        kr = K_from_A_r(A, r)
        ks = K_from_A_r(A, s)
        kt = K_from_A_r(A, t)

        residual = (
            (kr - ks) * (F(t) - F(r))
            - (kr - kt) * (F(s) - F(r))
        )

        if residual == 0:
            discoveries += 1

        if discoveries <= 12:
            print(
                f"{r:2d},{s:2d},{t:2d}"
                f"{'':4s}{residual}"
            )

    print()
    print(f"Exact cancellations found: {discoveries}")
    print()
    print(
        "These are expected eliminations of A."
    )
    print(
        "The important question is whether anything involving n or u "
        "survives the elimination."
    )


# ============================================================================
# 3. SEARCH FOR INTEGER LINEAR ELIMINATION
# ============================================================================

def search_linear_K_relations(A, rs):
    print()
    print("=" * 78)
    print("3. LINEAR K-RELATION SEARCH")
    print("=" * 78)

    print()
    print(
        "We search small integer coefficient combinations"
    )
    print(
        "    c1*K_r1 + c2*K_r2 + ... + c0"
    )
    print(
        "whose A-dependent term cancels."
    )

    Ks = {r: K_from_A_r(A, r) for r in rs}
    Fs = {r: F(r) for r in rs}

    found = []

    # Two-term relations.
    for r, s in combinations(rs, 2):

        # F(s) * K_r - F(r) * K_s
        lhs = Fs[s] * Ks[r] - Fs[r] * Ks[s]
        rhs = Fs[s] - Fs[r]

        if lhs == rhs:
            found.append(
                (
                    f"F({s})K_{r} - F({r})K_{s}",
                    rhs
                )
            )

    print()
    print("Exact A-free relations:")
    line()

    for expression, value in found[:20]:
        print(
            f"{expression:<35} = {value}"
        )

    print()
    print(f"Relations found: {len(found)}")

    print()
    print("IMPORTANT:")
    print(
        "Every relation above is independent of A,"
    )
    print(
        "but also independent of p,q once the auxiliary r-values are fixed."
    )
    print(
        "Therefore it is not automatically factorization information."
    )


# ============================================================================
# 4. SEARCH FOR n/u RESIDUALS
# ============================================================================

def search_n_u_residuals(p, q, rs):
    n = p * q
    u = p + q
    A = A_from_pq(p, q)

    print()
    print("=" * 78)
    print("4. SEARCH FOR RESIDUAL n/u INFORMATION")
    print("=" * 78)

    print()
    print(f"p={p}")
    print(f"q={q}")
    print(f"n={n}")
    print(f"u={u}")
    print(f"A={A}")

    print()
    print(
        "Testing whether auxiliary F(r) values interact with n and u "
        "through simple exact identities."
    )

    tests = []

    for r in rs:
        fr = F(r)

        tests.append(
            (
                f"u-r mod F(r)",
                (u - r) % fr
            )
        )

        tests.append(
            (
                f"u-(n+1) mod F(r)",
                (u - (n + 1)) % fr
            )
        )

        tests.append(
            (
                f"u^2-nu mod F(r)",
                (u * u - n * u) % fr
            )
        )

        tests.append(
            (
                f"n mod F(r)",
                n % fr
            )
        )

        tests.append(
            (
                f"n+1 mod F(r)",
                (n + 1) % fr
            )
        )

    print()
    print("r     F(r)      expression                  residue")
    line()

    for r in rs:
        fr = F(r)

        vals = [
            ("u-r", (u - r) % fr),
            ("u-(n+1)", (u - (n + 1)) % fr),
            ("u²-nu", (u * u - n * u) % fr),
            ("n", n % fr),
            ("n+1", (n + 1) % fr),
        ]

        print(f"{r:2d} {fr:8d}")

        for name, residue in vals:
            print(
                f"      {name:<20s} {residue}"
            )


# ============================================================================
# 5. THE IMPORTANT TEST:
#    CAN MULTIPLE AUXILIARIES REDUCE u WITHOUT A?
# ============================================================================

def u_space_test(p, q, rs):
    n = p * q
    u_true = p + q

    print()
    print("=" * 78)
    print("5. MULTI-AUXILIARY u-SPACE COLLAPSE")
    print("=" * 78)

    print()
    print(
        "NO A OR K VALUES ARE USED TO SELECT u."
    )

    print()
    print(f"n       = {n}")
    print(f"true u  = {u_true}")
    print()

    # Mathematical bounds for semiprime p,q >= 2.
    lo = 2 + (n // 2)
    hi = n + 1

    print("Raw mathematical u interval:")
    print(f"    {lo} <= u <= {hi}")

    # For each F(r), inspect gcd interactions with n.
    print()
    print("AUXILIARY MODULAR GEOMETRY")
    line()

    print(
        "r    m=F(r)    gcd(n,m)   gcd(n+1,m)   "
        "possible u residues"
    )
    line()

    for r in rs:
        m = F(r)

        # Number of possible u residues modulo m.
        residues = list(range(m))

        # n alone cannot determine A, so every residue remains possible.
        possible = len(residues)

        print(
            f"{r:2d} {m:10d} "
            f"{gcd(n,m):10d} "
            f"{gcd(n + 1,m):13d} "
            f"{possible:18d}"
        )

    print()
    print(
        "Now combine auxiliary moduli."
    )

    M = 1

    for r in rs:
        m = F(r)

        # Only multiply if coprime, otherwise retain lcm-like growth.
        M = M // gcd(M, m) * m

        possible_classes = M

        print(
            f"r={r:2d} "
            f"M_bits={bits(M):3d} "
            f"M={M} "
            f"possible u classes={possible_classes}"
        )

        # Stop before enormous output.
        if M.bit_length() > 80:
            break

    print()
    print(
        "RESULT:"
    )
    print(
        "Without an independent residue for A, the auxiliary moduli"
    )
    print(
        "do not select one of these u classes."
    )


# ============================================================================
# 6. CONTROL: PUT A BACK IN
# ============================================================================

def oracle_u_collapse(p, q, rs):
    n = p * q
    u_true = p + q
    A = A_from_pq(p, q)

    print()
    print("=" * 78)
    print("6. CONTROL: ORACLE A -> u COLLAPSE")
    print("=" * 78)

    print()
    print(
        "This section intentionally supplies A."
    )

    print()
    print(f"n       = {n}")
    print(f"A       = {A}")
    print(f"true u  = {u_true}")

    print()
    print("r     modulus F(r)    surviving u residues")
    line()

    for r in rs:

        m = F(r)
        target = A % m

        survivors = []

        for residue in range(m):
            value = A_from_nu(n % m, residue) % m

            if value == target:
                survivors.append(residue)

        print(
            f"{r:2d} {m:15d} "
            f"{str(survivors):>30s}"
        )

    print()
    print(
        "This is the contrast we need:"
    )
    print()
    print(
        "    n alone       -> many u residues"
    )
    print(
        "    n + A mod m   -> usually 1 or 2 u residues"
    )


# ============================================================================
# 7. FACTOR-PAIR COLLISION TEST
# ============================================================================

def collision_test(limit=5000):
    print()
    print("=" * 78)
    print("7. FIXED-n MULTI-AUXILIARY COLLISION TEST")
    print("=" * 78)

    collisions = []

    for n in range(4, limit + 1):

        pairs = factor_pairs(n)

        if len(pairs) < 2:
            continue

        signatures = {}

        for p, q in pairs:

            A = A_from_pq(p, q)

            # Full auxiliary signature.
            signature = tuple(
                K_from_A_r(A, r)
                for r in PRIMES[:6]
            )

            signatures.setdefault(signature, []).append((p, q))

        if len(signatures) < len(pairs):
            collisions.append(
                (n, pairs, signatures)
            )

    print()
    print(f"range checked      : 4 .. {limit}")
    print(f"collisions found   : {len(collisions)}")

    for n, pairs, signatures in collisions[:10]:

        print()
        print(f"n={n}")
        print(f"pairs={pairs}")

        for sig, grouped in signatures.items():
            if len(grouped) > 1:
                print(
                    f"  COLLISION: {grouped}"
                )

    if not collisions:
        print()
        print(
            "No multi-auxiliary collisions found."
        )

    print()
    print(
        "This measures how much information the full K-vector"
    )
    print(
        "would contain IF it were available."
    )


# ============================================================================
# 8. SEARCH FOR A POSSIBLE NEW INVARIANT
# ============================================================================

def invariant_search(p, q, rs):
    n = p * q
    u = p + q
    A = A_from_pq(p, q)

    print()
    print("=" * 78)
    print("8. SEARCH FOR NEW n,u,AUXILIARY INVARIANTS")
    print("=" * 78)

    print()
    print(
        "We test simple expressions for repeated zero/residue behaviour."
    )

    expressions = {
        "u^2 - (n+1)u": lambda r: (
            u * u - (n + 1) * u
        ),

        "u^2 - nu": lambda r: (
            u * u - n * u
        ),

        "A": lambda r: A,

        "A-F(r)": lambda r: (
            A - F(r)
        ),

        "A*F(r)-1": lambda r: (
            A * F(r) - 1
        ),

        "A-F(r)*n": lambda r: (
            A - F(r) * n
        ),

        "A-F(r)*(n+1)": lambda r: (
            A - F(r) * (n + 1)
        ),
    }

    print()
    print(
        "expression                         gcd over all auxiliary values"
    )
    line()

    for name, fn in expressions.items():

        values = [
            fn(r)
            for r in rs
        ]

        g = 0

        for value in values:
            g = gcd(g, abs(value))

        print(
            f"{name:<35s} {g}"
        )

    print()
    print(
        "A surprisingly large or structured gcd recurring across"
    )
    print(
        "many unrelated targets would be worth investigating."
    )


# ============================================================================
# 9. MULTI-TARGET EXPERIMENT
# ============================================================================

def multi_target_test(targets, rs):
    print()
    print("=" * 78)
    print("9. MULTI-TARGET STRUCTURAL TEST")
    print("=" * 78)

    print()
    print(
        "target   n_bits   u_bits   n                  u"
    )
    line()

    for i, (p, q) in enumerate(targets, 1):

        n = p * q
        u = p + q

        print(
            f"{i:6d} "
            f"{n.bit_length():7d} "
            f"{u.bit_length():7d} "
            f"{n:<18d} "
            f"{u}"
        )

    print()
    print(
        "Searching for auxiliary expressions that vanish"
    )
    print(
        "for EVERY target."
    )

    candidates = {
        "u-r": lambda n, u, r: u - r,
        "u-(n+1)": lambda n, u, r: u - (n + 1),
        "u^2-nu": lambda n, u, r: u * u - n * u,
        "u^2-(n+1)u": lambda n, u, r: (
            u * u - (n + 1) * u
        ),
    }

    print()
    print(
        "expression          zero count / total"
    )
    line()

    total = len(targets) * len(rs)

    for name, fn in candidates.items():

        zeros = 0

        for p, q in targets:
            n = p * q
            u = p + q

            for r in rs:
                if fn(n, u, r) == 0:
                    zeros += 1

        print(
            f"{name:<20s} {zeros:5d} / {total}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    random.seed(20260812)

    print()
    print("=" * 78)
    print("KAPPA MULTI-AUXILIARY / ELIMINATION / u-COLLAPSE SEARCH")
    print("=" * 78)

    print()
    print("Goal:")
    print()
    print(
        "Search directly for a non-oracle relation generated by"
    )
    print(
        "multiple auxiliary primes."
    )

    print()
    print(
        "The central distinction is:"
    )
    print()
    print(
        "    n + K_r values -> A -> u"
    )
    print()
    print(
        "versus"
    )
    print()
    print(
        "    n alone -> auxiliary structure -> u"
    )

    # Auxiliary primes.
    sequences = build_auxiliary_sequences()

    print()
    print("=" * 78)
    print("AUXILIARY SEQUENCES")
    print("=" * 78)

    for name, rs in sequences.items():

        B = 1

        for r in rs:
            B *= F(r)

        print()
        print(
            f"{name:<16s} count={len(rs):2d} "
            f"B_bits={B.bit_length():2d} "
            f"B={B}"
        )
        print(
            f"  r={rs}"
        )

    rs = sequences["first"]

    # Representative target.
    p = 22612043
    q = 26706517

    n = p * q
    u = p + q
    A = A_from_pq(p, q)

    print()
    print("=" * 78)
    print("REPRESENTATIVE TARGET")
    print("=" * 78)

    print(f"p       = {p}")
    print(f"q       = {q}")
    print(f"n       = {n}")
    print(f"n bits  = {n.bit_length()}")
    print(f"u       = {u}")
    print(f"u bits  = {u.bit_length()}")
    print(f"A bits  = {A.bit_length()}")

    test_auxiliary_algebra(A, rs)

    test_three_auxiliary_elimination(A, rs)

    search_linear_K_relations(A, rs)

    search_n_u_residuals(p, q, rs)

    u_space_test(p, q, rs)

    oracle_u_collapse(p, q, rs)

    invariant_search(p, q, rs)

    collision_test(5000)

    # Generate additional targets.
    targets = []

    for bits_target in [20, 30, 40, 50]:

        for _ in range(3):

            try:
                pp, qq, nn = make_semiprime(bits_target)

                targets.append((pp, qq))

            except RuntimeError:
                break

    multi_target_test(targets, rs)

    print()
    print("=" * 78)
    print("FINAL OUTPUT")
    print("=" * 78)

    print()
    print(
        "The experiment has tested four possibilities:"
    )

    print()
    print("1. Pairwise K elimination")
    print("2. Three-auxiliary A elimination")
    print("3. n/u modular interaction")
    print("4. Multi-target invariant behaviour")

    print()
    print(
        "The strongest result to look for is a relation that:"
    )

    print()
    print("    • uses multiple auxiliary primes")
    print("    • eliminates A")
    print("    • eliminates K")
    print("    • depends on n and u")
    print("    • and is NOT merely the identity A=A")
    print()

    print(
        "If such a relation appears repeatedly across targets,"
    )
    print(
        "the next step is to attempt algebraic elimination of u"
    )
    print(
        "or derive a direct congruence for p or q."
    )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

