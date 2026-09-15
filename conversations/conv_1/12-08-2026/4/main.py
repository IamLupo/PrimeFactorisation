#!/usr/bin/env python3
"""
KAPPA N-ONLY AUXILIARY STRUCTURE / U-COLLAPSE EXPERIMENT

Goal
----
Test whether auxiliary-prime identities can constrain

    u = p + q

from n alone, without supplying:

    p, q, A, K, or any oracle value.

Core definitions
----------------
    F(x) = x^2 - x + 1

For n = p*q and u = p+q:

    A = F(p)F(q)
      = u^2 -(n+1)u+n^2-n+1

The previous experiments showed:

    exact A -> u -> factors

but did NOT show:

    n -> A.

This experiment deliberately removes A and K from the
constraint-generation stage.

It searches for:
    1. n-only congruence structure
    2. auxiliary-prime induced residue structure
    3. CRT accumulation
    4. collapse of possible u values
    5. collisions between different factor pairs

The script prints results directly and does not create CSV files.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable


# ============================================================================
# CONFIGURATION
# ============================================================================

AUX_PRIMES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

MODULI = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47
]

TARGET_BITS = [20, 30, 40, 50]
TARGETS_PER_SIZE = 5

RANDOM_SEED = 19082026

# For the direct u-space experiment, enumerating all u is impossible
# at 50 bits. Instead we test residue classes modulo accumulated M.
MAX_CRT_BITS = 30

# Small semiprime collision search.
COLLISION_LIMIT = 5000


# ============================================================================
# BASIC POLYNOMIALS
# ============================================================================

def F(x: int) -> int:
    return x * x - x + 1


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    if n % 3 == 0:
        return n == 3

    d = 5
    step = 2

    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step

    return True


def primes_up_to(limit: int) -> list[int]:
    return [n for n in range(2, limit + 1) if is_prime(n)]


# ============================================================================
# SEMIPRIME GENERATION
# ============================================================================

def random_prime(bits: int) -> int:
    lo = 1 << (bits - 1)
    hi = (1 << bits) - 1

    while True:
        x = random.randint(lo, hi)
        x |= 1

        if is_prime(x):
            return x


@dataclass
class Target:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def u(self) -> int:
        return self.p + self.q

    @property
    def A(self) -> int:
        return F(self.p) * F(self.q)


def make_target(n_bits: int) -> Target:
    """
    Generate p,q whose product is approximately n_bits.
    """
    p_bits = n_bits // 2
    q_bits = n_bits - p_bits

    while True:
        p = random_prime(p_bits)
        q = random_prime(q_bits)

        if p == q:
            continue

        target = Target(min(p, q), max(p, q))

        if target.n.bit_length() in (n_bits - 1, n_bits):
            return target


# ============================================================================
# AUXILIARY PRODUCT STRUCTURE
# ============================================================================

def auxiliary_product(rs: Iterable[int]) -> int:
    result = 1

    for r in rs:
        result *= F(r)

    return result


def print_auxiliary_sequences() -> None:
    sequences = {
        "first": AUX_PRIMES[:9],
        "odd_first": [p for p in AUX_PRIMES if p != 2][:8],
        "every_other": AUX_PRIMES[::2],
        "larger_first": [11, 13, 17, 19, 23, 29],
    }

    print("=" * 78)
    print("AUXILIARY PRIME STRUCTURE")
    print("=" * 78)

    for name, rs in sequences.items():
        B = auxiliary_product(rs)

        print(
            f"{name:<16} count={len(rs):2d} "
            f"B_bits={B.bit_length():2d} "
            f"B={B}"
        )
        print(f"  r = {tuple(rs)}")

    print()


# ============================================================================
# THE u POLYNOMIAL
# ============================================================================

def polynomial_u_value(n: int, u: int) -> int:
    """
    This is A as a polynomial in n,u.

        A = u² -(n+1)u+n²-n+1
    """
    return (
        u * u
        - (n + 1) * u
        + n * n
        - n
        + 1
    )


def u_residues_from_polynomial(
    n: int,
    modulus: int,
    *,
    assumed_A_residue: int | None = None,
) -> list[int]:
    """
    Solve the polynomial congruence for u.

    Normally this requires A mod modulus.

    If A is not supplied, there is NO restriction.

    This function therefore makes the oracle boundary explicit.
    """
    if assumed_A_residue is None:
        return list(range(modulus))

    result = []

    for u in range(modulus):
        if polynomial_u_value(n, u) % modulus == assumed_A_residue:
            result.append(u)

    return result


# ============================================================================
# N-ONLY INVARIANTS
# ============================================================================

def n_only_residue_tests(n: int, modulus: int) -> dict[str, object]:
    """
    Search for identities involving n,u,F(r) that could potentially
    constrain u without knowing A.

    The crucial point is that u is unknown.

    We test whether the polynomial itself has structural factors or
    degeneracies modulo m.
    """

    values = []

    for u in range(modulus):
        a = polynomial_u_value(n, u) % modulus

        values.append(a)

    unique_values = len(set(values))

    collisions = modulus - unique_values

    return {
        "modulus": modulus,
        "unique_A_residues": unique_values,
        "collisions": collisions,
        "mapping": values,
    }


# ============================================================================
# RESULTANT / DIFFERENCE STRUCTURE
# ============================================================================

def difference_factor(n: int, u1: int, u2: int) -> int:
    """
    Difference of A(u1) and A(u2):

        A(u1)-A(u2)
        = (u1-u2)(u1+u2-(n+1))
    """
    return (
        polynomial_u_value(n, u1)
        - polynomial_u_value(n, u2)
    )


def verify_difference_identity(n: int, u1: int, u2: int) -> bool:
    lhs = difference_factor(n, u1, u2)

    rhs = (
        (u1 - u2)
        * (u1 + u2 - (n + 1))
    )

    return lhs == rhs


# ============================================================================
# DISCRIMINANT WITHOUT A
# ============================================================================

def discriminant_from_u(n: int, u: int) -> int:
    """
    D = (2u-(n+1))².
    """
    return (2 * u - (n + 1)) ** 2


def factor_pair_from_u(n: int, u: int) -> tuple[int, int] | None:
    """
    Recover p,q if u corresponds to an integer factor pair.
    """
    D = u * u - 4 * n

    if D < 0:
        return None

    root = math.isqrt(D)

    if root * root != D:
        return None

    if (u - root) % 2 != 0:
        return None

    p = (u - root) // 2
    q = (u + root) // 2

    if p * q != n:
        return None

    return (min(p, q), max(p, q))


# ============================================================================
# FACTOR PAIR ENUMERATION FOR SMALL TESTS
# ============================================================================

def factor_pairs(n: int) -> list[tuple[int, int]]:
    pairs = []

    d = 2

    while d * d <= n:
        if n % d == 0:
            pairs.append((d, n // d))
        d += 1

    return pairs


# ============================================================================
# COLLISION SEARCH
# ============================================================================

def collision_search(limit: int) -> None:
    print("=" * 78)
    print("FIXED-n FACTOR-PAIR / A COLLISION SEARCH")
    print("=" * 78)

    total_composite = 0
    multi_pair = 0
    collisions = []

    for n in range(4, limit + 1):
        pairs = factor_pairs(n)

        if not pairs:
            continue

        total_composite += 1

        if len(pairs) < 2:
            continue

        multi_pair += 1

        values = {}

        for p, q in pairs:
            a = F(p) * F(q)
            values.setdefault(a, []).append((p, q))

        if len(values) < len(pairs):
            collisions.append((n, pairs, values))

    print(f"range checked              : 4 .. {limit}")
    print(f"composite n                : {total_composite}")
    print(f"n with >=2 factor pairs    : {multi_pair}")
    print(f"A collisions               : {len(collisions)}")
    print()

    for n, pairs, values in collisions[:10]:
        print(f"n={n}")
        print(f"  pairs = {pairs}")
        print(f"  A map = {values}")
        print()

    if not collisions:
        print("No A collisions found.")
    else:
        print("Important: A can distinguish many factor pairs,")
        print("but A itself is still not known from n alone.")

    print()


# ============================================================================
# N-ONLY RANGE OF u
# ============================================================================

def u_range(n: int) -> tuple[int, int]:
    """
    For positive factors p,q:

        u=p+q >= 2sqrt(n)

    and trivially u <= n+1 for p=1,q=n.
    """
    lower = math.isqrt(4 * n)

    while lower * lower < 4 * n:
        lower += 1

    upper = n + 1

    return lower, upper


def print_u_geometry(target: Target) -> None:
    n = target.n
    u = target.u

    lo, hi = u_range(n)

    print("=" * 78)
    print("u GEOMETRY FROM n ALONE")
    print("=" * 78)

    print(f"n       = {n}")
    print(f"n bits  = {n.bit_length()}")
    print(f"true u  = {u}")
    print(f"u bits  = {u.bit_length()}")
    print()
    print(f"mathematical u range:")
    print(f"  {lo} <= u <= {hi}")
    print()
    print(
        "The raw interval contains approximately "
        f"{hi - lo + 1:,} integer u values."
    )
    print()


# ============================================================================
# MODULAR u-SPACE WITHOUT A
# ============================================================================

def modular_u_space_without_A(
    n: int,
    moduli: list[int],
) -> None:
    """
    Demonstrate the key limitation directly.

    Without A, every u residue is possible.

    We also inspect the polynomial map u -> A mod m and count
    how much information is lost by that map.
    """

    print("=" * 78)
    print("N-ONLY MODULAR u-SPACE TEST")
    print("=" * 78)

    print(
        "No A or K is supplied here.\n"
        "For each modulus we examine the map\n"
        "    u -> A(u) mod m\n"
        "and ask whether n alone selects one u residue."
    )
    print()

    print(
        f"{'m':>5} "
        f"{'u residues':>12} "
        f"{'A residues':>12} "
        f"{'compression':>14}"
    )
    print("-" * 50)

    for m in moduli:
        test = n_only_residue_tests(n, m)

        u_count = m
        a_count = test["unique_A_residues"]

        compression = u_count / a_count if a_count else 0

        print(
            f"{m:5d} "
            f"{u_count:12d} "
            f"{a_count:12d} "
            f"{compression:14.3f}"
        )

    print()
    print(
        "Interpretation: the polynomial can map several u residues to "
        "the same A residue, but without knowing which A residue is "
        "correct, n alone does not select among them."
    )
    print()


# ============================================================================
# AUXILIARY PRIME STRUCTURE
# ============================================================================

def auxiliary_structure_test(n: int) -> None:
    """
    Search for accidental relationships between n and F(r).

    We explicitly do NOT use p, q, u, A, or K.

    We record:
      gcd(n,F(r))
      gcd(n-1,F(r))
      gcd(n+1,F(r))
      gcd(n^2-1,F(r))

    These are legitimate n-only quantities.
    """

    print("=" * 78)
    print("N-ONLY AUXILIARY PRIME STRUCTURE")
    print("=" * 78)

    print(f"n = {n}")
    print()

    print(
        f"{'r':>4} "
        f"{'F(r)':>8} "
        f"{'gcd(n,F)':>10} "
        f"{'gcd(n-1,F)':>12} "
        f"{'gcd(n+1,F)':>12} "
        f"{'gcd(n²-1,F)':>14}"
    )
    print("-" * 70)

    for r in AUX_PRIMES:
        fr = F(r)

        print(
            f"{r:4d} "
            f"{fr:8d} "
            f"{math.gcd(n, fr):10d} "
            f"{math.gcd(n - 1, fr):12d} "
            f"{math.gcd(n + 1, fr):12d} "
            f"{math.gcd(n * n - 1, fr):14d}"
        )

    print()


# ============================================================================
# AUXILIARY EXTENSION ALGEBRA
# ============================================================================

def extension_identity_test() -> None:
    """
    Verify the exact extension law:

        K' - K = s(1-s)(1-K)

    and rewrite it using

        1-K = product F(x_i).

    The experiment searches for quantities where K cancels.
    """

    print("=" * 78)
    print("AUXILIARY EXTENSION / K-CANCELLATION TEST")
    print("=" * 78)

    print(
        "For an extension by s:"
    )
    print()
    print("    K' - K = s(1-s)(1-K)")
    print()
    print("Equivalently:")
    print()
    print("    (K'-K)/(1-K) = s(1-s)")
    print()
    print(
        "The right-hand side is known when s is known, but both "
        "K' and K remain oracle quantities."
    )
    print()

    for s in [2, 3, 5, 7, 11, 13]:
        lhs_coefficient = s * (1 - s)

        print(
            f"s={s:2d}  "
            f"coefficient s(1-s)={lhs_coefficient:5d}"
        )

    print()


# ============================================================================
# MULTI-TARGET TEST
# ============================================================================

def multi_target_test() -> list[Target]:
    print("=" * 78)
    print("MULTI-TARGET N-ONLY TEST")
    print("=" * 78)

    targets = []

    for bits in TARGET_BITS:
        for _ in range(TARGETS_PER_SIZE):
            target = make_target(bits)
            targets.append(target)

            print(
                f"n_bits={target.n.bit_length():2d} "
                f"n={target.n:<18d} "
                f"u_bits={target.u.bit_length():2d} "
                f"u={target.u}"
            )

    print()

    return targets


# ============================================================================
# IMPORTANT CONTROL: ORACLE VS N-ONLY
# ============================================================================

def oracle_control(target: Target) -> None:
    """
    Show the exact oracle route for comparison.

    This is intentionally separated from the n-only experiment.
    """

    print("=" * 78)
    print("CONTROL EXPERIMENT: ORACLE ROUTE")
    print("=" * 78)

    n = target.n
    A = target.A

    print(f"n = {n}")
    print(f"true A = {A}")
    print(f"true u = {target.u}")
    print()

    for r in AUX_PRIMES[:9]:
        B = F(r)

        # K = 1 - A*F(r)
        K = 1 - A * B

        recovered_A = (1 - K) // B

        print(
            f"r={r:2d} "
            f"F(r)={B:8d} "
            f"K_bits={abs(K).bit_length():3d} "
            f"A_recovered={recovered_A == A}"
        )

    print()


# ============================================================================
# DIFFERENCE IDENTITY SEARCH
# ============================================================================

def difference_identity_demo(target: Target) -> None:
    print("=" * 78)
    print("u-DIFFERENCE STRUCTURE")
    print("=" * 78)

    n = target.n
    true_u = target.u

    print(
        "For any u1,u2:"
    )
    print()
    print(
        "A(u1)-A(u2) = "
        "(u1-u2)(u1+u2-(n+1))"
    )
    print()

    tests = [
        (true_u, true_u + 1),
        (true_u, true_u + 2),
        (true_u, n + 1 - true_u),
    ]

    for u1, u2 in tests:
        lhs = difference_factor(n, u1, u2)
        rhs = (
            (u1 - u2)
            * (u1 + u2 - (n + 1))
        )

        print(
            f"u1={u1}  u2={u2}"
        )
        print(
            f"  difference = {lhs}"
        )
        print(
            f"  factored   = {rhs}"
        )
        print(
            f"  identity   = {lhs == rhs}"
        )
        print()

    print(
        "Notice the symmetry:"
    )
    print()
    print(
        "A(u) = A(n+1-u)."
    )
    print(
        "This is exactly why the discriminant produces two u candidates."
    )
    print()


# ============================================================================
# SEARCH FOR N-ONLY MODULAR ANNIHILATION
# ============================================================================

def search_n_only_annihilation(target: Target) -> None:
    """
    Try simple expressions involving n, r and u that might vanish
    modulo F(r), without using p,q,A,K.

    We search a deliberately modest expression family.

    If a relation survives for every target, it becomes interesting.
    """

    print("=" * 78)
    print("SEARCH: N-ONLY AUXILIARY ANNIHILATION")
    print("=" * 78)

    print(
        "Testing whether simple expressions in n,u,r vanish modulo F(r)."
    )
    print(
        "No A or K is used."
    )
    print()

    expressions = {
        "u-r": lambda n, u, r: u - r,
        "u+r": lambda n, u, r: u + r,
        "u-n": lambda n, u, r: u - n,
        "u-(n+1)": lambda n, u, r: u - (n + 1),
        "u^2-nu": lambda n, u, r: u * u - n * u,
        "u^2-(n+1)u": lambda n, u, r: u * u - (n + 1) * u,
        "n-u+r": lambda n, u, r: n - u + r,
        "n+u+r": lambda n, u, r: n + u + r,
    }

    for name, expr in expressions.items():
        hits = 0
        total = 0

        for r in AUX_PRIMES:
            m = F(r)

            value = expr(
                target.n,
                target.u,
                r,
            )

            total += 1

            if value % m == 0:
                hits += 1

        if hits:
            print(
                f"{name:<20} hits={hits:2d}/{total}"
            )

    print()
    print(
        "Only repeated/non-accidental identities across many targets "
        "would be interesting here."
    )
    print()


# ============================================================================
# CRT DEMONSTRATION
# ============================================================================

def crt_pair(a1: int, m1: int, a2: int, m2: int) -> tuple[int, int]:
    """
    Chinese remainder theorem for coprime moduli.
    """
    inv = pow(m1, -1, m2)

    t = ((a2 - a1) * inv) % m2

    x = a1 + m1 * t
    m = m1 * m2

    return x % m, m


def crt_without_oracle_demo(target: Target) -> None:
    """
    Demonstrate what CRT can and cannot do without A.

    Since A is not supplied, each modulus admits all possible u residues.
    Thus the combined CRT space remains the complete residue space.
    """

    print("=" * 78)
    print("CRT WITHOUT A/K ORACLE")
    print("=" * 78)

    print(
        "This is the crucial control."
    )
    print()
    print(
        "If A mod m is unknown, the congruence"
    )
    print()
    print(
        "u² -(n+1)u+n²-n+1 = A (mod m)"
    )
    print()
    print(
        "does not select a residue of u."
    )
    print()

    M = 1

    for m in MODULI:
        if math.gcd(M, m) != 1:
            continue

        new_M = M * m

        if new_M.bit_length() > MAX_CRT_BITS:
            break

        M = new_M

        print(
            f"adding m={m:2d} -> "
            f"CRT modulus bits={M.bit_length():2d} "
            f"possible u residues={M:,}"
        )

    print()

    print(
        f"Final modulus M={M:,} ({M.bit_length()} bits)"
    )
    print(
        "Without an A residue, every one of those residue classes "
        "remains possible."
    )
    print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    random.seed(RANDOM_SEED)

    print()
    print("=" * 78)
    print("KAPPA N-ONLY AUXILIARY / u-COLLAPSE EXPERIMENT")
    print("=" * 78)
    print()
    print("IMPORTANT:")
    print("  This experiment does NOT use A or K to generate constraints.")
    print("  A/K are shown only in separate control experiments.")
    print()
    print("Question:")
    print("  Can auxiliary-prime structure reduce the u-space from n alone?")
    print()

    # ------------------------------------------------------------------
    # 1. Auxiliary sequences
    # ------------------------------------------------------------------

    print_auxiliary_sequences()

    # ------------------------------------------------------------------
    # 2. Generate representative target
    # ------------------------------------------------------------------

    target = make_target(50)

    print("=" * 78)
    print("REPRESENTATIVE TARGET")
    print("=" * 78)

    print(f"p       = {target.p}")
    print(f"q       = {target.q}")
    print(f"n       = {target.n}")
    print(f"n bits  = {target.n.bit_length()}")
    print(f"u       = {target.u}")
    print(f"u bits  = {target.u.bit_length()}")
    print(f"A bits  = {target.A.bit_length()}")
    print()

    # ------------------------------------------------------------------
    # 3. u geometry
    # ------------------------------------------------------------------

    print_u_geometry(target)

    # ------------------------------------------------------------------
    # 4. N-only modular experiment
    # ------------------------------------------------------------------

    modular_u_space_without_A(
        target.n,
        MODULI,
    )

    # ------------------------------------------------------------------
    # 5. N-only auxiliary structure
    # ------------------------------------------------------------------

    auxiliary_structure_test(target.n)

    # ------------------------------------------------------------------
    # 6. Extension identity
    # ------------------------------------------------------------------

    extension_identity_test()

    # ------------------------------------------------------------------
    # 7. Difference structure
    # ------------------------------------------------------------------

    difference_identity_demo(target)

    # ------------------------------------------------------------------
    # 8. Search for simple n-only auxiliary relations
    # ------------------------------------------------------------------

    search_n_only_annihilation(target)

    # ------------------------------------------------------------------
    # 9. CRT control
    # ------------------------------------------------------------------

    crt_without_oracle_demo(target)

    # ------------------------------------------------------------------
    # 10. Oracle comparison
    # ------------------------------------------------------------------

    oracle_control(target)

    # ------------------------------------------------------------------
    # 11. Multi-target test
    # ------------------------------------------------------------------

    targets = multi_target_test()

    # ------------------------------------------------------------------
    # 12. Small collision search
    # ------------------------------------------------------------------

    collision_search(COLLISION_LIMIT)

    # ------------------------------------------------------------------
    # 13. Final report
    # ------------------------------------------------------------------

    print("=" * 78)
    print("FINAL REPORT")
    print("=" * 78)
    print()
    print("N-ONLY ROUTE")
    print("-" * 78)
    print()
    print("    n")
    print("     |")
    print("     +--> F(r) structure")
    print("     |")
    print("     +--> simple gcd/congruence tests")
    print("     |")
    print("     +--> polynomial A(u)")
    print("     |")
    print("     +--> possible u residues")
    print()
    print(
        "The experiment intentionally does not insert the hidden A "
        "residue."
    )
    print()
    print("ORACLE ROUTE")
    print("-" * 78)
    print()
    print("    n + A")
    print("       |")
    print("       v")
    print("    quadratic in u")
    print("       |")
    print("       v")
    print("    discriminant")
    print("       |")
    print("       v")
    print("    p,q")
    print()
    print(
        "The oracle route succeeds, but is not equivalent to an n-only "
        "factoring algorithm."
    )
    print()
    print("KEY OBSERVATION TO WATCH")
    print("-" * 78)
    print()
    print(
        "A genuinely interesting result would be a relation involving "
        "only n, auxiliary r-values and u that removes the unknown A/K "
        "quantity."
    )
    print()
    print(
        "In particular, look for repeated non-random congruences in the "
        "N-ONLY AUXILIARY ANNIHILATION section."
    )
    print()
    print(
        "If no such relation appears, the next experiment should attack "
        "the algebra of multiple auxiliary extensions directly rather "
        "than accumulating more CRT residues."
    )
    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()