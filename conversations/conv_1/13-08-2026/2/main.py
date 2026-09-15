#!/usr/bin/env python3

"""
==============================================================================
ADAPTIVE CYCLOTOMIC MODULUS / N-ONLY FACTOR-RECOVERY EXPERIMENT
==============================================================================

Goal
----
Attack the next question directly:

    Given only n=pq,

    can we choose auxiliary r so that

        m = F(r) = r^2-r+1

    interacts unusually strongly with n,

    and can that interaction produce information about

        A = F(p)F(q)

    or

        u = p+q

    without using p,q,A,K,u when constructing the constraint?

Core identities
---------------

    F(x) = x^2-x+1

    A = F(p)F(q)

    n = pq
    u = p+q

    A = u^2-(n+1)u+n^2-n+1

    A-F(n) = u(u-n-1)

    4A-3(n-1)^2 = (2u-(n+1))^2

The experiment DOES NOT use A, p, q, or u to create the
n-only modulus search.

It does use the true values separately as a control.

Main searches
-------------

1. Search r for unusually strong gcd relations between F(r) and
   n-1, n+1, n^2-1, n^2+n+1, etc.

2. For every n-only modulus m, determine what information about u
   follows from identities that are genuinely known from n.

3. Investigate whether divisibility conditions force A modulo m.

4. Search adaptive sequences of r and combine their moduli.

5. Measure the actual reduction in possible u residues.

6. Separately compare with the oracle route to see how much stronger
   the true A residue is.

No CSV files.
Everything is printed.
"""


from __future__ import annotations

from itertools import combinations
from math import gcd, isqrt
import random


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 500

AUX_LIMIT = 120

MODULUS_CANDIDATE_LIMIT = 30

TARGETS_PER_SIZE = 5

TARGET_BITS = [20, 30, 40, 50]

CRT_MODULUS_BIT_LIMIT = 60

RANDOM_SEED = 20260813


# =============================================================================
# BASIC ARITHMETIC
# =============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def primes_up_to(limit: int) -> list[int]:
    return [
        n
        for n in range(2, limit + 1)
        if is_prime(n)
    ]


PRIMES = primes_up_to(PRIME_LIMIT)


def F(x: int) -> int:
    return x * x - x + 1


def A_true(p: int, q: int) -> int:
    return F(p) * F(q)


def A_from_n_u(n: int, u: int) -> int:
    return (
        u * u
        - (n + 1) * u
        + n * n
        - n
        + 1
    )


def bitlen(x: int) -> int:
    if x == 0:
        return 1

    return abs(x).bit_length()


def integer_sqrt_exact(x: int):
    if x < 0:
        return None

    r = isqrt(x)

    if r * r == x:
        return r

    return None


# =============================================================================
# PRIME / SEMIPRIME GENERATION
# =============================================================================

def primes_in_bit_range(bits: int) -> list[int]:
    lo = 1 << (bits - 1)
    hi = 1 << bits

    return [
        p
        for p in PRIMES
        if lo <= p < hi
    ]


def make_semiprime(total_bits: int):
    factor_bits = max(3, total_bits // 2)

    candidates = primes_in_bit_range(factor_bits)

    if len(candidates) < 2:
        raise RuntimeError(
            f"Not enough primes for {factor_bits}-bit factors. "
            f"Increase PRIME_LIMIT."
        )

    p = random.choice(candidates)
    q = random.choice(candidates)

    while q == p:
        q = random.choice(candidates)

    return p, q, p * q


# =============================================================================
# FORMATTING
# =============================================================================

def section(title: str):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def line():
    print("-" * 78)


# =============================================================================
# N-ONLY CANDIDATE EXPRESSIONS
# =============================================================================

def n_only_expressions(n: int) -> dict[str, int]:
    """
    Expressions involving n only.

    The point is not to claim these determine A, but to see which
    factors of F(r) are exposed by n-only arithmetic.
    """
    return {
        "n-1": n - 1,
        "n+1": n + 1,
        "n^2-1": n * n - 1,
        "n^2+n+1": n * n + n + 1,
        "n^2-n+1": n * n - n + 1,
        "n^3-1": n**3 - 1,
        "n^3+1": n**3 + 1,
        "n^3-n": n**3 - n,
    }


# =============================================================================
# FACTOR OF m EXPOSED BY n
# =============================================================================

def modulus_profile(n: int, r: int):
    m = F(r)

    exprs = n_only_expressions(n)

    gcds = {
        name: gcd(abs(value), m)
        for name, value in exprs.items()
    }

    return m, gcds


# =============================================================================
# CRT HELPERS
# =============================================================================

def crt_pair(a1: int, m1: int, a2: int, m2: int):
    """
    Generalized CRT.

    Returns (a,m) for a mod m satisfying both congruences,
    or None if inconsistent.
    """
    g = gcd(m1, m2)

    if (a2 - a1) % g != 0:
        return None

    m1g = m1 // g
    m2g = m2 // g

    # Solve:
    # a1 + m1*t = a2 mod m2
    rhs = (a2 - a1) // g

    inv = pow(m1g, -1, m2g)

    t = (rhs * inv) % m2g

    modulus = m1 * m2g
    residue = (a1 + m1 * t) % modulus

    return residue, modulus


# =============================================================================
# U RANGE
# =============================================================================

def possible_u_range(n: int):
    """
    Very broad mathematical range for positive factors p,q >= 2:

        p+q >= 2 + n/2
        p+q <= n+1

    This is not a useful factoring interval by itself, but gives a
    baseline for counting candidate integers.
    """
    lo = n // 2 + 2
    hi = n + 1

    return lo, hi


def count_interval_residues(lo: int, hi: int, modulus: int, residues: set[int]):
    """
    Count integers u in [lo,hi] whose residue mod modulus lies in residues.
    """
    if lo > hi or not residues:
        return 0

    total = 0

    for residue in residues:
        first = lo + ((residue - lo) % modulus)

        if first <= hi:
            total += 1 + (hi - first) // modulus

    return total


# =============================================================================
# IMPORTANT CONTROL: WHAT A MODULUS WOULD REVEAL IF A WERE KNOWN?
# =============================================================================

def oracle_u_residues(n: int, A: int, m: int):
    residues = set()

    a_mod = A % m
    n_mod = n % m

    for ur in range(m):
        candidate = (
            ur * ur
            - (n_mod + 1) * ur
            + n_mod * n_mod
            - n_mod
            + 1
        ) % m

        if candidate == a_mod:
            residues.add(ur)

    return residues


# =============================================================================
# SECTION 1
# REPRESENTATIVE TARGET
# =============================================================================

def representative_target():
    section("1. REPRESENTATIVE 50-BIT TARGET")

    # Keep this deterministic.
    p = 22612043
    q = 26706517

    n = p * q
    u = p + q
    A = A_true(p, q)

    print()
    print(f"p       = {p}")
    print(f"q       = {q}")
    print(f"n       = {n}")
    print(f"n bits  = {bitlen(n)}")
    print(f"u       = {u}")
    print(f"u bits  = {bitlen(u)}")
    print(f"A       = {A}")
    print(f"A bits  = {bitlen(A)}")

    return p, q, n, u, A


# =============================================================================
# SECTION 2
# VERIFY THE BASIC STRUCTURE
# =============================================================================

def verify_basic_structure(p, q, n, u, A):
    section("2. BASIC ALGEBRAIC CONTROL")

    checks = {
        "A = F(p)F(q)": (
            A == F(p) * F(q)
        ),

        "A = A(n,u)": (
            A == A_from_n_u(n, u)
        ),

        "A-F(n) = u(u-n-1)": (
            A - F(n)
            == u * (u - n - 1)
        ),

        "discriminant square": (
            4 * A - 3 * (n - 1) ** 2
            == (2 * u - (n + 1)) ** 2
        ),
    }

    for name, result in checks.items():
        print(f"{name:<45s} {result}")


# =============================================================================
# SECTION 3
# ADAPTIVE AUXILIARY SEARCH
# =============================================================================

def adaptive_auxiliary_search(n: int):
    section("3. ADAPTIVE AUXILIARY-PRIME SEARCH")

    print()
    print(
        "Only n and r are used in this section."
    )

    print()
    print(
        "For each auxiliary r:"
    )
    print(
        "    m = F(r)"
    )
    print(
        "and we inspect gcd(m, n-only expressions)."
    )

    print()

    rows = []

    exprs = n_only_expressions(n)

    for r in PRIMES:

        if r > AUX_LIMIT:
            break

        m = F(r)

        gcds = {
            name: gcd(abs(value), m)
            for name, value in exprs.items()
        }

        # "strength" = largest nontrivial factor exposed.
        strength = max(gcds.values())

        rows.append(
            (
                r,
                m,
                bitlen(m),
                strength,
                gcds,
            )
        )

    rows.sort(
        key=lambda row: (
            row[3],
            row[2],
        ),
        reverse=True
    )

    print(
        "Top adaptive auxiliary moduli:"
    )
    line()

    print(
        "r    F(r)        bits   strongest gcd     expression"
    )
    line()

    shown = 0

    for r, m, m_bits, strength, gcds in rows:

        if strength <= 1:
            continue

        best_expr = None

        for name, g in gcds.items():
            if g == strength:
                best_expr = name
                break

        print(
            f"{r:3d}"
            f"{m:12d}"
            f"{m_bits:8d}"
            f"{strength:18d}"
            f"   {best_expr}"
        )

        shown += 1

        if shown >= 25:
            break

    if shown == 0:
        print(
            "No nontrivial gcd structure found."
        )


# =============================================================================
# SECTION 4
# FACTORIZATION OF F(r) AND N-ONLY RELATIONS
# =============================================================================

def section_factorization_profiles(n: int):
    section("4. FACTORS OF F(r) EXPOSED BY n")

    exprs = n_only_expressions(n)

    interesting = []

    for r in PRIMES:

        if r > AUX_LIMIT:
            break

        m = F(r)

        if m <= 1:
            continue

        gcds = {
            name: gcd(abs(value), m)
            for name, value in exprs.items()
        }

        nontrivial = [
            (name, g)
            for name, g in gcds.items()
            if 1 < g < m
        ]

        full_divisibility = [
            name
            for name, g in gcds.items()
            if g == m
        ]

        if nontrivial or full_divisibility:

            interesting.append(
                (
                    r,
                    m,
                    nontrivial,
                    full_divisibility
                )
            )

    if not interesting:
        print("No interesting profiles.")
        return

    for r, m, nontrivial, full_divisibility in interesting:

        print()
        print(
            f"r={r}, F(r)={m}, bits={bitlen(m)}"
        )

        if full_divisibility:
            print(
                "  FULL divisibility:"
            )

            for name in full_divisibility:
                print(
                    f"      F(r) | {name}"
                )

        if nontrivial:
            print(
                "  PARTIAL gcd factors:"
            )

            for name, g in nontrivial:
                print(
                    f"      gcd(F(r), {name}) = {g}"
                )


# =============================================================================
# SECTION 5
# DOES F(r) | n±1 GIVE A USEFUL u CONGRUENCE?
# =============================================================================

def divisibility_u_search(n: int, p: int, q: int, u: int, A: int):
    section("5. DIVISIBILITY-BASED u CONSTRAINT SEARCH")

    print()
    print(
        "This section searches for moduli m=F(r) satisfying"
    )
    print(
        "strong n-only divisibility conditions."
    )

    print(
        "It then asks what those conditions imply for the TRUE u"
    )
    print(
        "and, separately, what they would imply if A mod m were known."
    )

    print()
    print(
        "NO A IS USED TO SELECT THE AUXILIARY MODULI."
    )

    candidates = []

    exprs = n_only_expressions(n)

    for r in PRIMES:

        if r > AUX_LIMIT:
            break

        m = F(r)

        if m < 2:
            continue

        full = [
            name
            for name, value in exprs.items()
            if value % m == 0
        ]

        if full:

            candidates.append(
                (
                    r,
                    m,
                    full
                )
            )

    print()
    print(
        f"n-only divisible auxiliary moduli: {len(candidates)}"
    )

    if not candidates:
        print(
            "No r found with full divisibility."
        )
        return candidates

    print()
    print(
        "r    F(r)      condition(s) "
        "        true u mod F(r)"
    )
    line()

    for r, m, conditions in candidates:

        print(
            f"{r:3d}"
            f"{m:10d}"
            f"   {', '.join(conditions):<28s}"
            f"{u % m}"
        )

    return candidates


# =============================================================================
# SECTION 6
# CAN WE GET A WITHOUT A?
# =============================================================================

def search_n_only_a_residue(n: int, r: int):
    """
    Deliberately conservative.

    We examine possible congruence consequences from n mod m and
    algebraic identities, but we DO NOT use the hidden factorization.

    We ask whether n mod m itself uniquely determines a possible A mod m
    across all u residues.

    If every u gives the same A residue, then A mod m really is n-only.

    Otherwise n alone leaves ambiguity.
    """
    m = F(r)

    outputs = {}

    for ur in range(m):

        value = A_from_n_u(n % m, ur) % m

        outputs.setdefault(value, []).append(ur)

    return outputs


def section_n_only_A_residue(n: int):
    section("6. CAN n ALONE DETERMINE A mod F(r)?")

    print()
    print(
        "For each r we compute"
    )
    print(
        "    A(u) mod F(r)"
    )
    print(
        "for EVERY u residue."
    )

    print()
    print(
        "If there is only ONE possible A residue, then A mod F(r)"
    )
    print(
        "is genuinely determined by n alone."
    )

    print()
    print(
        "r    m=F(r)   distinct A residues   n-only?"
    )
    line()

    discoveries = []

    for r in PRIMES:

        if r > AUX_LIMIT:
            break

        m = F(r)

        # Don't attempt huge residue spaces.
        if m > 5000:
            continue

        mapping = search_n_only_a_residue(n, r)

        count = len(mapping)

        unique = count == 1

        print(
            f"{r:3d}"
            f"{m:10d}"
            f"{count:22d}"
            f"{str(unique):>12s}"
        )

        if unique:
            discoveries.append(
                (r, m, mapping)
            )

    print()

    if discoveries:
        print(
            "IMPORTANT: n-only A residue candidates found:"
        )

        for r, m, mapping in discoveries:
            print()
            print(
                f"r={r}, m={m}, mapping={mapping}"
            )

    else:
        print(
            "No modulus in the tested range made A mod F(r)"
        )
        print(
            "uniquely determined by n alone."
        )

    return discoveries


# =============================================================================
# SECTION 7
# POSSIBLE u RESIDUES FROM n-ONLY A STRUCTURE
# =============================================================================

def section_u_residue_candidates(
    n: int,
    u_true: int,
    A_true_value: int
):
    section("7. n-ONLY u-RESIDUE SPACE VS ORACLE u-RESIDUE SPACE")

    print()
    print(
        "For each modulus m:"
    )
    print(
        "  n-only: all u residues are possible unless an n-only"
    )
    print(
        "          relation reduces them."
    )
    print(
        "  oracle: A mod m is supplied."
    )

    print()
    print(
        "m       n-only u classes    oracle u classes    true u"
    )
    line()

    for r in PRIMES:

        if r > AUX_LIMIT:
            break

        m = F(r)

        if m > 5000:
            continue

        n_only_classes = m

        oracle_classes = oracle_u_residues(
            n,
            A_true_value,
            m
        )

        print(
            f"{m:5d}"
            f"{n_only_classes:22d}"
            f"{len(oracle_classes):21d}"
            f"{u_true % m:12d}"
        )


# =============================================================================
# SECTION 8
# ADAPTIVE CRT MODULUS BUILDING
# =============================================================================

def adaptive_crt_sequence(
    n: int,
    u_true: int,
    A_true_value: int
):
    section("8. ADAPTIVE CRT MODULUS GROWTH")

    print()
    print(
        "We build an n-only modulus from F(r) values selected because"
    )
    print(
        "they have strong gcd interaction with n."
    )

    print()
    print(
        "IMPORTANT:"
    )
    print(
        "The n-only branch NEVER uses A or u to choose r."
    )

    exprs = n_only_expressions(n)

    candidates = []

    for r in PRIMES:

        if r > AUX_LIMIT:
            break

        m = F(r)

        gcds = [
            gcd(abs(value), m)
            for value in exprs.values()
        ]

        score = max(gcds)

        candidates.append(
            (
                score,
                bitlen(m),
                r,
                m
            )
        )

    candidates.sort(
        key=lambda row: (
            row[0],
            row[1]
        ),
        reverse=True
    )

    selected = []

    M = 1

    print()
    print(
        "step   r     F(r)      gcd-score     CRT bits"
    )
    line()

    for score, m_bits, r, m in candidates:

        g = gcd(M, m)

        if m // g == 1:
            continue

        M_new = M // g * m

        selected.append(r)

        M = M_new

        print(
            f"{len(selected):4d}"
            f"{r:5d}"
            f"{m:12d}"
            f"{score:16d}"
            f"{bitlen(M):12d}"
        )

        if bitlen(M) >= CRT_MODULUS_BIT_LIMIT:
            break

    print()
    print(
        f"Selected r values: {selected}"
    )
    print(
        f"Final n-only CRT modulus: {M}"
    )
    print(
        f"Final modulus bits: {bitlen(M)}"
    )

    # Control: what happens if we actually know A mod each m?
    print()
    print(
        "ORACLE CONTROL:"
    )

    residues = []

    for r in selected:

        m = F(r)

        residues.append(
            (
                A_true_value % m,
                m
            )
        )

    crt_value = 0
    crt_modulus = 1
    consistent = True

    for residue, modulus in residues:

        result = crt_pair(
            crt_value,
            crt_modulus,
            residue,
            modulus
        )

        if result is None:
            consistent = False
            break

        crt_value, crt_modulus = result

    print(
        f"Oracle CRT A residue = {crt_value}"
    )
    print(
        f"Oracle CRT modulus   = {crt_modulus}"
    )
    print(
        f"True A                = {A_true_value}"
    )

    if consistent:

        exact = (
            crt_value == A_true_value % crt_modulus
        )

        print(
            f"Matches true A mod CRT modulus = {exact}"
        )

    # N-only branch.
    print()
    print(
        "N-ONLY interpretation:"
    )
    print(
        "The same CRT modulus exists, but there is no target residue."
    )
    print(
        "Therefore the modulus alone does not select A or u."
    )


# =============================================================================
# SECTION 9
# SEARCH SPECIAL DIVISIBILITY PATTERNS
# =============================================================================

def section_special_patterns(n: int):
    section("9. SPECIAL n-ONLY DIVISIBILITY PATTERNS")

    print()
    print(
        "Search for r where F(r) divides combinations naturally"
    )
    print(
        "associated with n."
    )

    expressions = n_only_expressions(n)

    for r in PRIMES:

        if r > AUX_LIMIT:
            break

        m = F(r)

        hits = []

        for name, value in expressions.items():

            if value % m == 0:
                hits.append(name)

        if hits:

            print()
            print(
                f"r={r:3d}, F(r)={m}"
            )

            for hit in hits:
                print(
                    f"    F(r) | {hit}"
                )


# =============================================================================
# SECTION 10
# MULTIPLE TARGETS
# =============================================================================

def section_multi_target_test(targets):
    section("10. MULTI-TARGET ADAPTIVE SEARCH")

    print()
    print(
        "We repeat the n-only adaptive search across independent"
    )
    print(
        "semiprimes and look for recurring modulus patterns."
    )

    recurring = {}

    for idx, (p, q) in enumerate(targets, 1):

        n = p * q

        exprs = n_only_expressions(n)

        local = []

        for r in PRIMES:

            if r > AUX_LIMIT:
                break

            m = F(r)

            score = max(
                gcd(abs(value), m)
                for value in exprs.values()
            )

            if score > 1:
                local.append(
                    (r, m, score)
                )

        local.sort(
            key=lambda row: (
                row[2],
                bitlen(row[1])
            ),
            reverse=True
        )

        best = tuple(
            r
            for r, _, _ in local[:10]
        )

        for r in best:
            recurring[r] = recurring.get(r, 0) + 1

        print()
        print(
            f"target {idx}: n={n}"
        )
        print(
            f"  p={p}"
        )
        print(
            f"  q={q}"
        )
        print(
            f"  best r={best}"
        )

    print()
    print(
        "RECURRING AUXILIARY PRIMES"
    )
    line()

    total = len(targets)

    for r, count in sorted(
        recurring.items(),
        key=lambda item: (
            item[1],
            item[0]
        ),
        reverse=True
    ):

        print(
            f"r={r:3d} "
            f"appears in {count}/{total} target top-sets"
        )


# =============================================================================
# SECTION 11
# FACTOR-RECOVERY CONTROL
# =============================================================================

def section_factor_recovery_control(
    p: int,
    q: int,
    n: int,
    u: int,
    A: int
):
    section("11. ORACLE FACTOR-RECOVERY CONTROL")

    D = (
        4 * A
        - 3 * (n - 1) ** 2
    )

    sqrt_D = integer_sqrt_exact(D)

    print()
    print(
        "Using exact A:"
    )

    print(
        f"D       = {D}"
    )

    print(
        f"sqrt(D) = {sqrt_D}"
    )

    if sqrt_D is None:
        print(
            "ERROR: discriminant is not a square."
        )
        return

    u_candidates = [
        (
            n + 1 + sqrt_D
        ) // 2,

        (
            n + 1 - sqrt_D
        ) // 2
    ]

    print()
    print(
        f"u candidates = {u_candidates}"
    )

    for uc in u_candidates:

        df = uc * uc - 4 * n

        sqrt_df = integer_sqrt_exact(df)

        if sqrt_df is None:
            continue

        p1 = (uc + sqrt_df) // 2
        q1 = (uc - sqrt_df) // 2

        print(
            f"u={uc} -> factors {(p1, q1)}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    random.seed(RANDOM_SEED)

    print()
    print("=" * 78)
    print(
        "KAPPA ADAPTIVE CYCLOTOMIC MODULUS / N-ONLY SEARCH"
    )
    print("=" * 78)

    print()
    print(
        "Central question:"
    )
    print()
    print(
        "Can n-only arithmetic choose cyclotomic moduli"
    )
    print(
        "that constrain A or u without using the hidden factors?"
    )

    # Representative.
    p, q, n, u, A = representative_target()

    verify_basic_structure(
        p,
        q,
        n,
        u,
        A
    )

    adaptive_auxiliary_search(n)

    section_factorization_profiles(n)

    divisibility_u_search(
        n,
        p,
        q,
        u,
        A
    )

    section_n_only_A_residue(n)

    section_u_residue_candidates(
        n,
        u,
        A
    )

    adaptive_crt_sequence(
        n,
        u,
        A
    )

    section_special_patterns(n)

    # Additional targets.
    targets = []

    for bits in TARGET_BITS:

        for _ in range(TARGETS_PER_SIZE):

            try:
                pv, qv, _ = make_semiprime(bits)
                targets.append((pv, qv))
            except RuntimeError:
                pass

    section_multi_target_test(targets)

    section_factor_recovery_control(
        p,
        q,
        n,
        u,
        A
    )

    # Final report.
    section("12. FINAL REPORT")

    print()
    print(
        "The experiment separates three levels:"
    )

    print()
    print(
        "LEVEL 1: n-only"
    )
    print(
        "    n -> F(r) -> gcd/divisibility structure"
    )

    print()
    print(
        "LEVEL 2: n + A residue"
    )
    print(
        "    n + A mod F(r) -> u residues"
    )

    print()
    print(
        "LEVEL 3: exact A"
    )
    print(
        "    n + A -> u -> p,q"
    )

    print()
    print(
        "The desired discovery is a LEVEL-1 mechanism that"
    )
    print(
        "produces a nontrivial A or u residue."
    )

    print()
    print(
        "Do NOT interpret a large CRT modulus by itself as"
    )
    print(
        "information about the hidden factorization."
    )

    print()
    print(
        "A successful experiment should show:"
    )
    print()
    print(
        "    n-only condition"
    )
    print(
        "          ↓"
    )
    print(
        "    specific A mod m"
    )
    print(
        "          ↓"
    )
    print(
        "    few u residues"
    )
    print(
        "          ↓"
    )
    print(
        "    recoverable p,q"
    )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

