#!/usr/bin/env python3

"""
==============================================================================
KAPPA NEXT EXPERIMENT
SIXTH-ROOT / F(r)=r^2-r+1 / N-ONLY ALGEBRAIC INTERACTION SEARCH
==============================================================================

Goal:

    Stop manipulating K-values.

    Search directly for algebraic structure connecting

        n
        u
        A(u) = u^2 -(n+1)u+n^2-n+1

    with the auxiliary polynomial

        F(r) = r^2-r+1.

Important:

    The discovery phase must NOT use p, q, u, A, or K.

    They are used only afterward as validation/oracle controls.

Core observation:

    F(r) = r^2-r+1

    and therefore modulo F(r):

        r^2 = r-1
        r^3 = -1
        r^6 = 1

This experiment asks whether that sixth-root structure creates
a nontrivial constraint involving n or u.

No CSV files are written.
All important results are printed.
"""


from math import gcd, isqrt
from itertools import combinations, product
from fractions import Fraction
import random


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814
random.seed(SEED)

TARGETS_PER_SIZE = 12
N_BITS_LIST = [20, 30, 40, 50]

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# A small set used for deeper symbolic searches.
R_DEEP = [2, 3, 5, 7, 11, 13, 17]

# Coefficients for small polynomial searches.
SMALL_COEFFS = range(-3, 4)


# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def next_prime(n):
    x = n + 1
    while not is_prime(x):
        x += 1
    return x


def random_prime_with_bits(bits):
    lo = 1 << (bits - 1)
    hi = (1 << bits) - 1

    while True:
        x = random.randrange(lo | 1, hi + 1, 2)
        if is_prime(x):
            return x


def F(r):
    return r * r - r + 1


def A_from_pq(p, q):
    n = p * q
    u = p + q
    return u * u - (n + 1) * u + n * n - n + 1


def A_from_nu(n, u):
    return u * u - (n + 1) * u + n * n - n + 1


def discriminant_from_A(n, A):
    return 4 * A - 3 * (n - 1) ** 2


def factor_pair_from_u(n, u):
    D = u * u - 4 * n
    if D < 0:
        return None

    s = isqrt(D)
    if s * s != D:
        return None

    p = (u - s) // 2
    q = (u + s) // 2

    if p * q == n:
        return p, q

    return None


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets():
    targets = []

    for bits in N_BITS_LIST:
        for _ in range(TARGETS_PER_SIZE):
            p_bits = bits // 2
            q_bits = bits - p_bits

            p = random_prime_with_bits(p_bits)
            q = random_prime_with_bits(q_bits)

            while p == q:
                q = random_prime_with_bits(q_bits)

            n = p * q
            u = p + q
            A = A_from_pq(p, q)

            targets.append({
                "p": p,
                "q": q,
                "n": n,
                "u": u,
                "A": A,
                "n_bits": n.bit_length(),
                "u_bits": u.bit_length(),
                "A_bits": A.bit_length(),
            })

    return targets


# ============================================================================
# MODULAR ROOT / IMAGE STRUCTURE
# ============================================================================

def roots_F_mod_m(m):
    """
    Find roots of r^2-r+1 modulo m.

    For the small m used here a direct scan is sufficient.
    """
    return [r for r in range(m) if F(r) % m == 0]


def polynomial_A_mod(n, u, m):
    return A_from_nu(n, u) % m


def quadratic_u_residues_from_A(n, a, m):
    return [
        u for u in range(m)
        if A_from_nu(n, u) % m == a % m
    ]


# ============================================================================
# SYMBOLIC REDUCTION MOD F(r)
# ============================================================================

def reduce_poly_in_r(coeffs):
    """
    Reduce a polynomial in r modulo

        r^2-r+1.

    coeffs[k] is coefficient of r^k.

    Since r^2 = r-1, every polynomial reduces to:

        c0 + c1*r.
    """

    coeffs = list(coeffs)

    while len(coeffs) > 2:
        k = len(coeffs) - 1
        c = coeffs.pop()

        # r^k = r^(k-2) * r^2
        #      = r^(k-1) - r^(k-2)
        if k - 1 >= len(coeffs):
            coeffs.extend([0] * (k - len(coeffs)))

        coeffs[k - 1] += c
        coeffs[k - 2] -= c

    while len(coeffs) < 2:
        coeffs.append(0)

    return coeffs[0], coeffs[1]


def sixth_root_power_reduction(k):
    """
    Return representation of r^k modulo F(r)
    in the form a + b*r.
    """

    if k == 0:
        return 1, 0

    # r^6 = 1
    k %= 6

    table = {
        0: (1, 0),
        1: (0, 1),
        2: (-1, 1),
        3: (-1, 0),
        4: (0, -1),
        5: (1, -1),
    }

    return table[k]


# ============================================================================
# N-ONLY AUXILIARY CONGRUENCE STRUCTURE
# ============================================================================

def n_only_auxiliary_table(n):
    rows = []

    for r in R_VALUES:
        m = F(r)

        rows.append({
            "r": r,
            "m": m,
            "g_n": gcd(n, m),
            "g_nm1": gcd(n - 1, m),
            "g_np1": gcd(n + 1, m),
            "n_mod": n % m,
            "n2_mod": (n * n) % m,
            "roots_F": roots_F_mod_m(m),
        })

    return rows


# ============================================================================
# TEST 1: DIRECT SIXTH-ROOT STRUCTURE
# ============================================================================

def test_sixth_root_structure():
    print()
    print("=" * 78)
    print("1. SIXTH-ROOT STRUCTURE OF F(r)=r^2-r+1")
    print("=" * 78)

    print()
    print("F(r) = r^2-r+1")
    print()
    print("Modulo F(r):")
    print("  r^2 = r-1")
    print("  r^3 = -1")
    print("  r^6 = 1")
    print()

    for r in R_DEEP:
        m = F(r)

        values = []
        x = r % m

        for k in range(1, 7):
            values.append(pow(x, k, m))

        print(
            f"r={r:2d} F={m:4d} "
            f"powers={values} "
            f"r^6 mod F={pow(r,6,m)}"
        )


# ============================================================================
# TEST 2: A(u) MOD F(r)
# ============================================================================

def test_A_mod_F(targets):
    print()
    print("=" * 78)
    print("2. A(u) MOD F(r): TARGET-DEPENDENT CONTROL")
    print("=" * 78)

    print()
    print("This section uses true u/A only as a diagnostic.")
    print("It is NOT part of the n-only discovery phase.")
    print()

    for idx, t in enumerate(targets[:5], 1):
        n = t["n"]
        u = t["u"]
        A = t["A"]

        print(
            f"target={idx:2d} "
            f"n_bits={t['n_bits']:2d} "
            f"u_bits={t['u_bits']:2d} "
            f"A_bits={t['A_bits']:3d}"
        )

        for r in R_DEEP:
            m = F(r)

            print(
                f"  r={r:2d} "
                f"F={m:4d} "
                f"A mod F={A % m:4d} "
                f"u mod F={u % m:4d}"
            )

        print()


# ============================================================================
# TEST 3: N-ONLY POLYNOMIAL REDUCTION
# ============================================================================

def test_n_only_reduction(targets):
    print()
    print("=" * 78)
    print("3. N-ONLY POLYNOMIAL REDUCTION MOD F(r)")
    print("=" * 78)

    print()
    print(
        "We reduce expressions involving n modulo F(r), "
        "without using p,q,u,A or K."
    )
    print()

    expressions = {
        "n": lambda n: n,
        "n+1": lambda n: n + 1,
        "n-1": lambda n: n - 1,
        "n^2-1": lambda n: n * n - 1,
        "n^2-n+1": lambda n: n * n - n + 1,
        "n^2+n+1": lambda n: n * n + n + 1,
    }

    for idx, t in enumerate(targets[:4], 1):
        n = t["n"]

        print(f"target {idx}: n={n}")

        for name, fn in expressions.items():
            residues = []

            for r in R_DEEP:
                m = F(r)
                residues.append(fn(n) % m)

            print(f"  {name:12s}: {residues}")

        print()


# ============================================================================
# TEST 4: CAN F(r) DIVIDE SIMPLE N-ONLY EXPRESSIONS?
# ============================================================================

def test_divisibility(targets):
    print()
    print("=" * 78)
    print("4. N-ONLY DIVISIBILITY SEARCH")
    print("=" * 78)

    expressions = {
        "n": lambda n: n,
        "n-1": lambda n: n - 1,
        "n+1": lambda n: n + 1,
        "n^2-1": lambda n: n * n - 1,
        "n^2-n+1": lambda n: n * n - n + 1,
        "n^2+n+1": lambda n: n * n + n + 1,
        "n^2+n-1": lambda n: n * n + n - 1,
        "n^2-n-1": lambda n: n * n - n - 1,
    }

    for name, fn in expressions.items():
        hits = 0
        total = 0
        examples = []

        for t in targets:
            n = t["n"]

            for r in R_VALUES:
                m = F(r)

                if fn(n) % m == 0:
                    hits += 1
                    if len(examples) < 6:
                        examples.append(
                            (n, r, m)
                        )

                total += 1

        print(
            f"{name:14s} "
            f"hits={hits:4d}/{total:4d} "
            f"fraction={hits/total:.4f}"
        )

        if examples:
            print(f"  examples={examples}")


# ============================================================================
# TEST 5: N-ONLY RESIDUE RELATIONS
# ============================================================================

def test_residue_relations(targets):
    print()
    print("=" * 78)
    print("5. N-ONLY RESIDUE RELATION SEARCH")
    print("=" * 78)

    """
    Search relations such as:

        n mod F(r) = constant
        n^2 mod F(r) = ...
        n(n-1) mod F(r) = ...
        n^2-n+1 mod F(r) = ...

    We are specifically looking for relations whose value depends
    on r but is independent of target.
    """

    expressions = {
        "n": lambda n: n,
        "n-1": lambda n: n - 1,
        "n+1": lambda n: n + 1,
        "n(n-1)": lambda n: n * (n - 1),
        "n(n+1)": lambda n: n * (n + 1),
        "n^2-1": lambda n: n * n - 1,
        "n^2-n+1": lambda n: n * n - n + 1,
        "n^2+n+1": lambda n: n * n + n + 1,
    }

    for name, fn in expressions.items():
        print()
        print(f"Expression: {name}")

        for r in R_DEEP:
            m = F(r)
            values = sorted({fn(t["n"]) % m for t in targets})

            print(
                f"  r={r:2d} F={m:4d} "
                f"distinct={len(values):2d} "
                f"values={values[:12]}"
            )


# ============================================================================
# TEST 6: SEARCH FOR A UNIVERSAL n MOD F(r) PATTERN
# ============================================================================

def test_universal_pattern(targets):
    print()
    print("=" * 78)
    print("6. UNIVERSAL POLYNOMIAL PATTERN SEARCH")
    print("=" * 78)

    """
    Search small expressions

        c0 + c1*n + c2*n^2

    and determine whether their residue modulo F(r)
    is universal across targets for each r.

    A universal zero would be particularly interesting.
    """

    candidates = []

    for c0 in SMALL_COEFFS:
        for c1 in SMALL_COEFFS:
            for c2 in SMALL_COEFFS:
                if c0 == c1 == c2 == 0:
                    continue

                candidates.append((c0, c1, c2))

    hits = []

    for c0, c1, c2 in candidates:
        all_zero = True
        all_universal = True

        for r in R_DEEP:
            m = F(r)

            vals = {
                (c0 + c1 * t["n"] + c2 * t["n"] * t["n"]) % m
                for t in targets
            }

            if len(vals) != 1:
                all_universal = False

            if vals != {0}:
                all_zero = False

        if all_zero:
            hits.append(
                ("ZERO", c0, c1, c2)
            )
        elif all_universal:
            hits.append(
                ("UNIVERSAL", c0, c1, c2)
            )

    print()
    print(f"candidate polynomials tested = {len(candidates)}")
    print(f"universal candidates         = {len(hits)}")

    for item in hits[:30]:
        print(
            f"  {item[0]:10s}: "
            f"{item[1]} + ({item[2]})n + ({item[3]})n^2"
        )

    if not hits:
        print("  No universal quadratic n-only relation found.")


# ============================================================================
# TEST 7: RESULTANT-STYLE SEARCH
# ============================================================================

def test_resultant_style(targets):
    print()
    print("=" * 78)
    print("7. RESULTANT-STYLE n/u ELIMINATION SEARCH")
    print("=" * 78)

    """
    A(u) is

        u^2 -(n+1)u+n^2-n+1.

    We ask whether simple candidate relations

        G(n,u,r)=0 mod F(r)

    can eliminate u.

    Discovery is intentionally restricted to expressions
    generated from n only after eliminating u symbolically.

    We test several natural combinations of the coefficients
    of A(u).
    """

    candidates = {
        "discriminant-coeff": lambda n: 3 * (n - 1) ** 2,
        "constant-term": lambda n: n * n - n + 1,
        "linear-coeff": lambda n: n + 1,
        "quadratic-coeff": lambda n: 1,
        "n^2-n": lambda n: n * n - n,
        "n^2-1": lambda n: n * n - 1,
        "(n-1)^2": lambda n: (n - 1) ** 2,
        "(n+1)^2": lambda n: (n + 1) ** 2,
    }

    for name, fn in candidates.items():
        universal = True
        value_map = {}

        for r in R_DEEP:
            m = F(r)

            vals = sorted({
                fn(t["n"]) % m
                for t in targets
            })

            value_map[r] = vals

            if len(vals) != 1:
                universal = False

        print(
            f"{name:20s} "
            f"universal_across_targets={universal}"
        )

        if universal:
            print(f"  values={value_map}")


# ============================================================================
# TEST 8: SIXTH-ROOT RESIDUE ORBITS
# ============================================================================

def test_orbits():
    print()
    print("=" * 78)
    print("8. SIXTH-ROOT ORBIT STRUCTURE")
    print("=" * 78)

    print()
    print(
        "For each r, examine the orbit generated by multiplication "
        "by r modulo F(r)."
    )
    print()

    for r in R_DEEP:
        m = F(r)

        x = 1
        orbit = []

        for _ in range(6):
            orbit.append(x)
            x = (x * r) % m

        print(
            f"r={r:2d} F={m:4d} "
            f"orbit={orbit} "
            f"returns_to_1={x == 1}"
        )


# ============================================================================
# TEST 9: N-ONLY CRT OF AUXILIARY MODULI
# ============================================================================

def test_auxiliary_crt(targets):
    print()
    print("=" * 78)
    print("9. N-ONLY AUXILIARY CRT STRUCTURE")
    print("=" * 78)

    """
    Combine F(r) moduli.

    We do NOT insert A residues.

    Instead we measure whether the resulting modulus creates
    any special gcd relationship with n, n-1, n+1, etc.
    """

    M = 1

    for r in R_VALUES:
        m = F(r)

        if gcd(M, m) == 1:
            M *= m

        print(
            f"r={r:2d} "
            f"F={m:5d} "
            f"M_bits={M.bit_length():3d} "
            f"gcd(M,n)={gcd(M, targets[0]['n']):3d} "
            f"gcd(M,n-1)={gcd(M, targets[0]['n']-1):3d} "
            f"gcd(M,n+1)={gcd(M, targets[0]['n']+1):3d}"
        )

    print()
    print(f"final M = {M}")
    print(f"final M bits = {M.bit_length()}")


# ============================================================================
# TEST 10: CAN n ALONE PREDICT ANY A MOD F(r)?
# ============================================================================

def test_prediction_boundary(targets):
    print()
    print("=" * 78)
    print("10. N-ONLY -> A MOD F(r) PREDICTION BOUNDARY")
    print("=" * 78)

    """
    This is a deliberately important control.

    For each r, determine whether A mod F(r) is uniquely determined
    by n mod F(r) across all generated targets.

    If multiple A residues occur for the same n residue,
    then n's residue alone does not determine A's residue.
    """

    for r in R_DEEP:
        m = F(r)

        mapping = {}

        for t in targets:
            n_res = t["n"] % m
            a_res = t["A"] % m

            mapping.setdefault(n_res, set()).add(a_res)

        ambiguous = sum(
            1 for vals in mapping.values()
            if len(vals) > 1
        )

        print(
            f"r={r:2d} F={m:4d} "
            f"n_residues={len(mapping):4d} "
            f"ambiguous_n_residues={ambiguous:4d}"
        )

        examples = [
            (nr, sorted(list(vals)))
            for nr, vals in mapping.items()
            if len(vals) > 1
        ]

        if examples:
            print(
                f"  sample ambiguity: {examples[:4]}"
            )


# ============================================================================
# TEST 11: VALIDATION — ORACLE A -> u
# ============================================================================

def validate_oracle(target):
    print()
    print("=" * 78)
    print("11. CONTROL: A -> u")
    print("=" * 78)

    n = target["n"]
    A = target["A"]
    true_u = target["u"]

    D = discriminant_from_A(n, A)

    print(f"n       = {n}")
    print(f"A       = {A}")
    print(f"true u  = {true_u}")
    print(f"D       = {D}")
    print(f"isqrt(D)^2 == D : {isqrt(D) ** 2 == D}")

    if isqrt(D) ** 2 == D:
        s = isqrt(D)

        candidates = [
            ((n + 1 + s) // 2),
            ((n + 1 - s) // 2),
        ]

        print(f"u candidates from A = {candidates}")
        print(f"true u recovered    = {true_u in candidates}")

        for u in candidates:
            pair = factor_pair_from_u(n, u)
            print(f"  u={u} -> pair={pair}")


# ============================================================================
# TEST 12: SEARCH FOR AN ACTUAL N-ONLY FACTOR SIGNAL
# ============================================================================

def test_n_only_factor_signal(targets):
    print()
    print("=" * 78)
    print("12. DIRECT N-ONLY FACTOR SIGNAL SEARCH")
    print("=" * 78)

    """
    We deliberately test only quantities computable from n.

    Candidate signals:

        gcd(n, F(r))
        gcd(n-1, F(r))
        gcd(n+1, F(r))
        gcd(n^2-1, F(r))
        n mod F(r)

    A useful signal would distinguish targets according to p/q
    substantially better than random residue behavior.
    """

    signal_names = [
        "gcd(n,F)",
        "gcd(n-1,F)",
        "gcd(n+1,F)",
        "gcd(n^2-1,F)",
    ]

    for name in signal_names:
        values = []

        for t in targets:
            row = []

            n = t["n"]

            for r in R_DEEP:
                m = F(r)

                if name == "gcd(n,F)":
                    x = gcd(n, m)
                elif name == "gcd(n-1,F)":
                    x = gcd(n - 1, m)
                elif name == "gcd(n+1,F)":
                    x = gcd(n + 1, m)
                else:
                    x = gcd(n * n - 1, m)

                row.append(x)

            values.append(tuple(row))

        distinct = len(set(values))

        print(
            f"{name:16s} "
            f"distinct target signatures={distinct}/{len(targets)}"
        )


# ============================================================================
# FINAL REPORT
# ============================================================================

def final_report(targets):
    print()
    print("=" * 78)
    print("FINAL REPORT")
    print("=" * 78)

    print()
    print("Experiment:")
    print("  F(r) = r^2-r+1")
    print("  sixth-root relation: r^6 = 1 mod F(r)")
    print()

    print("Discovery constraints:")
    print("  p,q,u,A,K are NOT used to discover n-only relations.")
    print()

    print("Targets:")
    print(f"  total = {len(targets)}")
    print(
        f"  n-bit range = "
        f"{min(t['n_bits'] for t in targets)}.."
        f"{max(t['n_bits'] for t in targets)}"
    )

    print()
    print("The experiment is looking specifically for:")
    print()
    print("    n")
    print("     |")
    print("     v")
    print("    F(r)=r^2-r+1")
    print("     |")
    print("     v")
    print("    nontrivial residue/algebraic relation")
    print("     |")
    print("     v")
    print("    constraint on u or A")
    print()

    print("A relation counts as a breakthrough only if:")
    print("  1. it is computable from n and r alone,")
    print("  2. it survives across independent targets,")
    print("  3. it is not a defining identity,")
    print("  4. it gives information about u/A/factors.")
    print()

    print("No CSV files are produced.")
    print("=" * 78)


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 78)
    print("KAPPA NEXT EXPERIMENT")
    print("SIXTH-ROOT / F(r)=r^2-r+1 / N-ONLY ALGEBRAIC INTERACTION SEARCH")
    print("=" * 78)

    print()
    print(f"random seed = {SEED}")
    print(f"R values    = {R_VALUES}")
    print(f"target sizes = {N_BITS_LIST}")
    print(f"targets/size = {TARGETS_PER_SIZE}")

    targets = generate_targets()

    representative = targets[-1]

    print()
    print("=" * 78)
    print("REPRESENTATIVE TARGET")
    print("=" * 78)

    print(f"p       = {representative['p']}")
    print(f"q       = {representative['q']}")
    print(f"n       = {representative['n']}")
    print(f"n bits  = {representative['n_bits']}")
    print(f"u       = {representative['u']}")
    print(f"u bits  = {representative['u_bits']}")
    print(f"A bits  = {representative['A_bits']}")

    test_sixth_root_structure()
    test_A_mod_F(targets)
    test_n_only_reduction(targets)
    test_divisibility(targets)
    test_residue_relations(targets)
    test_universal_pattern(targets)
    test_resultant_style(targets)
    test_orbits()
    test_auxiliary_crt(targets)
    test_prediction_boundary(targets)
    validate_oracle(representative)
    test_n_only_factor_signal(targets)
    final_report(targets)


if __name__ == "__main__":
    main()
