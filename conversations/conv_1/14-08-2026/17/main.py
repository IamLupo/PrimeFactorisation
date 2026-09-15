#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 59
CUBIC-CHARACTER ROOT-RATIO TEST
TRUE FACTOR SUM VS EXACT MATCHED DISCRIMINANT-QR FALSE NULL
CYCLOTOMIC VS RANDOM CONTROL
NO CSV OUTPUT
==============================================================================

Goal
----
For ell == 1 (mod 3), test the cubic residue character of the ratio of the
two roots of

    x^2 - s*x + n == 0 (mod ell)

For a genuine factorization n = p*q, the two roots are p,q modulo ell.

We test whether the true sum has an unusual cubic-character signature compared
with false sums drawn from the COMPLETE set of sums satisfying the same
cyclotomic discriminant-QR conditions.

This deliberately avoids the symmetric resultant R(n,s) branch.

Definitions
-----------
Let g be a primitive cubic-root generator:

    zeta^3 = 1, zeta != 1

For a nonzero a:

    a^((ell-1)/3) mod ell in {1, zeta, zeta^2}

This gives the cubic character class:

    0  -> chi_3(a) = 1
    1  -> chi_3(a) = zeta
    2  -> chi_3(a) = zeta^2

For candidate roots x,y, define:

    rho = x * y^(-1) mod ell

and test the cubic class of rho.

The unordered version uses:

    class(rho) in {0,1,2}
    nontrivial = class != 0

because swapping x,y maps class k -> -k mod 3.

We also retain an ordered canonical version for diagnostics by ordering the
two residue roots numerically.

Important
---------
This experiment does NOT assume that the cubic signature is a valid sieve.
It tests whether the TRUE sums have a distribution different from matched
FALSE discriminant-compatible sums.

A positive result requires:

    TRUE vs FALSE separation
        AND
    CYCLOTOMIC separation > RANDOM CONTROL separation

==============================================================================
"""

from __future__ import annotations

import math
import random
import time
from collections import Counter
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SEED = 20260814

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

TARGET_COUNT = 60

# Base cyclotomic family from the previous experiments.
CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

# Random control family, chosen with ell == 1 mod 3 so that a cubic character
# exists here as well.
CONTROL = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

# Number of false sums sampled per target. If the complete pool is smaller,
# use the complete pool.
FALSE_SAMPLE_CAP = 240

# For generating exact matched false pools.
# We require:
#     s even
#     D = s^2 - 4n is QR mod every CYCLOTOMIC ell
#
# The true sum is explicitly excluded.
#
# No prime-pair enumeration is performed.


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# ---------------------------------------------------------------------------
# Prime generation
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> list[int]:
    """Return primes in [lo, hi)."""
    if hi <= 2:
        return []

    size = hi + 1
    is_prime = bytearray(b"\x01") * size
    is_prime[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi)) + 1

    for p in range(2, limit):
        if is_prime[p]:
            start = p * p
            is_prime[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [x for x in range(lo, hi) if is_prime[x]]


def generate_targets(primes: list[int], count: int, rng: random.Random):
    """Generate distinct semiprime targets."""
    targets: list[Target] = []
    seen: set[int] = set()

    while len(targets) < count:
        p, q = rng.sample(primes, 2)
        if p > q:
            p, q = q, p

        n = p * q
        if n in seen:
            continue

        seen.add(n)
        targets.append(Target(
            p=p,
            q=q,
            n=n,
            s=p + q,
        ))

    return targets


# ---------------------------------------------------------------------------
# Quadratic-residue utilities
# ---------------------------------------------------------------------------

def qr_table(mod: int) -> set[int]:
    return {(x * x) % mod for x in range(mod)}


def is_qr_table(x: int, mod: int, table: set[int]) -> bool:
    return (x % mod) in table


# ---------------------------------------------------------------------------
# Cubic-character utilities
# ---------------------------------------------------------------------------

def primitive_root_prime(p: int) -> int:
    """
    Return a primitive root modulo an odd prime p.

    Straightforward trial implementation; all moduli here are small.
    """
    phi = p - 1

    factors = []
    x = phi
    d = 2

    while d * d <= x:
        if x % d == 0:
            factors.append(d)
            while x % d == 0:
                x //= d
        d += 1 if d == 2 else 2

    if x > 1:
        factors.append(x)

    for g in range(2, p):
        ok = True
        for f in factors:
            if pow(g, phi // f, p) == 1:
                ok = False
                break
        if ok:
            return g

    raise RuntimeError(f"No primitive root found for {p}")


def cubic_root_generator(mod: int) -> int:
    """
    Produce a nontrivial cube root of unity modulo mod.

    Since mod == 1 (mod 3), let g be primitive and set

        zeta = g^((mod-1)/3).

    Then zeta^3 = 1 and zeta != 1.
    """
    if mod % 3 != 1:
        raise ValueError(f"{mod} is not 1 mod 3")

    g = primitive_root_prime(mod)
    zeta = pow(g, (mod - 1) // 3, mod)

    assert zeta != 1
    assert pow(zeta, 3, mod) == 1

    return zeta


def build_cubic_character_table(mod: int) -> tuple[list[int], int]:
    """
    Return:
        character[a] in {0,1,2}
    where 0,1,2 represent powers of the chosen cubic root of unity.

    For a != 0:
        a^((ell-1)/3) = zeta^k
    """
    zeta = cubic_root_generator(mod)

    powers = {
        1: 0,
        zeta: 1,
        (zeta * zeta) % mod: 2,
    }

    character = [-1] * mod
    character[0] = -1

    exp = (mod - 1) // 3

    for a in range(1, mod):
        value = pow(a, exp, mod)

        if value not in powers:
            raise AssertionError(
                f"Unexpected cubic character value {value} mod {mod}"
            )

        character[a] = powers[value]

    return character, zeta


# ---------------------------------------------------------------------------
# Quadratic roots modulo ell
# ---------------------------------------------------------------------------

def modular_sqrt_prime(a: int, p: int) -> int | None:
    """Tonelli-Shanks for odd prime p."""
    a %= p

    if a == 0:
        return 0

    if pow(a, (p - 1) // 2, p) != 1:
        return None

    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    # Factor p-1 = q*2^s with q odd.
    q = p - 1
    s = 0

    while q % 2 == 0:
        s += 1
        q //= 2

    # Find quadratic non-residue z.
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)

    while t != 1:
        i = 1
        t2i = (t * t) % p

        while t2i != 1:
            t2i = (t2i * t2i) % p
            i += 1

            if i >= m:
                return None

        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p

    return r


def quadratic_roots_mod(s: int, n: int, mod: int) -> tuple[int, int] | None:
    """
    Solve x^2 - s*x + n == 0 mod mod.
    """
    D = (s * s - 4 * n) % mod
    root = modular_sqrt_prime(D, mod)

    if root is None:
        return None

    inv2 = pow(2, -1, mod)

    x1 = ((s + root) * inv2) % mod
    x2 = ((s - root) * inv2) % mod

    return x1, x2


# ---------------------------------------------------------------------------
# Cubic signature
# ---------------------------------------------------------------------------

def cubic_ratio_signature(
    s: int,
    n: int,
    mod: int,
    char_table: list[int],
):
    """
    Return a dictionary describing the cubic-character relationship of the
    two quadratic roots.

    None means:
        - discriminant is non-residue, or
        - one of the roots is 0 mod ell.

    unordered_class:
        0 => ratio is a cubic residue
        1/2 => nontrivial cubic class

    ordered_pair:
        (chi(x_small), chi(x_large))

    ratio_class:
        chi(x/y)
    """
    roots = quadratic_roots_mod(s, n, mod)
    if roots is None:
        return None

    x, y = roots

    # We need both roots nonzero for x/y.
    if x == 0 or y == 0:
        return None

    if x > y:
        x, y = y, x

    cx = char_table[x]
    cy = char_table[y]

    if cx < 0 or cy < 0:
        return None

    ratio = (x * pow(y, -1, mod)) % mod
    cr = char_table[ratio]

    # Character must satisfy:
    #     chi(x/y) = chi(x) - chi(y) mod 3
    check = (cx - cy) % 3

    if check != cr:
        raise AssertionError(
            f"Cubic character inconsistency mod {mod}: "
            f"x={x}, y={y}, cx={cx}, cy={cy}, cr={cr}"
        )

    # Unordered signature: 0 vs nonzero is swap-invariant.
    return {
        "ordered_pair": (cx, cy),
        "ratio_class": cr,
        "trivial_ratio": (cr == 0),
        "nontrivial_ratio": (cr != 0),
        "roots": (x, y),
    }


# ---------------------------------------------------------------------------
# False matched pools
# ---------------------------------------------------------------------------

def candidate_sums_for_n(
    n: int,
    s_min: int = 4_000_000,
    s_max: int = 8_399_998,
):
    if s_min % 2 != 0:
        s_min += 1

    if s_max % 2 != 0:
        s_max -= 1

    return range(s_min, s_max + 1, 2)


def build_exact_false_pool(
    n: int,
    true_s: int,
    qr_tables: dict[int, set[int]],
):
    """
    Complete set of even sums in the domain satisfying the same cyclotomic
    discriminant-QR conditions, excluding the true sum.
    """
    out = []

    cyclo = list(qr_tables)

    for s in candidate_sums_for_n(n):
        if s == true_s:
            continue

        D = s * s - 4 * n

        good = True

        for ell in cyclo:
            if (D % ell) not in qr_tables[ell]:
                good = False
                break

        if good:
            out.append(s)

    return out


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def entropy(counter: Counter) -> float:
    total = sum(counter.values())

    if total == 0:
        return 0.0

    h = 0.0

    for c in counter.values():
        p = c / total
        h -= p * math.log2(p)

    return h


def mutual_information_3x3(
    pairs: list[tuple[int, int]]
) -> float:
    if not pairs:
        return 0.0

    joint = Counter(pairs)
    xcount = Counter(x for x, _ in pairs)
    ycount = Counter(y for _, y in pairs)

    n = len(pairs)
    mi = 0.0

    for (x, y), c in joint.items():
        pxy = c / n
        px = xcount[x] / n
        py = ycount[y] / n
        mi += pxy * math.log2(pxy / (px * py))

    return mi


# ---------------------------------------------------------------------------
# Per-target measurement
# ---------------------------------------------------------------------------

def measure_target(
    target: Target,
    false_sums: list[int],
    cyclo: list[int],
    control: list[int],
    cyclo_chars: dict[int, list[int]],
    control_chars: dict[int, list[int]],
):
    data = {
        "cyclo_true": [],
        "cyclo_false": [],
        "control_true": [],
        "control_false": [],
    }

    # True signatures.
    for ell in cyclo:
        sig = cubic_ratio_signature(
            target.s,
            target.n,
            ell,
            cyclo_chars[ell],
        )

        data["cyclo_true"].append((ell, sig))

    for ell in control:
        sig = cubic_ratio_signature(
            target.s,
            target.n,
            ell,
            control_chars[ell],
        )

        data["control_true"].append((ell, sig))

    # False signatures.
    for s in false_sums:
        for ell in cyclo:
            sig = cubic_ratio_signature(
                s,
                target.n,
                ell,
                cyclo_chars[ell],
            )

            if sig is not None:
                data["cyclo_false"].append((ell, sig))

        for ell in control:
            sig = cubic_ratio_signature(
                s,
                target.n,
                ell,
                control_chars[ell],
            )

            if sig is not None:
                data["control_false"].append((ell, sig))

    return data


def summarize_family(
    true_data,
    false_data,
):
    true_ratio = []
    false_ratio = []

    true_ordered = []
    false_ordered = []

    for _, sig in true_data:
        if sig is None:
            continue

        true_ratio.append(sig["ratio_class"])
        true_ordered.append(sig["ordered_pair"])

    for _, sig in false_data:
        if sig is None:
            continue

        false_ratio.append(sig["ratio_class"])
        false_ordered.append(sig["ordered_pair"])

    true_nontrivial = sum(x != 0 for x in true_ratio)
    false_nontrivial = sum(x != 0 for x in false_ratio)

    return {
        "true_ratio": true_ratio,
        "false_ratio": false_ratio,
        "true_nontrivial_rate": (
            true_nontrivial / len(true_ratio)
            if true_ratio else 0.0
        ),
        "false_nontrivial_rate": (
            false_nontrivial / len(false_ratio)
            if false_ratio else 0.0
        ),
        "delta_nontrivial": (
            (true_nontrivial / len(true_ratio))
            - (false_nontrivial / len(false_ratio))
            if true_ratio and false_ratio else 0.0
        ),
        "true_ratio_entropy": entropy(Counter(true_ratio)),
        "false_ratio_entropy": entropy(Counter(false_ratio)),
        "true_ordered_mi": mutual_information_3x3(true_ordered),
        "false_ordered_mi": mutual_information_3x3(false_ordered),
        "true_count": len(true_ratio),
        "false_count": len(false_ratio),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 59")
    print("CUBIC-CHARACTER ROOT-RATIO TEST")
    print("TRUE FACTOR SUM VS EXACT MATCHED D-QR FALSE NULL")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. PRIME POPULATION
    # ----------------------------------------------------------------------

    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LO, PRIME_HI)
    prime_time = time.perf_counter() - t0

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_time:.6f}s")

    # ----------------------------------------------------------------------
    # 2. TARGETS
    # ----------------------------------------------------------------------

    targets = generate_targets(primes, TARGET_COUNT, rng)

    print()
    print("2. TARGETS")
    print("-" * 78)

    for i, t in enumerate(targets, 1):
        print(
            f"target {i:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    # ----------------------------------------------------------------------
    # 3. CUBIC MODULI
    # ----------------------------------------------------------------------

    cyclo = [x for x in CYCLOTOMIC if x % 3 == 1]
    control = [x for x in CONTROL if x % 3 == 1]

    print()
    print("3. CUBIC-CHARACTER MODULI")
    print("-" * 78)
    print(f"cyclotomic = {cyclo}")
    print(f"control    = {control}")

    # ----------------------------------------------------------------------
    # 4. BUILD CUBIC CHARACTER TABLES
    # ----------------------------------------------------------------------

    print()
    print("4. CUBIC CHARACTER TABLES")
    print("-" * 78)

    t0 = time.perf_counter()

    cyclo_chars = {}
    control_chars = {}

    for ell in cyclo:
        char, zeta = build_cubic_character_table(ell)
        cyclo_chars[ell] = char

        print(
            f"cyclo ell={ell:5d} "
            f"zeta={zeta:5d} "
            f"zeta^2={(zeta*zeta)%ell:5d}"
        )

    for ell in control:
        char, zeta = build_cubic_character_table(ell)
        control_chars[ell] = char

    print(f"table preparation = {time.perf_counter() - t0:.6f}s")

    # ----------------------------------------------------------------------
    # 5. ROOT-ALGEBRA VALIDATION
    # ----------------------------------------------------------------------

    print()
    print("5. TRUE-FACTOR CUBIC CHARACTER VALIDATION")
    print("-" * 78)

    validation_failures = 0
    true_cyclo_signatures = []

    for ti, target in enumerate(targets, 1):
        print(f"TARGET {ti:3d}")

        for ell in cyclo:
            sig = cubic_ratio_signature(
                target.s,
                target.n,
                ell,
                cyclo_chars[ell],
            )

            # Direct p,q check.
            pmod = target.p % ell
            qmod = target.q % ell

            if pmod == 0 or qmod == 0:
                expected = None
            else:
                cp = cyclo_chars[ell][pmod]
                cq = cyclo_chars[ell][qmod]
                expected = (cp - cq) % 3

            actual = None if sig is None else sig["ratio_class"]

            ok = (actual == expected)

            if not ok:
                validation_failures += 1

            print(
                f"  ell={ell:5d} "
                f"roots=({pmod:5d},{qmod:5d}) "
                f"ratio_class={str(actual):>2} "
                f"expected={str(expected):>2} "
                f"ok={ok}"
            )

            if sig is not None:
                true_cyclo_signatures.append(sig["ratio_class"])

    print(f"validation failures = {validation_failures}")
    print(f"status = {'PASS' if validation_failures == 0 else 'FAIL'}")

    # ----------------------------------------------------------------------
    # 6. EXACT MATCHED FALSE POOLS
    # ----------------------------------------------------------------------

    print()
    print("6. EXACT MATCHED FALSE-POOL CONSTRUCTION")
    print("-" * 78)
    print(
        "FALSE sums satisfy exactly the same cyclotomic "
        "D-QR conditions as the true sum."
    )

    # QR tables.
    qr_tables = {
        ell: qr_table(ell)
        for ell in cyclo
    }

    all_results = []

    t_false = time.perf_counter()

    for ti, target in enumerate(targets, 1):
        pool = build_exact_false_pool(
            target.n,
            target.s,
            qr_tables,
        )

        if len(pool) > FALSE_SAMPLE_CAP:
            false_sums = rng.sample(pool, FALSE_SAMPLE_CAP)
        else:
            false_sums = pool[:]

        print(
            f"target {ti:3d}: "
            f"full_false_pool={len(pool):5d} "
            f"sampled={len(false_sums):4d}"
        )

        all_results.append((target, false_sums))

    print(
        f"false-pool construction time = "
        f"{time.perf_counter() - t_false:.6f}s"
    )

    # ----------------------------------------------------------------------
    # 7. SIGNATURE MEASUREMENTS
    # ----------------------------------------------------------------------

    print()
    print("7. CUBIC ROOT-RATIO SIGNATURES")
    print("-" * 78)

    cyclo_true_all = []
    cyclo_false_all = []
    control_true_all = []
    control_false_all = []

    per_target = []

    for ti, (target, false_sums) in enumerate(all_results, 1):
        data = measure_target(
            target,
            false_sums,
            cyclo,
            control,
            cyclo_chars,
            control_chars,
        )

        csummary = summarize_family(
            data["cyclo_true"],
            data["cyclo_false"],
        )

        rsummary = summarize_family(
            data["control_true"],
            data["control_false"],
        )

        per_target.append((ti, csummary, rsummary))

        cyclo_true_all.extend(csummary["true_ratio"])
        cyclo_false_all.extend(csummary["false_ratio"])

        control_true_all.extend(rsummary["true_ratio"])
        control_false_all.extend(rsummary["false_ratio"])

        print()
        print(f"TARGET {ti:3d}")
        print(
            f"  cyclotomic "
            f"true_nontrivial={csummary['true_nontrivial_rate']:.6f} "
            f"false_nontrivial={csummary['false_nontrivial_rate']:.6f} "
            f"delta={csummary['delta_nontrivial']:+.6f}"
        )
        print(
            f"           "
            f"true_entropy={csummary['true_ratio_entropy']:.6f} "
            f"false_entropy={csummary['false_ratio_entropy']:.6f}"
        )
        print(
            f"  control    "
            f"true_nontrivial={rsummary['true_nontrivial_rate']:.6f} "
            f"false_nontrivial={rsummary['false_nontrivial_rate']:.6f} "
            f"delta={rsummary['delta_nontrivial']:+.6f}"
        )

    # ----------------------------------------------------------------------
    # 8. GLOBAL DISTRIBUTIONS
    # ----------------------------------------------------------------------

    print()
    print("8. GLOBAL CUBIC-CHARACTER DISTRIBUTIONS")
    print("-" * 78)

    for name, values in [
        ("CYCLOTOMIC TRUE", cyclo_true_all),
        ("CYCLOTOMIC FALSE", cyclo_false_all),
        ("CONTROL TRUE", control_true_all),
        ("CONTROL FALSE", control_false_all),
    ]:
        c = Counter(values)

        print()
        print(name)
        print(f"  class 0 = {c[0]:8d}")
        print(f"  class 1 = {c[1]:8d}")
        print(f"  class 2 = {c[2]:8d}")
        print(f"  total   = {len(values):8d}")

        if values:
            print(
                f"  nontrivial rate = "
                f"{sum(v != 0 for v in values)/len(values):.6f}"
            )
            print(
                f"  entropy = "
                f"{entropy(c):.6f} bits"
            )

    # ----------------------------------------------------------------------
    # 9. TARGET-LEVEL SEPARATION
    # ----------------------------------------------------------------------

    print()
    print("9. TARGET-LEVEL SEPARATION SUMMARY")
    print("-" * 78)

    c_deltas = [c["delta_nontrivial"] for _, c, _ in per_target]
    r_deltas = [r["delta_nontrivial"] for _, _, r in per_target]

    print(
        f"cyclotomic mean delta = {mean(c_deltas):+.6f}"
    )
    print(
        f"cyclotomic median delta = "
        f"{sorted(c_deltas)[len(c_deltas)//2]:+.6f}"
    )
    print(
        f"control mean delta = {mean(r_deltas):+.6f}"
    )
    print(
        f"control median delta = "
        f"{sorted(r_deltas)[len(r_deltas)//2]:+.6f}"
    )

    print()
    print("Per-target deltas:")
    for ti, c, r in per_target:
        print(
            f"  target {ti:3d}: "
            f"C={c['delta_nontrivial']:+.6f} "
            f"R={r['delta_nontrivial']:+.6f}"
        )

    # ----------------------------------------------------------------------
    # 10. POWER-3 SIGNATURE
    # ----------------------------------------------------------------------

    print()
    print("10. CUBIC CHARACTER CLASS HISTOGRAM")
    print("-" * 78)

    c_true = Counter(cyclo_true_all)
    c_false = Counter(cyclo_false_all)

    print("cyclotomic true:")
    print(
        f"  class 0: {c_true[0]:8d} "
        f"class 1: {c_true[1]:8d} "
        f"class 2: {c_true[2]:8d}"
    )

    print("cyclotomic false:")
    print(
        f"  class 0: {c_false[0]:8d} "
        f"class 1: {c_false[1]:8d} "
        f"class 2: {c_false[2]:8d}"
    )

    # ----------------------------------------------------------------------
    # 11. EXACT CONJUGATION CHECK
    # ----------------------------------------------------------------------

    print()
    print("11. SWAP / CONJUGATE CHECK")
    print("-" * 78)

    swap_failures = 0

    for target, _ in all_results:
        for ell in cyclo:
            sig1 = cubic_ratio_signature(
                target.s,
                target.n,
                ell,
                cyclo_chars[ell],
            )

            if sig1 is None:
                continue

            x, y = sig1["roots"]

            # Explicit swapped ratio.
            swapped = (y * pow(x, -1, ell)) % ell
            swapped_class = cyclo_chars[ell][swapped]

            expected = (-sig1["ratio_class"]) % 3

            if swapped_class != expected:
                swap_failures += 1

    print(f"swap failures = {swap_failures}")
    print(f"status = {'PASS' if swap_failures == 0 else 'FAIL'}")

    # ----------------------------------------------------------------------
    # 12. FINAL DIAGNOSTIC
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. FINAL DIAGNOSTIC")
    print("=" * 78)

    c_true_rate = (
        sum(v != 0 for v in cyclo_true_all) / len(cyclo_true_all)
        if cyclo_true_all else 0.0
    )

    c_false_rate = (
        sum(v != 0 for v in cyclo_false_all) / len(cyclo_false_all)
        if cyclo_false_all else 0.0
    )

    r_true_rate = (
        sum(v != 0 for v in control_true_all) / len(control_true_all)
        if control_true_all else 0.0
    )

    r_false_rate = (
        sum(v != 0 for v in control_false_all) / len(control_false_all)
        if control_false_all else 0.0
    )

    c_delta = c_true_rate - c_false_rate
    r_delta = r_true_rate - r_false_rate

    print()
    print("GLOBAL RATIO-CHARACTER RATES")
    print(
        f"cyclotomic true  nontrivial = {c_true_rate:.8f}"
    )
    print(
        f"cyclotomic false nontrivial = {c_false_rate:.8f}"
    )
    print(
        f"cyclotomic delta            = {c_delta:+.8f}"
    )
    print(
        f"control true  nontrivial     = {r_true_rate:.8f}"
    )
    print(
        f"control false nontrivial     = {r_false_rate:.8f}"
    )
    print(
        f"control delta                = {r_delta:+.8f}"
    )

    print()
    print("Interpretation")
    print("--------------")
    print("NO EFFECT:")
    print("  True and matched-false cubic ratio classes agree.")
    print()
    print("GENERIC EFFECT:")
    print("  True/false separation appears similarly in controls.")
    print()
    print("CYCLOTOMIC EFFECT:")
    print(
        "  True/false separation is materially stronger for "
        "cyclotomic ell than for matched random controls."
    )
    print()
    print("Important:")
    print(
        "  A difference in cubic-character distribution is not "
        "automatically a factorization algorithm."
    )
    print(
        "  It becomes genuinely interesting only if the effect "
        "survives larger target sets and changes the candidate "
        "set while preserving the true sum."
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 59 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

