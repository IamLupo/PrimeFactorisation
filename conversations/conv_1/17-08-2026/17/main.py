#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 87
DIRECT H6/H8 RESOLVENT FOR s^2

NEW SYMBOLIC DISCOVERY:
    N6 = 6 H6
    N8 = 24 H8
    N8 = N6 * (s^2 - 3n)

THEREFORE:
    H8 / H6 = (s^2 - 3n) / 4

AND:
    s^2 = 3n + 4 H8/H6

EXACT ORACLE RECONSTRUCTION
MODULAR RESOLVENT
CYCLOTOMIC VS CONTROL
CRT RECOVERY OF s
STRICT TARGET HOLDOUT

NO CSV OUTPUT
NO SKLEARN
==============================================================================

IMPORTANT:
    This experiment uses H6/H8 values computed from the known p,q only
    to measure the algebraic information available in the paper's
    coefficients.

    The reconstruction step uses n plus H6/H8 values as inputs.

    It does NOT yet compute H6/H8 from n alone.

==============================================================================
"""

from __future__ import annotations

import math
import statistics
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

NUM_TARGETS = 40
TRAIN_TARGETS = 30

P_MIN = 2_000_000
P_MAX = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_400_000

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67
]

CONTROL = [
    673, 4561, 4759, 6211,
    7879, 7951, 8689
]

RNG_SEED = 87087


# ============================================================================
# TARGET
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> List[int]:

    lo = max(2, lo)

    sieve = bytearray(
        b"\x01" * (hi + 1)
    )

    sieve[0:2] = b"\x00\x00"

    for p in range(
        2,
        math.isqrt(hi) + 1,
    ):
        if sieve[p]:
            start = p * p
            sieve[
                start:
                hi + 1:
                p
            ] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x
        for x in range(lo, hi + 1)
        if sieve[x]
    ]


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:

    state = RNG_SEED
    mask = (1 << 64) - 1

    out: List[Target] = []
    seen = set()

    N = len(primes)

    while len(out) < count:

        state = (
            6364136223846793005 * state
            + 1442695040888963407
        ) & mask

        i = state % N

        state = (
            6364136223846793005 * state
            + 1442695040888963407
        ) & mask

        j = state % N

        if i == j:
            continue

        p = primes[i]
        q = primes[j]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        out.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return out


# ============================================================================
# PAPER COEFFICIENTS
# ============================================================================

def sigma_semiprime(
    n: int,
    s: int,
    j: int,
) -> int:

    """
    sigma_j(pq)
      = (1+p^j)(1+q^j)
      = 1 + (p^j+q^j) + n^j

    Uses only n,s.
    """

    R = [0] * (j + 1)

    R[0] = 2

    if j >= 1:
        R[1] = s

    for k in range(
        2,
        j + 1,
    ):
        R[k] = (
            s * R[k - 1]
            - n * R[k - 2]
        )

    return (
        1
        + R[j]
        + n ** j
    )


def H6(
    n: int,
    s: int,
) -> int:

    numerator = (
        (n * n - n + 1)
        * sigma_semiprime(n, s, 1)
        - sigma_semiprime(n, s, 3)
    )

    if numerator % 6 != 0:
        raise ArithmeticError(
            "H6 integrality failure"
        )

    return numerator // 6


def H8(
    n: int,
    s: int,
) -> int:

    numerator = (
        -n * n
        * sigma_semiprime(n, s, 1)
        + (n * n + 1)
        * sigma_semiprime(n, s, 3)
        - sigma_semiprime(n, s, 5)
    )

    if numerator % 24 != 0:
        raise ArithmeticError(
            "H8 integrality failure"
        )

    return numerator // 24


# ============================================================================
# NEW RESOLVENT
# ============================================================================

def reconstructed_s2(
    n: int,
    h6: int,
    h8: int,
) -> int | None:

    """
    From:

        s^2 = 3n + 4 H8/H6

    Require exact divisibility.
    """

    if h6 == 0:
        return None

    numerator = (
        3 * n * h6
        + 4 * h8
    )

    if numerator % h6 != 0:
        return None

    return numerator // h6


def reconstructed_s(
    n: int,
    h6: int,
    h8: int,
) -> int | None:

    s2 = reconstructed_s2(
        n,
        h6,
        h8,
    )

    if s2 is None:
        return None

    if s2 < 0:
        return None

    root = math.isqrt(s2)

    if root * root != s2:
        return None

    return root


# ============================================================================
# MODULAR H6/H8
# ============================================================================

def h6_h8_mod_from_ns(
    n: int,
    s: int,
    ell: int,
) -> Tuple[int, int]:

    n %= ell
    s %= ell

    R = [0] * 6

    R[0] = 2 % ell
    R[1] = s

    for j in range(
        2,
        6,
    ):
        R[j] = (
            s * R[j - 1]
            - n * R[j - 2]
        ) % ell

    def sigma(j: int) -> int:
        return (
            1
            + R[j]
            + pow(n, j, ell)
        ) % ell

    inv6 = pow(
        6,
        ell - 2,
        ell,
    )

    inv24 = pow(
        24,
        ell - 2,
        ell,
    )

    h6 = (
        (
            (n * n - n + 1)
            * sigma(1)
            - sigma(3)
        )
        * inv6
    ) % ell

    h8 = (
        (
            -n * n * sigma(1)
            + (n * n + 1) * sigma(3)
            - sigma(5)
        )
        * inv24
    ) % ell

    return (
        h6,
        h8,
    )


def modular_s2_from_h6h8(
    n: int,
    h6_mod: int,
    h8_mod: int,
    ell: int,
) -> int | None:

    if h6_mod % ell == 0:
        return None

    inv_h6 = pow(
        h6_mod,
        ell - 2,
        ell,
    )

    return (
        3 * (n % ell)
        + 4 * h8_mod * inv_h6
    ) % ell


# ============================================================================
# MODULAR SQUARE ROOTS
# ============================================================================

def modular_square_roots(
    a: int,
    ell: int,
) -> List[int]:

    """
    Small/moderate experimental moduli, so direct enumeration is fine.
    """

    a %= ell

    return [
        x
        for x in range(ell)
        if (x * x) % ell == a
    ]


# ============================================================================
# CRT
# ============================================================================

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> Tuple[int, int]:

    """
    Combine:
        x = a1 mod m1
        x = a2 mod m2

    with gcd(m1,m2)=1.
    """

    inv = pow(
        m1,
        -1,
        m2,
    )

    t = (
        (a2 - a1)
        * inv
    ) % m2

    x = (
        a1
        + m1 * t
    )

    return (
        x % (m1 * m2),
        m1 * m2,
    )


def crt_residue_sets(
    residue_sets: List[List[int]],
    moduli: Sequence[int],
) -> List[Tuple[int, int]]:

    states = [
        (0, 1)
    ]

    for residues, modulus in zip(
        residue_sets,
        moduli,
    ):

        new_states = []

        for a, m in states:

            for r in residues:

                x, M = crt_pair(
                    a,
                    m,
                    r,
                    modulus,
                )

                new_states.append(
                    (x, M)
                )

        # Deduplicate.
        seen = set()

        states = []

        for x, M in new_states:

            key = (
                x % M,
                M,
            )

            if key in seen:
                continue

            seen.add(key)
            states.append(key)

    return states


# ============================================================================
# EVEN S DOMAIN FROM CRT
# ============================================================================

def lift_even_residues(
    crt_states: Sequence[Tuple[int, int]],
    lo: int,
    hi: int,
) -> List[int]:

    out = []

    for a, M in crt_states:

        first = a

        if first < lo:

            k = (
                (lo - first + M - 1)
                // M
            )

            first += k * M

        if first > hi:
            continue

        # Preserve even s.
        if first % 2:
            first += M

        if first > hi:
            continue

        for s in range(
            first,
            hi + 1,
            M,
        ):

            if s % 2 == 0:
                out.append(s)

    return sorted(
        set(out)
    )


# ============================================================================
# TARGET VALIDATION
# ============================================================================

def validate_target(
    t: Target,
) -> bool:

    h6_direct = H6(
        t.n,
        t.s,
    )

    h8_direct = H8(
        t.n,
        t.s,
    )

    # Direct factor formula.
    e1_direct = (
        t.s
        * (
            (t.n + 1) ** 2
            - t.s ** 2
        )
    )

    if e1_direct != 6 * h6_direct:
        return False

    predicted_s2 = reconstructed_s2(
        t.n,
        h6_direct,
        h8_direct,
    )

    if predicted_s2 != t.s * t.s:
        return False

    predicted_s = reconstructed_s(
        t.n,
        h6_direct,
        h8_direct,
    )

    return predicted_s == t.s


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 87")
    print("DIRECT H6/H8 RESOLVENT FOR s^2")
    print("H8/H6 = (s^2 - 3n)/4")
    print("EXACT + MODULAR + CRT")
    print("CYCLOTOMIC VS CONTROL")
    print("STRICT TARGET HOLDOUT")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)

    print(
        f"prime population = "
        f"{len(primes)}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    for i, t in enumerate(
        targets[:24],
        1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if NUM_TARGETS > 24:
        print(
            "... remaining generated targets omitted"
        )

    # ------------------------------------------------------------------
    # Identity validation
    # ------------------------------------------------------------------

    print("\n2. DIRECT H6/H8 RESOLVENT VALIDATION")
    print("-" * 78)

    failures = 0

    for i, t in enumerate(
        targets,
        1,
    ):

        ok = validate_target(t)

        if not ok:
            failures += 1

        if i <= 24:

            h6 = H6(
                t.n,
                t.s,
            )

            h8 = H8(
                t.n,
                t.s,
            )

            s2 = reconstructed_s2(
                t.n,
                h6,
                h8,
            )

            print(
                f"target {i:3d}: "
                f"ok={ok} "
                f"H6={h6} "
                f"H8={h8} "
                f"reconstructed_s2={s2 == t.s*t.s}"
            )

    print(
        f"identity failures = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "H6/H8 resolvent identity failed"
        )

    print("status = PASS")

    # ------------------------------------------------------------------
    # Symbolic factorization printed explicitly
    # ------------------------------------------------------------------

    print("\n3. ALGEBRAIC FACTORIZATION")
    print("-" * 78)

    print(
        "6 H6 = s * ((n+1)^2 - s^2)"
    )

    print(
        "24 H8 = "
        "s * (s^2 - 3n) * ((n+1)^2 - s^2)"
    )

    print(
        "Therefore:"
    )

    print(
        "H8 / H6 = (s^2 - 3n) / 4"
    )

    print(
        "and:"
    )

    print(
        "s^2 = 3n + 4 H8/H6"
    )

    # ------------------------------------------------------------------
    # Holdout
    # ------------------------------------------------------------------

    train = targets[
        :TRAIN_TARGETS
    ]

    test = targets[
        TRAIN_TARGETS:
    ]

    print("\n4. TARGET HOLDOUT")
    print("-" * 78)

    print(
        f"training targets = "
        f"{len(train)}"
    )

    print(
        f"test targets = "
        f"{len(test)}"
    )

    # ------------------------------------------------------------------
    # Exact oracle reconstruction
    # ------------------------------------------------------------------

    print("\n5. EXACT ORACLE RECONSTRUCTION")
    print("-" * 78)

    exact_train = 0
    exact_test = 0

    h6_zero = 0

    for t in train:

        h6 = H6(
            t.n,
            t.s,
        )

        h8 = H8(
            t.n,
            t.s,
        )

        if h6 == 0:
            h6_zero += 1
            continue

        recovered = reconstructed_s(
            t.n,
            h6,
            h8,
        )

        if recovered == t.s:
            exact_train += 1

    for t in test:

        h6 = H6(
            t.n,
            t.s,
        )

        h8 = H8(
            t.n,
            t.s,
        )

        if h6 == 0:
            h6_zero += 1
            continue

        recovered = reconstructed_s(
            t.n,
            h6,
            h8,
        )

        if recovered == t.s:
            exact_test += 1

    print(
        f"train exact recovery = "
        f"{exact_train}/{len(train)}"
    )

    print(
        f"test exact recovery = "
        f"{exact_test}/{len(test)}"
    )

    print(
        f"H6=0 targets = "
        f"{h6_zero}"
    )

    # ------------------------------------------------------------------
    # Modular experiment
    # ------------------------------------------------------------------

    print("\n6. MODULAR H6/H8 RESOLVENT")
    print("-" * 78)

    families = [
        (
            "C3",
            CYCLOTOMIC[:3],
        ),
        (
            "C5",
            CYCLOTOMIC[:5],
        ),
        (
            "C7",
            CYCLOTOMIC,
        ),
        (
            "R3",
            CONTROL[:3],
        ),
        (
            "R5",
            CONTROL[:5],
        ),
        (
            "R7",
            CONTROL,
        ),
    ]

    for name, moduli in families:

        print("\n")
        print(name)
        print(
            f"moduli = {list(moduli)}"
        )
        print("-" * 78)

        total_modular_classes = 0
        recovered_targets = 0
        unique_targets = 0
        degenerate_targets = 0

        for target in test:

            residue_sets = []

            degenerate = False

            for ell in moduli:

                h6m, h8m = (
                    h6_h8_mod_from_ns(
                        target.n,
                        target.s,
                        ell,
                    )
                )

                s2m = (
                    modular_s2_from_h6h8(
                        target.n,
                        h6m,
                        h8m,
                        ell,
                    )
                )

                if s2m is None:

                    degenerate = True
                    break

                roots = modular_square_roots(
                    s2m,
                    ell,
                )

                if not roots:

                    degenerate = True
                    break

                residue_sets.append(
                    roots
                )

            if degenerate:

                degenerate_targets += 1
                continue

            crt_states = crt_residue_sets(
                residue_sets,
                moduli,
            )

            candidates = lift_even_residues(
                crt_states,
                S_MIN,
                S_MAX,
            )

            total_modular_classes += (
                len(candidates)
            )

            if target.s in candidates:
                recovered_targets += 1

            if len(candidates) == 1:
                unique_targets += 1

        usable = (
            len(test)
            - degenerate_targets
        )

        if usable:

            mean_candidates = (
                total_modular_classes
                / usable
            )

        else:
            mean_candidates = float(
                "nan"
            )

        print(
            f"usable targets = "
            f"{usable}/{len(test)}"
        )

        print(
            f"mean s candidates = "
            f"{mean_candidates:.3f}"
        )

        print(
            f"true-s recovered = "
            f"{recovered_targets}/{len(test)}"
        )

        print(
            f"unique recovery = "
            f"{unique_targets}/{len(test)}"
        )

        print(
            f"H6 modular degeneracy = "
            f"{degenerate_targets}/{len(test)}"
        )

    # ------------------------------------------------------------------
    # Modular root diagnostics
    # ------------------------------------------------------------------

    print("\n7. MODULAR ROOT DIAGNOSTICS")
    print("-" * 78)

    for ell in CYCLOTOMIC:

        total_roots = 0
        zero_h6 = 0
        test_count = 0

        for target in test:

            h6m, h8m = (
                h6_h8_mod_from_ns(
                    target.n,
                    target.s,
                    ell,
                )
            )

            if h6m == 0:
                zero_h6 += 1
                continue

            s2m = modular_s2_from_h6h8(
                target.n,
                h6m,
                h8m,
                ell,
            )

            if s2m is None:
                continue

            roots = modular_square_roots(
                s2m,
                ell,
            )

            total_roots += len(
                roots
            )

            test_count += 1

        avg = (
            total_roots / test_count
            if test_count
            else float("nan")
        )

        print(
            f"ell={ell:5d} "
            f"H6_zero={zero_h6:2d} "
            f"mean square roots={avg:.3f}"
        )

    # ------------------------------------------------------------------
    # Special possibility: solve s^2 directly
    # ------------------------------------------------------------------

    print("\n8. DIRECT s^2 DIAGNOSTIC")
    print("-" * 78)

    errors = []

    for t in test:

        h6 = H6(
            t.n,
            t.s,
        )

        h8 = H8(
            t.n,
            t.s,
        )

        s2 = reconstructed_s2(
            t.n,
            h6,
            h8,
        )

        if s2 != t.s * t.s:
            errors.append(
                t.s
            )

    print(
        f"test s^2 errors = "
        f"{len(errors)}/{len(test)}"
    )

    if not errors:
        print(
            "Every held-out target satisfies the direct H6/H8 "
            "resolvent exactly."
        )

    # ------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "NEW ALGEBRAIC RESULT:"
    )

    print(
        "    H8/H6 = (s^2 - 3n)/4"
    )

    print(
        "Therefore the pair (H6,H8) determines s^2 directly."
    )

    print()
    print(
        "This is fundamentally different from the earlier "
        "N-only classification experiments."
    )

    print(
        "Those experiments searched for a predictor:"
    )

    print(
        "    n -> hidden sign"
    )

    print(
        "This experiment finds a direct algebraic resolvent:"
    )

    print(
        "    (n,H6,H8) -> s^2"
    )

    print()
    print(
        "CRITICAL BOTTLENECK:"
    )

    print(
        "The remaining question is now sharply isolated:"
    )

    print(
        "    Can H6 and H8, or their ratio, be obtained "
        "from n alone?"
    )

    print()
    print(
        "If yes, then:"
    )

    print(
        "    s^2 = 3n + 4 H8/H6"
    )

    print(
        "would immediately reveal s because s>0."
    )

    print()
    print(
        "That would bypass the enormous s-search completely."
    )

    print()
    print(
        "CAVEAT:"
    )

    print(
        "The current experiment still uses oracle H6/H8 values "
        "constructed from p,q."
    )

    runtime = (
        time.perf_counter()
        - start
    )

    print(
        f"\ntotal runtime = "
        f"{runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 87 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

