#!/usr/bin/env python3

"""
==============================================================================
KAPPA 50-BIT AUXILIARY / INFORMATION-GROWTH EXPERIMENT
==============================================================================

Core definitions:

    F(x) = x^2 - x + 1

    A = F(p)F(q)

For auxiliary primes r_i:

    B = product F(r_i)

and

    K = 1 - A*B.

Exact recovery:

    A = (1-K)/B

The experiment asks:

    1. How quickly does B reach ~50 bits?
    2. Does increasing B create new information about A?
    3. What happens if we observe K modulo small moduli?
    4. Can several auxiliary sequences accumulate information about A?
    5. Does CRT accumulation reduce the possible A/factor candidates?
    6. Is there any evidence that the auxiliary construction gives
       an n-only restriction, rather than merely re-encoding A?

IMPORTANT:

Exact K is allowed in the "oracle" experiments.

The n-only experiments explicitly distinguish this from information
that can only be obtained after knowing A/K.
"""

from __future__ import annotations

import csv
import math
import os
import random
from collections import Counter
from dataclasses import dataclass
from typing import Iterable


# ============================================================================
# CONFIGURATION
# ============================================================================

OUTDIR = "kappa_50_results"

AUX_TARGET_BITS = 50

TARGET_N_BITS = [20, 30, 40, 50]

# Number of semiprime examples for each target size.
TARGETS_PER_SIZE = 8

# Small moduli used for modular-information experiments.
MODULI = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47
]

# Maximum number of CRT steps shown per experiment.
MAX_CRT_STEPS = 14

RANDOM_SEED = 20260812


# ============================================================================
# BASIC MATH
# ============================================================================

def F(x: int) -> int:
    return x * x - x + 1


def bitlen(x: int) -> int:
    if x == 0:
        return 0
    return abs(x).bit_length()


def product(values: Iterable[int]) -> int:
    result = 1
    for x in values:
        result *= x
    return result


def gcd_many(a: int, b: int) -> int:
    return math.gcd(a, b)


def valuation_information(modulus: int) -> float:
    return math.log2(modulus)


# ============================================================================
# PRIME GENERATION
# ============================================================================

def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                ((limit - p * p) // p) + 1
            )

    return [i for i in range(2, limit + 1) if sieve[i]]


PRIMES = primes_up_to(10000)


# ============================================================================
# SEMIPRIME GENERATION
# ============================================================================

def is_prime(n: int) -> bool:
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


def choose_semiprime(bits: int, rng: random.Random) -> tuple[int, int, int]:
    """
    Choose p,q so that n=pq is approximately `bits` bits.
    """

    lower = 1 << (bits // 2 - 1)
    upper = 1 << (bits // 2 + 1)

    candidates = [
        p for p in PRIMES
        if lower <= p <= upper
    ]

    # For 50-bit experiments PRIMES is insufficient for p,q near 25 bits,
    # so fall back to generated probable primes.
    if not candidates:
        return random_semiprime(bits, rng)

    for _ in range(1000):
        p = rng.choice(candidates)
        q = rng.choice(candidates)

        if p == q:
            continue

        n = p * q

        if bitlen(n) == bits:
            return p, q, n

    return random_semiprime(bits, rng)


def random_odd_candidate(rng: random.Random, bits: int) -> int:
    x = rng.getrandbits(bits)
    x |= (1 << (bits - 1))
    x |= 1
    return x


def probable_prime(n: int, rounds: int = 12) -> bool:
    """
    Miller-Rabin for our experimental sizes.
    """

    if n < 2:
        return False

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # Deterministic enough for our experimental range with these bases.
    bases = [
        2, 325, 9375, 28178, 450775,
        9780504, 1795265022
    ]

    for a in bases[:rounds]:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x in (1, n - 1):
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(bits: int, rng: random.Random) -> int:
    while True:
        candidate = random_odd_candidate(rng, bits)

        if probable_prime(candidate):
            return candidate


def random_semiprime(bits: int, rng: random.Random) -> tuple[int, int, int]:
    """
    Generate p*q with approximately equal-size prime factors.
    """

    p_bits = bits // 2
    q_bits = bits - p_bits

    p = random_prime(p_bits, rng)
    q = random_prime(q_bits, rng)

    while q == p:
        q = random_prime(q_bits, rng)

    return p, q, p * q


# ============================================================================
# AUXILIARY SEQUENCES
# ============================================================================

@dataclass
class AuxiliarySequence:
    name: str
    rs: list[int]
    B: int


def build_auxiliary_sequence(
    name: str,
    target_bits: int,
    candidates: list[int],
) -> AuxiliarySequence:

    B = 1
    rs = []

    for r in candidates:
        rs.append(r)
        B *= F(r)

        if bitlen(B) >= target_bits:
            break

    return AuxiliarySequence(name, rs, B)


def build_sequences(target_bits: int) -> list[AuxiliarySequence]:

    sequences = []

    sequences.append(
        build_auxiliary_sequence(
            "first",
            target_bits,
            PRIMES,
        )
    )

    sequences.append(
        build_auxiliary_sequence(
            "odd_first",
            target_bits,
            [p for p in PRIMES if p != 2],
        )
    )

    sequences.append(
        build_auxiliary_sequence(
            "every_other",
            target_bits,
            PRIMES[::2],
        )
    )

    sequences.append(
        build_auxiliary_sequence(
            "larger_first",
            target_bits,
            [p for p in PRIMES if p >= 11],
        )
    )

    return sequences


# ============================================================================
# EXACT AUXILIARY EXPERIMENT
# ============================================================================

def exact_auxiliary_experiment(
    p: int,
    q: int,
    sequence: AuxiliarySequence,
) -> dict:

    A = F(p) * F(q)
    B = sequence.B

    K = 1 - A * B

    recovered_A = (1 - K) // B

    return {
        "p": p,
        "q": q,
        "n": p * q,
        "n_bits": bitlen(p * q),
        "A": A,
        "A_bits": bitlen(A),
        "B": B,
        "B_bits": bitlen(B),
        "K": K,
        "K_bits": bitlen(K),
        "recovered_A": recovered_A,
        "success": recovered_A == A,
        "sequence": sequence.name,
        "r_count": len(sequence.rs),
    }


# ============================================================================
# MODULAR INFORMATION
# ============================================================================

def modular_A_residue(A: int, B: int, modulus: int) -> int | None:
    """
    Given K = 1 - A*B mod m, recover A mod m when possible.

    This simulates the information obtained from K mod m.
    """

    g = math.gcd(B, modulus)

    if g != 1:
        return None

    # Simulated observed K.
    K_mod = (1 - A * B) % modulus

    inv_B = pow(B, -1, modulus)

    return ((1 - K_mod) * inv_B) % modulus


def candidate_residue_survival(
    n: int,
    A_true: int,
    modulus: int,
) -> tuple[int, int]:

    """
    Enumerate factor pairs d*(n/d)=n and count how many produce
    the same A modulo modulus.

    This is deliberately an experimental information measure.

    It is NOT intended as a factoring algorithm.
    """

    candidates = []

    d = 1

    while d * d <= n:
        if n % d == 0:
            q = n // d

            A = F(d) * F(q)

            candidates.append((d, q, A % modulus))

        d += 1

    true_residue = A_true % modulus

    survivors = sum(
        1 for _, _, residue in candidates
        if residue == true_residue
    )

    return len(candidates), survivors


# ============================================================================
# CRT
# ============================================================================

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> tuple[int, int] | None:

    g = math.gcd(m1, m2)

    if (a2 - a1) % g != 0:
        return None

    m1_g = m1 // g
    m2_g = m2 // g

    rhs = (a2 - a1) // g

    inv = pow(m1_g, -1, m2_g)

    t = (rhs * inv) % m2_g

    x = a1 + m1 * t
    modulus = m1 * m2_g

    return x % modulus, modulus


def crt_information_growth(
    A: int,
    B: int,
    moduli: list[int],
) -> list[dict]:

    residue = 0
    modulus = 1

    rows = []

    for i, m in enumerate(moduli):

        if i >= MAX_CRT_STEPS:
            break

        a_mod = modular_A_residue(A, B, m)

        if a_mod is None:
            rows.append({
                "modulus": m,
                "usable": False,
                "crt_modulus": modulus,
                "crt_bits": bitlen(modulus),
                "A_residue": "",
                "A_exact": False,
            })
            continue

        result = crt_pair(residue, modulus, a_mod, m)

        if result is None:
            continue

        residue, modulus = result

        exact = (
            modulus > A
            and residue == A
        )

        rows.append({
            "modulus": m,
            "usable": True,
            "crt_modulus": modulus,
            "crt_bits": bitlen(modulus),
            "A_residue": residue,
            "A_exact": exact,
        })

    return rows


# ============================================================================
# MULTI-SEQUENCE INFORMATION
# ============================================================================

def combined_residue_information(
    A: int,
    sequences: list[AuxiliarySequence],
    modulus: int,
) -> dict:

    residues = []

    for sequence in sequences:
        residue = modular_A_residue(A, sequence.B, modulus)

        if residue is not None:
            residues.append(residue)

    unique = sorted(set(residues))

    return {
        "modulus": modulus,
        "sequence_count": len(sequences),
        "usable_sequences": len(residues),
        "unique_A_residues": len(unique),
        "A_residue": A % modulus,
        "consistent": len(unique) == 1,
    }


# ============================================================================
# FIXED-n COLLISION ANALYSIS
# ============================================================================

def all_factor_pairs(n: int) -> list[tuple[int, int]]:
    pairs = []

    d = 2

    while d * d <= n:
        if n % d == 0:
            pairs.append((d, n // d))
        d += 1

    return pairs


def fixed_n_analysis(
    p: int,
    q: int,
    sequences: list[AuxiliarySequence],
    moduli: list[int],
) -> list[dict]:

    n = p * q
    A_true = F(p) * F(q)

    pairs = all_factor_pairs(n)

    rows = []

    for m in moduli:

        total = len(pairs)

        survivors = 0

        true_residue = A_true % m

        for a, b in pairs:
            candidate_A = F(a) * F(b)

            if candidate_A % m == true_residue:
                survivors += 1

        rows.append({
            "n": n,
            "n_bits": bitlen(n),
            "modulus": m,
            "factor_pairs": total,
            "survivors": survivors,
            "survival_fraction": (
                survivors / total if total else 1.0
            ),
        })

    return rows


# ============================================================================
# INFORMATION ESTIMATE
# ============================================================================

def information_bits(total: int, survivors: int) -> float:
    if survivors <= 0 or total <= 0:
        return 0.0

    return math.log2(total / survivors)


# ============================================================================
# CSV HELPERS
# ============================================================================

def ensure_output_dir():
    os.makedirs(OUTDIR, exist_ok=True)


def write_csv(filename: str, rows: list[dict]):

    path = os.path.join(OUTDIR, filename)

    if not rows:
        return

    keys = []

    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)

    with open(path, "w", newline="", encoding="utf-8") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=keys,
        )

        writer.writeheader()
        writer.writerows(rows)


# ============================================================================
# DISPLAY HELPERS
# ============================================================================

def print_rule(char="=", width=78):
    print(char * width)


def fmt_int(x: int) -> str:
    return f"{x:,}"


def print_sequence_summary(sequences):

    print()
    print("AUXILIARY SEQUENCES")
    print("-" * 78)

    print(
        f"{'sequence':<16}"
        f"{'count':>7}"
        f"{'B bits':>10}"
        f"{'first r':>10}"
        f"{'last r':>10}"
    )

    for seq in sequences:

        print(
            f"{seq.name:<16}"
            f"{len(seq.rs):>7}"
            f"{bitlen(seq.B):>10}"
            f"{seq.rs[0]:>10}"
            f"{seq.rs[-1]:>10}"
        )

    print()

    for seq in sequences:
        print(
            f"{seq.name:<16}: "
            f"r = {tuple(seq.rs)}"
        )


def print_exact_summary(rows):

    successes = sum(r["success"] for r in rows)

    print()
    print("EXACT A RECOVERY")
    print("-" * 78)

    print(
        f"tests={len(rows)} "
        f"successes={successes} "
        f"failures={len(rows)-successes}"
    )

    if rows:
        r = rows[0]

        print(
            f"example: n={r['n']} "
            f"A={r['A']} "
            f"B_bits={r['B_bits']} "
            f"K_bits={r['K_bits']}"
        )

    print(
        "RESULT: "
        "exact K plus known auxiliary B always recovers A."
    )


def print_crt_summary(A, B, rows):

    print()
    print("CRT INFORMATION GROWTH")
    print("-" * 78)

    print(
        f"A = {A:,} "
        f"(A bits={bitlen(A)})"
    )

    print(
        f"B bits={bitlen(B)}"
    )

    print()

    print(
        f"{'m':>6}"
        f"{'CRT bits':>12}"
        f"{'A residue':>24}"
        f"{'exact?':>10}"
    )

    for row in rows:

        residue = row["A_residue"]

        print(
            f"{row['modulus']:>6}"
            f"{row['crt_bits']:>12}"
            f"{str(residue):>24}"
            f"{str(row['A_exact']):>10}"
        )

    print()

    final = rows[-1] if rows else None

    if final:
        print(
            f"Final CRT modulus has "
            f"{final['crt_bits']} bits."
        )


def print_survival_summary(rows):

    print()
    print("FIXED-n MODULAR CANDIDATE SURVIVAL")
    print("-" * 78)

    print(
        f"{'m':>6}"
        f"{'pairs':>10}"
        f"{'survive':>10}"
        f"{'fraction':>12}"
        f"{'info bits':>12}"
    )

    for row in rows:

        info = information_bits(
            row["factor_pairs"],
            row["survivors"],
        )

        print(
            f"{row['modulus']:>6}"
            f"{row['factor_pairs']:>10}"
            f"{row['survivors']:>10}"
            f"{row['survival_fraction']:>12.4f}"
            f"{info:>12.4f}"
        )


def print_multi_sequence_summary(A, sequences):

    print()
    print("MULTI-AUXILIARY CONSISTENCY")
    print("-" * 78)

    print(
        f"{'mod':>6}"
        f"{'usable':>10}"
        f"{'unique residues':>18}"
        f"{'consistent?':>14}"
    )

    for m in MODULI:

        result = combined_residue_information(
            A,
            sequences,
            m,
        )

        print(
            f"{m:>6}"
            f"{result['usable_sequences']:>10}"
            f"{result['unique_A_residues']:>18}"
            f"{str(result['consistent']):>14}"
        )


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================

def main():

    rng = random.Random(RANDOM_SEED)

    ensure_output_dir()

    print_rule()
    print(
        "KAPPA 50-BIT AUXILIARY / INFORMATION-GROWTH EXPERIMENT"
    )
    print_rule()

    print()
    print("Configuration")
    print("-" * 78)
    print(f"Auxiliary target: ~{AUX_TARGET_BITS} bits")
    print(f"Target n sizes : {TARGET_N_BITS}")
    print(f"Targets/size   : {TARGETS_PER_SIZE}")
    print(f"Moduli         : {MODULI}")
    print()

    # ------------------------------------------------------------------
    # Auxiliary sequences
    # ------------------------------------------------------------------

    sequences = build_sequences(AUX_TARGET_BITS)

    print_sequence_summary(sequences)

    # ------------------------------------------------------------------
    # Generate target semiprimes
    # ------------------------------------------------------------------

    targets = []

    for bits in TARGET_N_BITS:

        for _ in range(TARGETS_PER_SIZE):

            p, q, n = choose_semiprime(bits, rng)

            targets.append({
                "p": p,
                "q": q,
                "n": n,
                "n_bits": bitlen(n),
                "A": F(p) * F(q),
                "A_bits": bitlen(F(p) * F(q)),
            })

    print()
    print("TARGETS")
    print("-" * 78)

    for target in targets:

        print(
            f"n_bits={target['n_bits']:>3} "
            f"n={target['n']} "
            f"A_bits={target['A_bits']:>3}"
        )

    # ------------------------------------------------------------------
    # Exact recovery
    # ------------------------------------------------------------------

    exact_rows = []

    for target in targets:

        for sequence in sequences:

            exact_rows.append(
                exact_auxiliary_experiment(
                    target["p"],
                    target["q"],
                    sequence,
                )
            )

    print_exact_summary(exact_rows)

    # ------------------------------------------------------------------
    # Select representative target for detailed experiments
    # ------------------------------------------------------------------

    representative = max(
        targets,
        key=lambda x: x["n_bits"]
    )

    p = representative["p"]
    q = representative["q"]
    n = representative["n"]
    A = representative["A"]

    print()
    print("REPRESENTATIVE TARGET")
    print("-" * 78)

    print(f"p = {p}")
    print(f"q = {q}")
    print(f"n = {n}")
    print(f"n bits = {bitlen(n)}")
    print(f"A = F(p)F(q) = {A}")
    print(f"A bits = {bitlen(A)}")

    # ------------------------------------------------------------------
    # CRT growth for each sequence
    # ------------------------------------------------------------------

    crt_rows = []

    print()
    print("CRT EXPERIMENTS")
    print("-" * 78)

    for sequence in sequences:

        rows = crt_information_growth(
            A,
            sequence.B,
            MODULI,
        )

        for row in rows:

            row["sequence"] = sequence.name
            row["B_bits"] = bitlen(sequence.B)
            row["A_bits"] = bitlen(A)

        crt_rows.extend(rows)

        print()
        print(
            f"Sequence: {sequence.name} "
            f"(B bits={bitlen(sequence.B)})"
        )

        print_crt_summary(
            A,
            sequence.B,
            rows,
        )

    # ------------------------------------------------------------------
    # Candidate survival
    #
    # Exact factor-pair enumeration is only practical for small n.
    # Therefore use the smaller targets for this experiment.
    # ------------------------------------------------------------------

    survival_rows = []

    small_targets = [
        t for t in targets
        if t["n_bits"] <= 30
    ]

    print()
    print("CANDIDATE SURVIVAL EXPERIMENTS")
    print("-" * 78)

    for target in small_targets[:4]:

        rows = fixed_n_analysis(
            target["p"],
            target["q"],
            sequences,
            MODULI,
        )

        survival_rows.extend(rows)

        print()
        print(
            f"n={target['n']} "
            f"(bits={target['n_bits']}) "
            f"p={target['p']} q={target['q']}"
        )

        print_survival_summary(rows)

    # ------------------------------------------------------------------
    # Multi-sequence consistency
    # ------------------------------------------------------------------

    multi_rows = []

    print_multi_sequence_summary(
        A,
        sequences,
    )

    for m in MODULI:

        result = combined_residue_information(
            A,
            sequences,
            m,
        )

        result["n"] = n
        result["A"] = A

        multi_rows.append(result)

    # ------------------------------------------------------------------
    # Auxiliary growth table
    # ------------------------------------------------------------------

    growth_rows = []

    for sequence in sequences:

        B = 1

        for i, r in enumerate(sequence.rs, start=1):

            B *= F(r)

            growth_rows.append({
                "sequence": sequence.name,
                "step": i,
                "r": r,
                "F_r": F(r),
                "B": B,
                "B_bits": bitlen(B),
            })

    print()
    print("AUXILIARY PRODUCT GROWTH")
    print("-" * 78)

    for sequence in sequences:

        rows = [
            r for r in growth_rows
            if r["sequence"] == sequence.name
        ]

        print(
            f"{sequence.name:<16}: "
            f"{len(rows)} primes -> "
            f"{rows[-1]['B_bits']} bits"
        )

    # ------------------------------------------------------------------
    # Target summary
    # ------------------------------------------------------------------

    target_summary_rows = []

    for target in targets:

        target_summary_rows.append({
            "n": target["n"],
            "n_bits": target["n_bits"],
            "A": target["A"],
            "A_bits": target["A_bits"],
        })

    # ------------------------------------------------------------------
    # Write CSVs
    # ------------------------------------------------------------------

    write_csv(
        "auxiliary_sequences.csv",
        [
            {
                "sequence": s.name,
                "r_count": len(s.rs),
                "rs": " ".join(map(str, s.rs)),
                "B": s.B,
                "B_bits": bitlen(s.B),
            }
            for s in sequences
        ],
    )

    write_csv(
        "auxiliary_growth.csv",
        growth_rows,
    )

    write_csv(
        "exact_recovery.csv",
        exact_rows,
    )

    write_csv(
        "crt_growth.csv",
        crt_rows,
    )

    write_csv(
        "candidate_survival.csv",
        survival_rows,
    )

    write_csv(
        "multi_sequence_information.csv",
        multi_rows,
    )

    write_csv(
        "target_summary.csv",
        target_summary_rows,
    )

    # ------------------------------------------------------------------
    # Final interpretation
    # ------------------------------------------------------------------

    print()
    print_rule()
    print("FINAL INTERPRETATION")
    print_rule()

    print()
    print("1. AUXILIARY SIZE")
    print("-" * 78)

    for sequence in sequences:

        print(
            f"{sequence.name:<16}: "
            f"{len(sequence.rs):>2} auxiliary primes, "
            f"B={bitlen(sequence.B):>3} bits"
        )

    print()
    print("2. EXACT INFORMATION")
    print("-" * 78)

    successes = sum(
        r["success"]
        for r in exact_rows
    )

    print(
        f"Exact A recovery: "
        f"{successes}/{len(exact_rows)}"
    )

    print(
        "This is expected: knowing exact K and B "
        "directly determines A."
    )

    print()
    print("3. CRT INFORMATION")
    print("-" * 78)

    print(
        f"A has {bitlen(A)} bits."
    )

    print(
        "CRT experiments measure how many residue bits "
        "can be accumulated about A."
    )

    print(
        "This is information about A, not automatically "
        "information obtainable from n."
    )

    print()
    print("4. CANDIDATE SURVIVAL")
    print("-" * 78)

    if survival_rows:

        for n_value in sorted(
            set(r["n"] for r in survival_rows)
        ):

            rows = [
                r for r in survival_rows
                if r["n"] == n_value
            ]

            total = rows[0]["factor_pairs"]

            best = min(
                rows,
                key=lambda r: r["survivors"]
            )

            info = information_bits(
                total,
                best["survivors"],
            )

            print(
                f"n={n_value}: "
                f"{total} factor pairs -> "
                f"best={best['survivors']} survivors "
                f"(~{info:.3f} bits candidate reduction)"
            )

    else:

        print(
            "No small-n candidate enumeration performed."
        )

    print()
    print("5. KEY QUESTION")
    print("-" * 78)

    print(
        "The experiment distinguishes:"
    )

    print(
        "    n + exact K + B  -> A"
    )

    print(
        "from:"
    )

    print(
        "    n alone          -> A ?"
    )

    print()
    print(
        "The first implication is algebraically certain."
    )

    print(
        "The second is the actual research question."
    )

    print()
    print(
        "If multiple auxiliary sequences only reproduce "
        "the same A residues, then they are multiple encodings "
        "of the same hidden quantity rather than independent "
        "sources of information."
    )

    print()
    print(
        f"Detailed CSV files written to: {OUTDIR}/"
    )

    print_rule()


if __name__ == "__main__":
    main()

