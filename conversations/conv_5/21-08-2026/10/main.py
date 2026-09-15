#!/usr/bin/env python3

"""
======================================================================
SAT FACTORIZATION / 2^z CONSTRAINT EXPERIMENT
======================================================================

Baseline:
    solve p*q = n

Structured:
    solve p*q = n
    plus

        p ≡  a       (mod 2^(z-1))
        q ≡ -a       (mod 2^(z-1))

    where

        n ≡ -a²      (mod 2^z)

The experiment compares:

    1. baseline factorization
    2. one structured branch for every admissible a
    3. number of SAT variables
    4. solving time
    5. recovered factor pair

Requires:
    pip install z3-solver
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass

from z3 import (
    Solver,
    BitVec,
    BitVecVal,
    sat,
)


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    5959,       # 59 * 101
    10403,      # 101 * 103
    18721,      # 113 *  ?? if not factorable, still test
    3233,       # 53 * 61
    100127 * 100129,
]

# z values to test.
Z_VALUES = [2, 3, 4, 5, 6, 7, 8, 10]

# Maximum number of z bits to expose.
# Increase gradually for larger experiments.
MAX_Z = 16

# Solver timeout in milliseconds.
TIMEOUT_MS = 30_000


# ======================================================================
# DATA
# ======================================================================

@dataclass
class Result:
    mode: str
    n: int
    z: int | None
    a: int | None
    elapsed: float
    sat: bool
    p: int | None
    q: int | None
    correct: bool


# ======================================================================
# HELPER: BIT WIDTH
# ======================================================================

def bit_width(n: int) -> int:
    return max(2, n.bit_length())


# ======================================================================
# FIND ADMISSIBLE a
# ======================================================================

def admissible_a_values(n: int, z: int) -> list[int]:
    """
    Find a such that

        n ≡ -a² mod 2^z.

    Only a modulo 2^(z-1) matters for the p/q constraints.
    """

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    target = n % modulus

    result = []

    for a in range(a_modulus):
        if (-a * a) % modulus == target:
            result.append(a)

    return result


# ======================================================================
# CREATE MULTIPLICATION CONSTRAINT
# ======================================================================

def multiplication_solver(n: int):
    """
    Build:

        p*q = n

    using fixed-width bit-vectors.
    """

    width = bit_width(n)

    p = BitVec("p", width)
    q = BitVec("q", width)

    solver = Solver()

    # Optional timeout.
    solver.set(timeout=TIMEOUT_MS)

    # Exact multiplication modulo 2^width.
    #
    # To avoid overflow ambiguity, extend by width bits.
    p_ext = p.zero_ext(width)
    q_ext = q.zero_ext(width)

    n_ext = BitVecVal(n, 2 * width)

    solver.add(p_ext * q_ext == n_ext)

    # Non-trivial factors.
    solver.add(p > 1)
    solver.add(q > 1)

    return solver, p, q, width


# ======================================================================
# BASELINE
# ======================================================================

def solve_baseline(n: int) -> Result:

    solver, p, q, width = multiplication_solver(n)

    start = time.perf_counter()

    status = solver.check()

    elapsed = time.perf_counter() - start

    if status != sat:
        return Result(
            mode="BASELINE",
            n=n,
            z=None,
            a=None,
            elapsed=elapsed,
            sat=False,
            p=None,
            q=None,
            correct=False,
        )

    model = solver.model()

    p_value = model[p].as_long()
    q_value = model[q].as_long()

    correct = (
        p_value * q_value == n
        and p_value > 1
        and q_value > 1
    )

    return Result(
        mode="BASELINE",
        n=n,
        z=None,
        a=None,
        elapsed=elapsed,
        sat=True,
        p=p_value,
        q=q_value,
        correct=correct,
    )


# ======================================================================
# STRUCTURED SAT
# ======================================================================

def solve_structured(
    n: int,
    z: int,
    a: int,
) -> Result:

    solver, p, q, width = multiplication_solver(n)

    low_modulus = 1 << (z - 1)

    p_residue = a % low_modulus
    q_residue = (-a) % low_modulus

    # Force low bits.
    #
    # Extracting z-1 bits directly is cleaner than writing
    # a large collection of Boolean clauses.
    p_low = p & (low_modulus - 1)
    q_low = q & (low_modulus - 1)

    solver.add(p_low == p_residue)
    solver.add(q_low == q_residue)

    start = time.perf_counter()

    status = solver.check()

    elapsed = time.perf_counter() - start

    if status != sat:
        return Result(
            mode="STRUCTURED",
            n=n,
            z=z,
            a=a,
            elapsed=elapsed,
            sat=False,
            p=None,
            q=None,
            correct=False,
        )

    model = solver.model()

    p_value = model[p].as_long()
    q_value = model[q].as_long()

    correct = (
        p_value * q_value == n
        and p_value > 1
        and q_value > 1
        and p_value % low_modulus == p_residue
        and q_value % low_modulus == q_residue
    )

    return Result(
        mode="STRUCTURED",
        n=n,
        z=z,
        a=a,
        elapsed=elapsed,
        sat=True,
        p=p_value,
        q=q_value,
        correct=correct,
    )


# ======================================================================
# PRINT RESULT
# ======================================================================

def print_result(r: Result):

    if not r.sat:
        print(
            f"{r.mode:10s} "
            f"n={r.n:<12d} "
            f"z={str(r.z):>3s} "
            f"a={str(r.a):>5s} "
            f"UNSAT/TIMEOUT "
            f"time={r.elapsed:.6f}s"
        )
        return

    print(
        f"{r.mode:10s} "
        f"n={r.n:<12d} "
        f"z={str(r.z):>3s} "
        f"a={str(r.a):>5s} "
        f"p={r.p:<10d} "
        f"q={r.q:<10d} "
        f"time={r.elapsed:.6f}s "
        f"correct={r.correct}"
    )


# ======================================================================
# SINGLE N EXPERIMENT
# ======================================================================

def experiment_n(n: int):

    print()
    print("=" * 100)
    print(f"N = {n}")
    print("=" * 100)

    # --------------------------------------------------------------
    # Baseline
    # --------------------------------------------------------------

    baseline = solve_baseline(n)

    print()
    print("BASELINE")
    print("-" * 100)

    print_result(baseline)

    if not baseline.sat:
        print("Baseline did not produce a factor pair.")
        return

    actual_p = baseline.p
    actual_q = baseline.q

    # Normalize only for display.
    true_small = min(actual_p, actual_q)
    true_large = max(actual_p, actual_q)

    print()
    print(
        f"Recovered factorization: "
        f"{true_small} × {true_large} = {n}"
    )

    # --------------------------------------------------------------
    # Structured experiments
    # --------------------------------------------------------------

    print()
    print("STRUCTURED 2^z EXPERIMENT")
    print("-" * 100)

    for z in Z_VALUES:

        if z >= n.bit_length():
            continue

        modulus = 1 << z
        target = n % modulus

        a_values = admissible_a_values(n, z)

        print()
        print(
            f"z={z:2d} "
            f"2^z={modulus:<6d} "
            f"n mod 2^z={target:<6d} "
            f"admissible a={a_values}"
        )

        if not a_values:
            print("  NO admissible a")
            continue

        for a in a_values:

            r = solve_structured(
                n=n,
                z=z,
                a=a,
            )

            print_result(r)


# ======================================================================
# VERIFY THE THEORY
# ======================================================================

def verify_congruence():

    print()
    print("=" * 100)
    print("THEORETICAL CONGRUENCE VERIFICATION")
    print("=" * 100)

    failures = 0

    for z in range(2, MAX_Z + 1):

        modulus = 1 << z
        k = 1 << (z - 1)

        for a in range(0, min(modulus, 256)):

            expected = (-a * a) % modulus

            for x in range(-5, 6):
                for y in range(-5, 6):

                    p = k * (y - x) + a
                    q = k * (y + x) - a

                    n = p * q

                    if n % modulus != expected:
                        failures += 1

                        print(
                            "FAILURE:",
                            z,
                            a,
                            x,
                            y,
                            p,
                            q,
                            n,
                            expected,
                        )

    print()
    print(f"Failures = {failures}")
    print(f"PASS     = {failures == 0}")


# ======================================================================
# RESIDUE TABLE
# ======================================================================

def residue_table(z: int):

    modulus = 1 << z

    print()
    print("=" * 100)
    print(f"REACHABLE RESIDUES FOR z={z}")
    print("=" * 100)

    mapping = {}

    for a in range(modulus):

        residue = (-a * a) % modulus

        mapping.setdefault(residue, []).append(a)

    print()

    for residue in sorted(mapping):

        print(
            f"n ≡ {residue:>5} (mod {modulus:>5})"
            f"  <- a = {mapping[residue]}"
        )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 100)
    print("2^z STRUCTURED SAT FACTORIZATION EXPERIMENT")
    print("=" * 100)

    verify_congruence()

    # Show reachable classes.
    for z in [3, 4, 5, 6, 7, 8]:
        residue_table(z)

    # Test actual factorization.
    for n in N_VALUES:
        experiment_n(n)


if __name__ == "__main__":
    main()
