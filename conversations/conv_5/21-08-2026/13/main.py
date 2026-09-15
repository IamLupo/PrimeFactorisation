#!/usr/bin/env python3

"""
======================================================================
2^z PARAMETERIZATION / PRIME FACTOR ENUMERATION EXPERIMENT
======================================================================

Family:

    k = 2^(z-1)

    p = k(y-x) + a
    q = k(y+x) - a

Therefore:

    pq ≡ -a^2 (mod 2^z)

Full parameterization constraints:

    p + q ≡ 0       (mod 2^z)
    q - p ≡ -2a     (mod 2^z)

This experiment does NOT stop at the first SAT solution.

For every (z, a) branch it:

    1. Enumerates SAT factor pairs.
    2. Tests p and q for primality.
    3. Blocks the discovered pair.
    4. Continues until UNSAT / timeout / candidate limit.
    5. Checks whether the known prime factor pair appears.
    6. Records x,y whenever the pair is compatible.

Requirements:

    pip install z3-solver sympy
======================================================================
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import sympy as sp

from z3 import (
    BitVec,
    BitVecVal,
    Solver,
    UGT,
    ULT,
    ZeroExt,
    Or,
    sat,
    unknown,
)


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    3233,                  # 53 * 61
    5959,                  # 59 * 101
    10403,                 # 101 * 103
    18721,                 # 97 * 193
    100127 * 100129,       # 100127 * 100129
]

Z_VALUES = [
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    10,
    12,
    14,
    16,
    20,
    24,
]

# Maximum SAT models examined for one (n,z,a) branch.
MAX_MODELS_PER_BRANCH = 100

# Maximum time spent enumerating one branch.
TIMEOUT_MS = 10_000

# Stop after finding this many prime pairs in a branch.
MAX_PRIME_PAIRS = 20


# ======================================================================
# RESULT
# ======================================================================

@dataclass
class Candidate:
    p: int
    q: int
    x: int | None
    y: int | None
    p_prime: bool
    q_prime: bool

    @property
    def both_prime(self) -> bool:
        return self.p_prime and self.q_prime


# ======================================================================
# BASIC MATH
# ======================================================================

def admissible_a(n: int, z: int) -> list[int]:
    """
    Solve

        n ≡ -a² (mod 2^z)

    for a modulo 2^(z-1).
    """

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    target = n % modulus

    return [
        a
        for a in range(a_modulus)
        if (-a * a) % modulus == target
    ]


def recover_xy(
    p: int,
    q: int,
    z: int,
    a: int,
):
    """
    Recover:

        y = (p+q)/2^z
        x = (q-p+2a)/2^z

    Return (x,y), or None if either is non-integral.
    """

    modulus = 1 << z

    y_num = p + q
    x_num = q - p + 2 * a

    if y_num % modulus != 0:
        return None

    if x_num % modulus != 0:
        return None

    x = x_num // modulus
    y = y_num // modulus

    return x, y


def verify_parameterization(
    p: int,
    q: int,
    z: int,
    a: int,
    x: int,
    y: int,
) -> bool:

    k = 1 << (z - 1)

    p2 = k * (y - x) + a
    q2 = k * (y + x) - a

    return (
        p == p2
        and q == q2
        and p * q == p2 * q2
    )


# ======================================================================
# SAT BUILDER
# ======================================================================

def build_solver(
    n: int,
    z: int,
    a: int,
):

    width = max(2, n.bit_length())

    p = BitVec("p", width)
    q = BitVec("q", width)

    solver = Solver()
    solver.set(timeout=TIMEOUT_MS)

    # --------------------------------------------------------------
    # p*q = n
    # --------------------------------------------------------------

    p_ext = ZeroExt(width, p)
    q_ext = ZeroExt(width, q)

    n_ext = BitVecVal(n, 2 * width)

    solver.add(
        p_ext * q_ext == n_ext
    )

    # --------------------------------------------------------------
    # non-trivial positive factors
    # --------------------------------------------------------------

    solver.add(
        UGT(p, BitVecVal(1, width))
    )

    solver.add(
        UGT(q, BitVecVal(1, width))
    )

    # Remove symmetry.
    #
    # We only enumerate p < q.
    #
    solver.add(
        ULT(p, q)
    )

    # --------------------------------------------------------------
    # FULL 2^z PARAMETERIZATION
    # --------------------------------------------------------------

    modulus = 1 << z
    mask = modulus - 1

    # p+q ≡ 0 mod 2^z

    solver.add(
        ((p + q) & mask)
        == BitVecVal(0, width)
    )

    # q-p ≡ -2a mod 2^z

    target = (-2 * a) % modulus

    solver.add(
        ((q - p) & mask)
        == BitVecVal(target, width)
    )

    return solver, p, q, width


# ======================================================================
# ENUMERATE ONE BRANCH
# ======================================================================

def enumerate_branch(
    n: int,
    z: int,
    a: int,
) -> list[Candidate]:

    solver, p, q, width = build_solver(
        n,
        z,
        a,
    )

    found: list[Candidate] = []

    start = time.perf_counter()

    models_checked = 0
    prime_pairs = 0

    while models_checked < MAX_MODELS_PER_BRANCH:

        elapsed_ms = (
            time.perf_counter() - start
        ) * 1000.0

        if elapsed_ms >= TIMEOUT_MS:
            break

        status = solver.check()

        if status == unknown:
            break

        if status != sat:
            break

        model = solver.model()

        pv = model.eval(p).as_long()
        qv = model.eval(q).as_long()

        xy = recover_xy(
            pv,
            qv,
            z,
            a,
        )

        if xy is None:
            # This should never happen for a correct FULL system.
            x = None
            y = None
        else:
            x, y = xy

        p_prime = sp.isprime(pv)
        q_prime = sp.isprime(qv)

        candidate = Candidate(
            p=pv,
            q=qv,
            x=x,
            y=y,
            p_prime=p_prime,
            q_prime=q_prime,
        )

        found.append(candidate)

        models_checked += 1

        if candidate.both_prime:
            prime_pairs += 1

            if prime_pairs >= MAX_PRIME_PAIRS:
                break

        # ----------------------------------------------------------
        # Block this exact factor pair.
        # ----------------------------------------------------------

        solver.add(
            Or(
                p != BitVecVal(pv, width),
                q != BitVecVal(qv, width),
            )
        )

    return found


# ======================================================================
# ACTUAL FACTORIZATION
# ======================================================================

def factor_n(n: int):

    """
    For benchmark purposes only.

    The experiment uses known factorization to determine whether
    the parameterized SAT search actually discovers the true prime pair.
    """

    factors = sp.factorint(n)

    if len(factors) != 2:
        return None

    expanded = []

    for prime, exponent in factors.items():
        expanded.extend(
            [int(prime)] * exponent
        )

    if len(expanded) != 2:
        return None

    p = min(expanded)
    q = max(expanded)

    if p * q != n:
        return None

    return p, q


# ======================================================================
# PRINT CANDIDATE
# ======================================================================

def print_candidate(
    candidate: Candidate,
    z: int,
    a: int,
    actual_pair,
):

    actual = False

    if actual_pair is not None:

        actual = (
            (candidate.p, candidate.q)
            == actual_pair
        )

    prime_marker = (
        "PRIME×PRIME"
        if candidate.both_prime
        else "composite"
    )

    actual_marker = (
        " <-- ACTUAL"
        if actual
        else ""
    )

    xy = ""

    if (
        candidate.x is not None
        and candidate.y is not None
    ):

        xy = (
            f"x={candidate.x} "
            f"y={candidate.y} "
        )

    print(
        f"      p={candidate.p:<12d} "
        f"q={candidate.q:<12d} "
        f"{prime_marker:<12s} "
        f"{xy}"
        f"{actual_marker}"
    )


# ======================================================================
# ONE N
# ======================================================================

def experiment_n(n: int):

    print()
    print("=" * 110)
    print(f"N = {n}")
    print(f"bits = {n.bit_length()}")
    print("=" * 110)

    actual_pair = factor_n(n)

    if actual_pair is not None:

        print()
        print(
            f"KNOWN PRIME PAIR: "
            f"{actual_pair[0]} × {actual_pair[1]}"
        )

    else:

        print()
        print(
            "KNOWN PRIME PAIR: unavailable"
        )

    # --------------------------------------------------------------
    # z scan
    # --------------------------------------------------------------

    for z in Z_VALUES:

        modulus = 1 << z

        if modulus > (1 << n.bit_length()):
            continue

        a_values = admissible_a(
            n,
            z,
        )

        print()
        print("-" * 110)
        print(
            f"z={z:2d} "
            f"2^z={modulus:<10d} "
            f"n mod 2^z={n % modulus:<10d} "
            f"a={a_values}"
        )

        if not a_values:

            print(
                "  no admissible a"
            )

            continue

        found_any = False
        found_actual = False

        for a in a_values:

            print()
            print(
                f"  BRANCH a={a}"
            )

            start = time.perf_counter()

            candidates = enumerate_branch(
                n,
                z,
                a,
            )

            elapsed = (
                time.perf_counter()
                - start
            )

            if not candidates:

                print(
                    f"    no models "
                    f"({elapsed:.6f}s)"
                )

                continue

            found_any = True

            print(
                f"    models={len(candidates)} "
                f"time={elapsed:.6f}s"
            )

            for candidate in candidates:

                print_candidate(
                    candidate,
                    z,
                    a,
                    actual_pair,
                )

                if actual_pair is not None:

                    if (
                        candidate.p,
                        candidate.q,
                    ) == actual_pair:

                        found_actual = True

            prime_count = sum(
                c.both_prime
                for c in candidates
            )

            print(
                f"    prime pairs in branch: "
                f"{prime_count}"
            )

        if found_actual:

            print()
            print(
                f"  >>> ACTUAL PRIME PAIR FOUND "
                f"AT z={z}"
            )

            # We intentionally continue scanning z.
            # This lets us observe how the actual pair survives
            # or disappears as z increases.


# ======================================================================
# LARGE EXAMPLE: SHOW ACTUAL a/x/y PROGRESSION
# ======================================================================

def analyze_actual_progression(
    n: int,
    p: int,
    q: int,
):

    print()
    print("=" * 110)
    print("ACTUAL FACTOR 2^z PROGRESSION")
    print("=" * 110)

    print(
        f"n={n}"
    )

    print(
        f"p={p}"
    )

    print(
        f"q={q}"
    )

    print()

    print(
        f"{'z':>3} "
        f"{'2^z':>10} "
        f"{'n mod 2^z':>12} "
        f"{'a':>12} "
        f"{'x':>12} "
        f"{'y':>12}"
    )

    print("-" * 70)

    for z in range(2, 25):

        modulus = 1 << z

        a_values = admissible_a(
            n,
            z,
        )

        matches = []

        for a in a_values:

            xy = recover_xy(
                p,
                q,
                z,
                a,
            )

            if xy is not None:

                x, y = xy

                matches.append(
                    (a, x, y)
                )

        if matches:

            text = "; ".join(
                f"a={a},x={x},y={y}"
                for a, x, y in matches
            )

        else:

            text = "-"

        print(
            f"{z:3d} "
            f"{modulus:10d} "
            f"{n % modulus:12d} "
            f"{text}"
        )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "2^z PARAMETERIZATION / PRIME ENUMERATION EXPERIMENT"
    )
    print("=" * 110)

    # --------------------------------------------------------------
    # Show the exact progression for your large example.
    # --------------------------------------------------------------

    large_n = 100127 * 100129

    analyze_actual_progression(
        large_n,
        100127,
        100129,
    )

    # --------------------------------------------------------------
    # Enumerate branches.
    # --------------------------------------------------------------

    for n in N_VALUES:

        experiment_n(n)


if __name__ == "__main__":
    main()
