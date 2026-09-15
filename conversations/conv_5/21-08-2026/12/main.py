#!/usr/bin/env python3

"""
======================================================================
FOCUSED 2^z STRUCTURED-SAT FACTORIZATION EXPERIMENT
======================================================================

Tests:

    p = 2^(z-1)(y-x) + a
    q = 2^(z-1)(y+x) - a

Therefore:

    pq ≡ -a^2 (mod 2^z)

Full parameterization requires:

    p + q ≡ 0       (mod 2^z)
    q - p ≡ -2a     (mod 2^z)

The experiment runs STRUCTURED SAT FIRST.

This is important for large n: the unconstrained p*q=n
problem can consume the entire timeout before the useful
2-adic constraints are tested.

Requires:

    pip install z3-solver
======================================================================
"""

from __future__ import annotations

import time

from z3 import (
    BitVec,
    BitVecVal,
    Solver,
    UGT,
    ULT,
    ZeroExt,
    sat,
    unknown,
)


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    5959,
    3233,
    10403,
    18721,

    # Your important large test:
    100127 * 100129,
]

Z_VALUES = [
    2, 3, 4, 5, 6, 7, 8,
    10, 12, 14, 16, 20, 24
]

# Timeout for each SAT query.
STRUCTURED_TIMEOUT_MS = 15_000
BASELINE_TIMEOUT_MS = 5_000

# Set False if you do not want the expensive baseline.
RUN_BASELINE = True


# ======================================================================
# HELPERS
# ======================================================================

def width_for(n: int) -> int:
    return max(2, n.bit_length())


def admissible_a(n: int, z: int) -> list[int]:
    """
    Find a in [0, 2^(z-1)-1] satisfying

        n ≡ -a^2 (mod 2^z)
    """

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    target = n % modulus

    return [
        a
        for a in range(a_modulus)
        if (-a * a) % modulus == target
    ]


def recover_xy(p: int, q: int, z: int, a: int):
    """
    From

        p = 2^(z-1)(y-x) + a
        q = 2^(z-1)(y+x) - a

    recover

        y = (p+q)/2^z
        x = (q-p+2a)/2^z
    """

    modulus = 1 << z

    y_num = p + q
    x_num = q - p + 2 * a

    if y_num % modulus:
        return None

    if x_num % modulus:
        return None

    return (
        x_num // modulus,
        y_num // modulus,
    )


# ======================================================================
# BUILD BASE SOLVER
# ======================================================================

def build_solver(n: int):
    width = width_for(n)

    p = BitVec("p", width)
    q = BitVec("q", width)

    solver = Solver()

    p_ext = ZeroExt(width, p)
    q_ext = ZeroExt(width, q)

    n_ext = BitVecVal(n, 2 * width)

    solver.add(
        p_ext * q_ext == n_ext
    )

    solver.add(
        UGT(p, BitVecVal(1, width))
    )

    solver.add(
        UGT(q, BitVecVal(1, width))
    )

    # Remove the p/q symmetry.
    solver.add(
        ULT(p, q)
    )

    return solver, p, q, width


# ======================================================================
# FULL PARAMETERIZED SAT
# ======================================================================

def solve_full(n: int, z: int, a: int):

    solver, p, q, width = build_solver(n)

    modulus = 1 << z
    mask = modulus - 1

    # --------------------------------------------------------------
    # p + q ≡ 0 (mod 2^z)
    # --------------------------------------------------------------

    solver.add(
        ((p + q) & mask)
        == BitVecVal(0, width)
    )

    # --------------------------------------------------------------
    # q - p ≡ -2a (mod 2^z)
    # --------------------------------------------------------------

    target = (-2 * a) % modulus

    solver.add(
        ((q - p) & mask)
        == BitVecVal(target, width)
    )

    solver.set(
        timeout=STRUCTURED_TIMEOUT_MS
    )

    start = time.perf_counter()

    status = solver.check()

    elapsed = time.perf_counter() - start

    if status == sat:

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
            return {
                "status": "BUG",
                "time": elapsed,
                "p": pv,
                "q": qv,
                "x": None,
                "y": None,
            }

        x, y = xy

        valid = (
            pv * qv == n
            and pv < qv
            and pv == (1 << (z - 1)) * (y - x) + a
            and qv == (1 << (z - 1)) * (y + x) - a
        )

        return {
            "status": "SAT",
            "time": elapsed,
            "p": pv,
            "q": qv,
            "x": x,
            "y": y,
            "valid": valid,
        }

    if status == unknown:

        return {
            "status": "TIMEOUT",
            "time": elapsed,
            "p": None,
            "q": None,
            "x": None,
            "y": None,
        }

    return {
        "status": "UNSAT",
        "time": elapsed,
        "p": None,
        "q": None,
        "x": None,
        "y": None,
    }


# ======================================================================
# BASELINE
# ======================================================================

def solve_baseline(n: int):

    solver, p, q, width = build_solver(n)

    solver.set(
        timeout=BASELINE_TIMEOUT_MS
    )

    start = time.perf_counter()

    status = solver.check()

    elapsed = time.perf_counter() - start

    if status == sat:

        model = solver.model()

        pv = model.eval(p).as_long()
        qv = model.eval(q).as_long()

        return (
            "SAT",
            elapsed,
            pv,
            qv,
        )

    if status == unknown:

        return (
            "TIMEOUT",
            elapsed,
            None,
            None,
        )

    return (
        "UNSAT",
        elapsed,
        None,
        None,
    )


# ======================================================================
# EXPERIMENT ONE N
# ======================================================================

def experiment(n: int):

    print()
    print("=" * 110)
    print(f"N = {n}")
    print(f"bits = {n.bit_length()}")
    print("=" * 110)

    print()
    print("LOW-BIT STRUCTURE")
    print("-" * 110)

    previous_target = None

    for z in Z_VALUES:

        modulus = 1 << z

        if modulus > (1 << n.bit_length()):
            break

        residue = n % modulus
        a_values = admissible_a(n, z)

        changed = (
            previous_target is None
            or residue != previous_target
        )

        if changed:

            print()
            print(
                f"z={z:2d} "
                f"2^z={modulus:<10d} "
                f"n mod 2^z={residue:<10d} "
                f"a={a_values}"
            )

        previous_target = residue

    # --------------------------------------------------------------
    # STRUCTURED SAT
    # --------------------------------------------------------------

    print()
    print("STRUCTURED SAT")
    print("-" * 110)

    solved = False

    for z in Z_VALUES:

        modulus = 1 << z

        a_values = admissible_a(
            n,
            z,
        )

        if not a_values:
            continue

        print()
        print(
            f"z={z:2d} "
            f"2^z={modulus:<8d} "
            f"n mod 2^z={n % modulus:<8d}"
        )

        for a in a_values:

            result = solve_full(
                n,
                z,
                a,
            )

            print(
                f"  a={a:<8d} "
                f"{result['status']:<8s} "
                f"time={result['time']:.6f}s",
                end="",
            )

            if result["status"] == "SAT":

                print(
                    f" "
                    f"p={result['p']} "
                    f"q={result['q']} "
                    f"x={result['x']} "
                    f"y={result['y']} "
                    f"valid={result['valid']}"
                )

                if result["valid"]:
                    solved = True
                    break

            else:

                print()

        if solved:
            break

    if not solved:
        print()
        print(
            "No structured factorization found "
            "within the tested z/a branches."
        )

    # --------------------------------------------------------------
    # BASELINE AFTER STRUCTURED
    # --------------------------------------------------------------

    if RUN_BASELINE:

        print()
        print("BASELINE SAT")
        print("-" * 110)

        status, elapsed, p, q = solve_baseline(n)

        if status == "SAT":

            print(
                f"BASELINE SAT "
                f"time={elapsed:.6f}s "
                f"p={p} "
                f"q={q}"
            )

        else:

            print(
                f"BASELINE {status} "
                f"time={elapsed:.6f}s"
            )


# ======================================================================
# EXACT KNOWN EXAMPLE
# ======================================================================

def verify_large_example():

    n = 100127 * 100129

    p = 100127
    q = 100129

    print()
    print("=" * 110)
    print("KNOWN LARGE EXAMPLE")
    print("=" * 110)

    print(f"n = {n}")
    print(f"p = {p}")
    print(f"q = {q}")

    for z in range(2, 13):

        modulus = 1 << z
        k = 1 << (z - 1)

        candidates = admissible_a(
            n,
            z,
        )

        matches = []

        for a in candidates:

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

        print(
            f"z={z:2d} "
            f"n%2^z={n % modulus:<5d} "
            f"a={candidates} "
            f"ACTUAL={matches}"
        )

        if z >= 8 and not candidates:
            pass


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print("2^z STRUCTURED FACTORIZATION SAT BENCHMARK")
    print("=" * 110)

    verify_large_example()

    for n in N_VALUES:
        experiment(n)


if __name__ == "__main__":
    main()
