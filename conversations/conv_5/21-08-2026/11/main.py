#!/usr/bin/env python3

"""
======================================================================
2^z PARAMETERIZATION SAT EXPERIMENT -- CORRECTED
======================================================================

GENERAL FORM

    k = 2^(z-1)

    p = k(y-x) + a
    q = k(y+x) - a

EXPANSION

    pq = k^2(y^2-x^2) + 2akx - a^2

MODULO 2^z

    pq ≡ -a^2 (mod 2^z)

RECONSTRUCTION

    p + q = 2^z y

    q - p = 2^z x - 2a

therefore

    y = (p+q) / 2^z

    x = (q-p+2a) / 2^z

and the SAT constraints are

    p + q ≡ 0       (mod 2^z)

    q - p ≡ -2a     (mod 2^z)


EXPERIMENTS

    BASELINE
        p*q = n

    RESIDUE
        p ≡  a  (mod 2^(z-1))
        q ≡ -a  (mod 2^(z-1))

    FULL
        p*q = n
        p+q ≡ 0   (mod 2^z)
        q-p ≡ -2a (mod 2^z)

The FULL system is the actual x,y parameterization.

======================================================================
"""

from __future__ import annotations

import time
from dataclasses import dataclass

from z3 import (
    BitVec,
    BitVecVal,
    Solver,
    UGT,
    ZeroExt,
    sat,
    unknown,
)


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    3233,                    # 53 * 61
    5959,                    # 59 * 101
    10403,                   # 101 * 103
    18721,                   # 97 * 193
    100127 * 100129,
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
]

TIMEOUT_MS = 10_000


# ======================================================================
# RESULT
# ======================================================================

@dataclass
class Result:
    mode: str
    n: int
    z: int | None
    a: int | None
    status: str
    elapsed: float
    p: int | None
    q: int | None
    x: int | None
    y: int | None
    correct: bool


# ======================================================================
# HELPERS
# ======================================================================

def bit_width(n: int) -> int:
    return max(2, n.bit_length())


def admissible_a_values(n: int, z: int) -> list[int]:
    """
    Find all a modulo 2^(z-1) satisfying

        n ≡ -a² (mod 2^z)
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
    Recover

        y = (p+q)/2^z

        x = (q-p+2a)/2^z

    Returns (x,y) or (None,None).
    """

    modulus = 1 << z

    y_num = p + q
    x_num = q - p + 2 * a

    if y_num % modulus != 0:
        return None, None

    if x_num % modulus != 0:
        return None, None

    y = y_num // modulus
    x = x_num // modulus

    return x, y


# ======================================================================
# BASE SAT
# ======================================================================

def build_base_solver(n: int):

    width = bit_width(n)

    p = BitVec("p", width)
    q = BitVec("q", width)

    solver = Solver()
    solver.set(timeout=TIMEOUT_MS)

    p_ext = ZeroExt(width, p)
    q_ext = ZeroExt(width, q)

    n_ext = BitVecVal(n, 2 * width)

    solver.add(
        p_ext * q_ext == n_ext
    )

    solver.add(
        UGT(
            p,
            BitVecVal(1, width)
        )
    )

    solver.add(
        UGT(
            q,
            BitVecVal(1, width)
        )
    )

    return solver, p, q, width


# ======================================================================
# BASELINE
# ======================================================================

def solve_baseline(n: int) -> Result:

    solver, p, q, _ = build_base_solver(n)

    start = time.perf_counter()
    status = solver.check()
    elapsed = time.perf_counter() - start

    if status == sat:

        model = solver.model()

        pv = model.eval(p).as_long()
        qv = model.eval(q).as_long()

        correct = (
            pv > 1
            and qv > 1
            and pv * qv == n
        )

        return Result(
            mode="BASELINE",
            n=n,
            z=None,
            a=None,
            status="SAT",
            elapsed=elapsed,
            p=pv,
            q=qv,
            x=None,
            y=None,
            correct=correct,
        )

    return Result(
        mode="BASELINE",
        n=n,
        z=None,
        a=None,
        status=(
            "UNKNOWN"
            if status == unknown
            else "UNSAT"
        ),
        elapsed=elapsed,
        p=None,
        q=None,
        x=None,
        y=None,
        correct=False,
    )


# ======================================================================
# RESIDUE ONLY
# ======================================================================

def solve_residue(
    n: int,
    z: int,
    a: int,
) -> Result:

    solver, p, q, width = build_base_solver(n)

    low_bits = z - 1
    low_modulus = 1 << low_bits
    low_mask = low_modulus - 1

    p_target = a % low_modulus
    q_target = (-a) % low_modulus

    solver.add(
        (p & low_mask)
        == BitVecVal(
            p_target,
            width,
        )
    )

    solver.add(
        (q & low_mask)
        == BitVecVal(
            q_target,
            width,
        )
    )

    start = time.perf_counter()
    status = solver.check()
    elapsed = time.perf_counter() - start

    if status != sat:

        return Result(
            mode="RESIDUE",
            n=n,
            z=z,
            a=a,
            status=(
                "UNKNOWN"
                if status == unknown
                else "UNSAT"
            ),
            elapsed=elapsed,
            p=None,
            q=None,
            x=None,
            y=None,
            correct=False,
        )

    model = solver.model()

    pv = model.eval(p).as_long()
    qv = model.eval(q).as_long()

    xv, yv = recover_xy(
        pv,
        qv,
        z,
        a,
    )

    return Result(
        mode="RESIDUE",
        n=n,
        z=z,
        a=a,
        status="SAT",
        elapsed=elapsed,
        p=pv,
        q=qv,
        x=xv,
        y=yv,
        correct=(
            pv * qv == n
            and pv > 1
            and qv > 1
        ),
    )


# ======================================================================
# FULL PARAMETERIZATION
# ======================================================================

def solve_full(
    n: int,
    z: int,
    a: int,
) -> Result:

    solver, p, q, width = build_base_solver(n)

    modulus = 1 << z
    mask = modulus - 1

    # --------------------------------------------------------------
    # p + q ≡ 0 (mod 2^z)
    # --------------------------------------------------------------

    solver.add(
        ((p + q) & mask)
        == BitVecVal(
            0,
            width,
        )
    )

    # --------------------------------------------------------------
    # q - p ≡ -2a (mod 2^z)
    #
    # THIS IS THE CORRECTED SIGN.
    # --------------------------------------------------------------

    target = (-2 * a) % modulus

    solver.add(
        ((q - p) & mask)
        == BitVecVal(
            target,
            width,
        )
    )

    start = time.perf_counter()
    status = solver.check()
    elapsed = time.perf_counter() - start

    if status != sat:

        return Result(
            mode="FULL",
            n=n,
            z=z,
            a=a,
            status=(
                "UNKNOWN"
                if status == unknown
                else "UNSAT"
            ),
            elapsed=elapsed,
            p=None,
            q=None,
            x=None,
            y=None,
            correct=False,
        )

    model = solver.model()

    pv = model.eval(p).as_long()
    qv = model.eval(q).as_long()

    xv, yv = recover_xy(
        pv,
        qv,
        z,
        a,
    )

    correct = (
        pv * qv == n
        and pv > 1
        and qv > 1
        and xv is not None
        and yv is not None
    )

    return Result(
        mode="FULL",
        n=n,
        z=z,
        a=a,
        status="SAT",
        elapsed=elapsed,
        p=pv,
        q=qv,
        x=xv,
        y=yv,
        correct=correct,
    )


# ======================================================================
# PRINT RESULT
# ======================================================================

def print_result(r: Result):

    z = "-" if r.z is None else str(r.z)
    a = "-" if r.a is None else str(r.a)

    if r.p is None:

        print(
            f"{r.mode:8s} "
            f"z={z:>3s} "
            f"a={a:>5s} "
            f"{r.status:12s} "
            f"time={r.elapsed:.6f}s"
        )

        return

    xy = ""

    if r.x is not None and r.y is not None:

        xy = (
            f" x={r.x:<8d}"
            f" y={r.y:<8d}"
        )

    print(
        f"{r.mode:8s} "
        f"z={z:>3s} "
        f"a={a:>5s} "
        f"{r.status:12s} "
        f"p={r.p:<10d} "
        f"q={r.q:<10d}"
        f"{xy} "
        f"time={r.elapsed:.6f}s "
        f"correct={r.correct}"
    )


# ======================================================================
# VERIFY PARAMETERIZATION
# ======================================================================

def verify_parameterization():

    print()
    print("=" * 100)
    print("PARAMETERIZATION VERIFICATION")
    print("=" * 100)

    identity_failures = 0
    modular_failures = 0
    reconstruction_failures = 0

    for z in range(2, 13):

        k = 1 << (z - 1)
        modulus = 1 << z

        for a in range(0, min(modulus, 64)):

            for x in range(-5, 6):
                for y in range(-5, 6):

                    p = k * (y - x) + a
                    q = k * (y + x) - a

                    n = p * q

                    # --------------------------------------------------
                    # Polynomial identity
                    # --------------------------------------------------

                    polynomial = (
                        k * k * (y * y - x * x)
                        + 2 * a * k * x
                        - a * a
                    )

                    if n != polynomial:
                        identity_failures += 1

                    # --------------------------------------------------
                    # Modular identity
                    # --------------------------------------------------

                    if n % modulus != (-a * a) % modulus:
                        modular_failures += 1

                    # --------------------------------------------------
                    # Reconstruction
                    # --------------------------------------------------

                    rx, ry = recover_xy(
                        p,
                        q,
                        z,
                        a,
                    )

                    if rx != x or ry != y:
                        reconstruction_failures += 1

    print()
    print(
        f"Identity failures       = {identity_failures}"
    )
    print(
        f"Modular failures        = {modular_failures}"
    )
    print(
        f"Reconstruction failures = {reconstruction_failures}"
    )

    print()

    ok = (
        identity_failures == 0
        and modular_failures == 0
        and reconstruction_failures == 0
    )

    print(f"PASS = {ok}")


# ======================================================================
# ANALYZE KNOWN FACTORS
# ======================================================================

def analyze_known_factors(
    n: int,
    p: int,
    q: int,
):

    print()
    print("=" * 100)
    print("KNOWN FACTOR PARAMETER ANALYSIS")
    print("=" * 100)

    print(f"n = {n}")
    print(f"p = {p}")
    print(f"q = {q}")

    for z in Z_VALUES:

        modulus = 1 << z

        if modulus > n * 2:
            break

        candidates = admissible_a_values(
            n,
            z,
        )

        print()
        print(
            f"z={z:2d} "
            f"2^z={modulus:<8d} "
            f"n mod 2^z={n % modulus:<8d}"
        )

        found = False

        for a in candidates:

            x, y = recover_xy(
                p,
                q,
                z,
                a,
            )

            if x is None:
                continue

            k = 1 << (z - 1)

            print(
                f"  a={a:<6d}"
                f" -> x={x:<8d}"
                f" y={y:<8d}"
                f" | "
                f"p={k}(y-x)+a"
                f" | "
                f"q={k}(y+x)-a"
            )

            found = True

        if not found:

            print(
                "  No integer x,y for this "
                "ordered factor pair."
            )


# ======================================================================
# SAT EXPERIMENT
# ======================================================================

def experiment_n(n: int):

    print()
    print("#" * 100)
    print(f"N = {n}")
    print("#" * 100)

    # --------------------------------------------------------------
    # Baseline only once.
    # --------------------------------------------------------------

    print()
    print("BASELINE")
    print("-" * 100)

    baseline = solve_baseline(n)

    print_result(baseline)

    if baseline.p is None:
        print("No factorization recovered.")
        return

    p = baseline.p
    q = baseline.q

    print()
    print(
        f"Recovered: "
        f"{min(p,q)} × {max(p,q)} = {n}"
    )

    analyze_known_factors(
        n,
        p,
        q,
    )

    # --------------------------------------------------------------
    # SAT comparison.
    # --------------------------------------------------------------

    print()
    print("=" * 100)
    print("SAT COMPARISON")
    print("=" * 100)

    for z in Z_VALUES:

        modulus = 1 << z

        if modulus > n.bit_length() * 2:
            continue

        candidates = admissible_a_values(
            n,
            z,
        )

        print()
        print(
            f"z={z:2d} "
            f"2^z={modulus:<8d} "
            f"n mod 2^z={n % modulus:<8d}"
        )

        if not candidates:

            print("  NO ADMISSIBLE a")
            continue

        for a in candidates:

            print()
            print(f"  a = {a}")

            residue = solve_residue(
                n,
                z,
                a,
            )

            full = solve_full(
                n,
                z,
                a,
            )

            print("    ", end="")
            print_result(residue)

            print("    ", end="")
            print_result(full)


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 100)
    print("CORRECTED 2^z PARAMETERIZATION SAT EXPERIMENT")
    print("=" * 100)

    print()
    print(
        "p = 2^(z-1)(y-x) + a"
    )
    print(
        "q = 2^(z-1)(y+x) - a"
    )

    print()
    print(
        "Therefore:"
    )
    print(
        "p + q ≡ 0       (mod 2^z)"
    )
    print(
        "q - p ≡ -2a     (mod 2^z)"
    )

    verify_parameterization()

    for n in N_VALUES:
        experiment_n(n)


if __name__ == "__main__":
    main()