#!/usr/bin/env python3

"""
======================================================================
SAFE 2^z BRANCH SAT EXPERIMENT
======================================================================

For

    k = 2^(z-1)

    p = k(y-x) + a
    q = k(y+x) - a

we have

    pq ≡ -a^2 (mod 2^z)

and the FULL parameterization requires

    p + q ≡ 0       (mod 2^z)
    q - p ≡ -2a     (mod 2^z)

This experiment deliberately performs ONLY ONE SAT QUERY
per (z,a) branch.

It does NOT enumerate all models.

It does NOT call factorint().

It is designed to avoid freezing the VM.

Requirements:

    pip install z3-solver sympy
"""

from __future__ import annotations

import time
import sympy as sp

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
    5959,                  # 59 * 101
    3233,                  # 53 * 61
    10403,                 # 101 * 103
    18721,                 # 97 * 193

    # IMPORTANT:
    # 10025616383 = 223 * 449 * 100129
    10025616383,
]

Z_VALUES = [
    2, 3, 4, 5, 6, 7, 8,
    10, 12, 14, 16, 20, 24
]

# Much shorter timeout to protect the machine.
TIMEOUT_MS = 1_000


# ======================================================================
# KNOWN FACTORS FOR BENCHMARKING
# ======================================================================

KNOWN_FACTORS = {
    5959: (59, 101),
    3233: (53, 61),
    10403: (101, 103),
    18721: (97, 193),

    # Full factorization:
    10025616383: (223, 449, 100129),
}


# ======================================================================
# HELPERS
# ======================================================================

def admissible_a(n: int, z: int) -> list[int]:
    """
    Solve

        n ≡ -a^2 (mod 2^z)

    for a modulo 2^(z-1).

    For the benchmark sizes this brute-force loop is tiny.
    """

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    # Avoid absurd searches for huge z.
    # For the current experiment we only use moderate z.
    if a_modulus > 1_000_000:
        return []

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
    y = (p+q)/2^z

    x = (q-p+2a)/2^z
    """

    modulus = 1 << z

    y_num = p + q
    x_num = q - p + 2 * a

    if y_num % modulus != 0:
        return None

    if x_num % modulus != 0:
        return None

    return (
        x_num // modulus,
        y_num // modulus,
    )


def is_known_pair(
    n: int,
    p: int,
    q: int,
) -> bool:

    known = KNOWN_FACTORS.get(n)

    if known is None:
        return False

    # Any pair of the known factors whose product is n.
    if p * q != n:
        return False

    return True


def factor_type(
    p: int,
    q: int,
) -> str:

    pp = sp.isprime(p)
    qq = sp.isprime(q)

    if pp and qq:
        return "PRIME x PRIME"

    if pp:
        return "PRIME x COMPOSITE"

    if qq:
        return "COMPOSITE x PRIME"

    return "COMPOSITE x COMPOSITE"


# ======================================================================
# SAT SOLVER
# ======================================================================

def solve_branch(
    n: int,
    z: int,
    a: int,
):

    width = max(2, n.bit_length())

    p = BitVec(
        f"p_{z}_{a}",
        width,
    )

    q = BitVec(
        f"q_{z}_{a}",
        width,
    )

    solver = Solver()

    solver.set(
        timeout=TIMEOUT_MS
    )

    # --------------------------------------------------------------
    # p*q=n
    # --------------------------------------------------------------

    p_ext = ZeroExt(width, p)
    q_ext = ZeroExt(width, q)

    n_ext = BitVecVal(
        n,
        2 * width,
    )

    solver.add(
        p_ext * q_ext == n_ext
    )

    # Positive nontrivial factors.

    solver.add(
        UGT(
            p,
            BitVecVal(1, width),
        )
    )

    solver.add(
        UGT(
            q,
            BitVecVal(1, width),
        )
    )

    # Remove symmetry.

    solver.add(
        ULT(p, q)
    )

    # --------------------------------------------------------------
    # FULL PARAMETERIZATION
    # --------------------------------------------------------------

    modulus = 1 << z
    mask = modulus - 1

    # p + q ≡ 0 mod 2^z

    solver.add(
        ((p + q) & mask)
        == BitVecVal(
            0,
            width,
        )
    )

    # q - p ≡ -2a mod 2^z

    target = (-2 * a) % modulus

    solver.add(
        ((q - p) & mask)
        == BitVecVal(
            target,
            width,
        )
    )

    # --------------------------------------------------------------
    # SOLVE ONCE
    # --------------------------------------------------------------

    start = time.perf_counter()

    status = solver.check()

    elapsed = time.perf_counter() - start

    if status == unknown:

        return {
            "status": "TIMEOUT",
            "time": elapsed,
        }

    if status != sat:

        return {
            "status": "UNSAT",
            "time": elapsed,
        }

    model = solver.model()

    p_value = model.eval(p).as_long()
    q_value = model.eval(q).as_long()

    xy = recover_xy(
        p_value,
        q_value,
        z,
        a,
    )

    if xy is None:

        return {
            "status": "BUG",
            "time": elapsed,
            "p": p_value,
            "q": q_value,
        }

    x, y = xy

    return {
        "status": "SAT",
        "time": elapsed,
        "p": p_value,
        "q": q_value,
        "x": x,
        "y": y,
        "type": factor_type(
            p_value,
            q_value,
        ),
        "known_pair": is_known_pair(
            n,
            p_value,
            q_value,
        ),
    }


# ======================================================================
# ONE N
# ======================================================================

def experiment_n(n: int):

    print()
    print("=" * 100)
    print(f"N = {n}")
    print(f"bits = {n.bit_length()}")
    print("=" * 100)

    known = KNOWN_FACTORS.get(n)

    if known is not None:

        print(
            "KNOWN FACTORS =",
            " × ".join(
                str(x)
                for x in known
            )
        )

    print()
    print("z / residue / admissible a")
    print("-" * 100)

    for z in Z_VALUES:

        modulus = 1 << z

        if modulus > (1 << n.bit_length()):
            continue

        a_values = admissible_a(
            n,
            z,
        )

        print(
            f"z={z:2d} "
            f"n mod 2^z={n % modulus:<10d} "
            f"a={a_values}"
        )

    print()
    print("ONE-SHOT SAT BRANCH TEST")
    print("-" * 100)

    total_sat = 0
    total_prime_prime = 0

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
            f"2^z={modulus}"
        )

        for a in a_values:

            result = solve_branch(
                n,
                z,
                a,
            )

            status = result["status"]

            print(
                f"  a={a:<8d} "
                f"{status:<8s} "
                f"{result['time']:.4f}s",
                end="",
            )

            if status == "SAT":

                total_sat += 1

                p = result["p"]
                q = result["q"]

                print(
                    f" "
                    f"p={p:<12d} "
                    f"q={q:<12d} "
                    f"{result['type']:<20s} "
                    f"x={result['x']:<8d} "
                    f"y={result['y']:<8d}",
                    end="",
                )

                if result["known_pair"]:

                    print(
                        " <-- KNOWN PAIR"
                    )

                else:

                    print()

                if result["type"] == "PRIME x PRIME":

                    total_prime_prime += 1

            else:

                print()

    print()
    print("-" * 100)
    print(
        f"SAT branches        = {total_sat}"
    )
    print(
        f"prime × prime hits  = {total_prime_prime}"
    )


# ======================================================================
# SPECIAL LARGE-N ANALYSIS
# ======================================================================

def large_example():

    n = 10025616383

    print()
    print("=" * 100)
    print("LARGE EXAMPLE FACTOR ANALYSIS")
    print("=" * 100)

    print()
    print(
        "Known:"
    )

    print(
        "  n = 223 × 449 × 100129"
    )

    print()
    print(
        "Your previous pair:"
    )

    print(
        "  100127 × 100129"
    )

    print()
    print(
        "but:"
    )

    print(
        "  100127 = 223 × 449"
    )

    print()

    for p in [
        223,
        449,
        100127,
        100129,
    ]:

        print(
            f"{p:8d} "
            f"prime={sp.isprime(p)}"
        )

    # --------------------------------------------------------------
    # Test known pairs.
    # --------------------------------------------------------------

    pairs = [
        (223, 449),
        (100127, 100129),
    ]

    print()
    print(
        "PARAMETERIZATION OF KNOWN FACTOR PAIRS"
    )
    print("-" * 100)

    for p, q in pairs:

        if p * q == n:

            print(
                f"{p} × {q} = n"
            )

        else:

            # Only the second pair is actually the cofactor pair.
            print(
                f"{p} × {q} = {p*q}"
            )

        for z in range(2, 13):

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
                        (
                            a,
                            x,
                            y,
                        )
                    )

            if matches:

                print(
                    f"  z={z:2d} "
                    f"a/x/y={matches}"
                )

        print()


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 100)
    print(
        "SAFE 2^z ONE-SHOT SAT EXPERIMENT"
    )
    print("=" * 100)

    large_example()

    for n in N_VALUES:

        experiment_n(n)


if __name__ == "__main__":
    main()
