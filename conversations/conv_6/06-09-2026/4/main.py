#!/usr/bin/env python3

# ============================================================
# START EXPERIMENT 162
# Product Collapse of the Cross-Radix Carry Correction
# ============================================================

from __future__ import annotations

import math
import random
from collections import defaultdict

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BIT_SIZES = [30, 36, 42, 48, 54]
SAMPLES_PER_SIZE = 2

SEED = 162

# ------------------------------------------------------------
# Basic number theory
# ------------------------------------------------------------

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37,
    )

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

    # Deterministic for the ranges used here.
    bases = [2, 3, 5, 7, 11, 13, 17]

    for a in bases:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
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
        x = rng.getrandbits(bits)
        x |= (1 << (bits - 1))
        x |= 1

        if is_probable_prime(x):
            return x


def make_semiprime(bits: int, rng: random.Random) -> tuple[int, int, int]:
    p_bits = bits // 2
    q_bits = bits - p_bits

    while True:
        p = random_prime(p_bits, rng)
        q = random_prime(q_bits, rng)

        if p != q:
            return p, q, p * q


# ------------------------------------------------------------
# CRT
# ------------------------------------------------------------

def crt_pair(a: int, m: int, b: int, n: int) -> int:
    """
    x == a (mod m)
    x == b (mod n)

    m,n coprime.
    """
    t = ((b - a) * pow(m, -1, n)) % n
    return a + m * t


# ------------------------------------------------------------
# Exact carry decomposition
# ------------------------------------------------------------

def exact_C_from_pq(
    p: int,
    q: int,
    r: int,
    s: int,
) -> tuple[int, int, int, int, int, int]:
    """
    Exact cell decomposition.

    Returns:
        a, b, alpha, beta, C, numerator
    """

    a = p % r
    b = q % s

    k = p // r
    ell = q // s

    alpha = (ell * a) % r
    beta = (k * b) % s

    numerator = r * beta + s * alpha + a * b
    C = numerator // (r * s)

    return a, b, alpha, beta, C, numerator


# ------------------------------------------------------------
# Direct residue-only reconstruction
# ------------------------------------------------------------

def exact_C_from_ab(
    n: int,
    a: int,
    b: int,
    r: int,
    s: int,
) -> tuple[int, int, int]:
    """
    Compute C from n,r,s and the local factor residues

        a = p mod r
        b = q mod s

    without knowing p or q.

    From

        pq = n

    we have

        p == n*b^{-1} (mod s)
        q == n*a^{-1} (mod r)

    and therefore recover k mod s and ell mod r.
    """

    # p mod s
    p_mod_s = (n * pow(b, -1, s)) % s

    # q mod r
    q_mod_r = (n * pow(a, -1, r)) % r

    # Reconstruct p mod rs.
    p_mod_rs = crt_pair(a, r, p_mod_s, s)

    # Reconstruct q mod rs.
    q_mod_rs = crt_pair(q_mod_r, r, b, s)

    # p = r*k + a
    k_mod_s = ((p_mod_rs - a) // r) % s

    # q = s*ell + b
    ell_mod_r = ((q_mod_rs - b) // s) % r

    alpha = (ell_mod_r * a) % r
    beta = (k_mod_s * b) % s

    numerator = r * beta + s * alpha + a * b
    C = numerator // (r * s)

    return alpha, beta, C


# ------------------------------------------------------------
# NEW simplified formula
# ------------------------------------------------------------

def simplified_C(
    n: int,
    a: int,
    b: int,
    r: int,
    s: int,
) -> tuple[int, int, int]:
    """
    New algebraic formula.

    Let

        D = n - ab.

    Since

        n-ab == r*beta (mod s),

    we have

        beta == (n-ab) * r^{-1} (mod s).

    Similarly,

        alpha == (n-ab) * s^{-1} (mod r).

    Thus:

        C =
        floor(
            [r*beta + s*alpha + ab] / (rs)
        )

    with alpha,beta obtained directly from n-ab.
    """

    z = a * b
    D = n - z

    beta = (D * pow(r, -1, s)) % s
    alpha = (D * pow(s, -1, r)) % r

    numerator = r * beta + s * alpha + z
    C = numerator // (r * s)

    return alpha, beta, C


# ------------------------------------------------------------
# Product-only function
# ------------------------------------------------------------

def product_only_C(
    n: int,
    z: int,
    r: int,
    s: int,
) -> int:
    """
    C as a function of z = ab alone.
    """

    D = n - z

    beta = (D * pow(r, -1, s)) % s
    alpha = (D * pow(s, -1, r)) % r

    return (r * beta + s * alpha + z) // (r * s)


# ------------------------------------------------------------
# Analyze one cell
# ------------------------------------------------------------

def analyze_cell(
    n: int,
    p: int,
    q: int,
    r: int,
    s: int,
) -> dict:

    raw_pairs = (r - 1) * (s - 1)

    product_to_C: dict[int, set[int]] = defaultdict(set)

    pair_counts = defaultdict(int)
    product_counts = defaultdict(set)

    formula_failures = 0
    product_formula_failures = 0

    true_a = p % r
    true_b = q % s
    true_z = true_a * true_b

    # --------------------------------------------------------
    # Exhaustive local state space
    # --------------------------------------------------------

    for a in range(1, r):
        for b in range(1, s):

            z = a * b

            alpha_direct, beta_direct, C_direct = \
                exact_C_from_ab(n, a, b, r, s)

            alpha_simple, beta_simple, C_simple = \
                simplified_C(n, a, b, r, s)

            if (
                alpha_direct != alpha_simple
                or beta_direct != beta_simple
                or C_direct != C_simple
            ):
                formula_failures += 1

            C_product = product_only_C(n, z, r, s)

            if C_product != C_direct:
                product_formula_failures += 1

            product_to_C[z].add(C_direct)

            pair_counts[C_direct] += 1
            product_counts[C_direct].add(z)

    # --------------------------------------------------------
    # Check product collapse
    # --------------------------------------------------------

    product_conflicts = sum(
        1
        for values in product_to_C.values()
        if len(values) != 1
    )

    distinct_products = len(product_to_C)

    # Number of distinct z values assigned to each C.
    c0_products = len(product_counts[0])
    c1_products = len(product_counts[1])
    c2_products = len(product_counts[2])

    c0_pairs = pair_counts[0]
    c1_pairs = pair_counts[1]
    c2_pairs = pair_counts[2]

    # Multiplicity:
    # how many (a,b) pairs correspond, on average, to one
    # distinct product z inside each C class?
    avg_mult = {}

    for c in (0, 1, 2):
        if product_counts[c]:
            avg_mult[c] = pair_counts[c] / len(product_counts[c])
        else:
            avg_mult[c] = 0.0

    true_C_direct = exact_C_from_pq(p, q, r, s)[4]
    true_C_simple = product_only_C(n, true_z, r, s)

    return {
        "r": r,
        "s": s,
        "raw_pairs": raw_pairs,
        "distinct_products": distinct_products,

        "formula_failures": formula_failures,
        "product_formula_failures": product_formula_failures,
        "product_conflicts": product_conflicts,

        "c0_pairs": c0_pairs,
        "c1_pairs": c1_pairs,
        "c2_pairs": c2_pairs,

        "c0_products": c0_products,
        "c1_products": c1_products,
        "c2_products": c2_products,

        "avg_mult_c0": avg_mult[0],
        "avg_mult_c1": avg_mult[1],
        "avg_mult_c2": avg_mult[2],

        "true_a": true_a,
        "true_b": true_b,
        "true_z": true_z,
        "true_C": true_C_direct,
        "true_C_product": true_C_simple,
    }


# ------------------------------------------------------------
# Display one cell
# ------------------------------------------------------------

def print_cell(result: dict) -> None:

    r = result["r"]
    s = result["s"]

    compression = (
        result["raw_pairs"] / result["distinct_products"]
        if result["distinct_products"]
        else 0.0
    )

    print(
        f"  ({r:2d},{s:2d}) "
        f"pairs={result['raw_pairs']:5d} "
        f"products={result['distinct_products']:5d} "
        f"collapse={compression:6.2f}x"
    )

    print(
        f"       C=0: "
        f"{result['c0_pairs']:5d} pairs / "
        f"{result['c0_products']:5d} products, "
        f"avg mult={result['avg_mult_c0']:.3f}"
    )

    print(
        f"       C=1: "
        f"{result['c1_pairs']:5d} pairs / "
        f"{result['c1_products']:5d} products, "
        f"avg mult={result['avg_mult_c1']:.3f}"
    )

    print(
        f"       C=2: "
        f"{result['c2_pairs']:5d} pairs / "
        f"{result['c2_products']:5d} products, "
        f"avg mult={result['avg_mult_c2']:.3f}"
    )

    print(
        f"       true: "
        f"a={result['true_a']:2d}, "
        f"b={result['true_b']:2d}, "
        f"ab={result['true_z']:4d}, "
        f"C={result['true_C']}"
    )

    print(
        f"       checks: "
        f"simplified_fail={result['formula_failures']}, "
        f"product_fail={result['product_formula_failures']}, "
        f"product_conflicts={result['product_conflicts']}"
    )


# ------------------------------------------------------------
# Main experiment
# ------------------------------------------------------------

def main() -> None:

    rng = random.Random(SEED)

    print("START EXPERIMENT 162")
    print("=" * 72)
    print("Product Collapse of the Cross-Radix Carry Correction")
    print()
    print("Testing:")
    print()
    print("    beta  == (n - a*b) * r^(-1) mod s")
    print("    alpha == (n - a*b) * s^(-1) mod r")
    print()
    print("and therefore:")
    print()
    print("    C = F(n,r,s,a*b)")
    print()
    print("where")
    print()
    print("    a = p mod r")
    print("    b = q mod s")
    print()
    print("The key question is whether the entire C-relation")
    print("collapses exactly onto the product z = a*b.")
    print("=" * 72)

    total_formula_failures = 0
    total_product_failures = 0
    total_product_conflicts = 0

    for bits in BIT_SIZES:

        print()
        print("-" * 72)
        print(f"BIT SIZE = {bits}")
        print("-" * 72)

        for sample in range(1, SAMPLES_PER_SIZE + 1):

            p, q, n = make_semiprime(bits, rng)

            print()
            print(
                f"Sample {sample}: "
                f"p={p}, q={q}"
            )
            print(f"n={n}")

            sample_formula_failures = 0
            sample_product_failures = 0
            sample_product_conflicts = 0

            results = []

            for r in R1:
                if math.gcd(n, r) != 1:
                    continue

                for s in R2:
                    if math.gcd(n, s) != 1:
                        continue

                    result = analyze_cell(n, p, q, r, s)

                    results.append(result)

                    print_cell(result)

                    sample_formula_failures += \
                        result["formula_failures"]

                    sample_product_failures += \
                        result["product_formula_failures"]

                    sample_product_conflicts += \
                        result["product_conflicts"]

            total_formula_failures += sample_formula_failures
            total_product_failures += sample_product_failures
            total_product_conflicts += sample_product_conflicts

            print()
            print("Sample summary:")
            print(
                f"  cells checked      = {len(results)}"
            )
            print(
                f"  simplified fails   = "
                f"{sample_formula_failures}"
            )
            print(
                f"  product fails      = "
                f"{sample_product_failures}"
            )
            print(
                f"  product conflicts  = "
                f"{sample_product_conflicts}"
            )

    # --------------------------------------------------------
    # Final theorem check
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("FINAL RESULT")
    print("=" * 72)

    print(
        f"Simplified-formula failures : "
        f"{total_formula_failures}"
    )

    print(
        f"Product-formula failures    : "
        f"{total_product_failures}"
    )

    print(
        f"Product conflicts           : "
        f"{total_product_conflicts}"
    )

    if (
        total_formula_failures == 0
        and total_product_failures == 0
        and total_product_conflicts == 0
    ):
        print()
        print("SUCCESS")
        print()
        print("Every tested C-cell is exactly a function of z = a*b.")
        print()
        print("Therefore the local carry correction has the exact form")
        print()
        print("    C_rs(a,b) = F_{n,r,s}(a*b)")
        print()
        print("The apparent 2-dimensional carry relation is therefore")
        print("a one-dimensional product-level relation in (a,b)-space.")

    else:
        print()
        print("FAILURE")
        print()
        print("The proposed product collapse was not universally valid")
        print("over the tested state spaces.")

    print()
    print("=" * 72)
    print("END OF EXPERIMENT 162")
    print("=" * 72)

    print("FINISHED EXPERIMENT 162")


if __name__ == "__main__":
    main()
