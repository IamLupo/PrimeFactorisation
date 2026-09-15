#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 440
========================================================================================================================

BOUNDED-K PULLBACK / X-SPACE COLLAPSE AUDIT

QUESTION
  Can a finite K interval be pulled back through

      K(x) = (r - x^2 - 2*d0*x) / 4

  so that the resulting admissible x interval becomes small enough
  to make K reconstruction practical?

IMPORTANT
  This experiment NEVER materializes a giant K set.

  It also does NOT construct the CRT product image.

  Instead:

    1. choose a hidden integer K_true
    2. choose x_true
    3. derive r so that the exact relation holds
    4. for a bounded K interval [-2^B, 2^B):
         solve the quadratic inequality for x
    5. count the integer x values in that interval
    6. apply the 2-adic K-image condition directly
    7. apply odd-prime K-image conditions directly
    8. count the surviving x values
    9. reconstruct K exactly from each surviving x
   10. verify the hidden K and x remain

  No giant residue set is constructed.

RULES
  exact integer arithmetic only
  no floating-point arithmetic for acceptance
  no resultants
  no Groebner basis
  no symbolic factorization
  modular conditions are necessary only
  exact K(x) reconstruction is authoritative
  exact quadratic relation is authoritative

CORE
  f(x,K) = x^2 + 2*d0*x + 4*K-r

  K(x) = (r-x^2-2*d0*x)/4

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import Iterable


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

K_BITS = [32, 40, 48, 56, 64, 72, 80]
U_VALUES = [8, 10, 12, 14, 16, 18, 20]

ODD_PRIMES = [3, 5, 7, 11, 13, 17, 19]

# Never enumerate a huge x interval.
MAX_X_ENUMERATION = 250_000

# For very large x intervals use only exact counting.
MATERIALIZE_X_THRESHOLD = 50_000

# Deterministic hidden-K multipliers for the audit.
# These are NOT used by the reconstruction logic.
HIDDEN_K_FACTORS = [17, 29, 43, 61, 89, 127, 191, 257]


# ----------------------------------------------------------------------------------------------------------------------
# INSTANCE DATA
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Instance:
    idx: int
    p: int
    q: int
    N: int
    S: int
    d: int
    x_true: int
    d0: int


INSTANCE_BASE = [
    (50411, 282599, -10),
    (1013, 10009, -17),
    (10009, 1000033, -24),
    (10009, 10037, -31),
    (50023, 50051, -38),
    (100019, 100043, -45),
    (200009, 200017, -52),
    (300017, 900007, -59),
]


def build_instances() -> list[Instance]:
    out: list[Instance] = []

    for idx, (p, q, x_true) in enumerate(INSTANCE_BASE, start=1):
        N = p * q
        S = p + q
        d = N - 2 * S + 1
        d0 = d - x_true

        # Basic consistency.
        assert d0 + x_true == d

        out.append(
            Instance(
                idx=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                x_true=x_true,
                d0=d0,
            )
        )

    return out


# ----------------------------------------------------------------------------------------------------------------------
# HIDDEN TEST GENERATION
# ----------------------------------------------------------------------------------------------------------------------

def hidden_K(inst: Instance) -> int:
    """
    Deterministic hidden K used only to generate a valid exact relation.

    The reconstruction itself never calls this function.

    We intentionally make |K| large enough that several B values are
    genuinely inside/outside the artificial bounded K range.
    """
    factor = HIDDEN_K_FACTORS[inst.idx - 1]

    scale = (
        inst.d0 * inst.d0
        + 13 * inst.N
        + 7 * inst.S * inst.S
        + factor * inst.d0
    )

    return -(scale * factor + inst.N * (factor + 3))


def build_r(inst: Instance, K_true: int) -> int:
    """
    Choose r so that the exact relation holds for the hidden x_true:

        K_true = (r - x^2 - 2*d0*x) / 4

    This is an algebraic test of the pullback mechanism.
    """
    x = inst.x_true
    return 4 * K_true + x * x + 2 * inst.d0 * x


# ----------------------------------------------------------------------------------------------------------------------
# EXACT K MAP
# ----------------------------------------------------------------------------------------------------------------------

def numerator_for_K(inst: Instance, r: int, x: int) -> int:
    return r - x * x - 2 * inst.d0 * x


def exact_K_from_x(inst: Instance, r: int, x: int) -> int | None:
    num = numerator_for_K(inst, r, x)

    if num % 4 != 0:
        return None

    return num // 4


def exact_relation_holds(inst: Instance, r: int, x: int, K: int) -> bool:
    return x * x + 2 * inst.d0 * x + 4 * K - r == 0


# ----------------------------------------------------------------------------------------------------------------------
# BOUNDED K INTERVAL
# ----------------------------------------------------------------------------------------------------------------------

def k_bounds(B: int) -> tuple[int, int]:
    lo = -(1 << B)
    hi = 1 << B
    return lo, hi


# ----------------------------------------------------------------------------------------------------------------------
# QUADRATIC PULLBACK
# ----------------------------------------------------------------------------------------------------------------------

def x_center(inst: Instance) -> int:
    """
    Center of the parabola

        K(x) = (r - x^2 - 2*d0*x)/4

    is at x = -d0.
    """
    return -inst.d0


def shifted_square_value(inst: Instance, r: int, x: int) -> int:
    """
    Complete the square:

      x^2 + 2*d0*x = (x+d0)^2 - d0^2

    Therefore

      K(x) = (r + d0^2 - (x+d0)^2) / 4.

    Let y = x+d0.
    """
    y = x + inst.d0
    return r + inst.d0 * inst.d0 - y * y


def pullback_y_radius(inst: Instance, r: int, B: int) -> tuple[int, int] | None:
    """
    K(x) in [-2^B, 2^B) means

      -2^B <= (r+d0^2-y^2)/4 < 2^B.

    Rearranging:

      r+d0^2 - 4*2^B < y^2 <= r+d0^2 + 4*2^B.

    This gives an exact annular interval in y^2.

    Return inclusive |y| bounds when non-empty:

      lower_abs <= |y| <= upper_abs.

    """
    k_lo, k_hi = k_bounds(B)

    A = r + inst.d0 * inst.d0

    # K >= k_lo:
    # A - y^2 >= 4*k_lo
    # y^2 <= A - 4*k_lo
    upper_sq = A - 4 * k_lo

    # K < k_hi:
    # A - y^2 < 4*k_hi
    # y^2 > A - 4*k_hi
    lower_sq_exclusive = A - 4 * k_hi

    if upper_sq < 0:
        return None

    upper_abs = isqrt(upper_sq)

    # Need y^2 > lower_sq_exclusive.
    if lower_sq_exclusive < 0:
        lower_abs = 0
    else:
        s = isqrt(lower_sq_exclusive)

        if s * s > lower_sq_exclusive:
            lower_abs = s
        else:
            lower_abs = s + 1

    if lower_abs > upper_abs:
        return None

    return lower_abs, upper_abs


def interval_count_for_radius(lower_abs: int, upper_abs: int) -> int:
    """
    Number of integer y with

        lower_abs <= |y| <= upper_abs.
    """
    if lower_abs > upper_abs:
        return 0

    return 2 * (upper_abs - lower_abs + 1) - (1 if lower_abs == 0 else 0)


def x_interval_count(inst: Instance, r: int, B: int) -> int:
    bounds = pullback_y_radius(inst, r, B)

    if bounds is None:
        return 0

    lo_abs, hi_abs = bounds

    return interval_count_for_radius(lo_abs, hi_abs)


# ----------------------------------------------------------------------------------------------------------------------
# RESIDUE TESTS
# ----------------------------------------------------------------------------------------------------------------------

def k_residue_from_x(inst: Instance, r: int, x: int, modulus: int) -> int | None:
    """
    Return K(x) modulo modulus without constructing the integer K.

    Need divisibility by 4 first because K is an integer.
    """
    num = numerator_for_K(inst, r, x)

    if num % 4 != 0:
        return None

    # Since num/4 is an integer, reduce after exact division.
    return (num // 4) % modulus


def local_k_image(p: int, r: int, d0: int, modulus: int) -> set[int]:
    """
    Generic local image of

        K(x) = (r - x^2 - 2*d0*x)/4

    modulo odd prime p.

    For odd p, 4 is invertible.
    """
    inv4 = pow(4, -1, p)
    out: set[int] = set()

    rr = r % p
    dd = d0 % p

    for x in range(p):
        value = (rr - x * x - 2 * dd * x) % p
        out.add((value * inv4) % p)

    return out


def local_2adic_image(
    inst: Instance,
    r: int,
    u: int,
) -> set[int]:
    """
    Build the small local K image modulo 2^(u-2).

    The x domain modulo 2^u is at most 2^20 here.
    """
    if u < 2:
        raise ValueError("u must be >= 2")

    x_mod = 1 << u
    k_mod = 1 << (u - 2)

    out: set[int] = set()

    for xr in range(x_mod):
        num = numerator_for_K(inst, r, xr)

        if num & 3:
            continue

        out.add((num // 4) % k_mod)

    return out


def odd_prime_image(
    inst: Instance,
    r: int,
    p: int,
) -> set[int]:
    return local_k_image(
        p=p,
        r=r,
        d0=inst.d0,
        modulus=p,
    )


# ----------------------------------------------------------------------------------------------------------------------
# DIRECT LOCAL TESTS
# ----------------------------------------------------------------------------------------------------------------------

def passes_2adic_x(inst: Instance, r: int, x: int, u: int) -> bool:
    k = exact_K_from_x(inst, r, x)

    if k is None:
        return False

    k_modulus = 1 << (u - 2)
    k_res = k % k_modulus

    image = local_2adic_image(inst, r, u)
    return k_res in image


def passes_odd_primes(inst: Instance, r: int, x: int, odd_images: dict[int, set[int]]) -> bool:
    k = exact_K_from_x(inst, r, x)

    if k is None:
        return False

    for p, image in odd_images.items():
        if (k % p) not in image:
            return False

    return True


# ----------------------------------------------------------------------------------------------------------------------
# SMART ENUMERATION OF PULLED-BACK X VALUES
# ----------------------------------------------------------------------------------------------------------------------

def candidate_x_values(
    inst: Instance,
    r: int,
    B: int,
) -> Iterable[int]:
    """
    Enumerate x only when the pulled-back interval is small enough.

    The y=x+d0 substitution makes the interval symmetric around zero.
    """
    bounds = pullback_y_radius(inst, r, B)

    if bounds is None:
        return

    lo_abs, hi_abs = bounds

    total = interval_count_for_radius(lo_abs, hi_abs)

    if total > MAX_X_ENUMERATION:
        return

    # y >= lo_abs
    if lo_abs == 0:
        for y in range(0, hi_abs + 1):
            yield y - inst.d0
    else:
        for y in range(lo_abs, hi_abs + 1):
            yield y - inst.d0
            yield -y - inst.d0


# ----------------------------------------------------------------------------------------------------------------------
# CASE ANALYSIS
# ----------------------------------------------------------------------------------------------------------------------

def run_case(inst: Instance, u: int, B: int, r: int, K_true: int) -> dict:
    k_lo, k_hi = k_bounds(B)

    true_in_range = k_lo <= K_true < k_hi

    bounds = pullback_y_radius(inst, r, B)

    if bounds is None:
        return {
            "interval_count": 0,
            "materialized": False,
            "x_candidates": 0,
            "x_2adic": 0,
            "x_crt": 0,
            "valid_k": 0,
            "true_survives": False,
            "unique_k": False,
            "unique_x": False,
            "true_in_range": true_in_range,
        }

    lo_abs, hi_abs = bounds
    interval_count = interval_count_for_radius(lo_abs, hi_abs)

    materialized = interval_count <= MAX_X_ENUMERATION

    if not materialized:
        return {
            "interval_count": interval_count,
            "materialized": False,
            "x_candidates": 0,
            "x_2adic": 0,
            "x_crt": 0,
            "valid_k": 0,
            "true_survives": None,
            "unique_k": None,
            "unique_x": None,
            "true_in_range": true_in_range,
        }

    # Local 2-adic image is small enough to be constructed exactly.
    image_2 = local_2adic_image(inst, r, u)

    odd_images = {
        p: odd_prime_image(inst, r, p)
        for p in ODD_PRIMES
    }

    x_candidates = []
    x_2adic = []
    x_crt = []
    valid_k_values: set[int] = set()
    valid_x_values: set[int] = set()

    for x in candidate_x_values(inst, r, B):
        x_candidates.append(x)

        Kx = exact_K_from_x(inst, r, x)
        if Kx is None:
            continue

        k_mod_2 = Kx % (1 << (u - 2))

        if k_mod_2 not in image_2:
            continue

        x_2adic.append(x)

        if not passes_odd_primes(inst, r, x, odd_images):
            continue

        x_crt.append(x)

        # Authoritative exact bounded-K test.
        if not (k_lo <= Kx < k_hi):
            continue

        valid_k_values.add(Kx)
        valid_x_values.add(x)

    true_survives = False

    if true_in_range:
        true_survives = (
            inst.x_true in valid_x_values
            and K_true in valid_k_values
            and exact_relation_holds(inst, r, inst.x_true, K_true)
        )

    return {
        "interval_count": interval_count,
        "materialized": True,
        "x_candidates": len(x_candidates),
        "x_2adic": len(x_2adic),
        "x_crt": len(x_crt),
        "valid_k": len(valid_k_values),
        "valid_x": len(valid_x_values),
        "true_survives": true_survives,
        "unique_k": len(valid_k_values) == 1,
        "unique_x": len(valid_x_values) == 1,
        "true_in_range": true_in_range,
        "k_values": valid_k_values,
        "x_values": valid_x_values,
        "lo_abs": lo_abs,
        "hi_abs": hi_abs,
        "two_adic_image_size": len(image_2),
    }


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 440")
    print("=" * 120)
    print()
    print("BOUNDED-K PULLBACK / X-SPACE COLLAPSE AUDIT")
    print()
    print("QUESTION")
    print("  Can a finite K interval be pulled back through K(x)")
    print("  so that the resulting admissible x-space becomes small?")
    print()
    print("K_BITS =", K_BITS)
    print("U_VALUES =", U_VALUES)
    print("ODD_PRIMES =", ODD_PRIMES)
    print("MAX_X_ENUMERATION =", MAX_X_ENUMERATION)
    print()

    instances = build_instances()

    global_cases = 0
    materialized_cases = 0
    true_survival_cases = 0
    unique_k_cases = 0
    unique_x_cases = 0

    aggregate = {
        B: {
            "cases": 0,
            "materialized": 0,
            "interval": 0,
            "x2": 0,
            "xcrt": 0,
            "unique_k": 0,
            "unique_x": 0,
            "true_survive": 0,
        }
        for B in K_BITS
    }

    for inst in instances:
        K_true = hidden_K(inst)
        r = build_r(inst, K_true)

        print("-" * 120)
        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} S={inst.S} d={inst.d} "
            f"d0={inst.d0} x_true={inst.x_true}"
        )

        print("  hidden K is used only to generate r; reconstruction never reads it")
        print()

        for B in K_BITS:
            global_cases += 1
            aggregate[B]["cases"] += 1

            result = run_case(
                inst=inst,
                u=max(U_VALUES),
                B=B,
                r=r,
                K_true=K_true,
            )

            global_cases += 0

            if not result["materialized"]:
                print(
                    f"  B={B:2d}"
                    f" interval={result['interval_count']:>20}"
                    f" materialized=NO"
                    f" true_in_range={str(result['true_in_range']):>5}"
                )
                continue

            materialized_cases += 1
            aggregate[B]["materialized"] += 1
            aggregate[B]["interval"] += result["interval_count"]
            aggregate[B]["x2"] += result["x_2adic"]
            aggregate[B]["xcrt"] += result["x_crt"]

            if result["unique_k"]:
                unique_k_cases += 1
                aggregate[B]["unique_k"] += 1

            if result["unique_x"]:
                unique_x_cases += 1
                aggregate[B]["unique_x"] += 1

            if result["true_survives"]:
                true_survival_cases += 1
                aggregate[B]["true_survive"] += 1

            print(
                f"  B={B:2d}"
                f" interval={result['interval_count']:>9}"
                f" 2adic={result['x_2adic']:>8}"
                f" crt={result['x_crt']:>8}"
                f" K={result['valid_k']:>8}"
                f" X={result['valid_x']:>8}"
                f" true={str(result['true_in_range']):>5}"
                f" survive={str(result['true_survives']):>5}"
                f" uniqueK={str(result['unique_k']):>5}"
                f" uniqueX={str(result['unique_x']):>5}"
            )

            # For a unique result, print it explicitly.
            if result["unique_k"] and result["k_values"]:
                only_k = next(iter(result["k_values"]))
                print(f"      UNIQUE K = {only_k}")

            if result["unique_x"] and result["x_values"]:
                only_x = next(iter(result["x_values"]))
                print(f"      UNIQUE x = {only_x}")

        print()

    # ------------------------------------------------------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL PULLBACK SUMMARY")
    print("=" * 120)
    print(
        f"  total cases                  = {global_cases}"
    )
    print(
        f"  materialized x intervals    = {materialized_cases}/{global_cases}"
    )
    print(
        f"  true K/x survival cases      = {true_survival_cases}"
    )
    print(
        f"  unique K reconstructions     = {unique_k_cases}"
    )
    print(
        f"  unique x reconstructions     = {unique_x_cases}"
    )
    print()

    print("BY K BIT BOUND")
    print(
        "  B   cases materialized avg-interval avg-2adic-x avg-crt-x "
        "true-survive uniqueK uniqueX"
    )
    print("-" * 120)

    for B in K_BITS:
        a = aggregate[B]

        mat = a["materialized"]

        avg_interval = (
            a["interval"] / mat
            if mat
            else 0.0
        )

        avg_x2 = (
            a["x2"] / mat
            if mat
            else 0.0
        )

        avg_xcrt = (
            a["xcrt"] / mat
            if mat
            else 0.0
        )

        print(
            f"{B:3d}"
            f" {a['cases']:5d}"
            f" {a['materialized']:11d}"
            f" {avg_interval:12.2f}"
            f" {avg_x2:12.2f}"
            f" {avg_xcrt:11.2f}"
            f" {a['true_survive']:12d}"
            f" {a['unique_k']:7d}"
            f" {a['unique_x']:7d}"
        )

    # ------------------------------------------------------------------------------------------------------------------
    # STRUCTURAL CHECKS
    # ------------------------------------------------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    all_true_relation = True
    all_true_integrality = True

    for inst in instances:
        K_true = hidden_K(inst)
        r = build_r(inst, K_true)

        recovered = exact_K_from_x(
            inst,
            r,
            inst.x_true,
        )

        if recovered != K_true:
            all_true_relation = False

        if numerator_for_K(inst, r, inst.x_true) % 4 != 0:
            all_true_integrality = False

    print(
        f"  hidden x_true -> exact K relation = {all_true_relation}"
    )
    print(
        f"  true K integrality checks          = {all_true_integrality}"
    )
    print(
        "  no giant K interval was materialized = True"
    )
    print(
        "  no CRT Cartesian product was built   = True"
    )

    # ------------------------------------------------------------------------------------------------------------------
    # INTERPRETATION
    # ------------------------------------------------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)
    print()
    print(
        "  The key quantity is the size of the pulled-back x population."
    )
    print()
    print(
        "  A bounded K interval is translated into an exact quadratic"
    )
    print(
        "  inequality in y = x+d0."
    )
    print()
    print(
        "  The experiment then measures whether:"
    )
    print(
        "      bounded K interval"
    )
    print(
        "          -> x interval"
    )
    print(
        "          -> 2-adic K residue filter"
    )
    print(
        "          -> odd-prime K residue filter"
    )
    print(
        "          -> exact K(x)"
    )
    print(
        "  produces a small population."
    )
    print()
    print(
        "  IMPORTANT:"
    )
    print(
        "  A small x population is useful only if the bound on K is"
    )
    print(
        "  independently justified. The hidden K generator in this"
    )
    print(
        "  experiment is only a controlled structural test."
    )
    print()
    print(
        "  This remains a candidate-isolation experiment, not a"
    )
    print(
        "  factorization theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 440 FINAL STATUS")
    print("=" * 120)
    print(
        f"  CASES = {global_cases}"
    )
    print(
        f"  MATERIALIZED CASES = {materialized_cases}"
    )
    print(
        f"  TRUE SURVIVAL CASES = {true_survival_cases}"
    )
    print(
        f"  UNIQUE K CASES = {unique_k_cases}"
    )
    print(
        f"  UNIQUE X CASES = {unique_x_cases}"
    )
    print(
        f"  ALL BASIC EXACT CHECKS = {all_true_relation and all_true_integrality}"
    )
    print("=" * 120)
    print("EXPERIMENT 440 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
