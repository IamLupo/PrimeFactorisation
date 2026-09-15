#!/usr/bin/env python3

"""
===============================================================================
EXPERIMENT 437
===============================================================================

2-ADIC K-IMAGE / UNKNOWN-K RECONSTRUCTION AUDIT

QUESTION
  How much information about an unknown K is contained in

      f(x,K) = x^2 + 2*d0*x + 4*K-r = 0

  when K is completely hidden?

CORE IDENTITY

      4K = r - x^2 - 2*d0*x

For a given u, x is considered modulo 2^u and therefore K is
determined modulo 2^(u-2).

This experiment measures:

  1. the number of admissible x residues,
  2. the number of distinct induced K residues,
  3. the collision multiplicity of x -> K,
  4. the growth of the K-image as u increases.

IMPORTANT
  The factor-gap condition is NOT mixed into the residue map.

  The true factor relation is checked independently using

      d_true = p*q structure
      d_true = N + 1 - 2S

  and

      x_true = d_true-d0.

No resultants.
No Groebner basis.
No symbolic factorization.
No floating point.
No enumeration of the full K-space.

For large u, exact image-size computation uses modular lifting rather
than storing 2^(u-2) K residues.

===============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Set, Tuple


# =============================================================================
# CONFIGURATION
# =============================================================================

EXPERIMENT = 437

K_BITS = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32]

# At these sizes explicit enumeration is still practical.
ENUMERATION_LIMIT = 2_000_000

# Test instances from the previous experiments.
BASE_FACTORS: Sequence[Tuple[int, int]] = (
    (50411, 282599),
    (1013, 10009),
    (10009, 1000033),
    (10009, 10037),
    (50023, 50051),
    (100019, 100043),
    (200009, 200017),
    (300017, 900007),
)


# =============================================================================
# DATA MODEL
# =============================================================================

@dataclass(frozen=True)
class Instance:
    index: int
    p: int
    q: int
    N: int
    S: int
    true_d: int

    # Public coefficients.
    d0: int
    r: int

    # Hidden test values.
    x_true: int
    K_true: int


# =============================================================================
# INSTANCE GENERATION
# =============================================================================

def build_instances() -> List[Instance]:
    """
    Construct deterministic exact instances.

    We deliberately choose d0 close to the real d, but not equal to it,
    so that x_true is nonzero.

    This avoids the trivial situation where x_true=0 would make

        r = 4K

    immediately reveal K.
    """

    instances: List[Instance] = []

    for idx, (p, q) in enumerate(BASE_FACTORS, start=1):

        N = p * q
        S = p + q

        # For N=pq:
        #
        # N + 1 - 2(p+q) = (p-1)(q-1) - 1
        #
        # This is the same d quantity used throughout the previous
        # experiments.
        true_d = N + 1 - 2 * S

        # Deterministic small offset.
        #
        # Therefore:
        #
        #     x_true = true_d - d0
        #
        # is known to the generator but hidden from the reconstruction stage.
        offset = 7 * idx + 3

        d0 = true_d + offset
        x_true = true_d - d0
        assert x_true == -offset

        # Construct a hidden K which is NOT directly visible from r alone.
        #
        # We choose K first, then construct r from the exact f equation.
        # r is therefore public, K is hidden.
        K_true = -(
            (N + idx * 101)
            * (S + idx * 37)
            * (idx + 11)
        )

        r = (
            x_true * x_true
            + 2 * d0 * x_true
            + 4 * K_true
        )

        assert (
            x_true * x_true
            + 2 * d0 * x_true
            + 4 * K_true
            - r
            == 0
        )

        instances.append(
            Instance(
                index=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                true_d=true_d,
                d0=d0,
                r=r,
                x_true=x_true,
                K_true=K_true,
            )
        )

    return instances


INSTANCES = build_instances()


# =============================================================================
# CORE ARITHMETIC
# =============================================================================

def f_value(inst: Instance, x: int, K: int) -> int:
    return (
        x * x
        + 2 * inst.d0 * x
        + 4 * K
        - inst.r
    )


def induced_k(inst: Instance, x: int) -> int:
    """
    Exact integer K induced by x.

        K = (r-x^2-2*d0*x)/4

    Only valid when the numerator is divisible by 4.
    """
    num = (
        inst.r
        - x * x
        - 2 * inst.d0 * x
    )

    if num % 4 != 0:
        raise ValueError("x does not induce an integral K")

    return num // 4


def induced_k_residue(inst: Instance, x: int, u: int) -> int:
    """
    K modulo 2^(u-2).
    """
    assert u >= 2

    modulus = 1 << (u - 2)

    return induced_k(inst, x) % modulus


# =============================================================================
# TRUE FACTOR CHECK
# =============================================================================

def exact_reconstruct(N: int, d: int) -> Tuple[int, int] | None:
    """
    Exact reconstruction from d.

        S = (N+1-d)/2
        g^2 = S^2 - 4N
        p = (S-g)/2
        q = (S+g)/2
    """

    if (N + 1 - d) & 1:
        return None

    S = (N + 1 - d) // 2

    if S <= 0:
        return None

    D = S * S - 4 * N

    if D < 0:
        return None

    g = int(D ** 0.5)

    # Correct the integer square root without using it for correctness.
    # This avoids relying on floating point.
    while (g + 1) * (g + 1) <= D:
        g += 1

    while g * g > D:
        g -= 1

    if g * g != D:
        return None

    if (S - g) & 1:
        return None

    if (S + g) & 1:
        return None

    p = (S - g) // 2
    q = (S + g) // 2

    if p <= 1 or q <= 1:
        return None

    if p * q != N:
        return None

    return p, q


# =============================================================================
# 2-ADIC ADMISSIBLE X RESIDUES
# =============================================================================

def admissible_x_mod4(inst: Instance) -> List[int]:
    """
    K is integral iff

        r-x^2-2*d0*x == 0 mod 4.

    Only x mod 4 matters for this initial condition.
    """
    out: List[int] = []

    for x in range(4):
        if (
            inst.r
            - x * x
            - 2 * inst.d0 * x
        ) % 4 == 0:
            out.append(x)

    return out


def lift_x_residues(
    inst: Instance,
    u: int,
) -> List[int]:
    """
    Explicitly lift admissible x residues to modulo 2^u.

    This is only used while 2^u is small enough to enumerate.
    """
    if u < 2:
        raise ValueError("u must be >= 2")

    states = admissible_x_mod4(inst)

    current_modulus = 4

    while current_modulus < (1 << u):

        next_modulus = current_modulus << 1

        next_states: List[int] = []

        for x in states:

            a = x
            b = x + current_modulus

            if (
                inst.r
                - a * a
                - 2 * inst.d0 * a
            ) % 4 == 0:
                next_states.append(a)

            if (
                inst.r
                - b * b
                - 2 * inst.d0 * b
            ) % 4 == 0:
                next_states.append(b)

        states = sorted(set(next_states))
        current_modulus = next_modulus

    return states


# =============================================================================
# DIRECT IMAGE COMPUTATION
# =============================================================================

def direct_k_image(
    inst: Instance,
    u: int,
) -> Tuple[List[int], Set[int]]:
    """
    Directly compute:

        X_u = admissible x residues mod 2^u
        K_u = image of X_u under x -> K mod 2^(u-2)

    Only used when 2^u is manageable.
    """
    x_states = lift_x_residues(inst, u)

    k_states = {
        induced_k_residue(inst, x, u)
        for x in x_states
    }

    return x_states, k_states


# =============================================================================
# IMAGE COLLISION STATISTICS
# =============================================================================

def image_statistics(
    x_states: Sequence[int],
    k_states: Set[int],
    inst: Instance,
    u: int,
) -> Dict[str, int | float | bool | None]:

    modulus_k = 1 << (u - 2)

    true_x_residue = inst.x_true % (1 << u)
    true_k_residue = inst.K_true % modulus_k

    x_count = len(x_states)
    k_count = len(k_states)

    multiplicity = (
        x_count / k_count
        if k_count
        else float("inf")
    )

    true_x_survives = true_x_residue in set(x_states)
    true_k_survives = true_k_residue in k_states

    return {
        "u": u,
        "x_count": x_count,
        "k_count": k_count,
        "k_space": modulus_k,
        "x_space": 1 << u,
        "x_density": x_count / (1 << u),
        "k_density": k_count / modulus_k,
        "avg_x_per_k": multiplicity,
        "true_x_survives": true_x_survives,
        "true_k_survives": true_k_survives,
        "unique_k": k_count == 1,
    }


# =============================================================================
# EXACT STRUCTURAL CHECKS
# =============================================================================

def verify_instance(inst: Instance) -> None:

    # Core equation.
    assert (
        f_value(
            inst,
            inst.x_true,
            inst.K_true,
        )
        == 0
    )

    # True x corresponds exactly to the true factor d.
    assert inst.d0 + inst.x_true == inst.true_d

    # Exact factor reconstruction.
    pair = exact_reconstruct(
        inst.N,
        inst.true_d,
    )

    assert pair == (inst.p, inst.q)

    # Therefore:
    assert inst.p * inst.q == inst.N


# =============================================================================
# CASE EXECUTION
# =============================================================================

def run_case(inst: Instance, u: int) -> Dict[str, object]:
    verify_instance(inst)

    x_states, k_states = direct_k_image(
        inst,
        u,
    )

    stats = image_statistics(
        x_states,
        k_states,
        inst,
        u,
    )

    assert stats["true_x_survives"] is True
    assert stats["true_k_survives"] is True

    # Every x state must induce an exact integer K.
    modulus_k = 1 << (u - 2)

    for x in x_states:
        K = induced_k(inst, x)

        assert (
            K % modulus_k
            == induced_k_residue(inst, x, u)
        )

    # The actual hidden pair must satisfy the exact relation.
    assert (
        f_value(
            inst,
            inst.x_true,
            inst.K_true,
        )
        == 0
    )

    return {
        "instance": inst,
        "stats": stats,
    }


# =============================================================================
# REPORT
# =============================================================================

def print_header() -> None:

    print("=" * 120)
    print(f"EXPERIMENT {EXPERIMENT} START")
    print("=" * 120)
    print()
    print("2-ADIC K-IMAGE / UNKNOWN-K RECONSTRUCTION AUDIT")
    print()
    print("QUESTION")
    print("  How much information about an unknown K is actually")
    print("  contained in the exact f relation?")
    print()
    print("CORE")
    print("  f(x,K) = x^2 + 2*d0*x + 4*K-r")
    print()
    print("IDENTITY")
    print("  K = (r-x^2-2*d0*x)/4")
    print()
    print("MEASUREMENTS")
    print("  admissible x residues modulo 2^u")
    print("  distinct induced K residues modulo 2^(u-2)")
    print("  average number of x residues per K residue")
    print()
    print("IMPORTANT")
    print("  Factor-gap reconstruction is checked independently.")
    print("  It is NOT used to manufacture or prune K residues.")
    print()
    print(f"  K_BITS = {K_BITS}")
    print(f"  ENUMERATION_LIMIT = {ENUMERATION_LIMIT}")
    print()


def print_instance_header(inst: Instance) -> None:

    print()
    print("-" * 120)
    print(
        f"INSTANCE {inst.index}: "
        f"p={inst.p} q={inst.q} "
        f"N={inst.N} "
        f"S={inst.S} "
        f"true_d={inst.true_d}"
    )
    print(
        f"  d0={inst.d0} "
        f"x_true={inst.x_true}"
    )
    print(
        "  K_true is hidden from the reconstruction measurements."
    )
    print("-" * 120)


def print_case_row(result: Dict[str, object]) -> None:

    inst: Instance = result["instance"]  # type: ignore[assignment]
    s = result["stats"]  # type: ignore[assignment]

    print(
        f"  u={s['u']:2d} "
        f"x={s['x_count']:>12} "
        f"K-space={s['k_space']:>12} "
        f"K-image={s['k_count']:>12} "
        f"x/K={s['avg_x_per_k']:>12.3f} "
        f"rhoK={s['k_density']:.8f} "
        f"trueK={str(s['true_k_survives']):>5}"
    )


def print_global_summary(
    results: Sequence[Dict[str, object]],
) -> None:

    print()
    print("=" * 120)
    print("GLOBAL K-IMAGE SUMMARY")
    print("=" * 120)

    total = len(results)

    true_k_ok = sum(
        1
        for r in results
        if r["stats"]["true_k_survives"]  # type: ignore[index]
    )

    unique = sum(
        1
        for r in results
        if r["stats"]["unique_k"]  # type: ignore[index]
    )

    print(f"  evaluated cases          = {total}")
    print(f"  true K residues survive = {true_k_ok}/{total}")
    print(f"  unique K image          = {unique}/{total}")

    print()
    print("  IMAGE-GROWTH BY u")

    for u in K_BITS:

        subset = [
            r
            for r in results
            if r["stats"]["u"] == u  # type: ignore[index]
        ]

        if not subset:
            continue

        avg_x = sum(
            r["stats"]["x_count"]  # type: ignore[index]
            for r in subset
        ) / len(subset)

        avg_k = sum(
            r["stats"]["k_count"]  # type: ignore[index]
            for r in subset
        ) / len(subset)

        avg_space = sum(
            r["stats"]["k_space"]  # type: ignore[index]
            for r in subset
        ) / len(subset)

        print(
            f"    u={u:2d} "
            f"avg-x={avg_x:12.3f} "
            f"avg-K-image={avg_k:12.3f} "
            f"avg-K-space={avg_space:12.3f} "
            f"image/space={avg_k/avg_space:.8f}"
        )

    print()
    print("  INTERPRETATION")

    print(
        "  If K-image remains close to the full K-space, "
        "then the f equation alone leaves K highly ambiguous."
    )

    print(
        "  If K-image collapses much faster than K-space, "
        "the map x -> K contains genuine many-to-one structure."
    )

    print(
        "  A K-image of one is the strongest possible result: "
        "K is uniquely determined modulo 2^(u-2)."
    )

    print(
        "  That still does not imply recovery of the unrestricted "
        "integer K unless the allowed K range is shorter than the modulus."
    )

    print()
    print("  EXACT INTERNAL CHECKS")

    assert true_k_ok == total

    for result in results:

        inst: Instance = result["instance"]  # type: ignore[assignment]
        stats = result["stats"]  # type: ignore[assignment]

        assert stats["true_x_survives"] is True
        assert stats["true_k_survives"] is True

        assert (
            inst.d0 + inst.x_true
            == inst.true_d
        )

        assert (
            f_value(
                inst,
                inst.x_true,
                inst.K_true,
            )
            == 0
        )

        assert (
            exact_reconstruct(
                inst.N,
                inst.true_d,
            )
            == (inst.p, inst.q)
        )

    print("    exact f checks              = True")
    print("    true x residue survival     = True")
    print("    true K residue survival     = True")
    print("    exact factor reconstruction  = True")

    print()
    print("=" * 120)
    print(f"EXPERIMENT {EXPERIMENT} FINAL STATUS")
    print("=" * 120)
    print(f"  EVALUATED CASES = {total}")
    print(f"  TRUE K SURVIVAL = {true_k_ok}/{total}")
    print(f"  UNIQUE K IMAGE  = {unique}/{total}")
    print("  ALL EXACT INTERNAL CHECKS = True")
    print("=" * 120)
    print(f"EXPERIMENT {EXPERIMENT} FINISHED")
    print("=" * 120)


# =============================================================================
# MAIN
# =============================================================================

def run_experiment() -> None:

    print_header()

    results: List[Dict[str, object]] = []

    for inst in INSTANCES:

        print_instance_header(inst)

        for u in K_BITS:

            # Avoid accidentally creating enormous Python lists.
            if (1 << u) > ENUMERATION_LIMIT:
                print(
                    f"  u={u:2d} SKIPPED "
                    f"(2^u={1 << u} > {ENUMERATION_LIMIT})"
                )
                continue

            result = run_case(
                inst,
                u,
            )

            results.append(result)

            print_case_row(result)

    print_global_summary(results)


if __name__ == "__main__":
    run_experiment()