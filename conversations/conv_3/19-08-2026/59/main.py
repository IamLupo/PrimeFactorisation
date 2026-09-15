#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 436
==============================================================================

JOINT 2^u / RESIDUE-WHEEL DENSITY AND INDEPENDENCE AUDIT

QUESTION
  Are the exact 2^u f-congruence condition and the quadratic-residue
  wheel approximately independent, or is there strong arithmetic
  correlation between them?

CORE MEASUREMENTS

  interval density:
      D_interval = number of d-values / interval size

  f density:
      D_f = # {d : f-condition holds} / interval

  wheel density:
      D_w = # {d : wheel condition holds} / interval

  joint density:
      D_fw = # {d : f-condition AND wheel} / interval

  independence prediction:
      D_expected = D_f * D_w

  correlation ratio:
      C = D_fw / (D_f * D_w)

  C ~ 1:
      approximately independent.

  C >> 1:
      positive correlation.

  C << 1:
      negative correlation.

IMPORTANT
  The wheel is only a necessary condition.
  Exact gap-square reconstruction remains authoritative.

RULES
  exact integer arithmetic only
  no floating point in arithmetic predicates
  no resultants
  no Groebner basis
  no symbolic factorization
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import Iterable


# ============================================================================
# CONFIGURATION
# ============================================================================

UNKNOWN_BITS = [
    20, 24, 28, 32, 36, 40, 44, 48, 52, 56
]

WHEEL_PREFIXES = [
    [],
    [3],
    [3, 5],
    [3, 5, 7],
    [3, 5, 7, 11],
    [3, 5, 7, 11, 13],
    [3, 5, 7, 11, 13, 17],
    [3, 5, 7, 11, 13, 17, 19],
]

MAX_INTERVAL = 2_000_000


# ============================================================================
# TEST INSTANCES
# ============================================================================

@dataclass(frozen=True)
class Instance:
    index: int
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q

    @property
    def gap(self) -> int:
        return self.q - self.p

    @property
    def d(self) -> int:
        return self.n + 1 - 2 * self.s


INSTANCES = [
    Instance(1, 50411, 282599),
    Instance(2, 1013, 10009),
    Instance(3, 10009, 1000033),
    Instance(4, 10009, 10037),
    Instance(5, 50023, 50051),
    Instance(6, 100019, 100043),
    Instance(7, 200009, 200017),
    Instance(8, 300017, 900007),
]


# ============================================================================
# EXACT K
# ============================================================================

def exact_K(inst: Instance) -> int:
    """
    Exact invariant used by the preceding experiments.
    """
    n = inst.n
    s = inst.s
    return -n * (s * s - 2 * n)


def centered_residue(value: int, modulus: int) -> int:
    r = value % modulus
    if r >= modulus // 2:
        r -= modulus
    return r


# ============================================================================
# INTERVAL MODEL
# ============================================================================

def candidate_interval(inst: Instance, u: int) -> tuple[int, int]:
    """
    Reconstruct the same transformed d uncertainty interval used by
    Experiment 435.

    This routine is intentionally isolated so that the interval definition
    can be changed without changing any of the density tests.
    """

    K = exact_K(inst)
    M = 1 << u

    K0 = centered_residue(K, M)
    kq = (K - K0) // M

    # Preserve the same empirical interval law used in the preceding
    # corrected experiment.
    radius = min(abs(kq), MAX_INTERVAL // 2)

    d_true = inst.d

    lo = d_true - radius
    hi = d_true + radius

    return lo, hi


# ============================================================================
# EXACT f CONDITION
# ============================================================================

def f_condition(inst: Instance, d: int, u: int) -> bool:
    """
    Exact 2^u congruence condition.

    Candidate d determines

        S = (N+1-d)/2

    and therefore

        K(d) = -N(S^2 - 2N).

    The candidate satisfies the known K information iff

        K(d) == K(true) mod 2^u.
    """

    n = inst.n

    if ((n + 1 - d) & 1) != 0:
        return False

    s = (n + 1 - d) // 2

    Kd = -n * (s * s - 2 * n)
    Ktrue = exact_K(inst)

    return (Kd - Ktrue) % (1 << u) == 0


# ============================================================================
# QUADRATIC-RESIDUE WHEEL
# ============================================================================

def quadratic_residues(p: int) -> set[int]:
    return {(x * x) % p for x in range(p)}


def wheel_modulus(primes: Iterable[int]) -> int:
    m = 1
    for p in primes:
        m *= p
    return m


def wheel_accepts(inst: Instance, d: int, primes: list[int]) -> bool:
    """
    Necessary condition for

        (N+1-d)^2 - 16N = 4g^2.
    """

    if not primes:
        return True

    n = inst.n

    for p in primes:
        qr = quadratic_residues(p)
        disc = ((n + 1 - d) * (n + 1 - d) - 16 * n) % p

        if disc not in qr:
            return False

    return True


# ============================================================================
# EXACT RECONSTRUCTION
# ============================================================================

def reconstruct(inst: Instance, d: int) -> tuple[int, int] | None:
    n = inst.n

    if ((n + 1 - d) & 1) != 0:
        return None

    s = (n + 1 - d) // 2

    disc = s * s - 4 * n

    if disc < 0:
        return None

    g = isqrt(disc)

    if g * g != disc:
        return None

    if (s - g) & 1:
        return None

    p = (s - g) // 2
    q = (s + g) // 2

    if p <= 0 or q <= 0:
        return None

    if p * q != n:
        return None

    return p, q


# ============================================================================
# ONE CASE
# ============================================================================

def run_case(
    inst: Instance,
    u: int,
) -> dict | None:

    lo, hi = candidate_interval(inst, u)
    interval = hi - lo + 1

    if interval > MAX_INTERVAL:
        return None

    # Precompute the f set once.
    f_set: set[int] = set()

    for d in range(lo, hi + 1):
        if f_condition(inst, d, u):
            f_set.add(d)

    # The true d must always survive.
    assert inst.d in f_set

    # Exact baseline reconstruction.
    exact_set = {
        d for d in f_set
        if reconstruct(inst, d) is not None
    }

    assert exact_set == {inst.d}

    wheel_data = []

    for primes in WHEEL_PREFIXES:

        wheel_set: set[int] = set()

        for d in range(lo, hi + 1):
            if wheel_accepts(inst, d, primes):
                wheel_set.add(d)

        joint_set = f_set & wheel_set

        exact_joint = {
            d for d in joint_set
            if reconstruct(inst, d) is not None
        }

        assert inst.d in wheel_set
        assert inst.d in joint_set
        assert exact_joint == {inst.d}

        wheel_data.append(
            {
                "primes": primes,
                "wheel_hits": len(wheel_set),
                "joint_hits": len(joint_set),
                "wheel_density_num": len(wheel_set),
                "joint_density_num": len(joint_set),
            }
        )

    return {
        "interval": interval,
        "f_count": len(f_set),
        "exact_count": len(exact_set),
        "wheels": wheel_data,
    }


# ============================================================================
# REPORT
# ============================================================================

def print_header() -> None:
    print("=" * 120)
    print("EXPERIMENT 436 START")
    print("=" * 120)
    print()
    print("JOINT 2^u / RESIDUE-WHEEL DENSITY AND INDEPENDENCE AUDIT")
    print()
    print("QUESTION")
    print("  Are the exact f-congruence and quadratic-residue wheel")
    print("  approximately independent?")
    print()


def main() -> None:

    print_header()

    print("CONFIGURATION")
    print(f"  instances       = {len(INSTANCES)}")
    print(f"  K bits          = {UNKNOWN_BITS}")
    print(f"  wheel prefixes  = {WHEEL_PREFIXES}")
    print(f"  max interval    = {MAX_INTERVAL}")
    print()

    results = []

    # ------------------------------------------------------------------------
    # Compact table
    # ------------------------------------------------------------------------

    print("=" * 120)
    print("COMPACT JOINT-DENSITY SUMMARY")
    print("=" * 120)

    print(
        f"{'i':>2} {'u':>3} {'gap':>8} {'interval':>10} "
        f"{'f':>10} {'f-density':>12}"
    )

    for inst in INSTANCES:
        for u in UNKNOWN_BITS:

            result = run_case(inst, u)

            if result is None:
                print(
                    f"{inst.index:2d} {u:3d} {inst.gap:8d} "
                    f"{'SKIPPED':>10}"
                )
                continue

            results.append((inst, u, result))

            interval = result["interval"]
            f_count = result["f_count"]

            print(
                f"{inst.index:2d} {u:3d} "
                f"{inst.gap:8d} {interval:10d} "
                f"{f_count:10d} "
                f"{f_count / interval:12.8f}"
            )

    # ------------------------------------------------------------------------
    # Joint density tables
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("JOINT DENSITY / INDEPENDENCE RATIO")
    print("=" * 120)

    for wi, primes in enumerate(WHEEL_PREFIXES):

        print()
        print(f"WHEEL = {primes}")

        print(
            f"{'i':>2} {'u':>3} {'interval':>10} "
            f"{'f':>8} {'wheel':>10} {'joint':>10} "
            f"{'expected':>12} {'corr':>12}"
        )

        for inst, u, result in results:

            interval = result["interval"]
            f_count = result["f_count"]

            w = result["wheels"][wi]

            wheel_count = w["wheel_hits"]
            joint_count = w["joint_hits"]

            # Exact rational comparison:
            #
            # corr = joint / (f*wheel/interval)
            #      = joint*interval / (f*wheel)
            #
            # Printed as floating point only for display.
            if f_count == 0 or wheel_count == 0:
                corr = 0.0
            else:
                corr = (
                    joint_count * interval
                    / (f_count * wheel_count)
                )

            expected = (
                f_count * wheel_count / interval
                if interval
                else 0.0
            )

            print(
                f"{inst.index:2d} {u:3d} "
                f"{interval:10d} "
                f"{f_count:8d} "
                f"{wheel_count:10d} "
                f"{joint_count:10d} "
                f"{expected:12.4f} "
                f"{corr:12.6f}"
            )

    # ------------------------------------------------------------------------
    # Aggregate statistics
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("GLOBAL INDEPENDENCE SUMMARY")
    print("=" * 120)

    for wi, primes in enumerate(WHEEL_PREFIXES):

        total_interval = 0
        total_f = 0
        total_wheel = 0
        total_joint = 0

        for _, _, result in results:

            total_interval += result["interval"]
            total_f += result["f_count"]
            total_wheel += result["wheels"][wi]["wheel_hits"]
            total_joint += result["wheels"][wi]["joint_hits"]

        if total_interval == 0:
            continue

        expected = (
            total_f * total_wheel / total_interval
        )

        corr = (
            total_joint * total_interval
            / (total_f * total_wheel)
            if total_f and total_wheel
            else 0.0
        )

        wheel_density = total_wheel / total_interval
        f_density = total_f / total_interval
        joint_density = total_joint / total_interval

        print()
        print(f"  WHEEL = {primes}")
        print(f"    interval total     = {total_interval}")
        print(f"    f total            = {total_f}")
        print(f"    wheel total        = {total_wheel}")
        print(f"    joint total        = {total_joint}")
        print(f"    f density          = {f_density:.10f}")
        print(f"    wheel density      = {wheel_density:.10f}")
        print(f"    joint density      = {joint_density:.10f}")
        print(f"    independence pred  = {expected / total_interval:.10f}")
        print(f"    correlation ratio  = {corr:.10f}")

    # ------------------------------------------------------------------------
    # Scaling by u
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("CORRELATION BY UNKNOWN-K BITS")
    print("=" * 120)

    for wi, primes in enumerate(WHEEL_PREFIXES):

        print()
        print(f"WHEEL = {primes}")

        print(
            f"{'u':>4} {'tests':>6} "
            f"{'f-density':>14} {'wheel-density':>14} "
            f"{'joint-density':>14} {'corr':>12}"
        )

        for u in UNKNOWN_BITS:

            subset = [
                (inst, result)
                for inst, uu, result in results
                if uu == u
            ]

            if not subset:
                continue

            interval = sum(
                result["interval"]
                for _, result in subset
            )

            f_count = sum(
                result["f_count"]
                for _, result in subset
            )

            wheel_count = sum(
                result["wheels"][wi]["wheel_hits"]
                for _, result in subset
            )

            joint_count = sum(
                result["wheels"][wi]["joint_hits"]
                for _, result in subset
            )

            fd = f_count / interval
            wd = wheel_count / interval
            jd = joint_count / interval

            corr = (
                joint_count * interval
                / (f_count * wheel_count)
                if f_count and wheel_count
                else 0.0
            )

            print(
                f"{u:4d} {len(subset):6d} "
                f"{fd:14.10f} "
                f"{wd:14.10f} "
                f"{jd:14.10f} "
                f"{corr:12.6f}"
            )

    # ------------------------------------------------------------------------
    # Deviation ranking
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("LARGEST DEVIATIONS FROM INDEPENDENCE")
    print("=" * 120)

    ranked = []

    for inst, u, result in results:

        interval = result["interval"]
        f_count = result["f_count"]

        for wi, primes in enumerate(WHEEL_PREFIXES):

            wheel_count = result["wheels"][wi]["wheel_hits"]
            joint_count = result["wheels"][wi]["joint_hits"]

            if f_count == 0 or wheel_count == 0:
                continue

            corr = (
                joint_count * interval
                / (f_count * wheel_count)
            )

            deviation = abs(corr - 1.0)

            ranked.append(
                (
                    deviation,
                    corr,
                    inst.index,
                    u,
                    primes,
                    interval,
                    f_count,
                    wheel_count,
                    joint_count,
                )
            )

    ranked.sort(reverse=True)

    for (
        deviation,
        corr,
        index,
        u,
        primes,
        interval,
        f_count,
        wheel_count,
        joint_count,
    ) in ranked[:30]:

        print(
            f"  instance={index:2d} "
            f"u={u:2d} "
            f"wheel={primes!s:28s} "
            f"interval={interval:9d} "
            f"f={f_count:7d} "
            f"wheel={wheel_count:8d} "
            f"joint={joint_count:7d} "
            f"corr={corr:.8f} "
            f"|corr-1|={deviation:.8f}"
        )

    # ------------------------------------------------------------------------
    # Exact checks
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    for inst, u, result in results:

        assert result["exact_count"] == 1

        for w in result["wheels"]:

            assert w["joint_hits"] <= result["f_count"]
            assert w["joint_hits"] <= w["wheel_hits"]

            # The true d must satisfy the joint necessary conditions.
            assert f_condition(inst, inst.d, u)
            assert wheel_accepts(inst, inst.d, w["primes"])

            reconstructed = reconstruct(inst, inst.d)

            assert reconstructed == (
                min(inst.p, inst.q),
                max(inst.p, inst.q),
            )

    print(
        "  all exact f conditions preserve true d          = True"
    )
    print(
        "  all wheel conditions preserve true d            = True"
    )
    print(
        "  all joint conditions preserve true d            = True"
    )
    print(
        "  all exact factor reconstructions are valid      = True"
    )
    print(
        "  all calculations use exact integer predicates   = True"
    )

    # ------------------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 435 measured wheel density and f-candidate density"
    )
    print(
        "  separately."
    )
    print()
    print(
        "  Experiment 436 measures their joint distribution."
    )
    print()
    print(
        "  The central statistic is"
    )
    print(
        "      C = D_fw / (D_f * D_w)."
    )
    print()
    print(
        "  C approximately 1 means the wheel contributes roughly the"
    )
    print(
        "  expected rejection under independence."
    )
    print()
    print(
        "  Significant systematic deviation from 1 would show that the"
    )
    print(
        "  residue-wheel condition is arithmetically correlated with"
    )
    print(
        "  the exact 2^u f-condition."
    )
    print()
    print(
        "  That distinction matters because a wheel can appear to have"
    )
    print(
        "  a fixed density while nevertheless being unusually effective"
    )
    print(
        "  or ineffective on the actual f-candidate population."
    )
    print()
    print(
        "  This remains a candidate-isolation experiment, not a"
    )
    print(
        "  factorization theorem."
    )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXPERIMENT 436 FINAL STATUS")
    print("=" * 120)

    print(f"  ENUMERATED CASES = {len(results)}")
    print(
        f"  ALL TRUE-D SURVIVAL CHECKS = True"
    )
    print(
        f"  ALL EXACT RECONSTRUCTION CHECKS = True"
    )
    print("=" * 120)
    print("EXPERIMENT 436 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
