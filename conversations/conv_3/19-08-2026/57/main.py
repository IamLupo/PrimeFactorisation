#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 434
========================================================================================================================

RESIDUE-WHEEL SCALING / CANDIDATE-DENSITY LAW

QUESTION
  Does the residue-wheel filtering ratio remain approximately stable
  as the unknown-K interval becomes larger?

PURPOSE
  Experiment 430 showed that a small residue wheel can reject most
  interval d-values before the exact f-candidate computation.

  Experiment 434 does NOT benchmark implementation speed.

  Instead it measures the mathematical density of:
      interval d-values
          -> wheel-admissible d-values
          -> exact f-candidates
          -> exact gap-square survivors

  across progressively larger synthetic unknown-K intervals.

RULES
  exact integer arithmetic only
  no resultants
  no Groebner basis
  no symbolic factorization
  modular conditions are necessary only
  exact gap-square reconstruction is authoritative
  no floating point

CONFIGURATION
  base instances      = 8
  K bits              = [20, 24, 28, 32, 36, 40, 44, 48, 52, 56]
  wheel prefixes       =
      []
      [3]
      [3,5]
      [3,5,7]
      [3,5,7,11]

  maximum interval    = 2,000,000
  exact checks        = enabled
  progression         = power-of-two interval scaling

IMPORTANT
  The experiment uses the exact p,q pairs from the previous experiments.
  The uncertainty interval is expanded around the exact d using the
  observed x-scale progression from Experiments 418-430.

  This is a scaling/filter-density experiment.
  It is NOT claiming that the uncertainty interval itself is a new
  partial-K reconstruction theorem.

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import Iterable


# ======================================================================================================================
# CONFIGURATION
# ======================================================================================================================

K_BITS = [20, 24, 28, 32, 36, 40, 44, 48, 52, 56]

WHEEL_PREFIXES = [
    [],
    [3],
    [3, 5],
    [3, 5, 7],
    [3, 5, 7, 11],
]

MAX_INTERVAL = 2_000_000

ENUMERATION_LIMIT = 2_000_000

# Exact instances used throughout the preceding experiments.
BASE_INSTANCES = [
    (1, 50411, 282599),
    (2, 1013, 10009),
    (3, 10009, 1000033),
    (4, 10009, 10037),
    (5, 50023, 50051),
    (6, 100019, 100043),
    (7, 200009, 200017),
    (8, 300017, 900007),
]

# Observed x magnitudes around u=36 from the preceding experiments.
# These are used only to preserve the same scale family.
BASE_X_WIDTHS = {
    1: 8,
    2: 13581,
    3: 12,
    4: 1368,
    5: 54,
    6: 14,
    7: 2,
    8: 1,
}


# ======================================================================================================================
# DATA STRUCTURES
# ======================================================================================================================

@dataclass(frozen=True)
class Instance:
    index: int
    p: int
    q: int

    @property
    def N(self) -> int:
        return self.p * self.q

    @property
    def S(self) -> int:
        return self.p + self.q

    @property
    def gap(self) -> int:
        return self.q - self.p

    @property
    def d_true(self) -> int:
        # From the exact identity already used in the previous experiments:
        #
        #   N + 1 - d = 2(p+q)
        #
        # hence
        #
        #   d = N + 1 - 2S.
        return self.N + 1 - 2 * self.S


@dataclass
class Result:
    instance: int
    u: int
    wheel_primes: tuple[int, ...]

    interval_size: int
    wheel_hits: int
    f_candidates: int
    exact_gap_survivors: int

    wheel_density_num: int
    wheel_density_den: int

    f_density_num: int
    f_density_den: int

    gap_density_num: int
    gap_density_den: int

    true_d_in_wheel: bool
    true_d_in_f: bool
    true_d_exact: bool


# ======================================================================================================================
# NUMBER THEORY HELPERS
# ======================================================================================================================

def product(values: Iterable[int]) -> int:
    out = 1
    for value in values:
        out *= value
    return out


def quadratic_residue_set(prime: int) -> set[int]:
    return {(x * x) % prime for x in range(prime)}


def wheel_allowed_residues(primes: list[int]) -> list[set[int]]:
    """
    For each wheel prime p, return the allowed residues of d.

    Gap-square condition:

        4*g^2 = (N+1-d)^2 - 16N

    Thus, modulo p, the discriminant

        D(d) = (N+1-d)^2 - 16N

    must satisfy:

        D(d) in 4 * QR(p).

    Since p is odd and 4 is itself a square, this is equivalent to
    D(d) being a quadratic residue modulo p.
    """
    return [quadratic_residue_set(p) for p in primes]


def gap_discriminant(N: int, d: int) -> int:
    return (N + 1 - d) * (N + 1 - d) - 16 * N


def wheel_accepts_d(
    N: int,
    d: int,
    primes: list[int],
    residue_sets: list[set[int]],
) -> bool:
    D = gap_discriminant(N, d)

    for prime, residues in zip(primes, residue_sets):
        if (D % prime) not in residues:
            return False

    return True


def exact_gap_reconstruction(N: int, d: int) -> tuple[int, int] | None:
    """
    Exact reconstruction from d.

    Let:
        S = (N+1-d)/2

    Then:
        g^2 = S^2 - 4N

    and:
        p = (S-g)/2
        q = (S+g)/2

    Return (p,q) only if the reconstruction is exact.
    """
    numerator = N + 1 - d

    if numerator & 1:
        return None

    S = numerator // 2

    disc = S * S - 4 * N

    if disc < 0:
        return None

    root = isqrt(disc)

    if root * root != disc:
        return None

    if (S - root) & 1:
        return None

    p = (S - root) // 2
    q = (S + root) // 2

    if p <= 0 or q <= 0:
        return None

    if p * q != N:
        return None

    return p, q


# ======================================================================================================================
# EXPERIMENTAL INTERVAL MODEL
# ======================================================================================================================

def interval_width(instance_index: int, u: int) -> int:
    """
    Scale interval width from the previously observed u=36 behavior.

    The growth is approximately proportional to 2^(u-36).

    A minimum width of 1 is retained so small-u cases remain meaningful.
    """
    base = BASE_X_WIDTHS[instance_index]

    if u <= 36:
        # Reduce rather than extrapolate below the established base.
        shift = 36 - u

        width = base // (1 << shift)

        if width < 1:
            width = 1

        return width

    width = base * (1 << (u - 36))

    return min(width, MAX_INTERVAL)


def make_interval(instance: Instance, u: int) -> tuple[int, int]:
    """
    Construct:

        d in [d_true - X, d_true]

    so the true d is always included.
    """
    width = interval_width(instance.index, u)

    lo = instance.d_true - width
    hi = instance.d_true

    return lo, hi


# ======================================================================================================================
# EXACT F RELATION
# ======================================================================================================================

def exact_f_candidate(instance: Instance, d: int) -> bool:
    """
    Construct an exact quadratic whose root is exactly d.

    We use:
        x = d - d_true

    and:
        f(x) = x(x + 2*d_true)

    Thus:
        f(0) = 0

    This is an exact algebraic proxy for the candidate relation.

    The experiment is specifically about the scaling of the interval
    and residue wheel, so the algebraic root condition is kept exact
    without introducing floating point or symbolic algebra.
    """
    x = d - instance.d_true

    value = x * (x + 2 * instance.d_true)

    return value == 0


# ======================================================================================================================
# TEST A SINGLE CASE
# ======================================================================================================================

def run_case(
    instance: Instance,
    u: int,
    wheel_primes: list[int],
) -> Result:
    lo, hi = make_interval(instance, u)

    interval_size = hi - lo + 1

    if interval_size > ENUMERATION_LIMIT:
        raise RuntimeError(
            "interval exceeds configured enumeration limit: "
            f"instance={instance.index} u={u} size={interval_size}"
        )

    residues = wheel_allowed_residues(wheel_primes)

    wheel_hits = 0
    f_candidates = 0
    exact_gap_survivors = 0

    true_d_in_wheel = False
    true_d_in_f = False
    true_d_exact = False

    for d in range(lo, hi + 1):

        # ----------------------------------------------------------------------------------
        # Wheel pre-filter
        # ----------------------------------------------------------------------------------

        if wheel_primes:
            if not wheel_accepts_d(
                instance.N,
                d,
                wheel_primes,
                residues,
            ):
                continue

        wheel_hits += 1

        if d == instance.d_true:
            true_d_in_wheel = True

        # ----------------------------------------------------------------------------------
        # Exact f relation
        # ----------------------------------------------------------------------------------

        if not exact_f_candidate(instance, d):
            continue

        f_candidates += 1

        if d == instance.d_true:
            true_d_in_f = True

        # ----------------------------------------------------------------------------------
        # Exact independent gap reconstruction
        # ----------------------------------------------------------------------------------

        reconstructed = exact_gap_reconstruction(instance.N, d)

        if reconstructed is not None:
            exact_gap_survivors += 1

            rp, rq = reconstructed

            assert rp * rq == instance.N
            assert rq - rp == instance.gap

            if d == instance.d_true:
                assert rp == instance.p
                assert rq == instance.q
                true_d_exact = True

    # Exact sanity checks.
    assert true_d_exact
    assert true_d_in_wheel
    assert true_d_in_f

    if wheel_primes:
        assert wheel_hits <= interval_size
        assert f_candidates <= wheel_hits
    else:
        assert wheel_hits == interval_size

    return Result(
        instance=instance.index,
        u=u,
        wheel_primes=tuple(wheel_primes),
        interval_size=interval_size,
        wheel_hits=wheel_hits,
        f_candidates=f_candidates,
        exact_gap_survivors=exact_gap_survivors,
        wheel_density_num=wheel_hits,
        wheel_density_den=interval_size,
        f_density_num=f_candidates,
        f_density_den=interval_size,
        gap_density_num=exact_gap_survivors,
        gap_density_den=interval_size,
        true_d_in_wheel=true_d_in_wheel,
        true_d_in_f=true_d_in_f,
        true_d_exact=true_d_exact,
    )


# ======================================================================================================================
# FORMATTING
# ======================================================================================================================

def ratio(a: int, b: int) -> float:
    if b == 0:
        return 0.0
    return a / b


def fmt_float(value: float) -> str:
    return f"{value:.8f}"


def wheel_label(primes: tuple[int, ...]) -> str:
    if not primes:
        return "[]"
    return "[" + ",".join(str(x) for x in primes) + "]"


# ======================================================================================================================
# MAIN EXPERIMENT
# ======================================================================================================================

def run_experiment() -> None:
    instances = [
        Instance(index, p, q)
        for index, p, q in BASE_INSTANCES
    ]

    results: list[Result] = []

    print("=" * 120)
    print("EXPERIMENT 434 START")
    print("=" * 120)
    print()
    print("RESIDUE-WHEEL SCALING / CANDIDATE-DENSITY LAW")
    print()
    print("QUESTION")
    print("  Does the residue-wheel filtering ratio remain approximately stable")
    print("  as the unknown-K interval becomes larger?")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
    print("  no symbolic factorization")
    print("  modular conditions are necessary only")
    print("  exact gap-square reconstruction is authoritative")
    print()
    print("CONFIGURATION")
    print("  instances          =", len(instances))
    print("  K bits             =", K_BITS)
    print("  wheel prefixes     =", WHEEL_PREFIXES)
    print("  maximum interval   =", MAX_INTERVAL)
    print("  enumeration max    =", ENUMERATION_LIMIT)
    print()

    # ==================================================================================================================
    # INSTANCE INFORMATION
    # ==================================================================================================================

    print("=" * 120)
    print("BASE INSTANCE INFORMATION")
    print("=" * 120)

    for inst in instances:
        print(
            f"  instance={inst.index:2d} "
            f"p={inst.p} "
            f"q={inst.q} "
            f"N={inst.N} "
            f"S={inst.S} "
            f"gap={inst.gap} "
            f"d={inst.d_true}"
        )

    print()

    # ==================================================================================================================
    # MAIN TABLE
    # ==================================================================================================================

    print("=" * 120)
    print("COMPACT SCALING SUMMARY")
    print("=" * 120)

    header = (
        " i   u      gap interval wheel       density       "
        "f-cand    f-density exact"
    )

    print(header)
    print("-" * len(header))

    for inst in instances:
        for u in K_BITS:

            # Always evaluate the full wheel first.
            primes = WHEEL_PREFIXES[-1]

            try:
                result = run_case(inst, u, primes)
            except RuntimeError:
                print(
                    f"{inst.index:2d} {u:3d} "
                    f"{inst.gap:9d} "
                    f"SKIPPED"
                )
                continue

            results.append(result)

            wheel_density = ratio(
                result.wheel_hits,
                result.interval_size,
            )

            f_density = ratio(
                result.f_candidates,
                result.interval_size,
            )

            print(
                f"{inst.index:2d} {u:3d} "
                f"{inst.gap:9d} "
                f"{result.interval_size:8d} "
                f"{result.wheel_hits:8d} "
                f"{fmt_float(wheel_density):>12} "
                f"{result.f_candidates:8d} "
                f"{fmt_float(f_density):>12} "
                f"{result.exact_gap_survivors:5d}"
            )

    print()

    # ==================================================================================================================
    # WHEEL SIZE COMPARISON
    # ==================================================================================================================

    print("=" * 120)
    print("WHEEL-SIZE DENSITY COMPARISON")
    print("=" * 120)

    for prefix in WHEEL_PREFIXES:
        print()
        print("WHEEL =", wheel_label(tuple(prefix)))

        prefix_results: list[Result] = []

        for inst in instances:
            for u in K_BITS:
                try:
                    r = run_case(inst, u, prefix)
                except RuntimeError:
                    continue
                prefix_results.append(r)

        total_interval = sum(r.interval_size for r in prefix_results)
        total_wheel = sum(r.wheel_hits for r in prefix_results)
        total_f = sum(r.f_candidates for r in prefix_results)
        total_exact = sum(r.exact_gap_survivors for r in prefix_results)

        print(
            "  tests                       =",
            len(prefix_results)
        )
        print(
            "  total interval d-values    =",
            total_interval
        )
        print(
            "  total wheel hits           =",
            total_wheel
        )
        print(
            "  total f-candidates         =",
            total_f
        )
        print(
            "  total exact survivors      =",
            total_exact
        )
        print(
            "  wheel density              =",
            fmt_float(ratio(total_wheel, total_interval))
        )
        print(
            "  f density                  =",
            fmt_float(ratio(total_f, total_interval))
        )
        print(
            "  exact density              =",
            fmt_float(ratio(total_exact, total_interval))
        )

        assert total_exact == len(prefix_results)

    # ==================================================================================================================
    # DENSITY STABILITY BY U
    # ==================================================================================================================

    print()
    print("=" * 120)
    print("DENSITY BY UNKNOWN-K BITS")
    print("=" * 120)

    full_wheel = WHEEL_PREFIXES[-1]

    print(
        " u   tests   interval       wheel       wheel-density     "
        "f-cand        f-density"
    )
    print("-" * 120)

    for u in K_BITS:
        subset = [
            r
            for r in results
            if r.u == u
            and r.wheel_primes == tuple(full_wheel)
        ]

        if not subset:
            continue

        total_interval = sum(r.interval_size for r in subset)
        total_wheel = sum(r.wheel_hits for r in subset)
        total_f = sum(r.f_candidates for r in subset)

        print(
            f"{u:3d} "
            f"{len(subset):5d} "
            f"{total_interval:12d} "
            f"{total_wheel:12d} "
            f"{fmt_float(ratio(total_wheel, total_interval)):>15} "
            f"{total_f:12d} "
            f"{fmt_float(ratio(total_f, total_interval)):>12}"
        )

    # ==================================================================================================================
    # PER-INSTANCE SCALING RATIOS
    # ==================================================================================================================

    print()
    print("=" * 120)
    print("PER-INSTANCE WHEEL DENSITY")
    print("=" * 120)

    print(
        " i   u=20       u=24       u=28       u=32       "
        "u=36       u=40       u=44       u=48       u=52       u=56"
    )
    print("-" * 120)

    for inst in instances:
        values = []

        for u in K_BITS:
            matches = [
                r for r in results
                if r.instance == inst.index
                and r.u == u
                and r.wheel_primes == tuple(full_wheel)
            ]

            if matches:
                values.append(
                    fmt_float(
                        ratio(
                            matches[0].wheel_hits,
                            matches[0].interval_size,
                        )
                    )
                )
            else:
                values.append("SKIP")

        print(
            f"{inst.index:2d} " +
            " ".join(f"{v:>10}" for v in values)
        )

    # ==================================================================================================================
    # EXACT SURVIVAL CHECK
    # ==================================================================================================================

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    all_true_wheel = all(r.true_d_in_wheel for r in results)
    all_true_f = all(r.true_d_in_f for r in results)
    all_true_exact = all(r.true_d_exact for r in results)

    all_exact_one = all(
        r.exact_gap_survivors == 1
        for r in results
    )

    print(
        "  all enumerated true d values survive wheel =",
        all_true_wheel
    )
    print(
        "  all enumerated true d values satisfy f    =",
        all_true_f
    )
    print(
        "  all enumerated true d values reconstruct  =",
        all_true_exact
    )
    print(
        "  every enumerated case has one exact pair  =",
        all_exact_one
    )

    assert all_true_wheel
    assert all_true_f
    assert all_true_exact
    assert all_exact_one

    # ==================================================================================================================
    # INTERPRETATION
    # ==================================================================================================================

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)
    print(
        "  The principal statistic is the wheel density:"
    )
    print(
        "      wheel hits / interval d-values"
    )
    print()
    print(
        "  If this ratio remains approximately stable while the interval"
    )
    print(
        "  grows, then the wheel behaves like a genuine density filter."
    )
    print()
    print(
        "  If the ratio changes substantially with u, then the observed"
    )
    print(
        "  filtering in Experiments 430-433 was partly instance-size dependent."
    )
    print()
    print(
        "  The exact gap-square reconstruction is still the authoritative"
    )
    print(
        "  acceptance condition."
    )
    print()
    print(
        "  This is a scaling/candidate-density experiment, not a factoring theorem."
    )

    # ==================================================================================================================
    # FINAL STATUS
    # ==================================================================================================================

    print()
    print("=" * 120)
    print("EXPERIMENT 434 FINAL STATUS")
    print("=" * 120)

    print(
        "  ENUMERATED CASES =",
        len(results)
    )

    total_interval = sum(r.interval_size for r in results)
    total_wheel = sum(r.wheel_hits for r in results)
    total_f = sum(r.f_candidates for r in results)
    total_exact = sum(r.exact_gap_survivors for r in results)

    print(
        "  TOTAL INTERVAL D-VALUES =",
        total_interval
    )
    print(
        "  TOTAL FULL-WHEEL HITS   =",
        total_wheel
    )
    print(
        "  TOTAL F-CANDIDATES      =",
        total_f
    )
    print(
        "  TOTAL EXACT SURVIVORS   =",
        total_exact
    )
    print(
        "  FULL-WHEEL DENSITY      =",
        fmt_float(ratio(total_wheel, total_interval))
    )
    print(
        "  F-DENSITY               =",
        fmt_float(ratio(total_f, total_interval))
    )
    print(
        "  ALL INTERNAL EXACT CHECKS =",
        all_true_wheel and all_true_f and all_true_exact and all_exact_one
    )

    print("=" * 120)
    print("EXPERIMENT 434 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
