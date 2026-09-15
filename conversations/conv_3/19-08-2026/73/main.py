#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ======================================================================================================================
# EXPERIMENT 449
# ======================================================================================================================
#
# INDEPENDENT K-BOUND -> EXACT d0 POPULATION AUDIT
#
# QUESTION
#
#   Can an independently justified finite K interval reduce the
#   d0 family to a small exact population?
#
# CORE EQUATION
#
#       4K = r + d0^2 - d^2
#
# hence
#
#       d0^2 = d^2 + 4K - r
#
# For
#
#       K_low <= K <= K_high
#
# we obtain
#
#       d^2 + 4*K_low - r <= d0^2 <= d^2 + 4*K_high - r.
#
# The number of integer d0 values is therefore obtained from
# TWO exact integer square roots.
#
# The reconstruction logic NEVER uses:
#   - d0_true
#   - x_true
#   - K_true
#
# Those values occur only in post-hoc validation.
#
# NO GIANT K ENUMERATION
# NO GIANT d0 ENUMERATION
# NO CRT CARTESIAN PRODUCT
# INTEGER-EXACT
#
# ======================================================================================================================


@dataclass(frozen=True)
class Instance:
    idx: int

    # Construction-only quantities.
    p: int
    q: int

    # Public / reconstruction quantities.
    N: int
    S: int
    d: int
    r: int

    # Hidden values: ONLY for post-hoc checks.
    d0_true: int
    x_true: int
    K_true: int


# ----------------------------------------------------------------------------------------------------------------------
# Same controlled instances as experiments 447/448.
# ----------------------------------------------------------------------------------------------------------------------

RAW = [
    (50411, 282599, -10, -3449849766661475506269),
    (1013, 10009, -17, -2968347695488785),
    (10009, 1000033, -24, -4306289434216722346403),
    (10009, 10037, -31, -615138736337991015),
    (50023, 50051, -38, -557809093589551261113),
    (100019, 100043, -45, -12714738365215418549091),
    (200009, 200017, -52, -305667239831581101061223),
    (300017, 900007, -59, -18737381797403407039327539),
]


def build_instances() -> list[Instance]:
    result: list[Instance] = []

    for idx, (p, q, x_true, K_true) in enumerate(RAW, start=1):

        N = p * q
        S = p + q
        d = N - 2 * S + 1
        d0_true = d - x_true

        # Generated only from the hidden construction.
        # Reconstruction never uses d0_true.
        r = 4 * K_true + d * d - d0_true * d0_true

        result.append(
            Instance(
                idx=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                r=r,
                d0_true=d0_true,
                x_true=x_true,
                K_true=K_true,
            )
        )

    return result


# ======================================================================================================================
# BASIC EXACT HELPERS
# ======================================================================================================================

def candidate_k(inst: Instance, d0: int) -> tuple[bool, int]:
    """
    Exact evaluation of

        K(d0) = (r + d0^2 - d^2) / 4.
    """

    numerator = inst.r + d0 * d0 - inst.d * inst.d

    if numerator % 4 != 0:
        return False, 0

    return True, numerator // 4


def required_d0_parity(inst: Instance) -> int:
    """
    Determine which parity of d0 makes K integral.

    Since only mod 4 matters, testing d0=0 and d0=1 is sufficient.
    """

    for parity in (0, 1):

        numerator = (
            inst.r
            + parity * parity
            - inst.d * inst.d
        )

        if numerator % 4 == 0:
            return parity

    raise AssertionError(
        f"no integral d0 parity exists for instance={inst.idx}"
    )


def ceil_sqrt(n: int) -> int:
    """
    Exact ceil(sqrt(n)) for arbitrary Python integers.

    For n <= 0 the square condition requires special handling,
    so callers should not use this as a direct population count
    without checking the domain.
    """

    if n <= 0:
        return 0

    r = isqrt(n)

    if r * r == n:
        return r

    return r + 1


# ======================================================================================================================
# EXACT d0 POPULATION FROM A K INTERVAL
# ======================================================================================================================

@dataclass(frozen=True)
class PopulationResult:
    count: int
    d0_low: int | None
    d0_high: int | None
    parity: int
    square_low: int
    square_high: int
    true_in: bool


def exact_d0_population(
    inst: Instance,
    K_low: int,
    K_high: int,
) -> PopulationResult:

    if K_low > K_high:
        raise ValueError("K_low must not exceed K_high")

    parity = required_d0_parity(inst)

    # Exact pullback:
    #
    #   d0^2 = d^2 + 4K - r
    #
    square_low = inst.d * inst.d + 4 * K_low - inst.r
    square_high = inst.d * inst.d + 4 * K_high - inst.r

    if square_high < 0:
        return PopulationResult(
            count=0,
            d0_low=None,
            d0_high=None,
            parity=parity,
            square_low=square_low,
            square_high=square_high,
            true_in=False,
        )

    # No nonnegative d0^2 can be below zero.
    effective_low = max(square_low, 0)

    root_low = ceil_sqrt(effective_low)
    root_high = isqrt(square_high)

    if root_low > root_high:
        return PopulationResult(
            count=0,
            d0_low=None,
            d0_high=None,
            parity=parity,
            square_low=square_low,
            square_high=square_high,
            true_in=False,
        )

    # Positive branch first.
    first = root_low

    if (first & 1) != parity:
        first += 1

    last = root_high

    if (last & 1) != parity:
        last -= 1

    if first > last:
        return PopulationResult(
            count=0,
            d0_low=None,
            d0_high=None,
            parity=parity,
            square_low=square_low,
            square_high=square_high,
            true_in=False,
        )

    positive_count = ((last - first) // 2) + 1

    # Negative branch is symmetric.
    total_count = 2 * positive_count

    # d0 = 0 must not be double-counted.
    if parity == 0 and first == 0:
        total_count -= 1

    # True value is a post-hoc check only.
    true_in = False

    if square_low <= inst.d0_true * inst.d0_true <= square_high:
        true_in = ((inst.d0_true & 1) == parity)

    return PopulationResult(
        count=total_count,
        d0_low=-last,
        d0_high=last,
        parity=parity,
        square_low=square_low,
        square_high=square_high,
        true_in=true_in,
    )


# ======================================================================================================================
# BOUND MODELS
# ======================================================================================================================
#
# Every model below depends only on N,S,d,r.
#
# K_TRUE is NEVER used here.
#
# The bounds are diagnostic hypotheses, not claims of validity.
# ======================================================================================================================


@dataclass(frozen=True)
class BoundModel:
    name: str

    def interval(self, inst: Instance) -> tuple[int, int]:
        raise NotImplementedError


@dataclass(frozen=True)
class SymmetricBits(BoundModel):
    bits: int

    def interval(self, inst: Instance) -> tuple[int, int]:
        R = 1 << self.bits
        return -R, R


@dataclass(frozen=True)
class SymmetricN(BoundModel):

    def interval(self, inst: Instance) -> tuple[int, int]:
        R = inst.N
        return -R, R


@dataclass(frozen=True)
class SymmetricS(BoundModel):

    def interval(self, inst: Instance) -> tuple[int, int]:
        R = inst.S
        return -R, R


@dataclass(frozen=True)
class SymmetricNS(BoundModel):

    def interval(self, inst: Instance) -> tuple[int, int]:
        R = inst.N * inst.S
        return -R, R


@dataclass(frozen=True)
class SymmetricN2(BoundModel):

    def interval(self, inst: Instance) -> tuple[int, int]:
        R = inst.N * inst.N
        return -R, R


@dataclass(frozen=True)
class SymmetricS2(BoundModel):

    def interval(self, inst: Instance) -> tuple[int, int]:
        R = inst.S * inst.S
        return -R, R


@dataclass(frozen=True)
class IntervalDataShift(BoundModel):
    """
    K interval derived from zero-centered arithmetic scale:

        K ∈ [-(scale), +(scale)]

    where scale is a known-data expression.
    """

    factor: int
    expression: str

    def interval(self, inst: Instance) -> tuple[int, int]:

        if self.expression == "N":
            R = self.factor * inst.N
        elif self.expression == "S":
            R = self.factor * inst.S
        elif self.expression == "NS":
            R = self.factor * inst.N * inst.S
        elif self.expression == "N2":
            R = self.factor * inst.N * inst.N
        elif self.expression == "S2":
            R = self.factor * inst.S * inst.S
        else:
            raise ValueError(
                f"unknown expression={self.expression}"
            )

        return -R, R


def make_models() -> list[BoundModel]:

    models: list[BoundModel] = []

    for bits in (
        32,
        40,
        48,
        56,
        64,
        68,
        72,
        76,
        80,
        84,
        88,
        96,
    ):
        models.append(
            SymmetricBits(
                name=f"+/-2^{bits}",
                bits=bits,
            )
        )

    models.extend(
        [
            SymmetricS(name="+/-S"),
            SymmetricN(name="+/-N"),
            SymmetricNS(name="+/-N*S"),
            SymmetricS2(name="+/-S^2"),
            SymmetricN2(name="+/-N^2"),
            IntervalDataShift(
                name="+/-2S",
                factor=2,
                expression="S",
            ),
            IntervalDataShift(
                name="+/-4S",
                factor=4,
                expression="S",
            ),
            IntervalDataShift(
                name="+/-2N",
                factor=2,
                expression="N",
            ),
            IntervalDataShift(
                name="+/-4N",
                factor=4,
                expression="N",
            ),
        ]
    )

    return models


# ======================================================================================================================
# OPTIONAL ONE-SIDED K MODELS
# ======================================================================================================================

@dataclass(frozen=True)
class NegativeScale(BoundModel):
    """
    K in [-R, 0].
    """

    expression: str

    def interval(self, inst: Instance) -> tuple[int, int]:

        if self.expression == "S":
            R = inst.S
        elif self.expression == "N":
            R = inst.N
        elif self.expression == "NS":
            R = inst.N * inst.S
        elif self.expression == "N2":
            R = inst.N * inst.N
        else:
            raise ValueError(self.expression)

        return -R, 0


def make_one_sided_models() -> list[BoundModel]:

    return [
        NegativeScale(
            name="[-S,0]",
            expression="S",
        ),
        NegativeScale(
            name="[-N,0]",
            expression="N",
        ),
        NegativeScale(
            name="[-N*S,0]",
            expression="NS",
        ),
        NegativeScale(
            name="[-N^2,0]",
            expression="N2",
        ),
    ]


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    instances = build_instances()

    models = make_models() + make_one_sided_models()

    total_cases = 0
    true_contained = 0
    empty_pullbacks = 0
    unique_d0_cases = 0
    small_d0_cases = 0

    print("=" * 120)
    print("EXPERIMENT 449")
    print("=" * 120)
    print()
    print("INDEPENDENT K-BOUND -> EXACT d0 POPULATION AUDIT")
    print()
    print("QUESTION")
    print("  Can an independently specified finite K interval reduce")
    print("  the d0 family to a small exact population?")
    print()
    print("CORE IDENTITY")
    print("  4K = r + d0^2 - d^2")
    print()
    print("EXACT PULLBACK")
    print("  d0^2 = d^2 + 4K - r")
    print()
    print("RULES")
    print("  K_true is NEVER used to construct a bound")
    print("  d0_true is NEVER used to construct a bound")
    print("  x_true is NEVER used to construct a bound")
    print("  no giant K enumeration")
    print("  no giant d0 enumeration")
    print("  no CRT Cartesian product")
    print("  integer-exact")
    print()

    for inst in instances:

        print("-" * 120)
        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} "
            f"S={inst.S} "
            f"d={inst.d} "
            f"r={inst.r}"
        )

        print()
        print("  POST-HOC TRUE VALUES")
        print(f"    true d0 = {inst.d0_true}")
        print(f"    true x  = {inst.x_true}")
        print(f"    true K  = {inst.K_true}")

        parity = required_d0_parity(inst)

        print()
        print("  REQUIRED d0 PARITY")
        print(f"    d0 mod 2 = {parity}")

        print()
        print(
            "  MODEL                        K-low"
            "                         K-high"
            "       d0-count"
            "    true-in"
            "    uniqueness"
        )
        print("  " + "-" * 112)

        instance_cases = 0

        for model in models:

            K_low, K_high = model.interval(inst)

            pop = exact_d0_population(
                inst,
                K_low,
                K_high,
            )

            total_cases += 1
            instance_cases += 1

            if pop.true_in:
                true_contained += 1

            if pop.count == 0:
                empty_pullbacks += 1

            if pop.count == 1:
                unique_d0_cases += 1

            if 1 < pop.count <= 10:
                small_d0_cases += 1

            if pop.count == 0:
                uniqueness = "EMPTY"
            elif pop.count == 1:
                uniqueness = "UNIQUE"
            elif pop.count <= 10:
                uniqueness = "SMALL"
            else:
                uniqueness = "LARGE"

            print(
                f"  {model.name:<20}"
                f"{K_low:>26}"
                f"{K_high:>26}"
                f"{pop.count:>16}"
                f"{str(pop.true_in):>12}"
                f"{uniqueness:>14}"
            )

        print()
        print(
            f"  instance interval tests = {instance_cases}"
        )

    # ==============================================================================================================
    # GLOBAL SUMMARY
    # ==============================================================================================================

    print()
    print("=" * 120)
    print("GLOBAL EXPERIMENT 449 SUMMARY")
    print("=" * 120)

    print(
        f"  independent K interval tests = {total_cases}"
    )

    print(
        f"  true-containing intervals    = {true_contained}"
        f"/{total_cases}"
    )

    print(
        f"  empty d0 pullbacks           = {empty_pullbacks}"
    )

    print(
        f"  unique d0 cases              = {unique_d0_cases}"
    )

    print(
        f"  small d0 cases (2..10)       = {small_d0_cases}"
    )

    print()
    print("EXACT PULLBACK FORMULA")
    print()
    print("  K_low <= K <= K_high")
    print()
    print("      <=>")
    print()
    print("  d^2 + 4K_low - r")
    print("       <= d0^2 <=")
    print("  d^2 + 4K_high - r")
    print()
    print("  Therefore the exact candidate count is obtained from")
    print("  integer square roots plus the required parity.")
    print()
    print("  No candidate K values are enumerated.")
    print("  No candidate d0 values are enumerated.")

    print()
    print("INTERPRETATION")
    print()
    print("  Experiment 448 showed that the bare quadratic relation")
    print("  leaves a parity-compatible d0 family.")
    print()
    print("  Experiment 449 asks whether an EXTERNAL K bound can")
    print("  collapse that family.")
    print()
    print("  A useful result is not merely a small population.")
    print("  The K interval itself must be independently justified.")
    print()
    print("  A diagnostic interval such as")
    print("      |K| < 2^72")
    print("  measures computational narrowing but is not a proof")
    print("  unless that bound follows independently from the")
    print("  available public data.")
    print()
    print("  The same applies to scales such as")
    print("      |K| <= N")
    print("      |K| <= N*S")
    print("      |K| <= N^2.")
    print()

    print("IMPORTANT")
    print()
    print("  If one of these genuinely independent K bounds produces")
    print("  a tiny d0 population, that is the next object to inspect.")
    print()
    print("  If every independently justified bound still produces")
    print("  a large d0 population, then the obstruction is an")
    print("  information bound rather than an enumeration problem.")

    print()
    print("=" * 120)
    print("EXPERIMENT 449 FINAL STATUS")
    print("=" * 120)

    print(
        "  O(1) K->d0 inversion       = True"
    )

    print(
        "  GIANT K ENUMERATION        = False"
    )

    print(
        "  GIANT d0 ENUMERATION       = False"
    )

    print(
        "  CRT CARTESIAN PRODUCT      = False"
    )

    print(
        "  K_TRUE USED IN BOUNDS      = False"
    )

    print(
        "  d0_TRUE USED IN BOUNDS     = False"
    )

    print(
        "  X_TRUE USED IN BOUNDS      = False"
    )

    print(
        "  INTEGER-EXACT              = True"
    )

    print()
    print(
        "  CONCLUSION:"
    )
    print(
        "    This experiment isolates whether an independently"
    )
    print(
        "    justified K bound can actually identify d0."
    )

    print("=" * 120)
    print("EXPERIMENT 449 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()

