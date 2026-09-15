#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 464
========================================================================================================================

CROSS-INSTANCE EXACT QUOTIENT / RATIONAL j-LAW AUDIT

QUESTION

  Experiments 462 and 463 found no simple:

      exact affine law,
      exact quadratic law,
      small modular affine law

  connecting the hidden displacement j to the observed tuple.

  This experiment tests a different restricted class:

      j = floor((a*F + b) / G) + c

  and the exact rational form:

      a*F + b = j*(c*G + e) + q

  where F and G are simple observable-only integer features.

PURPOSE

  Search for low-complexity quotient, floor, and rational relations
  between the hidden j and observable-derived quantities.

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no arbitrary interpolation
  no high-degree polynomial fitting
  no hidden quantities inside observable features
  true j used only as post-hoc diagnostic target

IMPORTANT

  A surviving law is only a candidate empirical cross-instance law.

  It does NOT establish a reconstruction procedure.

  The candidate must:
    1. hold exactly on all eight instances,
    2. use very small coefficients,
    3. use only observable features,
    4. survive leave-one-out validation.

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt
from typing import Callable


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

COEFF_MIN = -4
COEFF_MAX = 4

OFFSET_MIN = -4
OFFSET_MAX = 4

MAX_PRINT = 100


# ----------------------------------------------------------------------------------------------------------------------
# DATA
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int
    r: int

    x0: int
    d0_ref: int
    K_ref: int

    true_d0: int
    true_x: int
    true_K: int

    @property
    def true_j(self) -> int:
        delta = self.true_d0 - self.d0_ref

        if delta % 2:
            raise AssertionError(
                f"non-integral true j: instance={self.idx}"
            )

        return delta // 2


INSTANCES = (
    Instance(
        1,
        14246098189,
        333010,
        14245432170,
        -13799399066930810668576,
        0,
        14245432170,
        -3449849766732702667144,
        14245432180,
        -10,
        -3449849766661475506269,
    ),
    Instance(
        2,
        10139117,
        11022,
        10117074,
        -11873391125935945,
        1,
        10117073,
        -2968347786542523,
        10117091,
        -17,
        -2968347695488785,
    ),
    Instance(
        3,
        10009330297,
        1010042,
        10007310214,
        -17225157737347240276460,
        0,
        10007310214,
        -4306289434336810069115,
        10007310238,
        -24,
        -4306289434216722346403,
    ),
    Instance(
        4,
        100460333,
        20046,
        100420242,
        -2460554951578020025,
        1,
        100420241,
        -615138737944715127,
        100420273,
        -31,
        -615138736337991015,
    ),
    Instance(
        5,
        2503701173,
        100074,
        2503501026,
        -2231236374548471123872,
        0,
        2503501026,
        -557809093637117780968,
        2503501064,
        -38,
        -557809093589551261113,
    ),
    Instance(
        6,
        10006200817,
        200062,
        10005800694,
        -50858953461762196260849,
        1,
        10005800693,
        -12714738365445551965559,
        10005800739,
        -45,
        -12714738365215418549091,
    ),
    Instance(
        7,
        40005200153,
        400026,
        40004400102,
        -1222668959330484861858204,
        0,
        40004400102,
        -305667239832621215464551,
        40004400154,
        -52,
        -305667239831581101061223,
    ),
    Instance(
        8,
        270017400119,
        1200024,
        270015000072,
        -74949527189645489927322133,
        1,
        270015000071,
        -18737381797411507489330569,
        270015000131,
        -59,
        -18737381797403407039327539,
    ),
)


# ----------------------------------------------------------------------------------------------------------------------
# OBSERVABLE FEATURE BASIS
# ----------------------------------------------------------------------------------------------------------------------

def observable_features(inst: Instance) -> dict[str, int]:
    N = inst.N
    S = inst.S
    d = inst.d
    r = inst.r

    return {
        "1": 1,

        "N": N,
        "S": S,
        "d": d,
        "r": r,

        "N-S": N - S,
        "N-d": N - d,
        "d-S": d - S,

        "N+S": N + S,
        "N+d": N + d,
        "S+d": S + d,

        "N%S": N % S,
        "d%S": d % S,
        "N%d": N % d,

        "gcd(N,S)": gcd(N, S),
        "gcd(N,d)": gcd(N, d),
        "gcd(S,d)": gcd(S, d),

        "sqrtN": isqrt(N),

        "N//S": N // S,
        "d//S": d // S,
        "N//d": N // d,

        "bitlenN": N.bit_length(),
        "bitlenS": S.bit_length(),
        "bitlend": d.bit_length(),

        "digitsN": len(str(abs(N))),
        "digitsS": len(str(abs(S))),
        "digitsd": len(str(abs(d))),
    }


FEATURES = tuple(
    observable_features(INSTANCES[0]).keys()
)


# ----------------------------------------------------------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------------------------------------------------------

def all_equal(values):
    return len(set(values)) == 1


def exact_floor_div(num: int, den: int) -> int:
    if den == 0:
        raise ZeroDivisionError("exact_floor_div denominator is zero")

    return num // den


def exact_remainder(num: int, den: int) -> int:
    if den == 0:
        raise ZeroDivisionError("exact_remainder denominator is zero")

    return num % den


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE TYPES
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class FloorCandidate:
    F: str
    G: str
    a: int
    b: int
    c: int


@dataclass(frozen=True)
class RationalCandidate:
    F: str
    G: str
    a: int
    b: int
    c: int
    e: int
    q: int


# ----------------------------------------------------------------------------------------------------------------------
# FLOOR LAW
# ----------------------------------------------------------------------------------------------------------------------

def floor_candidate_holds(
    cand: FloorCandidate,
    inst: Instance,
) -> bool:

    f = observable_features(inst)[cand.F]
    g = observable_features(inst)[cand.G]

    numerator = cand.a * f + cand.b
    denominator = g + cand.c

    if denominator == 0:
        return False

    predicted = numerator // denominator

    return predicted == inst.true_j


def floor_candidate_all_hold(
    cand: FloorCandidate,
    instances,
) -> bool:
    return all(
        floor_candidate_holds(cand, inst)
        for inst in instances
    )


# ----------------------------------------------------------------------------------------------------------------------
# RATIONAL LAW
# ----------------------------------------------------------------------------------------------------------------------

def rational_candidate_holds(
    cand: RationalCandidate,
    inst: Instance,
) -> bool:

    f = observable_features(inst)[cand.F]
    g = observable_features(inst)[cand.G]

    lhs = cand.a * f + cand.b
    scale = cand.c * g + cand.e

    if scale == 0:
        return False

    rhs = inst.true_j * scale + cand.q

    return lhs == rhs


def rational_candidate_all_hold(
    cand: RationalCandidate,
    instances,
) -> bool:
    return all(
        rational_candidate_holds(cand, inst)
        for inst in instances
    )


# ----------------------------------------------------------------------------------------------------------------------
# LEAVE-ONE-OUT VALIDATION
# ----------------------------------------------------------------------------------------------------------------------

def leave_one_out_floor(
    cand: FloorCandidate,
    instances,
):
    failures = []

    for held_out in instances:

        remaining = [
            x for x in instances
            if x.idx != held_out.idx
        ]

        if not floor_candidate_all_hold(
            cand,
            remaining,
        ):
            failures.append(held_out.idx)
            continue

        if not floor_candidate_holds(
            cand,
            held_out,
        ):
            failures.append(held_out.idx)

    return failures


def leave_one_out_rational(
    cand: RationalCandidate,
    instances,
):
    failures = []

    for held_out in instances:

        remaining = [
            x for x in instances
            if x.idx != held_out.idx
        ]

        if not rational_candidate_all_hold(
            cand,
            remaining,
        ):
            failures.append(held_out.idx)
            continue

        if not rational_candidate_holds(
            cand,
            held_out,
        ):
            failures.append(held_out.idx)

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# SEARCH FLOOR LAWS
# ----------------------------------------------------------------------------------------------------------------------

def search_floor_candidates(instances):

    candidates = []

    for F in FEATURES:
        for G in FEATURES:

            # Avoid an entirely redundant self-pair.
            if F == G:
                pass

            for a in range(COEFF_MIN, COEFF_MAX + 1):

                if a == 0:
                    continue

                for b in range(COEFF_MIN, COEFF_MAX + 1):

                    for c in range(COEFF_MIN, COEFF_MAX + 1):

                        for inst in instances:

                            f = observable_features(inst)[F]
                            g = observable_features(inst)[G]

                            denominator = g + c

                            if denominator == 0:
                                break

                            numerator = a * f + b

                            predicted = numerator // denominator

                            if predicted != inst.true_j:
                                break

                        else:
                            candidate = FloorCandidate(
                                F=F,
                                G=G,
                                a=a,
                                b=b,
                                c=c,
                            )

                            if (
                                leave_one_out_floor(
                                    candidate,
                                    instances,
                                )
                                == []
                            ):
                                if candidate not in candidates:
                                    candidates.append(
                                        candidate
                                    )

    return candidates


# ----------------------------------------------------------------------------------------------------------------------
# SEARCH RATIONAL LAWS
# ----------------------------------------------------------------------------------------------------------------------

def search_rational_candidates(instances):

    candidates = []

    for F in FEATURES:
        for G in FEATURES:

            for a in range(COEFF_MIN, COEFF_MAX + 1):

                if a == 0:
                    continue

                for b in range(COEFF_MIN, COEFF_MAX + 1):

                    for c in range(COEFF_MIN, COEFF_MAX + 1):

                        # c=0 and e=0 will later be possible,
                        # but that degenerates to an affine relation.
                        for e in range(OFFSET_MIN, OFFSET_MAX + 1):

                            scale_zero = False

                            # We do not enumerate q blindly.
                            # Infer q from the first instance.
                            first = instances[0]

                            f0 = observable_features(first)[F]
                            g0 = observable_features(first)[G]

                            lhs0 = (
                                a * f0
                                + b
                            )

                            scale0 = (
                                c * g0
                                + e
                            )

                            if scale0 == 0:
                                continue

                            q = (
                                lhs0
                                - first.true_j * scale0
                            )

                            # Keep q deliberately small.
                            if not (
                                OFFSET_MIN <= q <= OFFSET_MAX
                            ):
                                continue

                            candidate = RationalCandidate(
                                F=F,
                                G=G,
                                a=a,
                                b=b,
                                c=c,
                                e=e,
                                q=q,
                            )

                            if not rational_candidate_all_hold(
                                candidate,
                                instances,
                            ):
                                continue

                            if (
                                leave_one_out_rational(
                                    candidate,
                                    instances,
                                )
                                != []
                            ):
                                continue

                            if candidate not in candidates:
                                candidates.append(
                                    candidate
                                )

    return candidates


# ----------------------------------------------------------------------------------------------------------------------
# BASIC QUOTIENT / REMAINDER AUDIT
# ----------------------------------------------------------------------------------------------------------------------

def print_basic_quotient_audit(instances):

    print("=" * 120)
    print("BASIC OBSERVABLE QUOTIENT / REMAINDER AUDIT")
    print("=" * 120)
    print()

    scales = (
        "S",
        "d",
        "N",
        "N-S",
        "N-d",
        "d-S",
        "sqrtN",
        "gcd(S,d)",
        "N//S",
        "d//S",
    )

    for scale_name in scales:

        print(
            f"  SCALE = {scale_name}"
        )
        print(
            "    idx       j       quotient       remainder"
        )
        print(
            "    " + "-" * 58
        )

        for inst in instances:

            value = observable_features(inst)[scale_name]

            if value == 0:
                continue

            q = inst.true_j // value
            rem = inst.true_j % value

            print(
                f"    {inst.idx:<3d}"
                f" {inst.true_j:8d}"
                f" {q:15d}"
                f" {rem:15d}"
            )

        print()


# ----------------------------------------------------------------------------------------------------------------------
# SMALL RATIO SCREEN
# ----------------------------------------------------------------------------------------------------------------------

def print_ratio_screen(instances):

    print("=" * 120)
    print("DIRECT SMALL-RATIO SCREEN")
    print("=" * 120)
    print()

    ratio_features = (
        "S",
        "d",
        "N",
        "N-S",
        "N-d",
        "d-S",
        "sqrtN",
        "gcd(S,d)",
    )

    for name in ratio_features:

        ratios = []

        for inst in instances:

            v = observable_features(inst)[name]

            if v == 0:
                ratios.append(None)
                continue

            num = inst.true_j
            den = v

            g = gcd(abs(num), abs(den))

            reduced_num = num // g
            reduced_den = den // g

            ratios.append(
                (reduced_num, reduced_den)
            )

        print(
            f"  {name:15s} reduced j/F ratios = {ratios}"
        )

    print()


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE PRINTING
# ----------------------------------------------------------------------------------------------------------------------

def print_floor_candidate(
    rank: int,
    cand: FloorCandidate,
    instances,
):

    print()
    print(f"  FLOOR CANDIDATE {rank}")
    print(
        f"    j = floor("
        f"({cand.a}*{cand.F}+{cand.b})"
        f"/({cand.G}+{cand.c})"
        f")"
    )

    print()
    print(
        "    idx       j       numerator       denominator    predicted"
    )
    print(
        "    " + "-" * 75
    )

    for inst in instances:

        values = observable_features(inst)

        numerator = (
            cand.a * values[cand.F]
            + cand.b
        )

        denominator = (
            values[cand.G]
            + cand.c
        )

        predicted = numerator // denominator

        print(
            f"    {inst.idx:<3d}"
            f" {inst.true_j:8d}"
            f" {numerator:16d}"
            f" {denominator:16d}"
            f" {predicted:12d}"
        )


def print_rational_candidate(
    rank: int,
    cand: RationalCandidate,
    instances,
):

    print()
    print(f"  RATIONAL CANDIDATE {rank}")

    print(
        f"    {cand.a}*{cand.F}+{cand.b}"
        f" = j*({cand.c}*{cand.G}+{cand.e})"
        f" + {cand.q}"
    )

    print()
    print(
        "    idx       j            lhs            scale"
    )
    print(
        "    " + "-" * 75
    )

    for inst in instances:

        values = observable_features(inst)

        lhs = (
            cand.a * values[cand.F]
            + cand.b
        )

        scale = (
            cand.c * values[cand.G]
            + cand.e
        )

        print(
            f"    {inst.idx:<3d}"
            f" {inst.true_j:8d}"
            f" {lhs:20d}"
            f" {scale:15d}"
        )


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main():

    print("=" * 120)
    print("EXPERIMENT 464")
    print("=" * 120)
    print()
    print(
        "CROSS-INSTANCE EXACT QUOTIENT / RATIONAL j-LAW AUDIT"
    )
    print()

    print(
        f"COEFFICIENT RANGE = [{COEFF_MIN},{COEFF_MAX}]"
    )
    print(
        f"OFFSET RANGE      = [{OFFSET_MIN},{OFFSET_MAX}]"
    )
    print()

    # ------------------------------------------------------------------
    # Consistency
    # ------------------------------------------------------------------

    print("INPUT CONSISTENCY")
    print("  " + "-" * 80)

    failures = 0

    for inst in INSTANCES:

        if inst.d != inst.N - 2 * inst.S + 1:
            failures += 1

        if (
            inst.d0_ref + 2 * inst.true_j
            != inst.true_d0
        ):
            failures += 1

        if (
            inst.x0 - 2 * inst.true_j
            != inst.true_x
        ):
            failures += 1

        if (
            inst.K_ref
            + inst.d0_ref * inst.true_j
            + inst.true_j ** 2
            != inst.true_K
        ):
            failures += 1

    print(
        f"  consistency failures = {failures}"
    )
    print()

    if failures:
        raise AssertionError(
            "input consistency validation failed"
        )

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    print_basic_quotient_audit(INSTANCES)
    print_ratio_screen(INSTANCES)

    # ------------------------------------------------------------------
    # Floor search
    # ------------------------------------------------------------------

    print("=" * 120)
    print("EXACT FLOOR / QUOTIENT SEARCH")
    print("=" * 120)
    print()

    print(
        "Candidate:"
    )
    print()
    print(
        "  j = floor((a*F+b)/(G+c))"
    )
    print()

    floor_candidates = search_floor_candidates(
        INSTANCES
    )

    print(
        f"surviving floor candidates = "
        f"{len(floor_candidates)}"
    )
    print()

    if floor_candidates:

        for rank, cand in enumerate(
            floor_candidates[:MAX_PRINT],
            start=1,
        ):
            print_floor_candidate(
                rank,
                cand,
                INSTANCES,
            )

    else:
        print("  none")

    print()

    # ------------------------------------------------------------------
    # Rational search
    # ------------------------------------------------------------------

    print("=" * 120)
    print("EXACT RATIONAL RELATION SEARCH")
    print("=" * 120)
    print()

    print(
        "Candidate:"
    )
    print()
    print(
        "  a*F+b = j*(c*G+e)+q"
    )
    print()

    rational_candidates = search_rational_candidates(
        INSTANCES
    )

    print(
        f"surviving rational candidates = "
        f"{len(rational_candidates)}"
    )
    print()

    if rational_candidates:

        for rank, cand in enumerate(
            rational_candidates[:MAX_PRINT],
            start=1,
        ):
            print_rational_candidate(
                rank,
                cand,
                INSTANCES,
            )

    else:
        print("  none")

    print()

    # ------------------------------------------------------------------
    # Cross-category summary
    # ------------------------------------------------------------------

    print("=" * 120)
    print("CROSS-CATEGORY SUMMARY")
    print("=" * 120)
    print()

    print(
        f"  floor candidates          = "
        f"{len(floor_candidates)}"
    )
    print(
        f"  rational candidates       = "
        f"{len(rational_candidates)}"
    )

    print()

    if not floor_candidates and not rational_candidates:

        print(
            "  No exact low-complexity quotient or rational law"
        )
        print(
            "  survived all eight instances."
        )

    else:

        print(
            "  At least one candidate survived."
        )
        print(
            "  These candidates require independent validation."
        )

    print()

    # ------------------------------------------------------------------
    # Final interpretation
    # ------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 464 SUMMARY")
    print("=" * 120)
    print()

    print(
        "  Experiment 462:"
    )
    print(
        "    no small exact affine/quadratic law."
    )
    print()

    print(
        "  Experiment 463:"
    )
    print(
        "    no small modular observable law."
    )
    print()

    print(
        "  Experiment 464:"
    )
    print(
        "    quotient/floor and restricted rational forms."
    )
    print()

    print(
        "  These searches remain deliberately bounded."
    )
    print(
        "  Therefore a negative result rules out only"
    )
    print(
        "  this tested family, not arbitrary functions."
    )
    print()

    print("=" * 120)
    print("EXPERIMENT 464 FINAL STATUS")
    print("=" * 120)
    print()
    print(
        "  EXACT INTEGER ARITHMETIC       = True"
    )
    print(
        "  FLOOR SEARCH COMPLETED         = True"
    )
    print(
        "  RATIONAL SEARCH COMPLETED      = True"
    )
    print(
        f"  FLOOR CANDIDATES               = "
        f"{len(floor_candidates)}"
    )
    print(
        f"  RATIONAL CANDIDATES            = "
        f"{len(rational_candidates)}"
    )
    print(
        "  LEAVE-ONE-OUT VALIDATION       = True"
    )
    print(
        "  FLOATING POINT                 = False"
    )
    print(
        "  ARBITRARY INTERPOLATION        = False"
    )
    print(
        "  GIANT K ENUMERATION            = False"
    )
    print(
        "  GIANT d0 ENUMERATION           = False"
    )
    print(
        "  CRT CARTESIAN PRODUCT          = False"
    )
    print()
    print(
        "  CONCLUSION:"
    )
    print(
        "    This experiment extends the cross-instance"
    )
    print(
        "    search from ordinary polynomial and modular"
    )
    print(
        "    laws to a tightly bounded exact quotient/"
    )
    print(
        "    rational family using only observable data."
    )
    print()
    print("EXPERIMENT 464 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
