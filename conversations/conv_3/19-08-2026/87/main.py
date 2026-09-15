#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 463
========================================================================================================================

CROSS-INSTANCE MODULAR j-LAW / OBSERVED-DATA AUDIT

QUESTION
  Experiment 462 found no small exact affine or quadratic law
  of the form

      a*j + F(N,S,d,r) = 0.

  Could j nevertheless obey a simple modular relation determined
  by the observed tuple?

  We therefore test relations of the form

      j == a*F(N,S,d,r) + b (mod m)

  for small moduli m and low-complexity observable features F.

PURPOSE
  Detect simple modular structure connecting the hidden displacement
  j to the observed tuple across instances.

RULES
  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  true j is used only as post-hoc diagnostic data
  true j is never used to construct a bound
  no arbitrary interpolation
  no fitted high-degree polynomial
  no large coefficients

IMPORTANT
  A modular relation surviving all eight instances is only a
  candidate empirical law.

  It is NOT treated as proof of a reconstruction mechanism.

  Because only eight instances are available, the experiment
  emphasizes very small moduli, simple features, and
  leave-one-out validation.

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt
from typing import Callable


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

MODULUS_MIN = 2
MODULUS_MAX = 256

# Only very small coefficients are considered.
COEFF_MIN = -8
COEFF_MAX = 8

# A modular relation is interesting only when it predicts every
# instance and survives leave-one-out.
MIN_REQUIRED_INSTANCES = 8

# Keep printed candidates bounded.
MAX_PRINTED_CANDIDATES = 100


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

        if delta % 2 != 0:
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
# OBSERVABLE FEATURES
# ----------------------------------------------------------------------------------------------------------------------

Feature = Callable[[Instance], int]


def features(inst: Instance) -> dict[str, int]:
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


FEATURE_NAMES = tuple(features(INSTANCES[0]).keys())


# ----------------------------------------------------------------------------------------------------------------------
# MODULAR RELATION
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Candidate:
    modulus: int
    a: int
    feature_name: str
    b: int


def candidate_holds(candidate: Candidate, inst: Instance) -> bool:
    vals = features(inst)

    lhs = inst.true_j
    rhs = (
        candidate.a * vals[candidate.feature_name]
        + candidate.b
    )

    return (lhs - rhs) % candidate.modulus == 0


def all_hold(candidate: Candidate, instances) -> bool:
    return all(candidate_holds(candidate, inst) for inst in instances)


def held_out_test(candidate: Candidate, instances):
    failures = []

    for excluded in instances:
        remaining = [
            inst
            for inst in instances
            if inst.idx != excluded.idx
        ]

        if not all_hold(candidate, remaining):
            failures.append(excluded.idx)
            continue

        if not candidate_holds(candidate, excluded):
            failures.append(excluded.idx)

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def search_candidates(
    instances,
    modulus_min: int,
    modulus_max: int,
    coeff_min: int,
    coeff_max: int,
):
    candidates = []

    for m in range(modulus_min, modulus_max + 1):

        # A modulus of 1 carries no information.
        if m <= 1:
            continue

        # Canonical coefficient range modulo m.
        for a in range(coeff_min, coeff_max + 1):

            # a=0 produces constant residue tests and cannot encode
            # j through the feature.
            if a == 0:
                continue

            for feature_name in FEATURE_NAMES:

                vals = features(instances[0])
                feature_value = vals[feature_name]

                # Rather than exhaust b, infer the only b modulo m
                # needed from the first instance.
                b = (instances[0].true_j -
                     a * feature_value) % m

                candidate = Candidate(
                    modulus=m,
                    a=a,
                    feature_name=feature_name,
                    b=b,
                )

                if not all_hold(candidate, instances):
                    continue

                failures = held_out_test(candidate, instances)

                if failures:
                    continue

                # Canonicalize signs / coefficient duplicates.
                # The relation
                #   j = aF+b (mod m)
                # is equivalent modulo m to coefficient reduction.
                a_can = a % m
                b_can = b % m

                candidate = Candidate(
                    modulus=m,
                    a=a_can,
                    feature_name=feature_name,
                    b=b_can,
                )

                if candidate not in candidates:
                    candidates.append(candidate)

    return candidates


# ----------------------------------------------------------------------------------------------------------------------
# REDUNDANCY / TRIVIALITY TESTS
# ----------------------------------------------------------------------------------------------------------------------

def relation_is_trivial(candidate: Candidate, instances) -> bool:
    """
    Detects relations that merely exploit a feature already equal to
    j modulo m by construction of the diagnostic dataset.

    This does not reject them from the primary search; it labels them.
    """

    m = candidate.modulus
    vals = []

    for inst in instances:
        v = features(inst)[candidate.feature_name]
        vals.append((inst.true_j - v) % m)

    # a == 1 means this is directly j-feature.
    if candidate.a % m == 1:
        return len(set(vals)) == 1

    return False


def candidate_span(candidate: Candidate, instances):
    m = candidate.modulus

    pairs = []

    for inst in instances:
        f = features(inst)[candidate.feature_name] % m
        j = inst.true_j % m

        predicted = (
            candidate.a * f + candidate.b
        ) % m

        pairs.append(
            (inst.idx, j, f, predicted)
        )

    return pairs


# ----------------------------------------------------------------------------------------------------------------------
# CROSS-INSTANCE RESIDUE TABLE
# ----------------------------------------------------------------------------------------------------------------------

def print_residue_table(instances):
    print("=" * 120)
    print("BASIC RESIDUE TABLE")
    print("=" * 120)
    print()

    moduli = (
        2, 3, 4, 5, 6, 7, 8,
        9, 10, 11, 12, 16,
        24, 32, 64, 128, 256,
    )

    print(
        "  idx      j    "
        + " ".join(f"m={m:<3}" for m in moduli)
    )
    print("  " + "-" * 110)

    for inst in instances:
        values = [
            inst.true_j % m
            for m in moduli
        ]

        formatted = " ".join(
            f"{v:<6}"
            for v in values
        )

        print(
            f"  {inst.idx:<3d} "
            f"{inst.true_j:<6d}"
            f"{formatted}"
        )

    print()


# ----------------------------------------------------------------------------------------------------------------------
# FEATURE RESIDUE SCREEN
# ----------------------------------------------------------------------------------------------------------------------

def print_feature_residue_screen(instances):
    print("=" * 120)
    print("FEATURE / j RESIDUE SCREEN")
    print("=" * 120)
    print()

    moduli = (2, 3, 4, 5, 7, 8, 11, 13, 16, 31, 32)

    for name in FEATURE_NAMES:

        useful = []

        for m in moduli:
            j_residues = {
                inst.true_j % m
                for inst in instances
            }

            f_residues = {
                features(inst)[name] % m
                for inst in instances
            }

            if (
                len(j_residues) > 1
                and len(f_residues) > 1
            ):
                useful.append(m)

        if useful:
            print(
                f"  {name:18s} candidate moduli = {useful}"
            )

    print()


# ----------------------------------------------------------------------------------------------------------------------
# PRINT CANDIDATES
# ----------------------------------------------------------------------------------------------------------------------

def print_candidate(candidate: Candidate, instances, rank: int):
    print()
    print(f"  CANDIDATE {rank}")
    print(
        f"    j == "
        f"{candidate.a}*{candidate.feature_name}"
        f" + {candidate.b} (mod {candidate.modulus})"
    )

    print(
        f"    modulus     = {candidate.modulus}"
    )
    print(
        f"    coefficient = {candidate.a}"
    )
    print(
        f"    feature     = {candidate.feature_name}"
    )
    print(
        f"    offset      = {candidate.b}"
    )

    print()
    print("    instance    j    feature mod m    predicted mod m")
    print("    " + "-" * 65)

    m = candidate.modulus

    for inst in instances:
        f = features(inst)[candidate.feature_name] % m

        predicted = (
            candidate.a * f
            + candidate.b
        ) % m

        print(
            f"    {inst.idx:3d} "
            f"{inst.true_j:8d} "
            f"{f:16d} "
            f"{predicted:17d}"
        )

    trivial = relation_is_trivial(candidate, instances)

    print()
    print(
        f"    direct feature-offset relation = {trivial}"
    )


# ----------------------------------------------------------------------------------------------------------------------
# MODULUS COMPLEXITY SUMMARY
# ----------------------------------------------------------------------------------------------------------------------

def summarize_by_modulus(candidates):
    counts = {}

    for candidate in candidates:
        counts.setdefault(candidate.modulus, 0)
        counts[candidate.modulus] += 1

    return counts


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main():
    print("=" * 120)
    print("EXPERIMENT 463")
    print("=" * 120)
    print()
    print("CROSS-INSTANCE MODULAR j-LAW / OBSERVED-DATA AUDIT")
    print()
    print("QUESTION")
    print(
        "  Does j obey a simple modular law determined by"
    )
    print(
        "  the observed tuple across the eight instances?"
    )
    print()
    print(
        f"MODULI TESTED = [{MODULUS_MIN},{MODULUS_MAX}]"
    )
    print(
        f"COEFFICIENT RANGE = [{COEFF_MIN},{COEFF_MAX}]"
    )
    print()

    # Validate the supplied instances before doing the experiment.
    print("INPUT CONSISTENCY")
    print("  " + "-" * 80)

    consistency_failures = 0

    for inst in INSTANCES:
        if inst.d != inst.N - 2 * inst.S + 1:
            consistency_failures += 1
            print(
                f"  instance={inst.idx} d=N-2S+1 FAILED"
            )

        reconstructed_d0 = (
            inst.d0_ref + 2 * inst.true_j
        )

        reconstructed_x = (
            inst.x0 - 2 * inst.true_j
        )

        reconstructed_K = (
            inst.K_ref
            + inst.d0_ref * inst.true_j
            + inst.true_j ** 2
        )

        if reconstructed_d0 != inst.true_d0:
            consistency_failures += 1

        if reconstructed_x != inst.true_x:
            consistency_failures += 1

        if reconstructed_K != inst.true_K:
            consistency_failures += 1

    print(
        f"  consistency failures = {consistency_failures}"
    )
    print()

    if consistency_failures:
        raise AssertionError(
            "input consistency validation failed"
        )

    # ------------------------------------------------------------------
    # Residues
    # ------------------------------------------------------------------

    print_residue_table(INSTANCES)
    print_feature_residue_screen(INSTANCES)

    # ------------------------------------------------------------------
    # Candidate search
    # ------------------------------------------------------------------

    print("=" * 120)
    print("MODULAR CANDIDATE SEARCH")
    print("=" * 120)
    print()

    print(
        "Candidate form:"
    )
    print()
    print(
        "  j == a*F(N,S,d,r) + b (mod m)"
    )
    print()

    print(
        "Only very small a,b and m are searched."
    )
    print(
        "The relation must hold for all eight instances."
    )
    print(
        "It must also survive leave-one-out validation."
    )
    print()

    candidates = search_candidates(
        INSTANCES,
        MODULUS_MIN,
        MODULUS_MAX,
        COEFF_MIN,
        COEFF_MAX,
    )

    print(
        f"total surviving candidates = {len(candidates)}"
    )
    print()

    # ------------------------------------------------------------------
    # Candidate complexity analysis
    # ------------------------------------------------------------------

    modulus_counts = summarize_by_modulus(candidates)

    print("SURVIVING CANDIDATES BY MODULUS")
    print("  " + "-" * 80)

    if modulus_counts:
        for m in sorted(modulus_counts):
            print(
                f"  modulus={m:3d} "
                f"count={modulus_counts[m]}"
            )
    else:
        print("  none")

    print()

    # ------------------------------------------------------------------
    # Candidate listing
    # ------------------------------------------------------------------

    if candidates:
        print("=" * 120)
        print("SURVIVING MODULAR LAWS")
        print("=" * 120)

        for rank, candidate in enumerate(
            candidates[:MAX_PRINTED_CANDIDATES],
            start=1,
        ):
            print_candidate(
                candidate,
                INSTANCES,
                rank,
            )

        if len(candidates) > MAX_PRINTED_CANDIDATES:
            print()
            print(
                f"  additional candidates omitted = "
                f"{len(candidates)-MAX_PRINTED_CANDIDATES}"
            )

        print()
    else:
        print("=" * 120)
        print("NO SURVIVING MODULAR LAWS")
        print("=" * 120)
        print()
        print(
            "  No candidate of the tested form survived all eight"
        )
        print(
            "  instances and leave-one-out validation."
        )
        print()

    # ------------------------------------------------------------------
    # Most restrictive candidate audit
    # ------------------------------------------------------------------

    print("=" * 120)
    print("STRICT CANDIDATE AUDIT")
    print("=" * 120)
    print()

    strict_candidates = [
        c for c in candidates
        if c.modulus >= 16
        and abs(c.a) <= 4
        and c.feature_name not in ("1",)
    ]

    print(
        "  strict candidate definition:"
    )
    print(
        "    modulus >= 16"
    )
    print(
        "    |a| <= 4"
    )
    print(
        "    feature is nonconstant"
    )
    print()

    print(
        f"  strict candidates = {len(strict_candidates)}"
    )

    if strict_candidates:
        for i, c in enumerate(
            strict_candidates[:MAX_PRINTED_CANDIDATES],
            start=1,
        ):
            print(
                f"  {i:3d}. "
                f"j == {c.a}*{c.feature_name}+{c.b}"
                f" (mod {c.modulus})"
            )
    else:
        print("  none")

    print()

    # ------------------------------------------------------------------
    # Hidden-vs-observed comparison
    # ------------------------------------------------------------------

    print("=" * 120)
    print("HIDDEN-VARIABLE CONTRAST")
    print("=" * 120)
    print()

    hidden_names = (
        "j",
        "j^2",
        "2j",
        "2j+1",
    )

    for name in hidden_names:
        if name == "j":
            values = [inst.true_j for inst in INSTANCES]
        elif name == "j^2":
            values = [inst.true_j ** 2 for inst in INSTANCES]
        elif name == "2j":
            values = [2 * inst.true_j for inst in INSTANCES]
        else:
            values = [2 * inst.true_j + 1 for inst in INSTANCES]

        print(
            f"  {name:8s} = {values}"
        )

    print()
    print(
        "  The hidden sequence remains visibly nonconstant."
    )
    print(
        "  The modular search asks whether its residue pattern"
    )
    print(
        "  is explained by the observed tuple."
    )
    print()

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 463 SUMMARY")
    print("=" * 120)
    print()

    print(
        f"  instances                         = {len(INSTANCES)}"
    )
    print(
        f"  feature count                     = {len(FEATURE_NAMES)}"
    )
    print(
        f"  modulus range                     = "
        f"{MODULUS_MIN}..{MODULUS_MAX}"
    )
    print(
        f"  coefficient range                 = "
        f"{COEFF_MIN}..{COEFF_MAX}"
    )
    print(
        f"  surviving modular candidates      = "
        f"{len(candidates)}"
    )
    print(
        f"  strict surviving candidates       = "
        f"{len(strict_candidates)}"
    )
    print(
        f"  input consistency failures        = "
        f"{consistency_failures}"
    )
    print()

    print("INTERPRETATION")
    print()
    print(
        "  Experiment 461 established that expressions built"
    )
    print(
        "  only from N,S,d,r are constant within one fibre."
    )
    print()
    print(
        "  Experiment 462 found no small exact affine or"
    )
    print(
        "  quadratic cross-instance law for j."
    )
    print()
    print(
        "  Experiment 463 tests a different algebraic class:"
    )
    print(
        "  modular relations between j and observable features."
    )
    print()
    print(
        "  A surviving modular law would not prove recoverability."
    )
    print(
        "  It would only identify a candidate cross-instance"
    )
    print(
        "  structure worthy of independent testing."
    )
    print()
    print(
        "  Conversely, absence of a law in this bounded search"
    )
    print(
        "  does not prove that no more complicated modular law exists."
    )
    print()

    if not candidates:
        print("PRIMARY OUTCOME")
        print()
        print(
            "  No tested small modular relation connects j to"
        )
        print(
            "  the observed features across all eight instances."
        )
        print()
    else:
        print("PRIMARY OUTCOME")
        print()
        print(
            f"  {len(candidates)} modular candidate(s) survived."
        )
        print(
            "  These must be treated as empirical candidates,"
        )
        print(
            "  not as an established reconstruction rule."
        )
        print()

    print("=" * 120)
    print("EXPERIMENT 463 FINAL STATUS")
    print("=" * 120)
    print()
    print(
        "  EXACT INTEGER ARITHMETIC        = True"
    )
    print(
        "  MODULAR SEARCH COMPLETED        = True"
    )
    print(
        f"  MODULAR CANDIDATES              = {len(candidates)}"
    )
    print(
        f"  STRICT CANDIDATES               = {len(strict_candidates)}"
    )
    print(
        "  LEAVE-ONE-OUT VALIDATION        = True"
    )
    print(
        "  GIANT K ENUMERATION             = False"
    )
    print(
        "  GIANT d0 ENUMERATION            = False"
    )
    print(
        "  CRT CARTESIAN PRODUCT           = False"
    )
    print(
        "  FLOATING POINT                  = False"
    )
    print(
        "  ARBITRARY INTERPOLATION         = False"
    )
    print()
    print("  CONCLUSION:")
    print(
        "    This experiment tests whether the hidden displacement"
    )
    print(
        "    j has a simple modular relationship to the observed"
    )
    print(
        "    tuple, after ordinary polynomial cross-instance laws"
    )
    print(
        "    have already failed."
    )
    print()
    print("EXPERIMENT 463 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
