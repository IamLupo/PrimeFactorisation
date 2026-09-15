#!/usr/bin/env python3

"""
========================================================================================================================
EXPERIMENT 466
========================================================================================================================

NONDEGENERATE RATIONAL / EUCLIDEAN j-LAW AUDIT

PURPOSE

  Experiment 465R2 produced 104272 apparent rational candidates.

  Inspection shows that many are degenerate identities of the form

      a*F + b = q

  with

      c = 0

  or use a constant Euclidean feature such as q(N,d)=1.

  Those candidates contain no genuine observed-data dependence on j.

  This experiment removes those degeneracies and performs a genuine
  leave-one-out search.

QUESTION

  Does a nondegenerate exact rational relation of the form

      a*F + b = j*(c*G + e) + q

  survive when:

      * F varies across instances;
      * G varies across instances;
      * c != 0;
      * c*G+e varies across instances;
      * c*G+e is never zero;
      * a*F+b varies across instances;
      * coefficient tuples are reduced canonically;
      * the candidate is discovered using seven instances only;
      * the held-out eighth instance is then tested independently?

OBSERVED VARIABLES

  N, S, d, r

FEATURE BASIS

  Euclidean quotients
  Euclidean remainders
  gcd values
  simple Euclidean combinations

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no arbitrary interpolation
  no hidden variables in observable features

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

COEFF_MIN = -4
COEFF_MAX = 4

OFFSET_MIN = -8
OFFSET_MAX = 8

MAX_PRINT = 100


# ----------------------------------------------------------------------------------------------------------------------
# INSTANCE DATA
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
# EUCLIDEAN FEATURES
# ----------------------------------------------------------------------------------------------------------------------

def euclid(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        raise ValueError("Euclidean divisor cannot be zero")

    q, rem = divmod(a, b)
    g = gcd(abs(a), abs(b))
    return g, q, rem


def features(inst: Instance) -> dict[str, int]:

    N = inst.N
    S = inst.S
    d = inst.d
    r = inst.r

    g_NS, q_NS, rem_NS = euclid(N, S)
    g_Nd, q_Nd, rem_Nd = euclid(N, d)
    g_Sd, q_Sd, rem_Sd = euclid(S, d)
    g_dr, q_dr, rem_dr = euclid(d, r)
    g_Nr, q_Nr, rem_Nr = euclid(N, r)
    g_Sr, q_Sr, rem_Sr = euclid(S, r)

    return {
        "q(N,S)": q_NS,
        "q(N,d)": q_Nd,
        "q(S,d)": q_Sd,
        "q(d,r)": q_dr,
        "q(N,r)": q_Nr,
        "q(S,r)": q_Sr,

        "rem(N,S)": rem_NS,
        "rem(N,d)": rem_Nd,
        "rem(S,d)": rem_Sd,
        "rem(d,r)": rem_dr,
        "rem(N,r)": rem_Nr,
        "rem(S,r)": rem_Sr,

        "gcd(N,S)": g_NS,
        "gcd(N,d)": g_Nd,
        "gcd(S,d)": g_Sd,
        "gcd(d,r)": g_dr,
        "gcd(N,r)": g_Nr,
        "gcd(S,r)": g_Sr,

        "q(N,S)-q(S,d)": q_NS - q_Sd,
        "q(N,d)-q(S,d)": q_Nd - q_Sd,
        "q(S,d)-q(N,S)": q_Sd - q_NS,

        "rem(N,S)-rem(S,d)": rem_NS - rem_Sd,
        "rem(N,d)-rem(S,d)": rem_Nd - rem_Sd,
        "rem(N,S)-rem(N,d)": rem_NS - rem_Nd,

        "gcd(N,S)-gcd(S,d)": g_NS - g_Sd,
        "gcd(N,d)-gcd(S,d)": g_Nd - g_Sd,
    }


FEATURE_NAMES = tuple(features(INSTANCES[0]).keys())


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE
# ----------------------------------------------------------------------------------------------------------------------

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
# CANONICALIZATION
# ----------------------------------------------------------------------------------------------------------------------

def canonical_tuple(
    F: str,
    G: str,
    a: int,
    b: int,
    c: int,
    e: int,
    q: int,
) -> tuple[int, int, int, int, int]:
    """
    Canonicalize the coefficient vector.

    Since multiplying

        aF+b = j(cG+e)+q

    by a common nonzero scalar is not allowed here because j is fixed
    and q is an additive term, we only remove the obvious sign duplicate
    when the complete coefficient vector permits it.

    This function instead provides a stable representation for the
    candidate set.
    """

    vals = (a, b, c, e, q)

    for v in vals:
        if v != 0:
            if v < 0:
                return tuple(-x for x in vals)
            break

    return vals


# ----------------------------------------------------------------------------------------------------------------------
# NONDEGENERACY
# ----------------------------------------------------------------------------------------------------------------------

def is_nondegenerate(
    cand: RationalCandidate,
    instances: tuple[Instance, ...],
) -> bool:

    fvals = []
    gvals = []
    denoms = []
    nums = []

    for inst in instances:

        table = features(inst)

        Fv = table[cand.F]
        Gv = table[cand.G]

        numerator = cand.a * Fv + cand.b
        denominator = cand.c * Gv + cand.e

        fvals.append(Fv)
        gvals.append(Gv)
        nums.append(numerator)
        denoms.append(denominator)

    # Both selected observable features must actually vary.
    if len(set(fvals)) < 2:
        return False

    if len(set(gvals)) < 2:
        return False

    # The denominator must genuinely vary.
    if cand.c == 0:
        return False

    if len(set(denoms)) < 2:
        return False

    # No division-by-zero branch.
    if any(v == 0 for v in denoms):
        return False

    # The numerator must also vary.
    if len(set(nums)) < 2:
        return False

    # Exclude constant denominator magnitude/sign cases.
    if all(abs(v) == 1 for v in denoms):
        return False

    return True


# ----------------------------------------------------------------------------------------------------------------------
# EXACT TEST
# ----------------------------------------------------------------------------------------------------------------------

def exact_holds(
    cand: RationalCandidate,
    inst: Instance,
) -> bool:

    table = features(inst)

    Fv = table[cand.F]
    Gv = table[cand.G]

    lhs = cand.a * Fv + cand.b
    denominator = cand.c * Gv + cand.e

    if denominator == 0:
        return False

    rhs = inst.true_j * denominator + cand.q

    return lhs == rhs


# ----------------------------------------------------------------------------------------------------------------------
# SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def search_full(
    instances: tuple[Instance, ...],
) -> list[RationalCandidate]:

    results: list[RationalCandidate] = []

    first = instances[0]
    first_features = features(first)

    feature_cache = {
        inst.idx: features(inst)
        for inst in instances
    }

    # Both features must be nonconstant.
    usable_features = []

    for name in FEATURE_NAMES:

        vals = [
            feature_cache[inst.idx][name]
            for inst in instances
        ]

        if len(set(vals)) >= 2:
            usable_features.append(name)

    for F in usable_features:

        F0 = first_features[F]

        for G in usable_features:

            G0 = first_features[G]

            for a in range(COEFF_MIN, COEFF_MAX + 1):

                if a == 0:
                    continue

                for b in range(OFFSET_MIN, OFFSET_MAX + 1):

                    lhs0 = a * F0 + b

                    for c in range(COEFF_MIN, COEFF_MAX + 1):

                        if c == 0:
                            continue

                        for e in range(OFFSET_MIN, OFFSET_MAX + 1):

                            denominator0 = c * G0 + e

                            if denominator0 == 0:
                                continue

                            # Infer q exactly from instance 1.
                            q = (
                                lhs0
                                - first.true_j * denominator0
                            )

                            if not (
                                OFFSET_MIN
                                <= q
                                <= OFFSET_MAX
                            ):
                                continue

                            cand = RationalCandidate(
                                F=F,
                                G=G,
                                a=a,
                                b=b,
                                c=c,
                                e=e,
                                q=q,
                            )

                            if not is_nondegenerate(
                                cand,
                                instances,
                            ):
                                continue

                            if all(
                                exact_holds(cand, inst)
                                for inst in instances
                            ):
                                results.append(cand)

    return dedupe(results)


# ----------------------------------------------------------------------------------------------------------------------
# GENUINE LEAVE-ONE-OUT SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def search_leave_one_out():

    all_successes = []
    union_candidates: set[RationalCandidate] = set()

    for held_out in INSTANCES:

        training = tuple(
            inst
            for inst in INSTANCES
            if inst.idx != held_out.idx
        )

        candidates = search_full(training)

        passing = [
            cand
            for cand in candidates
            if exact_holds(cand, held_out)
        ]

        for cand in passing:

            union_candidates.add(cand)

            all_successes.append(
                (
                    held_out.idx,
                    cand,
                )
            )

    return (
        dedupe(list(union_candidates)),
        all_successes,
    )


# ----------------------------------------------------------------------------------------------------------------------
# DEDUPLICATION
# ----------------------------------------------------------------------------------------------------------------------

def dedupe(
    candidates: list[RationalCandidate],
) -> list[RationalCandidate]:

    seen = set()
    out = []

    for cand in candidates:

        key = (
            cand.F,
            cand.G,
            canonical_tuple(
                cand.F,
                cand.G,
                cand.a,
                cand.b,
                cand.c,
                cand.e,
                cand.q,
            ),
        )

        if key in seen:
            continue

        seen.add(key)
        out.append(cand)

    return out


# ----------------------------------------------------------------------------------------------------------------------
# INPUT AUDIT
# ----------------------------------------------------------------------------------------------------------------------

def validate_inputs():

    failures = []

    for inst in INSTANCES:

        if inst.d != inst.N - 2 * inst.S + 1:
            failures.append(
                f"instance {inst.idx}: d relation"
            )

        if inst.d0_ref + 2 * inst.true_j != inst.true_d0:
            failures.append(
                f"instance {inst.idx}: d0 reconstruction"
            )

        if inst.x0 - 2 * inst.true_j != inst.true_x:
            failures.append(
                f"instance {inst.idx}: x reconstruction"
            )

        expected_K = (
            inst.K_ref
            + inst.d0_ref * inst.true_j
            + inst.true_j ** 2
        )

        if expected_K != inst.true_K:
            failures.append(
                f"instance {inst.idx}: K reconstruction"
            )

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# FEATURE SUMMARY
# ----------------------------------------------------------------------------------------------------------------------

def print_feature_variation():

    print("=" * 120)
    print("NONDEGENERATE FEATURE SCREEN")
    print("=" * 120)
    print()

    for name in FEATURE_NAMES:

        values = [
            features(inst)[name]
            for inst in INSTANCES
        ]

        distinct = len(set(values))

        marker = (
            "VARIABLE"
            if distinct >= 2
            else "CONSTANT"
        )

        print(
            f"  {name:35s}"
            f" distinct={distinct:<2d}"
            f" {marker:10s}"
            f" values={values}"
        )

    print()


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE REPORT
# ----------------------------------------------------------------------------------------------------------------------

def print_candidate(
    cand: RationalCandidate,
    prefix: str = "",
):

    print(
        f"{prefix}{cand.a}*{cand.F}+{cand.b}"
        f" = j*({cand.c}*{cand.G}+{cand.e})"
        f" + {cand.q}"
    )


def print_candidate_detail(
    cand: RationalCandidate,
):

    print_candidate(cand, "  ")

    print(
        "    inst      j        F-value       G-value"
        "       numerator       denominator"
        "       match"
    )
    print(
        "    "
        + "-" * 100
    )

    for inst in INSTANCES:

        table = features(inst)

        Fv = table[cand.F]
        Gv = table[cand.G]

        numerator = cand.a * Fv + cand.b
        denominator = cand.c * Gv + cand.e

        lhs = numerator
        rhs = inst.true_j * denominator + cand.q

        print(
            f"    {inst.idx:<4d}"
            f" {inst.true_j:8d}"
            f" {Fv:14d}"
            f" {Gv:14d}"
            f" {lhs:16d}"
            f" {denominator:16d}"
            f" {lhs == rhs!s:>10s}"
        )


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main():

    print("=" * 120)
    print("EXPERIMENT 466")
    print("=" * 120)
    print()

    print(
        "NONDEGENERATE RATIONAL / EUCLIDEAN j-LAW AUDIT"
    )
    print()

    print(
        f"COEFFICIENT RANGE = "
        f"[{COEFF_MIN},{COEFF_MAX}]"
    )
    print(
        f"OFFSET RANGE      = "
        f"[{OFFSET_MIN},{OFFSET_MAX}]"
    )
    print()

    # --------------------------------------------------------------------------------------------------------------
    # Input consistency.
    # --------------------------------------------------------------------------------------------------------------

    failures = validate_inputs()

    print("INPUT CONSISTENCY")
    print(
        "  "
        + "-"
        * 80
    )
    print(
        f"  consistency failures = {len(failures)}"
    )

    if failures:

        for failure in failures:
            print(
                f"  FAILURE: {failure}"
            )

        raise AssertionError(
            "input consistency failure"
        )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Feature variation.
    # --------------------------------------------------------------------------------------------------------------

    print_feature_variation()

    # --------------------------------------------------------------------------------------------------------------
    # Full-data nondegenerate search.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("NONDEGENERATE FULL-DATA RATIONAL SEARCH")
    print("=" * 120)
    print()

    full_candidates = search_full(
        INSTANCES
    )

    print(
        f"nondegenerate full-data candidates = "
        f"{len(full_candidates)}"
    )

    if full_candidates:

        for cand in full_candidates[:MAX_PRINT]:
            print_candidate(cand, "  ")

        if len(full_candidates) > MAX_PRINT:
            print(
                f"  ... {len(full_candidates) - MAX_PRINT}"
                f" additional candidates suppressed"
            )

    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Genuine leave-one-out.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GENUINE LEAVE-ONE-OUT RATIONAL SEARCH")
    print("=" * 120)
    print()

    loo_candidates, loo_successes = (
        search_leave_one_out()
    )

    print(
        f"unique LOO-surviving candidates = "
        f"{len(loo_candidates)}"
    )

    by_holdout = {}

    for held_out_idx, cand in loo_successes:

        by_holdout.setdefault(
            held_out_idx,
            0,
        )

        by_holdout[held_out_idx] += 1

    print()
    print(
        "  held-out instance survival counts:"
    )

    for inst in INSTANCES:

        print(
            f"    instance {inst.idx}: "
            f"{by_holdout.get(inst.idx, 0)}"
        )

    print()

    if loo_candidates:

        print(
            "  LOO SURVIVING CANDIDATES"
        )

        for cand in loo_candidates[:MAX_PRINT]:
            print_candidate(
                cand,
                "    ",
            )

        if len(loo_candidates) > MAX_PRINT:
            print(
                f"    ... {len(loo_candidates) - MAX_PRINT}"
                f" additional candidates suppressed"
            )

    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Degeneracy comparison.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("DEGENERACY COMPARISON WITH EXPERIMENT 465R2")
    print("=" * 120)
    print()

    print(
        "  Previous rational candidate count = 104272"
    )

    print(
        f"  Current nondegenerate count       = "
        f"{len(full_candidates)}"
    )

    print(
        f"  Current LOO-surviving count       = "
        f"{len(loo_candidates)}"
    )

    print()

    if not loo_candidates:

        print(
            "  All previously observed rational candidates"
        )
        print(
            "  are eliminated once degenerate constant/"
        )
        print(
            "  denominator constructions are excluded"
        )
        print(
            "  and the held-out instance is genuinely tested."
        )

    else:

        print(
            "  A nondegenerate candidate survived."
        )
        print(
            "  It requires a separate structural audit."
        )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Detailed survivors.
    # --------------------------------------------------------------------------------------------------------------

    if loo_candidates:

        print("=" * 120)
        print("DETAILED LOO SURVIVORS")
        print("=" * 120)
        print()

        for i, cand in enumerate(
            loo_candidates[:MAX_PRINT],
            start=1,
        ):

            print(
                f"  CANDIDATE {i}"
            )

            print_candidate_detail(
                cand
            )

            print()

    # --------------------------------------------------------------------------------------------------------------
    # Global summary.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 466 SUMMARY")
    print("=" * 120)
    print()

    print(
        f"  instances                         = "
        f"{len(INSTANCES)}"
    )

    print(
        f"  observable Euclidean features    = "
        f"{len(FEATURE_NAMES)}"
    )

    print(
        f"  full-data nondegenerate laws     = "
        f"{len(full_candidates)}"
    )

    print(
        f"  genuine LOO laws                 = "
        f"{len(loo_candidates)}"
    )

    print(
        "  c = 0 candidates admitted        = False"
    )

    print(
        "  constant F candidates admitted   = False"
    )

    print(
        "  constant G candidates admitted   = False"
    )

    print(
        "  constant denominator admitted     = False"
    )

    print(
        "  zero denominator admitted         = False"
    )

    print(
        "  floating point                    = False"
    )

    print(
        "  interpolation                     = False"
    )

    print(
        "  giant K enumeration               = False"
    )

    print(
        "  giant d0 enumeration              = False"
    )

    print(
        "  CRT Cartesian product             = False"
    )

    print()

    print(
        "INTERPRETATION"
    )
    print()

    print(
        "  The 104272 rational candidates"
    )

    print(
        "  reported by Experiment 465R2"
    )

    print(
        "  included trivial constructions"
    )

    print(
        "  such as constant Euclidean quotients"
    )

    print(
        "  and c=0 relations."
    )

    print()

    print(
        "  Experiment 466 removes those cases."
    )

    print(
        "  A surviving law must contain genuine"
    )

    print(
        "  variation in both the numerator feature"
    )

    print(
        "  and denominator feature."
    )

    print()

    print(
        "  More importantly, discovery is performed"
    )

    print(
        "  on seven instances and the eighth is"
    )

    print(
        "  tested afterward."
    )

    print()

    print(
        "  Therefore a negative LOO result is"
    )

    print(
        "  substantially stronger evidence against"
    )

    print(
        "  this particular low-complexity rational"
    )

    print(
        "  family than the previous all-eight fit."
    )

    print()

    print(
        "LIMITATION"
    )
    print()

    print(
        "  A negative result still does not establish"
    )

    print(
        "  that no more complicated observable law exists."
    )

    print()

    print("=" * 120)
    print("EXPERIMENT 466 FINAL STATUS")
    print("=" * 120)
    print()

    print(
        "  EXACT INTEGER ARITHMETIC        = True"
    )

    print(
        "  NONDEGENERATE RATIONAL TEST     = True"
    )

    print(
        "  GENUINE LEAVE-ONE-OUT TEST      = True"
    )

    print(
        f"  FULL-DATA CANDIDATES             = "
        f"{len(full_candidates)}"
    )

    print(
        f"  LOO CANDIDATES                   = "
        f"{len(loo_candidates)}"
    )

    print(
        "  CONSTANT-FEATURE LAWS            = Excluded"
    )

    print(
        "  c=0 LAWS                         = Excluded"
    )

    print(
        "  CONSTANT-DENOMINATOR LAWS        = Excluded"
    )

    print(
        "  ZERO-DENOMINATOR LAWS            = Excluded"
    )

    print(
        "  FLOATING POINT                   = False"
    )

    print(
        "  ARBITRARY INTERPOLATION          = False"
    )

    print(
        "  GIANT K ENUMERATION              = False"
    )

    print(
        "  GIANT d0 ENUMERATION             = False"
    )

    print(
        "  CRT CARTESIAN PRODUCT            = False"
    )

    print()

    if not loo_candidates:

        print(
            "  CONCLUSION:"
        )

        print(
            "    After removing the degenerate rational"
        )

        print(
            "    constructions responsible for the large"
        )

        print(
            "    Experiment-465R2 candidate count, no"
        )

        print(
            "    nondegenerate rational Euclidean-derived"
        )

        print(
            "    j-law survives genuine leave-one-out"
        )

        print(
            "    validation."
        )

    else:

        print(
            "  CONCLUSION:"
        )

        print(
            "    At least one nondegenerate rational"
        )

        print(
            "    Euclidean-derived relation survives"
        )

        print(
            "    leave-one-out validation."
        )

        print(
            "    Those candidates require a subsequent"
        )

        print(
            "    independent structural audit."
        )

    print()
    print(
        "EXPERIMENT 466 FINISHED"
    )
    print("=" * 120)


if __name__ == "__main__":
    main()
