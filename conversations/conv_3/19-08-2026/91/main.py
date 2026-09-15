#!/usr/bin/env python3

"""
========================================================================================================================
EXPERIMENT 467
========================================================================================================================

TWO-FEATURE RATIONAL / EUCLIDEAN j-LAW AUDIT

PURPOSE

  Experiment 466 tested

      a*F+b = j*(c*G+e)+q

  using one observable feature in the numerator and one in the
  denominator.

  Experiment 466 found no nondegenerate candidate.

  The next escalation permits TWO independently varying observable
  features in both numerator and denominator:

      a*F + b*G + c
        =
      j*(d*H + e*I + f) + q

  where all coefficients are small bounded integers.

  The objective is to test whether j can be represented by a
  genuinely multivariate exact rational law in the Euclidean-derived
  observable feature space.

NONDEGENERACY REQUIREMENTS

  1. F and G must vary across the training data.
  2. H and I must vary across the training data.
  3. The numerator must vary.
  4. The denominator must vary.
  5. The denominator must never vanish.
  6. At least one coefficient of each feature pair must be nonzero.
  7. Constant-only constructions are rejected.
  8. Scalar/sign duplicates are canonicalized.
  9. A candidate is discovered on seven instances and tested on
     the held-out eighth instance.

OBSERVED VARIABLES

  N, S, d, r

EUCLIDEAN FEATURE BASIS

  Euclidean quotients
  Euclidean remainders
  gcd values
  previously constructed Euclidean differences

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no arbitrary interpolation
  no hidden variables in feature construction

COEFFICIENT RANGE = [-3,3]
OFFSET RANGE      = [-6,6]

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from itertools import product


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

COEFF_MIN = -3
COEFF_MAX = 3

OFFSET_MIN = -6
OFFSET_MAX = 6

MAX_PRINT = 50


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
    return gcd(abs(a), abs(b)), q, rem


def make_features(inst: Instance) -> dict[str, int]:

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


FEATURE_NAMES = tuple(
    make_features(INSTANCES[0]).keys()
)


# ----------------------------------------------------------------------------------------------------------------------
# CACHE
# ----------------------------------------------------------------------------------------------------------------------

FEATURE_CACHE: dict[int, dict[str, int]] = {
    inst.idx: make_features(inst)
    for inst in INSTANCES
}


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class TwoFeatureRationalCandidate:
    F: str
    G: str
    H: str
    I: str

    a: int
    b: int
    c: int

    d: int
    e: int
    f: int

    q: int


# ----------------------------------------------------------------------------------------------------------------------
# CANONICALIZATION
# ----------------------------------------------------------------------------------------------------------------------

def canonical_coefficients(
    a: int,
    b: int,
    c: int,
    d: int,
    e: int,
    f: int,
    q: int,
) -> tuple[int, int, int, int, int, int, int]:

    vals = (a, b, c, d, e, f, q)

    for value in vals:

        if value == 0:
            continue

        if value < 0:
            return tuple(-x for x in vals)

        break

    return vals


# ----------------------------------------------------------------------------------------------------------------------
# TRAINING FEATURE LIST
# ----------------------------------------------------------------------------------------------------------------------

def varying_features(
    instances: tuple[Instance, ...],
) -> list[str]:

    result = []

    for name in FEATURE_NAMES:

        values = [
            FEATURE_CACHE[inst.idx][name]
            for inst in instances
        ]

        if len(set(values)) >= 2:
            result.append(name)

    return result


# ----------------------------------------------------------------------------------------------------------------------
# EXACT EVALUATION
# ----------------------------------------------------------------------------------------------------------------------

def evaluate_candidate(
    cand: TwoFeatureRationalCandidate,
    inst: Instance,
) -> tuple[int, int, int, int]:

    table = FEATURE_CACHE[inst.idx]

    Fv = table[cand.F]
    Gv = table[cand.G]
    Hv = table[cand.H]
    Iv = table[cand.I]

    numerator = (
        cand.a * Fv
        + cand.b * Gv
        + cand.c
    )

    denominator = (
        cand.d * Hv
        + cand.e * Iv
        + cand.f
    )

    lhs = numerator
    rhs = inst.true_j * denominator + cand.q

    return (
        lhs,
        rhs,
        numerator,
        denominator,
    )


def exact_holds(
    cand: TwoFeatureRationalCandidate,
    inst: Instance,
) -> bool:

    lhs, rhs, _, denominator = evaluate_candidate(
        cand,
        inst,
    )

    if denominator == 0:
        return False

    return lhs == rhs


# ----------------------------------------------------------------------------------------------------------------------
# NONDEGENERACY
# ----------------------------------------------------------------------------------------------------------------------

def is_nondegenerate(
    cand: TwoFeatureRationalCandidate,
    instances: tuple[Instance, ...],
) -> bool:

    numerator_values = []
    denominator_values = []

    F_values = []
    G_values = []
    H_values = []
    I_values = []

    for inst in instances:

        table = FEATURE_CACHE[inst.idx]

        Fv = table[cand.F]
        Gv = table[cand.G]
        Hv = table[cand.H]
        Iv = table[cand.I]

        F_values.append(Fv)
        G_values.append(Gv)
        H_values.append(Hv)
        I_values.append(Iv)

        numerator_values.append(
            cand.a * Fv
            + cand.b * Gv
            + cand.c
        )

        denominator_values.append(
            cand.d * Hv
            + cand.e * Iv
            + cand.f
        )

    # Each selected input feature must vary.
    if len(set(F_values)) < 2:
        return False

    if len(set(G_values)) < 2:
        return False

    if len(set(H_values)) < 2:
        return False

    if len(set(I_values)) < 2:
        return False

    # Numerator must genuinely vary.
    if len(set(numerator_values)) < 2:
        return False

    # Denominator must genuinely vary.
    if len(set(denominator_values)) < 2:
        return False

    # Denominator may never vanish.
    if any(value == 0 for value in denominator_values):
        return False

    # Reject numerator with no actual feature contribution.
    if cand.a == 0 and cand.b == 0:
        return False

    # Reject denominator with no actual feature contribution.
    if cand.d == 0 and cand.e == 0:
        return False

    # Require two-feature structure on both sides.
    if cand.a == 0 or cand.b == 0:
        return False

    if cand.d == 0 or cand.e == 0:
        return False

    return True


# ----------------------------------------------------------------------------------------------------------------------
# SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def search_training(
    training: tuple[Instance, ...],
) -> list[TwoFeatureRationalCandidate]:

    usable = varying_features(training)

    # Pairings are ordered because numerator and denominator roles differ.
    feature_pairs = [
        (F, G)
        for F in usable
        for G in usable
        if F != G
    ]

    candidates = []

    first = training[0]
    first_table = FEATURE_CACHE[first.idx]

    # Search feature pairs.
    for F, G in feature_pairs:

        F0 = first_table[F]
        G0 = first_table[G]

        # Numerator coefficient tuples.
        for a, b in product(
            range(COEFF_MIN, COEFF_MAX + 1),
            repeat=2,
        ):

            if a == 0 or b == 0:
                continue

            for c in range(
                OFFSET_MIN,
                OFFSET_MAX + 1,
            ):

                lhs0 = (
                    a * F0
                    + b * G0
                    + c
                )

                # Denominator feature pair.
                for H, I in feature_pairs:

                    H0 = first_table[H]
                    I0 = first_table[I]

                    for d, e in product(
                        range(COEFF_MIN, COEFF_MAX + 1),
                        repeat=2,
                    ):

                        if d == 0 or e == 0:
                            continue

                        for f in range(
                            OFFSET_MIN,
                            OFFSET_MAX + 1,
                        ):

                            den0 = (
                                d * H0
                                + e * I0
                                + f
                            )

                            if den0 == 0:
                                continue

                            # Infer q exactly from the first
                            # training point.
                            q = (
                                lhs0
                                - first.true_j * den0
                            )

                            if not (
                                OFFSET_MIN
                                <= q
                                <= OFFSET_MAX
                            ):
                                continue

                            cand = TwoFeatureRationalCandidate(
                                F=F,
                                G=G,
                                H=H,
                                I=I,
                                a=a,
                                b=b,
                                c=c,
                                d=d,
                                e=e,
                                f=f,
                                q=q,
                            )

                            if not is_nondegenerate(
                                cand,
                                training,
                            ):
                                continue

                            valid = True

                            for inst in training:

                                if not exact_holds(
                                    cand,
                                    inst,
                                ):
                                    valid = False
                                    break

                            if valid:
                                candidates.append(cand)

    return dedupe(candidates)


# ----------------------------------------------------------------------------------------------------------------------
# DEDUPLICATION
# ----------------------------------------------------------------------------------------------------------------------

def dedupe(
    candidates: list[TwoFeatureRationalCandidate],
) -> list[TwoFeatureRationalCandidate]:

    seen = set()
    result = []

    for cand in candidates:

        coeffs = canonical_coefficients(
            cand.a,
            cand.b,
            cand.c,
            cand.d,
            cand.e,
            cand.f,
            cand.q,
        )

        key = (
            cand.F,
            cand.G,
            cand.H,
            cand.I,
            coeffs,
        )

        if key in seen:
            continue

        seen.add(key)
        result.append(cand)

    return result


# ----------------------------------------------------------------------------------------------------------------------
# LEAVE-ONE-OUT SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def leave_one_out():

    successful_candidates = set()

    held_out_counts = {
        inst.idx: 0
        for inst in INSTANCES
    }

    for held_out in INSTANCES:

        training = tuple(
            inst
            for inst in INSTANCES
            if inst.idx != held_out.idx
        )

        training_candidates = search_training(
            training
        )

        passing = []

        for cand in training_candidates:

            if exact_holds(
                cand,
                held_out,
            ):

                passing.append(cand)
                successful_candidates.add(cand)

        held_out_counts[
            held_out.idx
        ] = len(passing)

    return (
        sorted(
            successful_candidates,
            key=lambda c: (
                c.F,
                c.G,
                c.H,
                c.I,
                canonical_coefficients(
                    c.a,
                    c.b,
                    c.c,
                    c.d,
                    c.e,
                    c.f,
                    c.q,
                ),
            ),
        ),
        held_out_counts,
    )


# ----------------------------------------------------------------------------------------------------------------------
# FULL-DATA SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def full_data_search():

    return dedupe(
        search_training(INSTANCES)
    )


# ----------------------------------------------------------------------------------------------------------------------
# INPUT VALIDATION
# ----------------------------------------------------------------------------------------------------------------------

def validate_inputs():

    failures = []

    for inst in INSTANCES:

        if inst.d != inst.N - 2 * inst.S + 1:
            failures.append(
                f"instance {inst.idx}: d != N-2S+1"
            )

        j = inst.true_j

        d0 = inst.d0_ref + 2 * j
        x = inst.x0 - 2 * j
        K = (
            inst.K_ref
            + inst.d0_ref * j
            + j * j
        )

        if d0 != inst.true_d0:
            failures.append(
                f"instance {inst.idx}: d0 reconstruction"
            )

        if x != inst.true_x:
            failures.append(
                f"instance {inst.idx}: x reconstruction"
            )

        if K != inst.true_K:
            failures.append(
                f"instance {inst.idx}: K reconstruction"
            )

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# FEATURE REPORT
# ----------------------------------------------------------------------------------------------------------------------

def print_feature_screen():

    print("=" * 120)
    print(
        "TWO-FEATURE OBSERVABLE SCREEN"
    )
    print("=" * 120)
    print()

    for name in FEATURE_NAMES:

        values = [
            FEATURE_CACHE[inst.idx][name]
            for inst in INSTANCES
        ]

        distinct = len(set(values))

        status = (
            "VARIABLE"
            if distinct >= 2
            else "CONSTANT"
        )

        print(
            f"  {name:35s}"
            f" distinct={distinct:<2d}"
            f" {status:10s}"
        )

    print()


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE PRINTING
# ----------------------------------------------------------------------------------------------------------------------

def print_candidate(
    cand: TwoFeatureRationalCandidate,
    prefix: str = "",
):

    numerator = (
        f"{cand.a}*{cand.F}"
        f" + {cand.b}*{cand.G}"
        f" + {cand.c}"
    )

    denominator = (
        f"{cand.d}*{cand.H}"
        f" + {cand.e}*{cand.I}"
        f" + {cand.f}"
    )

    print(
        f"{prefix}{numerator} = "
        f"j*({denominator}) + {cand.q}"
    )


def print_candidate_detail(
    cand: TwoFeatureRationalCandidate,
):

    print_candidate(cand, "  ")

    print()
    print(
        "    inst      j"
        "          numerator"
        "          denominator"
        "          lhs"
        "          rhs"
        "        match"
    )

    print(
        "    "
        + "-" * 100
    )

    for inst in INSTANCES:

        lhs, rhs, numerator, denominator = (
            evaluate_candidate(
                cand,
                inst,
            )
        )

        print(
            f"    {inst.idx:4d}"
            f" {inst.true_j:7d}"
            f" {numerator:18d}"
            f" {denominator:18d}"
            f" {lhs:18d}"
            f" {rhs:18d}"
            f" {str(lhs == rhs):>10s}"
        )


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main():

    print("=" * 120)
    print("EXPERIMENT 467")
    print("=" * 120)
    print()

    print(
        "TWO-FEATURE RATIONAL / EUCLIDEAN j-LAW AUDIT"
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
    # INPUT CONSISTENCY
    # --------------------------------------------------------------------------------------------------------------

    failures = validate_inputs()

    print("INPUT CONSISTENCY")
    print(
        "  "
        + "-" * 80
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
    # FEATURE SCREEN
    # --------------------------------------------------------------------------------------------------------------

    print_feature_screen()

    # --------------------------------------------------------------------------------------------------------------
    # FULL DATA
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "TWO-FEATURE FULL-DATA RATIONAL SEARCH"
    )
    print("=" * 120)
    print()

    full_candidates = full_data_search()

    print(
        "nondegenerate two-feature "
        f"full-data candidates = "
        f"{len(full_candidates)}"
    )

    if full_candidates:

        for cand in full_candidates[:MAX_PRINT]:
            print_candidate(
                cand,
                "  ",
            )

        if len(full_candidates) > MAX_PRINT:
            print(
                f"  ... "
                f"{len(full_candidates) - MAX_PRINT}"
                f" additional candidates suppressed"
            )

    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # LEAVE ONE OUT
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "GENUINE TWO-FEATURE LEAVE-ONE-OUT SEARCH"
    )
    print("=" * 120)
    print()

    loo_candidates, held_out_counts = (
        leave_one_out()
    )

    print(
        f"unique LOO-surviving candidates = "
        f"{len(loo_candidates)}"
    )

    print()

    print(
        "  held-out instance survival counts:"
    )

    for inst in INSTANCES:

        print(
            f"    instance {inst.idx}: "
            f"{held_out_counts[inst.idx]}"
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
                f"    ... "
                f"{len(loo_candidates) - MAX_PRINT}"
                f" additional candidates suppressed"
            )

    else:

        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # DETAILED SURVIVORS
    # --------------------------------------------------------------------------------------------------------------

    if loo_candidates:

        print("=" * 120)
        print(
            "DETAILED LOO SURVIVOR AUDIT"
        )
        print("=" * 120)
        print()

        for index, cand in enumerate(
            loo_candidates[:MAX_PRINT],
            start=1,
        ):

            print(
                f"  CANDIDATE {index}"
            )

            print_candidate_detail(
                cand
            )

            print()

    # --------------------------------------------------------------------------------------------------------------
    # COMPARISON WITH PRIOR EXPERIMENTS
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "CROSS-EXPERIMENT COMPARISON"
    )
    print("=" * 120)
    print()

    print(
        "  Experiment 462:"
    )

    print(
        "    small affine/quadratic laws = 0"
    )

    print(
        "  Experiment 463:"
    )

    print(
        "    small modular laws = 0"
    )

    print(
        "  Experiment 464:"
    )

    print(
        "    restricted quotient/rational laws = 0"
    )

    print(
        "  Experiment 465R2:"
    )

    print(
        "    raw rational candidates = 104272"
    )

    print(
        "  Experiment 466:"
    )

    print(
        "    nondegenerate one-feature rational laws = 0"
    )

    print(
        "  Experiment 467:"
    )

    print(
        f"    nondegenerate two-feature rational laws = "
        f"{len(full_candidates)}"
    )

    print(
        f"    genuine LOO survivors = "
        f"{len(loo_candidates)}"
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "GLOBAL EXPERIMENT 467 SUMMARY"
    )
    print("=" * 120)
    print()

    print(
        f"  instances                         = "
        f"{len(INSTANCES)}"
    )

    print(
        f"  Euclidean observable features    = "
        f"{len(FEATURE_NAMES)}"
    )

    print(
        "  numerator features per candidate = 2"
    )

    print(
        "  denominator features per candidate = 2"
    )

    print(
        f"  coefficient range                 = "
        f"[{COEFF_MIN},{COEFF_MAX}]"
    )

    print(
        f"  offset range                      = "
        f"[{OFFSET_MIN},{OFFSET_MAX}]"
    )

    print(
        f"  full-data candidates              = "
        f"{len(full_candidates)}"
    )

    print(
        f"  genuine LOO candidates            = "
        f"{len(loo_candidates)}"
    )

    print(
        "  constant-feature candidates      = excluded"
    )

    print(
        "  one-feature numerator             = excluded"
    )

    print(
        "  one-feature denominator           = excluded"
    )

    print(
        "  constant denominator              = excluded"
    )

    print(
        "  zero denominator                  = excluded"
    )

    print(
        "  c=0-style degeneracy              = excluded"
    )

    print(
        "  floating point                    = False"
    )

    print(
        "  arbitrary interpolation           = False"
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
        "  Experiment 466 eliminated the large"
    )

    print(
        "  rational candidate population by"
    )

    print(
        "  requiring genuine numerator and"
    )

    print(
        "  denominator variation."
    )

    print()

    print(
        "  Experiment 467 enlarges the tested"
    )

    print(
        "  function family by allowing two"
    )

    print(
        "  independent observable Euclidean"
    )

    print(
        "  features on each side."
    )

    print()

    print(
        "  The tested form is:"
    )

    print(
        "      aF + bG + c"
    )

    print(
        "        ="
    )

    print(
        "      j(dH + eI + f) + q"
    )

    print()

    print(
        "  A negative result only excludes this"
    )

    print(
        "  bounded two-feature rational family."
    )

    print(
        "  It does not prove that no more complex"
    )

    print(
        "  cross-instance relation exists."
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # FINAL STATUS
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "EXPERIMENT 467 FINAL STATUS"
    )
    print("=" * 120)
    print()

    print(
        "  EXACT INTEGER ARITHMETIC        = True"
    )

    print(
        "  TWO-FEATURE RATIONAL TEST       = True"
    )

    print(
        "  NONDEGENERATE FILTERING         = True"
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
        "  CONSTANT FEATURE LAWS            = Excluded"
    )

    print(
        "  ONE-FEATURE LAWS                 = Excluded"
    )

    print(
        "  CONSTANT DENOMINATOR LAWS        = Excluded"
    )

    print(
        "  ZERO DENOMINATOR LAWS            = Excluded"
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
            "    No nondegenerate two-feature"
        )

        print(
            "    rational Euclidean-derived j-law"
        )

        print(
            "    survives genuine leave-one-out"
        )

        print(
            "    validation within the tested"
        )

        print(
            "    coefficient and feature family."
        )

    else:

        print(
            "  CONCLUSION:"
        )

        print(
            "    A nondegenerate two-feature"
        )

        print(
            "    rational Euclidean-derived"
        )

        print(
            "    relation survived leave-one-out."
        )

        print(
            "    This candidate requires an"
        )

        print(
            "    independent structural audit."
        )

    print()
    print(
        "EXPERIMENT 467 FINISHED"
    )
    print("=" * 120)


if __name__ == "__main__":
    main()
