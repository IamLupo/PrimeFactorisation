#!/usr/bin/env python3

"""
========================================================================================================================
EXPERIMENT 468
========================================================================================================================

FAST TWO-FEATURE RATIONAL / EUCLIDEAN j-LAW AUDIT
MEET-IN-THE-MIDDLE EXACT SEARCH

PURPOSE

  Experiment 467 searched the two-feature rational family

      aF + bG + c = j(dH + eI + f) + q

  directly over numerator feature pairs and denominator feature
  pairs. That creates a large Cartesian product.

  Experiment 468 uses the equivalent form

      aF + bG - j(dH + eI + f) = q-c.

  Therefore, across the training instances, the left side must
  be CONSTANT.

  We compare normalized difference vectors instead of comparing
  every numerator candidate against every denominator candidate.

  This changes the search from an explicit pair-of-pairs
  Cartesian product into a hash/meet-in-the-middle search.

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no arbitrary interpolation
  no hidden variables in observable features

NONDEGENERACY

  numerator:
      a != 0
      b != 0
      F and G variable
      aF+bG must vary

  denominator:
      dcoef != 0
      ecoef != 0
      H and I variable
      dcoef*H+ecoef*I+f must vary
      denominator must never vanish

  scalar/sign duplicates are removed

SEARCH

  coefficient range = [-3,3]
  c and q range      = [-6,6]

IMPORTANT OPTIMIZATION

  c disappears from the cross-instance difference signature.

  Once

      Delta = aF0+bG0-j0(dH0+eI0+f)

  is known, we only require

      q-c = Delta.

  Since c,q are both in [-6,6],

      q-c must lie in [-12,12].

  Therefore c and q are NOT enumerated during the expensive
  feature search.

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from itertools import product
from collections import defaultdict


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

COEFF_MIN = -3
COEFF_MAX = 3

OFFSET_MIN = -6
OFFSET_MAX = 6

DIFF_MIN = OFFSET_MIN - OFFSET_MAX
DIFF_MAX = OFFSET_MAX - OFFSET_MIN

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

FEATURE_CACHE = {
    inst.idx: make_features(inst)
    for inst in INSTANCES
}


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class FastCandidate:
    F: str
    G: str

    H: str
    I: str

    a: int
    b: int

    dcoef: int
    ecoef: int
    f: int

    delta: int


# ----------------------------------------------------------------------------------------------------------------------
# VALIDATION
# ----------------------------------------------------------------------------------------------------------------------

def validate_inputs() -> list[str]:

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
                f"instance {inst.idx}: d0 mismatch"
            )

        if x != inst.true_x:
            failures.append(
                f"instance {inst.idx}: x mismatch"
            )

        if K != inst.true_K:
            failures.append(
                f"instance {inst.idx}: K mismatch"
            )

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# VARYING FEATURES
# ----------------------------------------------------------------------------------------------------------------------

def get_varying_features(
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
# NONDEGENERATE FEATURE PAIRS
# ----------------------------------------------------------------------------------------------------------------------

def get_feature_pairs(
    instances: tuple[Instance, ...],
) -> list[tuple[str, str]]:

    varying = get_varying_features(instances)

    return [
        (F, G)
        for F in varying
        for G in varying
        if F != G
    ]


# ----------------------------------------------------------------------------------------------------------------------
# SIGNATURE HELPERS
# ----------------------------------------------------------------------------------------------------------------------

def vector_difference(
    values: list[int],
) -> tuple[int, ...]:

    base = values[0]

    return tuple(
        value - base
        for value in values[1:]
    )


def canonical_sign(
    values: tuple[int, ...],
) -> tuple[int, ...]:

    for value in values:

        if value == 0:
            continue

        if value < 0:
            return tuple(
                -x
                for x in values
            )

        return values

    return values


# ----------------------------------------------------------------------------------------------------------------------
# NUMERATOR SIGNATURES
#
# Numerator:
#
#     P_i = aF_i + bG_i
#
# Since c is constant across all instances, it disappears from
# the difference signature.
# ----------------------------------------------------------------------------------------------------------------------

def build_numerator_signatures(
    instances: tuple[Instance, ...],
    feature_pairs: list[tuple[str, str]],
) -> dict[
    tuple[int, ...],
    list[tuple[str, str, int, int]],
]:

    buckets = defaultdict(list)

    for F, G in feature_pairs:

        F_values = [
            FEATURE_CACHE[inst.idx][F]
            for inst in instances
        ]

        G_values = [
            FEATURE_CACHE[inst.idx][G]
            for inst in instances
        ]

        for a, b in product(
            range(COEFF_MIN, COEFF_MAX + 1),
            repeat=2,
        ):

            if a == 0 or b == 0:
                continue

            values = [
                a * F_values[i]
                + b * G_values[i]
                for i in range(len(instances))
            ]

            if len(set(values)) < 2:
                continue

            sig = canonical_sign(
                vector_difference(values)
            )

            buckets[sig].append(
                (F, G, a, b)
            )

    return buckets


# ----------------------------------------------------------------------------------------------------------------------
# DENOMINATOR SIGNATURES
#
# Denominator:
#
#     D_i = j_i(dH_i + eI_i + f)
#
# The constant f DOES NOT disappear because it is multiplied by j_i.
# ----------------------------------------------------------------------------------------------------------------------

def build_denominator_signatures(
    instances: tuple[Instance, ...],
    feature_pairs: list[tuple[str, str]],
) -> dict[
    tuple[int, ...],
    list[tuple[str, str, int, int, int]],
]:

    buckets = defaultdict(list)

    j_values = [
        inst.true_j
        for inst in instances
    ]

    for H, I in feature_pairs:

        H_values = [
            FEATURE_CACHE[inst.idx][H]
            for inst in instances
        ]

        I_values = [
            FEATURE_CACHE[inst.idx][I]
            for inst in instances
        ]

        for dcoef, ecoef in product(
            range(COEFF_MIN, COEFF_MAX + 1),
            repeat=2,
        ):

            if dcoef == 0 or ecoef == 0:
                continue

            for f in range(
                OFFSET_MIN,
                OFFSET_MAX + 1,
            ):

                denom_values = [
                    dcoef * H_values[i]
                    + ecoef * I_values[i]
                    + f
                    for i in range(len(instances))
                ]

                if len(set(denom_values)) < 2:
                    continue

                if any(
                    value == 0
                    for value in denom_values
                ):
                    continue

                values = [
                    j_values[i] * denom_values[i]
                    for i in range(len(instances))
                ]

                sig = vector_difference(values)

                buckets[sig].append(
                    (
                        H,
                        I,
                        dcoef,
                        ecoef,
                        f,
                    )
                )

    return buckets


# ----------------------------------------------------------------------------------------------------------------------
# MATCHING
# ----------------------------------------------------------------------------------------------------------------------

def candidate_matches(
    numerator_entry,
    denominator_entry,
    instances: tuple[Instance, ...],
) -> FastCandidate | None:

    F, G, a, b = numerator_entry

    H, I, dcoef, ecoef, f = denominator_entry

    numerator0 = (
        a * FEATURE_CACHE[instances[0].idx][F]
        + b * FEATURE_CACHE[instances[0].idx][G]
    )

    denominator0 = (
        dcoef * FEATURE_CACHE[instances[0].idx][H]
        + ecoef * FEATURE_CACHE[instances[0].idx][I]
        + f
    )

    delta = (
        numerator0
        - instances[0].true_j * denominator0
    )

    if not (
        DIFF_MIN
        <= delta
        <= DIFF_MAX
    ):
        return None

    return FastCandidate(
        F=F,
        G=G,
        H=H,
        I=I,
        a=a,
        b=b,
        dcoef=dcoef,
        ecoef=ecoef,
        f=f,
        delta=delta,
    )


# ----------------------------------------------------------------------------------------------------------------------
# VALIDATE CANDIDATE
# ----------------------------------------------------------------------------------------------------------------------

def candidate_is_valid(
    cand: FastCandidate,
    instances: tuple[Instance, ...],
) -> bool:

    lhs_constants = []

    for inst in instances:

        table = FEATURE_CACHE[inst.idx]

        numerator = (
            cand.a * table[cand.F]
            + cand.b * table[cand.G]
        )

        denominator = (
            cand.dcoef * table[cand.H]
            + cand.ecoef * table[cand.I]
            + cand.f
        )

        if denominator == 0:
            return False

        value = (
            numerator
            - inst.true_j * denominator
        )

        if value != cand.delta:
            return False

        lhs_constants.append(value)

    if len(set(lhs_constants)) != 1:
        return False

    return True


# ----------------------------------------------------------------------------------------------------------------------
# CHECK THAT c,q CAN REALIZE delta
# ----------------------------------------------------------------------------------------------------------------------

def possible_offsets(
    delta: int,
) -> list[tuple[int, int]]:

    result = []

    for c in range(
        OFFSET_MIN,
        OFFSET_MAX + 1,
    ):

        q = delta + c

        if (
            OFFSET_MIN
            <= q
            <= OFFSET_MAX
        ):

            result.append(
                (c, q)
            )

    return result


# ----------------------------------------------------------------------------------------------------------------------
# CANONICALIZE
# ----------------------------------------------------------------------------------------------------------------------

def canonical_candidate_key(
    cand: FastCandidate,
    c: int,
    q: int,
):
    values = (
        cand.a,
        cand.b,
        cand.dcoef,
        cand.ecoef,
        cand.f,
        cand.delta,
        c,
        q,
    )

    first_nonzero = next(
        (
            x
            for x in values
            if x != 0
        ),
        0,
    )

    if first_nonzero < 0:
        values = tuple(
            -x
            for x in values
        )

    return (
        cand.F,
        cand.G,
        cand.H,
        cand.I,
        values,
    )


# ----------------------------------------------------------------------------------------------------------------------
# FAST SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def search_fast(
    instances: tuple[Instance, ...],
) -> list[tuple[FastCandidate, int, int]]:

    feature_pairs = get_feature_pairs(
        instances
    )

    numerator_buckets = (
        build_numerator_signatures(
            instances,
            feature_pairs,
        )
    )

    denominator_buckets = (
        build_denominator_signatures(
            instances,
            feature_pairs,
        )
    )

    matches = []

    for sig, numerator_entries in (
        numerator_buckets.items()
    ):

        denominators = (
            denominator_buckets.get(sig)
        )

        if not denominators:
            # A sign-canonical numerator signature
            # needs a correspondingly canonical denominator
            # signature. Try the opposite orientation too.
            opposite = tuple(
                -x
                for x in sig
            )
            denominators = (
                denominator_buckets.get(
                    opposite
                )
            )

            if not denominators:
                continue

        for numerator_entry in (
            numerator_entries
        ):

            for denominator_entry in (
                denominators
            ):

                cand = candidate_matches(
                    numerator_entry,
                    denominator_entry,
                    instances,
                )

                if cand is None:
                    continue

                if not candidate_is_valid(
                    cand,
                    instances,
                ):
                    continue

                for c, q in possible_offsets(
                    cand.delta
                ):

                    matches.append(
                        (
                            cand,
                            c,
                            q,
                        )
                    )

    # Deduplicate.
    unique = {}

    for cand, c, q in matches:

        key = canonical_candidate_key(
            cand,
            c,
            q,
        )

        unique[key] = (
            cand,
            c,
            q,
        )

    return list(
        unique.values()
    )


# ----------------------------------------------------------------------------------------------------------------------
# EXACT FULL FORM
# ----------------------------------------------------------------------------------------------------------------------

def evaluate_full(
    cand: FastCandidate,
    c: int,
    q: int,
    inst: Instance,
) -> tuple[int, int, int, int]:

    table = FEATURE_CACHE[inst.idx]

    numerator = (
        cand.a * table[cand.F]
        + cand.b * table[cand.G]
        + c
    )

    denominator = (
        cand.dcoef * table[cand.H]
        + cand.ecoef * table[cand.I]
        + cand.f
    )

    lhs = numerator

    rhs = (
        inst.true_j * denominator
        + q
    )

    return (
        numerator,
        denominator,
        lhs,
        rhs,
    )


# ----------------------------------------------------------------------------------------------------------------------
# FULL VALIDATION
# ----------------------------------------------------------------------------------------------------------------------

def exact_full_validation(
    candidate,
) -> bool:

    cand, c, q = candidate

    for inst in INSTANCES:

        numerator, denominator, lhs, rhs = (
            evaluate_full(
                cand,
                c,
                q,
                inst,
            )
        )

        if denominator == 0:
            return False

        if lhs != rhs:
            return False

    return True


# ----------------------------------------------------------------------------------------------------------------------
# LEAVE ONE OUT
# ----------------------------------------------------------------------------------------------------------------------

def genuine_loo_search():

    survivors = []
    survival_counts = {
        inst.idx: 0
        for inst in INSTANCES
    }

    for held_out in INSTANCES:

        training = tuple(
            inst
            for inst in INSTANCES
            if inst.idx != held_out.idx
        )

        candidates = search_fast(
            training
        )

        local = []

        for candidate in candidates:

            cand, c, q = candidate

            numerator, denominator, lhs, rhs = (
                evaluate_full(
                    cand,
                    c,
                    q,
                    held_out,
                )
            )

            if denominator == 0:
                continue

            if lhs == rhs:

                local.append(
                    candidate
                )

        survival_counts[
            held_out.idx
        ] = len(local)

        survivors.extend(local)

    # Unique survivors.
    unique = {}

    for candidate in survivors:

        cand, c, q = candidate

        key = canonical_candidate_key(
            cand,
            c,
            q,
        )

        unique[key] = candidate

    return (
        list(unique.values()),
        survival_counts,
    )


# ----------------------------------------------------------------------------------------------------------------------
# FEATURE STATISTICS
# ----------------------------------------------------------------------------------------------------------------------

def print_feature_screen():

    print("=" * 120)
    print(
        "NONDEGENERATE EUCLIDEAN FEATURE SCREEN"
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
            f" {status}"
        )

    print()


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE PRINT
# ----------------------------------------------------------------------------------------------------------------------

def print_candidate(
    candidate,
    prefix: str = "",
):

    cand, c, q = candidate

    numerator = (
        f"{cand.a}*{cand.F}"
        f" + {cand.b}*{cand.G}"
        f" + {c}"
    )

    denominator = (
        f"{cand.dcoef}*{cand.H}"
        f" + {cand.ecoef}*{cand.I}"
        f" + {cand.f}"
    )

    print(
        f"{prefix}{numerator}"
        f" = j*({denominator})"
        f" + {q}"
    )


def print_candidate_details(
    candidate,
):

    cand, c, q = candidate

    print_candidate(
        candidate,
        "  ",
    )

    print()

    print(
        "    idx"
        "       j"
        "        numerator"
        "        denominator"
        "        lhs"
        "        rhs"
        "      match"
    )

    print(
        "    "
        + "-" * 105
    )

    for inst in INSTANCES:

        numerator, denominator, lhs, rhs = (
            evaluate_full(
                cand,
                c,
                q,
                inst,
            )
        )

        print(
            f"    {inst.idx:3d}"
            f" {inst.true_j:8d}"
            f" {numerator:18d}"
            f" {denominator:18d}"
            f" {lhs:18d}"
            f" {rhs:18d}"
            f" {str(lhs == rhs):>8s}"
        )


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main():

    print("=" * 120)
    print(
        "EXPERIMENT 468"
    )
    print("=" * 120)
    print()

    print(
        "FAST TWO-FEATURE RATIONAL / EUCLIDEAN j-LAW AUDIT"
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

    print(
        "SEARCH METHOD     = meet-in-the-middle signatures"
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # INPUT
    # --------------------------------------------------------------------------------------------------------------

    failures = validate_inputs()

    print(
        "INPUT CONSISTENCY"
    )

    print(
        "  "
        + "-" * 80
    )

    print(
        f"  consistency failures = "
        f"{len(failures)}"
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
        "FAST TWO-FEATURE FULL-DATA SEARCH"
    )
    print("=" * 120)
    print()

    full_candidates = search_fast(
        INSTANCES
    )

    print(
        f"full-data candidates = "
        f"{len(full_candidates)}"
    )

    if full_candidates:

        for candidate in full_candidates[
            :MAX_PRINT
        ]:

            print_candidate(
                candidate,
                "  ",
            )

        if len(full_candidates) > MAX_PRINT:

            print(
                f"  ... "
                f"{len(full_candidates) - MAX_PRINT}"
                f" additional candidates suppressed"
            )

    else:

        print(
            "  none"
        )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # GENUINE LOO
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "GENUINE FAST LEAVE-ONE-OUT SEARCH"
    )
    print("=" * 120)
    print()

    loo_candidates, survival_counts = (
        genuine_loo_search()
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
            f"{survival_counts[inst.idx]}"
        )

    print()

    if loo_candidates:

        for candidate in loo_candidates[
            :MAX_PRINT
        ]:

            print_candidate(
                candidate,
                "  ",
            )

        if len(loo_candidates) > MAX_PRINT:

            print(
                f"  ... "
                f"{len(loo_candidates) - MAX_PRINT}"
                f" additional candidates suppressed"
            )

    else:

        print(
            "  none"
        )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # OPTIONAL DETAIL
    # --------------------------------------------------------------------------------------------------------------

    if loo_candidates:

        print("=" * 120)
        print(
            "LOO SURVIVOR DETAIL"
        )
        print("=" * 120)
        print()

        for index, candidate in enumerate(
            loo_candidates[
                :MAX_PRINT
            ],
            start=1,
        ):

            print(
                f"  CANDIDATE {index}"
            )

            print_candidate_details(
                candidate
            )

            print()

    # --------------------------------------------------------------------------------------------------------------
    # COMPLEXITY EXPLANATION
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "SEARCH COMPLEXITY COMPARISON"
    )
    print("=" * 120)
    print()

    print(
        "  Previous direct approach:"
    )

    print(
        "    numerator feature pair"
        " x denominator feature pair"
        " x coefficients"
    )

    print(
        "    This creates a large Cartesian product."
    )

    print()

    print(
        "  Current approach:"
    )

    print(
        "    1. enumerate numerator signatures"
    )

    print(
        "    2. enumerate denominator signatures"
    )

    print(
        "    3. hash equal difference signatures"
    )

    print(
        "    4. recover q-c as a small constant"
    )

    print(
        "  No numerator/denominator Cartesian product is built."
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "GLOBAL EXPERIMENT 468 SUMMARY"
    )
    print("=" * 120)
    print()

    print(
        f"  instances                         = "
        f"{len(INSTANCES)}"
    )

    print(
        f"  Euclidean features               = "
        f"{len(FEATURE_NAMES)}"
    )

    variable_count = len(
        get_varying_features(
            INSTANCES
        )
    )

    pair_count = variable_count * (
        variable_count - 1
    )

    print(
        f"  variable features                = "
        f"{variable_count}"
    )

    print(
        f"  ordered feature pairs            = "
        f"{pair_count}"
    )

    print(
        f"  coefficient range                = "
        f"[{COEFF_MIN},{COEFF_MAX}]"
    )

    print(
        f"  offset range                     = "
        f"[{OFFSET_MIN},{OFFSET_MAX}]"
    )

    print(
        f"  full-data candidates             = "
        f"{len(full_candidates)}"
    )

    print(
        f"  genuine LOO candidates           = "
        f"{len(loo_candidates)}"
    )

    print(
        "  constant features                = excluded"
    )

    print(
        "  zero denominator                 = excluded"
    )

    print(
        "  one-feature numerator            = excluded"
    )

    print(
        "  one-feature denominator          = excluded"
    )

    print(
        "  arbitrary interpolation           = False"
    )

    print(
        "  floating point                    = False"
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
        "  This experiment searches the same"
    )

    print(
        "  nondegenerate two-feature rational"
    )

    print(
        "  family as Experiment 467."
    )

    print()

    print(
        "  The mathematical search condition is:"
    )

    print(
        "      aF + bG - j(dH + eI + f)"
    )

    print(
        "  must be constant across the training"
    )

    print(
        "  instances."
    )

    print()

    print(
        "  Because c and q enter only through"
    )

    print(
        "      q-c,"
    )

    print(
        "  they are recovered after the expensive"
    )

    print(
        "  signature match rather than enumerated"
    )

    print(
        "  during feature-pair matching."
    )

    print()

    print(
        "  A negative result still excludes only"
    )

    print(
        "  this bounded feature/coefficient family."
    )

    print(
        "  It does not establish impossibility of"
    )

    print(
        "  arbitrary cross-instance reconstruction."
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # FINAL STATUS
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "EXPERIMENT 468 FINAL STATUS"
    )
    print("=" * 120)
    print()

    print(
        "  EXACT INTEGER ARITHMETIC        = True"
    )

    print(
        "  TWO-FEATURE RATIONAL FAMILY     = True"
    )

    print(
        "  MEET-IN-THE-MIDDLE SEARCH       = True"
    )

    print(
        "  NONDEGENERATE FILTER             = True"
    )

    print(
        "  GENUINE LEAVE-ONE-OUT            = True"
    )

    print(
        f"  FULL-DATA CANDIDATES              = "
        f"{len(full_candidates)}"
    )

    print(
        f"  LOO CANDIDATES                    = "
        f"{len(loo_candidates)}"
    )

    print(
        "  CONSTANT FEATURES                = Excluded"
    )

    print(
        "  ZERO DENOMINATOR                 = Excluded"
    )

    print(
        "  ONE-FEATURE NUMERATOR            = Excluded"
    )

    print(
        "  ONE-FEATURE DENOMINATOR          = Excluded"
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
            "    rational Euclidean-derived"
        )

        print(
            "    j-law survives genuine"
        )

        print(
            "    leave-one-out validation"
        )

        print(
            "    in the tested bounded family."
        )

    else:

        print(
            "  CONCLUSION:"
        )

        print(
            "    At least one nondegenerate"
        )

        print(
            "    two-feature rational relation"
        )

        print(
            "    survives leave-one-out."
        )

        print(
            "    It requires an independent"
        )

        print(
            "    structural audit."
        )

    print()

    print(
        "EXPERIMENT 468 FINISHED"
    )
    print("=" * 120)


if __name__ == "__main__":
    main()
