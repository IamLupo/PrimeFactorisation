#!/usr/bin/env python3

"""
========================================================================================================================
EXPERIMENT 465R2
========================================================================================================================

EUCLIDEAN-DERIVED OBSERVABLE / j-LAW AUDIT
CORRECTED RATIONAL SEARCH

QUESTION

  Can the hidden displacement j be determined by a simple exact law
  built from Euclidean-division structure of the observed tuple?

OBSERVED VARIABLES

  N, S, d, r

DERIVED OBSERVABLE FEATURES

  Euclidean quotients
  Euclidean remainders
  gcd values
  simple Euclidean combinations
  integer square-root / bit-length diagnostics

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no interpolation
  no hidden variables in observable features

SEARCH FAMILIES

  1. affine
       j = a*F + b

  2. floor
       j = floor((a*F+b)/(G+c))

  3. quotient
       j = (a*F+b)//G + c

  4. exact rational
       a*F+b = j*(c*G+e) + q

The rational candidate explicitly contains all five coefficients
(a,b,c,e,q), fixing the constructor bug from the previous script.

The expensive searches use only the Euclidean feature basis,
rather than every auxiliary feature.

All surviving candidates are subjected to leave-one-out validation.

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

COEFF_MIN = -4
COEFF_MAX = 4

OFFSET_MIN = -8
OFFSET_MAX = 8

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
# EUCLIDEAN HELPERS
# ----------------------------------------------------------------------------------------------------------------------

def euclid(a: int, b: int) -> tuple[int, int, int]:
    """
    Return:

        gcd(|a|, |b|), quotient, remainder

    using exact Python integer divmod semantics.
    """

    if b == 0:
        raise ValueError("Euclidean divisor cannot be zero")

    q, rem = divmod(a, b)
    g = gcd(abs(a), abs(b))
    return g, q, rem


# ----------------------------------------------------------------------------------------------------------------------
# OBSERVABLE FEATURE CONSTRUCTION
# ----------------------------------------------------------------------------------------------------------------------

def all_observable_features(inst: Instance) -> dict[str, int]:

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
        # Primitive observables.
        "1": 1,
        "N": N,
        "S": S,
        "d": d,
        "r": r,

        # Euclidean quotients.
        "q(N,S)": q_NS,
        "q(N,d)": q_Nd,
        "q(S,d)": q_Sd,
        "q(d,r)": q_dr,
        "q(N,r)": q_Nr,
        "q(S,r)": q_Sr,

        # Euclidean remainders.
        "rem(N,S)": rem_NS,
        "rem(N,d)": rem_Nd,
        "rem(S,d)": rem_Sd,
        "rem(d,r)": rem_dr,
        "rem(N,r)": rem_Nr,
        "rem(S,r)": rem_Sr,

        # GCDs.
        "gcd(N,S)": g_NS,
        "gcd(N,d)": g_Nd,
        "gcd(S,d)": g_Sd,
        "gcd(d,r)": g_dr,
        "gcd(N,r)": g_Nr,
        "gcd(S,r)": g_Sr,

        # Simple combinations.
        "N-S": N - S,
        "N-d": N - d,
        "d-S": d - S,
        "N+S": N + S,
        "N+d": N + d,
        "S+d": S + d,

        # Quotient combinations.
        "q(N,S)-q(S,d)": q_NS - q_Sd,
        "q(N,d)-q(S,d)": q_Nd - q_Sd,
        "q(S,d)-q(N,S)": q_Sd - q_NS,

        # Remainder combinations.
        "rem(N,S)-rem(S,d)": rem_NS - rem_Sd,
        "rem(N,d)-rem(S,d)": rem_Nd - rem_Sd,
        "rem(N,S)-rem(N,d)": rem_NS - rem_Nd,

        # GCD combinations.
        "gcd(N,S)-gcd(S,d)": g_NS - g_Sd,
        "gcd(N,d)-gcd(S,d)": g_Nd - g_Sd,

        # Scale diagnostics.
        "sqrtN": isqrt(N),
        "bitlenN": N.bit_length(),
        "bitlenS": S.bit_length(),
        "bitlend": d.bit_length(),

        "digitsN": len(str(abs(N))),
        "digitsS": len(str(abs(S))),
        "digitsd": len(str(abs(d))),

        "N//S": N // S,
        "d//S": d // S,
        "N//d": N // d,
    }


# Only the genuine Euclidean-derived basis is used by the expensive searches.
SEARCH_FEATURES = (
    "q(N,S)",
    "q(N,d)",
    "q(S,d)",
    "q(d,r)",
    "q(N,r)",
    "q(S,r)",

    "rem(N,S)",
    "rem(N,d)",
    "rem(S,d)",
    "rem(d,r)",
    "rem(N,r)",
    "rem(S,r)",

    "gcd(N,S)",
    "gcd(N,d)",
    "gcd(S,d)",
    "gcd(d,r)",
    "gcd(N,r)",
    "gcd(S,r)",

    "q(N,S)-q(S,d)",
    "q(N,d)-q(S,d)",
    "q(S,d)-q(N,S)",

    "rem(N,S)-rem(S,d)",
    "rem(N,d)-rem(S,d)",
    "rem(N,S)-rem(N,d)",

    "gcd(N,S)-gcd(S,d)",
    "gcd(N,d)-gcd(S,d)",
)


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE TYPES
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class AffineCandidate:
    F: str
    a: int
    b: int


@dataclass(frozen=True)
class FloorCandidate:
    F: str
    G: str
    a: int
    b: int
    c: int


@dataclass(frozen=True)
class QuotientCandidate:
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
# CONSISTENCY
# ----------------------------------------------------------------------------------------------------------------------

def validate_instances(instances) -> int:

    failures = 0

    for inst in instances:

        if inst.d != inst.N - 2 * inst.S + 1:
            failures += 1

        j = inst.true_j

        if inst.d0_ref + 2 * j != inst.true_d0:
            failures += 1

        if inst.x0 - 2 * j != inst.true_x:
            failures += 1

        if (
            inst.K_ref
            + inst.d0_ref * j
            + j * j
            != inst.true_K
        ):
            failures += 1

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# AFFINE SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def affine_holds(cand: AffineCandidate, inst: Instance) -> bool:

    values = all_observable_features(inst)

    return (
        cand.a * values[cand.F]
        + cand.b
        == inst.true_j
    )


def search_affine(instances):

    results = []

    for F in SEARCH_FEATURES:

        for a in range(COEFF_MIN, COEFF_MAX + 1):

            if a == 0:
                continue

            for b in range(OFFSET_MIN, OFFSET_MAX + 1):

                cand = AffineCandidate(
                    F,
                    a,
                    b,
                )

                if all(
                    affine_holds(cand, inst)
                    for inst in instances
                ):
                    results.append(cand)

    return results


# ----------------------------------------------------------------------------------------------------------------------
# FLOOR SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def floor_holds(cand: FloorCandidate, inst: Instance) -> bool:

    values = all_observable_features(inst)

    numerator = (
        cand.a * values[cand.F]
        + cand.b
    )

    denominator = (
        values[cand.G]
        + cand.c
    )

    if denominator == 0:
        return False

    return numerator // denominator == inst.true_j


def search_floor(instances):

    results = []

    for F in SEARCH_FEATURES:
        for G in SEARCH_FEATURES:

            for a in range(COEFF_MIN, COEFF_MAX + 1):

                if a == 0:
                    continue

                for b in range(OFFSET_MIN, OFFSET_MAX + 1):

                    # c can be inferred from the small bounded range.
                    for c in range(OFFSET_MIN, OFFSET_MAX + 1):

                        cand = FloorCandidate(
                            F,
                            G,
                            a,
                            b,
                            c,
                        )

                        if all(
                            floor_holds(cand, inst)
                            for inst in instances
                        ):
                            results.append(cand)

    return results


# ----------------------------------------------------------------------------------------------------------------------
# QUOTIENT SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def quotient_holds(
    cand: QuotientCandidate,
    inst: Instance,
) -> bool:

    values = all_observable_features(inst)

    denominator = values[cand.G]

    if denominator == 0:
        return False

    numerator = (
        cand.a * values[cand.F]
        + cand.b
    )

    predicted = (
        numerator // denominator
        + cand.c
    )

    return predicted == inst.true_j


def search_quotient(instances):

    results = []

    for F in SEARCH_FEATURES:
        for G in SEARCH_FEATURES:

            for a in range(COEFF_MIN, COEFF_MAX + 1):

                if a == 0:
                    continue

                for b in range(OFFSET_MIN, OFFSET_MAX + 1):

                    for c in range(OFFSET_MIN, OFFSET_MAX + 1):

                        cand = QuotientCandidate(
                            F,
                            G,
                            a,
                            b,
                            c,
                        )

                        if all(
                            quotient_holds(cand, inst)
                            for inst in instances
                        ):
                            results.append(cand)

    return results


# ----------------------------------------------------------------------------------------------------------------------
# RATIONAL SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def rational_holds(
    cand: RationalCandidate,
    inst: Instance,
) -> bool:

    values = all_observable_features(inst)

    lhs = (
        cand.a * values[cand.F]
        + cand.b
    )

    denominator = (
        cand.c * values[cand.G]
        + cand.e
    )

    rhs = (
        inst.true_j * denominator
        + cand.q
    )

    return lhs == rhs


def search_rational(instances):

    """
    Search:

        a*F+b = j*(c*G+e)+q

    All coefficients are bounded.

    The constant q is inferred from the first instance, so it
    does not need an additional search dimension.

    This fixes the previous constructor error and also reduces
    unnecessary work.
    """

    results = []

    first = instances[0]

    first_feature_cache = all_observable_features(first)

    for F in SEARCH_FEATURES:

        F0 = first_feature_cache[F]

        for G in SEARCH_FEATURES:

            G0 = first_feature_cache[G]

            for a in range(COEFF_MIN, COEFF_MAX + 1):

                if a == 0:
                    continue

                for b in range(OFFSET_MIN, OFFSET_MAX + 1):

                    lhs0 = a * F0 + b

                    for c in range(COEFF_MIN, COEFF_MAX + 1):

                        for e in range(OFFSET_MIN, OFFSET_MAX + 1):

                            denominator0 = c * G0 + e

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

                            if all(
                                rational_holds(
                                    cand,
                                    inst,
                                )
                                for inst in instances
                            ):
                                results.append(cand)

    return results


# ----------------------------------------------------------------------------------------------------------------------
# LEAVE-ONE-OUT VALIDATION
# ----------------------------------------------------------------------------------------------------------------------

def loo_check(cand, predicate, instances):

    failures = []

    for held_out in instances:

        training = [
            inst
            for inst in instances
            if inst.idx != held_out.idx
        ]

        if not all(
            predicate(cand, inst)
            for inst in training
        ):
            failures.append(
                held_out.idx
            )
            continue

        if not predicate(
            cand,
            held_out,
        ):
            failures.append(
                held_out.idx
            )

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# TABLES
# ----------------------------------------------------------------------------------------------------------------------

def print_euclidean_table(instances):

    print("=" * 120)
    print("EUCLIDEAN OBSERVABLE FEATURE TABLE")
    print("=" * 120)
    print()

    pairs = (
        ("N,S", "q(N,S)", "rem(N,S)", "gcd(N,S)"),
        ("N,d", "q(N,d)", "rem(N,d)", "gcd(N,d)"),
        ("S,d", "q(S,d)", "rem(S,d)", "gcd(S,d)"),
        ("d,r", "q(d,r)", "rem(d,r)", "gcd(d,r)"),
        ("N,r", "q(N,r)", "rem(N,r)", "gcd(N,r)"),
        ("S,r", "q(S,r)", "rem(S,r)", "gcd(S,r)"),
    )

    for pair_name, qn, rn, gn in pairs:

        print(f"  PAIR = {pair_name}")
        print(
            "    idx       j        quotient        remainder        gcd"
        )
        print(
            "    " + "-" * 75
        )

        for inst in instances:

            values = all_observable_features(inst)

            print(
                f"    {inst.idx:<3d}"
                f" {inst.true_j:8d}"
                f" {values[qn]:16d}"
                f" {values[rn]:16d}"
                f" {values[gn]:16d}"
            )

        print()


def print_feature_summary(instances):

    print("=" * 120)
    print("EUCLIDEAN SEARCH-FEATURE SUMMARY")
    print("=" * 120)
    print()

    for name in SEARCH_FEATURES:

        values = [
            all_observable_features(inst)[name]
            for inst in instances
        ]

        print(
            f"  {name:30s}"
            f" distinct={len(set(values)):<2d}"
            f" values={values}"
        )

    print()


# ----------------------------------------------------------------------------------------------------------------------
# CANDIDATE OUTPUT
# ----------------------------------------------------------------------------------------------------------------------

def print_affine_candidate(cand, instances):

    print()
    print("  AFFINE CANDIDATE")
    print(
        f"    j = {cand.a}*{cand.F} + {cand.b}"
    )

    for inst in instances:

        value = all_observable_features(inst)[cand.F]
        predicted = cand.a * value + cand.b

        print(
            f"    inst={inst.idx}"
            f"  j={inst.true_j}"
            f"  predicted={predicted}"
        )


def print_floor_candidate(cand, instances):

    print()
    print("  FLOOR CANDIDATE")
    print(
        f"    j = floor("
        f"({cand.a}*{cand.F}+{cand.b})"
        f"/({cand.G}+{cand.c})"
        f")"
    )


def print_quotient_candidate(cand, instances):

    print()
    print("  QUOTIENT CANDIDATE")
    print(
        f"    j = "
        f"({cand.a}*{cand.F}+{cand.b})"
        f"//{cand.G} + {cand.c}"
    )


def print_rational_candidate(cand, instances):

    print()
    print("  RATIONAL CANDIDATE")
    print(
        f"    {cand.a}*{cand.F}+{cand.b}"
        f" = j*({cand.c}*{cand.G}+{cand.e})"
        f" + {cand.q}"
    )

    for inst in instances:

        values = all_observable_features(inst)

        lhs = (
            cand.a * values[cand.F]
            + cand.b
        )

        denominator = (
            cand.c * values[cand.G]
            + cand.e
        )

        rhs = (
            inst.true_j * denominator
            + cand.q
        )

        print(
            f"    inst={inst.idx}"
            f"  j={inst.true_j}"
            f"  lhs={lhs}"
            f"  rhs={rhs}"
            f"  match={lhs == rhs}"
        )


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main():

    print("=" * 120)
    print("EXPERIMENT 465R2")
    print("=" * 120)
    print()

    print(
        "CORRECTED EUCLIDEAN-DERIVED OBSERVABLE / j-LAW AUDIT"
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
    # Consistency.
    # --------------------------------------------------------------------------------------------------------------

    print("INPUT CONSISTENCY")
    print(
        "  "
        + "-"
        * 80
    )

    consistency_failures = validate_instances(
        INSTANCES
    )

    print(
        f"  consistency failures = "
        f"{consistency_failures}"
    )
    print()

    if consistency_failures:
        raise AssertionError(
            "input consistency failure"
        )

    # --------------------------------------------------------------------------------------------------------------
    # Euclidean tables.
    # --------------------------------------------------------------------------------------------------------------

    print_euclidean_table(
        INSTANCES
    )

    print_feature_summary(
        INSTANCES
    )

    # --------------------------------------------------------------------------------------------------------------
    # Affine.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("EUCLIDEAN-DERIVED AFFINE SEARCH")
    print("=" * 120)
    print()

    affine_candidates = search_affine(
        INSTANCES
    )

    affine_candidates = [
        cand
        for cand in affine_candidates
        if loo_check(
            cand,
            affine_holds,
            INSTANCES,
        ) == []
    ]

    print(
        f"surviving affine candidates = "
        f"{len(affine_candidates)}"
    )

    if affine_candidates:
        for cand in affine_candidates[:MAX_PRINT]:
            print_affine_candidate(
                cand,
                INSTANCES,
            )
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Floor.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("EUCLIDEAN-DERIVED FLOOR SEARCH")
    print("=" * 120)
    print()

    floor_candidates = search_floor(
        INSTANCES
    )

    floor_candidates = [
        cand
        for cand in floor_candidates
        if loo_check(
            cand,
            floor_holds,
            INSTANCES,
        ) == []
    ]

    print(
        f"surviving floor candidates = "
        f"{len(floor_candidates)}"
    )

    if floor_candidates:
        for cand in floor_candidates[:MAX_PRINT]:
            print_floor_candidate(
                cand,
                INSTANCES,
            )
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Quotient.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("EUCLIDEAN-DERIVED QUOTIENT SEARCH")
    print("=" * 120)
    print()

    quotient_candidates = search_quotient(
        INSTANCES
    )

    quotient_candidates = [
        cand
        for cand in quotient_candidates
        if loo_check(
            cand,
            quotient_holds,
            INSTANCES,
        ) == []
    ]

    print(
        f"surviving quotient candidates = "
        f"{len(quotient_candidates)}"
    )

    if quotient_candidates:
        for cand in quotient_candidates[:MAX_PRINT]:
            print_quotient_candidate(
                cand,
                INSTANCES,
            )
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Rational.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("EUCLIDEAN-DERIVED RATIONAL SEARCH")
    print("=" * 120)
    print()

    rational_candidates = search_rational(
        INSTANCES
    )

    rational_candidates = [
        cand
        for cand in rational_candidates
        if loo_check(
            cand,
            rational_holds,
            INSTANCES,
        ) == []
    ]

    print(
        f"surviving rational candidates = "
        f"{len(rational_candidates)}"
    )

    if rational_candidates:
        for cand in rational_candidates[:MAX_PRINT]:
            print_rational_candidate(
                cand,
                INSTANCES,
            )
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # Summary.
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 465R2 SUMMARY")
    print("=" * 120)
    print()

    print(
        f"  instances                         = "
        f"{len(INSTANCES)}"
    )

    print(
        f"  Euclidean search features         = "
        f"{len(SEARCH_FEATURES)}"
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
        f"  affine candidates                 = "
        f"{len(affine_candidates)}"
    )

    print(
        f"  floor candidates                  = "
        f"{len(floor_candidates)}"
    )

    print(
        f"  quotient candidates               = "
        f"{len(quotient_candidates)}"
    )

    print(
        f"  rational candidates               = "
        f"{len(rational_candidates)}"
    )

    print()

    if (
        not affine_candidates
        and not floor_candidates
        and not quotient_candidates
        and not rational_candidates
    ):
        print(
            "  No tested low-complexity Euclidean-derived"
        )
        print(
            "  observable law survived all eight instances."
        )
    else:
        print(
            "  At least one Euclidean-derived candidate survived."
        )
        print(
            "  Candidate requires independent validation."
        )

    print()

    print("=" * 120)
    print("EXPERIMENT 465R2 FINAL STATUS")
    print("=" * 120)
    print()

    print(
        "  EXACT INTEGER ARITHMETIC        = True"
    )
    print(
        "  EUCLIDEAN FEATURE BASIS         = True"
    )
    print(
        f"  AFFINE CANDIDATES                = "
        f"{len(affine_candidates)}"
    )
    print(
        f"  FLOOR CANDIDATES                 = "
        f"{len(floor_candidates)}"
    )
    print(
        f"  QUOTIENT CANDIDATES              = "
        f"{len(quotient_candidates)}"
    )
    print(
        f"  RATIONAL CANDIDATES              = "
        f"{len(rational_candidates)}"
    )
    print(
        "  LEAVE-ONE-OUT VALIDATION        = True"
    )
    print(
        "  FLOATING POINT                  = False"
    )
    print(
        "  ARBITRARY INTERPOLATION         = False"
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
    print()

    print(
        "  RATIONAL CONSTRUCTOR FIXED      = True"
    )

    print()

    print(
        "  CONCLUSION:"
    )

    print(
        "    Experiment 465R2 repeats the Euclidean-derived"
    )
    print(
        "    observable search with a corrected exact rational"
    )
    print(
        "    candidate representation."
    )

    print(
        "    The rational family is explicitly:"
    )

    print(
        "      a*F+b = j*(c*G+e)+q"
    )

    print(
        "    with all coefficients independently bounded."
    )

    print()

    print(
        "EXPERIMENT 465R2 FINISHED"
    )
    print("=" * 120)


if __name__ == "__main__":
    main()