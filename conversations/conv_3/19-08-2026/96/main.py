#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 473
========================================================================================================================

SPARSE MIXED-QUADRATIC / EUCLIDEAN j-LAW AUDIT

QUESTION

  After failure of:

      Experiment 471:
          a*F + b*G + c = j*(d*H + e)

      Experiment 472:
          a*F^2 + b*G + c = j*(d*H + e)

  does a mixed quadratic observable interaction survive?

SEARCH FORM

      a*F*G + b*H + c = j*(d*Q + e)

where F,G,H,Q are independently varying observable Euclidean
features.

COEFFICIENT RANGE = [-2,-1,1,2]
OFFSET RANGE      = [-4,-3,-2,-1,0,1,2,3,4]

SEARCH METHOD
  meet-in-the-middle exact integer signatures

VALIDATION
  full-data search
  genuine leave-one-out search

RULES
  exact integer arithmetic
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no arbitrary interpolation
========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations, product
from math import gcd
from typing import Dict, List, Sequence, Tuple


# ======================================================================================================================
# CONFIGURATION
# ======================================================================================================================

COEFFS = (-2, -1, 1, 2)
OFFSETS = tuple(range(-4, 5))


# ======================================================================================================================
# INSTANCES
# ======================================================================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int
    r: int
    true_d0: int
    true_x: int
    true_K: int
    true_j: int


INSTANCES: Tuple[Instance, ...] = (
    Instance(
        1,
        14246098189,
        333010,
        14245432170,
        -13799399066930810668576,
        14245432180,
        -10,
        -3449849766661475506269,
        5,
    ),
    Instance(
        2,
        10139117,
        11022,
        10117074,
        -11873391125935945,
        10117091,
        -17,
        -2968347695488785,
        9,
    ),
    Instance(
        3,
        10009330297,
        1010042,
        10007310214,
        -17225157737347240276460,
        10007310238,
        -24,
        -4306289434216722346403,
        12,
    ),
    Instance(
        4,
        100460333,
        20046,
        100420242,
        -2460554951578020025,
        100420273,
        -31,
        -615138736337991015,
        16,
    ),
    Instance(
        5,
        2503701173,
        100074,
        2503501026,
        -2231236374548471123872,
        2503501064,
        -38,
        -557809093589551261113,
        19,
    ),
    Instance(
        6,
        10006200817,
        200062,
        10005800694,
        -50858953461762196260849,
        10005800739,
        -45,
        -12714738365215418549091,
        23,
    ),
    Instance(
        7,
        40005200153,
        400026,
        40004400102,
        -1222668959330484861858204,
        40004400154,
        -52,
        -305667239831581101061223,
        26,
    ),
    Instance(
        8,
        270017400119,
        1200024,
        270015000072,
        -74949527189645489927322133,
        270015000131,
        -59,
        -18737381797403407039327539,
        30,
    ),
)


# ======================================================================================================================
# PARITY-ROBUST REFERENCE DERIVATION
# ======================================================================================================================

@dataclass(frozen=True)
class Reference:
    x0: int
    d0_ref: int
    K_ref: int


def derive_reference_candidates(inst: Instance) -> List[Reference]:
    candidates: List[Reference] = []

    target = inst.d * inst.d - inst.r

    for x0 in (0, 1):
        d0_ref = inst.d - x0

        # Exact square compatibility.
        if (d0_ref * d0_ref - target) % 4 != 0:
            continue

        numerator = (
            inst.r
            - 2 * inst.d * x0
            + x0 * x0
        )

        if numerator % 4 != 0:
            continue

        K_ref = numerator // 4

        # Direct exact verification.
        if (
            4 * K_ref
            != inst.r
            + d0_ref * d0_ref
            - inst.d * inst.d
        ):
            continue

        candidates.append(
            Reference(
                x0=x0,
                d0_ref=d0_ref,
                K_ref=K_ref,
            )
        )

    return candidates


def derive_reference(inst: Instance) -> Reference:
    refs = derive_reference_candidates(inst)

    if len(refs) != 1:
        raise AssertionError(
            f"expected exactly one reference: "
            f"instance={inst.idx}, refs={refs}"
        )

    return refs[0]


# ======================================================================================================================
# FIBRE
# ======================================================================================================================

def fibre_values(
    inst: Instance,
    ref: Reference,
    j: int,
) -> Tuple[int, int, int]:

    d0 = ref.d0_ref + 2 * j
    x = ref.x0 - 2 * j
    K = ref.K_ref + ref.d0_ref * j + j * j

    return d0, x, K


def direct_K(inst: Instance, d0: int) -> int:
    numerator = (
        inst.r
        + d0 * d0
        - inst.d * inst.d
    )

    if numerator % 4 != 0:
        raise AssertionError(
            f"non-integral direct K: instance={inst.idx}"
        )

    return numerator // 4


def reference_consistency_failures() -> int:
    failures = 0

    for inst in INSTANCES:

        refs = derive_reference_candidates(inst)

        if len(refs) != 1:
            failures += 1
            continue

        ref = refs[0]

        if ref.x0 + ref.d0_ref != inst.d:
            failures += 1

        d0, x, K = fibre_values(
            inst,
            ref,
            inst.true_j,
        )

        if d0 != inst.true_d0:
            failures += 1

        if x != inst.true_x:
            failures += 1

        if K != inst.true_K:
            failures += 1

        try:
            if direct_K(inst, d0) != K:
                failures += 1
        except AssertionError:
            failures += 1

    return failures


# ======================================================================================================================
# OBSERVABLE EUCLIDEAN FEATURES
# ======================================================================================================================

def build_features(
    instances: Sequence[Instance],
) -> Dict[str, List[int]]:

    features: Dict[str, List[int]] = {}

    def add(name: str, fn) -> None:
        features[name] = [
            int(fn(inst))
            for inst in instances
        ]

    add("N", lambda i: i.N)
    add("S", lambda i: i.S)
    add("d", lambda i: i.d)
    add("r", lambda i: i.r)

    add("q(N,S)", lambda i: divmod(i.N, i.S)[0])
    add("q(N,d)", lambda i: divmod(i.N, i.d)[0])
    add("q(S,d)", lambda i: divmod(i.S, i.d)[0])
    add("q(d,r)", lambda i: divmod(i.d, i.r)[0])
    add("q(N,r)", lambda i: divmod(i.N, i.r)[0])
    add("q(S,r)", lambda i: divmod(i.S, i.r)[0])

    add("rem(N,S)", lambda i: divmod(i.N, i.S)[1])
    add("rem(N,d)", lambda i: divmod(i.N, i.d)[1])
    add("rem(S,d)", lambda i: divmod(i.S, i.d)[1])
    add("rem(d,r)", lambda i: divmod(i.d, i.r)[1])
    add("rem(N,r)", lambda i: divmod(i.N, i.r)[1])
    add("rem(S,r)", lambda i: divmod(i.S, i.r)[1])

    add("gcd(N,S)", lambda i: gcd(i.N, i.S))
    add("gcd(N,d)", lambda i: gcd(i.N, i.d))
    add("gcd(S,d)", lambda i: gcd(i.S, i.d))
    add("gcd(d,r)", lambda i: gcd(i.d, i.r))
    add("gcd(N,r)", lambda i: gcd(i.N, i.r))
    add("gcd(S,r)", lambda i: gcd(i.S, i.r))

    add(
        "q(N,S)-q(S,d)",
        lambda i:
            divmod(i.N, i.S)[0]
            - divmod(i.S, i.d)[0],
    )

    add(
        "q(N,d)-q(S,d)",
        lambda i:
            divmod(i.N, i.d)[0]
            - divmod(i.S, i.d)[0],
    )

    add(
        "q(S,d)-q(N,S)",
        lambda i:
            divmod(i.S, i.d)[0]
            - divmod(i.N, i.S)[0],
    )

    add(
        "rem(N,S)-rem(S,d)",
        lambda i:
            divmod(i.N, i.S)[1]
            - divmod(i.S, i.d)[1],
    )

    add(
        "rem(N,d)-rem(S,d)",
        lambda i:
            divmod(i.N, i.d)[1]
            - divmod(i.S, i.d)[1],
    )

    add(
        "rem(N,S)-rem(N,d)",
        lambda i:
            divmod(i.N, i.S)[1]
            - divmod(i.N, i.d)[1],
    )

    add(
        "gcd(N,S)-gcd(S,d)",
        lambda i:
            gcd(i.N, i.S)
            - gcd(i.S, i.d),
    )

    add(
        "gcd(N,d)-gcd(S,d)",
        lambda i:
            gcd(i.N, i.d)
            - gcd(i.S, i.d),
    )

    return features


def variable_features(
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> List[str]:

    return [
        name
        for name, values in features.items()
        if len({values[i] for i in indices}) > 1
    ]


# ======================================================================================================================
# CANDIDATE REPRESENTATION
# ======================================================================================================================

@dataclass(frozen=True)
class NumeratorPartial:
    F: str
    G: str
    a: int
    b: int
    c: int


@dataclass(frozen=True)
class DenominatorPartial:
    H: str
    d: int
    e: int


@dataclass(frozen=True)
class Candidate:
    F: str
    G: str
    H: str
    Q: str
    a: int
    b: int
    c: int
    d: int
    e: int

    def numerator(
        self,
        features: Dict[str, List[int]],
        idx: int,
    ) -> int:

        F = features[self.F][idx]
        G = features[self.G][idx]
        H = features[self.H][idx]

        return (
            self.a * F * G
            + self.b * H
            + self.c
        )

    def denominator(
        self,
        features: Dict[str, List[int]],
        idx: int,
    ) -> int:

        Q = features[self.Q][idx]

        return (
            self.d * Q
            + self.e
        )

    def holds(
        self,
        features: Dict[str, List[int]],
        idx: int,
    ) -> bool:

        den = self.denominator(
            features,
            idx,
        )

        if den == 0:
            return False

        return (
            self.numerator(features, idx)
            == INSTANCES[idx].true_j * den
        )

    def key(self) -> Tuple:
        return (
            self.F,
            self.G,
            self.H,
            self.Q,
            self.a,
            self.b,
            self.c,
            self.d,
            self.e,
        )

    def __str__(self) -> str:
        return (
            f"{self.a}*{self.F}*{self.G} + "
            f"{self.b}*{self.H} + "
            f"{self.c} = "
            f"j*({self.d}*{self.Q} + {self.e})"
        )


# ======================================================================================================================
# NONDEGENERACY
# ======================================================================================================================

def nondegenerate(
    cand: Candidate,
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> bool:

    # Require three distinct numerator/denominator features.
    if len({
        cand.F,
        cand.G,
        cand.H,
        cand.Q,
    }) != 4:
        return False

    # All four features must actually vary.
    for name in (
        cand.F,
        cand.G,
        cand.H,
        cand.Q,
    ):
        if len({
            features[name][i]
            for i in indices
        }) <= 1:
            return False

    # No coefficient degeneracy.
    if cand.a == 0:
        return False

    if cand.b == 0:
        return False

    if cand.d == 0:
        return False

    # Denominator must stay nonzero.
    for i in indices:
        if cand.denominator(
            features,
            i,
        ) == 0:
            return False

    return True


# ======================================================================================================================
# NUMERATOR SIGNATURES
# ======================================================================================================================

def build_numerator_signatures(
    names: Sequence[str],
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> Dict[
    Tuple[int, ...],
    List[NumeratorPartial],
]:

    buckets: Dict[
        Tuple[int, ...],
        List[NumeratorPartial],
    ] = {}

    for F, G in permutations(names, 2):

        for H in names:

            if H == F or H == G:
                continue

            for a, b in product(
                COEFFS,
                repeat=2,
            ):

                for c in OFFSETS:

                    signature = tuple(
                        a * features[F][i]
                        * features[G][i]
                        + b * features[H][i]
                        + c
                        for i in indices
                    )

                    item = NumeratorPartial(
                        F=F,
                        G=G,
                        a=a,
                        b=b,
                        c=c,
                    )

                    buckets.setdefault(
                        signature,
                        [],
                    ).append(item)

    return buckets


# ======================================================================================================================
# DENOMINATOR SIGNATURES
# ======================================================================================================================

def build_denominator_signatures(
    names: Sequence[str],
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> Dict[
    Tuple[int, ...],
    List[DenominatorPartial],
]:

    buckets: Dict[
        Tuple[int, ...],
        List[DenominatorPartial],
    ] = {}

    for Q in names:

        for d, e in product(
            COEFFS,
            OFFSETS,
        ):

            signature = []
            valid = True

            for i in indices:

                den = (
                    d * features[Q][i]
                    + e
                )

                if den == 0:
                    valid = False
                    break

                signature.append(
                    INSTANCES[i].true_j * den
                )

            if not valid:
                continue

            item = DenominatorPartial(
                H="",
                d=d,
                e=e,
            )

            buckets.setdefault(
                tuple(signature),
                [],
            ).append(
                (Q, item)
            )

    return buckets


# ======================================================================================================================
# SEARCH
# ======================================================================================================================

def search_training_set(
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> List[Candidate]:

    names = variable_features(
        features,
        indices,
    )

    numerator = build_numerator_signatures(
        names,
        features,
        indices,
    )

    denominator = build_denominator_signatures(
        names,
        features,
        indices,
    )

    result: Dict[
        Tuple,
        Candidate,
    ] = {}

    for signature, nitems in numerator.items():

        ditem = denominator.get(
            signature
        )

        if ditem is None:
            continue

        for n in nitems:

            for Q, dpart in ditem:

                cand = Candidate(
                    F=n.F,
                    G=n.G,
                    H=n.H,
                    Q=Q,
                    a=n.a,
                    b=n.b,
                    c=n.c,
                    d=dpart.d,
                    e=dpart.e,
                )

                if not nondegenerate(
                    cand,
                    features,
                    indices,
                ):
                    continue

                result[cand.key()] = cand

    return list(result.values())


# ======================================================================================================================
# LEAVE-ONE-OUT
# ======================================================================================================================

def genuine_loo_search(
    features: Dict[str, List[int]],
) -> Tuple[List[Candidate], List[int]]:

    n = len(INSTANCES)

    fold_sets: List[
        Dict[Tuple, Candidate]
    ] = []

    held_out_counts = [0] * n

    for held_out in range(n):

        training = [
            i
            for i in range(n)
            if i != held_out
        ]

        candidates = search_training_set(
            features,
            training,
        )

        survivors: Dict[
            Tuple,
            Candidate,
        ] = {}

        for cand in candidates:

            if cand.holds(
                features,
                held_out,
            ):
                survivors[cand.key()] = cand

        fold_sets.append(
            survivors
        )

        held_out_counts[
            held_out
        ] = len(
            survivors
        )

    common = set(
        fold_sets[0]
    )

    for fold in fold_sets[1:]:
        common &= set(fold)

    final = [
        fold_sets[0][key]
        for key in sorted(common)
    ]

    return final, held_out_counts


# ======================================================================================================================
# DISPLAY
# ======================================================================================================================

def print_candidate(
    cand: Candidate,
    features: Dict[str, List[int]],
) -> None:

    print(
        "  MIXED-QUADRATIC CANDIDATE"
    )

    print(
        f"    {cand}"
    )

    print()

    for idx, inst in enumerate(
        INSTANCES
    ):

        numerator = cand.numerator(
            features,
            idx,
        )

        denominator = cand.denominator(
            features,
            idx,
        )

        rhs = (
            inst.true_j
            * denominator
        )

        print(
            f"    inst={inst.idx} "
            f"j={inst.true_j} "
            f"numerator={numerator} "
            f"den={denominator} "
            f"rhs={rhs} "
            f"match={numerator == rhs}"
        )

    print()


# ======================================================================================================================
# REFERENCE DISPLAY
# ======================================================================================================================

def print_reference_audit() -> None:

    print("=" * 120)
    print(
        "REFERENCE / PARITY AUDIT"
    )
    print("=" * 120)
    print()

    for inst in INSTANCES:

        refs = derive_reference_candidates(
            inst
        )

        print(
            f"INSTANCE {inst.idx}"
        )

        print(
            f"  candidate parity references = "
            f"{len(refs)}"
        )

        for ref in refs:
            print(
                f"    x0={ref.x0} "
                f"d0_ref={ref.d0_ref} "
                f"K_ref={ref.K_ref}"
            )

        print()


# ======================================================================================================================
# TRUE BRANCH CHECK
# ======================================================================================================================

def true_branch_check() -> int:

    failures = 0

    print("=" * 120)
    print(
        "TRUE-BRANCH EXACT RECONSTRUCTION"
    )
    print("=" * 120)
    print()

    for inst in INSTANCES:

        ref = derive_reference(
            inst
        )

        d0, x, K = fibre_values(
            inst,
            ref,
            inst.true_j,
        )

        K_direct = direct_K(
            inst,
            d0,
        )

        exact = (
            d0 == inst.true_d0
            and x == inst.true_x
            and K == inst.true_K
            and K == K_direct
        )

        if not exact:
            failures += 1

        print(
            f"  instance={inst.idx} "
            f"j={inst.true_j:<3} "
            f"x0={ref.x0} "
            f"d0_ref={ref.d0_ref} "
            f"d0={d0} "
            f"x={x} "
            f"K={K} "
            f"K_direct={K_direct} "
            f"exact={exact}"
        )

    print()

    print(
        f"  reconstruction failures = "
        f"{failures}"
    )

    print()

    return failures


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 473")
    print("=" * 120)
    print()

    print(
        "SPARSE MIXED-QUADRATIC / EUCLIDEAN j-LAW AUDIT"
    )
    print()

    print(
        "SEARCH FORM"
    )
    print(
        "  a*F*G + b*H + c = j*(d*Q + e)"
    )
    print()

    print(
        f"COEFFICIENT RANGE = {list(COEFFS)}"
    )

    print(
        f"OFFSET RANGE      = {list(OFFSETS)}"
    )

    print(
        "SEARCH METHOD     = meet-in-the-middle signatures"
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # INPUT CONSISTENCY
    # --------------------------------------------------------------------------------------------------------------

    print(
        "INPUT CONSISTENCY"
    )

    consistency = (
        reference_consistency_failures()
    )

    print(
        f"  consistency failures = {consistency}"
    )

    print()

    if consistency:
        raise AssertionError(
            "Reference consistency failed."
        )

    print_reference_audit()

    # --------------------------------------------------------------------------------------------------------------
    # FEATURES
    # --------------------------------------------------------------------------------------------------------------

    features = build_features(
        INSTANCES
    )

    indices = list(
        range(len(INSTANCES))
    )

    variable = variable_features(
        features,
        indices,
    )

    print("=" * 120)
    print(
        "NONDEGENERATE MIXED-QUADRATIC FEATURE SCREEN"
    )
    print("=" * 120)
    print()

    for name, values in features.items():

        distinct = len(set(values))

        print(
            f"  {name:<35} "
            f"distinct={distinct:<3} "
            f"{'VARIABLE' if distinct > 1 else 'CONSTANT'}"
        )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # SEARCH SIZE
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "SEARCH SIZE"
    )
    print("=" * 120)
    print()

    print(
        f"  variable features                  = "
        f"{len(variable)}"
    )

    print(
        f"  ordered F,G feature pairs          = "
        f"{len(variable) * (len(variable) - 1)}"
    )

    print(
        f"  third numerator feature H choices  = "
        f"{max(0, len(variable) - 2)}"
    )

    print(
        f"  numerator coefficient tuples       = "
        f"{len(COEFFS) ** 2}"
    )

    print(
        f"  offsets                            = "
        f"{len(OFFSETS)}"
    )

    print(
        f"  denominator feature choices       = "
        f"{len(variable)}"
    )

    print(
        f"  denominator coefficient choices   = "
        f"{len(COEFFS) * len(OFFSETS)}"
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # FULL SEARCH
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "MIXED-QUADRATIC FULL-DATA SEARCH"
    )
    print("=" * 120)
    print()

    full_candidates = search_training_set(
        features,
        indices,
    )

    print(
        f"full-data nondegenerate candidates = "
        f"{len(full_candidates)}"
    )

    if full_candidates:

        for cand in full_candidates[:20]:
            print_candidate(
                cand,
                features,
            )

        if len(full_candidates) > 20:
            print(
                "  ... output truncated after 20 candidates"
            )

    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # LOO
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "GENUINE LEAVE-ONE-OUT SEARCH"
    )
    print("=" * 120)
    print()

    loo_candidates, held_out_counts = (
        genuine_loo_search(features)
    )

    print(
        "  held-out instance survival counts:"
    )

    for idx, count in enumerate(
        held_out_counts,
        start=1,
    ):
        print(
            f"    instance {idx}: {count}"
        )

    print()

    print(
        f"unique globally LOO-surviving candidates = "
        f"{len(loo_candidates)}"
    )

    if loo_candidates:

        for cand in loo_candidates[:20]:
            print_candidate(
                cand,
                features,
            )

        if len(loo_candidates) > 20:
            print(
                "  ... output truncated after 20 candidates"
            )

    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # TRUE BRANCH
    # --------------------------------------------------------------------------------------------------------------

    reconstruction_failures = (
        true_branch_check()
    )

    # --------------------------------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "GLOBAL EXPERIMENT 473 SUMMARY"
    )
    print("=" * 120)
    print()

    print(
        f"  instances                         = "
        f"{len(INSTANCES)}"
    )

    print(
        f"  variable observable features     = "
        f"{len(variable)}"
    )

    print(
        f"  coefficient range                = "
        f"{list(COEFFS)}"
    )

    print(
        f"  offset range                     = "
        f"{list(OFFSETS)}"
    )

    print(
        f"  full-data nondegenerate laws     = "
        f"{len(full_candidates)}"
    )

    print(
        f"  globally LOO laws                = "
        f"{len(loo_candidates)}"
    )

    print(
        f"  reconstruction failures          = "
        f"{reconstruction_failures}"
    )

    print()

    print(
        "INTERPRETATION"
    )

    print()

    print(
        "  Experiment 472 tested sparse pure"
    )

    print(
        "  quadratic numerator terms:"
    )

    print(
        "      aF^2 + bG + c = j(dH + e)"
    )

    print()

    print(
        "  Experiment 473 instead tests a mixed"
    )

    print(
        "  interaction between two independent"
    )

    print(
        "  observable features:"
    )

    print(
        "      aFG + bH + c = j(dQ + e)"
    )

    print()

    print(
        "  The mixed term can expose relationships"
    )

    print(
        "  which are invisible to a basis containing"
    )

    print(
        "  only single-feature powers."
    )

    print()

    print(
        "  The family remains deliberately bounded"
    )

    print(
        "  and requires genuine leave-one-out"
    )

    print(
        "  validation."
    )

    print()

    # --------------------------------------------------------------------------------------------------------------
    # FINAL
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print(
        "EXPERIMENT 473 FINAL STATUS"
    )
    print("=" * 120)
    print()

    print(
        "  PARITY-ROBUST REFERENCE DERIVATION = True"
    )

    print(
        "  EXACT INTEGER ARITHMETIC            = True"
    )

    print(
        "  SPARSE MIXED-QUADRATIC FAMILY      = True"
    )

    print(
        "  MEET-IN-THE-MIDDLE SEARCH           = True"
    )

    print(
        "  NONDEGENERATE FILTER                = True"
    )

    print(
        "  GENUINE LEAVE-ONE-OUT               = True"
    )

    print(
        "  FLOATING POINT                       = False"
    )

    print(
        "  ARBITRARY INTERPOLATION             = False"
    )

    print(
        "  GIANT K ENUMERATION                 = False"
    )

    print(
        "  GIANT d0 ENUMERATION                = False"
    )

    print(
        "  CRT CARTESIAN PRODUCT               = False"
    )

    print()

    if reconstruction_failures:
        print(
            "  CONCLUSION:"
        )
        print(
            "    Reference reconstruction failed."
        )

    elif loo_candidates:
        print(
            "  CONCLUSION:"
        )
        print(
            "    A bounded mixed-quadratic"
        )
        print(
            "    Euclidean-derived j-law survived"
        )
        print(
            "    genuine leave-one-out validation."
        )

    else:
        print(
            "  CONCLUSION:"
        )
        print(
            "    No bounded sparse mixed-quadratic"
        )
        print(
            "    Euclidean-derived rational j-law"
        )
        print(
            "    survived genuine leave-one-out"
        )
        print(
            "    validation."
        )

    print()

    print("=" * 120)
    print(
        "EXPERIMENT 473 FINISHED"
    )
    print("=" * 120)


if __name__ == "__main__":
    main()
