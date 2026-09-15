#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 474
========================================================================================================================

SPARSE CUBIC / EUCLIDEAN j-LAW AUDIT

QUESTION

  Does a low-complexity cubic interaction between three
  independently varying observable features determine j
  across the eight instances?

SEARCH FORM

      a*F*G*H + c = j*(d*Q + e)

where:

  F, G, H, Q are distinct observable Euclidean-derived
  features.

COEFFICIENT RANGE = [-2,-1,1,2]
OFFSET RANGE      = [-4,-3,-2,-1,0,1,2,3,4]

SEARCH METHOD
  exact signature hashing / meet-in-the-middle

VALIDATION
  full-data search
  genuine leave-one-out search

RULES
  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no arbitrary interpolation
  no hidden-variable features

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
# REFERENCE / PARITY
# ======================================================================================================================

@dataclass(frozen=True)
class Reference:
    x0: int
    d0_ref: int
    K_ref: int


def derive_reference_candidates(
    inst: Instance,
) -> List[Reference]:

    candidates: List[Reference] = []

    # d0 parity follows from d0^2 == d^2-r (mod 4).
    target = inst.d * inst.d - inst.r

    for x0 in (0, 1):
        d0_ref = inst.d - x0

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


def derive_reference(
    inst: Instance,
) -> Reference:

    refs = derive_reference_candidates(inst)

    if len(refs) != 1:
        raise AssertionError(
            f"expected exactly one parity-compatible reference: "
            f"instance={inst.idx}, candidates={len(refs)}"
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


def direct_K(
    inst: Instance,
    d0: int,
) -> int:

    numerator = (
        inst.r
        + d0 * d0
        - inst.d * inst.d
    )

    if numerator % 4 != 0:
        raise AssertionError(
            f"direct K not integral: instance={inst.idx}"
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
# OBSERVABLE EUCLIDEAN FEATURE BASIS
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
class CubicCandidate:
    F: str
    G: str
    H: str
    Q: str
    a: int
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
            self.a * F * G * H
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

        denominator = self.denominator(
            features,
            idx,
        )

        if denominator == 0:
            return False

        return (
            self.numerator(features, idx)
            == INSTANCES[idx].true_j * denominator
        )

    def key(self) -> Tuple:
        return (
            self.F,
            self.G,
            self.H,
            self.Q,
            self.a,
            self.c,
            self.d,
            self.e,
        )

    def __str__(self) -> str:
        return (
            f"{self.a}*{self.F}*{self.G}*{self.H}"
            f" + {self.c}"
            f" = j*({self.d}*{self.Q} + {self.e})"
        )


# ======================================================================================================================
# NONDEGENERACY
# ======================================================================================================================

def nondegenerate(
    cand: CubicCandidate,
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> bool:

    names = (
        cand.F,
        cand.G,
        cand.H,
        cand.Q,
    )

    # All features must be distinct.
    if len(set(names)) != 4:
        return False

    # All features must genuinely vary in the training set.
    for name in names:

        if len({
            features[name][i]
            for i in indices
        }) <= 1:
            return False

    # Nonzero leading / denominator coefficient.
    if cand.a == 0:
        return False

    if cand.d == 0:
        return False

    # Denominator must remain nonzero.
    for i in indices:

        if (
            cand.d * features[cand.Q][i]
            + cand.e
        ) == 0:
            return False

    return True


# ======================================================================================================================
# NUMERATOR SIGNATURE BUILD
# ======================================================================================================================

def build_numerator_signatures(
    names: Sequence[str],
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> Dict[
    Tuple[int, ...],
    List[Tuple[str, str, str, int, int]],
]:

    buckets: Dict[
        Tuple[int, ...],
        List[Tuple[str, str, str, int, int]],
    ] = {}

    # Ordered triples.
    for F, G, H in permutations(names, 3):

        for a in COEFFS:

            for c in OFFSETS:

                signature = tuple(
                    a
                    * features[F][i]
                    * features[G][i]
                    * features[H][i]
                    + c
                    for i in indices
                )

                item = (
                    F,
                    G,
                    H,
                    a,
                    c,
                )

                buckets.setdefault(
                    signature,
                    [],
                ).append(item)

    return buckets


# ======================================================================================================================
# DENOMINATOR SIGNATURE BUILD
# ======================================================================================================================

def build_denominator_signatures(
    names: Sequence[str],
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> Dict[
    Tuple[int, ...],
    List[Tuple[str, int, int]],
]:

    buckets: Dict[
        Tuple[int, ...],
        List[Tuple[str, int, int]],
    ] = {}

    for Q in names:

        for d, e in product(
            COEFFS,
            OFFSETS,
        ):

            signature = []
            valid = True

            for i in indices:

                denominator = (
                    d * features[Q][i]
                    + e
                )

                if denominator == 0:
                    valid = False
                    break

                signature.append(
                    INSTANCES[i].true_j
                    * denominator
                )

            if not valid:
                continue

            key = tuple(signature)

            buckets.setdefault(
                key,
                [],
            ).append(
                (Q, d, e)
            )

    return buckets


# ======================================================================================================================
# SEARCH
# ======================================================================================================================

def search_training_set(
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> List[CubicCandidate]:

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
        CubicCandidate,
    ] = {}

    for signature, numerator_items in numerator.items():

        denominator_items = denominator.get(
            signature
        )

        if denominator_items is None:
            continue

        for (
            F,
            G,
            H,
            a,
            c,
        ) in numerator_items:

            for (
                Q,
                d,
                e,
            ) in denominator_items:

                cand = CubicCandidate(
                    F=F,
                    G=G,
                    H=H,
                    Q=Q,
                    a=a,
                    c=c,
                    d=d,
                    e=e,
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
) -> Tuple[List[CubicCandidate], List[int]]:

    n = len(INSTANCES)

    fold_survivors: List[
        Dict[Tuple, CubicCandidate]
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
            CubicCandidate
        ] = {}

        for cand in candidates:

            if cand.holds(
                features,
                held_out,
            ):
                survivors[cand.key()] = cand

        fold_survivors.append(
            survivors
        )

        held_out_counts[
            held_out
        ] = len(
            survivors
        )

    common = set(
        fold_survivors[0]
    )

    for fold in fold_survivors[1:]:

        common &= set(fold)

    final = [
        fold_survivors[0][key]
        for key in sorted(common)
    ]

    return final, held_out_counts


# ======================================================================================================================
# CANDIDATE DISPLAY
# ======================================================================================================================

def print_candidate(
    cand: CubicCandidate,
    features: Dict[str, List[int]],
) -> None:

    print(
        "  CUBIC CANDIDATE"
    )

    print(
        f"    {cand}"
    )

    print()

    for idx, inst in enumerate(INSTANCES):

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
# REFERENCE AUDIT
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
        f"  reconstruction failures = {failures}"
    )

    print()

    return failures


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    print("=" * 120)
    print(
        "EXPERIMENT 474"
    )
    print("=" * 120)
    print()

    print(
        "SPARSE CUBIC / EUCLIDEAN j-LAW AUDIT"
    )
    print()

    print(
        "SEARCH FORM"
    )
    print(
        "  a*F*G*H + c = j*(d*Q + e)"
    )
    print()

    print(
        f"COEFFICIENT RANGE = {list(COEFFS)}"
    )

    print(
        f"OFFSET RANGE      = {list(OFFSETS)}"
    )

    print(
        "SEARCH METHOD     = exact signature hashing"
    )

    print()

    # ==============================================================================================================
    # CONSISTENCY
    # ==============================================================================================================

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

    # ==============================================================================================================
    # FEATURES
    # ==============================================================================================================

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
        "NONDEGENERATE CUBIC FEATURE SCREEN"
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

    # ==============================================================================================================
    # SEARCH SIZE
    # ==============================================================================================================

    print("=" * 120)
    print(
        "SEARCH SIZE"
    )
    print("=" * 120)
    print()

    print(
        f"  variable features                  = {len(variable)}"
    )

    print(
        f"  ordered F,G,H feature triples      = "
        f"{len(variable) * (len(variable) - 1) * (len(variable) - 2)}"
    )

    print(
        f"  coefficient choices                = {len(COEFFS)}"
    )

    print(
        f"  offsets                            = {len(OFFSETS)}"
    )

    print(
        f"  denominator feature choices       = {len(variable)}"
    )

    print(
        f"  denominator coefficient choices   = "
        f"{len(COEFFS) * len(OFFSETS)}"
    )

    print()

    # ==============================================================================================================
    # FULL SEARCH
    # ==============================================================================================================

    print("=" * 120)
    print(
        "SPARSE CUBIC FULL-DATA SEARCH"
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
        print(
            "  none"
        )

    print()

    # ==============================================================================================================
    # LOO
    # ==============================================================================================================

    print("=" * 120)
    print(
        "GENUINE LEAVE-ONE-OUT SEARCH"
    )
    print("=" * 120)
    print()

    loo_candidates, held_out_counts = (
        genuine_loo_search(
            features
        )
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
        print(
            "  none"
        )

    print()

    # ==============================================================================================================
    # TRUE BRANCH
    # ==============================================================================================================

    reconstruction_failures = (
        true_branch_check()
    )

    # ==============================================================================================================
    # SUMMARY
    # ==============================================================================================================

    print("=" * 120)
    print(
        "GLOBAL EXPERIMENT 474 SUMMARY"
    )
    print("=" * 120)
    print()

    print(
        f"  instances                         = {len(INSTANCES)}"
    )

    print(
        f"  variable observable features     = {len(variable)}"
    )

    print(
        f"  coefficient range                = {list(COEFFS)}"
    )

    print(
        f"  offset range                     = {list(OFFSETS)}"
    )

    print(
        f"  full-data nondegenerate laws     = {len(full_candidates)}"
    )

    print(
        f"  globally LOO laws                = {len(loo_candidates)}"
    )

    print(
        f"  reconstruction failures          = {reconstruction_failures}"
    )

    print()

    print(
        "INTERPRETATION"
    )

    print()

    print(
        "  Experiment 473 tested the mixed"
    )

    print(
        "  quadratic interaction:"
    )

    print(
        "      aFG + bH + c = j(dQ + e)"
    )

    print()

    print(
        "  Experiment 474 raises the nonlinear"
    )

    print(
        "  interaction to three observable"
    )

    print(
        "  factors:"
    )

    print(
        "      aFGH + c = j(dQ + e)"
    )

    print()

    print(
        "  This tests whether a sparse cubic"
    )

    print(
        "  interaction among independent"
    )

    print(
        "  observable Euclidean features can"
    )

    print(
        "  encode the hidden displacement."
    )

    print()

    print(
        "  The search remains bounded and"
    )

    print(
        "  requires genuine leave-one-out"
    )

    print(
        "  validation."
    )

    print()

    print("=" * 120)
    print(
        "EXPERIMENT 474 FINAL STATUS"
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
        "  SPARSE CUBIC RATIONAL FAMILY       = True"
    )

    print(
        "  SIGNATURE-HASH SEARCH               = True"
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
            "    A bounded sparse cubic"
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
            "    No bounded sparse cubic"
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
        "EXPERIMENT 474 FINISHED"
    )
    print("=" * 120)


if __name__ == "__main__":
    main()
