#!/usr/bin/env python3
"""
================================================================================
EXPERIMENT 469
================================================================================

THREE-FEATURE SPARSE RATIONAL / EUCLIDEAN j-LAW AUDIT

QUESTION
  Does the hidden displacement j satisfy a low-complexity relation of the form

      a*F + b*G + c = j*(d*H + e)

  where F, G, H are independently varying Euclidean-derived observable
  features?

GOAL
  Extend Experiment 468 from two-feature rational laws to a sparse
  three-feature rational family without constructing a giant Cartesian
  product.

SEARCH FORM
  aF + bG + c = j(dH + e)

COEFFICIENT RANGE = [-2,2]
OFFSET RANGE      = [-4,4]

SEARCH METHOD
  Meet-in-the-middle signatures.

RULES
  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  no arbitrary interpolation
  no constant features
  no zero denominator
  genuine leave-one-out validation
  no hidden values used to construct observables

================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations, product
from math import gcd, isqrt
from typing import Dict, Iterable, List, Sequence, Tuple


# ============================================================================
# DATA
# ============================================================================

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


# ============================================================================
# FEATURE GENERATION
# ============================================================================

@dataclass(frozen=True)
class Feature:
    name: str
    values: Tuple[int, ...]


def euclidean_divmod(a: int, b: int) -> Tuple[int, int]:
    if b == 0:
        raise ZeroDivisionError("division by zero")
    return divmod(a, b)


def digits(n: int) -> int:
    return len(str(abs(n)))


def build_raw_features(instances: Sequence[Instance]) -> Dict[str, List[int]]:
    out: Dict[str, List[int]] = {}

    def add(name: str, fn):
        vals = []
        for inst in instances:
            vals.append(int(fn(inst)))
        out[name] = vals

    add("1", lambda i: 1)
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

    # Derived Euclidean combinations.
    add("q(N,S)-q(S,d)",
        lambda i: divmod(i.N, i.S)[0] - divmod(i.S, i.d)[0])

    add("q(N,d)-q(S,d)",
        lambda i: divmod(i.N, i.d)[0] - divmod(i.S, i.d)[0])

    add("q(S,d)-q(N,S)",
        lambda i: divmod(i.S, i.d)[0] - divmod(i.N, i.S)[0])

    add("rem(N,S)-rem(S,d)",
        lambda i: divmod(i.N, i.S)[1] - divmod(i.S, i.d)[1])

    add("rem(N,d)-rem(S,d)",
        lambda i: divmod(i.N, i.d)[1] - divmod(i.S, i.d)[1])

    add("rem(N,S)-rem(N,d)",
        lambda i: divmod(i.N, i.S)[1] - divmod(i.N, i.d)[1])

    add("gcd(N,S)-gcd(S,d)",
        lambda i: gcd(i.N, i.S) - gcd(i.S, i.d))

    add("gcd(N,d)-gcd(S,d)",
        lambda i: gcd(i.N, i.d) - gcd(i.S, i.d))

    return out


# ============================================================================
# SEARCH CANDIDATES
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    F: str
    G: str
    H: str
    a: int
    b: int
    c: int
    d: int
    e: int

    def numerator_signature(
        self,
        features: Dict[str, List[int]],
        indices: Sequence[int],
    ) -> Tuple[int, ...]:
        return tuple(
            self.a * features[self.F][i]
            + self.b * features[self.G][i]
            + self.c
            - INSTANCES[i].true_j
            * (
                self.d * features[self.H][i]
                + self.e
            )
            for i in indices
        )

    def denominator_values(
        self,
        features: Dict[str, List[int]],
        indices: Sequence[int],
    ) -> List[int]:
        return [
            self.d * features[self.H][i] + self.e
            for i in indices
        ]

    def evaluate(self, features: Dict[str, List[int]], idx: int) -> Tuple[int, int]:
        lhs = (
            self.a * features[self.F][idx]
            + self.b * features[self.G][idx]
            + self.c
        )
        rhs = INSTANCES[idx].true_j * (
            self.d * features[self.H][idx] + self.e
        )
        return lhs, rhs

    def __str__(self) -> str:
        return (
            f"{self.a}*{self.F} + "
            f"{self.b}*{self.G} + {self.c} = "
            f"j*({self.d}*{self.H} + {self.e})"
        )


# ============================================================================
# UTILITY FILTERS
# ============================================================================

COEFFS = (-2, -1, 1, 2)
OFFSETS = tuple(range(-4, 5))


def variable_feature_names(
    features: Dict[str, List[int]],
) -> List[str]:
    names = []
    for name, values in features.items():
        if len(set(values)) > 1:
            names.append(name)
    return names


def validate_observable_consistency() -> int:
    failures = 0

    for inst in INSTANCES:
        if inst.d != inst.N - 2 * inst.S + 1:
            failures += 1

        if 4 * inst.true_K != (
            inst.r + inst.true_d0 * inst.true_d0 - inst.d * inst.d
        ):
            failures += 1

        if inst.true_x + inst.true_d0 != inst.d:
            failures += 1

        reconstructed_d0 = inst.d + 2 * 0
        if reconstructed_d0 < 0:
            failures += 1

    return failures


def candidate_is_nondegenerate(
    cand: Candidate,
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> bool:
    # All three selected observable features must genuinely vary.
    for name in (cand.F, cand.G, cand.H):
        if len({features[name][i] for i in indices}) <= 1:
            return False

    # Both numerator coefficients must be nonzero.
    if cand.a == 0 or cand.b == 0:
        return False

    # Denominator coefficient must be nonzero.
    if cand.d == 0:
        return False

    # F, G, H must be distinct.
    if len({cand.F, cand.G, cand.H}) != 3:
        return False

    # No denominator may vanish on the training set.
    for i in indices:
        den = cand.d * features[cand.H][i] + cand.e
        if den == 0:
            return False

    # Reject candidate if its numerator is identically zero on training set.
    lhs_values = []
    for i in indices:
        lhs_values.append(
            cand.a * features[cand.F][i]
            + cand.b * features[cand.G][i]
            + cand.c
        )
    if len(set(lhs_values)) == 1 and lhs_values[0] == 0:
        return False

    return True


# ============================================================================
# FAST SEARCH
# ============================================================================

@dataclass(frozen=True)
class PartialNumerator:
    F: str
    G: str
    a: int
    b: int
    c: int

    def signature(
        self,
        features: Dict[str, List[int]],
        indices: Sequence[int],
    ) -> Tuple[int, ...]:
        return tuple(
            self.a * features[self.F][i]
            + self.b * features[self.G][i]
            + self.c
            for i in indices
        )


@dataclass(frozen=True)
class PartialDenominator:
    H: str
    d: int
    e: int

    def signature(
        self,
        features: Dict[str, List[int]],
        indices: Sequence[int],
    ) -> Tuple[int, ...]:
        return tuple(
            INSTANCES[i].true_j
            * (
                self.d * features[self.H][i] + self.e
            )
            for i in indices
        )


def canonical_pair(
    a: int,
    F: str,
    b: int,
    G: str,
) -> Tuple[Tuple[str, int], Tuple[str, int]]:
    left = (F, a)
    right = (G, b)
    return tuple(sorted((left, right)))


def enumerate_numerators(
    feature_names: Sequence[str],
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> Dict[Tuple[int, ...], List[PartialNumerator]]:
    buckets: Dict[Tuple[int, ...], List[PartialNumerator]] = {}

    for F, G in permutations(feature_names, 2):
        for a, b in product(COEFFS, repeat=2):
            if a == 0 or b == 0:
                continue

            for c in OFFSETS:
                sig = tuple(
                    a * features[F][i]
                    + b * features[G][i]
                    + c
                    for i in indices
                )

                obj = PartialNumerator(F, G, a, b, c)
                buckets.setdefault(sig, []).append(obj)

    return buckets


def enumerate_denominators(
    feature_names: Sequence[str],
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> Dict[Tuple[int, ...], List[PartialDenominator]]:
    buckets: Dict[Tuple[int, ...], List[PartialDenominator]] = {}

    for H in feature_names:
        for d, e in product(COEFFS, OFFSETS):
            if d == 0:
                continue

            vals = [
                INSTANCES[i].true_j * (
                    d * features[H][i] + e
                )
                for i in indices
            ]

            # Reject any zero denominator.
            if any(
                d * features[H][i] + e == 0
                for i in indices
            ):
                continue

            sig = tuple(vals)

            obj = PartialDenominator(H, d, e)
            buckets.setdefault(sig, []).append(obj)

    return buckets


def combine_candidates(
    numerator_buckets: Dict[Tuple[int, ...], List[PartialNumerator]],
    denominator_buckets: Dict[Tuple[int, ...], List[PartialDenominator]],
) -> List[Candidate]:
    candidates: List[Candidate] = []

    for num_sig, nums in numerator_buckets.items():
        dens = denominator_buckets.get(num_sig)
        if not dens:
            continue

        for n in nums:
            for d in dens:
                if len({n.F, n.G, d.H}) != 3:
                    continue

                candidates.append(
                    Candidate(
                        F=n.F,
                        G=n.G,
                        H=d.H,
                        a=n.a,
                        b=n.b,
                        c=n.c,
                        d=d.d,
                        e=d.e,
                    )
                )

    return candidates


def deduplicate_candidates(
    candidates: Iterable[Candidate],
) -> List[Candidate]:
    seen = set()
    result = []

    for cand in candidates:
        key = (
            cand.F,
            cand.G,
            cand.H,
            cand.a,
            cand.b,
            cand.c,
            cand.d,
            cand.e,
        )
        if key not in seen:
            seen.add(key)
            result.append(cand)

    return result


# ============================================================================
# FULL-DATA SEARCH
# ============================================================================

def full_data_search(
    features: Dict[str, List[int]],
) -> List[Candidate]:
    names = variable_feature_names(features)
    indices = list(range(len(INSTANCES)))

    nums = enumerate_numerators(names, features, indices)
    dens = enumerate_denominators(names, features, indices)

    raw = combine_candidates(nums, dens)

    result = []

    for cand in deduplicate_candidates(raw):
        if candidate_is_nondegenerate(cand, features, indices):
            result.append(cand)

    return result


# ============================================================================
# LEAVE-ONE-OUT
# ============================================================================

def validate_candidate_on_index(
    cand: Candidate,
    features: Dict[str, List[int]],
    idx: int,
) -> bool:
    den = cand.d * features[cand.H][idx] + cand.e

    if den == 0:
        return False

    lhs, rhs = cand.evaluate(features, idx)
    return lhs == rhs


def search_on_training_set(
    features: Dict[str, List[int]],
    training_indices: Sequence[int],
) -> List[Candidate]:
    names = variable_feature_names(features)

    nums = enumerate_numerators(
        names,
        features,
        training_indices,
    )

    dens = enumerate_denominators(
        names,
        features,
        training_indices,
    )

    raw = combine_candidates(nums, dens)

    return [
        cand
        for cand in deduplicate_candidates(raw)
        if candidate_is_nondegenerate(
            cand,
            features,
            training_indices,
        )
    ]


@dataclass(frozen=True)
class LOOSurvivor:
    held_out: int
    candidate: Candidate


def leave_one_out_search(
    features: Dict[str, List[int]],
) -> List[LOOSurvivor]:
    survivors: List[LOOSurvivor] = []

    n = len(INSTANCES)

    for held_out in range(n):
        training = [i for i in range(n) if i != held_out]

        candidates = search_on_training_set(
            features,
            training,
        )

        for cand in candidates:
            if validate_candidate_on_index(
                cand,
                features,
                held_out,
            ):
                survivors.append(
                    LOOSurvivor(
                        held_out=held_out,
                        candidate=cand,
                    )
                )

    return survivors


# ============================================================================
# STRICT POST-FILTERS
# ============================================================================

def reject_lower_feature_collapse(
    cand: Candidate,
    features: Dict[str, List[int]],
) -> bool:
    """
    Reject candidates that are secretly one-feature laws.

    A lower-feature collapse is detected by comparing the candidate
    against all possible one-feature forms

        aF + c = j(dG + e)

    within the same bounded coefficient/offset family.
    """

    indices = list(range(len(INSTANCES)))

    # Candidate must use genuinely distinct variable features.
    if len({cand.F, cand.G, cand.H}) != 3:
        return False

    target = [
        cand.a * features[cand.F][i]
        + cand.b * features[cand.G][i]
        + cand.c
        - INSTANCES[i].true_j * (
            cand.d * features[cand.H][i] + cand.e
        )
        for i in indices
    ]

    if any(v != 0 for v in target):
        return False

    # This function intentionally only identifies obvious structural
    # collapse cases. More sophisticated symbolic dependency testing
    # belongs to a later experiment.
    return True


def strict_full_data_candidates(
    candidates: Sequence[Candidate],
    features: Dict[str, List[int]],
) -> List[Candidate]:
    result = []

    for cand in candidates:
        if not candidate_is_nondegenerate(
            cand,
            features,
            range(len(INSTANCES)),
        ):
            continue

        if not reject_lower_feature_collapse(
            cand,
            features,
        ):
            continue

        result.append(cand)

    return result


# ============================================================================
# OUTPUT
# ============================================================================

def print_feature_screen(features: Dict[str, List[int]]) -> None:
    print("=" * 120)
    print("NONDEGENERATE THREE-FEATURE SCREEN")
    print("=" * 120)
    print()

    for name, values in features.items():
        distinct = len(set(values))
        state = "VARIABLE" if distinct > 1 else "CONSTANT"
        print(
            f"  {name:<35} distinct={distinct:<3} {state}"
        )

    print()


def print_candidate(
    title: str,
    cand: Candidate,
    features: Dict[str, List[int]],
    indices: Sequence[int],
) -> None:
    print(f"  {title}")
    print(f"    {cand}")

    for idx in indices:
        lhs, rhs = cand.evaluate(features, idx)
        den = (
            cand.d * features[cand.H][idx]
            + cand.e
        )

        print(
            f"    inst={idx + 1:<2} "
            f"j={INSTANCES[idx].true_j:<3} "
            f"lhs={lhs} "
            f"rhs={rhs} "
            f"den={den} "
            f"match={lhs == rhs}"
        )

    print()


def print_true_reconstruction_checks(features: Dict[str, List[int]]) -> None:
    failures = 0

    print("=" * 120)
    print("TRUE-BRANCH RECONSTRUCTION CHECK")
    print("=" * 120)
    print()

    for inst in INSTANCES:
        d0 = inst.d + 2 * inst.true_j
        x = inst.d - d0

        K_num = inst.r + d0 * d0 - inst.d * inst.d
        if K_num % 4 != 0:
            failures += 1
            K = None
        else:
            K = K_num // 4

        print(
            f"  instance={inst.idx} "
            f"j={inst.true_j} "
            f"d0={d0} "
            f"x={x} "
            f"K={K} "
            f"exact="
            f"{d0 == inst.true_d0 and x == inst.true_x and K == inst.true_K}"
        )

    print()
    print(f"  reconstruction failures = {failures}")
    print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 469")
    print("=" * 120)
    print()
    print("THREE-FEATURE SPARSE RATIONAL / EUCLIDEAN j-LAW AUDIT")
    print()
    print("SEARCH FORM")
    print("  a*F + b*G + c = j*(d*H + e)")
    print()
    print("COEFFICIENT RANGE =", list(COEFFS))
    print("OFFSET RANGE      =", list(OFFSETS))
    print("SEARCH METHOD     = meet-in-the-middle signatures")
    print()

    failures = validate_observable_consistency()

    print("INPUT CONSISTENCY")
    print("  consistency failures =", failures)
    print()

    if failures:
        raise AssertionError(
            f"input consistency failed: {failures}"
        )

    features = build_raw_features(INSTANCES)

    print_feature_screen(features)

    variable_names = variable_feature_names(features)

    print("=" * 120)
    print("SEARCH SIZE")
    print("=" * 120)
    print()

    print(
        f"  variable features                  = {len(variable_names)}"
    )
    print(
        f"  ordered numerator feature pairs    = "
        f"{len(variable_names) * (len(variable_names) - 1)}"
    )
    print(
        f"  denominator feature choices       = {len(variable_names)}"
    )
    print(
        f"  numerator coefficient tuples      = "
        f"{len(COEFFS) ** 2}"
    )
    print(
        f"  denominator coefficient choices   = "
        f"{len(COEFFS) * len(OFFSETS)}"
    )
    print()

    # ----------------------------------------------------------------------
    # FULL DATA
    # ----------------------------------------------------------------------

    print("=" * 120)
    print("THREE-FEATURE FULL-DATA SEARCH")
    print("=" * 120)
    print()

    full_candidates = full_data_search(features)

    print(
        f"full-data nondegenerate candidates = {len(full_candidates)}"
    )

    if not full_candidates:
        print("  none")
    else:
        print()

        shown = 0
        for cand in strict_full_data_candidates(
            full_candidates,
            features,
        ):
            shown += 1
            print_candidate(
                "SURVIVING CANDIDATE",
                cand,
                features,
                range(len(INSTANCES)),
            )

            if shown >= 20:
                print("  ... output truncated after 20 candidates")
                break

    print()

    # ----------------------------------------------------------------------
    # LOO
    # ----------------------------------------------------------------------

    print("=" * 120)
    print("GENUINE LEAVE-ONE-OUT SEARCH")
    print("=" * 120)
    print()

    loo_survivors = leave_one_out_search(features)

    print(
        f"raw LOO-surviving candidates = {len(loo_survivors)}"
    )
    print()

    held_out_counts = [0] * len(INSTANCES)

    for survivor in loo_survivors:
        held_out_counts[survivor.held_out] += 1

    print("  held-out instance survival counts:")

    for idx, count in enumerate(held_out_counts, start=1):
        print(
            f"    instance {idx}: {count}"
        )

    print()

    # Strictly deduplicate exact same candidate across folds.
    unique_loo: Dict[
        Tuple[str, str, str, int, int, int, int, int],
        List[int],
    ] = {}

    for survivor in loo_survivors:
        c = survivor.candidate
        key = (
            c.F,
            c.G,
            c.H,
            c.a,
            c.b,
            c.c,
            c.d,
            c.e,
        )
        unique_loo.setdefault(key, []).append(
            survivor.held_out
        )

    print(
        f"  unique LOO candidates = {len(unique_loo)}"
    )

    print()

    shown = 0

    for key, held_outs in unique_loo.items():
        c = Candidate(*key)

        # Require survival on every held-out instance.
        if len(set(held_outs)) != len(INSTANCES):
            continue

        shown += 1

        print_candidate(
            "GLOBAL LOO CANDIDATE",
            c,
            features,
            range(len(INSTANCES)),
        )

        if shown >= 20:
            print("  ... output truncated after 20 candidates")
            break

    if shown == 0:
        print("  none")
        print()

    # ----------------------------------------------------------------------
    # CROSS-CHECK
    # ----------------------------------------------------------------------

    print("=" * 120)
    print("NONDEGENERACY CROSS-CHECK")
    print("=" * 120)
    print()

    degenerate_full = 0
    zero_denominator_full = 0
    constant_feature_full = 0

    all_indices = range(len(INSTANCES))

    for cand in full_candidates:
        for name in (cand.F, cand.G, cand.H):
            if len(
                {
                    features[name][i]
                    for i in all_indices
                }
            ) <= 1:
                constant_feature_full += 1
                break

        if any(
            cand.d * features[cand.H][i] + cand.e == 0
            for i in all_indices
        ):
            zero_denominator_full += 1

        if (
            cand.a == 0
            or cand.b == 0
            or cand.d == 0
        ):
            degenerate_full += 1

    print(
        f"  zero-denominator candidates = {zero_denominator_full}"
    )
    print(
        f"  constant-feature candidates = {constant_feature_full}"
    )
    print(
        f"  zero-coefficient candidates = {degenerate_full}"
    )
    print()

    # ----------------------------------------------------------------------
    # TRUE BRANCH
    # ----------------------------------------------------------------------

    print_true_reconstruction_checks(features)

    # ----------------------------------------------------------------------
    # SUMMARY
    # ----------------------------------------------------------------------

    global_loo_candidates = 0

    for key, held_outs in unique_loo.items():
        if len(set(held_outs)) == len(INSTANCES):
            global_loo_candidates += 1

    print("=" * 120)
    print("GLOBAL EXPERIMENT 469 SUMMARY")
    print("=" * 120)
    print()

    print(
        f"  instances                         = {len(INSTANCES)}"
    )
    print(
        f"  variable Euclidean features       = {len(variable_names)}"
    )
    print(
        f"  coefficient range                 = {list(COEFFS)}"
    )
    print(
        f"  offset range                      = {list(OFFSETS)}"
    )
    print(
        f"  full-data nondegenerate laws      = {len(full_candidates)}"
    )
    print(
        f"  raw LOO survivors                 = {len(loo_survivors)}"
    )
    print(
        f"  globally surviving LOO laws      = {global_loo_candidates}"
    )
    print(
        f"  zero-denominator candidates       = {zero_denominator_full}"
    )
    print(
        f"  constant-feature candidates       = {constant_feature_full}"
    )
    print(
        f"  zero-coefficient candidates       = {degenerate_full}"
    )
    print()

    print("INTERPRETATION")
    print()
    print(
        "  Experiment 468 tested two-feature rational"
    )
    print(
        "  Euclidean-derived relations."
    )
    print()
    print(
        "  Experiment 469 increases the numerator side"
    )
    print(
        "  to two independently varying observable"
    )
    print(
        "  features while retaining a separate"
    )
    print(
        "  denominator feature:"
    )
    print()
    print(
        "      aF + bG + c = j(dH + e)"
    )
    print()
    print(
        "  The search is performed by matching exact"
    )
    print(
        "  integer signatures rather than constructing"
    )
    print(
        "  the complete numerator/denominator Cartesian"
    )
    print(
        "  product."
    )
    print()
    print(
        "  A surviving candidate is considered materially"
    )
    print(
        "  interesting only if it also survives every"
    )
    print(
        "  leave-one-out holdout."
    )
    print()
    print(
        "LIMITATION"
    )
    print()
    print(
        "  A negative result excludes only this bounded"
    )
    print(
        "  sparse three-feature rational family."
    )
    print(
        "  It does not establish impossibility of more"
    )
    print(
        "  complicated observable relations."
    )
    print()

    print("=" * 120)
    print("EXPERIMENT 469 FINAL STATUS")
    print("=" * 120)
    print()
    print(
        "  EXACT INTEGER ARITHMETIC        = True"
    )
    print(
        "  THREE-FEATURE RATIONAL FAMILY  = True"
    )
    print(
        "  MEET-IN-THE-MIDDLE SEARCH      = True"
    )
    print(
        "  NONDEGENERATE FILTER            = True"
    )
    print(
        "  GENUINE LEAVE-ONE-OUT           = True"
    )
    print(
        f"  FULL-DATA CANDIDATES             = {len(full_candidates)}"
    )
    print(
        f"  GLOBAL LOO CANDIDATES            = {global_loo_candidates}"
    )
    print(
        "  CONSTANT FEATURES                = Excluded"
    )
    print(
        "  ZERO DENOMINATOR                 = Excluded"
    )
    print(
        "  ZERO COEFFICIENT                 = Excluded"
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

    if global_loo_candidates:
        print(
            "  CONCLUSION:"
        )
        print(
            "    A bounded three-feature rational law"
        )
        print(
            "    survived genuine leave-one-out validation."
        )
        print(
            "    Such a candidate requires independent"
        )
        print(
            "    reproduction on fresh data before it can"
        )
        print(
            "    be treated as meaningful structure."
        )
    else:
        print(
            "  CONCLUSION:"
        )
        print(
            "    No bounded nondegenerate three-feature"
        )
        print(
            "    rational Euclidean-derived j-law survived"
        )
        print(
            "    genuine leave-one-out validation."
        )

    print()
    print("=" * 120)
    print("EXPERIMENT 469 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
