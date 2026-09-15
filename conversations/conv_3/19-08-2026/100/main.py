#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 477
========================================================================================================================

FULL EUCLIDEAN / CONTINUED-FRACTION OBSERVABLE j-LAW AUDIT

QUESTION

  Does the hidden displacement j obey a simple exact law derived from
  the COMPLETE Euclidean algorithm data of the observed tuple?

  Observable input:
      (N, S, d, r)

  Euclidean pairs:
      (N,S)
      (N,d)
      (S,d)
      (d,r)
      (N,r)
      (S,r)

  Unlike earlier experiments, this does not use only the first quotient,
  first remainder, or gcd.  The complete Euclidean quotient chain and
  remainder chain are extracted.

SEARCH FAMILY

  1. Exact affine laws
       j = a*F + b*G + c

  2. Exact sparse quadratic laws
       j = a*F*G + b*H + c

  3. Exact rational laws
       a*F + b = j*(c*G + e)

  4. Continued-fraction / Euclidean summary features:
       quotient entries
       remainder entries
       quotient sums
       weighted quotient sums
       remainder sums
       Euclidean step count
       continued-fraction numerator/denominator data
       gcd
       first/last quotients
       alternating quotient sums
       mixed pair summaries

NONDEGENERACY

  Constant features are excluded.
  Zero denominators are excluded.
  Zero coefficients are excluded where appropriate.
  Duplicate feature signatures are removed.

VALIDATION

  Full eight-instance exact fit.
  Genuine leave-one-out validation:
      discover on 7 instances,
      test on the held-out 8th.

RULES

  exact integer arithmetic only
  no floating point
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no arbitrary interpolation

FINAL EXPERIMENT POLICY

  This is intended as the final observable-only audit.
  If no genuine law survives, the script reports that the tested
  observable-information route is exhausted.

========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt
from itertools import combinations
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple


# ======================================================================================================================
# DATA
# ======================================================================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int
    r: int
    true_j: int
    true_d0: int
    true_x: int
    true_K: int


INSTANCES: Tuple[Instance, ...] = (
    Instance(
        1,
        14246098189,
        333010,
        14245432170,
        -13799399066930810668576,
        5,
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
        9,
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
        12,
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
        16,
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
        19,
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
        23,
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
        26,
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
        30,
        270015000131,
        -59,
        -18737381797403407039327539,
    ),
)


# ======================================================================================================================
# EXACT REFERENCE / PARITY
# ======================================================================================================================

def valid_d0_parities(d: int, r: int) -> List[int]:
    """
    d0^2 == d^2-r (mod 4)
    """
    target = (d * d - r) % 4
    out: List[int] = []
    for p in (0, 1):
        if (p * p) % 4 == target:
            out.append(p)
    return out


def derive_reference(inst: Instance) -> Tuple[int, int, int]:
    """
    Canonical known-data reference:
        x0 = parity-compatible representative closest to 0
        d0_ref = d - x0
        K_ref = (r + d0_ref^2 - d^2)/4
    """
    candidates = valid_d0_parities(inst.d, inst.r)
    if not candidates:
        raise AssertionError(f"no valid d0 parity: instance={inst.idx}")

    x0_candidates = [
        x for x in (0, 1)
        if ((inst.d - x) & 1) in candidates
    ]

    if not x0_candidates:
        raise AssertionError(f"no x0 candidate: instance={inst.idx}")

    # Canonical choice used by previous experiments:
    # smallest absolute x0, then smallest x0.
    x0 = min(x0_candidates, key=lambda z: (abs(z), z))
    d0_ref = inst.d - x0

    numerator = inst.r + d0_ref * d0_ref - inst.d * inst.d
    if numerator % 4 != 0:
        raise AssertionError(
            f"reference K not integral: instance={inst.idx}"
        )

    K_ref = numerator // 4
    return x0, d0_ref, K_ref


def reconstruct_true(inst: Instance) -> Tuple[int, int, int]:
    x0, d0_ref, K_ref = derive_reference(inst)

    d0 = d0_ref + 2 * inst.true_j
    x = x0 - 2 * inst.true_j
    K = K_ref + d0_ref * inst.true_j + inst.true_j * inst.true_j

    return d0, x, K


def consistency_failures() -> int:
    failures = 0

    for inst in INSTANCES:
        if inst.d != inst.N - 2 * inst.S + 1:
            failures += 1

        x0, d0_ref, K_ref = derive_reference(inst)

        d0, x, K = reconstruct_true(inst)

        if d0 != inst.true_d0:
            failures += 1
        if x != inst.true_x:
            failures += 1
        if K != inst.true_K:
            failures += 1

        if x + d0 != inst.d:
            failures += 1

        if 4 * K != inst.r + d0 * d0 - inst.d * inst.d:
            failures += 1

        if x0 + d0_ref != inst.d:
            failures += 1

        if 4 * K_ref != inst.r + d0_ref * d0_ref - inst.d * inst.d:
            failures += 1

    return failures


# ======================================================================================================================
# EUCLIDEAN ALGORITHM
# ======================================================================================================================

def euclidean_chain(a: int, b: int) -> Tuple[List[int], List[int]]:
    """
    Exact Euclidean algorithm.

    Returns:
        quotients
        remainders

    Uses floor-style Euclidean division with nonnegative remainder
    for positive b.  For negative b, the signs are normalized first.
    """
    if b == 0:
        return [], [a]

    # Normalize sign of divisor so remainder is nonnegative.
    if b < 0:
        a = -a
        b = -b

    q: List[int] = []
    rem: List[int] = []

    x, y = a, b

    while y != 0:
        quotient = x // y
        remainder = x - quotient * y

        q.append(quotient)
        rem.append(remainder)

        x, y = y, remainder

    return q, rem


def continued_fraction_data(a: int, b: int) -> Dict[str, int]:
    q, rem = euclidean_chain(a, b)

    gcd_value = abs(gcd(a, b))

    q_sum = sum(q)
    q_abs_sum = sum(abs(v) for v in q)
    q_sq_sum = sum(v * v for v in q)

    weighted_q_sum = sum((i + 1) * v for i, v in enumerate(q))
    alternating_q_sum = sum(
        (1 if i % 2 == 0 else -1) * v
        for i, v in enumerate(q)
    )

    rem_sum = sum(rem)
    rem_abs_sum = sum(abs(v) for v in rem)

    weighted_rem_sum = sum((i + 1) * v for i, v in enumerate(rem))
    alternating_rem_sum = sum(
        (1 if i % 2 == 0 else -1) * v
        for i, v in enumerate(rem)
    )

    return {
        "steps": len(q),
        "gcd": gcd_value,
        "q_sum": q_sum,
        "q_abs_sum": q_abs_sum,
        "q_sq_sum": q_sq_sum,
        "q_weighted_sum": weighted_q_sum,
        "q_alternating_sum": alternating_q_sum,
        "rem_sum": rem_sum,
        "rem_abs_sum": rem_abs_sum,
        "rem_weighted_sum": weighted_rem_sum,
        "rem_alternating_sum": alternating_rem_sum,
        "q_first": q[0] if q else 0,
        "q_last": q[-1] if q else 0,
        "rem_first": rem[0] if rem else 0,
        "rem_last_nonzero": next(
            (v for v in reversed(rem) if v != 0),
            0,
        ),
        "q_len": len(q),
        "rem_len": len(rem),
    }


# ======================================================================================================================
# FULL OBSERVABLE FEATURE CONSTRUCTION
# ======================================================================================================================

FeatureValues = Dict[str, List[int]]


def build_features(instances: Sequence[Instance]) -> FeatureValues:
    data: FeatureValues = {}

    def add(name: str, values: List[int]) -> None:
        if name not in data:
            data[name] = values

    # Primitive observables.
    add("N", [i.N for i in instances])
    add("S", [i.S for i in instances])
    add("d", [i.d for i in instances])
    add("r", [i.r for i in instances])

    pairs = {
        "NS": lambda i: (i.N, i.S),
        "Nd": lambda i: (i.N, i.d),
        "Sd": lambda i: (i.S, i.d),
        "dr": lambda i: (i.d, i.r),
        "Nr": lambda i: (i.N, i.r),
        "Sr": lambda i: (i.S, i.r),
    }

    for pair_name, getter in pairs.items():
        qs: List[int] = []
        rs: List[int] = []
        gs: List[int] = []

        for inst in instances:
            a, b = getter(inst)

            # Explicit Euclidean division.
            if b == 0:
                raise AssertionError(
                    f"zero divisor in Euclidean pair {pair_name}, "
                    f"instance={inst.idx}"
                )

            q, rem = euclidean_chain(a, b)

            # First-level q/r.
            qs.append(q[0] if q else 0)
            rs.append(rem[0] if rem else 0)
            gs.append(abs(gcd(a, b)))

        add(f"q({pair_name})", qs)
        add(f"rem({pair_name})", rs)
        add(f"gcd({pair_name})", gs)

        summary_names = [
            "steps",
            "q_sum",
            "q_abs_sum",
            "q_sq_sum",
            "q_weighted_sum",
            "q_alternating_sum",
            "rem_sum",
            "rem_abs_sum",
            "rem_weighted_sum",
            "rem_alternating_sum",
            "q_first",
            "q_last",
            "rem_first",
            "rem_last_nonzero",
            "q_len",
            "rem_len",
        ]

        summaries: Dict[str, List[int]] = {name: [] for name in summary_names}

        for inst in instances:
            a, b = getter(inst)
            cd = continued_fraction_data(a, b)
            for name in summary_names:
                summaries[name].append(cd[name])

        for name in summary_names:
            add(f"{name}({pair_name})", summaries[name])

    # Cross-pair exact arithmetic features.
    def combine(
        name: str,
        f: Sequence[int],
        g: Sequence[int],
        op: Callable[[int, int], int],
    ) -> None:
        add(name, [op(a, b) for a, b in zip(f, g)])

    combine(
        "q(NS)-q(Sd)",
        data["q(NS)"],
        data["q(Sd)"],
        lambda a, b: a - b,
    )

    combine(
        "rem(NS)-rem(Nd)",
        data["rem(NS)"],
        data["rem(Nd)"],
        lambda a, b: a - b,
    )

    combine(
        "rem(Nd)-rem(Sd)",
        data["rem(Nd)"],
        data["rem(Sd)"],
        lambda a, b: a - b,
    )

    combine(
        "gcd(Sd)-gcd(dr)",
        data["gcd(Sd)"],
        data["gcd(dr)"],
        lambda a, b: a - b,
    )

    combine(
        "steps(NS)-steps(Nd)",
        data["steps(NS)"],
        data["steps(Nd)"],
        lambda a, b: a - b,
    )

    combine(
        "qsum(NS)-qsum(Nd)",
        data["q_sum(NS)"],
        data["q_sum(Nd)"],
        lambda a, b: a - b,
    )

    combine(
        "remabs(NS)-remabs(Nd)",
        data["rem_abs_sum(NS)"],
        data["rem_abs_sum(Nd)"],
        lambda a, b: a - b,
    )

    return data


# ======================================================================================================================
# FEATURE SCREEN
# ======================================================================================================================

def print_feature_screen(features: FeatureValues) -> List[str]:
    names = sorted(features)

    variable: List[str] = []

    for name in names:
        vals = features[name]
        distinct = len(set(vals))
        cls = "VARIABLE" if distinct > 1 else "CONSTANT"

        print(
            f"  {name:<38} distinct={distinct:<3} {cls}"
        )

        if distinct > 1:
            variable.append(name)

    return variable


# ======================================================================================================================
# SIGNATURE UTILITIES
# ======================================================================================================================

def signature(
    values: Sequence[int],
    indices: Sequence[int],
) -> Tuple[int, ...]:
    return tuple(values[i] for i in indices)


def is_constant_on(values: Sequence[int], indices: Sequence[int]) -> bool:
    return len({values[i] for i in indices}) <= 1


# ======================================================================================================================
# SEARCH HELPERS
# ======================================================================================================================

COEFFS = (-2, -1, 1, 2)
OFFSETS = tuple(range(-4, 5))


@dataclass(frozen=True)
class AffineCandidate:
    F: str
    G: str
    a: int
    b: int
    c: int

    def text(self) -> str:
        return (
            f"{self.a}*{self.F} + "
            f"{self.b}*{self.G} + {self.c} = j"
        )


@dataclass(frozen=True)
class QuadraticCandidate:
    F: str
    G: str
    H: str
    a: int
    b: int
    c: int
    d: int
    e: int

    def text(self) -> str:
        return (
            f"{self.a}*{self.F}*{self.G} + "
            f"{self.b}*{self.H} + {self.c} = "
            f"j*({self.d}*{self.H} + {self.e})"
        )


@dataclass(frozen=True)
class RationalCandidate:
    F: str
    G: str
    a: int
    b: int
    c: int
    e: int

    def text(self) -> str:
        return (
            f"{self.a}*{self.F} + {self.b} = "
            f"j*({self.c}*{self.G} + {self.e})"
        )


# ======================================================================================================================
# AFFINE SEARCH
# ======================================================================================================================

def search_affine(
    features: FeatureValues,
    variable_features: Sequence[str],
    train_indices: Sequence[int],
) -> List[AffineCandidate]:
    found: List[AffineCandidate] = []

    # Hash exact signatures of aF+bG.
    target = tuple(
        INSTANCES[i].true_j
        for i in train_indices
    )

    seen = set()

    for F, G in combinations(variable_features, 2):
        fvals = features[F]
        gvals = features[G]

        for a in COEFFS:
            for b in COEFFS:
                sig = tuple(
                    a * fvals[i] + b * gvals[i]
                    for i in train_indices
                )

                # j = expression + c
                delta = tuple(
                    target[k] - sig[k]
                    for k in range(len(train_indices))
                )

                if len(set(delta)) == 1:
                    c = delta[0]
                    if c in OFFSETS:
                        cand = AffineCandidate(F, G, a, b, c)
                        if cand not in seen:
                            seen.add(cand)
                            found.append(cand)

    return found


# ======================================================================================================================
# RATIONAL SEARCH
# ======================================================================================================================

def search_rational(
    features: FeatureValues,
    variable_features: Sequence[str],
    train_indices: Sequence[int],
) -> List[RationalCandidate]:
    """
    Searches:
        aF+b = j*(cG+e)

    Nondegenerate:
      F variable on training set
      G variable on training set
      c != 0
      denominator != 0 at every training point
    """
    found: List[RationalCandidate] = []
    seen = set()

    jvals = [INSTANCES[i].true_j for i in train_indices]

    for F, G in combinations(variable_features, 2):
        fvals = features[F]
        gvals = features[G]

        for a in COEFFS:
            for b in OFFSETS:
                lhs = tuple(
                    a * fvals[i] + b
                    for i in train_indices
                )

                for c in COEFFS:
                    if c == 0:
                        continue

                    for e in OFFSETS:
                        den = tuple(
                            c * gvals[i] + e
                            for i in train_indices
                        )

                        if any(v == 0 for v in den):
                            continue

                        ok = all(
                            lhs[k] == jvals[k] * den[k]
                            for k in range(len(train_indices))
                        )

                        if ok:
                            cand = RationalCandidate(
                                F, G, a, b, c, e
                            )
                            if cand not in seen:
                                seen.add(cand)
                                found.append(cand)

    return found


# ======================================================================================================================
# QUADRATIC SEARCH
# ======================================================================================================================

def search_quadratic(
    features: FeatureValues,
    variable_features: Sequence[str],
    train_indices: Sequence[int],
) -> List[QuadraticCandidate]:
    """
    Searches:
        a*F*G + b*H + c = j*(d*Q + e)

    This is deliberately sparse.
    """
    found: List[QuadraticCandidate] = []
    seen = set()

    jvals = [INSTANCES[i].true_j for i in train_indices]

    # Keep the feature pair ordered so F*G and G*F
    # are deduplicated.
    for F_pos in range(len(variable_features)):
        F = variable_features[F_pos]

        for G_pos in range(F_pos + 1, len(variable_features)):
            G = variable_features[G_pos]

            FG = [
                features[F][i] * features[G][i]
                for i in range(len(INSTANCES))
            ]

            for H in variable_features:
                if H in (F, G):
                    continue

                hvals = features[H]

                for Q in variable_features:
                    if Q in (F, G, H):
                        continue

                    qvals = features[Q]

                    for a in COEFFS:
                        for b in COEFFS:
                            for c in OFFSETS:
                                lhs = tuple(
                                    a * FG[i]
                                    + b * hvals[i]
                                    + c
                                    for i in train_indices
                                )

                                for d in COEFFS:
                                    den0 = tuple(
                                        d * qvals[i]
                                        for i in train_indices
                                    )

                                    for e in OFFSETS:
                                        den = tuple(
                                            den0[k] + e
                                            for k in range(len(train_indices))
                                        )

                                        if any(v == 0 for v in den):
                                            continue

                                        if all(
                                            lhs[k] == jvals[k] * den[k]
                                            for k in range(len(train_indices))
                                        ):
                                            cand = QuadraticCandidate(
                                                F,
                                                G,
                                                H,
                                                a,
                                                b,
                                                c,
                                                d,
                                                e,
                                            )

                                            if cand not in seen:
                                                seen.add(cand)
                                                found.append(cand)

    return found


# ======================================================================================================================
# CANDIDATE VALIDATION
# ======================================================================================================================

def validate_affine(
    cand: AffineCandidate,
    features: FeatureValues,
    indices: Sequence[int],
) -> bool:
    vals = [
        cand.a * features[cand.F][i]
        + cand.b * features[cand.G][i]
        + cand.c
        for i in indices
    ]

    target = [INSTANCES[i].true_j for i in indices]
    return vals == target


def validate_rational(
    cand: RationalCandidate,
    features: FeatureValues,
    indices: Sequence[int],
) -> bool:
    for i in indices:
        lhs = (
            cand.a * features[cand.F][i]
            + cand.b
        )

        den = (
            cand.c * features[cand.G][i]
            + cand.e
        )

        if den == 0:
            return False

        if lhs != INSTANCES[i].true_j * den:
            return False

    return True


def validate_quadratic(
    cand: QuadraticCandidate,
    features: FeatureValues,
    indices: Sequence[int],
) -> bool:
    for i in indices:
        lhs = (
            cand.a
            * features[cand.F][i]
            * features[cand.G][i]
            + cand.b * features[cand.H][i]
            + cand.c
        )

        den = (
            cand.d * features[cand.H][i]
            + cand.e
        )

        if den == 0:
            return False

        if lhs != INSTANCES[i].true_j * den:
            return False

    return True


# ======================================================================================================================
# LEAVE-ONE-OUT
# ======================================================================================================================

def loo_affine(
    features: FeatureValues,
    variable_features: Sequence[str],
) -> List[Tuple[AffineCandidate, int]]:
    out: List[Tuple[AffineCandidate, int]] = []

    for held_out in range(len(INSTANCES)):
        train = [
            i for i in range(len(INSTANCES))
            if i != held_out
        ]

        candidates = search_affine(
            features,
            variable_features,
            train,
        )

        for cand in candidates:
            if validate_affine(
                cand,
                features,
                [held_out],
            ):
                out.append((cand, held_out))

    return out


def loo_rational(
    features: FeatureValues,
    variable_features: Sequence[str],
) -> List[Tuple[RationalCandidate, int]]:
    out: List[Tuple[RationalCandidate, int]] = []

    for held_out in range(len(INSTANCES)):
        train = [
            i for i in range(len(INSTANCES))
            if i != held_out
        ]

        candidates = search_rational(
            features,
            variable_features,
            train,
        )

        for cand in candidates:
            if validate_rational(
                cand,
                features,
                [held_out],
            ):
                out.append((cand, held_out))

    return out


def loo_quadratic(
    features: FeatureValues,
    variable_features: Sequence[str],
) -> List[Tuple[QuadraticCandidate, int]]:
    out: List[Tuple[QuadraticCandidate, int]] = []

    for held_out in range(len(INSTANCES)):
        train = [
            i for i in range(len(INSTANCES))
            if i != held_out
        ]

        candidates = search_quadratic(
            features,
            variable_features,
            train,
        )

        for cand in candidates:
            if validate_quadratic(
                cand,
                features,
                [held_out],
            ):
                out.append((cand, held_out))

    return out


# ======================================================================================================================
# STRONGER FULL-CHAIN SIGNATURE SEARCH
# ======================================================================================================================

@dataclass(frozen=True)
class ChainSignatureCandidate:
    """
    Linear relation in a handful of complete Euclidean summary features:

        a1*F1 + a2*F2 + a3*F3 + a4*F4 + c = j

    """
    names: Tuple[str, str, str, str]
    coeffs: Tuple[int, int, int, int]
    offset: int

    def text(self) -> str:
        terms = []
        for coeff, name in zip(self.coeffs, self.names):
            terms.append(f"{coeff}*{name}")

        return " + ".join(terms) + f" + {self.offset} = j"


def search_chain_sparse(
    features: FeatureValues,
    variable_features: Sequence[str],
    train_indices: Sequence[int],
) -> List[ChainSignatureCandidate]:
    """
    Final broad-but-bounded search.

    Uses four distinct variable observable Euclidean features with
    coefficients in [-2,-1,1,2].

    This attacks interactions between full-chain summaries without
    multiplying four-feature combinations against a huge coefficient
    Cartesian product blindly.
    """
    found: List[ChainSignatureCandidate] = []
    seen = set()

    target = tuple(
        INSTANCES[i].true_j
        for i in train_indices
    )

    # Search in signature space.
    # Build all two-feature coefficient signatures first.
    pair_signatures: Dict[
        Tuple[str, str, int, int],
        Tuple[int, ...]
    ] = {}

    for A_idx in range(len(variable_features)):
        A = variable_features[A_idx]
        for B_idx in range(A_idx + 1, len(variable_features)):
            B = variable_features[B_idx]

            for ca in COEFFS:
                for cb in COEFFS:
                    sig = tuple(
                        ca * features[A][i]
                        + cb * features[B][i]
                        for i in train_indices
                    )

                    pair_signatures[(A, B, ca, cb)] = sig

    # Pair + pair = four-feature relation.
    keys = list(pair_signatures)

    for pidx in range(len(keys)):
        A, B, ca, cb = keys[pidx]
        sig1 = pair_signatures[keys[pidx]]

        used1 = {A, B}

        for qidx in range(pidx + 1, len(keys)):
            C, D, cc, cd = keys[qidx]

            if C in used1 or D in used1:
                continue

            sig2 = pair_signatures[keys[qidx]]

            # We require:
            # sig1 + sig2 + offset = j
            delta = tuple(
                target[k] - sig1[k] - sig2[k]
                for k in range(len(train_indices))
            )

            if len(set(delta)) != 1:
                continue

            offset = delta[0]

            if offset not in OFFSETS:
                continue

            cand = ChainSignatureCandidate(
                (A, B, C, D),
                (ca, cb, cc, cd),
                offset,
            )

            if cand not in seen:
                seen.add(cand)
                found.append(cand)

    return found


def validate_chain_sparse(
    cand: ChainSignatureCandidate,
    features: FeatureValues,
    indices: Sequence[int],
) -> bool:
    for i in indices:
        lhs = cand.offset

        for coeff, name in zip(cand.coeffs, cand.names):
            lhs += coeff * features[name][i]

        if lhs != INSTANCES[i].true_j:
            return False

    return True


def loo_chain_sparse(
    features: FeatureValues,
    variable_features: Sequence[str],
) -> List[Tuple[ChainSignatureCandidate, int]]:
    out: List[Tuple[ChainSignatureCandidate, int]] = []

    for held_out in range(len(INSTANCES)):
        train = [
            i for i in range(len(INSTANCES))
            if i != held_out
        ]

        candidates = search_chain_sparse(
            features,
            variable_features,
            train,
        )

        for cand in candidates:
            if validate_chain_sparse(
                cand,
                features,
                [held_out],
            ):
                out.append((cand, held_out))

    return out


# ======================================================================================================================
# RECONSTRUCTION CHECK
# ======================================================================================================================

def print_true_reconstruction() -> int:
    print("========================================================================================================================")
    print("TRUE-BRANCH EXACT RECONSTRUCTION")
    print("========================================================================================================================")

    failures = 0

    for inst in INSTANCES:
        x0, d0_ref, K_ref = derive_reference(inst)

        d0, x, K = reconstruct_true(inst)

        ok = (
            d0 == inst.true_d0
            and x == inst.true_x
            and K == inst.true_K
        )

        print(
            f"  instance={inst.idx} "
            f"j={inst.true_j:<3} "
            f"x0={x0:<2} "
            f"d0_ref={d0_ref} "
            f"d0={d0} "
            f"x={x} "
            f"K={K} "
            f"exact={ok}"
        )

        if not ok:
            failures += 1

    print()
    print(f"  reconstruction failures = {failures}")
    print()

    return failures


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 477")
    print("=" * 120)
    print()
    print("FULL EUCLIDEAN / CONTINUED-FRACTION OBSERVABLE j-LAW AUDIT")
    print()
    print("SEARCH FAMILY")
    print("  complete Euclidean quotient/remainder chains")
    print("  continued-fraction summary features")
    print("  exact affine / rational / sparse nonlinear laws")
    print("  genuine leave-one-out validation")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no factorization")
    print("  no giant K enumeration")
    print("  no giant d0 enumeration")
    print("  no CRT Cartesian product")
    print("  no arbitrary interpolation")
    print()

    failures = consistency_failures()

    print("INPUT CONSISTENCY")
    print("  --------------------------------------------------------------------------------")
    print(f"  consistency failures = {failures}")
    print()

    if failures:
        raise AssertionError(
            f"input consistency failed: {failures}"
        )

    # --------------------------------------------------------------------------------------------------------------
    # REFERENCE AUDIT
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("REFERENCE / PARITY AUDIT")
    print("=" * 120)
    print()

    for inst in INSTANCES:
        parities = valid_d0_parities(inst.d, inst.r)
        x0, d0_ref, K_ref = derive_reference(inst)

        print(f"INSTANCE {inst.idx}")
        print(f"  candidate d0 parities = {parities}")
        print(
            f"  x0={x0} "
            f"d0_ref={d0_ref} "
            f"K_ref={K_ref}"
        )
        print()

    # --------------------------------------------------------------------------------------------------------------
    # FEATURE BUILD
    # --------------------------------------------------------------------------------------------------------------

    features = build_features(INSTANCES)

    print("=" * 120)
    print("FULL EUCLIDEAN FEATURE SCREEN")
    print("=" * 120)
    print()

    variable_features = print_feature_screen(features)

    print()
    print(f"  total features    = {len(features)}")
    print(f"  variable features = {len(variable_features)}")
    print()

    if len(variable_features) < 4:
        raise AssertionError(
            "Not enough variable features for final sparse search"
        )

    # --------------------------------------------------------------------------------------------------------------
    # AFFINE SEARCH
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("FULL EUCLIDEAN AFFINE SEARCH")
    print("=" * 120)
    print()

    all_indices = list(range(len(INSTANCES)))

    affine = search_affine(
        features,
        variable_features,
        all_indices,
    )

    print(f"full-data affine candidates = {len(affine)}")

    if affine:
        for cand in affine[:20]:
            print(f"  {cand.text()}")
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # RATIONAL SEARCH
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("FULL EUCLIDEAN RATIONAL SEARCH")
    print("=" * 120)
    print()

    rational = search_rational(
        features,
        variable_features,
        all_indices,
    )

    print(f"full-data rational candidates = {len(rational)}")

    if rational:
        for cand in rational[:20]:
            print(f"  {cand.text()}")
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # QUADRATIC SEARCH
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("FULL EUCLIDEAN SPARSE-QUADRATIC SEARCH")
    print("=" * 120)
    print()

    quadratic = search_quadratic(
        features,
        variable_features,
        all_indices,
    )

    print(
        f"full-data sparse-quadratic candidates = {len(quadratic)}"
    )

    if quadratic:
        for cand in quadratic[:20]:
            print(f"  {cand.text()}")
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # FOUR-FEATURE FULL-CHAIN SEARCH
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("FULL EUCLIDEAN FOUR-FEATURE SIGNATURE SEARCH")
    print("=" * 120)
    print()

    chain = search_chain_sparse(
        features,
        variable_features,
        all_indices,
    )

    print(
        f"full-data four-feature candidates = {len(chain)}"
    )

    if chain:
        for cand in chain[:20]:
            print(f"  {cand.text()}")
    else:
        print("  none")

    print()

    # --------------------------------------------------------------------------------------------------------------
    # GENUINE LOO
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GENUINE LEAVE-ONE-OUT VALIDATION")
    print("=" * 120)
    print()

    loo_aff = loo_affine(
        features,
        variable_features,
    )

    loo_rat = loo_rational(
        features,
        variable_features,
    )

    loo_quad = loo_quadratic(
        features,
        variable_features,
    )

    loo_chain = loo_chain_sparse(
        features,
        variable_features,
    )

    print(
        f"  affine LOO matches      = {len(loo_aff)}"
    )
    print(
        f"  rational LOO matches    = {len(loo_rat)}"
    )
    print(
        f"  quadratic LOO matches   = {len(loo_quad)}"
    )
    print(
        f"  four-feature LOO matches= {len(loo_chain)}"
    )
    print()

    # Unique candidates that predict all held-out instances.
    def unique_all(
        entries,
    ):
        by_candidate: Dict[object, set] = {}

        for cand, held_out in entries:
            by_candidate.setdefault(cand, set()).add(held_out)

        return [
            cand
            for cand, heldouts in by_candidate.items()
            if heldouts == set(range(len(INSTANCES)))
        ]

    loo_aff_unique = unique_all(loo_aff)
    loo_rat_unique = unique_all(loo_rat)
    loo_quad_unique = unique_all(loo_quad)
    loo_chain_unique = unique_all(loo_chain)

    print(
        f"  unique globally LOO affine       = {len(loo_aff_unique)}"
    )
    print(
        f"  unique globally LOO rational     = {len(loo_rat_unique)}"
    )
    print(
        f"  unique globally LOO quadratic    = {len(loo_quad_unique)}"
    )
    print(
        f"  unique globally LOO four-feature = {len(loo_chain_unique)}"
    )
    print()

    for title, candidates in (
        ("AFFINE", loo_aff_unique),
        ("RATIONAL", loo_rat_unique),
        ("QUADRATIC", loo_quad_unique),
        ("FOUR-FEATURE", loo_chain_unique),
    ):
        if candidates:
            print(f"{title} SURVIVORS")
            for cand in candidates[:20]:
                print(f"  {cand.text()}")
            print()

    # --------------------------------------------------------------------------------------------------------------
    # RECONSTRUCTION
    # --------------------------------------------------------------------------------------------------------------

    reconstruction_failures = print_true_reconstruction()

    # --------------------------------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 477 SUMMARY")
    print("=" * 120)
    print()
    print(f"  instances                         = {len(INSTANCES)}")
    print(f"  total Euclidean features          = {len(features)}")
    print(f"  variable Euclidean features       = {len(variable_features)}")
    print(f"  affine full-data candidates       = {len(affine)}")
    print(f"  rational full-data candidates     = {len(rational)}")
    print(f"  quadratic full-data candidates    = {len(quadratic)}")
    print(f"  four-feature full-data candidates = {len(chain)}")
    print(f"  globally LOO affine laws          = {len(loo_aff_unique)}")
    print(f"  globally LOO rational laws        = {len(loo_rat_unique)}")
    print(f"  globally LOO quadratic laws       = {len(loo_quad_unique)}")
    print(
        f"  globally LOO four-feature laws    = "
        f"{len(loo_chain_unique)}"
    )
    print(
        f"  reconstruction failures           = "
        f"{reconstruction_failures}"
    )
    print()

    any_survivor = any(
        (
            loo_aff_unique,
            loo_rat_unique,
            loo_quad_unique,
            loo_chain_unique,
        )
    )

    print("INTERPRETATION")
    print()
    print("  Previous experiments searched isolated Euclidean features.")
    print()
    print("  This experiment uses the COMPLETE Euclidean/continued-fraction")
    print("  structure of observable pairs:")
    print()
    print("      quotient chain")
    print("      remainder chain")
    print("      chain length")
    print("      weighted quotient sums")
    print("      alternating quotient sums")
    print("      remainder aggregates")
    print("      gcd")
    print("      first/last terms")
    print()
    print("  The final sparse search combines four independently varying")
    print("  Euclidean-derived observable features.")
    print()
    print("  A candidate is considered meaningful only if it predicts the")
    print("  held-out instance after being discovered without it.")
    print()

    if any_survivor:
        print("PRIMARY OUTCOME")
        print()
        print("  A candidate survived the final bounded observable audit.")
        print("  It should be investigated independently before drawing")
        print("  any conclusion about recoverability.")
    else:
        print("PRIMARY OUTCOME")
        print()
        print("  No candidate survived the final bounded observable audit.")
        print()
        print("  Combined with Experiments 461-476, this means the tested")
        print("  observable-only route has failed across:")
        print()
        print("      polynomial relations")
        print("      modular relations")
        print("      quotient/floor relations")
        print("      Euclidean-derived relations")
        print("      rational relations")
        print("      sparse quadratic relations")
        print("      mixed quadratic relations")
        print("      cubic relations")
        print("      repeated-feature cubic relations")
        print("      complete Euclidean-chain summaries")
        print()
        print("  This does NOT constitute a proof that every imaginable")
        print("  observable law is impossible.")
        print("  It does provide a strong stopping point for this research")
        print("  branch.")

    print()
    print("=" * 120)
    print("EXPERIMENT 477 FINAL STATUS")
    print("=" * 120)
    print()
    print(f"  EXACT INTEGER ARITHMETIC        = {True}")
    print(f"  COMPLETE EUCLIDEAN CHAINS       = {True}")
    print(f"  CONTINUED-FRACTION SUMMARIES    = {True}")
    print(f"  AFFINE SEARCH                   = {True}")
    print(f"  RATIONAL SEARCH                 = {True}")
    print(f"  SPARSE QUADRATIC SEARCH         = {True}")
    print(f"  FOUR-FEATURE SIGNATURE SEARCH   = {True}")
    print(f"  GENUINE LEAVE-ONE-OUT           = {True}")
    print(f"  FLOATING POINT                  = {False}")
    print(f"  ARBITRARY INTERPOLATION         = {False}")
    print(f"  GIANT K ENUMERATION             = {False}")
    print(f"  GIANT d0 ENUMERATION            = {False}")
    print(f"  CRT CARTESIAN PRODUCT           = {False}")
    print(f"  RECONSTRUCTION FAILURES         = {reconstruction_failures}")
    print()

    if any_survivor:
        print("  CONCLUSION:")
        print("    A bounded exact observable law survived the final")
        print("    leave-one-out audit and requires independent analysis.")
    else:
        print("  CONCLUSION:")
        print("    No bounded observable law based on the full Euclidean/")
        print("    continued-fraction structure survived genuine")
        print("    leave-one-out validation.")
        print()
        print("    This is a reasonable stopping point for the")
        print("    current j-recovery approach.")

    print()
    print("EXPERIMENT 477 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
