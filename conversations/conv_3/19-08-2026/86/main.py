#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 462
========================================================================================================================

CROSS-INSTANCE j-LAW / OBSERVED-DATA RELATION AUDIT

QUESTION
  Experiment 461 established that, within a single instance, every
  deterministic function of (N,S,d,r) is constant along the j-fibre.

  A remaining possibility is that the hidden displacement j itself
  follows a common cross-instance law determined by the observed data.

  This experiment therefore asks:

      Is there a simple exact relation

          F(N,S,d,r,j) = 0

  shared by all eight instances?

PURPOSE
  Search for genuinely low-complexity cross-instance laws for j.

  Candidate laws are restricted to small-coefficient exact integer
  relations. This deliberately avoids arbitrary interpolation and
  avoids constructing a relation with enormous fitted coefficients.

RULES
  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no floating point
  true j is used ONLY as post-hoc diagnostic data
  true j is NEVER used to construct a reconstruction bound
  no arbitrary polynomial interpolation
  no large fitted coefficients

CROSS-INSTANCE DATA

  For each instance:

      observed = (N,S,d,r)
      hidden diagnostic = j_true

  with

      j_true = (d0_true-d0_ref)/2.

TEST FAMILIES
  1. Exact affine laws
       a*j + b*N + c*S + d0*d + e*r + f = 0

  2. Exact low-degree monomial laws
       a*j + sum(c_m*m(N,S,d,r)) + c0 = 0

     using selected degree <= 2 monomials.

  3. Small-coefficient search
       coefficients restricted to [-BOUND,BOUND].

  4. Leave-one-out verification:
       any discovered relation must hold for every held-out
       instance, not merely the fitting subset.

  5. Direct candidate diagnostics:
       j versus simple exact scales such as
       N/S, S, sqrt(N), bit length, digit length, gcds.

IMPORTANT
  Finding a relation is NOT treated as proof that it is the
  generative mechanism. It only identifies a candidate law for
  further testing.

  Failure to find a small exact law strengthens the conclusion
  that no simple cross-instance relation for j is visible in
  the current eight-instance dataset.
========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import gcd, isqrt
from typing import Sequence


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

COEFF_BOUND = 4

# Keep the feature spaces deliberately small.
# This prevents arbitrary interpolation.
AFFINE_FEATURE_NAMES = (
    "1",
    "N",
    "S",
    "d",
    "r",
)

QUADRATIC_FEATURE_NAMES = (
    "1",
    "N",
    "S",
    "d",
    "r",
    "N^2",
    "S^2",
    "d^2",
    "N*S",
    "N*d",
    "S*d",
    "N*r",
    "S*r",
    "d*r",
)

# Leave-one-out is expensive for larger coefficient sets, but these
# small exact searches are intentionally bounded.
MAX_PRINT = 20


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
            raise AssertionError(f"true j non-integral: instance={self.idx}")
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
# FEATURE CONSTRUCTION
# ----------------------------------------------------------------------------------------------------------------------

def affine_features(inst: Instance) -> tuple[int, ...]:
    return (
        1,
        inst.N,
        inst.S,
        inst.d,
        inst.r,
    )


def quadratic_features(inst: Instance) -> tuple[int, ...]:
    N, S, d, r = inst.N, inst.S, inst.d, inst.r

    return (
        1,
        N,
        S,
        d,
        r,
        N * N,
        S * S,
        d * d,
        N * S,
        N * d,
        S * d,
        N * r,
        S * r,
        d * r,
    )


def named_features(inst: Instance) -> dict[str, int]:
    N, S, d, r = inst.N, inst.S, inst.d, inst.r

    return {
        "1": 1,
        "N": N,
        "S": S,
        "d": d,
        "r": r,
        "N^2": N * N,
        "S^2": S * S,
        "d^2": d * d,
        "N*S": N * S,
        "N*d": N * d,
        "S*d": S * d,
        "N*r": N * r,
        "S*r": S * r,
        "d*r": d * r,
    }


# ----------------------------------------------------------------------------------------------------------------------
# EXACT RELATION TESTING
# ----------------------------------------------------------------------------------------------------------------------

def relation_value(
    coeffs: Sequence[int],
    features: Sequence[int],
    j: int,
    j_coeff: int,
) -> int:
    total = j_coeff * j

    for c, f in zip(coeffs, features):
        total += c * f

    return total


def relation_holds_all(
    instances: Sequence[Instance],
    feature_fn,
    coeffs: Sequence[int],
    j_coeff: int,
) -> bool:
    for inst in instances:
        features = feature_fn(inst)
        value = relation_value(coeffs, features, inst.true_j, j_coeff)

        if value != 0:
            return False

    return True


def relation_holds_subset(
    instances: Sequence[Instance],
    feature_fn,
    coeffs: Sequence[int],
    j_coeff: int,
) -> bool:
    return relation_holds_all(instances, feature_fn, coeffs, j_coeff)


def leave_one_out_holds(
    instances: Sequence[Instance],
    feature_fn,
    coeffs: Sequence[int],
    j_coeff: int,
) -> tuple[bool, list[int]]:
    failures: list[int] = []

    for excluded in instances:
        subset = [x for x in instances if x.idx != excluded.idx]

        # The relation is checked directly on the held-out instance.
        train_ok = relation_holds_subset(
            subset,
            feature_fn,
            coeffs,
            j_coeff,
        )

        held_value = relation_value(
            coeffs,
            feature_fn(excluded),
            excluded.true_j,
            j_coeff,
        )

        # Both are required.
        if not train_ok or held_value != 0:
            failures.append(excluded.idx)

    return len(failures) == 0, failures


# ----------------------------------------------------------------------------------------------------------------------
# SMALL-COEFFICIENT SEARCH
# ----------------------------------------------------------------------------------------------------------------------

def search_small_relation(
    instances: Sequence[Instance],
    feature_fn,
    feature_names: Sequence[str],
    bound: int,
    max_results: int = MAX_PRINT,
):
    """
    Search

        j_coeff*j + c0*f0 + ... + ck*fk = 0

    with all coefficients in [-bound,bound].

    j_coeff must be nonzero so that the relation genuinely constrains j.
    """

    zero = 0
    results = []

    coefficient_range = range(-bound, bound + 1)

    # Prevent the trivial all-zero relation.
    for j_coeff in coefficient_range:
        if j_coeff == 0:
            continue

        # Search the feature coefficients.
        for coeffs in product(coefficient_range, repeat=len(feature_names)):
            if all(c == 0 for c in coeffs):
                continue

            # Normalize by gcd so equivalent scalar multiples are not
            # reported repeatedly.
            g = abs(j_coeff)
            for c in coeffs:
                g = gcd(g, abs(c))

            if g != 1:
                continue

            # Canonical sign.
            first_nonzero = j_coeff
            if first_nonzero == 0:
                for c in coeffs:
                    if c:
                        first_nonzero = c
                        break

            if first_nonzero < 0:
                continue

            if not relation_holds_all(
                instances,
                feature_fn,
                coeffs,
                j_coeff,
            ):
                continue

            loo_ok, failures = leave_one_out_holds(
                instances,
                feature_fn,
                coeffs,
                j_coeff,
            )

            results.append(
                (
                    j_coeff,
                    tuple(coeffs),
                    loo_ok,
                    tuple(failures),
                )
            )

            if len(results) >= max_results:
                return results

    return results


# ----------------------------------------------------------------------------------------------------------------------
# SPECIAL SIMPLE RELATION TESTS
# ----------------------------------------------------------------------------------------------------------------------

def test_simple_relations(instances: Sequence[Instance]):
    tests = []

    def all_equal(values):
        return len(set(values)) == 1

    # j itself and basic exact scales.
    tests.append(
        (
            "j",
            lambda x: x.true_j,
            False,
        )
    )

    tests.append(
        (
            "j-S",
            lambda x: x.true_j - x.S,
            False,
        )
    )

    tests.append(
        (
            "2j-S",
            lambda x: 2 * x.true_j - x.S,
            False,
        )
    )

    tests.append(
        (
            "j-(d-S)",
            lambda x: x.true_j - (x.d - x.S),
            False,
        )
    )

    tests.append(
        (
            "j-(N-d)",
            lambda x: x.true_j - (x.N - x.d),
            False,
        )
    )

    tests.append(
        (
            "j-floor_sqrt_N",
            lambda x: x.true_j - isqrt(x.N),
            False,
        )
    )

    tests.append(
        (
            "j-mod-S",
            lambda x: x.true_j % x.S,
            False,
        )
    )

    tests.append(
        (
            "j-mod-d",
            lambda x: x.true_j % x.d,
            False,
        )
    )

    tests.append(
        (
            "j-mod-N",
            lambda x: x.true_j % x.N,
            False,
        )
    )

    # Print all exact values for diagnostics.
    rows = []

    for name, fn, _ in tests:
        values = [fn(inst) for inst in instances]
        rows.append((name, values, all_equal(values)))

    return rows


# ----------------------------------------------------------------------------------------------------------------------
# BIT / DIGIT / INTEGER-SCALE DIAGNOSTICS
# ----------------------------------------------------------------------------------------------------------------------

def integer_scale_diagnostics(inst: Instance):
    return {
        "j": inst.true_j,
        "N": inst.N,
        "S": inst.S,
        "d": inst.d,
        "floor_sqrt_N": isqrt(inst.N),
        "bitlen_N": inst.N.bit_length(),
        "bitlen_S": inst.S.bit_length(),
        "bitlen_d": inst.d.bit_length(),
        "digits_N": len(str(abs(inst.N))),
        "digits_S": len(str(abs(inst.S))),
        "digits_d": len(str(abs(inst.d))),
        "gcd(N,S)": gcd(inst.N, inst.S),
        "gcd(N,d)": gcd(inst.N, inst.d),
        "gcd(S,d)": gcd(inst.S, inst.d),
        "gcd(N,r)": gcd(inst.N, abs(inst.r)),
        "gcd(S,r)": gcd(inst.S, abs(inst.r)),
        "gcd(d,r)": gcd(inst.d, abs(inst.r)),
    }


# ----------------------------------------------------------------------------------------------------------------------
# EXACT INTEGER TABLE PRINTING
# ----------------------------------------------------------------------------------------------------------------------

def print_instance_table(inst: Instance):
    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"N={inst.N} S={inst.S} d={inst.d} r={inst.r}"
    )
    print()

    print("POST-HOC j")
    print(f"  true j = {inst.true_j}")
    print()

    print("OBSERVED STRUCTURE")
    print(f"  d = N-2S+1                  : {inst.d == inst.N - 2 * inst.S + 1}")
    print(f"  N-2S+1-d                    : {inst.N - 2 * inst.S + 1 - inst.d}")
    print()

    print("INTEGER-SCALE DIAGNOSTICS")
    diag = integer_scale_diagnostics(inst)

    for name, value in diag.items():
        print(f"  {name:30s} = {value}")

    print()

    print("TRUE FIBRE RECONSTRUCTION")
    d0, x, K = (
        inst.d0_ref + 2 * inst.true_j,
        inst.x0 - 2 * inst.true_j,
        inst.K_ref + inst.d0_ref * inst.true_j + inst.true_j**2,
    )

    print(f"  d0 reconstructed = {d0}")
    print(f"  x reconstructed  = {x}")
    print(f"  K reconstructed  = {K}")
    print(f"  exact d0 match   = {d0 == inst.true_d0}")
    print(f"  exact x match    = {x == inst.true_x}")
    print(f"  exact K match    = {K == inst.true_K}")
    print()


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 462")
    print("=" * 120)
    print()
    print("CROSS-INSTANCE j-LAW / OBSERVED-DATA RELATION AUDIT")
    print()
    print("QUESTION")
    print("  Does the hidden displacement j follow a simple exact")
    print("  law determined by the observed tuple across instances?")
    print()
    print(f"COEFFICIENT SEARCH BOUND = {COEFF_BOUND}")
    print()

    # ------------------------------------------------------------------
    # Instance diagnostics
    # ------------------------------------------------------------------

    for inst in INSTANCES:
        print_instance_table(inst)

    # ------------------------------------------------------------------
    # Simple relations
    # ------------------------------------------------------------------

    print("=" * 120)
    print("SIMPLE EXACT SCALE AUDIT")
    print("=" * 120)
    print()

    simple_rows = test_simple_relations(INSTANCES)

    print("  expression                         values")
    print("  " + "-" * 110)

    for name, values, constant in simple_rows:
        print(f"  {name:32s} {values}")
        print(f"    constant across instances = {constant}")

    print()

    # ------------------------------------------------------------------
    # Affine relation search
    # ------------------------------------------------------------------

    print("=" * 120)
    print("SMALL-COEFFICIENT AFFINE SEARCH")
    print("=" * 120)
    print()

    print("Candidate:")
    print()
    print("  a*j + b*N + c*S + d0*d + e*r + f = 0")
    print()
    print(f"coefficient range = [-{COEFF_BOUND},{COEFF_BOUND}]")
    print()

    affine_results = search_small_relation(
        INSTANCES,
        affine_features,
        AFFINE_FEATURE_NAMES,
        COEFF_BOUND,
    )

    print(f"relations found = {len(affine_results)}")

    if affine_results:
        for idx, (j_coeff, coeffs, loo_ok, failures) in enumerate(
            affine_results,
            start=1,
        ):
            print()
            print(f"  RELATION {idx}")

            terms = [f"{j_coeff}*j"]

            for name, coeff in zip(AFFINE_FEATURE_NAMES, coeffs):
                if coeff:
                    terms.append(f"{coeff}*{name}")

            print("    " + " + ".join(terms) + " = 0")
            print(f"    leave-one-out = {loo_ok}")
            print(f"    held-out failures = {failures}")
    else:
        print("  No small-coefficient affine relation found.")

    print()

    # ------------------------------------------------------------------
    # Quadratic relation search
    # ------------------------------------------------------------------

    print("=" * 120)
    print("SMALL-COEFFICIENT QUADRATIC FEATURE SEARCH")
    print("=" * 120)
    print()

    print("Candidate:")
    print()
    print("  a*j + sum(c_i * observable_monomial_i) = 0")
    print()
    print("Observable monomial basis:")
    for name in QUADRATIC_FEATURE_NAMES:
        print(f"  {name}")
    print()
    print(
        "The coefficient search is deliberately bounded and rejects"
        " scalar multiples of the same relation."
    )
    print()

    quadratic_results = search_small_relation(
        INSTANCES,
        quadratic_features,
        QUADRATIC_FEATURE_NAMES,
        1,   # Keep quadratic search conservative.
    )

    print("quadratic coefficient bound = [-1,1]")
    print(f"relations found = {len(quadratic_results)}")

    if quadratic_results:
        for idx, (j_coeff, coeffs, loo_ok, failures) in enumerate(
            quadratic_results,
            start=1,
        ):
            print()
            print(f"  RELATION {idx}")

            terms = [f"{j_coeff}*j"]

            for name, coeff in zip(QUADRATIC_FEATURE_NAMES, coeffs):
                if coeff:
                    terms.append(f"{coeff}*{name}")

            print("    " + " + ".join(terms) + " = 0")
            print(f"    leave-one-out = {loo_ok}")
            print(f"    held-out failures = {failures}")
    else:
        print("  No small-coefficient quadratic relation found.")

    print()

    # ------------------------------------------------------------------
    # Direct cross-instance ratios
    # ------------------------------------------------------------------

    print("=" * 120)
    print("DIRECT j / OBSERVED-SCALE AUDIT")
    print("=" * 120)
    print()

    print(
        "All quantities below use exact integer quotient/remainder pairs."
    )
    print()

    scales = (
        ("S", lambda x: x.S),
        ("d", lambda x: x.d),
        ("N", lambda x: x.N),
        ("N-S", lambda x: x.N - x.S),
        ("N-d", lambda x: x.N - x.d),
        ("d-S", lambda x: x.d - x.S),
        ("floor_sqrt_N", lambda x: isqrt(x.N)),
    )

    for scale_name, scale_fn in scales:
        print(f"  SCALE = {scale_name}")
        print("    idx       j        quotient        remainder")
        print("    " + "-" * 55)

        for inst in INSTANCES:
            scale = scale_fn(inst)

            if scale == 0:
                quotient = "DIV0"
                remainder = "DIV0"
            else:
                quotient = inst.true_j // scale
                remainder = inst.true_j % scale

            print(
                f"    {inst.idx:3d} "
                f"{inst.true_j:9d} "
                f"{str(quotient):15s} "
                f"{str(remainder):15s}"
            )

        print()

    # ------------------------------------------------------------------
    # Cross-instance exact sequence audit
    # ------------------------------------------------------------------

    print("=" * 120)
    print("CROSS-INSTANCE j SEQUENCE AUDIT")
    print("=" * 120)
    print()

    js = [inst.true_j for inst in INSTANCES]
    first_differences = [
        js[i + 1] - js[i]
        for i in range(len(js) - 1)
    ]
    second_differences = [
        first_differences[i + 1] - first_differences[i]
        for i in range(len(first_differences) - 1)
    ]

    print(f"  j sequence              = {js}")
    print(f"  first differences       = {first_differences}")
    print(f"  second differences      = {second_differences}")
    print(
        f"  arithmetic progression  = "
        f"{len(set(first_differences)) == 1}"
    )
    print(
        f"  quadratic-in-index      = "
        f"{len(set(second_differences)) == 1}"
    )
    print()

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 462 SUMMARY")
    print("=" * 120)
    print()

    affine_found = len(affine_results)
    quadratic_found = len(quadratic_results)

    print(f"  instances                         = {len(INSTANCES)}")
    print(f"  j values                          = {js}")
    print()
    print(f"  affine relations found           = {affine_found}")
    print(f"  quadratic relations found        = {quadratic_found}")
    print()
    print(
        f"  affine law search bound          = {COEFF_BOUND}"
    )
    print("  quadratic law search bound       = 1")
    print()
    print("INTERPRETATION")
    print()
    print("  Experiment 461 proved that, inside one instance,")
    print("  every deterministic function of (N,S,d,r) is fixed")
    print("  across the j-fibre.")
    print()
    print("  Experiment 462 asks a different question:")
    print("  whether the hidden j values across the eight test")
    print("  instances obey a simple common law in the observed data.")
    print()
    print("  A relation is interesting only if it is:")
    print("    1. exact,")
    print("    2. low-complexity,")
    print("    3. small-coefficient, and")
    print("    4. survives leave-one-out validation.")
    print()
    print("  Arbitrary interpolation is explicitly excluded.")
    print()
    print("LIMITATION")
    print()
    print("  Failure to find a relation in these bounded families")
    print("  does NOT prove that no more complicated cross-instance")
    print("  relation exists.")
    print()
    print("NEXT INTERPRETATION")
    print()
    print("  If no stable low-complexity law exists, the eight")
    print("  observed tuples provide no evidence that j can be")
    print("  inferred from a simple shared formula.")
    print()
    print(
        "  If a small exact law survives leave-one-out, it becomes"
    )
    print(
        "  a concrete candidate for an independent-data audit."
    )
    print()
    print("=" * 120)
    print("EXPERIMENT 462 FINAL STATUS")
    print("=" * 120)
    print()
    print(f"  EXACT INTEGER ARITHMETIC       = True")
    print(f"  AFFINE SEARCH COMPLETED        = True")
    print(f"  QUADRATIC SEARCH COMPLETED     = True")
    print(f"  AFFINE RELATIONS FOUND         = {affine_found}")
    print(f"  QUADRATIC RELATIONS FOUND      = {quadratic_found}")
    print(f"  GIANT K ENUMERATION            = False")
    print(f"  GIANT d0 ENUMERATION           = False")
    print(f"  CRT CARTESIAN PRODUCT          = False")
    print(f"  FLOATING POINT                 = False")
    print(f"  ARBITRARY INTERPOLATION        = False")
    print()
    print("  CONCLUSION:")
    print("    This experiment tests whether the hidden displacement")
    print("    j has a simple common cross-instance law in the")
    print("    observed data, without treating an interpolating")
    print("    formula as evidence of an actual reconstruction rule.")
    print()
    print("EXPERIMENT 462 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
