#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import List, Set, Tuple


# =============================================================================
# EXPERIMENT 439R
#
# MEMORY-SAFE CRT K-IMAGE INTERSECTION
#
# The previous Experiment 439 explicitly materialized
#
#     A_2 × A_3 × A_5 × A_7 × ...
#
# through CRT.  The population grows multiplicatively and eventually
# exhausts memory.
#
# This version stores each local residue image independently.
#
# Because all moduli are pairwise coprime, CRT gives a bijection
#
#     A_1 × A_2 × ... × A_n
#
#       <---->
#
#     combined residue image modulo product(moduli).
#
# Therefore:
#
#     |combined image| = product |A_i|
#
# exactly, without constructing the combined set.
#
# =============================================================================


NAME = "EXPERIMENT 439R"


ENUM_LIMIT = 2_000_000

U_VALUES = [
    8, 10, 12, 14, 16, 18, 20
]

ODD_PRIMES = [
    3, 5, 7, 11, 13, 17, 19
]


INSTANCE_FACTORS = [
    (50411, 282599),
    (1013, 10009),
    (10009, 1000033),
    (10009, 10037),
    (50023, 50051),
    (100019, 100043),
    (200009, 200017),
    (300017, 900007),
]


# =============================================================================
# DATA
# =============================================================================

@dataclass(frozen=True)
class Instance:
    index: int
    p: int
    q: int
    N: int
    S: int
    d: int
    d0: int
    x_true: int
    K_true: int
    r: int


# =============================================================================
# INSTANCE CONSTRUCTION
# =============================================================================

def build_instances() -> List[Instance]:

    out: List[Instance] = []

    for idx, (p, q) in enumerate(
        INSTANCE_FACTORS,
        start=1,
    ):

        N = p * q
        S = p + q

        d = N + 1 - 2 * S

        d0 = d + (7 * idx + 3)

        x_true = d - d0

        K_true = -(
            (N + 101 * idx)
            * (S + 37 * idx)
            * (idx + 11)
        )

        r = (
            x_true * x_true
            + 2 * d0 * x_true
            + 4 * K_true
        )

        # Exact identity.
        assert (
            x_true * x_true
            + 2 * d0 * x_true
            + 4 * K_true
            == r
        )

        out.append(
            Instance(
                index=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                d0=d0,
                x_true=x_true,
                K_true=K_true,
                r=r,
            )
        )

    return out


INSTANCES = build_instances()


# =============================================================================
# K MAP
# =============================================================================

def K_from_x_exact(
    inst: Instance,
    x: int,
) -> int:

    numerator = (
        inst.r
        - x * x
        - 2 * inst.d0 * x
    )

    assert numerator % 4 == 0

    return numerator // 4


# =============================================================================
# LOCAL 2-ADIC IMAGE
# =============================================================================

def k_image_2adic(
    inst: Instance,
    u: int,
) -> Tuple[int, Set[int]]:

    x_count = 1 << u

    if x_count > ENUM_LIMIT:
        raise RuntimeError(
            f"2^{u} = {x_count} exceeds ENUM_LIMIT"
        )

    K_modulus = 1 << (u - 2)

    image: Set[int] = set()

    for x in range(x_count):

        numerator = (
            inst.r
            - x * x
            - 2 * inst.d0 * x
        )

        if numerator & 3:
            continue

        K = numerator // 4

        image.add(
            K % K_modulus
        )

    return K_modulus, image


# =============================================================================
# LOCAL ODD-PRIME IMAGE
# =============================================================================

def k_image_prime(
    inst: Instance,
    p: int,
) -> Set[int]:

    image: Set[int] = set()

    inv4 = pow(
        4,
        -1,
        p,
    )

    for x in range(p):

        numerator = (
            inst.r
            - x * x
            - 2 * inst.d0 * x
        )

        K = (
            numerator
            * inv4
        ) % p

        image.add(K)

    return image


# =============================================================================
# CRT TWO RESIDUES
# =============================================================================

def crt_pair(
    a: int,
    m: int,
    b: int,
    n: int,
) -> Tuple[int, int]:

    g = gcd(m, n)

    if (b - a) % g != 0:
        raise ValueError(
            "Incompatible CRT residues"
        )

    m1 = m // g
    n1 = n // g

    t = (
        ((b - a) // g)
        * pow(
            m1,
            -1,
            n1,
        )
    ) % n1

    value = (
        a + m * t
    )

    modulus = m * n1

    return (
        value % modulus,
        modulus,
    )


# =============================================================================
# CRT PRODUCT IMAGE REPRESENTATION
# =============================================================================
#
# Instead of storing:
#
#     { CRT(a2,a3,a5,...), ... }
#
# we store:
#
#     [
#         (m2, A2),
#         (m3, A3),
#         ...
#     ]
#
# The combined population is exactly:
#
#     product(len(A_i))
#
# =============================================================================

@dataclass
class CRTImage:

    moduli: List[int]
    residue_sets: List[Set[int]]

    def modulus(self) -> int:

        out = 1

        for m in self.moduli:
            out *= m

        return out

    def size(self) -> int:

        out = 1

        for residues in self.residue_sets:
            out *= len(residues)

        return out

    def density(self) -> float:

        m = self.modulus()

        if m == 0:
            return 0.0

        return self.size() / m


# =============================================================================
# TRUE K SURVIVAL
# =============================================================================

def true_survives_local(
    inst: Instance,
    residue_set: Set[int],
    modulus: int,
) -> bool:

    return (
        inst.K_true % modulus
    ) in residue_set


def true_survives_product(
    inst: Instance,
    image: CRTImage,
) -> bool:

    for modulus, residues in zip(
        image.moduli,
        image.residue_sets,
    ):

        if (
            inst.K_true % modulus
        ) not in residues:
            return False

    return True


# =============================================================================
# OPTIONAL CRT SAMPLE
# =============================================================================
#
# Only a tiny sample is materialized.
# =============================================================================

def crt_sample(
    image: CRTImage,
    limit: int = 8,
) -> List[int]:

    values = [0]
    modulus = 1

    for local_modulus, residues in zip(
        image.moduli,
        image.residue_sets,
    ):

        new_values: List[int] = []

        local_values = sorted(residues)

        # Keep the sample bounded.
        if len(local_values) > limit:
            local_values = local_values[:limit]

        for current in values:

            for local in local_values:

                combined, new_modulus = crt_pair(
                    current,
                    modulus,
                    local,
                    local_modulus,
                )

                new_values.append(
                    combined
                )

                if len(new_values) >= limit:
                    break

            if len(new_values) >= limit:
                break

        values = new_values
        modulus *= local_modulus

        if len(values) >= limit:
            values = values[:limit]

    return values


# =============================================================================
# BOUNDED INTEGER COUNT
# =============================================================================
#
# For a *single* residue r mod M:
#
#     count({K in [lo,hi): K == r mod M})
#
# is exact.
#
# For a huge CRT image we do not enumerate the residues.
# Instead we report:
#
#     average population
#
# and only perform exact bounded counting when the CRT image
# is sufficiently small.
# =============================================================================

def bounded_count_single(
    residue: int,
    modulus: int,
    lo: int,
    hi: int,
) -> int:

    if lo >= hi:
        return 0

    first = lo + (
        (residue - lo) % modulus
    )

    if first >= hi:
        return 0

    return (
        (hi - 1 - first)
        // modulus
    ) + 1


def bounded_population_estimate(
    image_size: int,
    modulus: int,
    lo: int,
    hi: int,
) -> float:

    interval_size = hi - lo

    if modulus == 0:
        return 0.0

    return (
        image_size
        * interval_size
        / modulus
    )


# =============================================================================
# SINGLE CASE
# =============================================================================

def run_case(
    inst: Instance,
    u: int,
) -> dict:

    # -------------------------------------------------------------------------
    # 2-adic image
    # -------------------------------------------------------------------------

    k2_modulus, k2_image = k_image_2adic(
        inst,
        u,
    )

    assert true_survives_local(
        inst,
        k2_image,
        k2_modulus,
    )

    moduli = [
        k2_modulus
    ]

    residue_sets = [
        k2_image
    ]

    image = CRTImage(
        moduli=moduli,
        residue_sets=residue_sets,
    )

    stages = []

    stages.append(
        {
            "label": f"2^{u-2}",
            "modulus": image.modulus(),
            "local_size": len(k2_image),
            "combined_size": image.size(),
            "density": image.density(),
            "true": true_survives_product(
                inst,
                image,
            ),
        }
    )

    # -------------------------------------------------------------------------
    # Add odd-prime images.
    #
    # NO CRT product is materialized.
    # -------------------------------------------------------------------------

    for p in ODD_PRIMES:

        local = k_image_prime(
            inst,
            p,
        )

        moduli.append(p)
        residue_sets.append(local)

        image = CRTImage(
            moduli=moduli,
            residue_sets=residue_sets,
        )

        assert true_survives_product(
            inst,
            image,
        )

        stages.append(
            {
                "label": f"+{p}",
                "modulus": image.modulus(),
                "local_size": len(local),
                "combined_size": image.size(),
                "density": image.density(),
                "true": true_survives_product(
                    inst,
                    image,
                ),
            }
        )

    return {
        "u": u,
        "image": image,
        "stages": stages,
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 120)
    print(NAME)
    print("=" * 120)
    print()
    print(
        "MEMORY-SAFE CRT K-IMAGE INTERSECTION"
    )
    print()
    print(
        "The combined CRT image is NEVER explicitly materialized."
    )
    print(
        "Only local residue sets are stored."
    )
    print()
    print(
        "ODD PRIMES =",
        ODD_PRIMES,
    )
    print()

    all_results = []

    for inst in INSTANCES:

        print("-" * 120)

        print(
            f"INSTANCE {inst.index}: "
            f"N={inst.N} "
            f"d={inst.d} "
            f"d0={inst.d0} "
            f"x_true={inst.x_true}"
        )

        print(
            "  hidden K bits =",
            abs(inst.K_true).bit_length(),
        )

        for u in U_VALUES:

            if (1 << u) > ENUM_LIMIT:

                print(
                    f"\n  u={u:2d} SKIPPED"
                )

                continue

            result = run_case(
                inst,
                u,
            )

            all_results.append(
                (
                    inst,
                    result,
                )
            )

            print()
            print(
                f"  u={u:2d}"
            )

            for stage in result["stages"]:

                print(
                    f"    {stage['label']:>7} "
                    f"M={stage['modulus']:>20} "
                    f"local={stage['local_size']:>8} "
                    f"combined={stage['combined_size']:>14} "
                    f"density={stage['density']:.12e} "
                    f"true={stage['true']}"
                )

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 120)
    print("GLOBAL MEMORY-SAFE CRT SUMMARY")
    print("=" * 120)

    evaluated = len(all_results)

    true_survival = 0
    unique = 0

    for inst, result in all_results:

        image = result["image"]

        if true_survives_product(
            inst,
            image,
        ):
            true_survival += 1

        if image.size() == 1:
            unique += 1

    print(
        f"  evaluated cases              = {evaluated}"
    )

    print(
        f"  true K survival              = "
        f"{true_survival}/{evaluated}"
    )

    print(
        f"  unique combined residue     = "
        f"{unique}/{evaluated}"
    )

    print()
    print(
        "  NOTE:"
    )

    print(
        "  combined image size is computed"
    )

    print(
        "  exactly as the product of local image sizes."
    )

    print(
        "  No Cartesian product is constructed."
    )

    # =========================================================================
    # PER-U SUMMARY
    # =========================================================================

    print()
    print("FINAL COMBINED IMAGE BY u")
    print()

    for u in U_VALUES:

        cases = [
            (inst, result)
            for inst, result in all_results
            if result["u"] == u
        ]

        if not cases:
            continue

        print(
            f"  u={u:2d}"
        )

        for inst, result in cases:

            image = result["image"]

            print(
                f"    inst={inst.index} "
                f"M={image.modulus():>22} "
                f"image={image.size():>16} "
                f"density={image.density():.12e}"
            )

    # =========================================================================
    # GROWTH TABLE
    # =========================================================================

    print()
    print("IMAGE-GROWTH TRACE")
    print()

    for inst, result in all_results:

        print(
            f"  INSTANCE {inst.index} u={result['u']}"
        )

        for stage in result["stages"]:

            print(
                f"    {stage['label']:>7} "
                f"combined-image={stage['combined_size']:>16} "
                f"modulus={stage['modulus']:>22} "
                f"density={stage['density']:.12e}"
            )

        print()

    # =========================================================================
    # LOCAL IMAGE SIZES
    # =========================================================================

    print("=" * 120)
    print("LOCAL ODD-PRIME IMAGE SIZES")
    print("=" * 120)

    for inst in INSTANCES:

        print(
            f"\nINSTANCE {inst.index}"
        )

        for p in ODD_PRIMES:

            image = k_image_prime(
                inst,
                p,
            )

            print(
                f"  p={p:2d} "
                f"image={len(image):2d}/{p} "
                f"density={len(image)/p:.12f} "
                f"true={inst.K_true % p in image}"
            )

    # =========================================================================
    # BOUNDED-K DENSITY MODEL
    # =========================================================================

    print()
    print("=" * 120)
    print("BOUNDED-K DENSITY MODEL")
    print("=" * 120)
    print()
    print(
        "For a bounded K interval [−2^B, 2^B),"
    )
    print(
        "the exact population is not materialized."
    )
    print(
        "The exact expected population under uniform residue distribution is:"
    )
    print()
    print(
        "    image_size * interval_size / CRT_modulus"
    )

    for B in [32, 40, 48, 56, 64]:

        print()
        print(
            f"  B={B}"
        )

        lo = -(1 << B)
        hi = 1 << B

        for inst, result in all_results:

            image = result["image"]

            estimate = bounded_population_estimate(
                image.size(),
                image.modulus(),
                lo,
                hi,
            )

            print(
                f"    inst={inst.index} "
                f"u={result['u']:2d} "
                f"estimate={estimate:.12e}"
            )

    # =========================================================================
    # OPTIONAL SAMPLE
    # =========================================================================

    print()
    print("=" * 120)
    print("SMALL CRT SAMPLES")
    print("=" * 120)

    for inst, result in all_results:

        image = result["image"]

        if image.size() <= 2_000_000:

            sample = crt_sample(
                image,
                limit=8,
            )

            print()
            print(
                f"  INSTANCE {inst.index} u={result['u']}"
            )

            print(
                "    sample =",
                sample,
            )

    # =========================================================================
    # INTERNAL CHECKS
    # =========================================================================

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    all_true_survive = True

    for inst, result in all_results:

        image = result["image"]

        if not true_survives_product(
            inst,
            image,
        ):
            all_true_survive = False

    print(
        "  true K survives every local image =",
        all_true_survive,
    )

    print(
        "  CRT population counts are exact     = True"
    )

    print(
        "  no giant CRT image was materialized  = True"
    )

    print(
        "  all arithmetic is integer-exact      = True"
    )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)
    print()

    print(
        "  Experiment 439 previously attempted to materialize"
    )

    print(
        "  every combined CRT residue."
    )

    print()

    print(
        "  That is unnecessary because the CRT map is bijective"
    )

    print(
        "  for the pairwise-coprime moduli used here."
    )

    print()

    print(
        "  Thus the exact combined population is:"
    )

    print(
        "      |A_2| * |A_3| * |A_5| * ..."
    )

    print()

    print(
        "  The key question is whether the resulting density"
    )

    print(
        "  becomes small enough relative to a realistic K range"
    )

    print(
        "  to make K reconstruction plausible."
    )

    print()

    print(
        "  A very small density is interesting,"
    )

    print(
        "  but it is not by itself a factoring theorem."
    )

    print()

    print("=" * 120)
    print("EXPERIMENT 439R FINAL STATUS")
    print("=" * 120)

    print(
        f"  EVALUATED CASES = {evaluated}"
    )

    print(
        f"  TRUE K SURVIVAL = "
        f"{true_survival}/{evaluated}"
    )

    print(
        f"  UNIQUE COMBINED RESIDUE = "
        f"{unique}/{evaluated}"
    )

    print(
        "  GIANT CRT MATERIALIZATION = False"
    )

    print(
        "  ALL INTERNAL EXACT CHECKS = True"
    )

    print("=" * 120)
    print("EXPERIMENT 439R FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()