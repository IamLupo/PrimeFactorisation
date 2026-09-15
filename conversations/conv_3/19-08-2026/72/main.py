#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import gcd


# ======================================================================================================================
# EXPERIMENT 448
# ======================================================================================================================
#
# EXACT d0 AMBIGUITY / CONGRUENCE-IDENTIFIABILITY AUDIT
#
# QUESTION
#
#   Given only the observed quantities N,S,d,r, does the equation
#
#       4K = r + d0^2 - d^2
#
#   determine d0, or does every parity-compatible d0 produce an
#   equally valid integer K?
#
# CORE EQUATION
#
#       A = r + d0^2
#       4K = A - d^2
#
# therefore
#
#       K(d0) = (r + d0^2 - d^2)/4
#
# For two candidate d0 values:
#
#       d0' = d0 + 2t
#
# the K difference is
#
#       K(d0') - K(d0)
#           = ((d0+2t)^2 - d0^2)/4
#           = t*d0 + t^2
#
# which is always integral.
#
# This means parity-compatible d0 values may form an entire affine
# family of exact integer K solutions.
#
# NO hidden K/d0/x is used during reconstruction.
#
# RULES
#   exact integer arithmetic only
#   no factorization during reconstruction
#   no giant K enumeration
#   no giant X enumeration
#   no CRT Cartesian product
#
# ======================================================================================================================


@dataclass(frozen=True)
class Instance:
    idx: int

    # Construction-only factors.
    p: int
    q: int

    # Observed/reconstruction quantities.
    N: int
    S: int
    d: int
    r: int

    # Hidden construction values.
    d0_true: int
    x_true: int
    K_true: int


RAW = [
    # p, q, x_true, hidden K_true
    (50411, 282599, -10, -3449849766661475506269),
    (1013, 10009, -17, -2968347695488785),
    (10009, 1000033, -24, -4306289434216722346403),
    (10009, 10037, -31, -615138736337991015),
    (50023, 50051, -38, -557809093589551261113),
    (100019, 100043, -45, -12714738365215418549091),
    (200009, 200017, -52, -305667239831581101061223),
    (300017, 900007, -59, -18737381797403407039327539),
]


# ======================================================================================================================
# CONSTRUCTION
# ======================================================================================================================

def build_instances() -> list[Instance]:

    result: list[Instance] = []

    for idx, (p, q, x_true, K_true) in enumerate(RAW, start=1):

        N = p * q
        S = p + q

        d = N - 2 * S + 1

        # Controlled hidden d0.
        d0_true = d - x_true

        # From:
        #
        #   4K = r + d0^2 - d^2
        #
        # solve for the observed r used by this experiment.
        r = 4 * K_true + d * d - d0_true * d0_true

        result.append(
            Instance(
                idx=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                r=r,
                d0_true=d0_true,
                x_true=x_true,
                K_true=K_true,
            )
        )

    return result


# ======================================================================================================================
# OBSERVABLE FUNCTIONS
# ======================================================================================================================

def candidate_k(inst: Instance, d0: int) -> tuple[bool, int]:
    """
    Compute

        K(d0) = (r + d0^2 - d^2)/4.

    The returned boolean tells whether K is integral.
    """

    numerator = inst.r + d0 * d0 - inst.d * inst.d

    if numerator % 4 != 0:
        return False, 0

    return True, numerator // 4


def candidate_x(inst: Instance, d0: int) -> int:
    return inst.d - d0


def parity_required(inst: Instance) -> int:
    """
    Because

        r + d0^2 - d^2 ≡ 0 (mod 4),

    the parity of d0 is constrained by r-d^2.

    This function returns the parity that works.
    """

    for parity in (0, 1):

        test_d0 = parity

        if (inst.r + test_d0 * test_d0 - inst.d * inst.d) % 4 == 0:
            return parity

    raise AssertionError("No parity class satisfies the mod-4 condition")


# ======================================================================================================================
# LOCAL MODULUS AUDIT
# ======================================================================================================================

def valid_residue_classes(
    inst: Instance,
    modulus: int,
) -> list[int]:
    """
    Enumerate d0 residues modulo a SMALL modulus.

    This is not a giant search. It is only used for small local
    moduli in order to determine whether stronger congruences
    distinguish d0.

    A residue a is considered locally valid when

        r + a^2 - d^2 ≡ 0 (mod 4).

    and a is compatible with the requested modulus.

    The integrality condition itself only depends on mod 4,
    but we also record the resulting K residue modulo modulus.
    """

    result = []

    for a in range(modulus):
        numerator = inst.r + a * a - inst.d * inst.d

        if numerator % 4 == 0:
            result.append(a)

    return result


def k_residue_for_d0(
    inst: Instance,
    d0: int,
    modulus: int,
) -> int | None:
    """
    Return K(d0) modulo modulus when K(d0) is integral.
    """

    ok, k = candidate_k(inst, d0)

    if not ok:
        return None

    return k % modulus


# ======================================================================================================================
# WINDOW AUDIT
# ======================================================================================================================

def scan_d0_window(
    inst: Instance,
    radius: int,
) -> list[tuple[int, int, int]]:
    """
    Enumerate only a small diagnostic d0 window around the true hidden
    value.

    The true d0 is used ONLY by the test harness to place the diagnostic
    window. The reconstruction logic itself does not use it to construct
    K.

    Each result is:

        (d0_candidate, x_candidate, K_candidate)
    """

    out = []

    # Preserve only parity-compatible candidates.
    start = inst.d0_true - radius
    stop = inst.d0_true + radius

    for d0 in range(start, stop + 1):

        ok, K = candidate_k(inst, d0)

        if not ok:
            continue

        x = candidate_x(inst, d0)

        out.append((d0, x, K))

    return out


# ======================================================================================================================
# AFFINE FAMILY CHECK
# ======================================================================================================================

def affine_family_check(
    inst: Instance,
    t_values: list[int],
) -> bool:
    """
    Check

        d0' = d0 + 2t

        K' = K + t*d0 + t^2

    exactly.

    The true values appear only as a post-hoc algebraic check.
    """

    for t in t_values:

        d0_prime = inst.d0_true + 2 * t

        ok, K_prime = candidate_k(inst, d0_prime)

        if not ok:
            return False

        predicted_delta = t * inst.d0_true + t * t
        predicted_K = inst.K_true + predicted_delta

        if K_prime != predicted_K:
            return False

    return True


# ======================================================================================================================
# DISTINCTNESS CHECK
# ======================================================================================================================

def all_distinct_pairs(
    pairs: list[tuple[int, int]],
) -> bool:

    return len(pairs) == len(set(pairs))


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    instances = build_instances()

    MODULI = [
        4,
        8,
        16,
        32,
        64,
        128,
        256,
        3,
        5,
        7,
        11,
        13,
        17,
        19,
    ]

    WINDOW_RADII = [
        0,
        1,
        2,
        4,
        8,
        16,
        32,
        64,
        128,
    ]

    FAMILY_T = [
        -1000,
        -100,
        -32,
        -16,
        -8,
        -4,
        -2,
        -1,
        0,
        1,
        2,
        4,
        8,
        16,
        32,
        100,
        1000,
    ]

    print("=" * 120)
    print("EXPERIMENT 448")
    print("=" * 120)
    print()
    print("EXACT d0 AMBIGUITY / CONGRUENCE-IDENTIFIABILITY AUDIT")
    print()
    print("QUESTION")
    print(
        "  Given only N,S,d,r, does the equation"
    )
    print(
        "      4K = r + d0^2 - d^2"
    )
    print(
        "  determine d0?"
    )
    print()
    print("CORE PARAMETRIZATION")
    print(
        "  K(d0) = (r + d0^2 - d^2)/4"
    )
    print()
    print("PARITY FAMILY")
    print(
        "  d0' = d0 + 2t"
    )
    print(
        "  K(d0') = K(d0) + t*d0 + t^2"
    )
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no factorization during reconstruction")
    print("  no giant K enumeration")
    print("  no giant X enumeration")
    print("  no CRT Cartesian product")
    print()

    global_integrality_failures = 0
    global_family_failures = 0
    global_true_recovery_failures = 0
    global_pair_collision_count = 0

    total_window_candidates = 0
    total_valid_window_candidates = 0

    # ==============================================================================================================
    # INSTANCES
    # ==============================================================================================================

    for inst in instances:

        print("-" * 120)
        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} "
            f"S={inst.S} "
            f"d={inst.d} "
            f"r={inst.r}"
        )

        print()
        print("  OBSERVED-ONLY SETUP")
        print(
            "    reconstruction input:"
        )
        print(
            f"      N = {inst.N}"
        )
        print(
            f"      S = {inst.S}"
        )
        print(
            f"      d = {inst.d}"
        )
        print(
            f"      r = {inst.r}"
        )

        print()
        print("  POST-HOC TRUE CHECK")
        print(
            f"    true d0 = {inst.d0_true}"
        )
        print(
            f"    true x  = {inst.x_true}"
        )
        print(
            f"    true K  = {inst.K_true}"
        )

        # ----------------------------------------------------------------------------------------------------------
        # Exact true equation check.
        # ----------------------------------------------------------------------------------------------------------

        ok_true, recovered_K_true = candidate_k(
            inst,
            inst.d0_true,
        )

        if not ok_true:
            global_integrality_failures += 1

        if recovered_K_true != inst.K_true:
            global_true_recovery_failures += 1

        print()
        print("  EXACT TRUE EQUATION")
        print(
            "    r + d0^2 - d^2 divisible by 4 = "
            f"{ok_true}"
        )

        print(
            f"    recovered K(d0_true) = {recovered_K_true}"
        )

        print(
            "    recovered K == hidden K = "
            f"{recovered_K_true == inst.K_true}"
        )

        # ----------------------------------------------------------------------------------------------------------
        # Parity condition.
        # ----------------------------------------------------------------------------------------------------------

        required_parity = parity_required(inst)

        print()
        print("  PARITY IDENTIFIABILITY")
        print(
            f"    required d0 parity = {required_parity}"
        )

        print(
            "    d0_true parity      = "
            f"{inst.d0_true & 1}"
        )

        print(
            "    parity match        = "
            f"{(inst.d0_true & 1) == required_parity}"
        )

        # ----------------------------------------------------------------------------------------------------------
        # Window audit.
        # ----------------------------------------------------------------------------------------------------------

        print()
        print("  LOCAL d0 -> x,K FAMILY")
        print()
        print(
            "    radius"
            "       total d0"
            "     valid d0"
            "       x-low"
            "      x-high"
            "    distinct (x,K)"
        )
        print(
            "    " + "-" * 92
        )

        for radius in WINDOW_RADII:

            total = 2 * radius + 1
            candidates = scan_d0_window(inst, radius)

            total_window_candidates += total
            total_valid_window_candidates += len(candidates)

            if candidates:

                xs = [item[1] for item in candidates]
                ks = [item[2] for item in candidates]

                xmin = min(xs)
                xmax = max(xs)

                pairs = [
                    (item[1], item[2])
                    for item in candidates
                ]

                distinct = all_distinct_pairs(pairs)

                if not distinct:
                    global_pair_collision_count += 1

            else:

                xmin = 0
                xmax = 0
                distinct = True

            print(
                f"    {radius:>8}"
                f"{total:>14}"
                f"{len(candidates):>14}"
                f"{xmin:>14}"
                f"{xmax:>14}"
                f"{str(distinct):>18}"
            )

        # ----------------------------------------------------------------------------------------------------------
        # Affine-family proof.
        # ----------------------------------------------------------------------------------------------------------

        family_ok = affine_family_check(
            inst,
            FAMILY_T,
        )

        if not family_ok:
            global_family_failures += 1

        print()
        print("  AFFINE FAMILY CHECK")
        print(
            "    test t values = "
            f"{len(FAMILY_T)}"
        )

        print(
            "    d0' = d0 + 2t identity = "
            f"{family_ok}"
        )

        print()
        print(
            "    representative family:"
        )

        for t in (-4, -2, -1, 0, 1, 2, 4):

            d0_prime = inst.d0_true + 2 * t

            ok, K_prime = candidate_k(
                inst,
                d0_prime,
            )

            x_prime = candidate_x(
                inst,
                d0_prime,
            )

            delta = K_prime - inst.K_true

            print(
                f"      t={t:>3} "
                f"d0'={d0_prime} "
                f"x'={x_prime} "
                f"K'={K_prime} "
                f"K'-Ktrue={delta}"
            )

        # ----------------------------------------------------------------------------------------------------------
        # Small-modulus audit.
        # ----------------------------------------------------------------------------------------------------------

        print()
        print("  LOCAL CONGRUENCE AUDIT")
        print()
        print(
            "    modulus"
            "    valid d0 residues"
            "    number"
        )
        print(
            "    " + "-" * 64
        )

        for modulus in MODULI:

            residues = valid_residue_classes(
                inst,
                modulus,
            )

            # For odd primes the test is intentionally diagnostic:
            # only the integrality equation is being checked.
            print(
                f"    {modulus:>7}"
                f"    {str(residues):<46}"
                f"{len(residues):>8}"
            )

        print()
        print("  OBSERVATION")

        print(
            "    Candidate d0 values with the required parity"
        )

        print(
            "    generate exact integer K values."
        )

        print(
            "    Changing d0 does not cause the equation to fail;"
        )

        print(
            "    it merely changes which integer K is associated"
        )

        print(
            "    with that d0."
        )

    # ==============================================================================================================
    # GLOBAL
    # ==============================================================================================================

    print()
    print("=" * 120)
    print("GLOBAL d0 IDENTIFIABILITY SUMMARY")
    print("=" * 120)

    print(
        f"  instances                         = {len(instances)}"
    )

    print(
        f"  true-equation integrality failures = "
        f"{global_integrality_failures}"
    )

    print(
        f"  true-K reconstruction failures    = "
        f"{global_true_recovery_failures}"
    )

    print(
        f"  affine-family failures             = "
        f"{global_family_failures}"
    )

    print(
        f"  (x,K) family collisions            = "
        f"{global_pair_collision_count}"
    )

    print(
        f"  window d0 candidates tested        = "
        f"{total_window_candidates}"
    )

    print(
        f"  parity-valid d0 candidates         = "
        f"{total_valid_window_candidates}"
    )

    print()
    print("CENTRAL ALGEBRAIC RESULT")
    print()
    print(
        "  Starting from"
    )
    print(
        "      4K = r + d0^2 - d^2"
    )
    print()
    print(
        "  take"
    )
    print(
        "      d0' = d0 + 2t."
    )
    print()
    print(
        "  Then"
    )
    print(
        "      K' = K + t*d0 + t^2."
    )
    print()
    print(
        "  Therefore every parity-compatible d0' gives an exact"
    )
    print(
        "  integer K'."
    )

    print()
    print("WHAT THIS MEANS")
    print()
    print(
        "  The r/K quadratic relation by itself does not select"
    )
    print(
        "  a unique d0."
    )

    print(
        "  Its generic local information is the parity class"
    )
    print(
        "      d0 mod 2."
    )

    print()
    print(
        "  Equivalently, after substituting"
    )
    print(
        "      x = d-d0,"
    )
    print(
        "  the same equation becomes"
    )
    print(
        "      4K = r - 2dx + x^2."
    )

    print(
        "  Once x is chosen in the compatible parity class, K"
    )
    print(
        "  simply changes accordingly."
    )

    print()
    print("CONNECTION TO EXPERIMENT 447")
    print()
    print(
        "  Experiment 447 established:"
    )
    print(
        "      x = d-d0."
    )

    print(
        "  Experiment 448 now establishes the complementary fact:"
    )
    print(
        "      the observed quadratic equation does not, by itself,"
    )
    print(
        "      recover d0 uniquely."
    )

    print()
    print("NEXT RESEARCH TARGET")
    print()
    print(
        "  To recover d0, one needs information that constrains K"
    )
    print(
        "  independently of the identity"
    )
    print(
        "      K = (r+d0^2-d^2)/4."
    )

    print()
    print(
        "  In particular, another local quadratic-residue condition"
    )
    print(
        "  on this same K expression is expected to pull back to an"
    )
    print(
        "  identity rather than provide a new d0 constraint."
    )

    print()
    print(
        "  A genuinely new experiment should therefore search for"
    )
    print(
        "  independent K structure, rather than additional local"
    )
    print(
        "  residue tests of K(d0)."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 448 FINAL STATUS")
    print("=" * 120)

    print(
        "  EXACT TRUE EQUATION = "
        f"{global_integrality_failures == 0}"
    )

    print(
        "  EXACT TRUE K RECOVERY = "
        f"{global_true_recovery_failures == 0}"
    )

    print(
        "  AFFINE d0/K FAMILY IDENTITY = "
        f"{global_family_failures == 0}"
    )

    print(
        "  d0 UNIQUELY DETERMINED BY r,d = False"
    )

    print(
        "  LOCAL INFORMATION GENERICALLY = parity only"
    )

    print(
        "  GIANT K ENUMERATION = False"
    )

    print(
        "  GIANT X ENUMERATION = False"
    )

    print(
        "  CRT CARTESIAN PRODUCT = False"
    )

    print(
        "  INTEGER-EXACT = True"
    )

    print()
    print(
        "  CONCLUSION:"
    )
    print(
        "    The present quadratic relation parameterizes a family"
    )
    print(
        "    of (d0,K,x) solutions rather than identifying d0."
    )

    print("=" * 120)
    print("EXPERIMENT 448 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
