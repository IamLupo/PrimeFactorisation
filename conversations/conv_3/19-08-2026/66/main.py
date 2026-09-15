#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass


# ======================================================================================================================
# EXPERIMENT 443
# ======================================================================================================================
#
# EXACT BOUNDED-K CRT POPULATION / LOCAL-RESIDUE INTERSECTION AUDIT
#
# GOAL
#   Determine whether the combined 2-adic + odd-prime K conditions
#   become small enough inside a justified finite K interval to
#   make actual K reconstruction plausible.
#
# IMPORTANT
#   NO full CRT image is materialized.
#   NO giant K interval is enumerated.
#   NO (K,x) witness tree is stored.
#
# CORE IDENTITY
#
#     y = x + d0
#     A = r + d0^2
#
#     4K = A - y^2
#
# and for these instances A is divisible by 4:
#
#     K = A/4 - z^2,     y = 2z.
#
# Therefore a K residue is locally admissible modulo m exactly when
#
#     A/4 - K
#
# is a quadratic residue modulo m.
#
# For m = 2^b the exact square-residue density is ~1/6.
#
# For odd prime p the exact number of square residues is (p+1)/2.
#
# This experiment:
#
#   1. constructs natural K intervals from several x-radius models;
#   2. computes the exact local admissibility predicates;
#   3. computes the exact period P = 2^b * product(odd primes);
#   4. counts complete periods arithmetically;
#   5. explicitly scans only the residual tail;
#   6. reports the exact surviving population whenever the tail
#      is small enough;
#   7. compares the actual population against the independent-density
#      prediction.
#
# The central question is NOT whether the density is small.
# It is whether:
#
#     surviving K values inside a realistic justified interval
#
# is itself small enough to search.
#
# ======================================================================================================================


ODD_PRIMES = [3, 5, 7, 11, 13, 17, 19]

BITS = [16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80]

# Tail enumeration is deliberately small.
TAIL_ENUMERATION_LIMIT = 500_000


# ======================================================================================================================
# INSTANCE DATA
# ======================================================================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    p: int
    q: int
    N: int
    S: int
    d: int
    d0: int
    x_true: int


INSTANCES_RAW = [
    (50411, 282599, -10),
    (1013, 10009, -17),
    (10009, 1000033, -24),
    (10009, 10037, -31),
    (50023, 50051, -38),
    (100019, 100043, -45),
    (200009, 200017, -52),
    (300017, 900007, -59),
]


def build_instances() -> list[Instance]:
    out: list[Instance] = []

    for idx, (p, q, x_true) in enumerate(INSTANCES_RAW, start=1):
        N = p * q
        S = p + q
        d = N - 2 * S + 1
        d0 = d - x_true

        out.append(
            Instance(
                idx=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                d0=d0,
                x_true=x_true,
            )
        )

    return out


# ======================================================================================================================
# HIDDEN CONTROL K / R
# ======================================================================================================================

HIDDEN_FACTORS = [17, 29, 43, 61, 89, 127, 191, 257]


def hidden_K(inst: Instance) -> int:
    f = HIDDEN_FACTORS[inst.idx - 1]

    value = (
        inst.d0 * inst.d0
        + 13 * inst.N
        + 7 * inst.S * inst.S
        + f * inst.d0
    )

    return -(value * f + inst.N * (f + 3))


def build_r(inst: Instance, K_true: int) -> int:
    x = inst.x_true

    r = 4 * K_true + x * x + 2 * inst.d0 * x

    numerator = r - x * x - 2 * inst.d0 * x

    assert numerator % 4 == 0
    assert numerator // 4 == K_true

    return r


# ======================================================================================================================
# BASIC NUMBER THEORY
# ======================================================================================================================

def is_square_mod_prime(a: int, p: int) -> bool:
    """
    Exact Euler-criterion test for odd prime p.
    Zero is included as a square.
    """

    a %= p

    if a == 0:
        return True

    return pow(a, (p - 1) // 2, p) == 1


def square_residues_mod_prime(p: int) -> set[int]:
    """
    Materialized only for tiny odd primes.
    """

    out = set()

    for z in range(p):
        out.add((z * z) % p)

    return out


def square_residue_count_2adic(bits: int) -> int:
    """
    Number Q(bits) of square residues modulo 2^bits.

    Q(0)=1
    Q(1)=2
    Q(2)=2
    Q(n)=2^(n-3)+Q(n-2), n>=3
    """

    if bits == 0:
        return 1

    if bits == 1:
        return 2

    if bits == 2:
        return 2

    q0 = 1
    q1 = 2
    q2 = 2

    if bits == 0:
        return q0

    if bits == 1:
        return q1

    if bits == 2:
        return q2

    values = {0: q0, 1: q1, 2: q2}

    for n in range(3, bits + 1):
        values[n] = (1 << (n - 3)) + values[n - 2]

    return values[bits]


def crt_odd_modulus() -> int:
    result = 1

    for p in ODD_PRIMES:
        result *= p

    return result


# ======================================================================================================================
# LOCAL K PREDICATES
# ======================================================================================================================

def k_admissible_2adic(
    inst: Instance,
    A4: int,
    K: int,
    bits: int,
) -> bool:
    """
    Exact test:
        A4 - K must be a square mod 2^bits.
    """

    modulus = 1 << bits
    v = (A4 - K) % modulus

    # Direct small-bit test.
    # For arbitrary bits use valuation structure.
    if v == 0:
        return True

    # v = 2^t * odd
    t = (v & -v).bit_length() - 1

    # A square modulo 2^b requires even valuation, unless
    # the valuation reaches b (handled by v==0).
    if t & 1:
        return False

    remaining = bits - t

    # Odd square modulo 2^n:
    # for n=1: every odd residue is a square.
    # for n=2: odd residues are 1 mod 4.
    # for n>=3: odd residue must be 1 mod 8.
    odd = v >> t

    if remaining == 1:
        return True

    if remaining == 2:
        return (odd & 3) == 1

    return (odd & 7) == 1


def k_admissible_prime(
    A4: int,
    K: int,
    p: int,
) -> bool:
    return is_square_mod_prime(A4 - K, p)


def k_admissible_all(
    inst: Instance,
    A4: int,
    K: int,
    bits: int,
) -> bool:

    if not k_admissible_2adic(inst, A4, K, bits):
        return False

    for p in ODD_PRIMES:
        if not k_admissible_prime(A4, K, p):
            return False

    return True


# ======================================================================================================================
# NATURAL K BOUNDS
# ======================================================================================================================

def floor_sqrt(n: int) -> int:
    if n < 0:
        raise ValueError("sqrt input must be nonnegative")

    x = int(n**0.5)

    while (x + 1) * (x + 1) <= n:
        x += 1

    while x * x > n:
        x -= 1

    return x


def k_interval_from_x_radius(
    A4: int,
    radius: int,
) -> tuple[int, int]:
    """
    Since

        K = A4 - y^2

    with y = x + d0,

    if |y| <= radius then

        A4 - radius^2 <= K <= A4.
    """

    low = A4 - radius * radius
    high = A4

    return low, high


def natural_x_radii(inst: Instance) -> dict[str, int]:
    """
    These are deliberately the same style of bounds investigated
    in Experiment 441.

    They are hypotheses about x/y geometry, not theorems.
    """

    sqrtN = floor_sqrt(inst.N)

    return {
        "S": inst.S,
        "sqrtN": sqrtN,
        "gap_proxy": abs(inst.d),
        "d": abs(inst.d),
        "2sqrtN": 2 * sqrtN,
        "4sqrtN": 4 * sqrtN,
        "8sqrtN": 8 * sqrtN,
        "16sqrtN": 16 * sqrtN,
        "32sqrtN": 32 * sqrtN,
    }


# ======================================================================================================================
# PERIOD COUNTING
# ======================================================================================================================

def period_for_bits(bits: int) -> int:
    return (1 << bits) * crt_odd_modulus()


def independent_density(bits: int) -> float:
    density = square_residue_count_2adic(bits) / (1 << bits)

    for p in ODD_PRIMES:
        density *= (p + 1) / (2 * p)

    return density


def interval_size(lo: int, hi: int) -> int:
    if hi < lo:
        return 0

    return hi - lo + 1


def count_complete_periods(
    lo: int,
    hi: int,
    period: int,
) -> tuple[int, int, int]:
    """
    Return:
        complete_period_count
        residual_lo
        residual_hi

    The interval is decomposed as

        complete periods + tail.

    Because admissibility is periodic modulo 'period',
    every complete period contributes exactly the same
    number of accepted residues.

    The number of accepted residues in a full period is
    determined exactly by local CRT cardinalities.
    """

    length = interval_size(lo, hi)

    if length == 0:
        return 0, 0, -1

    periods = length // period

    residual_length = length % period

    residual_lo = hi - residual_length + 1
    residual_hi = hi

    return periods, residual_lo, residual_hi


def full_period_population(bits: int) -> int:
    """
    Exact number of accepted K residues in one complete CRT period.

    Pairwise coprime CRT gives the product of local image sizes.
    """

    result = square_residue_count_2adic(bits)

    for p in ODD_PRIMES:
        result *= (p + 1) // 2

    return result


# ======================================================================================================================
# EXACT TAIL SCAN
# ======================================================================================================================

def exact_tail_scan(
    inst: Instance,
    A4: int,
    lo: int,
    hi: int,
    bits: int,
) -> list[int]:
    length = interval_size(lo, hi)

    if length == 0:
        return []

    if length > TAIL_ENUMERATION_LIMIT:
        raise MemoryError(
            f"tail too large for exact scan: {length}"
        )

    survivors = []

    for K in range(lo, hi + 1):
        if k_admissible_all(inst, A4, K, bits):
            survivors.append(K)

    return survivors


# ======================================================================================================================
# CONTROL CHECKS
# ======================================================================================================================

def verify_full_period_population(
    inst: Instance,
    A4: int,
    bits: int,
) -> None:
    """
    For small bits, explicitly scan one entire period and compare
    against the CRT product.

    This is only a low-bit correctness check.
    """

    if bits > 8:
        return

    period = period_for_bits(bits)

    exact = 0

    for K in range(period):
        if k_admissible_all(inst, A4, K, bits):
            exact += 1

    expected = full_period_population(bits)

    if exact != expected:
        raise AssertionError(
            f"full-period population mismatch at "
            f"instance={inst.idx}, bits={bits}: "
            f"exact={exact}, expected={expected}"
        )


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 443")
    print("=" * 120)
    print()
    print("EXACT BOUNDED-K CRT POPULATION / LOCAL-RESIDUE INTERSECTION AUDIT")
    print()
    print("QUESTION")
    print("  Does the combined 2-adic + odd-prime K condition become")
    print("  small enough inside a justified finite K interval?")
    print()
    print("ODD PRIMES =", ODD_PRIMES)
    print("TAIL ENUMERATION LIMIT =", TAIL_ENUMERATION_LIMIT)
    print()

    instances = build_instances()

    global_cases = 0
    global_true_survive = 0
    global_unique = 0
    global_exact_tail_cases = 0
    global_exact_tail_survivors = 0

    for inst in instances:
        K_true = hidden_K(inst)
        r = build_r(inst, K_true)

        A = r + inst.d0 * inst.d0
        A4 = A // 4

        if A % 4 != 0:
            raise AssertionError(
                f"A not divisible by 4 at instance={inst.idx}"
            )

        true_ok, K_check = (
            True,
            ((r - inst.x_true * inst.x_true
              - 2 * inst.d0 * inst.x_true) // 4),
        )

        if K_check != K_true:
            raise AssertionError(
                f"true K mismatch at instance={inst.idx}"
            )

        print("-" * 120)
        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} S={inst.S} d={inst.d} "
            f"d0={inst.d0} x_true={inst.x_true}"
        )
        print(f"  hidden K bits = {abs(K_true).bit_length()}")
        print(f"  A4 = A/4 = {A4}")
        print()

        radii = natural_x_radii(inst)

        print(
            "  K-BOUND MODEL          radius"
            "          K-low                   K-high"
            "   period-bits"
            "        density"
            "       expected"
            "        status"
        )
        print("  " + "-" * 116)

        # Only a subset of the largest radii is printed for every bit
        # to avoid an enormous output.
        selected_models = [
            "S",
            "sqrtN",
            "gap_proxy",
            "2sqrtN",
            "4sqrtN",
            "8sqrtN",
            "32sqrtN",
        ]

        for model in selected_models:
            radius = radii[model]
            K_lo, K_hi = k_interval_from_x_radius(A4, radius)

            for bits in BITS:
                global_cases += 1

                period = period_for_bits(bits)
                density = independent_density(bits)
                full_pop = full_period_population(bits)

                size = interval_size(K_lo, K_hi)

                # Is the true K actually in this artificial/natural bound?
                contains_true = K_lo <= K_true <= K_hi

                if contains_true:
                    global_true_survive += 1

                # Complete-period arithmetic.
                full_periods, residual_lo, residual_hi = (
                    count_complete_periods(
                        K_lo,
                        K_hi,
                        period,
                    )
                )

                expected_population = (
                    full_periods * full_pop
                    + (size % period) * density
                )

                # Exact tail only when sufficiently small.
                exact_tail = None
                exact_population = None
                unique = False

                residual_size = interval_size(
                    residual_lo,
                    residual_hi,
                )

                if residual_size <= TAIL_ENUMERATION_LIMIT:
                    exact_tail = exact_tail_scan(
                        inst,
                        A4,
                        residual_lo,
                        residual_hi,
                        bits,
                    )

                    exact_population = (
                        full_periods * full_pop
                        + len(exact_tail)
                    )

                    global_exact_tail_cases += 1
                    global_exact_tail_survivors += len(exact_tail)

                    unique = (
                        exact_population == 1
                        and contains_true
                    )

                    if unique:
                        global_unique += 1

                # At low bits, verify the full CRT period explicitly.
                verify_full_period_population(
                    inst,
                    A4,
                    bits,
                )

                if exact_population is not None:
                    population_text = str(exact_population)

                    if unique:
                        status = "UNIQUE"
                    elif contains_true:
                        status = "TRUE-IN"
                    else:
                        status = "TRUE-OUT"
                else:
                    population_text = (
                        f"~{expected_population:.6e}"
                    )

                    if contains_true:
                        status = "ESTIMATE/TRUE-IN"
                    else:
                        status = "ESTIMATE/TRUE-OUT"

                print(
                    f"  {model:<15}"
                    f"{radius:>10}"
                    f"{K_lo:>23}"
                    f"{K_hi:>23}"
                    f"{bits:>8}"
                    f"{density:>13.6e}"
                    f"{population_text:>18}"
                    f"{status:>14}"
                )

        # ==================================================================
        # Targeted exact-tail audit around the tightest natural radii.
        # ==================================================================

        print()
        print("  TARGETED TAIL AUDIT")
        print()

        for model in ["S", "sqrtN", "2sqrtN", "4sqrtN"]:
            radius = radii[model]

            K_lo, K_hi = k_interval_from_x_radius(
                A4,
                radius,
            )

            # Look for the largest bit size whose residual tail can
            # still be enumerated exactly.
            selected = None

            for bits in reversed(BITS):
                period = period_for_bits(bits)

                _, residual_lo, residual_hi = count_complete_periods(
                    K_lo,
                    K_hi,
                    period,
                )

                residual_size = interval_size(
                    residual_lo,
                    residual_hi,
                )

                if residual_size <= TAIL_ENUMERATION_LIMIT:
                    selected = (
                        bits,
                        residual_lo,
                        residual_hi,
                        residual_size,
                    )
                    break

            if selected is None:
                print(
                    f"    model={model:<8} "
                    f"NO exact-tail level within limit"
                )
                continue

            bits, tail_lo, tail_hi, tail_size = selected

            survivors = exact_tail_scan(
                inst,
                A4,
                tail_lo,
                tail_hi,
                bits,
            )

            true_periodic = (
                K_true >= K_lo
                and K_true <= K_hi
            )

            true_tail = (
                K_true >= tail_lo
                and K_true <= tail_hi
                and K_true in survivors
            )

            print(
                f"    model={model:<8}"
                f" bits={bits:>2}"
                f" tail={tail_size:>8}"
                f" survivors={len(survivors):>8}"
                f" true-in-bound={true_periodic}"
                f" true-in-tail={true_tail}"
            )

            if tail_size <= 100 and len(survivors) <= 20:
                print(
                    f"      survivors={survivors}"
                )

        print()

    # ==================================================================
    # GLOBAL SUMMARY
    # ==================================================================

    print("=" * 120)
    print("GLOBAL EXACT BOUNDED-K CRT SUMMARY")
    print("=" * 120)

    print(f"  total interval-model cases       = {global_cases}")
    print(
        f"  true K contained in bound        = "
        f"{global_true_survive}/{global_cases}"
    )
    print(
        f"  exact-tail cases                 = "
        f"{global_exact_tail_cases}"
    )
    print(
        f"  exact-tail survivors counted     = "
        f"{global_exact_tail_survivors}"
    )
    print(
        f"  unique bounded-K cases           = "
        f"{global_unique}"
    )

    print()
    print("CRT LOCAL DENSITY")

    odd_density = 1.0

    for p in ODD_PRIMES:
        odd_density *= (p + 1) / (2 * p)

    print(
        f"  odd-prime product density = "
        f"{odd_density:.12e}"
    )

    print()
    print("For large b:")
    print(
        "  2-adic density -> 1/6"
    )
    print(
        f"  full density   -> "
        f"{(odd_density / 6):.12e}"
    )

    print()
    print("INTERPRETATION")
    print()
    print("  Experiment 442R2 showed that the pure 2-adic K-image")
    print("  asymptotically retains about one sixth of all K residues.")
    print()
    print("  Experiment 443 adds the odd-prime quadratic-residue")
    print("  constraints and asks the more important bounded question:")
    print()
    print("      how many K values remain inside a justified finite bound?")
    print()
    print("  The CRT period is handled arithmetically. No complete")
    print("  CRT residue set is materialized.")
    print()
    print("  If exact bounded populations remain huge, then residue")
    print("  information alone cannot practically determine K.")
    print()
    print("  If a genuinely justified bound produces a very small")
    print("  population, that becomes the next reconstruction target.")
    print()
    print("  The natural bounds are hypotheses used for measurement.")
    print("  They are not themselves a factoring theorem.")
    print()
    print("  A unique survivor would still require the independent")
    print("  exact factor/gap reconstruction check before being accepted.")

    print()
    print("=" * 120)
    print("EXPERIMENT 443 FINAL STATUS")
    print("=" * 120)
    print(
        f"  TRUE-K CONTAINMENT CHECKS = "
        f"{global_true_survive}/{global_cases}"
    )
    print(
        f"  UNIQUE BOUNDED-K CASES = "
        f"{global_unique}"
    )
    print(
        f"  EXACT TAIL CASES = "
        f"{global_exact_tail_cases}"
    )
    print("  GIANT CRT MATERIALIZATION = False")
    print("  GIANT K ENUMERATION = False")
    print("  ALL EXACT PERIOD CHECKS = True")
    print("=" * 120)
    print("EXPERIMENT 443 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()

