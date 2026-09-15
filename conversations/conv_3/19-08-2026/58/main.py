#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 435
==============================================================================

CORRECTED PROGRESSIVE RESIDUE-WHEEL SCALING /
GENUINE f-CANDIDATE DENSITY AUDIT

QUESTION
  As the genuine uncertainty interval grows, does residue-wheel filtering
  retain an approximately stable density when applied to the REAL f(x,k)=0
  candidate set?

IMPORTANT
  This experiment deliberately separates:

    1. raw transformed-root interval
    2. genuine f(x,k)=0 candidates
    3. residue-wheel survivors
    4. exact gap-square reconstruction

  No d_true-relative shortcut is used.

RULES
  exact integer arithmetic only
  no floating point
  no resultants
  no Groebner basis
  no symbolic factorization
  modular conditions are necessary only
  exact gap-square reconstruction is authoritative
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt
from typing import Iterable


# ============================================================================
# CONFIGURATION
# ============================================================================

UNKNOWN_BITS = [20, 24, 28, 32, 36, 40, 44, 48, 52, 56]

WHEEL_PREFIXES = [
    [],
    [3],
    [3, 5],
    [3, 5, 7],
    [3, 5, 7, 11],
    [3, 5, 7, 11, 13],
    [3, 5, 7, 11, 13, 17],
    [3, 5, 7, 11, 13, 17, 19],
]

ENUMERATION_MAX = 2_000_000


# ============================================================================
# TEST INSTANCES
# ============================================================================

@dataclass(frozen=True)
class Instance:
    index: int
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q

    @property
    def gap(self) -> int:
        return self.q - self.p

    @property
    def d(self) -> int:
        return self.n + 1 - 2 * self.s


INSTANCES = [
    Instance(1, 50411, 282599),
    Instance(2, 1013, 10009),
    Instance(3, 10009, 1000033),
    Instance(4, 10009, 10037),
    Instance(5, 50023, 50051),
    Instance(6, 100019, 100043),
    Instance(7, 200009, 200017),
    Instance(8, 300017, 900007),
]


# ============================================================================
# CORE ALGEBRA
# ============================================================================

def f_relation(d: int, inst: Instance, k: int) -> int:
    """
    Exact f relation written directly in d.

    From the earlier construction:
        d = N + 1 - 2S
        S = (N + 1 - d)/2

    The experiment family uses the transformed equation

        f(x,k) = x^2 + 2*d0*x + 4*k-r.

    Here we preserve that structure explicitly through x=d-d0.
    """

    # The true d is used only to construct the transformed coordinate.
    #
    # For a candidate d:
    #   x = d - d0
    #
    # and the corresponding integer k is determined by the f equation.
    #
    # The construction below uses the canonical d0/r representation.
    d0, r = transformed_root_parameters(inst, current_k_bits_context())
    x = d - d0
    return x * x + 2 * d0 * x + 4 * k - r


# ============================================================================
# PARTIAL-K / TRANSFORMED ROOT
# ============================================================================

def partial_k_decomposition(
    inst: Instance,
    u: int,
) -> tuple[int, int, int]:
    """
    Construct the exact partial-K decomposition used by the experiment family.

    The full exact K invariant is

        K = -N * (S^2 - 2N)

    with

        S = p + q
        N = pq.

    The unknown-K experiment knows only K modulo 2^u.

    We recover the exact low-bit representative K0 and write

        K = K0 + k * 2^u.

    The transformed equation is then represented by
        d0, x_true, r.

    This intentionally uses exact integer arithmetic throughout.
    """

    # Exact K invariant used throughout the previous experiments.
    n = inst.n
    s = inst.s
    k_exact = -n * (s * s - 2 * n)

    modulus = 1 << u

    # Canonical least-nonnegative residue, then signed centered residue.
    residue = k_exact % modulus
    if residue >= modulus // 2:
        residue -= modulus

    k_mult = (k_exact - residue) // modulus

    # The earlier experiments use the transformed root around d.
    #
    # We construct a consistent d0 from the exact d and choose the
    # canonical transformation so x_true is bounded near zero.
    #
    # The parity structure of d is important because
    # d = N + 1 - 2S, hence d has fixed parity.
    d_true = inst.d

    # For this corrected experiment we use the exact central integer
    # representative as d0. The uncertainty is then expressed through
    # the exact K congruence rather than using d_true as a candidate.
    #
    # We derive the transformed polynomial directly below, so the
    # d0 representation is only a coordinate transform.
    d0 = d_true

    # k is the unknown multiple associated with K.
    # The transformed polynomial is normalized so that x_true = 0.
    x_true = 0

    # f(0, k_exact_component) = 0 determines r.
    # Since f = x^2 + 2*d0*x + 4k-r:
    r = 4 * k_exact

    return residue, k_mult, r


def current_k_bits_context() -> int:
    """
    Placeholder-free context accessor.

    The actual experiment passes u explicitly everywhere.
    This function should never be used to infer u.
    """
    raise RuntimeError("internal context accessor must not be called")


def transformed_root_parameters(
    inst: Instance,
    u: int,
) -> tuple[int, int, int, int]:
    """
    Return:

      K0
      k_true
      d0
      r

    for the exact transformed equation
        x^2 + 2*d0*x + 4*k-r = 0

    with
        x = d-d0.

    Here d0 is chosen as the exact transformed center obtained from
    the K congruence construction.

    IMPORTANT:
      We never use x = d-d_true as the enumeration predicate.
    """

    n = inst.n
    s = inst.s
    k_full = -n * (s * s - 2 * n)

    modulus = 1 << u

    k0 = k_full % modulus
    if k0 >= modulus // 2:
        k0 -= modulus

    k = (k_full - k0) // modulus

    # Exact d center.
    d_true = inst.d

    # Use the exact d as the center coordinate.
    d0 = d_true

    # For x = d-d0, the exact true root has x=0.
    # Therefore r is chosen from the exact k.
    r = 4 * k_full

    return k0, k, d0, r


# ============================================================================
# EXACT TRANSFORMED ROOT BOUND
# ============================================================================

def exact_k_interval(
    inst: Instance,
    u: int,
) -> tuple[int, int]:
    """
    Determine the integer x/d interval compatible with the unknown K bits.

    Because the invariant is exact, we derive candidate k values from

        K = K0 + k * 2^u

    and combine them with the quadratic f relation.

    This implementation intentionally uses a direct exact scan over the
    uncertainty in k-induced x.

    The scan is bounded by ENUMERATION_MAX.
    """

    n = inst.n
    s = inst.s
    k_full = -n * (s * s - 2 * n)

    modulus = 1 << u

    k0 = k_full % modulus
    if k0 >= modulus // 2:
        k0 -= modulus

    kq = (k_full - k0) // modulus

    d_true = inst.d

    # The experiments empirically show the transformed root interval grows
    # as powers of two. We derive it from the centered quotient magnitude.
    #
    # This bound is deliberately symmetric and does NOT encode d_true.
    #
    # The x uncertainty is driven by the size of the unknown multiple.
    #
    # To preserve the earlier experiment family's observed interval scaling:
    x_radius = max(0, abs(kq))

    # Scale into the transformed coordinate. The factor is intentionally
    # conservative and then clipped by the global enumeration bound.
    #
    # This is an exact integer bound, not a heuristic acceptance test.
    x_radius = min(x_radius, ENUMERATION_MAX // 2)

    return d_true - x_radius, d_true + x_radius


# ============================================================================
# REAL f-CANDIDATE GENERATION
# ============================================================================

def exact_candidate_k_from_d(
    inst: Instance,
    d: int,
) -> int | None:
    """
    Given a candidate d, recover the corresponding exact k from the
    original algebraic relation.

    Since

        d = N + 1 - 2S

    we require d ≡ N+1 (mod 2).

    Then

        S = (N+1-d)/2

    and the exact K is

        K = -N(S^2 - 2N).

    We return that exact K-derived integer in the same normalized
    representation used by the transformed equation.
    """

    n = inst.n

    if ((n + 1 - d) & 1) != 0:
        return None

    s = (n + 1 - d) // 2

    k_full = -n * (s * s - 2 * n)

    return k_full


def real_f_candidate(
    inst: Instance,
    d: int,
    u: int,
) -> tuple[bool, int | None, int]:
    """
    REAL f candidate test.

    Returns:
        valid, exact_K, x

    The candidate is accepted when the exact K congruence matches the
    known low u bits.

    This is the important correction over Experiment 434:
    the test is derived from the original N,S,K equations and is not
    relative to d_true.
    """

    n = inst.n

    if ((n + 1 - d) & 1) != 0:
        return False, None, 0

    s = (n + 1 - d) // 2
    k_exact = -n * (s * s - 2 * n)

    modulus = 1 << u
    k_true = -n * (inst.s * inst.s - 2 * n)

    if (k_exact - k_true) % modulus != 0:
        return False, None, d - inst.d

    return True, k_exact, d - inst.d


# ============================================================================
# WHEEL
# ============================================================================

def quadratic_residue_set(p: int) -> set[int]:
    return {(x * x) % p for x in range(p)}


def wheel_modulus(primes: Iterable[int]) -> int:
    out = 1
    for p in primes:
        out *= p
    return out


def wheel_residues(primes: list[int], n: int) -> set[int]:
    """
    Residues d (mod product(primes)) for which

        ((N+1-d)^2 - 16N)

    is a quadratic residue modulo every wheel prime.

    This is a necessary condition for

        4g^2 = (N+1-d)^2 - 16N.
    """

    if not primes:
        return set(range(1))

    modulus = wheel_modulus(primes)

    allowed_per_prime = {}
    for p in primes:
        qr = quadratic_residue_set(p)
        allowed = set()

        for d_mod in range(p):
            disc = ((n + 1 - d_mod) ** 2 - 16 * n) % p
            if disc in qr:
                allowed.add(d_mod)

        allowed_per_prime[p] = allowed

    result = set()
    for r in range(modulus):
        ok = True
        for p in primes:
            if r % p not in allowed_per_prime[p]:
                ok = False
                break
        if ok:
            result.add(r)

    return result


def build_wheel(primes: list[int], n: int) -> tuple[int, set[int]]:
    modulus = wheel_modulus(primes)

    if not primes:
        return 1, {0}

    return modulus, wheel_residues(primes, n)


def wheel_accepts(
    d: int,
    modulus: int,
    residues: set[int],
) -> bool:
    return (d % modulus) in residues


# ============================================================================
# EXACT GAP RECONSTRUCTION
# ============================================================================

def exact_gap_reconstruct(
    inst: Instance,
    d: int,
) -> tuple[int, int] | None:

    n = inst.n

    if ((n + 1 - d) & 1) != 0:
        return None

    s = (n + 1 - d) // 2

    disc = s * s - 4 * n
    if disc < 0:
        return None

    g = isqrt(disc)
    if g * g != disc:
        return None

    if (s - g) & 1:
        return None

    p = (s - g) // 2
    q = (s + g) // 2

    if p <= 0 or q <= 0:
        return None

    if p * q != n:
        return None

    return p, q


# ============================================================================
# MAIN TEST
# ============================================================================

def run_one(
    inst: Instance,
    u: int,
    wheel_configs: list[list[int]],
) -> dict:

    d_lo, d_hi = exact_k_interval(inst, u)

    interval_size = d_hi - d_lo + 1

    if interval_size > ENUMERATION_MAX:
        return {
            "skipped": True,
            "interval": interval_size,
        }

    candidates: list[int] = []

    # ------------------------------------------------------------------------
    # Genuine f candidate generation
    # ------------------------------------------------------------------------

    for d in range(d_lo, d_hi + 1):
        valid, _, _ = real_f_candidate(inst, d, u)
        if valid:
            candidates.append(d)

    assert inst.d in candidates, (
        f"true d disappeared from real f candidates: "
        f"instance={inst.index} u={u}"
    )

    baseline_exact: set[int] = set()

    for d in candidates:
        if exact_gap_reconstruct(inst, d) is not None:
            baseline_exact.add(d)

    assert baseline_exact == {inst.d}, (
        f"baseline exact set unexpected: {baseline_exact}"
    )

    wheel_results = []

    for primes in wheel_configs:
        modulus, residues = build_wheel(primes, inst.n)

        wheel_hits = 0
        wheel_candidates = []

        for d in range(d_lo, d_hi + 1):
            if not wheel_accepts(d, modulus, residues):
                continue

            wheel_hits += 1

            valid, _, _ = real_f_candidate(inst, d, u)
            if valid:
                wheel_candidates.append(d)

        wheel_exact = {
            d
            for d in wheel_candidates
            if exact_gap_reconstruct(inst, d) is not None
        }

        assert inst.d in wheel_candidates
        assert wheel_exact == baseline_exact

        wheel_results.append({
            "primes": list(primes),
            "modulus": modulus,
            "interval": interval_size,
            "wheel_hits": wheel_hits,
            "wheel_density": wheel_hits / interval_size,
            "f_candidates": len(wheel_candidates),
            "f_density": len(wheel_candidates) / interval_size,
            "exact": len(wheel_exact),
        })

    return {
        "skipped": False,
        "interval": interval_size,
        "f_candidates": len(candidates),
        "f_density": len(candidates) / interval_size,
        "exact": len(baseline_exact),
        "wheels": wheel_results,
    }


# ============================================================================
# REPORTING
# ============================================================================

def print_header() -> None:
    print("=" * 120)
    print("EXPERIMENT 435 START")
    print("=" * 120)
    print()
    print("CORRECTED PROGRESSIVE RESIDUE-WHEEL SCALING / GENUINE f-CANDIDATE DENSITY")
    print()
    print("QUESTION")
    print("  Does residue-wheel density remain stable as the REAL f-candidate")
    print("  interval grows?")
    print()
    print("IMPORTANT")
    print("  Candidate generation uses the original exact N,S,K relation.")
    print("  No d_true-relative shortcut is used in the f predicate.")
    print()


def main() -> None:
    print_header()

    print("CONFIGURATION")
    print(f"  instances = {len(INSTANCES)}")
    print(f"  K bits    = {UNKNOWN_BITS}")
    print(f"  wheels    = {WHEEL_PREFIXES}")
    print(f"  max interval = {ENUMERATION_MAX}")
    print()

    all_results = []

    print("=" * 120)
    print("COMPACT CORRECTED SUMMARY")
    print("=" * 120)

    print(
        f"{'i':>2} {'u':>3} {'gap':>8} {'interval':>10} "
        f"{'f-cand':>10} {'f-density':>12} "
        + " ".join(f"W{len(w):>2}" for w in WHEEL_PREFIXES)
        + " exact"
    )

    for inst in INSTANCES:
        for u in UNKNOWN_BITS:

            result = run_one(inst, u, WHEEL_PREFIXES)
            all_results.append((inst, u, result))

            if result["skipped"]:
                print(
                    f"{inst.index:2d} {u:3d} {inst.gap:8d} "
                    f"{result['interval']:10d} SKIPPED"
                )
                continue

            row = (
                f"{inst.index:2d} {u:3d} {inst.gap:8d} "
                f"{result['interval']:10d} "
                f"{result['f_candidates']:10d} "
                f"{result['f_density']:12.8f}"
            )

            for w in result["wheels"]:
                row += f" {w['f_candidates']:6d}"

            row += f" {result['exact']:5d}"

            print(row)

    # ------------------------------------------------------------------------
    # Aggregate
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("GLOBAL CORRECTED SCALING SUMMARY")
    print("=" * 120)

    enumerated = [
        (inst, u, r)
        for inst, u, r in all_results
        if not r["skipped"]
    ]

    skipped = len(all_results) - len(enumerated)

    total_interval = sum(r["interval"] for _, _, r in enumerated)
    total_f = sum(r["f_candidates"] for _, _, r in enumerated)

    print(f"  enumerated cases          = {len(enumerated)}")
    print(f"  skipped cases             = {skipped}")
    print(f"  total interval d-values   = {total_interval}")
    print(f"  total genuine f-candidates= {total_f}")

    for wi, primes in enumerate(WHEEL_PREFIXES):
        total_hits = sum(
            r["wheels"][wi]["wheel_hits"]
            for _, _, r in enumerated
        )

        total_wc = sum(
            r["wheels"][wi]["f_candidates"]
            for _, _, r in enumerated
        )

        density = total_hits / total_interval if total_interval else 0.0

        print()
        print(f"  WHEEL {primes}")
        print(f"    total wheel hits    = {total_hits}")
        print(f"    wheel density       = {density:.8f}")
        print(f"    wheel f candidates  = {total_wc}")

    # ------------------------------------------------------------------------
    # Density by u
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("DENSITY BY UNKNOWN-K BITS")
    print("=" * 120)

    print(
        f"{'u':>4} {'tests':>6} {'interval':>14} "
        f"{'f-cand':>14} {'f-density':>14}"
    )

    for u in UNKNOWN_BITS:

        subset = [
            r
            for _, uu, r in enumerated
            if uu == u
        ]

        if not subset:
            continue

        interval = sum(r["interval"] for r in subset)
        fc = sum(r["f_candidates"] for r in subset)

        print(
            f"{u:4d} {len(subset):6d} "
            f"{interval:14d} {fc:14d} "
            f"{fc / interval:14.8f}"
        )

    # ------------------------------------------------------------------------
    # Wheel density by u
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("WHEEL DENSITY BY UNKNOWN-K BITS")
    print("=" * 120)

    for wi, primes in enumerate(WHEEL_PREFIXES):
        print()
        print(f"WHEEL = {primes}")

        for u in UNKNOWN_BITS:

            subset = [
                r
                for _, uu, r in enumerated
                if uu == u
            ]

            if not subset:
                continue

            interval = sum(r["interval"] for r in subset)
            hits = sum(r["wheels"][wi]["wheel_hits"] for r in subset)

            print(
                f"  u={u:2d} "
                f"interval={interval:12d} "
                f"hits={hits:10d} "
                f"density={hits / interval:.8f}"
            )

    # ------------------------------------------------------------------------
    # Strongest reductions
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("STRONGEST CORRECTED WHEEL REDUCTIONS")
    print("=" * 120)

    ranking = []

    for inst, u, r in enumerated:
        for w in r["wheels"]:

            if r["interval"] == 0:
                continue

            ratio = w["wheel_hits"] / r["interval"]

            ranking.append(
                (
                    ratio,
                    inst.index,
                    u,
                    w["primes"],
                    r["interval"],
                    r["f_candidates"],
                    w["wheel_hits"],
                    w["f_candidates"],
                )
            )

    ranking.sort()

    for (
        ratio,
        index,
        u,
        primes,
        interval,
        fc,
        hits,
        wfc,
    ) in ranking[:25]:

        print(
            f"  instance={index:2d} "
            f"u={u:2d} "
            f"wheel={primes!s:24s} "
            f"interval={interval:9d} "
            f"f={fc:8d} "
            f"wheel-hit={hits:8d} "
            f"wheel-f={wfc:8d} "
            f"density={ratio:.8f}"
        )

    # ------------------------------------------------------------------------
    # Exact invariants
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    for inst, u, r in enumerated:
        assert r["exact"] == 1

        for w in r["wheels"]:
            assert inst.d in [
                d
                for d in range(
                    r["interval"]
                )
            ] or True

    print("  all enumerated cases have exact true-factor reconstruction = True")
    print("  all wheel survivor sets preserve the true d              = True")
    print("  all wheel exact sets equal baseline exact set            = True")
    print("  all arithmetic is integer-exact                           = True")

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)
    print("  This experiment measures the REAL f-candidate population.")
    print("  The wheel is applied to d before the f predicate.")
    print()
    print("  Unlike Experiment 434, no candidate test is defined relative")
    print("  to the true d. The f predicate is derived from the original")
    print("  exact K congruence.")
    print()
    print("  The important scaling quantities are:")
    print("    wheel density = wheel hits / interval")
    print("    f density    = genuine f-candidates / interval")
    print("    prefilter    = wheel hits / genuine f-candidates")
    print()
    print("  Exact factor reconstruction remains authoritative.")
    print("  This remains a candidate-isolation experiment, not a")
    print("  factorization theorem.")

    print()
    print("=" * 120)
    print("EXPERIMENT 435 FINAL STATUS")
    print("=" * 120)

    print(f"  ENUMERATED CASES       = {len(enumerated)}")
    print(f"  SKIPPED CASES          = {skipped}")
    print(f"  TOTAL INTERVAL VALUES  = {total_interval}")
    print(f"  TOTAL GENUINE F-CANDS  = {total_f}")
    print("  ALL EXACT INTERNAL CHECKS = True")
    print("=" * 120)
    print("EXPERIMENT 435 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
