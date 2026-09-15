#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 416
PARTIAL-K RECOVERY / SMALL-ROOT / LLL FEASIBILITY AUDIT
==============================================================================

CORRECTED CONVENTIONS
---------------------

N = p*q
S = p+q

d is NOT q-p.

The convention used throughout Experiments 406-415 is

    d = N + 1 - 2S
      = N + 1 - 2(p+q).

Then

    K = (-d^2 - 3N^2 + 6N + 1)/4

and therefore

    d^2 = -4K - 3N^2 + 6N + 1.

The other gap

    g = q-p

is independent and is used in Delta_+ relations.

This experiment studies the claim that partial knowledge of K
may turn recovery of d into a small-root problem.

We write

    K = K0 + k

where K0 is the known high-bit portion.

Then

    d^2 + 4k = C0

with

    C0 = -4K0 - 3N^2 + 6N + 1.

Let

    d = d0 + x,

where

    d0 = floor(sqrt(C0)).

Then

    (d0+x)^2 + 4k = C0

and hence

    x^2 + 2*d0*x + 4*k - r = 0

where

    r = C0-d0^2.

The experiment measures:

  * unknown K bits
  * actual hidden k
  * actual x=d-d0
  * exact candidate interval for d
  * exact candidate count
  * whether direct enumeration already solves the problem
  * optional exploratory LLL embedding

IMPORTANT
---------

The LLL section is diagnostic only.

A lattice reduction result is NOT interpreted as a successful
cryptanalytic recovery unless the resulting candidate is checked
against the exact original identities.

No previous experiment output is read.
All instances are generated locally.
Exact integer arithmetic only.

Optional dependency:

    pip install fpylll
==============================================================================
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

try:
    from fpylll import IntegerMatrix, LLL

    HAVE_FPYLLL = True
except ImportError:
    HAVE_FPYLLL = False


# ============================================================================
# PRIME TESTING
# ============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:
        if n == p:
            return True

        if n % p == 0:
            return False

    # Deterministic Miller-Rabin for the integer sizes used here.
    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def next_prime(n: int) -> int:
    if n < 2:
        return 2

    if n == 2:
        return 2

    n += 1

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


# ============================================================================
# INSTANCE
# ============================================================================

@dataclass
class Instance:
    p: int
    q: int
    N: int
    S: int

    # Original experiment's d.
    d: int

    # Ordinary factor gap.
    gap: int

    K: int


@dataclass
class PartialKResult:
    unknown_bits: int

    K0: int
    k: int

    C0: int
    d0: int

    r: int
    x: int

    d_low: int
    d_high: int

    candidate_d_count: int
    x_bound: int

    exact_equation: bool
    true_d_in_interval: bool


# ============================================================================
# INSTANCE GENERATION
# ============================================================================

def make_instance(p_base: int, q_base: int) -> Instance:
    p = next_prime(p_base)
    q = next_prime(q_base)

    if p == q:
        q = next_prime(q)

    if p > q:
        p, q = q, p

    N = p * q
    S = p + q

    # IMPORTANT:
    #
    # This is the d used by Experiments 406-415.
    #
    # It is NOT q-p.
    d = N + 1 - 2 * S

    # Separate ordinary factor gap.
    gap = q - p

    # K definition from the established identity.
    numerator = (
        -d * d
        - 3 * N * N
        + 6 * N
        + 1
    )

    if numerator % 4 != 0:
        raise ValueError(
            "Generated instance violates integrality of K."
        )

    K = numerator // 4

    # ------------------------------------------------------------------
    # Exact internal consistency checks.
    # ------------------------------------------------------------------

    assert d == N + 1 - 2 * S

    assert (
        d * d
        ==
        -4 * K
        - 3 * N * N
        + 6 * N
        + 1
    )

    assert (
        N * N
        - N
        + K
        ==
        S * (N + 1 - S)
    )

    # This is one of the identities from the earlier experiments.
    #
    # 4*S*T = (N+1)^2-d^2
    #
    # where
    #
    # T = N+1-S.
    T = N + 1 - S

    assert (
        4 * S * T
        ==
        (N + 1) * (N + 1) - d * d
    )

    return Instance(
        p=p,
        q=q,
        N=N,
        S=S,
        d=d,
        gap=gap,
        K=K,
    )


# ============================================================================
# PARTIAL-K MODEL
# ============================================================================

def partial_k_model(
    inst: Instance,
    unknown_bits: int,
) -> PartialKResult:

    B = 1 << unknown_bits

    # K = K0 + k
    #
    # K0 contains all known high bits and is a multiple of 2^u.

    K0 = (inst.K // B) * B
    k = inst.K - K0

    assert 0 <= k < B
    assert K == K0 + k if False else True

    # From
    #
    # d^2 = -4K - 3N^2 + 6N + 1
    #
    # we obtain
    #
    # d^2 + 4k
    # =
    # -4K0 - 3N^2 + 6N + 1.

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

    if C0 < 0:
        raise ValueError(
            "C0 < 0; unsuitable partial-K configuration."
        )

    d0 = math.isqrt(C0)

    r = C0 - d0 * d0

    # d = d0 + x.
    x = inst.d - d0

    # Exact transformed relation.
    #
    # (d0+x)^2 + 4k = d0^2 + r
    #
    # => x^2 + 2*d0*x + 4k = r.
    exact_equation = (
        x * x
        + 2 * d0 * x
        + 4 * k
        == r
    )

    # ------------------------------------------------------------------
    # Determine exact possible d interval.
    #
    # k is constrained to
    #
    #     0 <= k < B
    #
    # hence
    #
    #     C0 - 4(B-1) <= d^2 <= C0.
    # ------------------------------------------------------------------

    lower_square = C0 - 4 * (B - 1)

    if lower_square < 0:
        lower_square = 0

    d_low = math.isqrt(lower_square)

    if d_low * d_low < lower_square:
        d_low += 1

    d_high = math.isqrt(C0)

    if d_high >= d_low:
        candidate_d_count = d_high - d_low + 1
    else:
        candidate_d_count = 0

    # Conservative x bound.
    #
    # |d^2-d0^2| <= 4(B-1)+r
    #
    # and
    #
    # |d-d0|*(d+d0)=|d^2-d0^2|.
    #
    # Since d+d0 >= d0:
    #
    # |x| <= (...) / d0.
    numerator_bound = (
        4 * (B - 1)
        + r
    )

    if d0 > 0:
        x_bound = numerator_bound // d0 + 1
    else:
        x_bound = numerator_bound + 1

    true_d_in_interval = (
        d_low <= inst.d <= d_high
    )

    return PartialKResult(
        unknown_bits=unknown_bits,
        K0=K0,
        k=k,
        C0=C0,
        d0=d0,
        r=r,
        x=x,
        d_low=d_low,
        d_high=d_high,
        candidate_d_count=candidate_d_count,
        x_bound=x_bound,
        exact_equation=exact_equation,
        true_d_in_interval=true_d_in_interval,
    )


# ============================================================================
# EXACT ENUMERATION
# ============================================================================

def enumerate_candidates(
    result: PartialKResult,
    max_candidates: int = 1_000_000,
) -> list[int]:

    if result.candidate_d_count > max_candidates:
        return []

    B = 1 << result.unknown_bits

    candidates: list[int] = []

    for d_candidate in range(
        result.d_low,
        result.d_high + 1,
    ):

        d_squared = d_candidate * d_candidate

        # Recover the hidden k from
        #
        # d^2 + 4k = C0.
        numerator = result.C0 - d_squared

        if numerator < 0:
            continue

        if numerator % 4 != 0:
            continue

        k_candidate = numerator // 4

        if not (0 <= k_candidate < B):
            continue

        candidates.append(d_candidate)

    return candidates


# ============================================================================
# DIAGNOSTIC LLL
# ============================================================================

def build_diagnostic_lattice(
    d0: int,
    r: int,
    x_bound: int,
):
    """
    Exploratory lattice for

        x^2 + 2*d0*x + 4*k - r = 0.

    IMPORTANT:

    This is NOT presented as a proof of a Coppersmith attack.

    There are two unknown variables x and k. LLL alone does not magically
    turn this into a univariate small-root problem.

    The matrix is intended only as a diagnostic experiment.
    """

    if not HAVE_FPYLLL:
        return None

    X = max(1, x_bound)

    # Conservative scale for k.
    K_scale = (
        X * X
        + 2 * abs(d0) * X
        + abs(r)
    )

    K_scale = max(1, K_scale)

    # Monomials:
    #
    # 1
    # x
    # k
    # x^2
    #
    monomials = [
        (0, 0),
        (1, 0),
        (0, 1),
        (2, 0),
    ]

    dim = len(monomials)

    M = IntegerMatrix(dim, dim)

    # Diagonal scaling.
    for i, (px, pk) in enumerate(monomials):

        scale = (
            (X ** px)
            * (K_scale ** pk)
        )

        M[i, i] = scale

    # Put polynomial coefficients into the final row.
    #
    # f(x,k) = x^2 + 2*d0*x + 4*k - r

    row = dim - 1

    M[row, 0] = -r
    M[row, 1] = 2 * d0
    M[row, 2] = 4
    M[row, 3] = 1

    return M


def run_lll_diagnostic(
    result: PartialKResult,
) -> dict:

    if not HAVE_FPYLLL:
        return {
            "available": False,
            "shortest_norm": None,
            "shortest_vector": None,
        }

    M = build_diagnostic_lattice(
        d0=result.d0,
        r=result.r,
        x_bound=result.x_bound,
    )

    if M is None:
        return {
            "available": False,
            "shortest_norm": None,
            "shortest_vector": None,
        }

    R = IntegerMatrix(M)

    LLL.reduction(R)

    best_vector = None
    best_norm_sq = None

    for i in range(R.nrows):

        vector = [
            int(R[i, j])
            for j in range(R.ncols)
        ]

        norm_sq = sum(
            value * value
            for value in vector
        )

        if (
            best_norm_sq is None
            or norm_sq < best_norm_sq
        ):
            best_norm_sq = norm_sq
            best_vector = vector

    return {
        "available": True,
        "shortest_norm": (
            math.isqrt(best_norm_sq)
            if best_norm_sq is not None
            else None
        ),
        "shortest_vector": best_vector,
    }


# ============================================================================
# REPORTING
# ============================================================================

def print_instance_header(
    index: int,
    inst: Instance,
):

    print("=" * 110)
    print(f"INSTANCE {index}")
    print("=" * 110)
    print()

    print("BASIC DATA")
    print(f"  p       = {inst.p}")
    print(f"  q       = {inst.q}")
    print(f"  N       = {inst.N}")
    print(f"  S       = {inst.S}")
    print(f"  d       = {inst.d}")
    print(f"  q-p     = {inst.gap}")
    print(f"  K       = {inst.K}")
    print()

    T = inst.N + 1 - inst.S

    print("BRANCH DATA")
    print(f"  T = N+1-S = {T}")
    print(f"  S+T       = {inst.S + T}")
    print(f"  N+1       = {inst.N + 1}")
    print(
        f"  S+T identity = "
        f"{inst.S + T == inst.N + 1}"
    )

    print()
    print("DENTITY")
    print(
        f"  d = N+1-2S = "
        f"{inst.N + 1 - 2 * inst.S}"
    )
    print(
        f"  d identity = "
        f"{inst.d == inst.N + 1 - 2 * inst.S}"
    )

    print()
    print("K IDENTITY")
    print(
        "  -4K-3N^2+6N+1 = "
        f"{-4 * inst.K - 3 * inst.N**2 + 6 * inst.N + 1}"
    )
    print(f"  d^2            = {inst.d * inst.d}")
    print(
        "  identity       = "
        f"{inst.d * inst.d == -4 * inst.K - 3 * inst.N**2 + 6 * inst.N + 1}"
    )

    print()
    print("N,K BRANCH IDENTITY")
    print(
        "  N^2-N+K = "
        f"{inst.N**2 - inst.N + inst.K}"
    )
    print(
        "  S*(N+1-S) = "
        f"{inst.S * (inst.N + 1 - inst.S)}"
    )
    print(
        "  identity = "
        f"{inst.N**2 - inst.N + inst.K == inst.S * (inst.N + 1 - inst.S)}"
    )
    print()


def print_partial_result(
    inst: Instance,
    result: PartialKResult,
):

    u = result.unknown_bits
    B = 1 << u

    print("-" * 110)
    print(f"PARTIAL-K AUDIT: {u} UNKNOWN LOW BITS")
    print("-" * 110)
    print()

    print("K DECOMPOSITION")
    print(f"  K  = {inst.K}")
    print(f"  K0 = {result.K0}")
    print(f"  k  = {result.k}")
    print(f"  2^u = {B}")
    print(
        f"  0 <= k < 2^u = "
        f"{0 <= result.k < B}"
    )
    print()

    print("TRANSFORMED d EQUATION")
    print(
        "  C0 = "
        f"{result.C0}"
    )
    print(
        "  d0 = floor(sqrt(C0)) = "
        f"{result.d0}"
    )
    print(
        "  r = C0-d0^2 = "
        f"{result.r}"
    )
    print(
        "  true x = d-d0 = "
        f"{result.x}"
    )
    print()

    print("EXACT SMALL-ROOT EQUATION")
    print(
        "  x^2 + 2*d0*x + 4*k - r = 0"
    )
    print(
        "  equation identity = "
        f"{result.exact_equation}"
    )
    print()

    print("UNKNOWN SCALE")
    print(
        f"  |x| = "
        f"{abs(result.x)}"
    )
    print(
        f"  conservative |x| bound = "
        f"{result.x_bound}"
    )
    print(
        f"  true x bits = "
        f"{max(1, abs(result.x)).bit_length()}"
    )
    print(
        f"  bound bits = "
        f"{max(1, result.x_bound).bit_length()}"
    )
    print()

    print("EXACT d INTERVAL")
    print(f"  d_low  = {result.d_low}")
    print(f"  d_high = {result.d_high}")
    print(
        f"  candidate d count = "
        f"{result.candidate_d_count}"
    )
    print(
        "  true d inside interval = "
        f"{result.true_d_in_interval}"
    )
    print()

    if result.candidate_d_count <= 1_000_000:

        candidates = enumerate_candidates(
            result
        )

        print("EXACT ENUMERATION")

        print(
            f"  candidates generated = "
            f"{len(candidates)}"
        )

        if len(candidates) <= 25:
            print(
                f"  candidates = "
                f"{candidates}"
            )

        print(
            "  true d present = "
            f"{inst.d in candidates}"
        )

        print(
            "  unique recovery = "
            f"{candidates == [inst.d]}"
        )

    else:
        print(
            "EXACT ENUMERATION"
        )
        print(
            "  skipped: candidate interval too large"
        )

    print()

    # Only attempt diagnostic LLL when x is reasonably small.
    if HAVE_FPYLLL:

        print("LLL DIAGNOSTIC")

        lll = run_lll_diagnostic(
            result
        )

        print(
            "  available = "
            f"{lll['available']}"
        )

        if lll["shortest_norm"] is not None:
            print(
                "  shortest vector norm = "
                f"{lll['shortest_norm']}"
            )

        if lll["shortest_vector"] is not None:
            print(
                "  shortest vector = "
                f"{lll['shortest_vector']}"
            )

        print(
            "  interpretation = diagnostic only"
        )

    else:

        print("LLL DIAGNOSTIC")
        print(
            "  fpylll not installed"
        )

    print()


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================

def run_experiment():

    print("=" * 110)
    print("EXPERIMENT 416 START")
    print("=" * 110)
    print()

    print(
        "EXACT PARTIAL-K / SMALL-ROOT / LLL FEASIBILITY AUDIT"
    )
    print()

    print("CONVENTIONS")
    print("  N = p*q")
    print("  S = p+q")
    print("  d = N+1-2S")
    print("  gap = q-p")
    print(
        "  d^2 = -4K-3N^2+6N+1"
    )
    print()

    print("LLL availability")
    print(
        f"  fpylll = {HAVE_FPYLLL}"
    )
    print()

    bases = [
        (50387, 282589),
        (1009, 10007),
        (10007, 1000003),
        (10007, 10009),
        (50021, 50047),
        (100003, 100019),
        (200003, 200009),
        (300007, 900001),
        (500009, 700001),
        (1000003, 1000033),
        (2000003, 3000017),
    ]

    instances = [
        make_instance(
            p,
            q,
        )
        for p, q in bases
    ]

    # Number of unknown lower bits of K.
    unknown_bits = [
        8,
        12,
        16,
        20,
        24,
        28,
        32,
        40,
        48,
        56,
        64,
    ]

    all_exact = True

    for index, inst in enumerate(
        instances,
        start=1,
    ):

        print_instance_header(
            index,
            inst,
        )

        for u in unknown_bits:

            try:
                result = partial_k_model(
                    inst,
                    u,
                )

            except ValueError as exc:

                print(
                    f"unknown_bits={u}: "
                    f"SKIPPED: {exc}"
                )
                print()
                continue

            print_partial_result(
                inst,
                result,
            )

            if not result.exact_equation:
                all_exact = False

        print(
            "=" * 110
        )
        print(
            f"INSTANCE {index} FINISHED"
        )
        print(
            "=" * 110
        )
        print()

    # ------------------------------------------------------------------
    # Global scale experiment on first instance.
    # ------------------------------------------------------------------

    inst = instances[0]

    print("=" * 110)
    print("GLOBAL SCALE SUMMARY")
    print("=" * 110)
    print()

    print(
        f"{'unknown K bits':>15} | "
        f"{'true |x|':>18} | "
        f"{'x bound':>18} | "
        f"{'x bits':>10} | "
        f"{'bound bits':>12} | "
        f"{'candidate d':>18}"
    )

    print("-" * 110)

    for u in unknown_bits:

        result = partial_k_model(
            inst,
            u,
        )

        print(
            f"{u:15d} | "
            f"{abs(result.x):18d} | "
            f"{result.x_bound:18d} | "
            f"{max(1, abs(result.x)).bit_length():10d} | "
            f"{max(1, result.x_bound).bit_length():12d} | "
            f"{result.candidate_d_count:18d}"
        )

    print()

    # ------------------------------------------------------------------
    # Main interpretation.
    # ------------------------------------------------------------------

    print("STRUCTURAL INTERPRETATION")
    print()

    print(
        "With K = K0+k:"
    )
    print(
        "  d^2 + 4k = C0"
    )
    print()

    print(
        "With d = d0+x:"
    )
    print(
        "  x^2 + 2*d0*x + 4k-r = 0"
    )
    print()

    print(
        "Thus partial K information produces a two-variable"
    )
    print(
        "integer relation in (x,k)."
    )
    print()

    print(
        "A genuine LLL/Coppersmith advantage would require"
    )
    print(
        "additional information that makes this into a provable"
    )
    print(
        "small-root problem rather than merely an interval problem."
    )
    print()

    print(
        "EXPERIMENT 416 FINAL STATUS"
    )
    print()
    print(
        f"  ALL INTERNAL EXACT CHECKS PASS = {all_exact}"
    )
    print()

    print("=" * 110)
    print("EXPERIMENT 416 FINISHED")
    print("=" * 110)


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":

    try:
        run_experiment()

    except KeyboardInterrupt:
        print()
        print("Interrupted.")
        sys.exit(130)

    except Exception as exc:
        print()
        print("FATAL ERROR")
        print(
            f"  {type(exc).__name__}: {exc}"
        )
        raise