from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# =============================================================================
# EXPERIMENT 418
# COMPACT EXACT BIVARIATE (x,k) LATTICE / LLL AUDIT
# =============================================================================
#
# Purpose:
#
#   Given
#
#       N = p*q
#       S = p+q
#       d = N+1-2S
#       K = (-d^2 - 3N^2 + 6N + 1)/4
#
#   and partial information
#
#       K = K0 + k,       0 <= k < 2^u,
#
#   define
#
#       C0 = -4K0 - 3N^2 + 6N + 1
#       d0 = floor(sqrt(C0))
#       r  = C0 - d0^2
#
#   with
#
#       d = d0 + x
#
#   giving the exact integer relation
#
#       f(x,k) = x^2 + 2*d0*x + 4*k - r = 0.
#
# This experiment is deliberately conservative:
#
#   * no resultants
#   * no symbolic multivariate factorization
#   * no floating-point arithmetic in the exact identities
#   * LLL is only treated as a lattice diagnostic
#
# Experiment 418 is a COMPACT reporting version of Experiment 417.
#
# Instead of printing every intermediate value for every (instance, u),
# it prints:
#
#   1. one summary row per test
#   2. LLL metrics for each test
#   3. detailed data only for interesting cases
#   4. a global summary
#
# A reduced vector is never called a "recovery" merely because it is short.
# A relation is considered root-vanishing only when its reconstructed
# integral polynomial evaluates exactly to zero at (x_true, k_true).
# =============================================================================


# -----------------------------------------------------------------------------
# Optional LLL support
# -----------------------------------------------------------------------------

try:
    from fpylll import IntegerMatrix, LLL

    HAVE_FPYLLL = True
    FPYLLL_IMPORT_ERROR = None

except Exception as exc:
    HAVE_FPYLLL = False
    FPYLLL_IMPORT_ERROR = repr(exc)


# -----------------------------------------------------------------------------
# Exact polynomial infrastructure
# -----------------------------------------------------------------------------

Monomial = Tuple[int, int]  # (x_degree, k_degree)
Poly = Dict[Monomial, int]


def poly_add(a: Poly, b: Poly) -> Poly:
    out = dict(a)

    for monomial, coeff in b.items():
        out[monomial] = out.get(monomial, 0) + coeff

        if out[monomial] == 0:
            del out[monomial]

    return out


def poly_scale(a: Poly, c: int) -> Poly:
    if c == 0:
        return {}

    return {
        monomial: coeff * c
        for monomial, coeff in a.items()
        if coeff * c != 0
    }


def poly_mul_monomial(a: Poly, dx: int, dk: int) -> Poly:
    return {
        (i + dx, j + dk): coeff
        for (i, j), coeff in a.items()
    }


def poly_eval(a: Poly, x: int, k: int) -> int:
    total = 0

    for (i, j), coeff in a.items():
        total += coeff * (x ** i) * (k ** j)

    return total


def poly_to_vector(
    a: Poly,
    monomials: List[Monomial],
    X: int,
    Y: int,
) -> List[int]:
    """
    Bound-scaled coefficient embedding:

        x -> X*Xvar
        k -> Y*Kvar

    Therefore

        c*x^i*k^j
            -> c*X^i*Y^j.
    """

    return [
        a.get((i, j), 0) * (X ** i) * (Y ** j)
        for i, j in monomials
    ]


def vector_to_poly(
    vec: List[int],
    monomials: List[Monomial],
    X: int,
    Y: int,
) -> Tuple[Poly, bool]:
    """
    Attempt exact inversion of the bound-scaled embedding.

    Returns:
        (polynomial, reconstructible)

    reconstructible=True only when every nonzero coordinate maps back to
    an integer coefficient.
    """

    out: Poly = {}

    for value, (i, j) in zip(vec, monomials):
        scale = (X ** i) * (Y ** j)

        if scale == 0 or value % scale != 0:
            return {}, False

        coeff = value // scale

        if coeff:
            out[(i, j)] = coeff

    return out, True


# -----------------------------------------------------------------------------
# Exact instance generation
# -----------------------------------------------------------------------------

@dataclass
class Instance:
    p: int
    q: int
    N: int
    S: int
    d: int
    gap: int
    K: int


def make_instance(p: int, q: int) -> Instance:
    N = p * q
    S = p + q
    d = N + 1 - 2 * S
    gap = q - p

    numerator = (
        -d * d
        - 3 * N * N
        + 6 * N
        + 1
    )

    if numerator % 4 != 0:
        raise AssertionError(
            f"K numerator is not divisible by 4 for p={p}, q={q}"
        )

    K = numerator // 4

    # Fundamental exact identities.
    assert -4 * K - 3 * N * N + 6 * N + 1 == d * d
    assert N * N - N + K == S * (N + 1 - S)

    return Instance(
        p=p,
        q=q,
        N=N,
        S=S,
        d=d,
        gap=gap,
        K=K,
    )


# -----------------------------------------------------------------------------
# Partial-K data
# -----------------------------------------------------------------------------

@dataclass
class PartialK:
    u: int
    K0: int
    k: int
    C0: int
    d0: int
    r: int
    x_true: int


def partial_k_data(inst: Instance, u: int) -> PartialK:
    modulus = 1 << u

    # Floor-aligned decomposition:
    #
    #   K0 = floor(K / 2^u) * 2^u
    #   k  = K-K0.
    #
    # This works for negative K in exact integer arithmetic.

    K0 = (inst.K // modulus) * modulus
    k = inst.K - K0

    if not (0 <= k < modulus):
        raise AssertionError("invalid partial-K decomposition")

    if K0 + k != inst.K:
        raise AssertionError("K decomposition does not reconstruct K")

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

    if C0 < 0:
        raise AssertionError("C0 must be nonnegative")

    d0 = math.isqrt(C0)
    r = C0 - d0 * d0
    x_true = inst.d - d0

    residual = (
        x_true * x_true
        + 2 * d0 * x_true
        + 4 * k
        - r
    )

    if residual != 0:
        raise AssertionError(
            f"exact transformed equation failed: {residual}"
        )

    return PartialK(
        u=u,
        K0=K0,
        k=k,
        C0=C0,
        d0=d0,
        r=r,
        x_true=x_true,
    )


# -----------------------------------------------------------------------------
# Exact interval
# -----------------------------------------------------------------------------

@dataclass
class IntervalData:
    d_low: int
    d_high: int
    x_low: int
    x_high: int
    candidate_count: int


def exact_d_interval(pk: PartialK) -> IntervalData:
    modulus = 1 << pk.u

    # Since 0 <= k <= 2^u-1:
    #
    #   C0 - 4(2^u-1) <= d^2 <= C0.

    lo_sq = pk.C0 - 4 * (modulus - 1)
    hi_sq = pk.C0

    if lo_sq <= 0:
        d_low = 0
    else:
        d_low = math.isqrt(lo_sq)

        while d_low * d_low < lo_sq:
            d_low += 1

    d_high = math.isqrt(hi_sq)

    x_low = d_low - pk.d0
    x_high = d_high - pk.d0

    candidate_count = max(0, d_high - d_low + 1)

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        x_low=x_low,
        x_high=x_high,
        candidate_count=candidate_count,
    )


# -----------------------------------------------------------------------------
# Shifted polynomial family
# -----------------------------------------------------------------------------

@dataclass
class ShiftedPolynomial:
    i: int
    j: int
    poly: Poly


def build_base_polynomial(pk: PartialK) -> Poly:
    return {
        (2, 0): 1,
        (1, 0): 2 * pk.d0,
        (0, 1): 4,
        (0, 0): -pk.r,
    }


def build_shifted_family(
    base: Poly,
    max_i: int,
    max_j: int,
) -> List[ShiftedPolynomial]:

    out: List[ShiftedPolynomial] = []

    for i in range(max_i + 1):
        for j in range(max_j + 1):
            out.append(
                ShiftedPolynomial(
                    i=i,
                    j=j,
                    poly=poly_mul_monomial(base, i, j),
                )
            )

    return out


def monomial_basis(
    polys: List[ShiftedPolynomial],
) -> List[Monomial]:

    mons = set()

    for item in polys:
        mons.update(item.poly.keys())

    return sorted(mons)


# -----------------------------------------------------------------------------
# Lattice
# -----------------------------------------------------------------------------

@dataclass
class LatticeData:
    monomials: List[Monomial]
    shifted: List[ShiftedPolynomial]
    rows: List[List[int]]
    X: int
    Y: int


def build_lattice(
    pk: PartialK,
    x_bound: int,
    k_bound: int,
    max_i: int,
    max_j: int,
) -> LatticeData:

    X = max(1, abs(x_bound))
    Y = max(1, abs(k_bound))

    base = build_base_polynomial(pk)

    shifted = build_shifted_family(
        base,
        max_i=max_i,
        max_j=max_j,
    )

    monomials = monomial_basis(shifted)

    rows = [
        poly_to_vector(
            item.poly,
            monomials,
            X,
            Y,
        )
        for item in shifted
    ]

    return LatticeData(
        monomials=monomials,
        shifted=shifted,
        rows=rows,
        X=X,
        Y=Y,
    )


# -----------------------------------------------------------------------------
# Numeric exact diagnostics
# -----------------------------------------------------------------------------

def vector_inf_norm(v: List[int]) -> int:
    return max((abs(x) for x in v), default=0)


def vector_l2_squared(v: List[int]) -> int:
    return sum(x * x for x in v)


# -----------------------------------------------------------------------------
# LLL
# -----------------------------------------------------------------------------

@dataclass
class LLLResult:
    rows: List[List[int]]
    success: bool
    error: Optional[str]


def run_lll(rows: List[List[int]]) -> LLLResult:
    if not HAVE_FPYLLL:
        return LLLResult(
            rows=[],
            success=False,
            error="fpylll unavailable",
        )

    if not rows:
        return LLLResult(
            rows=[],
            success=True,
            error=None,
        )

    nr = len(rows)
    nc = len(rows[0])

    try:
        M = IntegerMatrix(nr, nc)

        for i in range(nr):
            for j in range(nc):
                M[i, j] = rows[i][j]

        LLL.reduction(M)

        reduced = [
            [int(M[i, j]) for j in range(nc)]
            for i in range(nr)
        ]

        return LLLResult(
            rows=reduced,
            success=True,
            error=None,
        )

    except Exception as exc:
        return LLLResult(
            rows=[],
            success=False,
            error=repr(exc),
        )


# -----------------------------------------------------------------------------
# Reduced-vector audit
# -----------------------------------------------------------------------------

@dataclass
class VectorAudit:
    index: int
    inf_norm: int
    l2_squared: int
    reconstructible: bool
    nonzero_terms: int
    root_value: Optional[int]
    vanishes: bool


def audit_reduced_vectors(
    reduced_rows: List[List[int]],
    lattice: LatticeData,
    pk: PartialK,
) -> List[VectorAudit]:

    audits: List[VectorAudit] = []

    for index, row in enumerate(reduced_rows):

        poly, reconstructible = vector_to_poly(
            row,
            lattice.monomials,
            lattice.X,
            lattice.Y,
        )

        if reconstructible:
            root_value = poly_eval(
                poly,
                pk.x_true,
                pk.k,
            )
            vanishes = root_value == 0
            nonzero_terms = len(poly)
        else:
            root_value = None
            vanishes = False
            nonzero_terms = 0

        audits.append(
            VectorAudit(
                index=index,
                inf_norm=vector_inf_norm(row),
                l2_squared=vector_l2_squared(row),
                reconstructible=reconstructible,
                nonzero_terms=nonzero_terms,
                root_value=root_value,
                vanishes=vanishes,
            )
        )

    return audits


# -----------------------------------------------------------------------------
# Exact candidate enumeration
# -----------------------------------------------------------------------------

@dataclass
class EnumerationResult:
    enumerated: bool
    candidates: List[int]


def enumerate_exact_candidates(
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> EnumerationResult:

    if interval.candidate_count > limit:
        return EnumerationResult(
            enumerated=False,
            candidates=[],
        )

    candidates: List[int] = []

    modulus = 1 << pk.u

    for d in range(interval.d_low, interval.d_high + 1):

        delta = pk.C0 - d * d

        if delta < 0 or delta % 4 != 0:
            continue

        k = delta // 4

        if not (0 <= k < modulus):
            continue

        x = d - pk.d0

        residual = (
            x * x
            + 2 * pk.d0 * x
            + 4 * k
            - pk.r
        )

        if residual == 0:
            candidates.append(d)

    return EnumerationResult(
        enumerated=True,
        candidates=candidates,
    )


# -----------------------------------------------------------------------------
# Per-test summary record
# -----------------------------------------------------------------------------

@dataclass
class TestSummary:
    instance_id: int
    bits: int
    p: int
    q: int
    gap: int
    true_x: int
    X: int
    Y: int
    candidate_count: int
    enumerated: bool
    candidates_found: int
    true_present: bool
    unique_recovery: bool
    equation_ok: bool
    lattice_rows: int
    lattice_cols: int
    raw_min_inf: int
    raw_max_inf: int
    lll_success: bool
    lll_min_inf: Optional[int]
    lll_max_inf: Optional[int]
    vanishing_rows: int
    integral_rows: int
    interesting: bool
    lll_error: Optional[str]


# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------

def compact_bool(value: bool) -> str:
    return "YES" if value else "NO"


def compact_signed(value: int) -> str:
    return f"{value:d}"


def print_test_table(results: List[TestSummary]) -> None:
    print()
    print("=" * 118)
    print("COMPACT TEST SUMMARY")
    print("=" * 118)

    header = (
        f"{'i':>2} "
        f"{'u':>3} "
        f"{'gap':>9} "
        f"{'x':>10} "
        f"{'X':>10} "
        f"{'cand':>9} "
        f"{'unique':>6} "
        f"{'LLL':>5} "
        f"{'vanish':>6} "
        f"{'shortest':>12}"
    )

    print(header)
    print("-" * len(header))

    for r in results:

        shortest = (
            str(r.lll_min_inf)
            if r.lll_min_inf is not None
            else "-"
        )

        print(
            f"{r.instance_id:2d} "
            f"{r.bits:3d} "
            f"{r.gap:9d} "
            f"{r.true_x:10d} "
            f"{r.X:10d} "
            f"{r.candidate_count:9d} "
            f"{compact_bool(r.unique_recovery):>6} "
            f"{compact_bool(r.lll_success):>5} "
            f"{r.vanishing_rows:6d} "
            f"{shortest:>12}"
        )


def print_lll_table(results: List[TestSummary]) -> None:
    print()
    print("=" * 118)
    print("LLL DIAGNOSTIC SUMMARY")
    print("=" * 118)

    header = (
        f"{'i':>2} "
        f"{'u':>3} "
        f"{'dim':>9} "
        f"{'raw min':>14} "
        f"{'raw max':>18} "
        f"{'LLL min':>18} "
        f"{'LLL max':>18} "
        f"{'vanish':>6}"
    )

    print(header)
    print("-" * len(header))

    for r in results:

        lll_min = (
            str(r.lll_min_inf)
            if r.lll_min_inf is not None
            else "-"
        )

        lll_max = (
            str(r.lll_max_inf)
            if r.lll_max_inf is not None
            else "-"
        )

        print(
            f"{r.instance_id:2d} "
            f"{r.bits:3d} "
            f"{r.lattice_rows:4d}x{r.lattice_cols:<4d} "
            f"{r.raw_min_inf:14d} "
            f"{r.raw_max_inf:18d} "
            f"{lll_min:>18} "
            f"{lll_max:>18} "
            f"{r.vanishing_rows:6d}"
        )


# -----------------------------------------------------------------------------
# Detailed output only for interesting cases
# -----------------------------------------------------------------------------

def print_interesting_case(
    summary: TestSummary,
    inst: Instance,
    pk: PartialK,
    interval: IntervalData,
    enumeration: EnumerationResult,
    audits: List[VectorAudit],
) -> None:

    print()
    print("=" * 118)
    print(
        f"DETAIL: INSTANCE {summary.instance_id}, "
        f"UNKNOWN K BITS = {summary.bits}"
    )
    print("=" * 118)

    print(
        f"p={inst.p}  q={inst.q}  N={inst.N}  "
        f"S={inst.S}  gap={inst.gap}"
    )

    print(
        f"K={inst.K}"
    )

    print()
    print("PARTIAL-K")
    print(f"  K0        = {pk.K0}")
    print(f"  k         = {pk.k}")
    print(f"  2^u       = {1 << pk.u}")

    print()
    print("TRANSFORMED ROOT")
    print(f"  d0        = {pk.d0}")
    print(f"  x_true    = {pk.x_true}")
    print(f"  r         = {pk.r}")

    print()
    print("INTERVAL")
    print(
        f"  d = [{interval.d_low}, {interval.d_high}]"
    )
    print(
        f"  x = [{interval.x_low}, {interval.x_high}]"
    )
    print(
        f"  candidates in interval = {interval.candidate_count}"
    )

    if enumeration.enumerated:
        print()
        print("ENUMERATION")

        if len(enumeration.candidates) <= 30:
            print(
                f"  candidates = {enumeration.candidates}"
            )
        else:
            print(
                f"  candidates found = "
                f"{len(enumeration.candidates)}"
            )

        print(
            f"  true present = "
            f"{inst.d in enumeration.candidates}"
        )

        print(
            f"  unique = "
            f"{len(enumeration.candidates) == 1 and enumeration.candidates[0] == inst.d}"
        )

    else:
        print()
        print("ENUMERATION")
        print("  skipped: interval exceeds limit")

    print()
    print("LLL")

    if not summary.lll_success:
        print(
            f"  status = FAILED"
        )
        print(
            f"  error  = {summary.lll_error}"
        )
        return

    print(
        f"  shortest inf-norm = {summary.lll_min_inf}"
    )
    print(
        f"  longest inf-norm  = {summary.lll_max_inf}"
    )
    print(
        f"  integral reduced rows = {summary.integral_rows}"
    )
    print(
        f"  exact root-vanishing rows = {summary.vanishing_rows}"
    )

    interesting_audits = [
        a
        for a in audits
        if a.vanishes or a.reconstructible
    ]

    if interesting_audits:
        print()
        print("REDUCED-VECTOR DETAILS")

        for audit in interesting_audits[:10]:
            print(
                f"  row={audit.index:2d} "
                f"inf={audit.inf_norm} "
                f"terms={audit.nonzero_terms:2d} "
                f"integral={compact_bool(audit.reconstructible):3s} "
                f"value={audit.root_value} "
                f"vanishes={compact_bool(audit.vanishes)}"
            )


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 118)
    print("EXPERIMENT 418 START")
    print("=" * 118)

    print()
    print("COMPACT EXACT BIVARIATE (x,k) LATTICE / LLL AUDIT")

    print()
    print("CORE RELATION")
    print("  f(x,k) = x^2 + 2*d0*x + 4*k - r = 0")

    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no symbolic multivariate factorization")
    print("  LLL = diagnostic only")
    print("  short vector != recovery")
    print("  exact polynomial vanishing is checked separately")

    unknown_bits = [
        20,
        24,
        28,
        32,
        36,
        40,
        44,
        48,
    ]

    max_i = 2
    max_j = 2
    enumerate_limit = 200_000

    instances = [
        (50411, 282599),
        (1013, 10009),
        (10009, 1000033),
        (10009, 10037),
        (50023, 50051),
        (100019, 100043),
        (200009, 200017),
        (300017, 900007),
    ]

    print()
    print("CONFIGURATION")
    print(f"  instances      = {len(instances)}")
    print(f"  K bits         = {unknown_bits}")
    print(f"  x shifts       = 0..{max_i}")
    print(f"  k shifts       = 0..{max_j}")
    print(f"  lattice rows   = {(max_i + 1) * (max_j + 1)}")
    print(f"  enumeration max= {enumerate_limit}")
    print(f"  fpylll         = {HAVE_FPYLLL}")

    if not HAVE_FPYLLL:
        print(f"  import error   = {FPYLLL_IMPORT_ERROR}")

    # Configuration is constant across all tests.
    results: List[TestSummary] = []
    details = []

    for instance_id, (p, q) in enumerate(instances, start=1):

        inst = make_instance(p, q)

        for u in unknown_bits:

            pk = partial_k_data(inst, u)

            equation_ok = (
                pk.x_true * pk.x_true
                + 2 * pk.d0 * pk.x_true
                + 4 * pk.k
                - pk.r
                == 0
            )

            interval = exact_d_interval(pk)

            enumeration = enumerate_exact_candidates(
                pk,
                interval,
                enumerate_limit,
            )

            candidates_found = len(enumeration.candidates)

            true_present = (
                enumeration.enumerated
                and inst.d in enumeration.candidates
            )

            unique_recovery = (
                enumeration.enumerated
                and len(enumeration.candidates) == 1
                and enumeration.candidates[0] == inst.d
            )

            X = max(
                abs(interval.x_low),
                abs(interval.x_high),
                1,
            )

            Y = max((1 << u) - 1, 1)

            lattice = build_lattice(
                pk,
                x_bound=X,
                k_bound=Y,
                max_i=max_i,
                max_j=max_j,
            )

            raw_norms = [
                vector_inf_norm(row)
                for row in lattice.rows
            ]

            raw_min_inf = min(raw_norms)
            raw_max_inf = max(raw_norms)

            lll_result = run_lll(lattice.rows)

            lll_min_inf: Optional[int] = None
            lll_max_inf: Optional[int] = None
            vanishing_rows = 0
            integral_rows = 0
            audits: List[VectorAudit] = []

            if lll_result.success:

                reduced = lll_result.rows

                reduced_norms = [
                    vector_inf_norm(row)
                    for row in reduced
                ]

                if reduced_norms:
                    lll_min_inf = min(reduced_norms)
                    lll_max_inf = max(reduced_norms)

                audits = audit_reduced_vectors(
                    reduced,
                    lattice,
                    pk,
                )

                integral_rows = sum(
                    1
                    for audit in audits
                    if audit.reconstructible
                )

                vanishing_rows = sum(
                    1
                    for audit in audits
                    if audit.vanishes
                )

            interesting = (
                pk.x_true != 0
                or interval.candidate_count > 1
                or not unique_recovery
                or vanishing_rows > 0
                or (
                    lll_result.success
                    and lll_min_inf is not None
                    and lll_min_inf < raw_min_inf
                )
            )

            summary = TestSummary(
                instance_id=instance_id,
                bits=u,
                p=inst.p,
                q=inst.q,
                gap=inst.gap,
                true_x=pk.x_true,
                X=X,
                Y=Y,
                candidate_count=interval.candidate_count,
                enumerated=enumeration.enumerated,
                candidates_found=candidates_found,
                true_present=true_present,
                unique_recovery=unique_recovery,
                equation_ok=equation_ok,
                lattice_rows=len(lattice.rows),
                lattice_cols=len(lattice.monomials),
                raw_min_inf=raw_min_inf,
                raw_max_inf=raw_max_inf,
                lll_success=lll_result.success,
                lll_min_inf=lll_min_inf,
                lll_max_inf=lll_max_inf,
                vanishing_rows=vanishing_rows,
                integral_rows=integral_rows,
                interesting=interesting,
                lll_error=lll_result.error,
            )

            results.append(summary)

            if interesting:
                details.append(
                    (
                        summary,
                        inst,
                        pk,
                        interval,
                        enumeration,
                        audits,
                    )
                )

    # -------------------------------------------------------------------------
    # Compact global output
    # -------------------------------------------------------------------------

    print_test_table(results)
    print_lll_table(results)

    # -------------------------------------------------------------------------
    # Interesting cases
    # -------------------------------------------------------------------------

    print()
    print("=" * 118)
    print("INTERESTING CASES")
    print("=" * 118)

    print(
        "Detailed output is shown only when x != 0, the interval has "
        "multiple candidates, unique enumeration fails, or LLL finds an "
        "interesting reduced relation."
    )

    if not details:
        print("  none")
    else:
        for (
            summary,
            inst,
            pk,
            interval,
            enumeration,
            audits,
        ) in details:

            print_interesting_case(
                summary,
                inst,
                pk,
                interval,
                enumeration,
                audits,
            )

    # -------------------------------------------------------------------------
    # Global checks
    # -------------------------------------------------------------------------

    print()
    print("=" * 118)
    print("GLOBAL EXACT-CHECK SUMMARY")
    print("=" * 118)

    all_equations = all(r.equation_ok for r in results)

    all_true_inside = all(
        r.true_present
        for r in results
        if r.enumerated
    )

    print(
        f"  transformed equation identities = "
        f"{all_equations}"
    )

    print(
        f"  all enumerated intervals contain true d = "
        f"{all_true_inside}"
    )

    print(
        f"  total tests = {len(results)}"
    )

    print(
        f"  interesting tests = "
        f"{sum(r.interesting for r in results)}"
    )

    lll_success_count = sum(
        r.lll_success
        for r in results
    )

    print(
        f"  successful LLL reductions = "
        f"{lll_success_count}/{len(results)}"
    )

    total_vanishing = sum(
        r.vanishing_rows
        for r in results
    )

    print(
        f"  total exact root-vanishing reduced rows = "
        f"{total_vanishing}"
    )

    unique_count = sum(
        r.unique_recovery
        for r in results
    )

    print(
        f"  unique exact interval recoveries = "
        f"{unique_count}/{len(results)}"
    )

    # -------------------------------------------------------------------------
    # Aggregate by unknown K bits
    # -------------------------------------------------------------------------

    print()
    print("=" * 118)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 118)

    header = (
        f"{'u':>3} "
        f"{'tests':>5} "
        f"{'unique':>7} "
        f"{'multi':>7} "
        f"{'nonzero x':>10} "
        f"{'LLL':>7} "
        f"{'vanish':>8}"
    )

    print(header)
    print("-" * len(header))

    for u in unknown_bits:

        group = [
            r
            for r in results
            if r.bits == u
        ]

        print(
            f"{u:3d} "
            f"{len(group):5d} "
            f"{sum(r.unique_recovery for r in group):7d} "
            f"{sum(r.candidate_count > 1 for r in group):7d} "
            f"{sum(r.true_x != 0 for r in group):10d} "
            f"{sum(r.lll_success for r in group):7d} "
            f"{sum(r.vanishing_rows for r in group):8d}"
        )

    # -------------------------------------------------------------------------
    # Interpretation
    # -------------------------------------------------------------------------

    print()
    print("=" * 118)
    print("INTERPRETATION")
    print("=" * 118)

    print(
        "  The relation f(x,k)=0 is an exact integer equation."
    )

    print(
        "  LLL is therefore being evaluated as a lattice heuristic, "
        "not as a demonstrated Coppersmith theorem."
    )

    print(
        "  A reduced vector is interesting only when it also produces "
        "an exact integral polynomial relation."
    )

    print(
        "  Even multiple root-vanishing relations do not by themselves "
        "constitute a proof of factor recovery."
    )

    print()
    print("=" * 118)
    print("EXPERIMENT 418 FINAL STATUS")
    print("=" * 118)

    if all_equations:
        print("  ALL INTERNAL EXACT EQUATION CHECKS PASS = True")
    else:
        print("  ALL INTERNAL EXACT EQUATION CHECKS PASS = False")

    if lll_success_count == len(results):
        print("  ALL LLL REDUCTIONS SUCCESSFUL = True")
    else:
        print(
            "  ALL LLL REDUCTIONS SUCCESSFUL = "
            f"{lll_success_count == len(results)}"
        )

    print("=" * 118)
    print("EXPERIMENT 418 FINISHED")
    print("=" * 118)


if __name__ == "__main__":
    run_experiment()