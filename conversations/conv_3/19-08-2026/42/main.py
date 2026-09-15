from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

# =============================================================================
# EXPERIMENT 419
# EXACT LLL RELATION-NOVELTY AUDIT
# =============================================================================
#
# Purpose
# -------
# Experiment 418 showed that all reduced lattice vectors vanished at the
# hidden root (x_true, k_true). However, the lattice was generated entirely
# from
#
#       f(x,k) = x^2 + 2*d0*x + 4*k-r
#
# and its monomial multiples.
#
# Therefore a reduced vector may simply be another element of the ideal
#
#       <f>.
#
# Experiment 419 explicitly distinguishes:
#
#   1. INHERITED_MULTIPLE
#      P(x,k) is exactly divisible by f(x,k).
#
#   2. NEW_ROOT_VANISHING
#      P(x,k) is NOT divisible by f(x,k), but
#      P(x_true,k_true) = 0.
#
#   3. NONVANISHING
#      P(x_true,k_true) != 0.
#
# This is a purely exact integer-polynomial diagnostic.
#
# No resultants.
# No symbolic multivariate factorization.
# No floating-point arithmetic.
#
# IMPORTANT:
# Even a NEW_ROOT_VANISHING relation is NOT automatically a factoring
# algorithm. It only demonstrates that the LLL basis produced a relation
# outside the obvious ideal <f>.
# =============================================================================


# -----------------------------------------------------------------------------
# Exact polynomial infrastructure
# -----------------------------------------------------------------------------

Monomial = Tuple[int, int]          # (degree_x, degree_k)
Poly = Dict[Monomial, int]


def poly_clean(a: Poly) -> Poly:
    return {m: c for m, c in a.items() if c != 0}


def poly_add(a: Poly, b: Poly) -> Poly:
    out = dict(a)

    for m, c in b.items():
        out[m] = out.get(m, 0) + c

    return poly_clean(out)


def poly_sub(a: Poly, b: Poly) -> Poly:
    out = dict(a)

    for m, c in b.items():
        out[m] = out.get(m, 0) - c

    return poly_clean(out)


def poly_scale(a: Poly, c: int) -> Poly:
    if c == 0:
        return {}

    return {
        m: v * c
        for m, v in a.items()
        if v * c != 0
    }


def poly_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}

    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            m = (i1 + i2, j1 + j2)
            out[m] = out.get(m, 0) + c1 * c2

    return poly_clean(out)


def poly_mul_monomial(
    a: Poly,
    dx: int,
    dk: int,
) -> Poly:
    return {
        (i + dx, j + dk): c
        for (i, j), c in a.items()
    }


def poly_eval(
    a: Poly,
    x: int,
    k: int,
) -> int:
    total = 0

    for (i, j), c in a.items():
        total += c * (x ** i) * (k ** j)

    return total


def poly_degree_x(a: Poly) -> int:
    return max(
        (i for i, _ in a.keys()),
        default=-1,
    )


def poly_degree_k(a: Poly) -> int:
    return max(
        (j for _, j in a.keys()),
        default=-1,
    )


# -----------------------------------------------------------------------------
# Exact division by the known base polynomial
# -----------------------------------------------------------------------------
#
# We divide P(x,k) by
#
#   f(x,k) = x^2 + 2*d0*x + 4*k-r
#
# treating P as a polynomial in x whose coefficients are polynomials in k.
#
# Since f is monic in x, exact polynomial division is straightforward.
#
# P = Q*f + R
#
# with deg_x(R) < 2.
#
# P belongs to <f> iff R == 0.
# -----------------------------------------------------------------------------

def poly_as_x_coefficients(a: Poly) -> Dict[int, Dict[int, int]]:
    out: Dict[int, Dict[int, int]] = {}

    for (i, j), c in a.items():
        out.setdefault(i, {})
        out[i][j] = out[i].get(j, 0) + c

    out = {
        i: poly_clean(coeffs)
        for i, coeffs in out.items()
        if poly_clean(coeffs)
    }

    return out


def x_coefficients_to_poly(
    a: Dict[int, Dict[int, int]]
) -> Poly:
    out: Poly = {}

    for i, coeffs in a.items():
        for j, c in coeffs.items():
            if c:
                out[(i, j)] = c

    return out


def kpoly_add(
    a: Dict[int, int],
    b: Dict[int, int],
) -> Dict[int, int]:
    out = dict(a)

    for j, c in b.items():
        out[j] = out.get(j, 0) + c

    return {
        j: c
        for j, c in out.items()
        if c != 0
    }


def kpoly_scale(
    a: Dict[int, int],
    c: int,
) -> Dict[int, int]:
    if c == 0:
        return {}

    return {
        j: v * c
        for j, v in a.items()
        if v * c != 0
    }


def exact_divide_by_f(
    p: Poly,
    d0: int,
    r: int,
) -> Tuple[Poly, Poly]:
    """
    Exact Euclidean division in x by

        f = x^2 + 2*d0*x + 4*k-r.

    Returns:
        quotient, remainder
    """

    coeffs = poly_as_x_coefficients(p)

    if not coeffs:
        return {}, {}

    max_deg = max(coeffs.keys())

    # Mutable coefficient dictionary.
    work: Dict[int, Dict[int, int]] = {
        i: dict(v)
        for i, v in coeffs.items()
    }

    quotient: Dict[int, Dict[int, int]] = {}

    for deg in range(max_deg, 1, -1):

        lead = work.get(deg, {})

        if not lead:
            continue

        q = dict(lead)

        quotient[deg - 2] = kpoly_add(
            quotient.get(deg - 2, {}),
            q,
        )

        # Subtract q * x^2.
        work[deg] = kpoly_add(
            work.get(deg, {}),
            kpoly_scale(q, -1),
        )

        # Subtract q * 2*d0*x.
        work[deg - 1] = kpoly_add(
            work.get(deg - 1, {}),
            kpoly_scale(q, -2 * d0),
        )

        # Subtract q * (4k-r).
        #
        # 4k-r means coefficients:
        #   k^1 -> 4
        #   k^0 -> -r
        constant_part: Dict[int, int] = {}

        for j, c in q.items():
            constant_part[j + 1] = (
                constant_part.get(j + 1, 0) + 4 * c
            )
            constant_part[j] = (
                constant_part.get(j, 0) - r * c
            )

        work[deg - 2] = kpoly_add(
            work.get(deg - 2, {}),
            kpoly_scale(constant_part, -1),
        )

        # Clean zero entries.
        for d in list(work.keys()):
            if not work[d]:
                del work[d]

    remainder = {
        i: coeffs
        for i, coeffs in work.items()
        if i <= 1 and coeffs
    }

    quotient_poly = x_coefficients_to_poly(quotient)
    remainder_poly = x_coefficients_to_poly(remainder)

    return quotient_poly, remainder_poly


def is_exact_multiple_of_f(
    p: Poly,
    d0: int,
    r: int,
) -> Tuple[bool, Poly, Poly]:
    quotient, remainder = exact_divide_by_f(
        p,
        d0,
        r,
    )

    return (
        len(remainder) == 0,
        quotient,
        remainder,
    )


# -----------------------------------------------------------------------------
# Number theoretic instances
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


def make_instance(
    p: int,
    q: int,
) -> Instance:

    N = p * q
    S = p + q
    d = N + 1 - 2 * S
    gap = q - p

    K = (
        -d * d
        - 3 * N * N
        + 6 * N
        + 1
    ) // 4

    assert (
        -4 * K
        - 3 * N * N
        + 6 * N
        + 1
        == d * d
    )

    assert (
        N * N
        - N
        + K
        == S * (N + 1 - S)
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


# -----------------------------------------------------------------------------
# Partial K
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


def partial_k_data(
    inst: Instance,
    u: int,
) -> PartialK:

    modulus = 1 << u

    K0 = (inst.K // modulus) * modulus
    k = inst.K - K0

    assert 0 <= k < modulus
    assert K0 + k == inst.K

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

    d0 = math.isqrt(C0)
    r = C0 - d0 * d0

    x_true = inst.d - d0

    residual = (
        x_true * x_true
        + 2 * d0 * x_true
        + 4 * k
        - r
    )

    assert residual == 0

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


def exact_d_interval(
    pk: PartialK,
) -> IntervalData:

    lo_sq = pk.C0 - 4 * ((1 << pk.u) - 1)
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

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        x_low=x_low,
        x_high=x_high,
        candidate_count=max(
            0,
            d_high - d_low + 1,
        ),
    )


# -----------------------------------------------------------------------------
# Base polynomial
# -----------------------------------------------------------------------------

def build_base_polynomial(
    pk: PartialK,
) -> Poly:

    return {
        (2, 0): 1,
        (1, 0): 2 * pk.d0,
        (0, 1): 4,
        (0, 0): -pk.r,
    }


# -----------------------------------------------------------------------------
# Shifted polynomial family
# -----------------------------------------------------------------------------

@dataclass
class ShiftedPolynomial:
    i: int
    j: int
    poly: Poly


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
                    poly=poly_mul_monomial(
                        base,
                        i,
                        j,
                    ),
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
# Scaled coefficient lattice
# -----------------------------------------------------------------------------

def poly_to_vector(
    poly: Poly,
    monomials: List[Monomial],
    X: int,
    Y: int,
) -> List[int]:

    return [
        poly.get((i, j), 0)
        * (X ** i)
        * (Y ** j)
        for i, j in monomials
    ]


def vector_to_poly(
    vector: List[int],
    monomials: List[Monomial],
    X: int,
    Y: int,
) -> Tuple[Poly, bool]:

    out: Poly = {}
    integral = True

    for value, (i, j) in zip(
        vector,
        monomials,
    ):

        scale = (
            (X ** i)
            * (Y ** j)
        )

        if scale == 0:
            integral = False
            continue

        if value % scale != 0:
            integral = False
            continue

        coeff = value // scale

        if coeff:
            out[(i, j)] = coeff

    return out, integral


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
    X: int,
    Y: int,
    max_i: int,
    max_j: int,
) -> LatticeData:

    X = max(1, X)
    Y = max(1, Y)

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
# LLL
# -----------------------------------------------------------------------------

try:
    from fpylll import IntegerMatrix, LLL

    HAVE_FPYLLL = True

except Exception:
    HAVE_FPYLLL = False


def run_lll(
    rows: List[List[int]],
) -> Optional[List[List[int]]]:

    if not HAVE_FPYLLL:
        return None

    if not rows:
        return []

    nr = len(rows)
    nc = len(rows[0])

    M = IntegerMatrix(
        nr,
        nc,
    )

    for i in range(nr):
        for j in range(nc):
            M[i, j] = rows[i][j]

    LLL.reduction(M)

    return [
        [
            int(M[i, j])
            for j in range(nc)
        ]
        for i in range(nr)
    ]


# -----------------------------------------------------------------------------
# Vector statistics
# -----------------------------------------------------------------------------

def inf_norm(
    v: List[int],
) -> int:
    return max(
        (abs(x) for x in v),
        default=0,
    )


def nonzero_count(
    v: List[int],
) -> int:
    return sum(
        1
        for x in v
        if x != 0
    )


# -----------------------------------------------------------------------------
# Exact candidate enumeration
# -----------------------------------------------------------------------------

def enumerate_candidates(
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> Tuple[List[int], bool]:

    if interval.candidate_count > limit:
        return [], False

    candidates: List[int] = []

    for d in range(
        interval.d_low,
        interval.d_high + 1,
    ):

        delta = pk.C0 - d * d

        if delta < 0:
            continue

        if delta % 4:
            continue

        k = delta // 4

        if not (
            0
            <= k
            < (1 << pk.u)
        ):
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

    return candidates, True


# -----------------------------------------------------------------------------
# Reduced relation classification
# -----------------------------------------------------------------------------

@dataclass
class RelationAudit:
    index: int
    inf_norm: int
    terms: int
    integral: bool

    value_at_root: Optional[int]

    inherited_multiple: bool
    quotient_terms: int

    remainder_terms: int

    classification: str


def classify_reduced_row(
    index: int,
    row: List[int],
    lattice: LatticeData,
    pk: PartialK,
) -> RelationAudit:

    norm = inf_norm(row)
    terms = nonzero_count(row)

    poly, integral = vector_to_poly(
        row,
        lattice.monomials,
        lattice.X,
        lattice.Y,
    )

    if not integral:
        return RelationAudit(
            index=index,
            inf_norm=norm,
            terms=terms,
            integral=False,
            value_at_root=None,
            inherited_multiple=False,
            quotient_terms=0,
            remainder_terms=0,
            classification="NONINTEGRAL_EMBEDDING",
        )

    value = poly_eval(
        poly,
        pk.x_true,
        pk.k,
    )

    is_multiple, quotient, remainder = (
        is_exact_multiple_of_f(
            poly,
            pk.d0,
            pk.r,
        )
    )

    if is_multiple:
        classification = "INHERITED_MULTIPLE"

    elif value == 0:
        classification = "NEW_ROOT_VANISHING"

    else:
        classification = "NONVANISHING"

    return RelationAudit(
        index=index,
        inf_norm=norm,
        terms=terms,
        integral=True,
        value_at_root=value,
        inherited_multiple=is_multiple,
        quotient_terms=len(quotient),
        remainder_terms=len(remainder),
        classification=classification,
    )


# -----------------------------------------------------------------------------
# Compact per-test record
# -----------------------------------------------------------------------------

@dataclass
class TestResult:
    instance: int
    u: int
    gap: int
    x_true: int
    X: int
    candidate_count: int
    unique: bool

    successful_lll: bool
    reduced_rows: int

    inherited: int
    new_vanishing: int
    nonvanishing: int
    nonintegral: int

    shortest: Optional[int]
    longest: Optional[int]


# -----------------------------------------------------------------------------
# One test
# -----------------------------------------------------------------------------

def run_test(
    instance_index: int,
    inst: Instance,
    u: int,
    max_i: int,
    max_j: int,
    enumeration_limit: int,
) -> TestResult:

    pk = partial_k_data(
        inst,
        u,
    )

    interval = exact_d_interval(pk)

    candidates, enumerated = enumerate_candidates(
        pk,
        interval,
        enumeration_limit,
    )

    unique = (
        enumerated
        and len(candidates) == 1
        and candidates[0] == inst.d
    )

    X = max(
        abs(interval.x_low),
        abs(interval.x_high),
        1,
    )

    Y = max(
        (1 << u) - 1,
        1,
    )

    lattice = build_lattice(
        pk,
        X=X,
        Y=Y,
        max_i=max_i,
        max_j=max_j,
    )

    reduced = run_lll(
        lattice.rows,
    )

    if reduced is None:
        return TestResult(
            instance=instance_index,
            u=u,
            gap=inst.gap,
            x_true=pk.x_true,
            X=X,
            candidate_count=interval.candidate_count,
            unique=unique,
            successful_lll=False,
            reduced_rows=0,
            inherited=0,
            new_vanishing=0,
            nonvanishing=0,
            nonintegral=0,
            shortest=None,
            longest=None,
        )

    audits = [
        classify_reduced_row(
            index=i,
            row=row,
            lattice=lattice,
            pk=pk,
        )
        for i, row in enumerate(reduced)
    ]

    inherited = sum(
        a.classification == "INHERITED_MULTIPLE"
        for a in audits
    )

    new_vanishing = sum(
        a.classification == "NEW_ROOT_VANISHING"
        for a in audits
    )

    nonvanishing = sum(
        a.classification == "NONVANISHING"
        for a in audits
    )

    nonintegral = sum(
        a.classification == "NONINTEGRAL_EMBEDDING"
        for a in audits
    )

    norms = [
        a.inf_norm
        for a in audits
    ]

    return TestResult(
        instance=instance_index,
        u=u,
        gap=inst.gap,
        x_true=pk.x_true,
        X=X,
        candidate_count=interval.candidate_count,
        unique=unique,
        successful_lll=True,
        reduced_rows=len(reduced),
        inherited=inherited,
        new_vanishing=new_vanishing,
        nonvanishing=nonvanishing,
        nonintegral=nonintegral,
        shortest=min(norms) if norms else None,
        longest=max(norms) if norms else None,
    )


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 118)
    print("EXPERIMENT 419 START")
    print("=" * 118)
    print()
    print("EXACT LLL RELATION-NOVELTY AUDIT")
    print()
    print("CORE RELATION")
    print("  f(x,k) = x^2 + 2*d0*x + 4*k-r = 0")
    print()
    print("QUESTION")
    print("  Does LLL produce genuinely new root-vanishing polynomials,")
    print("  or only exact multiples of the known polynomial f?")
    print()
    print("CLASSIFICATION")
    print("  INHERITED_MULTIPLE  = P is exactly divisible by f")
    print("  NEW_ROOT_VANISHING  = P is not divisible by f, but P(root)=0")
    print("  NONVANISHING        = P(root)!=0")
    print("  NONINTEGRAL         = scaled vector is not an integral polynomial")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no symbolic multivariate factorization")
    print("  LLL = diagnostic only")
    print()

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
    enumeration_limit = 200000

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

    print("CONFIGURATION")
    print(f"  instances       = {len(instances)}")
    print(f"  K bits          = {unknown_bits}")
    print(f"  x shifts        = 0..{max_i}")
    print(f"  k shifts        = 0..{max_j}")
    print("  lattice rows    = 9")
    print(f"  enumeration max = {enumeration_limit}")
    print(f"  fpylll          = {HAVE_FPYLLL}")

    if not HAVE_FPYLLL:
        print()
        print("ERROR")
        print("  fpylll could not be imported.")
        print("  Install with:")
        print("    pip install fpylll cysignals")
        return

    # -------------------------------------------------------------------------
    # Generate instances
    # -------------------------------------------------------------------------

    generated: List[Instance] = [
        make_instance(p, q)
        for p, q in instances
    ]

    # -------------------------------------------------------------------------
    # Compact summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 118)
    print("COMPACT RELATION-NOVELTY SUMMARY")
    print("=" * 118)

    print(
        " i   u       gap          x        X      cand unique   "
        "inherit new   nonzero  shortest"
    )
    print("-" * 118)

    results: List[TestResult] = []

    for idx, inst in enumerate(
        generated,
        start=1,
    ):

        for u in unknown_bits:

            result = run_test(
                instance_index=idx,
                inst=inst,
                u=u,
                max_i=max_i,
                max_j=max_j,
                enumeration_limit=enumeration_limit,
            )

            results.append(result)

            print(
                f"{idx:2d} "
                f"{u:3d} "
                f"{inst.gap:10d} "
                f"{result.x_true:10d} "
                f"{result.X:8d} "
                f"{result.candidate_count:9d} "
                f"{'YES' if result.unique else 'NO ':>6} "
                f"{result.inherited:7d} "
                f"{result.new_vanishing:4d} "
                f"{result.nonvanishing:8d} "
                f"{result.shortest if result.shortest is not None else 0}"
            )

    # -------------------------------------------------------------------------
    # Interesting cases
    # -------------------------------------------------------------------------

    interesting = [
        r
        for r in results
        if (
            r.new_vanishing > 0
            or r.nonvanishing > 0
            or r.nonintegral > 0
            or (
                r.candidate_count > 1
                and not r.unique
            )
        )
    ]

    print()
    print("=" * 118)
    print("INTERESTING CASES")
    print("=" * 118)

    if not interesting:
        print("  none")
    else:

        for r in interesting:

            print()
            print(
                f"INSTANCE {r.instance}, "
                f"UNKNOWN K BITS = {r.u}"
            )
            print(
                f"  gap              = {r.gap}"
            )
            print(
                f"  true x           = {r.x_true}"
            )
            print(
                f"  X                = {r.X}"
            )
            print(
                f"  candidate count  = {r.candidate_count}"
            )
            print(
                f"  unique recovery  = {r.unique}"
            )
            print()
            print(
                f"  inherited        = {r.inherited}"
            )
            print(
                f"  new vanishing    = {r.new_vanishing}"
            )
            print(
                f"  nonvanishing     = {r.nonvanishing}"
            )
            print(
                f"  nonintegral      = {r.nonintegral}"
            )
            print(
                f"  shortest         = {r.shortest}"
            )
            print(
                f"  longest          = {r.longest}"
            )

            if r.new_vanishing > 0:
                print()
                print(
                    "  *** GENUINELY NEW ROOT-VANISHING RELATION "
                    "DETECTED ***"
                )

    # -------------------------------------------------------------------------
    # Aggregate classification
    # -------------------------------------------------------------------------

    total = len(results)

    total_inherited = sum(
        r.inherited
        for r in results
    )

    total_new = sum(
        r.new_vanishing
        for r in results
    )

    total_nonzero = sum(
        r.nonvanishing
        for r in results
    )

    total_nonintegral = sum(
        r.nonintegral
        for r in results
    )

    successful = sum(
        r.successful_lll
        for r in results
    )

    unique_count = sum(
        r.unique
        for r in results
    )

    print()
    print("=" * 118)
    print("GLOBAL RELATION CLASSIFICATION")
    print("=" * 118)

    print(f"  total tests                   = {total}")
    print(f"  successful LLL reductions     = {successful}/{total}")
    print(f"  total reduced vectors        = {successful * 9}")
    print(f"  inherited multiples           = {total_inherited}")
    print(f"  genuinely new vanishing      = {total_new}")
    print(f"  nonvanishing                  = {total_nonzero}")
    print(f"  nonintegral                  = {total_nonintegral}")
    print(f"  unique interval recoveries   = {unique_count}/{total}")

    # -------------------------------------------------------------------------
    # Aggregate by K-bit count
    # -------------------------------------------------------------------------

    print()
    print("=" * 118)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 118)

    print(
        " u   tests  unique   multi   inherited   new   nonzero   shortest-min"
    )
    print("-" * 118)

    for u in unknown_bits:

        subset = [
            r
            for r in results
            if r.u == u
        ]

        unique_u = sum(
            r.unique
            for r in subset
        )

        multi_u = len(subset) - unique_u

        inherited_u = sum(
            r.inherited
            for r in subset
        )

        new_u = sum(
            r.new_vanishing
            for r in subset
        )

        nonzero_u = sum(
            r.nonvanishing
            for r in subset
        )

        shortest_values = [
            r.shortest
            for r in subset
            if r.shortest is not None
        ]

        print(
            f"{u:2d} "
            f"{len(subset):5d} "
            f"{unique_u:7d} "
            f"{multi_u:7d} "
            f"{inherited_u:10d} "
            f"{new_u:5d} "
            f"{nonzero_u:9d} "
            f"{min(shortest_values) if shortest_values else 0}"
        )

    # -------------------------------------------------------------------------
    # Final interpretation
    # -------------------------------------------------------------------------

    print()
    print("=" * 118)
    print("INTERPRETATION")
    print("=" * 118)

    print("  The lattice is generated from monomial multiples of f.")
    print("  Therefore inherited root-vanishing relations are expected.")
    print()
    print("  The decisive statistic is NEW_ROOT_VANISHING.")
    print()
    print(
        "  If new vanishing = 0, LLL did not produce a relation outside"
    )
    print(
        "  the exact ideal <f> in this construction."
    )
    print()
    print(
        "  If new vanishing > 0, this is a genuinely new algebraic"
    )
    print(
        "  relation at the hidden root and deserves further study."
    )
    print()
    print(
        "  Even a new relation does NOT by itself prove factor recovery."
    )
    print(
        "  A recovery theorem would still require root isolation or"
    )
    print(
        "  another exact method connecting the relation to p and q."
    )

    print()
    print("=" * 118)
    print("EXPERIMENT 419 FINAL STATUS")
    print("=" * 118)

    if total_new == 0:
        print(
            "  NO GENUINELY NEW ROOT-VANISHING RELATIONS DETECTED = TRUE"
        )
    else:
        print(
            "  GENUINELY NEW ROOT-VANISHING RELATIONS DETECTED = TRUE"
        )
        print(
            f"  count = {total_new}"
        )

    print(
        f"  ALL LLL REDUCTIONS SUCCESSFUL = "
        f"{successful == total}"
    )

    print(
        "  ALL EXACT INTERNAL CHECKS = TRUE"
    )

    print("=" * 118)
    print("EXPERIMENT 419 FINISHED")
    print("=" * 118)


if __name__ == "__main__":
    run_experiment()
