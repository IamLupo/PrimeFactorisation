from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# =============================================================================
# EXPERIMENT 421
# THREE-VARIABLE (x,k,g) / DUAL-EQUATION LLL AUDIT
# =============================================================================
#
# Previous experiments:
#
#   419:
#       ideal-only lattice
#       -> 576/576 vectors were inherited multiples of f.
#
#   420:
#       added arbitrary free monomials
#       -> many "new" root-vanishing relations appeared,
#          but many were trivial low-degree point relations.
#
# Experiment 421 introduces a genuine independent algebraic relation.
#
# Variables:
#
#   x = d-d0
#   k = K-K0
#   g = q-p
#
# Existing equation:
#
#   f(x,k)
#       = x^2 + 2*d0*x + 4*k-r
#       = 0
#
# Independent gap equation:
#
#   4*g^2
#       = (N+1-d)^2 - 16N
#
# with d=d0+x.
#
# Therefore:
#
#   h(x,g)
#       = 4*g^2
#         - (N+1-d0-x)^2
#         + 16N
#       = 0.
#
# The two equations are generated independently.
#
# The lattice contains:
#
#   x^i k^j g^l f
#   x^i k^j g^l h
#
# for configured shift ranges.
#
# The main diagnostic is NOT merely whether vectors vanish at the hidden
# root. That is expected because the basis generators vanish.
#
# Instead we ask whether LLL produces a relation:
#
#   P(x,k,g)
#
# that
#
#   1. vanishes exactly at the hidden root,
#   2. is not an obvious multiple of f,
#   3. is not an obvious multiple of h,
#   4. contains no g variable.
#
# A g-free relation is especially interesting because it would be a relation
# in (x,k) obtained from the joint system without explicitly constructing
# a resultant.
#
# This experiment remains strictly a computational diagnostic.
#
# NO:
#   resultants
#   Groebner-basis computation
#   symbolic multivariate factorization
#   floating-point arithmetic in identities
#   claim of a Coppersmith theorem
#
# Exact integer arithmetic only.
# =============================================================================


# -----------------------------------------------------------------------------
# Polynomial infrastructure
# -----------------------------------------------------------------------------

# (x-degree, k-degree, g-degree)
Monomial = Tuple[int, int, int]
Poly = Dict[Monomial, int]


def poly_clean(p: Poly) -> Poly:
    return {
        m: c
        for m, c in p.items()
        if c != 0
    }


def poly_add(a: Poly, b: Poly) -> Poly:
    out = dict(a)

    for m, c in b.items():
        out[m] = out.get(m, 0) + c

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

    for (i1, j1, l1), c1 in a.items():
        for (i2, j2, l2), c2 in b.items():

            m = (
                i1 + i2,
                j1 + j2,
                l1 + l2,
            )

            out[m] = (
                out.get(m, 0)
                + c1 * c2
            )

    return poly_clean(out)


def poly_mul_monomial(
    p: Poly,
    dx: int,
    dk: int,
    dg: int,
) -> Poly:

    return {
        (
            i + dx,
            j + dk,
            l + dg,
        ): c
        for (i, j, l), c in p.items()
    }


def poly_eval(
    p: Poly,
    x: int,
    k: int,
    g: int,
) -> int:

    total = 0

    for (i, j, l), c in p.items():
        total += (
            c
            * (x ** i)
            * (k ** j)
            * (g ** l)
        )

    return total


def poly_total_degree(p: Poly) -> int:
    return max(
        (
            i + j + l
            for i, j, l in p
        ),
        default=0,
    )


def poly_is_g_free(p: Poly) -> bool:
    return all(
        l == 0
        for _, _, l in p
    )


# -----------------------------------------------------------------------------
# Number-theoretic instance
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

    d = (
        N
        + 1
        - 2 * S
    )

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
        ==
        S * (N + 1 - S)
    )

    # Independent gap identity.
    assert (
        4 * gap * gap
        ==
        (N + 1 - d) ** 2
        - 16 * N
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

    K0 = (
        inst.K // modulus
    ) * modulus

    k = (
        inst.K
        - K0
    )

    assert (
        0 <= k < modulus
    )

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

    d0 = math.isqrt(
        C0
    )

    r = (
        C0
        - d0 * d0
    )

    x_true = (
        inst.d
        - d0
    )

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
# Exact equations
# -----------------------------------------------------------------------------

def build_f(
    pk: PartialK,
) -> Poly:

    # x^2 + 2*d0*x + 4*k-r
    return {
        (2, 0, 0): 1,
        (1, 0, 0): 2 * pk.d0,
        (0, 1, 0): 4,
        (0, 0, 0): -pk.r,
    }


def build_h(
    inst: Instance,
    pk: PartialK,
) -> Poly:

    # h =
    #
    #   4*g^2
    #   - (N+1-d0-x)^2
    #   + 16N
    #
    #
    # Let A = N+1-d0.
    #
    # Then
    #
    #   h = 4g^2 - (A-x)^2 + 16N
    #     = 4g^2 - x^2 + 2A*x - A^2 + 16N.

    A = (
        inst.N
        + 1
        - pk.d0
    )

    return {
        (0, 0, 2): 4,
        (2, 0, 0): -1,
        (1, 0, 0): 2 * A,
        (0, 0, 0): (
            -A * A
            + 16 * inst.N
        ),
    }


# -----------------------------------------------------------------------------
# Exact interval for x
# -----------------------------------------------------------------------------

@dataclass
class IntervalData:
    x_low: int
    x_high: int
    candidate_count: int


def exact_x_interval(
    pk: PartialK,
) -> IntervalData:

    lo_sq = (
        pk.C0
        - 4 * ((1 << pk.u) - 1)
    )

    hi_sq = pk.C0

    if lo_sq <= 0:
        d_low = 0

    else:
        d_low = math.isqrt(
            lo_sq
        )

        while (
            d_low * d_low
            < lo_sq
        ):
            d_low += 1

    d_high = math.isqrt(
        hi_sq
    )

    x_low = (
        d_low
        - pk.d0
    )

    x_high = (
        d_high
        - pk.d0
    )

    return IntervalData(
        x_low=x_low,
        x_high=x_high,
        candidate_count=max(
            0,
            d_high - d_low + 1,
        ),
    )


# -----------------------------------------------------------------------------
# Lattice representation
# -----------------------------------------------------------------------------

@dataclass
class Lattice:
    monomials: List[Monomial]
    rows: List[List[int]]
    labels: List[str]
    X: int
    Y: int
    G: int


def monomial_basis(
    polys: List[Poly],
) -> List[Monomial]:

    mons = set()

    for p in polys:
        mons.update(
            p.keys()
        )

    return sorted(mons)


def poly_to_scaled_vector(
    p: Poly,
    monomials: List[Monomial],
    X: int,
    Y: int,
    G: int,
) -> List[int]:

    out = []

    for i, j, l in monomials:

        out.append(
            p.get(
                (i, j, l),
                0,
            )
            * (X ** i)
            * (Y ** j)
            * (G ** l)
        )

    return out


def vector_to_poly(
    row: List[int],
    monomials: List[Monomial],
    X: int,
    Y: int,
    G: int,
) -> Tuple[Poly, bool]:

    p: Poly = {}
    integral = True

    for value, (i, j, l) in zip(
        row,
        monomials,
    ):

        scale = (
            (X ** i)
            * (Y ** j)
            * (G ** l)
        )

        if scale == 0:
            integral = False
            continue

        if value % scale != 0:
            integral = False
            continue

        c = (
            value // scale
        )

        if c:
            p[(i, j, l)] = c

    return p, integral


# -----------------------------------------------------------------------------
# Build dual-equation lattice
# -----------------------------------------------------------------------------

def build_lattice(
    inst: Instance,
    pk: PartialK,
    X: int,
    Y: int,
    G: int,
    f_x: int,
    f_k: int,
    f_g: int,
    h_x: int,
    h_k: int,
    h_g: int,
) -> Lattice:

    f = build_f(pk)
    h = build_h(inst, pk)

    polys: List[Poly] = []
    labels: List[str] = []

    # Multiples of f.
    for i in range(f_x + 1):
        for j in range(f_k + 1):
            for l in range(f_g + 1):

                polys.append(
                    poly_mul_monomial(
                        f,
                        i,
                        j,
                        l,
                    )
                )

                labels.append(
                    "F"
                )

    # Multiples of h.
    for i in range(h_x + 1):
        for j in range(h_k + 1):
            for l in range(h_g + 1):

                polys.append(
                    poly_mul_monomial(
                        h,
                        i,
                        j,
                        l,
                    )
                )

                labels.append(
                    "H"
                )

    mons = monomial_basis(
        polys
    )

    rows = [
        poly_to_scaled_vector(
            p,
            mons,
            X,
            Y,
            G,
        )
        for p in polys
    ]

    return Lattice(
        monomials=mons,
        rows=rows,
        labels=labels,
        X=X,
        Y=Y,
        G=G,
    )


# -----------------------------------------------------------------------------
# LLL
# -----------------------------------------------------------------------------

try:
    from fpylll import (
        IntegerMatrix,
        LLL,
    )

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
# Relation diagnostics
# -----------------------------------------------------------------------------

@dataclass
class RelationAudit:
    index: int
    norm: int
    terms: int
    degree: int
    integral: bool
    value: Optional[int]
    g_free: bool
    classification: str


def classify_relation(
    p: Poly,
    value: int,
    integral: bool,
) -> str:

    if not integral:
        return "NONINTEGRAL"

    if value != 0:
        return "NONVANISHING"

    # Root-vanishing.
    #
    # We do NOT attempt symbolic ideal membership here.
    #
    # Instead:
    #
    #   g-containing + vanishing
    #       -> JOINT_VANISHING
    #
    #   g-free + vanishing
    #       -> GAP_ELIMINATED
    #
    # A GAP_ELIMINATED relation is the key result.
    #

    if poly_is_g_free(p):
        return "GAP_ELIMINATED"

    return "JOINT_VANISHING"


def audit_reduced_rows(
    reduced: List[List[int]],
    lattice: Lattice,
    pk: PartialK,
    inst: Instance,
) -> List[RelationAudit]:

    audits = []

    for index, row in enumerate(
        reduced
    ):

        norm = max(
            (
                abs(v)
                for v in row
            ),
            default=0,
        )

        terms = sum(
            v != 0
            for v in row
        )

        p, integral = vector_to_poly(
            row,
            lattice.monomials,
            lattice.X,
            lattice.Y,
            lattice.G,
        )

        if integral:

            value = poly_eval(
                p,
                pk.x_true,
                pk.k,
                inst.gap,
            )

        else:
            value = None

        audits.append(
            RelationAudit(
                index=index,
                norm=norm,
                terms=terms,
                degree=poly_total_degree(p),
                integral=integral,
                value=value,
                g_free=(
                    integral
                    and poly_is_g_free(p)
                ),
                classification=classify_relation(
                    p,
                    value if value is not None else 1,
                    integral,
                ),
            )
        )

    return audits


# -----------------------------------------------------------------------------
# Exact candidate filtering
# -----------------------------------------------------------------------------

def enumerate_candidates(
    inst: Instance,
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> Tuple[List[int], bool]:

    if interval.candidate_count > limit:
        return [], False

    candidates = []

    for x in range(
        interval.x_low,
        interval.x_high + 1,
    ):

        delta = (
            pk.r
            - x * x
            - 2 * pk.d0 * x
        )

        if delta % 4 != 0:
            continue

        k = delta // 4

        if not (
            0 <= k < (1 << pk.u)
        ):
            continue

        # Independent gap equation.
        h_value = (
            4 * inst.gap * inst.gap
            - (
                inst.N
                + 1
                - pk.d0
                - x
            ) ** 2
            + 16 * inst.N
        )

        if h_value != 0:
            continue

        d = (
            pk.d0
            + x
        )

        candidates.append(d)

    return candidates, True


# -----------------------------------------------------------------------------
# Result
# -----------------------------------------------------------------------------

@dataclass
class Result:
    instance: int
    u: int
    gap: int
    x_true: int

    X: int
    candidates: int
    unique: bool

    joint_vanishing: int
    gap_eliminated: int
    nonzero: int
    nonintegral: int

    shortest: int
    longest: int


# -----------------------------------------------------------------------------
# Single test
# -----------------------------------------------------------------------------

def run_test(
    idx: int,
    inst: Instance,
    u: int,
    f_x: int,
    f_k: int,
    f_g: int,
    h_x: int,
    h_k: int,
    h_g: int,
    enumeration_limit: int,
) -> Result:

    pk = partial_k_data(
        inst,
        u,
    )

    interval = exact_x_interval(
        pk
    )

    candidates, enumerated = (
        enumerate_candidates(
            inst,
            pk,
            interval,
            enumeration_limit,
        )
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

    # Conservative gap bound:
    #
    # g = q-p < sqrt(N)
    #
    # use isqrt(N)+1.
    G = max(
        math.isqrt(inst.N) + 1,
        1,
    )

    lattice = build_lattice(
        inst,
        pk,
        X,
        Y,
        G,
        f_x,
        f_k,
        f_g,
        h_x,
        h_k,
        h_g,
    )

    reduced = run_lll(
        lattice.rows
    )

    if reduced is None:
        raise RuntimeError(
            "fpylll unavailable"
        )

    audits = audit_reduced_rows(
        reduced,
        lattice,
        pk,
        inst,
    )

    joint = sum(
        a.classification
        == "JOINT_VANISHING"
        for a in audits
    )

    gap_eliminated = sum(
        a.classification
        == "GAP_ELIMINATED"
        for a in audits
    )

    nonzero = sum(
        a.classification
        == "NONVANISHING"
        for a in audits
    )

    nonintegral = sum(
        a.classification
        == "NONINTEGRAL"
        for a in audits
    )

    norms = [
        a.norm
        for a in audits
    ]

    return Result(
        instance=idx,
        u=u,
        gap=inst.gap,
        x_true=pk.x_true,
        X=X,
        candidates=interval.candidate_count,
        unique=unique,
        joint_vanishing=joint,
        gap_eliminated=gap_eliminated,
        nonzero=nonzero,
        nonintegral=nonintegral,
        shortest=min(norms),
        longest=max(norms),
    )


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 421 START")
    print("=" * 120)
    print()
    print("THREE-VARIABLE (x,k,g) / DUAL-EQUATION LLL AUDIT")
    print()
    print("CORE EQUATIONS")
    print("  f(x,k) = x^2 + 2*d0*x + 4*k-r")
    print(
        "  h(x,g) = 4*g^2-(N+1-d0-x)^2+16N"
    )
    print()
    print("NEW VARIABLE")
    print("  g = q-p")
    print()
    print("QUESTION")
    print(
        "  Can LLL find a new g-free relation from the joint system?"
    )
    print()
    print("KEY DIAGNOSTIC")
    print(
        "  GAP_ELIMINATED = exact root-vanishing polynomial with no g."
    )
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
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

    # Moderate shift ranges.
    f_x = 1
    f_k = 1
    f_g = 1

    h_x = 1
    h_k = 1
    h_g = 1

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
    print(
        f"  instances       = {len(instances)}"
    )
    print(
        f"  K bits          = {unknown_bits}"
    )
    print(
        f"  f shifts        = ({f_x},{f_k},{f_g})"
    )
    print(
        f"  h shifts        = ({h_x},{h_k},{h_g})"
    )
    print(
        f"  enumeration max = {enumeration_limit}"
    )
    print(
        f"  fpylll          = {HAVE_FPYLLL}"
    )

    if not HAVE_FPYLLL:
        print()
        print("ERROR")
        print(
            "  fpylll import failed."
        )
        print(
            "  Install with:"
        )
        print(
            "    pip install fpylll cysignals"
        )
        return

    generated = [
        make_instance(
            p,
            q,
        )
        for p, q in instances
    ]

    # -------------------------------------------------------------------------
    # Compact summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT DUAL-EQUATION SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap          x        X      cand unique "
        "joint g-free nonzero shortest"
    )

    print("-" * 120)

    results: List[Result] = []

    for idx, inst in enumerate(
        generated,
        start=1,
    ):

        for u in unknown_bits:

            result = run_test(
                idx,
                inst,
                u,
                f_x,
                f_k,
                f_g,
                h_x,
                h_k,
                h_g,
                enumeration_limit,
            )

            results.append(
                result
            )

            print(
                f"{idx:2d} "
                f"{u:3d} "
                f"{inst.gap:10d} "
                f"{result.x_true:10d} "
                f"{result.X:8d} "
                f"{result.candidates:9d} "
                f"{'YES' if result.unique else 'NO ':>6} "
                f"{result.joint_vanishing:4d} "
                f"{result.gap_eliminated:4d} "
                f"{result.nonzero:7d} "
                f"{result.shortest}"
            )

    # -------------------------------------------------------------------------
    # Global counts
    # -------------------------------------------------------------------------

    total = len(results)

    joint_total = sum(
        r.joint_vanishing
        for r in results
    )

    eliminated_total = sum(
        r.gap_eliminated
        for r in results
    )

    nonzero_total = sum(
        r.nonzero
        for r in results
    )

    nonintegral_total = sum(
        r.nonintegral
        for r in results
    )

    unique_total = sum(
        r.unique
        for r in results
    )

    print()
    print("=" * 120)
    print("GLOBAL DUAL-EQUATION CLASSIFICATION")
    print("=" * 120)

    print(
        f"  total tests                = {total}"
    )
    print(
        f"  reduced vectors            = {total * 16}"
    )
    print(
        f"  joint root-vanishing       = {joint_total}"
    )
    print(
        f"  g-free root relations      = {eliminated_total}"
    )
    print(
        f"  nonvanishing               = {nonzero_total}"
    )
    print(
        f"  nonintegral                = {nonintegral_total}"
    )
    print(
        f"  unique interval recoveries = {unique_total}/{total}"
    )

    # -------------------------------------------------------------------------
    # Aggregate
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests unique   joint   g-free   nonzero   shortest-min"
    )

    print("-" * 120)

    for u in unknown_bits:

        subset = [
            r
            for r in results
            if r.u == u
        ]

        print(
            f"{u:2d} "
            f"{len(subset):5d} "
            f"{sum(r.unique for r in subset):6d} "
            f"{sum(r.joint_vanishing for r in subset):7d} "
            f"{sum(r.gap_eliminated for r in subset):8d} "
            f"{sum(r.nonzero for r in subset):8d} "
            f"{min(r.shortest for r in subset)}"
        )

    # -------------------------------------------------------------------------
    # Interesting cases
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("INTERESTING CASES")
    print("=" * 120)

    interesting = [
        r
        for r in results
        if r.gap_eliminated > 0
    ]

    if not interesting:
        print(
            "  No g-free root-vanishing relation detected."
        )

    else:

        for r in interesting:

            print(
                f"  INSTANCE {r.instance}, u={r.u}: "
                f"g-free relations={r.gap_eliminated}, "
                f"x={r.x_true}, X={r.X}, "
                f"candidates={r.candidates}"
            )

    # -------------------------------------------------------------------------
    # Interpretation
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 420 used arbitrary free monomials."
    )

    print(
        "  Experiment 421 instead introduces the independent exact"
    )

    print(
        "  relation coming from the prime gap g=q-p."
    )

    print()
    print(
        "  A JOINT_VANISHING relation still contains g and is therefore"
    )

    print(
        "  not an eliminated relation."
    )

    print()
    print(
        "  GAP_ELIMINATED is the important statistic."
    )

    print(
        "  Such a polynomial depends only on (x,k) while being obtained"
    )

    print(
        "  from the joint (x,k,g) system."
    )

    print()
    print(
        "  This does NOT construct a resultant and does NOT establish"
    )

    print(
        "  a factorization theorem."
    )

    print()
    print(
        "  If g-free relations appear, the next experiment should inspect"
    )

    print(
        "  their exact coefficients and determine whether they reduce to"
    )

    print(
        "  the existing f relation or provide genuinely new information."
    )

    # -------------------------------------------------------------------------
    # Final
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXPERIMENT 421 FINAL STATUS")
    print("=" * 120)

    print(
        "  ALL EXACT INTERNAL CHECKS = TRUE"
    )

    print(
        f"  G-FREE ROOT RELATIONS = {eliminated_total}"
    )

    if eliminated_total:
        print(
            "  NEW GAP-ELIMINATED RELATION DETECTED = TRUE"
        )
    else:
        print(
            "  NO GAP-ELIMINATED RELATION DETECTED = TRUE"
        )

    print("=" * 120)
    print("EXPERIMENT 421 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
