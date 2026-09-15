from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# =============================================================================
# EXPERIMENT 422
# EXACT G-FREE RELATION NOVELTY / IDEAL-MEMBERSHIP AUDIT
# =============================================================================
#
# Experiment 421 found g-free relations in the joint (x,k,g) lattice.
#
# This experiment asks:
#
#   Are the g-free relations genuinely new modulo <f>?
#
# CORE RELATION
#
#   f(x,k) = x^2 + 2*d0*x + 4*k-r
#
# JOINT RELATION
#
#   h(x,g) = 4*g^2-(N+1-d0-x)^2+16N
#
# CLASSIFICATION
#
#   JOINT_VANISHING
#       exact polynomial vanishes at the hidden (x,k,g) root
#       and contains g
#
#   GFREE_INHERITED
#       g-free, vanishes, and is exactly divisible by f
#
#   GFREE_NEW_K_DEPENDENT
#       g-free, vanishes, not divisible by f, still contains k
#
#   GFREE_NEW_X_ONLY
#       g-free, vanishes, not divisible by f, contains only x
#
#   NONVANISHING
#       exact polynomial does not vanish
#
#   NONINTEGRAL
#       LLL vector cannot be reconstructed as an integral polynomial
#
# RULES
#
#   exact integer arithmetic only
#   no resultants
#   no Groebner basis
#   no symbolic factorization
#   LLL = diagnostic only
#
# IMPORTANT
#
#   The exact divisibility test exploits the fact that f is linear in k:
#
#       f = 4k + A(x)
#
#   and performs exact polynomial division in k.
#
# =============================================================================


# -----------------------------------------------------------------------------
# Polynomial representation
# -----------------------------------------------------------------------------

# (x_degree, k_degree, g_degree)
Monomial = Tuple[int, int, int]
Poly = Dict[Monomial, int]


def clean_poly(p: Poly) -> Poly:
    return {
        m: c
        for m, c in p.items()
        if c != 0
    }


def poly_add(a: Poly, b: Poly) -> Poly:
    out = dict(a)

    for m, c in b.items():
        out[m] = out.get(m, 0) + c

    return clean_poly(out)


def poly_scale(a: Poly, c: int) -> Poly:
    if c == 0:
        return {}

    return {
        m: v * c
        for m, v in a.items()
        if v * c != 0
    }


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


def is_g_free(p: Poly) -> bool:
    return all(
        l == 0
        for _, _, l in p
    )


def is_x_only(p: Poly) -> bool:
    return all(
        j == 0 and l == 0
        for _, j, l in p
    )


def has_k(p: Poly) -> bool:
    return any(
        j != 0
        for _, j, _ in p
    )


# -----------------------------------------------------------------------------
# Primitive / canonical polynomial form
# -----------------------------------------------------------------------------

def coefficient_gcd(p: Poly) -> int:
    g = 0

    for c in p.values():
        g = math.gcd(
            g,
            abs(c),
        )

    return g


def canonicalize_poly(
    p: Poly,
) -> Optional[Tuple[Tuple[Monomial, int], ...]]:

    if not p:
        return None

    g = coefficient_gcd(p)

    if g == 0:
        return None

    normalized = {
        m: c // g
        for m, c in p.items()
    }

    first_monomial = sorted(
        normalized
    )[0]

    if normalized[first_monomial] < 0:
        normalized = {
            m: -c
            for m, c in normalized.items()
        }

    return tuple(
        sorted(
            normalized.items()
        )
    )


def canonical_to_poly(
    key: Optional[
        Tuple[
            Tuple[Monomial, int],
            ...
        ]
    ],
) -> Poly:
    if key is None:
        return {}

    return dict(key)


def format_poly(
    p: Poly,
    max_terms: int = 10,
) -> str:

    if not p:
        return "0"

    parts: List[str] = []

    for (i, j, l), c in sorted(p.items()):

        if len(parts) >= max_terms:
            parts.append("...")
            break

        names: List[str] = []

        if i:
            names.append(
                "x" if i == 1 else f"x^{i}"
            )

        if j:
            names.append(
                "k" if j == 1 else f"k^{j}"
            )

        if l:
            names.append(
                "g" if l == 1 else f"g^{l}"
            )

        monomial = "*".join(names)

        if not monomial:
            parts.append(str(c))
        elif c == 1:
            parts.append(monomial)
        elif c == -1:
            parts.append(f"-{monomial}")
        else:
            parts.append(
                f"{c}*{monomial}"
            )

    return " + ".join(parts)


# -----------------------------------------------------------------------------
# Univariate-in-x polynomial helpers
# -----------------------------------------------------------------------------

XPoly = Dict[int, int]


def xpoly_clean(p: XPoly) -> XPoly:
    return {
        d: c
        for d, c in p.items()
        if c != 0
    }


def xpoly_add(
    a: XPoly,
    b: XPoly,
    scale_b: int = 1,
) -> XPoly:

    out = dict(a)

    for degree, coeff in b.items():
        out[degree] = (
            out.get(degree, 0)
            + scale_b * coeff
        )

        if out[degree] == 0:
            del out[degree]

    return out


def xpoly_scale(
    a: XPoly,
    c: int,
) -> XPoly:
    return {
        degree: coeff * c
        for degree, coeff in a.items()
        if coeff * c != 0
    }


# -----------------------------------------------------------------------------
# Exact divisibility by f
# -----------------------------------------------------------------------------

def f_x_part(
    pk,
) -> XPoly:
    # f(x,k) = 4k + A(x)
    #
    # A(x) = x^2 + 2*d0*x-r
    return {
        2: 1,
        1: 2 * pk.d0,
        0: -pk.r,
    }


def exact_divisible_by_f(
    p: Poly,
    pk,
) -> Tuple[bool, Optional[Poly]]:

    if not p:
        return True, {}

    if not is_g_free(p):
        return False, None

    # Group P by powers of k:
    #
    #   P(x,k) = sum_j A_j(x) k^j
    #
    coeffs: Dict[int, XPoly] = {}

    for (
        degree_x,
        degree_k,
        degree_g,
    ), coeff in p.items():

        if degree_g != 0:
            return False, None

        coeffs.setdefault(
            degree_k,
            {},
        )

        coeffs[degree_k][degree_x] = (
            coeffs[degree_k].get(
                degree_x,
                0,
            )
            + coeff
        )

    coeffs = {
        j: xpoly_clean(v)
        for j, v in coeffs.items()
        if xpoly_clean(v)
    }

    if not coeffs:
        return True, {}

    # f = 4k + A(x)
    A = f_x_part(pk)

    max_degree = max(coeffs)

    # Work on a mutable copy.
    remainder: Dict[int, XPoly] = {
        j: dict(v)
        for j, v in coeffs.items()
    }

    quotient_by_k: Dict[int, XPoly] = {}

    #
    # Polynomial long division in k.
    #
    # Since the leading coefficient of f in k is 4,
    # every coefficient of the current leading X-polynomial
    # must be divisible by 4 for an integral quotient.
    #
    for degree in range(
        max_degree,
        0,
        -1,
    ):

        current = xpoly_clean(
            remainder.get(
                degree,
                {},
            )
        )

        if not current:
            continue

        if any(
            coeff % 4 != 0
            for coeff in current.values()
        ):
            return False, None

        quotient_piece = {
            deg: coeff // 4
            for deg, coeff in current.items()
        }

        quotient_by_k[
            degree - 1
        ] = xpoly_add(
            quotient_by_k.get(
                degree - 1,
                {},
            ),
            quotient_piece,
        )

        # Subtract quotient_piece * A(x) from the k^(degree-1) coefficient.
        correction: XPoly = {}

        for deg_q, coeff_q in quotient_piece.items():
            for deg_a, coeff_a in A.items():

                target_degree = (
                    deg_q
                    + deg_a
                )

                correction[target_degree] = (
                    correction.get(
                        target_degree,
                        0,
                    )
                    + coeff_q * coeff_a
                )

        remainder[
            degree - 1
        ] = xpoly_add(
            remainder.get(
                degree - 1,
                {},
            ),
            correction,
            scale_b=-1,
        )

    # Exact divisibility requires the final k^0 remainder to be zero.
    final_remainder = xpoly_clean(
        remainder.get(
            0,
            {},
        )
    )

    if final_remainder:
        return False, None

    quotient: Poly = {}

    for degree_k, xpoly in quotient_by_k.items():
        for degree_x, coeff in xpoly.items():

            if coeff != 0:
                quotient[
                    (
                        degree_x,
                        degree_k,
                        0,
                    )
                ] = coeff

    return True, clean_poly(
        quotient
    )


# -----------------------------------------------------------------------------
# Number-theoretic instances
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

    # Existing exact identities.
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
        == S * (
            N + 1 - S
        )
    )

    # Gap identity:
    #
    #   (q-p)^2 = (p+q)^2 - 4pq = S^2 - 4N
    #
    # and
    #
    #   S = (N+1-d)/2
    #
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
# Partial K decomposition
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
# Core equations
# -----------------------------------------------------------------------------

def build_f(
    pk: PartialK,
) -> Poly:

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

    A = (
        inst.N
        + 1
        - pk.d0
    )

    return {
        (0, 0, 2): 4,
        (2, 0, 0): -1,
        (1, 0, 0): 2 * A,
        (0, 0, 0):
            -A * A
            + 16 * inst.N,
    }


# -----------------------------------------------------------------------------
# Exact interval
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
        - 4 * (
            (1 << pk.u) - 1
        )
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

    return IntervalData(
        x_low=d_low - pk.d0,
        x_high=d_high - pk.d0,
        candidate_count=max(
            0,
            d_high - d_low + 1,
        ),
    )


# -----------------------------------------------------------------------------
# Lattice
# -----------------------------------------------------------------------------

@dataclass
class Lattice:
    monomials: List[Monomial]
    rows: List[List[int]]
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


def poly_to_vector(
    p: Poly,
    mons: List[Monomial],
    X: int,
    Y: int,
    G: int,
) -> List[int]:

    # FIX:
    # The previous script referenced an undefined variable `m`.
    # Here we explicitly unpack each monomial as (i, j, l).

    vec: List[int] = []

    for i, j, l in mons:

        coefficient = p.get(
            (i, j, l),
            0,
        )

        scale = (
            (X ** i)
            * (Y ** j)
            * (G ** l)
        )

        vec.append(
            coefficient * scale
        )

    return vec


def vector_to_poly(
    row: List[int],
    mons: List[Monomial],
    X: int,
    Y: int,
    G: int,
) -> Tuple[Poly, bool]:

    p: Poly = {}
    integral = True

    for value, (
        i,
        j,
        l,
    ) in zip(
        row,
        mons,
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

        coeff = value // scale

        if coeff != 0:
            p[
                (
                    i,
                    j,
                    l,
                )
            ] = coeff

    return (
        clean_poly(p),
        integral,
    )


def build_lattice(
    inst: Instance,
    pk: PartialK,
    X: int,
    Y: int,
    G: int,
) -> Lattice:

    f = build_f(pk)
    h = build_h(
        inst,
        pk,
    )

    polys: List[Poly] = []

    #
    # Experiment 421 configuration:
    #
    #   f * x^i * k^j * g^l
    #   h * x^i * k^j * g^l
    #
    # with i,j,l in {0,1}
    #

    for i in range(2):
        for j in range(2):
            for l in range(2):

                polys.append(
                    poly_mul_monomial(
                        f,
                        i,
                        j,
                        l,
                    )
                )

    for i in range(2):
        for j in range(2):
            for l in range(2):

                polys.append(
                    poly_mul_monomial(
                        h,
                        i,
                        j,
                        l,
                    )
                )

    mons = monomial_basis(
        polys
    )

    rows = [
        poly_to_vector(
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

    matrix = IntegerMatrix(
        nr,
        nc,
    )

    for i in range(nr):
        for j in range(nc):
            matrix[i, j] = rows[i][j]

    LLL.reduction(
        matrix
    )

    return [
        [
            int(matrix[i, j])
            for j in range(nc)
        ]
        for i in range(nr)
    ]


# -----------------------------------------------------------------------------
# Classification
# -----------------------------------------------------------------------------

@dataclass
class Audit:
    index: int
    norm: int
    terms: int
    degree: int

    integral: bool
    value: Optional[int]

    g_free: bool
    inherited: bool
    x_only: bool

    classification: str

    canonical: Optional[
        Tuple[
            Tuple[Monomial, int],
            ...
        ]
    ]


def classify_polynomial(
    p: Poly,
    integral: bool,
    value: Optional[int],
    pk: PartialK,
) -> Tuple[
    str,
    bool,
    bool,
    bool,
    Optional[
        Tuple[
            Tuple[Monomial, int],
            ...
        ]
    ],
]:

    if not integral:
        return (
            "NONINTEGRAL",
            False,
            False,
            False,
            None,
        )

    if value != 0:
        return (
            "NONVANISHING",
            is_g_free(p),
            False,
            is_x_only(p),
            None,
        )

    # Vanishing and still contains g.
    if not is_g_free(p):
        return (
            "JOINT_VANISHING",
            False,
            False,
            False,
            None,
        )

    # g-free: now test exact membership in <f>.
    divisible, _quotient = (
        exact_divisible_by_f(
            p,
            pk,
        )
    )

    if divisible:
        return (
            "GFREE_INHERITED",
            True,
            True,
            is_x_only(p),
            None,
        )

    canonical = canonicalize_poly(
        p
    )

    if is_x_only(p):
        return (
            "GFREE_NEW_X_ONLY",
            True,
            False,
            True,
            canonical,
        )

    return (
        "GFREE_NEW_K_DEPENDENT",
        True,
        False,
        False,
        canonical,
    )


def audit_rows(
    reduced: List[List[int]],
    lattice: Lattice,
    inst: Instance,
    pk: PartialK,
) -> List[Audit]:

    audits: List[Audit] = []

    for index, row in enumerate(
        reduced
    ):

        norm = max(
            (
                abs(value)
                for value in row
            ),
            default=0,
        )

        terms = sum(
            value != 0
            for value in row
        )

        polynomial, integral = (
            vector_to_poly(
                row,
                lattice.monomials,
                lattice.X,
                lattice.Y,
                lattice.G,
            )
        )

        value: Optional[int] = None

        if integral:
            value = poly_eval(
                polynomial,
                pk.x_true,
                pk.k,
                inst.gap,
            )

        (
            classification,
            g_free,
            inherited,
            x_only,
            canonical,
        ) = classify_polynomial(
            polynomial,
            integral,
            value,
            pk,
        )

        audits.append(
            Audit(
                index=index,
                norm=norm,
                terms=terms,
                degree=poly_total_degree(
                    polynomial
                ),
                integral=integral,
                value=value,
                g_free=g_free,
                inherited=inherited,
                x_only=x_only,
                classification=classification,
                canonical=canonical,
            )
        )

    return audits


# -----------------------------------------------------------------------------
# Exact candidate verification
# -----------------------------------------------------------------------------

def enumerate_candidates(
    inst: Instance,
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> Tuple[List[int], bool]:

    if interval.candidate_count > limit:
        return [], False

    candidates: List[int] = []

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

        # Verify h exactly at the corresponding x and true g.
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

        candidates.append(
            pk.d0 + x
        )

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
    candidate_count: int
    unique: bool

    inherited: int
    new_k: int
    new_x: int
    joint: int
    nonzero: int
    nonintegral: int

    primitive_new_k: int
    primitive_new_x: int

    shortest: int
    longest: int


# -----------------------------------------------------------------------------
# One test
# -----------------------------------------------------------------------------

def run_test(
    instance_index: int,
    inst: Instance,
    u: int,
    enumeration_limit: int,
) -> Tuple[
    Result,
    List[Audit],
]:

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

    # A conservative exact scaling bound for g.
    #
    # The true g satisfies g^2 = S^2 - 4N.
    # Using N^(1/2)+1 is sufficient as a finite lattice scale.
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
    )

    reduced = run_lll(
        lattice.rows
    )

    if reduced is None:
        raise RuntimeError(
            "fpylll is unavailable"
        )

    audits = audit_rows(
        reduced,
        lattice,
        inst,
        pk,
    )

    inherited = sum(
        audit.classification
        == "GFREE_INHERITED"
        for audit in audits
    )

    new_k = sum(
        audit.classification
        == "GFREE_NEW_K_DEPENDENT"
        for audit in audits
    )

    new_x = sum(
        audit.classification
        == "GFREE_NEW_X_ONLY"
        for audit in audits
    )

    joint = sum(
        audit.classification
        == "JOINT_VANISHING"
        for audit in audits
    )

    nonzero = sum(
        audit.classification
        == "NONVANISHING"
        for audit in audits
    )

    nonintegral = sum(
        audit.classification
        == "NONINTEGRAL"
        for audit in audits
    )

    primitive_new_k = len({
        audit.canonical
        for audit in audits
        if (
            audit.classification
            == "GFREE_NEW_K_DEPENDENT"
            and audit.canonical is not None
        )
    })

    primitive_new_x = len({
        audit.canonical
        for audit in audits
        if (
            audit.classification
            == "GFREE_NEW_X_ONLY"
            and audit.canonical is not None
        )
    })

    norms = [
        audit.norm
        for audit in audits
    ]

    result = Result(
        instance=instance_index,
        u=u,
        gap=inst.gap,
        x_true=pk.x_true,
        X=X,
        candidate_count=interval.candidate_count,
        unique=unique,
        inherited=inherited,
        new_k=new_k,
        new_x=new_x,
        joint=joint,
        nonzero=nonzero,
        nonintegral=nonintegral,
        primitive_new_k=primitive_new_k,
        primitive_new_x=primitive_new_x,
        shortest=min(norms),
        longest=max(norms),
    )

    return (
        result,
        audits,
    )


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 422 START")
    print("=" * 120)

    print()
    print(
        "EXACT G-FREE RELATION NOVELTY / IDEAL-MEMBERSHIP AUDIT"
    )

    print()
    print("CORE RELATION")
    print(
        "  f(x,k) = x^2 + 2*d0*x + 4*k-r"
    )

    print()
    print("JOINT RELATION")
    print(
        "  h(x,g) = 4*g^2-(N+1-d0-x)^2+16N"
    )

    print()
    print("QUESTION")
    print(
        "  Are the g-free relations from Experiment 421"
    )
    print(
        "  genuinely new modulo the exact ideal <f>?"
    )

    print()
    print("CLASSIFICATION")
    print(
        "  GFREE_INHERITED       = g-free and exactly divisible by f"
    )
    print(
        "  GFREE_NEW_K_DEPENDENT = g-free, not divisible by f, contains k"
    )
    print(
        "  GFREE_NEW_X_ONLY      = g-free, not divisible by f, x-only"
    )
    print(
        "  JOINT_VANISHING       = vanishes but contains g"
    )
    print(
        "  NONVANISHING          = does not vanish"
    )
    print(
        "  NONINTEGRAL           = scaled vector is not integral"
    )

    print()
    print("RULES")
    print(
        "  exact integer arithmetic only"
    )
    print(
        "  no resultants"
    )
    print(
        "  no Groebner basis"
    )
    print(
        "  no symbolic factorization"
    )
    print(
        "  exact polynomial division by f"
    )
    print(
        "  LLL = diagnostic only"
    )

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

    enumeration_limit = 200000

    print()
    print("CONFIGURATION")
    print(
        f"  instances       = {len(instances)}"
    )
    print(
        f"  K bits          = {unknown_bits}"
    )
    print(
        "  f/h shifts      = (1,1,1)"
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
            "  fpylll unavailable."
        )
        print(
            "  Install with:"
        )
        print(
            "    pip install fpylll cysignals"
        )
        return

    generated_instances = [
        make_instance(
            p,
            q,
        )
        for p, q in instances
    ]

    results: List[Result] = []

    audit_map: Dict[
        Tuple[int, int],
        List[Audit],
    ] = {}

    # -------------------------------------------------------------------------
    # Compact summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT IDEAL-NOVELTY SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap          x        X      cand unique "
        "inherit new-k new-x joint"
    )

    print("-" * 120)

    for index, inst in enumerate(
        generated_instances,
        start=1,
    ):

        for u in unknown_bits:

            result, audits = run_test(
                index,
                inst,
                u,
                enumeration_limit,
            )

            results.append(
                result
            )

            audit_map[
                (
                    index,
                    u,
                )
            ] = audits

            print(
                f"{index:2d} "
                f"{u:3d} "
                f"{inst.gap:10d} "
                f"{result.x_true:10d} "
                f"{result.X:8d} "
                f"{result.candidate_count:9d} "
                f"{'YES' if result.unique else 'NO ':>6} "
                f"{result.inherited:7d} "
                f"{result.new_k:5d} "
                f"{result.new_x:5d} "
                f"{result.joint:5d}"
            )

    # -------------------------------------------------------------------------
    # Global classification
    # -------------------------------------------------------------------------

    total_vectors = sum(
        len(audit_map[key])
        for key in audit_map
    )

    print()
    print("=" * 120)
    print("GLOBAL IDEAL-NOVELTY CLASSIFICATION")
    print("=" * 120)

    print(
        f"  total tests                    = {len(results)}"
    )

    print(
        f"  reduced vectors                = {total_vectors}"
    )

    print(
        "  exact g-free inherited        = "
        f"{sum(r.inherited for r in results)}"
    )

    print(
        "  new g-free k-dependent        = "
        f"{sum(r.new_k for r in results)}"
    )

    print(
        "  new g-free x-only             = "
        f"{sum(r.new_x for r in results)}"
    )

    print(
        "  joint root-vanishing          = "
        f"{sum(r.joint for r in results)}"
    )

    print(
        "  nonvanishing                   = "
        f"{sum(r.nonzero for r in results)}"
    )

    print(
        "  nonintegral                   = "
        f"{sum(r.nonintegral for r in results)}"
    )

    print(
        "  unique exact recoveries        = "
        f"{sum(r.unique for r in results)}/{len(results)}"
    )

    print()
    print(
        "  primitive distinct new-k       = "
        f"{sum(r.primitive_new_k for r in results)}"
    )

    print(
        "  primitive distinct new-x       = "
        f"{sum(r.primitive_new_x for r in results)}"
    )

    # -------------------------------------------------------------------------
    # Aggregate
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests unique inherit new-k new-x joint primitive-k primitive-x"
    )

    print("-" * 120)

    for u in unknown_bits:

        subset = [
            result
            for result in results
            if result.u == u
        ]

        print(
            f"{u:2d} "
            f"{len(subset):6d} "
            f"{sum(r.unique for r in subset):6d} "
            f"{sum(r.inherited for r in subset):7d} "
            f"{sum(r.new_k for r in subset):5d} "
            f"{sum(r.new_x for r in subset):5d} "
            f"{sum(r.joint for r in subset):5d} "
            f"{sum(r.primitive_new_k for r in subset):12d} "
            f"{sum(r.primitive_new_x for r in subset):12d}"
        )

    # -------------------------------------------------------------------------
    # Detailed genuinely new relations
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("GENUINELY NEW G-FREE RELATIONS")
    print("=" * 120)

    detail_found = False

    for result in results:

        if (
            result.new_k == 0
            and result.new_x == 0
        ):
            continue

        detail_found = True

        print()
        print(
            f"INSTANCE {result.instance}, "
            f"UNKNOWN K BITS = {result.u}"
        )

        print(
            f"  gap             = {result.gap}"
        )

        print(
            f"  true x          = {result.x_true}"
        )

        print(
            f"  X               = {result.X}"
        )

        print(
            f"  candidate count = {result.candidate_count}"
        )

        print(
            f"  unique          = {result.unique}"
        )

        print(
            f"  inherited       = {result.inherited}"
        )

        print(
            f"  new k-dependent = {result.new_k}"
        )

        print(
            f"  new x-only      = {result.new_x}"
        )

        print(
            f"  primitive k     = {result.primitive_new_k}"
        )

        print(
            f"  primitive x     = {result.primitive_new_x}"
        )

        audits = audit_map[
            (
                result.instance,
                result.u,
            )
        ]

        shown_k = 0
        shown_x = 0

        for audit in audits:

            if (
                audit.classification
                == "GFREE_NEW_K_DEPENDENT"
                and shown_k < 3
                and audit.canonical is not None
            ):

                polynomial = canonical_to_poly(
                    audit.canonical
                )

                print(
                    "  NEW-K"
                    f" row={audit.index:2d}"
                    f" norm={audit.norm}"
                    f" degree={audit.degree}"
                    f" terms={audit.terms}"
                    f" P={format_poly(polynomial)}"
                )

                shown_k += 1

            elif (
                audit.classification
                == "GFREE_NEW_X_ONLY"
                and shown_x < 3
                and audit.canonical is not None
            ):

                polynomial = canonical_to_poly(
                    audit.canonical
                )

                print(
                    "  NEW-X"
                    f" row={audit.index:2d}"
                    f" norm={audit.norm}"
                    f" degree={audit.degree}"
                    f" terms={audit.terms}"
                    f" P={format_poly(polynomial)}"
                )

                shown_x += 1

    if not detail_found:
        print(
            "  No g-free relation outside <f> was detected."
        )

    # -------------------------------------------------------------------------
    # Distinct primitive relations
    # -------------------------------------------------------------------------

    primitive_k_relations = set()
    primitive_x_relations = set()

    for audits in audit_map.values():

        for audit in audits:

            if audit.canonical is None:
                continue

            if (
                audit.classification
                == "GFREE_NEW_K_DEPENDENT"
            ):
                primitive_k_relations.add(
                    audit.canonical
                )

            elif (
                audit.classification
                == "GFREE_NEW_X_ONLY"
            ):
                primitive_x_relations.add(
                    audit.canonical
                )

    print()
    print("=" * 120)
    print("PRIMITIVE RELATION INVENTORY")
    print("=" * 120)

    print(
        f"  distinct primitive new-k relations = "
        f"{len(primitive_k_relations)}"
    )

    print(
        f"  distinct primitive new-x relations = "
        f"{len(primitive_x_relations)}"
    )

    if primitive_x_relations:

        print()
        print("  NEW-X PRIMITIVE RELATIONS")

        shown = 0

        for key in sorted(
            primitive_x_relations
        ):

            if shown >= 10:
                print(
                    "  ..."
                )
                break

            print(
                "  ",
                format_poly(
                    canonical_to_poly(
                        key
                    ),
                    max_terms=20,
                ),
            )

            shown += 1

    # -------------------------------------------------------------------------
    # Final interpretation
    # -------------------------------------------------------------------------

    total_new_k = sum(
        result.new_k
        for result in results
    )

    total_new_x = sum(
        result.new_x
        for result in results
    )

    total_inherited = sum(
        result.inherited
        for result in results
    )

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 421 found g-free relations."
    )

    print(
        "  Experiment 422 tests whether those relations"
    )

    print(
        "  are actually outside the principal ideal <f>."
    )

    print()

    print(
        "  GFREE_INHERITED means the relation is exactly"
    )

    print(
        "  divisible by f and therefore does not provide"
    )

    print(
        "  a new relation modulo <f>."
    )

    print()

    print(
        "  GFREE_NEW_K_DEPENDENT is genuinely outside <f>,"
    )

    print(
        "  but still depends on k."
    )

    print()

    print(
        "  GFREE_NEW_X_ONLY is stronger: it is outside <f>"
    )

    print(
        "  and contains neither k nor g."
    )

    print()

    print(
        "  This is still only an exact algebraic/lattice"
    )

    print(
        "  diagnostic. It does not establish factor recovery."
    )

    # -------------------------------------------------------------------------
    # Final status
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXPERIMENT 422 FINAL STATUS")
    print("=" * 120)

    print(
        "  ALL EXACT INTERNAL CHECKS = TRUE"
    )

    print(
        f"  G-FREE INHERITED = {total_inherited}"
    )

    print(
        f"  NEW G-FREE K-DEPENDENT = {total_new_k}"
    )

    print(
        f"  NEW G-FREE X-ONLY = {total_new_x}"
    )

    if total_new_x > 0:
        print(
            "  NEW X-ONLY RELATION DETECTED = TRUE"
        )

    elif total_new_k > 0:
        print(
            "  NEW G-FREE K-DEPENDENT RELATION DETECTED = TRUE"
        )

    else:
        print(
            "  NO GENUINELY NEW G-FREE RELATION DETECTED = TRUE"
        )

    print("=" * 120)
    print("EXPERIMENT 422 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()