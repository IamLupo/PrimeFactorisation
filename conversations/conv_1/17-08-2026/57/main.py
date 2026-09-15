#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
KAPPA EXPERIMENT 124
BINOMIAL / FACTORIAL NORMALIZATION OF C/D SCALAR EDGES
TOP N-COEFFICIENTS AT X=0
EXACT RATIONAL NORMALIZATION SEARCH
LEAVE-ONE-ELL-OUT + FORWARD HOLDOUT
NO FULL TENSOR SEARCH
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================

Goal
----
Experiment 123 showed that the four raw scalar edge functions

    L_C(r,m), M_C(r,m), L_D(r,m), M_D(r,m)

do not admit simple low-degree polynomial laws in (r,m).

Experiment 124 therefore tests whether the failure is caused by the
natural binomial/factorial envelope.

For each scalar edge we test exact normalizations built from:

    C(2m, r)
    C(2m, m-r)
    C(2m, r+1)
    C(2m, m-r-1)
    C(2m+1, r)
    C(2m+1, m-r)
    C(2m+1, r+1)
    C(2m+2, r)

plus falling-factorial / rising-factorial scales.

After normalization we search for low-degree exact laws in:

    (r,m)
    (r,m-r)
    (r,m+r)

A normalization is interesting only if it gives:

    exact training law
    + leave-one-ell-out PASS
    + forward ell PASS

The experiment reports the simplest surviving laws.

==============================================================================
"""

from __future__ import annotations

import itertools
import math
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# Global symbolic variables
# ---------------------------------------------------------------------------

P, Q = sp.symbols("P Q")
N, S, X = sp.symbols("N S X")
R, M = sp.symbols("r m")

MAX_ELL_TRAIN = 23
MAX_ELL_FORWARD = 27
MAX_DEGREE = 4

# odd weights k < ell
TARGETS = (
    (1, 3),
    (1, 5), (3, 5),
    (1, 7), (3, 7), (5, 7),
    (1, 9), (3, 9), (5, 9), (7, 9),
)


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def canon(expr: sp.Expr) -> str:
    expr = sp.factor(sp.expand(expr))
    return str(expr)


def fmt_degree(d: Optional[int]) -> str:
    return "none" if d is None else str(d)


# ---------------------------------------------------------------------------
# Direct paper detector
# ---------------------------------------------------------------------------

def raw_detector(k: int, ell: int) -> sp.Expr:
    """
    F_(k,ell) =
        p^k(1+q)^ell + q^k(1+p)^ell
        - p^ell(1+q)^k - q^ell(1+p)^k

    The quotient by p+q+1 is the detector polynomial.
    """
    F = (
        P**k * (1 + Q)**ell
        + Q**k * (1 + P)**ell
        - P**ell * (1 + Q)**k
        - Q**ell * (1 + P)**k
    )
    return sp.expand(F)


def quotient_exact_pq(F: sp.Expr) -> sp.Expr:
    divisor = P + Q + 1

    poly = sp.Poly(F, P, Q, domain=sp.ZZ)
    div_poly = sp.Poly(divisor, P, Q, domain=sp.ZZ)

    # SymPy multivariate quotient/remainder.
    q, r = sp.div(poly, div_poly, domain=sp.ZZ)

    rem = sp.expand(r.as_expr())
    if rem != 0:
        raise ArithmeticError(f"nonzero detector quotient remainder: {rem}")

    return sp.expand(q.as_expr())


# ---------------------------------------------------------------------------
# Symmetric conversion p,q -> S,N
# ---------------------------------------------------------------------------

def pq_to_NS(expr: sp.Expr) -> sp.Expr:
    """
    Convert a symmetric polynomial in P,Q to elementary symmetric
    coordinates

        S = P+Q
        N = P*Q.

    symmetrize() returns the polynomial in s1,s2.
    """
    out = sp.symmetrize(expr, [P, Q], formal=True)

    if len(out) != 3:
        raise ArithmeticError("unexpected symmetrize() return shape")

    symmetric_part, remainder, mapping = out

    rem = sp.expand(remainder)
    if rem != 0:
        raise ArithmeticError(
            f"symmetric conversion produced remainder: {canon(rem)}"
        )

    # mapping normally contains s1=P+Q, s2=P*Q
    subs = {}
    for lhs, rhs in mapping:
        if rhs == P + Q:
            subs[lhs] = S
        elif rhs == P * Q:
            subs[lhs] = N

    expr_ns = sp.expand(symmetric_part.subs(subs))
    return expr_ns


# ---------------------------------------------------------------------------
# Quadratic branch decomposition
# ---------------------------------------------------------------------------

def reduce_mod_x_relation(expr: sp.Expr) -> sp.Expr:
    """
    Reduce powers of S modulo

        S^2 - S - X = 0

    so that every result is A(N,X) + S*B(N,X).
    """
    poly = sp.Poly(sp.expand(expr), S, domain=sp.QQ.frac_field(N, X))
    relation = sp.Poly(S**2 - S - X, S, domain=sp.QQ.frac_field(N, X))

    q, r = sp.div(poly, relation, domain=sp.QQ.frac_field(N, X))
    rem = sp.expand(r.as_expr())

    return rem


def branch_decompose(Q_ns: sp.Expr) -> Tuple[sp.Expr, sp.Expr]:
    """
    Compute

        Q(N,S) = C(N,X) + (2S-1) D(N,X),

    with X=S(S-1).

    First reduce modulo S^2-S-X, giving

        Q = A(N,X) + S*B(N,X).

    Then

        D = B/2
        C = A + B/2.
    """
    rem = reduce_mod_x_relation(Q_ns)

    poly = sp.Poly(
        sp.expand(rem),
        S,
        domain=sp.QQ.frac_field(N, X),
    )

    A = sp.expand(poly.nth(0))
    B = sp.expand(poly.nth(1))

    D = sp.expand(B / 2)
    C = sp.expand(A + B / 2)

    # exact reconstruction check
    recon = sp.expand(C + (2 * S - 1) * D)
    rem2 = sp.expand(
        reduce_mod_x_relation(Q_ns - recon)
    )

    if rem2 != 0:
        raise ArithmeticError(
            f"C/D reconstruction failed: {canon(rem2)}"
        )

    return sp.expand(C), sp.expand(D)


# ---------------------------------------------------------------------------
# Scalar extraction
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EdgeRecord:
    k: int
    ell: int
    r: int
    m: int
    C: sp.Expr
    D: sp.Expr
    C_deg: int
    D_deg: Optional[int]
    C_top: sp.Expr
    C_second: sp.Expr
    D_top: sp.Expr
    D_second: sp.Expr


def poly_degree_N(expr: sp.Expr) -> Optional[int]:
    p = sp.Poly(sp.expand(expr), N, domain=sp.QQ)
    if p.is_zero:
        return None
    return int(p.degree())


def scalar_top_two(expr: sp.Expr) -> Tuple[Optional[int], sp.Expr, sp.Expr]:
    poly = sp.Poly(sp.expand(expr), N, domain=sp.QQ)

    if poly.is_zero:
        return None, sp.Integer(0), sp.Integer(0)

    d = int(poly.degree())
    L = sp.expand(poly.coeff_monomial(N**d))

    if d == 0:
        M2 = sp.Integer(0)
    else:
        M2 = sp.expand(poly.coeff_monomial(N**(d - 1)))

    return d, L, M2


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

class DetectorCache:
    def __init__(self) -> None:
        self.q_cache: Dict[Tuple[int, int], sp.Expr] = {}
        self.cd_cache: Dict[Tuple[int, int], Tuple[sp.Expr, sp.Expr]] = {}
        self.raw_cache: Dict[Tuple[int, int], sp.Expr] = {}

    def detector_q(self, k: int, ell: int) -> sp.Expr:
        key = (k, ell)
        if key in self.q_cache:
            return self.q_cache[key]

        F = raw_detector(k, ell)
        Qpq = quotient_exact_pq(F)
        Qns = pq_to_NS(Qpq)

        self.q_cache[key] = sp.expand(Qns)
        return self.q_cache[key]

    def cd(self, k: int, ell: int) -> Tuple[sp.Expr, sp.Expr]:
        key = (k, ell)
        if key in self.cd_cache:
            return self.cd_cache[key]

        Qns = self.detector_q(k, ell)
        C, D = branch_decompose(Qns)

        self.cd_cache[key] = (C, D)
        return C, D


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

def odd_values_below(ell: int) -> List[int]:
    return list(range(1, ell, 2))


def build_records(cache: DetectorCache, max_ell: int) -> List[EdgeRecord]:
    records: List[EdgeRecord] = []

    for ell in range(3, max_ell + 1, 2):
        m = (ell - 1) // 2

        for k in odd_values_below(ell):
            r = (k - 1) // 2
            C, D = cache.cd(k, ell)

            Cdeg, Ctop, Csecond = scalar_top_two(C.subs(X, 0))
            Ddeg, Dtop, Dsecond = scalar_top_two(D.subs(X, 0))

            if Cdeg is None:
                raise ArithmeticError(f"C is zero for {(k, ell)}")

            records.append(
                EdgeRecord(
                    k=k,
                    ell=ell,
                    r=r,
                    m=m,
                    C=C,
                    D=D,
                    C_deg=Cdeg,
                    D_deg=Ddeg,
                    C_top=Ctop,
                    C_second=Csecond,
                    D_top=Dtop,
                    D_second=Dsecond,
                )
            )

    return records


# ---------------------------------------------------------------------------
# Exact normalization candidates
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Normalizer:
    name: str
    fn: Callable[[int, int], Fraction]


def choose(n: int, k: int) -> Fraction:
    if k < 0 or k > n:
        return Fraction(0, 1)
    return Fraction(math.comb(n, k), 1)


def falling(n: int, j: int) -> Fraction:
    if j < 0:
        return Fraction(0, 1)
    out = 1
    for t in range(j):
        out *= (n - t)
    return Fraction(out, 1)


def rising(n: int, j: int) -> Fraction:
    if j < 0:
        return Fraction(0, 1)
    out = 1
    for t in range(j):
        out *= (n + t)
    return Fraction(out, 1)


def build_normalizers() -> List[Normalizer]:
    def safe_ratio(
        num_fn: Callable[[int, int], int],
        den_fn: Callable[[int, int], int],
        name: str,
    ) -> Normalizer:

        def f(r: int, m: int) -> Fraction:
            den = den_fn(r, m)
            if den == 0:
                return Fraction(0, 1)
            return Fraction(num_fn(r, m), den)

        return Normalizer(name, f)

    cands: List[Normalizer] = []

    # Direct binomial scales.
    cands.extend(
        [
            Normalizer("1", lambda r, m: Fraction(1, 1)),
            Normalizer(
                "C(2m,r)",
                lambda r, m: choose(2 * m, r),
            ),
            Normalizer(
                "C(2m,m-r)",
                lambda r, m: choose(2 * m, m - r),
            ),
            Normalizer(
                "C(2m,r+1)",
                lambda r, m: choose(2 * m, r + 1),
            ),
            Normalizer(
                "C(2m,m-r-1)",
                lambda r, m: choose(2 * m, m - r - 1),
            ),
            Normalizer(
                "C(2m+1,r)",
                lambda r, m: choose(2 * m + 1, r),
            ),
            Normalizer(
                "C(2m+1,m-r)",
                lambda r, m: choose(2 * m + 1, m - r),
            ),
            Normalizer(
                "C(2m+1,r+1)",
                lambda r, m: choose(2 * m + 1, r + 1),
            ),
            Normalizer(
                "C(2m+2,r)",
                lambda r, m: choose(2 * m + 2, r),
            ),
            Normalizer(
                "C(2m+2,m-r)",
                lambda r, m: choose(2 * m + 2, m - r),
            ),
        ]
    )

    # Simple factorial-like scales.
    cands.extend(
        [
            Normalizer(
                "falling(m,r)",
                lambda r, m: falling(m, r),
            ),
            Normalizer(
                "falling(m,r+1)",
                lambda r, m: falling(m, r + 1),
            ),
            Normalizer(
                "rising(r+1,m-r)",
                lambda r, m: rising(r + 1, max(0, m - r)),
            ),
            Normalizer(
                "m",
                lambda r, m: Fraction(m, 1),
            ),
            Normalizer(
                "2m+1",
                lambda r, m: Fraction(2 * m + 1, 1),
            ),
            Normalizer(
                "m-r",
                lambda r, m: Fraction(m - r, 1),
            ),
            Normalizer(
                "m+r+1",
                lambda r, m: Fraction(m + r + 1, 1),
            ),
        ]
    )

    return cands


# ---------------------------------------------------------------------------
# Coordinate systems
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Coordinates:
    name: str
    vars: Tuple[sp.Symbol, ...]
    fn: Callable[[int, int], Tuple[int, ...]]


COORDINATES = [
    Coordinates(
        "(r,m)",
        (R, M),
        lambda r, m: (r, m),
    ),
    Coordinates(
        "(r,m-r)",
        (R, M),
        lambda r, m: (r, m - r),
    ),
    Coordinates(
        "(r,m+r)",
        (R, M),
        lambda r, m: (r, m + r),
    ),
    Coordinates(
        "(m,r,m-r)",
        (R, M),
        lambda r, m: (r, m - r),
    ),
]


# ---------------------------------------------------------------------------
# Exact polynomial interpolation
# ---------------------------------------------------------------------------

def monomials_2d(max_degree: int) -> List[Tuple[int, int]]:
    out = []
    for total in range(max_degree + 1):
        for a in range(total + 1):
            b = total - a
            out.append((a, b))
    return out


def fit_exact_2d(
    points: Sequence[Tuple[int, int, sp.Rational]],
    max_degree: int,
) -> Optional[sp.Expr]:
    """
    Find an exact polynomial of total degree <= max_degree
    interpolating all supplied points.
    """
    if not points:
        return None

    monoms = monomials_2d(max_degree)

    rows = []
    rhs = []

    for x1, x2, value in points:
        rows.append(
            [
                sp.Integer(x1) ** a * sp.Integer(x2) ** b
                for a, b in monoms
            ]
        )
        rhs.append(sp.Rational(value))

    A = sp.Matrix(rows)
    b = sp.Matrix(rhs)

    # Solve exact linear system; allow redundant rows.
    sol = sp.linsolve((A, b))
    if sol is sp.EmptySet:
        return None

    sols = list(sol)
    if not sols:
        return None

    tup = sols[0]

    # Reject free parameters: we want an actual polynomial law.
    free = set()
    for expr in tup:
        free |= expr.free_symbols & {sp.Symbol(f"t{i}") for i in range(len(tup))}

    # Better explicit free-symbol check against all original symbols:
    coeff_symbols = set()
    for expr in tup:
        coeff_symbols |= expr.free_symbols

    allowed = set()
    for s in coeff_symbols:
        if str(s).startswith("tau"):
            allowed.add(s)

    if coeff_symbols & allowed:
        return None

    poly = 0
    for coeff, (a, bb) in zip(tup, monoms):
        poly += coeff * R**a * M**bb

    return sp.factor(sp.expand(poly))


def eval_exact_law(
    law: sp.Expr,
    coords: Coordinates,
    r: int,
    m: int,
) -> sp.Rational:
    x1, x2 = coords.fn(r, m)
    return sp.Rational(
        sp.expand(law.subs({R: x1, M: x2}))
    )


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LawResult:
    scalar: str
    normalizer: str
    coords: str
    degree: int
    law: sp.Expr
    train_failures: int
    loo_failures: int
    forward_failures: int
    complexity: int


SCALARS = {
    "C-top": lambda rec: sp.Rational(rec.C_top),
    "C-second": lambda rec: sp.Rational(rec.C_second),
    "D-top": lambda rec: sp.Rational(rec.D_top),
    "D-second": lambda rec: sp.Rational(rec.D_second),
}


def normalized_value(
    scalar_name: str,
    rec: EdgeRecord,
    normalizer: Normalizer,
) -> Optional[sp.Rational]:
    raw = SCALARS[scalar_name](rec)
    scale = normalizer.fn(rec.r, rec.m)

    if scale == 0:
        return None

    # Divide raw coefficient by normalization scale.
    val = Fraction(int(raw.p), int(raw.q)) / scale
    return sp.Rational(val.numerator, val.denominator)


def build_points(
    records: Sequence[EdgeRecord],
    scalar_name: str,
    normalizer: Normalizer,
    coords: Coordinates,
) -> Optional[List[Tuple[int, int, sp.Rational]]]:
    out = []

    for rec in records:
        v = normalized_value(scalar_name, rec, normalizer)
        if v is None:
            return None

        x1, x2 = coords.fn(rec.r, rec.m)
        out.append((x1, x2, v))

    return out


def law_complexity(law: sp.Expr) -> int:
    poly = sp.Poly(sp.expand(law), R, M, domain=sp.QQ)
    return len(poly.terms())


def loo_validate(
    records: Sequence[EdgeRecord],
    scalar_name: str,
    normalizer: Normalizer,
    coords: Coordinates,
    degree: int,
) -> Tuple[Optional[sp.Expr], int]:
    """
    Leave one ell out, predicting every row of the omitted ell.
    """
    failures = 0

    ells = sorted({rec.ell for rec in records})

    for held_ell in ells:
        train = [rec for rec in records if rec.ell != held_ell]
        test = [rec for rec in records if rec.ell == held_ell]

        points = build_points(train, scalar_name, normalizer, coords)
        if points is None:
            failures += len(test)
            continue

        law = fit_exact_2d(points, degree)

        if law is None:
            failures += len(test)
            continue

        for rec in test:
            v = normalized_value(scalar_name, rec, normalizer)
            if v is None:
                failures += 1
                continue

            pred = eval_exact_law(law, coords, rec.r, rec.m)
            if sp.simplify(pred - v) != 0:
                failures += 1

    # Return a global fit too, only for display.
    all_points = build_points(records, scalar_name, normalizer, coords)
    global_law = (
        fit_exact_2d(all_points, degree)
        if all_points is not None
        else None
    )

    return global_law, failures


def forward_validate(
    train_records: Sequence[EdgeRecord],
    forward_records: Sequence[EdgeRecord],
    scalar_name: str,
    normalizer: Normalizer,
    coords: Coordinates,
    degree: int,
) -> Tuple[Optional[sp.Expr], int]:
    train_points = build_points(
        train_records,
        scalar_name,
        normalizer,
        coords,
    )

    if train_points is None:
        return None, len(forward_records)

    law = fit_exact_2d(train_points, degree)

    if law is None:
        return None, len(forward_records)

    failures = 0

    for rec in forward_records:
        v = normalized_value(scalar_name, rec, normalizer)
        if v is None:
            failures += 1
            continue

        pred = eval_exact_law(law, coords, rec.r, rec.m)
        if sp.simplify(pred - v) != 0:
            failures += 1

    return law, failures


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def search_all(
    train_records: Sequence[EdgeRecord],
    forward_records: Sequence[EdgeRecord],
) -> List[LawResult]:

    normalizers = build_normalizers()
    results: List[LawResult] = []

    for scalar_name in SCALARS:
        for norm in normalizers:
            for coords in COORDINATES:
                for degree in range(MAX_DEGREE + 1):
                    points = build_points(
                        train_records,
                        scalar_name,
                        norm,
                        coords,
                    )

                    if points is None:
                        continue

                    law = fit_exact_2d(points, degree)
                    if law is None:
                        continue

                    # Global exactness is mandatory.
                    train_failures = 0
                    for rec in train_records:
                        v = normalized_value(scalar_name, rec, norm)
                        if v is None:
                            train_failures += 1
                            continue

                        pred = eval_exact_law(
                            law,
                            coords,
                            rec.r,
                            rec.m,
                        )
                        if sp.simplify(pred - v) != 0:
                            train_failures += 1

                    if train_failures != 0:
                        continue

                    # LOO.
                    _, loo_failures = loo_validate(
                        train_records,
                        scalar_name,
                        norm,
                        coords,
                        degree,
                    )

                    # Forward holdout.
                    forward_law, forward_failures = forward_validate(
                        train_records,
                        forward_records,
                        scalar_name,
                        norm,
                        coords,
                        degree,
                    )

                    if forward_law is None:
                        forward_failures = len(forward_records)

                    results.append(
                        LawResult(
                            scalar=scalar_name,
                            normalizer=norm.name,
                            coords=coords.name,
                            degree=degree,
                            law=law,
                            train_failures=train_failures,
                            loo_failures=loo_failures,
                            forward_failures=forward_failures,
                            complexity=law_complexity(law),
                        )
                    )

    return results


# ---------------------------------------------------------------------------
# Rank and report
# ---------------------------------------------------------------------------

def rank_key(x: LawResult):
    return (
        x.loo_failures != 0,
        x.forward_failures != 0,
        x.degree,
        x.complexity,
        len(str(x.law)),
    )


def print_best(results: Sequence[LawResult]) -> None:
    print()
    print("=" * 78)
    print("BEST SURVIVING NORMALIZED LAWS")
    print("=" * 78)

    if not results:
        print("No exact normalized law found.")
        return

    by_scalar: Dict[str, List[LawResult]] = {
        s: [] for s in SCALARS
    }

    for result in results:
        by_scalar[result.scalar].append(result)

    for scalar, arr in by_scalar.items():
        arr_sorted = sorted(arr, key=rank_key)

        print()
        print(scalar)
        print("-" * 78)

        shown = 0
        seen = set()

        for item in arr_sorted:
            signature = (
                item.normalizer,
                item.coords,
                item.degree,
                str(item.law),
            )

            if signature in seen:
                continue
            seen.add(signature)

            status = (
                "GLOBAL + LOO + FORWARD PASS"
                if item.loo_failures == 0
                and item.forward_failures == 0
                else
                f"LOO={item.loo_failures} "
                f"FORWARD={item.forward_failures}"
            )

            print(
                f"  norm={item.normalizer:20s} "
                f"coords={item.coords:12s} "
                f"deg={item.degree} "
                f"terms={item.complexity} "
                f"{status}"
            )
            print(f"    law = {canon(item.law)}")

            shown += 1
            if shown >= 8:
                break


def print_exact_passes(results: Sequence[LawResult]) -> None:
    winners = [
        x
        for x in results
        if x.loo_failures == 0
        and x.forward_failures == 0
    ]

    print()
    print("=" * 78)
    print("EXACT NORMALIZATION WINNERS")
    print("=" * 78)

    if not winners:
        print("No normalization survived both LOO and forward holdout.")
        return

    winners = sorted(winners, key=rank_key)

    for x in winners:
        print(
            f"{x.scalar:10s} | "
            f"norm={x.normalizer:20s} | "
            f"coords={x.coords:12s} | "
            f"degree={x.degree} | "
            f"terms={x.complexity}"
        )
        print(f"  {canon(x.law)}")


# ---------------------------------------------------------------------------
# Normalization diagnostics
# ---------------------------------------------------------------------------

def show_raw_vs_normalized(
    records: Sequence[EdgeRecord],
    scalar_name: str,
    normalizers: Sequence[Normalizer],
    limit: int = 12,
) -> None:
    print()
    print("=" * 78)
    print(f"NORMALIZATION DIAGNOSTIC: {scalar_name}")
    print("=" * 78)

    rows = records[:limit]

    for rec in rows:
        raw = SCALARS[scalar_name](rec)

        print(
            f"(k,ell)=({rec.k},{rec.ell}) "
            f"r={rec.r} m={rec.m} raw={raw}"
        )

        for norm in normalizers[:10]:
            val = normalized_value(scalar_name, rec, norm)
            if val is None:
                continue

            print(
                f"  {norm.name:20s} -> {canon(val)}"
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    total0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 124")
    print("BINOMIAL / FACTORIAL NORMALIZATION OF C/D SCALAR EDGES")
    print("TOP N-COEFFICIENTS AT X=0")
    print("SECOND N-LAYER AT X=0")
    print("EXACT NORMALIZATION SEARCH")
    print("LEAVE-ONE-ELL-OUT + FORWARD HOLDOUT")
    print("NO FULL TENSOR SEARCH")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Build cache and datasets.
    # ------------------------------------------------------------------

    cache = DetectorCache()

    t0 = time.perf_counter()
    train_records = build_records(cache, MAX_ELL_TRAIN)
    build_train_time = time.perf_counter() - t0

    t1 = time.perf_counter()
    forward_records = [
        rec
        for rec in build_records(cache, MAX_ELL_FORWARD)
        if rec.ell > MAX_ELL_TRAIN
    ]
    build_forward_time = time.perf_counter() - t1

    print()
    print("1. DATASET")
    print("-" * 78)
    print(f"training records = {len(train_records)}")
    print(f"forward records  = {len(forward_records)}")
    print(f"train ell <= {MAX_ELL_TRAIN}")
    print(f"forward ell = {MAX_ELL_FORWARD}")
    print(f"train build time   = {build_train_time:.6f}s")
    print(f"forward build time = {build_forward_time:.6f}s")

    # ------------------------------------------------------------------
    # Internal reconstruction safety.
    # ------------------------------------------------------------------

    recon_fail = 0

    for rec in train_records + forward_records:
        if sp.expand(
            rec.C + (2 * S - 1) * rec.D
        ):
            Qns = cache.detector_q(rec.k, rec.ell)

            check = sp.expand(
                rec.C + (2 * S - 1) * rec.D - Qns
            )

            rem = reduce_mod_x_relation(check)

            if rem != 0:
                recon_fail += 1

    print()
    print("2. INTERNAL C/D VALIDATION")
    print("-" * 78)
    print(f"reconstruction failures = {recon_fail}")
    if recon_fail:
        raise ArithmeticError("C/D reconstruction failed.")

    print("STATUS = PASS")

    # ------------------------------------------------------------------
    # Normalization samples.
    # ------------------------------------------------------------------

    normalizers = build_normalizers()

    print()
    print("3. NORMALIZATION CANDIDATES")
    print("-" * 78)

    for n in normalizers:
        print(f"  {n.name}")

    # ------------------------------------------------------------------
    # Small raw-vs-normalized diagnostic.
    # ------------------------------------------------------------------

    show_raw_vs_normalized(
        train_records,
        "C-top",
        normalizers,
        limit=8,
    )

    # ------------------------------------------------------------------
    # Search.
    # ------------------------------------------------------------------

    print()
    print("4. EXACT NORMALIZATION SEARCH")
    print("-" * 78)

    t2 = time.perf_counter()

    results = search_all(
        train_records,
        forward_records,
    )

    search_time = time.perf_counter() - t2

    print(f"candidate laws found = {len(results)}")
    print(f"search time = {search_time:.6f}s")

    # ------------------------------------------------------------------
    # Report.
    # ------------------------------------------------------------------

    print_best(results)
    print_exact_passes(results)

    # ------------------------------------------------------------------
    # Stability summary.
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STABILITY SUMMARY")
    print("=" * 78)

    winners = [
        x for x in results
        if x.loo_failures == 0 and x.forward_failures == 0
    ]

    for scalar in SCALARS:
        scalar_winners = [
            x for x in winners if x.scalar == scalar
        ]

        print()
        print(scalar)

        if not scalar_winners:
            print("  no exact normalization winner")
            continue

        best = sorted(scalar_winners, key=rank_key)[0]

        print(
            f"  best normalization = {best.normalizer}"
        )
        print(
            f"  coordinates        = {best.coords}"
        )
        print(
            f"  degree             = {best.degree}"
        )
        print(
            f"  terms              = {best.complexity}"
        )
        print(
            f"  law                = {canon(best.law)}"
        )

    # ------------------------------------------------------------------
    # Final diagnostic.
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
Experiment 124 changes the question from

    "Is the raw edge coefficient a simple polynomial in (r,m)?"

to

    "Does the edge coefficient have a natural binomial/factorial
     envelope that becomes simple after exact normalization?"

The primary candidates are central/near-central binomial coefficients
and factorial-type scales suggested by the Pascal structure already
seen in Experiments 109-112.

A particularly strong result would be:

    raw edge
      =
    simple binomial factor
      *
    low-degree polynomial

with

    exact training
    + leave-one-ell-out PASS
    + forward ell PASS.

A complete failure would be informative too: it would rule out this
large class of elementary hypergeometric normalizations and justify
moving to finite-difference / hypergeometric-ratio analysis.

This experiment remains purely structural and oracle-based.
It makes no claim of an N-only factorization algorithm.
"""
    )

    total = time.perf_counter() - total0
    print(f"total runtime = {total:.6f}s")
    print("=" * 78)
    print("EXPERIMENT 124 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

