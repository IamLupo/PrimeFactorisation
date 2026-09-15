#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
KAPPA EXPERIMENT 121
C/D KERNEL COMPRESSION
CACHED QUADRATIC-BRANCH DATASET
LOW-DEGREE (r,m) COEFFICIENT LAW SEARCH
EXACT SYMBOLIC ARITHMETIC
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================

Goal
----
Experiment 120 established the exact decomposition

    Q(N,S) = C(N,X) + (2S-1) D(N,X),
    X = S(S-1),

but repeated symbolic recomputation made the experiment too slow.

This experiment changes only the computational strategy:

    1. Construct each Q exactly once.
    2. Construct C,D exactly once.
    3. Cache everything.
    4. Never recompute quotient or decomposition during analysis.
    5. Inspect the coefficient arrays of C and D as functions of

           r = (k-1)/2
           m = (ell-1)/2.

The main test is whether individual C/D coefficients are low-degree
polynomials in (r,m), with exact leave-one-ell-out validation.

This is a structural experiment only.
It does NOT claim an N-only factoring algorithm.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

MAX_ELL = 23

# To keep the first run fast, use all odd k<ell.
# The dataset is still small:
#   ell=3..23 odd
#   total records = 66
#
# Set to False if you want to restrict the detector family later.
USE_ALL_ODD_K = True

# Polynomial degree search in r,m.
# Total degree <= MAX_BIV_DEG.
MAX_BIV_DEG = 4

# We only attempt a coefficient law when enough observations exist.
MIN_POINTS_FOR_FIT = 8


# ---------------------------------------------------------------------------
# SYMBOLS
# ---------------------------------------------------------------------------

N, S, X = sp.symbols("N S X")
r_sym, m_sym = sp.symbols("r m")


# ---------------------------------------------------------------------------
# DATA STRUCTURE
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CDRecord:
    k: int
    ell: int
    r: int
    m: int
    Q: sp.Expr
    C: sp.Expr
    D: sp.Expr


# ---------------------------------------------------------------------------
# BASIC SYMMETRIC POWER SUMS
# ---------------------------------------------------------------------------

def power_sum_table(max_power: int) -> List[sp.Expr]:
    """
    P_j = p^j + q^j
    with
        P_0 = 2
        P_1 = S
        P_j = S P_{j-1} - N P_{j-2}.
    """
    P = [sp.Integer(0)] * (max_power + 1)
    P[0] = sp.Integer(2)

    if max_power >= 1:
        P[1] = S

    for j in range(2, max_power + 1):
        P[j] = sp.expand(S * P[j - 1] - N * P[j - 2])

    return P


# ---------------------------------------------------------------------------
# SYMMETRIC MONOMIAL
# ---------------------------------------------------------------------------

def symmetric_monomial(a: int, b: int, P: Sequence[sp.Expr]) -> sp.Expr:
    """
    Return

        p^a q^b + q^a p^b

    in the symmetric variables (N,S), where

        N = pq
        S = p+q.

    If a>=b:

        N^b * (p^(a-b)+q^(a-b)).
    """
    if a < b:
        a, b = b, a

    return sp.expand(N**b * P[a - b])


# ---------------------------------------------------------------------------
# RAW NUMERATOR
# ---------------------------------------------------------------------------

def raw_detector_numerator(k: int, ell: int) -> sp.Expr:
    """
    Exact numerator:

        F =
          p^k(1+q)^ell + q^k(1+p)^ell
        - p^ell(1+q)^k - q^ell(1+p)^k.
    """

    max_power = max(k, ell)
    P = power_sum_table(max_power)

    F = sp.Integer(0)

    # First pair:
    # p^k(1+q)^ell + q^k(1+p)^ell
    for t in range(ell + 1):
        coeff = sp.binomial(ell, t)
        F += coeff * symmetric_monomial(k, t, P)

    # Second pair:
    # p^ell(1+q)^k + q^ell(1+p)^k
    for t in range(k + 1):
        coeff = sp.binomial(k, t)
        F -= coeff * symmetric_monomial(ell, t, P)

    return sp.expand(F)


# ---------------------------------------------------------------------------
# EXACT QUOTIENT
# ---------------------------------------------------------------------------

def quotient_polynomial(k: int, ell: int) -> sp.Expr:
    """
    Q = F/(S+1), exact in ZZ[N,S].

    Because the symmetric numerator is known to be divisible by S+1,
    we perform polynomial division explicitly and check the remainder.
    """

    F = raw_detector_numerator(k, ell)

    poly = sp.Poly(F, S, domain=sp.ZZ.poly_ring(N))
    divisor = sp.Poly(S + 1, S, domain=sp.ZZ.poly_ring(N))

    quotient, remainder = sp.div(poly, divisor)

    rem_expr = sp.expand(remainder.as_expr())

    if rem_expr != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder for ({k},{ell}): "
            f"{sp.factor(rem_expr)}"
        )

    return sp.expand(quotient.as_expr())


# ---------------------------------------------------------------------------
# QUADRATIC-BRANCH REDUCTION
# ---------------------------------------------------------------------------

def reduce_mod_quadratic(expr: sp.Expr) -> Tuple[sp.Expr, sp.Expr]:
    """
    Reduce expr modulo

        S^2 - S - X = 0

    to

        A(N,X) + B(N,X) S.

    Then

        Q = C + (2S-1)D

    with

        D = B/2
        C = A + B/2.
    """

    # SymPy domain with X,N as coefficient variables.
    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.poly_ring(N, X),
    )

    modulus = sp.Poly(
        S**2 - S - X,
        S,
        domain=sp.QQ.poly_ring(N, X),
    )

    remainder = poly.rem(modulus)
    rem_expr = sp.expand(remainder.as_expr())

    # Extract coefficients of S^0 and S^1.
    rem_poly = sp.Poly(
        rem_expr,
        S,
        domain=sp.QQ.poly_ring(N, X),
    )

    A = sp.expand(rem_poly.coeff_monomial(S**0))
    B = sp.expand(rem_poly.coeff_monomial(S**1))

    return A, B


def branch_decompose(Q: sp.Expr) -> Tuple[sp.Expr, sp.Expr]:
    """
    Return C,D satisfying

        Q = C(N,X) + (2S-1)D(N,X).

    Since

        C - D + 2SD = A + BS,

    we have

        D = B/2
        C = A + B/2.
    """

    A, B = reduce_mod_quadratic(Q)

    D = sp.expand(B / 2)
    C = sp.expand(A + B / 2)

    reconstructed = sp.expand(
        C.subs(X, S * (S - 1))
        + (2 * S - 1) * D.subs(X, S * (S - 1))
        - Q
    )

    if reconstructed != 0:
        raise ArithmeticError(
            "Quadratic-branch reconstruction failed: "
            f"{sp.factor(reconstructed)}"
        )

    return C, D


# ---------------------------------------------------------------------------
# CACHE BUILD
# ---------------------------------------------------------------------------

def build_dataset(max_ell: int) -> List[CDRecord]:
    records: List[CDRecord] = []

    for ell in range(3, max_ell + 1, 2):
        if USE_ALL_ODD_K:
            ks = list(range(1, ell, 2))
        else:
            # Small default family if desired later.
            ks = [k for k in (1, 3, 5, 7, 9, 11) if k < ell]

        for k in ks:
            Q = quotient_polynomial(k, ell)
            C, D = branch_decompose(Q)

            records.append(
                CDRecord(
                    k=k,
                    ell=ell,
                    r=(k - 1) // 2,
                    m=(ell - 1) // 2,
                    Q=Q,
                    C=C,
                    D=D,
                )
            )

    return records


# ---------------------------------------------------------------------------
# BASIC VALIDATION
# ---------------------------------------------------------------------------

def validate_dataset(records: Sequence[CDRecord]) -> None:
    reconstruction_failures = 0
    involution_failures = 0

    for rec in records:
        q_rebuilt = sp.expand(
            rec.C.subs(X, S * (S - 1))
            + (2 * S - 1) * rec.D.subs(X, S * (S - 1))
        )

        if sp.expand(q_rebuilt - rec.Q) != 0:
            reconstruction_failures += 1

        involuted = sp.expand(rec.Q.subs(S, 1 - S))
        if rec.Q == rec.Q.subs(S, 1 - S):
            # In this case D must vanish.
            if sp.expand(rec.D) != 0:
                involution_failures += 1
        else:
            # General branch form should transform as
            # C - (2S-1)D.
            expected = sp.expand(
                rec.C.subs(X, S * (S - 1))
                - (2 * S - 1) * rec.D.subs(X, S * (S - 1))
            )

            if sp.expand(involuted - expected) != 0:
                involution_failures += 1

    print()
    print("VALIDATION")
    print("----------")
    print(f"reconstruction failures = {reconstruction_failures}/{len(records)}")
    print(f"involution failures      = {involution_failures}/{len(records)}")

    if reconstruction_failures or involution_failures:
        raise ArithmeticError("Cached C/D validation failed.")


# ---------------------------------------------------------------------------
# POLYNOMIAL COEFFICIENT EXTRACTION
# ---------------------------------------------------------------------------

def coeff_map(expr: sp.Expr) -> Dict[Tuple[int, int], sp.Rational]:
    """
    Return coefficients

        expr = sum c[a,b] N^a X^b
    """
    if expr == 0:
        return {}

    poly = sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.QQ,
    )

    out: Dict[Tuple[int, int], sp.Rational] = {}

    for (n_deg, x_deg), coeff in poly.terms():
        out[(n_deg, x_deg)] = sp.Rational(coeff)

    return out


# ---------------------------------------------------------------------------
# SUPPORT PROFILE
# ---------------------------------------------------------------------------

def complexity_profile(records: Sequence[CDRecord]) -> None:
    max_c_x = -1
    max_c_n = -1
    max_d_x = -1
    max_d_n = -1
    max_c_terms = 0
    max_d_terms = 0

    for rec in records:
        c_map = coeff_map(rec.C)
        d_map = coeff_map(rec.D)

        if c_map:
            max_c_n = max(max_c_n, max(a for a, _ in c_map))
            max_c_x = max(max_c_x, max(b for _, b in c_map))
            max_c_terms = max(max_c_terms, len(c_map))

        if d_map:
            max_d_n = max(max_d_n, max(a for a, _ in d_map))
            max_d_x = max(max_d_x, max(b for _, b in d_map))
            max_d_terms = max(max_d_terms, len(d_map))

    print()
    print("C/D COMPLEXITY")
    print("--------------")
    print(f"C max degree_N = {max_c_n}")
    print(f"C max degree_X = {max_c_x}")
    print(f"C max terms    = {max_c_terms}")
    print(f"D max degree_N = {max_d_n}")
    print(f"D max degree_X = {max_d_x}")
    print(f"D max terms    = {max_d_terms}")


# ---------------------------------------------------------------------------
# BINOMIAL / PARITY NORMALIZATION
# ---------------------------------------------------------------------------

def try_simple_factor(expr: sp.Expr) -> sp.Expr:
    """
    Factoring is intentionally conservative.
    """
    return sp.factor(expr)


def normalized_profile(records: Sequence[CDRecord]) -> None:
    """
    For each coefficient c(N,X), inspect its N-degree and determine whether
    the coefficient's numerical multiplier has simple common factors such as
    2, (2r+1), (2m+1), etc.

    This section deliberately only prints candidates; it does not assume
    a formula.
    """

    print()
    print("NORMALIZED COEFFICIENT SAMPLES")
    print("-------------------------------")

    # Keep output small.
    samples = 0

    for rec in records:
        if samples >= 20:
            break

        cmap = coeff_map(rec.C)
        dmap = coeff_map(rec.D)

        if cmap:
            key = max(cmap.keys(), key=lambda t: (t[0], t[1]))
            coeff = cmap[key]
            print(
                f"C ({rec.k},{rec.ell}) "
                f"N^{key[0]} X^{key[1]}: "
                f"{try_simple_factor(coeff)}"
            )
            samples += 1

        if samples >= 20:
            break

        if dmap:
            key = max(dmap.keys(), key=lambda t: (t[0], t[1]))
            coeff = dmap[key]
            print(
                f"D ({rec.k},{rec.ell}) "
                f"N^{key[0]} X^{key[1]}: "
                f"{try_simple_factor(coeff)}"
            )
            samples += 1


# ---------------------------------------------------------------------------
# BIVARIATE MONOMIAL BASIS
# ---------------------------------------------------------------------------

def monomial_pairs(total_degree: int) -> List[Tuple[int, int]]:
    basis: List[Tuple[int, int]] = []

    for i in range(total_degree + 1):
        for j in range(total_degree + 1 - i):
            basis.append((i, j))

    return basis


def design_matrix(points: Sequence[Tuple[int, int]],
                  degree: int) -> sp.Matrix:
    basis = monomial_pairs(degree)

    rows = []
    for r, m in points:
        rows.append([sp.Integer(r) ** i * sp.Integer(m) ** j
                     for i, j in basis])

    return sp.Matrix(rows)


def fit_bivariate_polynomial(
    points: Sequence[Tuple[int, int]],
    values: Sequence[sp.Rational],
    degree: int,
) -> Optional[sp.Expr]:
    """
    Exact interpolation/fit in QQ using total degree <= degree.

    Returns a polynomial if the linear system is exactly consistent
    and has a unique solution.
    """

    basis = monomial_pairs(degree)

    if len(points) < len(basis):
        return None

    A = design_matrix(points, degree)
    b = sp.Matrix([sp.Rational(v) for v in values])

    try:
        sol = sp.linsolve((A, b))
    except Exception:
        return None

    if sol == sp.EmptySet:
        return None

    tuples = list(sol)
    if len(tuples) != 1:
        return None

    coeff_tuple = tuples[0]

    if any(v.free_symbols for v in coeff_tuple):
        return None

    poly = sp.Integer(0)
    for coeff, (i, j) in zip(coeff_tuple, basis):
        poly += coeff * r_sym**i * m_sym**j

    return sp.expand(poly)


def evaluate_bivariate(
    expr: sp.Expr,
    r_value: int,
    m_value: int,
) -> sp.Rational:
    return sp.Rational(
        sp.expand(
            expr.subs(
                {
                    r_sym: r_value,
                    m_sym: m_value,
                }
            )
        )
    )


# ---------------------------------------------------------------------------
# COEFFICIENT TABLE
# ---------------------------------------------------------------------------

def coefficient_table(
    records: Sequence[CDRecord],
    which: str,
) -> Dict[Tuple[int, int], Dict[Tuple[int, int], sp.Rational]]:
    """
    Return:

        (N_degree, X_degree)
            -> {(r,m): coefficient}

    Missing monomials are treated as zero only after the global support
    has been determined.
    """

    result: Dict[
        Tuple[int, int],
        Dict[Tuple[int, int], sp.Rational]
    ] = {}

    all_keys = set()

    raw_maps = {}

    for rec in records:
        expr = rec.C if which == "C" else rec.D
        cmap = coeff_map(expr)
        raw_maps[(rec.r, rec.m)] = cmap
        all_keys.update(cmap.keys())

    for key in sorted(all_keys):
        result[key] = {}

        for rec in records:
            result[key][(rec.r, rec.m)] = raw_maps[(rec.r, rec.m)].get(
                key,
                sp.Rational(0),
            )

    return result


# ---------------------------------------------------------------------------
# LEAVE-ONE-ELL-OUT LAW SEARCH
# ---------------------------------------------------------------------------

@dataclass
class LawResult:
    which: str
    key: Tuple[int, int]
    degree: int
    polynomial: sp.Expr
    tested: int
    failures: int


def loo_search_for_key(
    records: Sequence[CDRecord],
    which: str,
    key: Tuple[int, int],
    max_degree: int,
) -> Optional[LawResult]:
    table = coefficient_table(records, which).get(key, {})

    if not table:
        return None

    unique_m = sorted(set(m for _, m in table.keys()))

    best: Optional[LawResult] = None

    for degree in range(max_degree + 1):
        failures_total = 0
        tests_total = 0
        all_success = True
        representative_poly: Optional[sp.Expr] = None

        # Leave out an entire ell-level, i.e. all records with same m.
        for holdout_m in unique_m:
            train_points = [
                p for p in table.keys()
                if p[1] != holdout_m
            ]

            test_points = [
                p for p in table.keys()
                if p[1] == holdout_m
            ]

            if len(train_points) < len(monomial_pairs(degree)):
                all_success = False
                break

            train_values = [table[p] for p in train_points]

            poly = fit_bivariate_polynomial(
                train_points,
                train_values,
                degree,
            )

            if poly is None:
                all_success = False
                break

            if representative_poly is None:
                representative_poly = poly

            for r_val, m_val in test_points:
                predicted = evaluate_bivariate(poly, r_val, m_val)
                actual = table[(r_val, m_val)]

                tests_total += 1

                if sp.expand(predicted - actual) != 0:
                    failures_total += 1

        if all_success and representative_poly is not None:
            best = LawResult(
                which=which,
                key=key,
                degree=degree,
                polynomial=representative_poly,
                tested=tests_total,
                failures=failures_total,
            )

            if failures_total == 0:
                return best

    return best


# ---------------------------------------------------------------------------
# RUN LAW SEARCH
# ---------------------------------------------------------------------------

def search_kernel_laws(
    records: Sequence[CDRecord],
    which: str,
    max_degree: int,
) -> List[LawResult]:
    table = coefficient_table(records, which)

    results: List[LawResult] = []

    for key in sorted(table.keys()):
        result = loo_search_for_key(
            records,
            which,
            key,
            max_degree,
        )

        if result is not None:
            results.append(result)

    return results


# ---------------------------------------------------------------------------
# SUMMARY OF LAW SEARCH
# ---------------------------------------------------------------------------

def print_law_summary(results: Sequence[LawResult], which: str) -> None:
    exact = [r for r in results if r.failures == 0]

    print()
    print(f"{which} KERNEL LAW SEARCH")
    print("----------------------")
    print(f"supported coefficient slots = {len(results)}")
    print(f"exact LOO polynomial laws    = {len(exact)}")

    if not results:
        print("No usable coefficient law found.")
        return

    degree_counts: Dict[int, int] = {}
    for r in exact:
        degree_counts[r.degree] = degree_counts.get(r.degree, 0) + 1

    print("degree distribution:")
    for degree in sorted(degree_counts):
        print(f"  degree <= {degree}: {degree_counts[degree]}")

    print()
    print("sample exact laws:")

    shown = 0
    for result in exact:
        n_deg, x_deg = result.key

        print(
            f"  {which}[N^{n_deg} X^{x_deg}] "
            f"(degree={result.degree}) = "
            f"{sp.factor(result.polynomial)}"
        )

        shown += 1
        if shown >= 20:
            break


# ---------------------------------------------------------------------------
# NORMALIZATION BY NATURAL r/m FACTORS
# ---------------------------------------------------------------------------

def inspect_top_rows(
    records: Sequence[CDRecord],
    which: str,
) -> None:
    """
    Look only at numerically simplest coefficient slots and compare them
    against elementary combinatorial quantities.

    This is exploratory and deliberately does not fit a theorem.
    """

    table = coefficient_table(records, which)

    print()
    print(f"{which} TOP-ROW STRUCTURE")
    print("-----------------------")

    printed = 0

    # Prefer highest X-degree, then highest N-degree.
    keys = sorted(
        table.keys(),
        key=lambda key: (-key[1], -key[0])
    )

    for key in keys:
        values = table[key]

        # Print only slots which are nonzero on at least one point.
        if not any(v != 0 for v in values.values()):
            continue

        # First few data points.
        examples = []
        for point in sorted(values.keys())[:8]:
            examples.append((point, values[point]))

        print(f"{which}[N^{key[0]} X^{key[1]}]: {examples}")

        printed += 1
        if printed >= 8:
            break


# ---------------------------------------------------------------------------
# HOLDOUT BY ELL
# ---------------------------------------------------------------------------

def extended_ell_validation(
    records: Sequence[CDRecord],
    max_degree: int,
) -> None:
    """
    Conservative check:
    discover coefficient laws using ell<=15 and test ell=17,19,21,23.

    Because the laws are already LOO-tested, this provides an additional
    forward-weight check.
    """

    train = [rec for rec in records if rec.m <= 8]
    future = [rec for rec in records if rec.m > 8]

    print()
    print("FORWARD ELL HOLDOUT")
    print("-------------------")
    print(f"training records = {len(train)}")
    print(f"future records   = {len(future)}")

    total_slots = 0
    exact_slots = 0
    future_failures = 0

    for which in ("C", "D"):
        train_table = coefficient_table(train, which)
        full_table = coefficient_table(records, which)

        for key, train_values in train_table.items():
            # Need enough distinct points.
            if len(train_values) < MIN_POINTS_FOR_FIT:
                continue

            found = None

            for degree in range(max_degree + 1):
                points = list(train_values.keys())
                values = list(train_values.values())

                poly = fit_bivariate_polynomial(
                    points,
                    values,
                    degree,
                )

                if poly is not None:
                    # Exact training fit.
                    ok_train = True
                    for p in points:
                        if evaluate_bivariate(poly, *p) != train_values[p]:
                            ok_train = False
                            break

                    if ok_train:
                        found = poly
                        break

            if found is None:
                continue

            total_slots += 1

            slot_ok = True

            for rec in future:
                actual = full_table[key].get(
                    (rec.r, rec.m),
                    sp.Rational(0),
                )

                predicted = evaluate_bivariate(
                    found,
                    rec.r,
                    rec.m,
                )

                if predicted != actual:
                    slot_ok = False
                    future_failures += 1

            if slot_ok:
                exact_slots += 1

    print(f"forward-test coefficient slots = {total_slots}")
    print(f"slots passing all future ell   = {exact_slots}")
    print(f"future coefficient failures     = {future_failures}")


# ---------------------------------------------------------------------------
# CACHE CONSISTENCY
# ---------------------------------------------------------------------------

def cache_consistency_check(records: Sequence[CDRecord]) -> None:
    """
    This deliberately checks only the cached objects.
    No symbolic quotient recomputation.
    """

    failures = 0

    for rec in records:
        C2 = sp.expand(rec.C)
        D2 = sp.expand(rec.D)

        if C2 != rec.C or D2 != rec.D:
            failures += 1

    print()
    print("CACHE CONSISTENCY")
    print("-----------------")
    print(f"cache failures = {failures}/{len(records)}")

    if failures:
        raise ArithmeticError("Cached dataset corruption detected.")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("KAPPA EXPERIMENT 121")
    print("C/D KERNEL COMPRESSION")
    print("CACHED QUADRATIC-BRANCH DATASET")
    print("LOW-DEGREE (r,m) COEFFICIENT LAW SEARCH")
    print("EXACT SYMBOLIC ARITHMETIC")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    total_start = time.perf_counter()

    # -----------------------------------------------------------------------
    # CACHE BUILD
    # -----------------------------------------------------------------------

    print()
    print("1. BUILDING EXACT C/D CACHE")
    print("----------------------------")

    t0 = time.perf_counter()
    records = build_dataset(MAX_ELL)
    build_time = time.perf_counter() - t0

    print(f"records = {len(records)}")
    print(f"build time = {build_time:.6f}s")

    # -----------------------------------------------------------------------
    # VALIDATION
    # -----------------------------------------------------------------------

    t0 = time.perf_counter()
    validate_dataset(records)
    print(f"validation time = {time.perf_counter() - t0:.6f}s")

    # -----------------------------------------------------------------------
    # CACHE CONSISTENCY
    # -----------------------------------------------------------------------

    cache_consistency_check(records)

    # -----------------------------------------------------------------------
    # COMPLEXITY
    # -----------------------------------------------------------------------

    complexity_profile(records)

    # -----------------------------------------------------------------------
    # NORMALIZED SAMPLES
    # -----------------------------------------------------------------------

    normalized_profile(records)

    # -----------------------------------------------------------------------
    # C KERNEL LAW SEARCH
    # -----------------------------------------------------------------------

    t0 = time.perf_counter()
    c_results = search_kernel_laws(
        records,
        "C",
        MAX_BIV_DEG,
    )
    c_time = time.perf_counter() - t0

    print()
    print(f"C law-search time = {c_time:.6f}s")
    print_law_summary(c_results, "C")

    # -----------------------------------------------------------------------
    # D KERNEL LAW SEARCH
    # -----------------------------------------------------------------------

    t0 = time.perf_counter()
    d_results = search_kernel_laws(
        records,
        "D",
        MAX_BIV_DEG,
    )
    d_time = time.perf_counter() - t0

    print()
    print(f"D law-search time = {d_time:.6f}s")
    print_law_summary(d_results, "D")

    # -----------------------------------------------------------------------
    # STRUCTURAL SAMPLES
    # -----------------------------------------------------------------------

    inspect_top_rows(records, "C")
    inspect_top_rows(records, "D")

    # -----------------------------------------------------------------------
    # FORWARD HOLDOUT
    # -----------------------------------------------------------------------

    t0 = time.perf_counter()
    extended_ell_validation(records, MAX_BIV_DEG)
    holdout_time = time.perf_counter() - t0

    print()
    print(f"holdout time = {holdout_time:.6f}s")

    # -----------------------------------------------------------------------
    # FINAL STATUS
    # -----------------------------------------------------------------------

    exact_c = sum(r.failures == 0 for r in c_results)
    exact_d = sum(r.failures == 0 for r in d_results)

    total_time = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print("The experiment does NOT regenerate the symbolic detector during")
    print("analysis. Each Q,C,D triple is constructed exactly once and then")
    print("reused from the cache.")

    print()
    print("The structural target is:")
    print()
    print("    C(N,X) coefficient = f(r,m)")
    print("    D(N,X) coefficient = g(r,m)")
    print()
    print("where")
    print("    r=(k-1)/2")
    print("    m=(ell-1)/2")
    print()
    print(f"exact C coefficient laws found = {exact_c}")
    print(f"exact D coefficient laws found = {exact_d}")

    print()
    print("Interpretation:")
    print("  * many exact low-degree laws -> genuine kernel compression candidate")
    print("  * few/no exact laws -> change basis before doing more symbolic work")
    print("  * forward-ell success -> substantially stronger evidence than")
    print("    interpolation on one fixed ell range")

    print()
    print(f"total runtime = {total_time:.6f}s")
    print("=" * 78)
    print("EXPERIMENT 121 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

