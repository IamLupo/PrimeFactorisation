#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 111
POST-BOUNDARY HYPERGEOMETRIC STRIP
EXACT KERNEL -> NEWTON COEFFICIENTS -> r>k CLASSIFICATION

GOAL
----
Experiment 110 established the terminal boundary:

    r < k:
        leading(c) = -C(k,r)

    r = k:
        leading(c) = C(k+1,2)

The unexplained region is:

    r > k.

This experiment isolates

    r = k + d, d = 1,2,3,4

and asks whether the resulting coefficients admit exact formulas in

    a = (k-1)/2
    b = (ell-1)/2.

The experiment is symbolic/oracle-only.
It performs NO factor search and does not use p,q in the analysis stage.

NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")
a, b = sp.symbols("a b")

MAX_ELL = 17
MAX_D = 4


# ============================================================================
# BASIC UTILITIES
# ============================================================================

def odd_weights(max_ell: int) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            out.append((k, ell))

    return out


def factor_clean(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.expand(expr))


# ============================================================================
# DIRECT PAPER DETECTOR
# ============================================================================

def detector_raw(k: int, ell: int) -> sp.Expr:
    """
    F_(k,ell) =
        p^k (1+q)^ell + q^k (1+p)^ell
        - p^ell (1+q)^k - q^ell (1+p)^k.
    """

    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def symmetric_in_NS(expr: sp.Expr) -> sp.Expr:
    """
    Convert a symmetric polynomial in p,q into S=p+q, N=pq.

    SymPy's symmetrize() is used only as an exact symbolic transformation.
    """

    result = sp.symmetrize(
        sp.expand(expr),
        [p, q],
        formal=True,
    )

    symmetric_part, remainder, mapping = result

    if remainder != 0:
        raise ArithmeticError(
            f"Symmetric reduction left remainder: {remainder}"
        )

    # mapping is [(s1, p+q), (s2, p*q)] in current SymPy versions.
    if len(mapping) != 2:
        raise ArithmeticError(
            f"Unexpected symmetric mapping: {mapping}"
        )

    s1, s2 = mapping[0][0], mapping[1][0]

    converted = sp.expand(
        symmetric_part.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )

    return converted


def quotient_polynomial(k: int, ell: int) -> sp.Expr:
    """
    Divide F(p,q) by p+q+1, then convert the symmetric quotient to N,S.
    """

    raw = detector_raw(k, ell)

    divisor = p + q + 1

    quotient, remainder = sp.div(
        raw,
        divisor,
        domain=sp.QQ,
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            f"Detector {(k,ell)} is not divisible by p+q+1."
        )

    return sp.expand(symmetric_in_NS(quotient))


# ============================================================================
# NEWTON POWER-SUM BASIS
# ============================================================================

def newton_basis(max_j: int) -> List[sp.Expr]:
    """
    P_0=2
    P_1=S
    P_j=S P_{j-1} - N P_{j-2}.
    """

    P = [sp.Integer(2)]

    if max_j >= 1:
        P.append(S)

    for j in range(2, max_j + 1):
        P.append(
            sp.expand(S * P[-1] - N * P[-2])
        )

    return P


def newton_coefficients(poly: sp.Expr) -> Dict[int, sp.Expr]:
    """
    Write

        poly(N,S) = C0(N) + sum_j c_j(N) P_j.

    The conversion is triangular because P_j has leading S-degree j
    with coefficient 1.
    """

    poly = sp.Poly(
        sp.expand(poly),
        S,
    )

    if poly.is_zero:
        return {}

    max_j = poly.degree()
    P = newton_basis(max_j)

    residual = sp.expand(poly.as_expr())
    coeffs: Dict[int, sp.Expr] = {}

    for j in range(max_j, 0, -1):
        residual_poly = sp.Poly(residual, S)
        c = sp.expand(residual_poly.coeff_monomial(S**j))

        if c != 0:
            coeffs[j] = sp.factor(c)
            residual = sp.expand(
                residual - c * P[j]
            )

    return coeffs


# ============================================================================
# DATA OBJECT
# ============================================================================

@dataclass(frozen=True)
class StripDatum:
    k: int
    ell: int
    d: int
    r: int
    j: int
    coefficient: sp.Expr
    degree_N: int
    leading: sp.Expr
    next_coefficient: sp.Expr


# ============================================================================
# EXACT COEFFICIENT EXTRACTION
# ============================================================================

def extract_strip_datum(
    k: int,
    ell: int,
    d: int,
) -> Optional[StripDatum]:

    r = k + d
    j = ell - 1 - r

    # Newton moments start at P1. So j must be >= 1.
    if j < 1:
        return None

    q_poly = quotient_polynomial(k, ell)
    coeffs = newton_coefficients(q_poly)

    coefficient = sp.expand(coeffs.get(j, 0))

    poly_N = sp.Poly(
        coefficient,
        N,
    )

    if poly_N.is_zero:
        degree_N = -1
        leading = sp.Integer(0)
        next_coefficient = sp.Integer(0)
    else:
        degree_N = int(poly_N.degree())
        leading = sp.factor(poly_N.LC())

        if degree_N >= 1:
            next_coefficient = sp.factor(
                poly_N.coeff_monomial(N**(degree_N - 1))
            )
        else:
            next_coefficient = sp.Integer(0)

    return StripDatum(
        k=k,
        ell=ell,
        d=d,
        r=r,
        j=j,
        coefficient=coefficient,
        degree_N=degree_N,
        leading=leading,
        next_coefficient=next_coefficient,
    )


# ============================================================================
# DATA GENERATION
# ============================================================================

def build_dataset(max_ell: int, max_d: int) -> List[StripDatum]:
    data: List[StripDatum] = []

    for k, ell in odd_weights(max_ell):
        for d in range(1, max_d + 1):
            item = extract_strip_datum(k, ell, d)

            if item is not None:
                data.append(item)

    return data


# ============================================================================
# EXACT POLYNOMIAL FITTING IN HALF-WEIGHTS
# ============================================================================

def monomial_basis(total_degree: int) -> List[Tuple[int, int]]:
    basis: List[Tuple[int, int]] = []

    for da in range(total_degree + 1):
        for db in range(total_degree + 1 - da):
            basis.append((da, db))

    return basis


def polynomial_fit_2var(
    samples: Sequence[Tuple[int, int, sp.Expr]],
    degree: int,
) -> Optional[sp.Expr]:
    """
    Exact QQ interpolation for y = P(a,b).

    Returns None unless the system is solvable exactly.
    """

    basis = monomial_basis(degree)

    if len(samples) < len(basis):
        return None

    rows = []
    rhs = []

    for av, bv, y in samples:
        rows.append([
            sp.Integer(av) ** da * sp.Integer(bv) ** db
            for da, db in basis
        ])
        rhs.append(sp.sympify(y))

    M = sp.Matrix(rows)
    Y = sp.Matrix(rhs)

    # Require consistency and unique coefficients.
    rank_M = M.rank()

    if rank_M < len(basis):
        return None

    if M.row_join(Y).rank() != rank_M:
        return None

    solution = M.gauss_jordan_solve(Y)[0]

    expr = sp.Integer(0)

    for coeff, (da, db) in zip(solution, basis):
        expr += coeff * a**da * b**db

    return sp.factor(sp.expand(expr))


def fit_leading_law(
    data: Sequence[StripDatum],
    d: int,
    degree: int,
) -> Optional[sp.Expr]:

    subset = [
        item
        for item in data
        if item.d == d
    ]

    samples = []

    for item in subset:
        av = (item.k - 1) // 2
        bv = (item.ell - 1) // 2

        samples.append(
            (
                av,
                bv,
                item.leading,
            )
        )

    return polynomial_fit_2var(
        samples,
        degree,
    )


# ============================================================================
# LEAVE-ONE-OUT STRUCTURAL TEST
# ============================================================================

def leave_one_out_test(
    data: Sequence[StripDatum],
    d: int,
    degree: int,
) -> Tuple[bool, int, Optional[sp.Expr]]:

    subset = [
        item
        for item in data
        if item.d == d
    ]

    if len(subset) < 2:
        return False, 0, None

    failures = 0
    discovered: Optional[sp.Expr] = None

    for holdout_index in range(len(subset)):

        train = (
            subset[:holdout_index]
            + subset[holdout_index + 1:]
        )

        samples = []

        for item in train:
            samples.append(
                (
                    (item.k - 1) // 2,
                    (item.ell - 1) // 2,
                    item.leading,
                )
            )

        law = polynomial_fit_2var(
            samples,
            degree,
        )

        if law is None:
            continue

        discovered = law

        held = subset[holdout_index]

        av = sp.Integer((held.k - 1) // 2)
        bv = sp.Integer((held.ell - 1) // 2)

        predicted = sp.expand(
            law.subs(
                {
                    a: av,
                    b: bv,
                }
            )
        )

        if sp.simplify(predicted - held.leading) != 0:
            failures += 1

    return failures == 0, failures, discovered


# ============================================================================
# NORMALIZATION SEARCH
# ============================================================================

def candidate_normalizations(
    d: int,
) -> List[Tuple[str, sp.Expr]]:

    """
    Candidate factors suggested by the Pascal boundary.

    These are deliberately modest: the experiment should discover
    whether simple binomial factors explain the post-boundary regime.
    """

    K = a
    L = b

    candidates = [
        ("1", sp.Integer(1)),
        ("a+1", K + 1),
        ("b+1", L + 1),
        ("a+b+1", K + L + 1),
        ("a+b+2", K + L + 2),
        ("2*a+1", 2 * K + 1),
        ("2*b+1", 2 * L + 1),
        (
            "binom(a+b+2,2)",
            (K + L + 2) * (K + L + 1) / 2,
        ),
        (
            "binom(a+b+3,2)",
            (K + L + 3) * (K + L + 2) / 2,
        ),
    ]

    # For higher d include products of consecutive linear factors.
    if d >= 2:
        candidates.extend(
            [
                (
                    "(a+b+1)(a+b+2)",
                    (K + L + 1) * (K + L + 2),
                ),
                (
                    "(a+b+2)(a+b+3)",
                    (K + L + 2) * (K + L + 3),
                ),
                (
                    "(b+1)(a+b+1)",
                    (L + 1) * (K + L + 1),
                ),
                (
                    "(a+1)(a+b+1)",
                    (K + 1) * (K + L + 1),
                ),
            ]
        )

    return candidates


def normalized_law_search(
    data: Sequence[StripDatum],
    d: int,
    max_degree: int,
) -> None:

    subset = [
        item
        for item in data
        if item.d == d
    ]

    if not subset:
        return

    print()
    print(f"d = {d}")
    print("-" * 78)

    for name, factor in candidate_normalizations(d):

        samples = []
        usable = True

        for item in subset:

            av = sp.Integer((item.k - 1) // 2)
            bv = sp.Integer((item.ell - 1) // 2)

            denom = sp.expand(
                factor.subs(
                    {
                        a: av,
                        b: bv,
                    }
                )
            )

            if denom == 0:
                usable = False
                break

            value = sp.factor(
                sp.Rational(item.leading) / denom
            )

            samples.append(
                (
                    int(av),
                    int(bv),
                    value,
                )
            )

        if not usable:
            continue

        for degree in range(max_degree + 1):

            law = polynomial_fit_2var(
                samples,
                degree,
            )

            if law is None:
                continue

            print(
                f"normalization={name}"
            )
            print(
                f"  quotient law degree={degree}"
            )
            print(
                f"  quotient = {factor_clean(law)}"
            )
            break


# ============================================================================
# FULL COEFFICIENT FACTORIZATION
# ============================================================================

def show_sample_factorizations(
    data: Sequence[StripDatum],
) -> None:

    print()
    print("SAMPLE POST-BOUNDARY FACTORIZATIONS")
    print("-" * 78)

    seen = set()

    for item in data:

        key = (item.k, item.ell, item.d)

        if key in seen:
            continue

        seen.add(key)

        print()
        print(
            f"(k,ell)=({item.k},{item.ell}), "
            f"d={item.d}, r={item.r}, j={item.j}"
        )

        print(
            f"  coefficient = "
            f"{factor_clean(item.coefficient)}"
        )

        if item.degree_N >= 0:
            print(
                f"  degree_N   = {item.degree_N}"
            )

            print(
                f"  leading    = "
                f"{factor_clean(item.leading)}"
            )

            print(
                f"  next       = "
                f"{factor_clean(item.next_coefficient)}"
            )

        if len(seen) >= 20:
            break


# ============================================================================
# DEGREE LAW SEARCH
# ============================================================================

def degree_pattern_test(
    data: Sequence[StripDatum],
) -> None:

    print()
    print("POST-BOUNDARY DEGREE PROFILE")
    print("-" * 78)

    for d in range(1, MAX_D + 1):

        subset = [
            item for item in data if item.d == d
        ]

        if not subset:
            continue

        values = sorted(
            {
                item.degree_N
                for item in subset
            }
        )

        print(
            f"d={d}: degrees={values}"
        )

        # Check whether degree depends only on k and d.
        groups: Dict[Tuple[int, int], set] = {}

        for item in subset:
            key = (item.k, item.d)

            groups.setdefault(
                key,
                set(),
            ).add(item.degree_N)

        stable = all(
            len(v) == 1
            for v in groups.values()
        )

        print(
            f"  degree independent of ell = {stable}"
        )


# ============================================================================
# NEW-WEIGHT EXTENSION
# ============================================================================

def extended_weight_check() -> None:

    print()
    print("EXTENDED WEIGHT CHECK")
    print("-" * 78)

    checks = [
        (1, 17),
        (3, 17),
        (5, 17),
        (7, 17),
        (9, 17),
        (11, 17),
        (13, 17),
        (15, 17),
    ]

    failures = 0

    for k, ell in checks:

        q_poly = quotient_polynomial(k, ell)
        coeffs = newton_coefficients(q_poly)

        for d in range(1, MAX_D + 1):

            r = k + d
            j = ell - 1 - r

            if j < 1:
                continue

            c = sp.expand(
                coeffs.get(j, 0)
            )

            if c == 0:
                print(
                    f"FAIL ({k},{ell}) d={d}: "
                    f"zero coefficient"
                )
                failures += 1

    print(
        f"extended failures = {failures}"
    )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print("KAPPA EXPERIMENT 111")
    print("POST-BOUNDARY HYPERGEOMETRIC STRIP")
    print("EXACT KERNEL -> NEWTON COEFFICIENTS -> r>k CLASSIFICATION")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    weights = odd_weights(MAX_ELL)

    print()
    print("1. PARAMETERS")
    print("-" * 78)
    print(f"max ell = {MAX_ELL}")
    print(f"d values = 1..{MAX_D}")
    print(f"weights = {len(weights)}")

    print()
    print("2. EXACT POST-BOUNDARY DATASET")
    print("-" * 78)

    data = build_dataset(
        MAX_ELL,
        MAX_D,
    )

    print(
        f"dataset rows = {len(data)}"
    )

    # Exact sanity pass.
    for item in data:

        if item.r != item.k + item.d:
            raise ArithmeticError(
                "r/d consistency failure"
            )

        if item.j != item.ell - 1 - item.r:
            raise ArithmeticError(
                "j/r/ell consistency failure"
            )

    print("structural consistency = PASS")

    degree_pattern_test(data)

    show_sample_factorizations(data)

    print()
    print("3. EXACT LEADING-COEFFICIENT SEARCH")
    print("-" * 78)

    for d in range(1, MAX_D + 1):

        print()
        print(f"d={d}")

        for degree in range(0, 7):

            law = fit_leading_law(
                data,
                d,
                degree,
            )

            if law is None:
                continue

            ok, failures, loo_law = leave_one_out_test(
                data,
                d,
                degree,
            )

            print(
                f"  degree={degree}"
            )
            print(
                f"  global law = "
                f"{factor_clean(law)}"
            )
            print(
                f"  LOO ok={ok} failures={failures}"
            )

            if loo_law is not None:
                print(
                    f"  LOO law = "
                    f"{factor_clean(loo_law)}"
                )

            break

    print()
    print("4. NORMALIZED LEADING-COEFFICIENT SEARCH")
    print("-" * 78)

    for d in range(1, MAX_D + 1):
        normalized_law_search(
            data,
            d,
            max_degree=4,
        )

    print()
    print("5. EXTENDED WEIGHT VALIDATION")

    extended_weight_check()

    print()
    print("6. FINAL DIAGNOSTIC")
    print("=" * 78)
    print(
        """
Experiment 111 isolates the region that Experiment 110 could not
explain:

    r > k.

The main object is

    c_{ell-1-r}(N),   r = k+d.

For d=1,2,3,4 the script:

  1. derives every coefficient exactly;
  2. measures its N-degree;
  3. factors the coefficient;
  4. extracts the leading N coefficient;
  5. searches for exact polynomial laws in
         a=(k-1)/2, b=(ell-1)/2;
  6. performs leave-one-out verification;
  7. tests simple binomial normalizations;
  8. validates new ell=17 weights.

The most valuable result is not a successful interpolation.

It is an exact factorization law of the form

    lead_{d}(k,ell)
       = product/binomial expression in k,ell

that survives leave-one-out and the new ell=17 cases.

The first regime to watch is d=1.

Experiment 110 already suggests:

    lead_{d=1} = ell + k.

If confirmed exactly, the next target is d=2, followed by d=3 and d=4.

A successful sequence would provide the missing post-boundary kernel
and potentially identify a hypergeometric structure.

This remains a structural/oracle experiment and makes no claim of
an N-only factorization algorithm.
"""
    )


if __name__ == "__main__":
    main()

