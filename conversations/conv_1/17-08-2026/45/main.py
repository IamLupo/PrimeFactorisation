#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 113R
SECOND-COEFFICIENT POST-BOUNDARY LAW
ACTUAL-DEGREE / CANCELLATION SAFE
EXACT BINOMIAL KERNEL -> NEWTON COEFFICIENTS
STRICT LEAVE-ONE-ELL-OUT VALIDATION
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================

Fix over Experiment 113
------------------------
Experiment 113 incorrectly assumed the Experiment 112 leading law was
universally valid for every admissible (k,ell,d).

That is false at exceptional cancellation points.

This revision:

  1. computes the coefficient exactly;
  2. measures its ACTUAL degree;
  3. records the actual leading coefficient;
  4. classifies rows where the expected degree drops;
  5. searches the second coefficient relative to the ACTUAL degree;
  6. never aborts merely because an expected leading law cancels.

For

    c(N) = L*N^D + M*N^(D-1) + ...

we search M as a function of

    a=(k-1)/2
    b=(ell-1)/2
    d=r-k.

The leading-law from Experiment 112 is retained only as a diagnostic.
It is NOT used as a hard validation condition.
"""

from __future__ import annotations

import sys
from functools import lru_cache

import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

N, S = sp.symbols("N S")
a, b = sp.symbols("a b")


# ============================================================================
# PARAMETERS
# ============================================================================

MAX_ELL_TRAIN = 27
MAX_ELL_EXTENDED = 31
MAX_D = 10

POLY_DEGREES = range(0, 9)


# ============================================================================
# POWER SUMS
# ============================================================================

@lru_cache(maxsize=None)
def power_sum(m: int) -> sp.Expr:
    """
    P_m = p^m + q^m

    P_0 = 2
    P_1 = S
    P_m = S P_(m-1) - N P_(m-2)
    """
    if m == 0:
        return sp.Integer(2)

    if m == 1:
        return S

    p0 = sp.Integer(2)
    p1 = S

    for _ in range(2, m + 1):
        p0, p1 = p1, sp.expand(S * p1 - N * p0)

    return sp.expand(p1)


@lru_cache(maxsize=None)
def symmetric_monomial_sum(x: int, y: int) -> sp.Expr:
    """
    p^x q^y + q^x p^y
    """
    m = min(x, y)
    delta = abs(x - y)

    return sp.expand(N**m * power_sum(delta))


# ============================================================================
# DIRECT PAPER DETECTOR
# ============================================================================

@lru_cache(maxsize=None)
def detector_symmetric(k: int, ell: int) -> sp.Expr:
    """
    F_(k,l) =
      p^k(1+q)^l + q^k(1+p)^l
      - p^l(1+q)^k - q^l(1+p)^k
    """

    if not (k < ell and k % 2 == 1 and ell % 2 == 1):
        raise ValueError("Expected odd k < odd ell.")

    expr = sp.Integer(0)

    for t in range(ell + 1):
        expr += sp.binomial(ell, t) * symmetric_monomial_sum(k, t)

    for t in range(k + 1):
        expr -= sp.binomial(k, t) * symmetric_monomial_sum(ell, t)

    return sp.expand(expr)


# ============================================================================
# EXACT QUOTIENT BY S+1
# ============================================================================

@lru_cache(maxsize=None)
def quotient_polynomial(k: int, ell: int) -> sp.Poly:
    """
    Q_(k,l)(N,S) = F_(k,l)/(S+1)

    SymPy-safe implementation:
    S is the explicit polynomial generator and N belongs to the coefficient
    field QQ(N).
    """

    F = detector_symmetric(k, ell)

    field = sp.QQ.frac_field(N)

    dividend = sp.Poly(F, S, domain=field)
    divisor = sp.Poly(S + 1, S, domain=field)

    quotient, remainder = sp.div(dividend, divisor)

    if not remainder.is_zero:
        raise ArithmeticError(
            f"Division failure for ({k},{ell}): "
            f"remainder={remainder.as_expr()}"
        )

    return quotient


def coefficient_S_power(k: int, ell: int, j: int) -> sp.Expr:
    return sp.expand(quotient_polynomial(k, ell).nth(j))


# ============================================================================
# BASIC VALIDATION
# ============================================================================

def validate_quotients(max_ell: int) -> None:
    print("1. EXACT QUOTIENT VALIDATION")

    failures = 0
    total = 0

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):

            F = detector_symmetric(k, ell)
            Q = quotient_polynomial(k, ell).as_expr()

            total += 1

            if sp.expand(F - (S + 1) * Q) != 0:
                failures += 1
                print(f"  FAIL ({k},{ell})")

    print(f"quotient failures = {failures}/{total}")

    if failures:
        raise ArithmeticError("Exact quotient validation failed.")

    print("STATUS = PASS\n")


# ============================================================================
# EXPERIMENT 112 LAW -- DIAGNOSTIC ONLY
# ============================================================================

def nominal_degree(k: int, d: int) -> int:
    """
    The degree law suggested by Experiment 112.

    IMPORTANT:
    This is only a nominal prediction. Exact cancellation can lower the
    actual degree.
    """
    return k + (d - 1) // 2


def nominal_leading_law(k: int, ell: int, d: int) -> sp.Expr:
    """
    Nominal Experiment-112 leading coefficient law.
    Diagnostic only.
    """
    if d % 2 == 1:
        return sp.Integer(k + ell)

    return -sp.Rational(k + ell, 2) * (k - ell + d + 1)


# ============================================================================
# COEFFICIENT RECORD
# ============================================================================

def coefficient_record(k: int, ell: int, d: int) -> dict:
    """
    r = k+d
    j = ell-1-r.

    Extract the EXACT coefficient polynomial c_j(N).
    """
    r = k + d
    j = ell - 1 - r

    if j < 0:
        raise ValueError(
            f"Invalid post-boundary position: "
            f"k={k}, ell={ell}, d={d}, r={r}, j={j}"
        )

    polynomial = sp.Poly(
        coefficient_S_power(k, ell, j),
        N,
        domain=sp.QQ,
    )

    if polynomial.is_zero:
        degree = -1
        lead = sp.Integer(0)
        second = sp.Integer(0)
        lower = sp.Integer(0)
    else:
        degree = polynomial.degree()

        lead = sp.expand(polynomial.LC())

        if degree >= 1:
            second = sp.expand(polynomial.nth(degree - 1))
        else:
            second = sp.Integer(0)

        lower = sp.expand(
            polynomial.nth(degree - 2)
        ) if degree >= 2 else sp.Integer(0)

    D_nom = nominal_degree(k, d)
    L_nom = sp.expand(nominal_leading_law(k, ell, d))

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "a": (k - 1) // 2,
        "b": (ell - 1) // 2,
        "r": r,
        "j": j,

        "poly": polynomial.as_expr(),

        "degree": degree,
        "lead": lead,
        "second": second,
        "third": lower,

        "nom_degree": D_nom,
        "nom_lead": L_nom,

        "degree_drop": degree < D_nom,
        "leading_match": sp.expand(lead - L_nom) == 0,
    }


# ============================================================================
# DATASET
# ============================================================================

def build_dataset(max_ell: int, max_d: int) -> list[dict]:
    rows = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            for d in range(1, max_d + 1):

                if k + d > ell - 1:
                    continue

                rows.append(
                    coefficient_record(k, ell, d)
                )

    return rows


# ============================================================================
# CANCELLATION REPORT
# ============================================================================

def report_actual_degree(rows: list[dict]) -> None:
    print("2. ACTUAL DEGREE / LEADING COEFFICIENT")
    print("=" * 78)

    degree_drops = [
        row for row in rows
        if row["degree_drop"]
    ]

    leading_mismatches = [
        row for row in rows
        if not row["leading_match"]
    ]

    print(
        f"actual degree drops = "
        f"{len(degree_drops)}/{len(rows)}"
    )

    print(
        f"nominal leading mismatches = "
        f"{len(leading_mismatches)}/{len(rows)}"
    )

    if degree_drops:
        print("\nCANCELLATION CASES")

        for row in degree_drops[:100]:
            print(
                f"  ({row['k']},{row['ell']}), "
                f"d={row['d']}: "
                f"actual_degree={row['degree']} "
                f"nominal_degree={row['nom_degree']}"
            )

    print()


# ============================================================================
# SPLIT NORMAL / EXCEPTIONAL ROWS
# ============================================================================

def split_rows(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """
    Normal rows:
        actual degree == nominal degree
        AND leading coefficient matches Experiment 112.

    Exceptional rows:
        anything else.
    """
    normal = []
    exceptional = []

    for row in rows:
        if (
            row["degree"] == row["nom_degree"]
            and row["leading_match"]
        ):
            normal.append(row)
        else:
            exceptional.append(row)

    return normal, exceptional


# ============================================================================
# POLYNOMIAL FIT
# ============================================================================

def monomials_total_degree(max_degree: int) -> list[tuple[int, int]]:
    terms = []

    for total in range(max_degree + 1):
        for ia in range(total + 1):
            ib = total - ia
            terms.append((ia, ib))

    return terms


def exact_polynomial_fit(
    samples: list[tuple[int, int, sp.Expr]],
    max_degree: int,
) -> sp.Expr | None:

    if not samples:
        return None

    terms = monomials_total_degree(max_degree)

    matrix = []
    rhs = []

    for av, bv, value in samples:

        matrix.append([
            sp.Integer(av) ** ia *
            sp.Integer(bv) ** ib
            for ia, ib in terms
        ])

        rhs.append(sp.sympify(value))

    M = sp.Matrix(matrix)
    Y = sp.Matrix(rhs)

    try:
        solution_set = sp.linsolve((M, Y))
    except Exception:
        return None

    solutions = list(solution_set)

    if not solutions:
        return None

    solution = solutions[0]

    free = set()

    for value in solution:
        free |= value.free_symbols

    if free:
        return None

    candidate = sp.Integer(0)

    for coeff, (ia, ib) in zip(solution, terms):
        candidate += coeff * a**ia * b**ib

    return sp.expand(candidate)


def evaluate_law(expr: sp.Expr, av: int, bv: int) -> sp.Expr:
    return sp.expand(
        expr.subs({
            a: av,
            b: bv,
        })
    )


# ============================================================================
# LEAVE-ONE-ELL-OUT
# ============================================================================

def loo_validate(
    rows: list[dict],
    law: sp.Expr,
) -> tuple[int, int]:

    ells = sorted(
        set(row["ell"] for row in rows)
    )

    degree = sp.Poly(
        law,
        a,
        b,
    ).total_degree()

    failures = 0
    total = 0

    for held_out in ells:

        train = [
            row
            for row in rows
            if row["ell"] != held_out
        ]

        test = [
            row
            for row in rows
            if row["ell"] == held_out
        ]

        training_samples = [
            (
                row["a"],
                row["b"],
                row["second"],
            )
            for row in train
        ]

        fitted = exact_polynomial_fit(
            training_samples,
            degree,
        )

        for row in test:

            total += 1

            if fitted is None:
                failures += 1
                continue

            prediction = evaluate_law(
                fitted,
                row["a"],
                row["b"],
            )

            if sp.expand(
                prediction - row["second"]
            ) != 0:
                failures += 1

    return failures, total


# ============================================================================
# SEARCH NORMAL SECOND-COEFFICIENT LAW
# ============================================================================

def search_second_laws(
    rows: list[dict],
) -> dict[int, sp.Expr]:

    print("3. SECOND-COEFFICIENT LAW SEARCH")
    print("=" * 78)

    laws = {}

    for d in sorted(
        set(row["d"] for row in rows)
    ):

        subset = [
            row
            for row in rows
            if row["d"] == d
        ]

        # Only use nondegenerate rows for the primary law.
        subset = [
            row
            for row in subset
            if (
                not row["degree_drop"]
                and row["leading_match"]
                and row["degree"] >= 1
            )
        ]

        samples = [
            (
                row["a"],
                row["b"],
                row["second"],
            )
            for row in subset
        ]

        found = None
        found_degree = None

        for degree in POLY_DEGREES:

            candidate = exact_polynomial_fit(
                samples,
                degree,
            )

            if candidate is None:
                continue

            exact = True

            for av, bv, value in samples:

                if sp.expand(
                    evaluate_law(
                        candidate,
                        av,
                        bv,
                    ) - value
                ) != 0:
                    exact = False
                    break

            if exact:
                found = sp.factor(candidate)
                found_degree = degree
                break

        print(f"d={d}:")

        if found is None:
            print("  normal-row law = NONE")
            continue

        print(f"  degree = {found_degree}")
        print(f"  law    = {found}")

        failures, total = loo_validate(
            subset,
            found,
        )

        print(
            f"  LOO = {failures}/{total}"
        )

        laws[d] = found

    print()

    return laws


# ============================================================================
# EXTENDED WEIGHTS
# ============================================================================

def validate_extended(
    laws: dict[int, sp.Expr],
    max_ell: int,
) -> None:

    print("4. EXTENDED-WEIGHT VALIDATION")
    print("=" * 78)

    failures = 0
    total = 0

    for ell in range(5, max_ell + 1, 2):

        for k in range(1, ell, 2):

            for d, law in sorted(laws.items()):

                if k + d > ell - 1:
                    continue

                row = coefficient_record(
                    k,
                    ell,
                    d,
                )

                # Do not use exceptional cancellation rows as tests of
                # the normal law.
                if (
                    row["degree_drop"]
                    or not row["leading_match"]
                    or row["degree"] < 1
                ):
                    continue

                total += 1

                prediction = evaluate_law(
                    law,
                    row["a"],
                    row["b"],
                )

                if sp.expand(
                    prediction - row["second"]
                ) != 0:

                    failures += 1

                    if failures <= 25:
                        print(
                            f"  FAIL "
                            f"(k,ell)=({k},{ell}), "
                            f"d={d}: "
                            f"pred={prediction}, "
                            f"actual={row['second']}"
                        )

    print(
        f"extended failures = {failures}/{total}"
    )

    print(
        "STATUS = "
        + ("PASS" if failures == 0 else "FAIL")
    )

    print()


# ============================================================================
# NORMALIZED M/L SEARCH
# ============================================================================

def normalized_search(
    rows: list[dict],
) -> None:

    print("5. NORMALIZED SECOND/LEADING SEARCH")
    print("=" * 78)

    for d in sorted(
        set(row["d"] for row in rows)
    ):

        subset = [
            row
            for row in rows
            if (
                row["d"] == d
                and not row["degree_drop"]
                and row["leading_match"]
                and row["lead"] != 0
            )
        ]

        samples = []

        for row in subset:

            normalized = sp.cancel(
                row["second"] / row["lead"]
            )

            samples.append(
                (
                    row["a"],
                    row["b"],
                    normalized,
                )
            )

        found = None

        for degree in range(0, 7):

            candidate = exact_polynomial_fit(
                samples,
                degree,
            )

            if candidate is None:
                continue

            if all(
                sp.expand(
                    evaluate_law(
                        candidate,
                        av,
                        bv,
                    ) - value
                ) == 0
                for av, bv, value in samples
            ):
                found = sp.factor(candidate)
                break

        print(
            f"d={d}: M/L = "
            + (
                str(found)
                if found is not None
                else "NONE"
            )
        )

    print()


# ============================================================================
# EXCEPTIONAL CANCELLATION ANALYSIS
# ============================================================================

def cancellation_analysis(
    rows: list[dict],
) -> None:

    print("6. EXCEPTIONAL CANCELLATION ANALYSIS")
    print("=" * 78)

    exceptional = [
        row for row in rows
        if row["degree_drop"]
        or not row["leading_match"]
    ]

    print(
        f"exceptional rows = "
        f"{len(exceptional)}/{len(rows)}"
    )

    by_d = {}

    for row in exceptional:
        by_d.setdefault(row["d"], []).append(row)

    for d in sorted(by_d):

        group = by_d[d]

        print(
            f"d={d}: "
            f"{len(group)} exceptional"
        )

        for row in group[:15]:

            print(
                f"  (k,ell)=({row['k']},{row['ell']}): "
                f"D={row['degree']} "
                f"D_nom={row['nom_degree']} "
                f"L={row['lead']} "
                f"L_nom={row['nom_lead']}"
            )

    print()


# ============================================================================
# REPRESENTATIVE TABLE
# ============================================================================

def representative_table(
    rows: list[dict],
) -> None:

    print("7. REPRESENTATIVE ACTUAL COEFFICIENTS")
    print("=" * 78)

    wanted = {
        (1, 5, 1),
        (1, 5, 2),
        (1, 7, 1),
        (1, 7, 2),
        (3, 7, 1),
        (3, 9, 2),
        (5, 9, 1),
        (7, 15, 6),
    }

    for row in rows:

        key = (
            row["k"],
            row["ell"],
            row["d"],
        )

        if key not in wanted:
            continue

        print(
            f"(k,ell)=({row['k']},{row['ell']}), "
            f"d={row['d']}, j={row['j']}"
        )

        print(
            f"  c(N) = {row['poly']}"
        )

        print(
            f"  actual degree = {row['degree']}"
        )

        print(
            f"  actual lead   = {row['lead']}"
        )

        print(
            f"  actual second = {row['second']}"
        )

        print(
            f"  nominal D     = {row['nom_degree']}"
        )

        print(
            f"  nominal lead  = {row['nom_lead']}"
        )

        print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print("KAPPA EXPERIMENT 113R")
    print("SECOND-COEFFICIENT POST-BOUNDARY LAW")
    print("ACTUAL-DEGREE / CANCELLATION SAFE")
    print("EXACT BINOMIAL KERNEL -> NEWTON COEFFICIENTS")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    print("PARAMETERS")
    print(f"training max ell  = {MAX_ELL_TRAIN}")
    print(f"extended max ell  = {MAX_ELL_EXTENDED}")
    print(f"d values          = 1..{MAX_D}")
    print()

    validate_quotients(MAX_ELL_TRAIN)

    rows = build_dataset(
        MAX_ELL_TRAIN,
        MAX_D,
    )

    print(
        f"DATASET rows = {len(rows)}"
    )
    print()

    report_actual_degree(rows)

    normal_rows, exceptional_rows = split_rows(rows)

    print("NORMAL / EXCEPTIONAL SPLIT")
    print(
        f"normal rows      = {len(normal_rows)}"
    )
    print(
        f"exceptional rows = {len(exceptional_rows)}"
    )
    print()

    representative_table(rows)

    laws = search_second_laws(
        normal_rows
    )

    normalized_search(
        normal_rows
    )

    cancellation_analysis(
        rows
    )

    validate_extended(
        laws,
        MAX_ELL_EXTENDED,
    )

    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "Experiment 113R removes a false assumption from 113."
    )
    print()
    print(
        "The Experiment-112 leading law is now treated as a"
    )
    print(
        "generic/nondegenerate law rather than an unconditional one."
    )
    print()
    print(
        "Exact cancellations are recorded separately."
    )
    print()
    print(
        "The primary target remains:"
    )
    print()
    print(
        "    c(N) = L*N^D + M*N^(D-1) + ..."
    )
    print()
    print(
        "with M searched as an exact function of"
    )
    print(
        "    a=(k-1)/2, b=(ell-1)/2, d."
    )
    print()
    print(
        "The strongest positive outcome is:"
    )
    print(
        "  exact M-law + LOO + extended-weight validation."
    )
    print()
    print(
        "A second useful outcome is a classification of the"
    )
    print(
        "exceptional cancellation locus."
    )
    print()
    print(
        "No numerical fitting is used."
    )
    print(
        "No factor-pair search is performed."
    )
    print(
        "No N-only factoring claim is made."
    )
    print()
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(
            "\nInterrupted.",
            file=sys.stderr,
        )
        raise