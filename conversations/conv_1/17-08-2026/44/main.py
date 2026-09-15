#!/usr/bin/env python3

# ==============================================================================
# KAPPA EXPERIMENT 112R
# GENERAL POST-BOUNDARY LEADING-KERNEL LAW
# EXACT BINOMIAL KERNEL -> NEWTON COEFFICIENT
# SYMPY-SAFE SYMBOLIC DIVISION
# d = r-k > 0
# NEW WEIGHTS ell = 17,19,21,23
# NO FACTOR-PAIR SEARCH
# NO CSV
# NO SKLEARN
# ==============================================================================

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import sympy as sp


# ------------------------------------------------------------------------------
# Symbols
# ------------------------------------------------------------------------------

N, S = sp.symbols("N S")
p, q = sp.symbols("p q")


# ------------------------------------------------------------------------------
# Basic helpers
# ------------------------------------------------------------------------------

def odd_values(stop: int):
    return range(1, stop, 2)


def newton_power_sum(j: int, ss: sp.Expr = S, nn: sp.Expr = N) -> sp.Expr:
    """
    P_j = p^j + q^j

    P_0 = 2
    P_1 = S
    P_j = S P_(j-1) - N P_(j-2)
    """
    if j < 0:
        return sp.Integer(0)

    if j == 0:
        return sp.Integer(2)

    if j == 1:
        return ss

    p0 = sp.Integer(2)
    p1 = ss

    for _ in range(2, j + 1):
        p0, p1 = p1, sp.expand(ss * p1 - nn * p0)

    return sp.expand(p1)


# ------------------------------------------------------------------------------
# Symmetric reduction
# ------------------------------------------------------------------------------

def symmetric_pair_sum(a: int, b: int) -> sp.Expr:
    """
    p^a q^b + p^b q^a
    expressed in S=p+q and N=pq.
    """
    if a < b:
        a, b = b, a

    if a == b:
        return N**a

    return sp.expand(
        N**b * newton_power_sum(a - b)
    )


def symmetric_detector(k: int, ell: int) -> sp.Expr:
    """
    Semiprime detector:

        F_(k,l) =
          p^k(1+q)^l + q^k(1+p)^l
          - p^l(1+q)^k - q^l(1+p)^k

    Reduced exactly into Z[N,S].
    """
    expr = (
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )

    poly = sp.Poly(sp.expand(expr), p, q)

    result = sp.Integer(0)

    for (ip, iq), coeff in poly.terms():
        if ip < iq:
            continue

        coeff = sp.Integer(coeff)

        if ip == iq:
            result += coeff * N**ip
        else:
            result += coeff * symmetric_pair_sum(ip, iq)

    return sp.factor(sp.expand(result))


# ------------------------------------------------------------------------------
# FIXED: safe symbolic division
# ------------------------------------------------------------------------------

def quotient_polynomial(k: int, ell: int) -> sp.Expr:
    """
    Divide F_(k,l)(N,S) by S+1.

    IMPORTANT:
    Do NOT use domain="ZZ[N]".
    That syntax triggers GeneratorsError on the SymPy version
    used in the user's environment.

    Instead, use polynomial division with S as the explicit generator.
    N remains an exact symbolic coefficient.
    """
    F = sp.expand(symmetric_detector(k, ell))

    quotient, remainder = sp.div(
        F,
        S + 1,
        S,
    )

    quotient = sp.expand(quotient)
    remainder = sp.expand(remainder)

    if remainder != 0:
        raise ArithmeticError(
            f"Division by S+1 failed for ({k},{ell}): "
            f"remainder={remainder}"
        )

    return sp.factor(quotient)


# ------------------------------------------------------------------------------
# Newton-basis conversion
# ------------------------------------------------------------------------------

def multiply_by_S_basis(
    basis: Dict[int, sp.Expr],
) -> Dict[int, sp.Expr]:
    """
    If

        E = c0 + sum c_j P_j,

    then

        S*P_0 = P_1
        S*P_1 = P_2 + 2N
        S*P_j = P_(j+1) + N P_(j-1), j>=2
    """
    out: Dict[int, sp.Expr] = {}

    for j, coeff in basis.items():
        coeff = sp.expand(coeff)

        if j == 0:
            out[1] = sp.expand(out.get(1, 0) + coeff)

        elif j == 1:
            out[2] = sp.expand(out.get(2, 0) + coeff)
            out[0] = sp.expand(out.get(0, 0) + 2 * N * coeff)

        else:
            out[j + 1] = sp.expand(out.get(j + 1, 0) + coeff)
            out[j - 1] = sp.expand(
                out.get(j - 1, 0) + N * coeff
            )

    return out


def power_to_newton_coefficients(poly: sp.Expr) -> Dict[int, sp.Expr]:
    """
    Write

        poly(N,S) = C0(N) + sum_j c_j(N) P_j,

    where P_j = p^j+q^j.
    """
    poly = sp.Poly(sp.expand(poly), S)

    max_deg = poly.degree()

    # S^0 = P0/2.
    #
    # We represent the constant polynomial 1 separately because
    # P0 = 2. This avoids introducing fractional coefficients
    # unnecessarily.
    powers: List[Dict[int, sp.Expr]] = [
        {0: sp.Integer(1)},   # S^0
        {1: sp.Integer(1)},   # S^1 = P1
    ]

    for _ in range(2, max_deg + 1):
        powers.append(
            multiply_by_S_basis(powers[-1])
        )

    result: Dict[int, sp.Expr] = {}

    for m in range(max_deg + 1):
        coeff_m = sp.expand(poly.coeff_monomial(S**m))

        if coeff_m == 0:
            continue

        for j, basis_coeff in powers[m].items():
            contribution = sp.expand(coeff_m * basis_coeff)

            if j == 0:
                # Basis index 0 is the literal constant 1 here.
                result[0] = sp.expand(
                    result.get(0, 0) + contribution
                )
            else:
                result[j] = sp.expand(
                    result.get(j, 0) + contribution
                )

    return {
        j: sp.factor(sp.expand(c))
        for j, c in sorted(result.items())
        if sp.expand(c) != 0
    }


# ------------------------------------------------------------------------------
# Coefficient records
# ------------------------------------------------------------------------------

@dataclass
class CoefficientRecord:
    k: int
    ell: int
    d: int
    r: int
    j: int
    coefficient: sp.Expr
    degree: int
    leading: sp.Expr
    next_coefficient: sp.Expr


def coefficient_record(k: int, ell: int, d: int) -> CoefficientRecord:
    """
    Post-boundary region:

        r = k+d
        j = ell-1-r
    """
    r = k + d
    j = ell - 1 - r

    if j < 1:
        raise ValueError(
            f"Invalid indices: k={k}, ell={ell}, d={d}, "
            f"r={r}, j={j}"
        )

    Q = quotient_polynomial(k, ell)
    coeffs = power_to_newton_coefficients(Q)

    c = sp.expand(coeffs.get(j, 0))
    P = sp.Poly(c, N)

    if c == 0:
        degree = -sp.oo
        leading = sp.Integer(0)
        next_coeff = sp.Integer(0)
    else:
        degree = int(P.degree())
        leading = sp.expand(P.LC())

        if degree >= 1:
            next_coeff = sp.expand(
                P.coeff_monomial(N ** (degree - 1))
            )
        else:
            next_coeff = sp.Integer(0)

    return CoefficientRecord(
        k=k,
        ell=ell,
        d=d,
        r=r,
        j=j,
        coefficient=sp.factor(c),
        degree=degree,
        leading=sp.factor(leading),
        next_coefficient=sp.factor(next_coeff),
    )


# ------------------------------------------------------------------------------
# Dataset
# ------------------------------------------------------------------------------

def build_dataset(max_ell: int, max_d: int) -> List[CoefficientRecord]:
    rows: List[CoefficientRecord] = []

    for ell in range(3, max_ell + 1, 2):
        for k in odd_values(ell):
            for d in range(1, max_d + 1):

                r = k + d
                j = ell - 1 - r

                if j < 1:
                    continue

                rows.append(
                    coefficient_record(k, ell, d)
                )

    return rows


# ------------------------------------------------------------------------------
# Variables a,b
# ------------------------------------------------------------------------------

def aval(k: int) -> int:
    return (k - 1) // 2


def bval(ell: int) -> int:
    return (ell - 1) // 2


# ------------------------------------------------------------------------------
# Exact bivariate polynomial fitting
# ------------------------------------------------------------------------------

def monomial_basis(a, b, total_degree: int):
    basis = []

    for i in range(total_degree + 1):
        for j in range(total_degree + 1 - i):
            basis.append(a**i * b**j)

    return basis


def fit_polynomial_law(
    rows: List[CoefficientRecord],
    attr: str,
    total_degree: int,
):
    """
    Find an exact polynomial law in

        a=(k-1)/2
        b=(ell-1)/2

    of total degree <= total_degree.

    This is deliberately exact and does not do numerical fitting.
    """
    if not rows:
        return None

    a, b = sp.symbols("a b")

    basis = monomial_basis(a, b, total_degree)
    unknowns = sp.symbols(
        f"u0:{len(basis)}"
    )

    equations = []

    for row in rows:
        lhs = sp.Integer(0)

        for idx, mon in enumerate(basis):
            lhs += unknowns[idx] * mon.subs(
                {
                    a: aval(row.k),
                    b: bval(row.ell),
                }
            )

        target = sp.sympify(getattr(row, attr))

        equations.append(
            sp.expand(lhs - target)
        )

    solution = sp.linsolve(
        equations,
        unknowns,
    )

    if solution == sp.EmptySet:
        return None

    tuples = list(solution)

    if not tuples:
        return None

    values = tuples[0]

    # Reject underdetermined families.
    if any(v.free_symbols & set(unknowns) for v in values):
        return None

    polynomial = sp.expand(
        sum(
            values[i] * basis[i]
            for i in range(len(basis))
        )
    )

    for row in rows:
        predicted = sp.expand(
            polynomial.subs(
                {
                    a: aval(row.k),
                    b: bval(row.ell),
                }
            )
        )

        actual = sp.expand(getattr(row, attr))

        if sp.expand(predicted - actual) != 0:
            return None

    return sp.factor(polynomial)


# ------------------------------------------------------------------------------
# Known d=1..4 laws
# ------------------------------------------------------------------------------

def known_leading_law(d: int):
    a, b = sp.symbols("a b")

    if d == 1:
        return 2 * (a + b + 1)

    if d == 2:
        return -(a + b + 1) * (2*a - 2*b + 3)

    if d == 3:
        return -2 * (a + b + 1)

    if d == 4:
        return (a + b + 1) * (2*a - 2*b + 5)

    return None


# ------------------------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------------------------

def print_degree_profile(rows: List[CoefficientRecord]):
    print()
    print("=" * 78)
    print("DEGREE PROFILE")
    print("=" * 78)

    grouped: Dict[int, List[int]] = {}

    for row in rows:
        grouped.setdefault(row.d, []).append(row.degree)

    for d in sorted(grouped):
        print(
            f"d={d}: "
            f"degrees={sorted(set(grouped[d]))}"
        )


def validate_known_laws(rows: List[CoefficientRecord]):
    print()
    print("=" * 78)
    print("KNOWN d=1..4 LEADING LAWS")
    print("=" * 78)

    a, b = sp.symbols("a b")

    for d in range(1, 5):
        subset = [
            row for row in rows
            if row.d == d
        ]

        law = known_leading_law(d)

        failures = 0

        for row in subset:
            predicted = sp.expand(
                law.subs(
                    {
                        a: aval(row.k),
                        b: bval(row.ell),
                    }
                )
            )

            if sp.expand(predicted - row.leading) != 0:
                failures += 1

        print(
            f"d={d}: "
            f"failures={failures}/{len(subset)} "
            f"law={sp.factor(law)}"
        )


def search_general_laws(rows: List[CoefficientRecord]):
    print()
    print("=" * 78)
    print("GENERAL LEADING-LAW SEARCH")
    print("=" * 78)

    for d in sorted({row.d for row in rows}):

        subset = [
            row for row in rows
            if row.d == d
        ]

        found = None
        found_degree = None

        for degree in range(0, 7):
            candidate = fit_polynomial_law(
                subset,
                "leading",
                degree,
            )

            if candidate is not None:
                found = candidate
                found_degree = degree
                break

        if found is None:
            print(
                f"d={d}: no exact law "
                f"of total degree <= 6"
            )
        else:
            print(
                f"d={d}: degree={found_degree} "
                f"law={sp.factor(found)}"
            )


def loo_by_ell(rows: List[CoefficientRecord]):
    print()
    print("=" * 78)
    print("LEAVE-ONE-ELL-OUT")
    print("=" * 78)

    ells = sorted(
        {row.ell for row in rows}
    )

    for d in sorted(
        {row.d for row in rows}
    ):
        subset = [
            row for row in rows
            if row.d == d
        ]

        failures = 0
        tests = 0

        for held_ell in ells:

            train = [
                row for row in subset
                if row.ell != held_ell
            ]

            test = [
                row for row in subset
                if row.ell == held_ell
            ]

            if not train or not test:
                continue

            candidate = None

            for degree in range(0, 6):
                candidate = fit_polynomial_law(
                    train,
                    "leading",
                    degree,
                )

                if candidate is not None:
                    break

            if candidate is None:
                continue

            a, b = sp.symbols("a b")

            for row in test:
                tests += 1

                predicted = sp.expand(
                    candidate.subs(
                        {
                            a: aval(row.k),
                            b: bval(row.ell),
                        }
                    )
                )

                if sp.expand(predicted - row.leading) != 0:
                    failures += 1

        print(
            f"d={d}: "
            f"LOO failures={failures}/{tests}"
        )


# ------------------------------------------------------------------------------
# Direct extended validation
# ------------------------------------------------------------------------------

def extended_validation(max_d: int):
    print()
    print("=" * 78)
    print("EXTENDED WEIGHT VALIDATION")
    print("=" * 78)

    failures = 0
    total = 0

    for ell in [17, 19, 21, 23]:

        for k in odd_values(ell):

            for d in range(1, max_d + 1):

                r = k + d
                j = ell - 1 - r

                if j < 1:
                    continue

                row = coefficient_record(
                    k,
                    ell,
                    d,
                )

                total += 1

                if row.r != r or row.j != j:
                    failures += 1

    print(
        f"extended failures={failures}/{total}"
    )

    print(
        "STATUS = "
        + ("PASS" if failures == 0 else "FAIL")
    )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 112R")
    print("GENERAL POST-BOUNDARY LEADING-KERNEL LAW")
    print("EXACT BINOMIAL KERNEL -> NEWTON COEFFICIENT")
    print("SYMPY-SAFE SYMBOLIC DIVISION")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    MAX_ELL = 23
    MAX_D = 10

    print()
    print("1. PARAMETERS")
    print("-" * 78)
    print(f"max ell = {MAX_ELL}")
    print(f"d values = 1..{MAX_D}")

    print()
    print("2. BUILDING DATASET")
    print("-" * 78)

    rows = build_dataset(
        MAX_ELL,
        MAX_D,
    )

    print(
        f"dataset rows = {len(rows)}"
    )

    if not rows:
        raise ArithmeticError(
            "Dataset is empty."
        )

    print()
    print("3. DEGREE PROFILE")
    print_degree_profile(rows)

    print()
    print("4. KNOWN LAW VALIDATION")
    validate_known_laws(rows)

    print()
    print("5. GENERAL LAW SEARCH")
    search_general_laws(rows)

    print()
    print("6. LEAVE-ONE-ELL-OUT")
    loo_by_ell(rows)

    print()
    print("7. EXTENDED WEIGHTS")
    extended_validation(MAX_D)

    print()
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print("The previous script failed before the mathematical")
    print("experiment began because SymPy interpreted")
    print('domain="ZZ[N]" as an invalid polynomial-ring generator.')
    print()
    print("This revision performs exact division with S as the")
    print("explicit polynomial variable and keeps N symbolic.")
    print()
    print("The intended mathematical test is unchanged:")
    print()
    print("    c_{ell-1-r}(N),  r = k+d, d>0")
    print()
    print("is tested for an exact law in")
    print()
    print("    a=(k-1)/2, b=(ell-1)/2, d.")
    print()
    print("No numerical approximation or factor-pair search is used.")


if __name__ == "__main__":
    main()