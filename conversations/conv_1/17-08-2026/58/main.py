#!/usr/bin/env python3

# ==============================================================================
# KAPPA EXPERIMENT 125R4
# EXACT C/D COORDINATE CONVERSION
# ROBUST QUADRATIC-BRANCH DECOMPOSITION
#
# Q(N,S) = C(N,X) + (2S-1) D(N,X)
# X = S(S-1)
#
# IMPORTANT:
#   No quotient-ring reconstruction is used.
#   C(X) and D(X) are recovered by exact coefficient matching in S.
#
# NO FACTOR-PAIR SEARCH
# NO CSV
# NO SKLEARN
# SAFE SYMPY DOMAINS
# ==============================================================================

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional, Sequence

import sympy as sp


# ==============================================================================
# SYMBOLS
# ==============================================================================

N, S, X = sp.symbols("N S X")
R, M = sp.symbols("r m")


# ==============================================================================
# PARAMETERS
# ==============================================================================

TRAIN_MAX_ELL = 23
FORWARD_ELL = 27

MAX_RATIO_DEGREE = 2


# ==============================================================================
# HELPERS
# ==============================================================================

def canon(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.expand(expr))


def safe_degree(expr: sp.Expr, var: sp.Symbol) -> Optional[int]:
    expr = sp.expand(expr)

    if expr == 0:
        return None

    p = sp.Poly(
        expr,
        var,
        domain=sp.QQ.frac_field(N, X),
    )

    return int(p.degree())


# ==============================================================================
# SYMMETRIC POWERS
#
# U_n = p^n + q^n
# U_0 = 2
# U_1 = S
# U_n = S U_(n-1) - N U_(n-2)
# ==============================================================================

def symmetric_power_table(max_n: int):

    if max_n < 0:
        return []

    U = [sp.Integer(2)]

    if max_n == 0:
        return U

    U.append(S)

    for n in range(2, max_n + 1):
        U.append(
            sp.expand(
                S * U[n - 1]
                - N * U[n - 2]
            )
        )

    return U


def mixed_power_sum(a: int, b: int, U) -> sp.Expr:
    """
    p^a q^b + q^a p^b
    """

    if a == b:
        return sp.expand(2 * N**a)

    if a > b:
        return sp.expand(
            N**b * U[a - b]
        )

    return sp.expand(
        N**a * U[b - a]
    )


# ==============================================================================
# PAPER DETECTOR
# ==============================================================================

def detector_polynomial(k: int, ell: int) -> sp.Expr:

    U = symmetric_power_table(
        max(k, ell)
    )

    A = sp.Integer(0)

    for j in range(ell + 1):
        A += (
            sp.binomial(ell, j)
            * mixed_power_sum(k, j, U)
        )

    B = sp.Integer(0)

    for j in range(k + 1):
        B += (
            sp.binomial(k, j)
            * mixed_power_sum(ell, j, U)
        )

    F = sp.expand(A - B)

    numerator = sp.Poly(
        F,
        S,
        domain=sp.QQ.frac_field(N),
    )

    denominator = sp.Poly(
        S + 1,
        S,
        domain=sp.QQ.frac_field(N),
    )

    Q_poly, remainder = sp.div(
        numerator,
        denominator,
    )

    rem = sp.expand(
        remainder.as_expr()
    )

    if rem != 0:
        raise ArithmeticError(
            f"Nonzero detector quotient remainder "
            f"for ({k},{ell}): {canon(rem)}"
        )

    return canon(
        Q_poly.as_expr()
    )


# ==============================================================================
# EXACT X-BASIS CONVERSION
#
# Given an S-polynomial f(S) which is known to be invariant under
#
#       S -> 1-S
#
# recover the unique polynomial F(X) such that
#
#       f(S) = F(S(S-1)).
#
# This deliberately avoids polynomial-ring quotient arithmetic.
# ==============================================================================

def polynomial_in_X_from_S(
    fS: sp.Expr,
    max_degree: Optional[int] = None,
) -> sp.Expr:
    """
    Recover F(X) from

        f(S) = F(S(S-1)).

    We deliberately avoid constructing a Poly over QQ(N) while
    temporary coefficient symbols are present.
    """

    fS = sp.expand(fS)

    if fS == 0:
        return sp.Integer(0)

    poly_f = sp.Poly(
        fS,
        S,
        domain=sp.QQ.frac_field(N),
    )

    degree_S = int(poly_f.degree())

    if degree_S % 2 != 0:
        raise ArithmeticError(
            "Invariant polynomial has odd S-degree: "
            f"degree={degree_S}"
        )

    if max_degree is None:
        max_degree = degree_S // 2

    # Unknown coefficients of F(X).
    z = sp.symbols(
        f"z0:{max_degree + 1}"
    )

    candidate = sp.expand(
        sum(
            z[j] * (S * (S - 1))**j
            for j in range(max_degree + 1)
        )
    )

    difference = sp.expand(
        candidate - fS
    )

    # IMPORTANT:
    # Do NOT use Poly(difference, S, domain=QQ(N))
    # because z0,z1,... are not elements of QQ(N).
    equations = []

    for power in range(degree_S + 1):
        coeff = sp.expand(
            difference.coeff(S, power)
        )

        if coeff != 0:
            equations.append(coeff)

    if not equations:
        return sp.Integer(0)

    solution = sp.solve(
        equations,
        z,
        dict=True,
        simplify=False,
    )

    if len(solution) != 1:
        raise ArithmeticError(
            "X-basis coefficient system is not uniquely solvable."
        )

    sol = solution[0]

    if any(
        zi not in sol
        for zi in z
    ):
        missing = [
            zi
            for zi in z
            if zi not in sol
        ]
        raise ArithmeticError(
            "Underdetermined X-basis conversion; "
            f"missing={missing}"
        )

    result = sp.expand(
        sum(
            sol[z[j]] * X**j
            for j in range(
                max_degree + 1
            )
        )
    )

    # Final exact identity check.
    check = sp.expand(
        result.subs(
            X,
            S * (S - 1),
        )
        - fS
    )

    if check != 0:
        raise ArithmeticError(
            "X-basis conversion failed exact identity check: "
            f"{canon(check)}"
        )

    return canon(result)


# ==============================================================================
# EXACT C/D DECOMPOSITION
#
# C = (Q(S)+Q(1-S))/2
#
# Q(S)-Q(1-S) is divisible by 2S-1.
#
# D = [Q(S)-Q(1-S)]/(2S-1)
#
# Both C and D are invariant under S -> 1-S and therefore belong to
# QQ(N)[S(S-1)].
# ==============================================================================

def branch_decompose(
    Q: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:

    Q = sp.expand(Q)

    Q_inv = sp.expand(
        Q.subs(
            S,
            1 - S,
        )
    )

    C_S = sp.cancel(
        sp.expand(
            (Q + Q_inv) / 2
        )
    )

    odd = sp.expand(
        Q - Q_inv
    )

    if odd == 0:
        D_S = sp.Integer(0)
    else:

        numerator = sp.Poly(
            odd,
            S,
            domain=sp.QQ.frac_field(N),
        )

        denominator = sp.Poly(
            2 * S - 1,
            S,
            domain=sp.QQ.frac_field(N),
        )

        quotient, remainder = sp.div(
            numerator,
            denominator,
        )

        rem = sp.expand(
            remainder.as_expr()
        )

        if rem != 0:
            raise ArithmeticError(
                "Odd part not divisible by "
                f"(2S-1): {canon(rem)}"
            )

        D_S = sp.expand(
            quotient.as_expr()
        )

    # Both branches must be invariant.
    if sp.expand(
        C_S.subs(S, 1 - S) - C_S
    ) != 0:
        raise ArithmeticError(
            "C branch is not invariant."
        )

    if sp.expand(
        D_S.subs(S, 1 - S) - D_S
    ) != 0:
        raise ArithmeticError(
            "D branch is not invariant."
        )

    C_deg_S = safe_degree(
        C_S,
        S,
    )

    D_deg_S = safe_degree(
        D_S,
        S,
    )

    C_X = polynomial_in_X_from_S(
        C_S,
        None
        if C_deg_S is None
        else C_deg_S // 2,
    )

    D_X = polynomial_in_X_from_S(
        D_S,
        None
        if D_deg_S is None
        else D_deg_S // 2,
    )

    return (
        canon(C_X),
        canon(D_X),
    )


# ==============================================================================
# DIRECT RECONSTRUCTION
# ==============================================================================

def reconstruct_from_CD(
    C: sp.Expr,
    D: sp.Expr,
) -> sp.Expr:

    C_S = sp.expand(
        C.subs(
            X,
            S * (S - 1),
        )
    )

    D_S = sp.expand(
        D.subs(
            X,
            S * (S - 1),
        )
    )

    return canon(
        C_S
        + (2 * S - 1)
        * D_S
    )


def validate_cd(
    Q: sp.Expr,
    C: sp.Expr,
    D: sp.Expr,
) -> bool:

    reconstructed = reconstruct_from_CD(
        C,
        D,
    )

    return sp.expand(
        reconstructed - Q
    ) == 0


# ==============================================================================
# SCALAR EDGES
# ==============================================================================

def scalar_edge(
    expr: sp.Expr,
) -> tuple[
    Optional[int],
    sp.Expr,
    sp.Expr,
]:

    expr0 = sp.expand(
        expr.subs(
            X,
            0,
        )
    )

    if expr0 == 0:
        return (
            None,
            sp.Integer(0),
            sp.Integer(0),
        )

    poly = sp.Poly(
        expr0,
        N,
        domain=sp.QQ,
    )

    degree = int(
        poly.degree()
    )

    top = sp.expand(
        poly.coeff_monomial(
            N**degree
        )
    )

    if degree == 0:
        second = sp.Integer(0)
    else:
        second = sp.expand(
            poly.coeff_monomial(
                N**(degree - 1)
            )
        )

    return (
        degree,
        canon(top),
        canon(second),
    )


# ==============================================================================
# RECORD
# ==============================================================================

@dataclass
class Record:

    k: int
    ell: int
    r: int
    m: int

    Q: sp.Expr
    C: sp.Expr
    D: sp.Expr

    C_degree: Optional[int]
    C_top: sp.Expr
    C_second: sp.Expr

    D_degree: Optional[int]
    D_top: sp.Expr
    D_second: sp.Expr


# ==============================================================================
# DATASET
# ==============================================================================

def build_records(
    max_ell: int,
) -> list[Record]:

    records = []

    for ell in range(
        3,
        max_ell + 1,
        2,
    ):

        m = (ell - 1) // 2

        for k in range(
            1,
            ell,
            2,
        ):

            r = (k - 1) // 2

            Q = detector_polynomial(
                k,
                ell,
            )

            C, D = branch_decompose(
                Q
            )

            C_degree, C_top, C_second = \
                scalar_edge(C)

            D_degree, D_top, D_second = \
                scalar_edge(D)

            records.append(
                Record(
                    k=k,
                    ell=ell,
                    r=r,
                    m=m,
                    Q=Q,
                    C=C,
                    D=D,
                    C_degree=C_degree,
                    C_top=C_top,
                    C_second=C_second,
                    D_degree=D_degree,
                    D_top=D_top,
                    D_second=D_second,
                )
            )

    return records


# ==============================================================================
# VALIDATION
# ==============================================================================

def validate_records(
    records: Sequence[Record],
) -> None:

    reconstruction_failures = 0
    involution_failures = 0

    for rec in records:

        if not validate_cd(
            rec.Q,
            rec.C,
            rec.D,
        ):

            reconstruction_failures += 1

            if reconstruction_failures <= 5:
                print()
                print(
                    "RECONSTRUCTION FAILURE"
                )
                print(
                    f"(k,ell)=({rec.k},{rec.ell})"
                )

                reconstructed = \
                    reconstruct_from_CD(
                        rec.C,
                        rec.D,
                    )

                print(
                    "  residual = "
                    f"{canon(reconstructed-rec.Q)}"
                )

        Q_inv = sp.expand(
            rec.Q.subs(
                S,
                1 - S,
            )
        )

        reconstructed_inv = sp.expand(
            reconstruct_from_CD(
                rec.C,
                rec.D,
            ).subs(
                S,
                1 - S,
            )
        )

        if sp.expand(
            reconstructed_inv - Q_inv
        ) != 0:

            involution_failures += 1

    print(
        f"reconstruction failures = "
        f"{reconstruction_failures}/{len(records)}"
    )

    print(
        f"involution failures      = "
        f"{involution_failures}/{len(records)}"
    )

    if reconstruction_failures:
        raise ArithmeticError(
            "C/D reconstruction failed."
        )

    if involution_failures:
        raise ArithmeticError(
            "C/D involution validation failed."
        )

    print(
        "STATUS = PASS"
    )


# ==============================================================================
# FINITE DIFFERENCES
# ==============================================================================

def finite_difference_order(
    values: Sequence[sp.Expr],
) -> Optional[int]:

    if len(values) < 2:
        return None

    current = [
        sp.expand(v)
        for v in values
    ]

    for order in range(
        len(current)
    ):

        if len(current) == 1:
            return order

        if all(
            sp.expand(v - current[0]) == 0
            for v in current
        ):
            return order

        current = [
            sp.expand(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

    return None


# ==============================================================================
# SMALL RATIONAL RATIO SEARCH
# ==============================================================================

def triangular_basis(
    degree: int,
):
    basis = []

    for total in range(
        degree + 1
    ):
        for a in range(
            total + 1
        ):
            b = total - a
            basis.append(
                (a, b)
            )

    return basis


def basis_poly(
    coeffs,
    basis,
):
    return sp.expand(
        sum(
            coeffs[i]
            * R**a
            * M**b
            for i, (a, b) in enumerate(
                basis
            )
        )
    )


def fit_ratio_law(
    records: Sequence[Record],
    field: str,
    direction: str,
    degree: int,
):

    values = {}

    for rec in records:

        value = getattr(
            rec,
            field,
        )

        if value is not None:
            values[
                (rec.r, rec.m)
            ] = sp.expand(value)

    basis = triangular_basis(
        degree
    )

    width = len(basis)

    pcoef = sp.symbols(
        f"p0:{width}"
    )

    qcoef = sp.symbols(
        f"q0:{width}"
    )

    P = basis_poly(
        pcoef,
        basis,
    )

    Q = basis_poly(
        qcoef,
        basis,
    )

    equations = []

    for rec in records:

        source = getattr(
            rec,
            field,
        )

        if source is None or source == 0:
            continue

        if direction == "m":
            target_key = (
                rec.r,
                rec.m + 1,
            )
        else:
            target_key = (
                rec.r + 1,
                rec.m,
            )

        if target_key not in values:
            continue

        target = values[target_key]

        eq = sp.expand(
            (
                target * Q
                - source * P
            ).subs(
                {
                    R: rec.r,
                    M: rec.m,
                }
            )
        )

        equations.append(eq)

    if not equations:
        return None

    unknowns = list(pcoef) + list(qcoef)

    A, b = sp.linear_eq_to_matrix(
        equations,
        unknowns,
    )

    ns = A.nullspace()

    if len(ns) != 1:
        return None

    vector = ns[0]

    qpart = vector[width:]

    pivot = next(
        (
            v
            for v in qpart
            if v != 0
        ),
        None,
    )

    if pivot is None:
        return None

    vector = [
        sp.simplify(
            v / pivot
        )
        for v in vector
    ]

    P_expr = basis_poly(
        vector[:width],
        basis,
    )

    Q_expr = basis_poly(
        vector[width:],
        basis,
    )

    if Q_expr == 0:
        return None

    return canon(
        sp.cancel(
            P_expr / Q_expr
        )
    )


def verify_ratio_law(
    records: Sequence[Record],
    field: str,
    direction: str,
    law: sp.Expr,
) -> int:

    values = {}

    for rec in records:

        value = getattr(
            rec,
            field,
        )

        if value is not None:
            values[
                (rec.r, rec.m)
            ] = sp.expand(value)

    failures = 0

    for rec in records:

        source = getattr(
            rec,
            field,
        )

        if source is None or source == 0:
            continue

        if direction == "m":
            target_key = (
                rec.r,
                rec.m + 1,
            )
        else:
            target_key = (
                rec.r + 1,
                rec.m,
            )

        if target_key not in values:
            continue

        predicted = sp.expand(
            source
            * law.subs(
                {
                    R: rec.r,
                    M: rec.m,
                }
            )
        )

        if sp.expand(
            predicted
            - values[target_key]
        ) != 0:

            failures += 1

    return failures


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 125R4")
    print("EXACT C/D COORDINATE CONVERSION")
    print("ROBUST QUADRATIC-BRANCH DECOMPOSITION")
    print("X = S(S-1)")
    print("DIRECT IDENTITY VALIDATION")
    print("NO FULL TENSOR SEARCH")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    # --------------------------------------------------------------------------
    # BUILD TRAINING SET
    # --------------------------------------------------------------------------

    print()
    print("1. BUILDING TRAINING DATASET")
    print("-" * 78)

    t_build = time.perf_counter()

    records = build_records(
        TRAIN_MAX_ELL
    )

    print(
        f"records = {len(records)}"
    )

    print(
        f"build time = "
        f"{time.perf_counter()-t_build:.6f}s"
    )

    # --------------------------------------------------------------------------
    # VALIDATION
    # --------------------------------------------------------------------------

    print()
    print("2. EXACT C/D VALIDATION")
    print("-" * 78)

    validate_records(
        records
    )

    # --------------------------------------------------------------------------
    # SCALAR EDGE TABLE
    # --------------------------------------------------------------------------

    print()
    print("3. SCALAR EDGE TABLE")
    print("-" * 78)

    for rec in records:

        if rec.ell > 13:
            break

        print(
            f"({rec.k},{rec.ell}) "
            f"r={rec.r} m={rec.m} "
            f"| "
            f"C deg={rec.C_degree} "
            f"L={rec.C_top} "
            f"M={rec.C_second} "
            f"| "
            f"D deg={rec.D_degree} "
            f"L={rec.D_top} "
            f"M={rec.D_second}"
        )

    # --------------------------------------------------------------------------
    # FINITE DIFFERENCES
    # --------------------------------------------------------------------------

    print()
    print("4. FINITE-DIFFERENCE TEST")
    print("-" * 78)

    for field in (
        "C_top",
        "C_second",
        "D_top",
        "D_second",
    ):

        print()
        print(field)

        max_r = max(
            rec.r
            for rec in records
        )

        for r0 in range(
            max_r + 1
        ):

            data = sorted(
                (
                    rec.m,
                    getattr(
                        rec,
                        field,
                    ),
                )
                for rec in records
                if rec.r == r0
                and getattr(
                    rec,
                    field,
                ) is not None
            )

            if len(data) < 4:
                continue

            values = [
                value
                for _, value in data
            ]

            order = finite_difference_order(
                values
            )

            print(
                f"  r={r0}: "
                f"order={order}"
            )

    # --------------------------------------------------------------------------
    # RATIO LAWS
    # --------------------------------------------------------------------------

    print()
    print("5. SMALL RATIONAL-RATIO SEARCH")
    print("-" * 78)

    found_laws = []

    for field in (
        "C_top",
        "C_second",
        "D_top",
        "D_second",
    ):

        for direction in (
            "m",
            "r",
        ):

            result = None

            for degree in range(
                MAX_RATIO_DEGREE + 1
            ):

                law = fit_ratio_law(
                    records,
                    field,
                    direction,
                    degree,
                )

                if law is None:
                    continue

                failures = verify_ratio_law(
                    records,
                    field,
                    direction,
                    law,
                )

                if failures == 0:

                    result = (
                        degree,
                        law,
                    )

                    break

            if result is None:

                print(
                    f"{field}, "
                    f"direction={direction}: "
                    "no exact law"
                )

            else:

                degree, law = result

                print()
                print(
                    f"{field}, "
                    f"direction={direction}"
                )
                print(
                    f"  degree={degree}"
                )
                print(
                    f"  law={canon(law)}"
                )
                print(
                    "  STATUS = PASS"
                )

                found_laws.append(
                    (
                        field,
                        direction,
                        degree,
                        law,
                    )
                )

    # --------------------------------------------------------------------------
    # FORWARD HOLDOUT
    # --------------------------------------------------------------------------

    print()
    print("6. FORWARD HOLDOUT")
    print("-" * 78)

    future_all = build_records(
        FORWARD_ELL
    )

    future = [
        rec
        for rec in future_all
        if rec.ell == FORWARD_ELL
    ]

    print(
        f"future ell = {FORWARD_ELL}"
    )

    print(
        f"future records = "
        f"{len(future)}"
    )

    training_values = {
        field: {
            (rec.r, rec.m):
            getattr(rec, field)
            for rec in records
            if getattr(
                rec,
                field,
            ) is not None
        }
        for field in (
            "C_top",
            "C_second",
            "D_top",
            "D_second",
        )
    }

    for (
        field,
        direction,
        degree,
        law,
    ) in found_laws:

        failures = 0
        tested = 0

        values = training_values[field]

        for rec in future:

            if direction == "m":
                source = (
                    rec.r,
                    rec.m - 1,
                )
            else:
                source = (
                    rec.r - 1,
                    rec.m,
                )

            if source not in values:
                continue

            source_value = values[source]

            if source_value is None:
                continue

            if source_value == 0:
                continue

            predicted = sp.expand(
                source_value
                * law.subs(
                    {
                        R: source[0],
                        M: source[1],
                    }
                )
            )

            actual = getattr(
                rec,
                field,
            )

            tested += 1

            if sp.expand(
                predicted - actual
            ) != 0:

                failures += 1

        print()
        print(
            f"{field}, "
            f"direction={direction}, "
            f"degree={degree}"
        )
        print(
            f"  law={canon(law)}"
        )
        print(
            f"  forward failures="
            f"{failures}/{tested}"
        )

    # --------------------------------------------------------------------------
    # FINAL
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The previous failure was caused by the "
        "S^2-S-X representation layer."
    )

    print()
    print(
        "This version does not perform quotient-ring "
        "reconstruction."
    )

    print()
    print(
        "Instead it uses the exact identities:"
    )
    print()
    print(
        "  C(S) = (Q(S)+Q(1-S))/2"
    )
    print(
        "  D(S) = (Q(S)-Q(1-S))/(2S-1)"
    )
    print()
    print(
        "and then solves exactly for C(X), D(X) "
        "with X=S(S-1)."
    )

    print()
    print(
        "The only subsequent search is for small "
        "rational recurrences of the four scalar "
        "upper-edge sequences."
    )

    print()
    print(
        "If this version passes reconstruction, "
        "I would treat the C/D representation as "
        "settled and stop spending experiments on "
        "representation bugs."
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter()-t0:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 125R4 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()