#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 104
QUADRATIC TRACE / NORM CLASSIFICATION
PAPER DETECTORS IN SHIFTED ROOT COORDINATES

a = p + 1
b = q + 1

tau = a + b = S + 2
nu  = a*b     = N + S + 1

SHIFTED POWER SUMS:
    U_j = a^j + b^j
    U_0 = 2
    U_1 = tau
    U_j = tau*U_(j-1) - nu*U_(j-2)

QUESTIONS:
    1. Do paper quotients simplify in (tau, nu)?
    2. Do they factor more strongly there?
    3. Do shifted-power-sum coordinates have lower complexity?
    4. Is there a common trace/norm generator structure?
    5. Does a larger detector family stabilize around a small algebra?

STRICT TARGET HOLDOUT
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
import random
import time

import sympy as sp


# =============================================================================
# PARAMETERS
# =============================================================================

P_MIN = 2_000_000
P_MAX = 4_200_000

TARGET_COUNT = 50
TRAIN_COUNT = 35

RANDOM_SEED = 104

# Original six detectors used throughout the previous experiments.
BASE_DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

# Additional low-weight members of the same general family.
EXTENDED_DETECTORS = (
    (1, 9),
    (3, 9),
    (5, 9),
    (7, 9),
)

ALL_DETECTORS = BASE_DETECTORS + EXTENDED_DETECTORS

MAX_MOMENT = 9


# =============================================================================
# TARGET
# =============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> list[int]:

    if hi < 2 or lo > hi:
        return []

    flags = bytearray(b"\x01") * (hi + 1)
    flags[0] = 0
    flags[1] = 0

    root = isqrt(hi)

    for p in range(2, root + 1):
        if flags[p]:
            start = p * p
            count = (hi - start) // p + 1
            flags[start:hi + 1:p] = b"\x00" * count

    return [
        n
        for n in range(max(2, lo), hi + 1)
        if flags[n]
    ]


# =============================================================================
# TARGET GENERATION
# =============================================================================

def generate_targets(
    primes: list[int],
    count: int,
) -> list[Target]:

    rng = random.Random(RANDOM_SEED)

    out: list[Target] = []
    seen: set[int] = set()

    while len(out) < count:

        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        out.append(
            Target(
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

    return out


# =============================================================================
# PAPER COEFFICIENT
# =============================================================================

def detector_direct(
    k: int,
    ell: int,
    p: int,
    q: int,
) -> int:

    return (
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )


# =============================================================================
# SYMBOLIC PAPER QUOTIENT
# =============================================================================

def quotient_polynomial(
    k: int,
    ell: int,
) -> tuple[
    sp.Symbol,
    sp.Symbol,
    sp.Expr,
]:

    N, S = sp.symbols("N S")
    p, q = sp.symbols("p q")

    f = (
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )

    symmetric, remainder, mapping = sp.symmetrize(
        sp.expand(f),
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            f"Symmetric reduction failed for ({k},{ell})."
        )

    mapping_dict = dict(mapping)

    if len(mapping_dict) != 2:
        raise ArithmeticError("Unexpected symmetric mapping.")

    mapping_items = list(mapping_dict.items())

    s1_symbol = None
    s2_symbol = None

    for symbol, expression in mapping_items:
        if sp.expand(expression - (p + q)) == 0:
            s1_symbol = symbol
        elif sp.expand(expression - p * q) == 0:
            s2_symbol = symbol

    if s1_symbol is None or s2_symbol is None:
        raise ArithmeticError("Could not identify elementary symmetric variables.")

    symmetric = sp.expand(
        symmetric.subs(
            {
                s1_symbol: S,
                s2_symbol: N,
            }
        )
    )

    field = sp.QQ.frac_field(N)

    numerator = sp.Poly(
        symmetric,
        S,
        domain=field,
    )

    divisor = sp.Poly(
        S + 1,
        S,
        domain=field,
    )

    quotient, rem = numerator.div(divisor)

    if sp.cancel(rem.as_expr()) != 0:
        raise ArithmeticError(
            f"Detector ({k},{ell}) is not divisible by S+1."
        )

    Q = sp.factor(
        sp.cancel(
            quotient.as_expr()
        )
    )

    check = sp.cancel(
        (S + 1) * Q - symmetric
    )

    if sp.simplify(check) != 0:
        raise ArithmeticError(
            f"Quotient reconstruction failed for ({k},{ell})."
        )

    return N, S, Q


# =============================================================================
# SHIFTED TRACE / NORM COORDINATES
# =============================================================================

def shifted_coordinates(
    Q: sp.Expr,
    N: sp.Symbol,
    S: sp.Symbol,
) -> tuple[
    sp.Symbol,
    sp.Symbol,
    sp.Expr,
]:

    tau, nu = sp.symbols("tau nu")

    # tau = S + 2
    # nu  = N + S + 1
    #
    # Hence:
    # S = tau - 2
    # N = nu - tau + 1

    transformed = sp.expand(
        Q.subs(
            {
                S: tau - 2,
                N: nu - tau + 1,
            }
        )
    )

    transformed = sp.factor(
        sp.cancel(transformed)
    )

    return tau, nu, transformed


# =============================================================================
# SHIFTED POWER SUMS U_j
# =============================================================================

def shifted_power_sum_polynomials(
    tau: sp.Symbol,
    nu: sp.Symbol,
    max_j: int,
) -> list[sp.Expr]:

    U = [
        sp.Integer(0)
    ] * (max_j + 1)

    U[0] = sp.Integer(2)

    if max_j >= 1:
        U[1] = tau

    for j in range(2, max_j + 1):
        U[j] = sp.expand(
            tau * U[j - 1]
            - nu * U[j - 2]
        )

    return U


# =============================================================================
# DECOMPOSE POLYNOMIAL IN SHIFTED POWER SUMS
# =============================================================================

def decompose_shifted_moments(
    Q_shifted: sp.Expr,
    tau: sp.Symbol,
    nu: sp.Symbol,
    max_moment: int,
) -> dict[int, sp.Expr]:

    field = sp.QQ.frac_field(nu)

    poly = sp.Poly(
        sp.expand(Q_shifted),
        tau,
        domain=field,
    )

    degree = poly.degree()

    if degree > max_moment:
        raise ArithmeticError(
            f"Shifted degree {degree} exceeds max_moment={max_moment}."
        )

    U = shifted_power_sum_polynomials(
        tau,
        nu,
        max_moment,
    )

    remainder = sp.expand(Q_shifted)
    coefficients: dict[int, sp.Expr] = {}

    for d in range(degree, 0, -1):

        current = sp.Poly(
            sp.expand(remainder),
            tau,
            domain=field,
        )

        lead = sp.factor(
            current.coeff_monomial(
                tau ** d
            )
        )

        coefficients[d] = lead

        remainder = sp.expand(
            remainder - lead * U[d]
        )

    remainder = sp.expand(remainder)

    if sp.degree(
        remainder,
        tau,
    ) > 0:

        raise ArithmeticError(
            "Residual tau dependence remained."
        )

    if remainder != 0:
        coefficients[0] = sp.factor(remainder)

    return {
        j: sp.factor(c)
        for j, c in coefficients.items()
    }


# =============================================================================
# COMPLEXITY
# =============================================================================

def polynomial_complexity(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> dict[str, int]:

    expanded = sp.expand(expr)

    terms = sp.Add.make_args(expanded)

    degrees = []

    for term in terms:
        powers = term.as_powers_dict()
        total = 0

        for var in variables:
            total += int(
                powers.get(var, 0)
            )

        degrees.append(total)

    return {
        "terms": len(terms),
        "total_degree": max(degrees) if degrees else 0,
    }


# =============================================================================
# FACTORIZATION COMPLEXITY
# =============================================================================

def factor_signature(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> list[sp.Expr]:

    factored = sp.factor(expr)

    if isinstance(factored, sp.Mul):
        return list(factored.args)

    return [factored]


# =============================================================================
# NUMERIC VALIDATION
# =============================================================================

def validate_detector(
    k: int,
    ell: int,
    N: sp.Symbol,
    S: sp.Symbol,
    Q: sp.Expr,
    targets: list[Target],
) -> int:

    failures = 0

    for t in targets:

        value = Q.subs(
            {
                N: t.n,
                S: t.s,
            }
        )

        expected = sp.Rational(
            detector_direct(
                k,
                ell,
                t.p,
                t.q,
            ),
            t.s + 1,
        )

        if sp.simplify(
            value - expected
        ) != 0:
            failures += 1

    return failures


# =============================================================================
# SHIFTED COORDINATE NUMERIC VALIDATION
# =============================================================================

def validate_shifted(
    Q_shifted: sp.Expr,
    tau: sp.Symbol,
    nu: sp.Symbol,
    detector: tuple[int, int],
    targets: list[Target],
) -> int:

    k, ell = detector
    failures = 0

    for t in targets:

        tau_value = t.p + t.q + 2
        nu_value = (
            (t.p + 1)
            * (t.q + 1)
        )

        lhs = sp.simplify(
            Q_shifted.subs(
                {
                    tau: tau_value,
                    nu: nu_value,
                }
            )
        )

        rhs = sp.Rational(
            detector_direct(
                k,
                ell,
                t.p,
                t.q,
            ),
            t.s + 1,
        )

        if sp.simplify(
            lhs - rhs
        ) != 0:
            failures += 1

    return failures


# =============================================================================
# COMMON GCD IN SHIFTED ALGEBRA
# =============================================================================

def pairwise_shifted_gcds(
    transformed: dict[
        tuple[int, int],
        sp.Expr,
    ],
    tau: sp.Symbol,
    nu: sp.Symbol,
) -> list[tuple]:

    rows = []

    items = list(
        transformed.items()
    )

    for i in range(len(items)):
        d1, q1 = items[i]

        for j in range(i + 1, len(items)):
            d2, q2 = items[j]

            g = sp.factor(
                sp.Poly(
                    sp.expand(q1),
                    tau,
                    nu,
                    domain=sp.QQ,
                ).gcd(
                    sp.Poly(
                        sp.expand(q2),
                        tau,
                        nu,
                        domain=sp.QQ,
                    )
                ).as_expr()
            )

            if g != 1:
                rows.append(
                    (
                        d1,
                        d2,
                        g,
                    )
                )

    return rows


# =============================================================================
# MAIN
# =============================================================================

def main():

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 104")
    print("QUADRATIC TRACE / NORM CLASSIFICATION")
    print("PAPER DETECTORS IN SHIFTED ROOT COORDINATES")
    print("a=p+1, b=q+1")
    print("tau=a+b=S+2")
    print("nu=a*b=N+S+1")
    print("STRICT TARGET HOLDOUT")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # =========================================================================
    # PRIME POPULATION
    # =========================================================================

    start = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)

    print(
        f"prime population = {len(primes)}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter()-start:.6f}s"
    )

    targets = generate_targets(
        primes,
        TARGET_COUNT,
    )

    for i, t in enumerate(
        targets[:24],
        start=1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if TARGET_COUNT > 24:
        print(
            "... remaining targets omitted"
        )

    train = targets[:TRAIN_COUNT]
    test = targets[TRAIN_COUNT:]

    print()
    print("2. TARGET HOLDOUT")
    print("-" * 78)

    print(
        f"training targets = {len(train)}"
    )

    print(
        f"test targets     = {len(test)}"
    )

    # =========================================================================
    # SYMBOLIC BUILD
    # =========================================================================

    N, S = sp.symbols("N S")
    tau, nu = sp.symbols("tau nu")

    transformed = {}
    shifted_moments = {}

    print()
    print("3. PAPER QUOTIENT STRUCTURE")
    print("-" * 78)

    for detector in ALL_DETECTORS:

        k, ell = detector

        N2, S2, Q = quotient_polynomial(
            k,
            ell,
        )

        failures = validate_detector(
            k,
            ell,
            N2,
            S2,
            Q,
            targets,
        )

        if failures:
            raise ArithmeticError(
                f"Detector validation failed "
                f"for {detector}: {failures}"
            )

        tau2, nu2, Qt = shifted_coordinates(
            Q,
            N2,
            S2,
        )

        shifted_failures = validate_shifted(
            Qt,
            tau2,
            nu2,
            detector,
            targets,
        )

        if shifted_failures:
            raise ArithmeticError(
                f"Shifted validation failed "
                f"for {detector}: {shifted_failures}"
            )

        transformed[detector] = Qt

        # tau degree / nu degree
        poly = sp.Poly(
            sp.expand(Qt),
            tau2,
            nu2,
            domain=sp.QQ,
        )

        degree_tau = poly.degree(tau2)
        degree_nu = poly.degree(nu2)

        comp_old = polynomial_complexity(
            Q,
            (N2, S2),
        )

        comp_new = polynomial_complexity(
            Qt,
            (tau2, nu2),
        )

        print()
        print(
            f"Q_({k},{ell})"
        )

        print(
            f"  original: "
            f"terms={comp_old['terms']} "
            f"degree={comp_old['total_degree']}"
        )

        print(
            f"  shifted:  "
            f"terms={comp_new['terms']} "
            f"degree={comp_new['total_degree']} "
            f"deg_tau={degree_tau} "
            f"deg_nu={degree_nu}"
        )

        print(
            f"  shifted expression = {Qt}"
        )

        shifted_moments[detector] = (
            decompose_shifted_moments(
                Qt,
                tau2,
                nu2,
                MAX_MOMENT,
            )
        )

    # =========================================================================
    # FACTOR STRUCTURE
    # =========================================================================

    print()
    print("4. SHIFTED TRACE / NORM FACTORIZATION")
    print("-" * 78)

    for detector in ALL_DETECTORS:

        Qt = transformed[detector]

        print()
        print(
            f"Q_({detector[0]},{detector[1]}) factors:"
        )

        factors = factor_signature(
            Qt,
            (tau, nu),
        )

        for factor in factors:
            print(
                f"  {factor}"
            )

    # =========================================================================
    # SHIFTED POWER-SUM BASIS
    # =========================================================================

    print()
    print("5. SHIFTED POWER-SUM DECOMPOSITION")
    print("-" * 78)

    for detector in ALL_DETECTORS:

        print()
        print(
            f"Q_({detector[0]},{detector[1]})"
        )

        d = shifted_moments[detector]

        for j in sorted(d):
            print(
                f"  U{j}: {d[j]}"
            )

    # =========================================================================
    # COMPLEXITY COMPARISON
    # =========================================================================

    print()
    print("6. COMPLEXITY COMPARISON")
    print("-" * 78)

    print(
        "detector | old_terms old_degree | "
        "shift_terms shift_degree"
    )

    for detector in ALL_DETECTORS:

        N2, S2, Q = quotient_polynomial(
            *detector
        )

        Qt = transformed[detector]

        old = polynomial_complexity(
            Q,
            (N2, S2),
        )

        new = polynomial_complexity(
            Qt,
            (tau, nu),
        )

        print(
            f"{detector!s:12s} "
            f"{old['terms']:9d} "
            f"{old['total_degree']:10d} | "
            f"{new['terms']:10d} "
            f"{new['total_degree']:11d}"
        )

    # =========================================================================
    # COMMON GCD
    # =========================================================================

    print()
    print("7. COMMON FACTORS IN SHIFTED ALGEBRA")
    print("-" * 78)

    gcd_rows = pairwise_shifted_gcds(
        transformed,
        tau,
        nu,
    )

    if not gcd_rows:
        print(
            "No nonconstant pairwise gcds."
        )
    else:
        for d1, d2, g in gcd_rows:
            print(
                f"{d1} / {d2}: gcd = {g}"
            )

    # =========================================================================
    # TRACE GENERATOR SEARCH
    # =========================================================================

    print()
    print("8. TRACE GENERATOR SEARCH")
    print("-" * 78)

    print(
        "Testing whether low-weight quotients can be represented "
        "using only U1,U2 plus polynomial coefficients in nu."
    )

    U = shifted_power_sum_polynomials(
        tau,
        nu,
        MAX_MOMENT,
    )

    for detector in ALL_DETECTORS:

        d = shifted_moments[detector]

        max_j = max(
            j
            for j in d
            if j > 0
        )

        high = [
            j
            for j in range(
                3,
                max_j + 1,
            )
            if j in d
            and sp.simplify(
                d[j]
            ) != 0
        ]

        if not high:
            status = "U1/U2 ONLY"
        else:
            status = (
                "requires "
                + ",".join(
                    f"U{j}"
                    for j in high
                )
            )

        print(
            f"Q_{detector}: {status}"
        )

    # =========================================================================
    # SEMIPRIME TRACE/NORM IDENTITIES
    # =========================================================================

    print()
    print("9. TRACE / NORM IDENTITIES")
    print("-" * 78)

    print(
        "tau = S + 2"
    )

    print(
        "nu  = N + S + 1"
    )

    print(
        "discriminant = tau^2 - 4*nu"
    )

    discriminant = sp.expand(
        tau**2 - 4 * nu
    )

    # Back-substitute.
    discriminant_NS = sp.factor(
        discriminant.subs(
            {
                tau: S + 2,
                nu: N + S + 1,
            }
        )
    )

    print(
        f"tau^2 - 4 nu = "
        f"{discriminant_NS}"
    )

    if sp.simplify(
        discriminant_NS
        - (
            S**2
            - 4*N
        )
    ) != 0:
        raise ArithmeticError(
            "Quadratic discriminant identity failed."
        )

    print(
        "PASS: shifted discriminant = (p-q)^2"
    )

    print()
    print(
        "norm shift:"
    )

    norm_shift = sp.factor(
        (
            N + S + 1
        )
    )

    print(
        f"(p+1)(q+1) = {norm_shift}"
    )

    # =========================================================================
    # ORACLE HOLDOUT SANITY
    # =========================================================================

    print()
    print("10. HELD-OUT STRUCTURAL SANITY")
    print("-" * 78)

    for i, t in enumerate(
        test[:10],
        start=1,
    ):

        tau_value = t.s + 2

        nu_value = (
            (t.p + 1)
            * (t.q + 1)
        )

        discr = (
            tau_value**2
            - 4 * nu_value
        )

        expected = (
            t.p - t.q
        ) ** 2

        print(
            f"target {TRAIN_COUNT+i:2d}: "
            f"tau={tau_value} "
            f"nu={nu_value} "
            f"discriminant_ok={discr == expected}"
        )

    # =========================================================================
    # FINAL DIAGNOSTIC
    # =========================================================================

    print()
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
This experiment asks whether the paper detector family becomes
structurally simpler after moving from

    (N,S)

to the shifted quadratic coordinates

    tau = (p+1)+(q+1) = S+2
    nu  = (p+1)(q+1) = N+S+1.

These are the trace and norm of the shifted conjugates p+1 and q+1.

A meaningful positive result is not merely another recovery of S.

The interesting outcomes are:

A. LOWER COMPLEXITY
   Multiple paper quotients have substantially fewer terms or
   lower effective degree in (tau,nu).

B. COMMON TRACE/NORM GENERATOR
   Many quotients reduce to a small set of shifted power sums U_j.

C. NEW FACTORIZATION STRUCTURE
   A common nontrivial factor appears in the shifted algebra.

D. STABLE GENERAL FAMILY
   The simplification persists as (k,l) increases.

The key quadratic identity is

    tau^2 - 4 nu = (p-q)^2

so the discriminant of the shifted quadratic algebra is exactly
the classical factor-producing square.

This experiment is intentionally structural. It does not claim an
N-only factoring algorithm.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter()-total_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 104 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

