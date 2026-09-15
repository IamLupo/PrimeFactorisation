#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 102R4-P
PATCHED EXACT RATIONAL MOMENT INVERSION

FIXES
------------------------------------------------------------------------------
1. SymPy Matrix.column() incompatibility:
       M.inv().column(j)
   replaced by ordinary Matrix slicing.

2. M.inv() is computed ONCE, not once per column.

3. Exact affine system is retained:
       Q - C0(N) = M(N) P

   where
       P = (P1,...,P6)^T.

4. Exact validation remains on unseen targets.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
STRICT TARGET HOLDOUT
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

RANDOM_SEED = 1024

DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

MAX_MOMENT = 6


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

def sieve_primes(lo: int, hi: int) -> list[int]:

    if hi < 2 or lo > hi:
        return []

    flags = bytearray(b"\x01") * (hi + 1)
    flags[0] = 0
    flags[1] = 0

    root = isqrt(hi)

    for p in range(2, root + 1):
        if flags[p]:

            start = p * p
            count = ((hi - start) // p) + 1

            flags[
                start:hi + 1:p
            ] = b"\x00" * count

    return [
        n
        for n in range(
            max(2, lo),
            hi + 1,
        )
        if flags[n]
    ]


# =============================================================================
# TARGET GENERATION
# =============================================================================

def generate_targets(
    primes: list[int],
    count: int,
) -> list[Target]:

    rng = random.Random(
        RANDOM_SEED
    )

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
# POWER SUM POLYNOMIALS
# =============================================================================

def power_sum_polynomials(
    N: sp.Symbol,
    S: sp.Symbol,
    max_j: int,
) -> list[sp.Expr]:

    P = [
        sp.Integer(0)
    ] * (max_j + 1)

    P[0] = sp.Integer(2)

    if max_j >= 1:
        P[1] = S

    for j in range(
        2,
        max_j + 1,
    ):
        P[j] = sp.expand(
            S * P[j - 1]
            - N * P[j - 2]
        )

    return P


def power_sums_numeric(
    p: int,
    q: int,
    max_j: int,
) -> list[int]:

    return [
        p ** j + q ** j
        for j in range(
            max_j + 1
        )
    ]


# =============================================================================
# PAPER DETECTOR
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
# SYMMETRIC QUOTIENT
# =============================================================================

def quotient_polynomial(
    k: int,
    ell: int,
) -> tuple[
    sp.Symbol,
    sp.Symbol,
    sp.Expr,
]:

    N, S = sp.symbols(
        "N S"
    )

    p, q = sp.symbols(
        "p q"
    )

    f = (
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )

    symmetric, remainder, mapping = (
        sp.symmetrize(
            sp.expand(f),
            [p, q],
            formal=True,
        )
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            f"Symmetric reduction failed "
            f"for f_({k},{ell})."
        )

    symbols = list(
        dict(mapping).keys()
    )

    if len(symbols) != 2:
        raise ArithmeticError(
            "Unexpected symmetric-variable mapping."
        )

    s1, s2 = symbols

    symmetric = sp.expand(
        symmetric.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )

    # Safe SymPy domain:
    field = sp.QQ.frac_field(N)

    numerator_poly = sp.Poly(
        symmetric,
        S,
        domain=field,
    )

    divisor_poly = sp.Poly(
        S + 1,
        S,
        domain=field,
    )

    quotient_poly, remainder_poly = (
        numerator_poly.div(
            divisor_poly
        )
    )

    remainder_expr = sp.expand(
        remainder_poly.as_expr()
    )

    if remainder_expr != 0:
        raise ArithmeticError(
            f"f_({k},{ell}) is not divisible "
            f"by S+1."
        )

    Q = sp.cancel(
        sp.expand(
            quotient_poly.as_expr()
        )
    )

    if sp.expand(
        (S + 1) * Q
        - symmetric
    ) != 0:
        raise ArithmeticError(
            f"Quotient reconstruction failed "
            f"for f_({k},{ell})."
        )

    return (
        N,
        S,
        sp.factor(Q),
    )


# =============================================================================
# POWER-SUM DECOMPOSITION
# =============================================================================

def decompose_in_power_sums(
    Q: sp.Expr,
    N: sp.Symbol,
    S: sp.Symbol,
    max_moment: int,
) -> dict[int, sp.Expr]:

    field = sp.QQ.frac_field(N)

    poly = sp.Poly(
        sp.expand(Q),
        S,
        domain=field,
    )

    if poly.is_zero:
        return {}

    degree = poly.degree()

    if degree > max_moment:
        raise ArithmeticError(
            f"S-degree {degree} exceeds "
            f"MAX_MOMENT={max_moment}."
        )

    P = power_sum_polynomials(
        N,
        S,
        max_moment,
    )

    remainder = sp.expand(Q)

    coefficients: dict[
        int,
        sp.Expr,
    ] = {}

    for d in range(
        degree,
        0,
        -1,
    ):

        current = sp.Poly(
            sp.expand(remainder),
            S,
            domain=field,
        )

        current_degree = current.degree()

        if current_degree < d:
            continue

        lead = sp.expand(
            current.coeff_monomial(
                S ** d
            )
        )

        coefficients[d] = (
            coefficients.get(
                d,
                sp.Integer(0),
            )
            + lead
        )

        remainder = sp.expand(
            remainder
            - lead * P[d]
        )

    remainder = sp.expand(
        remainder
    )

    if sp.degree(
        remainder,
        S,
    ) > 0:

        raise ArithmeticError(
            "Power-sum decomposition "
            "left an S-dependent remainder."
        )

    if remainder != 0:
        coefficients[0] = sp.factor(
            remainder
        )

    return {
        j: sp.factor(c)
        for j, c in coefficients.items()
        if sp.expand(c) != 0
    }


# =============================================================================
# VALIDATION
# =============================================================================

def validate_decomposition(
    Q: sp.Expr,
    decomposition: dict[int, sp.Expr],
    targets: list[Target],
    N: sp.Symbol,
    S: sp.Symbol,
) -> int:

    failures = 0

    highest = max(
        decomposition.keys(),
        default=0,
    )

    for t in targets:

        expected = int(
            sp.expand(Q).subs(
                {
                    N: t.n,
                    S: t.s,
                }
            )
        )

        moments = power_sums_numeric(
            t.p,
            t.q,
            highest,
        )

        recovered = 0

        for j, coeff in decomposition.items():

            coeff_value = int(
                sp.expand(coeff).subs(
                    N,
                    t.n,
                )
            )

            if j == 0:
                recovered += coeff_value

            else:
                recovered += (
                    coeff_value
                    * moments[j]
                )

        if recovered != expected:
            failures += 1

    return failures


# =============================================================================
# AFFINE MOMENT SYSTEM
# =============================================================================

def build_affine_system(
    decompositions: dict[
        tuple[int, int],
        dict[int, sp.Expr],
    ],
) -> tuple[
    sp.Symbol,
    sp.Matrix,
    sp.Matrix,
]:

    N = sp.Symbol("N")

    M_rows = []
    C0_values = []

    for detector in DETECTORS:

        decomposition = (
            decompositions[
                detector
            ]
        )

        C0_values.append(
            sp.factor(
                decomposition.get(
                    0,
                    sp.Integer(0),
                )
            )
        )

        M_rows.append(
            [
                sp.factor(
                    decomposition.get(
                        j,
                        sp.Integer(0),
                    )
                )
                for j in range(
                    1,
                    MAX_MOMENT + 1,
                )
            ]
        )

    M = sp.Matrix(
        M_rows
    )

    C0 = sp.Matrix(
        C0_values
    )

    return (
        N,
        M,
        C0,
    )


# =============================================================================
# MATRIX INVERSION
# =============================================================================

def exact_inverse(
    M: sp.Matrix,
) -> sp.Matrix:

    """
    FIX:
    Do not use:

        M.inv().column(j)

    SymPy versions differ in support for Matrix.column().

    More importantly, M.inv() is computed only once.

    We simply use the inverse matrix directly.
    """

    inverse = M.inv()

    # Verify immediately.
    identity_error = (
        M * inverse
        - sp.eye(M.rows)
    )

    for value in identity_error:
        if sp.simplify(value) != 0:
            raise ArithmeticError(
                "M * M^-1 != I."
            )

    return inverse.applyfunc(
        lambda x: sp.factor(
            sp.cancel(x)
        )
    )


# =============================================================================
# NUMERICAL ORACLE Q VECTOR
# =============================================================================

def oracle_q_vector(
    t: Target,
) -> sp.Matrix:

    values = []

    for k, ell in DETECTORS:

        f = detector_direct(
            k,
            ell,
            t.p,
            t.q,
        )

        values.append(
            sp.Rational(
                f,
                t.s + 1,
            )
        )

    return sp.Matrix(
        values
    )


# =============================================================================
# NEWTON FACTORIZATION
# =============================================================================

def newton_factor_from_p1(
    n: int,
    p1: int,
) -> tuple[int, int]:

    D = (
        p1 * p1
        - 4 * n
    )

    if D < 0:
        raise ArithmeticError(
            "Negative discriminant."
        )

    root = isqrt(D)

    if root * root != D:
        raise ArithmeticError(
            "Discriminant is not square."
        )

    if (p1 - root) % 2 != 0:
        raise ArithmeticError(
            "Parity failure."
        )

    p = (p1 - root) // 2
    q = (p1 + root) // 2

    if p * q != n:
        raise ArithmeticError(
            "Recovered roots do not multiply to N."
        )

    return (
        min(p, q),
        max(p, q),
    )


# =============================================================================
# RATIONAL COMPLEXITY
# =============================================================================

def rational_complexity(
    expr: sp.Expr,
    N: sp.Symbol,
) -> tuple[int, int, int, int]:

    num, den = sp.fraction(
        sp.together(expr)
    )

    num_poly = sp.Poly(
        sp.expand(num),
        N,
    )

    den_poly = sp.Poly(
        sp.expand(den),
        N,
    )

    return (
        num_poly.degree(),
        den_poly.degree(),
        len(num_poly.terms()),
        len(den_poly.terms()),
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 102R4-P")
    print("PATCHED EXACT RATIONAL MOMENT INVERSION")
    print("PAPER QUOTIENT FAMILY -> P1")
    print("STRICT TARGET HOLDOUT")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # PRIME POPULATION
    # -------------------------------------------------------------------------

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
        f"{time.perf_counter() - start:.6f}s"
    )

    if not primes:
        raise RuntimeError(
            "Prime generation failed."
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
        print("... remaining targets omitted")

    train = targets[
        :TRAIN_COUNT
    ]

    test = targets[
        TRAIN_COUNT:
    ]

    print()
    print("2. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = {len(train)}"
    )
    print(
        f"test targets     = {len(test)}"
    )

    # -------------------------------------------------------------------------
    # BUILD QUOTIENTS
    # -------------------------------------------------------------------------

    print()
    print("3. QUOTIENT -> MOMENT DECOMPOSITION")
    print("-" * 78)

    decompositions = {}

    for detector in DETECTORS:

        k, ell = detector

        N, S, Q = quotient_polynomial(
            k,
            ell,
        )

        decomposition = (
            decompose_in_power_sums(
                Q,
                N,
                S,
                MAX_MOMENT,
            )
        )

        failures = validate_decomposition(
            Q,
            decomposition,
            targets,
            N,
            S,
        )

        if failures:
            raise ArithmeticError(
                f"Moment decomposition failed "
                f"for Q_({k},{ell})."
            )

        decompositions[
            detector
        ] = decomposition

        print()
        print(
            f"Q_({k},{ell}) = "
            f"{sp.factor(Q)}"
        )

        for j in sorted(
            decomposition
        ):
            print(
                f"  P{j}: "
                f"{sp.factor(decomposition[j])}"
            )

        print(
            f"  validation = "
            f"{len(targets)}/{len(targets)}"
        )

    print()
    print("status = PASS")

    # -------------------------------------------------------------------------
    # AFFINE MATRIX
    # -------------------------------------------------------------------------

    print()
    print("4. AFFINE MOMENT SYSTEM")
    print("-" * 78)

    N, M, C0 = build_affine_system(
        decompositions
    )

    print(
        "Q - C0(N) = M(N) * P"
    )

    print()
    print(
        f"M shape = "
        f"{M.rows} x {M.cols}"
    )

    for detector, row in zip(
        DETECTORS,
        M.tolist(),
    ):

        print(
            f"{detector}: "
            + str(
                [
                    sp.factor(x)
                    for x in row
                ]
            )
        )

    print()
    print("C0:")

    for detector, value in zip(
        DETECTORS,
        C0,
    ):
        print(
            f"{detector}: "
            f"{sp.factor(value)}"
        )

    # -------------------------------------------------------------------------
    # DETERMINANT
    # -------------------------------------------------------------------------

    print()
    print("5. EXACT DETERMINANT")
    print("-" * 78)

    start = time.perf_counter()

    determinant = sp.factor(
        M.det()
    )

    print(
        f"determinant time = "
        f"{time.perf_counter()-start:.6f}s"
    )

    print(
        "det(M) ="
    )

    print(
        determinant
    )

    # -------------------------------------------------------------------------
    # EXACT INVERSE
    # -------------------------------------------------------------------------

    print()
    print("6. EXACT RATIONAL INVERSE")
    print("-" * 78)

    start = time.perf_counter()

    # FIX:
    #
    # OLD:
    #     solved = M.inv().column(j)
    #
    # This failed because Matrix.column() is unavailable in this SymPy build.
    #
    # NEW:
    #     inverse = M.inv()
    #
    # computed exactly once.
    inverse = exact_inverse(
        M
    )

    print(
        f"inverse time = "
        f"{time.perf_counter()-start:.6f}s"
    )

    # -------------------------------------------------------------------------
    # P1 RESOLVENT
    # -------------------------------------------------------------------------

    print()
    print("7. P1 RATIONAL RESOLVENT")
    print("-" * 78)

    Q_symbols = sp.symbols(
        "Q13 Q15 Q17 Q35 Q37 Q57"
    )

    p1_expr = sp.Integer(0)

    for j in range(
        MAX_MOMENT
    ):

        p1_expr += (
            inverse[0, j]
            * (
                Q_symbols[j]
                - C0[j]
            )
        )

    p1_expr = sp.factor(
        sp.cancel(
            p1_expr
        )
    )

    print(
        "P1 ="
    )

    print(
        p1_expr
    )

    # -------------------------------------------------------------------------
    # P1 COEFFICIENTS
    # -------------------------------------------------------------------------

    print()
    print("8. P1 COEFFICIENT COMPLEXITY")
    print("-" * 78)

    for j, detector in enumerate(
        DETECTORS
    ):

        coeff = sp.factor(
            sp.cancel(
                inverse[0, j]
            )
        )

        stats = rational_complexity(
            coeff,
            N,
        )

        print()
        print(
            f"Q_({detector}) coefficient:"
        )
        print(
            f"  degree numerator   = {stats[0]}"
        )
        print(
            f"  degree denominator = {stats[1]}"
        )
        print(
            f"  numerator terms    = {stats[2]}"
        )
        print(
            f"  denominator terms  = {stats[3]}"
        )
        print(
            f"  coefficient        = {coeff}"
        )

    # -------------------------------------------------------------------------
    # FULL SYMBOLIC IDENTITY
    # -------------------------------------------------------------------------

    print()
    print("9. AFFINE INVERSE IDENTITY CHECK")
    print("-" * 78)

    identity_error = (
        M * inverse
        - sp.eye(MAX_MOMENT)
    )

    failures = 0

    for value in identity_error:

        if sp.simplify(value) != 0:
            failures += 1

    print(
        f"M*M^-1 failures = "
        f"{failures}/{MAX_MOMENT**2}"
    )

    if failures:
        raise ArithmeticError(
            "Exact inverse identity failed."
        )

    print("status = PASS")

    # -------------------------------------------------------------------------
    # TARGET VALIDATION
    # -------------------------------------------------------------------------

    print()
    print("10. TARGET VALIDATION")
    print("-" * 78)

    all_exact = 0
    test_exact = 0
    test_p1 = 0
    test_factor = 0
    singular = 0

    for index, t in enumerate(
        targets
    ):

        inverse_at_n = inverse.subs(
            N,
            t.n,
        )

        # Detect singular points robustly.
        singular_target = False

        for entry in inverse_at_n:

            _, den = sp.fraction(
                sp.together(entry)
            )

            if sp.simplify(den) == 0:
                singular_target = True
                break

        if singular_target:
            singular += 1
            continue

        Q = oracle_q_vector(
            t
        )

        C0_at_n = C0.subs(
            N,
            t.n,
        )

        recovered = (
            inverse_at_n
            * (
                Q
                - C0_at_n
            )
        )

        expected = sp.Matrix(
            power_sums_numeric(
                t.p,
                t.q,
                MAX_MOMENT,
            )[1:]
        )

        exact_vector = (
            recovered == expected
        )

        if exact_vector:
            all_exact += 1

            if index >= TRAIN_COUNT:
                test_exact += 1

        p1_ok = (
            sp.simplify(
                recovered[0]
                - t.s
            ) == 0
        )

        if p1_ok and index >= TRAIN_COUNT:

            test_p1 += 1

            try:
                p, q = (
                    newton_factor_from_p1(
                        t.n,
                        int(recovered[0]),
                    )
                )

                if {
                    p,
                    q,
                } == {
                    t.p,
                    t.q,
                }:
                    test_factor += 1

            except Exception:
                pass

        if index < 3:
            print(
                f"target {index+1}: "
                f"vector_exact={exact_vector} "
                f"P1_exact={p1_ok}"
            )

    print()
    print(
        f"all-target exact vectors = "
        f"{all_exact}/{TARGET_COUNT}"
    )

    print(
        f"test exact vectors = "
        f"{test_exact}/{len(test)}"
    )

    print(
        f"test P1 recovery = "
        f"{test_p1}/{len(test)}"
    )

    print(
        f"test factor recovery = "
        f"{test_factor}/{len(test)}"
    )

    print(
        f"singular targets = "
        f"{singular}/{TARGET_COUNT}"
    )

    # -------------------------------------------------------------------------
    # DETERMINANT SINGULARITY STUDY
    # -------------------------------------------------------------------------

    print()
    print("11. DETERMINANT SINGULARITY STUDY")
    print("-" * 78)

    train_singular = []
    test_singular = []

    for i, t in enumerate(
        targets
    ):

        d = sp.simplify(
            determinant.subs(
                N,
                t.n,
            )
        )

        if d == 0:

            if i < TRAIN_COUNT:
                train_singular.append(
                    t.n
                )
            else:
                test_singular.append(
                    t.n
                )

    print(
        f"train singular = "
        f"{len(train_singular)}/{len(train)}"
    )

    print(
        f"test singular = "
        f"{len(test_singular)}/{len(test)}"
    )

    # -------------------------------------------------------------------------
    # NEWTON BRIDGE
    # -------------------------------------------------------------------------

    print()
    print("12. NEWTON FACTORIZATION BRIDGE")
    print("-" * 78)

    oracle_factor_success = 0

    for t in test:

        p, q = newton_factor_from_p1(
            t.n,
            t.s,
        )

        if {
            p,
            q,
        } == {
            t.p,
            t.q,
        }:
            oracle_factor_success += 1

    print(
        f"oracle P1 -> factors = "
        f"{oracle_factor_success}/{len(test)}"
    )

    # -------------------------------------------------------------------------
    # FINAL
    # -------------------------------------------------------------------------

    print()
    print("13. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
102R3 showed that no low-degree POLYNOMIAL N-coefficient combination
recovered P1.

102R4-P now asks the natural module-theoretic question:

    Is the six-dimensional detector system invertible over Q(N)?

If yes, then:

    Q - C0(N) = M(N) P

and

    P = M(N)^(-1)(Q-C0(N)).

The first component gives an exact rational expression for

    P1 = p+q.

This is much stronger than a fitted classifier.

A strong result is:

    det(M) != 0 generically
    + manageable rational coefficients
    + exact held-out P1 recovery.

The result would still be conditional on access to the detector
quotients Q. It would establish an exact algebraic coordinate system,
not yet an N-only factorization algorithm.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter()-total_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 102R4-P COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()