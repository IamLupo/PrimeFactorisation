#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 102R3
PAPER DIVISOR-MOMENT BASIS / NEWTON RECONSTRUCTION

PATCHED SYMpy POLYNOMIAL DIVISION
------------------------------------------------------------------------------
The previous failure came from:

    sp.div(..., domain="ZZ[N]")

which is not accepted by the installed SymPy version.

This version performs division in S over QQ(N):

    Poly(expr, S, domain=QQ.frac_field(N))

and therefore avoids the domain-parser problem.

STRUCTURAL OBJECT
------------------------------------------------------------------------------
For each paper detector

    f_(k,l)

we first compute the exact symmetric semiprime specialization

    f_(k,l)(N,S)

then remove the universal factor

    f_(k,l) = (S+1) Q_(k,l).

The quotient is decomposed in the Newton power-sum basis

    P_j = p^j + q^j

with

    P_0 = 2
    P_1 = S
    P_j = S P_(j-1) - N P_(j-2).

MAIN QUESTION
------------------------------------------------------------------------------
Does the paper detector family span the first divisor moment

    P_1 = p+q

using low-degree coefficients in N?

If yes:

    P_1 = S

and

    x^2 - P_1 x + N = 0

recovers p,q.

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

RANDOM_SEED = 102

DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

MAX_MOMENT = 10
MAX_COEFF_N_DEGREE = 3


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

    if len(primes) < 100:
        raise RuntimeError("Prime population unexpectedly small.")

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

    P = [sp.Integer(0)] * (max_j + 1)

    P[0] = sp.Integer(2)

    if max_j >= 1:
        P[1] = S

    for j in range(2, max_j + 1):
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
        for j in range(max_j + 1)
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
# QUOTIENT POLYNOMIAL
# =============================================================================

def quotient_polynomial(
    k: int,
    ell: int,
) -> tuple[sp.Symbol, sp.Symbol, sp.Expr]:

    """
    Construct

        Q_(k,l)(N,S) = f_(k,l)(N,S)/(S+1)

    using polynomial division in S over QQ(N).

    This is the patched part of Experiment 102R3.
    """

    N, S = sp.symbols("N S")
    p, q = sp.symbols("p q")

    f = (
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )

    f = sp.expand(f)

    symmetric, remainder_sym, mapping = sp.symmetrize(
        f,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder_sym) != 0:
        raise ArithmeticError(
            f"Symmetric reduction failed for f_({k},{ell})."
        )

    mapping_dict = dict(mapping)

    if len(mapping_dict) != 2:
        raise ArithmeticError(
            "Unexpected symmetric-variable mapping."
        )

    formal_symbols = list(mapping_dict.keys())

    s1, s2 = formal_symbols

    symmetric = sp.expand(
        symmetric.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )

    # -------------------------------------------------------------------------
    # PATCH:
    #
    # Perform division as a polynomial in S over QQ(N).
    # -------------------------------------------------------------------------

    coeff_domain = sp.QQ.frac_field(N)

    numerator_poly = sp.Poly(
        symmetric,
        S,
        domain=coeff_domain,
    )

    divisor_poly = sp.Poly(
        S + 1,
        S,
        domain=coeff_domain,
    )

    quotient_poly, remainder_poly = (
        numerator_poly.div(divisor_poly)
    )

    remainder = sp.expand(
        remainder_poly.as_expr()
    )

    if remainder != 0:
        raise ArithmeticError(
            f"f_({k},{ell}) is not divisible by S+1; "
            f"remainder={remainder}"
        )

    quotient = sp.cancel(
        sp.expand(
            quotient_poly.as_expr()
        )
    )

    # Exact reconstruction check.
    reconstruction = sp.expand(
        quotient * (S + 1)
    )

    if sp.expand(
        reconstruction - symmetric
    ) != 0:
        raise ArithmeticError(
            f"Quotient reconstruction failed "
            f"for f_({k},{ell})."
        )

    return N, S, sp.factor(quotient)


# =============================================================================
# POWER-SUM DECOMPOSITION
# =============================================================================

def decompose_in_power_sums(
    Q: sp.Expr,
    N: sp.Symbol,
    S: sp.Symbol,
    max_moment: int,
) -> dict[int, sp.Expr]:

    """
    Express Q(N,S) uniquely as

        C0(N) + C1(N)P1 + ... + Cd(N)Pd.

    Since P_d is monic in S of degree d, the reduction is triangular.
    """

    polynomial = sp.Poly(
        sp.expand(Q),
        S,
        domain=sp.QQ.frac_field(N),
    )

    if polynomial.is_zero:
        return {}

    degree = polynomial.degree()

    if degree > max_moment:
        raise ArithmeticError(
            f"S-degree {degree} exceeds MAX_MOMENT={max_moment}."
        )

    P = power_sum_polynomials(
        N,
        S,
        max_moment,
    )

    remainder = sp.expand(Q)
    coefficients: dict[int, sp.Expr] = {}

    for d in range(degree, 0, -1):

        current = sp.Poly(
            sp.expand(remainder),
            S,
            domain=sp.QQ.frac_field(N),
        )

        current_degree = current.degree()

        if current_degree < d:
            continue

        if current_degree != d:
            d = current_degree

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

    remainder_poly = sp.Poly(
        sp.expand(remainder),
        S,
        domain=sp.QQ.frac_field(N),
    )

    if remainder_poly.degree() > 0:
        raise ArithmeticError(
            "Power-sum triangular reduction left "
            "a positive-degree S remainder."
        )

    constant = sp.factor(
        sp.expand(
            remainder_poly.as_expr()
        )
    )

    if constant != 0:
        coefficients[0] = constant

    return {
        j: sp.factor(c)
        for j, c in coefficients.items()
        if sp.expand(c) != 0
    }


# =============================================================================
# DECOMPOSITION VALIDATION
# =============================================================================

def validate_decomposition(
    k: int,
    ell: int,
    Q: sp.Expr,
    decomposition: dict[int, sp.Expr],
    targets: list[Target],
    N: sp.Symbol,
    S: sp.Symbol,
) -> int:

    failures = 0

    highest_moment = max(
        decomposition.keys(),
        default=0,
    )

    for t in targets:

        expected_Q = int(
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
            highest_moment,
        )

        reconstructed = 0

        for j, coefficient in decomposition.items():

            coeff_value = int(
                sp.expand(
                    coefficient
                ).subs(
                    N,
                    t.n,
                )
            )

            if j == 0:
                reconstructed += coeff_value
            else:
                reconstructed += (
                    coeff_value
                    * moments[j]
                )

        if reconstructed != expected_Q:
            failures += 1

        # Independent direct detector check.
        oracle_f = detector_direct(
            k,
            ell,
            t.p,
            t.q,
        )

        if oracle_f != expected_Q * (t.s + 1):
            failures += 1

    return failures


# =============================================================================
# MOMENT MATRIX
# =============================================================================

def build_moment_matrix(
    decompositions: dict[
        tuple[int, int],
        dict[int, sp.Expr],
    ],
) -> tuple[sp.Symbol, sp.Matrix]:

    N = sp.Symbol("N")

    rows = []

    for detector in DETECTORS:

        decomposition = (
            decompositions[
                detector
            ]
        )

        rows.append(
            [
                sp.expand(
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

    return N, sp.Matrix(rows)


# =============================================================================
# RANK
# =============================================================================

def stacked_sample_rank(
    matrix: sp.Matrix,
    N: sp.Symbol,
    targets: list[Target],
) -> int:

    blocks = [
        matrix.subs(
            N,
            t.n,
        )
        for t in targets
    ]

    if not blocks:
        return 0

    stacked = sp.Matrix.vstack(
        *blocks
    )

    return int(
        stacked.rank()
    )


# =============================================================================
# SYMBOLIC MOMENT RECOVERY
# =============================================================================

def search_symbolic_moment_recovery(
    target_moment: int,
    max_n_degree: int,
) -> tuple[bool, sp.Expr | None]:

    """
    Search for

        P_target = sum_i A_i(N) Q_i(N,S)

    where deg A_i <= max_n_degree.

    A result counts only if it is an exact symbolic identity.
    """

    N, S = sp.symbols("N S")

    P = power_sum_polynomials(
        N,
        S,
        target_moment,
    )[target_moment]

    quotient_data = []

    for k, ell in DETECTORS:

        _, _, Q = quotient_polynomial(
            k,
            ell,
        )

        quotient_data.append(
            sp.expand(Q)
        )

    unknowns = []

    rhs = sp.Integer(0)

    for i, Q in enumerate(
        quotient_data
    ):

        A_i = sp.Integer(0)

        for d in range(
            max_n_degree + 1
        ):

            u = sp.Symbol(
                f"c_{i}_{d}"
            )

            unknowns.append(u)

            A_i += (
                u
                * N ** d
            )

        rhs += (
            A_i * Q
        )

    polynomial = sp.Poly(
        sp.expand(P - rhs),
        N,
        S,
    )

    equations = polynomial.coeffs()

    if not equations:
        return False, None

    try:
        A, b = sp.linear_eq_to_matrix(
            equations,
            unknowns,
        )

        solutions = list(
            sp.linsolve(
                (
                    A,
                    b,
                )
            )
        )

    except Exception:
        return False, None

    if not solutions:
        return False, None

    solution = solutions[0]

    free_symbols = set()

    for value in solution:
        free_symbols.update(
            value.free_symbols.intersection(
                set(unknowns)
            )
        )

    if free_symbols:
        return False, None

    if all(
        sp.expand(value) == 0
        for value in solution
    ):
        return False, None

    reconstructed = sp.Integer(0)

    index = 0

    for Q in quotient_data:

        A_i = sp.Integer(0)

        for d in range(
            max_n_degree + 1
        ):

            A_i += (
                solution[index]
                * N ** d
            )

            index += 1

        reconstructed += (
            A_i * Q
        )

    reconstructed = sp.expand(
        reconstructed
    )

    if sp.expand(
        reconstructed - P
    ) != 0:
        return False, None

    return True, sp.factor(
        reconstructed
    )


# =============================================================================
# TARGET-LEVEL NUMERICAL FIT
# =============================================================================

def exact_target_fit(
    targets: list[Target],
    target_moment: int,
    degree_n: int,
) -> tuple[int, int]:

    train = targets[:TRAIN_COUNT]
    test = targets[TRAIN_COUNT:]

    feature_count = len(DETECTORS) * (
        degree_n + 1
    )

    if len(train) < feature_count:
        return 0, 0

    A_rows = []
    b_rows = []

    for t in train:

        row = []

        for k, ell in DETECTORS:

            detector = detector_direct(
                k,
                ell,
                t.p,
                t.q,
            )

            q_value = sp.Rational(
                detector,
                t.s + 1,
            )

            for d in range(
                degree_n + 1
            ):
                row.append(
                    q_value
                    * (t.n ** d)
                )

        A_rows.append(row)

        b_rows.append(
            sp.Integer(
                t.p ** target_moment
                + t.q ** target_moment
            )
        )

    A = sp.Matrix(A_rows)
    b = sp.Matrix(b_rows)

    try:
        solutions = list(
            sp.linsolve(
                (
                    A,
                    b,
                )
            )
        )
    except Exception:
        return 0, 0

    if not solutions:
        return 0, 0

    solution = solutions[0]

    free = set()

    for value in solution:
        free.update(
            value.free_symbols
        )

    if free:
        return 0, 0

    def predict(t: Target):
        index = 0
        total = sp.Integer(0)

        for k, ell in DETECTORS:

            detector = detector_direct(
                k,
                ell,
                t.p,
                t.q,
            )

            q_value = sp.Rational(
                detector,
                t.s + 1,
            )

            for d in range(
                degree_n + 1
            ):

                total += (
                    solution[index]
                    * q_value
                    * (t.n ** d)
                )

                index += 1

        return sp.expand(total)

    train_ok = 0
    test_ok = 0

    for t in train:

        expected = (
            t.p ** target_moment
            + t.q ** target_moment
        )

        if predict(t) == expected:
            train_ok += 1

    for t in test:

        expected = (
            t.p ** target_moment
            + t.q ** target_moment
        )

        if predict(t) == expected:
            test_ok += 1

    return train_ok, test_ok


# =============================================================================
# NEWTON FACTORIZATION
# =============================================================================

def newton_factor_from_p1(
    n: int,
    p1: int,
) -> tuple[int, int]:

    discriminant = (
        p1 * p1
        - 4 * n
    )

    if discriminant < 0:
        raise ArithmeticError(
            "Negative discriminant."
        )

    root = isqrt(
        discriminant
    )

    if root * root != discriminant:
        raise ArithmeticError(
            "Discriminant is not a square."
        )

    if (p1 - root) % 2 != 0:
        raise ArithmeticError(
            "Parity failure."
        )

    p = (p1 - root) // 2
    q = (p1 + root) // 2

    if p * q != n:
        raise ArithmeticError(
            "Newton reconstruction failed."
        )

    return min(p, q), max(p, q)


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 102R3")
    print("PATCHED PAPER DIVISOR-MOMENT BASIS / NEWTON RECONSTRUCTION")
    print("SAFE SYMPY DOMAIN FOR POLYNOMIAL DIVISION")
    print("STRICT TARGET HOLDOUT")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # Prime population
    # -------------------------------------------------------------------------

    t0 = time.perf_counter()

    # IMPORTANT FIX:
    # The sieve takes both lo and hi.
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
        f"{time.perf_counter() - t0:.6f}s"
    )

    if not primes:
        raise RuntimeError(
            "Prime population is empty."
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

    # -------------------------------------------------------------------------
    # Quotient construction and decomposition
    # -------------------------------------------------------------------------

    print()
    print("3. PAPER QUOTIENT -> POWER-SUM BASIS")
    print("-" * 78)

    decompositions = {}

    for k, ell in DETECTORS:

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
            k,
            ell,
            Q,
            decomposition,
            targets,
            N,
            S,
        )

        decompositions[
            (k, ell)
        ] = decomposition

        print()
        print(
            f"Q_({k},{ell}) = "
            f"{sp.factor(Q)}"
        )

        print(
            "  moments:"
        )

        for j in sorted(
            decomposition.keys()
        ):
            print(
                f"    P{j}: "
                f"{sp.factor(decomposition[j])}"
            )

        print(
            f"  validation failures = "
            f"{failures}/{len(targets)}"
        )

        if failures:
            raise ArithmeticError(
                f"Moment decomposition failed "
                f"for Q_({k},{ell})."
            )

    print()
    print("status = PASS")

    # -------------------------------------------------------------------------
    # Moment matrix
    # -------------------------------------------------------------------------

    print()
    print("4. DETECTOR / MOMENT MATRIX")
    print("-" * 78)

    N, matrix = build_moment_matrix(
        decompositions
    )

    print(
        "columns = P1 ... P{}".format(
            MAX_MOMENT
        )
    )

    for i, detector in enumerate(
        DETECTORS
    ):
        print(
            f"{detector}: "
            f"{[sp.factor(matrix[i, j]) for j in range(matrix.cols)]}"
        )

    print()
    print(
        f"training stacked rank = "
        f"{stacked_sample_rank(matrix, N, train)}"
    )

    print(
        f"test stacked rank = "
        f"{stacked_sample_rank(matrix, N, test)}"
    )

    # -------------------------------------------------------------------------
    # Symbolic recovery
    # -------------------------------------------------------------------------

    print()
    print("5. SYMBOLIC P1 ... P6 RECOVERY")
    print("-" * 78)

    for j in range(1, 7):

        print()
        print(
            f"P{j}"
        )

        found = False

        for degree in range(
            0,
            MAX_COEFF_N_DEGREE + 1,
        ):

            ok, identity = (
                search_symbolic_moment_recovery(
                    j,
                    degree,
                )
            )

            if ok:
                print(
                    f"  degree_N={degree}: FOUND"
                )
                print(
                    f"  identity = {identity}"
                )
                found = True
                break

            print(
                f"  degree_N={degree}: none"
            )

        if not found:
            print(
                "  NO LOW-COMPLEXITY "
                "SYMBOLIC RECOVERY"
            )

    # -------------------------------------------------------------------------
    # Numerical holdout fits
    # -------------------------------------------------------------------------

    print()
    print("6. HELD-OUT NUMERICAL FIT DIAGNOSTIC")
    print("-" * 78)

    for j in range(1, 5):

        print()
        print(
            f"P{j}"
        )

        found = False

        for degree in range(
            0,
            MAX_COEFF_N_DEGREE + 1,
        ):

            train_ok, test_ok = (
                exact_target_fit(
                    targets,
                    j,
                    degree,
                )
            )

            print(
                f"  degree_N={degree}: "
                f"train={train_ok}/{TRAIN_COUNT} "
                f"test={test_ok}/{len(test)}"
            )

            if (
                train_ok == TRAIN_COUNT
                and test_ok == len(test)
            ):
                found = True
                break

        if not found:
            print(
                "  no exact held-out fit"
            )

    # -------------------------------------------------------------------------
    # Newton bridge
    # -------------------------------------------------------------------------

    print()
    print("7. NEWTON BRIDGE")
    print("-" * 78)

    success = 0

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
            success += 1

    print(
        f"oracle P1 -> factor recovery = "
        f"{success}/{len(test)}"
    )

    # -------------------------------------------------------------------------
    # Final diagnostic
    # -------------------------------------------------------------------------

    print()
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The corrected computation is:

    f_(k,l) -> symmetric f(N,S)
            -> divide by (S+1)
            -> Q_(k,l)(N,S)
            -> triangular Newton moment basis.

The basis is:

    P1 = p+q
    P2 = p^2+q^2
    P3 = p^3+q^3
    ...

with

    Pj = S P(j-1) - N P(j-2).

The highest-value outcome is an exact symbolic identity

    P1 = sum_i A_i(N) Q_i.

That would establish that the paper's quotient detector family contains
the ordinary first divisor moment in a low-complexity module.

A training interpolation without symbolic identity is explicitly
classified as non-evidence.

The experiment remains an information-content experiment:
Q_(k,l) values are oracle values constructed from p,q.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - total_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 102R3 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()