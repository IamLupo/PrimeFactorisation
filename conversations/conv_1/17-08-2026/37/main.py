#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 105
PAPER DETECTOR GENERATING-FUNCTION / RECURRENCE CLASSIFICATION

GOAL
----
Explain the INTERNAL STRUCTURE of the paper detector family.

Instead of trying to recover S=p+q, this experiment studies

    Q_(k,l)(N,S)

after removing the universal (S+1) factor and expressing the quotient in
the Newton power-sum basis

    P_j = p^j + q^j.

For a semiprime:

    P_0 = 2
    P_1 = S
    P_j = S*P_(j-1) - N*P_(j-2).

The experiment then:

1. Builds a grid of low-weight paper detectors.
2. Converts each detector to the P_j basis.
3. Extracts the coefficient tensor c[k,l,j].
4. Searches exact recurrences in j.
5. Searches exact recurrences in k/l.
6. Tests binomial-transform structure.
7. Builds generating-polynomial slices.
8. Looks for a small symbolic generator.
9. Validates all structural identities independently on unseen (p,q).

IMPORTANT
---------
This is NOT a factorization experiment.

It asks whether the detector family itself has an internal recurrence /
generating-function structure that has not yet been identified.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import isqrt
import random
import time

import sympy as sp


# =============================================================================
# PARAMETERS
# =============================================================================

P_MIN = 2_000_000
P_MAX = 4_200_000

TARGET_COUNT = 32
TRAIN_COUNT = 24

RANDOM_SEED = 105

# Detector grid.
K_VALUES = (1, 3, 5, 7, 9)
L_VALUES = (3, 5, 7, 9, 11)

# Only k < l are used.
DETECTORS = tuple(
    (k, ell)
    for k in K_VALUES
    for ell in L_VALUES
    if k < ell
)

MAX_DETECTOR_WEIGHT = max(
    k + ell
    for k, ell in DETECTORS
)

# The quotient after division by S+1 has S-degree at most k+l-2.
MAX_MOMENT = MAX_DETECTOR_WEIGHT - 2


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

    flags = bytearray(
        b"\x01"
    ) * (hi + 1)

    flags[0] = 0
    flags[1] = 0

    root = isqrt(hi)

    for p in range(
        2,
        root + 1,
    ):
        if flags[p]:

            start = p * p

            count = (
                (hi - start) // p
                + 1
            )

            flags[
                start:hi + 1:p
            ] = (
                b"\x00"
            ) * count

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
# PAPER DETECTOR
# =============================================================================

def detector_symbolic(
    k: int,
    ell: int,
    p: sp.Symbol,
    q: sp.Symbol,
) -> sp.Expr:

    return sp.expand(
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )


# =============================================================================
# SYMMETRIC REDUCTION
# =============================================================================

def symmetric_detector(
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

    f = detector_symbolic(
        k,
        ell,
        p,
        q,
    )

    symmetric, remainder, mapping = (
        sp.symmetrize(
            f,
            [p, q],
            formal=True,
        )
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            f"Symmetrization failed for "
            f"({k},{ell})."
        )

    mapping_dict = dict(mapping)

    s1 = None
    s2 = None

    for symbol, expr in mapping_dict.items():

        if sp.expand(
            expr - (p + q)
        ) == 0:
            s1 = symbol

        elif sp.expand(
            expr - p * q
        ) == 0:
            s2 = symbol

    if s1 is None or s2 is None:
        raise ArithmeticError(
            "Could not identify elementary symmetric variables."
        )

    symmetric = sp.expand(
        symmetric.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )

    return (
        N,
        S,
        symmetric,
    )


# =============================================================================
# REMOVE UNIVERSAL (S+1)
# =============================================================================

def quotient_detector(
    k: int,
    ell: int,
) -> tuple[
    sp.Symbol,
    sp.Symbol,
    sp.Expr,
]:

    N, S, f = symmetric_detector(
        k,
        ell,
    )

    field = sp.QQ.frac_field(N)

    poly = sp.Poly(
        f,
        S,
        domain=field,
    )

    divisor = sp.Poly(
        S + 1,
        S,
        domain=field,
    )

    quotient, remainder = poly.div(
        divisor
    )

    if sp.simplify(
        remainder.as_expr()
    ) != 0:

        raise ArithmeticError(
            f"Detector ({k},{ell}) "
            f"is not divisible by S+1."
        )

    Q = sp.factor(
        sp.cancel(
            quotient.as_expr()
        )
    )

    reconstruction = sp.cancel(
        (S + 1) * Q - f
    )

    if sp.simplify(
        reconstruction
    ) != 0:

        raise ArithmeticError(
            f"Quotient reconstruction "
            f"failed for ({k},{ell})."
        )

    return (
        N,
        S,
        Q,
    )


# =============================================================================
# NEWTON POWER SUMS
# =============================================================================

def newton_power_sums(
    N: sp.Symbol,
    S: sp.Symbol,
    max_j: int,
) -> list[sp.Expr]:

    P = [
        sp.Integer(0)
    ] * (
        max_j + 1
    )

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


# =============================================================================
# MOMENT DECOMPOSITION
# =============================================================================

def decompose_newton(
    Q: sp.Expr,
    N: sp.Symbol,
    S: sp.Symbol,
    max_j: int,
) -> dict[int, sp.Expr]:

    P = newton_power_sums(
        N,
        S,
        max_j,
    )

    field = sp.QQ.frac_field(N)

    poly = sp.Poly(
        Q,
        S,
        domain=field,
    )

    degree = poly.degree()

    if degree > max_j:
        raise ArithmeticError(
            f"S-degree {degree} > max_j {max_j}."
        )

    remainder = sp.expand(Q)

    result: dict[int, sp.Expr] = {}

    # Since P_j has leading term S^j,
    # this is triangular elimination.
    for d in range(
        degree,
        0,
        -1,
    ):

        current = sp.Poly(
            remainder,
            S,
            domain=field,
        )

        coefficient = sp.factor(
            current.coeff_monomial(
                S ** d
            )
        )

        result[d] = coefficient

        remainder = sp.expand(
            remainder
            - coefficient * P[d]
        )

    remainder = sp.expand(
        remainder
    )

    if sp.degree(
        remainder,
        S,
    ) > 0:

        raise ArithmeticError(
            "Newton decomposition left S dependence."
        )

    if remainder != 0:

        result[0] = sp.factor(
            remainder
        )

    return {
        j: sp.factor(c)
        for j, c in result.items()
    }


# =============================================================================
# NUMERIC DETECTOR VALIDATION
# =============================================================================

def validate_quotient(
    k: int,
    ell: int,
    N: sp.Symbol,
    S: sp.Symbol,
    Q: sp.Expr,
    targets: list[Target],
) -> int:

    failures = 0

    for t in targets:

        symbolic_value = sp.expand(
            Q.subs(
                {
                    N: t.n,
                    S: t.s,
                }
            )
        )

        direct = (
            (
                (1 + t.q) ** ell * t.p ** k
                - (1 + t.q) ** k * t.p ** ell
                + (1 + t.p) ** ell * t.q ** k
                - (1 + t.p) ** k * t.q ** ell
            )
            // (t.s + 1)
        )

        if symbolic_value != direct:
            failures += 1

    return failures


# =============================================================================
# BASIC COMPLEXITY
# =============================================================================

def term_count(
    expr: sp.Expr,
) -> int:

    expanded = sp.expand(
        expr
    )

    return len(
        sp.Add.make_args(
            expanded
        )
    )


def total_degree(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> int:

    poly = sp.Poly(
        sp.expand(expr),
        *variables,
        domain=sp.QQ,
    )

    return poly.total_degree()


# =============================================================================
# COEFFICIENT ARRAY
# =============================================================================

def coefficient_tensor(
    decompositions: dict[
        tuple[int, int],
        dict[int, sp.Expr],
    ],
) -> dict[
    tuple[int, int, int],
    sp.Expr,
]:

    tensor = {}

    for detector, decomposition in decompositions.items():

        k, ell = detector

        for j in range(
            1,
            MAX_MOMENT + 1,
        ):

            tensor[
                k,
                ell,
                j,
            ] = sp.factor(
                decomposition.get(
                    j,
                    sp.Integer(0),
                )
            )

    return tensor


# =============================================================================
# RECURRENCE SOLVER
# =============================================================================

def solve_exact_recurrence(
    sequence: list[sp.Expr],
    order: int,
) -> tuple[bool, list[sp.Expr] | None]:

    if len(sequence) <= order:
        return False, None

    unknowns = sp.symbols(
        f"c0:{order}"
    )

    equations = []

    for i in range(
        order,
        len(sequence),
    ):

        expr = (
            sequence[i]
            - sum(
                unknowns[j]
                * sequence[
                    i - 1 - j
                ]
                for j in range(order)
            )
        )

        equations.append(
            sp.Eq(
                sp.expand(expr),
                0,
            )
        )

    solution = sp.solve(
        equations,
        unknowns,
        dict=True,
    )

    if not solution:
        return False, None

    sol = solution[0]

    coefficients = [
        sp.factor(
            sol.get(
                u,
                0,
            )
        )
        for u in unknowns
    ]

    # Independent verification.
    for i in range(
        order,
        len(sequence),
    ):

        lhs = sp.expand(
            sequence[i]
        )

        rhs = sp.expand(
            sum(
                coefficients[j]
                * sequence[
                    i - 1 - j
                ]
                for j in range(order)
            )
        )

        if sp.simplify(
            lhs - rhs
        ) != 0:

            return False, None

    return True, coefficients


# =============================================================================
# RECURRENCE SEARCH ALONG j
# =============================================================================

def recurrence_search_j(
    tensor: dict[
        tuple[int, int, int],
        sp.Expr,
    ],
    detector: tuple[int, int],
) -> list[
    tuple[int, list[sp.Expr]]
]:

    k, ell = detector

    sequence = [
        tensor.get(
            (
                k,
                ell,
                j,
            ),
            sp.Integer(0),
        )
        for j in range(
            1,
            MAX_MOMENT + 1,
        )
    ]

    found = []

    for order in range(
        1,
        min(
            5,
            len(sequence) // 2 + 1,
        ),
    ):

        ok, coeffs = solve_exact_recurrence(
            sequence,
            order,
        )

        if ok:
            found.append(
                (
                    order,
                    coeffs,
                )
            )

    return found


# =============================================================================
# BINOMIAL TRANSFORM TEST
# =============================================================================

def binomial_transform_candidate(
    sequence: list[sp.Expr],
    sign: int,
) -> list[sp.Expr]:

    out = []

    for n in range(
        len(sequence)
    ):

        value = sp.Integer(0)

        for r in range(
            n + 1
        ):

            value += (
                sign ** (
                    n - r
                )
                * sp.binomial(
                    n,
                    r,
                )
                * sequence[r]
            )

        out.append(
            sp.factor(value)
        )

    return out


def detect_zero_tail(
    sequence: list[sp.Expr],
) -> int | None:

    for start in range(
        len(sequence)
    ):

        if all(
            sp.simplify(
                sequence[j]
            ) == 0
            for j in range(
                start,
                len(sequence)
            )
        ):

            return start

    return None


# =============================================================================
# GENERATING POLYNOMIAL
# =============================================================================

def generating_polynomial(
    coefficients: list[sp.Expr],
    z: sp.Symbol,
) -> sp.Expr:

    return sp.factor(
        sum(
            coefficients[j - 1]
            * z ** j
            for j in range(
                1,
                len(coefficients) + 1,
            )
        )
    )


# =============================================================================
# CROSS-DETECTOR PATTERN SEARCH
# =============================================================================

def cross_detector_difference(
    tensor: dict[
        tuple[int, int, int],
        sp.Expr,
    ],
    d1: tuple[int, int],
    d2: tuple[int, int],
    max_j: int,
) -> list[sp.Expr]:

    return [
        sp.factor(
            tensor[
                d1[0],
                d1[1],
                j,
            ]
            - tensor[
                d2[0],
                d2[1],
                j,
            ]
        )
        for j in range(
            1,
            max_j + 1,
        )
    ]


# =============================================================================
# MAIN
# =============================================================================

def main():

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 105")
    print("PAPER DETECTOR GENERATING-FUNCTION / RECURRENCE CLASSIFICATION")
    print("NEWTON MOMENT COEFFICIENT ARRAY")
    print("EXACT SYMBOLIC STRUCTURE")
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
        f"prime population = "
        f"{len(primes)}"
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
        targets[:20],
        start=1,
    ):

        print(
            f"target {i:2d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"N={t.n} "
            f"S={t.s}"
        )

    if TARGET_COUNT > 20:
        print(
            "... remaining targets omitted"
        )

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
        f"training targets = "
        f"{len(train)}"
    )

    print(
        f"test targets     = "
        f"{len(test)}"
    )

    # -------------------------------------------------------------------------
    # BUILD DETECTORS
    # -------------------------------------------------------------------------

    N, S = sp.symbols(
        "N S"
    )

    decompositions = {}

    print()
    print("3. BUILDING NEWTON COEFFICIENT ARRAY")
    print("-" * 78)

    all_validation_failures = 0

    for detector in DETECTORS:

        k, ell = detector

        N2, S2, Q = quotient_detector(
            k,
            ell,
        )

        failures = validate_quotient(
            k,
            ell,
            N2,
            S2,
            Q,
            targets,
        )

        all_validation_failures += failures

        decomposition = decompose_newton(
            Q,
            N2,
            S2,
            MAX_MOMENT,
        )

        decompositions[
            detector
        ] = decomposition

        print()
        print(
            f"Q_({k},{ell}) "
            f"S-degree={sp.degree(Q, S2)} "
            f"terms={term_count(Q)}"
        )

        print(
            "  moments:",
            {
                j: decomposition[j]
                for j in sorted(
                    decomposition
                )
                if j > 0
            },
        )

        print(
            f"  validation failures = "
            f"{failures}/{len(targets)}"
        )

    if all_validation_failures != 0:

        raise ArithmeticError(
            "Detector validation failed."
        )

    print()
    print(
        "identity validation = PASS"
    )

    # -------------------------------------------------------------------------
    # TENSOR
    # -------------------------------------------------------------------------

    tensor = coefficient_tensor(
        decompositions
    )

    print()
    print("4. COEFFICIENT TENSOR")
    print("-" * 78)

    print(
        "Rows are detector (k,l). "
        "Columns are Newton moments P_j."
    )

    header = [
        "detector"
    ] + [
        f"P{j}"
        for j in range(
            1,
            MAX_MOMENT + 1,
        )
    ]

    print(
        " | ".join(
            f"{x:>10s}"
            for x in header
        )
    )

    for detector in DETECTORS:

        row = [
            str(detector)
        ]

        for j in range(
            1,
            MAX_MOMENT + 1,
        ):

            row.append(
                str(
                    tensor[
                        detector[0],
                        detector[1],
                        j,
                    ]
                )
            )

        print(
            " | ".join(
                f"{x:>10s}"
                for x in row
            )
        )

    # -------------------------------------------------------------------------
    # RECURRENCE SEARCH IN MOMENT INDEX
    # -------------------------------------------------------------------------

    print()
    print("5. RECURRENCE SEARCH IN j")
    print("-" * 78)

    any_j_recurrence = False

    for detector in DETECTORS:

        found = recurrence_search_j(
            tensor,
            detector,
        )

        if found:

            any_j_recurrence = True

            print()
            print(
                f"detector={detector}"
            )

            for order, coeffs in found:

                print(
                    f"  order={order}"
                )

                print(
                    f"  coefficients={coeffs}"
                )

        else:

            print(
                f"{detector}: "
                f"no exact low-order recurrence"
            )

    # -------------------------------------------------------------------------
    # SEQUENCE NORMALIZATION BY N
    # -------------------------------------------------------------------------

    print()
    print("6. NORMALIZED MOMENT COEFFICIENTS")
    print("-" * 78)

    print(
        """
The raw coefficients mix detector weight with powers of N.

For each detector, divide the P_j coefficient by its highest
available N-power and inspect whether a stable coefficient pattern
appears as (k,l) changes.
"""
    )

    for detector in DETECTORS:

        k, ell = detector

        print()
        print(
            f"detector={detector}"
        )

        degree = max(
            (
                sp.degree(
                    tensor[
                        k,
                        ell,
                        j,
                    ],
                    N,
                )
                for j in range(
                    1,
                    MAX_MOMENT + 1,
                )
                if tensor[
                    k,
                    ell,
                    j,
                ] != 0
            ),
            default=0,
        )

        for j in range(
            1,
            MAX_MOMENT + 1,
        ):

            c = tensor[
                k,
                ell,
                j,
            ]

            if c == 0:
                continue

            print(
                f"  P{j}: "
                f"degree_N={sp.degree(c,N)} "
                f"coefficient={c}"
            )

    # -------------------------------------------------------------------------
    # GENERATING POLYNOMIALS
    # -------------------------------------------------------------------------

    print()
    print("7. MOMENT GENERATING POLYNOMIALS")
    print("-" * 78)

    z = sp.Symbol("z")

    generating = {}

    for detector in DETECTORS:

        coefficients = [
            tensor[
                detector[0],
                detector[1],
                j,
            ]
            for j in range(
                1,
                MAX_MOMENT + 1,
            )
        ]

        G = generating_polynomial(
            coefficients,
            z,
        )

        generating[
            detector
        ] = G

        print()
        print(
            f"G_{detector}(z) ="
        )

        print(
            G
        )

    # -------------------------------------------------------------------------
    # BINOMIAL TRANSFORM TEST
    # -------------------------------------------------------------------------

    print()
    print("8. BINOMIAL-TRANSFORM TEST")
    print("-" * 78)

    for detector in DETECTORS:

        sequence = [
            tensor[
                detector[0],
                detector[1],
                j,
            ]
            for j in range(
                1,
                MAX_MOMENT + 1,
            )
        ]

        print()
        print(
            f"detector={detector}"
        )

        for sign in (
            1,
            -1,
        ):

            transformed = (
                binomial_transform_candidate(
                    sequence,
                    sign,
                )
            )

            tail = detect_zero_tail(
                transformed
            )

            if tail is not None:
                print(
                    f"  sign={sign}: "
                    f"zero-tail starts at index {tail}"
                )
            else:
                print(
                    f"  sign={sign}: "
                    f"no exact zero-tail"
                )

    # -------------------------------------------------------------------------
    # CROSS DETECTOR DIFFERENCES
    # -------------------------------------------------------------------------

    print()
    print("9. CROSS-DETECTOR DIFFERENCE STRUCTURE")
    print("-" * 78)

    candidate_pairs = [
        (
            (1, 5),
            (3, 5),
        ),
        (
            (1, 7),
            (3, 7),
        ),
        (
            (3, 5),
            (5, 7),
        ),
        (
            (1, 9),
            (3, 9),
        ),
        (
            (3, 9),
            (5, 9),
        ),
        (
            (5, 9),
            (7, 9),
        ),
    ]

    for d1, d2 in candidate_pairs:

        if (
            d1 not in DETECTORS
            or d2 not in DETECTORS
        ):
            continue

        diff = cross_detector_difference(
            tensor,
            d1,
            d2,
            MAX_MOMENT,
        )

        nonzero = [
            (j + 1, value)
            for j, value in enumerate(
                diff
            )
            if sp.simplify(
                value
            ) != 0
        ]

        print()
        print(
            f"{d1} - {d2}"
        )

        for j, value in nonzero:
            print(
                f"  P{j}: {value}"
            )

    # -------------------------------------------------------------------------
    # WEIGHT-LINE SEARCH
    # -------------------------------------------------------------------------

    print()
    print("10. WEIGHT-LINE PATTERN SEARCH")
    print("-" * 78)

    # For each ell, compare detectors (k,ell).
    for ell in L_VALUES:

        row_detectors = [
            d
            for d in DETECTORS
            if d[1] == ell
        ]

        if len(row_detectors) < 2:
            continue

        print()
        print(
            f"ell={ell}"
        )

        for d in row_detectors:

            k = d[0]

            highest = max(
                j
                for j in range(
                    1,
                    MAX_MOMENT + 1,
                )
                if tensor[
                    k,
                    ell,
                    j,
                ] != 0
            )

            print(
                f"  k={k}: "
                f"highest_moment={highest}"
            )

    # -------------------------------------------------------------------------
    # TRAIN / TEST STRUCTURAL VALIDATION
    # -------------------------------------------------------------------------

    print()
    print("11. STRUCTURAL HOLDOUT VALIDATION")
    print("-" * 78)

    # The symbolic identities are global. Here we independently evaluate
    # the Newton expansion on held-out targets.
    total_failures = 0

    for detector in DETECTORS:

        k, ell = detector

        N2, S2, Q = quotient_detector(
            k,
            ell,
        )

        decomposition = (
            decompositions[
                detector
            ]
        )

        P = newton_power_sums(
            N2,
            S2,
            MAX_MOMENT,
        )

        reconstructed = sp.expand(
            decomposition.get(
                0,
                0,
            )
            + sum(
                decomposition.get(
                    j,
                    0,
                )
                * P[j]
                for j in range(
                    1,
                    MAX_MOMENT + 1,
                )
            )
        )

        identity = sp.simplify(
            reconstructed - Q
        )

        if identity != 0:
            raise ArithmeticError(
                f"Symbolic reconstruction failed "
                f"for {detector}."
            )

        failures = 0

        for t in test:

            lhs = sp.simplify(
                Q.subs(
                    {
                        N2: t.n,
                        S2: t.s,
                    }
                )
            )

            rhs = sp.simplify(
                reconstructed.subs(
                    {
                        N2: t.n,
                        S2: t.s,
                    }
                )
            )

            if lhs != rhs:
                failures += 1

        total_failures += failures

        print(
            f"{detector}: "
            f"test failures={failures}/{len(test)}"
        )

    if total_failures != 0:
        raise ArithmeticError(
            "Held-out reconstruction failed."
        )

    # -------------------------------------------------------------------------
    # GENERAL GENERATING-FUNCTION DIAGNOSTIC
    # -------------------------------------------------------------------------

    print()
    print("12. GENERATING-FUNCTION DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The decisive question is now:

    Are the detector coefficients c_(k,l,j)(N)
    generated by a small recurrence or a standard
    triangular transform?

A strong positive result would be:

    detector family
        ->
    one compact recurrence / generating function
        ->
    all Q_(k,l).

That would explain the internal mathematics of the paper
without using the hidden factors p and q.

A negative result would mean the observed triangular moment
structure is real but does not collapse to a small recurrence
at the tested weights.
"""
    )

    print()
    print(
        f"j-recurrence found = "
        f"{any_j_recurrence}"
    )

    print(
        f"total detectors = "
        f"{len(DETECTORS)}"
    )

    print(
        f"max moment = "
        f"{MAX_MOMENT}"
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter()-total_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 105 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()
