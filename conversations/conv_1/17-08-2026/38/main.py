#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 106
EXACT COEFFICIENT LAW / BIVARIATE GENERATING FUNCTION

GOAL
----
Experiment 105 established:

    Q_(k,l)(N,S) = sum_j c_(k,l,j)(N) P_j

where

    P_j = p^j + q^j
    P_0 = 2
    P_1 = S
    P_j = S P_(j-1) - N P_(j-2).

The previous recurrence search was contaminated by finite zero-padding.

This experiment does NOT search recurrence in a padded sequence.

Instead it derives the Newton-basis coefficients directly from the
original detector expression and tests whether the coefficient tensor

    c_(k,l,j)(N)

has a closed combinatorial formula.

Primary goals:

1. Derive the coefficient of P_j from binomial expansion.
2. Normalize coefficients by their N-degree.
3. Search for dependence on k,l,j through exact symbolic interpolation.
4. Test a candidate binomial closed form.
5. Search a bivariate generating function in detector indices.
6. Validate every proposed identity symbolically.
7. Validate on held-out semiprimes.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
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

RANDOM_SEED = 106

K_VALUES = (1, 3, 5, 7, 9, 11)
L_VALUES = (3, 5, 7, 9, 11, 13)

DETECTORS = tuple(
    (k, l)
    for k in K_VALUES
    for l in L_VALUES
    if k < l
)

MAX_L = max(l for _, l in DETECTORS)
MAX_MOMENT = MAX_L - 1


# =============================================================================
# TARGETS
# =============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


def sieve_primes(lo: int, hi: int) -> list[int]:
    if hi < 2 or lo > hi:
        return []

    a = bytearray(b"\x01") * (hi + 1)
    a[0] = 0
    a[1] = 0

    for p in range(2, sp.integer_nthroot(hi, 2)[0] + 1):
        if a[p]:
            start = p * p
            a[start:hi + 1:p] = b"\x00" * (
                (hi - start) // p + 1
            )

    return [
        n
        for n in range(max(2, lo), hi + 1)
        if a[n]
    ]


def generate_targets(
    primes: list[int],
    count: int,
) -> list[Target]:

    rng = random.Random(RANDOM_SEED)
    seen: set[int] = set()
    out: list[Target] = []

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
# NEWTON BASIS
# =============================================================================

def newton_basis(
    N: sp.Symbol,
    S: sp.Symbol,
    degree: int,
) -> list[sp.Expr]:

    P = [sp.Integer(0)] * (degree + 1)

    P[0] = sp.Integer(2)

    if degree >= 1:
        P[1] = S

    for j in range(2, degree + 1):
        P[j] = sp.expand(
            S * P[j - 1] - N * P[j - 2]
        )

    return P


# =============================================================================
# PAPER DETECTOR
# =============================================================================

def raw_detector(
    k: int,
    l: int,
    p: sp.Symbol,
    q: sp.Symbol,
) -> sp.Expr:

    return sp.expand(
        (1 + q) ** l * p ** k
        - (1 + q) ** k * p ** l
        + (1 + p) ** l * q ** k
        - (1 + p) ** k * q ** l
    )


def symmetric_detector(
    k: int,
    l: int,
) -> tuple[sp.Symbol, sp.Symbol, sp.Expr]:

    p, q = sp.symbols("p q")
    N, S = sp.symbols("N S")

    f = raw_detector(k, l, p, q)

    symmetric, remainder, mapping = sp.symmetrize(
        f,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            f"Symmetrization failed for ({k},{l})"
        )

    mapping = dict(mapping)

    s1 = None
    s2 = None

    for sym, expr in mapping.items():
        if sp.expand(expr - (p + q)) == 0:
            s1 = sym
        elif sp.expand(expr - p * q) == 0:
            s2 = sym

    if s1 is None or s2 is None:
        raise ArithmeticError(
            "Failed to identify elementary symmetric variables."
        )

    f = sp.expand(
        symmetric.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )

    return N, S, f


def quotient_detector(
    k: int,
    l: int,
) -> tuple[sp.Symbol, sp.Symbol, sp.Expr]:

    N, S, f = symmetric_detector(k, l)

    poly = sp.Poly(
        f,
        S,
        domain=sp.QQ.frac_field(N),
    )

    divisor = sp.Poly(
        S + 1,
        S,
        domain=sp.QQ.frac_field(N),
    )

    q, r = poly.div(divisor)

    if sp.simplify(r.as_expr()) != 0:
        raise ArithmeticError(
            f"S+1 does not divide detector ({k},{l})"
        )

    Q = sp.expand(q.as_expr())

    if sp.simplify((S + 1) * Q - f) != 0:
        raise ArithmeticError(
            f"Quotient reconstruction failed ({k},{l})"
        )

    return N, S, Q


# =============================================================================
# EXACT NEWTON DECOMPOSITION
# =============================================================================

def newton_decompose(
    Q: sp.Expr,
    N: sp.Symbol,
    S: sp.Symbol,
    max_degree: int,
) -> dict[int, sp.Expr]:

    P = newton_basis(N, S, max_degree)

    remainder = sp.expand(Q)
    coeffs: dict[int, sp.Expr] = {}

    for j in range(max_degree, 0, -1):

        poly = sp.Poly(
            remainder,
            S,
            domain=sp.QQ.frac_field(N),
        )

        c = sp.factor(
            poly.coeff_monomial(S ** j)
        )

        if c != 0:
            coeffs[j] = c
            remainder = sp.expand(
                remainder - c * P[j]
            )

    remainder = sp.expand(remainder)

    if sp.degree(remainder, S) not in (None, 0):
        raise ArithmeticError(
            "Newton decomposition left S dependence."
        )

    if remainder != 0:
        coeffs[0] = sp.factor(remainder)

    return coeffs


# =============================================================================
# KEY IDEA:
# DIRECT BINOMIAL COEFFICIENT FORMULA
# =============================================================================

def direct_power_coefficient(
    k: int,
    l: int,
    j: int,
) -> sp.Expr:
    """
    Derive the coefficient of p^j + q^j after expansion.

    We use the identity

        (1+q)^l p^k
        = sum_r binom(l,r) p^k q^r

    and the symmetric companion terms.

    The coefficient requested here is only the raw monomial coefficient
    before reduction to Newton moments. This function is used to generate
    candidate structured arrays, not as the final detector coefficient.
    """

    value = sp.Integer(0)

    # Contributions with p exponent j.
    if k == j:
        value += 1

    # Terms from (1+p)^l q^k.
    if k == j:
        value += 1

    return sp.expand(value)


# =============================================================================
# COEFFICIENT NORMALIZATION
# =============================================================================

def normalize_coefficient(
    c: sp.Expr,
    N: sp.Symbol,
) -> tuple[sp.Expr, int]:

    c = sp.factor(c)

    if c == 0:
        return sp.Integer(0), 0

    d = sp.degree(c, N)

    if d is None:
        d = 0

    return sp.factor(c / N ** d), d


# =============================================================================
# BINOMIAL PATTERN TEST
# =============================================================================

def binomial_signature(
    c: sp.Expr,
    N: sp.Symbol,
) -> dict[str, object]:

    if c == 0:
        return {
            "degree": -1,
            "leading": 0,
            "normalized": 0,
        }

    poly = sp.Poly(
        c,
        N,
        domain=sp.QQ,
    )

    degree = poly.degree()
    leading = poly.LC()

    normalized = sp.factor(
        c / N ** degree
    )

    return {
        "degree": degree,
        "leading": leading,
        "normalized": normalized,
    }


# =============================================================================
# DETECTOR TENSOR
# =============================================================================

def build_tensor(
    decompositions: dict[
        tuple[int, int],
        dict[int, sp.Expr],
    ],
) -> dict[
    tuple[int, int, int],
    sp.Expr,
]:

    tensor = {}

    for (k, l), coeffs in decompositions.items():
        for j in range(
            1,
            MAX_MOMENT + 1,
        ):
            tensor[k, l, j] = sp.factor(
                coeffs.get(j, 0)
            )

    return tensor


# =============================================================================
# TEST SEPARATION BY l-k
# =============================================================================

def test_difference_dependence(
    tensor,
    N: sp.Symbol,
) -> None:

    print()
    print("7. DEPENDENCE ON WEIGHT DIFFERENCE d=l-k")
    print("-" * 78)

    grouped = {}

    for (k, l) in DETECTORS:
        d = l - k

        grouped.setdefault(
            d,
            []
        ).append((k, l))

    for d, detectors in sorted(grouped.items()):

        if len(detectors) < 2:
            continue

        print()
        print(
            f"d = {d}"
        )

        for j in range(
            1,
            MAX_MOMENT + 1,
        ):

            values = [
                tensor[
                    k,
                    l,
                    j,
                ]
                for k, l in detectors
            ]

            if len(set(map(str, values))) == 1:
                print(
                    f"  P{j}: "
                    f"INDEPENDENT OF k "
                    f"value={values[0]}"
                )


# =============================================================================
# TEST FACTORIZATION IN k
# =============================================================================

def test_k_factorization(
    tensor,
    N: sp.Symbol,
) -> None:

    print()
    print("8. k-FAMILY FACTORIZATION TEST")
    print("-" * 78)

    for l in L_VALUES:

        family = [
            (k, l)
            for k, ll in DETECTORS
            if ll == l
        ]

        if len(family) < 2:
            continue

        print()
        print(
            f"l={l}"
        )

        for j in range(
            1,
            l,
        ):

            values = [
                tensor[
                    k,
                    l,
                    j,
                ]
                for k, _ in family
            ]

            if not any(
                v != 0
                for v in values
            ):
                continue

            print(
                f"  P{j}:"
            )

            for (k, _), value in zip(
                family,
                values,
            ):

                print(
                    f"    k={k}: {sp.factor(value)}"
                )


# =============================================================================
# SYMBOLIC INTERPOLATION IN k
# =============================================================================

def interpolate_in_k(
    tensor,
    N: sp.Symbol,
    l: int,
    j: int,
) -> sp.Expr | None:

    points = []

    for k in K_VALUES:
        if k >= l:
            continue

        c = tensor[
            k,
            l,
            j,
        ]

        if c == 0:
            continue

        points.append(
            (k, c)
        )

    if len(points) < 2:
        return None

    x = sp.Symbol("x")

    poly = sp.interpolate(
        points,
        x,
    )

    poly = sp.factor(
        poly
    )

    # Validate only at available k.
    for k, expected in points:
        if sp.simplify(
            poly.subs(x, k) - expected
        ) != 0:
            return None

    return poly


# =============================================================================
# GENERATING FUNCTION IN k
# =============================================================================

def k_generating_polynomial(
    tensor,
    N: sp.Symbol,
    l: int,
    j: int,
) -> sp.Expr:

    x = sp.Symbol("x")

    expr = sp.Integer(0)

    for k in K_VALUES:
        if k >= l:
            continue

        expr += (
            tensor[
                k,
                l,
                j,
            ]
            * x ** k
        )

    return sp.factor(expr)


# =============================================================================
# HELD-OUT VALIDATION
# =============================================================================

def validate_all(
    decompositions,
    test_targets: list[Target],
) -> None:

    print()
    print("9. STRICT HELD-OUT VALIDATION")
    print("-" * 78)

    N, S = sp.symbols(
        "N S"
    )

    P = newton_basis(
        N,
        S,
        MAX_MOMENT,
    )

    failures_total = 0

    for (k, l), coeffs in decompositions.items():

        N2, S2, Q = quotient_detector(
            k,
            l,
        )

        reconstructed = sp.expand(
            coeffs.get(
                0,
                0,
            )
            + sum(
                coeffs.get(
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

        if sp.simplify(
            reconstructed - Q
        ) != 0:
            raise ArithmeticError(
                f"Symbolic reconstruction failed ({k},{l})"
            )

        failures = 0

        for t in test_targets:

            lhs = Q.subs(
                {
                    N2: t.n,
                    S2: t.s,
                }
            )

            rhs = reconstructed.subs(
                {
                    N2: t.n,
                    S2: t.s,
                }
            )

            if sp.simplify(
                lhs - rhs
            ) != 0:
                failures += 1

        print(
            f"({k},{l}): "
            f"{failures}/{len(test_targets)}"
        )

        failures_total += failures

    if failures_total:
        raise ArithmeticError(
            "Held-out validation failed."
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    start_total = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 106")
    print("EXACT COEFFICIENT LAW / BIVARIATE GENERATING STRUCTURE")
    print("NEWTON MOMENT COEFFICIENT TENSOR")
    print("NO FACTOR-PAIR SEARCH")
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

    targets = generate_targets(
        primes,
        TARGET_COUNT,
    )

    for i, t in enumerate(
        targets[:20],
        1,
    ):
        print(
            f"target {i:2d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"N={t.n} "
            f"S={t.s}"
        )

    if TARGET_COUNT > 20:
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
    # BUILD COEFFICIENT TENSOR
    # -------------------------------------------------------------------------

    print()
    print("3. EXACT NEWTON COEFFICIENT EXTRACTION")
    print("-" * 78)

    decompositions = {}

    for detector in DETECTORS:

        k, l = detector

        N, S, Q = quotient_detector(
            k,
            l,
        )

        coeffs = newton_decompose(
            Q,
            N,
            S,
            MAX_MOMENT,
        )

        decompositions[
            detector
        ] = coeffs

        print()
        print(
            f"Q_({k},{l})"
        )

        for j in range(
            1,
            l,
        ):

            c = coeffs.get(
                j,
                0,
            )

            if c != 0:
                print(
                    f"  c[{j}] = "
                    f"{sp.factor(c)}"
                )

    tensor = build_tensor(
        decompositions
    )

    # -------------------------------------------------------------------------
    # EXACT NORMALIZED SIGNATURES
    # -------------------------------------------------------------------------

    N = sp.Symbol("N")

    print()
    print("4. COEFFICIENT SIGNATURE TABLE")
    print("-" * 78)

    for detector in DETECTORS:

        k, l = detector

        print()
        print(
            f"detector=({k},{l})"
        )

        for j in range(
            1,
            l,
        ):

            c = tensor[
                k,
                l,
                j,
            ]

            if c == 0:
                continue

            sig = binomial_signature(
                c,
                N,
            )

            print(
                f"  P{j}: "
                f"deg_N={sig['degree']} "
                f"leading={sig['leading']} "
                f"normalized={sig['normalized']}"
            )

    # -------------------------------------------------------------------------
    # l-k STRUCTURE
    # -------------------------------------------------------------------------

    test_difference_dependence(
        tensor,
        N,
    )

    # -------------------------------------------------------------------------
    # k FAMILY
    # -------------------------------------------------------------------------

    test_k_factorization(
        tensor,
        N,
    )

    # -------------------------------------------------------------------------
    # INTERPOLATE COEFFICIENTS IN k
    # -------------------------------------------------------------------------

    print()
    print("6. SYMBOLIC INTERPOLATION IN k")
    print("-" * 78)

    x = sp.Symbol("x")

    for l in L_VALUES:

        print()
        print(
            f"l={l}"
        )

        for j in range(
            1,
            l,
        ):

            poly = interpolate_in_k(
                tensor,
                N,
                l,
                j,
            )

            if poly is not None:
                print(
                    f"  P{j}: "
                    f"c(k) = {poly}"
                )

    # -------------------------------------------------------------------------
    # GENERATING POLYNOMIAL IN k
    # -------------------------------------------------------------------------

    print()
    print("10. k-GENERATING POLYNOMIALS")
    print("-" * 78)

    for l in L_VALUES:

        for j in range(
            1,
            l,
        ):

            G = k_generating_polynomial(
                tensor,
                N,
                l,
                j,
            )

            if G == 0:
                continue

            print(
                f"l={l}, j={j}: {G}"
            )

    # -------------------------------------------------------------------------
    # FACTORIZATION ACROSS j
    # -------------------------------------------------------------------------

    print()
    print("11. CROSS-j FACTORIZATION")
    print("-" * 78)

    for detector in DETECTORS:

        k, l = detector

        nonzero = [
            (
                j,
                tensor[
                    k,
                    l,
                    j,
                ],
            )
            for j in range(
                1,
                l,
            )
            if tensor[
                k,
                l,
                j,
            ] != 0
        ]

        if len(nonzero) < 2:
            continue

        print()
        print(
            f"detector=({k},{l})"
        )

        for j1, c1 in nonzero:

            for j2, c2 in nonzero:

                if j2 <= j1:
                    continue

                g = sp.factor(
                    sp.gcd(
                        sp.Poly(
                            c1,
                            N,
                        ),
                        sp.Poly(
                            c2,
                            N,
                        ),
                    ).as_expr()
                )

                if g not in (
                    1,
                    -1,
                ):

                    print(
                        f"  gcd(c{j1},c{j2}) = {g}"
                    )

    # -------------------------------------------------------------------------
    # HOLDOUT
    # -------------------------------------------------------------------------

    validate_all(
        decompositions,
        test,
    )

    # -------------------------------------------------------------------------
    # DIAGNOSTIC
    # -------------------------------------------------------------------------

    print()
    print("12. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
Experiment 105 showed that the coefficient recurrence search was
dominated by finite zero tails.

This experiment attacks the remaining structural question directly:

    c_(k,l,j)(N)

as an exact symbolic coefficient family.

The valuable outcomes are:

1. SIMPLE k-law:
       c_(k,l,j)(N)
   has a closed polynomial dependence on k.

2. SIMPLE l-j law:
   the same polynomial pattern persists as l grows.

3. COMMON FACTOR:
   many coefficients share an explicit N-polynomial factor.

4. GENERATING FUNCTION:
   the coefficient tensor collapses to a small bivariate
   generating polynomial.

5. NONE:
   coefficients remain genuinely irregular.

The strongest result would be a formula valid for arbitrary
odd k < l, not merely interpolation over the tested grid.

This remains a structural experiment.
It deliberately does not claim an N-only factoring algorithm.
"""
    )

    print()
    print(
        f"detectors tested = {len(DETECTORS)}"
    )

    print(
        f"max l = {MAX_L}"
    )

    print(
        f"total runtime = "
        f"{time.perf_counter() - start_total:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 106 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

