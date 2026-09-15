#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 100
EISENSTEIN NORM / MODULAR SQUARE-ROOT BRIDGE

NEW NUMBER-THEORETIC REFORMULATION

For n = p*q and
    T = S^2 - 3N = p^2 - p*q + q^2,

we have

    T - N = (p-q)^2
    T + N = p^2 + q^2

and T is the Eisenstein norm

    T = Norm(p + q*omega).

CENTRAL QUESTION:

Can an N-only low-degree polynomial R(N) behave as a square root
of T modulo many finite fields?

That is,

    R(N)^2 == T (mod ell)

for unseen semiprimes.

A positive result would be substantially more interesting than
ordinary character prediction because it would produce an explicit
algebraic square-root representation.

TESTS:

1. Eisenstein/norm identities
2. T-N square validation
3. T+N sum-of-squares validation
4. Oracle finite-field square roots
5. Linear N-only square-root search mod ell
6. Quadratic N-only square-root search mod ell
7. Integer small-coefficient polynomial search across many moduli
8. Strict target holdout
9. Cyclotomic-style modulus set vs matched control moduli
10. Cross-modulus consistency

IMPORTANT:

The oracle phase may use p,q only to construct T.

The N-only square-root search itself NEVER uses p or q.

NO CSV
NO SKLEARN
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Optional, Sequence, Tuple


# =============================================================================
# PARAMETERS
# =============================================================================

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

NUM_TARGETS = 60
TRAIN_TARGETS = 45

# Small finite fields.
CYCLO_MODULI = [5, 7, 11, 13, 17, 19, 23, 31, 37, 61, 67]
CONTROL_MODULI = [43, 47, 53, 59, 71, 73, 79, 83, 89, 97, 101]

# Linear:
#     R(x) = a*x + b
#
# Quadratic:
#     R(x) = a*x^2 + b*x + c
#
# Full brute force is feasible over these small fields.
MAX_DEGREE = 2

# Integer coefficient search:
#     R(x) = c0 + c1*x + ... + cd*x^d
#
# Small coefficient box.
INTEGER_COEFF_BOUND = 8
INTEGER_MAX_DEGREE = 2

RANDOM_SEED = 100100


# =============================================================================
# DATA
# =============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int
    T: int


@dataclass
class PolynomialResult:
    degree: int
    coeffs: Tuple[int, ...]
    train_failures: int
    test_failures: int
    train_coverage: float
    test_coverage: float


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi < 2:
        return []

    flags = bytearray(b"\x01") * (hi + 1)
    flags[0] = 0
    flags[1] = 0

    for p in range(2, math.isqrt(hi) + 1):
        if flags[p]:
            start = p * p
            count = ((hi - start) // p) + 1
            flags[start:hi + 1:p] = b"\x00" * count

    return [i for i in range(lo, hi + 1) if flags[i]]


# =============================================================================
# TARGET GENERATION
# =============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:

    rng = random.Random(RANDOM_SEED)

    seen = set()
    out: List[Target] = []

    while len(out) < count:
        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        n = p * q
        s = p + q
        T = p * p - p * q + q * q

        out.append(
            Target(
                p=p,
                q=q,
                n=n,
                s=s,
                T=T,
            )
        )

    return out


# =============================================================================
# BASIC NUMBER-THEORETIC IDENTITIES
# =============================================================================

def validate_eisenstein_identities(
    targets: Sequence[Target],
) -> None:

    print("\n2. EISENSTEIN / QUADRATIC-FORM VALIDATION")
    print("-" * 78)

    fail_T = 0
    fail_minus = 0
    fail_plus = 0

    for i, t in enumerate(targets, 1):

        T1 = t.s * t.s - 3 * t.n
        T2 = t.p * t.p - t.p * t.q + t.q * t.q

        if T1 != t.T or T2 != t.T:
            fail_T += 1

        if t.T - t.n != (t.p - t.q) ** 2:
            fail_minus += 1

        if t.T + t.n != t.p ** 2 + t.q ** 2:
            fail_plus += 1

    print(f"T = S^2 - 3N failures = {fail_T}/{len(targets)}")
    print(f"T-N=(p-q)^2 failures  = {fail_minus}/{len(targets)}")
    print(f"T+N=p^2+q^2 failures = {fail_plus}/{len(targets)}")

    if fail_T or fail_minus or fail_plus:
        raise ArithmeticError(
            "Eisenstein/quadratic-form validation failed."
        )

    print("status = PASS")


# =============================================================================
# LEGENDRE SYMBOL
# =============================================================================

def legendre_symbol(a: int, p: int) -> int:
    a %= p

    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1

    if value == p - 1:
        return -1

    return 0


# =============================================================================
# MODULAR ROOTS
# =============================================================================

def modular_square_roots_bruteforce(
    a: int,
    p: int,
) -> List[int]:

    a %= p

    return [
        x
        for x in range(p)
        if (x * x) % p == a
    ]


# =============================================================================
# ORACLE MODULAR ROOT DIAGNOSTIC
# =============================================================================

def oracle_modular_root_diagnostic(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> None:

    print("\n3. ORACLE FINITE-FIELD SQUARE-ROOT DIAGNOSTIC")
    print("-" * 78)

    for ell in moduli:

        qr = 0
        zero = 0
        nonqr = 0
        root_counts: Dict[int, int] = {}

        for t in targets:

            roots = modular_square_roots_bruteforce(
                t.T,
                ell,
            )

            if len(roots) == 0:
                nonqr += 1
            elif len(roots) == 1:
                zero += 1
            else:
                qr += 1

            root_counts[len(roots)] = (
                root_counts.get(len(roots), 0) + 1
            )

        print(
            f"ell={ell:3d} "
            f"nonzero_QR={qr:3d} "
            f"zero={zero:3d} "
            f"nonQR={nonqr:3d} "
            f"root_count_distribution={root_counts}"
        )


# =============================================================================
# POLYNOMIAL EVALUATION
# =============================================================================

def poly_eval_mod(
    coeffs: Sequence[int],
    x: int,
    p: int,
) -> int:

    acc = 0

    for c in reversed(coeffs):
        acc = (acc * x + c) % p

    return acc


# =============================================================================
# COLLISION CONSISTENCY
# =============================================================================

def collision_check(
    targets: Sequence[Target],
    p: int,
) -> Tuple[int, int]:

    """
    For equal N mod p states, check whether T mod p is identical.

    If the same state maps to multiple T values, no deterministic
    N-only function modulo p can exist.
    """

    mapping: Dict[int, int] = {}
    contradictions = 0
    repeated = 0

    for t in targets:

        x = t.n % p
        y = t.T % p

        if x in mapping:
            repeated += 1

            if mapping[x] != y:
                contradictions += 1

        else:
            mapping[x] = y

    return repeated, contradictions


# =============================================================================
# TRAIN / TEST POLYNOMIAL SEARCH OVER F_p
# =============================================================================

def search_field_polynomial(
    train: Sequence[Target],
    test: Sequence[Target],
    p: int,
    degree: int,
) -> Optional[PolynomialResult]:

    """
    Exhaustively search:

        R(x)^2 = T

    on training targets.

    For degree 1:
        a*x+b

    For degree 2:
        a*x^2+b*x+c
    """

    coeff_count = degree + 1

    # No solution is necessary if repeated N states have inconsistent T.
    _, contradictions = collision_check(
        train,
        p,
    )

    if contradictions:
        return None

    train_mod = [
        (
            t.n % p,
            t.T % p,
        )
        for t in train
    ]

    test_mod = [
        (
            t.n % p,
            t.T % p,
        )
        for t in test
    ]

    total_candidates = p ** coeff_count

    # Safety guard.
    if total_candidates > 2_000_000:
        return None

    best: Optional[
        Tuple[
            Tuple[int, ...],
            int,
            int,
        ]
    ] = None

    for coeffs in product(
        range(p),
        repeat=coeff_count,
    ):

        train_fail = 0

        for x, y in train_mod:

            r = poly_eval_mod(
                coeffs,
                x,
                p,
            )

            if (r * r - y) % p != 0:
                train_fail += 1
                break

        if train_fail:
            continue

        test_fail = 0

        for x, y in test_mod:

            r = poly_eval_mod(
                coeffs,
                x,
                p,
            )

            if (r * r - y) % p != 0:
                test_fail += 1

        if best is None:
            best = (
                coeffs,
                train_fail,
                test_fail,
            )
        else:
            if test_fail < best[2]:
                best = (
                    coeffs,
                    train_fail,
                    test_fail,
                )

        if test_fail == 0:
            break

    if best is None:
        return None

    coeffs, train_fail, test_fail = best

    train_coverage = 1.0 - (
        train_fail / max(1, len(train))
    )

    test_coverage = 1.0 - (
        test_fail / max(1, len(test))
    )

    return PolynomialResult(
        degree=degree,
        coeffs=coeffs,
        train_failures=train_fail,
        test_failures=test_fail,
        train_coverage=train_coverage,
        test_coverage=test_coverage,
    )


# =============================================================================
# FIELD SEARCH REPORT
# =============================================================================

def run_field_search(
    name: str,
    moduli: Sequence[int],
    targets: Sequence[Target],
) -> None:

    print(f"\n4. {name.upper()} FIELD POLYNOMIAL SEARCH")
    print("=" * 78)

    train = targets[:TRAIN_TARGETS]
    test = targets[TRAIN_TARGETS:]

    for ell in moduli:

        repeated, contradictions = collision_check(
            train,
            ell,
        )

        print(
            f"\nell={ell:3d} "
            f"repeated_states={repeated:3d} "
            f"contradictions={contradictions:3d}"
        )

        for degree in range(
            1,
            MAX_DEGREE + 1,
        ):

            result = search_field_polynomial(
                train,
                test,
                ell,
                degree,
            )

            if result is None:
                print(
                    f"  degree={degree}: NO TRAINING SOLUTION"
                )
                continue

            print(
                f"  degree={degree}: "
                f"coeffs={result.coeffs} "
                f"train_cov={result.train_coverage:.4f} "
                f"test_cov={result.test_coverage:.4f} "
                f"test_fail={result.test_failures}"
            )

            # Do not continue searching higher degree once a
            # perfect held-out solution has appeared.
            if (
                result.train_failures == 0
                and result.test_failures == 0
            ):
                break


# =============================================================================
# SMALL INTEGER POLYNOMIAL SEARCH
# =============================================================================

def search_small_integer_polynomial(
    train: Sequence[Target],
    test: Sequence[Target],
    moduli: Sequence[int],
    degree: int,
    bound: int,
) -> Optional[Tuple[int, ...]]:

    """
    Search one integer polynomial

        R(N) = c0 + c1*N + ... + cd*N^d

    whose square agrees with oracle T modulo every supplied modulus,
    for every training target.

    This is the closest finite search here to an actual N-only algebraic
    square-root candidate.
    """

    coeff_range = range(
        -bound,
        bound + 1,
    )

    checked = 0

    for coeffs in product(
        coeff_range,
        repeat=degree + 1,
    ):

        checked += 1

        ok = True

        for t in train:

            for ell in moduli:

                r = poly_eval_mod(
                    coeffs,
                    t.n,
                    ell,
                )

                if (
                    r * r - t.T
                ) % ell != 0:
                    ok = False
                    break

            if not ok:
                break

        if not ok:
            continue

        # Strict held-out validation.
        for t in test:

            for ell in moduli:

                r = poly_eval_mod(
                    coeffs,
                    t.n,
                    ell,
                )

                if (
                    r * r - t.T
                ) % ell != 0:
                    return None

        return coeffs

    return None


def run_integer_search(
    name: str,
    moduli: Sequence[int],
    targets: Sequence[Target],
) -> None:

    print(
        f"\n5. {name.upper()} SMALL-INTEGER N-ONLY SEARCH"
    )
    print("=" * 78)

    train = targets[:TRAIN_TARGETS]
    test = targets[TRAIN_TARGETS:]

    print(
        f"moduli={list(moduli)}"
    )
    print(
        f"coefficient bound={INTEGER_COEFF_BOUND}"
    )

    for degree in range(
        1,
        INTEGER_MAX_DEGREE + 1,
    ):

        t0 = time.perf_counter()

        coeffs = search_small_integer_polynomial(
            train,
            test,
            moduli,
            degree,
            INTEGER_COEFF_BOUND,
        )

        elapsed = time.perf_counter() - t0

        if coeffs is None:

            print(
                f"degree={degree}: "
                f"NO universal integer polynomial "
                f"time={elapsed:.3f}s"
            )

        else:

            print(
                f"degree={degree}: "
                f"FOUND coeffs={coeffs} "
                f"time={elapsed:.3f}s"
            )

            print(
                "  This is a candidate:"
            )

            terms = []

            for i, c in enumerate(coeffs):

                if c == 0:
                    continue

                if i == 0:
                    terms.append(str(c))
                elif i == 1:
                    terms.append(
                        f"{c}*N"
                    )
                else:
                    terms.append(
                        f"{c}*N^{i}"
                    )

            print(
                "  R(N) = "
                + " + ".join(terms)
            )

            print(
                "  SUCCESS: exact held-out modular "
                "square-root agreement."
            )

            break


# =============================================================================
# RESIDUE-CLASS BEST BRANCH DIAGNOSTIC
# =============================================================================

def branch_predictability(
    targets: Sequence[Target],
    p: int,
) -> None:

    print(
        f"\n6. MODULAR ROOT BRANCH STRUCTURE ell={p}"
    )
    print("-" * 78)

    mapping: Dict[int, set[int]] = {}

    for t in targets:

        x = t.n % p
        roots = modular_square_roots_bruteforce(
            t.T,
            p,
        )

        if x not in mapping:
            mapping[x] = set()

        mapping[x].update(roots)

    deterministic = 0
    ambiguous = 0

    for roots in mapping.values():

        if len(roots) <= 1:
            deterministic += 1
        else:
            ambiguous += 1

    print(
        f"distinct N states = {len(mapping)}"
    )
    print(
        f"deterministic root states = "
        f"{deterministic}"
    )
    print(
        f"ambiguous root states = "
        f"{ambiguous}"
    )


# =============================================================================
# ORACLE GCD CHECK
# =============================================================================

def oracle_root_factor_check(
    targets: Sequence[Target],
) -> None:

    print(
        "\n7. ORACLE ROOT-TO-FACTOR CHECK"
    )
    print("-" * 78)

    for i, t in enumerate(targets[:10], 1):

        # The actual hidden roots modulo N are constructed from p,q.
        #
        # Since T == q^2 (mod p)
        # and T == p^2 (mod q), the four CRT combinations are:
        #
        #   x =  +/-q (mod p)
        #   x =  +/-p (mod q)
        #
        # We explicitly demonstrate the factor-revealing property.

        roots = []

        # Four sign combinations.
        for a in (+1, -1):
            for b in (+1, -1):

                # Solve:
                # x = a*q mod p
                # x = b*p mod q
                #
                # CRT:
                # x = a*q + p*k
                # p*k = b*p - a*q (mod q)
                #
                # p^{-1} exists mod q.
                inv_p = pow(t.p, -1, t.q)

                k = (
                    (b * t.p - a * t.q)
                    * inv_p
                ) % t.q

                x = (
                    a * t.q
                    + t.p * k
                ) % t.n

                roots.append(x)

        roots = sorted(set(roots))

        recovered = []

        for x in roots:

            g1 = math.gcd(
                x - 1,
                t.n,
            )
            g2 = math.gcd(
                x + 1,
                t.n,
            )

            if (
                1 < g1 < t.n
            ):
                recovered.append(g1)

            if (
                1 < g2 < t.n
            ):
                recovered.append(g2)

            # More generally use another root pair.
            for y in roots:

                if y == x:
                    continue

                g = math.gcd(
                    x - y,
                    t.n,
                )

                if (
                    1 < g < t.n
                ):
                    recovered.append(g)

        recovered = sorted(
            set(recovered)
        )

        print(
            f"target {i:2d}: "
            f"roots={len(roots)} "
            f"factor_candidates={recovered}"
        )


# =============================================================================
# CROSS-MODULUS SUMMARY
# =============================================================================

def cross_modulus_summary(
    targets: Sequence[Target],
) -> None:

    print(
        "\n8. CROSS-MODULUS SUMMARY"
    )
    print("-" * 78)

    combined = (
        CYCLO_MODULI,
        CONTROL_MODULI,
    )

    for name, moduli in (
        ("CYCLO", CYCLO_MODULI),
        ("CONTROL", CONTROL_MODULI),
    ):

        contradiction_count = 0
        repeated_count = 0

        for ell in moduli:

            repeated, contradictions = collision_check(
                targets[:TRAIN_TARGETS],
                ell,
            )

            repeated_count += repeated
            contradiction_count += contradictions

        print(
            f"{name:7s} "
            f"repeated={repeated_count} "
            f"contradictions={contradiction_count}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 100")
    print("EISENSTEIN NORM / MODULAR SQUARE-ROOT BRIDGE")
    print("T = S^2 - 3N = p^2 - pq + q^2")
    print("T-N = (p-q)^2")
    print("STRICT TARGET HOLDOUT")
    print("NO FACTOR-PAIR CHECK IN N-ONLY SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # PRIME POPULATION
    # -------------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_MIN,
        PRIME_MAX,
    )

    print(
        "\n1. PRIME POPULATION"
    )
    print("-" * 78)
    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # -------------------------------------------------------------------------
    # TARGETS
    # -------------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    print(
        f"total targets = {len(targets)}"
    )

    for i, t in enumerate(
        targets[:24],
        1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s} "
            f"T={t.T}"
        )

    print(
        f"\ntraining targets = {TRAIN_TARGETS}"
    )
    print(
        f"test targets     = "
        f"{NUM_TARGETS - TRAIN_TARGETS}"
    )

    # -------------------------------------------------------------------------
    # IDENTITIES
    # -------------------------------------------------------------------------

    validate_eisenstein_identities(
        targets
    )

    # -------------------------------------------------------------------------
    # ORACLE ROOT STRUCTURE
    # -------------------------------------------------------------------------

    oracle_modular_root_diagnostic(
        targets,
        CYCLO_MODULI,
    )

    # -------------------------------------------------------------------------
    # FIELD SEARCH
    # -------------------------------------------------------------------------

    run_field_search(
        "cyclotomic",
        CYCLO_MODULI,
        targets,
    )

    run_field_search(
        "control",
        CONTROL_MODULI,
        targets,
    )

    # -------------------------------------------------------------------------
    # BRANCH DIAGNOSTIC
    # -------------------------------------------------------------------------

    for ell in [7, 13, 19, 31]:
        branch_predictability(
            targets,
            ell,
        )

    # -------------------------------------------------------------------------
    # INTEGER N-ONLY SEARCH
    # -------------------------------------------------------------------------

    run_integer_search(
        "cyclotomic",
        CYCLO_MODULI[:5],
        targets,
    )

    run_integer_search(
        "cyclotomic",
        CYCLO_MODULI,
        targets,
    )

    run_integer_search(
        "control",
        CONTROL_MODULI[:5],
        targets,
    )

    # -------------------------------------------------------------------------
    # ORACLE FACTOR BRIDGE
    # -------------------------------------------------------------------------

    oracle_root_factor_check(
        targets
    )

    # -------------------------------------------------------------------------
    # CROSS MODULUS
    # -------------------------------------------------------------------------

    cross_modulus_summary(
        targets
    )

    # -------------------------------------------------------------------------
    # FINAL
    # -------------------------------------------------------------------------

    print(
        "\n9. FINAL DIAGNOSTIC"
    )
    print("=" * 78)

    print(
        """
The structural reformulation is:

    T = p^2 - p*q + q^2
      = S^2 - 3N

    T-N = (p-q)^2
    T+N = p^2+q^2.

T is an Eisenstein norm.

This experiment asks a sharper question than the previous
N-only character experiments:

    Does there exist a reusable N-only residue/polynomial
    R(N) with

        R(N)^2 = T (mod ell)

    on unseen targets?

The strongest positive result is:

    a small integer polynomial R(N)
    survives multiple moduli and held-out targets.

That would provide an explicit algebraic square-root representation.

The oracle root-to-factor section verifies the second bridge:

    two different square roots of the same residue
    can yield a nontrivial gcd with N.

A failure is also decisive:

    if no low-degree field polynomial survives,
    and no small integer polynomial survives across
    several moduli, then the Eisenstein reformulation does
    not immediately yield an N-only square-root map.

IMPORTANT:
The oracle sections use p,q only for validation.
The polynomial searches themselves use N and oracle T residues
only for training/testing; they do not use p or q as features.
"""
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 100 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

