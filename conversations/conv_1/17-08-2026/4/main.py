#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 74
E1-FROM-N MODULAR STRUCTURE SEARCH
UNIVERSAL N-ONLY POLYNOMIALS VS CONTROLS
STRICT TARGET HOLDOUT
OPTIONAL CUBIC/CRT RECONSTRUCTION FROM PREDICTED E1
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Research question
-----------------
Experiment 73R showed that knowing

    E1 = (p^2 - 1)(q^2 - 1)(p + q)
       = s * ((n + 1)^2 - s^2)

is enormously informative.

The unresolved problem is therefore:

    Can E1 itself be predicted from N alone?

This experiment searches for SMALL UNIVERSAL integer polynomials

    F(n) = a0 + a1*n + ... + ad*n^d

with the SAME coefficients across all moduli.

Important:
    - coefficients are selected on training targets only;
    - test targets are never used for fitting;
    - cyclotomic and control families are compared;
    - no true s is used in prediction;
    - E1 is used only to define the oracle label during evaluation.

A second diagnostic takes the best N-only E1 formula and asks:

    F(n) mod ell
        ->
    roots of x^3 - (n+1)^2 x + F(n) = 0 mod ell
        ->
    CRT classes
        ->
    even-s-domain candidates.

This is the first experiment where the actual missing quantity E1,
rather than a sign proxy, is the object being attacked.

==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from itertools import product
from typing import Iterable, Sequence


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SEED = 74001
NUM_TARGETS = 120

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998

TRAIN_FRACTION = 0.67

# Cyclotomic family used throughout the Kappa experiments.
CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]

# Only controls with ell == 1 mod 3 have two nontrivial cube roots.
CONTROL = [
    673, 4561, 4759, 6211, 7879, 7951, 8689, 9781
]

# Search only small universal integer polynomials.
MAX_DEGREE = 3
COEF_RANGE = range(-4, 5)

# Reconstruction is deliberately conservative.
MAX_COMBINED_CLASSES = 100_000
MAX_DOMAIN_HITS = 50_000

# Number of best formulas whose CRT behavior gets inspected.
TOP_RECONSTRUCTION_FORMULAS = 5


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int
    e1: int


@dataclass(frozen=True)
class Polynomial:
    coeffs: tuple[int, ...]

    @property
    def degree(self) -> int:
        return len(self.coeffs) - 1

    def __str__(self) -> str:
        terms: list[str] = []

        for i, c in enumerate(self.coeffs):
            if c == 0:
                continue

            if i == 0:
                terms.append(str(c))
                continue

            if i == 1:
                base = "n"
            else:
                base = f"n^{i}"

            if c == 1:
                terms.append(base)
            elif c == -1:
                terms.append(f"-{base}")
            else:
                terms.append(f"{c}*{base}")

        if not terms:
            return "0"

        out = terms[0]
        for term in terms[1:]:
            if term.startswith("-"):
                out += " - " + term[1:]
            else:
                out += " + " + term

        return out


@dataclass
class Score:
    poly: Polynomial
    train_acc_cyclo: float
    test_acc_cyclo: float
    train_acc_control: float
    test_acc_control: float
    train_mae_cyclo: float
    test_mae_cyclo: float
    train_mae_control: float
    test_mae_control: float


# ---------------------------------------------------------------------------
# Prime generation
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> list[int]:
    if hi < 2:
        return []

    size = hi + 1
    is_prime = bytearray(b"\x01") * size
    is_prime[0:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(hi)) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [p for p in range(max(2, lo), hi + 1) if is_prime[p]]


# ---------------------------------------------------------------------------
# Target generation
# ---------------------------------------------------------------------------

def generate_targets(
    primes: Sequence[int],
    count: int,
    rng: random.Random,
) -> list[Target]:
    out: list[Target] = []
    seen: set[int] = set()

    while len(out) < count:
        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q
        if n in seen:
            continue

        s = p + q
        e1 = (p * p - 1) * (q * q - 1) * s

        # Independent identity checks.
        e1_alt = s * ((n + 1) * (n + 1) - s * s)

        if e1 != e1_alt:
            raise RuntimeError(
                f"E1 identity failure for target {len(out) + 1}"
            )

        cubic = s * s * s - (n + 1) * (n + 1) * s + e1
        if cubic != 0:
            raise RuntimeError(
                f"Cubic identity failure for target {len(out) + 1}"
            )

        out.append(
            Target(
                idx=len(out) + 1,
                p=p,
                q=q,
                n=n,
                s=s,
                e1=e1,
            )
        )
        seen.add(n)

    return out


# ---------------------------------------------------------------------------
# Basic modular helpers
# ---------------------------------------------------------------------------

def legendre_symbol(a: int, p: int) -> int:
    """
    Returns -1, 0, +1.
    p must be prime.
    """
    a %= p
    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1
    if value == p - 1:
        return -1

    raise RuntimeError(
        f"Unexpected Legendre value for a={a}, p={p}: {value}"
    )


def roots_of_cubic_mod_prime(
    a3: int,
    a1: int,
    a0: int,
    p: int,
) -> list[int]:
    """
    Brute-force cubic root enumeration.

    This is intentionally used only for the relatively small experimental
    primes above. It is simple and much less error-prone than introducing
    another root algorithm into this experiment.
    """
    roots = []

    for x in range(p):
        v = (
            x * x * x
            + a3 * x * x
            + a1 * x
            + a0
        ) % p

        if v == 0:
            roots.append(x)

    return roots


def crt_pair(a: int, m: int, b: int, n: int) -> tuple[int, int]:
    """
    Combine:
        x == a mod m
        x == b mod n

    Assumes gcd(m,n) = 1.
    """
    g = math.gcd(m, n)
    if g != 1:
        raise ValueError("CRT moduli are not coprime")

    t = ((b - a) * pow(m, -1, n)) % n
    x = a + m * t
    mod = m * n
    return x % mod, mod


def merge_residue_classes(
    classes: Sequence[int],
    old_mod: int,
    roots: Sequence[int],
    new_mod: int,
    limit: int,
) -> tuple[list[int], int]:
    """
    Incrementally combine CRT classes.

    Stops growing if the class count exceeds 'limit'.
    """
    out: list[int] = []

    for a in classes:
        for b in roots:
            x, _ = crt_pair(a, old_mod, b, new_mod)
            out.append(x)

            if len(out) > limit:
                return out, old_mod * new_mod

    return out, old_mod * new_mod


def domain_hits(
    residue_classes: Sequence[int],
    modulus: int,
    lo: int,
    hi: int,
) -> list[int]:
    """
    Return all x in [lo, hi] satisfying x == r mod modulus
    for at least one r.
    """
    hits: set[int] = set()

    for r in residue_classes:
        r %= modulus

        if r > hi:
            continue

        # Smallest x >= lo congruent to r mod modulus.
        if r >= lo:
            x = r
        else:
            k = (lo - r + modulus - 1) // modulus
            x = r + k * modulus

        while x <= hi:
            hits.add(x)
            if len(hits) > MAX_DOMAIN_HITS:
                return sorted(hits)
            x += modulus

    return sorted(hits)


# ---------------------------------------------------------------------------
# Polynomial evaluation
# ---------------------------------------------------------------------------

def eval_poly(poly: Polynomial, n: int, mod: int | None = None) -> int:
    """
    Horner evaluation.
    """
    value = 0

    if mod is None:
        for c in reversed(poly.coeffs):
            value = value * n + c
        return value

    n %= mod

    for c in reversed(poly.coeffs):
        value = (value * n + c) % mod

    return value


def polynomial_library() -> list[Polynomial]:
    """
    Small universal integer polynomial library.

    Important:
        The same coefficients are used for every modulus.
        This deliberately avoids per-modulus interpolation.

    We include:
      - constants
      - degree 1
      - degree 2
      - degree 3

    with coefficients in [-4,4].

    Duplicate trailing-zero representations are removed.
    """
    polys: set[Polynomial] = set()

    for degree in range(MAX_DEGREE + 1):
        for coeffs in product(COEF_RANGE, repeat=degree + 1):
            if all(c == 0 for c in coeffs):
                polys.add(Polynomial(coeffs))
                continue

            if coeffs[-1] == 0:
                continue

            # Avoid needless expressions with a huge common gcd.
            g = 0
            for c in coeffs:
                g = math.gcd(g, abs(c))

            if g > 1:
                continue

            polys.add(Polynomial(coeffs))

    return sorted(
        polys,
        key=lambda p: (p.degree, p.coeffs),
    )


# ---------------------------------------------------------------------------
# Evaluation dataset
# ---------------------------------------------------------------------------

def exact_e1_match(
    poly: Polynomial,
    target: Target,
    ell: int,
) -> bool:
    predicted = eval_poly(poly, target.n, ell)
    actual = target.e1 % ell
    return predicted == actual


def residue_error(
    poly: Polynomial,
    target: Target,
    ell: int,
) -> int:
    """
    Circular error in Z/ellZ.

    Exact prediction has error 0.
    """
    predicted = eval_poly(poly, target.n, ell)
    actual = target.e1 % ell

    d = abs(predicted - actual)
    return min(d, ell - d)


def score_poly(
    poly: Polynomial,
    train_targets: Sequence[Target],
    test_targets: Sequence[Target],
    moduli: Sequence[int],
) -> Score:
    def score_group(
        targets: Sequence[Target],
    ) -> tuple[float, float]:
        if not targets:
            return 0.0, 0.0

        total = len(targets) * len(moduli)

        correct = 0
        error_sum = 0

        for t in targets:
            for ell in moduli:
                correct += exact_e1_match(poly, t, ell)
                error_sum += residue_error(poly, t, ell)

        return correct / total, error_sum / total

    tc, te = score_group(train_targets)
    vc, ve = score_group(test_targets)

    return Score(
        poly=poly,
        train_acc_cyclo=tc,
        test_acc_cyclo=vc,
        train_acc_control=0.0,
        test_acc_control=0.0,
        train_mae_cyclo=te,
        test_mae_cyclo=ve,
        train_mae_control=0.0,
        test_mae_control=0.0,
    )


# ---------------------------------------------------------------------------
# Universal candidate search
# ---------------------------------------------------------------------------

def rank_polynomials(
    polys: Sequence[Polynomial],
    train_targets: Sequence[Target],
    test_targets: Sequence[Target],
    moduli: Sequence[int],
    top_k: int,
) -> list[Score]:
    scored: list[Score] = []

    for poly in polys:
        s = score_poly(
            poly,
            train_targets,
            test_targets,
            moduli,
        )
        scored.append(s)

    scored.sort(
        key=lambda s: (
            -s.train_acc_cyclo,
            s.train_mae_cyclo,
            s.poly.degree,
            s.poly.coeffs,
        )
    )

    return scored[:top_k]


def enrich_control_scores(
    scores: Sequence[Score],
    train_targets: Sequence[Target],
    test_targets: Sequence[Target],
    moduli: Sequence[int],
) -> list[Score]:
    out: list[Score] = []

    for s in scores:
        def group(
            targets: Sequence[Target],
        ) -> tuple[float, float]:
            if not targets:
                return 0.0, 0.0

            total = len(targets) * len(moduli)
            correct = 0
            error_sum = 0

            for t in targets:
                for ell in moduli:
                    predicted = eval_poly(s.poly, t.n, ell)
                    actual = t.e1 % ell

                    correct += predicted == actual

                    d = abs(predicted - actual)
                    error_sum += min(d, ell - d)

            return correct / total, error_sum / total

        tc, te = group(train_targets)
        vc, ve = group(test_targets)

        out.append(
            Score(
                poly=s.poly,
                train_acc_cyclo=s.train_acc_cyclo,
                test_acc_cyclo=s.test_acc_cyclo,
                train_acc_control=tc,
                test_acc_control=vc,
                train_mae_cyclo=s.train_mae_cyclo,
                test_mae_cyclo=s.test_mae_cyclo,
                train_mae_control=te,
                test_mae_control=ve,
            )
        )

    return out


# ---------------------------------------------------------------------------
# Fixed N-only baseline library
# ---------------------------------------------------------------------------

def baseline_polynomials() -> list[tuple[str, Polynomial]]:
    return [
        ("0", Polynomial((0,))),
        ("1", Polynomial((1,))),
        ("n", Polynomial((0, 1))),
        ("n+1", Polynomial((1, 1))),
        ("n-1", Polynomial((-1, 1))),
        ("n^2", Polynomial((0, 0, 1))),
        ("n^2+1", Polynomial((1, 0, 1))),
        ("n^2+n+1", Polynomial((1, 1, 1))),
        ("n^2+4n+1", Polynomial((1, 4, 1))),
        ("n^3", Polynomial((0, 0, 0, 1))),
        ("(n+1)^2", Polynomial((1, 2, 1))),
        ("n(n+1)", Polynomial((0, 1, 1))),
        ("n(n+1)^2", Polynomial((0, 1, 2, 1))),
    ]


def baseline_accuracy(
    poly: Polynomial,
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> tuple[float, float]:
    total = len(targets) * len(moduli)

    if total == 0:
        return 0.0, 0.0

    correct = 0
    error_sum = 0

    for t in targets:
        for ell in moduli:
            predicted = eval_poly(poly, t.n, ell)
            actual = t.e1 % ell

            correct += predicted == actual

            d = abs(predicted - actual)
            error_sum += min(d, ell - d)

    return correct / total, error_sum / total


# ---------------------------------------------------------------------------
# CRT reconstruction using a predicted E1
# ---------------------------------------------------------------------------

def reconstruct_from_predicted_e1(
    target: Target,
    poly: Polynomial,
    moduli: Sequence[int],
    max_classes: int = MAX_COMBINED_CLASSES,
) -> dict[str, int | bool | None]:
    """
    Use only:
        n
        predicted E1 = F(n)
        cubic modulo each ell

    No true s is used inside this function.
    """
    residue_classes = [0]
    modulus = 1

    predicted_root_counts: list[int] = []

    for ell in moduli:
        e1_pred = eval_poly(poly, target.n, ell)

        # cubic:
        # x^3 - (n+1)^2*x + E1 = 0
        a1 = -((target.n + 1) * (target.n + 1)) % ell
        a0 = e1_pred % ell

        roots = roots_of_cubic_mod_prime(
            0,
            a1,
            a0,
            ell,
        )

        predicted_root_counts.append(len(roots))

        if not roots:
            return {
                "predicted_classes": 0,
                "domain_hits": 0,
                "unique": False,
                "true_in": False,
            }

        residue_classes, modulus = merge_residue_classes(
            residue_classes,
            modulus,
            roots,
            ell,
            max_classes,
        )

        if len(residue_classes) > max_classes:
            return {
                "predicted_classes": len(residue_classes),
                "domain_hits": -1,
                "unique": False,
                "true_in": False,
            }

    hits = domain_hits(
        residue_classes,
        modulus,
        S_MIN,
        S_MAX,
    )

    return {
        "predicted_classes": len(residue_classes),
        "domain_hits": len(hits),
        "unique": len(hits) == 1,
        "true_in": target.s in hits,
    }


# ---------------------------------------------------------------------------
# Train/test reconstruction
# ---------------------------------------------------------------------------

def reconstruction_summary(
    targets: Sequence[Target],
    poly: Polynomial,
    moduli: Sequence[int],
) -> dict[str, float]:
    recovered = 0
    unique = 0
    hit_counts: list[int] = []

    for t in targets:
        result = reconstruct_from_predicted_e1(
            t,
            poly,
            moduli,
        )

        hits = int(result["domain_hits"] or 0)

        if hits >= 0:
            hit_counts.append(hits)

        recovered += int(bool(result["true_in"]))
        unique += int(bool(result["unique"]))

    mean_hits = (
        statistics.fmean(hit_counts)
        if hit_counts
        else float("nan")
    )

    median_hits = (
        statistics.median(hit_counts)
        if hit_counts
        else float("nan")
    )

    return {
        "recovered": recovered,
        "unique": unique,
        "recovery_rate": recovered / len(targets),
        "unique_rate": unique / len(targets),
        "mean_hits": mean_hits,
        "median_hits": median_hits,
    }


# ---------------------------------------------------------------------------
# Utility diagnostics
# ---------------------------------------------------------------------------

def print_header(text: str) -> None:
    print()
    print(text)
    print("-" * 78)


def mean(values: Sequence[float]) -> float:
    return statistics.fmean(values) if values else float("nan")


def median(values: Sequence[float]) -> float:
    return statistics.median(values) if values else float("nan")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    total_start = time.perf_counter()
    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 74")
    print("E1-FROM-N MODULAR STRUCTURE SEARCH")
    print("UNIVERSAL N-ONLY POLYNOMIALS VS CONTROLS")
    print("STRICT TARGET HOLDOUT")
    print("OPTIONAL CUBIC / CRT RECONSTRUCTION")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. Prime population
    # ----------------------------------------------------------------------
    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LO, PRIME_HI)
    generation_time = time.perf_counter() - t0

    print_header("1. PRIME POPULATION")
    print(f"prime population = {len(primes)}")
    print(f"generation time  = {generation_time:.6f}s")

    if len(primes) < 1000:
        raise RuntimeError("Prime population unexpectedly small")

    # ----------------------------------------------------------------------
    # 2. Targets
    # ----------------------------------------------------------------------
    targets = generate_targets(
        primes,
        NUM_TARGETS,
        rng,
    )

    print_header("2. TARGET SUMMARY")
    print(f"total targets = {len(targets)}")

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining generated targets omitted")

    # ----------------------------------------------------------------------
    # 3. Identity validation
    # ----------------------------------------------------------------------
    print_header("3. E1 IDENTITY VALIDATION")

    failures = 0

    for t in targets:
        e1_factor = (t.p * t.p - 1) * (t.q * t.q - 1) * t.s
        e1_sform = t.s * ((t.n + 1) ** 2 - t.s ** 2)
        cubic = t.s ** 3 - (t.n + 1) ** 2 * t.s + e1_factor

        ok_factor = e1_factor == t.e1
        ok_sform = e1_factor == e1_sform
        ok_cubic = cubic == 0

        if not (ok_factor and ok_sform and ok_cubic):
            failures += 1

        print(
            f"target {t.idx:3d}: "
            f"factor={ok_factor} "
            f"s-form={ok_sform} "
            f"cubic={ok_cubic}"
        )

    print(f"identity failures = {failures}")

    if failures:
        raise RuntimeError("E1 identity validation failed")

    print("status = PASS")

    # ----------------------------------------------------------------------
    # 4. Target split
    # ----------------------------------------------------------------------
    split = int(len(targets) * TRAIN_FRACTION)
    train_targets = targets[:split]
    test_targets = targets[split:]

    # ----------------------------------------------------------------------
    # 5. Baselines
    # ----------------------------------------------------------------------
    print_header("5. FIXED N-ONLY BASELINES")

    print(
        "Each formula is evaluated against the TRUE E1 residue, "
        "but the formula itself uses only n."
    )

    for name, poly in baseline_polynomials():
        c_train, e_train = baseline_accuracy(
            poly,
            train_targets,
            CYCLOTOMIC,
        )
        c_test, e_test = baseline_accuracy(
            poly,
            test_targets,
            CYCLOTOMIC,
        )

        r_train, _ = baseline_accuracy(
            poly,
            train_targets,
            CONTROL,
        )
        r_test, _ = baseline_accuracy(
            poly,
            test_targets,
            CONTROL,
        )

        print(
            f"{name:16s} "
            f"C_train={c_train:.6f} "
            f"C_test={c_test:.6f} "
            f"R_train={r_train:.6f} "
            f"R_test={r_test:.6f} "
            f"C_delta={c_test-r_test:+.6f}"
        )

    # ----------------------------------------------------------------------
    # 6. Universal polynomial search
    # ----------------------------------------------------------------------
    print_header("6. UNIVERSAL POLYNOMIAL SEARCH")

    polys = polynomial_library()
    print(f"candidate universal polynomials = {len(polys)}")
    print(
        "coefficient range = "
        f"[{min(COEF_RANGE)}, {max(COEF_RANGE)}]"
    )
    print(f"maximum degree = {MAX_DEGREE}")

    search_start = time.perf_counter()

    top = rank_polynomials(
        polys,
        train_targets,
        test_targets,
        CYCLOTOMIC,
        top_k=30,
    )

    top = enrich_control_scores(
        top,
        train_targets,
        test_targets,
        CONTROL,
    )

    search_time = time.perf_counter() - search_start
    print(f"search time = {search_time:.6f}s")

    print()
    print("TOP 20 BY CYCLOTOMIC TRAINING EXACT-MATCH RATE")
    print()

    for rank, score in enumerate(top[:20], 1):
        print(
            f"{rank:2d}. {str(score.poly):28s} "
            f"C_train={score.train_acc_cyclo:.6f} "
            f"C_test={score.test_acc_cyclo:.6f} "
            f"R_train={score.train_acc_control:.6f} "
            f"R_test={score.test_acc_control:.6f} "
            f"delta_test="
            f"{score.test_acc_cyclo-score.test_acc_control:+.6f}"
        )

    # ----------------------------------------------------------------------
    # 7. Best formula
    # ----------------------------------------------------------------------
    best = top[0]

    print_header("7. BEST UNIVERSAL N-ONLY E1 FORMULA")

    print(f"formula = {best.poly}")
    print(
        f"cyclotomic train exact = "
        f"{best.train_acc_cyclo:.6f}"
    )
    print(
        f"cyclotomic test exact  = "
        f"{best.test_acc_cyclo:.6f}"
    )
    print(
        f"control train exact    = "
        f"{best.train_acc_control:.6f}"
    )
    print(
        f"control test exact     = "
        f"{best.test_acc_control:.6f}"
    )
    print(
        f"test family delta       = "
        f"{best.test_acc_cyclo-best.test_acc_control:+.6f}"
    )
    print(
        f"cyclotomic train MAE    = "
        f"{best.train_mae_cyclo:.6f}"
    )
    print(
        f"cyclotomic test MAE     = "
        f"{best.test_mae_cyclo:.6f}"
    )

    # ----------------------------------------------------------------------
    # 8. Per-modulus performance of best formula
    # ----------------------------------------------------------------------
    print_header("8. BEST FORMULA BY MODULUS")

    for family_name, moduli in (
        ("CYCLOTOMIC", CYCLOTOMIC),
        ("CONTROL", CONTROL),
    ):
        print()
        print(family_name)

        for ell in moduli:
            train_acc, train_mae = baseline_accuracy(
                best.poly,
                train_targets,
                [ell],
            )
            test_acc, test_mae = baseline_accuracy(
                best.poly,
                test_targets,
                [ell],
            )

            print(
                f"  ell={ell:5d} "
                f"train={train_acc:.6f} "
                f"test={test_acc:.6f} "
                f"train_err={train_mae:.4f} "
                f"test_err={test_mae:.4f}"
            )

    # ----------------------------------------------------------------------
    # 9. Leave-one-modulus-out diagnostic
    # ----------------------------------------------------------------------
    print_header("9. LEAVE-ONE-MODULUS-OUT DIAGNOSTIC")

    # This is NOT a fit. The universal polynomial is already fixed from
    # target-level training. We simply report whether it works on each
    # modulus independently.
    c_loo: list[float] = []
    r_loo: list[float] = []

    for ell in CYCLOTOMIC:
        test_acc, _ = baseline_accuracy(
            best.poly,
            test_targets,
            [ell],
        )
        c_loo.append(test_acc)
        print(f"C ell={ell:5d} test={test_acc:.6f}")

    print(
        f"cyclotomic mean = {mean(c_loo):.6f}"
    )

    for ell in CONTROL:
        test_acc, _ = baseline_accuracy(
            best.poly,
            test_targets,
            [ell],
        )
        r_loo.append(test_acc)
        print(f"R ell={ell:5d} test={test_acc:.6f}")

    print(
        f"control mean    = {mean(r_loo):.6f}"
    )

    # ----------------------------------------------------------------------
    # 10. Information-theoretic residue compression
    # ----------------------------------------------------------------------
    print_header("10. MODULAR E1 COMPRESSION")

    print(
        "For each target and modulus, the oracle gives one E1 residue."
    )
    print(
        "The N-only formula predicts one residue."
    )
    print(
        "We measure exact prediction plus circular residue error."
    )

    for family_name, moduli in (
        ("CYCLOTOMIC", CYCLOTOMIC),
        ("CONTROL", CONTROL),
    ):
        exact_rates = []
        errors = []

        for t in test_targets:
            for ell in moduli:
                exact_rates.append(
                    int(
                        eval_poly(best.poly, t.n, ell)
                        == t.e1 % ell
                    )
                )
                errors.append(
                    residue_error(best.poly, t, ell)
                )

        print(
            f"{family_name:11s} "
            f"exact={mean(exact_rates):.6f} "
            f"mean_circular_error={mean(errors):.6f}"
        )

    # ----------------------------------------------------------------------
    # 11. Oracle E1 CRT upper-bound diagnostic
    # ----------------------------------------------------------------------
    print_header("11. ORACLE E1 CRT UPPER-BOUND")

    print(
        "This section deliberately uses TRUE E1."
    )
    print(
        "It quantifies the maximum information obtainable once E1 is known."
    )

    for family_name, moduli in (
        ("C3", CYCLOTOMIC[:3]),
        ("C5", CYCLOTOMIC[:5]),
        ("C7", CYCLOTOMIC[:7]),
    ):
        summary = reconstruction_summary(
            test_targets,
            Polynomial((0,)),  # placeholder, not used below
            moduli,
        )

        # Re-run with TRUE-E1 behavior implemented explicitly.
        recovered = 0
        unique = 0
        hit_counts: list[int] = []

        for t in test_targets:
            classes = [0]
            modulus = 1

            for ell in moduli:
                true_e1 = t.e1 % ell
                a1 = -((t.n + 1) ** 2) % ell
                a0 = true_e1

                roots = roots_of_cubic_mod_prime(
                    0,
                    a1,
                    a0,
                    ell,
                )

                if not roots:
                    classes = []
                    break

                classes, modulus = merge_residue_classes(
                    classes,
                    modulus,
                    roots,
                    ell,
                    MAX_COMBINED_CLASSES,
                )

                if len(classes) > MAX_COMBINED_CLASSES:
                    break

            hits = domain_hits(
                classes,
                modulus,
                S_MIN,
                S_MAX,
            )

            recovered += int(t.s in hits)
            unique += int(len(hits) == 1)
            hit_counts.append(len(hits))

        print(
            f"{family_name:4s} "
            f"recovered={recovered}/{len(test_targets)} "
            f"unique={unique}/{len(test_targets)} "
            f"mean_hits={mean(hit_counts):.3f} "
            f"median_hits={median(hit_counts):.3f}"
        )

    # ----------------------------------------------------------------------
    # 12. Predicted-E1 CRT reconstruction
    # ----------------------------------------------------------------------
    print_header("12. PREDICTED-E1 CRT RECONSTRUCTION")

    print(
        "Only the best N-only polynomial is used."
    )
    print(
        "TRUE E1 is NOT supplied to the reconstruction."
    )

    for name, moduli in (
        ("C3", CYCLOTOMIC[:3]),
        ("C5", CYCLOTOMIC[:5]),
        ("C7", CYCLOTOMIC[:7]),
    ):
        summary = reconstruction_summary(
            test_targets,
            best.poly,
            moduli,
        )

        print(
            f"{name:4s} "
            f"recovered={int(summary['recovered'])}/"
            f"{len(test_targets)} "
            f"recovery_rate={summary['recovery_rate']:.6f} "
            f"unique={int(summary['unique'])}/"
            f"{len(test_targets)} "
            f"unique_rate={summary['unique_rate']:.6f} "
            f"mean_hits={summary['mean_hits']:.3f} "
            f"median_hits={summary['median_hits']:.3f}"
        )

    # ----------------------------------------------------------------------
    # 13. Per-target examples
    # ----------------------------------------------------------------------
    print_header("13. SAMPLE TEST-TARGET DIAGNOSTICS")

    for t in test_targets[:12]:
        print(
            f"target {t.idx:3d}: n={t.n} s={t.s}"
        )

        for ell in CYCLOTOMIC[:7]:
            actual = t.e1 % ell
            predicted = eval_poly(best.poly, t.n, ell)

            print(
                f"  ell={ell:4d} "
                f"E1={actual:5d} "
                f"pred={predicted:5d} "
                f"ok={actual == predicted}"
            )

    # ----------------------------------------------------------------------
    # 14. Final diagnostic
    # ----------------------------------------------------------------------
    print_header("14. FINAL DIAGNOSTIC")

    c_test = best.test_acc_cyclo
    r_test = best.test_acc_control
    delta = c_test - r_test

    print(
        f"best universal formula = {best.poly}"
    )
    print(
        f"cyclotomic OOS exact E1 rate = {c_test:.6f}"
    )
    print(
        f"control OOS exact E1 rate    = {r_test:.6f}"
    )
    print(
        f"OOS family delta             = {delta:+.6f}"
    )

    print()
    print("Interpretation")
    print("--------------")

    if c_test < 0.10:
        print(
            "NO E1 PREDICTION: the searched universal polynomial family "
            "does not predict E1 residues."
        )
    elif c_test < 0.50:
        print(
            "WEAK E1 STRUCTURE: some N-only modular correlation exists, "
            "but exact reconstruction is poor."
        )
    elif c_test < 0.80:
        print(
            "INTERESTING E1 STRUCTURE: the universal polynomial predicts "
            "a substantial fraction of E1 residues out of sample."
        )
    else:
        print(
            "STRONG E1 STRUCTURE: a universal N-only polynomial predicts "
            "most E1 residues out of sample."
        )

    if delta > 0.05:
        print(
            "Cyclotomic advantage exceeds the control family by more than "
            "five percentage points."
        )
    else:
        print(
            "No strong cyclotomic advantage over the matched control family."
        )

    print()
    print(
        "CRITICAL TEST:"
    )
    print(
        "A genuinely important result requires BOTH:"
    )
    print(
        "  1. high out-of-sample E1 residue prediction, and"
    )
    print(
        "  2. actual s reconstruction from predicted E1."
    )
    print()
    print(
        "If E1 prediction remains near chance, the next useful direction "
        "is not more CRT work; it is returning to relation generation."
    )

    elapsed = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("EXPERIMENT 74 COMPLETE")
    print(f"total runtime = {elapsed:.6f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()

