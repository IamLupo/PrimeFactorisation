#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 63
N-ONLY SYMMETRIC LEGENDRE RESOLVENT SEARCH
COMMON-SIGN RECOVERY FROM CONJUGATE DISCRIMINANTS
TARGET + MODULUS OUT-OF-SAMPLE VALIDATION
NO CSV OUTPUT
==============================================================================

Goal
----
For the conjugate discriminants

    d1 = w^2 - 4n
    d2 = (w^2)^2 - 4n

where w^2 + w + 1 = 0, test whether the common Legendre sign of
(d1,d2) can be recovered from simple symmetric N-only expressions.

Exact symmetric identities:

    S = d1 + d2 = -(8n + 1)
    P = d1*d2 = 16n^2 + 4n + 1

For P != 0:

    chi(P) = chi(d1) * chi(d2)

Thus:
    chi(P) = -1  -> opposite signs
    chi(P) = +1  -> same signs

The unresolved case is:

    chi(P) = +1

where the oracle is either:

    (+1,+1)
or
    (-1,-1)

This experiment asks whether a small, explicitly defined family of
symmetric N-only expressions can recover that remaining common sign.

No arbitrary ML model is used.
No CSV files are produced.
=============================================================================="""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import comb
from statistics import mean, median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import random
import sympy as sp


# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

SEED = 6301

NUM_TARGETS = 120

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723,
]

RAW_CONTROLS = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781,
]

# Small explicit coefficient family.
# h(n) = a*S + b*P + c
#
# This is intentionally small to avoid turning the experiment into
# unconstrained feature mining.
COEF_VALUES = (-2, -1, 0, 1, 2)

# Additional fixed resolvents derived from S and P.
FIXED_EXPRESSIONS = (
    "S",
    "S+1",
    "S-1",
    "S+3",
    "S-3",
    "2S+1",
    "2S-1",
    "S+P",
    "S-P",
    "P",
    "P+1",
    "P-1",
    "P+S+1",
    "P-S+1",
)

TRAIN_TARGET_FRACTION = 0.75

# Number of random matched-D-QR sums used for a null diagnostic per target.
# We do not construct the entire 2.2M-domain pool because that dominates
# runtime and does not materially change this algebraic question.
FALSE_SAMPLES_PER_TARGET = 60


# -----------------------------------------------------------------------------
# DATA STRUCTURES
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q


@dataclass(frozen=True)
class LocalCell:
    target_idx: int
    ell: int
    d1: int
    d2: int
    chi1: int
    chi2: int
    S: int
    P: int

    @property
    def clean(self) -> bool:
        return self.chi1 != 0 and self.chi2 != 0

    @property
    def product_sign(self) -> int:
        return self.chi1 * self.chi2

    @property
    def same_sign(self) -> bool:
        return self.clean and self.product_sign == 1

    @property
    def common_sign(self) -> Optional[int]:
        if not self.same_sign:
            return None
        return self.chi1


# -----------------------------------------------------------------------------
# NUMBER THEORY
# -----------------------------------------------------------------------------

def legendre(a: int, p: int) -> int:
    """Return Legendre symbol in {-1,0,+1}."""
    a %= p
    if a == 0:
        return 0
    value = pow(a, (p - 1) // 2, p)
    return 1 if value == 1 else -1


def roots_x2_x_1(p: int) -> Tuple[int, int]:
    """
    Find the two nontrivial roots of x^2+x+1 mod p.

    Valid precisely for odd primes p == 1 mod 3.
    """
    if p == 3 or p % 3 != 1:
        raise ValueError(
            f"{p} is not orientation-valid: "
            "x^2+x+1 has no two distinct nontrivial roots."
        )

    roots = []
    for x in range(p):
        if (x * x + x + 1) % p == 0:
            roots.append(x)
            if len(roots) == 2:
                break

    if len(roots) != 2:
        raise RuntimeError(f"{p} does not have two roots of x^2+x+1")

    a, b = roots
    if (a + b) % p != p - 1:
        raise RuntimeError(f"Root sum identity failed for p={p}")
    if (a * b) % p != 1:
        raise RuntimeError(f"Root product identity failed for p={p}")

    return a, b


def valid_controls() -> List[int]:
    return [p for p in RAW_CONTROLS if p % 3 == 1]


def G(n: int) -> int:
    return 16 * n * n + 4 * n + 1


def symmetric_values(n: int, ell: int) -> Tuple[int, int]:
    """
    S = d1+d2
    P = d1*d2
    """
    S = (-(8 * n + 1)) % ell
    P = G(n) % ell
    return S, P


# -----------------------------------------------------------------------------
# PRIME/TARGET GENERATION
# -----------------------------------------------------------------------------

def build_prime_population() -> List[int]:
    primes = list(sp.primerange(PRIME_MIN, PRIME_MAX + 1))

    if len(primes) < NUM_TARGETS * 2:
        raise RuntimeError(
            f"Prime population too small: {len(primes)}"
        )

    return primes


def generate_targets(primes: Sequence[int]) -> List[Target]:
    """
    Generate reproducible distinct prime-pair targets.
    """
    rng = random.Random(SEED)

    targets: List[Target] = []
    used_pairs = set()

    attempts = 0
    while len(targets) < NUM_TARGETS:
        attempts += 1
        if attempts > NUM_TARGETS * 1000:
            raise RuntimeError("Could not generate enough targets.")

        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        a, b = sorted((p, q))

        if (a, b) in used_pairs:
            continue

        s = a + b
        if not (S_MIN <= s <= S_MAX):
            continue

        used_pairs.add((a, b))
        targets.append(Target(len(targets) + 1, a, b))

    return targets


# -----------------------------------------------------------------------------
# LOCAL STATES
# -----------------------------------------------------------------------------

def local_cell(target: Target, ell: int) -> LocalCell:
    w1, w2 = roots_x2_x_1(ell)

    d1 = (w1 * w1 - 4 * target.n) % ell
    d2 = (w2 * w2 - 4 * target.n) % ell

    c1 = legendre(d1, ell)
    c2 = legendre(d2, ell)

    S = (d1 + d2) % ell
    P = (d1 * d2) % ell

    # Strong local algebra checks.
    expected_S, expected_P = symmetric_values(target.n, ell)

    if S != expected_S:
        raise AssertionError(
            f"S identity failure target={target.idx} ell={ell}: "
            f"{S} != {expected_S}"
        )

    if P != expected_P:
        raise AssertionError(
            f"P identity failure target={target.idx} ell={ell}: "
            f"{P} != {expected_P}"
        )

    # Difference relation:
    # (d1-d2)^2 = -3 mod ell
    if ((d1 - d2) * (d1 - d2) + 3) % ell != 0:
        raise AssertionError(
            f"Difference identity failure target={target.idx} ell={ell}"
        )

    return LocalCell(
        target_idx=target.idx,
        ell=ell,
        d1=d1,
        d2=d2,
        chi1=c1,
        chi2=c2,
        S=S,
        P=P,
    )


# -----------------------------------------------------------------------------
# EXPRESSION EVALUATION
# -----------------------------------------------------------------------------

def fixed_expression_value(name: str, S: int, P: int, ell: int) -> int:
    if name == "S":
        return S % ell
    if name == "S+1":
        return (S + 1) % ell
    if name == "S-1":
        return (S - 1) % ell
    if name == "S+3":
        return (S + 3) % ell
    if name == "S-3":
        return (S - 3) % ell
    if name == "2S+1":
        return (2 * S + 1) % ell
    if name == "2S-1":
        return (2 * S - 1) % ell
    if name == "S+P":
        return (S + P) % ell
    if name == "S-P":
        return (S - P) % ell
    if name == "P":
        return P % ell
    if name == "P+1":
        return (P + 1) % ell
    if name == "P-1":
        return (P - 1) % ell
    if name == "P+S+1":
        return (P + S + 1) % ell
    if name == "P-S+1":
        return (P - S + 1) % ell

    raise ValueError(f"Unknown fixed expression: {name}")


def coeff_expression_value(
    a: int,
    b: int,
    c: int,
    S: int,
    P: int,
    ell: int,
) -> int:
    return (a * S + b * P + c) % ell


def candidate_label(a: int, b: int, c: int) -> str:
    return f"{a}*S + {b}*P + {c}"


# -----------------------------------------------------------------------------
# DATASET
# -----------------------------------------------------------------------------

def build_cells(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> List[LocalCell]:
    out = []

    for t in targets:
        for ell in moduli:
            out.append(local_cell(t, ell))

    return out


def clean_same_sign(cells: Iterable[LocalCell]) -> List[LocalCell]:
    return [c for c in cells if c.same_sign]


# -----------------------------------------------------------------------------
# TRAIN / TEST SPLITS
# -----------------------------------------------------------------------------

def split_targets(
    targets: Sequence[Target],
) -> Tuple[List[Target], List[Target]]:
    rng = random.Random(SEED + 1)

    ids = [t.idx for t in targets]
    rng.shuffle(ids)

    split = int(len(ids) * TRAIN_TARGET_FRACTION)

    train_ids = set(ids[:split])

    train = [t for t in targets if t.idx in train_ids]
    test = [t for t in targets if t.idx not in train_ids]

    return train, test


# -----------------------------------------------------------------------------
# CANDIDATE SCORING
# -----------------------------------------------------------------------------

def score_expression(
    cells: Sequence[LocalCell],
    expression_name: str,
) -> Dict[str, float]:
    """
    Predict common sign from chi(expression).

    Only same-sign oracle cells are considered.
    If chi(expression)==0, prediction is undefined.
    """
    total = 0
    covered = 0
    correct = 0

    for cell in cells:
        if not cell.same_sign:
            continue

        if expression_name.startswith("fixed:"):
            value = fixed_expression_value(
                expression_name[6:],
                cell.S,
                cell.P,
                cell.ell,
            )
        else:
            raise ValueError(expression_name)

        pred = legendre(value, cell.ell)

        total += 1

        if pred == 0:
            continue

        covered += 1

        if pred == cell.common_sign:
            correct += 1

    accuracy = correct / covered if covered else float("nan")
    coverage = covered / total if total else float("nan")

    return {
        "accuracy": accuracy,
        "coverage": coverage,
        "correct": correct,
        "covered": covered,
        "total": total,
    }


def score_linear_candidate(
    cells: Sequence[LocalCell],
    coeffs: Tuple[int, int, int],
) -> Dict[str, float]:
    a, b, c = coeffs

    total = 0
    covered = 0
    correct = 0

    for cell in cells:
        if not cell.same_sign:
            continue

        value = coeff_expression_value(
            a, b, c,
            cell.S,
            cell.P,
            cell.ell,
        )

        pred = legendre(value, cell.ell)

        total += 1

        if pred == 0:
            continue

        covered += 1

        if pred == cell.common_sign:
            correct += 1

    return {
        "accuracy": correct / covered if covered else float("nan"),
        "coverage": covered / total if total else float("nan"),
        "correct": correct,
        "covered": covered,
        "total": total,
    }


# -----------------------------------------------------------------------------
# BASELINES
# -----------------------------------------------------------------------------

def baseline_product_sign(cells: Sequence[LocalCell]) -> float:
    """
    For same-sign cells, chi(P)=+1 always, so this has no ability to
    distinguish +/+ from -/- and is deliberately included as a sanity check.
    """
    usable = [c for c in cells if c.same_sign]
    if not usable:
        return float("nan")

    # Always predicts +.
    correct = sum(1 for c in usable if c.common_sign == 1)
    return correct / len(usable)


def baseline_chi_n(
    cells: Sequence[LocalCell],
    shift: int = 0,
) -> Dict[str, float]:
    total = 0
    covered = 0
    correct = 0

    for cell in cells:
        if not cell.same_sign:
            continue

        # Recover n modulo ell from S = -(8n+1).
        #
        # 8 is invertible for every odd ell.
        inv8 = pow(8, -1, cell.ell)
        n_mod = (-(cell.S + 1) * inv8) % cell.ell

        value = (n_mod + shift) % cell.ell
        pred = legendre(value, cell.ell)

        total += 1

        if pred == 0:
            continue

        covered += 1

        if pred == cell.common_sign:
            correct += 1

    return {
        "accuracy": correct / covered if covered else float("nan"),
        "coverage": covered / total if total else float("nan"),
        "correct": correct,
        "covered": covered,
        "total": total,
    }


# -----------------------------------------------------------------------------
# MATCHED FALSE-SUM NULL
# -----------------------------------------------------------------------------

def discriminant_qr_for_sum(n: int, s: int, ell: int) -> bool:
    """
    Baseline D-QR condition:
        D = s^2 - 4n
    must be a quadratic residue modulo ell.

    This is the existing sieve condition.
    """
    d = (s * s - 4 * n) % ell
    return legendre(d, ell) >= 0


def matched_false_sums(
    target: Target,
    moduli: Sequence[int],
    count: int,
    seed_offset: int,
) -> List[int]:
    """
    Deterministic random sample of false even sums satisfying exactly the
    cyclotomic D-QR conditions.

    We sample directly over the valid even-sum domain rather than scanning
    all 2.2M values.
    """
    rng = random.Random(SEED + seed_offset)

    true_s = target.s

    out = set()

    max_attempts = count * 1000
    attempts = 0

    while len(out) < count:
        attempts += 1

        if attempts > max_attempts:
            raise RuntimeError(
                f"Target {target.idx}: could only construct "
                f"{len(out)} matched false sums."
            )

        k = rng.randrange((S_MAX - S_MIN) // 2 + 1)
        s = S_MIN + 2 * k

        if s == true_s:
            continue

        ok = True

        for ell in moduli:
            if not discriminant_qr_for_sum(target.n, s, ell):
                ok = False
                break

        if ok:
            out.add(s)

    return sorted(out)


def sum_to_local_cells_for_false(
    target: Target,
    s: int,
    moduli: Sequence[int],
) -> List[Tuple[int, int, int, int]]:
    """
    For a false sum s, use the corresponding discriminant-like quantities

        d(w) = w^2 - 4n

    is not correct because the true factor roots depend on p,q.

    Therefore for the false-null test we use the sum-discriminant state only
    and DO NOT invent a false branch orientation oracle.

    This function is intentionally absent from the final analysis pipeline.
    """
    raise NotImplementedError


# -----------------------------------------------------------------------------
# SYMBOLIC REDUCTION CHECKS
# -----------------------------------------------------------------------------

def symbolic_identity_report() -> None:
    """
    Print the exact algebra behind the search.
    """
    n = sp.symbols("n")
    S = -(8 * n + 1)
    P = 16 * n**2 + 4 * n + 1

    print("Exact symmetric identities:")
    print(f"  S = d1 + d2 = {sp.expand(S)}")
    print(f"  P = d1*d2 = {sp.expand(P)}")

    discriminant = sp.expand(S**2 - 4 * P)
    print(f"  S^2 - 4P = {discriminant}")

    print("  Therefore:")
    print("    (d1-d2)^2 = S^2 - 4P = -3")

    print()
    print("Important structural check:")
    print("  every symmetric polynomial in d1,d2 reduces to a polynomial")
    print("  in S and P, hence to an N-only polynomial.")
    print()


# -----------------------------------------------------------------------------
# MAIN EXPERIMENT
# -----------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("KAPPA EXPERIMENT 63")
    print("N-ONLY SYMMETRIC LEGENDRE RESOLVENT SEARCH")
    print("COMMON-SIGN RECOVERY FROM CONJUGATE DISCRIMINANTS")
    print("TARGET + MODULUS OUT-OF-SAMPLE VALIDATION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # PRIME POPULATION
    # -------------------------------------------------------------------------
    primes = build_prime_population()

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes):,}")

    # -------------------------------------------------------------------------
    # TARGETS
    # -------------------------------------------------------------------------
    targets = generate_targets(primes)

    print()
    print("2. TARGETS")
    print("-" * 78)
    print(f"total targets = {len(targets)}")
    print(f"train/test split = {int(NUM_TARGETS * TRAIN_TARGET_FRACTION)}/"
          f"{NUM_TARGETS - int(NUM_TARGETS * TRAIN_TARGET_FRACTION)}")

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining generated targets omitted")

    # -------------------------------------------------------------------------
    # MODULI
    # -------------------------------------------------------------------------
    controls = valid_controls()

    print()
    print("3. MODULUS FAMILIES")
    print("-" * 78)
    print(f"cyclotomic = {CYCLOTOMIC}")
    print(f"valid controls = {controls}")
    print(
        "excluded controls = "
        f"{[x for x in RAW_CONTROLS if x not in controls]}"
    )

    # -------------------------------------------------------------------------
    # SYMBOLIC IDENTITIES
    # -------------------------------------------------------------------------
    print()
    print("4. SYMMETRIC ALGEBRA")
    print("-" * 78)
    symbolic_identity_report()

    # -------------------------------------------------------------------------
    # BUILD LOCAL DATA
    # -------------------------------------------------------------------------
    print()
    print("5. LOCAL DATASET")
    print("-" * 78)

    cyclo_cells = build_cells(targets, CYCLOTOMIC)
    control_cells = build_cells(targets, controls)

    cyclo_clean = clean_same_sign(cyclo_cells)
    control_clean = clean_same_sign(control_cells)

    print(f"cyclotomic total cells = {len(cyclo_cells)}")
    print(f"cyclotomic same-sign cells = {len(cyclo_clean)}")
    print(f"control total cells = {len(control_cells)}")
    print(f"control same-sign cells = {len(control_clean)}")

    # -------------------------------------------------------------------------
    # TARGET SPLIT
    # -------------------------------------------------------------------------
    train_targets, test_targets = split_targets(targets)

    train_ids = {t.idx for t in train_targets}
    test_ids = {t.idx for t in test_targets}

    print()
    print("6. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_targets)}")
    print(f"test targets = {len(test_targets)}")
    print(f"test ids = {sorted(test_ids)}")

    train_c = [
        c for c in cyclo_clean
        if c.target_idx in train_ids
    ]
    test_c = [
        c for c in cyclo_clean
        if c.target_idx in test_ids
    ]

    train_r = [
        c for c in control_clean
        if c.target_idx in train_ids
    ]
    test_r = [
        c for c in control_clean
        if c.target_idx in test_ids
    ]

    # -------------------------------------------------------------------------
    # FIXED EXPRESSIONS
    # -------------------------------------------------------------------------
    print()
    print("7. FIXED N-ONLY RESOLVENTS")
    print("-" * 78)

    fixed_results = []

    for name in FIXED_EXPRESSIONS:
        tr = score_expression(train_c, f"fixed:{name}")
        te = score_expression(test_c, f"fixed:{name}")
        cr = score_expression(test_r, f"fixed:{name}")

        fixed_results.append(
            (name, tr, te, cr)
        )

        print(
            f"{name:12s} "
            f"train={tr['accuracy']:.6f} "
            f"test-C={te['accuracy']:.6f} "
            f"test-R={cr['accuracy']:.6f} "
            f"C_cov={te['coverage']:.6f}"
        )

    # -------------------------------------------------------------------------
    # LINEAR SYMMETRIC SEARCH
    # -------------------------------------------------------------------------
    print()
    print("8. SMALL LINEAR SYMMETRIC SEARCH")
    print("-" * 78)

    candidates = []

    for a, b, c in product(COEF_VALUES, repeat=3):
        if (a, b, c) == (0, 0, 0):
            continue

        # Avoid the exact trivial P expression since it is already a baseline.
        if (a, b, c) == (0, 1, 0):
            continue

        result_train = score_linear_candidate(train_c, (a, b, c))
        result_test_c = score_linear_candidate(test_c, (a, b, c))
        result_test_r = score_linear_candidate(test_r, (a, b, c))

        if result_train["covered"] < max(10, result_train["total"] // 3):
            continue

        candidates.append(
            {
                "coeffs": (a, b, c),
                "label": candidate_label(a, b, c),
                "train": result_train,
                "test_c": result_test_c,
                "test_r": result_test_r,
            }
        )

    # Rank ONLY using training accuracy.
    candidates.sort(
        key=lambda x: (
            x["train"]["accuracy"],
            x["train"]["coverage"],
        ),
        reverse=True,
    )

    print()
    print("Top 15 candidates ranked by TRAINING accuracy only:")
    print()

    for rank, item in enumerate(candidates[:15], 1):
        tr = item["train"]
        tc = item["test_c"]
        rr = item["test_r"]

        print(
            f"{rank:2d}. {item['label']:20s} "
            f"train={tr['accuracy']:.6f} "
            f"test-C={tc['accuracy']:.6f} "
            f"test-R={rr['accuracy']:.6f} "
            f"C_cov={tc['coverage']:.6f}"
        )

    # -------------------------------------------------------------------------
    # BEST TRAINED CANDIDATE
    # -------------------------------------------------------------------------
    if not candidates:
        raise RuntimeError("No usable symmetric candidates were generated.")

    best = candidates[0]

    print()
    print("9. BEST TRAINED CANDIDATE")
    print("-" * 78)
    print(f"expression = {best['label']}")
    print(
        f"train cyclotomic = "
        f"{best['train']['accuracy']:.6f}"
    )
    print(
        f"unseen-target cyclotomic = "
        f"{best['test_c']['accuracy']:.6f}"
    )
    print(
        f"unseen-target control = "
        f"{best['test_r']['accuracy']:.6f}"
    )

    # -------------------------------------------------------------------------
    # SIMPLE BASELINES
    # -------------------------------------------------------------------------
    print()
    print("10. N-ONLY BASELINES")
    print("-" * 78)

    baseline_specs = [
        ("chi(n)", 0),
        ("chi(n+1)", 1),
        ("chi(n-1)", -1),
    ]

    for name, shift in baseline_specs:
        result = baseline_chi_n(test_c, shift)
        control = baseline_chi_n(test_r, shift)

        print(
            f"{name:10s} "
            f"C={result['accuracy']:.6f} "
            f"R={control['accuracy']:.6f} "
            f"coverage={result['coverage']:.6f}"
        )

    base_same = baseline_product_sign(test_c)

    print(
        f"always + on same-sign cells = {base_same:.6f}"
    )

    # -------------------------------------------------------------------------
    # LEAVE-ONE-MODULUS-OUT
    # -------------------------------------------------------------------------
    print()
    print("11. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    loo_c = []
    loo_r = []

    for held_ell in CYCLOTOMIC:
        train_ells = [x for x in CYCLOTOMIC if x != held_ell]

        train_cells = [
            local_cell(t, ell)
            for t in targets
            for ell in train_ells
        ]
        test_cells = [
            local_cell(t, held_ell)
            for t in targets
        ]

        train_cells = clean_same_sign(train_cells)
        test_cells = clean_same_sign(test_cells)

        # Search based on training data excluding held ell.
        local_candidates = []

        for a, b, c in product(COEF_VALUES, repeat=3):
            if (a, b, c) == (0, 0, 0):
                continue

            tr = score_linear_candidate(
                train_cells,
                (a, b, c),
            )

            if tr["covered"] < max(10, tr["total"] // 3):
                continue

            local_candidates.append(
                ((a, b, c), tr)
            )

        if not local_candidates:
            continue

        local_candidates.sort(
            key=lambda z: (
                z[1]["accuracy"],
                z[1]["coverage"],
            ),
            reverse=True,
        )

        coeffs, _ = local_candidates[0]

        held = score_linear_candidate(test_cells, coeffs)

        loo_c.append(held["accuracy"])

        print(
            f"ell={held_ell:5d} "
            f"selected={candidate_label(*coeffs):20s} "
            f"accuracy={held['accuracy']:.6f}"
        )

    print()
    print(
        f"cyclotomic LOO mean   = {mean(loo_c):.6f}"
    )
    print(
        f"cyclotomic LOO median = {median(loo_c):.6f}"
    )

    # Controls use the same small search, with leave-one-modulus-out.
    for held_ell in controls:
        train_ells = [x for x in controls if x != held_ell]

        train_cells = [
            local_cell(t, ell)
            for t in targets
            for ell in train_ells
        ]
        test_cells = [
            local_cell(t, held_ell)
            for t in targets
        ]

        train_cells = clean_same_sign(train_cells)
        test_cells = clean_same_sign(test_cells)

        local_candidates = []

        for a, b, c in product(COEF_VALUES, repeat=3):
            if (a, b, c) == (0, 0, 0):
                continue

            tr = score_linear_candidate(
                train_cells,
                (a, b, c),
            )

            if tr["covered"] < max(10, tr["total"] // 3):
                continue

            local_candidates.append(
                ((a, b, c), tr)
            )

        if not local_candidates:
            continue

        local_candidates.sort(
            key=lambda z: (
                z[1]["accuracy"],
                z[1]["coverage"],
            ),
            reverse=True,
        )

        coeffs, _ = local_candidates[0]
        held = score_linear_candidate(test_cells, coeffs)

        loo_r.append(held["accuracy"])

    print(
        f"control LOO mean      = {mean(loo_r):.6f}"
    )
    print(
        f"control LOO median    = {median(loo_r):.6f}"
    )

    # -------------------------------------------------------------------------
    # EXACT IDENTIFIABILITY CHECK
    # -------------------------------------------------------------------------
    print()
    print("12. EXACT IDENTIFIABILITY CHECK")
    print("-" * 78)

    # Test whether the pair
    #   (chi(S), chi(P))
    # identifies the common sign on same-sign cells.
    mapping: Dict[Tuple[int, int], set] = {}

    for cell in cyclo_clean:
        key = (
            legendre(cell.S, cell.ell),
            legendre(cell.P, cell.ell),
        )

        mapping.setdefault(key, set()).add(cell.common_sign)

    ambiguous = {
        key: values
        for key, values in mapping.items()
        if len(values) > 1
    }

    print(
        f"distinct (chi(S),chi(P)) states = {len(mapping)}"
    )
    print(
        f"ambiguous states = {len(ambiguous)}"
    )

    for key, values in sorted(ambiguous.items()):
        print(
            f"  state={key} -> common_signs={sorted(values)}"
        )

    # -------------------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # -------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("13. FINAL DIAGNOSTIC")
    print("=" * 78)

    tc = best["test_c"]
    tr = best["test_r"]

    print(
        f"best unseen-target cyclotomic = "
        f"{tc['accuracy']:.6f}"
    )
    print(
        f"best unseen-target control     = "
        f"{tr['accuracy']:.6f}"
    )
    print(
        f"best target-holdout delta      = "
        f"{tc['accuracy'] - tr['accuracy']:+.6f}"
    )
    print(
        f"cyclotomic LOO mean            = "
        f"{mean(loo_c):.6f}"
    )
    print(
        f"control LOO mean               = "
        f"{mean(loo_r):.6f}"
    )

    print()
    print("INTERPRETATION")
    print("-" * 78)

    if (
        tc["accuracy"] >= 0.60
        and tc["accuracy"] - tr["accuracy"] >= 0.08
        and mean(loo_c) >= 0.55
        and mean(loo_c) - mean(loo_r) >= 0.05
    ):
        print("PROMISING:")
        print("A simple N-only symmetric expression shows")
        print("transferable cyclotomic predictive power.")
        print("A subsequent exact CRT reconstruction test is justified.")
    elif (
        tc["accuracy"] >= 0.55
        and mean(loo_c) >= 0.50
    ):
        print("WEAK SIGNAL:")
        print("Some N-only symmetric expression survives target")
        print("holdout, but the evidence is insufficient for an")
        print("algebraic factorization constraint.")
    else:
        print("NO ROBUST SYMMETRIC RESOLVENT SIGNAL:")
        print("The small algebraic family does not recover the")
        print("common branch sign in a transferable way.")
        print("Further ad-hoc polynomial searches should not be")
        print("treated as evidence unless independently justified.")

    print()
    print("Key point:")
    print("The search is deliberately restricted to symmetric")
    print("N-only expressions. It does not assume an orientation")
    print("label is directly observable from n.")
    print()
    print("=" * 78)
    print("EXPERIMENT 63 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

