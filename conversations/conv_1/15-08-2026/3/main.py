#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 61R
N-ONLY CONJUGATE-BRANCH ORIENTATION / CRT RECONSTRUCTION TEST
PATCHED CONTROL-FAMILY HANDLING
TARGET-LEVEL OUT-OF-SAMPLE VALIDATION
NO CSV OUTPUT
==============================================================================

PURPOSE
-------
Experiment 56:
    G(n) determines the PRODUCT of the two conjugate Legendre symbols,
    but not the orientation.

Experiment 60:
    CRT combines local congruences but does not create new information.

Experiment 61R:
    Test whether an N-only feature set can predict WHICH conjugate root
    carries the +1 Legendre symbol.

NO predictor feature uses:
    p, q, s, D(true_s), true local root choice, or factor information.

IMPORTANT CONTROL FIX
---------------------
The orientation label requires x^2+x+1 to have TWO roots modulo ell.

That occurs for primes ell == 1 (mod 3).

Therefore controls that are not 1 mod 3 are excluded from the control
orientation comparison rather than being treated as if they had roots.

For cubic character:
    if ell == 1 mod 3:
        three cubic classes exist;
    otherwise:
        every nonzero element has trivial cubic character,
        so class 1 is used and class 0 denotes zero.

The oracle CRT reconstruction remains a diagnostic only:
the unordered local pair is supplied there solely to test whether predicting
the missing orientation would be sufficient to recover the true CRT residue.
==============================================================================

"""

from __future__ import annotations

import math
import random
import statistics
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# ============================================================================
# CONFIG
# ============================================================================

SEED = 20260815

S_MIN = 4_000_000
S_MAX = 8_399_998
STEP = 2

TARGET_COUNT = 24

TRAIN_TARGETS = 16
TEST_TARGETS = TARGET_COUNT - TRAIN_TARGETS

CYCLO_MODULI = [
    7, 13, 19, 31, 37, 61, 67, 79, 127,
    307, 331, 631, 1723
]

CONTROL_MODULI_RAW = [
    673, 1109, 2309, 2351,
    4421, 4561, 4759, 6211,
    7879, 7951, 8273, 8689, 9781
]


# ============================================================================
# TARGETS
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q


TARGETS = [
    Target(3318013, 4042603),
    Target(2129167, 3402323),
    Target(2224517, 3978749),
    Target(3685051, 4020281),
    Target(2399627, 2452649),
    Target(2593039, 2996527),
    Target(2149859, 2772097),
    Target(2060543, 2514401),
    Target(2675423, 2883973),
    Target(2828887, 3960137),
    Target(3497381, 3793241),
    Target(2193509, 4011353),
    Target(3429689, 3983927),
    Target(2515871, 3050581),
    Target(2261297, 3374827),
    Target(2345537, 3673349),
    Target(2212039, 3452809),
    Target(2175373, 2476921),
    Target(2438509, 4188577),
    Target(2367553, 2399407),
    Target(2072897, 2087077),
    Target(3552023, 3707453),
    Target(3112909, 3210167),
    Target(2965819, 3239963),
]

assert len(TARGETS) == TARGET_COUNT


# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    a %= p

    if a == 0:
        return 0

    z = pow(a, (p - 1) // 2, p)

    if z == 1:
        return 1
    if z == p - 1:
        return -1

    raise RuntimeError(
        f"Invalid Legendre result: a={a}, p={p}, z={z}"
    )


def roots_cyclotomic(p: int) -> Tuple[int, int]:
    """
    Roots of x^2+x+1 modulo p.

    For prime p this requires:
        p == 1 (mod 3)
    except the special degenerate p=3 case, which is not used here.
    """
    if p == 3:
        raise RuntimeError("p=3 is degenerate for this experiment")

    if p % 3 != 1:
        raise RuntimeError(
            f"{p} is not 1 mod 3, so x^2+x+1 has no two nontrivial roots"
        )

    roots = []

    # Find a nontrivial cube root of unity:
    # x^3 = 1 and x != 1.
    zeta = pow(primitive_generator(p), (p - 1) // 3, p)

    if zeta != 1:
        roots = [zeta, (-zeta - 1) % p]
    else:
        # Extremely defensive fallback.
        for x in range(2, p):
            if (x * x + x + 1) % p == 0:
                roots.append(x)
                break

        if not roots:
            raise RuntimeError(
                f"Could not find root of x^2+x+1 modulo {p}"
            )

        roots.append((-roots[0] - 1) % p)

    roots = sorted(set(roots))

    if len(roots) != 2:
        raise RuntimeError(
            f"Expected exactly two roots modulo {p}, got {roots}"
        )

    return roots[0], roots[1]


def primitive_generator(p: int) -> int:
    """
    Find a primitive root modulo a prime p.

    The moduli here are small enough that trial search is inexpensive.
    """
    if p == 2:
        return 1

    phi = p - 1

    factors = []
    x = phi

    f = 2
    while f * f <= x:
        if x % f == 0:
            factors.append(f)
            while x % f == 0:
                x //= f
        f += 1

    if x > 1:
        factors.append(x)

    for g in range(2, p):
        ok = True

        for r in factors:
            if pow(g, phi // r, p) == 1:
                ok = False
                break

        if ok:
            return g

    raise RuntimeError(f"No primitive root found for {p}")


def has_cyclotomic_roots(p: int) -> bool:
    return p != 3 and p % 3 == 1


# ============================================================================
# CUBIC CHARACTER
# ============================================================================

def cubic_class(a: int, p: int) -> int:
    """
    Return cubic character class.

    For p == 1 mod 3:
        0 = zero
        1,2,3 = the three cubic classes

    For p != 1 mod 3:
        0 = zero
        1 = trivial nonzero class

    This allows the same N-only feature encoder to operate on controls
    without falsely requiring x^2+x+1 roots.
    """
    a %= p

    if a == 0:
        return 0

    if p % 3 != 1:
        return 1

    z = pow(a, (p - 1) // 3, p)

    if z == 1:
        return 1

    r1, r2 = roots_cyclotomic(p)

    if z == r1:
        return 2

    if z == r2:
        return 3

    raise RuntimeError(
        f"Unexpected cubic character: p={p}, a={a}, z={z}"
    )


# ============================================================================
# LOCAL STATE
# ============================================================================

@dataclass(frozen=True)
class LocalState:
    ell: int
    root1: int
    root2: int

    d1: int
    d2: int

    chi1: int
    chi2: int

    chi_n: int
    chi_n1: int
    chi_nminus1: int
    chi_4n1: int
    chi_G: int

    cubic_n: int
    cubic_n1: int
    cubic_4n1: int


def local_state_with_roots(n: int, ell: int) -> LocalState:
    """
    Full orientation-capable local state.

    Only call this for ell == 1 mod 3.
    """
    if not has_cyclotomic_roots(ell):
        raise ValueError(
            f"Orientation state requested for invalid ell={ell}"
        )

    r1, r2 = roots_cyclotomic(ell)

    d1 = (r1 * r1 - 4 * n) % ell
    d2 = (r2 * r2 - 4 * n) % ell

    G = 16 * n * n + 4 * n + 1

    return LocalState(
        ell=ell,
        root1=r1,
        root2=r2,
        d1=d1,
        d2=d2,
        chi1=legendre_symbol(d1, ell),
        chi2=legendre_symbol(d2, ell),

        chi_n=legendre_symbol(n, ell),
        chi_n1=legendre_symbol(n + 1, ell),
        chi_nminus1=legendre_symbol(n - 1, ell),
        chi_4n1=legendre_symbol(4 * n + 1, ell),
        chi_G=legendre_symbol(G, ell),

        cubic_n=cubic_class(n, ell),
        cubic_n1=cubic_class(n + 1, ell),
        cubic_4n1=cubic_class(4 * n + 1, ell),
    )


# ============================================================================
# FEATURES
# ============================================================================

FEATURE_NAMES = [
    "chi(n)",
    "chi(n+1)",
    "chi(n-1)",
    "chi(4n+1)",
    "chi(G)",
    "cubic_n==1",
    "cubic_n==2",
    "cubic_n==3",
    "cubic_n1==1",
    "cubic_n1==2",
    "cubic_n1==3",
    "cubic_4n1==1",
    "cubic_4n1==2",
    "cubic_4n1==3",
]


def encode_features(st: LocalState) -> List[float]:
    return [
        float(st.chi_n),
        float(st.chi_n1),
        float(st.chi_nminus1),
        float(st.chi_4n1),
        float(st.chi_G),

        float(st.cubic_n == 1),
        float(st.cubic_n == 2),
        float(st.cubic_n == 3),

        float(st.cubic_n1 == 1),
        float(st.cubic_n1 == 2),
        float(st.cubic_n1 == 3),

        float(st.cubic_4n1 == 1),
        float(st.cubic_4n1 == 2),
        float(st.cubic_4n1 == 3),
    ]


# ============================================================================
# ORIENTATION LABEL
# ============================================================================

def clean_orientation_label(st: LocalState) -> Optional[int]:
    """
    +1:
        root1 = QR, root2 = non-QR

    -1:
        root1 = non-QR, root2 = QR

    None:
        degenerate/non-oriented local state.
    """
    if st.chi1 == 1 and st.chi2 == -1:
        return 1

    if st.chi1 == -1 and st.chi2 == 1:
        return -1

    return None


# ============================================================================
# LOGISTIC MODEL
# ============================================================================

def sigmoid(z: float) -> float:
    if z >= 0:
        e = math.exp(-z)
        return 1.0 / (1.0 + e)

    e = math.exp(z)
    return e / (1.0 + e)


def train_logistic(
    X: List[List[float]],
    y: List[int],
    epochs: int = 2500,
    lr: float = 0.03,
    l2: float = 0.20,
) -> Tuple[List[float], float]:

    if not X:
        raise RuntimeError("No training samples")

    d = len(X[0])

    w = [0.0] * d
    b = 0.0

    n = len(X)

    for _ in range(epochs):
        gw = [0.0] * d
        gb = 0.0

        for xi, yi in zip(X, y):
            z = b

            for j in range(d):
                z += w[j] * xi[j]

            p = sigmoid(z)
            err = p - yi

            gb += err

            for j in range(d):
                gw[j] += err * xi[j]

        gb /= n

        for j in range(d):
            gw[j] /= n
            gw[j] += l2 * w[j]
            w[j] -= lr * gw[j]

        b -= lr * gb

    return w, b


def predict_probability(
    x: List[float],
    w: List[float],
    b: float,
) -> float:

    z = b

    for a, c in zip(x, w):
        z += a * c

    return sigmoid(z)


# ============================================================================
# FIXED N-ONLY RULES
# ============================================================================

def fixed_rule_predictions(st: LocalState) -> Dict[str, Optional[int]]:

    out: Dict[str, Optional[int]] = {}

    out["chi(n)"] = st.chi_n if st.chi_n != 0 else None
    out["chi(n+1)"] = st.chi_n1 if st.chi_n1 != 0 else None
    out["chi(4n+1)"] = st.chi_4n1 if st.chi_4n1 != 0 else None
    out["chi(G)"] = st.chi_G if st.chi_G != 0 else None

    if st.chi_n != 0 and st.chi_G != 0:
        out["chi(n)*chi(G)"] = st.chi_n * st.chi_G
    else:
        out["chi(n)*chi(G)"] = None

    if st.chi_n1 != 0 and st.chi_G != 0:
        out["chi(n+1)*chi(G)"] = st.chi_n1 * st.chi_G
    else:
        out["chi(n+1)*chi(G)"] = None

    return out


# ============================================================================
# CRT
# ============================================================================

def crt_pairwise(residues: List[Tuple[int, int]]) -> int:

    x = 0
    M = 1

    for ai, mi in residues:
        inv = pow(M, -1, mi)
        t = ((ai - x) * inv) % mi

        x += M * t
        M *= mi
        x %= M

    return x


# ============================================================================
# METRICS
# ============================================================================

def accuracy(labels: List[int], preds: List[int]) -> float:
    if not labels:
        return float("nan")

    return sum(
        int(a == b)
        for a, b in zip(labels, preds)
    ) / len(labels)


def median_or_nan(values: List[float]) -> float:
    if not values:
        return float("nan")

    return statistics.median(values)


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 61R")
    print("N-ONLY CONJUGATE-BRANCH ORIENTATION / CRT RECONSTRUCTION")
    print("PATCHED CONTROL-FAMILY HANDLING")
    print("TARGET-LEVEL OUT-OF-SAMPLE VALIDATION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. TARGETS
    # ----------------------------------------------------------------------

    print()
    print("1. TARGETS")
    print("-" * 78)

    for i, t in enumerate(TARGETS, 1):
        print(
            f"target {i:3d}: "
            f"p={t.p} q={t.q} "
            f"n={t.n} s={t.s}"
        )

    # ----------------------------------------------------------------------
    # 2. CONTROL VALIDITY
    # ----------------------------------------------------------------------

    control_valid = [
        ell
        for ell in CONTROL_MODULI_RAW
        if has_cyclotomic_roots(ell)
    ]

    control_excluded = [
        ell
        for ell in CONTROL_MODULI_RAW
        if not has_cyclotomic_roots(ell)
    ]

    print()
    print("2. MODULUS VALIDATION")
    print("-" * 78)

    print(f"cyclotomic moduli = {CYCLO_MODULI}")
    print(f"raw controls      = {CONTROL_MODULI_RAW}")

    print(
        "orientation-valid controls "
        f"(ell == 1 mod 3) = {control_valid}"
    )

    print(
        f"excluded controls  = {control_excluded}"
    )

    print(
        "Reason: x^2+x+1 has two nontrivial roots only for "
        "primes ell == 1 mod 3."
    )

    if not control_valid:
        raise RuntimeError(
            "No orientation-valid control moduli remain."
        )

    # ----------------------------------------------------------------------
    # 3. TRAIN/TEST TARGET SPLIT
    # ----------------------------------------------------------------------

    all_indices = list(range(TARGET_COUNT))
    rng.shuffle(all_indices)

    train_indices = sorted(all_indices[:TRAIN_TARGETS])
    test_indices = sorted(all_indices[TRAIN_TARGETS:])

    print()
    print("3. TARGET-LEVEL SPLIT")
    print("-" * 78)

    print(f"training targets = {TRAIN_TARGETS}")
    print(f"test targets     = {TEST_TARGETS}")

    print(
        f"train indices    = {[i + 1 for i in train_indices]}"
    )

    print(
        f"test indices     = {[i + 1 for i in test_indices]}"
    )

    # ----------------------------------------------------------------------
    # 4. LOCAL STATE
    # ----------------------------------------------------------------------

    print()
    print("4. CYCLOTOMIC LOCAL STATES")
    print("-" * 78)

    all_states: Dict[Tuple[int, int], LocalState] = {}

    for ti, t in enumerate(TARGETS):
        for ell in CYCLO_MODULI:
            all_states[(ti, ell)] = local_state_with_roots(
                t.n,
                ell,
            )

    print(
        f"prepared {TARGET_COUNT * len(CYCLO_MODULI)} "
        f"cyclotomic local states"
    )

    # ----------------------------------------------------------------------
    # 5. TRAINING DATA
    # ----------------------------------------------------------------------

    print()
    print("5. CLEAN ORIENTATION DATASET")
    print("-" * 78)

    train_X: List[List[float]] = []
    train_y: List[int] = []

    cyclo_clean_total = 0
    cyclo_degenerate = 0

    for ti in range(TARGET_COUNT):

        for ell in CYCLO_MODULI:

            st = all_states[(ti, ell)]
            label = clean_orientation_label(st)

            if label is None:
                cyclo_degenerate += 1
                continue

            cyclo_clean_total += 1

            if ti not in train_indices:
                continue

            train_X.append(
                encode_features(st)
            )

            train_y.append(
                1 if label == 1 else 0
            )

    print(
        f"clean cyclotomic samples = {cyclo_clean_total}"
    )

    print(
        f"degenerate cells         = {cyclo_degenerate}"
    )

    print(
        f"training samples          = {len(train_X)}"
    )

    # ----------------------------------------------------------------------
    # 6. TRAIN
    # ----------------------------------------------------------------------

    print()
    print("6. N-ONLY LOGISTIC MODEL")
    print("-" * 78)

    w, b = train_logistic(
        train_X,
        train_y,
    )

    print("coefficients:")

    for name, coeff in zip(FEATURE_NAMES, w):
        print(
            f"  {name:20s} {coeff:+.8f}"
        )

    print(
        f"  {'bias':20s} {b:+.8f}"
    )

    # ----------------------------------------------------------------------
    # 7. CYCLO HOLDOUT
    # ----------------------------------------------------------------------

    print()
    print("7. CYCLOTOMIC HOLDOUT")
    print("-" * 78)

    test_labels: List[int] = []
    test_preds: List[int] = []

    test_rows = []

    for ti in test_indices:

        for ell in CYCLO_MODULI:

            st = all_states[(ti, ell)]
            label = clean_orientation_label(st)

            if label is None:
                continue

            p_plus = predict_probability(
                encode_features(st),
                w,
                b,
            )

            pred = 1 if p_plus >= 0.5 else -1

            test_labels.append(label)
            test_preds.append(pred)

            test_rows.append(
                (
                    ti,
                    ell,
                    st,
                    label,
                    p_plus,
                    pred,
                )
            )

    cyclo_accuracy = accuracy(
        test_labels,
        test_preds,
    )

    print(
        f"test samples   = {len(test_labels)}"
    )

    print(
        f"accuracy       = {cyclo_accuracy:.6f}"
    )

    # ----------------------------------------------------------------------
    # 8. CONTROL HOLDOUT
    # ----------------------------------------------------------------------

    print()
    print("8. CONTROL HOLDOUT")
    print("-" * 78)

    control_labels: List[int] = []
    control_preds: List[int] = []

    for ti in test_indices:

        t = TARGETS[ti]

        for ell in control_valid:

            # Valid controls have the same two-root algebra,
            # so orientation labels are well-defined.
            st = local_state_with_roots(
                t.n,
                ell,
            )

            label = clean_orientation_label(st)

            if label is None:
                continue

            p_plus = predict_probability(
                encode_features(st),
                w,
                b,
            )

            pred = 1 if p_plus >= 0.5 else -1

            control_labels.append(label)
            control_preds.append(pred)

    control_accuracy = accuracy(
        control_labels,
        control_preds,
    )

    print(
        f"valid control moduli = {control_valid}"
    )

    print(
        f"control samples       = {len(control_labels)}"
    )

    print(
        f"accuracy              = {control_accuracy:.6f}"
    )

    # ----------------------------------------------------------------------
    # 9. PER-TARGET ACCURACY
    # ----------------------------------------------------------------------

    print()
    print("9. PER-TARGET HOLDOUT ACCURACY")
    print("-" * 78)

    target_accuracies = []

    for ti in test_indices:

        rows = [
            r
            for r in test_rows
            if r[0] == ti
        ]

        if not rows:
            print(
                f"target {ti + 1:3d}: "
                "no clean orientation cells"
            )
            continue

        ok = sum(
            1
            for _, _, _, label, _, pred in rows
            if label == pred
        )

        acc = ok / len(rows)

        target_accuracies.append(acc)

        print(
            f"target {ti + 1:3d}: "
            f"clean={len(rows):2d} "
            f"accuracy={acc:.6f}"
        )

    # ----------------------------------------------------------------------
    # 10. FIXED N-ONLY RULES
    # ----------------------------------------------------------------------

    print()
    print("10. FIXED N-ONLY RULES")
    print("-" * 78)

    rule_labels: Dict[str, List[int]] = {}
    rule_preds: Dict[str, List[int]] = {}

    for ti in test_indices:

        for ell in CYCLO_MODULI:

            st = all_states[(ti, ell)]
            label = clean_orientation_label(st)

            if label is None:
                continue

            rules = fixed_rule_predictions(st)

            for name, pred in rules.items():

                if pred is None:
                    continue

                rule_labels.setdefault(
                    name,
                    [],
                ).append(label)

                rule_preds.setdefault(
                    name,
                    [],
                ).append(pred)

    for name in rule_labels:

        print(
            f"{name:20s} "
            f"accuracy={accuracy(rule_labels[name], rule_preds[name]):.6f} "
            f"coverage={len(rule_labels[name]) / max(1, len(test_labels)):.6f}"
        )

    # ----------------------------------------------------------------------
    # 11. ORACLE CRT RECONSTRUCTION
    # ----------------------------------------------------------------------

    print()
    print("11. ORACLE-LOCAL-PAIR CRT RECONSTRUCTION")
    print("-" * 78)

    print(
        "The local root pair is used ONLY as a diagnostic."
    )

    print(
        "The orientation predictor still sees n-derived features only."
    )

    M_all = math.prod(CYCLO_MODULI)

    eligible = 0
    recovered = 0

    for ti in test_indices:

        t = TARGETS[ti]

        rows = [
            r
            for r in test_rows
            if r[0] == ti
        ]

        # Need every cyclotomic modulus to have a clean orientation.
        if len(rows) != len(CYCLO_MODULI):
            print(
                f"target {ti + 1:3d}: "
                f"not CRT-eligible "
                f"(clean={len(rows)}/{len(CYCLO_MODULI)})"
            )
            continue

        eligible += 1

        true_residues = []
        predicted_residues = []

        all_predictions_correct = True

        for (
            _ti,
            ell,
            st,
            label,
            _p_plus,
            pred,
        ) in rows:

            if label != pred:
                all_predictions_correct = False

            true_root = (
                st.root1
                if label == 1
                else st.root2
            )

            predicted_root = (
                st.root1
                if pred == 1
                else st.root2
            )

            true_residues.append(
                (true_root, ell)
            )

            predicted_residues.append(
                (predicted_root, ell)
            )

        true_crt = (
            crt_pairwise(true_residues)
            % M_all
        )

        predicted_crt = (
            crt_pairwise(predicted_residues)
            % M_all
        )

        true_matches = (
            true_crt == (t.s % M_all)
        )

        predicted_matches = (
            predicted_crt == (t.s % M_all)
        )

        if predicted_matches:
            recovered += 1

        print(
            f"target {ti + 1:3d}: "
            f"true_CRT={true_crt} "
            f"pred_CRT={predicted_crt} "
            f"orientation_all_correct={all_predictions_correct} "
            f"pred_recovered={predicted_matches}"
        )

    if eligible:
        crt_rate = recovered / eligible
    else:
        crt_rate = float("nan")

    print()
    print(
        f"CRT-eligible targets = {eligible}"
    )

    print(
        f"predicted recovered   = {recovered}"
    )

    print(
        f"CRT reconstruction rate = {crt_rate:.6f}"
    )

    # ----------------------------------------------------------------------
    # 12. FEATURE RANKING
    # ----------------------------------------------------------------------

    print()
    print("12. FEATURE COEFFICIENT RANKING")
    print("-" * 78)

    ranked = sorted(
        zip(FEATURE_NAMES, w),
        key=lambda z: abs(z[1]),
        reverse=True,
    )

    for name, coeff in ranked:
        print(
            f"{name:20s} {coeff:+.8f}"
        )

    # ----------------------------------------------------------------------
    # 13. FINAL DIAGNOSTIC
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        f"cyclotomic holdout accuracy = "
        f"{cyclo_accuracy:.6f}"
    )

    print(
        f"control holdout accuracy    = "
        f"{control_accuracy:.6f}"
    )

    if target_accuracies:
        print(
            f"median target accuracy      = "
            f"{median_or_nan(target_accuracies):.6f}"
        )

    print(
        f"CRT reconstruction rate     = "
        f"{crt_rate:.6f}"
    )

    print()
    print("INTERPRETATION")
    print("-" * 78)

    print(
        "NO SIGNAL:"
        " cyclotomic holdout ~ control holdout ~ 0.5."
    )

    print(
        "GENERIC SIGNAL:"
        " both families improve similarly."
    )

    print(
        "CYCLOTOMIC SIGNAL:"
        " cyclotomic holdout is materially stronger than controls."
    )

    print(
        "STRONG RESULT:"
        " N-only orientation predictions reconstruct the true CRT"
        " residue on unseen targets."
    )

    print()
    print("CONTROL NOTE")
    print("-" * 78)

    print(
        f"Raw controls excluded from orientation comparison: "
        f"{control_excluded}"
    )

    print(
        "They are excluded because the same two-root orientation label"
        " is not mathematically defined for ell != 1 mod 3."
    )

    print()
    print("CRITICAL LIMITATION")
    print("-" * 78)

    print(
        "The CRT section does not claim that the unordered local roots"
        " are already obtainable from n alone."
    )

    print(
        "It asks whether solving the missing orientation problem would"
        " be sufficient to turn the local information into an exact CRT"
        " reconstruction."
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 61R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()