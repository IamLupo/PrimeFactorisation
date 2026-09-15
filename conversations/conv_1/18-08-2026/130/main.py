#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 130 — EXACT (p,q)-KERNEL TO A/B-TRIANGLE PROVENANCE AUDIT
==============================================================================

Purpose
-------
Establish whether the A/B residual coefficient data studied in Experiments
95-129 is an exact reparameterization / basis transform of the original
symmetric p,q kernel from the paper.

Original kernel variables:
    p0, q0

Elementary symmetric variables:
    N = p0*q0
    S = p0 + q0
    X = S + 1

Original layered coefficients:
    G_{k,l}(N,X)
      = sum_a C[k,l,t,a] N^a X^(l-t-a)

Recent derived data:
    A/B channels -> centered residuals -> q[p,r]

This script searches only exact finite transformations.

IMPORTANT
---------
This file does NOT invent the original kernel.  Put the exact kernel adapter
in `kernel_data.py` next to this script.

Expected interface in kernel_data.py:

    def kernel(k, ell, p, q):
        return Fraction(...)

Optionally:

    def theorem_layer(k, ell, layer):
        # return dict {a: Fraction(...)}
        # where coefficient of N^a X^(ell-layer-a) is returned.

The A/B data are embedded below from the verified experiments.

Everything uses Fraction only.
No floating point.
No SymPy.
No extrapolation.
==============================================================================
"""

from fractions import Fraction
from math import gcd, factorial
from itertools import product
import importlib.util
import os
import sys


# ============================================================================
# 0. EXACT RECENT DATA
# ============================================================================

# q[p][r] = residual quotient in the corrected falling basis.
#
# These are the exact B/A residual-sector rows established in Experiments
# 115/126.  The present experiment primarily uses them as provenance targets.

Q_DATA = {
    "A-even": {
        0: [
            Fraction(-12879, 1),
            Fraction(-15362, 1),
            Fraction(8307, 1),
            Fraction(-748, 1),
            Fraction(-62305, 144),
            Fraction(34070797, 201600),
            Fraction(-102402481, 3628800),
        ],
        2: [
            Fraction(2797337, 11520),
            Fraction(7319861, 19200),
            Fraction(-266207819, 1612800),
            Fraction(-141551327, 14515200),
            Fraction(2461903867, 116121600),
            Fraction(-783039371, 116121600),
        ],
        4: [
            Fraction(-2083937, 921600),
            Fraction(-3930151, 921600),
            Fraction(218532413, 77414400),
            Fraction(238646881, 232243200),
            Fraction(-145404619, 92897280),
        ],
        6: [
            Fraction(85591, 22118400),
            Fraction(249013, 22118400),
            Fraction(-24042833, 619315200),
        ],
        8: [
            Fraction(-4913, 353894400),
        ],
    },

    "A-odd": {
        1: [
            Fraction(12143, 3360),
            Fraction(76427, 192),
            Fraction(1954873, 17920),
            Fraction(-61469491, 483840),
            Fraction(42852113, 1935360),
            Fraction(4384549, 645120),
        ],
        3: [
            Fraction(-989, 345600),
            Fraction(-36064769, 2764800),
            Fraction(-24988097, 2580480),
            Fraction(87300373, 19353600),
            Fraction(116226679, 77414400),
        ],
        5: [
            Fraction(-517, 2764800),
            Fraction(1189391, 5529600),
            Fraction(22183547, 103219200),
            Fraction(-69294643, 185794560),
        ],
        7: [
            Fraction(373, 928972800),
            Fraction(-616981, 371589120),
        ],
    },

    "B-even": {
        0: [
            Fraction(12980463, 1024),
            Fraction(12041869, 1024),
            Fraction(-594517923, 71680),
            Fraction(112271581, 71680),
            Fraction(123350021, 860160),
            Fraction(-1803389011, 12902400),
        ],
        2: [
            Fraction(-19344659, 76800),
            Fraction(-38450509, 115200),
            Fraction(490918171, 3225600),
            Fraction(-31781287, 1209600),
            Fraction(254480207, 12902400),
        ],
        4: [
            Fraction(25883, 12288),
            Fraction(682871, 184320),
            Fraction(-1505893, 2580480),
            Fraction(-25483559, 7741440),
        ],
        6: [
            Fraction(-2267, 614400),
            Fraction(-100657, 5529600),
        ],
    },

    "B-odd": {
        1: [
            Fraction(-584531, 35840),
            Fraction(-32224291, 21504),
            Fraction(74131151, 215040),
            Fraction(29406229, 129024),
            Fraction(-1338089411, 7741440),
            Fraction(495451247, 7741440),
        ],
        3: [
            Fraction(-59257, 230400),
            Fraction(7521137, 230400),
            Fraction(12697441, 3225600),
            Fraction(4595257, 1382400),
            Fraction(-140504813, 12902400),
        ],
        5: [
            Fraction(4457, 2764800),
            Fraction(-340837, 2764800),
            Fraction(-5342627, 12902400),
        ],
        7: [
            Fraction(-421, 38707200),
        ],
    },
}


CHANNEL_INFO = {
    "A-even": {"channel": 0, "parity": 0, "m": 8, "K": 6},
    "A-odd":  {"channel": 0, "parity": 1, "m": 8, "K": 6},
    "B-even": {"channel": 1, "parity": 0, "m": 7, "K": 5},
    "B-odd":  {"channel": 1, "parity": 1, "m": 7, "K": 5},
}


# ============================================================================
# 1. OPTIONAL ORIGINAL-KERNEL ADAPTER
# ============================================================================

def load_kernel_adapter():
    """
    Load kernel_data.py from the same directory.

    Required:
        def kernel(k, ell, p, q) -> Fraction/int

    Optional:
        def theorem_layer(k, ell, layer) -> dict[int, Fraction]

    No assumptions are made about implementation.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "kernel_data.py")

    if not os.path.exists(path):
        return None

    spec = importlib.util.spec_from_file_location("kernel_data", path)
    if spec is None or spec.loader is None:
        return None

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "kernel"):
        raise RuntimeError(
            "kernel_data.py exists but does not define kernel(k, ell, p, q)."
        )

    return module


# ============================================================================
# 2. SMALL EXACT POLYNOMIAL UTILITIES
# ============================================================================

def falling(n, r):
    """Exact falling factorial n_(r)."""
    if r < 0:
        return Fraction(0)
    result = 1
    for i in range(r):
        result *= n - i
    return result


def rising(n, r):
    """Exact rising factorial."""
    if r < 0:
        return Fraction(0)
    result = 1
    for i in range(r):
        result *= n + i
    return result


def binomial(n, r):
    if r < 0 or r > n:
        return 0
    if n < 0:
        raise ValueError("binomial requires n >= 0")
    return factorial(n) // (factorial(r) * factorial(n - r))


def primitive_signature(values):
    """
    Convert rational coefficients to a primitive integer vector up to sign.
    """
    nums = [v.numerator for v in values]
    dens = [v.denominator for v in values]

    if not values:
        return []

    lcm_den = 1
    for d in dens:
        lcm_den = lcm_den * d // gcd(lcm_den, d)

    ints = [n * (lcm_den // d) for n, d in zip(nums, dens)]

    g = 0
    for x in ints:
        g = gcd(g, abs(x))

    if g:
        ints = [x // g for x in ints]

    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-y for y in ints]
            break

    return ints


def coefficient_vector_from_samples(values, basis="falling"):
    """
    Reconstruct a polynomial in k from exact values k=0,...,K.

    Returned list is in the requested basis.

    Only the falling basis is needed by the current audit.
    """
    vals = [Fraction(v) for v in values]
    K = len(vals) - 1

    if basis != "falling":
        raise ValueError("Only falling basis is implemented here.")

    # Newton forward differences:
    # f(k) = sum_r Delta^r f(0)/r! * k_(r)
    rows = []
    current = vals[:]

    while current:
        rows.append(current[0])
        if len(current) == 1:
            break
        current = [
            current[i + 1] - current[i]
            for i in range(len(current) - 1)
        ]

    coeffs = [
        Fraction(rows[r], factorial(r))
        for r in range(K + 1)
    ]
    return coeffs


# ============================================================================
# 3. ORIGINAL KERNEL -> SYMMETRIC VARIABLES
# ============================================================================

def kernel_symmetry_test(kernel, k, ell, samples):
    """
    Check exact symmetry F(p,q)=F(q,p).
    """
    for p, q in samples:
        a = Fraction(kernel(k, ell, p, q))
        b = Fraction(kernel(k, ell, q, p))
        if a != b:
            return False, (p, q, a, b)
    return True, None


def reconstruct_symmetric_coefficients(kernel, k, ell, max_total_degree):
    """
    Recover G(N,X) coefficients from exact evaluations.

    Since X=p+q+1 and N=pq, we solve in the monomial basis

        N^a X^b,  a+b <= max_total_degree.

    We deliberately use a square exact linear system with Fraction Gaussian
    elimination rather than relying on floating-point interpolation.
    """

    monomials = []
    for total in range(max_total_degree + 1):
        for a in range(total + 1):
            b = total - a
            monomials.append((a, b))

    # Choose deterministic sample points.
    points = []
    limit = max(4, max_total_degree + 3)

    for p in range(0, limit + 1):
        for q in range(0, limit + 1):
            points.append((p, q))
            if len(points) >= len(monomials):
                break
        if len(points) >= len(monomials):
            break

    A = []
    bvec = []

    for p, q in points:
        N = p * q
        X = p + q + 1

        row = [
            Fraction(N ** a) * Fraction(X ** b)
            for a, b in monomials
        ]
        A.append(row)
        bvec.append(Fraction(kernel(k, ell, p, q)))

    solution = solve_exact(A, bvec)

    if solution is None:
        return None

    return {
        monomials[i]: solution[i]
        for i in range(len(monomials))
        if solution[i] != 0
    }


# ============================================================================
# 4. EXACT LINEAR ALGEBRA
# ============================================================================

def rref_solve(A, b):
    """
    Return one exact solution for A x = b if unique.
    Otherwise return None.
    """
    A = [[Fraction(x) for x in row] for row in A]
    b = [Fraction(x) for x in b]

    n = len(A)
    if n == 0:
        return []

    m = len(A[0])
    M = [A[i] + [b[i]] for i in range(n)]

    row = 0
    pivots = []

    for col in range(m):
        pivot = None
        for r in range(row, n):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[row], M[pivot] = M[pivot], M[row]

        pv = M[row][col]
        M[row] = [x / pv for x in M[row]]

        for r in range(n):
            if r == row:
                continue
            if M[r][col] == 0:
                continue

            f = M[r][col]
            M[r] = [
                M[r][c] - f * M[row][c]
                for c in range(m + 1)
            ]

        pivots.append(col)
        row += 1

        if row == n:
            break

    # inconsistency
    for r in range(n):
        if all(M[r][c] == 0 for c in range(m)) and M[r][m] != 0:
            return None

    if len(pivots) != m:
        return None

    x = [Fraction(0) for _ in range(m)]
    for r, c in enumerate(pivots):
        x[c] = M[r][m]
    return x


def solve_exact(A, b):
    return rref_solve(A, b)


def matrix_rank(A):
    if not A:
        return 0

    M = [[Fraction(x) for x in row] for row in A]
    rows = len(M)
    cols = len(M[0])

    rank = 0
    for col in range(cols):
        pivot = None
        for r in range(rank, rows):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = M[pivot], M[rank]

        pv = M[rank][col]
        M[rank] = [x / pv for x in M[rank]]

        for r in range(rows):
            if r != rank and M[r][col] != 0:
                f = M[r][col]
                M[r] = [
                    M[r][c] - f * M[rank][c]
                    for c in range(cols)
                ]

        rank += 1

    return rank


# ============================================================================
# 5. HOMOGENEOUS LAYER EXTRACTION
# ============================================================================

def extract_layer(coeffs, total_degree):
    """
    Given {(a,b): c}, return the homogeneous layer a+b=total_degree.
    """
    return {
        (a, b): c
        for (a, b), c in coeffs.items()
        if a + b == total_degree and c != 0
    }


def layer_vector(coeffs, total_degree):
    """
    Return coefficients by N-power:
        [N^0 X^D, N^1 X^(D-1), ...].
    """
    result = []
    for a in range(total_degree + 1):
        result.append(Fraction(coeffs.get((a, total_degree - a), 0)))
    return result


# ============================================================================
# 6. SEARCH FOR EXACT DICTIONARY TO Q-DATA
# ============================================================================

def evaluate_layer_vector_at(layer, degree):
    return layer_vector(layer, degree)


def candidate_transforms(layer_vec, channel_info, degree):
    """
    Generate only simple exact transforms.

    Candidate coordinate choices:
      * forward order
      * reverse order
      * signed reverse
      * falling transform
      * shifted index p=a+t

    No nonlinear fitting is used.
    """
    out = []

    # raw
    out.append(("raw", list(layer_vec)))

    # reversed
    out.append(("reverse", list(reversed(layer_vec))))

    # signed reverse
    out.append(("signed_reverse", [
        x if (i % 2 == 0) else -x
        for i, x in enumerate(reversed(layer_vec))
    ]))

    # falling transform of the coefficient list
    try:
        f = coefficient_vector_from_samples(layer_vec, basis="falling")
        out.append(("falling", f))
    except Exception:
        pass

    return out


def compare_exact(a, b):
    if len(a) != len(b):
        return False
    return all(Fraction(x) == Fraction(y) for x, y in zip(a, b))


def search_q_match(layer_vec, target, label):
    """
    Search exact structural matches between an original layer coefficient
    vector and a target Q-row.

    The target may have a lower degree because endpoint factors have already
    been removed.  We therefore also test exact truncation by leading or
    trailing support.
    """
    matches = []

    variants = candidate_transforms(layer_vec, CHANNEL_INFO[label],
                                     len(layer_vec) - 1)

    for name, vec in variants:
        if compare_exact(vec, target):
            matches.append(name)

        if len(vec) >= len(target):
            if compare_exact(vec[:len(target)], target):
                matches.append(name + "_prefix")

            if compare_exact(vec[-len(target):], target):
                matches.append(name + "_suffix")

    return matches


# ============================================================================
# 7. ORIGINAL-KERNEL PROVENANCE TEST
# ============================================================================

def provenance_test(kernel_module):
    """
    Attempt the following exact dictionary:

    original:
        F(k,ell,p,q)
             -> G(N,X)
             -> homogeneous layer
             -> coefficient vector

    derived:
        Q_DATA[channel][power]

    We test several small fixed values of k, ell supplied below.

    This function never declares success merely because dimensions agree.
    It requires exact coefficient equality.
    """
    if kernel_module is None:
        return {
            "status": "MISSING_KERNEL_ADAPTER",
            "matches": [],
        }

    kernel = kernel_module.kernel
    matches = []

    # These are deliberately small finite tests.
    # Change/add values once the kernel adapter is known.
    test_triples = [
        (1, 10),
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 19),
        (9, 21),
    ]

    for k, ell in test_triples:
        # We do not know the exact maximum total degree from this transcript.
        # Search a conservative finite range.
        for degree in range(0, ell + 1):
            coeffs = reconstruct_symmetric_coefficients(
                kernel,
                k,
                ell,
                degree,
            )

            if coeffs is None:
                continue

            layer = extract_layer(coeffs, degree)
            vec = layer_vector(layer, degree)

            for label, qrows in Q_DATA.items():
                for power, qrow in qrows.items():
                    target_degree = len(qrow) - 1

                    if target_degree != len(vec) - 1:
                        continue

                    names = search_q_match(
                        vec,
                        qrow,
                        label,
                    )

                    for name in names:
                        matches.append(
                            {
                                "k": k,
                                "ell": ell,
                                "degree": degree,
                                "label": label,
                                "power": power,
                                "transform": name,
                            }
                        )

    return {
        "status": "COMPLETE",
        "matches": matches,
    }


# ============================================================================
# 8. INDEX-DICTIONARY TEST WITHOUT CLAIMING A MATCH
# ============================================================================

def report_recent_structure():
    print("KNOWN DERIVED OBJECTS")
    print("-" * 78)

    for label, rows in Q_DATA.items():
        info = CHANNEL_INFO[label]

        print(
            f"  {label}: channel={info['channel']} "
            f"parity={info['parity']} m={info['m']} K={info['K']}"
        )

        for power, coeffs in rows.items():
            degree = len(coeffs) - 1
            support = list(range(len(coeffs)))

            print(
                f"    power={power}: "
                f"degree={degree} support={support}"
            )

        print()


# ============================================================================
# 9. SAFE OPTIONAL THEOREM-LAYER ADAPTER
# ============================================================================

def direct_theorem_layer_test(kernel_module):
    """
    If kernel_data.py exposes theorem_layer(), use it directly.

    This is much stronger than trying to reconstruct it indirectly.
    """
    if kernel_module is None or not hasattr(kernel_module, "theorem_layer"):
        return None

    results = []

    samples = [
        (1, 10),
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 19),
        (9, 21),
    ]

    for k, ell in samples:
        for layer in range(0, 5):
            try:
                d = kernel_module.theorem_layer(k, ell, layer)
            except Exception:
                continue

            if not isinstance(d, dict):
                continue

            items = sorted(
                (int(a), Fraction(v))
                for a, v in d.items()
                if Fraction(v) != 0
            )

            results.append((k, ell, layer, items))

    return results


# ============================================================================
# 10. MAIN
# ============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 130 — EXACT (p,q)-KERNEL TO A/B-TRIANGLE PROVENANCE AUDIT")
    print("=" * 78)
    print()

    kernel_module = load_kernel_adapter()

    # ------------------------------------------------------------------------
    # 1. Derived-data validation
    # ------------------------------------------------------------------------
    print("=" * 78)
    print("1. DERIVED A/B DATA VALIDATION")
    print("=" * 78)

    derived_ok = True

    for label, rows in Q_DATA.items():
        for power, row in rows.items():
            ok = (
                len(row) >= 1
                and all(isinstance(x, Fraction) for x in row)
            )
            print(
                f"  {label}, power={power}: "
                f"entries={len(row)} exact={ok}"
            )
            derived_ok = derived_ok and ok

    print()
    print(f"  derived_data_exact={derived_ok}")
    print()

    # ------------------------------------------------------------------------
    # 2. Original-kernel availability
    # ------------------------------------------------------------------------
    print("=" * 78)
    print("2. ORIGINAL KERNEL AVAILABILITY")
    print("=" * 78)

    kernel_available = kernel_module is not None

    print(f"  kernel_adapter_present={kernel_available}")

    if kernel_available:
        print("  kernel_data.py loaded successfully")
        print("  required interface: kernel(k, ell, p, q)")
    else:
        print(
            "  kernel_data.py NOT FOUND"
            "\n  Provenance cannot be claimed without the original kernel."
        )

    print()

    # ------------------------------------------------------------------------
    # 3. Symmetry test if available
    # ------------------------------------------------------------------------
    symmetry_ok = None

    print("=" * 78)
    print("3. ORIGINAL p,q SYMMETRY AUDIT")
    print("=" * 78)

    if kernel_available:
        kernel = kernel_module.kernel

        symmetry_ok = True
        failures = []

        for k, ell in [(1, 10), (3, 7), (3, 11), (5, 11)]:
            ok, detail = kernel_symmetry_test(
                kernel,
                k,
                ell,
                [(0, 1), (1, 2), (2, 3), (1, 4)],
            )

            print(
                f"  k={k}, ell={ell}: symmetric={ok}"
            )

            if not ok:
                symmetry_ok = False
                failures.append((k, ell, detail))

        if failures:
            print(f"  failures={failures}")
    else:
        print("  SKIPPED — original kernel unavailable")

    print()

    # ------------------------------------------------------------------------
    # 4. Direct theorem-layer adapter, if supplied
    # ------------------------------------------------------------------------
    print("=" * 78)
    print("4. DIRECT HOMOGENEOUS-LAYER ADAPTER")
    print("=" * 78)

    direct_layers = direct_theorem_layer_test(kernel_module)

    if direct_layers is None:
        print("  theorem_layer() not supplied")
    else:
        print(f"  recovered_exact_layer_records={len(direct_layers)}")

        for rec in direct_layers[:12]:
            k, ell, layer, items = rec
            print(
                f"    (k={k}, ell={ell}, layer={layer}) "
                f"support={[a for a, _ in items]}"
            )

    print()

    # ------------------------------------------------------------------------
    # 5. Provenance search
    # ------------------------------------------------------------------------
    print("=" * 78)
    print("5. EXACT PROVENANCE SEARCH")
    print("=" * 78)

    prov = provenance_test(kernel_module)

    print(f"  status={prov['status']}")

    if prov["matches"]:
        print(f"  exact_match_count={len(prov['matches'])}")

        for match in prov["matches"][:50]:
            print(
                "  MATCH "
                f"k={match['k']} "
                f"ell={match['ell']} "
                f"degree={match['degree']} "
                f"sector={match['label']} "
                f"power={match['power']} "
                f"transform={match['transform']}"
            )
    else:
        print("  no exact dictionary match found in tested family")

    print()

    # ------------------------------------------------------------------------
    # 6. No-false-positive audit
    # ------------------------------------------------------------------------
    print("=" * 78)
    print("6. NO-FALSE-POSITIVE POLICY")
    print("=" * 78)

    provenance_established = (
        kernel_available
        and prov["status"] == "COMPLETE"
        and len(prov["matches"]) > 0
    )

    print(
        "  provenance_established="
        + str(provenance_established)
    )

    if not provenance_established:
        print(
            "  IMPORTANT: the A/B residual operator is NOT identified "
            "with the original p,q kernel."
        )
    else:
        print(
            "  IMPORTANT: exact pointwise matches were found; "
            "inspect the dictionary above before merging theories."
        )

    print()

    # ------------------------------------------------------------------------
    # 7. Structural report
    # ------------------------------------------------------------------------
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The original paper begins with

    n = p q,
    N = p q,
    X = p + q + 1,

and the exact symmetric kernel

    G_{k,ell}(N,X).

Experiments 98-129 instead study derived coefficient objects
with endpoint factors, centered parity sectors, and residual
falling-basis arrays q[p,r].

This experiment tests the missing bridge:

    original p,q kernel
        |
        v
    symmetric (N,X) layers
        |
        v
    coefficient arrays
        |
        ? 
        v
    A/B residual q[p,r].

The bridge is accepted only when exact coefficient equality is
found after an explicitly stated finite transform.

Dimension agreement, similar support, large rational coefficients,
or a visually similar recurrence are NOT accepted as evidence.

Therefore the current output distinguishes:

    PROVENANCE ESTABLISHED
        from
    STRUCTURAL COMPATIBILITY ONLY.

This is the critical distinction needed before inserting
Experiments 95-129 into the main p,q-kernel paper.
"""
    )

    print()

    # ------------------------------------------------------------------------
    # 8. Final status
    # ------------------------------------------------------------------------
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = derived_ok

    if kernel_available:
        final_ok = final_ok and (
            prov["status"] == "COMPLETE"
        )

    print(f"  derived_data_exact={derived_ok}")
    print(f"  kernel_available={kernel_available}")
    print(f"  provenance_established={provenance_established}")
    print(f"  failures={0 if final_ok else 1}")

    if final_ok:
        print("  ALL BASIC CHECKS PASS=True")
    else:
        print("  ALL BASIC CHECKS PASS=False")

    print()
    if provenance_established:
        print("EXPERIMENT 130 COMPLETE — EXACT BRIDGE CANDIDATE FOUND")
    else:
        print("EXPERIMENT 130 COMPLETE — BRIDGE NOT YET ESTABLISHED")


if __name__ == "__main__":
    main()

