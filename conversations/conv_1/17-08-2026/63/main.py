#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 128
NEWTON-INDEX ALIGNMENT CERTIFICATE

REFERENCE:
    EXPERIMENT 126 DIRECT POST-BOUNDARY LEADING-LAW CERTIFICATE

GOAL
----

Experiment 127R failed because it assumed that the target Newton index was

    j = d + 1

for every post-boundary row.

This experiment does NOT assume an index mapping.

Instead:

  1. Build the exact quotient.
  2. Convert it exactly to S=p+q, N=pq.
  3. Compute the COMPLETE Newton coefficient tensor for small ell.
  4. For the exact coefficient polynomials reported by Experiment 126,
     search every Newton index j.
  5. Require exact polynomial equality, not merely matching degree/lead.
  6. Record every matching index.

The experiment therefore answers:

    "Which Newton coefficient did Experiment 126 actually certify?"

This is an INDEX-IDENTIFICATION experiment, not a new fitting experiment.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
S, N = sp.symbols("S N")


# ============================================================================
# Exact detector kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def quotient_Q_pq(k, ell):
    """
    Exact division by p+q+1.
    """

    F = kernel_F(k, ell)

    poly_F = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )

    poly_D = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    quotient, remainder = sp.div(poly_F, poly_D)

    rem = sp.expand(remainder.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder ({k},{ell}): "
            f"{sp.factor(rem)}"
        )

    return sp.expand(quotient.as_expr())


# ============================================================================
# Symmetric p,q -> S,N conversion
# ============================================================================

def symmetric_to_NS(expr):
    """
    Substitute q=S-p and reduce modulo

        p^2 - S*p + N.

    Because the input is symmetric, the final remainder must be
    independent of p.
    """

    substituted = sp.expand(expr.subs(q, S - p))

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    poly = sp.Poly(
        substituted,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    rem = sp.rem(poly, modulus).as_expr()
    rem = sp.expand(rem)

    p_coeff = sp.expand(rem.coeff(p, 1))

    if p_coeff != 0:
        raise ArithmeticError(
            "Symmetric reduction retained p-dependence: "
            + str(sp.factor(p_coeff))
        )

    return sp.expand(rem.coeff(p, 0))


# ============================================================================
# Newton basis
# ============================================================================

def newton_basis(ell):
    """
    Return P_0,...,P_{ell-1}

        P_0 = 2
        P_1 = S
        P_j = S P_{j-1} - N P_{j-2}.
    """

    P = [sp.Integer(2)]

    if ell >= 2:
        P.append(S)

    for j in range(2, ell):
        P.append(
            sp.expand(
                S * P[j - 1]
                - N * P[j - 2]
            )
        )

    return P


# ============================================================================
# Complete Newton tensor
# ============================================================================

def full_newton_tensor(Q_NS, ell):
    """
    Compute

        Q = sum_j c_j(N) P_j

    for j=0,...,ell-1.

    P_0=2 is handled explicitly.
    """

    P = newton_basis(ell)

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    coeffs = {}

    # P_j is monic in S for j>=1.
    for j in range(ell - 1, 0, -1):
        c = sp.expand(work.coeff_monomial(S**j))
        coeffs[j] = c

        if c != 0:
            sub_poly = sp.Poly(
                sp.expand(c * P[j]),
                S,
                domain=sp.QQ.frac_field(N),
            )
            work = work - sub_poly

    # Remaining part must be the P0 contribution.
    residual = sp.expand(work.as_expr())

    residual_poly = sp.Poly(
        residual,
        S,
        domain=sp.QQ.frac_field(N),
    )

    if residual_poly.degree() > 0:
        raise ArithmeticError(
            "Newton decomposition retained positive S-degree: "
            + str(sp.factor(residual))
        )

    coeffs[0] = sp.expand(residual / 2)

    # Exact reconstruction.
    recon = sp.expand(
        sum(
            coeffs[j] * P[j]
            for j in range(ell)
        )
    )

    resid = sp.expand(
        sp.cancel(Q_NS - recon)
    )

    if resid != 0:
        raise ArithmeticError(
            "Newton reconstruction failed: "
            + str(sp.factor(resid))
        )

    return coeffs


# ============================================================================
# Cache
# ============================================================================

def build_tensor_cache(max_ell):
    cache = {}

    print("=" * 78)
    print("BUILDING SMALL EXACT NEWTON TENSOR CACHE")
    print("=" * 78)

    for ell in range(3, max_ell + 1, 2):
        print(f"ell={ell}")

        for k in range(1, ell, 2):
            Qpq = quotient_Q_pq(k, ell)
            QNS = symmetric_to_NS(Qpq)
            coeffs = full_newton_tensor(QNS, ell)

            cache[(k, ell)] = {
                "Qpq": Qpq,
                "QNS": QNS,
                "coeffs": coeffs,
            }

    print()
    print(
        f"cache entries = {len(cache)}"
    )
    print()

    return cache


# ============================================================================
# Experiment-126 exact coefficient fingerprints
# ============================================================================
#
# These are copied from the successful Experiment-126 certificate.
#
# The purpose is NOT to refit them.
# We search the complete Newton tensor for exact polynomial equality.
# ============================================================================

REFERENCE = {
    (1, 5, 1):
        6*N - 1,

    (3, 7, 1):
        2*N*(5*N**2 - 5*N + 3),

    (5, 9, 1):
        7*N**3*(2*N**2 - 5*N + 8),

    (1, 7, 2):
        12*N + 1,

    (3, 9, 2):
        18*N**3 + 15*N**2 - 7*N + 1,

    (5, 11, 2):
        N**2*(24*N**3 + 70*N**2 - 84*N + 45),

    (1, 9, 3):
        -10*N**2 + 60*N - 1,

    (3, 11, 4):
        -21*N**4 + 175*N**3 + 28*N**2 - 9*N + 1,

    (5, 13, 5):
        N**3 * (
            18*N**4
            - 249*N**3
            + 1044*N**2
            - 330*N
            + 220
        ),

    (1, 13, 6):
        35*N**3 - 315*N**2 + 917*N + 1,

    (3, 17, 8):
        -50*N**6
        + 875*N**5
        - 4795*N**4
        + 12705*N**3
        + 66*N**2
        - 13*N
        + 1,

    (1, 23, 10):
        132*N**5
        - 3806*N**4
        + 38676*N**3
        - 203445*N**2
        + 646635*N
        + 1,
}


# ============================================================================
# Exact polynomial helper
# ============================================================================

def exact_equal(a, b):
    return sp.expand(a - b) == 0


def degree_N(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ,
    )

    d = poly.degree()

    if d == -sp.oo:
        return None

    return int(d)


def lead_N(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return sp.Integer(0)

    return sp.expand(poly.LC())


# ============================================================================
# Search reference coefficient through all j
# ============================================================================

def search_reference(cache, key, target):
    k, ell, d = key

    coeffs = cache[(k, ell)]["coeffs"]

    exact_matches = []

    for j in sorted(coeffs):
        c = sp.expand(coeffs[j])

        if exact_equal(c, target):
            exact_matches.append(j)

    return exact_matches


# ============================================================================
# Print complete local tensor for unmatched rows
# ============================================================================

def print_tensor_row(cache, key):
    k, ell, d = key

    print()
    print(
        f"LOCAL NEWTON TENSOR ({k},{ell}), d={d}"
    )
    print("-" * 78)

    coeffs = cache[(k, ell)]["coeffs"]

    for j in sorted(
        coeffs,
        reverse=True,
    ):
        c = sp.factor(coeffs[j])

        print(
            f"j={j:2d} "
            f"degN={degree_N(c)!s:>3} "
            f"lead={lead_N(c)!s:>8} "
            f"c_j={c}"
        )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 128")
    print("NEWTON-INDEX ALIGNMENT CERTIFICATE")
    print("EXPERIMENT 126 REFERENCE COEFFICIENTS")
    print("EXACT FULL TENSOR ON SMALL ELL")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    MAX_ELL = 23

    # ----------------------------------------------------------------------
    # 1. Build exact cache.
    # ----------------------------------------------------------------------

    cache = build_tensor_cache(MAX_ELL)

    # ----------------------------------------------------------------------
    # 2. Verify reference rows are available.
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("REFERENCE-ROW AVAILABILITY")
    print("=" * 78)

    missing = []

    for key in REFERENCE:
        k, ell, d = key

        if (k, ell) not in cache:
            missing.append(key)

    print(
        f"reference rows = {len(REFERENCE)}"
    )
    print(
        f"missing rows   = {len(missing)}"
    )

    if missing:
        for key in missing:
            print("MISSING:", key)

        raise ArithmeticError(
            "Reference rows missing from cache."
        )

    print("STATUS = PASS")
    print()

    # ----------------------------------------------------------------------
    # 3. Exact index search.
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("EXACT NEWTON-INDEX SEARCH")
    print("=" * 78)

    found = 0
    ambiguous = 0
    unmatched = 0

    mapping = {}

    for key, target in REFERENCE.items():

        matches = search_reference(
            cache,
            key,
            sp.expand(target),
        )

        k, ell, d = key

        print(
            f"({k},{ell}) d={d}"
        )
        print(
            f"  reference = {sp.factor(target)}"
        )
        print(
            f"  degree    = {degree_N(target)}"
        )
        print(
            f"  lead      = {lead_N(target)}"
        )
        print(
            f"  matching j = {matches}"
        )

        mapping[key] = matches

        if len(matches) == 1:
            found += 1
        elif len(matches) > 1:
            ambiguous += 1
        else:
            unmatched += 1

            print_tensor_row(
                cache,
                key,
            )

        print()

    print("=" * 78)
    print("INDEX SEARCH SUMMARY")
    print("=" * 78)

    print(
        f"unique matches   = {found}/{len(REFERENCE)}"
    )
    print(
        f"ambiguous matches = {ambiguous}"
    )
    print(
        f"unmatched         = {unmatched}"
    )
    print()

    # ----------------------------------------------------------------------
    # 4. Attempt to identify an index law.
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("INDEX-LAW CANDIDATE TABLE")
    print("=" * 78)

    for key in REFERENCE:
        matches = mapping[key]

        if len(matches) == 1:
            k, ell, d = key
            j = matches[0]

            print(
                f"(k,ell,d)=({k},{ell},{d})"
                f"  -> j={j}"
            )

    print()

    # ----------------------------------------------------------------------
    # Simple candidate forms.
    # ----------------------------------------------------------------------

    candidates = {
        "d": lambda k, ell, d: d,
        "d+1": lambda k, ell, d: d + 1,
        "ell-d": lambda k, ell, d: ell - d,
        "ell-d-1": lambda k, ell, d: ell - d - 1,
        "k+d": lambda k, ell, d: k + d,
        "k+d-1": lambda k, ell, d: k + d - 1,
        "ell-k-d": lambda k, ell, d: ell - k - d,
        "ell-k-d+1": lambda k, ell, d: ell - k - d + 1,
        "k+ell-d": lambda k, ell, d: k + ell - d,
        "k+ell-d-1": lambda k, ell, d: k + ell - d - 1,
    }

    print("=" * 78)
    print("SIMPLE INDEX-LAW CHECK")
    print("=" * 78)

    for name, fn in candidates.items():

        failures = 0
        used = 0

        for key, matches in mapping.items():

            if len(matches) != 1:
                continue

            k, ell, d = key
            predicted = fn(k, ell, d)

            used += 1

            if predicted != matches[0]:
                failures += 1

        print(
            f"{name:18s} "
            f"used={used:2d} "
            f"failures={failures:2d}"
        )

    print()

    # ----------------------------------------------------------------------
    # 5. Exact degree/leading diagnostic.
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("REFERENCE VS FOUND COEFFICIENT DIAGNOSTIC")
    print("=" * 78)

    for key, matches in mapping.items():

        if len(matches) != 1:
            continue

        k, ell, d = key
        j = matches[0]

        actual = cache[(k, ell)]["coeffs"][j]

        print(
            f"({k},{ell}) d={d} j={j}"
        )
        print(
            f"  degree = {degree_N(actual)}"
        )
        print(
            f"  lead   = {lead_N(actual)}"
        )
        print(
            f"  exact  = {sp.factor(actual)}"
        )
        print()

    # ----------------------------------------------------------------------
    # Final diagnostic.
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    if unmatched == 0:
        print(
            "Every Experiment-126 reference coefficient was found"
        )
        print(
            "exactly somewhere in the independently reconstructed"
        )
        print(
            "small-ell Newton tensor."
        )
        print()
        print(
            "The next step is now to derive the resulting index map."
        )
    else:
        print(
            "Some Experiment-126 reference coefficients were NOT"
        )
        print(
            "found in the reconstructed Newton tensor."
        )
        print()
        print(
            "That would indicate an indexing convention mismatch,"
        )
        print(
            "a coefficient normalization mismatch, or a difference"
        )
        print(
            "between the Experiment-126 extraction and the standard"
        )
        print(
            "Newton tensor."
        )

    print()
    print(
        "IMPORTANT:"
    )
    print(
        "This experiment does not assume j=d+1."
    )
    print(
        "It identifies j from exact polynomial equality."
    )
    print(
        "Therefore it isolates the indexing problem before any"
    )
    print(
        "new symbolic-law search is attempted."
    )
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print()
        print(
            "FATAL:",
            type(exc).__name__,
            str(exc),
        )
        sys.exit(1)

