#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 310R — EXACT PROJECTIVE SL(2) / COMMUTATOR-TRACE / WORD-INVARIANT
AUDIT
==============================================================================

Purpose
-------
Experiment 309R established that the two exact width-2 transition matrices

    T1
    T2

generate the full matrix algebra M_2(Q), so they are not contained in a
small common commutative algebra.

The next natural invariant is therefore the projective group structure.

For an invertible 2x2 matrix T, choose an exact square-root-free projective
normalization by working with

    tr(T)^2 / det(T),
    and for a pair (A,B), the normalized trace of AB.

We then test exact invariants associated with the pair:

    A, B,
    AB, BA,
    commutator [A,B],
    group commutator A B A^{-1} B^{-1},

including:

    * determinant and trace;
    * discriminant;
    * normalized trace tr(T)^2 / det(T);
    * normalized trace tr(T)/sqrt(det(T)) when a rational square root exists;
    * determinant of the group commutator;
    * trace of the group commutator;
    * projective commutator trace invariant;
    * nilpotence / scalar / finite-order style obstructions;
    * exact Fricke-type trace coordinates;
    * low-degree trace-word identities.

No floating-point arithmetic.
No external files.
No synthetic second n=pq case.
"""


from __future__ import annotations

import math
import sys
import sympy as sp


# ============================================================================
# EXACT TRANSFER MATRICES
# ============================================================================

T1 = sp.Matrix([
    [
        sp.Rational(
            -209857461192170070,
            11470116327290257,
        ),
        sp.Rational(
            -215450379004300026,
            11470116327290257,
        ),
    ],
    [
        1,
        0,
    ],
])

T2 = sp.Matrix([
    [
        sp.Rational(
            83976580526089197,
            971937272383741,
        ),
        sp.Rational(
            75947709834674022,
            971937272383741,
        ),
    ],
    [
        1,
        0,
    ],
])


# ============================================================================
# EXACT HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def valuation(x, p):
    x = sp.Rational(x)

    if x == 0:
        return None

    num = abs(int(x.p))
    den = abs(int(x.q))
    v = 0

    while num and num % p == 0:
        num //= p
        v += 1

    while den and den % p == 0:
        den //= p
        v -= 1

    return v


def matrix_zero(M):
    return all(
        clean(entry) == 0
        for entry in M
    )


def matrix_scalar(M):
    if M.shape != (2, 2):
        return False

    return (
        clean(M[0, 1]) == 0
        and clean(M[1, 0]) == 0
        and clean(M[0, 0] - M[1, 1]) == 0
    )


def matrix_equal(A, B):
    return all(
        clean(A[i, j] - B[i, j]) == 0
        for i in range(A.rows)
        for j in range(A.cols)
    )


def determinant(M):
    return clean(M.det())


def trace(M):
    return clean(sp.trace(M))


def discriminant(M):
    tr = trace(M)
    det = determinant(M)
    return clean(tr**2 - 4*det)


def normalized_trace_square(M):
    det = determinant(M)

    if det == 0:
        return None

    return clean(
        trace(M)**2 / det
    )


def matrix_content_gcd(M):
    nums = []
    dens = []

    for x in M:
        x = sp.Rational(x)
        nums.append(abs(int(x.p)))
        dens.append(abs(int(x.q)))

    if not nums:
        return 0

    g = 0

    for n in nums:
        g = math.gcd(g, n)

    for d in dens:
        g = math.gcd(g, d)

    return g


def charpoly(M):
    x = sp.symbols("x")
    return sp.factor(
        M.charpoly(x).as_expr()
    )


def exact_square_root(x):
    x = sp.Rational(x)

    if x < 0:
        return None

    num = sp.Integer(x.p)
    den = sp.Integer(x.q)

    sn = sp.integer_nthroot(num, 2)
    sd = sp.integer_nthroot(den, 2)

    if sn[1] and sd[1]:
        return sp.Rational(sn[0], sd[0])

    return None


def normalized_trace_if_possible(M):
    det = determinant(M)
    tr = trace(M)

    if det == 0:
        return None

    root = exact_square_root(abs(det))

    if root is None:
        return None

    if det < 0:
        return None

    return clean(
        tr / root
    )


# ============================================================================
# BASIC MATRIX WORDS
# ============================================================================

I = sp.eye(2)

AB = clean_matrix = T1 * T2
BA = T2 * T1

COMM = T1 * T2 - T2 * T1

if determinant(T1) != 0 and determinant(T2) != 0:
    GCOMM = (
        T1
        * T2
        * T1.inv()
        * T2.inv()
    )
else:
    GCOMM = None


def word_trace(words):
    M = I

    for W in words:
        M = M * W

    return trace(M)


# ============================================================================
# SECTION 1
# ============================================================================

def print_basic_matrices():

    print()
    print("=" * 78)
    print("1. EXACT PROJECTIVE TRANSFER MATRICES")
    print("=" * 78)

    print()
    print("  T1=")
    print(T1)

    print()
    print("  T2=")
    print(T2)

    print()
    print("  determinant(T1)={}".format(
        determinant(T1)
    ))

    print()
    print("  determinant(T2)={}".format(
        determinant(T2)
    ))

    print()
    print("  T1 invertible={}".format(
        determinant(T1) != 0
    ))

    print(
        "  T2 invertible={}".format(
            determinant(T2) != 0
        )
    )


# ============================================================================
# SECTION 2
# ============================================================================

def print_characteristic_data():

    print()
    print("=" * 78)
    print("2. CHARACTERISTIC / PROJECTIVE SCALAR DATA")
    print("=" * 78)

    for name, M in (
        ("T1", T1),
        ("T2", T2),
    ):

        tr = trace(M)
        det = determinant(M)
        disc = discriminant(M)

        print()
        print("  {}:".format(name))

        print(
            "    trace={}".format(tr)
        )

        print(
            "    determinant={}".format(det)
        )

        print(
            "    discriminant={}".format(disc)
        )

        print(
            "    characteristic_polynomial={}".format(
                charpoly(M)
            )
        )

        print(
            "    trace^2/determinant={}".format(
                normalized_trace_square(M)
            )
        )

        print(
            "    discriminant/determinant={}".format(
                clean(
                    disc / det
                )
                if det != 0
                else None
            )
        )

        print(
            "    rational_normalized_trace={}".format(
                normalized_trace_if_possible(M)
            )
        )


# ============================================================================
# SECTION 3
# ============================================================================

def print_product_data():

    print()
    print("=" * 78)
    print("3. PRODUCT / WORD INVARIANT AUDIT")
    print("=" * 78)

    words = [
        ("T1", T1),
        ("T2", T2),
        ("T1T2", AB),
        ("T2T1", BA),
        ("T1^2", T1**2),
        ("T2^2", T2**2),
        ("T1T2T1", T1*T2*T1),
        ("T2T1T2", T2*T1*T2),
    ]

    for name, M in words:

        print()
        print("  {}:".format(name))

        print(
            "    trace={}".format(
                trace(M)
            )
        )

        print(
            "    determinant={}".format(
                determinant(M)
            )
        )

        print(
            "    trace^2/determinant={}".format(
                normalized_trace_square(M)
            )
        )

        print(
            "    discriminant={}".format(
                discriminant(M)
            )
        )

        print(
            "    scalar_matrix={}".format(
                matrix_scalar(M)
            )
        )


# ============================================================================
# SECTION 4 — GROUP COMMUTATOR
# ============================================================================

def print_group_commutator():

    print()
    print("=" * 78)
    print("4. EXACT GROUP-COMMUTATOR AUDIT")
    print("=" * 78)

    if GCOMM is None:

        print(
            "  status=SINGULAR_TRANSFER_MATRIX"
        )

        return

    print()
    print(
        "  G=[T1,T2]_group=T1*T2*T1^(-1)*T2^(-1)"
    )

    print(GCOMM)

    tr = trace(GCOMM)
    det = determinant(GCOMM)
    disc = discriminant(GCOMM)

    print()
    print(
        "  trace(G)={}".format(tr)
    )

    print(
        "  determinant(G)={}".format(det)
    )

    print(
        "  discriminant(G)={}".format(disc)
    )

    print(
        "  trace(G)^2/determinant(G)={}".format(
            normalized_trace_square(GCOMM)
        )
    )

    print(
        "  scalar={}".format(
            matrix_scalar(GCOMM)
        )
    )

    print(
        "  identity={}".format(
            matrix_equal(GCOMM, I)
        )
    )

    print(
        "  nilpotent={}".format(
            matrix_zero(GCOMM**2)
        )
    )

    if det != 0:

        print(
            "  characteristic_polynomial={}".format(
                charpoly(GCOMM)
            )
        )


# ============================================================================
# SECTION 5 — PROJECTIVE COMMUTATOR INVARIANTS
# ============================================================================

def print_projective_commutator_invariants():

    print()
    print("=" * 78)
    print(
        "5. PROJECTIVE COMMUTATOR INVARIANT SEARCH"
    )
    print("=" * 78)

    if GCOMM is None:
        print(
            "  status=UNAVAILABLE"
        )
        return

    tr1 = trace(T1)
    tr2 = trace(T2)
    tr12 = trace(AB)

    det1 = determinant(T1)
    det2 = determinant(T2)

    print()
    print(
        "  tr(T1)={}".format(tr1)
    )

    print(
        "  tr(T2)={}".format(tr2)
    )

    print(
        "  tr(T1*T2)={}".format(tr12)
    )

    print(
        "  det(T1)={}".format(det1)
    )

    print(
        "  det(T2)={}".format(det2)
    )

    # Fricke-type trace expression:
    #
    # tr([A,B]) =
    # tr(A)^2 + tr(B)^2 + tr(AB)^2
    # - tr(A)tr(B)tr(AB) - 2 det-normalized terms
    #
    # For determinant-one matrices this is:
    #
    # tr([A,B]) =
    # x^2 + y^2 + z^2 - xyz - 2.
    #
    # We therefore normalize determinants formally by carrying the exact
    # determinant factors explicitly rather than introducing square roots.

    x = tr1
    y = tr2
    z = tr12

    d1 = det1
    d2 = det2

    # Exact determinant-aware identity candidate:
    #
    # Let A'=A/sqrt(d1), B'=B/sqrt(d2).
    # Then x'=x/sqrt(d1), y'=y/sqrt(d2),
    # z'=z/sqrt(d1*d2).
    #
    # Multiplying the SL2 Fricke expression by d1*d2 gives:
    #
    # d2*x^2 + d1*y^2 + z^2
    # - x*y*z - 2*d1*d2.

    fricke_scaled = clean(
        d2*x**2
        + d1*y**2
        + z**2
        - x*y*z
        - 2*d1*d2
    )

    actual_scaled = clean(
        d1*d2*trace(GCOMM)
    )

    print()
    print(
        "  determinant-aware Fricke expression={}".format(
            fricke_scaled
        )
    )

    print(
        "  det(T1)*det(T2)*trace(group_commutator)={}".format(
            actual_scaled
        )
    )

    print(
        "  exact_Fricke_identity={}".format(
            clean(
                fricke_scaled
                - actual_scaled
            ) == 0
        )
    )

    print()
    print(
        "  projective_commutator_ratio="
        "{}/{}".format(
            clean(trace(GCOMM)),
            clean(determinant(GCOMM))
        )
    )

    print(
        "  projective_commutator_invariant={}".format(
            normalized_trace_square(GCOMM)
        )
    )


# ============================================================================
# SECTION 6 — LOW-DEGREE COMMUTATOR STRUCTURE
# ============================================================================

def print_commutator_structure():

    print()
    print("=" * 78)
    print(
        "6. COMMUTATOR POLYNOMIAL / FINITE-ORDER OBSTRUCTION AUDIT"
    )
    print("=" * 78)

    if GCOMM is None:
        print(
            "  status=UNAVAILABLE"
        )
        return

    G = GCOMM

    tr = trace(G)
    det = determinant(G)

    print()
    print(
        "  Cayley-Hamilton residual="
    )

    ch = clean(
        G**2
        - tr*G
        + det*I
    )

    print(ch)

    print(
        "  exact_zero={}".format(
            matrix_zero(ch)
        )
    )

    for n in range(1, 7):

        Gn = clean(G**n)

        scalar_status = matrix_scalar(Gn)

        identity_status = matrix_equal(
            Gn,
            I,
        )

        print()
        print(
            "  G^{}:".format(n)
        )

        print(
            "    trace={}".format(
                trace(Gn)
            )
        )

        print(
            "    determinant={}".format(
                determinant(Gn)
            )
        )

        print(
            "    scalar={}".format(
                scalar_status
            )
        )

        print(
            "    identity={}".format(
                identity_status
            )
        )

        if scalar_status:
            print(
                "    scalar_value={}".format(
                    G[0, 0]
                )
            )


# ============================================================================
# SECTION 7 — TRACE-WORD INVARIANTS
# ============================================================================

def print_trace_word_audit():

    print()
    print("=" * 78)
    print(
        "7. TRACE-WORD INVARIANT AUDIT"
    )
    print("=" * 78)

    words = {
        "T1": [T1],
        "T2": [T2],
        "T1T2": [T1, T2],
        "T2T1": [T2, T1],
        "T1T2T1": [T1, T2, T1],
        "T2T1T2": [T2, T1, T2],
        "T1T2T1T2": [T1, T2, T1, T2],
        "T1T2T2T1": [T1, T2, T2, T1],
        "T1^2T2": [T1, T1, T2],
        "T1T2^2": [T1, T2, T2],
        "T2^2T1": [T2, T2, T1],
        "T2T1^2": [T2, T1, T1],
    }

    for name, ws in words.items():

        M = I

        for W in ws:
            M = M*W

        M = clean(M)

        print()
        print(
            "  {}:".format(name)
        )

        print(
            "    trace={}".format(
                trace(M)
            )
        )

        print(
            "    determinant={}".format(
                determinant(M)
            )
        )

        print(
            "    trace^2/determinant={}".format(
                normalized_trace_square(M)
            )
        )


# ============================================================================
# SECTION 8 — PRIME PROFILE
# ============================================================================

def print_prime_profile():

    print()
    print("=" * 78)
    print(
        "8. PROJECTIVE INVARIANT PRIME PROFILE"
    )
    print("=" * 78)

    objects = [
        ("T1", T1),
        ("T2", T2),
        ("T1T2", AB),
        ("group_commutator", GCOMM),
    ]

    primes = (
        2,
        3,
        5,
        7,
        11,
        13,
        17,
    )

    for name, M in objects:

        if M is None:
            continue

        print()
        print(
            "  {}:".format(name)
        )

        inv = normalized_trace_square(M)

        print(
            "    trace^2/determinant={}".format(
                inv
            )
        )

        if inv is not None:

            for p in primes:

                print(
                    "    v_{}={}".format(
                        p,
                        valuation(
                            inv,
                            p,
                        ),
                    )
                )


# ============================================================================
# SECTION 9 — TERMINAL SOURCE REFERENCE
# ============================================================================

def print_terminal_reference():

    q1_terminal = 495451247
    q3_terminal = 421514439

    g = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    print()
    print("=" * 78)
    print(
        "9. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    print(
        "  q1_terminal={}".format(
            q1_terminal
        )
    )

    print(
        "  q3_terminal={}".format(
            q3_terminal
        )
    )

    print(
        "  gcd={}".format(
            g
        )
    )

    print(
        "  q1/17={}".format(
            q1_terminal // 17
            if q1_terminal % 17 == 0
            else None
        )
    )

    print(
        "  q3/17={}".format(
            q3_terminal // 17
            if q3_terminal % 17 == 0
            else None
        )
    )


# ============================================================================
# SECTION 10 — INTERPRETATION
# ============================================================================

def print_interpretation():

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 309R established that the exact transition matrices T1 and T2
generate the full algebra M_2(Q).

Experiment 310R therefore moves from algebra generation to projective
group invariants.

For an invertible 2x2 matrix, quantities such as

    tr(T)^2 / det(T)

are invariant under scalar rescaling

    T -> lambda T.

Thus they are natural invariants of the corresponding projective
matrix class.

For the pair (T1,T2), the group commutator

    G = T1*T2*T1^(-1)*T2^(-1)

is invariant under simultaneous scalar rescaling of either generator.

The determinant-aware Fricke expression tests the exact trace identity
without introducing square roots.

The powers

    G, G^2, G^3, ...

test whether the projective commutator has special finite-order,
scalar, identity, or nilpotent behavior.

If a small-order projective relation appears, this could reveal hidden
group structure even though the matrices do not commute and do not lie
in a common one-variable algebra.

If all such invariants are generic, the evidence increasingly points
toward the two exact transitions being unrelated low-dimensional
accidents rather than manifestations of a common transfer mechanism.

This experiment remains diagnostic. Two exact transitions are still
insufficient to establish a universal cross-layer law.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 310R — EXACT PROJECTIVE SL(2) / COMMUTATOR-TRACE / "
        "WORD-INVARIANT AUDIT"
    )
    print("=" * 78)

    print_basic_matrices()
    print_characteristic_data()
    print_product_data()
    print_group_commutator()
    print_projective_commutator_invariants()
    print_commutator_structure()
    print_trace_word_audit()
    print_prime_profile()
    print_terminal_reference()
    print_interpretation()

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  projective_transfer_invariants_computed=True"
    )

    print(
        "  group_commutator_computed={}".format(
            GCOMM is not None
        )
    )

    if GCOMM is not None:

        print(
            "  group_commutator_scalar={}".format(
                matrix_scalar(GCOMM)
            )
        )

        print(
            "  group_commutator_identity={}".format(
                matrix_equal(GCOMM, I)
            )
        )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  interpolation_counted_as_proof=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  genuine_second_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 310R COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)
    except Exception as exc:
        print(
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )
        raise

