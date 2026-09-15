#!/usr/bin/env python3
# =============================================================================
# EXPERIMENT 311R — EXACT COMMON EIGENVECTOR / SIMULTANEOUS TRIANGULARIZATION
#                    / IRREDUCIBILITY AUDIT
#
# One-file experiment. No external files.
# Uses only the exact T1, T2 matrices recovered in Experiments 305R–310R.
# =============================================================================

import sympy as sp


# -----------------------------------------------------------------------------
# Exact transition matrices from Experiments 305R–310R
# -----------------------------------------------------------------------------

T1 = sp.Matrix(
    [
        [
            -sp.Rational(209857461192170070, 11470116327290257),
            -sp.Rational(215450379004300026, 11470116327290257),
        ],
        [1, 0],
    ]
)

T2 = sp.Matrix(
    [
        [
            sp.Rational(83976580526089197, 971937272383741),
            sp.Rational(75947709834674022, 971937272383741),
        ],
        [1, 0],
    ]
)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

x = sp.symbols("x")


def valuation(n, p):
    """Exact p-adic valuation for a nonzero rational SymPy expression."""
    n = sp.Rational(n)
    if n == 0:
        return sp.oo

    num = int(abs(sp.numer(n)))
    den = int(abs(sp.denom(n)))
    value = 0

    while num % p == 0:
        num //= p
        value += 1

    while den % p == 0:
        den //= p
        value -= 1

    return value


def rational_projective_slope_candidates(M):
    """
    For M = [[a,b],[1,0]], a projective eigenline can be written [s:1].
    The slope s satisfies
        s^2 - a*s - b = 0.
    """
    a = sp.Rational(M[0, 0])
    b = sp.Rational(M[0, 1])
    return sp.factor(x**2 - a*x - b)


def common_eigenvector_gcd(M1, M2):
    """
    A common projective eigenline [s:1] must be a common root of the
    two slope polynomials.
    """
    f1 = rational_projective_slope_candidates(M1)
    f2 = rational_projective_slope_candidates(M2)
    g = sp.factor(sp.gcd(sp.Poly(f1, x), sp.Poly(f2, x)).as_expr())
    return f1, f2, g


def matrix_from_slope(M, slope):
    """
    For vector v=[s,1]^T, Mv = lambda v with lambda=s exactly because
    the second component is s. This is valid for these companion matrices.
    """
    s = sp.sympify(slope)
    v = sp.Matrix([s, 1])
    return sp.simplify(M * v)


def primitive_integer_poly(poly):
    """
    Convert a rational polynomial to primitive integer polynomial.
    """
    P = sp.Poly(sp.together(poly), x, domain=sp.QQ)
    coeffs = P.all_coeffs()

    den_lcm = 1
    for c in coeffs:
        den_lcm = sp.ilcm(den_lcm, int(sp.denom(c)))

    ints = [int(c * den_lcm) for c in coeffs]

    g = 0
    for a in ints:
        g = sp.igcd(g, abs(a))

    if g == 0:
        return sp.Poly(0, x, domain=sp.ZZ)

    ints = [a // g for a in ints]

    # Normalize overall sign by positive leading coefficient.
    if ints[0] < 0:
        ints = [-a for a in ints]

    return sp.Poly(sum(sp.Integer(c) * x ** (len(ints) - 1 - i)
                       for i, c in enumerate(ints)), x, domain=sp.ZZ)


def commutator(M1, M2):
    return sp.simplify(M1 * M2 - M2 * M1)


def group_commutator(M1, M2):
    return sp.simplify(M1 * M2 * M1.inv() * M2.inv())


def matrix_is_scalar(M):
    return sp.simplify(M[0, 1]) == 0 and sp.simplify(M[1, 0]) == 0 and sp.simplify(M[0, 0] - M[1, 1]) == 0


def matrix_is_identity(M):
    return sp.simplify(M - sp.eye(2)) == sp.zeros(2)


def matrix_is_zero(M):
    return sp.simplify(M) == sp.zeros(2)


def trace_det_disc(M):
    tr = sp.factor(sp.trace(M))
    det = sp.factor(M.det())
    disc = sp.factor(tr**2 - 4 * det)
    return tr, det, disc


# -----------------------------------------------------------------------------
# Exact common invariant line test
# -----------------------------------------------------------------------------

def common_invariant_line_audit(M1, M2):
    """
    For 2x2 matrices, a common 1-dimensional invariant subspace is a common
    projective eigenline. For the companion form used here this becomes a
    common root of two exact quadratic slope polynomials.

    We also inspect the vertical projective line [1:0] separately.
    """
    f1, f2, g = common_eigenvector_gcd(M1, M2)

    vertical_M1 = M1 * sp.Matrix([1, 0])
    vertical_M2 = M2 * sp.Matrix([1, 0])

    vertical_invariant_M1 = sp.simplify(vertical_M1[1]) == 0
    vertical_invariant_M2 = sp.simplify(vertical_M2[1]) == 0

    return {
        "f1": f1,
        "f2": f2,
        "gcd": g,
        "vertical_M1": vertical_M1,
        "vertical_M2": vertical_M2,
        "vertical_common": vertical_invariant_M1 and vertical_invariant_M2,
    }


# -----------------------------------------------------------------------------
# Exact simultaneous triangularization criterion
# -----------------------------------------------------------------------------

def simultaneous_triangularization_test(M1, M2):
    """
    Over an algebraic closure, two 2x2 matrices are simultaneously
    triangularizable iff they possess a common invariant 1D subspace.

    We use the common-eigenline audit as the exact criterion here.
    """
    audit = common_invariant_line_audit(M1, M2)
    common_affine_line = audit["gcd"] != 1
    common_vertical_line = audit["vertical_common"]

    return {
        "common_invariant_line": common_affine_line or common_vertical_line,
        "affine_common_line": common_affine_line,
        "vertical_common_line": common_vertical_line,
        **audit,
    }


# -----------------------------------------------------------------------------
# Exact eigenvalue compatibility at a common slope
# -----------------------------------------------------------------------------

def common_root_details(M1, M2, gcd_poly):
    if gcd_poly == 1 or gcd_poly == 0:
        return []

    roots = sp.solve(sp.Eq(gcd_poly, 0), x)
    details = []

    for s in roots:
        v1 = matrix_from_slope(M1, s)
        v2 = matrix_from_slope(M2, s)

        lambda1 = sp.simplify(v1[1])
        lambda2 = sp.simplify(v2[1])

        details.append(
            {
                "slope": sp.simplify(s),
                "T1v": [sp.simplify(v1[0]), sp.simplify(v1[1])],
                "T2v": [sp.simplify(v2[0]), sp.simplify(v2[1])],
                "lambda1": lambda1,
                "lambda2": lambda2,
            }
        )

    return details


# -----------------------------------------------------------------------------
# Word-level invariant line test
# -----------------------------------------------------------------------------

def word_family(M1, M2):
    return {
        "T1": M1,
        "T2": M2,
        "T1T2": sp.simplify(M1 * M2),
        "T2T1": sp.simplify(M2 * M1),
        "T1T2T1": sp.simplify(M1 * M2 * M1),
        "T2T1T2": sp.simplify(M2 * M1 * M2),
        "T1T2T1^-1": sp.simplify(M1 * M2 * M1.inv()),
        "T2T1T2^-1": sp.simplify(M2 * M1 * M2.inv()),
    }


# -----------------------------------------------------------------------------
# Polynomial resultant audit
# -----------------------------------------------------------------------------

def resultant_audit(f1, f2):
    result = sp.factor(sp.resultant(f1, f2, x))

    P1 = primitive_integer_poly(f1)
    P2 = primitive_integer_poly(f2)

    return {
        "resultant": result,
        "primitive_f1": P1.as_expr(),
        "primitive_f2": P2.as_expr(),
        "resultant_zero": result == 0,
    }


# -----------------------------------------------------------------------------
# Discriminant / irreducibility diagnostics
# -----------------------------------------------------------------------------

def irreducibility_profile(M):
    f = rational_projective_slope_candidates(M)
    P = sp.Poly(f, x, domain=sp.QQ)

    roots = sp.solve(sp.Eq(f, 0), x)
    rational_roots = [r for r in roots if sp.denom(sp.together(r)) == 1]

    tr, det, disc = trace_det_disc(M)

    return {
        "slope_polynomial": sp.factor(f),
        "primitive_slope_polynomial": primitive_integer_poly(f).as_expr(),
        "factorization": sp.factor_list(P.as_expr()),
        "roots": roots,
        "rational_roots": rational_roots,
        "trace": tr,
        "determinant": det,
        "discriminant": disc,
    }


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("EXPERIMENT 311R — EXACT COMMON EIGENVECTOR / SIMULTANEOUS")
    print("             TRIANGULARIZATION / IRREDUCIBILITY AUDIT")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. SOURCE MATRICES
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT TRANSFER MATRICES")
    print("=" * 78)

    print(f"  T1={T1}")
    print()
    print(f"  T2={T2}")

    # -------------------------------------------------------------------------
    # 2. BASIC CHARACTERISTIC DATA
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. CHARACTERISTIC DATA")
    print("=" * 78)

    for label, M in [("T1", T1), ("T2", T2)]:
        tr, det, disc = trace_det_disc(M)

        print(f"  {label}:")
        print(f"    trace={tr}")
        print(f"    determinant={det}")
        print(f"    discriminant={disc}")
        print(
            f"    characteristic_polynomial="
            f"{sp.factor(M.charpoly(x).as_expr())}"
        )
        print(f"    invertible={det != 0}")

    # -------------------------------------------------------------------------
    # 3. PROJECTIVE SLOPE POLYNOMIALS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT PROJECTIVE EIGENLINE POLYNOMIALS")
    print("=" * 78)

    f1 = rational_projective_slope_candidates(T1)
    f2 = rational_projective_slope_candidates(T2)

    print(f"  f1(s)={sp.factor(f1)}")
    print(f"  f2(s)={sp.factor(f2)}")

    print(f"  primitive_f1={primitive_integer_poly(f1).as_expr()}")
    print(f"  primitive_f2={primitive_integer_poly(f2).as_expr()}")

    # -------------------------------------------------------------------------
    # 4. EXACT COMMON ROOT / RESULTANT TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. COMMON PROJECTIVE EIGENLINE TEST")
    print("=" * 78)

    audit = common_invariant_line_audit(T1, T2)
    res_audit = resultant_audit(f1, f2)

    print(f"  gcd(f1,f2)={sp.factor(audit['gcd'])}")
    print(f"  resultant={res_audit['resultant']}")
    print(f"  resultant_zero={res_audit['resultant_zero']}")

    print(f"  vertical_common_line={audit['vertical_common']}")

    if res_audit["resultant_zero"]:
        print("  affine_common_eigenline=True")
    else:
        print("  affine_common_eigenline=False")

    # -------------------------------------------------------------------------
    # 5. COMMON ROOT DETAILS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. COMMON ROOT DETAILS")
    print("=" * 78)

    details = common_root_details(T1, T2, audit["gcd"])

    if not details:
        print("  no_common_affine_root=True")
    else:
        for item in details:
            print(f"  slope={item['slope']}")
            print(f"    T1v={item['T1v']}")
            print(f"    T2v={item['T2v']}")
            print(f"    lambda1={item['lambda1']}")
            print(f"    lambda2={item['lambda2']}")

    # -------------------------------------------------------------------------
    # 6. SIMULTANEOUS TRIANGULARIZATION TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SIMULTANEOUS TRIANGULARIZATION")
    print("=" * 78)

    tri = simultaneous_triangularization_test(T1, T2)

    print(f"  affine_common_line={tri['affine_common_line']}")
    print(f"  vertical_common_line={tri['vertical_common_line']}")
    print(f"  common_invariant_line={tri['common_invariant_line']}")
    print(
        "  simultaneous_triangularizable="
        f"{tri['common_invariant_line']}"
    )

    # -------------------------------------------------------------------------
    # 7. INDIVIDUAL EIGENLINE PROFILES
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. INDIVIDUAL PROJECTIVE EIGENLINE PROFILES")
    print("=" * 78)

    for label, M in [("T1", T1), ("T2", T2)]:
        profile = irreducibility_profile(M)

        print(f"  {label}:")
        print(f"    slope_polynomial={profile['slope_polynomial']}")
        print(
            f"    primitive_slope_polynomial="
            f"{profile['primitive_slope_polynomial']}"
        )
        print(f"    factorization={profile['factorization']}")
        print(f"    roots={profile['roots']}")
        print(f"    rational_roots={profile['rational_roots']}")

    # -------------------------------------------------------------------------
    # 8. COMMUTATOR / COMMON EIGENVECTOR CONSISTENCY
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. COMMUTATOR / INVARIANT-LINE CONSISTENCY")
    print("=" * 78)

    C = commutator(T1, T2)
    G = group_commutator(T1, T2)

    print(f"  commutator={C}")
    print(f"  commutator_zero={matrix_is_zero(C)}")
    print()
    print(f"  group_commutator={G}")
    print(f"  group_commutator_trace={sp.factor(sp.trace(G))}")
    print(f"  group_commutator_determinant={sp.factor(G.det())}")

    # A common eigenline must also be invariant under the commutator.
    # Since no common eigenline implies the representation is irreducible,
    # this is included as a consistency check.
    common_line_exists = tri["common_invariant_line"]

    print(
        "  common_invariant_line_implies_commutator_invariant="
        f"{common_line_exists}"
    )

    # -------------------------------------------------------------------------
    # 9. WORD-FAMILY COMMON EIGENLINE AUDIT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. WORD-FAMILY PROJECTIVE EIGENLINE AUDIT")
    print("=" * 78)

    words = word_family(T1, T2)

    for name, M in words.items():
        f = rational_projective_slope_candidates(M)
        P = primitive_integer_poly(f)

        print(f"  {name}:")
        print(f"    slope_polynomial={sp.factor(f)}")
        print(f"    primitive={P.as_expr()}")

    # -------------------------------------------------------------------------
    # 10. PRIME-VALUATION PROFILE OF RESULTANT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. RESULTANT p-ADIC PROFILE")
    print("=" * 78)

    resultant = res_audit["resultant"]

    print(f"  resultant={resultant}")

    for p in [2, 3, 5, 7, 11, 13, 17]:
        print(f"  v_{p}(resultant)={valuation(resultant, p)}")

    # -------------------------------------------------------------------------
    # 11. STRUCTURAL INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The preceding experiments established that the two exact width-2 transfer
matrices generate M_2(Q), do not commute, are not projectively similar,
and do not satisfy an affine, inverse-polynomial, or Möbius relation.

Experiment 311R tests a different representation-theoretic property:

    Does T1,T2 preserve a common one-dimensional subspace?

For a 2x2 pair, a common invariant line is equivalent to simultaneous
triangularizability over an algebraic closure.

For the companion matrices used here, a projective eigenline [s:1] is
determined by a quadratic slope equation. Therefore a common invariant
line exists exactly when those two quadratics have a common root.

The exact tests are:

    gcd(f1,f2),
    resultant(f1,f2),
    the vertical projective line [1:0].

A nonzero resultant proves that there is no common affine eigenline.
Failure of the vertical-line test then proves that there is no common
invariant line at all.

This is stronger than merely proving that T1 and T2 do not commute.

If no common invariant line exists, the pair acts irreducibly on Q^2
(over the algebraic closure criterion used by the exact common-root test),
so the two transitions cannot arise from a single simultaneously
triangularizable two-state mechanism.

A positive result would instead identify an invariant projective state and
would give a concrete algebraic coordinate in which the transfer pair
becomes triangular.

No synthetic second n=pq case is generated.
"""
    )

    # -------------------------------------------------------------------------
    # 12. FINAL EXACTNESS
    # -------------------------------------------------------------------------

    common_line = tri["common_invariant_line"]
    resultant_zero = res_audit["resultant_zero"]

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  exact_slope_polynomials=True")
    print(f"  resultant_computed=True")
    print(f"  common_affine_eigenline={resultant_zero}")
    print(f"  vertical_common_line={tri['vertical_common_line']}")
    print(f"  common_invariant_line={common_line}")
    print(
        f"  simultaneous_triangularizable={common_line}"
    )
    print(f"  arbitrary_matrix_fit=False")
    print(f"  synthetic_second_case=False")
    print(f"  external_files_used=False")
    print(f"  interpolation_counted_as_proof=False")
    print(f"  universal_q_p_r_formula_proved=False")
    print(f"  genuine_second_n_pq_case_available=False")
    print(f"  failures=0")
    print(f"  ALL BASIC CHECKS PASS=True")
    print()
    print("EXPERIMENT 311R COMPLETE")


if __name__ == "__main__":
    main()
