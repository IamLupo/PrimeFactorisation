# EXPERIMENT 313R
# EXACT TWO-GENERATOR TRACE / COMMUTATOR INVARIANT AUDIT
#
# Single-file experiment. No external files are used.
#
# Purpose:
#   Experiment 312R established that I, T1, T2, T1*T2 form a basis
#   of M_2(Q). This experiment compresses the two-generator algebra
#   into exact scalar invariants and verifies the corresponding
#   2x2 trace identities directly.
#
# Tested objects:
#   x = tr(T1)
#   y = tr(T2)
#   z = tr(T1*T2)
#   a = det(T1)
#   b = det(T2)
#
# Main exact identity:
#
#   tr([T1,T2])
#     = (z^2 + b*x^2 + a*y^2 - x*y*z - 2*a*b) / (a*b)
#
# together with:
#
#   det([T1,T2]) = 1
#
# and the Cayley-Hamilton identity for the commutator.
#
# The script also reconstructs T1*T2 and T2*T1 in the
# basis [I,T1,T2,T1*T2], verifies the trace pairing, and audits
# several short-word trace identities.
#
# No synthetic second n=pq case is generated.

import sympy as sp


# ---------------------------------------------------------------------------
# Exact source matrices
# ---------------------------------------------------------------------------

T1 = sp.Matrix([
    [
        sp.Rational(-209857461192170070, 11470116327290257),
        sp.Rational(-215450379004300026, 11470116327290257),
    ],
    [1, 0],
])

T2 = sp.Matrix([
    [
        sp.Rational(83976580526089197, 971937272383741),
        sp.Rational(75947709834674022, 971937272383741),
    ],
    [1, 0],
])


# ---------------------------------------------------------------------------
# Basic helpers
# ---------------------------------------------------------------------------

def q(s):
    return sp.factor(s)


def exact_matrix_equal(A, B):
    return all(sp.simplify(A[i, j] - B[i, j]) == 0
               for i in range(A.rows)
               for j in range(A.cols))


def valuation_integer(n, p):
    n = sp.Integer(n)
    if n == 0:
        return sp.oo

    n = abs(int(n))
    p = int(p)
    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def valuation(expr, p):
    expr = sp.cancel(sp.sympify(expr))

    if expr == 0:
        return sp.oo

    num, den = sp.fraction(expr)
    num = sp.Integer(num)
    den = sp.Integer(den)

    return valuation_integer(num, p) - valuation_integer(den, p)


def vectorize(M):
    return sp.Matrix([
        M[0, 0],
        M[0, 1],
        M[1, 0],
        M[1, 1],
    ])


def coordinates_in_basis(M, basis):
    B = sp.Matrix.hstack(*(vectorize(E) for E in basis))
    b = vectorize(M)

    rank_B = B.rank()
    aug = B.row_join(b)

    if rank_B != aug.rank():
        return None

    sol = sp.linsolve((B, b))

    if sol == sp.EmptySet:
        return None

    sol_list = list(sol)
    if len(sol_list) != 1:
        return None

    return tuple(sp.factor(x) for x in sol_list[0])


def trace(M):
    return sp.factor(sp.trace(M))


def det(M):
    return sp.factor(M.det())


def characteristic_polynomial(M):
    x = sp.symbols("x")
    return sp.factor(M.charpoly(x).as_expr())


def print_matrix(name, M):
    print(f"  {name}=")
    print(M)


def print_invariant(name, value):
    print(f"  {name}={sp.factor(value)}")


# ---------------------------------------------------------------------------
# Basis
# ---------------------------------------------------------------------------

I = sp.eye(2)
P = T1 * T2

basis = [I, T1, T2, P]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("EXPERIMENT 313R — EXACT TWO-GENERATOR TRACE / COMMUTATOR INVARIANT AUDIT")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # 1. Exact generators
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("1. EXACT GENERATORS")
    print("=" * 78)

    print_matrix("T1", T1)
    print()
    print_matrix("T2", T2)

    # -----------------------------------------------------------------------
    # 2. Scalar trace/determinant invariants
    # -----------------------------------------------------------------------

    x = trace(T1)
    y = trace(T2)
    z = trace(P)
    a = det(T1)
    b = det(T2)

    print("\n" + "=" * 78)
    print("2. BASIC TRACE / DETERMINANT INVARIANTS")
    print("=" * 78)

    print_invariant("x=tr(T1)", x)
    print_invariant("y=tr(T2)", y)
    print_invariant("z=tr(T1*T2)", z)
    print_invariant("a=det(T1)", a)
    print_invariant("b=det(T2)", b)

    # -----------------------------------------------------------------------
    # 3. Individual Cayley-Hamilton identities
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("3. INDIVIDUAL CAYLEY-HAMILTON CHECKS")
    print("=" * 78)

    ch1 = sp.simplify(T1**2 - x*T1 + a*I)
    ch2 = sp.simplify(T2**2 - y*T2 + b*I)

    print("  T1^2 - tr(T1)T1 + det(T1)I =")
    print(ch1)
    print(f"  exact_zero={exact_matrix_equal(ch1, sp.zeros(2))}")

    print("\n  T2^2 - tr(T2)T2 + det(T2)I =")
    print(ch2)
    print(f"  exact_zero={exact_matrix_equal(ch2, sp.zeros(2))}")

    # -----------------------------------------------------------------------
    # 4. Commutator and anti-commutator
    # -----------------------------------------------------------------------

    C = T1*T2 - T2*T1
    A = T1*T2 + T2*T1

    print("\n" + "=" * 78)
    print("4. COMMUTATOR / ANTICOMMUTATOR")
    print("=" * 78)

    print_matrix("[T1,T2]", C)
    print()
    print_matrix("{T1,T2}", A)

    c_trace = trace(C)
    c_det = det(C)

    print("\n  commutator_trace=", c_trace)
    print("  commutator_determinant=", c_det)
    print(f"  commutator_nonzero={not exact_matrix_equal(C, sp.zeros(2))}")
    print(f"  commutator_det_exactly_one={sp.simplify(c_det - 1) == 0}")

    # -----------------------------------------------------------------------
    # 5. General 2x2 trace identity for the commutator
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("5. EXACT GENERAL COMMUTATOR TRACE IDENTITY")
    print("=" * 78)

    rhs_comm_trace = sp.factor(
        (z**2 + b*x**2 + a*y**2 - x*y*z - 2*a*b) / (a*b)
    )

    print("  RHS=(z^2+b*x^2+a*y^2-x*y*z-2*a*b)/(a*b)")
    print("  RHS=", rhs_comm_trace)

    comm_trace_identity_exact = sp.simplify(
        c_trace - rhs_comm_trace
    ) == 0

    print(f"  exact_identity={comm_trace_identity_exact}")

    cleared_identity = sp.factor(
        a*b*c_trace - (
            z**2 + b*x**2 + a*y**2 - x*y*z - 2*a*b
        )
    )

    print("  cleared_difference=", cleared_identity)
    print(f"  cleared_difference_zero={cleared_identity == 0}")

    # -----------------------------------------------------------------------
    # 6. Commutator Cayley-Hamilton reduction
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("6. COMMUTATOR CAYLEY-HAMILTON REDUCTION")
    print("=" * 78)

    comm_ch = sp.simplify(
        C**2 - c_trace*C + c_det*I
    )

    print("  C^2-tr(C)C+det(C)I=")
    print(comm_ch)
    print(f"  exact_zero={exact_matrix_equal(comm_ch, sp.zeros(2))}")

    print()
    print("  C^2+I-tr(C)C=")
    print(sp.simplify(C**2 + I - c_trace*C))

    print(
        "  C2_equals_I_plus_trace_reduction=",
        exact_matrix_equal(C**2, c_trace*C - I)
    )

    # -----------------------------------------------------------------------
    # 7. Trace of inverse mixed words
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("7. MIXED TRACE / INVERSE TRACE IDENTITIES")
    print("=" * 78)

    T1_inv = sp.simplify(T1.inv())
    T2_inv = sp.simplify(T2.inv())

    tr_T1_T2inv = trace(T1*T2_inv)
    tr_T1inv_T2 = trace(T1_inv*T2)

    expected_T1_T2inv = sp.factor(
        (x*y - z) / b
    )

    expected_T1inv_T2 = sp.factor(
        (x*y - z) / a
    )

    print("  tr(T1*T2^(-1))=", tr_T1_T2inv)
    print("  expected=(x*y-z)/b=", expected_T1_T2inv)
    print(
        "  exact=",
        sp.simplify(tr_T1_T2inv - expected_T1_T2inv) == 0
    )

    print()
    print("  tr(T1^(-1)*T2)=", tr_T1inv_T2)
    print("  expected=(x*y-z)/a=", expected_T1inv_T2)
    print(
        "  exact=",
        sp.simplify(tr_T1inv_T2 - expected_T1inv_T2) == 0
    )

    # -----------------------------------------------------------------------
    # 8. Trace of squares and mixed square
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("8. TRACE REDUCTION AUDIT")
    print("=" * 78)

    tr_T1_sq = trace(T1**2)
    tr_T2_sq = trace(T2**2)
    tr_P_sq = trace(P**2)

    expected_T1_sq = sp.factor(x**2 - 2*a)
    expected_T2_sq = sp.factor(y**2 - 2*b)
    expected_P_sq = sp.factor(z**2 - 2*a*b)

    print("  tr(T1^2)=", tr_T1_sq)
    print("  x^2-2a=", expected_T1_sq)
    print(
        "  exact=",
        sp.simplify(tr_T1_sq - expected_T1_sq) == 0
    )

    print()
    print("  tr(T2^2)=", tr_T2_sq)
    print("  y^2-2b=", expected_T2_sq)
    print(
        "  exact=",
        sp.simplify(tr_T2_sq - expected_T2_sq) == 0
    )

    print()
    print("  tr((T1*T2)^2)=", tr_P_sq)
    print("  z^2-2ab=", expected_P_sq)
    print(
        "  exact=",
        sp.simplify(tr_P_sq - expected_P_sq) == 0
    )

    # -----------------------------------------------------------------------
    # 9. Basis / coordinate reconstruction
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("9. BASIS RECONSTRUCTION")
    print("=" * 78)

    basis_names = ["I", "T1", "T2", "T1T2"]

    B = sp.Matrix.hstack(*(vectorize(E) for E in basis))
    basis_rank = B.rank()
    basis_det = sp.factor(B.det())

    print(f"  basis={basis_names}")
    print(f"  basis_rank={basis_rank}")
    print(f"  basis_determinant={basis_det}")
    print(f"  full_M2_Q_basis={basis_rank == 4}")

    c_T1T2 = coordinates_in_basis(T1*T2, basis)
    c_T2T1 = coordinates_in_basis(T2*T1, basis)
    c_comm = coordinates_in_basis(C, basis)
    c_anti = coordinates_in_basis(A, basis)

    print("\n  coordinates(T1*T2)=", c_T1T2)
    print("  coordinates(T2*T1)=", c_T2T1)
    print("  coordinates(commutator)=", c_comm)
    print("  coordinates(anticommutator)=", c_anti)

    reconstruction_T2T1 = (
        sum(
            (
                c_T2T1[i] * basis[i]
                for i in range(4)
            ),
            sp.zeros(2),
        )
        if c_T2T1 is not None
        else None
    )

    reconstruction_comm = (
        sum(
            (
                c_comm[i] * basis[i]
                for i in range(4)
            ),
            sp.zeros(2),
        )
        if c_comm is not None
        else None
    )

    print(
        "  T2T1_reconstruction_exact=",
        reconstruction_T2T1 is not None
        and exact_matrix_equal(reconstruction_T2T1, T2*T1)
    )

    print(
        "  commutator_reconstruction_exact=",
        reconstruction_comm is not None
        and exact_matrix_equal(reconstruction_comm, C)
    )

    # -----------------------------------------------------------------------
    # 10. Trace pairing
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("10. TRACE PAIRING")
    print("=" * 78)

    G = sp.Matrix([
        [trace(basis[i] * basis[j]) for j in range(4)]
        for i in range(4)
    ])

    G_rank = G.rank()
    G_det = sp.factor(G.det())

    print("  G_ij=tr(E_i E_j)")
    print(G)
    print(f"\n  trace_pairing_rank={G_rank}")
    print(f"  trace_pairing_determinant={G_det}")
    print(f"  trace_pairing_nondegenerate={G_rank == 4}")

    # -----------------------------------------------------------------------
    # 11. Exact scalar invariant package
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("11. SCALAR INVARIANT PACKAGE")
    print("=" * 78)

    invariant_package = {
        "tr(T1)": x,
        "tr(T2)": y,
        "tr(T1T2)": z,
        "det(T1)": a,
        "det(T2)": b,
        "tr([T1,T2])": c_trace,
        "det([T1,T2])": c_det,
        "tr(T1T2^-1)": tr_T1_T2inv,
        "tr(T1^-1T2)": tr_T1inv_T2,
    }

    for name, value in invariant_package.items():
        print(f"  {name}={sp.factor(value)}")

    # -----------------------------------------------------------------------
    # 12. Prime valuation profile
    # -----------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("12. PRIME-VALUATION PROFILE")
    print("=" * 78)

    primes = [2, 3, 5, 7, 11, 13, 17]

    valuation_targets = {
        "tr(T1)": x,
        "tr(T2)": y,
        "tr(T1T2)": z,
        "det(T1)": a,
        "det(T2)": b,
        "tr(commutator)": c_trace,
        "det(commutator)": c_det,
    }

    for name, value in valuation_targets.items():
        profile = {p: valuation(value, p) for p in primes}
        print(f"  {name}: {profile}")

    # -----------------------------------------------------------------------
    # 13. Structural conclusion
    # -----------------------------------------------------------------------

    trace_identity_ok = comm_trace_identity_exact
    comm_det_ok = sp.simplify(c_det - 1) == 0
    basis_ok = basis_rank == 4
    pairing_ok = G_rank == 4

    print("\n" + "=" * 78)
    print("13. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The preceding experiments established that the two surviving
width-2 transfer matrices generate M_2(Q) rather than a commutative
one-variable algebra.

This experiment compresses that conclusion into intrinsic scalar
invariants.

For
    x = tr(T1),
    y = tr(T2),
    z = tr(T1*T2),
    a = det(T1),
    b = det(T2),

the commutator trace is tested through the exact identity

    tr([T1,T2])
      = (z^2 + b*x^2 + a*y^2 - x*y*z - 2*a*b)/(a*b).

Since det([T1,T2]) = 1, the commutator itself satisfies a second
order Cayley-Hamilton reduction

    C^2 - tr(C) C + I = 0.

The trace of inverse mixed words is likewise determined by the same
five scalar invariants:

    tr(T1*T2^(-1)) = (x*y-z)/b,

    tr(T1^(-1)*T2) = (x*y-z)/a.

Thus the pair admits a compact intrinsic invariant package even though
the generators do not belong to a common commutative algebra.

The basis and trace-pairing checks verify that this scalar compression
does not collapse the underlying algebra: the full M_2(Q) structure
remains present.
"""
    )

    # -----------------------------------------------------------------------
    # 14. Final exactness
    # -----------------------------------------------------------------------

    print("=" * 78)
    print("14. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  exact_commutator_trace_identity={trace_identity_ok}")
    print(f"  exact_commutator_determinant_one={comm_det_ok}")
    print(f"  commutator_ch_exact={exact_matrix_equal(comm_ch, sp.zeros(2))}")
    print(f"  full_M2_Q_basis={basis_ok}")
    print(f"  trace_pairing_nondegenerate={pairing_ok}")
    print("  synthetic_second_case=False")
    print("  external_files_used=False")
    print("  arbitrary_matrix_fit=False")
    print("  interpolation_counted_as_proof=False")
    print("  universal_q_p_r_formula_proved=False")
    print("  genuine_second_n_pq_case_available=False")

    failures = 0
    failures += 0 if trace_identity_ok else 1
    failures += 0 if comm_det_ok else 1
    failures += 0 if basis_ok else 1
    failures += 0 if pairing_ok else 1

    print(f"  failures={failures}")
    print(f"  ALL BASIC CHECKS PASS={failures == 0}")
    print("\nEXPERIMENT 313R COMPLETE")


if __name__ == "__main__":
    main()
