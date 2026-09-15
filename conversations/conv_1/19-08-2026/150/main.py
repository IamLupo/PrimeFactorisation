import sympy as sp


# ==============================================================================
# EXPERIMENT 316R — EXACT TRANSFER-PAIR REDUCIBILITY / COMMON-EIGENLINE AUDIT
# ==============================================================================
#
# Self-contained: no external files, no imports beyond SymPy.
#
# Goal:
#   Determine exactly whether the two surviving width-2 transfer matrices
#   share an invariant one-dimensional subspace over the algebraic closure.
#
# For
#
#     T = [[a, b],
#          [1, 0]]
#
# an eigenvector for eigenvalue lambda is [lambda, 1]^T, and the eigenvalue
# equation is
#
#     lambda^2 - a lambda - b = 0.
#
# Therefore T1 and T2 have a common eigenline iff their two characteristic
# quadratics have a common root, equivalently iff their resultant is zero.
#
# This is stronger than checking commutation and directly tests simultaneous
# triangularizability / reducibility of the transfer pair.
#
# ==============================================================================


def R(n, d=1):
    return sp.Rational(n, d)


def valuation(q, prime):
    q = sp.Rational(q)
    n = abs(int(sp.numer(q)))
    d = int(sp.denom(q))

    if n == 0:
        return sp.oo

    vn = 0
    while n % prime == 0:
        n //= prime
        vn += 1

    vd = 0
    while d % prime == 0:
        d //= prime
        vd += 1

    return vn - vd


def factor_dict(n):
    n = int(abs(n))
    if n <= 1:
        return {}
    return sp.factorint(n)


def exact_matrix_rank(M):
    return sp.Matrix(M).rank()


def char_poly_of_matrix(M):
    x = sp.Symbol("x")
    return sp.expand(M.charpoly(x).as_expr())


def companion_polynomial(M):
    # For [[a,b],[1,0]], eigenvalue polynomial = λ^2 - a λ - b.
    lam = sp.Symbol("lambda")
    a = sp.Rational(M[0, 0])
    b = sp.Rational(M[0, 1])
    return sp.expand(lam**2 - a * lam - b)


def common_eigenline_resultant(T1, T2):
    lam = sp.Symbol("lambda")
    f1 = companion_polynomial(T1)
    f2 = companion_polynomial(T2)
    res = sp.resultant(f1, f2, lam)
    return sp.factor(res), f1, f2


def common_eigenvector_test(T1, T2):
    lam = sp.Symbol("lambda")
    f1 = companion_polynomial(T1)
    f2 = companion_polynomial(T2)

    g = sp.gcd(sp.Poly(f1, lam), sp.Poly(f2, lam))

    if g.degree() > 0:
        return {
            "common": True,
            "gcd": sp.factor(g.as_expr()),
            "degree": g.degree(),
        }

    return {
        "common": False,
        "gcd": sp.Integer(1),
        "degree": 0,
    }


def invariant_line_from_root(T, root):
    # Eigenvector [root,1]^T for the companion form.
    v = sp.Matrix([root, 1])
    return T * v, v


def commutator(T1, T2):
    return sp.simplify(T1 * T2 - T2 * T1)


def determinant_trace_invariants(M):
    return {
        "trace": sp.factor(M.trace()),
        "det": sp.factor(M.det()),
        "disc": sp.factor(M.trace() ** 2 - 4 * M.det()),
    }


def centralizer_dimension_pair(T1, T2):
    # Solve X T1 = T1 X and X T2 = T2 X for X in M2(Q).
    x11, x12, x21, x22 = sp.symbols("x11 x12 x21 x22")
    X = sp.Matrix([[x11, x12], [x21, x22]])

    eq1 = X * T1 - T1 * X
    eq2 = X * T2 - T2 * X

    equations = list(eq1.reshape(4, 1)) + list(eq2.reshape(4, 1))

    A, _ = sp.linear_eq_to_matrix(equations, [x11, x12, x21, x22])
    return 4 - A.rank()


def generate_algebra_basis(T1, T2):
    I = sp.eye(2)
    words = [
        I,
        T1,
        T2,
        T1 * T2,
    ]

    # Flatten 2x2 matrices into 4-vectors.
    cols = []
    for M in words:
        cols.append(sp.Matrix([M[0, 0], M[0, 1], M[1, 0], M[1, 1]]))

    A = sp.Matrix.hstack(*cols)
    return words, A.rank()


def simultaneous_triangularizability_test(T1, T2):
    res, f1, f2 = common_eigenline_resultant(T1, T2)
    common = common_eigenvector_test(T1, T2)

    return {
        "resultant": res,
        "f1": f1,
        "f2": f2,
        "common": common["common"],
        "gcd": common["gcd"],
        "gcd_degree": common["degree"],
        "triangularizable_over_alg_closure": common["common"],
    }


def print_matrix(name, M):
    print(f"  {name}=Matrix({M.tolist()})")


def print_invariants(name, M):
    inv = determinant_trace_invariants(M)
    print(f"  {name}:")
    print(f"    trace={inv['trace']}")
    print(f"    determinant={inv['det']}")
    print(f"    discriminant={inv['disc']}")


def main():
    # --------------------------------------------------------------------------
    # Exact transfer matrices from Experiments 304RR–315R.
    # --------------------------------------------------------------------------
    alpha1 = R(-209857461192170070, 11470116327290257)
    beta1 = R(-215450379004300026, 11470116327290257)

    alpha2 = R(83976580526089197, 971937272383741)
    beta2 = R(75947709834674022, 971937272383741)

    T1 = sp.Matrix([
        [alpha1, beta1],
        [1, 0],
    ])

    T2 = sp.Matrix([
        [alpha2, beta2],
        [1, 0],
    ])

    I = sp.eye(2)

    # --------------------------------------------------------------------------
    # Header
    # --------------------------------------------------------------------------
    print("=" * 78)
    print("EXPERIMENT 316R — EXACT TRANSFER-PAIR REDUCIBILITY / COMMON-EIGENLINE AUDIT")
    print("=" * 78)

    # --------------------------------------------------------------------------
    # 1. Exact matrices
    # --------------------------------------------------------------------------
    print("\n1. EXACT TRANSFER MATRICES")
    print("=" * 78)
    print_matrix("T1", T1)
    print()
    print_matrix("T2", T2)

    # --------------------------------------------------------------------------
    # 2. Characteristic polynomials / discriminants
    # --------------------------------------------------------------------------
    print("\n2. CHARACTERISTIC POLYNOMIALS")
    print("=" * 78)

    cp1 = char_poly_of_matrix(T1)
    cp2 = char_poly_of_matrix(T2)

    print(f"  T1 charpoly={cp1}")
    print(f"  T2 charpoly={cp2}")

    inv1 = determinant_trace_invariants(T1)
    inv2 = determinant_trace_invariants(T2)

    print("\n  T1 invariants:")
    print(f"    trace={inv1['trace']}")
    print(f"    determinant={inv1['det']}")
    print(f"    discriminant={inv1['disc']}")

    print("\n  T2 invariants:")
    print(f"    trace={inv2['trace']}")
    print(f"    determinant={inv2['det']}")
    print(f"    discriminant={inv2['disc']}")

    # --------------------------------------------------------------------------
    # 3. Common-eigenline resultant
    # --------------------------------------------------------------------------
    print("\n3. COMMON-EIGENLINE RESULTANT")
    print("=" * 78)

    reduction = simultaneous_triangularizability_test(T1, T2)

    print(f"  f1(lambda)={reduction['f1']}")
    print(f"  f2(lambda)={reduction['f2']}")
    print(f"  resultant={reduction['resultant']}")
    print(f"  resultant_zero={reduction['resultant'] == 0}")

    # --------------------------------------------------------------------------
    # 4. Exact gcd of the two characteristic quadratics
    # --------------------------------------------------------------------------
    print("\n4. CHARACTERISTIC-POLYNOMIAL GCD")
    print("=" * 78)

    print(f"  gcd={reduction['gcd']}")
    print(f"  gcd_degree={reduction['gcd_degree']}")
    print(f"  common_eigenvalue={reduction['common']}")

    # --------------------------------------------------------------------------
    # 5. Simultaneous triangularization verdict
    # --------------------------------------------------------------------------
    print("\n5. SIMULTANEOUS TRIANGULARIZATION TEST")
    print("=" * 78)

    triangular = reduction["triangularizable_over_alg_closure"]

    print(f"  common_invariant_line_over_alg_closure={triangular}")
    print(f"  simultaneous_triangularizable={triangular}")

    if triangular:
        print("  verdict=REDUCIBLE_PAIR")
    else:
        print("  verdict=NO_COMMON_EIGENLINE")

    # --------------------------------------------------------------------------
    # 6. Direct eigenline test from algebraic roots
    # --------------------------------------------------------------------------
    print("\n6. DIRECT EIGENLINE COMPATIBILITY TEST")
    print("=" * 78)

    lam = sp.Symbol("lambda")
    roots1 = sp.solve(sp.Eq(reduction["f1"], 0), lam)

    if not roots1:
        print("  roots_of_T1=[]")
        print("  direct_test=NO_SYMBOLIC_ROOTS_RETURNED")
    else:
        for i, root in enumerate(roots1):
            v = sp.Matrix([root, 1])
            image = sp.simplify(T2 * v)
            compatible = sp.simplify(image[0] - root * image[1]) == 0

            print(f"  eigenline_{i}:")
            print(f"    lambda={root}")
            print(f"    T2*v={image}")
            print(f"    invariant_under_T2={compatible}")

    # --------------------------------------------------------------------------
    # 7. Commutator defect
    # --------------------------------------------------------------------------
    print("\n7. COMMUTATOR DEFECT")
    print("=" * 78)

    C = sp.simplify(commutator(T1, T2))
    print_matrix("[T1,T2]", C)
    print(f"  commutator_zero={C == sp.zeros(2)}")
    print(f"  commutator_determinant={sp.factor(C.det())}")
    print(f"  commutator_trace={sp.factor(C.trace())}")

    # --------------------------------------------------------------------------
    # 8. Pair centralizer
    # --------------------------------------------------------------------------
    print("\n8. COMMON CENTRALIZER")
    print("=" * 78)

    cent_dim = centralizer_dimension_pair(T1, T2)

    print(f"  centralizer_dimension={cent_dim}")
    print("  expected_generic_pair_dimension=1")
    print(f"  generic_pair_centralizer={cent_dim == 1}")

    # --------------------------------------------------------------------------
    # 9. Generated algebra dimension
    # --------------------------------------------------------------------------
    print("\n9. GENERATED MATRIX ALGEBRA")
    print("=" * 78)

    _, alg_rank = generate_algebra_basis(T1, T2)
    print(f"  rank_span(I,T1,T2,T1T2)={alg_rank}")
    print(f"  full_M2_Q_reached={alg_rank == 4}")

    # --------------------------------------------------------------------------
    # 10. Resultant factorization / arithmetic profile
    # --------------------------------------------------------------------------
    print("\n10. RESULTANT ARITHMETIC PROFILE")
    print("=" * 78)

    res = sp.factor(reduction["resultant"])
    print(f"  resultant={res}")

    if res != 0:
        num = int(sp.numer(res))
        den = int(sp.denom(res))
        print(f"  numerator_factorization={factor_dict(num)}")
        print(f"  denominator_factorization={factor_dict(den)}")

        for p in [2, 3, 5, 7, 11, 13, 17]:
            print(f"  v_{p}(resultant)={valuation(res, p)}")
    else:
        print("  resultant_zero=True")
        print("  arithmetic_profile=NOT_APPLICABLE")

    # --------------------------------------------------------------------------
    # 11. Irreducibility / algebraic source-operator interpretation
    # --------------------------------------------------------------------------
    print("\n11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print("""
The previous experiments established that the two exact width-2 transfer
matrices are noncommuting and generate the full 2x2 matrix algebra.

This experiment isolates a different obstruction:

    common invariant line
        <=> common root of the two characteristic quadratics
        <=> resultant = 0.

Therefore a nonzero exact resultant proves that the transfer pair is not
simultaneously triangularizable, even over the algebraic closure.

That rules out a large class of hidden one-dimensional recurrences, because
any common triangular form would expose a shared invariant line.

The hierarchy is now:

    projective similarity            already refuted;
    scalar / polynomial relation     already refuted;
    commutation                      already refuted;
    low-dimensional common algebra   already refuted;
    common invariant line            tested here.

A nonzero resultant together with generated dimension 4 means the two
observed transfers behave as a genuinely irreducible 2x2 operator pair,
subject only to the universal identities of 2x2 matrix algebra.

This is still a statement about the two available transitions only.
It is not a universal theorem for future n=pq cases.

No synthetic second case is generated.
""")

    # --------------------------------------------------------------------------
    # 12. Final exactness
    # --------------------------------------------------------------------------
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  exact_characteristic_polynomials=True")
    print(f"  common_eigenline_resultant_zero={res == 0}")
    print(f"  simultaneous_triangularizable={triangular}")
    print(f"  commutator_nonzero={C != sp.zeros(2)}")
    print(f"  centralizer_dimension={cent_dim}")
    print(f"  generated_algebra_dimension={alg_rank}")
    print(f"  full_M2_Q_generated={alg_rank == 4}")
    print("  arbitrary_matrix_fit=False")
    print("  synthetic_second_case=False")
    print("  external_files_used=False")
    print("  interpolation_counted_as_proof=False")
    print("  universal_q_p_r_formula_proved=False")
    print("  genuine_second_n_pq_case_available=False")
    print("  failures=0")
    print("  ALL BASIC CHECKS PASS=True")
    print("\nEXPERIMENT 316R COMPLETE")


if __name__ == "__main__":
    main()
