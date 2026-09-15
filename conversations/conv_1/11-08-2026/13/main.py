import sympy as sp
from itertools import combinations


# =============================================================================
# KAPPA STRUCTURE — 3 -> 4 -> 5 EXACT TRANSITION ANALYSIS
# =============================================================================

sp.init_printing()

p, q, r, s, t = sp.symbols("p q r s t")
e1, e2, e3, e4, e5, n = sp.symbols(
    "e1 e2 e3 e4 e5 n"
)


# =============================================================================
# HELPERS
# =============================================================================

def elementary_symmetric(vars_):
    """
    Return elementary symmetric polynomials e1, e2, ..., en
    for the supplied variables.
    """
    result = []

    for k in range(1, len(vars_) + 1):
        result.append(
            sp.expand(
                sum(
                    sp.prod(combo)
                    for combo in combinations(vars_, k)
                )
            )
        )

    return result


def M1(vars_):
    """
    M1 = product(1 + x_i)
    """
    result = sp.Integer(1)

    for x in vars_:
        result *= (1 + x)

    return sp.expand(result)


def X(vars_):
    """
    X = product(1 + x_i^3)
    """
    result = sp.Integer(1)

    for x in vars_:
        result *= (1 + x**3)

    return sp.expand(result)


def exact_P(vars_):
    """
    P = cancel(X / M1)
    """
    return sp.cancel(X(vars_) / M1(vars_))


def exact_kappa(vars_):
    """
    kappa = 1 - P
    """
    return sp.expand(1 - exact_P(vars_))


def symmetrize_exact(expr, vars_):
    """
    SymPy's symmetrize() returns:

        (symmetric_polynomial, remainder, mapping)

    in current versions.

    This wrapper normalizes the interface and returns all three pieces.
    """
    result = sp.symmetrize(
        expr,
        vars_,
        formal=True
    )

    if len(result) != 3:
        raise RuntimeError(
            f"Unexpected symmetrize() return length: {len(result)}"
        )

    symmetric_part, remainder, mapping = result

    return (
        sp.expand(symmetric_part),
        sp.expand(remainder),
        mapping
    )


# =============================================================================
# ELEMENTARY-SYMMETRIC VARIABLES
# =============================================================================

E3 = elementary_symmetric([p, q, r])
E4 = elementary_symmetric([p, q, r, s])
E5 = elementary_symmetric([p, q, r, s, t])


# =============================================================================
# TEST 1: EXACT 3-BODY DEFINITIONS
# =============================================================================

print("=" * 80)
print("KAPPA STRUCTURE — 3 -> 4 -> 5 EXACT TRANSITION ANALYSIS")
print("=" * 80)

print()
print("=" * 80)
print("TEST 1: EXACT 3-BODY DEFINITIONS")
print("=" * 80)

vars3 = [p, q, r]

M1_3 = M1(vars3)
X3 = X(vars3)
P3 = exact_P(vars3)
kappa3 = exact_kappa(vars3)

print("M1_3 =")
print(M1_3)
print()

print("X3 =")
print(X3)
print()

print("P3 = cancel(X3/M1_3)")
print(P3)
print()

print("kappa3 = 1 - P3")
print(kappa3)
print()

check_norm3 = sp.simplify(
    kappa3 - (1 - P3)
)

print("Check:")
print(check_norm3)

assert check_norm3 == 0

print()
print("PASS: exact 3-body normalization confirmed.")


# =============================================================================
# TEST 2: EXACT 4-BODY DEFINITIONS
# =============================================================================

print()
print("=" * 80)
print("TEST 2: EXACT 4-BODY DEFINITIONS")
print("=" * 80)

vars4 = [p, q, r, s]

M1_4 = M1(vars4)
X4 = X(vars4)
P4 = exact_P(vars4)
kappa4 = exact_kappa(vars4)

print("M1_4 =")
print(M1_4)
print()

print("X4 =")
print(X4)
print()

print("P4 = cancel(X4/M1_4)")
print(P4)
print()

print("kappa4 = 1 - P4")
print(kappa4)
print()

check_norm4 = sp.simplify(
    kappa4 - (1 - P4)
)

assert check_norm4 == 0

print("PASS: exact 4-body normalization confirmed.")


# =============================================================================
# TEST 3: EXACT 5-BODY DEFINITIONS
# =============================================================================

print()
print("=" * 80)
print("TEST 3: EXACT 5-BODY DEFINITIONS")
print("=" * 80)

vars5 = [p, q, r, s, t]

M1_5 = M1(vars5)
X5 = X(vars5)
P5 = exact_P(vars5)
kappa5 = exact_kappa(vars5)

print("M1_5 =")
print(M1_5)
print()

print("X5 =")
print(X5)
print()

print("P5 = cancel(X5/M1_5)")
print(P5)
print()

print("kappa5 = 1 - P5")
print(kappa5)
print()

check_norm5 = sp.simplify(
    kappa5 - (1 - P5)
)

assert check_norm5 == 0

print("PASS: exact 5-body normalization confirmed.")


# =============================================================================
# TEST 4: EXACT ELEMENTARY-SYMMETRIC REDUCTION
# =============================================================================

print()
print("=" * 80)
print("TEST 4: EXACT ELEMENTARY-SYMMETRIC REDUCTION")
print("=" * 80)


kappa3_sym_raw, rem3, mapping3 = symmetrize_exact(
    kappa3,
    vars3
)

kappa4_sym_raw, rem4, mapping4 = symmetrize_exact(
    kappa4,
    vars4
)

kappa5_sym_raw, rem5, mapping5 = symmetrize_exact(
    kappa5,
    vars5
)


# -------------------------------------------------------------------------
# Replace SymPy's formal symbols with our e1, e2, ... symbols.
#
# With formal=True SymPy supplies a mapping such as:
#
#   [(s1, p + q + r), ...]
#
# We explicitly use the elementary-symmetric expressions already constructed
# above, which makes the result independent of the particular dummy names
# chosen by SymPy.
# -------------------------------------------------------------------------

formal3 = [sp.Symbol("s1"), sp.Symbol("s2"), sp.Symbol("s3")]
formal4 = [sp.Symbol("s1"), sp.Symbol("s2"), sp.Symbol("s3"), sp.Symbol("s4")]
formal5 = [
    sp.Symbol("s1"),
    sp.Symbol("s2"),
    sp.Symbol("s3"),
    sp.Symbol("s4"),
    sp.Symbol("s5"),
]


def formal_to_e(expr, count):
    replacements = {
        sp.Symbol("s1"): e1,
        sp.Symbol("s2"): e2,
        sp.Symbol("s3"): e3,
        sp.Symbol("s4"): e4,
        sp.Symbol("s5"): e5,
    }

    return sp.expand(expr.subs(replacements))


kappa3_sym = formal_to_e(kappa3_sym_raw, 3)
kappa4_sym = formal_to_e(kappa4_sym_raw, 4)
kappa5_sym = formal_to_e(kappa5_sym_raw, 5)


print()
print("kappa3 symmetric =")
print(kappa3_sym)

print()
print("3-body symmetrization remainder =")
print(rem3)

print()
print("kappa4 symmetric =")
print(kappa4_sym)

print()
print("4-body symmetrization remainder =")
print(rem4)

print()
print("kappa5 symmetric =")
print(kappa5_sym)

print()
print("5-body symmetrization remainder =")
print(rem5)

assert sp.expand(rem3) == 0
assert sp.expand(rem4) == 0
assert sp.expand(rem5) == 0

print()
print("PASS: exact elementary-symmetric reductions confirmed.")


# =============================================================================
# TEST 5: DIRECT KAPPA VS SYMMETRIC KAPPA
# =============================================================================

print()
print("=" * 80)
print("TEST 5: DIRECT KAPPA VS SYMMETRIC KAPPA")
print("=" * 80)


# Build the elementary-symmetric expressions directly.
#
# These are obtained from the exact direct polynomials and are checked by
# substituting:
#
#   e1 = sum x_i
#   e2 = sum x_i*x_j
#   ...
#
# back into the symmetric formulas.
# -------------------------------------------------------------------------

direct3_from_sym = sp.expand(
    kappa3_sym.subs({
        e1: E3[0],
        e2: E3[1],
        e3: E3[2],
    })
)

direct4_from_sym = sp.expand(
    kappa4_sym.subs({
        e1: E4[0],
        e2: E4[1],
        e3: E4[2],
        e4: E4[3],
    })
)

direct5_from_sym = sp.expand(
    kappa5_sym.subs({
        e1: E5[0],
        e2: E5[1],
        e3: E5[2],
        e4: E5[3],
        e5: E5[4],
    })
)


check_direct3 = sp.expand(direct3_from_sym - kappa3)
check_direct4 = sp.expand(direct4_from_sym - kappa4)
check_direct5 = sp.expand(direct5_from_sym - kappa5)


print()
print("3-body check:")
print(check_direct3)

print()
print("4-body check:")
print(check_direct4)

print()
print("5-body check:")
print(check_direct5)


assert check_direct3 == 0
assert check_direct4 == 0
assert check_direct5 == 0

print()
print("PASS: all direct kappa formulas exactly match")
print("their elementary-symmetric representations.")


# =============================================================================
# TEST 6: CLEAN EXACT KAPPA FORMULAS
# =============================================================================

print()
print("=" * 80)
print("TEST 6: CLEAN EXACT KAPPA FORMULAS")
print("=" * 80)

print()
print("kappa3(e1,e2,e3) =")
print(kappa3_sym)

print()
print("kappa4(e1,e2,e3,e4) =")
print(kappa4_sym)

print()
print("kappa5(e1,e2,e3,e4,e5) =")
print(kappa5_sym)


# =============================================================================
# TEST 7: EXACT 3 -> 4 TRANSITION
# =============================================================================

print()
print("=" * 80)
print("TEST 7: EXACT 3 -> 4 CORRECTION")
print("=" * 80)


# Important:
#
# The 4-body elementary symmetric variables are NOT obtained by simply
# taking the 3-body variables and setting e4 = n while keeping e1,e2,e3
# fixed.
#
# When the fourth variable is n:
#
#   E1(4) = e1 + n
#   E2(4) = e2 + n*e1
#   E3(4) = e3 + n*e2
#   E4(4) = n*e3
#
# Therefore the correct transition is obtained from:
#
#   kappa4(
#       e1+n,
#       e2+n*e1,
#       e3+n*e2,
#       n*e3
#   )
#
# minus kappa3(e1,e2,e3).
#
# This is the actual 3 -> 4 transition with one new variable n.
# -------------------------------------------------------------------------

E1_4_from_3 = e1 + n
E2_4_from_3 = e2 + n * e1
E3_4_from_3 = e3 + n * e2
E4_4_from_3 = n * e3


kappa4_extended_from_3 = sp.expand(
    kappa4_sym.subs({
        e1: E1_4_from_3,
        e2: E2_4_from_3,
        e3: E3_4_from_3,
        e4: E4_4_from_3,
    })
)

correction34 = sp.factor(
    sp.expand(
        kappa4_extended_from_3 - kappa3_sym
    )
)


print()
print("Correct 3 -> 4 elementary-symmetric update:")
print("e1 -> e1 + n")
print("e2 -> e2 + n*e1")
print("e3 -> e3 + n*e2")
print("e4 -> n*e3")

print()
print("kappa4(updated) - kappa3 =")
print(sp.expand(correction34))

print()
print("Factored:")
print(correction34)

print()
print("Collected in n:")
print(sp.collect(sp.expand(correction34), n))


# Verify against direct variable substitution.
direct4_with_n = sp.expand(
    kappa4.subs(s, n)
)

direct3 = kappa3

direct_correction34 = sp.expand(
    direct4_with_n - direct3
)

# Substitute p,q,r elementary symmetric representation into the direct
# correction and compare with the elementary-symmetric transition.
direct_correction34_sym = sp.symmetrize(
    direct_correction34,
    [p, q, r],
    formal=True
)

direct34_sym_part = direct_correction34_sym[0]
direct34_sym_remainder = direct_correction34_sym[1]

direct34_sym_part = formal_to_e(
    direct34_sym_part,
    3
)

# The direct correction contains the new variable n separately.
#
# Symmetrizing only in p,q,r therefore gives a polynomial in e1,e2,e3,n.
#
# Compare it with the transition derived above.
check34 = sp.expand(
    direct34_sym_part - correction34
)

print()
print("Direct 3 -> 4 correction:")
print(direct_correction34)

print()
print("Symmetric-form transition check:")
print(check34)

print()
print("Symmetrization remainder:")
print(direct34_sym_remainder)

assert sp.expand(direct34_sym_remainder) == 0
assert check34 == 0

print()
print("PASS: exact 3 -> 4 transition confirmed.")


# =============================================================================
# TEST 8: EXACT 4 -> 5 TRANSITION
# =============================================================================

print()
print("=" * 80)
print("TEST 8: EXACT 4 -> 5 CORRECTION")
print("=" * 80)


# Adding a new variable n to the 4-body system gives:
#
#   E1' = e1 + n
#   E2' = e2 + n*e1
#   E3' = e3 + n*e2
#   E4' = e4 + n*e3
#   E5' = n*e4
#
# The 5-body formula must be evaluated using ALL of these updates.
# -------------------------------------------------------------------------

E1_5_from_4 = e1 + n
E2_5_from_4 = e2 + n * e1
E3_5_from_4 = e3 + n * e2
E4_5_from_4 = e4 + n * e3
E5_5_from_4 = n * e4


kappa5_extended_from_4 = sp.expand(
    kappa5_sym.subs({
        e1: E1_5_from_4,
        e2: E2_5_from_4,
        e3: E3_5_from_4,
        e4: E4_5_from_4,
        e5: E5_5_from_4,
    })
)

correction45 = sp.factor(
    sp.expand(
        kappa5_extended_from_4 - kappa4_sym
    )
)


print()
print("Correct 4 -> 5 elementary-symmetric update:")
print("e1 -> e1 + n")
print("e2 -> e2 + n*e1")
print("e3 -> e3 + n*e2")
print("e4 -> e4 + n*e3")
print("e5 -> n*e4")

print()
print("kappa5(updated) - kappa4 =")
print(sp.expand(correction45))

print()
print("Factored:")
print(correction45)

print()
print("Collected in n:")
print(sp.collect(sp.expand(correction45), n))


# Direct 4 -> 5 check.
direct5_with_n = sp.expand(
    kappa5.subs(t, n)
)

direct4 = kappa4

direct_correction45 = sp.expand(
    direct5_with_n - direct4
)


# Symmetrize direct correction in p,q,r,s.
direct45_sym_result = sp.symmetrize(
    direct_correction45,
    [p, q, r, s],
    formal=True
)

direct45_sym_part = direct45_sym_result[0]
direct45_sym_remainder = direct45_sym_result[1]

direct45_sym_part = formal_to_e(
    direct45_sym_part,
    4
)


check45 = sp.expand(
    direct45_sym_part - correction45
)


print()
print("Direct 4 -> 5 correction:")
print(direct_correction45)

print()
print("Symmetric-form transition check:")
print(check45)

print()
print("Symmetrization remainder:")
print(direct45_sym_remainder)

assert sp.expand(direct45_sym_remainder) == 0
assert check45 == 0

print()
print("PASS: exact 4 -> 5 transition confirmed.")


# =============================================================================
# TEST 9: DEGREE IN THE NEW VARIABLE
# =============================================================================

print()
print("=" * 80)
print("TEST 9: DEGREE IN NEW VARIABLE")
print("=" * 80)


degree34 = sp.degree(
    sp.Poly(
        sp.expand(correction34),
        n
    ),
    n
)

degree45 = sp.degree(
    sp.Poly(
        sp.expand(correction45),
        n
    ),
    n
)


print()
print("degree of 3 -> 4 correction in n =", degree34)
print("degree of 4 -> 5 correction in n =", degree45)


# The exact transition is quadratic in the newly added variable.
assert degree34 == 2
assert degree45 == 2

print()
print("PASS: both exact transitions are quadratic in the new variable.")


# =============================================================================
# TEST 10: EXPLICIT QUADRATIC COEFFICIENTS
# =============================================================================

print()
print("=" * 80)
print("TEST 10: QUADRATIC TRANSITION STRUCTURE")
print("=" * 80)


poly34 = sp.Poly(
    sp.expand(correction34),
    n
)

poly45 = sp.Poly(
    sp.expand(correction45),
    n
)


a34 = sp.expand(poly34.coeff_monomial(n**2))
b34 = sp.expand(poly34.coeff_monomial(n))
c34 = sp.expand(poly34.coeff_monomial(1))

a45 = sp.expand(poly45.coeff_monomial(n**2))
b45 = sp.expand(poly45.coeff_monomial(n))
c45 = sp.expand(poly45.coeff_monomial(1))


print()
print("3 -> 4:")
print("a34 =", a34)
print("b34 =", b34)
print("c34 =", c34)

print()
print("4 -> 5:")
print("a45 =", a45)
print("b45 =", b45)
print("c45 =", c45)


# Because adding a new variable must reduce to the original system at n=0,
# the constant term must vanish.
assert c34 == 0
assert c45 == 0

print()
print("PASS: both transition corrections vanish identically at n = 0.")


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print()
print("=" * 80)
print("FINAL RESULT")
print("=" * 80)

print()
print("Exact kappa3:")
print(kappa3_sym)

print()
print("Exact kappa4:")
print(kappa4_sym)

print()
print("Exact kappa5:")
print(kappa5_sym)

print()
print("Exact 3 -> 4 correction:")
print(sp.factor(correction34))

print()
print("Exact 4 -> 5 correction:")
print(sp.factor(correction45))

print()
print("Degrees:")
print("3 -> 4:", degree34)
print("4 -> 5:", degree45)

print()
print("ALL TESTS PASSED.")