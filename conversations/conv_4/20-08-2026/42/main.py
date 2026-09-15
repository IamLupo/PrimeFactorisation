#!/usr/bin/env python3
# =============================================================================
# EXPERIMENT 518
# =============================================================================
# HOMOGENEOUS-LAYER CHARACTERISTIC POLYNOMIAL
# -> KAPPA RECURRENCE COEFFICIENT
#
# PURPOSE
# -------
# Move strictly upstream.
#
# We already know symbolically that the KAPPA sequence F_n has
#
#     chi_F(z) = z^4 - a1 z^3 + a2 z^2 - a3 z + a4
#
# with
#
#     a1 = 2(S+1)
#     a2 = 2N + S^2 + 3S + 1
#
# and therefore
#
#     S     = a1/2 - 1
#     Delta = 3*a1^2/4 - 2*a2 - 1.
#
# The remaining question is:
#
#     Can the ORIGINAL N-ONLY homogeneous-layer construction
#     generate recurrence invariants equivalent to a1 and a2
#     without explicitly introducing S, p, or q?
#
# NO DATA ANALYSIS
# ----------------
# This experiment is exact symbolic algebra only.
#
# HOW TO USE
# ----------
# Replace HOMOGENEOUS_LAYERS below with the exact symbolic
# expressions from the original homogeneous-layer construction.
#
# Example:
#
#     L0 = ...
#     L1 = ...
#     L2 = ...
#
# All expressions must depend only on N (and fixed constants).
#
# =============================================================================

import sympy as sp

N, z, t = sp.symbols("N z t")
S, Delta = sp.symbols("S Delta")
a1, a2, a3, a4 = sp.symbols("a1 a2 a3 a4")

# -----------------------------------------------------------------------------
# USER INPUT: ORIGINAL N-ONLY HOMOGENEOUS LAYER
# -----------------------------------------------------------------------------
#
# IMPORTANT:
# These MUST be the actual symbolic homogeneous-layer observables.
#
# Do not insert p, q, S, Delta, or any factor-dependent quantity here.
#
# Replace these example placeholders with your real construction.
#
# -----------------------------------------------------------------------------

HOMOGENEOUS_LAYERS = [
    # sp.expand(...),   # L_0
    # sp.expand(...),   # L_1
    # sp.expand(...),   # L_2
    # sp.expand(...),   # L_3
    # sp.expand(...),   # L_4
    # sp.expand(...),   # L_5
]

# -----------------------------------------------------------------------------
# GENERAL EXACT UTILITIES
# -----------------------------------------------------------------------------

def simp(expr):
    return sp.factor(sp.expand(sp.cancel(expr)))


def zero_certificate(name, expr):
    e = simp(expr)
    ok = (e == 0)
    print(f"  {name}")
    print(f"    difference = {e}")
    print(f"    PASS = {ok}")
    print()
    return ok


def depends_on(expr, symbol):
    return symbol in sp.sympify(expr).free_symbols


def polynomial_coefficients(expr, var):
    poly = sp.Poly(sp.expand(expr), var)
    return poly.all_coeffs()


# -----------------------------------------------------------------------------
# KAPPA TARGET SIDE
# -----------------------------------------------------------------------------

print("=" * 78)
print("EXPERIMENT 518 START")
print("=" * 78)
print("HOMOGENEOUS-LAYER CHARACTERISTIC POLYNOMIAL")
print("-> KAPPA RECURRENCE COEFFICIENT")
print()

print("[1] TARGET KAPPA RECURRENCE")
print("-" * 78)

kappa_chi = (
    z**4
    - 2 * (S + 1) * z**3
    + (2 * N + S**2 + 3 * S + 1) * z**2
    - (S + 1) * (2 * N + S) * z
    + N * (N + S + 1)
)

target_a1 = 2 * (S + 1)
target_a2 = 2 * N + S**2 + 3 * S + 1

print("  chi_F(z) =")
print("   ", sp.factor(kappa_chi))
print()
print("  a1_target =", target_a1)
print("  a2_target =", target_a2)
print()

print("[2] TARGET DOWNSTREAM RECOVERY")
print("-" * 78)

target_S = target_a1 / 2 - 1
target_Delta = sp.Rational(3, 4) * target_a1**2 - 2 * target_a2 - 1

print("  S_target     =", sp.expand(target_S))
print("  Delta_target =", sp.expand(target_Delta))
print()

zero_certificate(
    "S target identity",
    target_S - S,
)

zero_certificate(
    "Delta target identity",
    target_Delta - (S**2 - 4*N),
)

# -----------------------------------------------------------------------------
# INPUT VALIDATION
# -----------------------------------------------------------------------------

print("[3] HOMOGENEOUS-LAYER INPUT VALIDATION")
print("-" * 78)

if not HOMOGENEOUS_LAYERS:
    print("  NO HOMOGENEOUS_LAYERS HAVE BEEN SUPPLIED.")
    print()
    print("  This script is intentionally waiting for the exact")
    print("  N-only homogeneous-layer expressions from the earlier")
    print("  construction.")
    print()
    print("  No dummy data will be generated.")
    print()
    print("=" * 78)
    print("EXPERIMENT 518 FINISHED")
    print("=" * 78)
    raise SystemExit(0)

input_failures = 0

for i, L in enumerate(HOMOGENEOUS_LAYERS):
    L = simp(L)

    bad = sorted(
        str(s)
        for s in L.free_symbols
        if s not in {N}
    )

    print(f"  L_{i} = {L}")

    if bad:
        print(f"    INVALID symbols = {bad}")
        input_failures += 1
    else:
        print("    N-only = True")

print()

if input_failures:
    print("  INPUT VALIDATION FAILED.")
    print("  Factor-dependent symbols are not permitted.")
    print()
    print("=" * 78)
    print("EXPERIMENT 518 FINISHED")
    print("=" * 78)
    raise SystemExit(1)

# -----------------------------------------------------------------------------
# RECURRENCE SEARCH
# -----------------------------------------------------------------------------
#
# Search:
#
#   L_{n+r} = c_1 L_{n+r-1} + ... + c_r L_n
#
# with coefficients c_i depending ONLY on N.
#
# The smallest order with a nonzero exact solution is selected.
#
# -----------------------------------------------------------------------------

print("[4] MINIMAL HOMOGENEOUS-LAYER RECURRENCE SEARCH")
print("-" * 78)

layers = [simp(x) for x in HOMOGENEOUS_LAYERS]
m = len(layers)

found_recurrence = None

# We can only solve an order-r recurrence if enough symbolic terms exist.
for order in range(1, m):
    if m < 2 * order + 1:
        continue

    coeffs = sp.symbols(f"c0:{order}")

    equations = []

    # Use all windows available from the supplied symbolic sequence.
    for n0 in range(m - order):
        lhs = layers[n0 + order]

        rhs = sum(
            coeffs[j] * layers[n0 + order - 1 - j]
            for j in range(order)
        )

        equations.append(sp.expand(lhs - rhs))

    # Convert polynomial identities in N into coefficient equations.
    linear_eqs = []

    for eq in equations:
        poly = sp.Poly(sp.expand(eq), N)

        for coeff in poly.all_coeffs():
            linear_eqs.append(sp.expand(coeff))

    if not linear_eqs:
        continue

    M, rhs = sp.linear_eq_to_matrix(linear_eqs, coeffs)

    try:
        sol = sp.linsolve((M, rhs), coeffs)
    except Exception:
        continue

    if sol is sp.EmptySet:
        continue

    solutions = list(sol)

    if not solutions:
        continue

    candidate = solutions[0]

    # Reject unresolved free parameters.
    free = set()
    for value in candidate:
        free |= value.free_symbols

    free -= {N}

    if free:
        print(f"  order={order}: unresolved parameters = {sorted(free, key=str)}")
        continue

    found_recurrence = (
        order,
        candidate,
    )

    print(f"  EXACT RECURRENCE ORDER = {order}")

    for i, value in enumerate(candidate, 1):
        print(f"    c{i} = {simp(value)}")

    print()
    break

if found_recurrence is None:
    print("  NO EXACT N-ONLY RECURRENCE FOUND")
    print()
    print("  This is a genuine upstream negative result.")
    print()
    print("=" * 78)
    print("EXPERIMENT 518 FINISHED")
    print("=" * 78)
    raise SystemExit(0)

order, recurrence_coefficients = found_recurrence

# -----------------------------------------------------------------------------
# HOMOGENEOUS CHARACTERISTIC POLYNOMIAL
# -----------------------------------------------------------------------------

print("[5] HOMOGENEOUS CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

# Recurrence:
#
#   L_{n+r} = c1 L_{n+r-1} + ... + cr L_n
#
# Characteristic polynomial:
#
#   chi_L(t)
#     = t^r - c1 t^(r-1) - c2 t^(r-2) - ... - cr
#
# -----------------------------------------------------------------------------

chi_L = t**order

for i, c in enumerate(recurrence_coefficients):
    chi_L -= simp(c) * t**(order - 1 - i)

chi_L = sp.Poly(sp.expand(chi_L), t).as_expr()

print("  chi_L(t) =")
print("   ", sp.factor(chi_L))
print()

# -----------------------------------------------------------------------------
# HOMOGENEOUS CHARACTERISTIC COEFFICIENTS
# -----------------------------------------------------------------------------

print("[6] HOMOGENEOUS CHARACTERISTIC COEFFICIENTS")
print("-" * 78)

coeffs_L = polynomial_coefficients(chi_L, t)

for i, c in enumerate(coeffs_L):
    print(f"  coefficient[{i}] = {simp(c)}")

print()

# -----------------------------------------------------------------------------
# N-ONLY AUDIT
# -----------------------------------------------------------------------------

print("[7] N-ONLY CHARACTERISTIC AUDIT")
print("-" * 78)

for i, c in enumerate(coeffs_L):
    bad = c.free_symbols - {N}

    print(
        f"  coefficient[{i}] N-only = {not bad}"
        + (f"  bad={sorted(map(str, bad))}" if bad else "")
    )

print()

# -----------------------------------------------------------------------------
# DIRECT TARGET MORPHISM SEARCH
# -----------------------------------------------------------------------------
#
# Search simple exact combinations of homogeneous characteristic
# coefficients that could equal:
#
#   a1
#   a2
#
# without explicitly introducing S.
#
# Since a1/a2 themselves contain S, an EXACT bridge should instead
# appear as an observable expression which, after eliminating S,
# gives an invariant such as Delta.
#
# The most important target is therefore:
#
#   Delta = 3*a1^2/4 - 2*a2 - 1.
#
# We search low-degree polynomial combinations of homogeneous
# characteristic coefficients for an expression in N only.
# -----------------------------------------------------------------------------

print("[8] N-ONLY CHARACTERISTIC INVARIANT SEARCH")
print("-" * 78)

# Exclude leading coefficient 1.
observable_coeffs = coeffs_L[1:]

hits = []

# Test each individual coefficient.
for i, c in enumerate(observable_coeffs, start=1):
    if not depends_on(c, N):
        continue

    if c.free_symbols <= {N}:
        hits.append((f"c{i}", simp(c)))

if hits:
    for name, expr in hits:
        print(f"  N-only coefficient: {name} = {expr}")
else:
    print("  no individually N-only characteristic coefficient")

print()

# -----------------------------------------------------------------------------
# LOW-DEGREE LINEAR COMBINATION SEARCH
# -----------------------------------------------------------------------------

print("[9] LINEAR CHARACTERISTIC COMBINATION SEARCH")
print("-" * 78)

combination_hits = []

if len(observable_coeffs) >= 2:
    u = sp.symbols(f"u0:{len(observable_coeffs)}")

    combo = sum(
        u[i] * observable_coeffs[i]
        for i in range(len(observable_coeffs))
    )

    targets = [
        ("N", N),
        ("N^2", N**2),
        ("N*(N+1)", N * (N + 1)),
        ("N*(N+2)", N * (N + 2)),
        ("N^2+N+1", N**2 + N + 1),
    ]

    for target_name, target_expr in targets:
        equation = sp.Poly(
            sp.expand(combo - target_expr),
            N,
        )

        eqs = equation.all_coeffs()

        if not eqs:
            continue

        M, rhs = sp.linear_eq_to_matrix(eqs, u)

        sol = sp.linsolve((M, rhs), u)

        if sol is not sp.EmptySet:
            sols = list(sol)

            if sols:
                candidate = sols[0]

                free = set()
                for value in candidate:
                    free |= value.free_symbols

                free -= {N}

                if not free:
                    expression = simp(
                        combo.subs(
                            dict(zip(u, candidate))
                        )
                    )

                    if simp(expression - target_expr) == 0:
                        combination_hits.append(
                            (target_name, candidate)
                        )

for target_name, candidate in combination_hits:
    print(f"  HIT -> {target_name}")
    print(f"    coefficients = {candidate}")

if not combination_hits:
    print("  no exact N-only linear characteristic combination found")

print()

# -----------------------------------------------------------------------------
# RECURRENCE-COEFFICIENT MORPHISM
# -----------------------------------------------------------------------------
#
# Search whether two homogeneous recurrence coefficients can encode
# the KAPPA invariants under a rational transformation.
#
# We specifically test the invariant structure:
#
#   X = 3*A^2/4 - 2*B - 1
#
# and ask whether X is N-only.
# -----------------------------------------------------------------------------

print("[10] HOMOGENEOUS RECURRENCE -> GAP SIGNATURE")
print("-" * 78)

if len(recurrence_coefficients) >= 2:

    A = simp(recurrence_coefficients[0])
    B = simp(recurrence_coefficients[1])

    gap_candidate = simp(
        sp.Rational(3, 4) * A**2
        - 2 * B
        - 1
    )

    print("  A =", A)
    print("  B =", B)
    print()
    print("  G = 3*A^2/4 - 2*B - 1")
    print("    =", gap_candidate)

    bad = gap_candidate.free_symbols - {N}

    print("    N-only =", not bad)

    if not bad:
        print("    POSITIVE: homogeneous recurrence contains an N-only")
        print("    gap-like invariant.")
    else:
        print("    NEGATIVE: invariant still contains non-N symbols.")

else:
    print("  recurrence order is too small for the test")

print()

# -----------------------------------------------------------------------------
# RECURRENCE SHIFT / TRACE SEARCH
# -----------------------------------------------------------------------------
#
# The KAPPA spectrum has bases
#
#   P, Q, P+1, Q+1.
#
# Its trace is
#
#   a1 = 2(S+1).
#
# We therefore test whether simple transformations of the homogeneous
# characteristic polynomial expose a "trace-like" invariant.
# -----------------------------------------------------------------------------

print("[11] TRACE / CENTER SEARCH")
print("-" * 78)

leading = coeffs_L[0]

if leading != 0:
    trace_candidate = -coeffs_L[1] / leading

    print("  homogeneous trace candidate =", simp(trace_candidate))
    print("  N-only =", (trace_candidate.free_symbols <= {N}))

    if order == 4:
        center_candidate = simp(trace_candidate / 4)

        print("  order-4 spectral center candidate =",
              center_candidate)

else:
    print("  invalid leading coefficient")

print()

# -----------------------------------------------------------------------------
# CONCEPTUAL CENTERING TEST
# -----------------------------------------------------------------------------
#
# For a quartic:
#
#   z^4 + b3 z^3 + b2 z^2 + b1 z + b0
#
# translating z = xi - b3/4 eliminates the cubic term.
#
# We test whether the centered polynomial has a simple even structure.
# -----------------------------------------------------------------------------

print("[12] CENTERED HOMOGENEOUS CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

if order == 4:

    b3 = coeffs_L[1]
    center = simp(-b3 / 4)

    xi = sp.symbols("xi")

    centered = simp(
        chi_L.subs(t, xi + center)
    )

    print("  center =", center)
    print("  chi_L(xi+center) =")
    print("   ", sp.factor(centered))
    print()

    centered_poly = sp.Poly(
        sp.expand(centered),
        xi,
    )

    odd_terms = []

    for power in range(1, 5, 2):
        coeff = centered_poly.coeff_monomial(xi**power)

        if simp(coeff) != 0:
            odd_terms.append((power, simp(coeff)))

    print("  odd centered coefficients =", odd_terms)
    print("  EVEN CENTERED POLYNOMIAL =", not odd_terms)

else:
    print("  centered-quartic analysis requires order 4")

print()

# -----------------------------------------------------------------------------
# TARGET INTERFACE SUMMARY
# -----------------------------------------------------------------------------

print("[13] UPSTREAM INTERFACE SUMMARY")
print("-" * 78)

print("""
The KAPPA downstream interface is already known:

    F_0,...,F_7
        |
        v
    KAPPA recurrence
        |
        +--> a1
        |      |
        |      +--> S = a1/2 - 1
        |
        +--> a2
               |
               +--> Delta = 3*a1^2/4 - 2*a2 - 1
                              |
                              v
                       z^2 - S z + N
                              |
                              v
                             p,q

This experiment asks whether the original homogeneous-layer
sequence has an analogous recurrence whose invariants can be
mapped into the KAPPA recurrence invariants.

The desired bridge is:

    homogeneous layers
          |
          v
       chi_L(t)
          |
          v
    recurrence invariants
          |
          +----> a1
          |
          +----> a2
          |
          v
       S, Delta

No factor variables are permitted in the upstream construction.
""")

# -----------------------------------------------------------------------------
# FINAL STATUS
# -----------------------------------------------------------------------------

print("[14] EXPERIMENT STATUS")
print("-" * 78)

print("  homogeneous layers supplied      = True")
print("  exact recurrence found           = True")
print("  characteristic polynomial built = True")
print("  N-only input audit               =",
      input_failures == 0)

if order == 4:
    print("  centered quartic analysis        = True")
else:
    print("  centered quartic analysis        = False")

print()
print("PRIMARY QUESTION:")
print("""
    Does the original N-only homogeneous-layer construction
    possess a recurrence algebra whose characteristic invariants
    map exactly to the KAPPA recurrence coefficients?
""")

print("=" * 78)
print("EXPERIMENT 518 FINISHED")
print("=" * 78)
