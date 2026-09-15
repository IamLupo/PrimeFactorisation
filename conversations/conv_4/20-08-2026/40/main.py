#!/usr/bin/env python3

import sympy as sp
from math import gcd


print("=" * 78)
print("EXPERIMENT 516 START")
print("=" * 78)
print("KAPPA SEQUENCE -> UNCENTERED RECURRENCE -> GAP-SQUARE")
print("=" * 78)


# ============================================================================
# Helpers
# ============================================================================

def fact(expr):
    return sp.factor(sp.expand(expr))


def simp(expr):
    return sp.simplify(sp.expand(expr))


def is_zero(expr):
    return sp.simplify(expr) == 0


def cert(name, expr):
    diff = sp.factor(sp.simplify(expr))
    ok = diff == 0
    print(f"  {name}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")
    return ok


# ============================================================================
# KAPPA sequence
# ============================================================================
#
# F_n =
#   P(Q+1)^n + Q(P+1)^n
#   - (Q+1)P^n - (P+1)Q^n
#
# The four exponential bases are:
#
#   P, Q, P+1, Q+1
#
# Therefore the sequence satisfies an order-4 recurrence.
# ============================================================================

def F_symbolic(n, P, Q):
    return sp.expand(
        P * (Q + 1) ** n
        + Q * (P + 1) ** n
        - (Q + 1) * P ** n
        - (P + 1) * Q ** n
    )


def F_numeric(n, p, q):
    return (
        p * (q + 1) ** n
        + q * (p + 1) ** n
        - (q + 1) * p ** n
        - (p + 1) * q ** n
    )


# ============================================================================
# Prime-pair test set
# ============================================================================

pairs = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
    (3, 5),
    (5, 13),
    (13, 17),
    (3, 7),
    (7, 11),
]


# ============================================================================
# [1] SYMBOLIC FOUR-BASE POLYNOMIAL
# ============================================================================

print()
print("[1] FOUR-BASE CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

P, Q, z = sp.symbols("P Q z")

chi = sp.expand(
    (z - P)
    * (z - Q)
    * (z - P - 1)
    * (z - Q - 1)
)

poly = sp.Poly(chi, z)
coeffs = poly.all_coeffs()

a1 = sp.expand(-coeffs[1])
a2 = sp.expand(coeffs[2])
a3 = sp.expand(-coeffs[3])
a4 = sp.expand(coeffs[4])

print(f"  chi(z) = {sp.factor(chi)}")
print()
print(f"  a1 = {sp.factor(a1)}")
print(f"  a2 = {sp.factor(a2)}")
print(f"  a3 = {sp.factor(a3)}")
print(f"  a4 = {sp.factor(a4)}")
print()


# ============================================================================
# [2] CORRECT GAP-SQUARE IDENTITY
# ============================================================================

print()
print("[2] CORRECT UNCENTERED GAP-SQUARE IDENTITY")
print("-" * 78)

Delta_pq = sp.expand((P - Q) ** 2)

candidate_old = sp.expand(
    sp.Rational(3, 8) * a1 ** 2
    - a2
    - sp.Rational(1, 2)
)

candidate_correct = sp.expand(
    sp.Rational(3, 4) * a1 ** 2
    - 2 * a2
    - 1
)

print("  old candidate:")
print(f"    3*a1^2/8 - a2 - 1/2 = {sp.factor(candidate_old)}")
print()

cert(
    "old candidate = Delta/2",
    candidate_old - Delta_pq / 2,
)

print("  corrected candidate:")
print(f"    3*a1^2/4 - 2*a2 - 1 = {sp.factor(candidate_correct)}")
print()

symbolic_delta_pass = cert(
    "corrected candidate = Delta",
    candidate_correct - Delta_pq,
)


# ============================================================================
# [3] SAME IDENTITY IN N,S
# ============================================================================

print()
print("[3] GAP-SQUARE IN N,S")
print("-" * 78)

N, S = sp.symbols("N S")

a1_NS = 2 * (S + 1)
a2_NS = 2 * N + S ** 2 + 3 * S + 1

delta_NS = sp.factor(
    sp.Rational(3, 4) * a1_NS ** 2
    - 2 * a2_NS
    - 1
)

print(f"  a1 = {a1_NS}")
print(f"  a2 = {a2_NS}")
print(f"  Delta candidate = {delta_NS}")
print()

ns_pass = cert(
    "Delta = S^2 - 4N",
    delta_NS - (S ** 2 - 4 * N),
)


# ============================================================================
# [4] DIRECT RECURRENCE FORM
# ============================================================================

print()
print("[4] ORDER-4 RECURRENCE")
print("-" * 78)

print("""
For

    chi(z)
      = z^4 - a1 z^3 + a2 z^2 - a3 z + a4,

the KAPPA sequence obeys

    F_(n+4)
      - a1 F_(n+3)
      + a2 F_(n+2)
      - a3 F_(n+1)
      + a4 F_n
      = 0.
""")


# ============================================================================
# [5] SYMBOLIC RECURRENCE CHECK
# ============================================================================

print()
print("[5] SYMBOLIC RECURRENCE CHECK")
print("-" * 78)

symbolic_recurrence_failures = 0

for n in range(0, 9):
    expr = sp.expand(
        F_symbolic(n + 4, P, Q)
        - a1 * F_symbolic(n + 3, P, Q)
        + a2 * F_symbolic(n + 2, P, Q)
        - a3 * F_symbolic(n + 1, P, Q)
        + a4 * F_symbolic(n, P, Q)
    )

    ok = is_zero(expr)

    print(f"  n={n}: PASS={ok}")

    if not ok:
        symbolic_recurrence_failures += 1

print()
print(f"  SYMBOLIC RECURRENCE FAILURES = {symbolic_recurrence_failures}")


# ============================================================================
# [6] NUMERICAL RECURRENCE COEFFICIENT RECOVERY
# ============================================================================
#
# This is the important step.
#
# We are given only:
#
#   F_0, F_1, ..., F_7
#
# and solve for:
#
#   a1,a2,a3,a4
#
# by linear equations.
# ============================================================================

print()
print("[6] NUMERICAL RECURRENCE COEFFICIENT RECOVERY")
print("-" * 78)


def recover_recurrence_from_terms(values):
    """
    Recover a1,a2,a3,a4 from

      F[n+4] - a1 F[n+3] + a2 F[n+2]
             - a3 F[n+1] + a4 F[n] = 0.

    Uses exact rational arithmetic.
    """

    if len(values) < 8:
        raise ValueError("Need at least 8 sequence values.")

    rows = []
    rhs = []

    # Equation:
    #
    # -a1 F[n+3] + a2 F[n+2] - a3 F[n+1] + a4 F[n] = -F[n+4]
    #
    for n in range(len(values) - 4):
        rows.append([
            -values[n + 3],
             values[n + 2],
            -values[n + 1],
             values[n],
        ])
        rhs.append(-values[n + 4])

    A = sp.Matrix(rows)
    b = sp.Matrix(rhs)

    solutions = sp.linsolve((A, b))

    if solutions == sp.EmptySet:
        return None

    sol_list = list(solutions)

    if len(sol_list) != 1:
        return None

    sol = sol_list[0]

    # Reject solutions containing free symbolic parameters.
    for item in sol:
        if getattr(item, "free_symbols", set()):
            return None

    return tuple(sp.simplify(item) for item in sol)


recovery_failures = 0
numeric_results = []

for p, q in pairs:
    values = [
        sp.Integer(F_numeric(n, p, q))
        for n in range(12)
    ]

    recovered = recover_recurrence_from_terms(values)

    expected = (
        sp.Integer(2 * (p + q + 1)),
        sp.Integer(
            p * p
            + 4 * p * q
            + 3 * p
            + q * q
            + 3 * q
            + 1
        ),
        sp.Integer(
            (p + q + 1)
            * (2 * p * q + p + q)
        ),
        sp.Integer(
            p * q * (p + 1) * (q + 1)
        ),
    )

    ok = recovered == expected

    print(
        f"  ({p},{q})"
        f" recovered={recovered}"
        f" expected={expected}"
        f" PASS={ok}"
    )

    if not ok:
        recovery_failures += 1

    numeric_results.append(
        (p, q, values, recovered, expected)
    )

print()
print(f"  RECURRENCE RECOVERY FAILURES = {recovery_failures}")


# ============================================================================
# [7] GAP-SQUARE FROM OBSERVED F_n ONLY
# ============================================================================

print()
print("[7] GAP-SQUARE FROM OBSERVED F_n")
print("-" * 78)

delta_failures = 0

for p, q, values, recovered, expected in numeric_results:

    if recovered is None:
        print(f"  ({p},{q}) SKIP: recurrence unavailable")
        delta_failures += 1
        continue

    r1, r2, r3, r4 = recovered

    delta = sp.factor(
        sp.Rational(3, 4) * r1 ** 2
        - 2 * r2
        - 1
    )

    expected_delta = sp.Integer((p - q) ** 2)

    ok = delta == expected_delta

    print(
        f"  ({p},{q})"
        f" Delta={delta}"
        f" expected={expected_delta}"
        f" PASS={ok}"
    )

    if not ok:
        delta_failures += 1

print()
print(f"  DELTA RECOVERY FAILURES = {delta_failures}")


# ============================================================================
# [8] WINDOW STABILITY
# ============================================================================
#
# Important:
# recovering a recurrence from one oversized system is less interesting
# than showing that different short windows give the same coefficients.
#
# Four unknowns require four independent equations in the generic case.
# We test windows of four equations:
#
#   n=0..3
#   n=1..4
#   n=2..5
#   n=3..6
# ============================================================================

print()
print("[8] RECURRENCE WINDOW STABILITY")
print("-" * 78)


def recover_window(values, start):
    rows = []
    rhs = []

    for n in range(start, start + 4):
        rows.append([
            -values[n + 3],
             values[n + 2],
            -values[n + 1],
             values[n],
        ])
        rhs.append(-values[n + 4])

    A = sp.Matrix(rows)
    b = sp.Matrix(rhs)

    try:
        sol = sp.linsolve((A, b))
    except Exception:
        return None

    if sol == sp.EmptySet:
        return None

    sol_list = list(sol)

    if len(sol_list) != 1:
        return None

    candidate = sol_list[0]

    for item in candidate:
        if getattr(item, "free_symbols", set()):
            return None

    return tuple(sp.simplify(x) for x in candidate)


window_failures = 0

for p, q in pairs:

    values = [
        sp.Integer(F_numeric(n, p, q))
        for n in range(12)
    ]

    recovered_windows = []

    for start in range(4):
        recovered_windows.append(
            recover_window(values, start)
        )

    consistent = all(
        x == recovered_windows[0]
        for x in recovered_windows
    )

    print(
        f"  ({p},{q})"
        f" windows={recovered_windows}"
        f" consistent={consistent}"
    )

    if not consistent:
        window_failures += 1

print()
print(f"  WINDOW STABILITY FAILURES = {window_failures}")


# ============================================================================
# [9] CAN DELTA BE EXTRACTED WITHOUT S?
# ============================================================================

print()
print("[9] DIRECT INFORMATION SEPARATION")
print("-" * 78)

print("""
The operational chain tested here is:

    observed F_n
       |
       v
    a1,a2,a3,a4
       |
       v
    Delta = 3*a1^2/4 - 2*a2 - 1

No explicit S extraction is needed.

The only inputs to the final Delta formula are recurrence
coefficients reconstructed from the sequence itself.
""")


# ============================================================================
# [10] RECONSTRUCT THE FACTOR GAP
# ============================================================================

print()
print("[10] FACTOR-GAP RECOVERY")
print("-" * 78)

gap_failures = 0

for p, q, values, recovered, expected in numeric_results:

    if recovered is None:
        gap_failures += 1
        continue

    r1, r2, _, _ = recovered

    delta = sp.Integer(
        sp.Rational(3, 4) * r1 ** 2
        - 2 * r2
        - 1
    )

    gap = sp.sqrt(delta)
    expected_gap = abs(p - q)

    # SymPy can leave sqrt(k^2) as Abs(k), so compare squares.
    ok = sp.simplify(delta - expected_gap ** 2) == 0

    print(
        f"  ({p},{q})"
        f" Delta={delta}"
        f" |p-q|={expected_gap}"
        f" PASS={ok}"
    )

    if not ok:
        gap_failures += 1

print()
print(f"  GAP RECOVERY FAILURES = {gap_failures}")


# ============================================================================
# [11] OPTIONAL S RECOVERY AFTER DELTA
# ============================================================================
#
# Once a1 is known, S is actually available:
#
#   a1 = 2(S+1)
#   S  = a1/2 - 1
#
# This is useful because it shows the entire quadratic can be recovered.
# ============================================================================

print()
print("[11] S AFTER RECURRENCE RECOVERY")
print("-" * 78)

s_failures = 0
root_failures = 0

for p, q, values, recovered, expected in numeric_results:

    if recovered is None:
        s_failures += 1
        root_failures += 1
        continue

    r1, r2, _, _ = recovered

    s = sp.Integer(r1) / 2 - 1
    delta = sp.Integer(
        sp.Rational(3, 4) * r1 ** 2
        - 2 * r2
        - 1
    )

    expected_s = p + q

    s_ok = sp.simplify(s - expected_s) == 0

    # Reconstruct roots from:
    #
    #   z^2 - S z + N = 0
    #
    n = p * q
    disc = sp.simplify(s ** 2 - 4 * n)

    sqrt_disc = sp.sqrt(disc)

    roots = {
        sp.simplify((s + sqrt_disc) / 2),
        sp.simplify((s - sqrt_disc) / 2),
    }

    expected_roots = {sp.Integer(p), sp.Integer(q)}

    root_ok = roots == expected_roots

    print(
        f"  ({p},{q})"
        f" S={s} expected={expected_s}"
        f" S_PASS={s_ok}"
        f" roots={roots}"
        f" ROOT_PASS={root_ok}"
    )

    if not s_ok:
        s_failures += 1

    if not root_ok:
        root_failures += 1

print()
print(f"  S RECOVERY FAILURES = {s_failures}")
print(f"  ROOT RECOVERY FAILURES = {root_failures}")


# ============================================================================
# [12] SAME-N CONTROL
# ============================================================================

print()
print("[12] SAME-N CONTROL")
print("-" * 78)

same_n_pairs = [
    (2, 6),
    (3, 4),
    (2, 9),
    (3, 6),
    (4, 5),
    (2, 10),
    (5, 6),
    (3, 10),
]

same_n_failures = 0

for p, q in same_n_pairs:

    values = [
        sp.Integer(F_numeric(n, p, q))
        for n in range(12)
    ]

    recovered = recover_recurrence_from_terms(values)

    if recovered is None:
        print(f"  ({p},{q}) recurrence recovery FAILED")
        same_n_failures += 1
        continue

    r1, r2, _, _ = recovered

    delta = sp.Integer(
        sp.Rational(3, 4) * r1 ** 2
        - 2 * r2
        - 1
    )

    expected_delta = sp.Integer((p - q) ** 2)

    ok = delta == expected_delta

    print(
        f"  ({p},{q})"
        f" N={p*q}"
        f" S={p+q}"
        f" Delta={delta}"
        f" expected={expected_delta}"
        f" PASS={ok}"
    )

    if not ok:
        same_n_failures += 1

print()
print(f"  SAME-N CONTROL FAILURES = {same_n_failures}")


# ============================================================================
# [13] FINAL EXACT AUDIT
# ============================================================================

print()
print("[13] FINAL EXACT AUDIT")
print("-" * 78)

checks = [
    ("symbolic corrected Delta identity", symbolic_delta_pass),
    ("N,S Delta identity", ns_pass),
    ("symbolic recurrence", symbolic_recurrence_failures == 0),
    ("recurrence recovery", recovery_failures == 0),
    ("Delta recovery from F_n", delta_failures == 0),
    ("window stability", window_failures == 0),
    ("gap recovery", gap_failures == 0),
    ("S recovery", s_failures == 0),
    ("root recovery", root_failures == 0),
    ("same-N control", same_n_failures == 0),
]

for name, ok in checks:
    print(f"  {name:<34} = {ok}")

overall = all(ok for _, ok in checks)

print()
print(f"  OVERALL EXACT AUDIT = {overall}")


# ============================================================================
# [14] CONCLUSION
# ============================================================================

print()
print("[14] RESEARCH CONCLUSION")
print("-" * 78)

print(
r"""
The corrected identity is

    Delta = 3*a1^2/4 - 2*a2 - 1.

The important operational result is:

    F_n
      -> order-4 recurrence
      -> a1,a2
      -> Delta
      -> |p-q|.

Furthermore,

    S = a1/2 - 1,

so the same recurrence supplies both symmetric coordinates:

    S     = a1/2 - 1
    Delta = 3*a1^2/4 - 2*a2 - 1.

Then

    N = pq

and therefore

    z^2 - S z + N = 0

recovers p and q.

The crucial point is that this experiment reconstructs the
recurrence from exact F_n values rather than inserting P and Q
into the recurrence coefficients.

The next question is therefore no longer:

    "Can the four-base recurrence encode Delta?"

That is established.

The real upstream question is:

    Can the original N-only construction generate enough F_n
    values to make recurrence recovery possible without already
    knowing p and q?

If yes, the recurrence becomes a genuine bridge from the
observable kernel to the factor pair.
"""
)

print()
print("=" * 78)
print("EXPERIMENT 516 FINISHED")
print("=" * 78)
