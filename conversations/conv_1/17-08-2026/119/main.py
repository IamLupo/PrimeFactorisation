# ==============================================================
# KAPPA EXPERIMENT 183
# ODD-k C BRANCH: EXPLICIT SUPPORT-TRUNCATION THEOREM
# ==============================================================

import sympy as sp
from math import comb

# --------------------------------------------------------------
# Safe combinatorial binomial
# --------------------------------------------------------------

def Cbin(n, r):
    n = int(n)
    r = int(r)

    if r < 0:
        return 0
    if n < 0:
        return 0
    if r > n:
        return 0

    return comb(n, r)


# --------------------------------------------------------------
# Exact original C summand
# --------------------------------------------------------------

def C_summand(k, ell, s, a):
    q = ell - k - s - 1
    t = s - 2*k + 1

    first = Cbin(q, t + a)
    second = Cbin(q - 1, t + a - 1)

    return (-1)**s * Cbin(k, a) * (first + second)


# --------------------------------------------------------------
# Exact C
# --------------------------------------------------------------

def C_exact(k, ell, s):
    return sp.Integer(
        sum(C_summand(k, ell, s, a) for a in range(k + 1))
    )


# --------------------------------------------------------------
# Full Vandermonde bulk
# --------------------------------------------------------------

def C_bulk(k, ell, s):
    upper1 = ell - s - 1
    lower1 = s - k + 1

    upper2 = ell - s - 2
    lower2 = s - k

    return sp.Integer(
        (-1)**(ell + s + 1)
        * (
            Cbin(upper1, lower1)
            + Cbin(upper2, lower2)
        )
    )


# --------------------------------------------------------------
# Explicit support interval
#
# First piece:
#   0 <= t+a <= q
#
# Second piece:
#   0 <= t+a-1 <= q-1
#
# Their union gives
#   A = max(0, 2k-s-1)
#   B = min(k, ell+k-2s-2)
# --------------------------------------------------------------

def support_bounds(k, ell, s):
    A = max(0, 2*k - s - 1)
    B = min(k, ell + k - 2*s - 2)
    return A, B


# --------------------------------------------------------------
# Truncated exact support sum
# --------------------------------------------------------------

def C_support(k, ell, s):
    A, B = support_bounds(k, ell, s)

    if A > B:
        return sp.Integer(0)

    return sp.Integer(
        sum(C_summand(k, ell, s, a)
            for a in range(A, B + 1))
    )


# --------------------------------------------------------------
# Omitted-tail correction
# --------------------------------------------------------------

def C_tail(k, ell, s):
    A, B = support_bounds(k, ell, s)

    return sp.Integer(
        sum(
            C_summand(k, ell, s, a)
            for a in range(k + 1)
            if not (A <= a <= B)
        )
    )


# --------------------------------------------------------------
# 1. SUPPORT CERTIFICATE
# --------------------------------------------------------------

print("=" * 78)
print("1. EXPLICIT SUPPORT CERTIFICATE")
print("=" * 78)

support_failures = 0

for k in range(1, 12, 2):
    for ell in range(1, 52):
        for s in range(0, ell + 1):

            A, B = support_bounds(k, ell, s)

            # Determine actual active a's directly.
            active = [
                a for a in range(k + 1)
                if C_summand(k, ell, s, a) != 0
            ]

            predicted = (
                list(range(A, B + 1))
                if A <= B else []
            )

            if active != predicted:
                support_failures += 1
                print(
                    "SUPPORT FAIL:",
                    f"k={k}, ell={ell}, s={s}",
                    f"predicted={predicted}",
                    f"active={active}"
                )

print("support failures =", support_failures)


# --------------------------------------------------------------
# 2. FULL-SUPPORT THEOREM
# --------------------------------------------------------------

print()
print("=" * 78)
print("2. FULL-SUPPORT THEOREM")
print("=" * 78)

bulk_failures = 0
full_support_tests = 0

for k in range(1, 12, 2):
    for ell in range(1, 52):
        for s in range(0, ell + 1):

            A, B = support_bounds(k, ell, s)

            # Full-support region:
            # every a=0,...,k is active.
            if A == 0 and B == k:

                full_support_tests += 1

                exact = C_exact(k, ell, s)
                bulk = C_bulk(k, ell, s)

                if sp.simplify(exact - bulk) != 0:
                    bulk_failures += 1
                    print(
                        "BULK FAIL:",
                        f"k={k}, ell={ell}, s={s}",
                        f"exact={exact}",
                        f"bulk={bulk}"
                    )

print("full-support tests =", full_support_tests)
print("full-support failures =", bulk_failures)


# --------------------------------------------------------------
# 3. SUPPORT-TRUNCATION IDENTITY
# --------------------------------------------------------------

print()
print("=" * 78)
print("3. SUPPORT-TRUNCATION IDENTITY")
print("=" * 78)

tail_failures = 0

for k in range(1, 12, 2):
    for ell in range(1, 52):
        for s in range(0, ell + 1):

            exact = C_exact(k, ell, s)
            bulk = C_bulk(k, ell, s)
            tail = C_tail(k, ell, s)

            reconstructed = sp.simplify(bulk - tail)

            if sp.simplify(exact - reconstructed) != 0:
                tail_failures += 1

                A, B = support_bounds(k, ell, s)

                print(
                    "TAIL FAIL:",
                    f"k={k}, ell={ell}, s={s}",
                    f"support=({A},{B})",
                    f"exact={exact}",
                    f"bulk={bulk}",
                    f"tail={tail}",
                    f"reconstructed={reconstructed}"
                )

print("support-tail failures =", tail_failures)


# --------------------------------------------------------------
# 4. EMPTY-SUPPORT REGION
# --------------------------------------------------------------

print()
print("=" * 78)
print("4. EMPTY-SUPPORT REGION")
print("=" * 78)

empty_failures = 0

for k in range(1, 12, 2):
    for ell in range(1, 52):
        for s in range(0, ell + 1):

            A, B = support_bounds(k, ell, s)

            if A > B:
                exact = C_exact(k, ell, s)

                if exact != 0:
                    empty_failures += 1
                    print(
                        "EMPTY FAIL:",
                        f"k={k}, ell={ell}, s={s}",
                        f"support=({A},{B})",
                        f"exact={exact}"
                    )

print("empty-support failures =", empty_failures)


# --------------------------------------------------------------
# 5. FIRST UPPER-BOUNDARY LAYERS
# --------------------------------------------------------------

print()
print("=" * 78)
print("5. FIRST UPPER-BOUNDARY LAYERS")
print("=" * 78)

for k in range(1, 12, 2):
    print()
    print(f"k={k}")

    for ell in range(1, 30, 2):
        rows = []

        for s in range(0, ell + 1):
            A, B = support_bounds(k, ell, s)

            # Upper-support strip: A <= B < k
            if A <= B < k:

                exact = C_exact(k, ell, s)
                bulk = C_bulk(k, ell, s)
                tail = C_tail(k, ell, s)

                rows.append(
                    (
                        ell, s, A, B,
                        exact, bulk, tail
                    )
                )

        if rows:
            print(f"ell={ell}")
            for row in rows:
                ell0, s0, A, B, exact, bulk, tail = row

                print(
                    f"  s={s0}: "
                    f"support=[{A},{B}] "
                    f"exact={exact} "
                    f"bulk={bulk} "
                    f"tail={tail}"
                )


# --------------------------------------------------------------
# 6. RANDOM HOLDOUT
# --------------------------------------------------------------

print()
print("=" * 78)
print("6. HOLDOUT")
print("=" * 78)

holdout_failures = 0

holdouts = [
    (1, 47), (1, 48), (1, 49), (1, 50), (1, 51),
    (3, 47), (3, 48), (3, 49), (3, 50), (3, 51),
    (5, 47), (5, 48), (5, 49), (5, 50), (5, 51),
    (7, 47), (7, 48), (7, 49), (7, 50), (7, 51),
    (9, 47), (9, 48), (9, 49), (9, 50), (9, 51),
    (11, 47), (11, 48), (11, 49), (11, 50), (11, 51),
]

for k, ell in holdouts:

    failures = []

    for s in range(ell + 1):

        exact = C_exact(k, ell, s)
        reconstructed = sp.simplify(
            C_bulk(k, ell, s) - C_tail(k, ell, s)
        )

        if sp.simplify(exact - reconstructed) != 0:
            failures.append(
                (s, exact, reconstructed)
            )

    if failures:
        holdout_failures += 1
        print(
            f"(k={k},ell={ell}) status=FAIL"
        )

        for s, exact, closed in failures:
            print(
                f"  s={s}: exact={exact} "
                f"closed={closed}"
            )
    else:
        print(
            f"(k={k},ell={ell}) status=PASS"
        )


# --------------------------------------------------------------
# FINAL DIAGNOSTIC
# --------------------------------------------------------------

print()
print("=" * 78)
print("FINAL DIAGNOSTIC")
print("=" * 78)

print("support failures       =", support_failures)
print("full-support failures  =", bulk_failures)
print("tail reconstruction    =", tail_failures)
print("empty-support failures =", empty_failures)
print("holdout failures       =", holdout_failures)

if (
    support_failures == 0
    and bulk_failures == 0
    and tail_failures == 0
    and empty_failures == 0
    and holdout_failures == 0
):
    print()
    print("STATUS = COMPLETE SUPPORT-TRUNCATION PASS")
else:
    print()
    print("STATUS = INVESTIGATE")

