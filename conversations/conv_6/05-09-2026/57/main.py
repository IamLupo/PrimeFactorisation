#!/usr/bin/env python3

import math
import random
import time

import numpy as np


# ========================================================================
# START EXPERIMENT 140
# Rank-1 approximation + exact integer lifting
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 140")
print("Rank-1 approximation + exact integer lifting")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
# ------------------------------------------------------------------------

R1 = [17, 43, 59]
R2 = [19, 37, 61]

print()
print("R1 =", R1)
print("R2 =", R2)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:

        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = x * x % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ------------------------------------------------------------------------
# PRIME
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME
# ------------------------------------------------------------------------

def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        p, q = sorted((p, q))

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ------------------------------------------------------------------------
# CELL
# ------------------------------------------------------------------------

def cell(n, p, q, r, s):

    k, a = divmod(p, r)
    ell, b = divmod(q, s)

    c1 = (k * b) // s
    c2 = (ell * a) // r

    beta = (k * b) % s
    alpha = (ell * a) % r

    c3 = (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)

    E = c1 + c2 + c3

    Q = n // (r * s)

    K = Q - E

    assert K == k * ell

    return {
        "Q": Q,
        "K": K,
        "E": E,
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "alpha": alpha,
        "beta": beta,
    }


# ------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------

def build_grid(n, p, q):

    g = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            g[(i, j)] = cell(
                n,
                p,
                q,
                r,
                s,
            )

    return g


# ------------------------------------------------------------------------
# DIRECT FACTOR CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    for p in range(
        2,
        math.isqrt(n) + 1,
    ):

        if n % p == 0:

            return p, n // p

    return None, None


# ------------------------------------------------------------------------
# RANK-1 SVD
#
# Q = sigma * u * v^T + residual
#
# Because the true K = k*l^T is rank 1 and E is relatively small,
# the dominant singular vectors should approximate k and ell.
# ------------------------------------------------------------------------

def rank1_ratios(Q):

    A = np.asarray(
        Q,
        dtype=np.float64,
    )

    U, S, VT = np.linalg.svd(
        A,
        full_matrices=False,
    )

    u = U[:, 0]
    v = VT[0, :]

    # Make signs positive.
    if u[0] < 0:
        u = -u

    if v[0] < 0:
        v = -v

    rho1 = u[1] / u[0]
    rho2 = u[2] / u[0]

    # Approximate l ratios too.
    eta1 = v[1] / v[0]
    eta2 = v[2] / v[0]

    rank1 = (
        S[0]
        * np.outer(
            u,
            v,
        )
    )

    residual = A - rank1

    return {
        "u": u,
        "v": v,
        "singular_values": S,
        "rho1": float(rho1),
        "rho2": float(rho2),
        "eta1": float(eta1),
        "eta2": float(eta2),
        "rank1": rank1,
        "residual": residual,
    }


# ------------------------------------------------------------------------
# INTEGER LIFT
#
# Exact relations:
#
#   r0*k0 - r1*k1 = a1-a0
#
#   r0*k0 - r2*k2 = a2-a0
#
# Let:
#
#   rho1 ~= k1/k0
#   rho2 ~= k2/k0
#
# so:
#
#   d1*k0 ~= delta1
#
#   d2*k0 ~= delta2
#
# where:
#
#   d1 = r0-r1*rho1
#   d2 = r0-r2*rho2
#
# and:
#
#   delta1 in [-(r1-1), r1-1]
#   delta2 in [-(r2-1), r2-1].
#
# Enumerating delta values is tiny.
#
# We estimate k0 using both equations simultaneously.
# ------------------------------------------------------------------------

def lift_k_vector(
    rho1,
    rho2,
):

    r0 = R1[0]
    r1 = R1[1]
    r2 = R1[2]

    d1 = (
        r0
        - r1 * rho1
    )

    d2 = (
        r0
        - r2 * rho2
    )

    candidates = []

    # ------------------------------------------------------------
    # If both coefficients vanish numerically, the approximation
    # is insufficient to determine scale.
    # ------------------------------------------------------------

    denom = (
        d1 * d1
        + d2 * d2
    )

    if denom < 1e-30:
        return {
            "d1": d1,
            "d2": d2,
            "candidates": [],
        }

    # ------------------------------------------------------------
    # Small residue differences.
    # ------------------------------------------------------------

    delta1_min = -(r1 - 1)
    delta1_max = r1 - 1

    delta2_min = -(r2 - 1)
    delta2_max = r2 - 1

    # ------------------------------------------------------------
    # Include exact zero and small values.
    # ------------------------------------------------------------

    for delta1 in range(
        delta1_min,
        delta1_max + 1,
    ):

        for delta2 in range(
            delta2_min,
            delta2_max + 1,
        ):

            # Least-squares estimate:
            #
            # minimize
            #
            #   (d1*x-delta1)^2
            # + (d2*x-delta2)^2
            #
            k0_real = (
                d1 * delta1
                + d2 * delta2
            ) / denom

            if k0_real <= 0:
                continue

            k0 = int(
                round(k0_real)
            )

            if k0 <= 0:
                continue

            k1 = int(
                round(
                    rho1 * k0
                )
            )

            k2 = int(
                round(
                    rho2 * k0
                )
            )

            if k1 <= 0 or k2 <= 0:
                continue

            # Exact residue differences implied by the integers.
            actual_delta1 = (
                r0 * k0
                - r1 * k1
            )

            actual_delta2 = (
                r0 * k0
                - r2 * k2
            )

            # Must be valid residue differences.
            if not (
                -(r1 - 1)
                <= actual_delta1
                <= r1 - 1
            ):
                continue

            if not (
                -(r2 - 1)
                <= actual_delta2
                <= r2 - 1
            ):
                continue

            # Score the approximation.
            error = (
                abs(
                    rho1
                    - FractionLike(
                        k1,
                        k0,
                    )
                )
                +
                abs(
                    rho2
                    - FractionLike(
                        k2,
                        k0,
                    )
                )
            )

            candidates.append(
                {
                    "k": (
                        k0,
                        k1,
                        k2,
                    ),
                    "delta": (
                        actual_delta1,
                        actual_delta2,
                    ),
                    "error": error,
                    "k0_real": k0_real,
                }
            )

    # Deduplicate.
    unique = {}

    for c in candidates:

        key = c["k"]

        if (
            key not in unique
            or c["error"]
            < unique[key]["error"]
        ):

            unique[key] = c

    candidates = sorted(
        unique.values(),
        key=lambda x: x["error"],
    )

    return {
        "d1": d1,
        "d2": d2,
        "candidates": candidates,
    }


# ------------------------------------------------------------------------
# Small helper avoiding Fraction construction in the lift loop.
# ------------------------------------------------------------------------

def FractionLike(a, b):

    return a / b


# ------------------------------------------------------------------------
# EXACT RADIX RESIDUE TEST
#
# Given candidate k-vector:
#
#   p = r0*k0 + a0
#
# with unknown a0.
#
# We recover possible a0 by the exact relationships.
# ------------------------------------------------------------------------

def reconstruct_p_from_k(
    k,
):

    k0, k1, k2 = k

    r0, r1, r2 = R1

    # From:
    #
    #   r0*k0+a0 = r1*k1+a1
    #
    # a1-a0 is fixed.
    #
    d1 = (
        r0 * k0
        - r1 * k1
    )

    d2 = (
        r0 * k0
        - r2 * k2
    )

    # Need an a0 satisfying:
    #
    # a1 = a0+d1
    # a2 = a0+d2
    #
    # with:
    #
    # 0 <= a0 < r0
    # 0 <= a1 < r1
    # 0 <= a2 < r2.
    #
    for a0 in range(r0):

        a1 = a0 + d1
        a2 = a0 + d2

        if not (
            0 <= a1 < r1
        ):
            continue

        if not (
            0 <= a2 < r2
        ):
            continue

        p = (
            r0 * k0
            + a0
        )

        # Verify all representations.
        if (
            p % r0 == a0
            and p % r1 == a1
            and p % r2 == a2
        ):

            return p, (
                a0,
                a1,
                a2,
            )

    return None, None


# ------------------------------------------------------------------------
# FULL CARRY VERIFICATION
# ------------------------------------------------------------------------

def verify_factor(
    n,
    p,
    q,
):

    if p * q != n:
        return False

    for r in R1:

        for s in R2:

            C = cell(
                n,
                p,
                q,
                r,
                s,
            )

            if C["K"] != C["k"] * C["ell"]:
                return False

    return True


# ------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------

random.seed(140)

TEST_BITS = [
    30,
    36,
    42,
    48,
    54,
]


for bits in TEST_BITS:

    print()
    print("=" * 72)
    print(
        "GENERATING",
        bits,
        "-BIT SEMIPRIME"
    )
    print("=" * 72)

    n, true_p, true_q = random_semiprime(bits)

    g = build_grid(
        n,
        true_p,
        true_q,
    )

    Q = [
        [
            g[(i, j)]["Q"]
            for j in range(3)
        ]
        for i in range(3)
    ]

    true_k = [
        g[(i, 0)]["k"]
        for i in range(3)
    ]

    true_ell = [
        g[(0, j)]["ell"]
        for j in range(3)
    ]

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("TRUE k =", true_k)
    print("TRUE ell =", true_ell)

    # ------------------------------------------------------------
    # Rank-1 approximation.
    # ------------------------------------------------------------

    result = rank1_ratios(Q)

    print()
    print("SVD SINGULAR VALUES")

    for x in result["singular_values"]:
        print(
            f"{x:.12e}"
        )

    print()
    print("OBSERVED RANK-1 RATIOS")

    print(
        "rho1 =",
        f"{result['rho1']:.15f}"
    )

    print(
        "rho2 =",
        f"{result['rho2']:.15f}"
    )

    print()
    print("TRUE RATIOS")

    true_rho1 = (
        true_k[1]
        / true_k[0]
    )

    true_rho2 = (
        true_k[2]
        / true_k[0]
    )

    print(
        "true rho1 =",
        f"{true_rho1:.15f}"
    )

    print(
        "true rho2 =",
        f"{true_rho2:.15f}"
    )

    print()
    print("RATIO ERRORS")

    print(
        "rho1 error =",
        f"{abs(result['rho1'] - true_rho1):.15e}"
    )

    print(
        "rho2 error =",
        f"{abs(result['rho2'] - true_rho2):.15e}"
    )

    # ------------------------------------------------------------
    # Rank-1 residual.
    # ------------------------------------------------------------

    residual_norm = np.linalg.norm(
        result["residual"]
    )

    Q_norm = np.linalg.norm(
        np.asarray(
            Q,
            dtype=np.float64,
        )
    )

    print()
    print("RANK-1 RESIDUAL")

    print(
        "||Q-Q1|| =",
        f"{residual_norm:.12e}"
    )

    print(
        "||Q|| =",
        f"{Q_norm:.12e}"
    )

    print(
        "relative residual =",
        f"{residual_norm / Q_norm:.12e}"
    )

    # ------------------------------------------------------------
    # Integer lifting.
    # ------------------------------------------------------------

    print()
    print("INTEGER LIFT")

    t0 = time.perf_counter()

    lift = lift_k_vector(
        result["rho1"],
        result["rho2"],
    )

    lift_time = (
        time.perf_counter()
        - t0
    )

    print(
        "d1 =",
        f"{lift['d1']:.15e}"
    )

    print(
        "d2 =",
        f"{lift['d2']:.15e}"
    )

    print(
        "lift candidates =",
        len(lift["candidates"])
    )

    print(
        "runtime =",
        f"{lift_time:.9f}",
        "s"
    )

    # ------------------------------------------------------------
    # Try candidates.
    # ------------------------------------------------------------

    recovered = []

    for candidate in lift["candidates"]:

        k = candidate["k"]

        p_candidate, residues = (
            reconstruct_p_from_k(
                k
            )
        )

        if p_candidate is None:
            continue

        if n % p_candidate != 0:
            continue

        q_candidate = (
            n
            // p_candidate
        )

        if not verify_factor(
            n,
            p_candidate,
            q_candidate,
        ):
            continue

        recovered.append(
            (
                p_candidate,
                q_candidate,
                k,
                residues,
                candidate,
            )
        )

    # ------------------------------------------------------------
    # Deduplicate.
    # ------------------------------------------------------------

    unique = {}

    for item in recovered:

        key = (
            item[0],
            item[1],
        )

        unique[key] = item

    recovered = list(
        unique.values()
    )

    print()
    print("EXACT FACTOR VERIFICATION")

    print(
        "verified factors =",
        len(recovered)
    )

    correct = any(
        {
            item[0],
            item[1],
        }
        ==
        {
            true_p,
            true_q,
        }
        for item in recovered
    )

    print(
        "TRUE FACTORIZATION FOUND =",
        correct
    )

    if recovered:

        print()
        print("RECOVERED")

        for item in recovered[:10]:

            p2, q2, k2, residues2, metadata = item

            print(
                "p =",
                p2,
                "q =",
                q2,
            )

            print(
                "k =",
                k2,
            )

            print(
                "residues =",
                residues2,
            )

            print(
                "lift error =",
                f"{metadata['error']:.15e}",
            )

    # ------------------------------------------------------------
    # Direct comparison.
    # ------------------------------------------------------------

    direct_work = math.isqrt(n) - 1

    print()
    print("COMPARISON")

    print(
        "sqrt(n) search scale =",
        direct_work
    )

    print(
        "integer lift states =",
        len(lift["candidates"])
    )

    if len(lift["candidates"]) > 0:

        print(
            "search reduction =",
            f"{direct_work / len(lift['candidates']):.3f}x"
        )

    # ------------------------------------------------------------
    # Verify that the true vector is represented.
    # ------------------------------------------------------------

    true_in_lift = any(
        c["k"] == tuple(true_k)
        for c in lift["candidates"]
    )

    print(
        "TRUE k VECTOR IN LIFT CANDIDATES =",
        true_in_lift
    )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 140
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 140")
print("=" * 72)
