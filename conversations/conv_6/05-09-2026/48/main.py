#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 131
# Exact carry-floor Diophantine recovery
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 131")
print("Exact carry-floor Diophantine recovery")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
# ------------------------------------------------------------------------

R1 = [17, 43]
R2 = [19, 37]

r = R1[0]
s = R2[0]

R = r * s

print()
print("R1 =", R1)
print("R2 =", R2)
print()
print("r =", r)
print("s =", s)
print("R =", R)


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
    t = 0

    while d % 2 == 0:
        d //= 2
        t += 1

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

        for _ in range(t - 1):

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
# EXACT CARRY STATE
# ------------------------------------------------------------------------

def carry_state(p, q, r, s):

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

    Q = (p * q) // (r * s)

    return {
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "E": E,
        "Q": Q,
        "K": Q - E,
    }


# ------------------------------------------------------------------------
# DIRECT CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    tested = 0

    for p in range(2, math.isqrt(n) + 1):

        tested += 1

        if n % p == 0:
            return p, n // p, tested

    return None, None, tested


# ------------------------------------------------------------------------
# NEW DERIVATION
#
# We have:
#
#     n = (rk+a)(sl+b)
#
# therefore:
#
#     Q = floor(n/(rs))
#
# and exactly:
#
#     Q =
#       k*l
#       + floor(k*b/s)
#       + floor(l*a/r)
#       + c3
#
# Let
#
#     u = floor(k*b/s)
#     v = floor(l*a/r)
#
# Then:
#
#     Q = k*l + u + v + c3.
#
# Because u and v are directly tied to k and l:
#
#     u*s <= k*b < (u+1)*s
#
#     v*r <= l*a < (v+1)*r
#
# The experiment enumerates ONLY a,b,c3 and transforms these
# floor conditions into intervals for k and l.
#
# The question:
#
# Can the exact quotient equation determine k,l without scanning
# the sqrt(n)-scale factor interval?
# ------------------------------------------------------------------------


def solve_state(n, a, b, c3):

    Q = n // R

    solutions = []

    # ------------------------------------------------------------
    # Degenerate residues.
    #
    # If a=0 or b=0 the corresponding floor term vanishes and the
    # equations need a separate treatment.
    # ------------------------------------------------------------

    # For our random semiprimes these cases are unlikely, but handle
    # them explicitly.
    if a == 0 or b == 0:
        return solutions

    # ------------------------------------------------------------
    # Since:
    #
    #     u = floor(k*b/s)
    #
    # we have approximately:
    #
    #     u ~= k*b/s.
    #
    # Likewise:
    #
    #     v ~= l*a/r.
    #
    # We eliminate v using:
    #
    #     l*a = v*r + alpha
    #
    # for alpha in [0,r-1].
    #
    # This gives:
    #
    #     l = (v*r + alpha)/a.
    #
    # Instead of scanning l, scan the small alpha state.
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # We do the symmetric construction using:
    #
    #     beta = (k*b) mod s
    #
    #     alpha = (l*a) mod r
    #
    # where:
    #
    #     0 <= beta < s
    #     0 <= alpha < r.
    #
    # Thus:
    #
    #     k*b = u*s + beta
    #
    #     l*a = v*r + alpha
    #
    # and:
    #
    #     u = (k*b-beta)/s
    #
    #     v = (l*a-alpha)/r.
    #
    # Substitute into:
    #
    #     Q = kl + u + v + c3.
    #
    # Multiply by a*s*r etc. and solve for k,l.
    # ------------------------------------------------------------

    for alpha in range(r):

        # l*a == alpha mod r
        #
        # We don't directly know l. We will derive l from v.
        if math.gcd(a, r) != 1:
            continue

        for beta in range(s):

            if math.gcd(b, s) != 1:
                continue

            # ----------------------------------------------------
            # We have:
            #
            #   k*b = u*s + beta
            #
            # so
            #
            #   u = (k*b-beta)/s
            #
            # and:
            #
            #   l*a = v*r + alpha
            #
            #   v = (l*a-alpha)/r
            #
            # Substitute:
            #
            # Q-c3 =
            #     k*l
            #   + (k*b-beta)/s
            #   + (l*a-alpha)/r
            #
            # Multiply by r*s:
            #
            # r*s*(Q-c3) =
            #     r*s*k*l
            #   + r*k*b
            #   + s*l*a
            #   - r*beta
            #   - s*alpha
            #
            # Rearrange:
            #
            # (r*k + a)*(s*l + b)
            #     = r*s*Q + r*beta + s*alpha
            #
            # This is essentially the original multiplication,
            # but now alpha,beta are bounded by r,s.
            # ----------------------------------------------------

            target = (
                R * (Q - c3)
                + r * beta
                + s * alpha
            )

            if target <= 0:
                continue

            # ----------------------------------------------------
            # We now need:
            #
            #     (r*k+a)(s*l+b) = target
            #
            # Let:
            #
            #     P = r*k+a
            #     Q2 = s*l+b
            #
            # so P*Q2=target.
            #
            # We can use sqrt(target) factor enumeration.
            #
            # BUT the whole point is to see whether the restricted
            # residue conditions dramatically reduce the divisor
            # search.
            #
            # Therefore we enumerate divisors only of target.
            # ----------------------------------------------------

            limit = math.isqrt(target)

            d = 1

            while d <= limit:

                if target % d == 0:

                    p = d
                    q = target // d

                    if (
                        p % r == a
                        and q % s == b
                    ):

                        # Recover k,l.
                        k = (p - a) // r
                        ell = (q - b) // s

                        if k >= 0 and ell >= 0:

                            # Verify exact carry state.
                            c = carry_state(
                                p,
                                q,
                                r,
                                s,
                            )

                            if c["c3"] == c3:
                                if (
                                    c["alpha"]
                                    if "alpha" in c
                                    else (
                                        (ell * a) % r
                                    )
                                ) == alpha:

                                    if (
                                        (k * b) % s
                                    ) == beta:

                                        if (
                                            p * q
                                            == n
                                        ):
                                            solutions.append(
                                                (
                                                    p,
                                                    q,
                                                    a,
                                                    b,
                                                    alpha,
                                                    beta,
                                                    c3,
                                                )
                                            )

                d += 1

    return solutions


# ------------------------------------------------------------------------
# IMPORTANT OPTIMIZED VERSION
#
# The naive state solver above still factors "target".
#
# Instead of actually doing that, use the fact that target is very close
# to R*n/R = n and differs only by a bounded carry/residue correction.
#
# We therefore use gcd(target,n):
#
#     target = p*q-like quantity
#
# and see whether the correction causes a useful gcd leak.
# ------------------------------------------------------------------------

def gcd_state_test(n, a, b, alpha, beta, c3):

    Q = n // R

    target = (
        R * (Q - c3)
        + r * beta
        + s * alpha
    )

    g = math.gcd(
        n,
        target
    )

    return target, g


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(131)

TEST_BITS = [
    30,
    36,
    42,
]


for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, true_p, true_q = random_semiprime(bits)

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    true_c = carry_state(
        true_p,
        true_q,
        r,
        s,
    )

    print()
    print("TRUE STATE")

    for key in (
        "k",
        "ell",
        "a",
        "b",
        "c1",
        "c2",
        "c3",
        "E",
        "Q",
        "K",
    ):
        print(
            key,
            "=",
            true_c[key]
        )

    alpha_true = (
        true_c["ell"]
        * true_c["a"]
    ) % r

    beta_true = (
        true_c["k"]
        * true_c["b"]
    ) % s

    print(
        "alpha =",
        alpha_true
    )

    print(
        "beta  =",
        beta_true
    )

    # ------------------------------------------------------------
    # Direct control
    # ------------------------------------------------------------

    print()
    print("DIRECT SQRT(n) SCAN")

    t0 = time.perf_counter()

    pd, qd, tested = direct_factor(n)

    t1 = time.perf_counter()

    direct_time = t1 - t0

    print(
        "p tested =",
        tested
    )

    print(
        "runtime =",
        f"{direct_time:.6f}",
        "s"
    )

    # ------------------------------------------------------------
    # GCD LEAK TEST
    #
    # This is the actual new test.
    #
    # Enumerate only:
    #
    #   alpha < r
    #   beta < s
    #   c3
    #
    # and compute gcd(n,target).
    #
    # If a nontrivial factor repeatedly appears, we have discovered
    # a useful algebraic leakage mechanism.
    # ------------------------------------------------------------

    print()
    print("GCD CARRY-STATE TEST")

    t0 = time.perf_counter()

    tests = 0
    nontrivial = []
    true_state_hit = False

    # c3 is extremely small.  For these decompositions it is at most 1,
    # but use 0..2 conservatively.
    for alpha in range(r):

        for beta in range(s):

            for c3 in range(3):

                tests += 1

                target, g = gcd_state_test(
                    n,
                    true_c["a"],
                    true_c["b"],
                    alpha,
                    beta,
                    c3,
                )

                if 1 < g < n:

                    nontrivial.append(
                        (
                            alpha,
                            beta,
                            c3,
                            target,
                            g,
                        )
                    )

                if (
                    alpha == alpha_true
                    and beta == beta_true
                    and c3 == true_c["c3"]
                ):
                    true_state_hit = True

    t1 = time.perf_counter()

    gcd_time = t1 - t0

    print(
        "states tested =",
        tests
    )

    print(
        "nontrivial gcd hits =",
        len(nontrivial)
    )

    print(
        "true state enumerated =",
        true_state_hit
    )

    print(
        "runtime =",
        f"{gcd_time:.6f}",
        "s"
    )

    if nontrivial:

        print()
        print("NONTRIVIAL GCD HITS")

        for item in nontrivial[:20]:

            alpha, beta, c3, target, g = item

            print(
                "alpha =", alpha,
                "beta =", beta,
                "c3 =", c3,
                "gcd =", g
            )

            if (
                true_p % g == 0
                or true_q % g == 0
            ):
                print(
                    ">>> TRUE FACTOR DIVISOR <<<"
                )

    else:

        print()
        print(
            "NO NONTRIVIAL GCD LEAK"
        )

    # ------------------------------------------------------------
    # FULL STATE SOLVER ONLY FOR 30-BIT CONTROL
    #
    # This intentionally avoids the expensive solver on larger
    # instances.  It tells us whether the state equation is even
    # correctly reconstructing the factorization.
    # ------------------------------------------------------------

    if bits == 30:

        print()
        print("FULL STATE RECOVERY CONTROL")

        t0 = time.perf_counter()

        found = []

        # Enumerate actual residue states.
        for a in range(r):

            for b in range(s):

                for c3 in range(3):

                    sols = solve_state(
                        n,
                        a,
                        b,
                        c3,
                    )

                    found.extend(
                        sols
                    )

        t1 = time.perf_counter()

        print(
            "solutions =",
            len(found)
        )

        print(
            "runtime =",
            f"{t1 - t0:.6f}",
            "s"
        )

        exact = any(
            {x[0], x[1]}
            == {
                true_p,
                true_q,
            }
            for x in found
        )

        print(
            "TRUE FACTORIZATION FOUND =",
            exact
        )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 131
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 131")
print("=" * 72)
