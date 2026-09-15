#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 130
# Carry-conditioned residue prediction
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 130")
print("Carry-conditioned residue prediction")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
#
# Use several independent radix pairs.
# ------------------------------------------------------------------------

R1 = [17, 43, 47]
R2 = [19, 37, 53]

print()
print("R1 =", R1)
print("R2 =", R2)


# ------------------------------------------------------------------------
# PRIME TEST
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
# RANDOM PRIME
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
# CARRY STATE
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

    K = Q - E

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
        "K": K,
    }


# ------------------------------------------------------------------------
# DIRECT FACTOR CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    limit = math.isqrt(n)

    tested = 0

    for p in range(2, limit + 1):

        tested += 1

        if n % p == 0:
            return p, n // p, tested

    return None, None, tested


# ------------------------------------------------------------------------
# STATE TABLE FOR ONE RADIX PAIR
#
# For fixed (r,s) and residue pair (a,b), derive:
#
#     p = r*k + a
#     q = s*l + b
#
# and therefore
#
#     n = rs*k*l + rb*k + sa*l + ab.
#
# Rearranged:
#
#     Q = k*l + floor(
#             (rb*k + sa*l + ab)/(rs)
#         )
#
# We do NOT know k,l.
#
# But the carry split gives:
#
#     E = c1+c2+c3
#
# and
#
#     K = Q-E = k*l.
#
# The experiment asks whether consistency across several radix pairs
# can eliminate residue pairs without enumerating p.
# ------------------------------------------------------------------------

def build_observed_cells(n):

    cells = {}

    for r in R1:

        for s in R2:

            cells[(r, s)] = n // (r * s)

    return cells


# ------------------------------------------------------------------------
# RESIDUE PREDICTION
#
# For a chosen (a,b), derive congruence information about k,l from
#
#     r*k+a = p
#     s*l+b = q
#
# and
#
#     pq = n.
#
# Reducing n=(rk+a)(sl+b) modulo r gives:
#
#     n == a*(s*l+b) mod r
#
# hence, when gcd(a,r)=1:
#
#     q == n*a^{-1} mod r.
#
# Likewise:
#
#     p == n*b^{-1} mod s.
#
# This gives cross-radix consistency without scanning p.
# ------------------------------------------------------------------------

def residue_constraints_for_state(n, r, s, a, b):

    out = {
        "valid": True,
        "q_mod_r": None,
        "p_mod_s": None,
    }

    # ------------------------------------------------------------
    # From:
    #
    #     n = p*q
    #     p == a mod r
    #
    # if gcd(a,r)=1:
    #
    #     q == n*a^-1 mod r.
    #
    # ------------------------------------------------------------

    if math.gcd(a, r) == 1:

        out["q_mod_r"] = (
            n * pow(a, -1, r)
        ) % r

    # ------------------------------------------------------------
    # Likewise:
    #
    #     p == n*b^-1 mod s.
    # ------------------------------------------------------------

    if math.gcd(b, s) == 1:

        out["p_mod_s"] = (
            n * pow(b, -1, s)
        ) % s

    return out


# ------------------------------------------------------------------------
# BUILD A COMPLETE RESIDUE SIGNATURE
#
# Signature consists of the cross constraints generated by every radix
# pair.
# ------------------------------------------------------------------------

def residue_signature(n, residue_vector_p):

    signature = []

    for i, r in enumerate(R1):

        a = residue_vector_p[i]

        if math.gcd(a, r) != 1:

            signature.append(
                (r, a, None)
            )

            continue

        for s in R2:

            qmodr = (
                n * pow(a, -1, r)
            ) % r

            signature.append(
                (
                    r,
                    a,
                    qmodr
                )
            )

    return tuple(signature)


# ------------------------------------------------------------------------
# CANDIDATE RESIDUE VECTORS
#
# Enumerate residue vectors only.
#
# This is deliberately small and independent of sqrt(n).
# ------------------------------------------------------------------------

def enumerate_residue_vectors(n):

    states = []

    total = 0
    surviving = 0

    for a0 in range(R1[0]):

        for a1 in range(R1[1]):

            for a2 in range(R1[2]):

                total += 1

                a = (a0, a1, a2)

                # ------------------------------------------------
                # Each residue must actually be possible for a
                # divisor p of n.
                #
                # We derive q modulo every r_i.
                # ------------------------------------------------

                valid = True
                q_constraints = []

                for ai, r in zip(a, R1):

                    if ai == 0:

                        # p divisible by r.
                        #
                        # Then n divisible by r, which is extremely
                        # restrictive for our random semiprimes.
                        if n % r != 0:

                            valid = False
                            break

                        q_constraints.append(
                            (r, 0)
                        )

                    else:

                        if math.gcd(ai, r) != 1:

                            valid = False
                            break

                        qmod = (
                            n
                            * pow(ai, -1, r)
                        ) % r

                        q_constraints.append(
                            (r, qmod)
                        )

                if not valid:
                    continue

                surviving += 1

                states.append(
                    {
                        "a": a,
                        "q_constraints": tuple(
                            q_constraints
                        ),
                    }
                )

    return {
        "total": total,
        "surviving": surviving,
        "states": states,
    }


# ------------------------------------------------------------------------
# CROSS-RADIX CONSISTENCY
#
# A candidate p residue vector predicts q modulo every R1 radix.
#
# But q also has a residue vector modulo R2.
#
# We enumerate b-vectors and ask whether the predicted q residues and
# actual q CRT residue are mutually compatible.
#
# This is the key experiment:
#
# Can we determine (a0,a1,a2) without knowing p?
# ------------------------------------------------------------------------

def filter_p_residue_vectors(n):

    p_states = []

    total_a = 0
    total_ab_tests = 0

    for a0 in range(R1[0]):

        for a1 in range(R1[1]):

            for a2 in range(R1[2]):

                a = (
                    a0,
                    a1,
                    a2,
                )

                total_a += 1

                # ------------------------------------------------
                # Predict q modulo every r_i.
                # ------------------------------------------------

                q_from_r = []

                valid_a = True

                for ai, r in zip(a, R1):

                    if math.gcd(ai, r) != 1:

                        valid_a = False
                        break

                    qmod = (
                        n
                        * pow(ai, -1, r)
                    ) % r

                    q_from_r.append(
                        (r, qmod)
                    )

                if not valid_a:
                    continue

                # ------------------------------------------------
                # Enumerate q residues b_j.
                # ------------------------------------------------

                for b0 in range(R2[0]):

                    for b1 in range(R2[1]):

                        for b2 in range(R2[2]):

                            total_ab_tests += 1

                            b = (
                                b0,
                                b1,
                                b2,
                            )

                            # ------------------------------------------------
                            # b must be compatible with n modulo every
                            # R2 radix:
                            #
                            # p == n*b^-1 mod s.
                            #
                            # That predicted p residue must agree with
                            # the a-vector modulo s when the moduli are
                            # different but we can test p itself through
                            # CRT.
                            #
                            # For each s_j:
                            #
                            #     p == n*b_j^-1 mod s_j
                            #
                            # Build the resulting CRT class of p.
                            # ------------------------------------------------

                            valid_b = True

                            p_constraints = []

                            for bj, s in zip(
                                b,
                                R2
                            ):

                                if math.gcd(bj, s) != 1:

                                    valid_b = False
                                    break

                                pmod = (
                                    n
                                    * pow(bj, -1, s)
                                ) % s

                                p_constraints.append(
                                    (s, pmod)
                                )

                            if not valid_b:
                                continue

                            # ------------------------------------------------
                            # We now have:
                            #
                            # p mod R1_i = a_i
                            # p mod R2_j = derived pmod_j
                            #
                            # Combine these as a single consistency system.
                            #
                            # There is a contradiction if two moduli
                            # overlap nontrivially and the residues disagree.
                            # Here all radices are pairwise coprime, so we
                            # need both CRT classes to correspond to some p.
                            #
                            # Because everything is coprime, they always do.
                            #
                            # Therefore we go one step further:
                            #
                            # construct the complete p CRT class and q CRT
                            # class and check whether their product can equal
                            # n modulo the combined modulus.
                            # ------------------------------------------------

                            # p CRT:
                            # directly combine a-vector with p_constraints.
                            px = 0
                            pM = 1

                            for rr, aa in (
                                list(zip(R1, a))
                                + p_constraints
                            ):

                                inv = pow(
                                    pM,
                                    -1,
                                    rr
                                )

                                t = (
                                    (aa - px)
                                    * inv
                                ) % rr

                                px += pM * t
                                pM *= rr
                                px %= pM

                            # q CRT:
                            qx = 0
                            qM = 1

                            for rr, qq in (
                                q_from_r
                                + list(zip(R2, b))
                            ):

                                inv = pow(
                                    qM,
                                    -1,
                                    rr
                                )

                                t = (
                                    (qq - qx)
                                    * inv
                                ) % rr

                                qx += qM * t
                                qM *= rr
                                qx %= qM

                            # ------------------------------------------------
                            # Necessary condition:
                            #
                            #     px*qx == n (mod gcd(pM,qM))
                            #
                            # Since our moduli are pairwise coprime within
                            # each group and across groups here, gcd=1.
                            #
                            # Therefore this condition is automatically
                            # true.
                            #
                            # We retain the state so that Experiment 130
                            # can measure exactly how much this route buys.
                            # ------------------------------------------------

                            p_states.append(
                                {
                                    "a": a,
                                    "b": b,
                                    "p_class": px,
                                    "p_modulus": pM,
                                    "q_class": qx,
                                    "q_modulus": qM,
                                }
                            )

    return {
        "total_a": total_a,
        "total_ab_tests": total_ab_tests,
        "states": p_states,
    }


# ------------------------------------------------------------------------
# NEW TEST:
#
# Rather than falsely claiming the cross congruences are enough, compare
# their information content against the known factor.
#
# For each residue state we ask:
#
#     How many p < sqrt(n) fit the complete CRT class?
#
# We then determine whether the TRUE class is unique.
# ------------------------------------------------------------------------

def count_class_candidates(n, A, M):

    limit = math.isqrt(n)

    if A == 0:
        first = M
    else:
        first = A

    while first < 2:
        first += M

    if first > limit:
        return 0

    return (
        (limit - first)
        // M
    ) + 1


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(130)

TEST_BITS = [
    30,
    36,
    42,
    48,
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

    # ------------------------------------------------------------
    # True residue vector.
    # ------------------------------------------------------------

    true_a = [
        true_p % r
        for r in R1
    ]

    true_b = [
        true_q % s
        for s in R2
    ]

    print()
    print("TRUE p residues =", true_a)
    print("TRUE q residues =", true_b)

    # ------------------------------------------------------------
    # Direct control.
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    pd, qd, direct_tested = direct_factor(n)

    t1 = time.perf_counter()

    direct_time = t1 - t0

    print()
    print("DIRECT SQRT(n) SCAN")
    print("p tested =", direct_tested)
    print("runtime  =", f"{direct_time:.6f}", "s")

    # ------------------------------------------------------------
    # Build residue prediction.
    # ------------------------------------------------------------

    print()
    print("CARRY-CONDITIONED RESIDUE PREDICTION")

    t0 = time.perf_counter()

    result = enumerate_residue_vectors(n)

    t1 = time.perf_counter()

    prediction_time = t1 - t0

    print()
    print(
        "p residue vectors total =",
        result["total"]
    )

    print(
        "p residue vectors surviving gcd constraints =",
        result["surviving"]
    )

    print(
        "prediction runtime =",
        f"{prediction_time:.6f}",
        "s"
    )

    # ------------------------------------------------------------
    # Find the true residue vector inside the surviving states.
    # ------------------------------------------------------------

    true_survives = any(
        state["a"] == tuple(true_a)
        for state in result["states"]
    )

    print(
        "TRUE p RESIDUE VECTOR SURVIVES =",
        true_survives
    )

    # ------------------------------------------------------------
    # For each p residue vector, combine it with the inferred
    # q congruences and determine its candidate count below sqrt(n).
    # ------------------------------------------------------------

    print()
    print("CRT CLASS SIZE FOR p RESIDUE VECTORS")

    class_sizes = []

    for state in result["states"]:

        a = state["a"]

        A = 0
        M = 1

        for ai, r in zip(a, R1):

            inv = pow(
                M,
                -1,
                r
            )

            t = (
                (ai - A)
                * inv
            ) % r

            A += M * t
            M *= r
            A %= M

        count = count_class_candidates(
            n,
            A,
            M
        )

        class_sizes.append(
            (
                count,
                a,
                A,
                M,
            )
        )

    class_sizes.sort()

    print()
    print("minimum class size =", class_sizes[0][0])
    print("maximum class size =", class_sizes[-1][0])

    # ------------------------------------------------------------
    # TRUE CLASS.
    # ------------------------------------------------------------

    true_A = 0
    true_M = 1

    for ai, r in zip(true_a, R1):

        inv = pow(
            true_M,
            -1,
            r
        )

        t = (
            (ai - true_A)
            * inv
        ) % r

        true_A += true_M * t
        true_M *= r
        true_A %= true_M

    true_count = count_class_candidates(
        n,
        true_A,
        true_M
    )

    print()
    print("TRUE CRT CLASS")
    print("A =", true_A)
    print("M =", true_M)
    print("candidate count =", true_count)

    # ------------------------------------------------------------
    # Show smallest classes.
    # ------------------------------------------------------------

    print()
    print("SMALLEST CRT CLASSES")

    for count, a, A, M in class_sizes[:10]:

        exact = (
            tuple(true_a)
            == a
        )

        print(
            "a =", a,
            "A =", A,
            "M =", M,
            "count =", count,
            "TRUE =", exact
        )

    # ------------------------------------------------------------
    # FINAL INTERPRETATION.
    # ------------------------------------------------------------

    unique_classes = sum(
        1
        for count, _, _, _
        in class_sizes
        if count <= 1
    )

    true_unique = (
        true_count <= 1
    )

    print()
    print("SUMMARY")

    print(
        "candidate classes with <=1 p =",
        unique_classes
    )

    print(
        "TRUE class unique =",
        true_unique
    )

    print(
        "TRUE p recoverable from CRT class =",
        true_A == true_p
        if true_unique
        else (
            true_p % true_M
            == true_A
        )
    )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 130
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 130")
print("=" * 72)
